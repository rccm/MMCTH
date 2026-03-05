module pfaast

	use FASCODE_routines

	private
	
	public :: MODIS_fascode
	

      real ::	coefd(Fncd,Fnm,0:Fnd,Fnr),&
				coefo(Fnco,Fnm,0:Fnd,Fnr), &
                coefl(Fncl,Fnm,0:Fnd,Fnr),&
				coefs(Fncs,Fnm,0:Fnd,Fnr), &
                coefc(Fncc,Fnm,0:Fnd,Fnr)
	
contains


	subroutine MODIS_fascode(coeff_dir_path, year, jday, temp, wvmr, ozmr, theta, ang_2way, platform, &
							kban, jdet, taut, taut_2way, newang, newatm, new_2way, do_2way, iok, xxx, yyy ) 
	
! * MODIS band/detector 101-level fast transmittance routine
! .... version of 06.08.03

!    coeff_dir_path = directory that contains fast model coefficients
!        year = profile year
!        jday = profile day-of-year
!	 temp = temperature (Kelvin) profile
!	 wvmr = water-vapor mixing-ratio (g/kg) profile
!	 ozmr = ozone mixing-ratio (ppmv) profile
!	theta = local zenith angle
!	craft = TERRA, AQUA (either upper or lower case)
!	 kban = band number     (20...36)
!	 jdet = detector number (0...10) [ Product Order ]
!		 detector 0 is based on band-average response functions

!		 taut = total transmittance (see note below)
!		  iok = 0 if successful, 1 if I/O problem

! * NOTE: for kban = 26, return-arrays are filled with 1.0

! * PLOD/PFAAST regression model based on LBLRTM line-by-line transmittances.
! * Input temperatures, and water-vapor and ozone mixing ratios, must
! *	be defined at the pressure levels in array 'pstd'
! *    (see block data 'reference_atmosphere').
! * Units: temperature, deg-K; water vapor, g/kg; ozone, ppmv.
! * Logical units 31-35 are used for coefficient files.
! * Component taus are returned through common, product in 'taut'.
		character(*), intent(in) :: coeff_dir_path
		integer, intent(in) :: year, jday
      real, intent(inout) :: temp(*),wvmr(*),ozmr(*),taut(*), taut_2way(*)
		integer, intent(in) :: platform
		real, intent(in) :: theta, ang_2way
		integer, intent(in) :: kban, jdet, xxx, yyy
		integer, intent(inout) :: iok
		logical, intent(in) :: newang, newatm, new_2way, do_2way
		
      integer, parameter :: nd=10,nk=5,nl=101,nm=nl-1,koff=19,nr=17
	  integer, parameter :: nxc= 4,ncc=nxc+1,lencc=ncc*nm,lenccb=lencc*4
     integer, parameter :: nxd= 8,ncd=nxd+1,lencd=ncd*nm,lencdb=lencd*4
      integer, parameter :: nxo= 9,nco=nxo+1,lenco=nco*nm,lencob=lenco*4
      integer, parameter :: nxl= 2,ncl=nxl+1,lencl=ncl*nm,lenclb=lencl*4
      integer, parameter :: nxs=11,ncs=nxs+1,lencs=ncs*nm,lencsb=lencs*4
      integer, parameter :: ndt=nd+1,nrps=nr*ndt,nxw=nxl+nxs
      real, parameter :: slp=1.5/365.0
      real, parameter :: smag=3.0
      real, parameter :: pi=3.14159
	  real, parameter :: soff=0.41
      real, parameter :: coff=337.5

!      common/stdatm/pstd(nl),tstd(nl),wstd(nl),ostd(nl)
!     common/taudwo/taud(nl),tauw(nl),tauo(nl)
	  real :: taud(nl), tauw(nl), tauo(nl)
	  real :: taud_2way(nl), tauw_2way(nl), tauo_2way(nl)
	  
       real :: bufs(lencs)
!      real :: pavg(nm),tref(nm),wref(nm),oref(nm)
!      real ::  tavg(nm),wamt(nm),oamt(nm),secz(nm)
      real :: tauc(nl),tauc_2way(nl), tlas(nl),wlas(nl),olas(nl)

      real*4 x,rco2,ratio,tau_test
      character*28 xfile
      character*256 path
      character*256 cfile(nk),dfile
      character*6 craft
      character*6 cbt,cba
      character*6 cst,csa
      character*3 comp(nk)
      character*3 cbe,cle
      integer*4 lengcf(nk)
      integer*4 lengcx(nk)
      integer*4 iuc(nk)
      integer*2 jj
      logical big_endian

	  real ::  zlas, zlas2way
	  
	  integer :: ksat, iux, m, j, kk, ikrec, krec, krecx, k, lencx, l, i
	  integer :: lencf
	  integer*2 :: nsat
	  real :: dt, dw, fdo, datm, zsec, zsec_2way

!      data cinit/'zzzzzz'/
	

		tlas = nl*0.
		wlas = nl*0.
		olas = nl*0.
		zlas = -999.
		zlas2way = -999.
		xfile = '/modisdet.com.101.xxx_end.v3'
		cbe = 'big'
		cle = 'lit'
		
		cbt = 'TERRA'
		cba = 'AQUA'
		cst = 'terra'
		csa = 'aqua'

		comp = (/'dry','ozo','wts','wtl','wco'/)
		lengcf = (/lencdb,lencob,lencsb,lenclb,lenccb/)
		lengcx = (/lencd,lenco,lencs,lencl,lencc/)

		iok = 0
		
		if (platform == 0) then 
!       Shifted bands 34-36, LBL 11.3
			craft = 'TERRA'
			path = coeff_dir_path
		else if (platform == 1) then 
			craft = 'AQUA'
			path = coeff_dir_path
		else
			write(*,'(''tran_modisd101- unknown spacecraft '',i2)') platform
			iok = 1
			return 
		endif
				
		if (craft /= cinit) then 
		
			if (craft == "TERRA") then	
				ksat = 1
			else if (craft == "AQUA") then 
				ksat = 2
			else
				write(*,'(''tran_modisd101- unknown spacecraft '',a6)') craft
				iok = 1
				return
			endif
			
! select big- or little-endian libraries.
!			if (big_endian()) then 
				xfile(19:21) = cbe
!			else
!				xfile(19:21) = cle
!			endif
			
! *     define and open the files
			iux = 30
			do m=1, nk
			
				iux = iux + 1
				xfile(11:13) = comp(m)
! RAF     This code block modified to enable file name manipulation
!         when passed in from a C routine.

				dfile = path
				j = len(dfile)
				kk = len(xfile)
				i = 1
				do while( (i <= j) .and. (ichar(dfile(i:i)) .ne. 0) )
					i = i + 1
				enddo
				
				dfile(i:i+kk-1) = xfile(1:kk)
				jj = i + kk
				do while( jj <= j)
					dfile(jj:jj) = ' '
					jj = jj + 1
				enddo

				lencf=lengcf(m)
		
				open(iux,file=dfile,recl=lencf,access='direct', status='old', convert='big_endian')
				iuc(m)=iux
				cfile(m)=dfile
	      
			enddo
! *     first read each files fill-record for band 26/det 0
!		and verify satellite number stored in word 1
! *     note: number of levels is in word 2, creation date (yyyyddd) is in word 3
			ikrec=nrps*(ksat-1)
			krecx=ikrec+7
			do k=1,nk
				lencx=lengcx(k)
				read(iuc(k),rec=krecx) (bufs(j),j=1,lencx)
				nsat=bufs(1)
				if(nsat /= ksat) then
					dfile=cfile(k)
					write(*,'(''In tran_modisd101 ... requested data for '', &
& ''satellite '',i1/'' but read data for '', ''satellite '',i1,'' from file '',a80)') ksat,nsat,dfile
					iok = 1
					return
				endif
			enddo

! *     now read in the coefficients
			krec=ikrec
			do l=0,nd
				do k=1,nr
					krec=krec+1
					read(iuc(1),rec=krec) ((coefd(i,j,l,k),i=1,ncd),j=1,nm)
					read(iuc(2),rec=krec) ((coefo(i,j,l,k),i=1,nco),j=1,nm)
					read(iuc(3),rec=krec) ((coefs(i,j,l,k),i=1,ncs),j=1,nm)
					read(iuc(4),rec=krec) ((coefl(i,j,l,k),i=1,ncl),j=1,nm)
					read(iuc(5),rec=krec) ((coefc(i,j,l,k),i=1,ncc),j=1,nm)
				enddo
			enddo
			do k=1,nk
				close(iuc(k))
			enddo
		
			call conpir(pstd,tstd,wstd,ostd,nl,1,pavg,tref,wref,oref)

			cinit=craft
			iok=0
	  
!     End "first-time-through-only" block.	   
      endif
		

	  if(newatm) then
	   call conpir(pstd,temp,wvmr,ozmr,nl,1,pavg,tavg,wamt,oamt)
	  endif

! we have to enable doing two-way angles on an as-needed basis

	if(newang) then
		zsec = secant(theta)
		secz = zsec
	endif

	if (do_2way .and. new_2way) then 
		zsec_2way = secant(ang_2way)
		secz_2way = zsec_2way
	endif

	if(newang .or. newatm ) then
	   call calpir(tref,wref,oref,tavg,wamt,oamt,pavg,secz, &
     			 nm,nxd,nxw,nxo,nxc,xdry,xwet,xozo,xcon)
	endif

	if (do_2way .and. (new_2way .or. newatm)) then
		   call calpir(tref,wref,oref,tavg,wamt,oamt,pavg,secz_2way, &
     			 nm,nxd,nxw,nxo,nxc,xdry_2way,xwet_2way,xozo_2way,xcon_2way)
	endif


	if(kban == 26) then
	   do j=1,nl
	      taud(j)=1.0
	      tauo(j)=1.0
	      tauw(j)=1.0
	      taut(j)=1.0
	   enddo
	   return
	 endif

	  j=jdet
	  k=kban-koff
! *   dry
	  call taudoc(ncd,nxd,nm,coefd(:,:,j,k),xdry,taud)

!		do j=1, nl
!			print*, j, taud(j)
!		end do

!         Adjust dry tau for changes in CO2 concentration from model (R. Frey)

#ifdef CT_CODE
          x = (year - 1980) * 365.25 + jday
          rco2 = (slp*x - smag*sin(2*pi*(x/365.25 + soff))) + coff
          ratio=rco2/380.0
          do jj=1,nl
            tau_test = taud(jj)
            if(taud(jj) > 0.0 .and. taud(jj) < 1.0) then
              taud(jj)=taud(jj)**ratio
            endif
          enddo

! *   ozo
	  call taudoc(nco,nxo,nm,coefo(:,:,j,k),xozo,tauo)
#else
! we set the ozone to 0. if we're not doing CT, so transmittance due to ozone is 1.0
	   tauo = 1.0
#endif

! *   wet
	  call tauwtr(ncs,ncl,nxs,nxl,nxw,nm,coefs(:,:,j,k),coefl(:,:,j,k),xwet,tauw)
	  
	  call taudoc(ncc,nxc,nm,coefc(:,:,j,k),xcon,tauc)

	  do jj=1,nl
	   tauw(jj)=tauw(jj)*tauc(jj)
	  enddo
! * total
	  do jj=1,nl
	   taut(jj)=taud(jj)*tauo(jj)*tauw(jj)
	  enddo
	
	  if (do_2way) then 
	  
		j=jdet
	  
		call taudoc(ncd,nxd,nm,coefd(:,:,j,k),xdry_2way,taud_2way)
		do jj=1,nl
			tau_test = taud_2way(jj)
            if(taud_2way(jj) > 0.0 .and. taud_2way(jj) < 1.0) then
              taud_2way(jj)=taud_2way(jj)**ratio

            endif
		enddo
#ifdef CT_CODE 
		call taudoc(nco,nxo,nm,coefo(:,:,j,k),xozo_2way,tauo_2way)
#else
		tauo_2way = 1.0
#endif
		call tauwtr(ncs,ncl,nxs,nxl,nxw,nm,coefs(:,:,j,k),coefl(:,:,j,k),xwet_2way,tauw_2way)		
		call taudoc(ncc,nxc,nm,coefc(:,:,j,k),xcon_2way,tauc_2way)


		do jj=1,nl
			tauw_2way(jj)=tauw_2way(jj)*tauc_2way(jj)
		enddo
! * total
		do jj=1,nl
			taut_2way(jj)=taud_2way(jj)*tauo_2way(jj)*tauw_2way(jj)
		enddo
	  	  
	  
	  endif
	
	
	
	end subroutine MODIS_fascode	





end module pfaast
