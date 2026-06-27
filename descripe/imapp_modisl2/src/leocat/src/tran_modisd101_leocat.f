	  subroutine tran_modisd101(coeff_dir_path,year,jday,
     +                           temp,wvmr,ozmr,theta,platform,
     +                           kban,jdet,
     +				     taut,iok)
c * MODIS band/detector 101-level fast transmittance routine
c .... version of 06.08.03

c    coeff_dir_path = directory that contains fast model coefficients
c        year = profile year
c        jday = profile day-of-year
c	 temp = temperature (Kelvin) profile
c	 wvmr = water-vapor mixing-ratio (g/kg) profile
c	 ozmr = ozone mixing-ratio (ppmv) profile
c	theta = local zenith angle
c	craft = TERRA, AQUA (either upper or lower case)
c	 kban = band number     (20...36)
c	 jdet = detector number (0...10) [ Product Order ]
c		 detector 0 is based on band-average response functions

c		 taut = total transmittance (see note below)
c		  iok = 0 if successful, 1 if I/O problem

c * NOTE: for kban = 26, return-arrays are filled with 1.0

c * PLOD/PFAAST regression model based on LBLRTM line-by-line transmittances.
c * Input temperatures, and water-vapor and ozone mixing ratios, must
c *	be defined at the pressure levels in array 'pstd'
c *    (see block data 'reference_atmosphere').
c * Units: temperature, deg-K; water vapor, g/kg; ozone, ppmv.
c * Logical units 31-35 are used for coefficient files.
c * Component taus are returned through common, product in 'taut'.
 
      parameter (nd=10,nk=5,nl=101,nm=nl-1,koff=19,nr=17)
      parameter (nxc= 4,ncc=nxc+1,lencc=ncc*nm,lenccb=lencc*4)
      parameter (nxd= 8,ncd=nxd+1,lencd=ncd*nm,lencdb=lencd*4)
      parameter (nxo= 9,nco=nxo+1,lenco=nco*nm,lencob=lenco*4)
      parameter (nxl= 2,ncl=nxl+1,lencl=ncl*nm,lenclb=lencl*4)
      parameter (nxs=11,ncs=nxs+1,lencs=ncs*nm,lencsb=lencs*4)
      parameter (ndt=nd+1,nrps=nr*ndt,nxw=nxl+nxs)
      parameter (slp=1.5/365.0)
      parameter (smag=3.0)
      parameter (pi=3.14159)
      parameter (soff=0.41)
      parameter (coff=337.5)

      common/stdatm/pstd(nl),tstd(nl),wstd(nl),ostd(nl)
      common/taudwo/taud(nl),tauw(nl),tauo(nl)
	  
      dimension temp(*),wvmr(*),ozmr(*),taut(*)
      dimension coefd(ncd,nm,0:nd,nr),coefo(nco,nm,0:nd,nr),
     +           coefl(ncl,nm,0:nd,nr),coefs(ncs,nm,0:nd,nr),
     +           coefc(ncc,nm,0:nd,nr)
      dimension bufs(lencs)
      dimension pavg(nm),tref(nm),wref(nm),oref(nm)
      dimension tavg(nm),wamt(nm),oamt(nm),secz(nm)
      dimension tauc(nl),tlas(nl),wlas(nl),olas(nl)
      dimension xdry(nxd,nm),xozo(nxo,nm),xwet(nxw,nm),xcon(nxc,nm)

      real*4 secant,z
      real*4 x,rco2,ratio,tau_test
      character*28 xfile
      character*256 coeff_dir_path
      character*256 path
      character*256 cfile(nk),dfile
      character*6 craft,cinit
      character*6 cbt,cba
      character*6 cst,csa
      character*3 comp(nk)
      character*3 cbe,cle
      integer*4 platform
      integer*4 lengcf(nk)
      integer*4 lengcx(nk)
      integer*4 iuc(nk)
      integer*2 year,jday,jj
      logical big_endian,newang,newatm

      data tlas/nl*0./,wlas/nl*0./,olas/nl*0./,zlas/-999./
      data xfile/'/modisdet.com.101.xxx_end.v3'/
      data cinit/'zzzzzz'/
      data cbt/'TERRA'/,cba/'AQUA'/
      data cst/'terra'/,csa/'aqua'/
      data comp/'dry','ozo','wts','wtl','wco'/
      data cbe/'big'/,cle/'lit'/
      data lengcf/lencdb,lencob,lencsb,lenclb,lenccb/
      data lengcx/lencd,lenco,lencs,lencl,lencc/

      save

      secant(z)=1./cos(0.01745329*z)

      iok=0

c     write(*,'(''begin 0'')')
      	
      if (platform .eq. 0) then
	craft = 'TERRA'	    
c       Shifted bands 34-36, LBL 11.3
        path = coeff_dir_path      
      elseif (platform .eq. 1) then
	craft = 'AQUA'
        path = coeff_dir_path
c       write(*,'(''In platform if block '', a80)') path
c		path = '/data/raf/fm_coeff_lbl-11.3/srfs/shifted_band34'   
      else
        write(*,'(''tran_modisd101- unknown spacecraft '',i2)') platform
      endif
c     write(*,'(''End platform if block '', a80)') path

      if(craft.ne.cinit) then
	  
        if(craft.eq.cbt.or.craft.eq.cst) then
	  ksat=1
	elseif(craft.eq.cba.or.craft.eq.csa) then
	  ksat=2
	else
          write(*,'(''tran_modisd101- unknown spacecraft '',a6)') craft
	  go to 200
	endif

c *     determine which coefficient files to use
	if(big_endian()) then
	  xfile(19:21)=cbe
	else
	  xfile(19:21)=cle
	endif
	   
c *     define and open the files
	iux=30
	do m=1,nk

	  iux=iux+1
	  xfile(11:13)=comp(m)

c RAF     This code block modified to enable file name manipulation
c         when passed in from a C routine.

          dfile = path
          j = len(dfile)
          kk = len(xfile)
          i = 1
          do while( (i .le. j) .and. (ichar(dfile(i:i)) .ne. 0) )
            i = i + 1
          enddo

c         dfile(i:i+3-1) = '/..'
c         i = i + 3
          dfile(i:i+kk-1) = xfile(1:kk)
          jj = i + kk
          do while( jj .le. j)
            dfile(jj:jj) = ' '
            jj = jj + 1
          enddo

	  lencf=lengcf(m)
	      
c         write(*,'(a25)') xfile
          write(*,'(''file open '', 6i8,2x,a120)') j,kk,i,m,iux,lencf,
     +                   dfile
          open(iux,file=dfile,recl=lencf,access='direct',
     +		    status='old',err=200)
	  iuc(m)=iux
	  cfile(m)=dfile
	      
	enddo
        write(*,'(''cfiles '', a120)') cfile(1),cfile(2),cfile(3),
     +          cfile(4),cfile(5)

c *     first read each files fill-record for band 26/det 0
c		and verify satellite number stored in word 1
c *     note: number of levels is in word 2, creation date (yyyyddd) is in word 3
	ikrec=nrps*(ksat-1)
	krecx=ikrec+7
	do k=1,nk
	  lencx=lengcx(k)
	  read(iuc(k),rec=krecx) (bufs(j),j=1,lencx)
	  nsat=bufs(1)
	  if(nsat.ne.ksat) then
	    dfile=cfile(k)
	    go to 100
	  endif
	enddo

c *     now read in the coefficients
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
	  
c     End "first-time-through-only" block.	   
      endif
	  
c     write(*,'(''Begin'')')
      
	  dt=0.
	  dw=0.
	  do=0.
	  do j=1,nl
	   dt=dt+abs(temp(j)-tlas(j))
	   tlas(j)=temp(j)
	   dw=dw+abs(wvmr(j)-wlas(j))
	   wlas(j)=wvmr(j)
	   do=do+abs(ozmr(j)-olas(j))
	   olas(j)=ozmr(j)
	   taud(j)=1.0
	   tauw(j)=1.0
	   tauc(j)=1.0
	   tauo(j)=1.0
	   taut(j)=1.0
	  enddo
	  datm=dt+dw+do
	  newatm=datm.ne.0.
	  
c	  write(*,'(''newatm '', l4)') newatm
	  
	  if(newatm) then
	   call conpir(pstd,temp,wvmr,ozmr,nl,1,pavg,tavg,wamt,oamt)
c	   write(*,'(''done with conpir'')')
c	   call flush(6)
	  endif

	  newang=theta.ne.zlas
	  if(newang) then
	   zsec=secant(theta)
	   do l=1,nm
	      secz(l)=zsec
	   enddo
	   zlas=theta
	  endif

c    write(*,'(''newang '',l4)') newang

	  if(newang.or.newatm) then
	   call calpir(tref,wref,oref,tavg,wamt,oamt,pavg,secz,
     +			 nm,nxd,nxw,nxo,nxc,xdry,xwet,xozo,xcon)
c      write(*,'(''done with calpir'')')
c      call flush(6)
	  endif

	  if(kban.eq.26) then
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
c *   dry
	  call taudoc(ncd,nxd,nm,coefd(1,1,j,k),xdry,taud)

c         Adjust dry tau for changes in CO2 concentration from model (R. Frey)
          x = (year - 1980) * 365.25 + jday
          rco2 = (slp*x - smag*sin(2*pi*(x/365.25 + soff))) + coff
          ratio=rco2/380.0
          do jj=1,nl
            tau_test = taud(jj)
            if(taud(jj).gt.0.0 .and. taud(jj).lt.1.0) then
              taud(jj)=taud(jj)**ratio

c             write(*,'(''dry tau adj. '',2i4,5f12.4)') kban,jj,x,rco2,
c    *                 ratio,tau_test,taud(jj)

            endif
          enddo

c *   ozo
	  call taudoc(nco,nxo,nm,coefo(1,1,j,k),xozo,tauo)
c *   wet
	  call tauwtr(ncs,ncl,nxs,nxl,nxw,nm,coefs(1,1,j,k),
     +		     coefl(1,1,j,k),xwet,tauw)
	  call taudoc(ncc,nxc,nm,coefc(1,1,j,k),xcon,tauc)
	  do j=1,nl
	   tauw(j)=tauw(j)*tauc(j)
	  enddo
c * total
	  do j=1,nl
	   taut(j)=taud(j)*tauo(j)*tauw(j)
	  enddo
	  return
100	  write(*,'(''In tran_modisd101 ... requested data for '',
     +           ''satellite '',i1/'' but read data for '',
     +           ''satellite '',i1,'' from file '',a80)')
     +       ksat,nsat,dfile
200	  iok=1
	  return
	  end

c * PLOD/PFAAST subprograms for infrared transmittance
c	at 101-level AIRS/MODIS pressure coordinate
c .... version of 18.03.03

contents:

c     block data reference_atmosphere
c     subroutine calpir
c     subroutine conpir
c     subroutine gphite
c     subroutine taudoc
c     subroutine tauwtr

c *********************

	  block data reference_atmosphere
c *   Reference Atmosphere is 1976 U.S. Standard
	  parameter (nl=101)
	  common/stdatm/pstd(nl),tstd(nl),wstd(nl),ostd(nl)

      data pstd    / 0.0050,    0.0161,    0.0384,    0.0769,    0.1370,
     +    0.2244,    0.3454,    0.5064,    0.7140,    0.9753,    1.2972,
     +    1.6872,    2.1526,    2.7009,    3.3398,    4.0770,    4.9204,
     +    5.8776,    6.9567,    8.1655,    9.5119,   11.0038,   12.6492,
     +   14.4559,   16.4318,   18.5847,   20.9224,   23.4526,   26.1829,
     +   29.1210,   32.2744,   35.6505,   39.2566,   43.1001,   47.1882,
     +   51.5278,   56.1260,   60.9895,   66.1253,   71.5398,   77.2396,
     +   83.2310,   89.5204,   96.1138,  103.0172,  110.2366,  117.7775,
     +  125.6456,  133.8462,  142.3848,  151.2664,  160.4959,  170.0784,
     +  180.0183,  190.3203,  200.9887,  212.0277,  223.4415,  235.2338,
     +  247.4085,  259.9691,  272.9191,  286.2617,  300.0000,  314.1369,
     +  328.6753,  343.6176,  358.9665,  374.7241,  390.8926,  407.4738,
     +  424.4698,  441.8819,  459.7118,  477.9607,  496.6298,  515.7200,
     +  535.2322,  555.1669,  575.5248,  596.3062,  617.5112,  639.1398,
     +  661.1920,  683.6673,  706.5654,  729.8857,  753.6275,  777.7897,
     +  802.3714,  827.3713,  852.7880,  878.6201,  904.8659,  931.5236,
     +  958.5911,  986.0666, 1013.9476, 1042.2319, 1070.9170, 1100.0000/

      data tstd               / 190.19, 203.65, 215.30, 226.87, 237.83, 
     +  247.50, 256.03, 263.48, 267.09, 270.37, 266.42, 261.56, 256.40, 
     +  251.69, 247.32, 243.27, 239.56, 236.07, 232.76, 230.67, 228.71, 
     +  227.35, 226.29, 225.28, 224.41, 223.61, 222.85, 222.12, 221.42, 
     +  220.73, 220.07, 219.44, 218.82, 218.23, 217.65, 217.18, 216.91, 
     +  216.70, 216.70, 216.70, 216.70, 216.70, 216.70, 216.70, 216.70, 
     +  216.70, 216.70, 216.70, 216.70, 216.70, 216.70, 216.70, 216.71, 
     +  216.71, 216.72, 216.81, 217.80, 218.77, 219.72, 220.66, 222.51, 
     +  224.57, 226.59, 228.58, 230.61, 232.61, 234.57, 236.53, 238.48, 
     +  240.40, 242.31, 244.21, 246.09, 247.94, 249.78, 251.62, 253.45, 
     +  255.26, 257.04, 258.80, 260.55, 262.28, 264.02, 265.73, 267.42, 
     +  269.09, 270.77, 272.43, 274.06, 275.70, 277.32, 278.92, 280.51, 
     +  282.08, 283.64, 285.20, 286.74, 288.25, 289.75, 291.22, 292.68/ 

      data wstd               /  0.001,  0.001,  0.002,  0.003,  0.003, 
     +   0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003, 
     +   0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003, 
     +   0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003, 
     +   0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003, 
     +   0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003,  0.003, 
     +   0.003,  0.003,  0.004,  0.004,  0.005,  0.005,  0.007,  0.009, 
     +   0.011,  0.012,  0.014,  0.020,  0.025,  0.030,  0.035,  0.047, 
     +   0.061,  0.075,  0.089,  0.126,  0.162,  0.197,  0.235,  0.273, 
     +   0.310,  0.356,  0.410,  0.471,  0.535,  0.601,  0.684,  0.784, 
     +   0.886,  0.987,  1.094,  1.225,  1.353,  1.519,  1.686,  1.852, 
     +   2.036,  2.267,  2.496,  2.721,  2.947,  3.170,  3.391,  3.621, 
     +   3.848,  4.084,  4.333,  4.579,  4.822,  5.061,  5.296,  5.528/ 

      data ostd               /0.47330,0.27695,0.28678,0.51816,0.83229, 
     + 1.18466,1.69647,2.16633,3.00338,3.76287,4.75054,5.61330,6.33914, 
     + 7.03675,7.50525,7.75612,7.81607,7.69626,7.56605,7.28440,7.01002, 
     + 6.72722,6.44629,6.17714,5.92914,5.69481,5.47387,5.26813,5.01252, 
     + 4.68941,4.35141,4.01425,3.68771,3.37116,3.06407,2.77294,2.50321, 
     + 2.24098,1.98592,1.74840,1.54451,1.34582,1.17824,1.02513,0.89358, 
     + 0.78844,0.69683,0.62654,0.55781,0.50380,0.45515,0.42037,0.38632, 
     + 0.35297,0.32029,0.28832,0.25756,0.22739,0.19780,0.16877,0.14901, 
     + 0.13190,0.11511,0.09861,0.08818,0.07793,0.06786,0.06146,0.05768, 
     + 0.05396,0.05071,0.04803,0.04548,0.04301,0.04081,0.03983,0.03883, 
     + 0.03783,0.03685,0.03588,0.03491,0.03395,0.03368,0.03349,0.03331, 
     + 0.03313,0.03292,0.03271,0.03251,0.03190,0.03126,0.03062,0.02990, 
     + 0.02918,0.02850,0.02785,0.02721,0.02658,0.02596,0.02579,0.02579/ 

	  end
ccccccccc
      subroutine calpir(t_avg_ref,  amt_wet_ref, amt_ozo_ref,
     +                        t_avg,      amt_wet,     amt_ozo,
     +                        p_avg,      sec_theta,   n_layers,
     +                        n_dry_pred, n_wet_pred,  n_ozo_pred,
     +                        n_con_pred,
     +                        pred_dry,   pred_wet,    pred_ozo,
     +                        pred_con)
c ... version of 18.03.03

c  PURPOSE:

c    Routine to calculate the predictors for the dry (temperature), 
c      wet and ozone components of a fast transmittance model for a
c      scanning satellite based instrument.

c  REFERENCES:

c    AIRS FTC package science notes and software, S. Hannon and L. Strow,
c      Uni. of Maryland, Baltimore County (UMBC)

c  CREATED:

c    19-Sep-1996 HMW

c  ARGUMENTS:

c      Input
c    -----------
c     t_avg_ref  - REAL*4 reference layer average temperature array (K)

c    amt_wet_ref - REAL*4 reference water vapour amount array (k.mol)/cm^2

c    amt_ozo_ref - REAL*4 reference ozone amount array (k.mol)/cm^2

c      t_avg     - REAL*4 layer average temperature array (K)

c     amt_wet    - REAL*4 water vapour amount array (k.mol)/cm^2

c     amt_ozo    - REAL*4 ozone amount array (k.mol)/cm^2

c      p_avg     - REAL*4 layer average pressure array (mb)

c    sec_theta   - REAL*4 secant of the zenith angle array

c     n_layers   - INT*4 Number of atmospheric layers

c    n_dry_pred  - INT*4 number of dry (temperature) predictors

c    n_wet_pred  - INT*4 number of water vapour predictors

c    n_ozo_pred  - INT*4 number of ozone predictors

c    n_con_pred  - INT*4 number of water vapour continuum predictors

c      Output
c    -----------
c     pred_dry   - REAL*4 dry gas (temperature) predictor matrix

c     pred_wet   - REAL*4 water vapour predictor matrix

c     pred_ozo   - REAL*4 ozone predictor matrix

c     pred_con   - REAL*4 water vapour continuum predictor matrix

c  COMMENTS:

c    Levels or Layers?
c    -----------------
c      Profile data is input at a number of *LAYERS*.

c    Layer Numbering pt. A
c    ---------------------
c      Layer 1   => Atmosphere between LEVELs 1 & 2
c      Layer 2   => Atmosphere between LEVELs 2 & 3
c                        .
c                        .
c                        .
c      Layer L-1 => Atmosphere between LEVELs L-1 & L

c    Layer Numbering pt. B
c    ---------------------
c      For the HIS instrument, Layer 1 is at the top of the atmosphere
c        and Layer L-1 is at the surface.    

c    Layer Numbering pt. C
c    ---------------------
c      In this routine the number of *LAYERS* is passed in the argument
c        list, _not_ the number of LEVELS.  This was done to improve
c        the readability of this code, i.e. loop from 1->L(ayers) 
c        rather than from 1->L(evels)-1.

c=======================================================================

c-----------------------------------------------------------------------
c                 Turn off implicit type declaration
c-----------------------------------------------------------------------

       implicit none

c------------------------------------------------------------------------
c                             Arguments
c------------------------------------------------------------------------

c -- Input

      integer*4 n_layers,
     +          n_dry_pred, n_wet_pred, n_ozo_pred, n_con_pred

      real*4    t_avg_ref(*), amt_wet_ref(*), amt_ozo_ref(*),
     +          t_avg(*),     amt_wet(*),     amt_ozo(*),
     +          p_avg(*),     sec_theta(*)

c -- Output

      real*4    pred_dry(n_dry_pred, *),
     +          pred_wet(n_wet_pred, *),
     +          pred_ozo(n_ozo_pred, *),
     +          pred_con(n_con_pred, *)

c------------------------------------------------------------------------
c                           Local variables
c------------------------------------------------------------------------

c -- Parameters

      integer*4 max_layers
      parameter ( max_layers = 100 )

      integer*4 max_dry_pred, max_wet_pred, max_ozo_pred, max_con_pred
      parameter ( max_dry_pred = 8,
     +            max_wet_pred = 13,
     +            max_ozo_pred = 9,
     +            max_con_pred = 4 )

c -- Scalars

      integer*4 l

c -- Arrays

c     ....Pressure
      real*4    p_dp(max_layers),
     +          p_norm(max_layers)

c     ....Temperature
      real*4    delta_t(max_layers),
     +          t_ratio(max_layers),
     +          pw_t_ratio(max_layers)      ! Pressure weighted

c     ....Water vapour
      real*4    wet_ratio(max_layers),
     +          pw_wet(max_layers),         ! Pressure weighted
     +          pw_wet_ref(max_layers),     ! Pressure weighted
     +          pw_wet_ratio(max_layers)    ! Pressure weighted

c     ....Ozone
      real*4    ozo_ratio(max_layers), 
     +          pw_ozo_ratio(max_layers),   ! Pressure weighted
     +          pow_t_ratio(max_layers)     ! Pressure/ozone weighted

c************************************************************************
c                         ** Executable code **
c************************************************************************

c------------------------------------------------------------------------
c                   -- Check that n_layers is o.k. --
c------------------------------------------------------------------------

      if( n_layers .gt. max_layers )then
        write(*,'(/10x,''*** calpir : n_layers > max_layers'')')
        stop
      end if 

c------------------------------------------------------------------------
c         -- Check that numbers of predictors is consistent --
c------------------------------------------------------------------------

c     ---------------------------------
c     # of dry (temperature) predictors
c     ---------------------------------

      if( n_dry_pred .ne. max_dry_pred )then
        write(*,'(/10x,''*** calpir : invalid n_dry_pred'')')
        stop
      end if 

c     ----------------------------
c     # of water vapour predictors
c     ----------------------------

      if( n_wet_pred .ne. max_wet_pred )then
        write(*,'(/10x,''*** calpir : invalid n_wet_pred'')')
        stop
      end if 

c     ---------------------
c     # of ozone predictors
c     ---------------------

      if( n_ozo_pred .ne. max_ozo_pred )then
        write(*,'(/10x,''*** calpir : invalid n_ozo_pred'')')
        stop
      end if 

c     --------------------------------------
c     # of water vapour continuum predictors
c     --------------------------------------

      if( n_con_pred .ne. max_con_pred )then
        write(*,'(/10x,''*** calpir : invalid n_con_pred'')')
        stop
      end if 

c------------------------------------------------------------------------
c         -- Calculate ratios, offsets, etc, for top layer --
c------------------------------------------------------------------------

c     ------------------
c     Pressure variables
c     ------------------

      p_dp(1)   = p_avg(1) * ( p_avg(2) - p_avg(1) )
      p_norm(1) = 0.0

c     ---------------------
c     Temperature variables
c     ---------------------

      delta_t(1)    = t_avg(1) - t_avg_ref(1)
      t_ratio(1)    = t_avg(1) / t_avg_ref(1)
      pw_t_ratio(1) = 0.0

c     ----------------
c     Amount variables
c     ----------------

c     ....Water vapour
 
      wet_ratio(1)    = amt_wet(1) / amt_wet_ref(1)
      pw_wet(1)       = p_dp(1) * amt_wet(1)
      pw_wet_ref(1)   = p_dp(1) * amt_wet_ref(1)
      pw_wet_ratio(1) = wet_ratio(1)

c     ....Ozone

      ozo_ratio(1)    = amt_ozo(1) / amt_ozo_ref(1)
      pw_ozo_ratio(1) = 0.0
      pow_t_ratio(1)  = 0.0

c------------------------------------------------------------------------
c         -- Calculate ratios, offsets, etc, for all layers --
c------------------------------------------------------------------------

      do l = 2, n_layers

c       ------------------
c       Pressure variables
c       ------------------

        p_dp(l) = p_avg(l) * ( p_avg(l) - p_avg(l-1) )
        p_norm(l) = p_norm(l-1) + p_dp(l)

c       ---------------------
c       Temperature variables
c       ---------------------

        delta_t(l)    = t_avg(l) - t_avg_ref(l)
        t_ratio(l)    = t_avg(l) / t_avg_ref(l)
        pw_t_ratio(l) = pw_t_ratio(l-1) + ( p_dp(l) * t_ratio(l-1) )

c       ----------------
c       Amount variables
c       ----------------

c       ..Water vapour

        wet_ratio(l)  = amt_wet(l) / amt_wet_ref(l)
        pw_wet(l)     = pw_wet(l-1) + ( p_dp(l) * amt_wet(l) )
        pw_wet_ref(l) = pw_wet_ref(l-1) + ( p_dp(l) * amt_wet_ref(l) )
        
c       ..Ozone

        ozo_ratio(l)    = amt_ozo(l) / amt_ozo_ref(l)
        pw_ozo_ratio(l) = pw_ozo_ratio(l-1) +
     +                      ( p_dp(l) * ozo_ratio(l-1) )
        pow_t_ratio(l)  = pow_t_ratio(l-1) +
     +                      ( p_dp(l) * ozo_ratio(l-1) * delta_t(l-1) )

      end do

c------------------------------------------------------------------------
c              -- Scale the pressure dependent variables --
c------------------------------------------------------------------------

      do l = 2, n_layers

        pw_t_ratio(l)   = pw_t_ratio(l) / p_norm(l)
        pw_wet_ratio(l) = pw_wet(l) / pw_wet_ref(l)
        pw_ozo_ratio(l) = pw_ozo_ratio(l) / p_norm(l)
        pow_t_ratio(l)  = pow_t_ratio(l) / p_norm(l)
 
      end do

c------------------------------------------------------------------------
c                     -- Load up predictor arrays --
c------------------------------------------------------------------------

      do l = 1, n_layers

c       ----------------------
c       Temperature predictors
c       ----------------------

        pred_dry(1,l) = sec_theta(l)
        pred_dry(2,l) = sec_theta(l) * sec_theta(l)
        pred_dry(3,l) = sec_theta(l) * t_ratio(l)
        pred_dry(4,l) = pred_dry(3,l) * t_ratio(l)
        pred_dry(5,l) = t_ratio(l)
        pred_dry(6,l) = t_ratio(l) * t_ratio(l)
        pred_dry(7,l) = sec_theta(l) * pw_t_ratio(l)
        pred_dry(8,l) = pred_dry(7,l) / t_ratio(l) 

c       -----------------------
c       Water vapour predictors
c       -----------------------

        pred_wet(1,l)  = sec_theta(l) * wet_ratio(l)
        pred_wet(2,l)  = sqrt( pred_wet(1,l) )
        pred_wet(3,l)  = pred_wet(1,l) * delta_t(l)
        pred_wet(4,l)  = pred_wet(1,l) * pred_wet(1,l)
        pred_wet(5,l)  = abs( delta_t(l) ) * delta_t(l) * pred_wet(1,l)
        pred_wet(6,l)  = pred_wet(1,l) * pred_wet(4,l)
        pred_wet(7,l)  = sec_theta(l) * pw_wet_ratio(l)
        pred_wet(8,l)  = pred_wet(2,l) * delta_t(l)
        pred_wet(9,l)  = sqrt( pred_wet(2,l) )
        pred_wet(10,l) = pred_wet(7,l) * pred_wet(7,l)
        pred_wet(11,l) = sqrt( pred_wet(7,l) )
        pred_wet(12,l) = pred_wet(1,l)
c +++ old
c       pred_wet(13,l) = pred_wet(2,l)
c +++ new
        pred_wet(13,l) = pred_wet(1,l) / pred_wet(10,l)

c       ----------------
c       Ozone predictors
c       ----------------

        pred_ozo(1,l) = sec_theta(l) * ozo_ratio(l)
        pred_ozo(2,l) = sqrt( pred_ozo(1,l) )
        pred_ozo(3,l) = pred_ozo(1,l) * delta_t(l)
        pred_ozo(4,l) = pred_ozo(1,l) * pred_ozo(1,l)
        pred_ozo(5,l) = pred_ozo(2,l) * delta_t(l)
        pred_ozo(6,l) = sec_theta(l) * pw_ozo_ratio(l)
        pred_ozo(7,l) = sqrt( pred_ozo(6,l) ) * pred_ozo(1,l)
        pred_ozo(8,l) = pred_ozo(1,l) * pred_wet(1,l)
        pred_ozo(9,l) = sec_theta(l) * pow_t_ratio(l) * pred_ozo(1,l)

c       ---------------------------------
c       Water vapour continuum predictors
c       ---------------------------------

        pred_con(1,l) = sec_theta(l) * wet_ratio(l) /
     *	 	    ( t_ratio(l) * t_ratio(l) )  
        pred_con(2,l) = pred_con(1,l) * pred_con(1,l) / sec_theta(l)
        pred_con(3,l) = sec_theta(l) * wet_ratio(l) / t_ratio(l) 
        pred_con(4,l) = pred_con(3,l) * wet_ratio(l)

      end do
         
      return
      end
ccccccccc
      subroutine conpir( p, t, w, o, n_levels, i_dir,
     +                         p_avg, t_avg, w_amt, o_amt)
c ... version of 19.09.96

c  PURPOSE:

c    Function to convert atmospheric water vapour (g/kg) and ozone (ppmv)
c      profiles specified at n_levels layer BOUNDARIES to n_levels-1
c      integrated layer amounts of units (k.moles)/cm^2.  The average
c      LAYER pressure and temperature are also returned.

c  REFERENCES:

c    AIRS LAYERS package science notes, S. Hannon and L. Strow, Uni. of
c      Maryland, Baltimore County (UMBC)

c  CREATED:

c    19-Sep-1996 HMW

c  ARGUMENTS:

c     Input
c    --------
c       p     - REAL*4 pressure array (mb)

c       t     - REAL*4 temperature profile array (K)

c       w     - REAL*4 water vapour profile array (g/kg)

c       o     - REAL*4 ozone profile array (ppmv)

c    n_levels - INT*4 number of elements used in passed arrays

c     i_dir   - INT*4 direction of increasing layer number

c                 i_dir = +1, Level(1) == p(top)         } satellite/AC
c                             Level(n_levels) == p(sfc)  }    case

c                 i_dir = -1, Level(1) == p(sfc)         } ground-based
c                             Level(n_levels) == p(top)  }    case

c     Output
c    --------
c     p_avg   - REAL*4 average LAYER pressure array (mb)

c     t_avg   - REAL*4 average LAYER temperature (K)

c     w_amt   - REAL*4 integrated LAYER water vapour amount array (k.moles)/cm^2

c     o_amt   - REAL*4 integrated LAYER ozone amount array (k.moles)/cm^2

c  ROUTINES:

c    Subroutines:
c    ------------
c      gphite      - calculates geopotential height given profile data.

c    Functions:
c    ----------
c      NONE

c  COMMENTS:

c    Levels or Layers?
c    -----------------
c      Profile data is input at a number of *LEVELS*.  Number densitites
c        are calculated for *LAYERS* that are bounded by these levels.
c        So, for L levels there are L-1 layers.

c    Layer Numbering
c    ---------------
c      Layer 1   => Atmosphere between LEVELs 1 & 2
c      Layer 2   => Atmosphere between LEVELs 2 & 3
c                        .
c                        .
c                        .
c      Layer L-1 => Atmosphere between LEVELs L-1 & L

c=======================================================================

c-----------------------------------------------------------------------
c              -- Prevent implicit typing of variables --
c-----------------------------------------------------------------------

      implicit none

c-----------------------------------------------------------------------
c                           -- Arguments --
c-----------------------------------------------------------------------

c -- Arrays

      real*4    p(*), t(*), w(*), o(*), 
     +          p_avg(*), t_avg(*), w_amt(*), o_amt(*)

c -- Scalars

      integer*4 n_levels, i_dir

c-----------------------------------------------------------------------
c                         -- Local variables --
c-----------------------------------------------------------------------

c -- Parameters

      integer*4 max_levels
      parameter ( max_levels = 101 )         ! Maximum number of layers

      real*4    r_equator, r_polar, r_avg
      parameter ( r_equator = 6.378388e+06, ! Earth radius at equator
     +            r_polar   = 6.356911e+06, ! Earth radius at pole
     +            r_avg     = 0.5*(r_equator+r_polar) )

      real*4    g_sfc
      parameter ( g_sfc = 9.80665 )         ! Gravity at surface

      real*4    rho_ref
      parameter ( rho_ref = 1.2027e-12 )    ! Reference air "density"

      real*4    mw_dryair, mw_h2o, mw_o3
      parameter ( mw_dryair = 28.97,        ! Molec. wgt. of dry air (g/mol)
     +            mw_h2o    = 18.0152,      ! Molec. wgt. of water
     +            mw_o3     = 47.9982 )     ! Molec. wgt. of ozone

      real*4    R_gas, R_air
      parameter ( R_gas = 8.3143,           ! Ideal gas constant (J/mole/K)
     +            R_air = 0.9975*R_gas )    ! Gas constant for air (worst case) 

c -- Scalars

      integer*4 l, l_start, l_end, l_indx

      real*4    rho1, rho2, p1, p2, w1, w2, o1, o2, z1, z2,
     +          c_avg, g_avg, z_avg, w_avg, o_avg,
     +          dz, dp, r_hgt, A, B

c -- Arrays

      real*4    z(max_levels),              ! Pressure heights (m)
     +          g(max_levels),              ! Acc. due to gravity (m/s/s)
     +          mw_air(max_levels),         ! Molec. wgt. of air (g/mol)
     +          rho_air(max_levels),        ! air mass density (kg.mol)/m^3
     +          c(max_levels),              ! (kg.mol.K)/(N.m)
     +          w_ppmv(max_levels)          ! h2o LEVEL amount (ppmv)

c***********************************************************************
c                         ** Executable code **
c***********************************************************************

c-----------------------------------------------------------------------
c           -- Calculate initial values of pressure heights --
c-----------------------------------------------------------------------

      call gphite( p, t, w, 0.0, n_levels, i_dir, z)

c-----------------------------------------------------------------------
c      -- Set loop bounds for direction sensitive calculations --
c      -- so loop iterates from surface to the top             --
c-----------------------------------------------------------------------

      if( i_dir .gt. 0 )then

c       --------------------
c       Data stored top down
c       --------------------

        l_start = n_levels
        l_end   = 1

      else

c       ---------------------
c       Data stored bottom up
c       ---------------------

        l_start = 1
        l_end   = n_levels

      end if

c-----------------------------------------------------------------------
c          -- Air molecular mass and density, and gravity --
c          -- as a function of LEVEL                      --
c-----------------------------------------------------------------------

c     -----------------------
c     Loop from bottom to top
c     -----------------------

      do l = l_start, l_end, -1*i_dir

c       ---------------------------------
c       Convert water vapour g/kg -> ppmv
c       ---------------------------------

        w_ppmv(l) = 1.0e+03 * w(l) * mw_dryair / mw_h2o

c       -----------------------------------------
c       Calculate molecular weight of air (g/mol)
c       ----------------------------------------

        mw_air(l) = ( ( 1.0 - (w_ppmv(l)/1.0e+6) ) * mw_dryair ) +
     +              ( ( w_ppmv(l)/1.0e+06 ) * mw_h2o )

c       ----------------
c       Air mass density
c       ----------------

        c(l) = 0.001 * mw_air(l) / R_air    ! 0.001 factor for g -> kg
        rho_air(l) = c(l) * p(l) / t(l)

c       -------
c       Gravity
c       -------

        r_hgt = r_avg + z(l)                !  m
        g(l) = g_sfc -                      !  m/s^2
     +         g_sfc*( 1.0 - ( (r_avg*r_avg)/(r_hgt*r_hgt) ) )

      end do
 
c-----------------------------------------------------------------------
c                        -- LAYER quantities --
c-----------------------------------------------------------------------

c     -----------------------
c     Loop from bottom to top
c     -----------------------

      do l = l_start, l_end+i_dir, -1*i_dir

c       -------------------------------------------------------
c       Determine output array index.  This is done so that the
c       output data is always ordered from 1 -> L-1 regardless
c       of the orientation of the input data.  This is true by
c       default only for the bottom-up case.  For the top down
c       case no correction would give output layers from 2 -> L
c       -------------------------------------------------------

        if( i_dir .gt. 0 )then

          l_indx = l - 1

        else

          l_indx = l

        end if

c       ---------------------------------------
c       Assign current layer boundary densities
c       ---------------------------------------
 
        rho1 = rho_air(l)
        rho2 = rho_air(l-i_dir)
 
c       ---------
c       Average c
c       ---------

        c_avg = ( (rho1*c(l)) + (rho2*c(l-i_dir)) ) / ( rho1 + rho2 )

c       ---------
c       Average t
c       ---------

        t_avg(l_indx) = 
     +          ( (rho1*t(l)) + (rho2*t(l-i_dir)) ) / ( rho1 + rho2 )

c       ---------
c       Average p
c       ---------

        p1 = p(l)
        p2 = p(l-i_dir)

        z1 = z(l)
        z2 = z(l-i_dir)

        dp = p2 - p1

        A = log(p2/p1) / (z2-z1)
        B = p1 / exp(A*z1)

        p_avg(l_indx) = dp / log(p2/p1)

c       ------------------------------------------------
c       LAYER thickness (rather long-winded as it is not
c       assumed the layers are thin) in m. Includes
c       correction for altitude/gravity.
c       ------------------------------------------------

c       ...Initial values
        g_avg = g(l)
        dz = -1.0 * dp * t_avg(l_indx) / ( g_avg*c_avg*p_avg(l_indx) )

c       ...Calculate z_avg
        z_avg = z(l) + ( 0.5*dz )

c       ...Calculate new g_avg
        r_hgt = r_avg + z_avg 
        g_avg = g_sfc - g_sfc*( 1.0 - ( (r_avg*r_avg)/(r_hgt*r_hgt) ) )

c       ...Calculate new dz
        dz = -1.0 * dp * t_avg(l_indx) / ( g_avg*c_avg*p_avg(l_indx) )

c       ----------------------------------------
c       Calculate LAYER amounts for water vapour
c       ----------------------------------------

        w1 = w_ppmv(l)
        w2 = w_ppmv(l-i_dir)

        w_avg =  ( (rho1*w1) + (rho2*w2) ) / ( rho1+rho2 )

        w_amt(l_indx) =
     +       rho_ref * w_avg * dz * p_avg(l_indx) / t_avg(l_indx)

c       ---------------------------------
c       Calculate LAYER amounts for ozone
c       ---------------------------------

        o1 = o(l)
        o2 = o(l-i_dir)

        o_avg =  ( (rho1*o1) + (rho2*o2) ) / ( rho1+rho2 )

        o_amt(l_indx) = 
     +       rho_ref * o_avg * dz * p_avg(l_indx) / t_avg(l_indx)

      end do

      return
      end
ccccccccc
      subroutine gphite( p, t, w, z_sfc, n_levels, i_dir, z)
c .... version of 18.05.00

c PURPOSE:

c  Routine to compute geopotential height given the atmospheric state.
c    Includes virtual temperature adjustment.

c CREATED:

c  19-Sep-1996 Received from Hal Woolf, recoded by Paul van Delst
c  18-May-2000 Logic error related to z_sfc corrected by Hal Woolf

c  ARGUMENTS:

c     Input
c    --------
c       p     - REAL*4 pressure array (mb)

c       t     - REAL*4 temperature profile array (K)

c       w     - REAL*4 water vapour profile array (g/kg)

c     z_sfc   - REAL*4 surface height (m).  0.0 if not known.

c    n_levels - INT*4 number of elements used in passed arrays

c     i_dir   - INT*4 direction of increasing layer number

c                 i_dir = +1, Level(1) == p(top)         } satellite/AC
c                             Level(n_levels) == p(sfc)  }    case

c                 i_dir = -1, Level(1) == p(sfc)         } ground-based
c                             Level(n_levels) == p(top)  }    case

c     Output
c    --------
c       z     - REAL*4 pressure level height array (m)

c COMMENTS:

c   Dimension of height array may not not be the same as that of the
c     input profile data.

c=======================================================================

c-----------------------------------------------------------------------
c              -- Prevent implicit typing of variables --
c-----------------------------------------------------------------------

      implicit none

c-----------------------------------------------------------------------
c                           -- Arguments --
c-----------------------------------------------------------------------

c -- Arrays

      real*4    p(*), t(*), w(*), 
     +          z(*)

c -- Scalars

      integer*4 n_levels, i_dir

      real*4    z_sfc

c-----------------------------------------------------------------------
c                         -- Local variables --
c-----------------------------------------------------------------------

c -- Parameters

      real*4    rog, fac
      parameter ( rog = 29.2898, 
     +            fac = 0.5 * rog )

c -- Scalars

      integer*4 i_start, i_end, l

      real*4    v_lower, v_upper, algp_lower, algp_upper, hgt

c***********************************************************************
c                         ** Executable code **
c***********************************************************************

c-----------------------------------------------------------------------
c  -- Calculate virtual temperature adjustment and exponential       --
c  -- pressure height for level above surface.  Also set integration --
c  -- loop bounds                                                    --
c-----------------------------------------------------------------------

      if( i_dir .gt. 0 )then

c       --------------------
c       Data stored top down
c       --------------------

        v_lower = t(n_levels) * ( 1.0 + ( 0.00061 * w(n_levels) ) )

        algp_lower = alog( p(n_levels) )

        i_start = n_levels-1
        i_end   = 1

      else

c       ---------------------
c       Data stored bottom up
c       ---------------------

        v_lower = t(1) * ( 1.0 + ( 0.00061 * w(1) ) )

        algp_lower = alog( p(1) )

        i_start = 2
        i_end   = n_levels

      end if

c-----------------------------------------------------------------------
c                     -- Assign surface height --
c-----------------------------------------------------------------------

      hgt = z_sfc

c .. Following added 18 May 2000 ... previously, z(n_levels) for downward
c       (usual) case was not defined!

      if(i_dir.gt.0) then
         z(n_levels) = z_sfc
      else
         z(1) = z_sfc
      endif

c .. End of addition

c-----------------------------------------------------------------------
c             -- Loop over layers always from sfc -> top --
c-----------------------------------------------------------------------

      do l = i_start, i_end, -1*i_dir

c       ----------------------------------------------------
c       Apply virtual temperature adjustment for upper level
c       ----------------------------------------------------

        v_upper = t(l)
        if( p(l) .ge. 300.0 )
     +    v_upper = v_upper * ( 1.0 + ( 0.00061 * w(l) ) )

c       ----------------------------------------------------- 
c       Calculate exponential pressure height for upper layer
c       ----------------------------------------------------- 

        algp_upper = alog( p(l) )

c       ----------------
c       Calculate height
c       ----------------

        hgt = hgt + ( fac*(v_upper+v_lower)*(algp_lower-algp_upper) )

c       -------------------------------
c       Overwrite values for next layer
c       -------------------------------

        v_lower = v_upper
        algp_lower = algp_upper

c       ---------------------------------------------
c       Store heights in same direction as other data
c       ---------------------------------------------

        z(l) = hgt

      end do

      return
      end
ccccccccc
	subroutine taudoc(nc,nx,ny,cc,xx,tau)
c * Strow-Woolf model ... for dry, ozo(ne), and wco (water-vapor continuum)
c .... version of 05.09.02
	dimension cc(nc,ny),xx(nx,ny),tau(*)
	data trap/-999.99/

	last=ny
	taulyr=1.
	taulev=1.
	tau(1)=taulev
	do 110 j=1,ny
	if(taulev.eq.0.) go to 110
	if(j.gt.last) go to 100
	yy=cc(nc,j)
	if(yy.eq.trap) then
	   last=j
	   go to 100
	endif
	do i=1,nx
	   yy=yy+cc(i,j)*xx(i,j)
	enddo
	yy=amax1(yy,0.)
	if(yy.gt.0.) taulyr=exp(-yy)
100	taulev=taulev*taulyr
110	tau(j+1)=taulev
	return
	end
ccccccccc
	subroutine tauwtr(ncs,ncl,nxs,nxl,nxw,ny,ccs,ccl,xx,tau)
c * Strow-Woolf model ... for 'wet' (water-vapor other than continuum)
c .... version of 05.09.02
	dimension ccs(ncs,ny),ccl(ncl,ny),xx(nxw,ny),tau(*)

	odsum=0.
	taulyr=1.
	taulev=1.
	tau(1)=taulev
	do 110 j=1,ny
	if(odsum.lt.5.) then
	   yy=ccs(ncs,j)
	   do i=1,nxs
	      yy=yy+ccs(i,j)*xx(i,j)
	   enddo
	else
	   yy=ccl(ncl,j)
	   do i=1,nxl
	      yy=yy+ccl(i,j)*xx(i+11,j)
	   enddo
	endif
	yy=amax1(yy,0.)
	odsum=odsum+yy
	if(yy.gt.0.) taulyr=exp(-yy)
	taulev=taulev*taulyr
110	tau(j+1)=taulev
	return
	end
