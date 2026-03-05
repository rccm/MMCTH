      subroutine modis(fh, QAFlags, QAflags_O, iDTaot, dims)
!-----------------------------------------------------------------------
! !F90
!
! !DESCRIPTION:
!
!    The modis subroutine drives the processing of MODIS data through
!    the Deep Blue algorithm.  modis sequentially steps through the
!    intermediate data file, passing each datapoint first through
!    a reflectivity preprocessor, and then through the Deep Blue 
!    algorithm itself (find_v).  As the processed data is passed out 
!    from find_v, it is sorted and binned.  Once all the data has been 
!    processed, it is averaged and then written out to the MOD04 data 
!    product.
!
! !INPUT PARAMETERS:
!
!    Type       Name             Description
!    ====       ====             ===========
!    INTEGER*4  fh               file handle for the intermediate 
!                                 data file
!    BYTE       QAFlags          3-D array of bit flags from the 
!                                 MOD04 product
!    INTEGER*4  dims             array of SDS dimensions from the 
!                                 MOD04 product
!
! !OUTPUT PARAMETERS:  none
!
! !REVISION HISTORY:
!
!    Initial Version by Jeremy Warner   12/01/2006
!    Updated by Clare Salustro          05/02/2008
!
! !TEAM-UNIQUE HEADER:
!
!    This software is developed by the Deep Blue Science Team
!    for the National Aeronautics and Space Administration,
!    Goddard Space Flight Center, under contract NAS5-02041.
!
! !REFERENCES AND CREDITS
!
! !DESIGN NOTES:
!
!
!   Externals:
!
!     MODIS_W_GENERIC            (MODIS_39500.f)
!
!   Functions:
!
!     extract_data
!     read_data
!     total
!     find_v
!     write_db
!
! !END
!-----------------------------------------------------------------------
      use modis_surface, only: get_LER412,                    &
                               get_LER470,                    &
                               get_LER650,                    & 
                               get_LER865,                    &
                               terrain_flag_new
                               
      use dbdt_utils, only:   load_dbdt_region_table,         &
                              unload_dbdt_region_table,       &
                              create_dbdt_aot550
                              
      use core_arrays, only:  latitude,                       &
                              longitude,                      &
                              solar_zenith_angle,             &
                              sensor_zenith_angle
      
      implicit none
      INCLUDE 'PGS_MODIS_39500.f'
      INCLUDE 'PGS_SMF.f'
      INCLUDE 'PGS_PC.f'
      INCLUDE 'PGS_PC_9.f'
      INCLUDE 'PGS_IO.f'
      INCLUDE 'PGS_IO_1.f'

      integer fh
      INTEGER      No_byte,Fmax,Lmax, No_byte_O
      PARAMETER    (No_byte=6,Fmax=150,Lmax=550, No_byte_O=5)
      logical end_file, error_file
      integer(kind=4) :: len, dims(3)
      integer(kind=4) :: i, j, k, l, m
      integer(kind=4) :: LandSeaFlag
      integer(kind=4) :: idim, jdim
      integer(kind=4) :: iidx, jidx, maxiidx, maxjidx
      integer(kind=4) :: indsol,indscn,iofset
      integer(kind=4) :: intlpt(20)
      integer(kind=4) :: ierr, add_db_mod04
      integer(kind=4) :: isno,icat
      integer(kind=4) :: r412_count(Fmax,Lmax)
      integer(kind=4) :: t650_count(Fmax,Lmax)
      integer(kind=2) :: naot550_avg(Fmax,Lmax), nae_avg(Fmax,Lmax)
      integer(kind=2) :: naot_avg(Fmax,Lmax,3), nssa_avg(Fmax,Lmax,3)
      integer(kind=2) :: iaot_avg(Fmax,Lmax,3), issa_avg(Fmax,Lmax,3)
      integer(kind=2) :: nsr_avg(Fmax,Lmax,3), nref_avg(Fmax,Lmax,3)
      integer(kind=2) :: isr_avg(Fmax,Lmax,3), iref_avg(Fmax,Lmax,3)
      integer(kind=2) :: iaot550_avg(Fmax,Lmax), iae_avg(Fmax,Lmax)
      integer(kind=2) :: iaot550_best(Fmax,Lmax), iest_uncert(Fmax,Lmax)
      integer(kind=2) :: isd550(Fmax,Lmax)
      integer(kind=2) :: icld_frac(Fmax,Lmax)
      integer(kind=1) :: QAFlags(No_byte,Fmax,Lmax), QAFlags_O(No_byte_O,Fmax,Lmax)
      integer(kind=4) :: flags(4)
      
! New for combined SDS - 07 may 2010
      integer         :: nflag,n1,n2,n3,n4,n5,n6,n7,n8
      integer(kind=1) :: Quality_Dt(Fmax,Lmax), Quality_Dt_O(Fmax,Lmax)
      integer(kind=2) :: iDTaot(Fmax,Lmax), iDTDBaot(Fmax,Lmax), iDTDBflag(Fmax,Lmax)
      integer(kind=4) :: iDTDBqa(Fmax,Lmax)
! end new for combined SDS - 07 may 2010
      
      integer(kind=2) :: alg_cnt(Fmax,Lmax,2), alg_flag(Fmax,Lmax) 
      integer(kind=4) :: conf_flag(Fmax,Lmax), usefulness(Fmax,Lmax), aero_type(Fmax,Lmax), rcond(Fmax,Lmax)
      integer(kind=2) :: edgecnt(Fmax,Lmax)
      
      real(kind=4) :: realbuf, outbuf(20)
      real(kind=4) :: xdenom, xnvalm5(3), band26, Dstar1, tmpvg(6)
      real(kind=4) :: xzlog, xlog, densol, denscn, cthet0, ctheta, cofs, p1, pr
      real(kind=4) :: xxzlog(10), xxlog(8)
      real(kind=4) :: dum10, var
      real(kind=4) :: xlonlo, xlonhi
      real(kind=4) :: wlenth, fvalue
      real(kind=4) :: aot550_avg(Fmax,Lmax), ae_avg(Fmax,Lmax)
      real(kind=4) :: aot550_min(Fmax,Lmax), aot550_max(Fmax,Lmax)
      real(kind=4) :: aot550_list(Fmax,Lmax,100), sd550(Fmax,Lmax)
      real(kind=4) :: aot_avg(Fmax,Lmax,3), ssa_avg(Fmax,Lmax,3)
!      real(kind=4) :: aot_list(Fmax,Lmax,3,100)
      real(kind=4) :: sr_avg(Fmax,Lmax,3), ref_avg(Fmax,Lmax,3)
      real(kind=4) :: view(Fmax,Lmax), scat_ang(Fmax,Lmax)
      real(kind=4) :: sfc_typ(Fmax,Lmax)
      real(kind=4) :: alpha,  alpha_550, xxrat,  dd, aa550
      real(kind=4) :: alpha2, xxrat2, alpha2_550
      real(kind=4) :: refbuf(3), vazim, sazim
      real(kind=4) :: xsza(Fmax,Lmax), r412b_y40(Fmax,Lmax)
      real(kind=4) :: lat_sav(Fmax,Lmax), lon_sav(Fmax,Lmax)
      real(kind=4) :: sza_sav(Fmax,Lmax), vza_sav(Fmax,Lmax)
      real(kind=4) :: dstar_avg(Fmax, Lmax), amf(Fmax,Lmax)
      integer(kind=2) ::  nll_sav(Fmax,Lmax)
            
      real(kind=4) :: rbuff_11, rfbuf_11  ! MJ test20100929
      real(kind=4) :: rbuff_07, rfbuf_07  ! MJ test20100929
      real(kind=4) :: rbuff_09, rfbuf_09  ! MJ test20100929
      real(kind=4) :: cl_flag, outbufvg(20)  ! MJ test20110727
!     integer(kind=4) :: iimmod, jjmmod   ! MJ test20100917
      
      integer, dimension(Fmax,Lmax)      :: gzone
      character(len=255)  ::  dbdt_file
      character(256)  ::  msg
      integer         ::  status, doy, month, cnt
      integer         ::  i1, i2, j1, j2
      
      logical         ::  contains_dateline
      
      logical llprint(20)
      character*80 fname
      
      character(8)    ::  InstrumentMode, char_buf, Platform
      character(256)  ::  HDFAttrName,ECSParmName
      integer(kind=4) ::  rtn, Vrsn, idalg  ! "idalg" added by MJ 
      integer(kind=4) ::  pgs_met_getpcattr_s
      integer(kind=4) ::  LRN_Geo
      parameter          (LRN_Geo=600000)
      
!      INTEGER(kind=4) ::  pgs_io_gen_closef,   pgs_io_gen_openf, &
!                          pgs_pc_getreference, pgs_met_getpcattr_s, add_db_mod04, pgs_met_getpcattr
      INTEGER(kind=4) ::  pgs_pc_getreference
                              
      include 'sample90.inc'
      include 'table90.inc'
      include 'contrl.inc'
    
      common /bufout/ realbuf(26)
      common/lpoly/ xzlog(10),xlog(8),densol(4,7),denscn(4,5), &
                    cthet0(4),ctheta(4),cofs(16),indsol,indscn,iofset, &
                    p1,pr,dum10(10)
!     ---common parameter      
      common  /xday/ doy


      data xxzlog/0.0000000, 0.00977964, 0.0395086, 0.0904221, &
                 0.164818, 0.266515, 0.401776, 0.581261, &
                 0.824689, 1.17436/
!     -- these numbers correspond to satellite zenith angle
!     -- node points of 0.,16.,30.,40.,48.,54.,58.,60. degrees
      data xxlog /0.000000,0.0395086,0.143841,0.266515,0.401776, &
                 0.531394,0.635031,0.693147/
!
      data llprint/20*.false./
      
      real :: xtmp, cc, psi, sca
      integer:: itmp, ilat5, ilon5
      real, parameter ::  d2r = 3.14159/180.
      
!     -- initialize stuff
!		  open(1000, file="aot550.txt", form="formatted", status="new")
      end_file = .false.
      error_file = .false.
      maxiidx = 0
      maxjidx = 0
      do i = 1,20
        lprint(i) = llprint(i)
      enddo
      do i = 1,10
        xzlog(i) = xxzlog(i)
      enddo
      do i = 1,8
        xlog(i) = xxlog(i)
      enddo
      
      naot_avg(:,:,:) = 0
      iaot_avg(:,:,:) = -9999
      aot_avg(:,:,:)  = 0.0
      iref_avg(:,:,:) = -9999
      ref_avg(:,:,:)  = 0.0
      nssa_avg(:,:,:) = 0
      issa_avg(:,:,:) = -9999
      isr_avg(:,:,:)  = -9999
      ssa_avg(:,:,:)  = 0.0
      sr_avg(:,:,:)   = 0.0
      nsr_avg(:,:,:)  = 0
      nref_avg(:,:,:) = 0
      
      nae_avg(:,:) = 0
      ae_avg(:,:) = 0.0
      iae_avg(:,:) = -9999
      naot550_avg(:,:) = 0
      aot550_avg(:,:) = 0.0
      iaot550_avg(:,:) = -9999
      iaot550_best(:,:) = -9999
      r412_count(:,:) = 0
      t650_count(:,:) = 0
      aot550_min(:,:) = 999.
      aot550_max(:,:) = -999.
      sd550(:,:) = -999.
      isd550(:,:) = -9999
      iest_uncert(:,:) = -9999
      view(:,:) = -999.
      scat_ang(:,:) = -999.
      sfc_typ(:,:) = -999.
      xsza(:,:) = -999.
      r412b_y40(:,:) = -999.
      icld_frac(:,:) = 1000
      iDTDBaot(:,:) = -9999
      iDTDBqa(:,:) = 0
      iDTDBflag(:,:) = -999   
      Quality_Dt(:,:) = 0
      Quality_Dt_O(:,:) = 0   
      lat_sav(:,:) = 0.0
      lon_sav(:,:) = 0.0
      sza_sav(:,:) = 0.0
      vza_sav(:,:) = 0.0
      dstar_avg(:,:) = 0.0
      amf(:,:) = -999.0
      nll_sav(:,:) = 0
      gzone(:,:) = -999.0
      
!     -- initialize algorithm flag arrays
      alg_cnt(:,:,:)  = 0
      alg_flag(:,:)   = -999
      conf_flag(:,:)  = 0
      usefulness(:,:) = 0
      aero_type(:,:)  = 0
      rcond(:,:)      = 0
      
!     -- calculate values needed in table interpolation

      do j = 1, 7
        do k = j, j + 3
          xdenom = 1.0
          do i = j, j + 3
            if (i .ne. k) xdenom = xdenom * (xzlog(k) - xzlog(i))
          end do
          densol(k - j + 1, j) = xdenom
        end do
      end do
      do j = 1, 5
        do k = j, j + 3
          xdenom = 1.0
          do i = j, j + 3
            if (i .ne. k) xdenom = xdenom * (xlog(k) - xlog(i))
          end do
          denscn(k - j + 1, j) = xdenom
        end do
      end do

!     Determine platform
			HDFAttrName = 'CoreMetadata.0'
      ECSParmName = 'ASSOCIATEDPLATFORMSHORTNAME.1'
      Vrsn        = 1

      rtn=pgs_met_getpcattr_s(LRN_GEO,Vrsn,HDFAttrName,ECSParmName,char_buf)
      Platform 		= char_buf
!      print *, 'Platform is:', Platform

!     -- rewind measurement file
      call rewind_data(fh)

!     Here's the main processing loop
!     First read in a line of data
!      open(10000, file="ndvi.txt", form="FORMATTED")
10    call extract_data(fh,idim,jdim,xlat,xlong,sza,xthet,xphi,xnvalm5, &
                      cl_flag,vazim,sazim,stdv,LandSeaFlag,Dstar1,tmpvg,ierr)
!c.                   band26,vazim,sazim,LandSeaFlag,Dstar1,tmpvg,ierr)      
      if (ierr .lt. 0) goto 810
      outbuf(:) = -999.0
      
      call read_data(error_file,xnvalm5)
            
      iidx = idim/10 + 1
      jidx = jdim/10 + 1
      if (iidx .gt. maxiidx) maxiidx = iidx
      if (jidx .gt. maxjidx) maxjidx = jidx
      if (error_file) then
         go to 10
      endif
      
      if (LandSeaFlag.eq.0 .or. LandSeaFlag.gt.4) go to 10
      if (LandSeaFlag.eq.3) go to 10
      
      pcloud = 0.7
      isno   = 0
      icat   = 0

!     Check and set a snow/ice flag

      isnow = 0

!     Parse the input data into a form deepblue is comfortable with
 			itmp = xlat
      xtmp = xlat - itmp
      ilat5 = 10*(itmp+90) + 10*xtmp
      ilat5 = ilat5 + 1
      if (ilat5.gt.1800) ilat5 = 1800
      if (ilat5.lt.1)   ilat5 = 1

      itmp = xlong
      xtmp = xlong - itmp
      ilon5 = 10*(itmp+180) + 10*xtmp
      ilon5 = ilon5 + 1
      if (ilon5.gt.3600) ilon5 = 3600
      if (ilon5.lt.1)   ilon5 = 1
            
!     -- calc scattering angle
      cc     = 3.14159/180.
      psi    = acos(cos(sza*cc)*cos(xthet*cc) -       &
              sin(sza*cc)*sin(xthet*cc)*cos(xphi*cc))
      sca = 180. - psi/cc
!      print *, 'xlat, xlong, ilat5, ilon5: ', xlat, xlong, ilat5, ilon5
      sfref412 = get_LER412(ilat5, ilon5, realbuf(14), sca, xphi)/100.0
      sfref470 = get_LER470(ilat5, ilon5, realbuf(14), sca, xphi)/100.0
      sfref650 = get_LER650(ilat5, ilon5, realbuf(14), sca, xphi)/100.0
      call total
      
      realbuf(16) = Dstar1
!c.   realbuf(17) = band26  ! meant to be band26 or cl_flag? (MJ)
      realbuf(17) = cl_flag ! MJ modified (07/27/2011) @@@@@
      realbuf(14) = (tmpvg(5) - tmpvg(1)) / (tmpvg(5) + tmpvg(1))   ! TOA NDVI
!      write(10000,'3(F10.6,1X)') xlat, xlong, realbuf(14)

      if (realbuf(1) .lt. -900.0 .or. realbuf(2) .lt. -900.0) go to 10
      
!...+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!...  test by MJ (to write TOA mean refl. on MOD04_L2 hdf data), r20100929
      rbuff_11=-1.0*realbuf(11)/100.0
      rfbuf_11=3.14159*(10.0**rbuff_11)  ! B8, pi*L/F
      rbuff_07=-1.0*realbuf(7)/100.0
      rfbuf_07=3.14159*(10.0**rbuff_07)  ! B1?, pi*L/F
      rbuff_09=-1.0*realbuf(9)/100.0
      rfbuf_09=3.14159*(10.0**rbuff_09)  ! B3?, pi*L/F
!     iimmod=mod(idim,100)
!     jjmmod=mod(jdim,100)
!     if(idim.lt.1200.and.jdim.lt.1900.and.iimmod.eq.0.and.jjmmod.eq.0) then
!        print *, idim, jdim, rfbuf_11, rfbuf_11*3.14159/cos(sza*3.14159/180.0)
!     endif
!...+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!     Extract reflectance from realbuf (are ref*100)
!     refbuf(1) = realbuf(22)  ! original; LER at 412nm (MJ)
!     refbuf(2) = realbuf(23)  ! original; LER at 470nm (MJ)
!     refbuf(3) = realbuf(6)   ! original
      
      refbuf(1) = rfbuf_11*100.  ! MJ test TOA refl(20100929), B8
      refbuf(2) = rfbuf_09*100.  ! MJ test TOA refl(20101105), B3
      refbuf(3) = rfbuf_07*100.  ! MJ test TOA refl(20101105), B1

!c. ################################################################
!c. ################################################################
!c. Call over-vegetation retrieval (test version by MJ; 07/28/2011)
!c.   Position-A
!     call find_v_veg_test(idim,jdim,realbuf,tmpvg,outbufvg)
!c.   call find_v_veg_test(idim,jdim,realbuf,tmpvg,outbufvg,qa_flag)
!c. ################################################################
!c. ################################################################
!     Call the Deep Blue algorithm and dump output to file
      If (Platform .EQ. 'Terra') Then
				call find_v_terra(realbuf, outbuf, tmpvg, flags)
				
      Else If (Platform .EQ. 'Aqua') Then
				call find_v_aqua(realbuf, outbuf, tmpvg, flags)
      
      End If
!			if (abs(xlat-24.0) < 0.01 .AND. abs(xlong-25.0) < 0.01) then
 !     	print *, 'out: ', xlat, xlong, outbuf
!			end if
!c. ################################################################
!c. ################################################################
!c. Call over-vegetation retrieval (test version by MJ; 07/29/2011)
!c.   Position-B
!c    print *, 'position B calling find_v_veg_test'

!       use this one
!      call find_v_veg_test(idim,jdim,realbuf,tmpvg,outbufvg)

!c    print *, 'position C end of calling find_v_veg_test'
!c.   call find_v_veg_test(idim,jdim,realbuf,tmpvg,outbufvg,qa_flag)
!c. ################################################################
!c. ################################################################

!c. ################################################################
!c. ################################################################
!c. ################################################################
!c... A routine to combine standard DeepBlue and over-vegetation 
!c... retrievals may be placed here... (MJ; 08/05/2011) @@@@@

!      call cmb_stddb_veg(outbuf, outbufvg, tmpvg, idalg)
!c...  --> "outbuf" is input and output, "idalg" is output.
            
!     @TODO do we need to consider mixing the two?
!     -- record algorithm flag for this pixel
      select case (int(outbuf(18)))
        case (0)      ! Deep Blue
          alg_cnt(iidx,jidx,1) = alg_cnt(iidx,jidx,1) + 1
        case (1)      ! Vegetated
          alg_cnt(iidx,jidx,2) = alg_cnt(iidx,jidx,2) + 1
        case (-999)   ! no retrieval
        case default
          print *, "ERROR: Invalid algorithm flag: ", outbuf(18)
          go to 10
      end select
!      print *, 'iidx, jidx, outbuf(18): ', iidx, jidx, int(outbuf(18)), outbuf(18)
!     -->  check if gas correction is on in "DeepBlue.f90".
!c. ################################################################
!c. ################################################################
!c. ################################################################

!     Extract and average output

!      if (outbuf(1) .gt. -990 ) then 
!      	print *, 'outbuf(1,2,3,7) ',outbuf(1), outbuf(2),outbuf(3),outbuf(7)
!      endif
!      write(1000,'(3(F12.6,1x))'), xlat, xlong, outbuf(7)
      do i=1,3
         if (outbuf(i) .ge. 0.0 .and. outbuf(i) .le. 5.0) then
            aot_avg(iidx,jidx,i)  = aot_avg(iidx,jidx,i) + outbuf(i)
            naot_avg(iidx,jidx,i) = naot_avg(iidx,jidx,i) + 1
!            aot_list(iidx,jidx,i,naot_avg(iidx,jidx,i)) = outbuf(i)
            !if (i == 2 .AND. (iidx == 110 .AND. jidx == 188)) then
            !  print *, 'aot470, naot470: ', iidx, jidx, outbuf(i), aot_avg(iidx,jidx,i), naot_avg(iidx,jidx,i)
            !endif
         endif
         if (refbuf(i) .gt. 0.0) then 
					 ref_avg(iidx,jidx,i)  = ref_avg(iidx,jidx,i) + refbuf(i)
					 nref_avg(iidx,jidx,i) = nref_avg(iidx,jidx,i) + 1
		 		 endif
         if (outbuf(i+3) .ge. 0.0 .and. outbuf(i+3) .le. 5.0) then
            ssa_avg(iidx,jidx,i)  = ssa_avg(iidx,jidx,i) + outbuf(i+3)
            nssa_avg(iidx,jidx,i) = nssa_avg(iidx,jidx,i) + 1
         endif
     enddo
      
      if (outbuf(7) .gt. -900.0) then
         aot550_avg(iidx,jidx) = aot550_avg(iidx,jidx) + outbuf(7)
         naot550_avg(iidx,jidx) = naot550_avg(iidx,jidx) + 1
         if (aot550_min(iidx,jidx) .gt. outbuf(7)) aot550_min(iidx,jidx)=outbuf(7)
         if (aot550_max(iidx,jidx) .lt. outbuf(7)) aot550_max(iidx,jidx)=outbuf(7)
         aot550_list(iidx,jidx,naot550_avg(iidx,jidx)) = outbuf(7)
         
         dstar_avg(iidx,jidx) = dstar_avg(iidx,jidx) + Dstar1
      endif

      if (outbuf(8) .gt. -900.0) then
         ae_avg(iidx,jidx) = ae_avg(iidx,jidx) + outbuf(8)
         nae_avg(iidx,jidx) = nae_avg(iidx,jidx) + 1
      endif

      if (outbuf(9) .gt. 0.12) r412_count(iidx,jidx) = r412_count(iidx,jidx)+1
      if (outbuf(10) .gt. 0.) t650_count(iidx,jidx) = t650_count(iidx,jidx)+1

      if (outbuf(9) .gt. 0.0) then
        if (outbuf(18) == 0) then    ! SR412 is only defined for DB retrievals.
         sr_avg(iidx,jidx,1) = sr_avg(iidx,jidx,1) + outbuf(9)/100
         nsr_avg(iidx,jidx,1) = nsr_avg(iidx,jidx,1) + 1
        end if
      endif
      if (outbuf(11) .gt. 0.0) then
         sr_avg(iidx,jidx,2) = sr_avg(iidx,jidx,2) + outbuf(11)/100
         nsr_avg(iidx,jidx,2) = nsr_avg(iidx,jidx,2) + 1
      endif
      if (outbuf(12) .gt. 0.0) then
         sr_avg(iidx,jidx,3) = sr_avg(iidx,jidx,3) + outbuf(12)/100
         nsr_avg(iidx,jidx,3) = nsr_avg(iidx,jidx,3) + 1
      endif

!     Want values for central pixel (or nearby) for three variables

      if (mod(idim,10) .lt. 1. .or. mod(jdim,10) .lt. 1. .or. &
              mod(idim,10) .gt. 8.5 .or. mod(jdim,10) .gt. 8.5) go to 10

      if (mod(idim,10) .le. 5 .and. mod(jdim,10) .le. 5) then
         if (outbuf(13) .gt. -900.) then
            view(iidx,jidx) = outbuf(13)
         endif
  
         if (outbuf(14) .gt. -900.) then
            scat_ang(iidx,jidx) = outbuf(14)
         endif

         if (outbuf(15) .gt. -900.) then
            sfc_typ(iidx,jidx) = outbuf(15)
         endif

         if (outbuf(16) .gt. -900.) then
            xsza(iidx,jidx) = outbuf(16)
         endif

         if (outbuf(17) .gt. -900.) then
            r412b_y40(iidx,jidx) = outbuf(17)
         endif
      endif
      
      if (mod(idim,10) .gt. 5 .or. mod(jdim,10) .gt. 5) then
         if (outbuf(13) .gt. -900. .and. view(iidx,jidx) .lt. -900.) then
            view(iidx,jidx) = outbuf(13)
         endif
  
         if (outbuf(14) .gt. -900. .and. scat_ang(iidx,jidx) .lt. -900.) then
            scat_ang(iidx,jidx) = outbuf(14)
         endif

         if (outbuf(15) .gt. -900. .and. sfc_typ(iidx,jidx) .lt. -900.) then
            sfc_typ(iidx,jidx) = outbuf(15)
         endif

         if (outbuf(16) .gt. -900. .and. xsza(iidx,jidx) .lt. -900.) then
            xsza(iidx,jidx) = outbuf(16)
         endif

         if (outbuf(17) .gt. -900. .and. r412b_y40(iidx,jidx) .lt. -900.) then
            r412b_y40(iidx,jidx) = outbuf(17)
         endif
      endif

      go to 10
      

! End of the main loop - that's pretty short

810   continue

!close(10000)
print *, "End of pixel processing."
!print *, 'max i,j: ', maxiidx, maxjidx

!   -- calculate average lat/lon for each 10x10 pixel cell.
    contains_dateline = .false.
    if (maxval(longitude) > 179.0 .AND. minval(longitude) < -179.0) contains_dateline = .true.
    do i = 1, size(latitude,1)
      do j = 1, size(latitude,2)

!     -- skip undefined lat/lon values.
      if (latitude(i,j) < -900.0 .OR. longitude(i,j) < -900.0) cycle
      
      iidx = i/10 + 1
      jidx = j/10 + 1
    
      lat_sav(iidx,jidx) = lat_sav(iidx,jidx) + latitude(i,j)
      
!     -- convert longitudes to all positive values, otherwise values along dateline
!     -- will cancel each other out.
      if (contains_dateline .AND. longitude(i,j) < 0.0) then
        lon_sav(iidx,jidx) = lon_sav(iidx,jidx) + longitude(i,j) + 360.0
      else
        lon_sav(iidx,jidx) = lon_sav(iidx,jidx) + longitude(i,j)
      end if
          
      sza_sav(iidx,jidx) = sza_sav(iidx,jidx) + solar_zenith_angle(i,j)
      vza_sav(iidx,jidx) = vza_sav(iidx,jidx) + sensor_zenith_angle(i,j)
      
      nll_sav(iidx,jidx) = nll_sav(iidx,jidx) + 1

      end do
    end do
    where (nll_sav > 0) 
      lat_sav = lat_sav / nll_sav
      lon_sav = lon_sav / nll_sav
      sza_sav = sza_sav / nll_sav
      vza_sav = vza_sav / nll_sav
    elsewhere
      lat_sav = -999.0
      lon_sav = -999.0
      sza_sav = -999.0
      vza_sav = -999.0
    end where
    
!   -- back out our longitude conversion to positive values where needed.
    if (contains_dateline) then
      where (lon_sav > 180.0)  
        lon_sav = lon_sav - 360.0
      end where 
    end if
    
! Compute AOT and SSA averages
    do j = 1, maxjidx
      do i = 1, maxiidx
         
            do k = 1, 3
               
               if (naot_avg(i,j,k) .gt. 9) then
                  iaot_avg(i,j,k) = 1000*aot_avg(i,j,k)/naot_avg(i,j,k)
               else 
                  iaot_avg(i,j,k) = -9999
                  iref_avg(i,j,k) = -9999
               endif
               if (nref_avg(i,j,k) .gt. 9 .and. ref_avg(i,j,k) .gt. 0) then
                 iref_avg(i,j,k) = 100*ref_avg(i,j,k)/nref_avg(i,j,k)
               else
                 iref_avg(i,j,k) = -9999
               endif
               
               if (nssa_avg(i,j,k) .gt. 9) then
                  issa_avg(i,j,k) = 1000*ssa_avg(i,j,k)/nssa_avg(i,j,k)
               else 
                  issa_avg(i,j,k) = -9999
               endif
               if (nsr_avg(i,j,k) .gt. 0) then
                  isr_avg(i,j,k) = 1000*sr_avg(i,j,k)/nsr_avg(i,j,k)
               else
                  isr_avg(i,j,k) = -9999
               endif
            enddo
            
            if (nae_avg(i,j) .gt. 9) then
               ae_avg(i,j) = ae_avg(i,j)/nae_avg(i,j)
            else
               ae_avg(i,j) = -999.
            endif
            if (naot550_avg(i,j) .gt. 9) then
               aot550_avg(i,j) = aot550_avg(i,j)/naot550_avg(i,j)
               dstar_avg(i,j) = dstar_avg(i,j)/naot550_avg(i,j)
            else
               aot550_avg(i,j) = -999.
               naot550_avg(i,j) = 0
               dstar_avg(i,j) = -999.0
            endif

!           -- set cell algorithm flag to mode of alg_cnt(i,j).
            if (max(alg_cnt(i,j,1), alg_cnt(i,j,2)) .gt. 0) then
              if (alg_cnt(i,j,2) == 0) then       ! no vegetated retrievals, set to DB
                alg_flag(i,j) = 0
              elseif (alg_cnt(i,j,1) == 0) then   ! no DB retrievals, set to vegetated
                alg_flag(i,j) = 1
              else
                alg_flag(i,j) = 2
              end if
            else
              alg_flag(i,j) = -999        ! no retrievals in this cell
            end if
!            print *, 'alg flag: ', i, j, alg_cnt(i,j,1), alg_cnt(i,j,2), alg_flag(i,j)
            
          
          enddo          
      enddo
      
!     -- Calc standard deviation of tau550
      do j = 1, maxjidx
        do i = 1, maxiidx

          if (naot550_avg(i,j) .gt. 9 .and. aot550_avg(i,j) .gt. -900.0) then
            var = 0.0       
            do m=1,naot550_avg(i,j)
              var = var + (aot550_avg(i,j) - aot550_list(i,j,m))**2              
            enddo
            var = var/(naot550_avg(i,j)-1)
            sd550(i,j) = sqrt(var)
            isd550(i,j) = 1000*sd550(i,j)
          endif
            
         enddo
      enddo   
      
! Calculate tau550 and alpha
      do j = 1, maxjidx
        do i = 1, maxiidx   
!         -- get geozone
          ilat5 = (lat_sav(i,j) + 90.) *10.
          ilat5 = ilat5 + 1
          if (ilat5.gt.1800) ilat5 = 1800
          if (ilat5.lt.1)   ilat5 = 1

          ilon5 = (lon_sav(i,j) + 180.) *10.
          ilon5 = ilon5 + 1
          if (ilon5.gt.3600) ilon5 = 3600
          if (ilon5.lt.1)   ilon5 = 1
          
          gzone(i,j) = terrain_flag_new(ilon5,ilat5)
          
!           -- for vegetated retrievals, keep current AOT550 values -- CWB 09/28/11
!            if (alg_flag(i,j) == 1) then	
!              iaot550_avg(i,j) = aot550_avg(i,j) * 1000
!              cycle
!            end if
  
!           -- calculate AE for 550nm
            alpha_550  = -999.
            alpha2_550 = -999.
            xxrat  = (1.0*iaot_avg(i,j,1)) / (1.0*iaot_avg(i,j,2))
            xxrat2 = (1.0*iaot_avg(i,j,3)) / (1.0*iaot_avg(i,j,2))
            if (iaot_avg(i,j,1) .gt. 0 .and. xxrat .gt. 0.0) then
               dd     = alog(412./470.)
               alpha_550  = alog(xxrat)
               alpha_550  = -1.*alpha_550/dd
               if (iaot_avg(i,j,2) .lt. 200) alpha_550  = 1.0
               if (alpha_550 .gt. 1.8) alpha_550  = 1.8
            endif
            if (iaot_avg(i,j,3) .gt. 0 .and. xxrat2 .gt. 0.0) then
               dd     = alog(650./470.)
               alpha2_550 = alog(xxrat2)
               alpha2_550 = -1.*alpha2_550/dd
               if (iaot_avg(i,j,2) .lt. 200) alpha2_550  = 1.0
               if (alpha2_550 .gt. 1.8) alpha2_550  = 1.8
            endif
            if (alpha_550 .lt. -900.0 .and. alpha2_550 .gt. -900.0) alpha_550 = alpha2_550
            if (iaot_avg(i,j,1) .lt. 0 .and. iaot_avg(i,j,2) .lt. 0) then
               if (ae_avg(i,j) .gt. -900.0) alpha_550 = ae_avg(i,j)
            endif
            if (aot550_avg(i,j) .gt. 1.8 .and. ae_avg(i,j) .gt. -900.0) then
               alpha_550 = ae_avg(i,j)
            endif
            if (alpha_550 .lt. -0.4 .and. alpha_550 .gt. -900.0) alpha_550 = -0.4
            if (alpha_550 .gt. 1.8) alpha_550  = 1.8
            
!           -- case of thin dust over very bright surface. Reset AE to 1.0.
            if (iaot_avg(i,j,2) .lt. 700 .and. isr_avg(i,j,3) .gt. 250 .and. &
                alpha_550 .gt. 1.0) alpha_550  = 1.0
                
!           -- now calculate AE for the rest of the bands
            alpha  = -999.
            alpha2 = -999.
            xxrat  = (1.0*iaot_avg(i,j,1)) / (1.0*iaot_avg(i,j,2))
            xxrat2 = (1.0*iaot_avg(i,j,3)) / (1.0*iaot_avg(i,j,2))
            if (iaot_avg(i,j,1) .gt. 0 .and. xxrat .gt. 0.0) then
               dd     = alog(412./470.)
               alpha  = alog(xxrat)
               alpha  = -1.*alpha/dd
               if (iaot_avg(i,j,2) .lt. 200) alpha  = 1.5
               if (alpha .gt. 1.8) alpha  = 1.8
            endif
            if (iaot_avg(i,j,3) .gt. 0 .and. xxrat2 .gt. 0.0) then
               dd     = alog(650./470.)
               alpha2 = alog(xxrat2)
               alpha2 = -1.*alpha2/dd
               if (iaot_avg(i,j,2) .lt. 200) alpha2  = 1.5
               if (alpha2 .gt. 1.8) alpha2  = 1.8
            endif
            if (alpha .lt. -900.0 .and. alpha2 .gt. -900.0) alpha = alpha2
            if (iaot_avg(i,j,1) .lt. 0 .and. iaot_avg(i,j,2) .lt. 0) then
               if (ae_avg(i,j) .gt. -900.0) alpha = ae_avg(i,j)
            endif
            if (aot550_avg(i,j) .gt. 1.8 .and. ae_avg(i,j) .gt. -900.0) then
               alpha = ae_avg(i,j)
            endif
            if (alpha .lt. 0.0 .and. alpha .gt. -900.0) alpha = 0.0
            if (alpha .gt. 1.8) alpha  = 1.8

!           -- case of thin dust over very bright surface. Reset AE to 1.0.
            if (iaot_avg(i,j,2) .lt. 700 .and. isr_avg(i,j,3) .gt. 250 .and. &
                alpha .gt. 1.0) alpha  = 1.0
      
!           -- D* says dust. Set AE to dust. Only over N.Africa and Arabia.
            if ((alpha > -900.0 .AND. dstar_avg(i,j) > 1.02) .AND. &
            &   ((gzone(i,j) >= 1 .AND. gzone(i,j) <= 5) .OR. (gzone(i,j) == 26 .OR. gzone(i,j) == 27))) then
              alpha = 0.0
            end if
            
!           -- if AOT is maxed out for all channels, AE will be 0.0. But what if
!           -- it's smoke? Use D* to detect smoke (D* < 0.97) and set AE accordingly.
            if (iaot_avg(i,j,1) .gt. 3470 .AND. dstar_avg(i,j) .lt. 1.06) then
              alpha = 1.8
            endif
            
!           -- apply AE to AOT and get everything consistent.
            if (aot550_avg(i,j) .gt. 0.0 .and. alpha .gt. -900.) then 
               dd     = 550./500.
               dd     = dd**(-1.*alpha_550)
               iaot550_avg(i,j) = aot550_avg(i,j) * dd * 1000
               aot550_avg(i,j)  = aot550_avg(i,j) * dd
               
               dd     = 412./550.
               dd     = dd**(-1.*alpha)             
               iaot_avg(i,j,1)  = aot550_avg(i,j) * dd * 1000
               dd     = 470./550.
               dd     = dd**(-1.*alpha)
               iaot_avg(i,j,2)  = aot550_avg(i,j) * dd * 1000
               dd     = 650./550.
               dd     = dd**(-1.*alpha)
               iaot_avg(i,j,3)  = aot550_avg(i,j) * dd * 1000
            else
               do k = 1, 3
                iaot_avg(i,j,k)   = -9999
               enddo
               iaot550_avg(i,j)  = -9999
            endif

            if (alpha .lt. 0.0) then
               iae_avg(i,j) = -9999
            else
               iae_avg(i,j) = alpha*1000
            endif

         enddo
         
      enddo

!     -- reset high AOT's to 3.5
      where (iaot_avg > 3500)    iaot_avg    = 3500
      where (aot550_avg > 3.5)   aot550_avg  = 3.5
      where (iaot550_avg > 3500) iaot550_avg = 3500
      
! Cloud edge screening
      do k = 1, 3
        do j = 2, maxjidx - 1
          do i = 2, maxiidx - 1
               if (iaot_avg(i,j,k) .gt. 1200 .and. &
                   abs(iaot_avg(i-1,j,k)-iaot_avg(i,j,k)) .gt. 800 .and. &
                   abs(iaot_avg(i,j,k)-iaot_avg(i+1,j,k)) .gt. 800) then 
                  iaot_avg(i,j,:) = -9999
                  aot550_avg(i,j) = -999.0
                  iaot550_avg(i,j) = -9999
                  iae_avg(i,j) = -9999
               endif
               if (iaot_avg(i,j,k) .gt. 1200 .and. &
                   abs(iaot_avg(i,j+1,k)-iaot_avg(i,j,k)) .gt. 800 .and. &
                   abs(iaot_avg(i,j,k)-iaot_avg(i,j-1,k)) .gt. 800) then 
                  iaot_avg(i,j,:) = -9999
                  aot550_avg(i,j) = -999.0
                  iaot550_avg(i,j) = -9999
                  iae_avg(i,j) = -9999

               endif
            enddo
         enddo
      enddo
      
      do j = 2, maxjidx - 1
        do i = 2, maxiidx - 1
         
               if (aot550_avg(i,j) .gt. 1.2 .and. &
                   abs(aot550_avg(i-1,j)-aot550_avg(i,j)) .gt. 0.8 .and. &
                   abs(aot550_avg(i,j)-aot550_avg(i+1,j)) .gt. 0.8) then
                  aot550_avg(i,j) = -999.
                  iaot550_avg(i,j) = -9999
                  iaot_avg(i,j,:) = -9999
                  iae_avg(i,j) = -9999

               endif
               if (aot550_avg(i,j) .gt. 1.2 .and. &
                   abs(aot550_avg(i,j+1)-aot550_avg(i,j)) .gt. 0.8 .and. &
                   abs(aot550_avg(i,j)-aot550_avg(i,j-1)) .gt. 0.8) then
                  aot550_avg(i,j) = -999.
                  iaot550_avg(i,j) = -9999
                  iaot_avg(i,j,:) = -9999
                  iae_avg(i,j) = -9999
               endif
         enddo
      enddo
      
!     -- calculate and save a rudimentary cloud fraction.
      do j = 1, maxjidx
        do i = 1, maxiidx
       
          if (naot550_avg(i,j) .gt. -900) then
            icld_frac(i,j) = 1000*((100 - naot550_avg(i,j))/100.0)
          end if
          
        end do 
      end do 
      
! Replace all variables with fill value if aot550_avg undefined
      do j = 1, maxjidx
        do i = 1, maxiidx
        
        
					 if (iaot550_avg(i,j) .lt. -900) then
							do k = 1, 3
								 iaot_avg(i,j,k) = -9999
								 issa_avg(i,j,k) = -9999
								 isr_avg(i,j,k)  = -9999
								 iref_avg(i,j,k) = -9999
								 naot_avg(i,j,k) = 0
							enddo ! k
							iae_avg(i,j) = -9999
							isd550(i,j)  = -9999
							
							alg_flag(i,j) = -999
							icld_frac(i,j) = 1000
							dstar_avg(i,j) = -999.0
							aot550_avg(i,j) = -999.0
					 endif

         enddo ! i
      enddo ! j

    edgecnt(:,:) = 0
    do j = 2, maxjidx-1
      do i = 2, maxiidx-1
        edgecnt(i,j) = 0
        
        if (aot550_avg(i,j) .gt. -1 .AND. aot550_avg(i,j) .lt. 1.0) then

          if (abs(aot550_avg(i,j)-aot550_avg(i+1,j)) .gt. 0.4) then
            edgecnt(i,j)=edgecnt(i,j)+1
          endif
          if (abs(aot550_avg(i,j)-aot550_avg(i-1,j)) .gt. 0.4) then
            edgecnt(i,j)=edgecnt(i,j)+1
          endif
          if (abs(aot550_avg(i,j)-aot550_avg(i,j+1)) .gt. 0.4) then
            edgecnt(i,j)=edgecnt(i,j)+1
          endif
          if (abs(aot550_avg(i,j)-aot550_avg(i,j-1)) .gt. 0.4) then
            edgecnt(i,j)=edgecnt(i,j)+1
          endif
          
          if (abs(aot550_avg(i,j)-aot550_avg(i+1,j+1)) .gt. 0.4) then
            edgecnt(i,j)=edgecnt(i,j)+1
          endif
          if (abs(aot550_avg(i,j)-aot550_avg(i+1,j-1)) .gt. 0.4) then
            edgecnt(i,j)=edgecnt(i,j)+1
          endif
          if (abs(aot550_avg(i,j)-aot550_avg(i-1,j+1)) .gt. 0.4) then
            edgecnt(i,j)=edgecnt(i,j)+1
          endif
          if (abs(aot550_avg(i,j)-aot550_avg(i-1,j-1)) .gt. 0.4) then
            edgecnt(i,j)=edgecnt(i,j)+1
          endif
        else if (aot550_avg(i,j) .gt. -1.0 .AND. aot550_avg(i,j) .ge. 1.0) then
          if (abs(aot550_avg(i,j)-aot550_avg(i+1,j))/aot550_avg(i,j) .gt. 0.3) then
            edgecnt(i,j)=edgecnt(i,j)+1
          endif
          if (abs(aot550_avg(i,j)-aot550_avg(i-1,j))/aot550_avg(i,j) .gt. 0.3) then
            edgecnt(i,j)=edgecnt(i,j)+1
          endif
          if (abs(aot550_avg(i,j)-aot550_avg(i,j+1))/aot550_avg(i,j) .gt. 0.3) then
            edgecnt(i,j)=edgecnt(i,j)+1
          endif
          if (abs(aot550_avg(i,j)-aot550_avg(i,j-1))/aot550_avg(i,j) .gt. 0.3) then
            edgecnt(i,j)=edgecnt(i,j)+1
          endif
          
          if (abs(aot550_avg(i,j)-aot550_avg(i+1,j+1))/aot550_avg(i,j) .gt. 0.3) then
            edgecnt(i,j)=edgecnt(i,j)+1
          endif
          if (abs(aot550_avg(i,j)-aot550_avg(i+1,j-1))/aot550_avg(i,j) .gt. 0.3) then
            edgecnt(i,j)=edgecnt(i,j)+1
          endif
          if (abs(aot550_avg(i,j)-aot550_avg(i-1,j+1))/aot550_avg(i,j) .gt. 0.3) then
            edgecnt(i,j)=edgecnt(i,j)+1
          endif
          if (abs(aot550_avg(i,j)-aot550_avg(i-1,j-1))/aot550_avg(i,j) .gt. 0.3) then
            edgecnt(i,j)=edgecnt(i,j)+1
          endif
        end if
      end do
    end do
    
! Set QA Flags
      do j = 1, maxjidx
        do i = 1, maxiidx
            
            flags(1) = 0
            flags(2) = 0
            flags(3) = 0
            flags(4) = 0
            if (iaot550_avg(i,j).gt.-1) then

               flags(1) = 1
               flags(2) = 1

!             -- second to last parameter to set_qa is supposed to be geozone. currently disabled, fix later.
!              print *, 'in: ', 1,0.0,naot550_avg(i,j),0.0,view(i,j),aot550_avg(i,j),scat_ang(i,j),0,0,flags(2)
              call set_qa(1,0.0,naot550_avg(i,j),sd550(i,j),view(i,j), &
                      aot550_avg(i,j),scat_ang(i,j),edgecnt(i,j),0, &
                      gzone(i,j),dstar_avg(i,j),flags(2))
!							 If (Platform .EQ. 'Terra') Then
!							   call qa_flag_terra(sfc_typ(i,j),naot550_avg(i,j),naot550_avg(i,j),sd550(i,j),view(i,j), &
!                                realiaot550_avg(i,j),scat_ang(i,j),xsza(i,j), &
!                                r412b_y40(i,j),flags(2))
!					    
!							 Else If (Platform .EQ. 'Aqua') Then
!							   call qa_flag_aqua(sfc_typ(i,j),naot550_avg(i,j),naot550_avg(i,j),sd550(i,j), &
!							                view(i,j), realiaot550_avg(i,j),scat_ang(i,j),flags(2))
!               End If
               if (iae_avg(i,j).le.500) flags(3)= 1
               if (iae_avg(i,j).gt.1200) flags(3)= 2
               if (r412_count(i,j) .gt. 20) flags(4) = 1
               if (naot_avg(i,j,2) .lt. 50) flags(4) = 2
               if (t650_count(i,j) .gt. 50) flags(4) = 3

!               if (view(i,j) .lt. -900 .or. scat_ang(i,j) .lt. -900) then
!                  print *, "problem",naot_avg(i,j,1),view(i,j),scat_ang(i,j),sfc_typ(i,j)
!                  print *, flags
!               endif
               
            endif

!            call set_bits(flags(1),0,QAFlags(5,i,j))
!            call set_bits(flags(2),1,QAFlags(5,i,j))
!            call set_bits(flags(3),3,QAFlags(5,i,j))
!            call set_bits(flags(4),5,QAFlags(5,i,j))
            
            conf_flag(i,j)  = flags(2)
            usefulness(i,j) = flags(1)
            aero_type(i,j) = flags(3)
            rcond(i,j)    = flags(4)
            
         enddo
      enddo
      
!     -- reset QA to 1 if cell is adjacent to an empty (cloudy) pixel, i.e QA=0.
!     -- set usefulness flag to 0 if QA==0,1.  then save QA to MYD04 file.
      edgecnt(:,:) = 0
      do j = 1, maxjidx
        do i = 1, maxiidx

!         -- skip cells over N.Africa with Dstar > 1.2 -- indicate dust. Don't reset.
!         -- Otherwise, use 1.06 threshold elsewhere.
          if ((gzone(i,j) >= 1 .AND. gzone(i,j) <= 5) .OR. (gzone(i,j) == 26 .OR. gzone(i,j) == 27)) then
            if (dstar_avg(i,j) > 1.2) cycle
          else
            if (dstar_avg(i,j) >= 1.06) cycle
          endif
           
          i1 = max(i-1,1)
          i2 = min(i+1,maxiidx)
          j1 = max(j-1,1)
          j2 = min(j+1,maxjidx)
          if (conf_flag(i,j) /= 0) then
            cnt = count(conf_flag(i1:i2,j1:j2) == 0)
            if (cnt > 0) then
              conf_flag(i,j) = 1
            end if
          end if
        
	        !if (i>1 .AND. j>1 .AND. i<maxiidx .AND. j<maxjidx) then 
          !  if (conf_flag(i,j) == 0) then
          !    cnt = count(conf_flag(i-1:i+1,j-1:j+1) > 1)
          !    if (cnt >= 0) then
          !      do l = j-1, j+1
		      !        do k = i-1, i+1           
			    !          if (conf_flag(k,l) > 1) edgecnt(k,l) = 1
 			    !        end do
	        !      end do
  	      !    end if
          !  end if
	        !end if
	        
        end do
      end do
!      where (edgecnt == 1) conf_flag = 1

!     -- update usefulness flags
      where (conf_flag == 1) usefulness = 0
      
!     -- update QA flag array (QAFlags)
      do j = 1, Lmax
        do i = 1, Fmax
            call set_bits(usefulness(i,j),0,QAFlags(5,i,j))
            call set_bits(conf_flag(i,j),1,QAFlags(5,i,j))
            call set_bits(aero_type(i,j),3,QAFlags(5,i,j))
        end do
      end do
      
!     -- copy AOT550 value to best estimate array if QA=3 
      where (conf_flag == 3 .OR. conf_flag == 2) iaot550_best = iaot550_avg
      
!     -- calculate estimated uncertainty
      amf = (1.0/cos(d2r*sza_sav)) + (1.0/cos(d2r*vza_sav))
      where (conf_flag == 1 .AND. aot550_avg > -900)  
        iest_uncert = ((0.083 + 0.83*aot550_avg)/amf)*1000.0
      end where
      where (conf_flag == 2 .AND. aot550_avg > -900) 
        iest_uncert = ((0.10 + 0.60*aot550_avg)/amf)*1000.0
      end where
      where (conf_flag == 3 .AND. aot550_avg > -900) 
        iest_uncert = ((0.086 + 0.56*aot550_avg)/amf)*1000.0
      end where
      
! Create combined Dark Target/Deep Blue AOT550 SDS based on QA
! Extract Quality Flag - Code from Shana      
      do k=1,No_byte
        do i=1,Fmax
          do j=1,Lmax

            nflag=QAFlags(k,i,j)        
            call get_bits(nflag,n1,n2,n3,n4,n5,n6,n7,n8)
!            print *, nflag,n1,n2,n3,n4,n5,n6,n7,n8

      ! Find quality flag for Dtland(1) & Deep blue( 5th)             
            if( k.eq.1 ) then
              Quality_Dt(i,j)=(n2+2*n3+4*n4)
 !           print *,nflag,n1,(n2+2*n3+4*n4),n5,(n6+2*n7+4*n8)
            endif 
            
          end do
        end do
      end do

!     -- get ocean AOT quality flags      
      do k=1, No_byte_O
        do i= 1, Fmax
          do j=1, Lmax
!           -- get ocean QA flags
            nflag=QAFlags_O(k,i,j)        
            call get_bits(nflag,n1,n2,n3,n4,n5,n6,n7,n8)
            if( k.eq.1 ) then
              Quality_Dt_O(i,j)=(n2+2*n3+4*n4)
            endif       
          enddo
        enddo
      enddo

      Vrsn = 1
      rtn=pgs_pc_getreference(412008,Vrsn,dbdt_file)
      if (rtn .lt. 0) then
         msg = "Error retrieving LRN_DBDT lun from pcf."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      endif
      msg = "Will use DBDT region file "//dbdt_file
      call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC, msg, 'DeepBlue')
      
      if (doy < 31) then
        month = 1
      else if (doy <= 59) then
        month = 2
      else if (doy <= 90) then
        month = 3
      else if (doy <= 120) then
        month = 4
      else if (doy <= 151) then
        month = 5
      else if (doy <= 181) then
        month = 6
      else if (doy <= 212) then
        month = 7
      else if (doy <= 243) then
        month = 8
      else if (doy <= 273) then
        month = 9
      else if (doy <= 304) then
        month = 10
      else if (doy <= 334) then
        month = 11
      else if (doy <= 366) then
        month = 12
      end if

      
      status = load_dbdt_region_table(dbdt_file, month)
      if (status /= 0) then
        msg = "Error loading DBDT region input file."
        call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      end if
      
      status = create_dbdt_aot550(lat_sav, lon_sav, iaot550_best, conf_flag, &
      &                   iDTaot, Quality_Dt, iDTaot, Quality_Dt_O, iDTDBaot, &
      &                   iDTDBqa, iDTDBflag)
      if (status /= 0) then
        msg = "ERROR: Failed to combine DT and DB AOT550 values."
        print *, trim(msg), status
        call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      end if
      
      status = unload_dbdt_region_table()
      if (status /= 0) then
        msg = "Error unloading DBDT region input file."
        call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC, msg, 'DeepBlue')
      end if

!     -- update QA flag array (QAFlags) with DB/DT QA.
!     -- the usefulness flag (bit 0) is not currently used -- assumed that all 
!     -- points are useful, set to 1 everywhere.
      do j = 1, Lmax
        do i = 1, Fmax
        
          if (iDTDBqa(i,j) > 0) then
            call set_bits(1, 0, QAFlags(6,i,j))           ! usefulness flag
          else
            call set_bits(0, 0, QAFlags(6,i,j))
          end if
          
          call set_bits(iDTDBqa(i,j),1,QAFlags(6,i,j))  ! quality assurance flag
        end do
      end do
     
!      close(1000)
! Create new SDSs in MOD04 product
! THIS SHOULD BE COMMENTED OUT FOR OPERATIONAL VERSION
!     print *,'start add_db_mod04()...'
!     ierr = add_db_mod04()
!     print *, 'add_db_mod04 status: ', ierr
! Write averaged data and updated flags to MODO4 product
     call write_db(dims,QAFlags,iaot550_avg,iaot550_best,iae_avg,iaot_avg,issa_avg,isr_avg, &
                       isd550,iref_avg,naot550_avg,icld_frac,alg_flag,conf_flag,iDTDBaot,iDTDBqa,iDTDBflag, &
                       iest_uncert)
     return
     end

! Below subroutine is modified from Shana's code (via email) - 07 may 2010 -------------------------------             
subroutine get_bits(nflag,n8,n7,n6,n5,n4,n3,n2,n1)
    
    implicit none
    integer num,nflag,n1,n2,n3,n4,n5,n6,n7,n8
    integer ncheck,nbit ! added by CES when compilation failed

    if(nflag.lt.0)nflag=nflag+128
  
    num=nflag
    n1=0
    n2=0
    n3=0
    n4=0
    n5=0
    n6=0
    n7=0
    n8=0
    
    if (nflag .ge. 0 .and. nflag.le.255 ) then
      ncheck=128
    if(num.ge.ncheck) then
      nbit=1
      num=num-ncheck
    else
      nbit=0
    endif
    
    n1=nbit
    ncheck=64
    if(num.ge.ncheck) then
      nbit=1
      num=num-ncheck
    else
      nbit=0
    endif
    
    n2=nbit
    ncheck=32
    if(num.ge.ncheck) then
      nbit=1
      num=num-ncheck
    else
      nbit=0
    endif
    
    n3=nbit
    ncheck=16
    if(num.ge.ncheck) then
      nbit=1
      num=num-ncheck
    else
      nbit=0
    endif
      
    n4=nbit
    ncheck=8
    if(num.ge.ncheck) then
      nbit=1
      num=num-ncheck
    else
      nbit=0
    endif
     
    n5=nbit
    ncheck=4
    if(num.ge.ncheck) then
      nbit=1
      num=num-ncheck
    else
      nbit=0
    endif
    
    n6=nbit
    ncheck=2
    if(num.ge.ncheck) then
      nbit=1
      num=num-ncheck
    else
      nbit=0
    endif
      n7=nbit
      n8=num
    endif
      
    return
    end 
