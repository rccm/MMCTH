!-----------------------------------------------------------------------
! !F90
!
! !DESCRIPTION:
!
!    DeepBlue.f90 is the main driver for the MOD_PR04_DB process.
!    MOD_PR04_DB reads in MOD02, MOD03, and MOD06 data, processes
!    the data using the Deep Blue aerosol algorithm, and writes
!    the Deep Blue data fields into the MOD04 product.
!
! !INPUT PARAMETERS:  none
!
! !OUTPUT PARAMETERS:  none
!
! !REVISION HISTORY:
!
!    Initial Version by Jeremy Warner   12/01/2006
!    Updated by Clare Salustro     07/13/2009
!    Last Updated by Myeong-Jae Jeong (MJ)  10/26/2010
!
! !TEAM-UNIQUE HEADER:
!
!    This software is developed by the Deep Blue Science Team
!    for the National Aeronautics and Space Administration,
!    Goddard Space Flight Center, under contract NAS5-02041.
!
! !REFERENCES AND CREDITS
!
!    Modelled on code originally developed by Mark Gray and 
!    modified by J. Wei.
!
! !DESIGN NOTES:
!
!
!   Externals:
!
!        MODIS_F_GENERIC            (MODIS_39500.f)
!        MODIS_S_GENERIC            (MODIS_39500.f)
!        MODIS_S_SUCCESS            (MODIS_39500.f)
!        SR_band1                   (core_arrays.f90)
!        SR_band2                   (core_arrays.f90)
!        SR_band3                   (core_arrays.f90)
!        SR_band4                   (core_arrays.f90)
!        SR_band5                   (core_arrays.f90)
!        SR_band6                   (core_arrays.f90)
!        SR_band7                   (core_arrays.f90)
!        landsea                    (core_arrays.f90)
!        latitude                   (core_arrays.f90)
!        longitude                  (core_arrays.f90)
!        solar_zenith_angle         (core_arrays.f90)
!        sensor_zenith_angle        (core_arrays.f90)
!        relative_azimuth_angle     (core_arrays.f90)
!        sensor_azimuth_angle       (core_arrays.f90)
!        solar_azimuth_angle     		(core_arrays.f90)
!
!     Functions and Subroutines:
!
!        get_tables
!        allocate_arrays
!        get_modis_cube_data
!        GetCloudMask
!        getQAFlags
!        dump_data
!        modis
!        MODIS_SMF_SETDYNAMICMSG
!
! !END
!-----------------------------------------------------------------------
program DeepBlue

      use GeneralAuxType
      use core_arrays
      use primary_interface
      use lut_arrays   
      
      use landcover, only: load_landcover, unload_landcover
      use seawifs_surface_pressure, only: load_surfterr_table
      use modis_surface, only: 														&
      													load_brdf, 								&
      													unload_brdf, 							&
      													load_terrainflg_tables,		&
      													load_seasonal_desert, 		&
      													set_limits, 							&
      													load_hdfLER, 							&
      													get_LER412, 							&
      													get_LER470, 							&
      													get_LER650, 							&
      													get_LER865,               &
      													terrain_flag_new
      
!      use db_debug
      												      			      
      implicit none

      include 'PGS_MODIS_39500.f'
      INCLUDE 'PGS_SMF.f'
      INCLUDE 'PGS_PC.f'
      INCLUDE 'PGS_PC_9.f'
      INCLUDE 'PGS_IO.f'
      INCLUDE 'PGS_IO_1.f'
      
      include 'newaottbl90.inc'
      include 'sfc21tbl90.inc' ! @@@

      INTEGER        ierr
      INTEGER        No_byte,   Fmax,     Lmax,     Fbmax,      Lbmax
      PARAMETER     (No_byte=6, Fmax=150, Lmax=550, Fbmax=1500, Lbmax=5500)
      INTEGER        No_byte_O
      PARAMETER     (No_byte_O=5)
      CHARACTER*(*)  FUNCNAME
      PARAMETER     (FUNCNAME = 'DeepBlue')

      integer(kind=1) :: QAFlags(No_byte,Fmax,Lmax), QAFlags_O(No_byte_O,Fmax,Lmax)
      integer(kind=2) :: DTaot(Fmax,Lmax)
      integer(kind=4) :: getQAFlags, dims(3), get_1km_dims, dims_1km(3)
      integer(kind=4) :: GetCloudMask, LandSeaFlag(Fbmax,Lbmax), lsf
!      real(kind=4)    :: gref412(3600,1800),gref470(3600,1800), &
!      real(kind=4)		:: gref650(3600,1800)
!      real(kind=4)    :: ps(360,720)
      real(kind=4)    :: dd, rr, rr2, bt11, btd8, btd11, Dstar1 !changed for JH test 6/7/11
      real(kind=4)    :: sfc_refl_412, sfc_refl_470, sfc_refl_650
      real(kind=4)    :: delta
      parameter          (delta = -0.000001)

! --- include the HDF FORTRAN parameters
      real(kind=4)    :: tmp(12), tmpvg(6)  ! tmpvg added by MJ@@@ 07/26/11

! --- HDF integer functions (must be declared)
      integer, dimension (2)    :: start, edge, stride
      integer, dimension (2)    :: localstart, localedge, localstride

! --- HDF variables
      integer(kind=4)           :: status, i, number_of_bands, checkvariable, checkvariable26,checkvariable03
      integer(kind=4)           :: xdim, ydim, j, ilat5, ilon5
      integer(kind=4)           :: ii, ii2, jj, jj2, ix
      real(kind=4)              :: xmax, xmin, xx, var
      integer, dimension (11)   :: sp_bands

      character(256)  ::  geofiles, l1bfiles, tstfile, l1bfiles_base, polcorfile, lc_file, elev_file
      character(len=256)  ::  file_name
      character(256)  ::  polcorfile_aq, igbplcfile, nv21file
      character(256)  ::  msg, HDFAttrName,ECSParmName
      real(kind=4)    ::  cl_flag

      INTEGER(kind=4) ::  file_version, prtn, rtn, Vrsn 
      INTEGER(kind=4) ::  LRN_CLDMSK,        LRN_L1B_1km,    LRN_L1B_1km_RA,        LRN_Geo
      PARAMETER          (LRN_CLDMSK=422500, LRN_L1B_1km=700002, LRN_L1B_1km_RA=430001, LRN_Geo=600000)

      INTEGER(kind=4) ::  LRN_LCOVER,        LRN_ELEV,          LRN_650BRDF
      PARAMETER          (LRN_LCOVER=412307, LRN_ELEV=412001,   LRN_650BRDF=412400)
      
      INTEGER(kind=4) ::  LRN_GZFLG,         LRN_SDSRT
      PARAMETER          (LRN_GZFLG=412007,  LRN_SDSRT=412006)

    
      INTEGER(kind=4) ::  pgs_io_gen_closef,   pgs_io_gen_openf, &
                          pgs_pc_getreference, pgs_met_getpcattr_s, add_db_mod04, pgs_met_getpcattr

      integer         ::  reclen, fh, fv, fhcorr, EndDate
      LOGICAL         ::  modis_flag = .true., error_flag = .true. 
      REAL(kind=4)    ::  MinSolarZenithAngle, SolarZenithAngleZEPS
      CHARACTER(8)    ::  InstrumentMode, char_buf, Platform
      CHARACTER(512)  ::  usrlog
      PARAMETER          (SolarZenithAngleZEPS = 84.000001)


! --- POLCOR
! --- HDF variables
			integer			  		:: ipix, iframe, iscan, l1bfiles_len
			integer		        :: mside, it1, it2   
			real(kind=4)    	:: band26

! -- Local variables
			real(double)		       				 :: alpha  
			real(double),dimension(6)      :: xrvsc1, xrvsc2   
			real(double),dimension(6)      :: xm12c1, xm12c2   
			real(double),dimension(6)      :: xm13c1, xm13c2   
			real(double),dimension(3)      :: xam12, xam13, xrvs  ! m12, m13 
			integer                        :: iyear, ijday, ibset 
			real(double)                   :: sec_inp
			real(double)                   :: wt, pixnum, ap, ap2, ap3, ap4, ap5, x1, x2
			real(single)                   :: xlat,xlong,sza,xthet,xphi,sazim
			real(single)                   :: vazim, df_azi
			real(single), dimension(3)     :: xnvalm5, lx_m
			integer                        :: idim, jdim
			integer			       						 :: idet
			real(single), dimension(3)     :: polcor=1    !initialized with 1
      real(single), allocatable, dimension (:,:)  :: polcor2d, polcor2d_2
			real(single), allocatable, dimension (:,:)  :: xrvs_test2d, xam12_test2d
			real(single), allocatable, dimension (:,:)  :: xam13_test2d, alpha_test2d
			real(single), allocatable, dimension (:,:)  :: lr_q_test2d, lr_u_test2d, &
																						 lr_x_test2d
			real(single), allocatable, dimension (:,:)  :: xrvs_test2d_2, xam12_test2d_2, &
																						 xam13_test2d_2
			real(single), allocatable, dimension (:,:)  :: lr_q_test2d_2, lr_u_test2d_2, &
																						 lr_x_test2d_2
			real(single)		       							:: cossza 
			integer                        :: icof
			real(double),dimension(6)      :: xrvsc1_2, xrvsc2_2  ! for B03 
			real(double),dimension(6)      :: xm12c1_2, xm12c2_2  ! (using
			real(double),dimension(6)      :: xm13c1_2, xm13c2_2  ! B10 data)
			real(double),dimension(6)      :: xrvsc1_3, xrvsc2_3  ! for B01 
			real(double),dimension(6)      :: xm12c1_3, xm12c2_3  ! (using
			real(double),dimension(6)      :: xm13c1_3, xm13c2_3  ! B13 data)
			real(single), allocatable, dimension (:,:)  :: rfl_case1B08, rfl_case2B08, &
			                                       rfl_case3B08, rfl_case4B08
			real(single), allocatable, dimension (:,:)  :: rfl_case1B03, rfl_case2B03, &
			                                       rfl_case3B03, rfl_case4B03

      integer                         ::  check_subset
      integer, dimension(2)           ::  pixel_offsets
      logical                         ::  is_subset

      integer(kind=4)	:: LUN_LUT,LUN_NV21,LUN_LUT_aq, LUN_IGBPLC
      parameter(LUN_LUT=412301, LUN_LUT_aq=412305, LUN_IGBPLC=412304)                   !@MJ@
      parameter(LUN_NV21=412306)  ! MJ 2011/08/11
			
			real(kind=4), dimension(:,:), allocatable   ::  stdv02, stdv26, stdv03
			integer(kind=4)                             ::  m0, m1, n0, n1
		  integer(kind=4)                             ::  m, n
		  integer(kind=4)                             ::  cnt
		  real(kind=4)                                ::  mean, m2, sdelta
		  real(kind=4)                                ::  mean26, m26, sdelta26
		  real(kind=4)                                ::  mean03, m03, sdelta03
		  real(kind=4)                                ::  sfctmp,ugrd,vgrd,pwat,ozone
		  real(kind=4)                                ::  t1,t2,t3,t4
		  real(kind=4)                                ::  ndvi, sca, cc, psi
		  real(kind=4)                                ::  ndsi
		  
		  integer, dimension(2)                       ::  dims2
		  
! -- RT data for desert boundary condition on I,Q, U,V
  		real(single), dimension(3)     :: Lr_q, Lr_u, L_x, Lr_x   
  		real(single), dimension(3)     :: Lr_xtm  

!..... polcor new addition begins (Oct 2010)                                                              !@MJ@
            real(single)                     :: xrvs_tmp1, xrvs_tmp2, wtz, zaoi                           !@MJ@
	    integer(kind=4)		     :: nscan, nswath, iaoi1, iaoi2          !@MJ@
  	    real(double),dimension(6) :: xrvsc1_2a,xrvsc2_2a,xm12c1_2a,xm12c2_2a,xm13c1_2a,xm13c2_2a !B9  !@MJ@
  	    real(double),dimension(6) :: xrvsc1_2b,xrvsc2_2b,xm12c1_2b,xm12c2_2b,xm13c1_2b,xm13c2_2b !B10 !@MJ@
            real(single), dimension (5)      :: aoiref
            real(double), dimension (10,5,2) :: b8ym12ref, b8ym13ref
            real(double), dimension (10,5,2) :: b3ym12ref, b3ym13ref
            real(double), dimension (10,5,2) :: b3om12ref, b3om13ref
            common /prepol/ b8ym12ref,b8ym13ref,b3ym12ref,b3ym13ref,b3om12ref,b3om13ref
!..... polcor new addition ends (Oct 2010)                                                                !@MJ@
      
!      real, dimension(:,:), allocatable       ::  dstar_arr
      
! --- END POLCOR

            real(single)                     :: sza_in, vza_in, pwv_in, to3_in   !@MJ@ 08/04/2011
            real(double)                     :: srb1,srb2,srb4,srb9,srb10        !@MJ@ 08/04/2011
      
      integer                               ::  iidx, jidx


      integer tsti, tstj
      real    tstLER
      integer ::  gzflg


      integer  num_args
      integer  FlagRA
      character FlagBuff*10
      integer  iargc


      num_args = iargc ( )
      
      if(num_args .eq. 1) then
         call getarg ( 1, FlagBuff )
         read (FlagBuff,*) FlagRA
      else
         !     This is the default value
         FlagRA = 0
      endif
  
      msg = "Beginning Deep Blue processing of MODIS data."
      call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC, msg, 'DeepBlue')

      status = 0


!-------------------------------------------------------------------------------
!
! Retrieve Product Specific Attribute (PSA), MinSolarZenithAngle, from CloudMask 
! product.  If unable to retrieve MinSolarZenithAngle, exit fatal error 
!
!-------------------------------------------------------------------------------

      HDFAttrName = 'CoreMetadata.0'
      ECSParmName = 'MinSolarZenithAngle'
      Vrsn        = 1

      rtn=pgs_met_getpcattr_s(LRN_CldMsk,Vrsn,HDFAttrName,ECSParmName,char_buf)

      If (rtn .NE. PGS_S_SUCCESS) Then
         error_flag = .true.

         usrlog =                                                                            &
         'Halt process! pgs_met_getpcattr_s unable to retrieve MinSolarZenithAngle '         &
         // char(10) // 'from Cloud Mask product.  MOD_PR04DB exiting fail code 1. '         &
         // char(10) // 'Operator Action: Notify SDST.  Also check that MOD03 product '      &
         // char(10) // 'contains valid geolocation data including the solar zenith angle. ' 

         Call MODIS_SMF_SetDynamicMSG(MODIS_F_GENERIC,usrlog,FUNCNAME)
         Call Exit(0)

      End If

      read(char_buf,'(f8.2)') MinSolarZenithAngle

!-------------------------------------------------------------------------------
!
! Retrieve DayNightFlag from Geolocation product.  Exit fatal error if unable to
! to retrieve DayNightFlag.  
!
!-------------------------------------------------------------------------------

      HDFAttrName = 'CoreMetadata.0'
      ECSParmName = 'DAYNIGHTFLAG'
      Vrsn        = 1

      rtn=pgs_met_getpcattr_s(LRN_GEO,Vrsn,HDFAttrName,ECSParmName,char_buf)

      If (rtn .NE. PGS_S_SUCCESS) Then
         error_flag = .true.

         usrlog =                                                                              &
         'Halt process! pgs_met_getpcattr_s unable to retrieve DayNightFlag '                  &
         // char(10) // 'from Geolocation (MOD03) product.  MOD_PR04DB exiting fail code 1. '  &
         // char(10) // 'Operator Action: Notify SDST.  Also check for possibly '              &
         // char(10) // 'corrupt Geolocation product.'

         Call MODIS_SMF_SetDynamicMSG(MODIS_F_GENERIC,usrlog,FUNCNAME)
         Call Exit(0)
     End If

      InstrumentMode = char_buf
      modis_flag     = .true.


    
!-------------------------------------------------------------------------------
!
! Set logic to process the Aerosol Product 
! 
! Compare Cloud Mask MinSolarZenithAngle against Dark Target threshold of 84 
! degrees.  If MinSolarZenithAngle exceeds 84, skip granule.
!
! NOTE: We do our own SZA checks later at 72 degrees. Due to the combined DT/DB 
!     product, we have to "process" granules up to 84 degrees to any DT data will
!     be placed in the combined SDS.  
!
!-------------------------------------------------------------------------------

      If (MinSolarZenithAngle .GE. SolarZenithAngleZEPS) Then

         usrlog = 'Minimum Solar zenith angle is greater than 84 degrees. '  &
         // char(10) // 'Exit DeepBlue Algorithm '                        

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_M_GENERIC,usrlog,FUNCNAME)
         Call Exit(0)


      Else If (MinSolarZenithAngle .GE. 0.0) Then
! --- Minumum solar zenith angle is between 0 and 72 degrees

         If (InstrumentMode .EQ. 'Night') Then
! ------ instrument in night mode over illuminated earth - process as night

            usrlog = 'The earth scene is illuminated, but MODIS instrument is in night mode.' &
            // char(10) // 'Deep Blue Algorithm not run'                    

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_M_GENERIC,usrlog,FUNCNAME)
            Call Exit(0)

         Else If ( (InstrumentMode .EQ. 'Both') .OR. (InstrumentMode .EQ. 'Day') ) Then
! ------ process as day

            usrlog = 'Earth scene is illuminated and MODIS instrument is in '  &
            // char(10) // 'Day or Both mode.  Run Deep Blue Algorithm.'                    

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_M_GENERIC,usrlog,FUNCNAME)

         End If

      Else 
! --- Oh brother, minumum solar zenith angle is less than zero. Exit fail

         usrlog = 'Oh brother, minumum solar zenith angle is less than zero.  PGE04 '      &
         // char(10) // 'halting, no retrieval performed.  Process exiting fail, code 1. ' & 
         // char(10) // 'Operator Action: Notify SDST.  Also check for possibly '          &
         // char(10) // 'corrupt Geolocation and MOD35_L2 product.' 

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,usrlog,FUNCNAME)
         Call Exit(0)

      End If

!-------------------------------------------------------------------------------
!
! Retrieve ASSOCIATEDPLATFORMSHORTNAME from Geolocation product.  Exit fatal  
! error if unable to retrieve ASSOCIATEDPLATFORMSHORTNAME.  
!
!-------------------------------------------------------------------------------

      HDFAttrName = 'CoreMetadata.0'
      ECSParmName = 'ASSOCIATEDPLATFORMSHORTNAME.1'
      Vrsn        = 1

      rtn=pgs_met_getpcattr_s(LRN_GEO,Vrsn,HDFAttrName,ECSParmName,char_buf)
      
      If (rtn .NE. PGS_S_SUCCESS) Then
         error_flag = .true.

         usrlog =                                                                              &
         'Halt process! pgs_met_getpcattr_s unable to retrieve ASSOCIATEDPLATFORMSHORTNAME '   &
         // char(10) // 'from Geolocation (MOD03) product.  MOD_PR04DB exiting fail code 1. '  &
         // char(10) // 'Operator Action: Notify SDST.  Also check for possibly '              &
         // char(10) // 'corrupt Geolocation product.'

         Call MODIS_SMF_SetDynamicMSG(MODIS_F_GENERIC,usrlog,FUNCNAME)
         Call Exit(0)
      End If

      Platform = char_buf
      modis_flag     = .true.
! TESTING --- To turn off polcor for Terra:-------------------------------------
!	    Platform = 'TEST'
! TESTING ----------------------------------------------------------------------
			print *, "Platform is: ",Platform

!      If (Platform .EQ. 'Terra') Then
!! ------ MODIS-Terra - Need to check date, must be before 2008
!!     		print *, "* TERRA - Will need to perform polarization correction"

!				HDFAttrName = 'CoreMetadata.0'
!				ECSParmName = 'RANGEENDINGDATE'
!				Vrsn        = 1
	
!				rtn=pgs_met_getpcattr_s(LRN_GEO,Vrsn,HDFAttrName,ECSParmName,char_buf)
!				read(char_buf,'(I4)') EndDate

!!			If (EndDate .ge. 2008) Then  ! ori. for C051 Terra/MODIS
!				If (EndDate .gt. 2011) Then
!!				print *, '* Acquisition date after Dec 31, 2007! Deep Blue cannot proceed for Terra!'
!					print *, '* Acquisition date after Dec 31, 2010! Deep Blue cannot proceed for Terra!'
!					
!					usrlog =                                                                        &
!					              'Acquisition date is after Dec 31, 2010. Cannot proceeed with '   &
!					// char(10) //'Terra processing because RVS/polarization correction info is '   &
!					// char(10) //'not available. No Deep Blue SDSs will be filled in MOD04_L2. '   &
!					// char(10) //'Ending Deep Blue program.'
!					
!					Call MODIS_SMF_SetDynamicMSG(MODIS_M_GENERIC,usrlog,FUNCNAME)
!          Call Exit(0)
!          
!				ELSE				
!! ------ MODIS-Terra - Need to perform polarization correction
!
!          usrlog = 'Working with Terra Platform. Need to perform polarization '&
!          // char(10) // 'corection.'
!          
!          Call MODIS_SMF_SETDYNAMICMSG(MODIS_M_GENERIC,usrlog,FUNCNAME)
!				
!				End If
!				
!      Else If (Platform .EQ. 'Aqua') Then
!! ------ MODIS-Aqua - Now need to perform polarization correction
!!   	print *, "* AQUA - Will need to perform polarization correction"
!			HDFAttrName = 'CoreMetadata.0'
!			ECSParmName = 'RANGEENDINGDATE'
!			Vrsn        = 1
!			rtn=pgs_met_getpcattr_s(LRN_GEO,Vrsn,HDFAttrName,ECSParmName,char_buf)
!			read(char_buf,'(I4)') EndDate
!
!			        If (EndDate .gt. 2011) Then
!  		   	           print *, '* Acquisition date after Dec 31, 2010! Deep Blue cannot proceed for Aqua!'
!			
!					usrlog =                                                                        &
!					              'Acquisition date is after Dec 31, 2010. Cannot proceeed with '   &
!					// char(10) //'Aqua processing because RVS correction information is '   &
!					// char(10) //'not available. No Deep Blue SDSs will be filled in MYD04_L2. '   &
!					// char(10) //'Ending Deep Blue program.'
!					
!					Call MODIS_SMF_SetDynamicMSG(MODIS_M_GENERIC,usrlog,FUNCNAME)
!          Call Exit(0)
!          
!				ELSE				
!
!!@MJ@      usrlog = 'Working with Aqua Platform. No polarization correction '&
!!@MJ@     // char(10) // 'required.'
!           usrlog = 'Working with Aqua Platform. Need to perform polarization '&
!          // char(10) // 'correction.'
!
!          Call MODIS_SMF_SETDYNAMICMSG(MODIS_M_GENERIC,usrlog,FUNCNAME)
!
!				End If
!
!      End If  ! ... Else If (Platform .EQ. 'Aqua') Then
      
      status = get_1km_dims(dims_1km)
!      print *, 'dims_1km: ', dims_1km
      
!     Read in Flags from MOD04 and MOD35 products

      ierr = getQAFlags(QAflags, QAflags_O, DTaot, dims)
      nscan = dims(2)
      nswath = dims(1)
!      print *, 'nscan,nswath from getQAFlags: ',nscan,nswath
      
			if (nscan .ne. 203) then
         msg = "WARNING: Number of scan lines not equal to 203. Program will continue."
!         print *, msg
         call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC, msg, 'DeepBlue')
			endif

			if (nswath .ne. 135) then
         msg = "WARNING: Number of swaths not equal to 135. Program will exit."
!				 print *, msg
				 call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC, msg, 'DeepBlue')
			endif	
			
!     Read in tables

      call get_tables(fh)
      call get_newtables
!      call get_surf_tables
      
!     Get MOD02 and MOD03 filenames from PCF

      file_version = 1
      if( FlagRA .eq. 1) then
         prtn=pgs_pc_getreference(LRN_L1B_1km_RA,file_version,l1bfiles)
      else
         prtn=pgs_pc_getreference(LRN_L1B_1km,file_version,l1bfiles)
      endif
      if (prtn .lt. 0) then
         msg = "Error retrieving LRN_L1B_1km lun from pcf."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      endif
      msg = "Will use MOD02 L1B file "//l1bfiles
      call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC, msg, 'DeepBlue')

!			This doesn't always work - not sure why!
			l1bfiles_len = len_trim(l1bfiles)
      l1bfiles_base = l1bfiles(index(l1bfiles,"/",.true.)+1:l1bfiles_len)
      print *, trim(l1bfiles_base)//'.txt'
      
      file_version = 1
      prtn=pgs_pc_getreference(LRN_Geo,file_version,geofiles)
      if (prtn .lt. 0) then
         msg = "Error retrieving LRN_Geo lun from pcf."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      endif
      msg = "Will use MOD03 file "//geofiles
      call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC, msg, 'DeepBlue')
      
      is_subset = .false.
			pixel_offsets = (/-999, -999/)
			status = check_subset(geofiles, l1bfiles, is_subset, pixel_offsets)
			if (status /= 0) then
			  msg = 'ERROR: Could not determine if file is subset.'
			  print *, msg
			  call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
			end if
!   		print *, 'is_subet, offsets: ', is_subset, pixel_offsets


! -- load auxiliary input files
!     -- land cover
      file_version = 1
      prtn=pgs_pc_getreference(LRN_LCOVER,file_version,lc_file)
      if (prtn .lt. 0) then
         msg = "Error retrieving LRN_LCOVER lun from pcf."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      endif
      msg = "Will use land cover file "//lc_file
      call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC, msg, 'DeepBlue')
      
      status = load_landcover(lc_file)
      if (status /= 0) then
        msg = "Error loading land cover input file."
        call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      end if
      
!     -- elevation
      file_version = 1
      prtn=pgs_pc_getreference(LRN_ELEV,file_version,elev_file)
      if (prtn .lt. 0) then
         msg = "Error retrieving LRN_ELEV lun from pcf."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      endif
      msg = "Will use elevation file "//elev_file
      call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC, msg, 'DeepBlue')
      
      status = load_surfterr_table(elev_file)
      if (status /= 0) then
        msg = "Error loading elevation input file."
        call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      end if
      
!     -- geo zones
      file_version = 1
      prtn=pgs_pc_getreference(LRN_GZFLG,file_version,file_name)
      if (prtn .lt. 0) then
         msg = "Error retrieving LRN_GZFLG lun from pcf."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      endif
      msg = "Will use geozone file "//file_name
      call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC, msg, 'DeepBlue')
      
      status = load_terrainflg_tables(file_name)
      if (status /= 0) then
        msg = "Error loading geozone input file."
        call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      end if
      
!     -- seasonal desert data
      file_version = 1
      prtn=pgs_pc_getreference(LRN_SDSRT,file_version,file_name)
      if (prtn .lt. 0) then
         msg = "Error retrieving LRN_SDSRT lun from pcf."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      endif
      msg = "Will use seasonal desert file "//file_name
      call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC, msg, 'DeepBlue')
!      print *, msg
      status = load_seasonal_desert(file_name)
      if (status /= 0) then
        msg = "Error loading seasonal desert input file."
        call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      end if
      
!     -- 670nm BRDF data
      file_version = 1
      prtn=pgs_pc_getreference(LRN_650BRDF,file_version,file_name)
      if (prtn .lt. 0) then
         msg = "Error retrieving LRN_650BRDF lun from pcf."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      endif
      msg = "Will use 650brdf file "//file_name
      call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC, msg, 'DeepBlue')

      status = load_brdf(file_name)
      if (status /= 0) then
        msg = "Error loading elevation input file."
        call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      end if
      
! --- set up input arrays for "sfrdata" (hyperslab data read)
! need to fix this arbitrary input
!================================
        start = (/1,1/)
!        edge  = (/1354,nscan*10/)
        edge	= (/dims_1km(1),dims_1km(2)/)
        stride = (/1,1/)

        sp_bands = (/1,3,8,7,29,31,32,26,2,5,4/)

        xdim = edge(1)
        ydim = edge(2)	
                
!================================
        number_of_bands=size(sp_bands)
        allocate (bands(number_of_bands), stat = checkvariable)
        if ( checkvariable /= 0 ) then
           msg = "Error allocating bands."
           call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
        endif
       
        do i=1, number_of_bands
           bands(i)=sp_bands(i)
!           print *, 'band',i,'=',bands(i) 
        enddo

      localstart=(/start(1)-1, start(2)-1/)
      localedge=edge
      localstride=stride

! setup 
      call allocate_arrays (localedge, status )
      if (status .ne. 0 ) then
         msg = "Error detected in allocate_arrays."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      endif
!      print *, 'successfully passed through allocatearrays routine',status
  
      ! get data cube to be processed and ancillary data arrays
      call get_modis_data_cube(l1bfiles,geofiles,localstart, localedge, &
                                 localstride, sp_bands, status)
      if (status .ne. 0 ) then
         msg = "Error detected in get_modis_data_cube."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      endif
      
      ! put band_measurement into different variable
      allocate (SR_band1(xdim,ydim), stat = checkvariable)
      if ( checkvariable /= 0 ) then
         msg = "Error allocating SR_band1."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      endif
      SR_band1(:,:)=band_measurements(:,1,:)

      allocate (SR_band2(xdim,ydim), stat = checkvariable)
      if ( checkvariable /= 0 ) then
         msg = "Error allocating SR_band2."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      endif
      SR_band2(:,:)=band_measurements(:,2,:)

      allocate (SR_band3(xdim,ydim), stat = checkvariable)
      if ( checkvariable /= 0 ) then
         msg = "Error allocating SR_band3."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      endif
      SR_band3(:,:)=band_measurements(:,3,:)
      
      allocate (SR_band4(xdim,ydim), stat = checkvariable)
      if ( checkvariable /= 0 ) then
         msg = "Error allocating SR_band4."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      endif
      SR_band4(:,:)=band_measurements(:,4,:)

      allocate (SR_band5(xdim,ydim), stat = checkvariable)
      if ( checkvariable /= 0 ) then
         msg = "Error allocating SR_band5."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      endif
      SR_band5(:,:)=band_measurements(:,5,:)

      allocate (SR_band6(xdim,ydim), stat = checkvariable)
      if ( checkvariable /= 0 ) then
         msg = "Error allocating SR_band6."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      endif
      SR_band6(:,:)=band_measurements(:,6,:)

      allocate (SR_band7(xdim,ydim), stat = checkvariable)
      if ( checkvariable /= 0 ) then
         msg = "Error allocating SR_band7."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      endif
      SR_band7(:,:)=band_measurements(:,7,:)

      allocate (SR_band8(xdim,ydim), stat = checkvariable)
      if ( checkvariable /= 0 ) then
         msg = "Error allocating SR_band8."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      endif
      SR_band8(:,:)=band_measurements(:,8,:)

!c... B2 & B5 were added on 21 Jul 2011 for over-vegetation retrieval (MJ)
      allocate (SR_band9(xdim,ydim), stat = checkvariable)
      if ( checkvariable /= 0 ) then
         msg = "Error allocating SR_band9."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      endif
      SR_band9(:,:)=band_measurements(:,9,:)

      allocate (SR_band10(xdim,ydim), stat = checkvariable)
      if ( checkvariable /= 0 ) then
         msg = "Error allocating SR_band10."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      endif
      SR_band10(:,:)=band_measurements(:,10,:)
      
      allocate (SR_band11(xdim,ydim), stat = checkvariable)
      if ( checkvariable /= 0 ) then
         msg = "Error allocating SR_band11."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      endif
      SR_band11(:,:)=band_measurements(:,11,:)

!     -- reset bad 1.38um (band 26) data to small, non-zero value
      where (SR_band5 < 0.0) SR_band5 = 0.00001
       
!      print *,'SR_band7 =',SR_band7 
      
! Allocate for Terra correction --> also need for Aqua
!@MJ@		if (Platform .EQ. 'Terra') then

				allocate (polcor2d(xdim,ydim), stat = checkvariable)
				if ( checkvariable /= 0 ) then
					 msg = "Error allocating polcor2d."
					 call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
				endif

				allocate (polcor2d_2(xdim,ydim), stat = checkvariable)
				if ( checkvariable /= 0 ) then
					 msg = "Error allocating polcor2d_2."
					 call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
				endif
				
				allocate (xrvs_test2d(xdim,ydim), stat = checkvariable)
				if ( checkvariable /= 0 ) then
					 msg = "Error allocating xrvs_test2d."
					 call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
				endif

				allocate (xam12_test2d(xdim,ydim), stat = checkvariable)
				if ( checkvariable /= 0 ) then
					 msg = "Error allocating xam12_test2d."
					 call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
				endif

				allocate (xam13_test2d(xdim,ydim), stat = checkvariable)
				if ( checkvariable /= 0 ) then
					 msg = "Error allocating xam13_test2d."
					 call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
				endif

				allocate (alpha_test2d(xdim,ydim), stat = checkvariable)
				if ( checkvariable /= 0 ) then
					 msg = "Error allocating alpha_test2d."
					 call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
				endif

				allocate (lr_q_test2d(xdim,ydim), stat = checkvariable)
				if ( checkvariable /= 0 ) then
					 msg = "Error allocating lr_q_test2d."
					 call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
				endif

				allocate (lr_u_test2d(xdim,ydim), stat = checkvariable)
				if ( checkvariable /= 0 ) then
					 msg = "Error allocating lr_u_test2d."
					 call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
				endif

				allocate (lr_x_test2d(xdim,ydim), stat = checkvariable)
				if ( checkvariable /= 0 ) then
					 msg = "Error allocating lr_x_test2d."
					 call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
				endif

				allocate (xrvs_test2d_2(xdim,ydim), stat = checkvariable)
				if ( checkvariable /= 0 ) then
					 msg = "Error allocating xrvs_test2d_2."
					 call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
				endif

				allocate (xam12_test2d_2(xdim,ydim), stat = checkvariable)
				if ( checkvariable /= 0 ) then
					 msg = "Error allocating xam12_test2d_2."
					 call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
				endif

				allocate (xam13_test2d_2(xdim,ydim), stat = checkvariable)
				if ( checkvariable /= 0 ) then
					 msg = "Error allocating xam13_test2d_2."
					 call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
				endif

				allocate (lr_q_test2d_2(xdim,ydim), stat = checkvariable)
				if ( checkvariable /= 0 ) then
					 msg = "Error allocating lr_q_test2d_2."
					 call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
				endif

				allocate (lr_u_test2d_2(xdim,ydim), stat = checkvariable)
				if ( checkvariable /= 0 ) then
					 msg = "Error allocating lr_u_test2d_2."
					 call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
				endif

				allocate (lr_x_test2d_2(xdim,ydim), stat = checkvariable)
				if ( checkvariable /= 0 ) then
					 msg = "Error allocating lr_x_test2d_2."
					 call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
				endif

				allocate (rfl_case1B08(xdim,ydim), stat = checkvariable)
				if ( checkvariable /= 0 ) then
					 msg = "Error allocating rfl_case1B08."
					 call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
				endif
        rfl_case1B08(:,:) = -999.0
        
				allocate (rfl_case2B08(xdim,ydim), stat = checkvariable)
				if ( checkvariable /= 0 ) then
					 msg = "Error allocating rfl_case2B08."
					 call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
				endif
        rfl_case2B08(:,:) = -999.0
        
				allocate (rfl_case3B08(xdim,ydim), stat = checkvariable)
				if ( checkvariable /= 0 ) then
					 msg = "Error allocating rfl_case3B08."
					 call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
				endif
				rfl_case3B08(:,:) = -999.0
        
				allocate (rfl_case4B08(xdim,ydim), stat = checkvariable)
				if ( checkvariable /= 0 ) then
					 msg = "Error allocating rfl_case4B08."
					 call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
				endif
				rfl_case4B08(:,:) = -999.0
        
				allocate (rfl_case1B03(xdim,ydim), stat = checkvariable)
				if ( checkvariable /= 0 ) then
					 msg = "Error allocating rfl_case1B03."
					 call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
				endif
        rfl_case1B03(:,:) = -999.0
        
				allocate (rfl_case2B03(xdim,ydim), stat = checkvariable)
				if ( checkvariable /= 0 ) then
					 msg = "Error allocating rfl_case2B03."
					 call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
				endif
        rfl_case2B03(:,:) = -999.0
        
				allocate (rfl_case3B03(xdim,ydim), stat = checkvariable)
				if ( checkvariable /= 0 ) then
					 msg = "Error allocating rfl_case3B03."
					 call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
				endif
        rfl_case3B03(:,:) = -999.0
        
				allocate (rfl_case4B03(xdim,ydim), stat = checkvariable)
				if ( checkvariable /= 0 ) then
					 msg = "Error allocating rfl_case4B03."
					 call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
				endif
				rfl_case4B03(:,:) = -999.0
        
        
!        allocate (dstar_arr(xdim,ydim), stat=status)
!        if (status /= 0) then
!          print *, "Error allocating dstar_arr: ", status
!          stop
!        endif
        
        ! JAW input limited dimension of LER tables based on input lat/long range

!        if (Platform .EQ. 'Aqua') then
          ierr = set_limits(localedge, latitude, longitude)
          if (ierr /= 0) then
             msg = "Error getting lat/long limits."
             call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
          endif

          ierr = load_hdfLER()
          if (ierr /= 0) then
             msg = "Error getting lat/long limits."
             call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
!          endif

          !tsti = floor(10*(90.0 + latitude(1200,1200)))
          !tstj = floor(10*(180.0+longitude(1200,1200)))
          !print *,"tst coords = ",tsti,tstj
          !tstLER = get_LER865(tsti,tstj)
          !print *,"tst 865 = ",tstLER
          !tstLER = get_LER412(tsti,tstj,0.2,127.5,110.0)
          !print *,"tst 412 = ",tstLER
          !tstLER = get_LER412(tsti,tstj,0.2,127.5,60.0)
          !print *,"tst 412 = ",tstLER
          !tstLER = get_LER470(tsti,tstj,0.1,127.5,110.0)
          !print *,"tst 470 = ",tstLER
          !tstLER = get_LER470(tsti,tstj,0.1,127.5,60.0)
          !print *,"tst 470 = ",tstLER
          !tstLER = get_LER650(tsti,tstj,0.4,127.5,110.0)
          !print *,"tst 650 = ",tstLER
          !tstLER = get_LER650(tsti,tstj,0.4,127.5,60.0)
          !print *,"tst 650 = ",tstLER
          !stop
        endif
        ! end JAW mods
			
!     -- calculate stdv of 412nm reflectance within ~0.2 deg radius of each pixel
!     -- using the Welford single-pass algorithm.
!     -- See http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#On-line_algorithm)
      !print *, '---- start stdv02 ------'
      allocate (stdv02(xdim,ydim), stat = checkvariable)
      allocate (stdv26(xdim,ydim), stat = checkvariable26)
      allocate (stdv03(xdim,ydim), stat = checkvariable03)
      if ( checkvariable /= 0 ) then
         msg = "Error allocating stdv02 array."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      endif
     if ( checkvariable26 /= 0 ) then
         msg = "Error allocating stdv26 array."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      endif
     if ( checkvariable03 /= 0 ) then
         msg = "Error allocating stdv03 array."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      endif
      stdv02(:,:) = -999.0
      stdv26(:,:) = -999.0
      stdv03(:,:) = -999.0
      
!      print *, 'granule size: ', xdim, ydim
      do j = 1, ydim
        do i = 1, xdim
          m0 = i - 10
          m1 = i + 10
          n0 = j - 10
          n1 = j + 10
          if (m0 < 1) m0 = 1
          if (m1 > xdim) m1 = xdim
          if (n0 < 1) n0 = 1
          if (n1 > ydim) n1 = ydim
          
          cnt = 0
          mean  = 0.0
          mean26= 0.0
          mean03 = 0.0
          m2    = 0.0
          m26 = 0.0
          m03 = 0.0
          sdelta = -999.0
          sdelta26 = -999.0
          sdelta03 = -999.0
          do n = n0, n1
            do m = m0, m1
              if (SR_band3(m,n) < 0.0)  cycle
              cnt   = cnt + 1
              sdelta = SR_band3(m,n) - mean
              sdelta26 = SR_band8(m,n) - mean26
              sdelta03 = SR_band2(m,n) - mean03
              mean  = mean + sdelta/cnt
              mean26  = mean26 + sdelta26/cnt
              mean03  = mean03 + sdelta03/cnt
              m2    = m2 + sdelta*(SR_band3(m,n) - mean)
              m26    = m26 + sdelta26*(SR_band8(m,n) - mean26)
              m03    = m03 + sdelta03*(SR_band2(m,n) - mean03)
            end do
          end do
          stdv02(i,j) = sqrt(m2/(cnt-1))
          stdv26(i,j) = sqrt(m26/(cnt-1))
          stdv03(i,j) = sqrt(m03/(cnt-1))
        end do
      end do
      !print *, '---- end stdv02 ------'
      
!     Repackage the output so it can be handled by the deepblue program
!      do 12 j = 1, 2030
!        do 10 i = 1, 1354
      
! -- Read polarization sensitivity parameters for Terra/MODIS
! The parameters would better be read before accessing individual
! MODIS granules.
!       read parameters (m12, m13, year, doy, ...) from the hdf file
      file_version = 1
      prtn=pgs_pc_getreference(LUN_LUT,file_version,polcorfile)
      if (prtn .lt. 0) then
         msg = "Error retrieving POLCOR lun from pcf."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      endif
      msg = "Will use POLCOR file "//polcorfile
      call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC, msg, 'DeepBlue')

        status = 0
        if (status == 0 ) then
           call get_lut_xcal_modist_ts(polcorfile, status)  ! time series (max 117)
!           write(*,*)'get_lut_xcal_modist_ts status       >',status
        else
           write(*,*)'Failure detected before get_lut_xcal_modist_ts' 
        endif

!@MJ --below added for Aqua/MODIS RVS/polcor parameters
      file_version = 1
      prtn=pgs_pc_getreference(LUN_LUT_aq,file_version,polcorfile_aq)
      if (prtn .lt. 0) then
         msg = "Error retrieving POLCOR lun from pcf."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      endif
      msg = "Will use POLCOR file "//polcorfile_aq
      call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC, msg, 'DeepBlue')

        status = 0
        if (status == 0 ) then
           call get_lut_xcal_modisa_ts(polcorfile_aq, status)  ! time series (max 95)
!           write(*,*)'get_lut_xcal_modisa_ts status       >',status
        else
           write(*,*)'Failure detected before get_lut_xcal_modisa_ts' 
        endif
!@MJ --above added for Aqua/MODIS RVS/polcor parameters


!c...  Get pre-launch polarization parameters for Aqua/MODIS B8 !@MJ@
!c...  and B3, and Terra/MODIS B3
       call get_mxd_prepoltbls                                  !@MJ@
      
! -- Read IGBP Land Cover parameters
      file_version = 1
      prtn=pgs_pc_getreference(LUN_IGBPLC,file_version,igbplcfile)
      if (prtn .lt. 0) then
         msg = "Error retrieving IGBP lun from pcf."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      endif
      msg = "Will use IGBP file "//igbplcfile
      call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC, msg, 'DeepBlue')

      status = 0
      if (status == 0 ) then
         call get_lut_igbp_land_cover(igbplcfile, status)
      else
         write(*,*)'Failure detected before get_lut_igbp_land_cover' 
      endif

!c...  below added Aug 11, 2011 
      file_version = 1
      prtn=pgs_pc_getreference(LUN_NV21,file_version,nv21file)
      if (prtn .lt. 0) then
         msg = "Error retrieving NVAL21 lun from pcf."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      endif
      msg = "Will use NVALX21 file "//nv21file
      call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC, msg, 'DeepBlue')

        status = 0
        if (status == 0 ) then
           call get_lut_211sfc(nv21file, status)
        else
           write(*,*)'Failure detected before get_lut_211sfc' 
        endif
!c...  above added Aug 11, 2011 
      
! A SDS "EV center time" (1D array; size=203) in MOD03 hdf file(s) is added to read.
! This SDS will be used to interpolate m12 & m13 from the LUT (polarization 
! params from xcal file). It contains EV obs. time since 1993/01/01 in second unit.
!================================
!  now calculate alpha value
!  we also need to read in values for Lr_q, Lr_u, and L_x
!  Initialize Lr_q, Lr_u, and L_x (measured radiance)  
     Lr_q(:)=0.0   ! Q for Rayleigh atm.
     Lr_u(:)=0.0   ! U for Rayleigh atm.
     L_x(:)=0.0    ! I, measured radiance
  
!-------------------------------------------------------------------------------
!			Inserting new Polarization correction
!-------------------------------------------------------------------------------
!     call open_data(fh)
!     print *, 'This is a test for C006 RVS/polcor @@@@@###@@@@@' ! mj added
!     open(5000, file="dstar.txt", form="formatted", status="replace")
     do 13 iscan=1, nscan   ! original
!      print *, 'iscan, mirror_side(iscan) ', iscan, mirror_side(iscan)
!@MJ@  if (Platform .EQ. 'Terra') then
     	   if (mirror_side(iscan).lt.0 .or. mirror_side(iscan).gt.1) then 
     	 	   print *, 'Bad value: iscan, mirror_side(iscan): ', iscan, mirror_side(iscan)
     	 		 goto 13
         endif
!@MJ@  endif ! check mirror_side for Terra only --> now also for Aqua 
       
!--> interpolation of m12 & m13 (from xcal~) starts here!! 
!... -----------------------------------------------------------
        sec_inp=ev_cntr_time(iscan)                        
        if (Platform .EQ. 'Terra') then                            !@MJ@
           call set_indx_xcaltab(sectab, sec_inp, it1, it2)   
        elseif (Platform .EQ. 'Aqua') then                         !@MJ@
           call set_indx_xcaltab_aq(sectab_aq, sec_inp, it1, it2)  !@MJ@
        endif                                                      !@MJ@
!... -----------------------------------------------------------
        do 12 idet=1,10

      	   ipix=idet+(iscan-1)*10
					 j = ipix
           wt = -999.0
			 		 if (Platform .EQ. 'Terra') then

						 do icof=1,6  ! @chg@ --> till the end of the loop
								ibset=1 ! Band08(in the Pol.Corr. table from Bryan Franz) 
								xrvsc1(icof)=xxrvs(it1,ibset,mirror_side(iscan)+1,idet,icof)
								xrvsc2(icof)=xxrvs(it2,ibset,mirror_side(iscan)+1,idet,icof)
								xm12c1(icof)=xxam12(it1,ibset,mirror_side(iscan)+1,idet,icof)
								xm12c2(icof)=xxam12(it2,ibset,mirror_side(iscan)+1,idet,icof)
								xm13c1(icof)=xxam13(it1,ibset,mirror_side(iscan)+1,idet,icof)
								xm13c2(icof)=xxam13(it2,ibset,mirror_side(iscan)+1,idet,icof)
	
!...	          --> get B03 and B10 parameters to correct for B03 RVS !@MJ@
								ibset=2 ! Band09(in the Pol.Corr. table from Bryan Franz)
								xrvsc1_2a(icof)=xxrvs(it1,ibset,mirror_side(iscan)+1,idet,icof)
								xrvsc2_2a(icof)=xxrvs(it2,ibset,mirror_side(iscan)+1,idet,icof)
								xm12c1_2a(icof)=xxam12(it1,ibset,mirror_side(iscan)+1,idet,icof)
								xm12c2_2a(icof)=xxam12(it2,ibset,mirror_side(iscan)+1,idet,icof)
								xm13c1_2a(icof)=xxam13(it1,ibset,mirror_side(iscan)+1,idet,icof)
								xm13c2_2a(icof)=xxam13(it2,ibset,mirror_side(iscan)+1,idet,icof)

                                ibset=3 ! Band10(in the Pol.Corr. table from Bryan Franz)                                
								xrvsc1_2b(icof)=xxrvs(it1,ibset,mirror_side(iscan)+1,idet,icof)
								xrvsc2_2b(icof)=xxrvs(it2,ibset,mirror_side(iscan)+1,idet,icof)
								xm12c1_2b(icof)=xxam12(it1,ibset,mirror_side(iscan)+1,idet,icof)
								xm12c2_2b(icof)=xxam12(it2,ibset,mirror_side(iscan)+1,idet,icof)
								xm13c1_2b(icof)=xxam13(it1,ibset,mirror_side(iscan)+1,idet,icof)
								xm13c2_2b(icof)=xxam13(it2,ibset,mirror_side(iscan)+1,idet,icof)


!...	          --> below to be aplided to B01 just for a test 
!... 	            (not yet to be applied as of Jan 29, 2009)
								ibset=6 ! Band13(in the Pol.Corr. table from Bryan Franz)
								xrvsc1_3(icof)=xxrvs(it1,ibset,mirror_side(iscan)+1,idet,icof)
								xrvsc2_3(icof)=xxrvs(it2,ibset,mirror_side(iscan)+1,idet,icof)
								xm12c1_3(icof)=xxam12(it1,ibset,mirror_side(iscan)+1,idet,icof)
								xm12c2_3(icof)=xxam12(it2,ibset,mirror_side(iscan)+1,idet,icof)
								xm13c1_3(icof)=xxam13(it1,ibset,mirror_side(iscan)+1,idet,icof)
								xm13c2_3(icof)=xxam13(it2,ibset,mirror_side(iscan)+1,idet,icof)  
						 enddo  !  do icof=1,4

						 wt = (sec_inp-sectab(it1))/(sectab(it2)-sectab(it1))

!						 xrvsc1(1)=xxrvs(it1,ibset,mirror_side(iscan)+1,idet,1)
!						 xrvsc1(2)=xxrvs(it1,ibset,mirror_side(iscan)+1,idet,2)
!						 xrvsc1(3)=xxrvs(it1,ibset,mirror_side(iscan)+1,idet,3)
!						 xrvsc1(4)=xxrvs(it1,ibset,mirror_side(iscan)+1,idet,4)
!						 xm12c1(1)=xxam12(it1,ibset,mirror_side(iscan)+1,idet,1)
!						 xm12c1(2)=xxam12(it1,ibset,mirror_side(iscan)+1,idet,2)
!						 xm12c1(3)=xxam12(it1,ibset,mirror_side(iscan)+1,idet,3)
!						 xm12c1(4)=xxam12(it1,ibset,mirror_side(iscan)+1,idet,4)
!						 xm13c1(1)=xxam13(it1,ibset,mirror_side(iscan)+1,idet,1)
!						 xm13c1(2)=xxam13(it1,ibset,mirror_side(iscan)+1,idet,2)
!						 xm13c1(3)=xxam13(it1,ibset,mirror_side(iscan)+1,idet,3)
!						 xm13c1(4)=xxam13(it1,ibset,mirror_side(iscan)+1,idet,4)
!	
!						 xrvsc2(1)=xxrvs(it2,ibset,mirror_side(iscan)+1,idet,1)
!						 xrvsc2(2)=xxrvs(it2,ibset,mirror_side(iscan)+1,idet,2)
!						 xrvsc2(3)=xxrvs(it2,ibset,mirror_side(iscan)+1,idet,3)
!						 xrvsc2(4)=xxrvs(it2,ibset,mirror_side(iscan)+1,idet,4)
!						 xm12c2(1)=xxam12(it2,ibset,mirror_side(iscan)+1,idet,1)
!						 xm12c2(2)=xxam12(it2,ibset,mirror_side(iscan)+1,idet,2)
!						 xm12c2(3)=xxam12(it2,ibset,mirror_side(iscan)+1,idet,3)
!						 xm12c2(4)=xxam12(it2,ibset,mirror_side(iscan)+1,idet,4)
!						 xm13c2(1)=xxam13(it2,ibset,mirror_side(iscan)+1,idet,1)
!						 xm13c2(2)=xxam13(it2,ibset,mirror_side(iscan)+1,idet,2)
!						 xm13c2(3)=xxam13(it2,ibset,mirror_side(iscan)+1,idet,3)
!						 xm13c2(4)=xxam13(it2,ibset,mirror_side(iscan)+1,idet,4)

!@MJ@ below added for Aqua RVS and polarization corrections	
	                elseif (Platform .EQ. 'Aqua') then

						 do icof=1,6  ! @chg@ --> till the end of the loop
								ibset=1 ! Band08(in the Pol.Corr. table from Bryan Franz) 
								xrvsc1(icof)=xxrvs_aq(it1,ibset,mirror_side(iscan)+1,idet,icof)
								xrvsc2(icof)=xxrvs_aq(it2,ibset,mirror_side(iscan)+1,idet,icof)
								xm12c1(icof)=xxam12_aq(it1,ibset,mirror_side(iscan)+1,idet,icof)
								xm12c2(icof)=xxam12_aq(it2,ibset,mirror_side(iscan)+1,idet,icof)
								xm13c1(icof)=xxam13_aq(it1,ibset,mirror_side(iscan)+1,idet,icof)
								xm13c2(icof)=xxam13_aq(it2,ibset,mirror_side(iscan)+1,idet,icof)
	
!...	          --> get B03 and B10 parameters to correct for B03 RVS !@MJ@
								ibset=2 ! Band09(in the Pol.Corr. table from Bryan Franz)
								xrvsc1_2a(icof)=xxrvs_aq(it1,ibset,mirror_side(iscan)+1,idet,icof)
								xrvsc2_2a(icof)=xxrvs_aq(it2,ibset,mirror_side(iscan)+1,idet,icof)
								xm12c1_2a(icof)=xxam12_aq(it1,ibset,mirror_side(iscan)+1,idet,icof)
								xm12c2_2a(icof)=xxam12_aq(it2,ibset,mirror_side(iscan)+1,idet,icof)
								xm13c1_2a(icof)=xxam13_aq(it1,ibset,mirror_side(iscan)+1,idet,icof)
								xm13c2_2a(icof)=xxam13_aq(it2,ibset,mirror_side(iscan)+1,idet,icof)

                                ibset=3 ! Band10(in the Pol.Corr. table from Bryan Franz)                                
								xrvsc1_2b(icof)=xxrvs_aq(it1,ibset,mirror_side(iscan)+1,idet,icof)
								xrvsc2_2b(icof)=xxrvs_aq(it2,ibset,mirror_side(iscan)+1,idet,icof)
								xm12c1_2b(icof)=xxam12_aq(it1,ibset,mirror_side(iscan)+1,idet,icof)
								xm12c2_2b(icof)=xxam12_aq(it2,ibset,mirror_side(iscan)+1,idet,icof)
								xm13c1_2b(icof)=xxam13_aq(it1,ibset,mirror_side(iscan)+1,idet,icof)
								xm13c2_2b(icof)=xxam13_aq(it2,ibset,mirror_side(iscan)+1,idet,icof)
								
!...	          --> below to be aplided to B01 just for a test 
!... 	            (not yet to be applied as of Jan 29, 2009)->should not use!
								ibset=6 ! Band13(in the Pol.Corr. table from Bryan Franz)
								xrvsc1_3(icof)=xxrvs_aq(it1,ibset,mirror_side(iscan)+1,idet,icof)
								xrvsc2_3(icof)=xxrvs_aq(it2,ibset,mirror_side(iscan)+1,idet,icof)
								xm12c1_3(icof)=xxam12_aq(it1,ibset,mirror_side(iscan)+1,idet,icof)
								xm12c2_3(icof)=xxam12_aq(it2,ibset,mirror_side(iscan)+1,idet,icof)
								xm13c1_3(icof)=xxam13_aq(it1,ibset,mirror_side(iscan)+1,idet,icof)
								xm13c2_3(icof)=xxam13_aq(it2,ibset,mirror_side(iscan)+1,idet,icof)  
						 enddo  !  do icof=1,4

						 wt = (sec_inp-sectab_aq(it1))/(sectab_aq(it2)-sectab_aq(it1))
!@MJ@ above added for Aqua RVS and polarization corrections						 
					 endif
					 
					 do 10 iframe=1, edge(1)	           ! frame=1,1354
							i = iframe
							 
!             -- check lat/lon for fill values. Skip if found.
              if (latitude(i,j) < -900.0 .OR. longitude(i,j) < -900.0) go to 11

!             -- skip pixels at high SZA's
              if (solar_zenith_angle(i,j) > 72.0) go to 11
              
							! From DeepBlue.f90
				      lsf = landsea(iframe,ipix)
				      
				      ilat5 = (latitude(i,j) + 90.) *10.
              ilat5 = ilat5 + 1
              if (ilat5.gt.1800) ilat5 = 1800
              if (ilat5.lt.1)   ilat5 = 1
 
              ilon5 = (longitude(i,j) + 180.) *10.
              ilon5 = ilon5 + 1
              if (ilon5.gt.3600) ilon5 = 3600
              if (ilon5.lt.1)   ilon5 = 1
              
              sca    = -999.0
              cc     = 3.14159/180.
              psi    = acos(cos(solar_zenith_angle(i,j)*cc)*cos(sensor_zenith_angle(i,j)*cc) -  &
              &  sin(solar_zenith_angle(i,j)*cc)*sin(sensor_zenith_angle(i,j)*cc)*cos(relative_azimuth_angle(i,j)*cc))
              sca = 180. - psi/cc
           
              ndvi = (SR_band9(i,j) - SR_band1(i,j)) / (SR_band9(i,j) + SR_band1(i,j))
              sfc_refl_650 = get_LER650(ilat5, ilon5, ndvi, sca, relative_azimuth_angle(i,j))/100.0
              sfc_refl_412 = get_LER412(ilat5, ilon5, ndvi, sca, relative_azimuth_angle(i,j))/100.0
              sfc_refl_470 = get_LER470(ilat5, ilon5, ndvi, sca, relative_azimuth_angle(i,j))/100.0
          
!							if (lsf.eq.0 .or. lsf.gt.4) go to 11
!              if (lsf.eq.3) go to 11
            
            
!          Check for negative reflectivity in SR_band3 only in a very limited
!          geographic range.
           if (latitude(i,j) .le. 17.5 .and. latitude(i,j) .ge. 16.0 .and. &
               longitude(i,j) .ge. 14.5 .and. longitude(i,j) .le. 17.5 .and. &
               SR_band3(i,j) .lt. 0.0) go to 11

           cl_flag = -999.

           ix = 1
           xmin = 999.
           xmax = -999.
           do ii = 0, 2
              do jj = 0, 2
                 ii2 = ii -1
                 jj2 = jj -1
                 if (ii2+i .gt. 0 .and. ii2+i .lt. xdim+1 .and. &
                     jj2+j .gt. 0 .and. jj2+j .lt. ydim+1) then
                    xx = SR_band3(ii2+i,jj2+j)
                    if (xx.lt.xmin) xmin = xx
                    if (xx.gt.xmax) xmax = xx
                 endif
              enddo
           enddo
            
           if (xmax.lt.1.E-8.or.xmin.lt.1.E-8 ) go to 11
           var = xmax/xmin + delta           
                      
!         -- cirrus cloud tests (Jingfeng Huang)
!         ------------------------------------------------------------------------------------------
!         -- get precipitable water and total ozone from GDAS/NCEP data.
           pwv_in=3.0    ! @@@@@ PW(cm or g/cm^2), just for test (NCEP data should be fed)
           to3_in=300.0  ! @@@@@ O3(DU), just for test (NCEP data need be used)
           
           call GetAncData_PGE04(latitude(i,j),longitude(i,j),sfctmp,ugrd,vgrd,pwat,ozone)
           pwv_in=pwat/10.0  ! convert PW to cm from mm
           to3_in=ozone
           
           rr = SR_band8(i,j)/SR_band1(i,j)   
           rr2 = SR_band2(i,j)/SR_band3(i,j)
           
!          -- skip invalid 11um (band 31) data, used in cloud filters and D*.  
           if (SR_band6(i,j) < 0.0 .OR. SR_band7(i,j) < 0.0) go to 11
           
           bt11=0.0144*(1/(11.02e-6*log(7.3388e+2/SR_band6(i,j)+1)))
           
           btd8=0.0144*(1/(8.6e-6*log(25.354e+2/SR_band5(i,j)+1)))- &
           0.0144*(1/(11.02e-6*log(7.3388e+2/SR_band6(i,j)+1)))

           btd11=0.0144*(1/(11.02e-6*log(7.3388e+2/SR_band6(i,j)+1)))- &
           0.0144*(1/(12.02e-6*log(4.7534e+2/SR_band7(i,j)+1)))        

!          -- calculate D* parameter. Skip bad/cloudy pixels (signalled by Dstar1==-999)
           call Dstar(btd8,btd11,Dstar1)
           if (Dstar1 < -900.0) cycle   
!           dstar_arr(i,j) = Dstar1
!           if (Dstar1 > 1000.0) then
!              Dstar1 = 1000.0
!           endif
!           write(5000,'(21(F14.6,1x))') latitude(i,j), longitude(i,j), sfc_refl_650, sfc_refl_412, &
!           &  sfc_refl_470, pwv_in, bt11, btd8, btd11, Dstar1, ndvi, SR_band1(i,j), &
!           &  SR_band2(i,j), SR_band3(i,j), SR_band4(i,j), SR_band5(i,j), SR_band6(i,j), &
!           &  SR_band7(i,j), SR_band8(i,j), SR_band9(i,j), SR_band10(i,j) 
!           go to 11

!         -- @TODO - save gzflg and pass to modis() -> find_v rather than redoing lookup
!         --    in find_v().
          gzflg = terrain_flag_new(ilon5,ilat5)   ! intented cast to integer
          
!         -- if D* is high enough (indicates aerosols roughly), skip spatvar checks and
!         -- other cloud filters to avoid losing the aerosol plume.
!         -- Exclude N. Africa regions (gzone 1-5) though.
          if (Dstar1 >= 1.06 .AND. (gzflg < 1 .OR. (gzflg > 5 .AND. gzflg /= 26 .AND. gzflg /= 27))) go to 100 
           
!          -- use different spatial variability thresholds over bright and dark surfaces
           if (sfc_refl_650 >= 0.08 .AND. var > 1.1) go to 11
           if (sfc_refl_650 < 0.08  .AND. var > 1.15) go to 11

!           if (Dstar1 > 1000.0) then
!              Dstar1 = 1000.0
!           endif
!           write(1000,'(21(F14.6,1x))') latitude(i,j), longitude(i,j), sfc_refl_650, sfc_refl_412, &
!           &  sfc_refl_470, pwv_in, bt11, btd8, btd11, Dstar1, ndvi, SR_band1(i,j), &
!           &  SR_band2(i,j), SR_band3(i,j), SR_band4(i,j), SR_band5(i,j), SR_band6(i,j), &
!           &  SR_band7(i,j), SR_band8(i,j), SR_band9(i,j), SR_band10(i,j) 
!           go to 11
           
!          Over bright surfaces, apply RR1.38/0.66 and BTD11-12 threshold !JH
!          test addition - 6/7/2011
          if (sfc_refl_650 .ge. 0.08) then
             
             if (Dstar1 .lt. 1.12) then
             
              if (SR_band8(i,j).gt.0.018 .and. pwv_in.ge.0.9) go to 11 			! Thin Cirrus Screening, and Thin Cirrus Over-Correction Rectification
              if (SR_band8(i,j).gt.0.018 .and. pwv_in.lt.0.9 .and. btd11.gt.-1.0 .and. SR_band1(i,j).lt.0.55) go to 11 
              if (bt11 .lt. 270.0) go to 11 	! Cloud Screening
              if (sfc_refl_650 .ge. 0.16) then
                if (bt11 .ge. 270.0 .and. bt11 .lt. 281.0 .and. btd11 .gt. -0.5) go to 11    ! Cloud Edge Contamination Removal
              end if
             endif
             
           end if
          
!          Over dark surfaces

           if (sfc_refl_650.lt.0.08) then
            
            if (Dstar1 .lt. 1.12) then
           
              dd = SR_band1(i,j)/SR_band3(i,j)
           		if (SR_band4(i,j).gt.0.36 .and. dd.lt.0.95)     go to 11   ! Cloud Screening ?
            	if (SR_band8(i,j).gt.0.018 .and. pwv_in.gt.0.4) go to 11   ! Thin Cirrus Screening, and Thin Cirrus Over-Correction Rectification
              if (bt11 .lt. 270.0)                            go to 11 ! Cloud Screening
              if (bt11 .ge. 270.0 .and. bt11 .lt. 281.0 .and. SR_band4(i,j) .lt. 0.15) go to 11 		! Snow cover and cloud edge contamination removal
            	
            end if 
          
           end if

!         -- only apply these cloud tests over region 5 (N. Africa) and southern Africa.
!         -- developed using MYD021KM.A2011171.A12[25|30]*.hdf as test case.

!          -- over *really* bright surfaces, except Arabian Peninsula (6 <= gzflg <= 11)
          if (gzflg < 6 .OR. gzflg > 11) then 
            if (sfc_refl_650 > 0.25) then
              if (Dstar1 .lt. 0.85) go to 11
            end if
          end if
         
			    if ((gzflg == 5 .OR. gzflg == 1 .OR. gzflg == 26 .OR. gzflg == 27) .OR. (latitude(i,j) < 15.0 .AND. ((longitude(i,j) >-20.0 &
			    &   .AND. longitude(i,j) < 55.0) .AND. gzflg == -999))) then
		      
			      if (sfc_refl_650 < 0.08 .AND. Dstar1 .lt. 1.12) then 
			        if (SR_band2(i,j) > 0.15) then 
			          if (SR_band4(i,j) > 0.16) go to 11  ! used to be 0.16
            	  if (btd11 > 1.3) go to 11           ! used to be 1.3
						  end if
						end if
						if (sfc_refl_650 > 0.08 .AND. Dstar1 .lt. 1.12) then
						  if (SR_band2(i,j) > 0.25) then
						    if (btd11 > 3.0) go to 11
						  end if
           	end if
          end if
         
!           if (Dstar1 > 1000.0) then
!              Dstar1 = 1000.0
!           endif
!           write(1000,'(21(F14.6,1x))') latitude(i,j), longitude(i,j), sfc_refl_650, sfc_refl_412, &
!           &  sfc_refl_470, pwv_in, bt11, btd8, btd11, Dstar1, ndvi, SR_band1(i,j), &
!           &  SR_band2(i,j), SR_band3(i,j), SR_band4(i,j), SR_band5(i,j), SR_band6(i,j), &
!           &  SR_band7(i,j), SR_band8(i,j), SR_band9(i,j), SR_band10(i,j) 
!           go to 11
          
!         -- try to detect snow cover and skip pixel.
!         -- source: http://modis-snow-ice.gsfc.nasa.gov/?c=atbd&t=atbd
          ndsi = -999.0
          if (SR_band11(i,j) > -900.0 .AND. SR_band4(i,j) > -900.0) then
            ndsi = (SR_band11(i,j) - SR_band4(i,j)) / (SR_band11(i,j) + SR_band4(i,j))
          else
            go to 11
          end if
          if (ndsi > 0.35 .AND. SR_band9(i,j) > 0.11 .AND. SR_band11(i,j) > 0.1 .AND. bt11 < 283.0) then
            go to 11
          end if
          
!					-- only apply this to IBGP region 7 and 8 (in MJ's IGBP table), north of 45N.
!					-- inserted to deal with some snow/cloud issues in Canada (see 2008/276-279).
          if ((xlcvr_2(ilon5, ilat5) == 7 .OR. xlcvr_2(ilon5, ilat5) == 8) .AND. latitude(i,j) > 45.0) then
          	if ((ndsi > 0.05 .AND. ndsi < 0.35) .AND. SR_band9(i,j) > 0.25 .AND. SR_band11(i,j) > 0.25) then
            	go to 11
          	endif
					endif
!##################################################################
!##################################################################
 100       continue

!c...      Need gas absorption correction here ... !MJ 08/04/2011
           sza_in=solar_zenith_angle(i,j)
           vza_in=sensor_zenith_angle(i,j)

           srb1=SR_band1(i,j)
           srb2=SR_band2(i,j)
           srb4=SR_band4(i,j)
           srb9=SR_band9(i,j)
           srb10=SR_band10(i,j)

!c... Calling a routine to correct TOA refl for gas-absorption-effect
          call corr_trn_gas(sza_in,vza_in,pwv_in,to3_in,  &
                           srb1,srb2,srb4,srb9,srb10)     ! MJ 08/04/2011
!c... To turn off corrections for certain band(s), comment out corresponding
!c... line(s) below.
!           SR_band1(i,j)=srb1    ! MODIS B01
!           SR_band2(i,j)=srb2    ! MODIS B03
!           SR_band4(i,j)=srb4    ! MODIS B07
!           SR_band9(i,j)=srb9    ! MODIS B02
!           SR_band10(i,j)=srb10  ! MODIS B05
!##################################################################
!##################################################################

           tmp(1) = latitude(i,j)
           tmp(2) = longitude(i,j)
           tmp(3) = solar_zenith_angle(i,j)
           tmp(4) = sensor_zenith_angle(i,j)
           tmp(5) = relative_azimuth_angle(i,j)
           tmp(6) = SR_band1(i,j)
           tmp(7) = SR_band2(i,j)
           tmp(8) = SR_band3(i,j)
           tmp(9) = cl_flag
           tmp(10)= sensor_azimuth_angle(i,j)
           tmp(11)= solar_azimuth_angle(i,j) 
           tmp(12)= stdv02(i,j)
           
           tmpvg(1)=SR_band1(i,j)  ! B1  ! MJ@@@ 07/26/11
           tmpvg(2)=SR_band2(i,j)  ! B3
           tmpvg(3)=SR_band3(i,j)  ! B8
           tmpvg(4)=SR_band4(i,j)  ! B7
           tmpvg(5)=SR_band9(i,j)  ! B2
           tmpvg(6)=SR_band10(i,j) ! B5


! Use this version of dump_data for POLCOR
!           call dump_data (fh,i,j,tmp:1(1),tmp(2),tmp(3),tmp(4),tmp(5),tmp(6),&
!            tmp(7),tmp(8),tmp(10),tmp(11))
           					  
!@MJ@					 if (Platform .EQ. 'Terra') then ! also for Aqua  
							pixnum = float(iframe)-1.0D0 ! according to Bryan Franz 
							if (is_subset) then   ! have to recalculate original pixel number for subsets
 							  pixnum = float(iframe+pixel_offsets(1)-1)-1.0D0
							end if
							
							ap = pixnum                  ! 
							ap2 = ap*ap
							ap3 = ap2*ap
							ap4 = ap3*ap
							ap5 = ap4*ap
							x1 = xrvsc1(1)+ap*xrvsc1(2)+ap2*xrvsc1(3)+ap3*xrvsc1(4)+ap4*xrvsc1(5)+ap5*xrvsc1(6)
							x2 = xrvsc2(1)+ap*xrvsc2(2)+ap2*xrvsc2(3)+ap3*xrvsc2(4)+ap4*xrvsc2(5)+ap5*xrvsc2(6)
							xrvs(1) = x1 + wt*(x2-x1)
	            
							x1 = xm12c1(1)+ap*xm12c1(2)+ap2*xm12c1(3)+ap3*xm12c1(4)+ap4*xm12c1(5)+ap5*xm12c1(6)
							x2 = xm12c2(1)+ap*xm12c2(2)+ap2*xm12c2(3)+ap3*xm12c2(4)+ap4*xm12c2(5)+ap5*xm12c2(6)
							xam12(1) = x1 + wt*(x2-x1)
	
							x1 = xm13c1(1)+ap*xm13c1(2)+ap2*xm13c1(3)+ap3*xm13c1(4)+ap4*xm13c1(5)+ap5*xm13c1(6)
							x2 = xm13c2(1)+ap*xm13c2(2)+ap2*xm13c2(3)+ap3*xm13c2(4)+ap4*xm13c2(5)+ap5*xm13c2(6)
							xam13(1) = x1 + wt*(x2-x1)

!@MJ@ below added to use pre-launch polarization parameters
              zaoi=10.5 + (65.5 - 10.5)*pixnum/(1354.0) ! calculate AOI
!c...         get index with which interpolation will be made.
              call set_indx_prepol(aoiref,zaoi,iaoi1,iaoi2)
              wtz=(zaoi-aoiref(iaoi1))/(aoiref(iaoi2)-aoiref(iaoi1))

!c... For B8, Aqua/MODIS polcor will be computed using pre-launch characterization,
!c... whereas Terra/MODIS polcor will be don using OBPG's xcal data.
!c... (Ref. Kwiatkowska et al., 2008, Appl. Opt.; Jeong et al., 2010, IEEE TGRS)
              if (Platform .EQ. 'Aqua') then         ! use pre-launch polarization data
                  x1 = b8ym12ref(idet,iaoi1,mirror_side(iscan)+1) 
                  x2 = b8ym12ref(idet,iaoi2,mirror_side(iscan)+1) 
                  xam12(1) = x1 + wtz*(x2-x1)

                  x1 = b8ym13ref(idet,iaoi1,mirror_side(iscan)+1) 
                  x2 = b8ym13ref(idet,iaoi2,mirror_side(iscan)+1) 
                  xam13(1) = x1 + wtz*(x2-x1)
              endif ! if (Platform .EQ. 'Aqua') then

!... @chg@ below - to get xrvs values for the other bands (B03 [and B01])
!@MJ@         x1 = xrvsc1_2(1)+ap*xrvsc1_2(2)+ap2*xrvsc1_2(3)+ap3*xrvsc1_2(4)
!@MJ@         x2 = xrvsc2_2(1)+ap*xrvsc2_2(2)+ap2*xrvsc2_2(3)+ap3*xrvsc2_2(4)
!@MJ@         xrvs(2) = x1 + wt*(x2-x1)           

!@MJ@         x1 = xm12c1_2(1)+ap*xm12c1_2(2)+ap2*xm12c1_2(3)+ap3*xm12c1_2(4)
!@MJ@         x2 = xm12c2_2(1)+ap*xm12c2_2(2)+ap2*xm12c2_2(3)+ap3*xm12c2_2(4)
!@MJ@         xam12(2) = x1 + wt*(x2-x1)

!@MJ@         x1 = xm13c1_2(1)+ap*xm13c1_2(2)+ap2*xm13c1_2(3)+ap3*xm13c1_2(4)
!@MJ@         x2 = xm13c2_2(1)+ap*xm13c2_2(2)+ap2*xm13c2_2(3)+ap3*xm13c2_2(4)
!@MJ@         xam13(2) = x1 + wt*(x2-x1)

              x1 = xrvsc1_2a(1)+ap*xrvsc1_2a(2)+ap2*xrvsc1_2a(3)+ap3*xrvsc1_2a(4)+ap4*xrvsc1_2a(5)+ap5*xrvsc1_2a(6)
              x2 = xrvsc2_2a(1)+ap*xrvsc2_2a(2)+ap2*xrvsc2_2a(3)+ap3*xrvsc2_2a(4)+ap4*xrvsc2_2a(5)+ap5*xrvsc2_2a(6)
              xrvs_tmp1 = x1 + wt*(x2-x1)           

              x1 = xrvsc1_2b(1)+ap*xrvsc1_2b(2)+ap2*xrvsc1_2b(3)+ap3*xrvsc1_2b(4)+ap4*xrvsc1_2b(5)+ap5*xrvsc1_2b(6)
              x2 = xrvsc2_2b(1)+ap*xrvsc2_2b(2)+ap2*xrvsc2_2b(3)+ap3*xrvsc2_2b(4)+ap4*xrvsc2_2b(5)+ap5*xrvsc2_2b(6)
              xrvs_tmp2 = x1 + wt*(x2-x1)           

!@MJ@ RVS for B3 linearly interpolated w.r.t. wavelength using RVSs for B9 and B10.
!@MJ@ polcor for B3 are based on pre-launch characterizations.
              wt=(466.0-443.0)/(488.0-443.0)
              x1=xrvs_tmp1
              x2=xrvs_tmp2
              xrvs(2) = x1 + wt*(x2-x1)           

              if (Platform .EQ. 'Aqua') then
                  x1 = b3ym12ref(idet,iaoi1,mirror_side(iscan)+1) 
                  x2 = b3ym12ref(idet,iaoi2,mirror_side(iscan)+1) 
                  xam12(2) = x1 + wtz*(x2-x1)

                  x1 = b3ym13ref(idet,iaoi1,mirror_side(iscan)+1) 
                  x2 = b3ym13ref(idet,iaoi2,mirror_side(iscan)+1) 
                  xam13(2) = x1 + wtz*(x2-x1)
              elseif (Platform .EQ. 'Terra') then
                  x1 = b3om12ref(idet,iaoi1,mirror_side(iscan)+1) 
                  x2 = b3om12ref(idet,iaoi2,mirror_side(iscan)+1) 
                  xam12(2) = x1 + wtz*(x2-x1)

                  x1 = b3om13ref(idet,iaoi1,mirror_side(iscan)+1) 
                  x2 = b3om13ref(idet,iaoi2,mirror_side(iscan)+1) 
                  xam13(2) = x1 + wtz*(x2-x1)
              endif ! if (Platform .EQ. 'Aqua') then

!@MJ@ below 11 lines are place holder for B1 just in case needed in the future 
              x1 = xrvsc1_3(1)+ap*xrvsc1_3(2)+ap2*xrvsc1_3(3)+ap3*xrvsc1_3(4)+ap4*xrvsc1_3(5)+ap5*xrvsc1_3(6)
              x2 = xrvsc2_3(1)+ap*xrvsc2_3(2)+ap2*xrvsc2_3(3)+ap3*xrvsc2_3(4)+ap4*xrvsc2_3(5)+ap5*xrvsc2_3(6)
              xrvs(3) = x1 + wt*(x2-x1)           

              x1 = xm12c1_3(1)+ap*xm12c1_3(2)+ap2*xm12c1_3(3)+ap3*xm12c1_3(4)+ap4*xm12c1_3(5)+ap5*xm12c1_3(6)
              x2 = xm12c2_3(1)+ap*xm12c2_3(2)+ap2*xm12c2_3(3)+ap3*xm12c2_3(4)+ap4*xm12c2_3(5)+ap5*xm12c2_3(6)
              xam12(3) = x1 + wt*(x2-x1)

              x1 = xm13c1_3(1)+ap*xm13c1_3(2)+ap2*xm13c1_3(3)+ap3*xm13c1_3(4)+ap4*xm13c1_3(5)+ap5*xm13c1_3(6)
              x2 = xm13c2_3(1)+ap*xm13c2_3(2)+ap2*xm13c2_3(3)+ap3*xm13c2_3(4)+ap4*xm13c2_3(5)+ap5*xm13c2_3(6)
              xam13(3) = x1 + wt*(x2-x1)
!... @chg@ above - to get xrvs values for the other bands (B03 [and B01])

!             read data from temporary binary file
!             call extract_data(fh,idim,jdim,xlat,xlong,sza,xthet,xphi, &
!                            xnvalm5,vazim,sazim,ierr)
!							call extract_data(fh,idim,jdim,xlat,xlong,sza,xthet,xphi,xnvalm5, &
!														 band26,vazim,sazim,LandSeaFlag,ierr)
					  	xlat 	= tmp(1)
							xlong = tmp(2)
							sza		= tmp(3)
							xthet = tmp(4)
							xphi	= tmp(5)
							xnvalm5(1) = tmp(6)
							xnvalm5(2) = tmp(7)
							xnvalm5(3) = tmp(8)
!c.						band26 = tmp(9)
							vazim = tmp(10)
							sazim = tmp(11)							

!						  print *, xlat,xlong,sza,xthet,xphi,xnvalm5, &
!															 band26,vazim,sazim !,LandSeaFlag,ierr
!...          note: interpolation does not work when x3 (RAA) = 180.
							if(xphi.ge.179.9999) xphi=179.9999

!...   				df_azi = SAZ - VAZ (determines the sign of the stoke vector, U)
							df_azi=sazim-vazim
							if(df_azi.lt.0.) df_azi=360.0+df_azi
!
!...   ==================================================================
!...     compute alpha (angle; in degree)
							call calc_alpha(xlong,xlat,xthet,vazim,tt_ecr(:,iscan),alpha)
!...   ==================================================================
!...   ==================================================================
!...   get simulated I, Q, and U for Rayleigh atmosphere (at 412nm)										
						call get_rayl_table(nvalc, sza, xthet, xphi, Lr_q(1), Lr_u(1), Lr_x(1))
						call get_rayl_table(nvalc2,sza, xthet, xphi, Lr_q(2), Lr_u(2), Lr_x(2))
						call get_rayl_table(nvalc3,sza, xthet, xphi, Lr_q(3), Lr_u(3), Lr_x(3))
!...   I,Q, and U for the other two bands (B01 & B03) were set to zero.
!...   ==================================================================

!...        =========================
!...        =========================
!...          make units consistent with MODIS L1B refl. 
!...          to be consistent with MODIS L1B refl. @@@@@
              Lr_q(1)=Lr_q(1)*pi
              Lr_u(1)=Lr_u(1)*pi
              Lr_x(1)=Lr_x(1)*pi

              Lr_q(2)=Lr_q(2)*pi
              Lr_u(2)=Lr_u(2)*pi
              Lr_x(2)=Lr_x(2)*pi

              Lr_q(3)=Lr_q(3)*pi
              Lr_u(3)=Lr_u(3)*pi
              Lr_x(3)=Lr_x(3)*pi
!							print *, 'Lr_q',Lr_q
!							print *, 'Lr_u',Lr_u
!							print *, 'Lr_x',Lr_x
!... dependence of the sign of U on azimuthal angles
						if(df_azi.lt.180.0) Lr_u(:)=-1.0*Lr_u(:)
!...        =========================
!...        =========================

!... note: xrvs, xam12, and xam13 are now vectors! (Jan29,2009) @chg@
             cossza=cos(d2r*sza)
             lx_m(1)=xnvalm5(3)*cossza/xrvs(1)  ! 412nm --> OLD METHOD, PRE EQN CORR - Eqn A
! TESTING ----------------------------------------------------------------------
!             lx_m(1)=xnvalm5(3)*cossza!/xrvs(1)  ! 412nm --> OLD METHOD, PRE EQN CORR - Eqn A
! TESTING ----------------------------------------------------------------------
             lx_m(2)=xnvalm5(2)*cossza/xrvs(2)  ! 466nm
!...         lx_m(3)=xnvalm5(1)*cossza/xrvs(3)  ! 650nm --> not yet
             lx_m(3)=xnvalm5(1)*cossza          ! 650nm
!             lx_m(1)=xnvalm5(3)*cossza  				! 412nm --> for test06 - Eqns B,C,D
!             lx_m(2)=xnvalm5(2)*cossza  				! 466nm
!             lx_m(3)=xnvalm5(1)*cossza				  ! 650nm --> not yet
             Lr_xtm(1)=Lr_x(1)/xrvs(1)  				! pol.crr. only to Rayleigh
             Lr_xtm(2)=Lr_x(2)/xrvs(2)
             Lr_xtm(3)=Lr_x(3)/xrvs(3)

!... Note: xlat, xphi, xam12, xam13, xrvs, and alpha are scalar, while
!... 	  	   Lr_q, Lr_u, and polcor are vectors with three elements(3wvl).
!...  	     Currently, polcor(1)=1.0 and polcor(2)=1.0
!						 print *, 'First',xlat, xphi, xam12, xam13, xrvs, alpha
!						 print *, 'Lr_q',Lr_q
!						 print *, 'Lr_u',Lr_u
!						 print *, 'Lr_x',Lr_x
							
						 call polar_correct_xcal( xlat, xphi, xam12, xam13, xrvs, &		! EqnA
																			alpha, Lr_q, Lr_u, lx_m, polcor )
!						 call polar_correct_xcal( xlat, xphi, xam12, xam13, xrvs, &		! Eqns B,C,D
!																			alpha, Lr_q, Lr_u, Lr_xtm, polcor )

!...          a temporary setup to write "polcor" for 412nm
              polcor2d(iframe,ipix)=polcor(1)		! B08
              polcor2d_2(iframe,ipix)=polcor(2) ! @chg@ B03
!...          polcor2d_3(iframe,ipix)=polcor(3) ! @chg@ B01, not yet
!...    ******************************************************
!...    for test
              xrvs_test2d(iframe,ipix)=xrvs(1)   ! @chg@
              xam12_test2d(iframe,ipix)=xam12(1)
              xam13_test2d(iframe,ipix)=xam13(1)
              alpha_test2d(iframe,ipix)=alpha
              lr_q_test2d(iframe,ipix)=Lr_q(1)
              lr_u_test2d(iframe,ipix)=Lr_u(1)
              lr_x_test2d(iframe,ipix)=Lr_x(1)
              xrvs_test2d_2(iframe,ipix)=xrvs(2)      ! @chg@
              xam12_test2d_2(iframe,ipix)=xam12(2)   
              xam13_test2d_2(iframe,ipix)=xam13(2)   
              lr_q_test2d_2(iframe,ipix)=Lr_q(2)  
              lr_u_test2d_2(iframe,ipix)=Lr_u(2)  
!...          xrvs_test2d_3(iframe,ipix)=xrvs(3)      ! @chg@
!...          xam12_test2d_3(iframe,ipix)=xam12(3)    ! not yet
!...          xam13_test2d_3(iframe,ipix)=xam13(3)    ! to be 
!...          lr_q_test2d_3(iframe,ipix)=Lr_q(3)      ! used
!...          lr_u_test2d_3(iframe,ipix)=Lr_u(3)  
!...          lr_x_test2d_3(iframe,ipix)=Lr_x(3) 
!...    ******************************************************

!...   ==================================================================
!...   ==================================================================
!... --> procedures for correction of MODIS measured reflectance for 
!...     polarization effects need be here; e.g., 
!...     real(single), dimension(1354,2030)  :: rfl_pcr1,rfl_pcr2,rfl_pcr3
!...     rfl_pcr1(iframe,ipix)=xnvalm5(3)/xrvs(1)/polcor(1)  ! @chg@ (B08)
!...     rfl_pcr2(iframe,ipix)=xnvalm5(2)/xrvs(2)/polcor(2)  ! @chg@ (B03)
!...     rfl_pcr3(iframe,ipix)=xnvalm5(1)/xrvs(3)/polcor(3)  ! @chg@ (B01)
!...   ==================================================================
!...   ==================================================================
	
						! Haven't actually applied the correction yet... Just testing I/O
!  	         call dump_data (fh,idim,jdim,xlat,xlong,sza,xthet,xphi,xnvalm5, &
!                               band26,vazim,sazim,LandSeaFlag)

!       	For TERRA, apply polarization correction -----------------------------

!---			For 412 - SR_Band3 - tmp8  - B08 -------------------------------------
 							rfl_case1B08(iframe,ipix)=xnvalm5(3)
 							rfl_case2B08(iframe,ipix)=xnvalm5(3)/xrvs(1)
 							rfl_case3B08(iframe,ipix)=xnvalm5(3)/polcor(1)
              lr_x_test2d_2(iframe,ipix)=Lr_x(2) 
							rfl_case4B08(iframe,ipix)=xnvalm5(3)/xrvs(1)/polcor(1)
!              print *, xlat, xlong, xrvs(1), polcor(1), xnvalm5(3), rfl_case1B08(iframe,ipix), rfl_case4B08(iframe,ipix)
!							rfl_case4B08(iframe,ipix)=(Lr_x(1)/polcor(1)+xnvalm5(3)-Lr_x(1))/xrvs(1)
!							rfl_case4B08(iframe,ipix)=(Lr_x(1)/cossza/polcor(1)+xnvalm5(3)-Lr_x(1)/cossza)/xrvs(1)  ! jul 14
		        if (Platform .EQ. 'Terra') then 
		        tmp(8) = rfl_case4B08(iframe,ipix)  ! original, RVS+polcor(C006)
		        else
		          tmp(8) = rfl_case1B08(iframe,ipix)  ! original, RVS+polcor(C006)
		        end if
! 						tmp(8) = rfl_case1B08(iframe,ipix)  ! No correction (for B8)
!						tmp(8) = rfl_case2B08(iframe,ipix)  ! RVS corr. only (for B8), test
!						tmp(8) = rfl_case3B08(iframe,ipix)  ! pol. corr. only (for B8), test
							
							
!---			For 470 - SR_Band2 - tmp7 - B03 --------------------------------------
 							rfl_case1B03(iframe,ipix)=xnvalm5(2)
 							rfl_case2B03(iframe,ipix)=xnvalm5(2)/xrvs(2)
 							rfl_case3B03(iframe,ipix)=xnvalm5(2)/polcor(2)
 							rfl_case4B03(iframe,ipix)=xnvalm5(2)/xrvs(2)/polcor(2)  ! to be default !@MJ@
!
!							rfl_case4B03(iframe,ipix)=(Lr_x(2)/polcor(2)+xnvalm5(2)-Lr_x(2))/xrvs(2) ! test
              if (Platform .EQ. 'Terra') then 
 							  tmp(7) = rfl_case4B03(iframe,ipix)  ! RVS+polcor(C006) for B3
 							else
 							  tmp(7) = rfl_case1B03(iframe,ipix)  ! RVS+polcor(C006) for B3
 							end if
!							tmp(7) = rfl_case1B03(iframe,ipix)  ! No correction for B3
							
!@MJ@						endif ! for Terra only - Polarization correction -->also for Aqua
            if (Platform .EQ. 'Terra') then 
              tmpvg(2)=rfl_case4B03(iframe,ipix) ! No Crr.    SR_band2(i,j)  ! B3  ! MJ@@@ 07/26/11
              tmpvg(3)=rfl_case4B08(iframe,ipix) ! No Crr.    SR_band3(i,j)  ! B8  ! MJ@@@ 07/26/11
            else
              tmpvg(2)=rfl_case1B03(iframe,ipix) ! No Crr.    SR_band2(i,j)  ! B3  ! MJ@@@ 07/26/11
              tmpvg(3)=rfl_case1B08(iframe,ipix) ! No Crr.    SR_band3(i,j)  ! B8  ! MJ@@@ 07/26/11
            end if
!            tmpvg(2)=rfl_case4B03(iframe,ipix) 
!            tmpvg(3)=rfl_case4B08(iframe,ipix) ! TEMPORARY FOR TESTING, DELETE for C006!
	          
!            if (abs(xlat-10.5) > 0.05) cycle
!            if (abs(xlong+1.0) > 0.05) cycle
            
!            if (xlat < 29.9 .OR. xlat > 30.1) cycle
!            if (xlong < 114.9 .OR. xlong > 115.1) cycle
!c          call dump_data (fh,i-1,j-1,tmp,lsf,Dstar1)
            call dump_data (fh,i-1,j-1,tmp,lsf,Dstar1,tmpvg) !MJ@@@
 11        continue
 10     continue
 12   continue
 13   continue
!10				 enddo  !... do iframe=1, edge(1) !frame=1,1354
!12			 enddo    !... do idet=1,10
!13		 enddo      !... do iscan=1, 203
      msg = "Done repackaging MODIS data for Deep Blue."
      call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC, msg, 'DeepBlue')
			
!      print *, 'Starting debug output...'
!      file_name = trim(l1bfiles_base)//".dstar.hdf"
!      status = db_debug_output(trim(file_name), latitude, longitude, dstar_arr)
!      if (status /= 0) then
!        print *, "ERROR: Failed to create debug output file: ", status
!        stop
!      end if

!     Clean ups this mess 
      
      deallocate(bands)
      deallocate(SR_band1)
      deallocate(SR_band2)
      deallocate(SR_band3)
      deallocate(SR_band4)
      deallocate(SR_band5)
      deallocate(SR_band6)
      deallocate(SR_band7)
      deallocate(SR_band8)
      deallocate(SR_band9)
      deallocate(SR_band10)
      deallocate(SR_band11)
      deallocate(stdv02)
      deallocate(stdv26)
      deallocate(stdv03)

!     -- corrections-related arrays
      deallocate(polcor2d)
      deallocate(polcor2d_2)
      deallocate(xrvs_test2d)
      deallocate(xam12_test2d)
      deallocate(xam13_test2d)
      deallocate(alpha_test2d)
      deallocate(lr_q_test2d)
      deallocate(lr_u_test2d)
      deallocate(lr_x_test2d)
      deallocate(xrvs_test2d_2)
      deallocate(xam12_test2d_2)
      deallocate(xam13_test2d_2)
      deallocate(lr_q_test2d_2)
      deallocate(lr_u_test2d_2)
      deallocate(lr_x_test2d_2)
      deallocate(rfl_case1B08)
      deallocate(rfl_case2B08)
      deallocate(rfl_case3B08)
      deallocate(rfl_case4B08)
      deallocate(rfl_case1B03)
      deallocate(rfl_case2B03)
      deallocate(rfl_case3B03)
      deallocate(rfl_case4B03)
			
!     Now process all that data!!!
      call modis(fh, QAflags, QAflags_O, DTaot, dims)
      prtn = pgs_io_gen_closef(fh)

!      close(5000)

!     Update PGE04 metadata with Deep Blue attributes 
      call MOD04_MetaDataUpdate()
      
      msg = "Successfully completed Deep Blue processing of MODIS data."
      call MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS, msg, 'DeepBlue')

!     -- unload various input data
      call unload_landcover(status)
      if (status /= 0) then
        msg = "Failed to unload land cover database."
        call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC, msg, 'DeepBlue')
      end if
      
      call unload_brdf(status)
      if (status /= 0) then
        msg = "ERROR: Unable to unload BRDF data."
        call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC, msg, 'DeepBlue')
        stop
      end if
      
! Go home

end program DeepBlue

