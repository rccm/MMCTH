subroutine get_modis_data_cube (level1b_name, geolocation_name, start, edge, stride, spbands, status)
!-----------------------------------------------------------------------
! !F90
!
! !DESCRIPTION:
!
!    get_modis_data_cube extracts MODIS L1B data for use by Deep Blue
!    algorithm.  Core sds of interest in the MODIS level 1B file are 
!    the following:
!   latitude (MOD03)
!   longitude (MOD03)
!   solar_zenith_angle (MOD03)
!   sensor_zenith_angle (MOD03)
!   solar_azimuth_angle (MOD03)
!   sensor_azimuth_angle (MOD03)
!   band_measurements (MOD02)
!
! !INPUT PARAMETERS:
!
!    Type       Name             Description
!    ====       ====             ===========
!    CHARACTER  level1b_name     name of MODIS L1B file (MOD02)
!    CHARACTER  geolocation_name name of MODIS L1B geolocation file (MOD03)
!    INTEGER*4  start            array of start locations for read of SDS
!    INTEGER*4  edge             array of end locations for read of SDS
!    INTEGER*4  stride           array of strides for read of SDS
!    INTEGER*4  spbands          array of bands to read
!
! !OUTPUT PARAMETERS:
!
!    Type       Name             Description
!    ====       ====             ===========
!    INTEGER*4  status           return code
!
! !REVISION HISTORY:
!
!    Initial Version by Jeremy Warner   12/01/2006
!     (not really, but I don't know who wrote the original
!      that I modified, probably J. Wei).
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
!     MODIS_F_GENERIC         (MODIS_39500.f)
!     landsea                 (core_arrays.f90)
!     latitude                (core_arrays.f90)
!     longitude               (core_arrays.f90)
!     solar_zenith_angle      (core_arrays.f90)
!     sensor_zenith_angle     (core_arrays.f90)
!     solar_azimuth_angle     (core_arrays.f90)
!     sensor_azimuth_angle    (core_arrays.f90)
!     relative_azimuth_ang    (core_arrays.f90)
!     band_measurements       (core_arrays.f90)
!     bands                   (core_arrays.f90)
!     SR_band1                (core_arrays.f90)
!     SR_band2                (core_arrays.f90)
!     SR_band3                (core_arrays.f90)
!     SR_band4                (core_arrays.f90)
!     SR_band5                (core_arrays.f90)
!     SR_band6                (core_arrays.f90)
!     SR_band7                (core_arrays.f90)
!
!   Functions:
!
!    HDF functions:
!     SFstart
!     sffinfo
!     sfn2index
!     sfselect
!     sfginfo
!     sfrdata
!     sfendacc
!     sffattr
!     SFgainfo
!     sfrnatt
!     sfend
!
! !END
!-----------------------------------------------------------------------


   use GeneralAuxType
   use core_arrays
   use hdf

   implicit none

   include 'PGS_MODIS_39500.f'

	 integer, dimension (2), intent (in)       :: start, edge, stride
   integer, dimension (:), intent (in)       :: spbands
  
   character(*),           intent (in)       :: level1b_name, geolocation_name
   integer,                intent (out)      :: status

   real(double)                              :: d_scale_factor  
   real(single), allocatable                 :: scale_factor(:), offset_factor(:)
   integer                                   :: numberofbands ,checkvariable, i, j
   integer                                   :: scaleflag, offsetflag, index, value, count, array_index
   integer(integer_fourbyte)                 :: rank,  number_type, nattrs
   integer(integer_fourbyte)                 :: sds_id,sds_index,attr_index, hdfid
   integer(integer_fourbyte), dimension(2)   :: start2d, edge2d, stride2d, dimsizes_2d
   integer(integer_fourbyte), dimension(3)   :: start3d, edge3d, stride3d, dimsizes_3d
   integer(integer_fourbyte), dimension(1)   :: start1dmr,edge1dmr,stride1dmr, &
                                                dimsizes_1dmr
   integer(integer_fourbyte), dimension(1)   :: start1dmr2, edge1dmr2, &
                                                stride1dmr2, dimsizes_1dmr2
   integer(integer_twobyte),dimension(:,:), allocatable  :: temp_array
   integer                                   :: xdimension, ydimension, pixels_across_scan,temp_x, band_to_tweak
   character(len=64)                         :: sds_name, local_sds_name, attr_name, &
                                                scale_factor_name, offset_factor_name
   integer(integer_fourbyte), dimension(1)   :: nscan
	 character(256)                            ::  msg
   
   status=0
   numberofbands = size(spbands)
   scaleflag = 1
   offsetflag = 0
   start2d = start
   edge2d = edge
   stride2d = stride

   start3d  = (/ start2d(1), start2d(2), 0 /)
   edge3d   = (/ edge2d(1), edge2d(2),numberofbands /)
   stride3d = (/ stride2d(1), stride2d(2), 1 /)

   xdimension =   edge2d(1) 
   ydimension =   edge2d(2)
   
   nscan = edge(2)/10
   
!   print *, edge2d, nscan
!    print *,'spbands in cube SR =',spbands   
!==============================================================================
! start the geolocation file
!==============================================================================
   hdfid = SFstart(geolocation_name,DFACC_READ)
   status=sffinfo(hdfid, number_type, nattrs)
! Latitude
   sds_name = 'Latitude'
   sds_index = sfn2index(hdfid, sds_name)
   sds_id    = sfselect(hdfid , sds_index)

   status=sfginfo(sds_id,local_sds_name,rank, dimsizes_2d, number_type , nattrs)
   !print *,'SDS name',sds_name
   
! Use this variable to make sure that dimensions (ncans) is the same between MOD04, MOD03/MOD021KM files
   if (dimsizes_2d(2) /= ydimension) then
     msg = "Dimensions not equal between files. Program exiting. Notify SDST - input files may need reprocessing."
     print *, msg
     print *, 'From M?D04_L2: ',xdimension, ydimension
     print *, 'From M?D03:    ',dimsizes_2d
     call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
   endif
   
   status=sfrdata(sds_id,start2d,stride2d,edge2d,latitude)
   status = sfendacc(sds_id)
	 
! Longitude
   sds_name = 'Longitude'
   sds_index = sfn2index(hdfid, sds_name)
   sds_id    = sfselect(hdfid , sds_index)

   status=sfginfo(sds_id,local_sds_name,rank, dimsizes_2d, number_type , nattrs)

   status=sfrdata(sds_id,start2d,stride2d,edge2d,longitude)
   status = sfendacc(sds_id)

! landsea
   sds_name = 'Land/SeaMask'
   sds_index = sfn2index(hdfid, sds_name)
   sds_id    = sfselect(hdfid , sds_index)

   status=sfginfo(sds_id,local_sds_name,rank, dimsizes_2d, number_type , nattrs)

   status=sfrdata(sds_id,start2d,stride2d,edge2d,landsea)
   status = sfendacc(sds_id)

! SensorZenith 
   allocate(temp_array(edge2d(1), edge2d(2)))
   sds_name = 'SensorZenith'
   sds_index = sfn2index(hdfid, sds_name)
   sds_id    = sfselect(hdfid , sds_index)

   status=sfginfo(sds_id,local_sds_name,rank, dimsizes_2d, number_type , nattrs)
   attr_index= sffattr(sds_id, 'scale_factor')

   status = sfgainfo(sds_id, attr_index, attr_name, number_type, count)
   status = sfrnatt(sds_id, attr_index, d_scale_factor)

   status=sfrdata(sds_id,start2d,stride2d,edge2d,temp_array)
   status = sfendacc(sds_id)
   sensor_zenith_angle = temp_array * d_scale_factor

! SensorAzimuth
   sds_name = 'SensorAzimuth'
   sds_index = sfn2index(hdfid, sds_name)
   sds_id    = sfselect(hdfid , sds_index)

   status=sfginfo(sds_id,local_sds_name,rank, dimsizes_2d, number_type , nattrs)
   attr_index= sffattr(sds_id, 'scale_factor')

   status = sfrnatt(sds_id, attr_index, d_scale_factor)

   status=sfrdata(sds_id,start2d,stride2d,edge2d,temp_array)
   status = sfendacc(sds_id)

   sensor_azimuth_angle = temp_array * d_scale_factor

! SolarZenith
   sds_name = 'SolarZenith'
   sds_index = sfn2index(hdfid, sds_name)
   sds_id    = sfselect(hdfid , sds_index)

   status=sfginfo(sds_id,local_sds_name,rank, dimsizes_2d, number_type , nattrs)
   attr_index= sffattr(sds_id, 'scale_factor')

   status = SFgainfo(sds_id, attr_index, attr_name, number_type, count)
   status = sfrnatt(sds_id, attr_index, d_scale_factor)

   status=sfrdata(sds_id,start2d,stride2d,edge2d,temp_array)
   status = sfendacc(sds_id)

   solar_zenith_angle =  temp_array * d_scale_factor

! SolarAzimuth
   sds_name = 'SolarAzimuth'
   sds_index = sfn2index(hdfid, sds_name)
   sds_id    = sfselect(hdfid , sds_index)

   status=sfginfo(sds_id,local_sds_name,rank, dimsizes_2d, number_type , nattrs)
   attr_index= sffattr(sds_id, 'scale_factor')

   status = SFgainfo(sds_id, attr_index, attr_name, number_type, count)
   status = sfrnatt(sds_id, attr_index, d_scale_factor)

   status=sfrdata(sds_id,start2d,stride2d,edge2d,temp_array)
   status = sfendacc(sds_id)

   solar_azimuth_angle = temp_array * d_scale_factor

!   status = sfend(hdfid) 

   relative_azimuth_angle =  (sensor_azimuth_angle - solar_azimuth_angle) -180.
   !print *, 'relative azimuth angle',relative_azimuth_angle

   do i = 1, size(relative_azimuth_angle,1)
     do j = 1,  size(relative_azimuth_angle,2)
       if (relative_azimuth_angle(i,j) > 180.) relative_azimuth_angle(i,j) = relative_azimuth_angle(i,j) - 360. 
       if (relative_azimuth_angle(i,j) < -180.) relative_azimuth_angle(i,j) = relative_azimuth_angle(i,j) + 360.
       if (relative_azimuth_angle(i,j) < 0.) relative_azimuth_angle(i,j) = -relative_azimuth_angle(i,j) 
     enddo
   enddo

! Now needs to read a different type(size) of SDS - added by CES 05 dec 2008
! mirror_side
   sds_name = 'Mirror side'
   sds_index = sfn2index(hdfid, sds_name)
   sds_id    = sfselect(hdfid , sds_index)

   start1dmr  = (/ 0 /)     
   edge1dmr   = (/ nscan /)
   stride1dmr = (/ 1 /)
   status=sfrdata(sds_id,start1dmr,stride1dmr,edge1dmr,mirror_side)
   status = sfendacc(sds_id) 
!   if (nscan(1) /= 203) then
!   		print *, 'Checking mirror_side: ',mirror_side
!	 endif
	 
! T_inst2ECR
   sds_name = 'T_inst2ECR'
   sds_index = sfn2index(hdfid, sds_name)
   sds_id    = sfselect(hdfid , sds_index)

! Note, we only take first set of vectors (anti-velocity) direction
! check MOD03 ATBD, personal talk with Fred Patt
   start3d  = (/ 0, 0, 0 /)
   edge3d   = (/ 1, 3, nscan /)
   stride3d = (/ 1, 1, 1 /)

   status=sfginfo(sds_id,local_sds_name,rank, dimsizes_3d, number_type , nattrs)
   status=sfrdata(sds_id,start3d,stride3d,edge3d,tt_ecr)
   status = sfendacc(sds_id)

! ev_cntr_time    
   sds_name = 'EV center time'
   sds_index = sfn2index(hdfid, sds_name)
   sds_id    = sfselect(hdfid , sds_index)
   
   start1dmr2 = (/ 0 /)     
   edge1dmr2  = (/ nscan /)
   stride1dmr2= (/ 1 /)
   status=sfrdata(sds_id,0,1,nscan,ev_cntr_time)
   status = sfendacc(sds_id)
   !print *, 'SDS_NAME', sds_name                Reads ok at this point 

! close this geolocation HDF file
   status = sfend(hdfid) 
!   print *,'status=',status
   start3d  = (/ start2d(1), start2d(2), 0 /)
   edge3d   = (/ edge2d(1), edge2d(2),numberofbands /)
   stride3d = (/ stride2d(1), stride2d(2), 1 /)
!   print *,'after '  
!   print *,'var bands =',bands,size(bands) 
!==============================================================================
! start the level1b file
!==============================================================================
   array_index = -999
   hdfid = sfstart(level1b_name,DFACC_READ)
   do i = 1, size(bands)
     sds_name = ''
     if (bands(i) >=1 .and. bands(i) <=2) then
	sds_name = 'EV_250_Aggr1km_RefSB'
!	sds_name = 'BAND1'
	scale_factor_name = 'reflectance_scales'
	offset_factor_name = 'reflectance_offsets'
	array_index = bands(i)
     endif 
     if (bands(i) >=3 .and. bands(i) <=7) then
	sds_name = 'EV_500_Aggr1km_RefSB'
!	sds_name = 'BAND3'
	scale_factor_name = 'reflectance_scales'
	offset_factor_name = 'reflectance_offsets'
	array_index = bands(i) - 2
     endif
     if (bands(i) >=8 .and. bands(i) <=19 .or. bands(i) == 26) then
	sds_name = 'EV_1KM_RefSB'
!	sds_name = 'BAND8'
	scale_factor_name = 'reflectance_scales'
	offset_factor_name = 'reflectance_offsets'
	if (bands(i) <= 12) array_index = bands(i) - 7
	if (bands(i) == 13) array_index =  6
	if (bands(i) == 14) array_index =  8
	if (bands(i) > 14 .and. bands(i) <=19 ) array_index = bands(i) - 5
	if (bands(i) == 26)  array_index = 15
     endif
     if (bands(i) >=27 .and. bands(i) <=36) then
        sds_name = 'EV_1KM_Emissive'
!       sds_name = 'BAND29'
        scale_factor_name = 'radiance_scales'
        offset_factor_name = 'radiance_offsets'
        array_index = bands(i) - 20
     endif

!====
     sds_index = sfn2index(hdfid, sds_name)
     sds_id    = sfselect(hdfid , sds_index)

     status=sfginfo(sds_id,local_sds_name,rank, dimsizes_3d, number_type , nattrs)
     attr_index= sffattr(sds_id, scale_factor_name)
     status = SFgainfo(sds_id, attr_index, attr_name, number_type, count)
     allocate(scale_factor(count))
     allocate(offset_factor(count))

     status = sfrnatt(sds_id, attr_index, scale_factor)

     attr_index= sffattr(sds_id, offset_factor_name)
     status = sfrnatt(sds_id, attr_index, offset_factor)
     start3d = (/start2d(1), start2d(2),array_index-1/)
     stride3d = (/stride2d(1), stride2d(2), 1/)
     edge3d = (/edge3d(1), edge3d(2), 1/)

     status=sfrdata(sds_id,start3d,stride3d,edge3d,temp_array)
     band_measurements(:,i,:) = scale_factor(array_index) *( temp_array - offset_factor(array_index))
     if (i == 7) then
!        print *, 'band measurements =',band_measurements(:,i,:)
     endif 
     status = sfendacc(sds_id)
     deallocate(scale_factor)
     deallocate(offset_factor)
   enddo
  
!  scale reflectances by cosine of the solar zenith angle
   do i = 1, size(bands)
     if(bands(i) < 20 .or. bands(i) == 26) then
       band_measurements(:,i,:) = band_measurements(:,i,:) / cos(d2r*solar_zenith_angle) 
     endif
   enddo

   status = sfend(hdfid)
   deallocate(temp_array)

end subroutine get_modis_data_cube

