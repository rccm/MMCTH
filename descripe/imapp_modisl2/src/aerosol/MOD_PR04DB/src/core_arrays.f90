module core_arrays
!-----------------------------------------------------------------------
! !F90
!
! !DESCRIPTION:
!
!    core_arrays declares data arrays for MODIS level 1b data.
!
! !INPUT PARAMETERS:  none
!
! !OUTPUT PARAMETERS:  none
!
! !REVISION HISTORY:
!
!    Initial Version by Jennifer Wei, July 8, 2004
!    Jeremy Warner added prolog and several arrays, December 1, 2006
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
! !END
!-----------------------------------------------------------------------

   use GeneralAuxType

   implicit none

   !input data arrays (measurement)
   integer(kind=1),allocatable, dimension (:,:)   :: landsea
   real(single),   allocatable, dimension (:,:)   :: latitude
   real(single),   allocatable, dimension (:,:)   :: longitude
   real(single),   allocatable, dimension (:,:)   :: solar_zenith_angle
   real(single),   allocatable, dimension (:,:)   :: sensor_zenith_angle
   real(single),   allocatable, dimension (:,:)   :: solar_azimuth_angle
   real(single),   allocatable, dimension (:,:)   :: sensor_azimuth_angle
   real(single),   allocatable, dimension (:,:)   :: relative_azimuth_angle

   logical(kind=4),   allocatable, dimension (:,:)   :: cloud_flag, cirrus_flag
   real(single),   allocatable, dimension (:,:,:) :: band_measurements

   integer(integer_fourbyte),allocatable, dimension (:)     :: bands

   !output data arrays (after calculation)
   real(single),   allocatable, dimension (:,:)   :: SR_band1
   real(single),   allocatable, dimension (:,:)   :: SR_band2
   real(single),   allocatable, dimension (:,:)   :: SR_band3
   real(single),   allocatable, dimension (:,:)   :: SR_band4
   real(single),   allocatable, dimension (:,:)   :: SR_band5
   real(single),   allocatable, dimension (:,:)   :: SR_band6
   real(single),   allocatable, dimension (:,:)   :: SR_band7
   real(single),   allocatable, dimension (:,:)   :: SR_band8
   real(single),   allocatable, dimension (:,:)   :: SR_band9
   real(single),   allocatable, dimension (:,:)   :: SR_band10
   real(single),   allocatable, dimension (:,:)   :: SR_band11

	 !for polarization correction - added by CES 05 dec 2008
	 real(double),   allocatable, dimension (:,:,:) :: tt_inst2ecr
   real(double),   allocatable, dimension (:,:)   :: tt_ecr

   integer(kind=2),allocatable, dimension (:) :: mirror_side

   real(double),   allocatable, dimension (:)     :: ev_cntr_time  ! [MJadded]

end module core_arrays


