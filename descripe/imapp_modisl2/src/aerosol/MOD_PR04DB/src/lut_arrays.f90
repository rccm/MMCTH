module lut_arrays

!f90
!
!description:
!           data arrays declaration from MODIS level 1b & Geolocation
!           general parameters declaration
!
!input parameters: none
!
!
!output parameters: none
!revision history:
!team-unique header:
!             mark gray, lac
!             climate and radiation branch
!             nasa goddard space flight center
!             greenbelt, maryland, u.s.a.
!
!modified by M.J. Jeong
!             parameters and variables related to polarization  
!             sensitivity (stored as LUT in a hdf file) are added.
!modified by Myeong-Jae Jeong (Oct 2010)
!        value of "ntimer", related to C006 Terra/MODIS RVS corr. is changed
!        "ntimer_aq" related to C006 Aqua/MODIS RVS corr. is added

   use GeneralAuxType

   implicit none

   !Polarization LUT table arrays
    integer, parameter :: ndet=10	       ! number of detectors per channel
    integer, parameter :: naoi=6	       ! number of mirror incidence angles in LUT
    integer, parameter :: nside=2	       ! number of mirror sides
    integer, parameter :: nbandlut=11	       ! number of bands in LUT (incl. TDI bands)
    				   	       ! 11 bands=8,9,10,11,12,13L,13H,14L,14H,15,16
    integer, parameter :: nbands=3	       ! number of bands in polar. correction from MODIS, band 1,3,8
!   integer, parameter :: ntimer=200     ! number of times (~months) for which polarization coeff. were derived(C051), Terra
    integer, parameter :: ntimer=127     ! number of times for which Terra/MODIS polar. coeff. were derived(C006), 20100928
    integer, parameter :: ntimer_aq=121   ! number of times for which Aqua/MODIS polar. coeff. were derived(C006), 20101012
    integer, parameter :: nbandlutr=9    ! number of bands in LUT, 9 bands (B08-B16)  
    integer, parameter :: ncoeff=6       ! number of polynomial coefficients (3rd order polynomical) 
    integer, dimension(8) :: chindx=(/ 1,2,3,4,5,6,10,11 /) ! I don't understand this variable

    integer, dimension(nbands) :: bindx=(/ 1,3,8 /)    ! DeepBlue targetted bands in MODIS
   
    real(single),dimension (naoi)     :: angle_of_incident
    real(single),dimension (nbandlut,ndet,naoi,nside) :: am12, am13

!... xxrvs, xxam12, xxam13  ! pol. coeff. for ocn bands 
!... xxyear, xxday          ! year and DoY for which the pol. coeff. were derived
!... xxwave                 ! wavelengths for ocean bands
!... sectab                 ! reference time for LUT (seconds)
    real(double),dimension (ntimer,nbandlutr,nside,ndet,ncoeff) :: xxrvs, xxam12, xxam13 
    integer(integer_twobyte),dimension (ntimer)   :: xxyear, xxday   
    real(single),dimension (nbandlutr)            :: xxwave           
    real(double),dimension (ntimer)               :: sectab          

    integer(integer_twobyte),dimension (ntimer_aq)       :: xxyear_aq, xxday_aq   
    real(single),dimension (nbandlutr)            :: xxwave_aq        
    real(double),dimension (ntimer_aq,nbandlutr,nside,ndet,ncoeff) :: xxrvs_aq, xxam12_aq, xxam13_aq 
    real(double),dimension (ntimer_aq)                   :: sectab_aq          
  
end module lut_arrays


