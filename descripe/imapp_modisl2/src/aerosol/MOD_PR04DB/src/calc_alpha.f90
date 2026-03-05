!* -----------------------------------------------------------------------------
!* calc_alpha() - computes the frame rotation angle for polarization correction 
!*                                                                              
!* adapted from Seadas src in 'l1_modis_hdf.c'  (<--Seadas5.1)                                 
!*                                                                              
!* arguments:	lon: longitude
!* 		lat: latitude
!*		senz:sensor zenith angle
!*		sena:sensor azimuth angle
!*		mnorm:taken from MODIS geolocation T_inst2ECR
!* output: alpha
!* modified by JCW
!* version 1: keep same way as original file
!*
!* modified by MJ
!* The sign of alpha was changed to conform Gordon/Meister papers 
!* and to conform "l1_modis_hdf.c" in seadas5.2).
!* -----------------------------------------------------------------------------

  subroutine calc_alpha (lon, lat, senz, sena, mnorm, alpha)
         
      use GeneralAuxType  

      implicit none

      real(single),  intent(in):: lon, lat, senz,sena
      real(double), dimension(3), intent(in):: mnorm
      real(double),  intent(out):: alpha

! local variables 
      integer		:: i
      real(single)	:: slon, clon
      real(single)	:: slat, clat
      real(single)	:: szen, czen
      real(single)	:: sazi, cazi

      real(single), dimension(3) :: e, n, v, r, s
      real(single)    	 :: sdvcr, vdr, sdr, sdv 

!    /* invert mirror normal */
      do i = 1, 3
	 		r(i) = -mnorm(i)
      enddo

      slon = sin(d2r*lon) 
      clon = cos(d2r*lon) 
      slat = sin(d2r*lat) 
      clat = cos(d2r*lat) 
      szen = sin(d2r*senz) 
      czen = cos(d2r*senz) 
      sazi = sin(d2r*sena) 
      cazi = cos(d2r*sena) 

!     /* pixel coordinate system (north, east, vertical) in ECR */

       e(1) = -slon 
       e(2) =  clon 
       e(3) =  0.0 

       n(1) = -slat * clon 
       n(2) = -slat * slon 
       n(3) =  clat 

       v(1) =  clat * clon 
       v(2) =  clat * slon 
       v(3) =  slat 

!     /* sensor view vector in ECR */

       do  i = 1, 3
          s(i) = e(i) * szen * sazi         &
                + n(i) * szen * cazi        &
                + v(i) * czen 
       enddo

!     /* compute rotation angle (alpha) from pixel normal (v) to mirror */
!     /* normal (r) about sensor view vector (s)  (Wertz, p. 757)       */

       sdvcr = s(1) * (v(2)*r(3)-v(3)*r(2)) &
             + s(2) * (v(3)*r(1)-v(1)*r(3)) &
             + s(3) * (v(1)*r(2)-v(2)*r(1)) 

       vdr   = v(1)*r(1) + v(2)*r(2) + v(3)*r(3) 

       sdr   = s(1)*r(1) + s(2)*r(2) + s(3)*r(3) 

       sdv   = v(1)*s(1) + v(2)*s(2) + v(3)*s(3) 

!      /* negated to be consistent with Gordon et al. */
!...   alpha = (1/d2r_d) * atan2(sdvcr,(vdr-sdr*sdv)) - 90.0  ! <--Seadas5.1

       alpha = -1.0D0 * ((1/d2r_d) * atan2(sdvcr,(vdr-sdr*sdv)) - 90.0) 


      end subroutine calc_alpha

