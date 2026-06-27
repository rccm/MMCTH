      SUBROUTINE GET_SFCPRES_PWAT(LAT, LON, SFCPRES, GDAS_PWAT)
      
C-----------------------------------------------------------------------
C!F77
C
C!DESCRIPTION:
C     Gets bilinear interpolated surface pressure for a given lat/lon.
C
C!INPUT PARAMETERS:
C     LAT        Latitude (degrees, -90S to +90N)
C     LON        Longitude (degrees, -180W to +180E, Greenwich zero)
C
C!OUTPUT PARAMETERS:
C     SFCPRES    Interpolated surface pressure (hPa)
C     GDAS_PWAT  Interpolated precipitable water (mm)
C
C!REVISION HISTORY:
C     Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
C
C March 03
C     SWS Modified get_sfcpres.f to read PWAT
C
C!END
C-----------------------------------------------------------------------

      implicit none
      
c ... Input arguments     
      real lat, lon
      
c ... Output arguments
      real sfcpres
      real gdas_pwat
      
c ... Local variables
      real pres(0:15), temp(0:15), mixr(0:15),
     &  land, sfctmp, pwat, ugrd, vgrd, ozone, icec, sst
      integer nise
      integer met_date(4), ozn_date(4), ice_date(4), sst_date(4),
     &  nise_date(4)
      real lat1, lat2, lon1, lon2
      real pwat1, pwat2, pwat3, pwat4
      real p1, p2, p3, p4
      real t, u

c ... External functions
      integer floor, ceil
      external floor, ceil
            
c ... Set default return value
      sfcpres = -1.0

c ... Get coordinates of grid points surrounding this lat/lon
      lat1 = real(floor(lat))
      lat2 = real(ceil(lat))
      if (abs(lat2 - lat1) .lt. 1.0e-6) lat2 = lat1 + 1.0
      lon1 = real(floor(lon))
      lon2 = real(ceil(lon))
      if (abs(lon2 - lon1) .lt. 1.0e-6) lon2 = lon1 + 1.0

c ... Get the surface pressures at the 4 grid points
c ... surrounding this lat/lon
      call get_ancillary(lat1, lon1, pres, temp, mixr, land,
     &  sfctmp, p1, pwat1, ugrd, vgrd, ozone, icec, sst, nise,
     &  met_date, ozn_date, ice_date, sst_date, nise_date )
      if ((p1 .lt. 10.0) .or. (p1 .gt. 1200.0)) return
      if ((pwat1 .lt. 0.0) .or. (pwat1 .gt. 110.0)) return
      
      call get_ancillary(lat1, lon2, pres, temp, mixr, land,
     &  sfctmp, p2, pwat2, ugrd, vgrd, ozone, icec, sst, nise,
     &  met_date, ozn_date, ice_date, sst_date, nise_date )
      if ((p2 .lt. 10.0) .or. (p2 .gt. 1200.0)) return
      if ((pwat2 .lt. 0.0) .or. (pwat2 .gt. 110.0)) return

      call get_ancillary(lat2, lon2, pres, temp, mixr, land,
     &  sfctmp, p3, pwat3, ugrd, vgrd, ozone, icec, sst, nise,
     &  met_date, ozn_date, ice_date, sst_date, nise_date )
      if ((p3 .lt. 10.0) .or. (p3 .gt. 1200.0)) return
      if ((pwat3 .lt. 0.0) .or. (pwat3 .gt. 110.0)) return
     
      call get_ancillary(lat2, lon1, pres, temp, mixr, land,
     &  sfctmp, p4, pwat4, ugrd, vgrd, ozone, icec, sst, nise,
     &  met_date, ozn_date, ice_date, sst_date, nise_date )
      if ((p4 .lt. 10.0) .or. (p4 .gt. 1200.0)) return
      if ((pwat4 .lt. 0.0) .or. (pwat4 .gt. 110.0)) return

c ... Compute bilinear interpolated surface pressure
      u = (lat - lat1) / (lat2 - lat1)
      t = (lon - lon1) / (lon2 - lon1)
      sfcpres =
     &  (1.0 - t) * (1.0 - u) * p1 +
     &  t * (1.0 - u) * p2 +
     &  t * u * p3 +
     &  (1.0 - t) * u * p4

      gdas_pwat =
     &  (1.0 - t) * (1.0 - u) * pwat1 +
     &  t * (1.0 - u) * pwat2 +
     &  t * u * pwat3 +
     &  (1.0 - t) * u * pwat4

      END
