      subroutine get_angles(modfil_Geo,cube,scan_type,max_sol,
     +                      min_sol,buf_geo_size,rlat,rlon,elev,
     +                      view_zen,solar_zen,cube_topog,rel_angle)

      implicit none
      save

c     Common block used for writing debug output statements
      common / bug / debug, h_output

c----------------------------------------------------------------------
c!F77
c
c!Description:
c     Routine for reading and putting MODIS cloud mask geometry
c     values into arrays for processing.
c
c!Input Parameters:
c modfil_Geo    Geolocation file array
c cube          Scan cube counter
c scan_type     Day or night mode of scan
c max_sol       Maximum solar zenith angle of this granule
c min_sol       Minimum solar zenith anlge of this granule
c buf_geo_size  Line and element granule size
c
c!Output Parameters:
c
c rlat          Array containing scan cube of 1km pixel latitude
c               values
c rlon          Array containing scan cube of 1km pixel longitude
c               values
c elev          Array containing scan cube of 1km pixel elevation
c               values
c view_zen      Array containing scan cube of 1km pixel viewing
c               zenith angles
c solar_zen     Array containing scan cube of 1km pixel solar
c               azimulth angles
c cube_topog    Array containing scan cube of 1 km land/sea
c               information
c rel_angle     Array containing scan cube of 1km pixel relative
c               angles
c
c!Revision History:
c
c!Team-Unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c Revision 01.00  1996/2/13 08:20:37
c K. Strabala (kathys@ssec.wisc.edu)
c Original version.
c
c!END
c----------------------------------------------------------------------

        include 'mapi.inc'
        include 'global.inc'

        real pi
        parameter(pi = 3.14159)
        real dtr
        parameter(dtr = pi/180.0)
        real rtd
        parameter(rtd = 180.0/pi)

c       scalar arguments
        character*1 scan_type
        integer cube
        real max_sol,min_sol

c       scalar arrays
        integer buf_geo_size(2),modfil_Geo(MODFILLEN),
     +          cube_topog(npixel,scans_cube)
        real rlat(npixel,scans_cube),rlon(npixel,scans_cube),
     +       view_zen(npixel,scans_cube),solar_zen(npixel,scans_cube),
     +       rel_angle(npixel,scans_cube),elev(npixel,scans_cube)

c       local scalars
        integer rtn,ii,jj,mm,nn,ElementSize,LineSize,debug,h_output,i,j
        real vzar,szar,razr,cossna

c       local arrays
        integer data_size(2)
        real buf(npixel,scans_cube),
     +       sens_azim(npixel,scans_cube),solar_azim(npixel,scans_cube),
     +       rel_azim(npixel,scans_cube)

c       external subroutines
        external message

c       external functions
        integer Read_GEO_V2
        external Read_GEO_V2

        intrinsic abs,sin,cos,acos,nint

c ... debug statement ............................................
      if (debug .gt. 3) then
         write(h_output,'(10x/,''Within get_angles routine '',/)')
         write(h_output,'(10x/,''Scan mode is: '',A1,/)') scan_type
      endif
c ...............................................................

c ...   Grab the latitude information out of the geolocation file
c ...   First initialize output paramaters
         data_size(1) = 0
         data_size(2) = 0
         rtn = 0
         ElementSize = buf_geo_size(1)
         LineSize = buf_geo_size(2)

         do 134 mm = 1 , buf_geo_size(2)
           do 135 nn = 1 , buf_geo_size(1)
              buf(nn,mm) = bad_data
              rlat(nn,mm) = bad_data
              rlon(nn,mm) = bad_data
              elev(nn,mm) = bad_data
              view_zen(nn,mm) = bad_data
              solar_zen(nn,mm) = bad_data
              sens_azim(nn,mm) = bad_data
              solar_azim(nn,mm) = bad_data
              rel_azim(nn,mm) = bad_data
              cube_topog(nn,mm) = -1
              rel_angle(nn,mm) = 999.0
  135      continue
  134    continue

c ------------------------------------------------------------------
c ...    Read the latitude values
         rtn = Read_GEO_V2(modfil_Geo,1,cube,ElementSize,LineSize,
     +                     buf,data_size)

         if (rtn .ne. 0 )then
            call message ('get_angles',
     +      'Failed to extract latitude from Geo. file' //
     +      ' [OPERATOR ACTION: Check geolocation file, rerun PGE]',
     +      0, 2 )
         endif

         do 145 ii = 1 , data_size(2)
         do 145 jj = 1 , data_size(1)
          if ( buf(jj,ii) .ge. -90.0 .and. buf(jj,ii) .le.  90.0 )
     +      rlat(jj,ii) = buf(jj,ii)
 145     continue
c ------------------------------------------------------------------


c ------------------------------------------------------------------
c ...   Read the longitude values
         do i = 1 , buf_geo_size(2)
           do  j = 1 , buf_geo_size(1)
              buf(j,i) = bad_data
           enddo
         enddo

         rtn = Read_GEO_V2(modfil_Geo,2,cube,ElementSize,LineSize,
     +                     buf,data_size)

         if( rtn .ne. 0 )then
            call message ('get_angles',
     +      'Failed to extract longitude from Geo. file' //
     +      ' [OPERATOR ACTION: Check geolocation file, rerun PGE]',
     +      0, 2 )
         endif

         do 155 ii = 1 , data_size(2)
         do 155 jj = 1 , data_size(1)
          if ( buf(jj,ii) .ge. -180.0 .and. buf(jj,ii) .le.  180.0 )
     +      rlon(jj,ii) = buf(jj,ii)
 155     continue
c ------------------------------------------------------------------


c ------------------------------------------------------------------
c ...   Read elevation values
         do i = 1 , buf_geo_size(2)
           do  j = 1 , buf_geo_size(1)
              buf(j,i) = bad_data
           enddo
         enddo

         rtn = Read_GEO_V2(modfil_Geo,3,cube,ElementSize,LineSize,
     +                     buf,data_size)

         if( rtn .ne. 0 )then
            call message ('get_angles',
     +      'Failed to extract elevation from Geo. file' //
     +      ' [OPERATOR ACTION: Check geolocation file, rerun PGE]',
     +      0, 2 )
         endif

         do 156 ii = 1 , data_size(2)
         do 156 jj = 1 , data_size(1)
          if ( buf(jj,ii) .ge. -500.0 .and. buf(jj,ii) .le.  9000.0 )
     +      elev(jj,ii) = buf(jj,ii)
 156     continue
c ------------------------------------------------------------------


c ------------------------------------------------------------------
c ...   Read the solar zenith angle information
         do i = 1 , buf_geo_size(2)
           do  j = 1 , buf_geo_size(1)
              buf(j,i) = bad_data
           enddo
         enddo

         rtn = Read_GEO_V2(modfil_Geo,7,cube,ElementSize,LineSize,
     +                     buf,data_size)

         if( rtn .ne. 0 )then
            call message ('get_angles',
     +      'Failed to extract solar zenith from Geo. file' //
     +      ' [OPERATOR ACTION: Check geolocation file, rerun PGE]',
     +      0, 2 )
         endif

         do 165 ii = 1 , data_size(2)
         do 165 jj = 1 , data_size(1)
           if (buf(jj,ii) .ge. 0.0 .and. buf(jj,ii) .le.  180.0 ) then
            solar_zen(jj,ii) = buf(jj,ii)
c ...       Determine max and min solar zenith angles of granule
            if (nint(solar_zen(jj,ii)) .ne. -32767 .and.
     +          solar_zen(jj,ii).gt.max_sol)  max_sol = solar_zen(jj,ii)
            if (nint(solar_zen(jj,ii)) .ne. -32767 .and.
     +          solar_zen(jj,ii).lt.min_sol)  min_sol = solar_zen(jj,ii)
           endif
 165     continue
c ------------------------------------------------------------------


c ------------------------------------------------------------------
c ...   Read the sensor zenith angle information
         do i = 1 , buf_geo_size(2)
           do  j = 1 , buf_geo_size(1)
              buf(j,i) = bad_data
           enddo
         enddo

         rtn = Read_GEO_V2(modfil_Geo,4,cube,ElementSize,LineSize,
     +                     buf,data_size)

         if( rtn .ne. 0 )then
            call message ('get_angles',
     +      'Failed to extract sensor zenith from Geo. file' //
     +      ' [OPERATOR ACTION: Check geolocation file, rerun PGE]',
     +      0, 2 )
         endif

         do 175 ii = 1 , data_size(2)
         do 175 jj = 1 , data_size(1)
           if (buf(jj,ii) .ge. 0.0 .and. buf(jj,ii) .le.  180.0 )
     +       view_zen(jj,ii) = buf(jj,ii)
 175     continue
c ------------------------------------------------------------------

c ------------------------------------------------------------------
c ...   Read the Land Sea flag
         do i = 1 , buf_geo_size(2)
           do  j = 1 , buf_geo_size(1)
              buf(j,i) = bad_data
           enddo
         enddo

         rtn = Read_GEO_V2(modfil_Geo,9,cube,ElementSize,LineSize,
     +                     buf,data_size)

         if( rtn .ne. 0 )then
            call message ('get_angles',
     +      'Failed to extract land/sea mask info from Geo. file' //
     +      ' [OPERATOR ACTION: Check geolocation file, rerun PGE]',
     +      0, 2 )
         endif

         do 180 ii = 1 , data_size(2)
         do 180 jj = 1 , data_size(1)
            cube_topog(jj,ii) = nint(buf(jj,ii))
 180     continue
c ------------------------------------------------------------------


c     Added check for day or night scan mode.  Read the
c     azimuth data and determine a relative anlge only
c     if in day mode
      if (scan_type .eq. 'D') then

c ------------------------------------------------------------------
c ...   Read the solar azimuth angle information
         do i = 1 , buf_geo_size(2)
           do  j = 1 , buf_geo_size(1)
              buf(j,i) = bad_data
           enddo
         enddo

         rtn = Read_GEO_V2(modfil_Geo,8,cube,ElementSize,LineSize,
     +                     buf,data_size)

         if( rtn .ne. 0 )then
            call message ('get_angles',
     +      'Failed to extract solar zenith from Geo. file' //
     +      ' [OPERATOR ACTION: Check geolocation file, rerun PGE]',
     +      0, 2 )
         endif

         do 185 ii = 1 , data_size(2)
         do 185 jj = 1 , data_size(1)
           if (buf(jj,ii) .ge. -180.0 .and. buf(jj,ii) .le.  180.0 )
     +       solar_azim(jj,ii) = buf(jj,ii)
 185     continue
c ------------------------------------------------------------------


c ------------------------------------------------------------------
c ...   Read the sensor azimuth angle information
         do i = 1 , buf_geo_size(2)
           do  j = 1 , buf_geo_size(1)
              buf(j,i) = bad_data
           enddo
         enddo

         rtn = Read_GEO_V2(modfil_Geo,5,cube,ElementSize,LineSize,
     +                     buf,data_size)

         if( rtn .ne. 0 )then
            call message ('get_angles',
     +      'Failed to extract sensor azimuth from Geo. file' //
     +      ' [OPERATOR ACTION: Check geolocation file, rerun PGE]',
     +      0, 2 )
         endif

         do 195 ii = 1 , data_size(2)
         do 195 jj = 1 , data_size(1)
           if (buf(jj,ii) .ge. -180.0 .and. buf(jj,ii) .le.  180.0 )
     +       sens_azim(jj,ii) = buf(jj,ii)
 195     continue
c ------------------------------------------------------------------


c ------------------------------------------------------------------
c ...    get relative azimuth angle
         do 205 ii = 1 , data_size(2)
         do 205 jj = 1 , data_size(1)
            rel_azim(jj,ii) = abs(180.0-abs(sens_azim(jj,ii)-
     +                        solar_azim(jj,ii)))
 205     continue

c       Now calculate the relative angle (value that sun glint is based
c        upon.)
         do 215 ii = 1 , data_size(2)
         do 215 jj = 1 , data_size(1)
            if (nint(view_zen(jj,ii)) .ne. nint(bad_data) .and.
     +          nint(solar_azim(jj,ii)) .ne. nint(bad_data) .and.
     +          nint(solar_zen(jj,ii)) .ne. nint(bad_data) .and.
     +          nint(sens_azim(jj,ii)) .ne. nint(bad_data)) then
                vzar = view_zen(jj,ii) * dtr
                szar = solar_zen(jj,ii) * dtr
                razr = rel_azim(jj,ii) * dtr
                cossna = min( (sin(vzar) * sin(szar) * cos(razr) +
     +                   cos(vzar) * cos(szar)) , 1.0)
                rel_angle(jj,ii) = acos(cossna) * rtd
            endif
 215     continue
c ------------------------------------------------------------------
      endif

c ... debug statement ............................................
      if (debug .gt. 3) then
        write(h_output,'(15x,'' ANGLE VALUES: '',2i10)') data_size(1),
     *        data_size(2)
        write(h_output,'(2x,'' Latitude'',/,20(5f10.4/))')
     +       ((rlat(j,i),j=825,834),i=1,10)
        write(h_output,'(2x,'' Longitude'',/,20(5f10.4/))')
     +       ((rlon(j,i),j=825,834),i=1,10)
        write(h_output,'(2x,'' Elevation'',/,20(5f10.4/))')
     +       ((elev(j,i),j=825,834),i=1,10)
        write(h_output,'(2x,'' Solar Zenith'',/,20(5f10.4/))')
     +       ((solar_zen(j,i),j=825,834),i=1,10)
        write(h_output,'(2x,'' Viewing Zenith'',/,20(5f10.4/))')
     +       ((view_zen(j,i),j=825,834),i=1,10)
        write(h_output,'(2x,'' Solar Azimuth '',/,20(5f10.4/))')
     +       ((solar_azim(j,i),j=825,834),i=1,10)
        write(h_output,'(2x,'' Sensor Azimuth '',/,20(5f10.4/))')
     +       ((sens_azim(j,i),j=825,834),i=1,10)
        write(h_output,'(2x,'' Relative Angle'',/,20(5f10.4/))')
     +       ((rel_angle(j,i),j=825,834),i=1,10)
        write(h_output,'(2x,'' Land/Sea Mask'',/,20(5I10/))')
     +       ((cube_topog(j,i),j=825,834),i=1,10)
        write(h_output,'(15x,''MAX/MIN SOLAR ZENITH ANGLES:'',
     +        2f10.4,/)') max_sol,min_sol
      endif
c ................................................................

         return
         end
