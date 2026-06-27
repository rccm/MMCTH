      subroutine get_anc_data(h_eco1,h_eco2,cube,buf_anc_size,max_pixel,
     +                        rlat,rlon,TAI_end,cube_pw,cube_eco,
     +                        cube_ice,cube_snow,cube_sfctmp,cube_msl,
     +                        cube_ugrd,cube_vgrd,cube_sst,qa_bits)

      implicit none
      save

c ... Common statement for debug purposes
      common / bug / debug, h_output

c----------------------------------------------------------------------
C!F77
c
c!Description:
C     Routine for reading and putting MODIS cloud mask ancillary
c     data values into arrays for processing.
C
c!Input parameters:
c cube          Current scan cube number
c h_eco1        1 km Ecosystem file handle number
c h_eco2        10 minute Ecosystem file handle number
c cube          Scan cube number
c buf_anc_size  Line and element granule size
C max_pixel     Maximum number of pixels per scan in granule
c rlat          Array containing scan cube of 1km pixel latitude
c               values
c rlon          Array containing scan cube of 1km pixel longitude
c               values
c TAI_end       Ending TAI time of granule
c
c!Output Parameters:
c
c cube_pw       Array containing scan cube of 1km pixel precipitable
c               water values (g/cm2)
c cube_eco      Array containing scan cube of 1km pixel ecosystem
c               values (Olson World Ecosystem Map)
c cube_ice      Array containing scan cube of 1km pixel sea ice
c               concentration values (fraction)
c cube_snow     Array containing scan cube of 1km pixel snow cover
c               values
c cube_sfctmp   Array containing scan cube of 1km pixel surface
c               temperatures - from model (K)
c cube_msl      Array containing scan cube of 1km pixel mean sea level
c               pressure values from model (hPa)
c cube_ugrd     Array containing scan cube of 1km pixel surface wind
c               u component (m/s)
c cube_vgrd     Array containing scan cube of 1km pixel surface wind
c               v component (m/s)
c cube_sst      Array containing scan cube of 1km pixel surface
c               temperatures - from Reynolds Blended SSTs (K)
c qa_bits       10 byte array containing qa bit results
c
c!Revision History:
c
c!Team-Unique Header:
c This software is developed by the MODIS Science Data Support
c Team for the National Aeronautics and Space Administration,
c Goddard Space Flight Center, under contract NAS5-32373.
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c       message              Error messaging subroutine
c       get_ancillary        Subroutine to read ancillary data
c       date_to_ccsds        Converts date from CCSDS to format
c                            passed through get_ancillary
c       TAItoUTC             PGS routine to convert TAI time to CCSDS
c
c!END
c----------------------------------------------------------------------

        include 'global.inc'
        include 'PGS_MODIS_39500.f'

c       scalar arguements
        integer h_eco1,h_eco2,cube,max_pixel

c       array arguments
        integer buf_anc_size(2),cube_snow(npixel,scans_cube)
        double precision TAI_end
        real rlat(npixel,scans_cube),
     +       rlon(npixel,scans_cube),cube_pw(npixel,scans_cube),
     +       cube_ice(npixel,scans_cube),cube_sfctmp(npixel,scans_cube),
     +       cube_ugrd(npixel,scans_cube),cube_vgrd(npixel,scans_cube),
     +       cube_msl(npixel,scans_cube),cube_sst(npixel,scans_cube)
        byte cube_eco(npixel,scans_cube),qa_bits(10)

c       local scalars
        real b_lat,b_lon,land,prmsl,pwat,ugrd,vgrd,ozone,icec,
     +       sst,landtmp,lat_temp,lon_temp
        integer lat_indx,lon_indx,bytloc,ii,jj,mm,nn,h_output,debug,
     +          i,j,newlin_eco2,iocheck,newele_eco2,nise,nise_out
        logical init

c       local arrays
        real pres(0:25),temp(0:25),mixr(0:25),plat(1,1),plon(1,1),
     +       pxl_bias(5)
        double precision delta
        integer met_date(4), ozn_date(4), ice_date(4), sst_date(4),
     +          nise_date(4)
        byte map_eco(2),madison_eco,ecotest(1,1)

c       external subroutines
        external get_ancillary,message

c ...   external functions
        integer compare_times
        external compare_times

c ...   data statement - time window for use of nise data (7 days in
c ...    seconds)
        data delta  /604800d0/
        data init / .true. /

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(10x/,''Subroutine get_anc_data '',/)')
      endif
c ................................................................

c ...   Grab the latitude information out of the geolocation file
c ...   First initialize output paramaters
         iocheck = 0
         plat(1,1) = 0.0
         plon(1,1) = 0.0
         ecotest(1,1) = 0

         do 134 mm = 1 , buf_anc_size(2)
         do 135 nn = 1 , buf_anc_size(1)
            cube_pw(nn,mm) = 0.0
            cube_eco(nn,mm) = 0
            cube_ice(nn,mm) = 0.0
            cube_snow(nn,mm) = 0
            cube_sfctmp(nn,mm) = 0.0
            cube_ugrd(nn,mm) = 0.0
            cube_vgrd(nn,mm) = 0.0
            cube_msl(nn,mm) = 0.0
            cube_sst(nn,mm) = 0.0
  135    continue
  134    continue

c----------------------------------------------------------------------
         if (init) then

           nise_out = 0
c ...      Use Mr. Gumley's ancillary data subroutine to extract
c ...      all parameters needed for atmosphere group at the
c ...      given lat and lon.

c ...      Do first call using one lat and lon to check time and
c ...      date of data
           lat_temp = rlat(1,1)
           lon_temp = rlon(1,1)
           call get_ancillary(lat_temp,lon_temp,pres,temp,
     +          mixr,land,landtmp,prmsl,pwat,ugrd,vgrd,ozone,icec,sst,
     +          nise,met_date,ozn_date,ice_date,sst_date,nise_date,
     +          pxl_bias)

c ...      Use function to compare tai time of granule to that of
c ...       of nise data set - see if it is within delta
           nise_out = compare_times(nise_date,TAI_end,delta)
           if (nise_out .ne. 0) then
             call modis_smf_setdynamicmsg( MODIS_W_GENERIC,
     +       'NISE data is either missing or time ' //
     +       'of NISE data is outside acceptable range. ' //
     +       'Processing will continue without use of NISE. ',
     +       'get_anc_data')
           endif

         endif

c ...    Now extract a scan cube's worth of data
         do 140 ii = 1 , scans_cube
            do 145 jj = 1 , max_pixel
               call get_ancillary(rlat(jj,ii),rlon(jj,ii),pres,temp,
     +         mixr,land,landtmp,prmsl,pwat,ugrd,vgrd,ozone,icec,sst,
     +         nise,met_date,ozn_date,ice_date,sst_date,nise_date,
     +         pxl_bias)

              cube_ice(jj,ii) = icec
              if (nise_out .eq. 0) then
                cube_snow(jj,ii) = nise
              end if
              cube_pw(jj,ii) = pwat
              cube_sfctmp(jj,ii) = landtmp
              cube_sst(jj,ii) = sst
              cube_ugrd(jj,ii) = ugrd
              cube_vgrd(jj,ii) = vgrd
              cube_msl(jj,ii) = prmsl
 145       continue
 140      continue


c----------------------------------------------------------------------
c ...    Section added to set QA bits if data extracted correctly
c ...     using ancillary data subroutine
c ...     Set qa_bits for NCEP data used
         if (landtmp .lt. 0.0) then
c ...      Preliminary logic for backup, and setting of qa bits
c ...      All NMC variables - First land surface temp not used
           call set_qa_bit(qa_bits,58)
           call set_qa_bit(qa_bits,59)
c ...      Surface winds - not used
           call set_qa_bit(qa_bits,62)
           call set_qa_bit(qa_bits,63)
c ...      Precipitable water - not used
           call set_qa_bit(qa_bits,73)
           call set_qa_bit(qa_bits,74)
        endif

c ...   Ice ancillary data
        if (icec .lt. 0.0) then
c ...      Preliminary logic for backup, and setting of qa bits
           call set_qa_bit(qa_bits,68)
           call set_qa_bit(qa_bits,69)
        else
c ...      SSMI product
           call set_qa_bit(qa_bits,68)
        endif

c ...   Snow mask - if missing value, else set SSMI use bit
        if (nise .lt. 0 .or. nise_out .ne. 0) then
           call set_qa_bit(qa_bits,66)
           call set_qa_bit(qa_bits,67)
        else
           call set_qa_bit(qa_bits,66)
        endif

c ...   Reynolds blended sst data
        if (sst .lt. 0.0) then
c ...      Reynolds blended SST - not used
           call set_qa_bit(qa_bits,60)
           call set_qa_bit(qa_bits,61)
        endif
c----------------------------------------------------------------------

c----------------------------------------------------------------------
c        Okay, now read the ecosystem information

         if (init) then
c ...      First, let's make sure we are reading from the correct
c ...       file.  Extract Madison Wisconsin eco_type from file
           if (h_eco1 .ne. -5555) then
             plat(1,1) = 43.083
             plon(1,1) = -89.305
             call read_goode( 1,1,plat,plon,h_eco1,ecotest )
             if ( ecotest(1,1) .ne. 22 ) call message( 'get_anc_data',
     &         'Extracted incorrect ecosystem value from 1 km file. ' //
     &         '[OPERATOR ACTION: Verify size and format of file. ' //
     &         'If error persists, contact SDST.]', 0, 2 )
           else
c ...        Read from the 10 minute global Olson ecosystem file
             b_lat = 43.083
             b_lon = -89.305
             lat_indx = int((90.0 - b_lat) * 6.0 + 1.0)
             lon_indx = int((b_lon + 180.0) * 6.0 + 1.0)
             if(lat_indx .gt. 1080) lat_indx = 1080
             if(lon_indx .gt. 2160) lon_indx = 2160
             bytloc = ((lat_indx-1) * 2160) + lon_indx
             newlin_eco2 = bytloc / 2 + 1
             newele_eco2 = mod(bytloc, 2)
             if (mod( bytloc,2) .eq. 0) then
               newlin_eco2 = newlin_eco2 - 1
               newele_eco2 = 2
             endif

c ...        read value from 10minute file
             read (h_eco2, rec=newlin_eco2, iostat=iocheck) map_eco
             if (iocheck .ne. 0) then
               call message( 'get_anc_data',
     &         'Error reading ecosystem value from file. ' //
     &         'Make sure correct 10 minute file is loaded. ' //
     &         ' [OPERATOR ACTION: Notify SDST.]',
     &         0, 2 )
             endif

             madison_eco = map_eco(newele_eco2)
c ...        Compare value extracted with known correct value
             if (madison_eco .ne. 55) then
               call message( 'get_anc_data',
     &         'Extracted incorrect ecosystem value from file. ' //
     &         'Make sure correct 10 minute file is loaded. ' //
     &         ' [OPERATOR ACTION: Notify SDST.]',
     &         0, 2 )
             endif
           endif

           init = .false.
         endif
c -------------------------------------------------------------------

         if (h_eco1 .ne. -5555) then
c ...      Read all ecosystem values for this scan cube

           call read_goode( buf_anc_size( 1 ), buf_anc_size( 2 ),
     &       rlat, rlon, h_eco1, cube_eco )
         else

c ...      Read out of 10 minute global file
           do 170 ii = 1 , scans_cube
             do 175 jj = 1 , max_pixel
                b_lat = rlat(jj,ii)
                b_lon = rlon(jj,ii)
                lat_indx = int((90.0 - b_lat) * 6.0 + 1.0)
                lon_indx = int((b_lon + 180.0) * 6.0 + 1.0)
                if(lat_indx .gt. 1080) lat_indx = 1080
                if(lon_indx .gt. 2160) lon_indx = 2160
                bytloc = ((lat_indx-1) * 2160) + lon_indx
                newlin_eco2 = bytloc / 2 + 1
                newele_eco2 = mod(bytloc, 2)
                if (mod( bytloc,2) .eq. 0) then
                  newlin_eco2 = newlin_eco2 - 1
                  newele_eco2 = 2
                endif

c ...           read value from file
                read (h_eco2, rec=newlin_eco2, iostat=iocheck) map_eco
                if (iocheck .ne. 0) then
                  call message( 'get_anc_data',
     &            'Error reading ecosystem value from file. ' //
     &            ' [OPERATOR ACTION: Notify SDST.]',
     &            0, 2 )
                endif
                cube_eco(jj,ii) = map_eco(newele_eco2)
 175         continue
 170       continue

         endif

c----------------------------------------------------------------------

c ... debug statement ............................................
      if (debug .gt. 3) then
        write(h_output,'(15x,'' ANCILLARY DATA FOR CUBE: '',I5)') cube
        write(h_output,'(2x,'' Ice concentraion'',/,20(5f10.4/))')
     +       ((cube_ice(j,i),j=1,10),i=1,10)
        write(h_output,'(2x,'' Snow from NISE'',/,20(5I10/))')
     +       ((cube_snow(j,i),j=1,10),i=1,10)
        write(h_output,'(2x,''Precipitable Water (g/cm**2)'',
     +       /,20(5f10.4/))') ((cube_pw(j,i),j=1,10),i=1,10)
        write(h_output,'(2x,'' Ecosystem Type '',
     +       /,20(5I10/))') ((cube_eco(j,i),j=1,10),i=1,10)
        write(h_output,'(2x,'' Surface Temperature (K)'',
     +       /,20(5f10.4/))') ((cube_sfctmp(j,i),j=1,10),i=1,10)
        write(h_output,'(2x,'' U Wind Component (m/s) '',
     +       /,20(5f10.4/))') ((cube_ugrd(j,i),j=1,10),i=1,10)
        write(h_output,'(2x,'' V Wind Component (m/s) '',
     +       /,20(5f10.4/))') ((cube_vgrd(j,i),j=1,10),i=1,10)
        write(h_output,'(2x,'' Mean Sea Level Pressure (hPa) '',
     +       /,20(5f10.4/))') ((cube_msl(j,i),j=1,10),i=1,10)
      endif
c ................................................................

         return
         end
