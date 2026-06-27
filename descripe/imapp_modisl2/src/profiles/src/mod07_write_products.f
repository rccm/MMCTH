      SUBROUTINE MOD07_WRITE_PRODUCTS( OUT_HANDLE, SCAN, NBOXES )

c-----------------------------------------------------------------------
c !F77
c
c !DESCRIPTION:
c     Write MOD07 retrieval products and QA data to the output file.
c
c !INPUT PARAMETERS:
c     OUT_HANDLE    Integer array returned by OPMFIL for output file
c     SCAN          Scan number within L1B granule
c     NBOXES        Number of retrieval boxes along the scan
c                   (e.g. 1354 pixels, 5x5 sampling gives 1354/5 boxes)
c
c !OUTPUT PARAMETERS:
c     None
c
c     The following arrays in COMMON /MOD07_DATA/ are written to output:
c     BRIGHTNESS_TEMP       Brightness temperature
c     GUESS_TEMP_PROFILE    Guess temperature profile
c     GUESS_WVMR_PROFILE    Guess water vapour mixing ratio profile
c     RETR_TEMP_PROFILE     Retrieved temperature profile
c     RETR_DEWP_PROFILE     Retrieved dewpoint profile
c     RETR_WVMR_PROFILE     Retrieved water vapour mixing ratio profile
c     RETR_HITE_PROFILE     Retrieved geopotential height profile
c     RETR_OZONE_PROFILE     Retrieved ozone profile
c     SFC_TEMP              Surface temperature
c     SFC_PRES              Surface pressure
c     SFC_ELEV              Surface elevation
c     HEIGHT_TROPOPAUSE     Height of the tropopause
c     WATER_VAPOR           Water vapor for entire column, by integrating regression profiles 
c     WATER_VAPOR_DIRECT    Water vapor for entire column, by direct regression
c     WATER_VAPOR_LOW       Water vapor for surface to 900 hPa
c     WATER_VAPOR_HIGH      Water vapor for 700 hPa to 300 hPa
c     TOTAL_OZONE           Total column ozone
c     TOTAL_TOTALS          Total totals index
c     LIFTED_INDEX          Lifted index
c     K_INDEX               K index
c     PRODUCT_QA            Product run time QA
c     WATER_VAPOR_QA        IR water vapor product QA
c
c !REVISION HISTORY:
c    26 February 2002: SWS included Water_Vapor_Direct
c    May 2006: EvaBorbas added ozone profiles to the outputs
c    March 2007 EvaBorbas changed the Surface_Temperature name to Skin_Temperature
c    Jan 2009   mixing ration and dewpoint temperature are both output values
c    Jan 2010: EvaBorbas removed  ozone profiles to the outputs

c !TEAM-UNIQUE HEADER:
c     Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
c
c !END
c-----------------------------------------------------------------------

      IMPLICIT NONE

      INCLUDE 'mod07_data.inc'
      INCLUDE 'mapi.inc'

c ... arguments

      INTEGER out_handle( MODFILLEN ), scan, nboxes

c ... local variables

      INTEGER i

c ... Write each product array for this scan to output file

      i = 1
      call write_output( out_handle, scan, nboxes, i,
c     &  sfc_temp, 'Surface_Temperature' )
     &  sfc_temp, 'Skin_Temperature' )

      call write_output( out_handle, scan, nboxes, i,
     &  sfc_pres, 'Surface_Pressure' )

      call write_output( out_handle, scan, nboxes, i,
     &  sfc_elev, 'Surface_Elevation' )

      call write_output( out_handle, scan, nboxes, i,
     &  height_tropopause, 'Tropopause_Height' )

      call write_output( out_handle, scan, nboxes, i,
     &  total_ozone, 'Total_Ozone' )

      call write_output( out_handle, scan, nboxes, i,
     &  total_totals, 'Total_Totals' )

      call write_output( out_handle, scan, nboxes, i,
     &  lifted_index, 'Lifted_Index' )

      call write_output( out_handle, scan, nboxes, i,
     &  k_index, 'K_Index' )

      call write_output( out_handle, scan, nboxes, i,
     &  water_vapor, 'Water_Vapor' )

      call write_output( out_handle, scan, nboxes, i,
     &  water_vapor_direct, 'Water_Vapor_Direct' )

      call write_output( out_handle, scan, nboxes, i,
     &  water_vapor_low, 'Water_Vapor_Low' )

      call write_output( out_handle, scan, nboxes, i,
     &  water_vapor_high, 'Water_Vapor_High' )

      do i = 1, nbands

        call write_output( out_handle, scan, nboxes, i,
     &    brightness_temp(1,1,i), 'Brightness_Temperature' )

      end do

      do i = 1, outlevels

        call write_output( out_handle, scan, nboxes, i,
     &    guess_temp_profile(1,1,i), 'Guess_Temperature_Profile' )

        call write_output( out_handle, scan, nboxes, i,
     &    guess_wvmr_profile(1,1,i), 'Guess_Moisture_Profile' )

        call write_output( out_handle, scan, nboxes, i,
     &    retr_temp_profile(1,1,i), 'Retrieved_Temperature_Profile' )

        call write_output( out_handle, scan, nboxes, i,
     &    retr_dewp_profile(1,1,i), 'Retrieved_Moisture_Profile' )

       call write_output( out_handle, scan, nboxes, i,
     &    retr_wvmr_profile(1,1,i), 'Retrieved_WV_Mixing_Ratio_Profile' )


        call write_output( out_handle, scan, nboxes, i,
     &    retr_hite_profile(1,1,i), 'Retrieved_Height_Profile' )

c        call write_output( out_handle, scan, nboxes, i,
c     &    retr_ozone_profile(1,1,i), 'Retrieved_Ozone_Profile' )

      end do

c ... Write all QA arrays for this scan to output file

      call write_qa( out_handle, scan, nboxes )

      END

c-----------------------------------------------------------------------

      subroutine write_qa( out_handle, scan, nboxes )

c-----------------------------------------------------------------------
c !F77
c
c !DESCRIPTION:
c     Write QA data to the output file.
c
c !INPUT PARAMETERS:
c     OUT_HANDLE    Integer array returned by OPMFIL for output file
c     SCAN          Scan number within L1B granule
c     NBOXES        Number of retrieval boxes along the scan
c                   (e.g. 1354 pixels, 5x5 sampling gives 1354/5 boxes)
c
c !OUTPUT PARAMETERS:
c     None
c
c !REVISION HISTORY:
c
c !TEAM-UNIQUE HEADER:
c     Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
c
c !END
c-----------------------------------------------------------------------

      implicit none

      include 'mod07_data.inc'
      include 'mapi.inc'

c ... arguments

      integer out_handle( MODFILLEN ), scan, nboxes

c ... local variables

      integer out, i, j, k, rtn, start( 3 ), dims( 3 )

      byte byte_data( max_samp_pixel * max_samp_line * nproduct_qa )

      character*32 arrnm, grpnm, text

c ... Create the cloudmask output array

      do i = 1, max_samp_pixel * max_samp_line * nproduct_qa
        byte_data( i ) = 0
      end do

      out = 1
      do j = 1, max_samp_line
        do i = 1, nboxes
          byte_data( out ) =
     &      mask1( 1, 1+isamp/2+(i-1)*isamp, 1+isamp/2+(j-1)*isamp )
          out = out + 1
        end do
      end do

c ... Write the cloudmask array to the output file

      start( 1 ) = 0
      start( 2 ) = ( scan - 1 ) * max_samp_line
      dims( 1 ) = nboxes
      dims( 2 ) = 2
      arrnm = 'Cloud_Mask'
      grpnm = ' '
      rtn = -1
      rtn = pmar( out_handle, arrnm, grpnm, start, dims, byte_data )
      if ( rtn .ne. 0 ) then
        write( text, '(''Write failed on output array '',a)') arrnm
        call message ( 'mod07_write_products', text //
     &  ' [OPERATOR ACTION: Check available disk space.' //
     &  ' If problem persists, contact SDST]', 0, 2 )
      endif

c ... Create the product QA output array

      do i = 1, max_samp_pixel * max_samp_line * nproduct_qa
        byte_data( i ) = 0
      end do

      out = 1
      do j = 1, max_samp_line
        do i = 1, nboxes
          do k = 1, nproduct_qa
            byte_data( out ) = product_qa( k, i, j )
            out = out + 1
          end do
        end do
      end do

c ... Write the product QA array to the output file

      start( 1 ) = 0
      start( 2 ) = 0
      start( 3 ) = ( scan - 1 ) * max_samp_line
      dims( 1 ) = nproduct_qa
      dims( 2 ) = nboxes
      dims( 3 ) = max_samp_line
      arrnm = 'Quality_Assurance'
      grpnm = ' '
      rtn = -1
      rtn = pmar( out_handle, arrnm, grpnm, start, dims, byte_data )
      if ( rtn .ne. 0 ) then
        write( text, '(''Write failed on output array '',a)') arrnm
        call message ( 'mod07_write_products', text //
     &  ' [OPERATOR ACTION: Check available disk space.' //
     &  ' If problem persists, contact SDST]', 0, 2 )
      endif

c ... Create the IR water vapor QA output array

      do i = 1, max_samp_pixel * max_samp_line * nproduct_qa
        byte_data( i ) = 0
      end do

      out = 1
      do j = 1, max_samp_line
        do i = 1, nboxes
          do k = 1, nwater_vapor_qa
            byte_data( out ) = water_vapor_qa( k, i, j )
            out = out + 1
          end do
        end do
      end do

c ... Write the IR water vapor QA array to the output file

      start( 1 ) = 0
      start( 2 ) = 0
      start( 3 ) = ( scan - 1 ) * max_samp_line
      dims( 1 ) = nwater_vapor_qa
      dims( 2 ) = nboxes
      dims( 3 ) = max_samp_line
      arrnm = 'Quality_Assurance_Infrared'
      grpnm = ' '
      rtn = -1
      rtn = pmar( out_handle, arrnm, grpnm, start, dims, byte_data )
      if ( rtn .ne. 0 ) then
        write( text, '(''Write failed on output array '',a)') arrnm
        call message ( 'mod07_write_products', text //
     &  ' [OPERATOR ACTION: Check available disk space.' //
     &  ' If problem persists, contact SDST]', 0, 2 )
      endif

      end

c-----------------------------------------------------------------------

      subroutine write_output( out_handle, scan, nboxes, level, array,
     &  arrnm )

c Purpose:
c    Write a REAL*4 product array for one MODIS scan to the output file.
c
c Input:
c    OUT_HANDLE    MAPI handle array returned by OPMFIL for output file
c    SCAN          MODIS scan number
c    NBOXES        Number of 5x5 retrieval boxes for this scan
c    LEVEL         Level number in output array
c                  (=1 for 2D arrays)
c    ARRAY         Array of product values
c                  (Bad values should be set to BAD_VALUE in the calling
c                  routine: see mod06uw_data.inc for BAD_VALUE value)
c    ARRNM         Name of the array (SDS name) in the output file
c
c Output:
c    None
c
c Revised:
c    12-DEC-1997 Kathy Strabala
c                Changed BYTE to INTEGER*1 so that mapi would
c                like it
c    10-DEC-1997 Liam Gumley
c                Added BYTE data types
c                Now gets _FillValue directly from output file

      implicit none

      include 'mod07_data.inc'
      include 'mapi.inc'

c ... arguments

      integer out_handle( MODFILLEN ), scan, nboxes, level
      real array( max_samp_pixel, max_samp_line )
      character*(*) arrnm

c ... local variables

      character*20 grpnm, dtype, attr
      double precision scale_factor, add_offset
      integer rtn, rank, start( 3 ), dims( 3 ), nms

      integer i, j, out

      real float_data( max_samp_pixel * max_samp_line )
      integer*2 short_data( max_samp_pixel * max_samp_line )
      byte byte_data( max_samp_pixel * max_samp_line )

      character*72 text

      character*20 fill_type
      real fill_float
      integer*2 fill_short
      byte fill_byte

c ... external functions

      external gmardm, gmarin, pmar

c ... Get the scale factor

      grpnm = ' '
      attr = 'scale_factor'
      dtype = 'REAL*8'
      nms = 1
      rtn = -1
      rtn = gmarin( out_handle, arrnm, grpnm, attr, dtype,
     &  nms, scale_factor )
      if ( rtn .ne. 0 ) scale_factor = 1.0d0

c ... Get the offset

      grpnm = ' '
      attr = 'add_offset'
      dtype = 'REAL*8'
      nms = 1
      rtn = -1
      rtn = gmarin( out_handle, arrnm, grpnm, attr, dtype,
     &  nms, add_offset )
      if ( rtn .ne. 0 ) add_offset = 0.0d0

c ... Get the data type

      grpnm = ' '
      rank = 3
      rtn = -1
      rtn = gmardm( out_handle, arrnm, grpnm, dtype, rank, dims )
      if ( rtn .ne. 0 ) then
        write( text, '(''Could not find output array '',a)' ) arrnm
        call message ( 'mod07_write_products',
     &    text // ' [OPERATOR ACTION: Contact SDST]', 0, 2 )
      endif
      if ( dtype(1:6) .ne. 'REAL*4' .and.
     &     dtype(1:9) .ne. 'INTEGER*2' .and.
     &     dtype(1:9) .ne. 'INTEGER*1' ) then
        write( text, '(''Cannot handle type '',a,'' for array '',a)' )
     &     dtype, arrnm
        call message ( 'mod07_write_products',
     &    text // ' [OPERATOR ACTION: Contact SDST]', 0, 2 )
      endif

c ... Get the fill value

      fill_float = 0.0
      fill_short = 0
      fill_byte = 0
      grpnm = ' '
      attr = '_FillValue'
      fill_type = dtype
      nms = 1
      rtn = -1
      if ( fill_type(1:6) .eq. 'REAL*4' ) rtn = gmarin( out_handle,
     &  arrnm, grpnm, attr, dtype, nms, fill_float )
      if ( fill_type(1:9) .eq. 'INTEGER*2' ) rtn = gmarin( out_handle,
     &  arrnm, grpnm, attr, dtype, nms, fill_short )
      if ( fill_type(1:9) .eq. 'INTEGER*1' ) rtn = gmarin( out_handle,
     &  arrnm, grpnm, attr, dtype, nms, fill_byte )
      if ( rtn .ne. 0 ) then
        write( text, '(''_FillValue not found for array '',a)' ) arrnm
        call message ( 'mod07_write_products',
     &    text // ' [OPERATOR ACTION: Contact SDST]', 0, 2 )
      endif

c ... Create the scaled output array

      out = 1

      do j = 1, max_samp_line

        do i = 1, nboxes

c ...     Set default output value (fill value)

          float_data( out ) = fill_float
          short_data( out ) = fill_short
          byte_data( out ) = fill_byte

c ...     If array value is good, scale it and store it

          if ( abs( array( i, j ) - bad_value ) .gt. 1.0e-5 ) then

            if ( dtype(1:6) .eq. 'REAL*4' ) float_data( out ) =
     &        array( i, j ) / real( scale_factor ) + real( add_offset )

            if ( dtype(1:9) .eq. 'INTEGER*2' ) short_data( out ) =
     &        nint( array( i, j ) / real( scale_factor ) +
     &        real( add_offset ) )

            if ( dtype(1:9) .eq. 'INTEGER*1' ) byte_data( out ) =
     &        nint( array( i, j ) / real( scale_factor ) +
     &        real( add_offset ) )

          endif

          out = out + 1

        end do

      end do

c ... Write the scaled output array to the output file

      start( 1 ) = 0
      start( 2 ) = ( scan - 1 ) * max_samp_line
      start( 3 ) = level - 1
      dims( 1 ) = nboxes
      dims( 2 ) = 2
      dims( 3 ) = 1
      grpnm = ' '
      rtn = -1

      if ( dtype(1:6) .eq. 'REAL*4' )
     &  rtn = pmar( out_handle, arrnm, grpnm, start, dims, float_data )

      if ( dtype(1:9) .eq. 'INTEGER*2' )
     &  rtn = pmar( out_handle, arrnm, grpnm, start, dims, short_data )

      if ( dtype(1:9) .eq. 'INTEGER*1' )
     &  rtn = pmar( out_handle, arrnm, grpnm, start, dims, byte_data )

      if ( rtn .ne. 0 ) then
        write( text, '(''Write failed on output array '',a)') arrnm
        call message ( 'mod07_write_products', text //
     &  ' [OPERATOR ACTION: Check available disk space.' //
     &  ' If problem persists, contact SDST]', 0, 2 )
      endif

      end
