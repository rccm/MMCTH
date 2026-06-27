      subroutine mod06_write_products( out_handle, scan, nboxes )

c-----------------------------------------------------------------------
c
c!F77
c
c!Description:
c    Write MOD06 retrieval products and QA data to the output file.
c
c!Input Parameters:
c    OUT_HANDLE    Integer array returned by OPMFIL for output file
c    SCAN          Scan number within L1B granule
c    NBOXES        Number of retrieval boxes along the scan
c                  (e.g. 1354 pixels, 5x5 sampling gives 1354/5 boxes)
c
c!Output Parameters:
c    None
c
c    The following arrays in COMMON /MOD06_DATA/ are written to output:
c    BRIGHTNESS_TEMP       Brightness temperature
c    SFC_TEMP              Surface temperature
c    SFC_PRES              Surface pressure
c    PRODUCT_QA            Product run time QA
c
c!Revision History:
c $Id: mod06_write_products.f,v 1.7 1999/06/11 22:35:23 kis Exp $
c
c!Team-unique Header:
c
c!End
c
c-----------------------------------------------------------------------
      
      implicit none


      include 'mod06uw_data.inc'
      include 'mapi.inc'
      
c ... arguments

      integer out_handle( MODFILLEN ), scan, nboxes

c ... local variables

      integer i , j , vband

c whuang 05/07/01: Added save statement for common / mod06_data /
      save / mod06_data / 

c ... local arrays
      real ctp_night(max_samp_pixel,max_samp_line),
     +     ctp_night_nnad(max_samp_pixel,max_samp_line),
     +     ctp_day(max_samp_pixel,max_samp_line),
     +     ctp_day_nnad(max_samp_pixel,max_samp_line),
     +     cth_night_nnad(max_samp_pixel,max_samp_line),
     +     cth_day_nnad(max_samp_pixel,max_samp_line),
     +     ctt_night(max_samp_pixel,max_samp_line),
     +     ctt_night_nnad(max_samp_pixel,max_samp_line),
     +     ctt_day(max_samp_pixel,max_samp_line),
     +     ctt_day_nnad(max_samp_pixel,max_samp_line),
     +     cee_night(max_samp_pixel,max_samp_line),
     +     cee_night_nnad(max_samp_pixel,max_samp_line),
     +     cee_day(max_samp_pixel,max_samp_line),
     +     cee_day_nnad(max_samp_pixel,max_samp_line),
     +     cld_irp_night(max_samp_pixel,max_samp_line),
     +     cld_irp_day(max_samp_pixel,max_samp_line),
     +     fraction_night(max_samp_pixel,max_samp_line),
     +     fraction_night_nnad(max_samp_pixel,max_samp_line),
     +     fraction_day(max_samp_pixel,max_samp_line),
     +     fraction_day_nnad(max_samp_pixel,max_samp_line)
     

c ... Initialize day/night arrays
      do j = 1 , max_samp_line
         do i = 1 , max_samp_pixel
            ctp_night(i,j) = bad_value
            ctp_night_nnad(i,j) = bad_value
            ctt_night(i,j) = bad_value
            ctt_night_nnad(i,j) = bad_value
            cth_night_nnad(i,j) = bad_value
            cee_night(i,j) = bad_value
            cee_night_nnad(i,j) = bad_value
            cld_irp_night(i,j) = bad_value
            fraction_night(i,j) = bad_value
            fraction_night_nnad(i,j) = bad_value
            ctp_day(i,j) = bad_value
            ctp_day_nnad(i,j) = bad_value
            cth_day_nnad(i,j) = bad_value
            ctt_day(i,j) = bad_value 
            ctt_day_nnad(i,j) = bad_value 
            cee_day(i,j) =  bad_value
            cee_day_nnad(i,j) =  bad_value
            cld_irp_day(i,j) = bad_value
            fraction_day(i,j) = bad_value
            fraction_day_nnad(i,j) = bad_value
         enddo
      enddo

c ..  Separate 4 arrays into day and night products

      do j = 1 , max_samp_line
         do i = 1 , max_samp_pixel
            if (day_night_flag(i,j) .eq. 0) then
                ctp_night(i,j) = cloudtop_pres(i,j)
                ctp_night_nnad(i,j) = cloudtop_pres_nearnad(i,j)
                cth_night_nnad(i,j) = cloudtop_height_nearnad(i,j)
                ctt_night(i,j) = cloudtop_temp(i,j)
                ctt_night_nnad(i,j) = cloudtop_temp_nearnad(i,j)
                cee_night(i,j) = cloudtop_eff_emi(i,j)
                cee_night_nnad(i,j) = cloudtop_eff_emi_nearnad(i,j)
                cld_irp_night(i,j) = cloud_phase(i,j)
                fraction_night(i,j) = cloud_fraction(i,j)
                fraction_night_nnad(i,j) = cloud_fraction_nearnad(i,j)
            elseif (day_night_flag(i,j) .eq. 1) then
                ctp_day(i,j) = cloudtop_pres(i,j)
                ctp_day_nnad(i,j) = cloudtop_pres_nearnad(i,j)
                cth_day_nnad(i,j) = cloudtop_height_nearnad(i,j)
                ctt_day(i,j) = cloudtop_temp(i,j)
                ctt_day_nnad(i,j) = cloudtop_temp_nearnad(i,j)
                cee_day(i,j) = cloudtop_eff_emi(i,j)
                cee_day_nnad(i,j) = cloudtop_eff_emi_nearnad(i,j)
                cld_irp_day(i,j) = cloud_phase(i,j)
                fraction_day(i,j) = cloud_fraction(i,j)
                fraction_day_nnad(i,j) = cloud_fraction_nearnad(i,j)
            endif
         enddo
      enddo
c ... Write each product array for this scan to output file

      do i = 1 , n_output
        call write_output( out_handle, scan, nboxes, i,
     &    brightness_temp(1,1,i), 'Brightness_Temperature' )
      end do

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  sfc_temp, 'Surface_Temperature' )

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  sfc_pres, 'Surface_Pressure' )

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  cloud_h_method, 'Cloud_Height_Method' )

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  cloudtop_pres, 'Cloud_Top_Pressure' )

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  cloudtop_pres_nearnad, 'Cloud_Top_Pressure_Nadir' )

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  cloudtop_height, 'Cloud_Top_Height' )

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  cloudtop_height_nearnad, 'Cloud_Top_Height_Nadir' )

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  cth_day_nnad, 'Cloud_Top_Height_Nadir_Day' )

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  cth_night_nnad, 'Cloud_Top_Height_Nadir_Night' )

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  ctp_night, 'Cloud_Top_Pressure_Night' )

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  ctp_night_nnad, 'Cloud_Top_Pressure_Nadir_Night' )

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  ctp_day, 'Cloud_Top_Pressure_Day' )

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  ctp_day_nnad, 'Cloud_Top_Pressure_Nadir_Day' )

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  cloudtop_temp, 'Cloud_Top_Temperature' )

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  cloudtop_temp_nearnad, 'Cloud_Top_Temperature_Nadir' )

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  ctt_night, 'Cloud_Top_Temperature_Night' )

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  ctt_night_nnad, 'Cloud_Top_Temperature_Nadir_Night' )

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  ctt_day, 'Cloud_Top_Temperature_Day' )

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  ctt_day_nnad, 'Cloud_Top_Temperature_Nadir_Day' )

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  height_tropopause, 'Tropopause_Height' )

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  cloud_fraction, 'Cloud_Fraction' )

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  cloud_fraction_nearnad, 'Cloud_Fraction_Nadir' )

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  fraction_night, 'Cloud_Fraction_Night' )

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  fraction_night_nnad, 'Cloud_Fraction_Nadir_Night' )

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  fraction_day, 'Cloud_Fraction_Day' )

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  fraction_day_nnad, 'Cloud_Fraction_Nadir_Day' )

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  cloudtop_eff_emi, 'Cloud_Effective_Emissivity' )

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  cloudtop_eff_emi_nearnad, 'Cloud_Effective_Emissivity_Nadir' )

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  cee_night, 'Cloud_Effective_Emissivity_Night' )

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  cee_night_nnad, 'Cloud_Effective_Emissivity_Nadir_Night' )

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  cee_day, 'Cloud_Effective_Emissivity_Day' )

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  cee_day_nnad, 'Cloud_Effective_Emissivity_Nadir_Day' )

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  cloudtop_pre_ir, 'Cloud_Top_Pressure_Infrared' )

      do i = 1, max_solutns
        call write_output( out_handle, scan, nboxes, i,
     &    spec_cloud_forcing(1,1,i), 'Spectral_Cloud_Forcing' )
        call write_output( out_handle, scan, nboxes, i,
     &    cloudtop_pres_from_ratios(1,1,i), 
     &    'Cloud_Top_Pressure_From_Ratios' )
      end do
            
c     do i = 1, n_output
c       call write_output( out_handle, scan, nboxes, i,
c    &    radiance_variance(1,1,i), 'Radiance_Variance' )
c     end do      
c     Output only for band 31 (RAF).
      i = 1
      vband = 2
      call write_output( out_handle, scan, nboxes, i,
     &    radiance_variance(1,1,vband), 'Radiance_Variance' )

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  cloud_phase, 'Cloud_Phase_Infrared' )

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  cld_irp_night, 'Cloud_Phase_Infrared_Night' )

      i = 1
      call write_output( out_handle, scan, nboxes, i,
     &  cld_irp_day, 'Cloud_Phase_Infrared_Day' )

c ... Write all QA arrays for this scan to output file

      call write_qa( out_handle, scan, nboxes )
                  
      end

c-----------------------------------------------------------------------

      subroutine write_qa( out_handle, scan, nboxes )

c-----------------------------------------------------------------------
c
c!F77
c
c!Description:
c    Write MOD06 QA data to the output file.
c
c!Input Parameters:
c    OUT_HANDLE    Integer array returned by OPMFIL for output file
c    SCAN          Scan number within L1B granule
c    NBOXES        Number of retrieval boxes along the scan
c                  (e.g. 1354 pixels, 5x5 sampling gives 1354/5 boxes)
c
c!Output Parameters:
c    None
c
c    The following arrays in COMMON /MOD06_DATA/ are written to output:
c    PRODUCT_QA            Product run time QA
c
c!Revision History:
c    12-DEC-1997 Liam Gumley, CIMSS/SSEC
c                Created
c
c    11-MAY-2010 Rich Frey, CIMSS/SSEC
c                Added 2nd byte to the Cloud_Mask_5km SDS
c
c!Team-unique Header:
c
c!End
c
c-----------------------------------------------------------------------
      
      
      implicit none

      
      include 'mod06uw_data.inc'
      include 'mapi.inc'
      
c ... arguments

      integer out_handle( MODFILLEN ), scan, nboxes

c ... local variables

      integer out, i, j, k, rtn, jline, jbox, start( 3 ), dims( 3 ),
     &        mask_byte, mask_test, mask_scene
      
      byte byte_data( max_samp_pixel * max_samp_line * nproduct_qa ),
     &     mask_valid, mask_glint, mask_snow, mask_day,
     &     mt
      
      character*32 arrnm, grpnm

      character*72 text

      logical glint, snow, land, water, desert, coast, day

c whuang 05/07/01: Added save statement for common / mod06_data /
      save / mod06_data /
           
c-----------------------------------------------------------------------
c ... Create the cloudmask output array
      
      do i = 1, max_samp_pixel * max_samp_line * nproduct_qa
        byte_data( i ) = 0
      end do
      
      out = 1
      do j = 1, max_samp_line
        do i = 1, nboxes

          jline = 1+isamp/2+(j-1)*isamp
          jbox = 1+isamp/2+(i-1)*isamp

c-----------------------------------------------------------------------
c         First byte.

          byte_data( out ) =
c    &      mask1( 1, 1+isamp/2+(i-1)*isamp, 1+isamp/2+(j-1)*isamp )
     &      mask1( 1, jbox, jline )

          out = out + 1

c-----------------------------------------------------------------------
c         Second byte

          mask_test = 1
          mask_byte = mask1( 1, jbox, jline )
          mt = mask_byte
          if(mask_byte .lt. 0) mask_byte = mask_byte + 256
          mask_valid = iand( mask_byte, mask_test )

          if( mask_valid .eq. 1 .and. cloudtop_pres(i,j) .ne. bad_value) then

            glint = .false.
            snow = .false.
            water = .false.
            land = .false.
            coast = .false.
            desert = .false.
            day = .true.

            mask_test = 8
            mask_day = iand( mask_byte, mask_test )
            if ( mask_day .lt. 8 ) day = .false.

            mask_test = 16
            mask_glint = iand( mask_byte, mask_test )
            if ( mask_glint .lt. 16 ) glint = .true.

            mask_test = 32
            mask_snow = iand( mask_byte, mask_test )
            if ( mask_snow .lt. 32 ) snow = .true.

            mask_test = 192
            mask_scene = iand( mask_byte, mask_test )
            if ( mask_scene .lt. 64 ) then
              water = .true.
            else if ( mask_scene .lt. 128 ) then 
              coast = .true.
            else if ( mask_scene .lt. 192 ) then 
              desert = .true.
            else if ( mask_scene .ge. 192 ) then
              land = .true.
            end if

            if(glint) then
              call set_qa_bit( byte_data(out), 1 )
            else
              call set_qa_bit( byte_data(out), 0 )
            end if
            if(snow) then
              call set_qa_bit( byte_data(out), 3 )
            else
              call set_qa_bit( byte_data(out), 2 )
            end if
            if(water) then
              call set_qa_bit( byte_data(out), 4 )
            else if(coast) then
              call set_qa_bit( byte_data(out), 5 )
            else if(desert) then
              call set_qa_bit( byte_data(out), 4 )
              call set_qa_bit( byte_data(out), 5 )
            else if(land) then
              call set_qa_bit( byte_data(out), 6 )
            else
              write(*,'(''Invalid scene type in write_qa'')')
            end if
            if(day) then
              call set_qa_bit( byte_data(out), 7 )
            end if

          end if

c         if(mt .lt. 0) then
c           write(*,'(f12.2,6l3,2i10)') cloudtop_pres(i,j), glint, snow, water, land, coast, desert,
c    *               mask_byte, mask_scene
c           write(*,'(3i10)') mask1(1,jbox,jline),byte_data(out-1),byte_data(out)
c           write(*,'(1x)')
c           read(*,*)
c         end if

          out = out + 1

c-----------------------------------------------------------------------

        end do
      end do

c ... Write the cloudmask array to the output file

      start( 1 ) = 0
      start( 2 ) = 0
      start( 3 ) = ( scan - 1 ) * max_samp_line
      dims( 1 ) = 2
      dims( 2 ) = nboxes
      dims( 3 ) = max_samp_line
      arrnm = 'Cloud_Mask_5km'
      grpnm = ' '
      rtn = -1
      rtn = pmar( out_handle, arrnm, grpnm, start, dims, byte_data )
      if ( rtn .ne. 0 ) then
        write( text, '(''Write failed on output array '',a)') arrnm
        call message( 'mod06_write_products', text //
     &   ' [OPERATOR ACTION: Check system resources]' ,
     &   0, 2 )
      endif

c-----------------------------------------------------------------------

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
      arrnm = 'Quality_Assurance_5km'
      grpnm = ' '
      rtn = -1
      rtn = pmar( out_handle, arrnm, grpnm, start, dims, byte_data )
      if ( rtn .ne. 0 ) then
        write( text, '(''Write failed on output array '',a)') arrnm
        call message( 'mod06_write_products', text //
     &   ' [OPERATOR ACTION: Check system resources]' ,
     &   0, 2 )
      endif

      end
      
c-----------------------------------------------------------------------

      subroutine write_output( out_handle, scan, nboxes, level, array,
     &  arrnm )
c-----------------------------------------------------------------------
c!F77
c
c!Description:
c    Write MOD06 retrieval products to the output file.
c
c!Input Parameters:
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
c!Output Parameters:
c    None
c
c!Revision History:
c    12-DEC-1997 Kathy Strabala 
c                Changed BYTE to INTEGER*1 so that mapi would 
c                like it
c    10-DEC-1997 Liam Gumley
c                Added BYTE data types
c                Now gets _FillValue directly from output file
c
c!Team-unique Header:
c
c!End
c-----------------------------------------------------------------------

      implicit none

      
      include 'mod06uw_data.inc'
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

c whuang 05/07/01: Added save  satement for common / mod06_data /
      save / mod06_data /
      
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
        call message( 'mod06_write_products', text //
     &   ' [OPERATOR ACTION: Check system resources]' ,
     &   0, 2 )
      endif
      if ( dtype(1:6) .ne. 'REAL*4' .and.
     &     dtype(1:9) .ne. 'INTEGER*2' .and.
     &     dtype(1:9) .ne. 'INTEGER*1' ) then
        write( text, '(''Cannot handle type '',a,'' for array '',a)' )
     &     dtype, arrnm       
        call message( 'mod06_write_products', text //
     &   ' [OPERATOR ACTION: Check system resources]' ,
     &   0, 2 )
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
        call message( 'mod06_write_products', text //
     &   ' [OPERATOR ACTION: Check system resources]' ,
     &   0, 2 )
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
        call message( 'mod06_write_products', text //
     &   ' [OPERATOR ACTION: Check system resources]' ,
     &   0, 2 )
      endif
                
      end
