      program convert_cloudmask
      IMPLICIT NONE	        

C--------------------------------------------------------------------
C      This program will read in the 64 bit MODIS cloud mask and create 
C      an HDF4 product that includes SDS arrays for each of the first
C      8 bits.  It is hoped this will make the use and visualization 
C      of the cloud mask easier.
C
C      It will also copy the 5 km geolocation information into the 
C      file too.
C      
C      Written by Kathleen Strabala November 2009
c      Space Science and Engineering Center
C      University of Wisconsin - Madison
C      1225 West Dayton Street
C      Madison,  WI   53706
C
C     
C      USE:  convert_cloudmask.exe  <input_cloudmask.hdf> <outfile>
C
C--------------------------------------------------------------------

c ... Include files
      include 'hdf.inc'

c ... Variables
      character*255 input_file, output_file

      integer varid_imask, varid_lat, varid_lon, file_id, outlun, isamp
      parameter(isamp = 5)
      integer rtn
      integer dimsizes(3), data_type, num_attrs
      character*100 sdsname
      character*200 long_name
      integer num_pixels, num_lines, num_scans, mask_bytes, num_arg,
     +        iostat, imask, i, irt,  rank, sds_id
      integer*2 fill_value
      real fill_value_r

      character*40 dim_name(2)
      integer dims(2), valid_range(2)
      real valid_range_r(2)
      double precision offset, scale_fac

      integer lines_1km
      parameter (lines_1km = 10)

      integer start1(3), stride(3), edge1(3)
      integer start2(2), stride2(2), edge2(2)
      integer start3(2), stride3(2), edge3(2)

      integer max_pixel, mask_band, max_line
      parameter (max_pixel = 1354, mask_band = 6, max_line=8000)

      byte mask_buffer(max_pixel*max_line)
      integer mask_integer
      integer*2 retrieved(max_pixel*max_line), category(max_pixel*max_line), 
     +        day_night(max_pixel*max_line), sunglint(max_pixel*max_line), 
     +         snow_ice(max_pixel*max_line), land_cat(max_pixel*max_line)
      real lat_buffer((max_pixel/isamp)*(max_line/isamp))
      real lon_buffer((max_pixel/isamp)*(max_line/isamp))

c ... intrinsic functions
      intrinsic iand, ishft

c ... External functions
      integer sfstart, sfselect, sfn2index, sfginfo, sfrdata,
     &  sfendacc, sfend, sfcreate, sfsdmname, sfsfill, sfsrange,
     &  sfdimid, sfsnatt, sfwdata, sfscatt, copy_info
      external sfstart, sfselect, sfn2index, sfginfo, sfrdata,
     &  sfendacc, sfend, sfcreate, sfsdmname, sfsfill, sfsrange,
     &  sfdimid, sfsnatt, sfwdata, sfscatt, copy_info

      integer iargc

      integer string_length
      external string_length

c ... Read input arguments
      
      write( *, '(2x,''Initiated MODIS Cloud Mask Converter  '')')

c ... Get input arguments

c ... Check number of arguments
      num_arg = iargc()
      if (num_arg .ne. 2) then
          print *, 'Usage: convert_cloudmask.exe input_file output_file'
          print *, 'where'
          print *, 'input_file is the input MODIS cloudmask HDF file name'    
          print *, 'output_file is the output MODIS cloudmask HDF file name'
          call exit(-1)
      endif

c ... Extract arguments
      call getarg(1,input_file)
      call getarg(2,output_file)

c-----------------------------------------------------------------------
c     OPEN FILES
c-----------------------------------------------------------------------

c ... Print the input and output filenames
      print '(''Input Cloud Mask file: '', a)', input_file(1: string_length(input_file))
      print '(''Output file: '', a)', output_file(1 : string_length(output_file))

c ... Open Files
c ... (DFACC_READ is defined in hdf.inc)
      file_id = sfstart(input_file, DFACC_READ)
      if (file_id .eq. -1) then
        print *, 'Error: Could not open input MOD021KM HDF file for reading'
        call exit(-1)
      endif

c ... Get the SDS id for the cloud mask array
      varid_imask = sfselect(file_id, sfn2index(file_id, 'Cloud_Mask'))
      if (varid_imask .eq. -1) then
        print *, 'Error: Input file is not a Cloud Mask HDF file'
        call exit(-1)
      endif

c ... Get the number of pixels, lines, and scans
      rtn = sfginfo(varid_imask, sdsname, rank, dimsizes, data_type, num_attrs)
      num_pixels = dimsizes(1)
      num_lines  = dimsizes(2)
      mask_bytes  = dimsizes(3)
      num_scans  = num_lines / lines_1km
      print '(''Pixels = '', i6)', num_pixels
      print '(''Lines  = '', i6)', num_lines
      print '(''Mask Bytes  = '', i6)', mask_bytes
      print '(''Scans  = '', i6)', num_scans


c ... Set start vector for both the mask and maskqa array
      start1(1) = 0
      start1(2) = 0
      start1(3) = 0

c ... Set stride vector
      stride(1) = 1
      stride(2) = 1
      stride(3) = 1

c ... Set edge vector
      edge1(1) = num_pixels
      edge1(2) = num_lines
      edge1(3) = 1

c ... Read in the first byte of the cloud mask array for the whole granule
      rtn = sfrdata(varid_imask, start1, stride, edge1, mask_buffer)
      if (rtn .eq. -1) then
        print *, 'Error reading cloudmask data'
        call exit(-1)
      endif

c ... Now reconfigure the array and place in individual SDSs
      do i = 1 , num_pixels*num_lines
            retrieved(i)=9999
            category(i)=9999
            day_night(i)=9999
            sunglint(i)=9999
            snow_ice(i)=9999
            land_cat(i)=9999
      enddo

      
      do i = 1 , num_pixels*num_lines
            imask=1
            mask_integer=0
            mask_integer=mask_buffer(i)
            
c ...       Was mask retrieved?
            if (iand(imask,mask_integer) .eq. 1) then
                retrieved(i)=1
                category(i)=0
                day_night(i)=1
                sunglint(i)=2
                snow_ice(i)=2
                land_cat(i)=0

c ...       Four category cloud mask
                if ((iand(4,mask_integer) .eq. 4) .and. 
     +             (iand(2,mask_integer) .eq. 2)) category(i)=4
                if ((iand(4,mask_integer) .eq. 4) .and. 
     +             (iand(2,mask_integer) .eq. 0)) category(i)=3
                if ((iand(4,mask_integer) .eq. 0) .and. 
     +             (iand(2,mask_integer) .eq. 2)) category(i)=2
                if ((iand(4,mask_integer) .eq. 0) .and. 
     +             (iand(2,mask_integer) .eq. 0)) category(i)=1
c ...       Day/Night
                if (iand(ishft(mask_integer,-3),1) .eq. 1) day_night(i)=2
c ...       Sunglint
                if (iand(16,mask_integer) .ne. 16) sunglint(i)=1
c ...       Snow and ice
                if (iand(32,mask_integer) .ne. 32) snow_ice(i)=1
c ...       Four category land/sea flag
                if ((iand(128,mask_integer) .eq. 128) .and.
     +             (iand(64,mask_integer) .eq. 64)) land_cat(i)=4
                if ((iand(128,mask_integer) .eq. 128) .and. 
     +             (iand(64,mask_integer) .eq. 0)) land_cat(i)=3
                if ((iand(128,mask_integer) .eq. 0) .and. 
     +             (iand(64,mask_integer) .eq. 64)) land_cat(i)=2
                if ((iand(128,mask_integer) .eq. 0) .and. 
     +             (iand(64,mask_integer) .eq. 0)) land_cat(i)=1
            else
                continue
            endif
      enddo

c ... Create and open the output hdf file
      outlun=sfstart(output_file,DFACC_CREATE)
      if (outlun .lt. 0) then
        write(*,'( 2x,''Unable to open output file '',A)')
     +   output_file
        call exit(1)
      endif

c ------------------------------------------------------------------------
c     WRITE OUTPUT ARRAYS 
c ------------------------------------------------------------------------

c ... 1. CLOUD MASK
c ------------------------------------------------------------------------
c ... Write out arrays
c ... Create the Science Data Sets
c ... Create the output sds
      dims(1) = num_pixels
      dims(2) = num_lines
      rank = 2
      sdsname = 'MODIS_Cloud_Mask'
      sds_id = sfcreate(outlun,sdsname,DFNT_UINT16,rank,dims)
      if (sds_id .lt. 0) then
         write(*,'( 2x,''Unble to create MODIS cloud mask array '',
     +              I10)') sds_id
      endif

c ... name the dimensions
      dim_name(1) = 'Cell_Across_Swath_1km'
      dim_name(2) = 'Cell_Along_Swath_1km'
      irt = sfsdmname(sfdimid(sds_id,0),dim_name(1))
      irt = sfsdmname(sfdimid(sds_id,1),dim_name(2))

c ... Attach local attributes to the SDS
      long_name = 'MODIS Cloud Mask Categories: '
     +      //'9999 - No Mask Determined,  1 - Cloudy, 2 - Undetermined '
     +      //' 3 - Confident Clear, 4 - High Confident Clear'
      irt = sfscatt(sds_id,'long_name',DFNT_CHAR8,string_length(long_name),long_name)
      if( irt.lt.0 )
     +  write(*,'( 2x,''Unable to attach long name to MODIS_Cloud_Mask'')')

c ..  Initialize the SDS with fill values
      fill_value = 9999
      irt = sfsfill(sds_id,fill_value)
      if( irt.lt.0 )
     +  write(*,'( 2x,''Unable to put fill value into array '')')

c ..  Set the range flags
      valid_range(1) = 0
      valid_range(2) = 4
      irt = sfsrange(sds_id,valid_range(2),valid_range(1))
      if( irt.lt.0 )
     +  write(*,'( 2x,''Unable to put valid range value into array '')')

c ..  Write scale factor local attribute
      scale_fac = 1.0d0
      irt = sfsnatt(sds_id,'scale_factor',DFNT_FLOAT64,1,scale_fac)
      if( irt.lt.0 )
     +  write(*,'( 2x,''Unable to put scale factor into array '')')

c ..  Write offset local attribute
      offset = 0.0d0
      irt = sfsnatt(sds_id,'add_offset',DFNT_FLOAT64,1,offset)
      if( irt.lt.0 )
     +  write(*,'( 2x,''Unable to put offset into array '')')

c ... Write output data to file
      start2(1) = 0
      start2(2) = 0
      stride2(1) = 1
      stride2(2) = 1
      edge2(1) = num_pixels
      edge2(2) = num_lines

      irt = sfwdata(sds_id,start2,stride2,edge2,category)
      if (irt .lt. 0) then
          write(*,'( 2x,''Unable to write the Cloud Mask Output array''
     +         ,i5)') irt
      endif

c ... terminate access to the SDS
      irt = sfendacc(sds_id)


c ... 2. DAY / NIGHT FLAG
c -----------------------------------------------------------------------
c ... Create the output sds
      dims(1) = num_pixels
      dims(2) = num_lines
      rank = 2
      sdsname = 'MODIS_Day_Night_Flag'
      sds_id = sfcreate(outlun,sdsname,DFNT_UINT16,rank,dims)
      if (sds_id .lt. 0) then
         write(*,'( 2x,''Unble to create MODIS day/night array '',
     +              I10)') sds_id
      endif

c ... name the dimensions
      dim_name(1) = 'Cell_Across_Swath_1km'
      dim_name(2) = 'Cell_Along_Swath_1km'
      irt = sfsdmname(sfdimid(sds_id,0),dim_name(1))
      irt = sfsdmname(sfdimid(sds_id,1),dim_name(2))

c ... Attach local attributes to the SDS
      long_name = 'MODIS Day or Night Flag: '
     +      //'9999 - No Mask Determined,  1 - Night, 2 - Day '
      irt = sfscatt(sds_id,'long_name',DFNT_CHAR8,string_length(long_name),long_name)
      if( irt.lt.0 )
     +  write(*,'( 2x,''Unable to attach long name to MODIS_Day_Night_Flag'')')

c ..  Initialize the SDS with fill values
      fill_value = 9999
      irt = sfsfill(sds_id,fill_value)
      if( irt.lt.0 )
     +  write(*,'( 2x,''Unable to put fill value into array '')')

c ..  Set the range flags
      valid_range(1) = 0
      valid_range(2) = 2
      irt = sfsrange(sds_id,valid_range(2),valid_range(1))
      if( irt.lt.0 )
     +  write(*,'( 2x,''Unable to put valid range value into array '')')

c ..  Write scale factor local attribute
      scale_fac = 1.0d0
      irt = sfsnatt(sds_id,'scale_factor',DFNT_FLOAT64,1,scale_fac)
      if( irt.lt.0 )
     +  write(*,'( 2x,''Unable to put scale factor into array '')')

c ..  Write offset local attribute
      offset = 0.0d0
      irt = sfsnatt(sds_id,'add_offset',DFNT_FLOAT64,1,offset)
      if( irt.lt.0 )
     +  write(*,'( 2x,''Unable to put offset into array '')')

c ... Write output data to file
      start2(1) = 0
      start2(2) = 0
      stride2(1) = 1
      stride2(2) = 1
      edge2(1) = num_pixels
      edge2(2) = num_lines

      irt = sfwdata(sds_id,start2,stride2,edge2,day_night)
      if (irt .lt. 0) then
          write(*,'( 2x,''Unable to write the Day/Night Output array''
     +         ,i5)') irt
      endif

c ... terminate access to the SDS
      irt = sfendacc(sds_id)



c ... 3. SUNGLINT  FLAG
c -----------------------------------------------------------------------
c ... Create the output sds
      dims(1) = num_pixels
      dims(2) = num_lines
      rank = 2
      sdsname = 'MODIS_Sunglint_Flag'
      sds_id = sfcreate(outlun,sdsname,DFNT_UINT16,rank,dims)
      if (sds_id .lt. 0) then
         write(*,'( 2x,''Unble to create MODIS sunglint array '',
     +              I10)') sds_id
      endif

c ... name the dimensions
      dim_name(1) = 'Cell_Across_Swath_1km'
      dim_name(2) = 'Cell_Along_Swath_1km'
      irt = sfsdmname(sfdimid(sds_id,0),dim_name(1))
      irt = sfsdmname(sfdimid(sds_id,1),dim_name(2))

c ... Attach local attributes to the SDS
      long_name = 'MODIS Sunglint Flag: '
     +      //'9999 - No Mask Determined,  1 - Sunglint in Pixel, 2 - No Sunglint in Pixel '
      irt = sfscatt(sds_id,'long_name',DFNT_CHAR8,string_length(long_name),long_name)
      if( irt.lt.0 )
     +  write(*,'( 2x,''Unable to attach long name to MODIS_Sunglint_Flag'')')

c ..  Initialize the SDS with fill values
      fill_value = 9999
      irt = sfsfill(sds_id,fill_value)
      if( irt.lt.0 )
     +  write(*,'( 2x,''Unable to put fill value into array '')')

c ..  Set the range flags
      valid_range(1) = 0
      valid_range(2) = 2
      irt = sfsrange(sds_id,valid_range(2),valid_range(1))
      if( irt.lt.0 )
     +  write(*,'( 2x,''Unable to put valid range value into array '')')

c ..  Write scale factor local attribute
      scale_fac = 1.0d0
      irt = sfsnatt(sds_id,'scale_factor',DFNT_FLOAT64,1,scale_fac)
      if( irt.lt.0 )
     +  write(*,'( 2x,''Unable to put scale factor into array '')')

c ..  Write offset local attribute
      offset = 0.0d0
      irt = sfsnatt(sds_id,'add_offset',DFNT_FLOAT64,1,offset)
      if( irt.lt.0 )
     +  write(*,'( 2x,''Unable to put offset into array '')')

c ... Write output data to file
      start2(1) = 0
      start2(2) = 0
      stride2(1) = 1
      stride2(2) = 1
      edge2(1) = num_pixels
      edge2(2) = num_lines

      irt = sfwdata(sds_id,start2,stride2,edge2,sunglint)
      if (irt .lt. 0) then
          write(*,'( 2x,''Unable to write the Sunglint Output array''
     +         ,i5)') irt
      endif

c ... terminate access to the SDS
      irt = sfendacc(sds_id)

c ... 4. SNOW / ICE FLAG
c -----------------------------------------------------------------------
c ... Create the output sds
      dims(1) = num_pixels
      dims(2) = num_lines
      rank = 2
      sdsname = 'MODIS_Snow_Ice_Flag'
      sds_id = sfcreate(outlun,sdsname,DFNT_UINT16,rank,dims)
      if (sds_id .lt. 0) then
         write(*,'( 2x,''Unble to create MODIS snow and ice array '',
     +              I10)') sds_id
      endif

c ... name the dimensions
      dim_name(1) = 'Cell_Across_Swath_1km'
      dim_name(2) = 'Cell_Along_Swath_1km'
      irt = sfsdmname(sfdimid(sds_id,0),dim_name(1))
      irt = sfsdmname(sfdimid(sds_id,1),dim_name(2))

c ... Attach local attributes to the SDS
      long_name = 'MODIS Snow/Ice Flag: '
     +      //'9999 - No Mask Determined,  1 - Snow/Ice in Pixel, 2 - No Snow/Ice in Pixel '
      irt = sfscatt(sds_id,'long_name',DFNT_CHAR8,string_length(long_name),long_name)
      if( irt.lt.0 )
     +  write(*,'( 2x,''Unable to attach long name to MODIS_Snow_Ice_Flag'')')

c ..  Initialize the SDS with fill values
      fill_value = 9999
      irt = sfsfill(sds_id,fill_value)
      if( irt.lt.0 )
     +  write(*,'( 2x,''Unable to put fill value into array '')')

c ..  Set the range flags
      valid_range(1) = 0
      valid_range(2) = 2
      irt = sfsrange(sds_id,valid_range(2),valid_range(1))
      if( irt.lt.0 )
     +  write(*,'( 2x,''Unable to put valid range value into array '')')

c ..  Write scale factor local attribute
      scale_fac = 1.0d0
      irt = sfsnatt(sds_id,'scale_factor',DFNT_FLOAT64,1,scale_fac)
      if( irt.lt.0 )
     +  write(*,'( 2x,''Unable to put scale factor into array '')')

c ..  Write offset local attribute
      offset = 0.0d0
      irt = sfsnatt(sds_id,'add_offset',DFNT_FLOAT64,1,offset)
      if( irt.lt.0 )
     +  write(*,'( 2x,''Unable to put offset into array '')')

c ... Write output data to file
      start2(1) = 0
      start2(2) = 0
      stride2(1) = 1
      stride2(2) = 1
      edge2(1) = num_pixels
      edge2(2) = num_lines

      irt = sfwdata(sds_id,start2,stride2,edge2,snow_ice)
      if (irt .lt. 0) then
          write(*,'( 2x,''Unable to write the Snow and Ice Output array''
     +         ,i5)') irt
      endif

c ... terminate access to the SDS
      irt = sfendacc(sds_id)



c ... 1. MODIS CLOUDMASK LAND/SEA MASK
c ------------------------------------------------------------------------
c ... Write out arrays
c ... Create the Science Data Sets
c ... Create the output sds
      dims(1) = num_pixels
      dims(2) = num_lines
      rank = 2
      sdsname = 'MODIS_Simple_LandSea_Mask'
      sds_id = sfcreate(outlun,sdsname,DFNT_UINT16,rank,dims)
      if (sds_id .lt. 0) then
         write(*,'( 2x,''Unble to create MODIS land/sea mask array '',
     +              I10)') sds_id
      endif

c ... name the dimensions
      dim_name(1) = 'Cell_Across_Swath_1km'
      dim_name(2) = 'Cell_Along_Swath_1km'
      irt = sfsdmname(sfdimid(sds_id,0),dim_name(1))
      irt = sfsdmname(sfdimid(sds_id,1),dim_name(2))

c ... Attach local attributes to the SDS
      long_name = 'MODIS Simple Land and Sea Mask Categories: '
     +      //'9999 - No Mask Determined,  1 - Water, 2 - Coast '
     +      //' 3 - Desert, 4 - Land'
      irt = sfscatt(sds_id,'long_name',DFNT_CHAR8,string_length(long_name),long_name)
      if( irt.lt.0 )
     +  write(*,'( 2x,''Unable to attach long name to MODIS_Land_Sea_Mask'')')

c ..  Initialize the SDS with fill values
      fill_value = 9999
      irt = sfsfill(sds_id,fill_value)
      if( irt.lt.0 )
     +  write(*,'( 2x,''Unable to put fill value into array '')')

c ..  Set the range flags
      valid_range(1) = 0
      valid_range(2) = 4
      irt = sfsrange(sds_id,valid_range(2),valid_range(1))
      if( irt.lt.0 )
     +  write(*,'( 2x,''Unable to put valid range value into array '')')

c ..  Write scale factor local attribute
      scale_fac = 1.0d0
      irt = sfsnatt(sds_id,'scale_factor',DFNT_FLOAT64,1,scale_fac)
      if( irt.lt.0 )
     +  write(*,'( 2x,''Unable to put scale factor into array '')')

c ..  Write offset local attribute
      offset = 0.0d0
      irt = sfsnatt(sds_id,'add_offset',DFNT_FLOAT64,1,offset)
      if( irt.lt.0 )
     +  write(*,'( 2x,''Unable to put offset into array '')')

c ... Write output data to file
      start2(1) = 0
      start2(2) = 0
      stride2(1) = 1
      stride2(2) = 1
      edge2(1) = num_pixels
      edge2(2) = num_lines

      irt = sfwdata(sds_id,start2,stride2,edge2,land_cat)
      if (irt .lt. 0) then
          write(*,'( 2x,''Unable to write the MODIS land/sea mask array''
     +         ,i5)') irt
      endif

c ... terminate access to the SDS
      irt = sfendacc(sds_id)

      write( *, '(2x)' )
      write( *, '(2x,''Finished writing output cloud mask parameters. '')' )
      write( *, '(2x)' )
c --------------------------------------------------------------------
c     FINISHED WRITING CLOUD MASK ARRAYS TO OUTPUT FILE
c --------------------------------------------------------------------
c --------------------------------------------------------------------
c     COPY CLOUD MASK LATITUDES AND LONGITUDES TO OUTPUT ARRAY
c --------------------------------------------------------------------

c ... Now copy the lat/lon information into the output file
c ... Get the SDS id for the cloud mask array
      varid_lat = sfselect(file_id, sfn2index(file_id, 'Latitude'))
      if (varid_lat .eq. -1) then
        print *, 'Error: Could not find latitude array in Mask HDF file'
        call exit(-1)
      endif
      
      varid_lon = sfselect(file_id, sfn2index(file_id, 'Longitude'))
      if (varid_lon .eq. -1) then
        print *, 'Error: Could not find longitude in Mask HDF file'
        call exit(-1)
      endif

c ... Get the number of pixels, lines, and scans for lat/lon arrays
      rtn = sfginfo(varid_lat, sdsname, rank, dimsizes, data_type, num_attrs)
      num_pixels = dimsizes(1)
      num_lines  = dimsizes(2)
      print '(''GEO Pixels = '', i6)', num_pixels
      print '(''GEO lines  = '', i6)', num_lines

c ... Set start vector for both the lat and lon arrays
      start3(1) = 0
      start3(2) = 0

c ... Set stride vector
      stride3(1) = 1
      stride3(2) = 1

c ... Set edge vector
      edge3(1) = num_pixels
      edge3(2) = num_lines

c ... Read in the latitude from the geolocation file
      rtn = sfrdata(varid_lat, start3, stride3, edge3, lat_buffer)
      if (rtn .eq. -1) then
        print *, 'Error reading latitude data'
        call exit(-1)
      endif

c ... Read in the longitude from the geolocation file
      rtn = sfrdata(varid_lon, start3, stride3, edge3, lon_buffer)
      if (rtn .eq. -1) then
        print *, 'Error reading longitude data'
        call exit(-1)
      endif

c ... Now write it into our new output file
c ... Create the output sds
      dims(1) = num_pixels
      dims(2) = num_lines
      rank = 2
      sdsname = 'Latitude'
      sds_id = sfcreate(outlun,sdsname,DFNT_FLOAT32,rank,dims)
      if (sds_id .lt. 0) then
         write(*,'( 2x,''Unble to create MODIS latitude array '',
     +              I10)') sds_id
      endif

c ... name the dimensions
      dim_name(1) = 'Cell_Across_Swath_5km'
      dim_name(2) = 'Cell_Along_Swath_5km'
      irt = sfsdmname(sfdimid(sds_id,0),dim_name(1))
      irt = sfsdmname(sfdimid(sds_id,1),dim_name(2))

c ... Attach local attributes to the SDS
      long_name = 'Geodetic Latitude'
      irt = sfscatt(sds_id,'long_name',DFNT_CHAR8,string_length(long_name),long_name)
      if( irt.lt.0 )
     +  write(*,'( 2x,''Unable to attach long name to Latitude'')')

c ..  Initialize the SDS with fill values
      fill_value_r = -999.99
      irt = sfsfill(sds_id,fill_value_r)
      if( irt.lt.0 )
     +  write(*,'( 2x,''Unable to put fill value into array '')')

c ..  Set the range flags
      valid_range_r(1) = -90.0
      valid_range_r(2) = 90.0
      irt = sfsrange(sds_id,valid_range_r(2),valid_range_r(1))
      if( irt.lt.0 )
     +  write(*,'( 2x,''Unable to put valid range value into array '')')

c ..  Write scale factor local attribute
      scale_fac = 1.0d0
      irt = sfsnatt(sds_id,'scale_factor',DFNT_FLOAT64,1,scale_fac)
      if( irt.lt.0 )
     +  write(*,'( 2x,''Unable to put scale factor into array '')')

c ..  Write offset local attribute
      offset = 0.0d0
      irt = sfsnatt(sds_id,'add_offset',DFNT_FLOAT64,1,offset)
      if( irt.lt.0 )
     +  write(*,'( 2x,''Unable to put offset into array '')')

c ... Write output data to file
      irt = sfwdata(sds_id,start3,stride3,edge3,lat_buffer)
      if (irt .lt. 0) then
          write(*,'( 2x,''Unable to write the MODIS latitude array''
     +         ,i5)') irt
      endif
      
c ... terminate access to the latitude SDS
      irt = sfendacc(sds_id)

c ... Now write it into our new output file
c ... Create the output sds
      dims(1) = num_pixels
      dims(2) = num_lines
      rank = 2
      sdsname = 'Longitude'
      sds_id = sfcreate(outlun,sdsname,DFNT_FLOAT32,rank,dims)
      if (sds_id .lt. 0) then
         write(*,'( 2x,''Unble to create MODIS longitude array '',
     +              I10)') sds_id
      endif

c ... name the dimensions
      dim_name(1) = 'Cell_Across_Swath_5km'
      dim_name(2) = 'Cell_Along_Swath_5km'
      irt = sfsdmname(sfdimid(sds_id,0),dim_name(1))
      irt = sfsdmname(sfdimid(sds_id,1),dim_name(2))

c ... Attach local attributes to the SDS
      long_name = 'Geodetic Longitude'
      irt = sfscatt(sds_id,'long_name',DFNT_CHAR8,string_length(long_name),long_name)
      if( irt.lt.0 )
     +  write(*,'( 2x,''Unable to attach long name to Longitude'')')

c ..  Initialize the SDS with fill values
      fill_value_r = -999.99
      irt = sfsfill(sds_id,fill_value_r)
      if( irt.lt.0 )
     +  write(*,'( 2x,''Unable to put fill value into array '')')

c ..  Set the range flags
      valid_range_r(1) = -180.0
      valid_range_r(2) = 180.0
      irt = sfsrange(sds_id,valid_range_r(2),valid_range_r(1))
      if( irt.lt.0 )
     +  write(*,'( 2x,''Unable to put valid range value into array '')')

c ..  Write scale factor local attribute
      scale_fac = 1.0d0
      irt = sfsnatt(sds_id,'scale_factor',DFNT_FLOAT64,1,scale_fac)
      if( irt.lt.0 )
     +  write(*,'( 2x,''Unable to put scale factor into array '')')

c ..  Write offset local attribute
      offset = 0.0d0
      irt = sfsnatt(sds_id,'add_offset',DFNT_FLOAT64,1,offset)
      if( irt.lt.0 )
     +  write(*,'( 2x,''Unable to put offset into array '')')

c ... Write output data to file
      irt = sfwdata(sds_id,start3,stride3,edge3,lon_buffer)
      if (irt .lt. 0) then
          write(*,'( 2x,''Unable to write the MODIS longitude array''
     +         ,i5)') irt
      endif
      
c ... terminate access to the latitude SDS
      irt = sfendacc(sds_id)


c ... Close output file
      irt = sfend(outlun)

      write( *, '(2x)' )
      write( *, '(2x,''Finished copying lat/lon arrays to output file. '')' )

c --------------------------------------------------------------------
c     FINISHED
c --------------------------------------------------------------------

c ... 
      end


