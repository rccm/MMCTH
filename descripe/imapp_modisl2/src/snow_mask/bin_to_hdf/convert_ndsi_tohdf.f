C*----------------------------------------------------------------------
C    Copyright (C) 2011,  Space Science and Engineering Center,
C    University C  of Wisconsin-Madison, Madison WI.
C
C    This program is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    This program is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with this program.  If not, see <http://www.gnu.org/licenses/>.
C------------------------------------------------------------------------*/
C
C
              PROGRAM CONVERT_NDSI_TOHDF

c-----------------------------------------------------------------------
c Purpose:
c     To convert one of William Straka's MODIS Snow Mask binary files
c     into an HDF 4 file.
c
c     His file contains one 1 km parameter that matches the size
c     of the L1B file. The input L1B file is used to extract the
c     number of lines and elements and to copy the lat/lon values
c     from.
c
c     Usage: convert_ndsi_tohdf.exe L1B_1km_file Binary_file Output_file
c
c     Written by Kathleen Strabala 8 February 2008
c-----------------------------------------------------------------------

      implicit none

c ... Include files
      include 'hdf.inc'

c ... Local variables
      integer file_id, hdf_lun, outrec, line
      integer varid_lat, varid_lon, varid_1km, lat_id, lon_id, ndsi_id
      integer rtn, num_arg, irt, iostat
      integer rank, dimsizes(3), data_type, num_attrs, dims(2)
      character*100 sds_name 
      character*40 dim_name(2)
      character*100 arrayname
      integer num_pixels_5km, num_pixels_1km, num_lines_5km, 
     +        num_lines_1km, num_scans
      double precision scale_fac, offset
      
      integer lines_1km
      parameter (lines_1km = 10)
  
      integer start(2), stride(2), edge(2)

      integer max_5km_pixel, max_1km_pixel
      parameter (max_5km_pixel = 300, max_1km_pixel = 1354)

      integer inlun
      parameter (inlun = 31)

      real float_buffer(max_5km_pixel)
      integer*2  image(max_1km_pixel)
            
      integer dim

      character*200 in_file, in_bin, out_hdf

c ... External functions
      integer sfstart, sfselect, sfn2index, sfginfo, sfrdata,
     &  sfendacc, sfend, sfcreate, sfwdata, sfsfill, sfsdmname,
     &  sfdimid, sfsnatt, sfsrange, sfscatt
      external sfstart, sfselect, sfn2index, sfginfo, sfrdata,
     &  sfendacc, sfend, sfcreate, sfwdata, sfsfill, sfsdmname,
     &  sfdimid, sfsnatt, sfsrange, sfscatt

      integer iargc

      integer string_length
      external string_length

c-----------------------------------------------------------------------
c     GET INPUT ARGUMENTS
c-----------------------------------------------------------------------

c ... Check number of arguments
      num_arg = iargc()
      if (num_arg .ne. 3) then
        print *, 'Usage: convert_ndsi_tohdf.exe in_hdf in_bin out_hdf'
        print *, 'where'
        print *, 'in_hdf:  name of input L1B 1km HDF file'
        print *, 'in_bin:  name of input binary NDSI file'
        print *, 'out_hdf: name of output NDSI HDF file'
        call exit(-1)
      endif

c ... Extract arguments
      call getarg(1, in_file)
      call getarg(2, in_bin)
      call getarg(3, out_hdf)

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     OPEN FILES
c-----------------------------------------------------------------------

c ... Print the input and output filenames
      print '(''Input L1B 1km HDF file: '', a)', in_file(1: string_length(in_file))
      print '(''Input binary NDSI file: '', a)', in_bin(1 : string_length(in_bin))
      print '(''Output HDF NDSI  file: '', a)' ,out_hdf(1 : string_length(out_hdf))



c ... Open the input HDF file for read only
c ... (DFACC_READ is defined in hdf.inc)
      file_id = sfstart(in_file, DFACC_READ)
      if (file_id .eq. -1) then
        print *, 'Error: Could not open input MOD03 HDF file for reading'
        call exit(-1)
      endif

c ... Get the SDS id for each lat/lon array
      varid_lat  = sfselect(file_id, sfn2index(file_id, 'Latitude'))
      varid_lon  = sfselect(file_id, sfn2index(file_id, 'Longitude'))
      if (varid_lat  .eq. -1 .or.
     &    varid_lon  .eq. -1) then
        print *, 'Error: Unable to read LAT/LON values from MODIS 1KM HDF file' 
        call exit(-1) 
      endif

c ... Get the number of 5km pixels, lines, and scans
      rtn = sfginfo(varid_lat, sds_name, rank, dimsizes, data_type, num_attrs) 
      print '(''5KM DIMENSION SIZES '')'
      num_pixels_5km = dimsizes(1)
      num_lines_5km  = dimsizes(2)
      print '(''5 km Pixels = '', i6)', num_pixels_5km
      print '(''5 km Lines  = '', i6)', num_lines_5km

c ... Get 1km values
      varid_1km  = sfselect(file_id, sfn2index(file_id, 'EV_1KM_RefSB'))

c ... Get the number of 1km pixels, lines, and scans
      rtn = sfginfo(varid_1km, sds_name, rank, dimsizes, data_type, num_attrs) 
      num_pixels_1km = dimsizes(1)
      num_lines_1km  = dimsizes(2)
      num_scans = dimsizes(2) / 10
      print '(''1KM DIMENSION SIZES '')'
      print '(''1 km Pixels = '', i6)', num_pixels_1km
      print '(''1 km Lines  = '', i6)', num_lines_1km
      print '(''Number of Scans  = '', i6)', num_scans
      
c ... Open the input snow mask binary file
      open(inlun, file=in_bin, iostat=iostat, status='unknown',
     &  access='direct', form='unformatted', recl=(num_pixels_1km * 2))
      if (iostat .ne. 0) then
        print *, 'Error: Could not open input binary file for reading'
        call exit(-1)
      endif

c ... Open the output HDF file
      hdf_lun = sfstart(out_hdf, DFACC_CREATE)
      if (hdf_lun .eq. -1) then
         print *, 'Error: Could not open output HDF file'
         call exit(-1) 
      endif

c ... Copy the Latitude and Longitude fields
      dims(1) = num_pixels_5km
      dims(2) = num_lines_5km
      rank = 2
      lat_id = sfcreate(hdf_lun,'Latitude',DFNT_FLOAT32,rank,
     &             dims)
      if (lat_id .lt. 0) then
          print *, 'Error: Unable to copy LAT values from MODIS 1KM HDF file'
      endif
      dim_name(1) = 'Cell_Across_Swath_5km'
      dim_name(2) = 'Cell_Along_Swath_5km'
      irt = sfsdmname(sfdimid(lat_id,0),dim_name(1))
      irt = sfsdmname(sfdimid(lat_id,1),dim_name(2))
      irt = sfsfill(lat_id,-999.99)


c ... Copy the Longitude fields
      dims(1) = num_pixels_5km
      dims(2) = num_lines_5km
      rank = 2
      lon_id = sfcreate(hdf_lun,'Longitude',DFNT_FLOAT32,rank,
     &             dims)
      if (lon_id .lt. 0) then
          print *, 'Error: Unable to copy LON values from MODIS 1KM HDF file'
      endif
      dim_name(1) = 'Cell_Across_Swath_5km'
      dim_name(2) = 'Cell_Along_Swath_5km'
      irt = sfsdmname(sfdimid(lon_id,0),dim_name(1))
      irt = sfsdmname(sfdimid(lon_id,1),dim_name(2))
      irt = sfsfill(lon_id,-999.99)

      dims(1) = num_pixels_1km
      dims(2) = num_lines_1km
      rank = 2
      ndsi_id = sfcreate(hdf_lun,'Snow_Mask',DFNT_INT16,rank,dims)
      if (ndsi_id .lt. 0) then
         print *, 'Error: Unable to create snow mask HDF array'
      endif
      dim_name(1) = 'Cell_Across_Swath_1km'
      dim_name(2) = 'Cell_Along_Swath_1km'
      irt = sfsdmname(sfdimid(ndsi_id,0),dim_name(1))
      irt = sfsdmname(sfdimid(ndsi_id,1),dim_name(2))

c ... Attach local attributes to the SDS
      scale_fac = 1.0d0
      irt = sfsnatt(ndsi_id,'scale_factor',DFNT_FLOAT64,1,scale_fac)
      offset = 0.0d0
      irt = sfsnatt(ndsi_id,'add_offset',DFNT_FLOAT64,1,offset)
      irt = sfsrange(ndsi_id,1000,-1000)
      irt = sfsfill(ndsi_id,-32768)

      print '(''Copying Latitude and Longitude 5km arrays '')'
      print '(''----------------------------------------- '')'
c ... Loop over each line
      do line = 0, (num_lines_5km - 1), 1

c ...   Report progress 
        if (mod(line, 100) .eq. 0) print '(''Scan: '', i6)', line/2

c ...   Set start vector
        start(1) = 0
        start(2) = line

c ...   Set stride vector
        stride(1) = 1
        stride(2) = 1
      
c ...   Set edge vector
        edge(1) = num_pixels_5km
        edge(2) = 1

c ...   Read latitude data (float values, no scaling)
        rtn = sfrdata(varid_lat, start, stride, edge, float_buffer)
        if (rtn .eq. -1) then
          print *, 'Error reading latitude data'
          call exit(-1)
        endif

c ...   Write latitude data back to the file
        rtn = sfwdata(lat_id, start, stride, edge, float_buffer)
        if (rtn .eq. -1) then
          print *, 'Error writing latitude data'
          call exit(-1)
        endif

c ...   Read longitude data (float values, no scaling)
        rtn = sfrdata(varid_lon, start, stride, edge, float_buffer)
        if (rtn .eq. -1) then
          print *, 'Error reading latitude data'
          call exit(-1)
        endif

c ...   Write longitude data back to the file
        rtn = sfwdata(lon_id, start, stride, edge, float_buffer)
        if (rtn .eq. -1) then
          print *, 'Error writing latitude data'
          call exit(-1)
        endif
        
      enddo

c ... Read and write Snow Mask
c ... Initialize output record counter
      outrec = 1
      print '('' '')'
      print '(''Converting Snow Mask 1km array '')'
      print '(''------------------------------ '')'
      do line = 0, (num_lines_1km - 1), 1
 
c ...   Report progress on every 10th scan
        if (mod(line, 100) .eq. 0) print '(''Scan: '', i6)', line / 10

c ...   Read the snow mask data
        read(inlun, rec=outrec, iostat=iostat) image
        if (iostat .ne. 0) then
          write(*,*)'iostat', iostat
          print *, 'Error: Read of input snow mask file failed'
          call exit(-1)
        endif
    
c ...   Set start vector
        start(1) = 0
        start(2) = line

c ...   Set stride vector
        stride(1) = 1
        stride(2) = 1

c ...   Set edge vector
        edge(1) = num_pixels_1km
        edge(2) = 1

c ...   Write the snow mask data (float values, no scaling)
        rtn = sfwdata(ndsi_id, start, stride, edge, image)
        if (rtn .eq. -1) then
          print *, 'Error writing snow mask data'
          call exit(-1)
        endif

        outrec = outrec + 1

      enddo
      print '(''Snow Mask Conversion to HDF is complete '')'

c ... Close new HDF file
      irt = sfend(hdf_lun)

      end 


