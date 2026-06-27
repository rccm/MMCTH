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
              PROGRAM CONVERT_INVERSIONS_TOHDF

c-----------------------------------------------------------------------
c Purpose:
c     To convert one set of Yinghui Liu, William Straka and Jeff Key's 
c     MODIS Inversion Depth and Inversion Strength Binary files into
c     and HDF4 file.
c
c     The IMAPP Inversion Depth binary file contains a 1 km array
c     field in short integers scaled by 10 in units of meters.  The 
c     IMAPP Inversion Strength binary file contains a 1km array
c     field in short integers scaled by 10 in units of degrees C.  This
c     software will open both binary files, read out the values,
c     and write them into an HDF4 file.  It will also copy over the
c     5 km latitudes and longitudes out of the MODIS 1km L1B file
c     into the new HDF4 output file.
c
c     Usage: 
c      convert_inversions_tohdf.exe L1B_1km_file Inversion_Depth_binary_file 
c                Inversion_Strength_binary_file Output_HDF4_file
c
c     Written by Kathleen Strabala April 2011
c-----------------------------------------------------------------------

      implicit none

c ... Include files
      include 'hdf.inc'

c ... Local variables
      integer file_id, hdf_lun, outrec, line
      integer varid_lat, varid_lon, varid_1km, lat_id, lon_id
      integer invd_id, invs_id, rtn, num_arg, irt, iostat
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

      integer inlun1, inlun2
      parameter (inlun1 = 31,inlun2 = 32)

      real float_buffer(max_5km_pixel)
      integer*2  inv_dep(max_1km_pixel), inv_str(max_1km_pixel)
            
      integer dim

      character*132 in_file, in_bin1, in_bin2, out_hdf

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
      if (num_arg .ne. 4) then
        print *, 'Usage: convert_inversions_tohdf.exe in_hdf in_bin1 inbin2 out_hdf'
        print *, 'where'
        print *, 'in_hdf:  name of input L1B 1km HDF file'
        print *, 'in_bin1:  name of input binary Inversion Depth file'
        print *, 'in_bin2:  name of input binary Inversion Strength file'
        print *, 'out_hdf: name of output inversions HDF file'
        call exit(-1)
      endif

c ... Extract arguments
      call getarg(1, in_file)
      call getarg(2, in_bin1)
      call getarg(3, in_bin2)
      call getarg(4, out_hdf)

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     OPEN FILES
c-----------------------------------------------------------------------

c ... Print the input and output filenames
      print '(''Input L1B 1km HDF file: '', a)', in_file(1: string_length(in_file))
      print '(''Input binary Inversion Depth file: '', a)', in_bin1(1 : string_length(in_bin1))
      print '(''Input binary Inversion Strength file: '', a)', in_bin2(1 : string_length(in_bin2))
      print '(''Output HDF NDSI  file: '', a)' ,out_hdf(1 : string_length(out_hdf))



c ... Open the input HDF file for read only
c ... (DFACC_READ is defined in hdf.inc)
      file_id = sfstart(in_file, DFACC_READ)
      if (file_id .eq. -1) then
        print *, 'Error: Could not open input MOD1KM HDF file for reading'
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
      
c ... Open the input inversion depth binary file
      open(inlun1, file=in_bin1, iostat=iostat, status='unknown',
     &  access='direct', form='unformatted', recl=(num_pixels_1km * 2))
      if (iostat .ne. 0) then
        print *, 'Error: Could not open input binary inversion depth file for reading'
        call exit(-1)
      endif

c ... Open the input inversion strength binary file
      open(inlun2, file=in_bin2, iostat=iostat, status='unknown',
     &  access='direct', form='unformatted', recl=(num_pixels_1km * 2))
      if (iostat .ne. 0) then
        print *, 'Error: Could not open input binary inversion strength file for reading'
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

c ... Create the Inversion Depth SDS and attach attributes
      dims(1) = num_pixels_1km
      dims(2) = num_lines_1km
      rank = 2
      invd_id = sfcreate(hdf_lun,'Inversion_Depth',DFNT_INT16,rank,dims)
      if (invd_id .lt. 0) then
         print *, 'Error: Unable to create the Inversion Depth  HDF array'
      endif
      dim_name(1) = 'Cell_Across_Swath_1km'
      dim_name(2) = 'Cell_Along_Swath_1km'
      irt = sfsdmname(sfdimid(invd_id,0),dim_name(1))
      irt = sfsdmname(sfdimid(invd_id,1),dim_name(2))

c ... Attach local attributes to the SDS
      irt = sfsnatt(invd_id,"units",DFNT_CHAR8,1,"m")
      scale_fac = 0.1d0
      irt = sfsnatt(invd_id,'scale_factor',DFNT_FLOAT64,1,scale_fac)
      offset = 0.0d0
      irt = sfsnatt(invd_id,'add_offset',DFNT_FLOAT64,1,offset)
      irt = sfsrange(invd_id,20000,0)
      irt = sfsfill(invd_id,-10)

c ... Create the Inversion Strength SDS and attach attributes
      dims(1) = num_pixels_1km
      dims(2) = num_lines_1km
      rank = 2
      invs_id = sfcreate(hdf_lun,'Inversion_Strength',DFNT_INT16,rank,dims)
      if (invs_id .lt. 0) then
         print *, 'Error: Unable to create the Inversion Strength  HDF array'
      endif
      dim_name(1) = 'Cell_Across_Swath_1km'
      dim_name(2) = 'Cell_Along_Swath_1km'
      irt = sfsdmname(sfdimid(invs_id,0),dim_name(1))
      irt = sfsdmname(sfdimid(invs_id,1),dim_name(2))

c ... Attach local attributes to the SDS
      irt = sfsnatt(invs_id,"units",DFNT_CHAR8,1,"C")
      scale_fac = 0.1d0
      irt = sfsnatt(invs_id,'scale_factor',DFNT_FLOAT64,1,scale_fac)
      offset = 0.0d0
      irt = sfsnatt(invs_id,'add_offset',DFNT_FLOAT64,1,offset)
      irt = sfsrange(invs_id,1000,0)
      irt = sfsfill(invs_id,-10)


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

c ... Read and write Inversion fields
c ... Initialize output record counter
      outrec = 1
      print '('' '')'
      print '(''Converting Inversions 1km array '')'
      print '(''------------------------------------------------- '')'
      do line = 0, (num_lines_1km - 1), 1
 
c ...   Report progress on every 10th scan
        if (mod(line, 100) .eq. 0) print '(''Scan: '', i6)', line / 10

c ...   Read the Binary Inversions Depth data from the binary file
        read(inlun1, rec=outrec, iostat=iostat) inv_dep
        if (iostat .ne. 0) then
          write(*,*)'iostat', iostat
          print *, 'Error: Read of input Inversion Depth binary file failed'
          call exit(-1)
        endif

c ...   Read the Binary Inversion Strength data from the binary file
        read(inlun2, rec=outrec, iostat=iostat)  inv_str
        if (iostat .ne. 0) then
          write(*,*)'iostat', iostat
          print *, 'Error: Read of input Inversion Strength binary file failed'
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

c ...   Write the Inversion Depth data to the output HDF4 file (scaled integers)
        rtn = sfwdata(invd_id, start, stride, edge, inv_dep)
        if (rtn .eq. -1) then
          print *, 'Error writing Inversion Depth data to output HDF file'
          call exit(-1)
        endif

c ...   Write the Inversion Strength data to the output HDF4 file (scaled integers)
        rtn = sfwdata(invs_id, start, stride, edge, inv_str)
        if (rtn .eq. -1) then
          print *, 'Error writing Inversion Strength data to output HDF file'
          call exit(-1)
        endif

        outrec = outrec + 1

      enddo
      print '(''Inversion Files binary to HDF conversion is complete '')'

c ... Close new HDF file
      irt = sfend(hdf_lun)

      end 


