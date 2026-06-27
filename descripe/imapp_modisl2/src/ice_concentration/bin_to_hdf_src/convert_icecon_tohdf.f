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
              PROGRAM CONVERT_ICECON_TOHDF

c-----------------------------------------------------------------------
c Purpose:
c     To read in one of Yinghui Liu's MODIS Ice Concentration 
c     binary files and write out an HDF 4 file.
c
c     The IMAPP Ice Concentration binary file contains 4 fields
c     of 4 byte data:
c         Latitude, Longitude - 4 byte float
c         Ice Cover - 4 byte integer
c         Ice Concentration - 4 byte float
c     The data is stored as band sequential, meaning the Latitude
c     array for the entire pass will be first, followed by the Longitude,
c     the Ice Mask, etc.  
c 
c     This software will extract that information along with the 
c     latitude and longitudes fields at 5km from the L1B file
c     write them to an HDF4 file with the extension icecon.hdf.
c
c     Usage: convert_icecon_tohdf.exe L1B_1km_file Binary_Ice_Conc_file
c                                     Output_hdf_file
c
c     Written by Kathleen Strabala March 2011
c-----------------------------------------------------------------------

      implicit none

c ... Include files
      include 'hdf.inc'

c ... Local variables
      integer file_id, hdf_lun, line, length, rec_mask, rec_conc
      integer varid_lat, varid_lon, varid_1km, lat_id, lon_id, icemask_id,
     +        icec_id
      integer rtn, num_arg, irt, iostat
      integer rank, dimsizes(2), data_type, num_attrs, dims(2)
      character*100 sds_name 
      character*40 dim_name(2)
      character*100 arrayname
      integer num_pixels_5km, num_pixels_1km, num_lines_5km, 
     +        num_lines_1km, num_scans
      double precision scale_fac, offset
      
      integer start(2), stride(2), edge(2)

      integer max_5km_pixel, max_1km_pixel
      parameter (max_5km_pixel = 300, max_1km_pixel = 1354)

      integer inlun
      parameter (inlun = 31)

      real float_buffer(max_5km_pixel)
      integer ice_mask(max_1km_pixel)
      real ice_conc(max_1km_pixel)
            
      integer dim

      character*132 in_file, in_bin, out_hdf

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
        print *, 'Usage: convert_icecon_tohdf.exe in_hdf in_bin out_hdf'
        print *, 'where'
        print *, 'in_hdf:  name of input L1B 1km HDF file'
        print *, 'in_bin:  name of input binary ice concentration file'
        print *, 'out_hdf: name of output ice concentration HDF file'
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
      print '(''Input binary IST file: '', a)', in_bin(1 : string_length(in_bin))
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

c ... Now create the Ice Mask SDS Array
      dims(1) = num_pixels_1km
      dims(2) = num_lines_1km
      rank = 2
      icemask_id = sfcreate(hdf_lun,'Ice_Mask',DFNT_INT32,rank,dims)
      if (icemask_id .lt. 0) then
         print *, 'Error: Unable to create ice mask HDF array'
      endif
      dim_name(1) = 'Cell_Across_Swath_1km'
      dim_name(2) = 'Cell_Along_Swath_1km'
      irt = sfsdmname(sfdimid(icemask_id,0),dim_name(1))
      irt = sfsdmname(sfdimid(icemask_id,1),dim_name(2))

c ... Attach local attributes to the SDS
      irt = sfsnatt(icemask_id,"units",DFNT_CHAR8,4,"None")
      scale_fac = 1.0d0
      irt = sfsnatt(icemask_id,'scale_factor',DFNT_FLOAT64,1,scale_fac)
      offset = 0.0d0
      irt = sfsnatt(icemask_id,'add_offset',DFNT_FLOAT64,1,offset)
      irt = sfsrange(icemask_id,2,-2)
      irt = sfsfill(icemask_id,-999)

c ... Now create the Ice Concentration SDS Array
      dims(1) = num_pixels_1km
      dims(2) = num_lines_1km
      rank = 2
      icec_id = sfcreate(hdf_lun,'Ice_Concentration',DFNT_FLOAT32,rank,dims)
      if (icec_id .lt. 0) then
         print *, 'Error: Unable to create ice concentration HDF array'
      endif
      dim_name(1) = 'Cell_Across_Swath_1km'
      dim_name(2) = 'Cell_Along_Swath_1km'
      irt = sfsdmname(sfdimid(icec_id,0),dim_name(1))
      irt = sfsdmname(sfdimid(icec_id,1),dim_name(2))

c ... Attach local attributes to the SDS
      irt = sfsnatt(icec_id,"units",DFNT_CHAR8,1,"%")
      scale_fac = 1.0d0
      irt = sfsnatt(icec_id,'scale_factor',DFNT_FLOAT64,1,scale_fac)
      offset = 0.0d0
      irt = sfsnatt(icec_id,'add_offset',DFNT_FLOAT64,1,offset)
      irt = sfsrange(icec_id,100.0,0.0)
      irt = sfsfill(icec_id,-999.0)


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

c ... Open the input ice concentration binary file
c ... 
      rec_mask = num_lines_1km * 2 + 1
      rec_conc = num_lines_1km * 3 + 1
      length = num_pixels_1km * 4
      open(inlun, file=in_bin, iostat=iostat, status='unknown',
     &  access='direct', form='unformatted', recl=length)
      print '('' '')'
      print '(''Reading Ice Mask and Ice Concentration binary file '')'
      print '(''------------------------------------------------- '')'
      print '('' '')'
      do line = 0, (num_lines_1km - 1), 1

c ...   Report progress on every 10th scan
        if (mod(line, 100) .eq. 0) print '(''Scan: '', i6)', line / 10
  
c ...   Read the ice mask variable out of the file
        read(inlun, rec=rec_mask, iostat=iostat) ice_mask
        if (iostat .ne. 0) then
          print *, 'Error: Read of input Ice Mask failed'
          call exit(-1)
        endif
        rec_mask=rec_mask + 1

c ...   Read the ice concentration variable out of the file
        read(inlun,rec=rec_conc,iostat=iostat) ice_conc
        if (iostat .ne. 0) then
          print *, 'Error: Read of input Ice Conctration failed'
          call exit(-1)
        endif
        rec_conc = rec_conc + 1

c ...   Set start vector
        start(1) = 0
        start(2) = line

c ...   Set stride vector
        stride(1) = 1
        stride(2) = 1

c ...   Set edge vector
        edge(1) = num_pixels_1km
        edge(2) = 1

c ...   Write the Ice Mask data to the output HDF4 file (long integers)
        rtn = sfwdata(icemask_id, start, stride, edge, ice_mask)
        if (rtn .eq. -1) then
          print *, 'Error writing Ice Mask data to output HDF file'
          call exit(-1)
        endif

c ...   Write the Ice Concentration data to the output HDF4 file (float values)
        rtn = sfwdata(icec_id, start, stride, edge, ice_conc)
        if (rtn .eq. -1) then
          print *, 'Error writing Ice Concentration data to output HDF file'
          call exit(-1)
        endif
     
      enddo 

      close(inlun)

      print '(''Finished Writing Ice Mask and Ice Concentration 1km array '')'

c ... Close new HDF file
      irt = sfend(hdf_lun)

      end 


