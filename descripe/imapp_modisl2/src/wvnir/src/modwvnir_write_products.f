C--------------------------------------------------------------------
C  Copyright (C) 2002,  Space Science and Engineering Center, 
C  University C  of Wisconsin-Madison, Madison WI.
C      
C  This program is free software; you can redistribute it 
C  and/or modify it under the terms of the GNU General 
C  Public License as published by the Free Software Foundation; 
C  either version 2 of the License, or (at your option) any 
C  later version.
C
C  This program is distributed in the hope that it will be 
C  useful, but WITHOUT ANY WARRANTY; without even the implied 
C  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
C  See the  GNU General Public License for more details.
C
C  You should have received a copy of the GNU General Public 
C  License along with this program; if not, write to the Free 
C  Software Foundation, Inc., 59 Temple Place, Suite 330, 
C  Boston, MA  02111-1307 USA
C--------------------------------------------------------------------
C
C
      subroutine modwvnir_write_products(output_type,out_handle,hdr_lun,
     &                                   hdf_lun,scan,nscans,
     &                                   platform_name)

c-----------------------------------------------------------------------
c
c!F77
c
c!Description:
c    Write NIR water vapor retrieval product to the output file.
c
c!Input Parameters:
c    OUTPUT_TYPE Flag defining type of output:  
c                         1 = binary only
c                         2 = hdf only 
c                         3 = binary and hdf
c    OUT_HANDLE    LUN for NIR water vapor output binary file
c    HDR_LUN       LUN for product biary header output file
C    HDF_LUN       LUN for wvnir output HDF file
c    SCAN          Scan number within L1B granule
c    NSCANS        Total number of scans in this data set
c
c!Output Parameters:
c    None
c
c    The following array in COMMON /MODWVNIR_DATA/ is written to output:
c    WVNIR_OUT      NIR water vapor output product
c
c!Revision History:
c
c!Team-unique Header:
c
c!End
c
c-----------------------------------------------------------------------
      
      implicit none


      include 'modwvnir_data.inc'
      include 'hdf.inc'
      
c ... arguments

      integer out_handle, scan, nscans, nboxes, hdr_lun,
     &        output_type, hdf_lun
      character*5 platform_name

c ... local variables
      character*40 dim_name(2)
      character*100 arrayname

      integer start(2), stride(2), dims(2), rank, data_type, wv_id, ii, 
     +        irt, lat_id, lon_id, ipxl, ilin, nlins, neles, iii, jjj, jj
      real    lat1_stor(max_pixel/5,MAX_LINE/5),
     +        lon1_stor(max_pixel/5,MAX_LINE/5),
     +        wvnir_true(max_pixel,MAX_LINE)
      logical bin_file, hdf_file, init
      double precision scale_fac, offset

c ... functions      
      integer sfcreate, sfwdata, sfend, sfsfill, sfsdmname,
     +        sfdimid, sfsnatt, sfscatt
      external sfcreate, sfwdata, sfend, sfsfill, sfsdmname,
     +        sfdimid, sfsnatt, sfscatt

c ... data statements     
      data              init /.true./

c ... Determine which files to produce
      bin_file = .false.
      hdf_file = .false.
      if (output_type .eq. 1  .or.  output_type .eq. 3) then
          bin_file = .true.
      endif
      if (output_type .eq. 2  .or.  output_type .eq. 3) then
          hdf_file = .true.
      endif

      nlins = nscans * MAX_LINE
      neles = max_pixel 

c ... NIR water vapour field
c ... This subroutine becomes useful when more than one dataset
c ... is written to the output file, e.g. quality flags ...

      if (bin_file) then
         call write_output( out_handle, hdr_lun, scan, nscans, 1,
     &      WVNIR_OUT,'WVNIR_LAND','kg/m^2', 
     &      platform_name)
      endif

c ... Now for HDF file
      if (init) then
        if (hdf_file) then
          dims(1) = neles/5
          dims(2) = nlins/5
          rank = 2
          lat_id = sfcreate(hdf_lun,'Latitude',DFNT_FLOAT32,rank,
     &             dims)
          if (lat_id .lt. 0) then
              call message( 'wvnir_write_products',
     &          'Could not create the latitude hdf array  ',
     &          0, 3 )
          endif
          dim_name(1) = 'Cell_Across_Swath_5km'
          dim_name(2) = 'Cell_Along_Swath_5km'
          irt = sfsdmname(sfdimid(lat_id,0),dim_name(1))
          irt = sfsdmname(sfdimid(lat_id,1),dim_name(2))
          irt = sfsfill(lat_id,-999.99)

          dims(1) = neles/5
          dims(2) = nlins/5
          rank = 2
          lon_id = sfcreate(hdf_lun,'Longitude',DFNT_FLOAT32,rank,
     &             dims)
          if (lon_id .lt. 0) then
              call message( 'wvnir_write_products',
     &          'Could not create the longitude hdf array  ',
     &          0, 3 )
          endif
          dim_name(1) = 'Cell_Across_Swath_5km'
          dim_name(2) = 'Cell_Along_Swath_5km'
          irt = sfsdmname(sfdimid(lon_id,0),dim_name(1))
          irt = sfsdmname(sfdimid(lon_id,1),dim_name(2))
          irt = sfsfill(lon_id,-999.99)

          dims(1) = neles
          dims(2) = nlins
          rank = 2
          wv_id = sfcreate(hdf_lun,'Water_Vapor_Direct_Broadcast_Near_Infrared',
     &             DFNT_FLOAT32,rank,dims)
          if (wv_id .lt. 0) then
              call message( 'wvnir_write_products',
     &          'Could not create the output wvnir hdf array  ',
     &          0, 3 )
          endif
          dim_name(1) = 'Cell_Across_Swath_1km'
          dim_name(2) = 'Cell_Along_Swath_1km'
          irt = sfsdmname(sfdimid(wv_id,0),dim_name(1))
          irt = sfsdmname(sfdimid(wv_id,1),dim_name(2))

c ...     Attach local attributes to the SDS
          irt = sfscatt(wv_id,'units',DFNT_CHAR8,5,'kg/m2')
          scale_fac = 1.0d0
          irt = sfsnatt(wv_id,'scale_factor',DFNT_FLOAT64,1,scale_fac)
          offset = 0.0d0
          irt = sfsnatt(wv_id,'add_offset',DFNT_FLOAT64,1,offset)
          irt = sfsfill(wv_id,-1.0)

        endif
        init = .false.
      endif

      if (hdf_file) then

         do jj = 1, MAX_LINE
           do ii = 1, neles
              if (abs(WVNIR_OUT(ii,jj) - bad_value ) .gt. 1.0e-5 ) then
                 wvnir_true(ii,jj) = WVNIR_OUT(ii,jj)
              else
                 wvnir_true(ii,jj) = -1.0
              endif
           enddo
         enddo

c ...    Write the NIRWV values to the output file
         arrayname = 'Water_Vapor_Direct_Broadcast_Near_Infrared'
         rank = 2
         start(1) = 0
         start(2) = (scan-1) * MAX_LINE
         dims(1) = neles
         dims(2) = MAX_LINE
         stride(1) = 1
         stride(2) = 1

         irt = sfwdata(wv_id,start,stride,dims,wvnir_true)

         if (irt .ne. 0) then
           call message( 'wvnir_write_products',
     &     'Could not write wvnir array to output hdf file ',
     &     0, 3 )
         endif

c ...    Now the lat/lon files
         jjj = 1
         do ilin = 3 , MAX_LINE, 5
           iii = 1
           do ipxl = 3 , neles-5 , 5
             lat1_stor(iii,jjj) = lat1(ipxl,ilin)
             lon1_stor(iii,jjj) = lon1(ipxl,ilin)
             iii = iii+1
           enddo
           jjj = jjj + 1
         enddo

         arrayname = 'Latitude'
         rank = 2
         start(1) = 0
         start(2) = ((scan-1) * MAX_LINE) / 5
         dims(1) = neles/5
         dims(2) = MAX_LINE/5
         stride(1) = 1
         stride(2) = 1

         irt = sfwdata(lat_id,start,stride,dims,lat1_stor)

         if (irt .ne. 0) then
           call message( 'mod28_write_products',
     &     'Could not write latitude array to output hdf file ',
     &     0, 3 )
         endif

         arrayname = 'Longitude'
         rank = 2
         start(1) = 0
         start(2) = ((scan-1) * MAX_LINE) / 5
         dims(1) = neles/5
         dims(2) = MAX_LINE/5
         stride(1) = 1
         stride(2) = 1

         irt = sfwdata(lon_id,start,stride,dims,lon1_stor)

         if (irt .ne. 0) then
           call message( 'mod28_write_products',
     &     'Could not write longitude array to output hdf file ',
     &     0, 3 )
         endif
      endif

      END

c-----------------------------------------------------------------------

      subroutine write_output(out_handle,hdr_lun,scan,nscans,
     &                        req_band,array,arrnm,b_unit_name, 
     &                        platform_name)
c-----------------------------------------------------------------------
c!F77
c
c!Description:
c    Write NIR water vapor retrieval product to the output file.
c
c!Input Parameters:
c    OUT_HANDLE    LUN for MOD28 output file
c    HDR_LUN       LUN for product header output file
c    SCAN          MODIS scan number
c    NSCANS        Total number of scans in this data set
c    REQ_BAND      Output file band number
c    ARRAY         Array of product values
c                  (Bad values should be set to BAD_VALUE in the calling
c                  routine: see modwvnir_data.inc for BAD_VALUE value)
c    ARRNM         Name of the array (SDS name) in the output file
c    b_unit_name   textual name of band for header file
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
      save

      
      include 'modwvnir_data.inc'
      
c ... arguments

      integer out_handle, scan, nscans, req_band, hdr_lun
      real array(max_pixel,max_line)
      character*(*) arrnm, b_unit_name
      logical init
      character*5 platform_name

c ... local variables
            
      integer i, j, out, ii, out_samples, out_lines, rtn, fil_datatype,
     &        fil_lun, fil_resolution, fil_offset, fil_samples, 
     &        fil_lines, fil_bands, out_maxsamples, out_maxlines,
     &        req_lines, req_samples, hdr_flg
      

      real fill_float, fil_error

c ... Local arrays
      byte out_flag(max_pixel*max_line)
      character*4 fil_interleave
      character*72 text
      character*10 fil_bandunits(OUT_BAND)
      character*80 fil_bandnames(OUT_BAND)
      character*256 fil_desc
      real float_data( max_pixel * max_line )

c ... external routines
      external modwvnir_write_flat_file
      integer modwvnir_write_flat_file

      integer string_length
      external string_length

      data init /.true./, ii /0/

      if (init) then
c ...   Set output file parameters
        req_samples = max_pixel
        req_lines = max_line

        fil_lun = out_handle
        fil_datatype = 4
        fil_interleave = 'bil'
        fil_resolution = 1
        fil_desc = 'Direct Broadcast NIR Water Vapor Land Product'
        fil_error = -1
        fil_offset = 0
        fil_lines = nscans * max_line
        fil_bands = OUT_BAND
        fil_samples = max_pixel
        out_maxsamples = max_pixel
        out_maxlines = max_line
        init = .false.
      endif

c ... Fill in band and band unit information
      ii = ii + 1
      if (ii .le. OUT_BAND) then
        fil_bandnames(ii) = arrnm
        fil_bandunits(ii) = b_unit_name(1:string_length(b_unit_name))
      endif

c ... If we have gathered all the band names, then write header
      if (ii .eq. OUT_BAND) then
          hdr_flg = 1
      else
          hdr_flg = 0
      endif

      out = 1

c ... initialize the validity flag holder
      do j = 1 , max_line
        out_flag(out) = 0
        out = out + 1
      enddo

c ... Get the fill value
      
      fill_float = -1.0
      
c ... Create the scaled output array

      out = 1
      
      do j = 1, max_line

        do i = 1, max_pixel

c ...     Set default output value (fill value)
        
          float_data( out ) = fill_float

c ...     If array value is good, scale it and store it
        
          if ( abs( array( i, j ) - bad_value ) .gt. 1.0e-5 ) then

            float_data( out ) =  array( i, j ) 

          endif
          
          out = out + 1

        end do

      end do

c ... Write the output array and header to the output files

c ------ write the mask array to a flat output file
      rtn = modwvnir_write_flat_file(
     &                         scan,
     &                         req_samples,
     &                         req_lines,
     &                         req_band,
     &                         fil_lun,
     &                         hdr_flg,
     &                         hdr_lun,
     &                         platform_name,
     &                         fil_desc,
     &                         fil_datatype,
     &                         fil_interleave,
     &                         fil_resolution,
     &                         fil_error,
     &                         fil_offset,
     &                         fil_samples,
     &                         fil_lines,
     &                         fil_bands,
     &                         fil_bandnames,
     &                         fil_bandunits,
     &                         out_maxsamples,
     &                         out_maxlines,
     &                         float_data,
     &                         out_flag,
     &                         out_samples,
     &                         out_lines
     &                       )

      if ( rtn .ne. 0 ) then
        write(*, '(''Write failed on output array '',a)') arrnm
        write(*, '(5x,''rtn status = '',I9)') rtn
        write(text, '(''Write failed on output array '',a)') arrnm
        call message( 'modwvnir_write_products', text ,
     &   0, 2 )
      endif
                
      end


