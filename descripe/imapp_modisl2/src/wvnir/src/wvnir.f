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

      Program wvnir
      implicit none
C-----------------------------------------------------------------------
C
C!F77
C
C!Description:
C    This is the main program for the Direct Broadcast
C    version of the NIR water vapor product.
C
C!Input Parameters:
C
C The program reads a list of input files from the wvnir.cfg file.  
C   
C  Required 
C   Inputs Are:  MODIS 1 km flat file radiances and reflectances (MOD)21)
C                MODIS flat file observation geometry parameters (MOD03)
C
C   Inputs provided with the code are:
C                ASCII files specifying neural network coefficients:
C                   Minima and maxima of training data input
C                   Coefficients of two matrizes
C                   Scaling coefficients of input and output data
C
C!Output Parameters:
C
C  Output: Writes binary flat file of wvnir.
C
C!Usage:
C  wvnir.exe wvnir.cfg
C
C!References and Credits:
C   See Albert et. al 2004 (TBD) (modwvnir.pdf)
C   This code was created by modifying the source for IMAPP SST
C 
C!End
C
C-----------------------------------------------------------------------


C ... Include files
C ... External variable declarations

      include 'modwvnir_data.inc'

C ... Define neural network data

      include 'modwvnir_nn.inc'

C ... Internal variable declarations

C ... Constants

      character*18 routine_name
      parameter (routine_name = 'wvnir')
      
C ... Local arrays

      character*132 cfgname
      character*(PATH_MAX) l1b_1km_hdr,l1b_geo_hdr
      character*5 coutput_type
      character*4 interleave_1km,interleave_geo
      character*80 bandnames_1km(INBAND),bandnames_geo(20),
     &             bandunits_1km(INBAND),bandunits_geo(20)
      character*5 platform_name
      character*72 text

C ... Local scalars

      integer  num_arg, iargc, output_type
      integer  beg_lin, nlins, beg_ele, neles
      integer  l1b_1km_lun, geo_1km_lun,
     &         modwvnir_lun, hdr_lun, hdf_lun
      integer  nscans, npixels, datatype_1km, datatype_geo,
     &         resolution_1km, resolution_geo, 
     &         offset_1km, offset_geo, samples_1km, 
     &         samples_geo, bands_1km, bands_geo, 
     &         lines_1km, lines_geo
      integer  beg_scan, ibes
      integer  scan, iostat
      integer  line, pixel
      real     error_1km, error_geo

C#######################################################################
C     Initialization
C#######################################################################

C ... Required bands

      bands(1) = 2
      bands(2) = 17
      bands(3) = 18
      bands(4) = 19
      
C ... Main program

C ... Check number of arguments

      num_arg = iargc()
      if (num_arg .ne. 2) then
        print *, 'Usage: wvnir.exe cfgname output_type'
        print *, '       cfgname: name of WVNIR configuration file'
        print *, 'output_type: Format of output file'
        print *, '              1 = binary only '
        print *, '              2 = hdf only '
        print *, '              3 = binary and hdf '
        call exit(-1)
      endif

C ... Extract arguments:
C ... cfgname is the configuration file containing source
C ... input and output file information.

      call getarg(1, cfgname)
      call getarg(2, coutput_type)
      read(coutput_type,'(I1)',iostat=iostat) output_type

C ... Get debug value and processing interval information from PCF

      call modwvnir_get_rp(cfgname,debug,l1b_1km_hdr,l1b_geo_hdr,
     &                     beg_lin,nlins,beg_ele,neles)

C ... Open required files for this processing section

      call modwvnir_file_open(cfgname,neles,output_type,l1b_1km_lun,
     &                        geo_1km_lun,modwvnir_lun,hdr_lun,
     &                        hdf_lun,debug)

C ... Get number of scans and pixels for 1km bands

      call modwvnir_get_metadata(l1b_1km_hdr,l1b_geo_hdr,
     &     nscans,npixels,datatype_1km,datatype_geo,
     &     interleave_1km,interleave_geo,resolution_1km,
     &     resolution_geo,offset_1km,offset_geo,
     &     samples_1km,samples_geo,lines_1km,
     &     lines_geo,error_1km,error_geo,bands_1km,
     &     bands_geo,bandnames_1km,bandnames_geo,
     &     bandunits_1km,bandunits_geo,
     &     platform_name,
     &     debug)

C ... Check User defined processing intervals versus scan boundaries

      call modwvnir_chk_input(nscans,npixels,beg_lin,nlins,beg_ele,
     &                        neles,beg_scan,ibes, debug)

C#######################################################################
C     Scan processing
C#######################################################################
      
      do scan = beg_scan, nscans

        write(text, '(''INFO: Processing scan # '',i6)') scan
        call message( 'wvnir', text , 0, debug )

C ...   Extract radiance, geolocation, mask, and ancillary data

	call modwvnir_get_data(l1b_1km_lun,geo_1km_lun,datatype_1km,
     &       datatype_geo,interleave_1km,interleave_geo,
     &       resolution_1km,resolution_geo,
     &       offset_1km,offset_geo,samples_1km,
     &       samples_geo,lines_1km,lines_geo,
     &       error_1km,error_geo,bands_1km,bands_geo,
     &       bandnames_1km,bandnames_geo,bandunits_1km,
     &       bandunits_geo,scan, debug)

        call modwvnir_initialize_output()
 
C ...   Start loop over pixel boxes

        do line = 1 , max_line
          do pixel = ibes , npixels
	  
            call modwvnir_create_nn_input(line,pixel)
            call modwvnir_compute_product(line,pixel)
          end do
        end do

C ...   Write products for this scan

        call modwvnir_write_products( output_type,modwvnir_lun,hdr_lun,
     &                                hdf_lun,scan,nscans,platform_name)
C ... End of loop over all scans
      end do

C#######################################################################
C     End science processing
C#######################################################################

c ... close files for this section of code

      call modwvnir_file_close(l1b_1km_lun,geo_1km_lun,modwvnir_lun,
     &                         hdr_lun,output_type,hdf_lun)

c ... Write debug info and then close debug file

      call message( 'wvnir', 'WVNIR Processing completed' , 0, debug )

      end
      
