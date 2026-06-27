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
      Program sst
      implicit none

c-----------------------------------------------------------------------
c
c!F77
c
c!Description:
c    This is the main program for the Direct Broadcast
c    version of the Sea Surface Temperature product (MODIS product mod28).
c
c!Input Parameters:
c
c The program reads a list of input files from the sst.cfg file.  
c   
c  Required 
c   Inputs Are:  MODIS 1 km flat file radiances and reflectances
c                MODIS flat file geolocation parameters
c 
c  Optional 
c    Inputs Are: 1 degree OISST file of the week
c
c!Output Parameters:
c
c  Output: Writes binary flat file of sst.  For more specific
c          information on the arrays included in the output file,
c          see the SST_File_Description.txt file.
c
c!References and Credits:
c   See MODIS Infrared Sea Surface Temperature Algorithm ATBD-MOD-25.
c   This code was created by modifying the source for IMAPP cloudtop
c 
c!End
c
c-----------------------------------------------------------------------

C#######################################################################
C     DECLARATIONS
C#######################################################################

c ... UW include files

      include 'db_mod28uw_data.inc'
      include 'db_mod28uw_debug.inc'
      include 'platform_name.inc'


c ... output scalars

      real success_sst1, success_sst4, 
     &     error_1km, error_geo

c ... Local arrays

      character*5               coutput_type
      character*4 interleave_1km,interleave_geo
      character*5 sat_name, platform_1km
      character*80 bandnames_1km(INBAND),bandnames_geo(20),
     &             bandunits_1km(INBAND),bandunits_geo(20)
      character*132 cfgname
      
      character*(PATH_MAX) l1b_1km_hdr,l1b_geo_hdr
      integer icode( max_line*max_pixel )
      byte qual_flag_ir ( max_line*max_pixel ),
     &     conf_flag_ir ( max_line*max_pixel )


c ... Local scalars

      real     gfit, Chisq
      integer  nscans, npixels, scan, num_arg, iargc, 
     &         npix, ng_irp, ssctpr, slcd, smcd, shcd, sthncd,
     &         sthkcd, sopcd, scicd, ni, nw, nx, nu, qual_flag,
     &         conf_flag, maxl, clus, nmaxl, line, pixel,
     &         npix_scan, beg_lin, nlins, beg_ele, neles, ibes, 
     &         beg_scan, hdr_lun, l1b_1km_lun, geo_1km_lun,
     &         mod28_lun, datatype_1km, datatype_geo, 
     &         resol_mask, resolution_1km,resolution_geo,offset_1km,
     &         offset_geo,offset_mask,samples_1km,
     &         samples_geo,sampl_mask,bands_1km,hdf_lun,
     &         bands_geo,bands_mask,lines_1km,lines_geo,
     &         lines_mask,datat_mask,output_type,iostat

c ... Save Statements
      save /platform_name_common/
     
      save / mod28_data /, / mod28_debug / 

c ... External functions

      integer copy_info
      external copy_info
      
      
C#######################################################################
C     INITIALIZATION
C#######################################################################

c ... Data initialization

      bands(1) = 20
      bands(2) = 22
      bands(3) = 23
      bands(4) = 24
      bands(5) = 25
      bands(6) = 27
      bands(7) = 28
      bands(8) = 29
      bands(9) = 30
      bands(10) = 31
      bands(11) = 32
      bands(12) = 33
      bands(13) = 34
      bands(14) = 35
      bands(15) = 36
      bands(16) = 1
      bands(17) = 6
      bands(18) = 26


c ... list bands needed for SST algorithms     
      sstband(1) = 20
      sstband(2) = 22
      sstband(3) = 23
      sstband(4) = 31
      sstband(5) = 32
      
      sstbandptr(1) = 1
      sstbandptr(2) = 2
      sstbandptr(3) = 3
      sstbandptr(4) = 10
      sstbandptr(5) = 11

c ..  initialize the output metadata fields

      success_sst1 = 0.0
      success_sst4 = 0.0

      write( *, '(2x,''Initiated DB SST Processing  '')' ) 

c --- Get input argument of configuration file containing source
c ---  input and output file information.

c --- Check number of arguments
      num_arg = iargc()
      if (num_arg .ne. 3) then
        print *, 'Usage: sst.exe cfgname sat_name'
        print *, 'where'
        print *, 'cfgname: name of SST configuration file'
        print *, 'sat_name: MODIS satellite platform name (Aqua or Terra)'
        print *, 'output_type: Format of output file'
        print *, '              1 = binary only '
        print *, '              2 = hdf only '
        print *, '              3 = binary and hdf '
        call exit(-1)
      endif

c ... Extract arguments
      call getarg(1, cfgname)
      call getarg(2, sat_name)
      call getarg(3, coutput_type)
      read(coutput_type,'(I1)',iostat=iostat) output_type

c ... Check input platform name to make sure it is either Aqua or Terra
      platform_name = ' '
      if (sat_name(1:5) .eq. 'Terra'   .or.
     +    sat_name(1:5) .eq. 'terra'   .or.
     +    sat_name(1:5) .eq. 'TERRA') then
          platform_name = 'Terra'
      elseif (sat_name(1:4) .eq. 'Aqua'   .or.
     +    sat_name(1:4) .eq. 'aqua'   .or.
     +    sat_name(1:4) .eq. 'AQUA') then
          platform_name = 'Aqua'
      else
          print '(''Error: Incorrect Satellite Platform name entered '', a)',
     +              sat_name(1: len(sat_name))
          call exit(-1)
      endif

c ... Write out platform name
      print '('' For Satellite Platform '', a)',
     +              platform_name(1: len(platform_name))

c ... Get debug value and processing interval information from PCF
      call db_mod28_get_rp(cfgname,debug,l1b_1km_hdr,l1b_geo_hdr,
     +                     beg_lin,nlins,beg_ele,neles)

c ....................................................................

c *********************************************************************

c ... Open required files for this processing section

      call db_mod28_file_open(cfgname,neles,output_type,l1b_1km_lun,
     +                        geo_1km_lun,mod28_lun,hdr_lun,hdf_lun,
     +                        h_output)

c ... Get number of scans and pixels for 1km bands

      call db_mod28_Get_Metadata(l1b_1km_hdr,l1b_geo_hdr,
     &     nscans,npixels,datatype_1km,datatype_geo,
     &     interleave_1km,interleave_geo,resolution_1km,
     &     resolution_geo,offset_1km,offset_geo,
     &     samples_1km,samples_geo,lines_1km,
     &     lines_geo,error_1km,error_geo,bands_1km,
     &     bands_geo,bandnames_1km,bandnames_geo,
     &     bandunits_1km,bandunits_geo,platform_1km)


c ... Check User defined processing intervals versus scan boundaries

      call db_mod28_chk_input(nscans,npixels,beg_lin,nlins,beg_ele,
     +                        neles,beg_scan,ibes,platform_1km)

C#######################################################################
C     INITIALIZATION OF GRANULE BASED VARIABLES
C#######################################################################

      npix = 0

c ... Initialize metadata counters for SST product


C#######################################################################
C     START SCIENCE PROCESSING
C#######################################################################

c ... Begin scan loop

      do scan = beg_scan , nscans


c ...   Write debug info

        write( *, '(''Processing scan # '',i4)' ) scan

c ...   Extract radiance, geolocation, mask, and ancillary data

        call db_mod28_get_data(l1b_1km_lun,geo_1km_lun,datatype_1km,
     &       datatype_geo,interleave_1km,interleave_geo,
     &       resolution_1km,resolution_geo,
     &       offset_1km,offset_geo,samples_1km,
     &       samples_geo,lines_1km,lines_geo,
     &       error_1km,error_geo,bands_1km,bands_geo,
     &       bandnames_1km,bandnames_geo,bandunits_1km,
     &       bandunits_geo,scan)

        call db_mod28_initialize_output()

        npix_scan = 0

c ...   Start loop over pixel boxes

        do line = 1 , max_line

          do pixel = ibes , npixels

             npix = npix + 1
             npix_scan = npix_scan + 1

             call db_mod28_compute_products(line,pixel,qual_flag)

          end do

        end do

c ...   Write products for this scan

        call db_mod28_write_products(output_type,mod28_lun,hdr_lun,
     +                               hdf_lun,scan,nscans)

c ... End of loop over all scans

      end do

C#######################################################################
C     END SCIENCE PROCESSING
C#######################################################################

c ... close files for this section of code

      call db_mod28_file_close(l1b_1km_lun,geo_1km_lun,mod28_lun,
     +                         hdr_lun,output_type,hdf_lun)

c ... Write debug info and then close debug file

      write( *, '(2x,''DB SST Processing Completed. '')' ) 

      close( h_output )

      end
