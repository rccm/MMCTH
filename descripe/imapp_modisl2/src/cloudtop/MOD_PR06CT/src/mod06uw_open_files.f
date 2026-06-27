      subroutine mod06uw_open_files( cubefile, geofile, cldfile,
     &  outfile, bugfile, cube_handle, geo_handle, cld_handle, 
     &  out_handle, eco1_lun, bug_lun, error_flag )

c-----------------------------------------------------------------------
c
c!F77
c
c!Description:
c    Get input and output filenames;
c    Open input and output files;
c
c!Input Parameters:
c    None
c
c!Output Parameters:
c    CUBEFILE       Name of MODIS L1B 1km granule file
c    GEOFILE        Name of MODIS L1B geolocation file
c    CLDFILE        Name of MODIS L2 cloud mask file
c    OUTFILE        Name of MODIS L2 product file
c    BUGFILE        Name of runtime debug file
c    CUBE_HANDLE    Integer array returned by OPMFIL for L1B file
c    GEO_HANDLE     Integer array returned by OPMFIL for geolocation
c    CLD_HANDLE     Integer array returned by OPMFIL for cloudmask
c    OUT_HANDLE     Integer arrya returned by OPMFILE for output file
c    ECO1_LUN       FORTRAN logical unit number for ecosystem 1km file
c    BUG_LUN        FORTRAN logical unit number for runtime debug file
c    ERROR_FLAG     1 - Error, 0 - completed normally
c
c!Revision History:
c
c!Team-unique Header:
c
c!End
c
c-----------------------------------------------------------------------

      implicit none


      include 'PGS_SMF.f'
      include 'PGS_IO.f'
      include 'mapi.inc'
      include 'mod06uw_pcfnum.inc'
      include 'mod06uw_debug.inc'
      
c ... arguments

      character*(*) cubefile, geofile, cldfile, outfile, bugfile
      integer cube_handle(MODFILLEN),geo_handle(MODFILLEN),
     &        cld_handle(MODFILLEN),out_handle(MODFILLEN), 
     &        bug_lun, eco1_lun, error_flag

c ... local variables

      integer ier, file_id, file_version, file_length, file_access, i
      logical open_flag

c whuang 05/07/01: Added save statement for common / mod06_debug /
      save / mod06_debug / 

c ... set program name for error messaging
      character*32 routine_name
      parameter ( routine_name = 'mod06uw_open_files' )

c ... external functions

      integer pgs_pc_getreference, pgs_io_gen_openf, chk_input_l2
      external pgs_pc_getreference, pgs_io_gen_openf, chk_input_l2,
     &  opmfil      

      integer num_args
      integer FlagRA
      character FlagBuff*10
      integer iargc

      num_args = iargc ( )
      
      if(num_args .eq. 1) then
         call getarg ( 1, FlagBuff )
         read (FlagBuff,*) FlagRA
      else
C     This is the default value
         FlagRA = 0
      endif

c ... Initialize files which need mapi commands to open
      do i = 1 , MODFILLEN
        cube_handle(i) = 0
        geo_handle(i) = 0
        cld_handle(i) = 0
        out_handle(i) = 0
      enddo

c ... initialize successful open flag
      open_flag = .false.
            
c ... get Level 1B scan cube file name from process control file
      if( FlagRA .eq. 1) then
         file_id = dscube_pcfnum_RA
      else
         file_id = dscube_pcfnum
      endif
      file_version = 1
      cubefile = ' '
      ier = pgs_pc_getreference( file_id, file_version, cubefile )
      if( ier .ne. PGS_S_SUCCESS) then
         call message( routine_name,
     &  'Error getting scan cube filename from PCF.' //
     &  ' [OPERATOR ACTION: Check system resources]', 0, 2 )
        error_flag = 1
      endif

c ... get geolocation file name from process control file

      file_id = geo_pcfnum
      file_version = 1
      geofile = ' '
      ier = pgs_pc_getreference( file_id, file_version, geofile )
      if( ier .ne. PGS_S_SUCCESS ) then
        call message( routine_name,
     &  'Error getting geolocation filename from PCF.' //
     &  ' [OPERATOR ACTION: Check system resources]', 0, 2 )
        error_flag = 1
      endif

c ... get cloud mask file name from process control file
      file_id = cld_pcfnum
      file_version = 1
      cldfile = ' '
      ier = pgs_pc_getreference( file_id, file_version, cldfile )
      if( ier .ne. PGS_S_SUCCESS ) then
        call message( routine_name,
     &  'Error getting cloud mask filename from PCF.' //
     &  ' [OPERATOR ACTION: Check system resources]', 0, 2 )
        error_flag = 1
      endif

c ... get output file name from process control file

      file_id = out_pcfnum
      file_version = 1
      outfile = ' '
      ier = pgs_pc_getreference( file_id, file_version, outfile )
      if( ier .ne. PGS_S_SUCCESS ) then
        call message( routine_name,
     &  'Error getting output filename from PCF.' //
     &  ' [OPERATOR ACTION: Check system resources]', 0, 2 )
        error_flag = 1
      endif

c ... get debug file name from process control file

      file_id = bug_pcfnum
      file_version = 1
      bugfile = ' '
      ier = pgs_pc_getreference( file_id, file_version, bugfile )
      if( ier .ne. PGS_S_SUCCESS ) then
        call message( routine_name,
     &  'Error getting debug filename from PCF. ' //
     &  ' [OPERATOR ACTION: Check system resources]', 0, 1 )
        error_flag = 1
      endif

c ... Open ecosystem files

c ... Open the global 1 km file
      file_version = 1
      file_access = PGSd_IO_Gen_RDirUnf
      file_length = 512

      ier = pgs_io_gen_openf(eco1_pcfnum, file_access,
     &    file_length, eco1_lun, file_version)
      if( ier .eq. PGS_S_SUCCESS ) then
        call message( routine_name,
     &  'Success opening global 1km ecosystem file ',
     &  0, 3 )
        open_flag = .false.
        inquire(eco1_lun,opened=open_flag)
        if (.not. open_flag) then
          call message( routine_name,
     &    'Could not open open global 1km ecosystem file. ' //
     &    ' [OPERATOR ACTION: Check input file and PCF file]',
     &    0, 3 )
         endif
      else
        call message( routine_name,
     &  'Could not open open global 1km ecosystem file. ' //
     &  ' [OPERATOR ACTION: Check system resources]',
     &  0, 1 )
        error_flag = 1
      endif


c ... Open L1B file and get file handle

      if ( opmfil( cubefile, 'r', cube_handle ) .ne.
     &  MAPIOK ) then
        call message( routine_name,
     &  'Error opening MODIS L1B 1km input file.' //
     &  ' [OPERATOR ACTION: Check system resources]', 0, 2 )
        error_flag = 1
      endif

c ... Open geolocation file and get file handle

      if ( opmfil( geofile, 'r', geo_handle ) .ne.
     &  MAPIOK ) then
        call message( routine_name,
     &  'Error opening MODIS geolocation input file.' //
     &  ' [OPERATOR ACTION: Check system resources]', 0, 2 )
        error_flag = 1
      endif

c ... Open cloudmask file and get file handle

      if ( opmfil( cldfile, 'r', cld_handle ) .ne.
     &  MAPIOK ) then
        call message( routine_name,
     &  'Error opening cloudmask input file.' //
     &  ' [OPERATOR ACTION: Check system resources]', 0, 2 )
        error_flag = 1
      endif

c ... Open debug file and get logical unit number

      file_length = 1
      file_version = 1
      bug_lun = -1
      ier = pgs_io_gen_openf( bug_pcfnum, PGSd_IO_Gen_WSeqFrm,
     &  file_length, bug_lun, file_version )
      if ( ier .ne. PGS_S_SUCCESS ) then
        call message( routine_name,
     &  'Error opening output debug file.' //
     &  ' [OPERATOR ACTION: Check system resources]', 0, 1 )
        error_flag = 1
      endif
      open_flag = .false.
      inquire(h_output,opened=open_flag)
      if ( .not. open_flag) then
        call message( routine_name,
     &  'Could not open output debug file. ' //
     &  ' [OPERATOR ACTION: Check input file and PCF file]',
     &  0, 2 )
       endif


c ... Open output file and get file with update status 

      if ( opmfil( outfile, 'a', out_handle ) .ne.
     &  MAPIOK ) then
        call message( routine_name,
     &  'Error opening output product file.' //
     &  ' [OPERATOR ACTION: Check system resources]', 0, 2 )
        error_flag = 1
      endif

c ... Check time of L1B file against PCF collection time

      file_version = 1
      if( FlagRA .eq. 1) then
         ier = chk_input_l2( dscube_pcfnum_RA, file_version )
      else
         ier = chk_input_l2( dscube_pcfnum, file_version )
      endif
      if ( ier .eq. -1 ) then
        call message( routine_name,
     &  'Granule L1B 1km time does not match processing interval.' //
     &  ' [OPERATOR ACTION: Check system resources]',
     &  0, 2 )
        error_flag = 1
      endif

c ... Check time of geolocation file against PCF collection time

      file_version = 1
      ier = chk_input_l2( geo_pcfnum, file_version )
      if ( ier .eq. -1 ) then
        call message( routine_name,
     &  'Granule Geo. time does not match processing interval.' //
     &  ' [OPERATOR ACTION: Check system resources]',
     &  0, 2 )
        error_flag = 1
      endif

c ... Check time of cloudmask file against PCF collection time

      file_version = 1
      ier = chk_input_l2( cld_pcfnum, file_version )
      if ( ier .eq. -1 ) then
        call message( routine_name,
     &  'Granule C. Mask. time does not match processing interval.' //
     &  ' [OPERATOR ACTION: Check system resources]',
     &  0, 2 )
        error_flag = 1
      endif


c ... Write file information

      write( h_output, '(a,/)' ) '        MOD06CT QA Output File '
      write( h_output, '(72(''-''))' )
      write( h_output, '(a,/)' ) '             Input Files '
      write( h_output, '(72(''-''))' )
      write( h_output, '(a)' ) 'L1B granule file = ', cubefile
      write( h_output, '(a)' ) 'Geolocation file = ', geofile
      write( h_output, '(a)' ) 'Cloudmask file   = ', cldfile
      write( h_output, '(a,/)' ) '        MOD06CT QA Output File '
      write( h_output, '(a)' ) 'Output file      = ', outfile
c      write( h_output, '(a,/)' ) '        MOD06CT QA Output File '
c      write( h_output, '(a),/' ) 'Debug file       = ', bugfile

      end
