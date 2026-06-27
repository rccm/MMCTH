      SUBROUTINE MOD07_OPEN_FILES( CUBEFILE, GEOFILE, CLDFILE, OUTFILE,
     &  BUGFILE, CUBE_HANDLE, GEO_HANDLE, CLD_HANDLE, BUG_LUN )

c-----------------------------------------------------------------------
c !F77
c
c !DESCRIPTION:
c     Get input and output filenames;
c     Open input and output files;
c     Check file times against PCF collection time.
c
c !INPUT PARAMETERS:
c     None
c
c !OUTPUT PARAMETERS:
c     CUBEFILE       Name of MODIS L1B 1km granule file
c     GEOFILE        Name of MODIS L1B geolocation file
c     CLDFILE        Name of MODIS L2 cloud mask file
c     OUTFILE        Name of MODIS L2 product file
c     BUGFILE        Name of runtime debug file
c     CUBE_HANDLE    Integer array returned by OPMFIL for L1B file
c     GEO_HANDLE     Integer array returned by OPMFIL for geolocation
c     CLD_HANDLE     Integer array returned by OPMFIL for cloudmask
c     BUG_LUN        FORTRAN logical unit number for runtime debug file
c
c !REVISION HISTORY:
C 12 apr 2004: SWS changed file_length of regression coefficients 
c              for 335 elements in RTV (15 emis)
c
C 4 june 2005: SWS changed file_length of regression coefficients 
c              for 321 elements in RTV (8 emis)
C Nov 2006: EvaB changed file_length for 323 elements in RTV (10 emis)
c
c !TEAM-UNIQUE HEADER:
c     Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
c
c !END
c-----------------------------------------------------------------------

      IMPLICIT NONE
      
      INCLUDE 'debug.inc'
      INCLUDE 'mod07_pcfnum.inc'
      INCLUDE 'mapi.inc'
      INCLUDE 'PGS_SMF.f'
      INCLUDE 'PGS_IO.f'
      
c ... Arguments

      CHARACTER*(*) cubefile, geofile, cldfile, outfile, bugfile
      INTEGER cube_handle(*), geo_handle(*), cld_handle(*), bug_lun

c ... Local variables

      INTEGER ier, file_id, file_version, file_length
      CHARACTER*32 routine_name
      PARAMETER ( routine_name = 'mod07_open_files' )
      LOGICAL open_flag

c ... External functions

      INTEGER pgs_pc_getreference, pgs_io_gen_openf, chk_input_l2
      EXTERNAL pgs_pc_getreference, pgs_io_gen_openf, chk_input_l2,
     &  opmfil      

      integer num_args
      integer FlagRA
      character FlagBuff*10
      num_args = command_argument_count()
      
      if(num_args == 1) then
         call get_command_argument(1,FlagBuff)
         read (FlagBuff,*) FlagRA
      else
      !This is the default value
         FlagRA = 0
      endif      
      
c-----------------------------------------------------------------------
c    GET NAMES OF INPUT FILES
c-----------------------------------------------------------------------

c ... Get Level 1B scan cube file name from process control file
      if( FlagRA == 1) then
         file_id = cube_pcfnum_RA
      else
         file_id = cube_pcfnum
      endif
      file_version = 1
      cubefile = ' '
      ier = pgs_pc_getreference( file_id, file_version, cubefile )
      if( ier .ne. PGS_S_SUCCESS) call message( routine_name,
     &  'Error getting MODIS L1B filename from PCF' //
     &  ' [OPERATOR ACTION: Check PCF entry]', 0, 2 )

c ... Get geolocation file name from process control file

      file_id = geo_pcfnum
      file_version = 1
      geofile = ' '
      ier = pgs_pc_getreference( file_id, file_version, geofile )
      if( ier .ne. PGS_S_SUCCESS ) call message( routine_name,
     &  'Error getting MODIS geolocation filename from PCF' //
     &  ' [OPERATOR ACTION: Check PCF entry]', 0, 2 )

c ... Get cloud mask file name from process control file

      file_id = cld_pcfnum
      file_version = 1
      cldfile = ' '
      ier = pgs_pc_getreference( file_id, file_version, cldfile )
      if( ier .ne. PGS_S_SUCCESS ) call message( routine_name,
     &  'Error getting MODIS cloud mask filename from PCF' //
     &  ' [OPERATOR ACTION: Check PCF entry]', 0, 2 )

c ... Get output file name from process control file

      file_id = out_pcfnum
      file_version = 1
      outfile = ' '
      ier = pgs_pc_getreference( file_id, file_version, outfile )
      if( ier .ne. PGS_S_SUCCESS ) call message( routine_name,
     &  'Error getting MODIS output filename from PCF' //
     &  ' [OPERATOR ACTION: Check PCF entry]', 0, 2 )

c ... Get debug file name from process control file

      file_id = bug_pcfnum
      file_version = 1
      bugfile = ' '
      ier = pgs_pc_getreference( file_id, file_version, bugfile )
      if( ier .ne. PGS_S_SUCCESS ) call message( routine_name,
     &  'Error getting debug filename from PCF' //
     &  ' [OPERATOR ACTION: Check PCF entry]', 0, 2 )

c-----------------------------------------------------------------------
c    OPEN INPUT FILES
c-----------------------------------------------------------------------

c ... Open L1B file and get file handle

      if ( opmfil( cubefile, 'r', cube_handle ) .ne.
     &  MAPIOK ) call message( routine_name,
     &  'Error opening MODIS L1B input file' //
     &  ' [OPERATOR ACTION: Check that file exists]', 0, 2 )

c ... Open geolocation file and get file handle

      if ( opmfil( geofile, 'r', geo_handle ) .ne.
     &  MAPIOK ) call message( routine_name,
     &  'Error opening MODIS geolocation input file' //
     &  ' [OPERATOR ACTION: Check that file exists]', 0, 2 )

c ... Open cloudmask file and get file handle

      if ( opmfil( cldfile, 'r', cld_handle ) .ne.
     &  MAPIOK ) call message( routine_name,
     &  'Error opening MODIS cloudmask input file' //
     &  ' [OPERATOR ACTION: Check that file exists]', 0, 2 )

c ... Open debug file and get logical unit number

      file_length = 1
      file_version = 1
      bug_lun = -1
      ier = pgs_io_gen_openf( bug_pcfnum, PGSd_IO_Gen_WSeqFrm,
     &  file_length, bug_lun, file_version )
      if ( ier .ne. PGS_S_SUCCESS ) call message( routine_name,
     &  'Error opening new debug file' //
     &  ' [OPERATOR ACTION: Check that file creation is possible]',
     &  0, 2 )

c ... Open view angle and regression coefficient files
c ... (IUANG and IUCTW are declared in mod07_pcfnum.inc)

      file_length = 2720
      file_version = 1
      call modis_io_gen_openf( ang_pcfnum, PGSD_IO_GEN_RDIRUNF,
     &  file_length, iuang, file_version )
      inquire( iuang, opened = open_flag )
      if ( .not. open_flag )
     &  call message( 'mod07_compute_products',
     &  'Error opening view angle file' //
     &  ' [OPERATOR ACTION: Check that file exists]', 0, 2 )

C 4 june 2005: SWS changed file_length for 321 elements in RTV (8 emis)
C Nov 2006: EvaB changed file_length for 323 elements in RTV (10 emis)
c      file_length = 1300 
      file_length = 1308 
      file_version = 1
      call modis_io_gen_openf( trc_pcfnum, PGSD_IO_GEN_RDIRUNF,
     &  file_length, iuctw, file_version )
      inquire( iuctw, opened = open_flag )
      if ( .not. open_flag )
     &  call message( 'mod07_compute_products',
     &  'Error opening regression coefficient file' //
     &  ' [OPERATOR ACTION: Check that file exists]', 0, 2 )

c-----------------------------------------------------------------------
c    CHECK TIMES OF INPUT FILES
c-----------------------------------------------------------------------

c ... Check time of L1B file against PCF collection time

      file_version = 1
      if( FlagRA == 1) then
         ier = chk_input_l2( cube_pcfnum_RA, file_version )
      else
         ier = chk_input_l2( cube_pcfnum, file_version )
      endif
      if ( ier .eq. -1 ) call message( routine_name,
     &  'Input MODIS L1B file not within PCF collection time' //
     &  ' [OPERATOR ACTION: Check PCF entry, or stage correct file]',
     &  0, 2 )

c ... Check time of geolocation file against PCF collection time

      file_version = 1
      ier = chk_input_l2( geo_pcfnum, file_version )
      if ( ier .eq. -1 ) call message( routine_name,
     &  'Input MODIS geolocation file not within PCF collection time' //
     &  ' [OPERATOR ACTION: Check PCF entry, or stage correct file]',
     &  0, 2 )

c ... Check time of cloudmask file against PCF collection time

      file_version = 1
      ier = chk_input_l2( cld_pcfnum, file_version )
      if ( ier .eq. -1 ) call message( routine_name,
     &  'Input MODIS cloudmask file not within PCF collection time' //
     &  ' [OPERATOR ACTION: Check PCF entry, or stage correct file]',
     &  0, 2 )

c ... Write debug information

      if ( debug .eq. 1 ) then
        write( h_output, '(72(''-''))' )
        write( h_output, '(a,/)' ) 'mod07_open_files debug'
        write( h_output, '(a)' ) 'L1B granule file = ', cubefile
        write( h_output, '(a)' ) 'Geolocation file = ', geofile
        write( h_output, '(a)' ) 'Cloudmask file   = ', cldfile
        write( h_output, '(a)' ) 'Output file      = ', outfile
        write( h_output, '(a)' ) 'Debug file       = ', bugfile
      endif

      END
