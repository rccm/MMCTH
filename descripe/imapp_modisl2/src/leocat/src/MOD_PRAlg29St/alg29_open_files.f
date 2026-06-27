      SUBROUTINE ALG29_OPEN_FILES( LEOFILE, CLDFILE,  
     &  BUGFILE, LEO_HANDLE, CLD_HANDLE, BUG_LUN )

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
c     LEOFILE        Name of MODIS L2 LEOCAT output file
c     CLDFILE        Name of MODIS L2 cloud mask file
c     BUGFILE        Name of runtime debug file
c     LEO_HANDLE     Integer array returned by OPMFIL for cloudmask
c     CLD_HANDLE     Integer array returned by OPMFIL for cloudmask
c     BUG_LUN        FORTRAN logical unit number for runtime debug file
c
c !REVISION HISTORY:
c
c !TEAM-UNIQUE HEADER:
c     Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
c
c !END
c-----------------------------------------------------------------------

      IMPLICIT NONE
      
      INCLUDE 'debug.inc'
      INCLUDE 'mod35.inc'
      INCLUDE 'mapi.inc'
      INCLUDE 'PGS_SMF.f'
      INCLUDE 'PGS_IO.f'
      INCLUDE 'PGS_PC.f'
      
c ... Arguments

      CHARACTER*(*) leofile, cldfile, bugfile
      INTEGER  cld_handle(*), leo_handle(*), bug_lun

c ... Local variables

      INTEGER ier, file_id, file_version, file_length
      CHARACTER*32 routine_name
      PARAMETER ( routine_name = 'alg29_open_files' )
      LOGICAL open_flag

c ... External functions

      INTEGER pgs_pc_getreference, pgs_io_gen_openf 
      EXTERNAL pgs_pc_getreference, pgs_io_gen_openf, 
     &  opmfil      
            
      debug=1
c-----------------------------------------------------------------------
c    GET NAMES OF INPUT FILES
c-----------------------------------------------------------------------

c ... Get LEOCAT cloud mask file name from process control file

      file_id = LRN_mod35_LEOCAT
      file_version = 1
      leofile = ' '
      ier = pgs_pc_getreference( file_id, file_version, leofile )
      if( ier .ne. PGS_S_SUCCESS ) call message( routine_name,
     &  'Error getting MODIS cloud mask filename from PCF' //
     &  ' [OPERATOR ACTION: Check PCF entry]', 0, 2 )

c ... Get Collection 5 cloud mask file name from process control file
      file_id = LRN_mod35
      file_version = 1
      cldfile = ' '
      ier = pgs_pc_getreference( file_id, file_version, cldfile )
      if( ier .ne. PGS_S_SUCCESS ) call message( routine_name,
     &  'Error getting MODIS cloud mask filename from PCF' //
     &  ' [OPERATOR ACTION: Check PCF entry]', 0, 2 )

cc ... Get debug file name from process control file
c
c      file_id = LRN_out
c      file_version = 1
c      bugfile = ' '
c      ier = pgs_pc_getreference( file_id, file_version, bugfile )
c      if( ier .ne. PGS_S_SUCCESS ) call message( routine_name,
c     &  'Error getting debug filename from PCF' //
c     &  ' [OPERATOR ACTION: Check PCF entry]', 0, 2 )

c-----------------------------------------------------------------------
c    OPEN INPUT FILES
c-----------------------------------------------------------------------

c ... Open C6 cloudmask file and get file handle

      if ( opmfil( leofile, 'r', leo_handle ) .ne.
     &  MAPIOK ) call message( routine_name,
     &  'Error opening MODIS C6 cloudmask input file' //
     &  ' [OPERATOR ACTION: Check that file exists]', 0, 1 )

c ... Open C5 cloudmask file and get file handle

      if ( opmfil( cldfile, 'a', cld_handle ) .ne.
     &  MAPIOK ) call message( routine_name,
     &  'Error opening MODIS cloudmask input file' //
     &  ' [OPERATOR ACTION: Check that file exists]', 0, 2 )

cc ... Open debug file and get logical unit number
c
c      file_length = 1
c      file_version = 1
c      bug_lun = -1
c      ier = pgs_io_gen_openf( LRN_out, PGSd_IO_Gen_WSeqFrm,
c     &  file_length, bug_lun, file_version )
c      if ( ier .ne. PGS_S_SUCCESS ) call message( routine_name,
c     &  'Error opening new debug file' //
c     &  ' [OPERATOR ACTION: Check that file creation is possible]',
c     &  0, 2 )

c ... Write debug information

      if ( debug .eq. 1 ) then
        write( bug_lun, '(72(''-''))' )
        write( bug_lun, '(a,/)' ) 'alg29_open_files debug'
        write( bug_lun, '(a)' ) 'C6 Cloudmask file   = ', leofile
        write( bug_lun, '(a)' ) 'C5 Cloudmask file   = ', cldfile
c        write( bug_lun, '(a)' ) 'Debug file       = ', bugfile
      endif

      END
