        subroutine mod06_Create_Get_Metadata( nscans, npixels )

        IMPLICIT NONE
        SAVE

C!F77------------------------------------------------------------------
C
C!DESCRIPTION:  Extracts Global metadata from the L1B granule which
c               is needed for MOD06 product generation.
c               Please note that calling this subroutine will cause
c               a warning message to be inserted into the LogStatus
c               file.  Please see description of message on line 73.
C
C!INPUT PARAMETERS: NONE
c
C!OUTPUT PARAMETERS:
C       INTEGER  nscans     Number of scans for this granule
C       INTEGER  npixels    Maximum number of pixels per scan in
c                           this particular granule
C
C!REVISION HISTORY:
c 11/19/98 Kathy Strabala
c Added OPERATOR disregard messages following calls
c to GMFIN. MAPI returns an error message when there
c really is no problem.
c
c 10/28/98 Kathy Strabala
c Started with Get_Global_Metadata.f from MOD06_PR06UW.
c Stripped out debug portion for use in just the create
c module for MOD06.
C
C!TEAM-UNIQUE HEADER:
C
C
C!DESIGN NOTES:
C
C    Externals:
C       Functions:
C          GMFIN            (mapi.inc)
c       Subroutines:
c          message          Writes to LogStatus file
C!END-----------------------------------------------------------------

      include 'mod06_pcfnum.inc'
      include 'mapi.inc'
      include 'PGS_PC.f'
      include 'PGS_SMF.f'

c ... Arguments

      integer nscans, npixels

c ... Local scalars

      character*72 dtype, attrib
      integer nelmnt, value, rtn, cube_handle(MODFILLEN)

      integer file_id, file_version
      character*(PGSd_PC_FILE_PATH_MAX) cubefile

c ... Set program name for error messaging

      character*32 routine_name
      parameter ( routine_name = 'mod06_Create_Get_Metadata' )

c ... External routines

      integer pgs_pc_getreference
      external pgs_pc_getreference

      external opmfil, clmfil

      integer num_args
      integer FlagRA
      character FlagBuff*10
      integer iargc

      num_args = iargc ( )
      
      if(num_args .eq. 1) then
         call getarg ( 1, FlagBuff )
         read (FlagBuff,*) FlagRA
      else
C      This is the default value
         FlagRA = 0
      endif


c -----------------------------------------------------------------

c ... get Level 1B scan cube file name from process control file

      if( FlagRA .eq. 1) then
               file_id = cube_pcfnum_RA
      else
         file_id = cube_pcfnum
      endif

      file_version = 1
      cubefile = ' '
      rtn = -1
      rtn = pgs_pc_getreference( file_id, file_version, cubefile )
      if( rtn .ne. PGS_S_SUCCESS) call message( routine_name,
     &  'Error getting scan cube filename from PCF.' //
     &  ' [OPERATOR ACTION: Check PCF entry]', 0, 2 )

c ... Open L1B file and get file handle

      if ( opmfil( cubefile, 'r', cube_handle ) .ne. MAPIOK )then
        call message( routine_name,
     &  'Error opening MODIS L1B 1km input file.' //
     &  ' [OPERATOR ACTION: Check that file exists]', 0, 2 )
      else
        call message( routine_name,
     &  'Successfully opened MODIS L1B 1km input file.', 0, 3 )
      endif

c -----------------------------------------------------------------
c     Want to get 2 pieces of information out of global metadata file
c     First get number of scans for this particular granule

      dtype = '  '
      attrib = 'Number of Scans'
      nelmnt = 0
      value = 0
      rtn = 0

C     'MTYPEf2c():MAPI_E_ERR:324431360 (PID=26923)
C     ERROR: MTYPEf2c unable to use empty input parameter.'
C     when empty inputs for variables 'dtype' and 'nelmnt' are passed to
C     GMFIN(). This is an unfortunate message since the M-API library
c     was designed to accept empty inputs, and to return their actual
c     values.
c     Make inital call to get correct input parameters
c     without hardwiring
      rtn = GMFIN(cube_handle,attrib,dtype,nelmnt,value)
      if (rtn .ne. MAPIOK) then
         call message( routine_name,
     &  ' [OPERATOR: Disregard previous MAPI message]', 0, 3 )
      endif

c     Second call will replace variables with correct inputs
      rtn = GMFIN(cube_handle,attrib,dtype,nelmnt,value)
      if (rtn .ne. 0) then
        call message( routine_name,
     &    'Error getting nscan info from L1B file' //
     &    ' [OPERATOR ACTION: Stage correct L1B, rerun PGE]',
     &    0, 2 )
      endif

      nscans = value

c -----------------------------------------------------------------
c     Second parameter is Max Earth Veiw Frames

      dtype = '  '
      attrib = 'Max Earth View Frames'
      nelmnt = 0
      value = 0
      rtn = 0

c     Make inital call to get correct input parameters
c     without hardwiring
      rtn = GMFIN(cube_handle,attrib,dtype,nelmnt,value)
      if (rtn .ne. MAPIOK) then
         call message( routine_name,
     &  ' [OPERATOR: Disregard previous MAPI message]', 0, 3 )
      endif
c     Second call will replace variables with correct inputs
      rtn = GMFIN(cube_handle,attrib,dtype,nelmnt,value)
      if (rtn .ne. 0) then
        call message( routine_name,
     &    'Error getting max EV frames from L1B file' //
     &    ' [OPERATOR ACTION: Stage correct L1B, rerun PGE]',
     &    0, 2 )
      end if

c     Save value into variable which will indicate the maximum number
c     of pixels in a scan for this granule.

      npixels = value

c -----------------------------------------------------------------

c ... Close the 1km L1B file

      rtn = -1
      rtn = clmfil(cube_handle)
      if( rtn .eq. 0 ) then
        call message( routine_name,
     &  'Success closing scan cube 1km data file ', 0, 3 )
      else
        call message( routine_name,
     &  'Error closing scan cube 1km data file.' //
     &  ' [OPERATOR ACTION: Check system resources]',
     &  0, 2 )
      endif

      end
