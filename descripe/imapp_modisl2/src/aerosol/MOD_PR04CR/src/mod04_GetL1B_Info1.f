      subroutine mod04_GetL1B_Info( nscans, npixels )

      IMPLICIT NONE
      SAVE

      include 'mod04_pcfnum1.inc'
      include 'mapi.inc'
      include 'PGS_PC.f'
      include 'PGS_SMF.f'


C-----------------------------------------------------------------------
C !F77
C 
C !DESCRIPTION:  This routine retrieves the number of scans in the MODIS
C                L1B 1KM input product and the number of 1km pixels in a 
C                scan. 
C
C
C !INPUT PARAMETERS: NONE
C
C 
C !OUTPUT PARAMETERS:
C        INTEGER  nscans     number of scans for this granule
C        INTEGER  npixels    Maximum number of pixels per scan in
C                            this particular granule
C 
C
C !REVISION HISTORY:
C  11/19/98 Kathy Strabala
C  Added OPERATOR disregard messages following calls
C  to GMFIN. MAPI returns an error message when there
C  really is no problem.
C 
C  10/28/98 Kathy Strabala
C  Started with Get_Global_Metadata.f from MOD06_PR06UW.
C  Stripped out debug portion for use in just the create
C  module for MOD06.
C 
C !TEAM-UNIQUE HEADER:
C 
C    This software is developed by the MODIS Science Data Support Team
C    for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C    HDF portions developed at the National Center for Supercomputing
C    Applications at the University of Illinois at Urbana-Champaign.
C
C 
C !DESIGN NOTES:
C
C   Please note that calling this subroutine will cause the following
C   warning message to be inserted into the LogStatus file:  
C
C     'MTYPEf2c():MAPI_E_ERR:324431360 (PID=26923)
C      ERROR: MTYPEf2c unable to use empty input parameter.'
C
C   This occurs when empty inputs for variables 'dtype' and 'nelmnt' 
C   are passed to GMFIN(). This is an unfortunate message since the 
C   M-API library was designed to accept empty inputs, and to return 
C   their actual values.
C
C 
C    Externals:
C
C       Functions and Subroutines:
C          clmfil                    (mapi lib)
C          GMFIN                     (mapi lib)
C          message                   (src_UW)
C          mod04_GetL1B_Info         (science code)
C          opmfil                    (mapi lib)
C          pgs_pc_getreference       (PGSTK lib)
C
C
C       Named Constant:
C          MAPIOK                    (mapi.inc)
C          MODFILLEN                 (mapi.inc)
C          PGSd_PC_FILE_PATH_MAX     (PGS_PC.f)
C          PGS_S_SUCCESS             (PGS_SMF.f)
C 
C 
C !END-----------------------------------------------------------------

c ... Arguments
      integer nscans, npixels

c ... Local scalars
      character*72 dtype, attrib
      integer nelmnt, value, rtn, cube_handle(MODFILLEN)

      integer file_id, file_version
      character*(PGSd_PC_FILE_PATH_MAX) cubefile

c ... Set program name for error messaging
      character*32   ROUTINE_NAME 
      parameter    ( ROUTINE_NAME = 'mod04_GetL1B_Info' )

c ... External routines
c      integer  pgs_pc_getreference, opmfil, clmfil
      integer  pgs_pc_getreference
      external pgs_pc_getreference, opmfil, clmfil


c -----------------------------------------------------------------

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

      if( rtn .ne. PGS_S_SUCCESS) call message( ROUTINE_NAME,
     &  'Error getting scan cube filename from PCF.' //
     &  ' [OPERATOR ACTION: Check PCF entry]', 0, 2 )


c ... Open L1B file and get file handle

      if ( opmfil( cubefile, 'r', cube_handle ) .ne. MAPIOK )then
        call message( ROUTINE_NAME,
     &  'Error opening MODIS L1B 1km input file.' //
     &  ' [OPERATOR ACTION: Check that file exists]', 0, 2 )
      else
        call message( ROUTINE_NAME,
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

c     Make inital call to get correct input parameters
c     without hardwiring

      rtn = GMFIN(cube_handle,attrib,dtype,nelmnt,value)
      if (rtn .ne. MAPIOK) then
         call message( ROUTINE_NAME,
     &  ' [OPERATOR: Disregard previous MAPI message]', 0, 3 )
      endif

c     Second call will replace variables with correct inputs
      rtn = GMFIN(cube_handle,attrib,dtype,nelmnt,value)
      if (rtn .ne. 0) then
        call message( ROUTINE_NAME,
     &    'Error getting nscan info from L1B file' //
     &    ' [OPERATOR ACTION: Stage correct L1B, rerun PGE]',
     &    0, 2 )
      endif

      nscans = value


c -----------------------------------------------------------------
c     Second parameter is Max Earth View Frames

      dtype = '  '
      attrib = 'Max Earth View Frames'
      nelmnt = 0
      value = 0
      rtn = 0

c     Make inital call to get correct input parameters
c     without hardwiring
      rtn = GMFIN(cube_handle,attrib,dtype,nelmnt,value)

      if (rtn .ne. MAPIOK) then
         call message( ROUTINE_NAME,
     &  ' [OPERATOR: Disregard previous MAPI message]', 0, 3 )
      endif

c     Second call will replace variables with correct inputs
      rtn = GMFIN(cube_handle,attrib,dtype,nelmnt,value)

      if (rtn .ne. 0) then
        call message( ROUTINE_NAME,
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
        call message( ROUTINE_NAME,
     &  'Success closing scan cube 1km data file ', 0, 3 )
      else
        call message( ROUTINE_NAME,
     &  'Error closing scan cube 1km data file.' //
     &  ' [OPERATOR ACTION: Check system resources]',
     &  0, 2 )
      endif

      end
