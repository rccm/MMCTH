        subroutine mod06_Get_Global_Metadata( nscans, npixels )

        IMPLICIT NONE
        SAVE

        include 'mod06uw_debug.inc'      

C!F77------------------------------------------------------------------
C
C!DESCRIPTION:  Extracts Global metadata from the L1B granule which
c               is needed for UW product generation.
c               Please note that calling this subroutine will cause
c               a warning message to be inserted into the LogStatus
c               file.  Please see description of message on line 73.
C
C!INPUT PARAMETERS:   
C
C       INTEGER  modfil(3)  Passes vital stats about input L1B granule
c
C!OUTPUT PARAMETERS: 
C       INTEGER  nscans     Number of scans for this granule
C       INTEGER  npixels    Maximum number of pixels per scan in 
c                           this particular granule
C
C!REVISION HISTORY:  
C
C!TEAM-UNIQUE HEADER: 
c HDF portions developed at the National Center for Supercomputing
c Applications at the University of Illinois at Urbana-Champaign.
C
C!REFERENCES and CREDITS:
C
C!DESIGN NOTES:
C
C    Externals:
C       Functions:
C          GMFIN            (mapi.inc)
c       Subroutines:        
c          message          Writes to LogStatus file
C!END-----------------------------------------------------------------

      include 'mod06uw_pcfnum.inc'
      include 'mapi.inc'
      include 'hdf.inc'
      include 'PGS_IO.f'
      include 'PGS_PC.f'
      include 'PGS_SMF.f'
      include 'PGS_MODIS_39500.f'

c ... Arguments

      integer nscans, npixels

c ... Local scalars

      character*72 dtype, attrib
      integer nelmnt, value, rtn, cube_handle(MODFILLEN)

      integer file_id, file_version
      character*(PGSd_PC_FILE_PATH_MAX) cubefile

c ... Set program name for error messaging

      character*32 routine_name
      parameter ( routine_name = 'mod06_Get_Global_Metadata' )

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
C     This is the default value
         FlagRA = 0
      endif

c -----------------------------------------------------------------

c ... get Level 1B scan cube file name from process control file
      if( FlagRA .eq. 1) then
         file_id = dscube_pcfnum_RA
      else
         file_id = dscube_pcfnum
      endif
      file_version = 1
      cubefile = ' '
      rtn = -1
      rtn = pgs_pc_getreference( file_id, file_version, cubefile )
      if( rtn .ne. PGS_S_SUCCESS) call message( routine_name,
     &  'Error getting scan cube filename from PCF.' //
     &  ' [OPERATOR ACTION: Check PCF entry]', 0, 2 )

c ... Open L1B file and get file handle

      if ( opmfil( cubefile, 'r', cube_handle ) .ne. MAPIOK )
     &  call message( routine_name,
     &  'Error opening MODIS L1B 1km input file.' //
     &  ' [OPERATOR ACTION: Check that file exists]', 0, 2 )

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

c ... Debug statement

      if (debug .gt. 0) then
        write( h_output, '(72(''-''))' )
        write( h_output, '(a,/)' ) 'mod06_Get_Global_Metadata debug'
        write(h_output,'(10x,''Number of scans in this granule: '', 
     +  i10)') nscans
        write(h_output,'(10x,'' Max number of elements in this granule:
     +  '',i10,/)') npixels
      endif

      end
