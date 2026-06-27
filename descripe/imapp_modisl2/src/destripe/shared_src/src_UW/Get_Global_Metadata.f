      SUBROUTINE GET_GLOBAL_METADATA( MODFIL, NSCANS, MAX_PIXEL )

C-----------------------------------------------------------------------
C !F77
C
C !PURPOSE:
C    Extract number of scans and maximum pixel number from L1B granule.
C
C !DESCRIPTION:
C    Extracts Global metadata from the L1B granule which
C    is needed for UW product generation.
C    Please note that calling this subroutine will cause
C    a warning message to be inserted into the LogStatus
C    file.  Please see description of message MTYPEf2c().
C
C !INPUT PARAMETERS:
C    MODFIL       File information array containing
C                 HDF SD_ID (index 1)
C                 VS file identifier (index 2)
C                 File access mode (index 3)
C
C !OUTPUT PARAMETERS:
C    NSCANS       Number of earth scans in this granule
C    MAX_PIXEL    Maximum number of pixels per scan in this granule
C
C !REVISION HISTORY:
C
C !TEAM-UNIQUE HEADER:
C    Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
C
C !END
C-----------------------------------------------------------------------

      IMPLICIT NONE

      INCLUDE 'mapi.inc'
      INCLUDE 'PGS_SMF.f'

c ... Arguments

      integer modfil( MODFILLEN ), nscans, max_pixel

c ... Local variables

      character*72 dtype, attrib
      integer nelmnt, value, rtn

c ... Get number of scans for this granule

      dtype = '  '
      attrib = 'Number of Scans'
      nelmnt = 0
      value = 0
      rtn = 0

c ... The error message
c ... 'MTYPEf2c():MAPI_E_ERR:324431360 (PID=26923)
c ... ERROR: MTYPEf2c unable to use empty input parameter.'
c ... appears in LogStatus when empty inputs for variables 'dtype'
c ... and 'nelmnt' are passed to GMFIN().
c ... This is an unfortunate message since the M-API library
c ... was designed to accept empty inputs, and to return their actual
c ... values.

c ... Make inital call to get correct input parameters
c ... without hardwiring

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )
      
c ... Second call will replace variables with correct inputs

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )
      
      if (rtn .ne. 0) then
        call message( 'Get_Global_Metadata',
     &    'Error getting Number of Scans from L1B file' //
     &    ' [OPERATOR ACTION: Stage correct L1B, rerun PGE]', 0, 2 )
      endif

c ... Set return value

      nscans = value 

c ... Get number of earth scans for this granule

      dtype = '  '
      attrib = 'Max Earth View Frames'
      nelmnt = 0
      value = 0
      rtn = 0

c ... Make inital call to get correct input parameters 
c ... without hardwiring

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )

c ... Second call will replace variables with correct inputs

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )
      
      if (rtn .ne. 0) then
        call message( 'Get_Global_Metadata',
     &    'Error getting Max Earth View Frames from L1B file' //
     &    ' [OPERATOR ACTION: Stage correct L1B, rerun PGE]', 0, 2 )
      end if

c ... Set return value

      max_pixel = value 

      END
