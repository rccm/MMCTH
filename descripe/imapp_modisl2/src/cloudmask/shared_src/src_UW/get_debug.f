      SUBROUTINE GET_DEBUG(RP_DEBUG)
      
c!F77 ************************************************************
c ...
c!Description:
c     Routine which gets the runtime debug parameter out of the
c     pcf file.  
c     0 = no debug writes
c     1 = debug writes
c
c!Input parameters:
c None
c 
c!Output Parameters:
c rp_debug   Debug value - 0 or 1
c
c!Revision History:
c
c!Team-Unique Header:
c    Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c          external subroutines:  message.f
c          external functions:    pgs_pc_getconfigdata.f
c
c!Design Notes:
c
c!END****************************************************************

      implicit none

      save

      INCLUDE 'PGS_PC.f'

c     scalar argument
      integer rp_debug

c     character string version of debug value from pcf
      character*(PGSd_PC_VALUE_LENGTH_MAX) cbug

      integer LUN_debug_no,rtn
      parameter (LUN_debug_no = 10911)

c     toolkit runtime parameter extracter
      integer pgs_pc_getconfigdata

c     initialize
      rp_debug = 0
      cbug = ' '

c     Initial notice to operators 
        call message( 'get_debug',
     &  'OPERATOR NOTE:  If error message appears, please' //
     &  ' follow suggested action to solve problem and rerun PGE. '//
     &  ' If error remains, contact the SDST.',
     &  0, 3 )

      rtn = pgs_pc_getconfigdata(LUN_debug_no,cbug)
      if (rtn .ne. 0) then
        call message( 'get_debug',
     &  'Error reading debug value from pcf LUN 10911.' //
     &  ' [OPERATOR ACTION: Check PCF file, fix if corrupt.]',
     &  0, 1 )
      endif

c     convert to integer
      read(cbug,'(I5)') rp_debug

      END
