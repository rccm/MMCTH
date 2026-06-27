      subroutine mod06_get_rp(rp_debug,beg_lin,nlins,beg_ele,neles,
     +                        pcf_satid)
      
      implicit none
      save

c!F77 ************************************************************
c ...
c!Description:
c     Routine which gets the runtime debug parameter and user
c      defined processing intervals out of the pcf file.  
c      0 = no debug writes
c      1-4 = levels of debug output, from little to reams
c
c!Input parameters:
c None
c 
c!Output Parameters:
c rp_debug   Debug value - 0 or 1
c beg_lin    Beginning line number to process
c nlins      Number of lines to process
c beg_ele    Beginning element number to process
c neles      Number of elements to process
c pcf_satid  Platform code from pcf file
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

      INCLUDE 'PGS_PC.f'

c     scalar argument
      integer rp_debug,beg_lin,nlins,beg_ele,neles

c     character string version of debug value from pcf
      character*(PGSd_PC_VALUE_LENGTH_MAX) cbug, cblin,cnlins,
     &           cbele,cneles
      character*4 pcf_satid

      integer LUN_debug_no,LUN_blin_no,LUN_nlin_no,LUN_bele_no,
     &        LUN_nele_no,LUN_Sat_Instrument,rtn
      parameter (LUN_debug_no = 413150, LUN_blin_no = 413151,
     &           LUN_nlin_no = 413152, LUN_bele_no = 413153,
     &           LUN_nele_no = 413154, LUN_Sat_Instrument = 800510)
   
c     toolkit runtime parameter extracter
      integer pgs_pc_getconfigdata

c     initialize
      rp_debug = 0
      cbug = ' '
      beg_lin = 0
      cblin = ' '
      nlins = 0
      cnlins = ' '
      beg_ele = 0
      cbele = ' '
      neles = 0
      cneles = ' '
 
c     Initial notice to operators 
        call message( 'mod06_get_rp',
     &  'OPERATOR NOTE:  If error message appears, please' //
     &  ' follow suggested action to solve problem and rerun PGE. '//
     &  ' If error remains, contact the SDST.',
     &  0, 3 )

c     Get satellite instrument name.
      rtn = pgs_pc_getconfigdata(LUN_Sat_Instrument,pcf_satid)
      if (rtn .ne. 0) then
        call message( 'mod06_get_rp',
     &  'Error reading instrument name from pcf LUN 800510.' //
     &  ' [OPERATOR ACTION: Check PCF file, fix if corrupt.]',
     &  0, 1 )
      endif

      rtn = pgs_pc_getconfigdata(LUN_debug_no,cbug)
      if (rtn .ne. 0) then
        call message( 'mod06_get_rp',
     &  'Error reading debug value from pcf LUN 413150.' //
     &  ' [OPERATOR ACTION: Check PCF file, fix if corrupt.]',
     &  0, 1 )
      endif

c     convert bug level  to integer
      read(cbug,'(I5)') rp_debug
      if (rp_debug .lt. 0 .or. rp_debug .gt. 4) then
        call message( 'mod06_get_rp',
     &  'Debug runtime value is out of range [0-4], please ' //
     &  ' change in PCF file LUN 413150 and run again.]',
     &  0, 2 )
      endif



c ... Now retrieve beginning line number
      rtn = pgs_pc_getconfigdata(LUN_blin_no,cblin)
      if (rtn .ne. 0) then
        call message( 'mod06_get_rp',
     &  'Error reading beginning line value from pcf LUN 413151.' //
     &  ' [OPERATOR ACTION: Check PCF file, fix if corrupt.]',
     &  0, 1 )
      endif

c     convert beginning line number to integer
      read(cblin,'(I5)') beg_lin
      if (beg_lin .lt. 0 .or. beg_lin .gt. 2030) then
        call message( 'mod06_get_rp',
     &  'Begin line runtime value is out of range [0-2030], please ' //
     &  ' change in PCF file LUN 413151 and run again.]',
     &  0, 2 )
      endif


c ... Now retrieve number of lines to process
      rtn = pgs_pc_getconfigdata(LUN_nlin_no,cnlins)
      if (rtn .ne. 0) then
        call message( 'mod06_get_rp',
     &  'Error reading number of lines value from pcf LUN 413152.' //
     &  ' [OPERATOR ACTION: Check PCF file, fix if corrupt.]',
     &  0, 1 )
      endif

c     convert bug level  to integer
      read(cnlins,'(I5)') nlins
      if (nlins .lt. 0 .or. nlins .gt. 5500) then
        call message( 'mod06_get_rp',
     &  'Number of lines to process is out of range [0-5500],' //
     &  ' please change in PCF file LUN 413152 and run again.]',
     &  0, 2 )
      endif


c ... Now retrieve beginning element number
      rtn = pgs_pc_getconfigdata(LUN_bele_no,cbele)
      if (rtn .ne. 0) then
        call message( 'mod06_get_rp',
     &  'Error reading beginning element value from pcf LUN 413153.' //
     &  ' [OPERATOR ACTION: Check PCF file, fix if corrupt.]',
     &  0, 1 )
      endif

c     convert bug level  to integer
      read(cbele,'(I5)') beg_ele
      if (beg_ele .lt. 0 .or. beg_ele .gt. 1354) then
        call message( 'mod06_get_rp',
     &  'Beginning element number is out of range [0-1354], please ' //
     &  ' change in PCF file LUN 413153 and run again.]',
     &  0, 2 )
      endif


c ... Now retrieve number of elements to process
      rtn = pgs_pc_getconfigdata(LUN_nele_no,cneles)
      if (rtn .ne. 0) then
        call message( 'mod06_get_rp',
     &  'Error reading number of elements value from pcf LUN 413154.' //
     &  ' [OPERATOR ACTION: Check PCF file, fix if corrupt.]',
     &  0, 1 )
      endif

c     convert bug level  to integer
      read(cneles,'(I5)') neles
      if (neles .lt. 0 .or. neles .gt. 1354) then
        call message( 'mod06_get_rp',
     &  'Number of elements to process is out of range [0-1354],' //
     &  ' please change in PCF file LUN 413154 and run again.]',
     &  0, 2 )
      endif

      return
      end
