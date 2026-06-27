      subroutine file_close(h_eco1,h_eco2,no_250,modfil_Geo,
     +                      modfil_mod35,modfil_L1B_250,
     +                      modfil_L1B_1km,MET_Groups)

      implicit none
      save

      include 'mapi.inc'
      include 'PGS_SMF.f'

c     Common Block for debug statement
      common /bug/ debug, h_output

C-----------------------------------------------------------------------
C!F77
C
C!DESCRIPTION:
c Program which close all files used in the product science
c code by using PGS and mapi toolkit functions. Also attaches
c core, inventory and archive metadata to the file.
C
C!INPUT PARAMETERS:
C
C   integer h_*       File handles for generic files
C   integer modfil_*  HDF file id's
c   MET_Groups        Contains group names form metadata filing
c   no_250            Logical flag set to true if 250 data is used
C
C
C!OUTPUT PARAMETERS:  NONE
C
C!REVISION HISTORY:
c
C!TEAM-UNIQUE HEADER:
C
C!REFERENCES AND CREDITS
C
C!DESIGN NOTES:
C
C       The PCF file reference number must be accessed through
c       the product include file.
C
C EXTERNALS:
C
C       NAMED CONSTANTS:
C               PGS_S_SUCCESS                   (PGS_SMF.f)
C
C       FUNCTIONS AND SUBROUTINES:
c               pgs_io_gen_closef               (libPGSTK.a)
c               clmfil                          (mapi.a)
c               message.f
C
C!END-------------------------------------------------------------------

C     arguments
      character*(PGSd_MET_GROUP_NAME_L)
     *           MET_Groups(PGSd_MET_NUM_OF_GROUPS)
      integer modfil_mod35(MODFILLEN),modfil_Geo(MODFILLEN),
     *        modfil_L1B_250(MODFILLEN),modfil_L1B_1km(MODFILLEN),
     *        h_eco1,h_eco2
      logical no_250

c     local parameter values
c     Both are needed for correct metadata attachment
      integer ARCHIVEMETADATA, INVENTORYMETADATA
      parameter (INVENTORYMETADATA = 2,ARCHIVEMETADATA = 3)

c     Local scalars
      integer rtn,ier,debug,h_output,NUM_of_HDFAttrNms

c     Local arrays
      character*(MAX_ECS_NAME_L-1) HDFAttNms(PGSd_MET_NUM_OF_GROUPS)

C     Functions
      integer pgs_io_gen_closef

C ... First, close the generic files

c ... Close ecosystem files

c ... First the 1km NA map
      if (h_eco1 .ne. -5555) then
        ier = pgs_io_gen_closef(h_eco1)
        if( ier .eq. PGS_S_SUCCESS ) then
         call message( 'file_close',
     &    'Success closing 1km NA eco. file using pgs_io_gen_closef ',
     &    0, 3 )
        else
         call message( 'file_close',
     &    'Can-t close 1km NA eco. file using pgs_io_gen_closef'//
     &    ' [OPERATOR ACTION: Check system resources]',
     &    0, 2 )
        endif
      else

c ...   Next the global 10 minute map
        ier = pgs_io_gen_closef(h_eco2)
        if( ier .eq. PGS_S_SUCCESS ) then
         call message( 'file_close',
     &    'Success closing global eco. file using pgs_io_gen_closef ',
     &    0, 3 )
        else
         call message( 'file_close',
     &    'Can-t close global eco. file using pgs_io_gen_closef' //
     &    ' [OPERATOR ACTION: Check system resources]',
     &    0, 2 )
        endif
      endif

c ... Close runtime output debug file
      ier = pgs_io_gen_closef(h_output)
      if( ier .eq. PGS_S_SUCCESS ) then
        call message( 'file_close',
     &  'Success closing output debug file using pgs_io_gen_closef ',
     &  0, 3 )
      else
        call message( 'file_close',
     &  'Can-t close output debug file using pgs_io_gen_closef' //
     &  ' [OPERATOR ACTION: Check system resources]',
     &  0, 2 )
      endif


c ... Now for the hdf files

c ... Close the 250m L1B file - only close if file was opened correctly

      if (.not. no_250) then
        rtn = clmfil(modfil_L1B_250)
        if( rtn .eq. 0 ) then
          call message( 'file_close',
     &    'Success closing scan cube 250m data file ', 0, 3 )
        else
          call message( 'file_close',
     &    'Error closing scan cube 250m data file' //
     &    ' [OPERATOR ACTION: Check system resources]',
     &    0, 0 )
        endif
      endif

c ... Close the 1km L1B file

      rtn = clmfil(modfil_L1B_1km)
      if( rtn .eq. 0 ) then
        call message( 'file_close',
     &  'Success closing scan cube 1km data file ', 0, 3 )
      else
        call message( 'file_close',
     &  'Error closing scan cube 1km data file' //
     &  ' [OPERATOR ACTION: Check system resources]',
     &  0, 2 )
      endif

c ... Close the geolocation file

      rtn = clmfil(modfil_Geo)
      if( rtn .eq. 0 ) then
        call message( 'file_close',
     &  'Success closing geolocation data file ', 0, 3 )
      else
        call message( 'file_close',
     &  'Error closing geolocation data file' //
     &  ' [OPERATOR ACTION: Check system resources]',
     &  0, 2 )
      endif

c ... Close the output product file
c ... In this case, we need to use cpmfil to close the file in order
c ...  to attach the metadata fields to it.
      NUM_of_HDFAttrNms = 2
      HDFAttNms(INVENTORYMETADATA) = MECS_CORE
      HDFAttNms(ARCHIVEMETADATA) = MECS_ARCHIVE

      rtn = cpmfil(modfil_mod35,MET_Groups,HDFAttNms,NUM_of_HDFAttrNms)
      if( rtn .eq. 0 ) then
        call message( 'file_close',
     &  'Success closing output cloud mask file ', 0, 3 )
      else
        call message( 'file_close',
     &  'Error closing output cloud mask file' //
     &  ' [OPERATOR ACTION: Check system resources]',
     &  0, 2 )
      endif

      return
      end
