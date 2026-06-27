      subroutine Set_Runtime_Meta

      implicit none
      save

      include 'PGS_SMF.f'
      include 'mapi.inc'
      include 'mod35.inc'
      include 'mod35_meta.inc'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C
c    This simple program initializes and writes the
c    required core metadata into the runtime output
c    debug file.  No parameters are passed because
c    all it requires are the PCF LUN's of the
c    MCF_QA file and ASCII MCF file (which points to the
c    runtime output debuf file).  These numbers are
c    found in the product include file.
c
c!Input Parameters:
c None.
c
c !Output Parameters:
c None
c
C !REVISION HISTORY:
c 04/14/1999 fhliang fixed prolog.
C
C !Team-Unique Header:
c This software is developed by the MODIS Science Data Support
c Team for the National Aeronautics and Space Administration,
c Goddard Space Flight Center, under contract NAS5-32373.
C
C !REFERENCES AND CREDITS:
C
c    Written by Kathy Strabala with an outline from
c    Rich Hucek.
C
C !DESIGN NOTES:
C    Externals:
C        PGS_MET_Init    (libPGSTK.a)
C        PGS_MET_Write   (libPGSTK.a)
c        message         (atmos shared code)
C
C !END
C---------------------------------------------------------------------

c ... Runtime output number of groups
      character*(PGSd_MET_GROUP_NAME_L)
     *           groups(PGSd_MET_NUM_OF_GROUPS)

c fhliang 01/08/98: commented 2 lines
c ... Added for runtime output file
c     integer ODL_IN_MEMORY
c     parameter (ODL_IN_MEMORY = 1)

c gfireman 12/27/04: declared FileNameSuffix 
      character*3 FileNameSuffix

c fhliang 01/08/98: added declaration of Set_CoreMetadata_QC
c gfireman 12/27/04: changed to Write_ECSMetadata_nonHDF
      integer i,rtn, Write_ECSMetadata_nonHDF

c ... external subroutines
      external message, Write_ECSMetadata_nonHDF

c fhliang 01/08/98: commented 2 lines
c ... external functions
c     integer pgs_met_init,pgs_met_write
c     external pgs_met_init,pgs_met_write

c fhliang 01/08/98: Commented 12 lines and added call to function
c                   Set_CoreMetadata_QC, and changed status message.
cc ... Initialization of the output holder of the ODL metadata
c      do i = 1, PGSd_MET_NUM_OF_GROUPS
c         groups(i) = ' '
c      end do
c
cc ... initialize mcf for runtime output file
c      rtn = pgs_met_init(LRN_MCF_QC,groups)
c      if (rtn .ne. PGS_S_SUCCESS) then
c       call message( 'Set_Runtime_Meta',
c     &  'Error initializing runtime output metadata ', 0, 3 )
c      endif
c
cc ... Now write the metadata to an ascii file using PGS tool
c      rtn = pgs_met_write(groups(ODL_IN_MEMORY),' ',LRN_MCF_QC_met)
c      if (rtn .ne. PGS_S_SUCCESS) then
c       call message( 'Set_Runtime_Meta',
c     &  'Error writing runtime metadata to ascii file ', 0, 3 )
c      endif

c gfireman 12/28/04: changed to call Write_ECSMetadata_nonHDF
      FileNameSuffix = 'txt'
      rtn = Write_ECSMetadata_nonHDF( LRN_MCF_QC,
     1                                LRN_Of_QCFile_RP,
     2                                LUN_OF_NUM_INV_RP_PAIRS,
     3                                0,
     4                                FileNameSuffix )

      if (rtn .ne. PGS_S_SUCCESS) then
         call message( 'Set_Runtime_Meta',
     &   'Non-zero return detected from Write_ECSMetadata_nonHDF' //
     &   '[OPERATOR ACTION: Notify SDST]',
     &    0, 1 )
      else
         call message( 'Set_Runtime_Meta',
     2   'Call to Write_ECSMetadata_nonHDF returned successfully.', 0, 3 )
      end if


      return
      end
