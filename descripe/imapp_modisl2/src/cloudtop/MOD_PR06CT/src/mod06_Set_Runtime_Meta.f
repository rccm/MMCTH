      subroutine mod06_Set_Runtime_Meta

      implicit none
      save

      include 'PGS_SMF.f'
      include 'PGS_PC.f'
      include 'PGS_MODIS_39500.f'
      include 'PGS_PC_9.f'
      include 'PGS_IO.f'
      include 'PGS_IO_1.f'
      include 'mapi.inc'
      include 'mod06uw_pcfnum.inc'


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
c!Input parameters: None
c 
c !Output Parameters: None
c 
C !REVISION HISTORY:
c $Id: mod06_Set_Runtime_Meta.f,v 1.5 1999/04/16 22:56:05 kis Exp $
C
C !Team Unique Header:
C
C !REFERENCES AND CREDITS:
C
c    Written by Kathy Strabala with an outline from
c    Rich Hucek.
C
C !DESIGN NOTES:
C    Externals:
C        PGS_MET_Init              (libPGSTK.a)
C        PGS_MET_Write             (libPGSTK.a)
c        message                    
C
C !END
C---------------------------------------------------------------------

c ... scalar argument
c     integer mcf_qcnum, bug_pcfnum

      integer i,rtn, Set_CoreMetadata_QC

c ... external subroutines
      external message

c 02/24/98 fhliang: commented following block and added call to function
c                   Set_CoreMetadata_QC
cc ... Runtime output number of groups
c      character*(PGSd_MET_GROUP_NAME_L) 
c     *           groups(PGSd_MET_NUM_OF_GROUPS)
c
cc ... Added for runtime output file
c      integer ODL_IN_MEMORY
c      parameter (ODL_IN_MEMORY = 1)
c
cc ... external functions
c      integer pgs_met_init,pgs_met_write
c      external pgs_met_init,pgs_met_write
c
cc ... Initialization of the output holder of the ODL metadata
c      do i = 1, PGSd_MET_NUM_OF_GROUPS
c         groups(i) = ' '
c      end do
c
cc ... initialize mcf for runtime output file 
c      rtn = pgs_met_init(mcf_qcnum,groups)
c      if (rtn .ne. PGS_S_SUCCESS) then
c       call message( 'mod06_Set_Runtime_Meta',
c     &  'Error initializing runtime output metadata ', 0, 3 )
c      endif
c
cc ... Now write the metadata to an ascii file using PGS tool
c      rtn = pgs_met_write(groups(ODL_IN_MEMORY),' ',bug_pcfnum)
c      if (rtn .ne. PGS_S_SUCCESS) then
c       call message( 'mod06_Set_Runtime_Meta',
c     &  'Error writing runtime metadata to ascii file ', 0, 3 )
c      endif

      rtn = Set_CoreMetadata_QC(mcf_qcnum, bugRPref_pcfnum)

      if (rtn .ne. PGS_S_SUCCESS) then
      call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,
     2   'Non-zero return from Set_CoreMetadata_QC detected',
     3   'mod06_Set_Runtime_Meta')
      else
         call modis_smf_setdynamicmsg(MODIS_S_SUCCESS,
     2   'Create metadata file to non-HDF (mod06CT_QC) without error',
     3   'mod06_Set_Runtime_Meta')
      end if

      return
      end
