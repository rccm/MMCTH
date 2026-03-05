   Subroutine MOD04_MetaDataUpdate()

       implicit none
!       include 'mapi.f90.inc'
       include 'PGS_SMF.f'
       include 'PGS_PC.f'
       include 'PGS_PC_9.f'
       include 'PGS_MET.f' 
!       include 'hdf.f90'
       include 'dffunc.f90'
       include 'PGS_MODIS_39500.f'



!-----------------------------------------------------------------------
! Parameter declarations and definitions 
!-----------------------------------------------------------------------

      character*(*), PARAMETER :: BLANK = ' '
      integer, PARAMETER       :: INVENTORYMETADATA = 2,          &
                                  ARCHIVEMETADATA   = 3

      integer, PARAMETER       :: LUN_MOD04     = 405000,         &
                                  LUN_MOD06     = 412500,         &
                                  LUN_MOD35     = 422500,         &
                                  LUN_GEO       = 600000,         &
                                  LUN_L1B1KM    = 700002,         &
                                  LUN_L1B1KM_RA    = 430001

      integer, PARAMETER       :: MAX_LEN_AUTOQAFLAG_MP     = 64,  &
                                  MAX_LEN_AUTOQAFLAGEXPL_MP = 255, &
                                  MAX_LEN_NAME_MP           = 40,  &
                                  MAX_NUM_MP                = 2,   &
                                  MAX_NAME_LEN_PSA = 30

      integer, PARAMETER       :: MAX_ECS_NAME_L = 256,            &
                                  MODFILLEN      =   5

      integer, PARAMETER       :: DFACC_WRITE    = 2
                          
      character*(*), PARAMETER :: MECS_CORE    = 'CoreMetadata.0', &
                                  MECS_ARCHIVE = 'ArchiveMetadata.0'
      

      character*(*), PARAMETER :: FUNCNAME = 'MOD04_MetaDataUpdate'
     
      integer, PARAMETER       :: FAIL    = -1,   &
                                  SUCCESS =  1,   &
                                  MAPIOK  =  1


!-----------------------------------------------------------------------
! Initializations 
!-----------------------------------------------------------------------

! Measured Parameters
      character*(MAX_LEN_NAME_MP)           :: Name_MP           = BLANK
      integer                               :: NUM_MP            = 0
      integer                               :: QAPctMiss_MP      = 0
      character*(MAX_LEN_AUTOQAFLAG_MP)     :: AutoQAFlag_MP     = BLANK
      character*(MAX_LEN_AUTOQAFLAGEXPL_MP) :: AutoQAFlagExpl_MP = BLANK

! PSAs 
! On 08/09/2012 Changed the NUM_PSA to 15 in order to add the doi metadata into the output.
      integer                               :: NUM_PSA           = 15
      character*(MAX_NAME_LEN_PSA)          :: Name_PSA          = BLANK
      real                                  :: Value_PSA         = 0.0

      character*(PGSd_PC_FILE_PATH_MAX)     :: FileName
      character*2048                        :: msgbuf 

! Archive PSA
      integer                               :: NUM_ARCHPSA_SC   = 0
      character*(100)                       :: Name_ARCHPSA_SC  = BLANK
      real                                  :: Value_ARCHPSA_SC = 0. 

! Archive IP
      integer                               :: NUM_ARCH_IP      = 1,       & 
                                               LUN_ARCH_IP      = LUN_GEO, &
                                               VRSN_ARCH_IP     = 1
            
      
      integer :: NUM_IP
      integer, dimension(:), allocatable :: LUN_IP, VRSN_IP

      integer :: rtn, FileVrsn, SDID
      integer :: modfil_04(MODFILLEN)
      integer :: Num_of_HDfAttrNms = 2
      integer :: opmfil, cpmfil, Set_InvMet_MOD04, pgs_pc_getreference,   & 
                 sfstart, sfend, pgs_met_write, Set_ArchMet_MOD04

      character*(PGSd_MET_GROUP_NAME_L) Met_Handles(PGSd_MET_NUM_OF_GROUPS)
      character*(MAX_ECS_NAME_L-1)      HDFAttNms(PGSd_MET_NUM_OF_GROUPS)


      logical :: ExtGeoPntr_Flag = .false.
      logical :: error_flag 

      integer :: SUCCEED = 0

      integer :: Set_PSA, pgs_pc_getconfigdata
      character*(4) :: pcf_satid
      character*(26) :: doi_char
      integer :: LUN_Sat_Instrument = 800510


      integer  num_args
      integer  FlagRA
      character FlagBuff*10
      integer  iargc
      
      
      num_args = iargc ( )
      
      if(num_args .eq. 1) then
         call getarg ( 1, FlagBuff )
         read (FlagBuff,*) FlagRA
      else
         !     This is the default value
         FlagRA = 0
      endif

!-----------------------------------------------------------------------
! Executable code begins here
!-----------------------------------------------------------------------
      FileVrsn = 1
      rtn      = pgs_pc_getreference(412100, FileVrsn, FileName)
      
      If (FileVrsn .eq. 0) Then           !Aqua?
         NUM_IP = 22 
         allocate (LUN_IP(NUM_IP), VRSN_IP(NUM_IP))

         if( FlagRA .eq. 1) then
            LUN_IP(:)  = (/ 422500, 600000, 430001,                                         & 
                 412000, 412001, 412002, 412003, 412004, 412006, 412007, &
                 412008, 412100, 412101, 412300, 412301, 412302, 412303, & 
                 412304, 412307, 412305, 412306, 412400/)
         else
            LUN_IP(:)  = (/ 422500, 600000, 700002,                                         & 
                 412000, 412001, 412002, 412003, 412004, 412006, 412007, &
                 412008, 412100, 412101, 412300, 412301, 412302, 412303, & 
                 412304, 412307, 412305, 412306, 412400/)
         endif

         VRSN_IP(:) = (/ 1,1,1,            &
                         1,1,1,1,1,1,1,    &
                         1,1,1,1,1,1,1,    &
                         1,1,1,1,1  /)

      Else 
         NUM_IP = 36
         allocate (LUN_IP(NUM_IP), VRSN_IP(NUM_IP))

         if( FlagRA .eq. 1) then
            LUN_IP(:)  = (/ 422500, 600000, 430001,                                         & 
                 412000, 412001, 412002, 412003, 412004, 412005, 412006, 412007, & 
                 412008, 412009, 412010, 412011, 412012, 412013, 412014,         &
                 412100, 412100, 412101, 412101, 412102, 412102,                 &
                 412104, 412104, 412105, 412105, 412106, 412106,                 &
                 412108, 412108, 412109, 412109, 412110, 412110                  /)
         else
            LUN_IP(:)  = (/ 422500, 600000, 700002,                                         & 
                 412000, 412001, 412002, 412003, 412004, 412005, 412006, 412007, & 
                 412008, 412009, 412010, 412011, 412012, 412013, 412014,         &
                 412100, 412100, 412101, 412101, 412102, 412102,                 &
                 412104, 412104, 412105, 412105, 412106, 412106,                 &
                 412108, 412108, 412109, 412109, 412110, 412110                  /)
         endif

         VRSN_IP(:) = (/ 1,1,1,            &
                         1,1,1,1,1,1,1,1,  &
                         1,1,1,1,1,1,1,    &
                         1,2,1,2,1,2,      &
                         1,2,1,2,1,2,      &
                         1,2,1,2,1,2      /)
      End If

      
      rtn = Set_InvMet_MOD04( ExtGeoPntr_Flag,                                    &
                              NUM_IP,LUN_IP,VRSN_IP,                              &
                              NUM_MP,Name_MP,QAPctMiss_MP,                        &
                              AutoQAFlag_MP,AutoQAFlagExpl_MP,                    &
                              NUM_PSA, Name_PSA, Value_PSA,                       &
                              MET_Handles )

      If (rtn .eq. FAIL) Then
          msgbuf = 'Set InvMet_MOD04 failed. ' // char(10) // 'Operator Action:  Notify SDST' 
          Call MODIS_SMF_SetDynamicMsg( MODIS_E_GENERIC, msgbuf, FUNCNAME )
      Else
          msgbuf = 'Set InvMet_MOD04 succeeded. ' 
          Call MODIS_SMF_SetDynamicMsg( MODIS_S_GENERIC, msgbuf, FUNCNAME )
      EndIf


!-----------------------------------------------------------------------
!  Set the doi Attributes (PSAs).
!-----------------------------------------------------------------------
!           Get satellite instrument name.
      rtn = pgs_pc_getconfigdata(LUN_Sat_Instrument,pcf_satid)
      if (rtn .ne. 0) then
         call MODIS_SMF_SetDynamicMsg( MODIS_E_GENERIC, 'Metadata_MOD04_V2', &
         'Error reading instrument name from pcf LUN 800510 [OPERATOR ACTION: Check PCF file, fix if corrupt.]', FUNCNAME )
      endif

      if (pcf_satid .eq. 'AM1M') then
         doi_char = '10.5067/MODIS/MOD04_L2.006'
      else
         doi_char = '10.5067/MODIS/MYD04_L2.006'
      endif

      rtn = -1
      rtn = Set_PSA(MET_Handles, 'identifier_product_doi', NUM_PSA+1, doi_char)
      If (rtn .eq. FAIL) Then
         
         msgbuf = 'Set_PSA detected error setting DOI PSAs Operator Action:  ' // &
         'Refer to prior low level LogStatus error messages originating in function Set_InvPSA_Atmos'
         
         
         Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,'Metadata_MOD04_V2')
      EndIf

      rtn = -1
      rtn = Set_PSA(MET_Handles, 'identifier_product_doi_authority', NUM_PSA+2, 'http://dx.doi.org')
      If (rtn .eq. FAIL) Then
         
         msgbuf = 'Set_PSA detected error setting DOI PSAs Operator Action:'// &
         '  Refer to prior low level LogStatus error messages originating in function Set_InvPSA_Atmos.'
                  
         Call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,'Metadata_MOD04_V2')
      EndIf




      rtn = Set_ArchMet_MOD04( Met_Handles,                                       &
                               NUM_ARCH_IP,LUN_ARCH_IP,VRSN_ARCH_IP,              &
                               NUM_ARCHPSA_SC,                                    &
                               Name_ARCHPSA_SC,                                   &
                               Value_ARCHPSA_SC )

      If (rtn .eq. FAIL) Then
          msgbuf = 'Set ArchMet_MOD04 failed. ' // char(10) // 'Operator Action:  Notify SDST' 
          Call MODIS_SMF_SetDynamicMsg( MODIS_E_GENERIC, msgbuf, FUNCNAME )
      Else
          msgbuf = 'Set ArchMet_MOD04 succeeded. ' 
          Call MODIS_SMF_SetDynamicMsg( MODIS_S_GENERIC, msgbuf, FUNCNAME )
      EndIf



!---- Re-open HDF product file
      FileVrsn = 1
      rtn      = pgs_pc_getreference(LUN_MOD04,FileVrsn,FileName)
      
!     rtn      = opmfil(FileName,'a',modfil_04)
      SDID     = sfstart(FileName, DFACC_WRITE)
      
!.....write ECS metadata to aerosol product and close file
      Num_of_HDfAttrNms = 2
      HDFAttNms(INVENTORYMETADATA) = MECS_CORE
      HDFAttNms(ARCHIVEMETADATA)   = MECS_ARCHIVE

!...Write ECS inventory (Core) metadata to hdf product file
    rtn = pgs_met_write( MET_Handles(INVENTORYMETADATA), 'CoreMetadata.0', SDID)

    if (rtn .ne. 0) then
       msgbuf = 'pgs_met_write detected error writing CoreMetadata to Aerosol Product'   &
                // char(10) // 'HDF file.  Operator Action:  Notify SDST. '
       Call MODIS_SMF_SetDynamicMsg( MODIS_F_GENERIC, msgbuf, FUNCNAME )
    end if

!...Write ECS archive metadata to hdf product file
    rtn = pgs_met_write( MET_Handles(ARCHIVEMETADATA), 'ArchiveMetadata.0', SDID)

    if (rtn .ne. 0) then
       msgbuf = 'pgs_met_write detected error writing ArchiveMetadata to Aerosol Product'   &
                // char(10) // 'HDF file.  Operator Action:  Notify SDST. '
        Call MODIS_SMF_SetDynamicMsg( MODIS_F_GENERIC, msgbuf, FUNCNAME )
    end if

!...close product file
    rtn = sfend(SDID)

    if (rtn .ne. 0) then
       msgbuf = 'sfend detected error closing Aerosol Product HDF file. '   &
                // char(10) // 'Operator Action:  Notify SDST. '
       Call MODIS_SMF_SetDynamicMsg( MODIS_F_GENERIC, msgbuf, FUNCNAME )
    end if


END SUBROUTINE MOD04_MetaDataUpdate 
