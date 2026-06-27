      Function Set_CoreMetadata_QC(LUN_Of_QCFile_MCF,LUN_Of_QCFile_RP)

      implicit none
      include 'mapi.inc'
      include 'PGS_SMF.f'
      include 'PGS_PC.f'
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !PURPOSE:
C   Write EOSDIS Core System (ECS) inventory metadata associated with
C   MODIS atmosphere "diagnostic" QC files.  These metadata are written 
C   as component fields to separate ascii files as provided by the SDPTK
C   for non-HDF files.  They are intended to satisfy ECS's requirements 
C   for search and identification of "diagnostic" files. 
C
C
C !DESCRIPTION:
C
C   The individual ECS metadata fields written by function
C   Set_CoreMetadata_QC, the source of the field values, and the
C   "Data Location" parameter of the Metadata Configuration File (MCF)
C   are listed below.  Data sources are identified as follows:
C
C   MCF      - Metadata Configuration File
C   TK       - set by SDP toolkit
C   PGE      - Passed in from calling routine
C   PCF      - Process Control File
C   internal - computed internally within function Set_CoreMetadata_QC
C   static   - predefined value set within function Set_CoreMetadata_QC
C
C                    ECS Inventory Metadata Fields
C                    -----------------------------
C
C                                          Source of    Data Location
C     ECS Core Metadata Objects              Value      Value in MCF
C   -------------------------------        ---------    -------------
C      1    ShortName                         MCF           MCF
C      2    VersionID                         MCF           MCF
C      3    ProductionDateTime                TK            TK
C      5    RangeBeginningTime                PCF           PGE
C      6    RangeEndingTime                   PCF           PGE
C      7    RangeBeginningDate                PCF           PGE
C      8    RangeEndingDate                   PCF           PGE
C
C
C
C
C !INPUT PARAMETERS:
C                     
C
C  integer  LUN_Of_QCFile_MCF    LUN of QC file MCF 
C  integer  LUN_OF_QCFile_RP     LUN of RP that references QC file
C
C
C !OUTPUT PARAMETERS: N/A
C                            
C
C !REVISION HISTORY:
c Revision 1.6  1998/11/20  20:18:35  lipo
c Update error message associated with pgs_met_write call.
c
c Revision 1.5  1998/06/18  15:55:49  fhliang
c Update error messages with "Operator Action" strings.
c
c Revision 1.4  1997/12/17  18:59:38  rhucek
c Set error_flag=.true. if return from pgs_met_write not SUCCESS
c
c Revision 1.3  1997/12/17  12:13:35  rhucek
c Updated code comments.
c
c Revision 1.2  1997/12/11  01:32:53  rhucek
c removed unreferenced variable "i"
c
c Revision 1.1  1997/12/11  01:18:45  rhucek
c Initial revision
c
C
C
C !TEAM-UNIQUE HEADER:
C
C    This software is developed by the MODIS Science Data Support Team
C    for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C !REFERENCES AND CREDITS:
C
C    Written by         Liqun Ma   November 1997
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd,  Seabrook MD 20706
C    lma@ltpmail.gsfc.nasa.gov
C
C !DESIGN NOTES:
C  
C  Returns:     0 if successful, -1 if an error occurs
C
C  Externals:
C
C    Functions and Subroutines:
C        PGS_MET_Init              (libPGSTK.a)
C        PGS_MET_SetAttr_*         (libPGSTK.a)
C        PGS_MET_Write             (libPGSTK.a)
C        PGS_PC_GetConfigData      (libPGSTK.a)
C        MODIS_SMF_SetDynamicMsg   (science code)
C        String_Loc                (science code)
C
C
C    Named Constant:
C        MCORE_*               (mapi.inc)
C        MECS_CORE             (mapi.inc)
C        MODIS_E_GENERIC       (PGS_MODIS_39500.f)
C        PGSd_PC_*             (PGS_PC.f)
C        PGSd_MET_*            (PGS_MET.f: included in mapi.inc)
C        PGS_S_SUCCESS         (PGS_SMF.f)
C    
C !END
C----------------------------------------------------------------------

C PARAMETER declarations
      character*(*) FUNCNAME
      PARAMETER (FUNCNAME = 'Set_CoreMetadata_QC')

      character*1 BLANK
      parameter (BLANK = ' ')

      integer FAIL, SUCCEED
      parameter (FAIL = -1, SUCCEED = 0)
      
      integer LUN_COLLECTION_START
      integer LUN_COLLECTION_STOP
      parameter (LUN_COLLECTION_START=10258)
      parameter (LUN_COLLECTION_STOP=10259)

      integer ODL_IN_MEMORY 
      parameter ( ODL_IN_MEMORY = 1 )

      integer INVENTORYMETADATA 
      parameter ( INVENTORYMETADATA = 2 )


C Declaration of function arguments
      integer LUN_Of_QCFile_MCF, LUN_Of_QCFile_RP


C Local Variable Declarations
      character*25  msg25
      character*(PGSd_PC_VALUE_LENGTH_MAX) Col_StartDateTime
      character*(PGSd_PC_VALUE_LENGTH_MAX) Col_StopDateTime
      character*128 AttrN,AttrV_s
      character*512 msgbuf
      character*(PGSd_MET_GROUP_NAME_L) handles(PGSd_MET_NUM_OF_GROUPS)

      integer    rtn,iparm

      integer    pgs_met_init,pgs_met_setattr_s,
     1           PGS_PC_GetConfigData,Set_CoreMetadata_QC,
     4           pgs_met_write,String_Loc

      integer    fbyte,fbyte2,lbyte,lbyte2

      logical    error_flag
      
C------------------------
C Initialization
C------------------------

      Set_CoreMetadata_QC = FAIL
      error_flag = .false.


C-----------------------------------------------------------------------
C Perform input argument checks and return if not valid
C-----------------------------------------------------------------------

c-----Check for valid input file id 


c-----Check for valid LUN of QC file MCF  
      If (LUN_Of_QCFile_MCF .lt. 1) Then
         error_flag = .true.

         write(msg25,'(I25)') LUN_Of_QCFile_MCF 
         rtn = String_Loc(msg25,fbyte,lbyte)

         msgbuf = 'PCF logical unit number of QC file MCF "'//msg25(fbyte:lbyte)
     2   // '" is out of bounds. '
     3   // char(10) //'PCF LUN numbers must be greater than zero.' 
     4   // char(10) // 'Operator Action:  Notify SDST.'

         call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf

c-----Check for valid LUN of QC file RP
      If (LUN_Of_QCFile_RP .lt. 1) Then
         error_flag = .true.

         write(msg25,'(I25)') LUN_Of_QCFile_RP 
         rtn = String_Loc(msg25,fbyte,lbyte)

         msgbuf = 'PCF logical unit number of QC file RP reference "'
     2   // msg25(fbyte:lbyte) // '" is out of bounds. '
     3   // char(10) //'PCF LUN numbers must be greater than zero.' 
     4   // char(10) // 'Operator Action:  Notify SDST.'

         call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf

c-----Return if invalid inputs found
      If (error_flag) Return


C-------------------------------------------------------------------
C Initialize metadata tool defining array "Handles", and set 
C metadata objects whose values are hard coded in the MCF.
C-------------------------------------------------------------------

      rtn = pgs_met_init(LUN_Of_QCFile_MCF,Handles)

      If (rtn.ne.PGS_S_SUCCESS) Then
         write(msg25,'(I25)') LUN_Of_QCFile_MCF
         rtn = String_Loc(msg25,fbyte,lbyte)

         msgbuf = 'pgs_met_init() unable to initialize MCF on LUN '
     2   // msg25(fbyte:lbyte) // '.'
     3   // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     4   // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     5   // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

         call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

         Return
      End If


C-----------------------------------------------------------------------
C Read the collection start "datetime" from PCF and set 
C
C "RANGEBEGINNINGDATE"
C "RANGEBEGINNINGTIME"
C-----------------------------------------------------------------------
 
c --- Read collection start datetime from PCF

      write(msg25,'(I25)') LUN_COLLECTION_START 
      rtn = String_Loc(msg25,fbyte,lbyte)

      rtn = PGS_PC_GetConfigData(LUN_COLLECTION_START,Col_StartDateTime)

      If (rtn .ne. PGS_S_SUCCESS) Then
         error_flag = .true.

         msgbuf = 'PGS_PC_GetConfigData() unable to read the collection '
     2   // char(10) // 'start datetime on LUN '//msg25(fbyte:lbyte) // '.'
     3   // char(10) // 'Operator Action:  Check for valid PCF Runtime Parameter (RP) '
     4   // char(10) // 'entry.  If LUN is nonexistent or RP value is blank, stage '
     5   // char(10) // 'correct PCF and rerun PGE.  Otherwise, notify SDST.'

         call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)


c --- check for blank string
      Else If (Col_StartDateTime .eq. BLANK) Then
         error_flag = .true.
         
         msgbuf = 'Collection start datetime on RP LUN ' 
     2   // msg25(fbyte:lbyte) // ' is blank.' 
     3   // char(10) // 'PGS_PC_GetConfigData() unable to set attributes '
     4   // char(10) // 'RangeBeginningDate and RangeBeginningTime.' 
     5   // char(10) // 'Operator Action:  Check for valid PCF Runtime Parameter '
     6   // char(10) // 'value.  If blank, stage correct PCF and rerun PGE. '
     7   // char(10) // 'Otherwise, notify SDST.'

         call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)


c --- Collection start datetime non blank.  Set Range beginning date and time 
      Else
         rtn = String_Loc(Col_StartDateTime,fbyte2,lbyte2)
          
c --- iparm=1, set start "date";  iparm=2, set start "time"
         Do 200 iparm  = 1,2
     
            If (iparm .eq. 1) Then
               AttrN   = MCORE_RANGE_BEG_DATE
               AttrV_s = Col_StartDateTime(fbyte2:fbyte2+9)
            Else
               AttrN   = MCORE_RANGE_BEG_TIME
               AttrV_s = Col_StartDateTime(fbyte2+11:lbyte2)
            End If

            rtn = PGS_MET_SetAttr_s(handles(INVENTORYMETADATA),
     1            AttrN,AttrV_s)

            If (rtn .ne. PGS_S_SUCCESS) Then
               error_flag = .true.
               rtn = String_Loc(AttrN,fbyte,lbyte)

               msgbuf = 'PGS_MET_SetAttr_s() unable to set ' // AttrN(fbyte:lbyte) // '.'
     2         // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     3         // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     4         // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

               call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
            End If
  200    Continue        

      End If

C-----------------------------------------------------------------------
C Read the collection stop "datetime" from PCF and set 
C
C "RANGEENDINGDATE"
C "RANGEENDINGTIME"
C-----------------------------------------------------------------------

c --- Read collection stop datetime from PCF

      write(msg25,'(I25)') LUN_COLLECTION_STOP
      rtn = String_Loc(msg25,fbyte,lbyte)

      rtn = PGS_PC_GetConfigData(LUN_COLLECTION_STOP,Col_StopDateTime)

      If (rtn .ne. PGS_S_SUCCESS) Then
         error_flag = .true.

         msgbuf = 'PGS_PC_GetConfigData() unable to read the collection '
     2   //char(10) //'stop datetime on LUN '//msg25(fbyte:lbyte) // '.'
     3   // char(10) // 'Operator Action:  Check for valid PCF Runtime Parameter (RP) '
     4   // char(10) // 'entry.  If LUN is nonexistent or RP value is blank, stage '
     5   // char(10) // 'correct PCF and rerun PGE.  Otherwise, notify SDST.'

         call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)


c --- check for blank string
      Else If (Col_StopDateTime .eq. BLANK) Then
         error_flag = .true.
         
         msgbuf = 'Collection stop datetime on RP LUN ' 
     2   // msg25(fbyte:lbyte) // ' is blank.' 
     3   // char(10) // 'PGS_PC_GetConfigData() unable to set attributes '
     4   // char(10) // 'RangeEndingDate and RangeEndingTime.' 
     5   // char(10) // 'Operator Action:  Check for valid PCF Runtime Parameter '
     6   // char(10) // 'value.  If blank, stage correct PCF and rerun PGE. '
     7   // char(10) // 'Otherwise, notify SDST.'

          call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)


c --- Collection stop datetime non blank.  Set Range ending date and time 
      Else
         rtn = String_Loc(Col_StopDateTime,fbyte2,lbyte2)
          
c --- iparm=1, set stop "date";  iparm=2, set stop "time"
         Do 300 iparm  = 1,2
     
            If (iparm .eq. 1) Then
               AttrN   = MCORE_RANGE_ENDING_DATE
               AttrV_s = Col_StopDateTime(fbyte2:fbyte2+9)
            Else
               AttrN   = MCORE_RANGE_ENDING_TIME 
               AttrV_s = Col_StopDateTime(fbyte2+11:lbyte2)
            End If

            rtn = PGS_MET_SetAttr_s(handles(INVENTORYMETADATA),
     1            AttrN,AttrV_s)

            If (rtn .ne. PGS_S_SUCCESS) Then
               error_flag = .true.
               rtn = String_Loc(AttrN,fbyte,lbyte)

               msgbuf = 'PGS_MET_SetAttr_s() unable to set ' // AttrN(fbyte:lbyte) // '.'
     2         // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     3         // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     4         // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'

               call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
            End If
  300    Continue        

      End If



C----------------------------------------------------------------------
C Write QC file ECS metadata to separate ascii file. 
C-----------------------------------------------------------------------


C Write inventory Metadata to QC file.

      rtn = pgs_met_write(handles(ODL_IN_MEMORY),' ',LUN_Of_QCFile_RP)

      if (rtn.ne.PGS_S_SUCCESS) Then
         error_flag = .true.

         write(msg25,'(I25)') LUN_Of_QCFile_RP 
         rtn = String_Loc(msg25,fbyte,lbyte)

         msgbuf = 'pgs_met_write detected error setting ascii metadata to non-HDF file ' 
     1   // char(10) // 'referenced on Runtime Parameter LUN = ' // msg25(fbyte:lbyte) // '.'
     2   // char(10) // 'Operator Action: Check for wrong or corrupted PCF or input'
     3   // char(10) // 'file. If a fault is identified, stage correct PCF or input file'
     4   // char(10) // 'and rerun PGE. Otherwise, notify SDST.'

         call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

      end if

      If (.not.error_flag) Set_CoreMetadata_QC = SUCCEED 

      Return

      End
