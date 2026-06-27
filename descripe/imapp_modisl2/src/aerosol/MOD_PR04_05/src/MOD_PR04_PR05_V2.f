      PROGRAM Driver_Aerosol_WaterVapor

      IMPLICIT NONE

      INCLUDE 'mapi.inc'
      INCLUDE 'PGS_PC.f'
      INCLUDE 'PGS_SMF.f'
      INCLUDE 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !Description: This program is a driver MODIS MOD04 and MOD05 products
C
C
C !Input Parameters:
C      MODIS radiance data, solar and view zenith angles, cloud mask,
C      surface elevation, surface type (land or water).
C
C !Output Parameters:
C      MOD04 and MOD05 products (see file specs)
C
c !REVISION HISTORY:
c  01/29/98 fhliang
c  fixed prolog.
c
c !TEAM-UNIQUE HEADER:
c
c   Developed by MODIS Aerosol/Solar Water Vapor Retrieval Team
c   GSFC, Greenbelt, MD
c
c
C !REFERENCES and CREDITS:
C
C   Written by
C   Dr. Allen Chu                                    10/9/97
C   code 913
C   NASA Goddard Space Flight Center
C   Greenbelt, MD 20771
C
C !DESIGN NOTES:
C
C   At present, this program is set up to process one granule of MODIS
C   sensor data. A granule consists of 100 MODIS scan swathes, each
C   containing 1354 1-km pixels in the scan direction and 10 pixels
C   along the spacecraft flight direction.  This program also requires
C
C Externals:
C
C   Named Constants:
C
C Internals:
C
C   Functions and Subroutines:
C
C !END
C-----------------------------------------------------------------------

c ... lipo 30/10/98:  Added parameter declarations for temporary file LUNs
      integer    LUN_TEMP_GDAS,LUN_TEMP_OZN,LUN_TEMP_ICE
      parameter (LUN_TEMP_GDAS=497000, LUN_TEMP_OZN=497020, LUN_TEMP_ICE=497040)

      CHARACTER*(*) FUNCNAME
      PARAMETER    (FUNCNAME = 'Driver_Aerosol_WaterVapor')

      INTEGER    SUCCESS
      PARAMETER (SUCCESS = 0)

C miscellaneous local variables
      character*25 msg25
      character*1024 msgbuf

      INTEGER i,idebug,rtn,rtn_3,rtn_L1B,rtn_GEO,rtn_CldMsk,rtn_MOD07,
     1        rtn_met,rtn_ozn,rtn_sst,rtn_ice,RTN_NCEP
      integer    PGS_IO_Gen_Temp_Delete, rtn2,string_loc,fbyte,lbyte

C rhucek 11/22/02: added NumTempFiles and TempFile_LUNList declarations
      integer Delete_CommonTempFiles_Atmos, NumTempFiles,
     1        TempFile_LUNList(4) /4*0/

      LOGICAL error_flag, modis_flag

C For mod05 creation
      CHARACTER*512 usrlog
      INTEGER create_mod05,chk_input_L2,PGS_PC_GetConfigData

C For METADATA use
      INTEGER FILE_VERSION,No_PSA,MODFIL_MOD05(MODFILLEN),
     1        MODFIL_MOD06CD(MODFILLEN),
     2        MODFIL_MOD04(MODFILLEN),NumHandles

      PARAMETER(No_PSA=16,NumHandles=2)
      REAL QA_Metadata_MOD05(No_PSA)
      CHARACTER*(PGSd_MET_GROUP_NAME_L) groups(PGSd_MET_NUM_OF_GROUPS)
      CHARACTER*(MAX_ECS_NAME_L-1) HDFAttNms(PGSd_MET_NUM_OF_GROUPS)
      CHARACTER*(PGSd_PC_FILE_PATH_MAX) HDF_FILENAME
      CHARACTER*1024 MSG
      CHARACTER*1 work_buf
      INTEGER LRN_L1B,LRN_L1B_RA,LRN_GEO,LRN_CldMsk,LRN_MOD05,LRN_MOD07,LRN_DIAG,
     1        LRN_WISC_ANC_MET,LRN_WISC_ANC_OZN,LRN_WISC_ANC_SST,
     2        LRN_WISC_ANC_ICE
      INTEGER pgs_pc_getreference,pgs_met_getpcattr_s,Vrsn_No
      CHARACTER*8 InstrumentMode, char_buf 
      CHARACTER*256 HDFAttrName,ECSParmName
      REAL MinSolarZenithAngle,SolarZenithAngleZEPS
      REAL SLOPE_MEAN_LAND(3),SLOPE_MEAN_OCEAN(3)
      PARAMETER(SolarZenithAngleZEPS=84.000001)

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

c rhucek 01/19/98:  added status message
      usrlog = char(10) 
     1   // '-----------------------------------------------'
     2   // char(10)
     3   // 'Running Aerosol/Water Vapor Retrieval Algorithm'
     4   // char(10)
     5   // '-----------------------------------------------'

c      usrlog = 'Running MODIS Global Aerosol/Water Vapor Retrieval '
c     1         // 'Algorithm'

      Call MODIS_SMF_SETDYNAMICMSG(MODIS_M_GENERIC,usrlog,FUNCNAME)


C.....Initialization
      error_flag = .false.
      Vrsn_No = 1
      NumTempFiles = 0

      LRN_L1B     = 700002
      LRN_L1B_RA     = 430001
      LRN_GEO     = 600000
      LRN_CldMsk  = 422500
      LRN_MOD05   = 410000
      LRN_MOD07   = 420000
      LRN_DIAG    = 10911
      LRN_WISC_ANC_MET=900000
!      LRN_WISC_ANC_OZN=900020
      LRN_WISC_ANC_SST=900030
      LRN_WISC_ANC_ICE=900040


C
C Check for NCEP data (met,ozn,sst and ice) 8/2/99
C
      file_version = 1
      HDF_FILENAME = ' '
      rtn_met = -1
      rtn_met = pgs_pc_getreference(LRN_WISC_ANC_MET,file_version,
     &                             HDF_FILENAME)
      If ( rtn_met .eq. 0) Then
         NumTempFiles = NumTempFiles + 1
         TempFile_LUNList(NumTempFiles) = LUN_TEMP_GDAS 
      End If

!      file_version = 1
!      HDF_FILENAME = ' '
!      rtn_ozn = -1
!      rtn_ozn = pgs_pc_getreference(LRN_WISC_ANC_OZN,file_version,
!     &                             HDF_FILENAME)
!      If ( rtn_ozn .eq. 0) Then
!         NumTempFiles = NumTempFiles + 1
! !        TempFile_LUNList(NumTempFiles) = LUN_TEMP_OZN 
!      End If

C rhucek 11/22/02: NCEP sst and  sea ice concentration products no
C                  longer used by MOD_PR04_05 
C      file_version = 1
C      HDF_FILENAME = ' '
C      rtn_sst = -1
C      rtn_sst = pgs_pc_getreference(LRN_WISC_ANC_SST,file_version,
C     &                             HDF_FILENAME)
C      file_version = 1
C      HDF_FILENAME = ' '
C      rtn_ice = -1
C      rtn_ice = pgs_pc_getreference(LRN_WISC_ANC_ICE,file_version,
C     &                             HDF_FILENAME)

c      RTN_NCEP = rtn_met + rtn_ozn 
       RTN_NCEP = rtn_met 
C.....Check that temporal coverage of input ESDTs MOD03, MOD07_L2, and
C     MOD35_L2 overlap processing interval

      rtn_CldMsk=chk_input_L2(LRN_CldMsk,Vrsn_No)
      rtn_GEO=chk_input_L2(LRN_GEO,Vrsn_No)
      rtn_MOD07=chk_input_L2(LRN_MOD07,Vrsn_No)
C      rtn_MOD07=0

      rtn_3 = rtn_GEO + rtn_MOD07 + rtn_CldMsk

      If (rtn_3 .EQ. MAPIOK) Then
         usrlog = 'MOD03 and MODIS L2 inputs satisfy input file time check'
         Call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,usrlog,FUNCNAME)

C........Retrieve MinSolarZenithAngle from Cloud Mask product 
         HDFAttrName = 'CoreMetadata.0'
         ECSParmName = 'MinSolarZenithAngle'

         rtn=pgs_met_getpcattr_s(LRN_CldMsk,Vrsn_No,HDFAttrName,ECSParmName,char_buf)

         If (rtn .NE. PGS_S_SUCCESS) Then
            error_flag = .true.

            usrlog =
     1      'Halt process! pgs_met_getpcattr_s unable to retrieve MinSolarZenithAngle '
     2      // char(10) // 'from MOD35_L2 product.  MOD_PR04_05 exiting fail code 1. '
     3      // char(10) // 'Operator Action: Notify SDST.  Also check that MOD03 product ' 
     4      // char(10) // 'contains valid geolocation data including the solar zenith angle. '

            Call MODIS_SMF_SetDynamicMSG(MODIS_F_GENERIC,usrlog,FUNCNAME)
         End If

         read(char_buf,'(f8.2)') MinSolarZenithAngle


C........Retrieve DayNightFlag from Geolocation product
         HDFAttrName = 'CoreMetadata.0'
         ECSParmName = 'DAYNIGHTFLAG'

         rtn=pgs_met_getpcattr_s(LRN_GEO,Vrsn_No,HDFAttrName,ECSParmName,char_buf)

         If (rtn .NE. PGS_S_SUCCESS) Then
            error_flag = .true.

            usrlog =
     1      'Halt process! pgs_met_getpcattr_s unable to retrieve DayNightFlag '
     2      // char(10) // 'from Geolocation (MOD03) product.  MOD_PR04_05 exiting fail code 1. '
     3      // char(10) // 'Operator Action: Notify SDST.  Also check for possibly '
     4      // char(10) // 'corrupt Geolocation product.' 

            Call MODIS_SMF_SetDynamicMSG(MODIS_F_GENERIC,usrlog,FUNCNAME)
         End If

         InstrumentMode = char_buf
         modis_flag   = .true.

C........Read PCF for debug status (1:debugging; 0: no debugging)
         rtn=PGS_PC_GetConfigData(LRN_DIAG,work_buf)

         If (rtn .NE. PGS_S_SUCCESS) Then
            error_flag = .true.
            usrlog = 'PGS_PC_GetConfigData unable to retrieve Debug Mode Runtime Parameter'
            Call MODIS_SMF_SetDynamicMSG(MODIS_W_GENERIC,usrlog,FUNCNAME)
         Else
            READ(work_buf,'(I1)') idebug
         End If


C........Process as night - MOD05 IR only
         If (MinSolarZenithAngle .GE. SolarZenithAngleZEPS) Then

C...........Copy IR water vapor from MOD07_L2
            CALL MOD_PR05_V2(modis_flag,idebug,MODFIL_MOD05,MinSolarZenithAngle,RTN_NCEP)

            usrlog = 'Minimum Solar zenith angle is greater than 72 degrees. '
     1      // char(10) // 'MOD04 (Aerosol Product) not generated; only infrared parameters ' 
     2      // char(10) // 'included in MOD05 (Water Vapor Product).'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_M_GENERIC,usrlog,FUNCNAME)

C........Minumum solar zenith angle between 0 and 72 degrees
         Else If (MinSolarZenithAngle .GT. 0.0) Then

C...........instrument in night mode over illuminated earth - process as night
            If (InstrumentMode .EQ. 'Night') Then

C              Copy IR water vapor from MOD07_L2
                CALL MOD_PR05_V2(modis_flag,idebug,MODFIL_MOD05,MinSolarZenithAngle,RTN_NCEP)

               usrlog = 'The earth scene is illuminated, but MODIS instrument is in night mode. ' 
     1         // char(10) // 'MOD04 (Aerosol Product) not generated; only infrared parameters ' 
     2         // char(10) // 'included in MOD05 (Water Vapor Product).'

               Call MODIS_SMF_SETDYNAMICMSG(MODIS_M_GENERIC,usrlog,FUNCNAME)

C...........process as day
            Else If ( (InstrumentMode .EQ. 'Both') .OR. (InstrumentMode .EQ. 'Day') ) Then

C              Check that temporal coverage of MOD02 input overlaps processing interval.
               if( FlagRA .eq. 1) then
                  rtn_L1B=chk_input_L2(LRN_L1B_RA,Vrsn_No)
               else
                  rtn_L1B=chk_input_L2(LRN_L1B,Vrsn_No)
               endif

               If (rtn_L1B .ne. SUCCESS) Then
                  error_flag = .true.

                  usrlog = 'MODIS L1B 1KM product failed input file time check. '
     1            // char(10) // 'Process to halt.'
     2            // char(10) // 'Operator Action:  Stage MODIS L1B granule that matches '
     3            // char(10) // 'collection period, and rerun PGE.  If L1B granule and '
     4            // char(10) // 'collection period do match, contact SDST.'

                  Call MODIS_SMF_SetDynamicMSG(MODIS_E_GENERIC,usrlog,FUNCNAME)
               Else

C.................Perform Cirrus correction (add 7/20/99)
C                 CALL MOD_PR06CD_V2(modis_flag,idebug,MODFIL_MOD06CD,
C     1                    MinSolarZenithAngle,SLOPE_MEAN_LAND,SLOPE_MEAN_OCEAN)

C MOD05 processing commented out 18 June 2014 R. Cintineo
C.................Retrieve NIR water vapor and copy IR water vapor from MOD07_L2
C                   CALL MOD_PR05_V2(modis_flag,idebug,MODFIL_MOD05,
C     1                             MinSolarZenithAngle,RTN_NCEP)

C.................Process MOD04 and MOD05 correction
                  CALL MOD_PR04_V2(modis_flag,idebug,MODFIL_MOD04,RTN_NCEP,
     1                             SLOPE_MEAN_LAND,SLOPE_MEAN_OCEAN)

C MOD05 aerosol correction commented out 18 June 2014 R. Cintineo 
C.................Apply MOD05 aerosol correction; Correction to be applied post-launch
C                  CALL MOD_PR05_CORR_V2(modis_flag,idebug)
               End If   ! End L1B time check

C..............clean up sdptk temporary files, but only ones that exist 
               IF (RTN_NCEP.EQ.0) THEN
                  rtn = Delete_CommonTempFiles_Atmos(TempFile_LUNList, NumTempFiles)
               ENDIF

C...........Oh No, instrument mode is not applicable "NA".  Do not process granule.
            Else
               usrlog = 'Oh brother, instrument mode set to "NA".  PGE04 halting, ' 
     1         // char(10) // 'no retrieval performed.  Process exiting fail, code 1. ' 
     2         // char(10) // 'Operator Action: Notify SDST.  Also check for possibly '
     3         // char(10) // 'corrupt Geolocation product.' 

               Call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,usrlog,FUNCNAME)

            End If   ! End check on instrument mode 

C........Oh brother, minumum solar zenith angle is less than zero
         Else
            usrlog = 'Oh brother, minumum solar zenith angle is less than zero.  PGE04 '
     1      // char(10) // 'halting, no retrieval performed.  Process exiting fail, code 1. ' 
     2      // char(10) // 'Operator Action: Notify SDST.  Also check for possibly '
     3      // char(10) // 'corrupt Geolocation and MOD35_L2 product.' 

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,usrlog,FUNCNAME)
            
         End If   ! End check on solar zenith angle 

C.....MODIS input time checks failed
      Else
         error_flag = .true.

         IF(rtn_GEO.ne.MAPIOK) THEN

            usrlog = 'MODIS geolocation input failed input file time check. '
     1      // char(10) // 'Process to halt.'
     2      // char(10) // 'Operator Action:  Stage MODIS geolocation granule that matches '
     3      // char(10) // 'collection period, and rerun PGE.  If geolocation granule and '
     4      // char(10) // 'collection period do match, contact SDST.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,usrlog,FUNCNAME)

         ENDIF

         IF(rtn_CldMsk.ne.MAPIOK) THEN

            usrlog = 'MODIS cloud mask input failed input file time check. '
     1      // char(10) // 'Process to halt.'
     2      // char(10) // 'Operator Action:  Stage cloud mask granule that matches '
     3      // char(10) // 'collection period, and rerun PGE.  If cloud mask granule and '
     4      // char(10) // 'collection period do match, contact SDST.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,usrlog,FUNCNAME)

         ENDIF

         IF (rtn_MOD07 .ne. MAPIOK) THEN

            usrlog = 'MODIS MOD07 input failed input file time check. '
     1      // char(10) // 'Process to halt.'
     2      // char(10) // 'Operator Action:  Stage MOD07 granule that matches '
     3      // char(10) // 'collection period, and rerun PGE.  If MOD07 granule  and '
     4      // char(10) // 'collection period do match, contact SDST.'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,usrlog,FUNCNAME)

         ENDIF

      End If   ! End MOD03, MOD35_L2, and MOD07_L2 input file time checks

      If (error_flag) Then
         usrlog = 'MODIS process MOD_PR04_05 failed, exiting code 1'
     1   // char(10) // 'Operator Action: Action should be taken based upon '
     2   // char(10) // 'previous messages, such as due to solor zenith angle '
     3   // char(10) // 'not matched or inputs time interval not matched '

         Call MODIS_SMF_SetDynamicMSG(MODIS_E_GENERIC,usrlog,FUNCNAME)
         Call Exit(1)
      Else
         usrlog = 'MOD_PR04_05 completed successfully,'
     1   // char(10) // 'exiting code 0.'

         Call MODIS_SMF_SetDynamicMSG(MODIS_S_SUCCESS,usrlog,FUNCNAME)
         Call Exit(0)
      End If

      END
