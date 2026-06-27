      SUBROUTINE MOD_PR06CD_V2(modis_flag,idebug,Modfil_mod06cd,
     *           MinSolarZenithAngle,SLOPE_MEAN_LAND,SLOPE_MEAN_OCEAN)
C-----------------------------------------------------------------------
C !F77
C
C !Description: This subroutine is an interface for deriving cirrus
C               path radiance.
C
C              (This subroutine is much simpler compared to that in
C               deriving MOD06CD products. This has excludedthe inputs
C               time check, day/nigth flag, and metadata parameters
C               since they are done in Driver_Aerosol_WaterVapor.f)
C
C
C !Input Parameters:
C
C    modis_flag             Flag for using MODIS data
C    idebug                 Debugging index
C    Modfil_mod06cd         ID for opening MOD06CD HDF file (not used)
C    MinSolarZenithAngle    Minimum solar zenith angle within a granule
C
C !Output Parameters: none
C
C !Revision History:
C  Modified by D. A. Chu 6/11/99
C  Implementation of cirrus correction in PGE04 for aerosol retrieval
C
C  Modified by Liqun Ma   03/11/98
C  Updated the logical and some error msgs.
C
C  Modified by Liqun Ma   02/18/98
C  Some unreferenced variables and unused include files are moved out.
C  Some temporary output were moved out
C  Prolog was updated
C  Replaced calling Set_CoreMetadata and Set_ArchiveMetadata
C  to calling Update_InvMet_MOD06 and Update_ArchMet_MOD06
C  Added input file check section.
C
C !TEAM-UNIQUE HEADER:
C
C   This software was developed by the MODIS Science Team
C   for the National Aeronautics and Space Administration,
C   Goddard Space Flight Center.
C
C !REFERENCES and CREDITS:
C
C   Written by
C   Dr. Wei Han                                    11/24/97
C   SFA, Inc.
C   Naval Research Laboratory
C   Code 7212
C   Washington, DC 20375
C
C !DESIGN NOTES:
C
C   At present, this program is set up to process one granule of MODIS
C   data.
C
C   Functions and Subroutines:
c
C       Process_Mod06CD
C       MODIS_SMF_SETDYNAMICMSG
C
C !END
C-----------------------------------------------------------------------
C
      IMPLICIT NONE

      INCLUDE 'mapi.inc'
      INCLUDE 'PGS_MODIS_39500.f'

      CHARACTER*(*) FUNCNAME
      PARAMETER    (FUNCNAME = 'MOD_PR06CD_V2')
      CHARACTER*512 msgbuf
      REAL SLOPE_MEAN_LAND(3),SLOPE_MEAN_OCEAN(3)
      REAL MinSolarZenithAngle,SolarZenithAngleZEPS
      PARAMETER(SolarZenithAngleZEPS=72.000001)
      INTEGER Modfil_mod06cd(MODFILLEN),idebug
      LOGICAL modis_flag,error_flag,error_flag2

C
C Print message to LogStatus
C

      msgbuf =
     1   '----------------------------------------'
     2// char(10)
     3// 'Begin Cirrus cloud Detection Algorithm '
     4// char(10)
     5// '----------------------------------------'

      Call MODIS_SMF_SETDYNAMICMSG(MODIS_A_GENERIC,msgbuf,FUNCNAME)

C
C Initialization of error flags
C
      error_flag = .FALSE.
      error_flag2 = .FALSE.

C
C Cirrus path radiance is only derived when MinSolarZenithAngle < 72 degrees
C

      IF(modis_flag) THEN

        IF(MinSolarZenithAngle.LT.SolarZenithAngleZEPS) THEN

          CALL Process_Mod06CD(error_flag2,idebug,SLOPE_MEAN_LAND,SLOPE_MEAN_OCEAN)

          IF(error_flag2)  Then     !never happens; save for later use
            error_flag = .true.

            msgbuf =
     1      'Process_Mod06CD failed.'
     2      // char(10) // 'Operator Action: None'

            Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
          ENDIF

        END IF

        If (error_flag) Then
          msgbuf = 'MODIS process MOD_PR06CD Failed, exiting code 1'
          Call MODIS_SMF_SETDYNAMICMSG(MODIS_A_GENERIC,msgbuf,FUNCNAME)
C          Call exit(1)
        Else
          msgbuf = 'MODIS process MOD_PR06CD completed successfully,'
     2          // ' exiting code 0.'
          Call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,msgbuf,FUNCNAME)
C          Call exit(0)
        End If

      END IF

      RETURN
      END
