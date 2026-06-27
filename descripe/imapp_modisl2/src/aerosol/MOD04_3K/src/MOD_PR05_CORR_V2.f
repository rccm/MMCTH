      SUBROUTINE MOD_PR05_CORR_V2(modis_flag,idebug)

C-----------------------------------------------------------------------
C !F77
C
C !Description: This program is to correct column water vapor values
C                above clear land surfaces and above clouds from MODIS
C                channels near 1 micron.
C
C !Input Parameters:
C
C      MODIS radiance data, solar and view zenith angles, cloud mask,
C      surface elevation, surface type (land or water), MOD04 aerosol
C      optical thickness and MOD05 column water vapor.
C
C !Output Parameters:
C
C      Correction factors to column water vapor values retrieved from MODIS
C      data based on a atmospheric transmittance model.
C
C !Revision History:
c 10/08/1999 fhliang
c resolved 'too many continuation lines' problem for Tab_value.
c
c 01/29/98 rhucek
c added explanation in Input /Output Parameters sections in prolog.
c
c 01/29/98 fhliang
c added NCSA acknowledgement.
c
C !TEAM-UNIQUE HEADER:
C
C   This software was developed by the MODIS Atmosphere Science Team
C   for the National Aeronautics and Space Administration at
C   Goddard Space Flight Center.
C
c HDF portions developed at the National Center for Supercomputing
c Applications at the University of Illinois at Urbana-Champaign.
c
C !REFERENCES and CREDITS:
C
C   Written by
C   Dr. Allen Chu                 10/17/97
C   Code 913
C   NASA Goddard Space Flight Center
C   Greenbelt, MD 20771
C   achu@climate.gsfc.nasa.gov
C
C !DESIGN NOTES:
C
C   At present, this program is set up to process one granule of MODIS
C   sensor data. A granule consists of 100 MODIS scan swathes, each
C   containing 1354 1-km pixels in the scan direction and 10 pixels
C   along the spacecraft flight direction.
C
C Externals:
C   Named Constants:
C           MAPIOK                         (mapi.inc)
C           Max_Frames_Per_Line            (MOD05_CORR.inc)
C           No_Lines_Per_Scan              (MOD05_CORR.inc)
C
C !END
C-----------------------------------------------------------------------

      IMPLICIT NONE
      SAVE

      INCLUDE 'MOD05_CORR.inc'
      INCLUDE 'mapi.inc'
      INCLUDE 'PGS_PC.f'
      INCLUDE 'PGS_MODIS_39500.f'

C rhucek 01/09/98:  Added parameter statement for FUNCNAME
      character*(*) FUNCNAME
      PARAMETER (FUNCNAME = 'MOD_PR05_CORR_V2')


c  Parameters for cloud array

      BYTE Cloud(Buf_cldmsk,Max_Frames_Per_Line, No_Lines_Per_Scan),
     2  QA_Cloud(Buf_cldmsk_QA,Max_Frames_Per_Line,No_Lines_Per_Scan)
      CHARACTER*5   ioflag
      CHARACTER*40  att_N,dtype
      CHARACTER*512 msgbuf
      CHARACTER*(PGSd_PC_FILE_PATH_MAX) FN_L1B
      CHARACTER*(PGSd_MET_GROUP_NAME_L) groups(PGSd_MET_NUM_OF_GROUPS)

      INTEGER Max_EV_Frames

      CHARACTER*4  Scan_Type(Max_Scans_Per_Granule)

      INTEGER i1,j1
      INTEGER i, idebug, j, Rank, RTN, nms, Scans_Per_Granule,
     &        Total_No_Lines, BeginScan_No, EndScan_No, Scan_No,
     &        No_of_Scans, handle_QC,
     &        Buf_Size1, Buf_Size2,
     &        Modfil_Geo(MODFILLEN),Modfil_CldMsk(MODFILLEN),
     &        Modfil_MOD05(MODFILLEN),Modfil_MOD04(MODFILLEN),
     &        Modfil_L1B(MODFILLEN),
     &        handle_corr(2), Dim_Size(2), Data_Size(2),
     &        Data_Size_Anc(2),
     &        No_Frames(Max_Scans_Per_Granule),bsize,
     &        No_Grids
      INTEGER GetModisDat_MOD05_CORR,Set_CoreMetadata_QC
      INTEGER CldMsk      (Max_Frames_Per_Line, No_Lines_Per_Scan),
     2        LandSea_Flag(Max_Frames_Per_Line, No_Lines_Per_Scan)

      INTEGER  NumHandles,ODL_IN_MEMORY,pgs_met_init,pgs_met_write
      PARAMETER(ODL_IN_MEMORY=1)

      LOGICAL modis_flag

      REAL SatZen  (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     SolZen  (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     RelAz   (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Refl_2  (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Refl_19 (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Ratio_Refl_19_2(Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Height  (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Vapor   (Max_Frames_Per_Line, No_Lines_Per_Scan)

      BYTE  UN_2    (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &      UN_19   (Max_Frames_Per_Line, No_Lines_Per_Scan)

C      REAL Aerosol_OD(Max_Grids_Per_Line,No_Grids_Per_Scan)

      BYTE Vflag_2(Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Vflag_19(Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Buf_Un(Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Buf_Sa(Max_Frames_Per_Line, No_Lines_Per_Scan)

c rhucek 01/09/98:  added status message
      msgbuf = char(10) 
     1// '-----------------------------------------------------------'
     2// char(10)
     3// 'Begin Aerosol Correction Algorithm to Water Vapor Retrieval'
     4// char(10)
     5// '-----------------------------------------------------------'

      Call MODIS_SMF_SETDYNAMICMSG(MODIS_M_GENERIC,msgbuf,FUNCNAME)


C-----------------------------------------------------------------------
C If modis_flag = .true., process MODIS Data; otherwise, open
C "non-MODIS" data files.  Only MODIS option is presently available.
C-----------------------------------------------------------------------

      IF (modis_flag) THEN

C Initialize buffers to hold HDF VS and SD file Ids

         DO 20 i = 1, MODFILLEN
            Modfil_L1B(i) = 0
            Modfil_Geo(i) = 0
            Modfil_CldMsk(i) = 0
            Modfil_MOD04(i) = 0
            Modfil_MOD05(i) = 0
20       CONTINUE

         ioflag = 'OPEN'

         CALL FILEOC_CORR(ioflag,modis_flag,handle_corr,
     &               handle_QC,FN_L1B,
     &               modfil_L1B,modfil_Geo,modfil_CldMsk,
     &               modfil_MOD04,modfil_MOD05)
C         WRITE(*,*)'Pass FILEOC'

         att_N = 'Number of Scans'
         dtype = 'INTEGER*4'
         nms = 1
         RTN = GMFIN(modfil_L1B,att_N,dtype,nms,Scans_Per_Granule)

         IF (rtn.ne.mapiok) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'MAPI function 1st GMFIN failed',
     &    'MOD_PR05_CORR_V2.f')
C         WRITE(*,*)'Pass Scans_Per_Granule',Scans_Per_Granule

         att_N = 'Max Earth View Frames'
         dtype = 'INTEGER*4'
         nms = 1
         RTN = GMFIN(modfil_L1B,att_N,dtype,nms,Max_EV_Frames)

         IF (rtn.ne.mapiok) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'MAPI function 2nd GMFIN failed',
     &    'MOD_PR05_CORR_V2.f')
C         WRITE(*,*)'Pass Max_EV_Frames',Max_EV_Frames

         bsize = Max_Scans_Per_Granule * 4
         RTN = GMTBL(modfil_L1B,'Level 1B Swath Metadata',' ',
     &         'EV_Frames',0,Scans_Per_Granule,bsize,No_Frames)

         IF (rtn.ne.mapiok) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'MAPI function 1st GMTBL failed',
     &   'MOD_PR05_CORR_V2.f')
C         WRITE(*,*)'Pass No_Frames'

         bsize = Max_Scans_Per_Granule*4
         rtn = GMTBL(modfil_L1B,'Level 1B Swath Metadata',' ',
     &         'Scan Type',0,Scans_Per_Granule,bsize,Scan_Type)

         IF (rtn.ne.mapiok) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'MAPI function 2nd GMTBL failed',
     &   'MOD_PR05_CORR_V2.f')
C         WRITE(*,*)'Pass Scan_Type',Scan_Type

         Total_No_Lines = Scans_Per_Granule * 10

         Buf_Size1 = Max_Frames_Per_Line
         Buf_Size2 = No_Lines_Per_Scan

         BeginScan_No = 1
         EndScan_No   = Scans_Per_Granule

C
C Read look-up tables for correction
C

c         open(7, file='Res2.out', status='unknown')
         CALL READ_REFL_0P86_TABLE_AEROSOL(handle_corr)

         CALL READ_RATIO_TABLE_AEROSOL(handle_corr)

C         WRITE(*,*)'Pass reading look-up table'


C Loop over the number of scan swathes.  For each swath, get MODIS
C inputs including sensor reflectances, geolocation, and cloud mask.

         DO 999 Scan_No = BeginScan_No, EndScan_No

C         WRITE(*,*) 'Scan_No',Scan_No
C
C "IF - ELSE - ENDIF " added by vlin   06/21/96
 
c 04/26/2000 fhliang: following 2 lines were added.
         Data_Size(1) = No_Frames(Scan_No)
         Data_Size(2) = No_Lines_Per_Scan

         IF (Scan_Type(Scan_No).eq.'D') THEN

            RTN=GetModisDat_MOD05_CORR(
     &        Modfil_L1B,Modfil_Geo,Modfil_CldMsk,FN_L1B,
     &        Scan_No,Buf_Size1,Buf_Size2,Data_Size,
     &        SatZen,SolZen,RelAz,Height,
     &        Refl_2,Refl_19,Ratio_Refl_19_2,
     &        Un_2,Un_19,Vflag_2,Vflag_19,
     &        Buf_Un,Buf_Sa,CldMsk,LandSea_Flag,Cloud,QA_Cloud)

            IF (rtn.ne.0) CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     &         'Call to GetModisDat_MOD05_CORR failed',
     &          'MOD_PR05_CORR_V2.f')

            IF (idebug .EQ. 1) THEN

               write(handle_QC,'(/,"Scan Number ",i3)') Scan_No
               write(handle_QC,*)


               write(handle_QC,'(A)') 'first 11 reflectances for ' //
     &                       'band 2 for 10 lines'
               write(handle_QC,'(11f11.5)') ((Refl_2(i,j),i=1,11),
     &                                    j=1,10)

               write(handle_QC,*)
               write(handle_QC,'(A)') 'first 11 reflectances for ' //
     &                             'band 19 for 10 lines'
               write(handle_QC,'(11f11.5)') ((Refl_19(i,j),i=1,11),
     &                                    j=1,10)

               write(handle_QC,*)
               write(handle_QC,'(A)') 'first 11 ratios for band 19' //
     &                             ' and band 2 for 10 lines'
               write(handle_QC,'(11f11.5)') ((Ratio_Refl_19_2(i,j),i=1,
     &                              11),j=1,10)

               write(handle_QC,*)
               write(handle_QC,'(A)') 'first 11 sat. zenith angles ' //
     &                             'for 10 lines'
               write(handle_QC,'(11f11.5)') ((SatZen(i,j),i=1,11),
     &                                    j=1,10)

               write(handle_QC,*)
               write(handle_QC,'(A)') 'first 11 solar zenith angles ' //
     &                             'for 10 lines'
               write(handle_QC,'(11f11.5)') ((SolZen(i,j),i=1,11),
     &                                    j=1,10)

               write(handle_QC,*)
               write(handle_QC,'(A)') 'first 11 rel. azimuth angles ' //
     &                             'for 10 lines'
               write(handle_QC,'(11f11.5)') ((RelAz(i,j),i=1,11),
     &                                    j=1,10)

               write(handle_QC,*)
               write(handle_QC,'(A)') 'first 11 Cloud Mask Decisions ' //
     &                             'for 10 lines'
               write(handle_QC,'(10I4)') ((CldMsk(i,j),i=1,11),
     &                                 j=1,10)

               write(handle_QC,*)
               write(handle_QC,'(A)') 'first 11 Land/Sea Flags for ' //
     &                             '10 lines'
               write(handle_QC,'(11I4)') ((LandSea_Flag(i,j),i=1,11),
     &                                  j=1,10)

               write(handle_QC,*)
               write(handle_QC,'(A)') 'first 11 Heights for 10',
     &                             ' lines'
               write(handle_QC,'(11F7.1)') ((Height(i,j),i=1,11),
     &                                  j=1,10)

            END IF

            No_Grids=No_Frames(Scan_No)/10

            CALL READER_L2(Modfil_MOD04,Modfil_MOD05,Modfil_L1B,Scan_No,
     &                     VAPOR,AEROSOL_OD)

            DO J = 1, No_Lines_Per_Scan
            DO I = 1, No_Frames(Scan_No)

              REFL_0P86(I,J)         = Refl_2(i,j)
              REFL_0P94(I,J)         = Refl_19(i,j)
              RATIO_0P94_0P86(I,J)   = Ratio_Refl_19_2(i,j)
              ZENITH_SOL(I,J)        = SolZen(i,j)
              ZENITH_VIEW(I,J)       = SatZen(i,j)
              AZIMUTH_REL(I,J)       = RelAZ(i,j)
              I_LAND_WATER(I,J)      = LandSea_Flag(i,j)
              I_CLOUD(I,J)           = CldMsk(i,j)
              VAPOR_OLD(I,J)         = Vapor(i,j)
              VAPOR_CORRECTED_FACTOR(I,J)=1.0

              i1=int((I-1)/10)+1
              j1=int((J-1)/ 10)+1

              IF (i1 .eq. NO_FRAMES(Scan_NO)) THEN
                TAU_AER (I,J) = Aerosol_OD(NO_FRAMES(Scan_NO)-1,j1)
              ELSE
                TAU_AER (I,J) = Aerosol_OD(i1,j1)
              END IF

            END DO
            END DO

            CALL PROCESS_MOD05_CORR(Data_Size)

c 05/01/2000 fhliang: added the DO-loops J and I.
            DO J = 1, No_Lines_Per_Scan
            DO I = 1, No_Frames(Scan_No)
C           IF (FACTOR_FINAL(I,J) .LE. 0.0) FACTOR_FINAL(I,J)=1.0
               Vapor_Corrected_Factor(I,J)=FACTOR_FINAL(I,J)
            END DO
            END DO

         ELSE

c 04/27/2000 fhliang: added the DO-loops J and I.
            DO J= 1,No_Lines_Per_Scan
            DO I = 1,No_Frames(Scan_No)
C              Data_Size(1) = 0
               VAPOR_CORRECTED_FACTOR(I,J)=1.0
            END DO
            END DO

         END IF


         CALL MOD05_CORR_PutData(modfil_MOD05,Data_Size,
     &        Buf_Size1,Buf_Size2,Vapor_Corrected_Factor,
     &        Max_EV_Frames, Total_No_Lines)


999      CONTINUE

      END IF

      ioflag = 'CLOSE'

      CALL FILEOC_CORR(ioflag,modis_flag,handle_corr,
     &            handle_QC,FN_L1B,
     &            modfil_L1B,modfil_Geo,modfil_CldMsk,
     &            modfil_MOD04,modfil_MOD05)

c\* Adding metadata to mod05.qc file
c rhucek 12/17/97: commented 8 lines, added call to function
c Set_CoreMetadata_QC and changed error message.
c
c     rtn = pgs_met_init(LRN_MCFQC,groups)
c
c     rtn = pgs_met_write(groups(ODL_IN_MEMORY),' ',LRN_QCMET)

      rtn = Set_CoreMetadata_QC(LRN_MCFQC,LRN_QCMET)

      if (rtn.ne.mapiok) then
         msgbuf =
     1   'Set_CoreMetadata failed in MOD05 QC file'
     2   // char(10) // 'Operator Action:  Check for valid MCF file and PCF reference '
     3   // char(10) // 'to MCF file.  If MCF/PCF files are wrong or corrupted, stage '
     4   // char(10) // 'correct files and rerun PGE.  Otherwise, notify SDST.'
         call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,
     &   'MOD_PR05_CORR_V2.f')
      else
         call modis_smf_setdynamicmsg(MODIS_S_SUCCESS,
     2   'Create metadata file to non-HDF (mod5C_QC) without error',
     3   'MOD_PR05_CORR_V2')
      end if

      RETURN
      END



C--------------------------------------------------------------------------
        SUBROUTINE READ_REFL_0P86_TABLE_AEROSOL(handle_corr)
C---------------------------------------------------------------------------
C!F77
C
C!DESCRIPTION: This subroutine is to read in the table of the surface
C              reflectance at 0.86 micron.
C
c!INPUT PARAMETERS: handle_corr
c
c!OUTPUT PARAMETERS: REFL_0P86_T
c
c!REVISION HISTORY:
c 02/02/98 achu
c added arguments in INPUT / OUTPUT PARAMETERS sections in prolog.
c
c 01/28/98 fhliang
c fixed prolog.
c
c!TEAM-UNIQUE HEADER:
c
c Developed by MODIS Aerosol/Solar Water Vapor Retrieval Team
c GSFC, Greenbelt, MD
c
C !References and Credits:
C
C        Written by
C        Dr. Wei Han    06/20/97
C        han@neptune.nrl.navy.mil
C
C DESIGN NOTES: data and variables are passed with common statements
C               contained in 'MOD05_CORR.inc'.
C
c!END----------------------------------------------------------------------
C
        IMPLICIT NONE

        INCLUDE 'MOD05_CORR.inc'
          SAVE

        INTEGER handle_corr(2),I1,I2,I3,I4,I5,I6

        DO  I1 = 1, N_REFL_SUR
        DO  I2 = 1, N_VAPOR
        DO  I3 = 1, N_TAU
        DO  I4 = 1, N_ZENITH_SOL
        DO  I5 = 1, N_AZIMUTH_REL

        READ(handle_corr(1), *) (REFL_0P86_T(I1, I2, I3, I4, I5, I6),
     &                           I6=1,N_ZENITH_VIEW)

        END DO
        END DO
        END DO
        END DO
        END DO

        RETURN
        END


C---------------------------------------------------------------------------
        SUBROUTINE READ_RATIO_TABLE_AEROSOL(handle_corr)
C
C---------------------------------------------------------------------------
C!F77
C
C!DESCRIPTION: This subroutine is to read in the table of the ratio of 0.94
C              to 0.86 micron.
C
c!INPUT PARAMETERS: handle_corr
c
c!OUTPUT PARAMETERS: RATIO_T_AER
c
c!REVISION HISTORY:
c 02/02/98 achu
c added arguments in INPUT / OUTPUT PARAMETERS sections in prolog.
c
c 01/28/98 fhliang
c fixed prolog.
c
c!TEAM-UNIQUE HEADER:
c
c Developed by MODIS Aerosol/Solar Water Vapor Retrieval Team
c GSFC, Greenbelt, MD
c
C !References and Credits:
C
C        Written by
C        Dr. Wei Han    06/20/97
C        han@neptune.nrl.navy.mil
C
C DESIGN NOTES: data and variables are passed with common statements
C               contained in 'MOD05_CORR.inc'.
C
c!END-----------------------------------------------------------------------
C
        IMPLICIT NONE

        INCLUDE 'MOD05_CORR.inc'
         SAVE

        INTEGER handle_corr(2),I1,I2,I3,I4,I5,I6

        DO I1 = 1, N_VAPOR
        DO I2 = 1, N_REFL_SUR
        DO I3 = 1, N_TAU
        DO I4 = 1, N_ZENITH_SOL
        DO I5 = 1, N_AZIMUTH_REL

        READ(handle_corr(2), *) (RATIO_T_AER(I1, I2, I3, I4, I5, I6),
     &                          I6=1,N_ZENITH_VIEW)

        END DO
        END DO
        END DO
        END DO
        END DO

        RETURN
        END

C*************************************************************************

      SUBROUTINE PROCESS_MOD05_CORR(Data_Size)
C-----------------------------------------------------------------------
C!F77
C
C!Description:This subroutine is to correct vertical column water vapor
C              values using a lookup table procedure. Data from one MODIS
C              swath with 10 lines are processed.
C
C!Input Parameters: MODIS L1B data, MOD04 aerosol optical thickness,
C                   MOD05 total precipitable water, cloud mask,
C                   land and water masks, solar and view zenith angles,
C                   and relative azimuth angle
C
C                   Table of surface reflectance, table of ratio near 1 micron
C                   with transmission-only method, table of ratio near 1
C                   micron with aerosol effect, MODIS reflectance at 0.86 and
C                   0.94 micron, cloud mask, solar and view zenith angles,
C                   relative azimuth angles, and aerosol optical depth.
C
C!OUTPUT PARAMETERS: correction factor or initial correction factor
C                    at each pixel.
C
c!REVISION HISTORY:
c 01/28/98 fhliang
c 06/20/97 Wei Han
c fixed prolog.
c
c!TEAM-UNIQUE HEADER:
c
c Developed by MODIS Aerosol/Solar Water Vapor Retrieval Team
c GSFC, Greenbelt, MD
c
C !References and Credits:
C
C        Written by
C        Drs. Shaohua Shen and Bo-Cai Gao 09/07/99
C        shen@climate.gsfc.nasa.gov
C
C DESIGN NOTES: data and variables are passed with common statements
C               contained in 'MOD05_CORR.inc'.
C
C Internal subroutines called:
C		FIND_VAPOR_TRANS_ONLY
C		FIND_VAPOR_AER
C               GET_FACTOR
C
c!END-----------------------------------------------------------------------
C
        IMPLICIT NONE

        INCLUDE 'MOD05_CORR.inc'
         SAVE
        INTEGER I,J,Data_Size(2),handle_corr(2)

C
C Setting up initial values of factor and factor2 at full resolution
C (pixel-by-pixel)
C I_CLOUD(I,J)=0, it means there are cloud and  I_CLOUD(I,J)=1, it means
C clear atmosphere.
        izs=1
        izr=1
        izv=1
        itau=1
        ivp=1
        iref=1
        iref2=1
        DO J = 1, Data_Size(2)
        DO I = 1, Data_Size(1)


          FACTOR(I,J)  = 1.0
          FACTOR2(I,J) = 1.0

        if ((I_CLOUD(I,J) .EQ. 1) .and. ((REFL_0P86(I,J) .LE. 0.1) .or.
     $  ((REFL_0P86(I,J) .GE. 0.5) .and. (TAU_AER(I,J) .GE. 0.5)))) then
            FACTOR(I,J)=0.0

            CALL FIND_VAPOR_TRANS_ONLY(I,J)

            CALL FIND_VAPOR_AER(I,J)

            if (VAPOR_TRANS_ONLY(I,J) .GT. 0.0) then
              FACTOR2(I,J)=VAPOR_AER(I,J)/VAPOR_TRANS_ONLY(I,J)
            else
             FACTOR2(I,J)=1.0
            end if
            IF (FACTOR2(I,J) .GE. 2.0) FACTOR2(I,J)=1.0

          end if

        END DO
        END DO
        CALL GET_FACTOR(Data_size)

        RETURN
        END

C---------------------------------------------------------------------------
        SUBROUTINE FIND_VAPOR_TRANS_ONLY(I, J)
C---------------------------------------------------------------------------
C!F77
C
C!DESCRIPTION: This subroutine derives transmission-only water vapor amount
C              values for one pixel using a look-up table.
C
C!INPUT PARAMETERS: table of ratio near 1 micron with transmission-only
C                   method, solar and view zenith angles, MODIS reflectance
C                   at 0.86 and 0.94 micron.
C
C!OUTPUT PARAMETERS: transmission-only water vapor amount for one pixel.
C
c!REVISION HISTORY:
c 01/28/98 fhliang
c 06/20/97 Wei Han
c fixed prolog.
c
c!TEAM-UNIQUE HEADER:
c
c Developed by MODIS Aerosol/Solar Water Vapor Retrieval Team
c GSFC, Greenbelt, MD
c
C !References and Credits:
C
C        Written by
C        Drs. Shaohua Shen and Bo-Cai Gao 90/07/99
C        shen@climate.gsfc.nasa.gov
C
C
C DESIGN NOTES: data and variables are passed with common statements
C               contained in 'MOD05_CORR.inc'.
C
c!END-----------------------------------------------------------------------
C
        IMPLICIT NONE

        INCLUDE 'MOD05_CORR.inc'
        SAVE

        INTEGER I,J,K
        REAL RATIO,GEOM_FACTOR,DEG_TO_RAD
        REAL VAPOR_T(N_VAPOR), FIT
        REAL VAPOR_TOTAL(N_VAPOR)
        REAL RATIO_T_TRANS_ONLY(N_VAPOR)
        DATA VAPOR_T/0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0,
     $               7.0, 10.0, 15.0/
        DATA RATIO_T_TRANS_ONLY/0.8668932,0.8111669,0.7397158,0.6696488,
     $      0.6201043,0.5601289,0.5157256,0.4515401,0.4060718,0.3707590,
     $      0.3193024,0.2673745,0.2129472/
        DATA DEG_TO_RAD /0.0174533/



        RATIO = RATIO_0P94_0P86(I,J)

C
C Use two-way vertical water vapor amount when looking up the table
C

        DO K = 1, N_VAPOR

          VAPOR_TOTAL(K) = 2.0 * VAPOR_T(K)

        END DO

        DO K = 1, N_VAPOR-1

          IF ( RATIO .LT. RATIO_T_TRANS_ONLY(K) + ZEPS .AND.
     $         RATIO .GT. RATIO_T_TRANS_ONLY(K+1) -ZEPS ) THEN

          VAPOR_TRANS_ONLY(I,J) = FIT(RATIO_T_TRANS_ONLY(K),
     $            RATIO_T_TRANS_ONLY(K+1), VAPOR_TOTAL(K),
     $                           VAPOR_TOTAL(K+1), RATIO)

        ENDIF

        END DO

        IF ( RATIO .GE. RATIO_T_TRANS_ONLY(1) + ZEPS)
     $      VAPOR_TRANS_ONLY(I, J) = VAPOR_TOTAL(1)

        IF ( RATIO .LE. RATIO_T_TRANS_ONLY(N_VAPOR) - ZEPS)
     $      VAPOR_TRANS_ONLY(I, J) = VAPOR_TOTAL(N_VAPOR)


C Compute geometrical factor related to the solar and view zenith angles


        GEOM_FACTOR  =   1. / COS(ZENITH_SOL(I,J)  * DEG_TO_RAD)
     $                 + 1. / COS(ZENITH_VIEW(I,J) * DEG_TO_RAD)

C
C Derive one way vertical water vapor amount (in cm)
C

        VAPOR_TRANS_ONLY(I,J) = VAPOR_TRANS_ONLY(I,J) / GEOM_FACTOR

        RETURN
        END


C---------------------------------------------------------------------------
        REAL FUNCTION FIT(X1, X2, Y1, Y2, X)
C---------------------------------------------------------------------------
C!F77
C
C!DESCRIPTION: This subroutine does a one-dimensional linear interpolation job.
C
C!INPUT PARAMETERS: (x1,y1) and (x2,y2), x
C
C!OUTPUT PARAMETERS: y
C
c!REVISION HISTORY:
c 01/28/98 fhliang
c fixed prolog.
c
c!TEAM-UNIQUE HEADER:
c
c Developed by MODIS Aerosol/Solar Water Vapor Retrieval Team
c GSFC, Greenbelt, MD
c
C !References and Credits:
C
C        Written by
C        Dr. Wei Han    06/20/97
C        han@neptune.nrl.navy.mil
C
c!END-----------------------------------------------------------------------
C
        IMPLICIT NONE

        REAL X,X1,X2,Y1,Y2

        FIT = (Y2 - Y1) * (X - X1) / (X2 - X1) + Y1

        RETURN
        END

C---------------------------------------------------------------------------
        SUBROUTINE GET_FACTOR(Data_Size)
C---------------------------------------------------------------------------
C!F77
C
C!DESCRIPTION: This subroutine yields correction factors and corrected column
C              water vapor amount values at each pixel for one MODIS scan.
C
c!INPUT PARAMETERS: correction factor derived from subroutine INTERPOLATION,
C                   and initial correction factor at each pixel.
C
C!OUTPUT PARAMETERS: correction factor and corrected column water vapor
C                    amount at each pixel.
C
c!REVISION HISTORY:
c 01/28/98 fhliang
c 06/20/97 Wei Han
c fixed prolog.
c
c!TEAM-UNIQUE HEADER:
c
c Developed by MODIS Aerosol/Solar Water Vapor Retrieval Team
c GSFC, Greenbelt, MD
c
C !References and Credits:
C
C        Written by
C        Drs. Shaohua Shen and Bo-Cai Gao 09/07/99
C        shen@climate.gsfc.nasa.gov
C
C DESIGN NOTES: data and variables are passed with common statements
C               contained in 'MOD05_CORR.inc'.
C
c!END-----------------------------------------------------------------------
C
        IMPLICIT NONE

        INCLUDE 'MOD05_CORR.inc'

        INTEGER I,J,Data_Size(2)

        SAVE

        DO J = 1, Data_Size(2)
        DO I = 1, Data_Size(1)

        IF ((FACTOR(I,J)+ZEPS-1.0).GE.0.0) FACTOR_FINAL(I,J) = 1.0
        IF ((FACTOR(I,J)-ZEPS)/ZEPS .LT. 0.0) THEN
          FACTOR_FINAL(I,J) = FACTOR2(I,J)
        ELSE
          FACTOR_FINAL(I,J) = 1.0
        ENDIF

        IF (FACTOR_FINAL(I,J) .LT. 0.0) FACTOR_FINAL(I,J)=1.0
        END DO
        END DO

C        write(6,*) FACTOR_FINAL
        RETURN
        END


      Subroutine Find_VAPOR_AER(I,J)
C---------------------------------------------------------------------------
C!F77
C
C!DESCRIPTION: This subroutine derives surface reflectance values at 0.86 micron
C              for one pixel using a look up table procedure.
C
C!INPUT PARAMETERS: table of reflectance at 0.86 micron, MODIS reflectance
C                   at 0.86 micron, aerosol optical depth, solar and view
C                   zenith angles, relative azimuth angles, transmission
C                   -only water vapor amount.
C
C!OUTPUT PARAMETERS: surface reflectance at 0.86 micron for one pixel.
C
c!REVISION HISTORY:
c 01/28/98 fhliang
c fixed prolog.
c
c!TEAM-UNIQUE HEADER:
c
c Developed by MODIS Aerosol/Solar Water Vapor Retrieval Team
c GSFC, Greenbelt, MD
c
C !References and Credits:
C
C        Written by
C        Drs. S. Shen and Bo-Cai Gao, 09/07/99
C        shen@climate.gsfc.nasa.gov
C
C DESIGN NOTES: data and variables are passed with common statements
C               contained in 'MOD05_CORR.inc'.
C
C Functions called:
C             FIT
C             HUNT,
C             LOCATE2
C             RESUCE_DIM
C
c!END-----------------------------------------------------------------------
        IMPLICIT NONE
        INCLUDE 'MOD05_CORR.inc'

C        Real ts, tv, tr, tf, tp,  tt,tf2
C        integer izs, izr, izv, iref, ivp, itau,iref2
        SAVE

        Integer  k, Iflag,Scan_No,Rscan_No,i,j, locate2
        Real  z_sol, az_rel, z_view, t_aer, v_trans, r_0p86, r_sur,
     $          r_op9486, fit

        Real VAPOR_T(N_VAPOR), TAU_AER_T(N_TAU),
     $    REFL_SUR_0P86_T(N_REFL_SUR), ZENITH_SOL_T(N_ZENITH_SOL),
     $    ZENITH_VIEW_T(N_ZENITH_VIEW), AZIMUTH_REL_T(N_AZIMUTH_REL),
     $    ratio_vapor(N_VAPOR), refl_atmos(N_REFL_SUR)

        DATA VAPOR_T/0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0,
     $               7.0, 10.0, 15.0/
        DATA TAU_AER_T/0.0, 0.1, 0.25, 0.5, 1.0, 2.5, 3.5, 5.0/
        DATA REFL_SUR_0P86_T/0.01, 0.03, 0.05, 0.1, 0.2, 0.4, 0.6, 1.0/
        DATA ZENITH_SOL_T/0.0, 10.0, 20.0, 30.0,40.0, 50.0, 60.0, 70.0,
     $                   75.0, 80.0/
        DATA ZENITH_VIEW_T/0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0/
        DATA AZIMUTH_REL_T/0.0, 5.0, 10.0, 20.0, 30.0, 40.0, 60.0, 80.0,
     $                     90.0, 100.0, 110.0, 120.0, 130.0, 140.0,
     $                     150.0, 160.0, 170.0, 175.0, 180.0/


C The most important thing is to chose correct index for all six parameters.
C if the index we choose is correct, then we can directly get the
C surface reflectance  at 0.86 or water vapor using a lookup table and
C muti_dimension interpolation.

C First, choose the correct index for Zenith_solar.

        Iflag=1
C        z_sol=zenith_sol(i,j)
C        izs=2
C        izs=Locate1(zenith_sol_t, zenith_sol(i,j))
         call hunt1(zenith_sol_t, N_ZENITH_SOL, zenith_sol(i,j), izs)
c     t is the fractional distance from lower part of the interval to
c     in the particular interval that zenith_sol is in. If outside our
c     tabulated values, we just set it to the max/min.
       if (zenith_sol(i,j) .lt. zenith_sol_t(1)) then
           ts = 0
           izs = 1
       else if (zenith_sol(i,j) .ge.
     $     zenith_sol_t(N_ZENITH_SOL)) then
           ts = 1
          izs = N_ZENITH_SOL - 1
       endif
       if( (zenith_sol(i,j) .ge. zenith_sol_t(1)) .and.
     $     (zenith_sol(i,j) .lt. zenith_sol_t(N_ZENITH_SOL)) ) then
       ts =(zenith_sol(i,j) - zenith_sol_t(izs))/
     $     (zenith_sol_t(izs+1) -  zenith_sol_t(izs))
       endif

c     Now choose index for azimuth_rel

C      az_rel=azimuth_rel(i,j)
C      izr=1
C      izr=Locate1(azimuth_rel_t, azimuth_rel(i,j))
      call hunt1(azimuth_rel_t,N_AZIMUTH_REL, azimuth_rel(i,j), izr)
       if( izr .ge. N_AZIMUTH_REL) izr =  N_AZIMUTH_REL - 1
       if (azimuth_rel(i,j) .lt. azimuth_rel_t(1)) then
          tr = 0
         izr = 1
       else if (azimuth_rel(i,j) .ge.
     $    azimuth_rel_t(N_AZIMUTH_REL)) then
          tr = 1
         izr = N_AZIMUTH_REL - 1
       endif
       if( (azimuth_rel(i,j) .ge. azimuth_rel_t(1)) .and.
     $     (azimuth_rel(i,j) .lt. azimuth_rel_t(N_AZIMUTH_REL)) ) then
      tr=(azimuth_rel(i,j) - azimuth_rel_t(izr))/
     $   (azimuth_rel_t(izr+1) - azimuth_rel_t(izr))
       endif


C  Now choose an index for zenith zenith_view


C        z_view=zenith_view(i,j)
C        izv=1
C       izv=Locate1(zenith_view_t, zenith_view(i,j))
        call hunt1(zenith_view_t,N_ZENITH_VIEW, zenith_view(i,j), izv)
        if (zenith_view(i,j) .lt. zenith_view_t(1)) then
         tv = 0
        izv = 1
        else if (zenith_view(i,j) .ge.
     $          zenith_view_t(N_ZENITH_VIEW)) then
         tv = 1
        izv = N_ZENITH_VIEW - 1
        endif
        if( (zenith_view(i,j) .ge. zenith_view_t(1)) .and.
     $      (zenith_view(i,j) .lt. zenith_view_t(N_ZENITH_VIEW)) ) then
        tv = (zenith_view(i,j)-zenith_view_t(izv))/
     $     (zenith_view_t(izv+1) - zenith_view_t(izv))
        endif


C  Now chose an index for aerosol optical thickness

C        t_aer=tau_aer(i,j)
C        itau=1
C       itau=Locate1(tau_aer_t, tau_aer(i,j))
        call hunt1(tau_aer_t, N_TAU, tau_aer(i,j), itau)
        if(tau_aer(i,j) .lt.  tau_aer_t(1)) then
         tt   = 0
         itau = 1
        else if (tau_aer(i,j) .ge. tau_aer_t(N_tau)) then
         tt   = 1
         itau = N_TAU - 1
        endif
        if( (tau_aer(i,j) .ge. tau_aer_t(1)) .and.
     $      (tau_aer(i,j) .lt. tau_aer_t(N_tau)) ) then
        tt=(tau_aer(i,j) - tau_aer_t(itau))/(tau_aer_t(itau+1)
     $      - tau_aer_t(itau))
        endif

C  Now chose an index for vapor_trans_only

C        v_trans=vapor_trans_only(i,j)
C        ivp= 2
C       ivp=Locate1(vapor_t, vapor_trans_only(i,j))
         call hunt1(vapor_t,N_VAPOR, vapor_trans_only(i,j),ivp)
        if(vapor_trans_only(i,j) .lt. vapor_t(1)) then
            tp = 0
           ivp = 1
        else if (vapor_trans_only(i,j) .ge. vapor_t(N_VAPOR)) then
            tp = 1
           ivp = N_VAPOR - 1
        endif
        if( (vapor_trans_only(i,j) .ge. vapor_t(1)) .and.
     $      (vapor_trans_only(i,j) .lt. vapor_t(N_VAPOR)) ) then
        tp=(vapor_trans_only(i,j) -vapor_t(ivp))
     $      / (vapor_t(ivp+1) - vapor_t(ivp))
        endif

C Based on above 5 index, now we reduce the previous six dimensional array to a
C one dimensional array of the reflectance at atmospheric top at 0.86um.

      call Reduce_dim (refl_0p86_t, n_refl_sur, n_vapor, n_tau,
     $        n_zenith_sol, n_azimuth_rel, n_zenith_view, tp, tt,
     $        ts, tr, tv, ivp, itau, izs, izr, izv, refl_atmos)


C       write(7,*) 'refl_atmos=', (refl_atmos(k), k=1, N_REFL_SUR)


C  Now chose an index for reflectance of atmospheric top  at 0.86 micron,
C  if reflectance  at atmospheric top is monotonic.

C      r_0p86=refl_0p86(i,j)
      if (Iflag .eq. 1) then
C      iref=1
C      iref=Locate1(refl_atmos, refl_0p86(i,j))
      CALL hunt1(refl_atmos, N_REFL_SUR, refl_0p86(i,j), iref)
        if(refl_0p86(i,j) .lt. refl_atmos(1)) then
           tf = 0
         iref = 1
        else if (refl_0p86(i,j) .ge. refl_atmos(N_refl_sur)) then
           tf = 1
         iref = N_REFL_SUR - 1
        endif
        if( (refl_0p86(i,j) .ge. refl_atmos(1)) .and.
     $      (refl_0p86(i,j) .lt. refl_atmos(N_refl_sur)) ) then
        tf=(refl_0p86(i,j) - refl_atmos(iref))/(refl_atmos(iref+1) -
     $    refl_atmos(iref))
        endif

C    Computer the surface reflectance at 0.86 in term of lookup table.

      if (iref .lt. 1)  Refl_sur_0p86(i,j) = refl_sur_0p86_t(1)
      if (iref .ge. N_refl_sur)  Refl_sur_0p86(i,j) 
     $                         = refl_sur_0p86_t(N_refl_sur)
      if( (iref .ge. 1) .and. (iref .lt. N_refl_sur) ) then
       Refl_sur_0p86(i,j)=(refl_sur_0p86_t(iref+1) -
     $     refl_sur_0p86_t(iref))*tf + refl_sur_0p86_t(iref)
      endif

C if reflectance  at atmospheric top is not monotonic, using following methods:

       else
        do k=1, N_REFL_SUR-1
         if(refl_0p86(i,j) .gt. refl_atmos(k) -zeps
     $        .and. refl_0p86(i,j) .lt. refl_atmos(k) -zeps) then
           Refl_sur_0p86(i,j)=fit(refl_atmos(k), refl_atmos(k+1),
     $      refl_sur_0p86_t(k),  refl_sur_0p86_t(k+1),refl_0p86(i,j))
           endif
        end do
          if (refl_0p86(i,j) .lt. refl_atmos(1) +zeps) then
             Refl_sur_0p86(i,j)=refl_sur_0p86_t(1)
          else if (refl_0p86(i,j) .ge. refl_atmos(N_REFL_SUR)) then
             Refl_sur_0p86(i,j)=refl_sur_0p86_t(N_REFL_SUR)
          endif
        endif


C     Again chose an index for the new surface reflectance.

C        r_sur=Refl_sur_0p86(i,j)
C        iref2 = 2
C      iref2=Locate1(refl_sur_0p86_t, Refl_sur_0p86(i,j))
        call hunt1(refl_sur_0p86_t, N_REFL_SUR,
     &           Refl_sur_0p86(i,j), iref2)
        if (Refl_sur_0p86(i,j) .lt. refl_sur_0p86_t(1)) then
            tf2 = 0
          iref2 = 1
        else if (Refl_sur_0p86(i,j)  .ge.
     $          refl_sur_0p86_t(N_refl_sur)) then
            tf2 = 1
          iref2 =  N_REFL_SUR - 1
        endif
        if((Refl_sur_0p86(i,j) .ge. refl_sur_0p86_t(1)) .and.
     $     (Refl_sur_0p86(i,j) .lt. refl_sur_0p86_t(N_refl_sur))) then
        tf2=(Refl_sur_0p86(i,j) - refl_sur_0p86_t(iref2))/
     $  (refl_sur_0p86_t(iref2+1) - refl_sur_0p86_t(iref2))
        endif

C Based on above 5 index, now we reduce the previous six dimensional array to a
C one dimensional array of the ratio of water vapor.

      call Reduce_dim (ratio_t_aer, n_vapor, n_refl_sur, n_tau,
     $      n_zenith_sol, n_azimuth_rel, n_zenith_view,tf2,tt,
     $      ts, tr,tv, iref2, itau,izs,izr,izv, ratio_vapor)


C        write(7,*) 'AA=', (ratio_vapor(k), k=1, n_vapor)

C   Chose an index for the Ratio of reflectance.

C        r_op9486=ratio_0p94_0p86(i,j)
        ivp=Locate2(ratio_vapor,N_VAPOR, ratio_0p94_0p86(i,j))
        if(ratio_0p94_0p86(i,j) .gt. ratio_vapor(1)) then
            tp=0
            ivp=1
        else if (ratio_0p94_0p86(i,j)  .le.
     $           ratio_vapor(N_vapor)) then
            tp=1
            ivp=N_vapor-1
        endif
        if((ratio_0p94_0p86(i,j) .le. ratio_vapor(1)) .and.
     $     (ratio_0p94_0p86(i,j) .gt. ratio_vapor(N_vapor)) )then
        tp=(ratio_0p94_0p86(i,j) - ratio_vapor(ivp))
     $    /(ratio_vapor(ivp+1)   - ratio_vapor(ivp))
        endif


C  Computer the water Vapor using interpolation.

      if (ivp .lt. 1) Vapor_aer(i,j) = vapor_t(ivp)
      if (ivp .ge. N_vapor) Vapor_aer(i,j) = vapor_t(N_vapor)
      if((ivp .ge. 1) .and. (ivp .lt. N_vapor) ) then
       Vapor_aer(i,j)=(vapor_t(ivp+1) - vapor_t(ivp))*tp+vapor_t(ivp)
      endif

       RETURN
       End
C-----------------------------------------------------------------------

      Subroutine Reduce_dim (refl_table, n1,n2,n3,n4,n5,n6,
     &    tp, tt, ts, tr, tv, ivp, itau,izs,izr,izv, Tab_value)

C-----------------------------------------------------------------------
C!F77
C
C!Description:   This subroutine is used to drive a one dimensional
C                array of reflectance at the top of the atmospghere at
c                0.86 micron or the ratio of water vapor at  0.86 micron
C                from a six dimensional array.
C
C!Input Parameters: Table of surface reflectance, table of ratio near 1 micron
C                   with transmission-only method, table of ratio near 1
C                   micron with aerosol effect, MODIS reflectance at 0.86 and
C                   0.94 micron, the indexs of solar and view zenith angles,
C                   relative azimuth angles, and aerosol optical depth.
C
C!OUTPUT PARAMETERS: one dimensional array of reflectance at the top of the
C                    atmospghere 0.86 micron or the ratio of water vapor
C                    at 0.86 micron.
C
C
c!REVISION HISTORY:
c 09/07/99 Shaohua Shen
c fixed prolog.
c
c!TEAM-UNIQUE HEADER:
c
c Developed by MODIS Aerosol/Solar Water Vapor Retrieval Team
c GSFC, Greenbelt, MD
c
c
C !References and Credits:
C
C        Written by
C        Drs. Shaohua Shen and Bo-Cai Gao  09/07/99
C        Shen@climate.gsfc.nasa.gov
C
C DESIGN NOTES: data and variables are passed with common statements
C               contained in 'MOD05_CORR.inc'.
C
C Internal subroutines called:
C              NONE
C
c!END-----------------------------------------------------------------------

      Implicit None

      Integer n1,n2,n3,n4,n5,n6,ii
      Real refl_table(n1,n2,n3,n4,n5,n6)
      Real Tab_value(n1)
      Real tp,tt,ts,tr,tv,tp1,tt1,ts1,tr1,tv1
      Integer ivp,itau,izs,izr,izv, ivp1,itau1, izs1,izr1,izv1

C assume That we get all indexes in the six component. Now we calculate
C the surface_reflectance or Ratio  using multi_dimension interpolation and t
C the look up ! table.

      tp1=1-tp
      tt1=1-tt
      ts1=1-ts
      tr1=1-tr
      tv1=1-tv

      ivp1=ivp+1
      itau1=itau+1
      izs1=izs+1
      izr1=izr+1
      izv1=izv+1

      do ii=1,n1
      Tab_value(ii) = (tp1*tt1*ts1*tr1*tv1)*
     &      refl_table(ii,ivp,itau,izs,izr,izv)
     & + (tp1*tt*ts1*tr1*tv1)*refl_table(ii,ivp,itau1,izs,izr,izv)
     & + (tp1*tt1*ts*tr1*tv1)*refl_table(ii,ivp,itau,izs1,izr,izv)
     & + (tp1*tt1*ts1*tr*tv1)*refl_table(ii,ivp,itau,izs,izr1,izv)
     & + (tp1*tt1*ts1*tr1*tv)*refl_table(ii,ivp,itau,izs,izr,izv1)
     & + (tp*tt1*ts1*tr1*tv1)*refl_table(ii,ivp1,itau,izs,izr,izv)
     & + (tp*tt*ts1*tr1*tv1)*refl_table(ii,ivp1,itau1,izs,izr,izv)
     & + (tp*tt1*ts*tr1*tv1)*refl_table(ii,ivp1,itau,izs1,izr,izv)
     & + (tp*tt1*ts1*tr*tv1)*refl_table(ii,ivp1,itau,izs,izr1,izv)
     & + (tp*tt1*ts1*tr1*tv)*refl_table(ii,ivp1,itau,izs,izr,izv1)
     & + (tp1*tt*ts*tr1*tv1)*refl_table(ii,ivp,itau1,izs1,izr,izv)
     & + (tp1*tt*ts1*tr*tv1)*refl_table(ii,ivp,itau1,izs,izr1,izv)
     & + (tp1*tt*ts1*tr1*tv)*refl_table(ii,ivp,itau1,izs,izr,izv1)
     & + (tp1*tt1*ts*tr*tv1)*refl_table(ii,ivp,itau,izs1,izr1,izv)
     & + (tp1*tt1*ts*tr1*tv)*refl_table(ii,ivp,itau,izs1,izr,izv1)
     & + (tp1*tt1*ts1*tr*tv)*refl_table(ii,ivp,itau,izs,izr1,izv1)
     & + (tp*tt*ts*tr1*tv1)*refl_table(ii,ivp1,itau1,izs1,izr,izv)
     & + (tp*tt*ts1*tr*tv1)*refl_table(ii,ivp1,itau1,izs,izr1,izv)

      Tab_value(ii) = Tab_value(ii)
     & + (tp*tt*ts1*tr1*tv)*refl_table(ii,ivp1,itau1,izs,izr,izv1)
     & + (tp*tt1*ts*tr*tv1)*refl_table(ii,ivp1,itau,izs1,izr1,izv)
     & + (tp*tt1*ts*tr1*tv)*refl_table(ii,ivp1,itau,izs1,izr,izv1)
     & + (tp*tt1*ts1*tr*tv)*refl_table(ii,ivp1,itau,izs,izr1,izv1)
     & + (tp1*tt*ts*tr*tv1)*refl_table(ii,ivp,itau1,izs1,izr1,izv)
     & + (tp1*tt*ts*tr1*tv)*refl_table(ii,ivp,itau1,izs1,izr,izv1)
     & + (tp1*tt*ts1*tr*tv)*refl_table(ii,ivp,itau1,izs,izr1,izv1)
     & + (tp1*tt1*ts*tr*tv)*refl_table(ii,ivp,itau,izs1,izr1,izv1)
     & + (tp1*tt*ts*tr*tv)*refl_table(ii,ivp,itau1,izs1,izr1,izv1)
     & + (tp*tt1*ts*tr*tv)*refl_table(ii,ivp1,itau,izs1,izr1,izv1)
     & + (tp*tt*ts1*tr*tv)*refl_table(ii,ivp1,itau1,izs,izr1,izv1)
     & + (tp*tt*ts*tr1*tv)*refl_table(ii,ivp1,itau1,izs1,izr,izv1)
     & + (tp*tt*ts*tr*tv1)*refl_table(ii,ivp1,itau1,izs1,izr1,izv)
     & + (tp*tt*ts*tr*tv)*refl_table(ii,ivp1,itau1,izs1,izr1,izv1)


       end do
      return
      End


      INTEGER FUNCTION  Locate1(xx,n, x)
C-----------------------------------------------------------------------
C!F77
C
C!Description:This subroutine is to find the index of the given
C             solar and view zenith angles, relative azimuth angles,
C             and aerosol optical depth in the related look up table.
C             Here XX is an  array XX with length N.  Given a value X,
C             return a value J such that X is between XX(J) and XX(J+1). XX
C             must be monotonic and  increase. J=0 or J=N is return to
C             indicate that X is out of range.
C
C
C!Input Parameters:  solar angle (view zenith angle, relative azimuth angle,
C                    or aerosol optical depth, and the related
C                    look up table.
C
C!OUTPUT PARAMETERS: Index of the given solar angle (view zenith angle,
C                    relative azimuth angle, or aerosol optical depth)
C                    in the look up table
C
c!REVISION HISTORY:
c 09/07/99  Shaohua Shen
c fixed prolog.
c
c!TEAM-UNIQUE HEADER:
c
c Developed by MODIS Aerosol/Solar Water Vapor Retrieval Team
c GSFC, Greenbelt, MD
c
C !References and Credits:
C
C        Written by
C        Drs. Shaohua Shen and Bo-Cai Gao 09/07/99
C        shen@climate.gsfc.nasa.gov
C
C DESIGN NOTES: data and variables are passed with common statements
C               contained in 'MOD05_CORR.inc'.
C
C
c!END---------------------------------------------------------------------

      Implicit None

      Integer n,jl,ju,jm
      Real xx(n)
      Real x
      Logical  ascnd

      ascnd = (xx(n) .ge. xx(1))
      jl=1
      ju=n
  10  if(ju-jl .gt. 1) then
           jm=(ju+jl)/2
           if( ascnd .eqv. (x .ge.  xx(jm))) then
             jl=jm
           else
             ju=jm
           endif
       go to 10
       endif
      if ((x .le. xx(1)) .and. ascnd) then
           locate1=1
      else if ((x .gt. xx(n)) .and. ascnd) then
             locate1=n
      else
          locate1=jl
      end if
      return
      End


      INTEGER FUNCTION  Locate2(xx, n, x)
C-----------------------------------------------------------------------
C!F77
C
C!Description:This subroutine is to find the index of the given
C             ratio of reflectance at 0.94 and 0.86um in the look up table.
C             Here XX is an  array XX with length N.  Given a value X,
C             return a value J such that X is between XX(J) and XX(J+1). XX
C             must be monotonic and  decrease. J=0 or J=N is return to
C             inidicate that X is out of range.
C
C
C!Input Parameters:  ratio of reflectance of 0.94 to 0.86um, and the
C                   look up table for the ratio.
C
C!OUTPUT PARAMETERS: Index of the given ratio of reflectance at 0.94
C                    and 0.86um in the look up table
C
c!REVISION HISTORY:
c 09/07/99  Shaohua Shen
c fixed prolog.
c
c!TEAM-UNIQUE HEADER:
c
c Developed by MODIS Aerosol/Solar Water Vapor Retrieval Team
c GSFC, Greenbelt, MD
c
C !References and Credits:
C
C        Written by
C        Drs. Shaohua Shen and Bo-Cai Gao 09/07/99
C        shen@climate.gsfc.nasa.gov
C
C DESIGN NOTES: data and variables are passed with common statements
C               contained in 'MOD05_CORR.inc'.
C
C
c!END---------------------------------------------------------------------

      Implicit None

      Integer n,jl,ju,jm
      Real xx(n)
      Real x
      Logical dscnd

      dscnd = (xx(n) .lt. xx(1))
      jl=1
      ju=n
 10   if(ju-jl .gt. 1) then
           jm=(ju+jl)/2
           if ( dscnd .eqv. (x .le. xx(jm))) then
             jl=jm
           else
             ju=jm
           endif
      go to 10
      endif
      if  ((x .ge. xx(1)) .and. dscnd) then
           locate2=1
      else if ((x .le. xx(n)) .and. dscnd) then
             locate2=n
      else
          locate2=jl
      end if
      return
      End


	SUBROUTINE hunt1(xx,n, x, jlo)
C-----------------------------------------------------------------------
C!F77
C
C!Description:This subroutine is to find the index of the given
C             solar and view zenith angles, relative azimuth angles,
C             and aerosol optical depth in the related look up table.
C             Here XX is an  array XX with length N.  Given a value X,
C             return a value Jlo such that X is between XX(JLO) and
C             XX(JLO+1). XX must be monotonic and  increase or decrease.
C            J=0 or J=N is return to inidicate that X is out of range.
C
C
C!Input Parameters:  solar angle (view zenith angle, relative azimuth angle,
C                    or aerosol optical depth, and the related
C                    look up table.
C
C!OUTPUT PARAMETERS: Index of the given solar angle (view zenith angle,
C                    relative azimuth angle, or aerosol optical depth)
C                    in the look up table
C
c!REVISION HISTORY:
c 09/07/99  Shaohua Shen
c fixed prolog.
c
c!TEAM-UNIQUE HEADER:
c
c Developed by MODIS Aerosol/Solar Water Vapor Retrieval Team
c GSFC, Greenbelt, MD
c
C !References and Credits:
C
C        Written by
C        Drs. Shaohua Shen and Bo-Cai Gao 09/07/99
C        shen@climate.gsfc.nasa.gov
C
C DESIGN NOTES: data and variables are passed with common statements
C               contained in 'MOD05_CORR.inc'.
C Internal subroutines called:
C              NONE
C
c!END-----------------------------------------------------------------------

	IMPLICIT NONE

	INTEGER jlo
	INTEGER n,inc,jhi,jm
	REAL x
	REAL xx(n)
	LOGICAL  ascnd

	ascnd = (xx(n) .gt. xx(1))
	if (jlo .le. 0 .or. jlo .gt. n) then
 	     jlo=1
	     jhi=n
             go to 3
        endif
	inc=1
	if (x .ge. xx(jlo) .eqv. ascnd) then
   1           jhi=jlo+inc
               if (jhi .gt. n) then
		  jhi=n
	       else if (x .ge. xx(jhi) .eqv. ascnd) then
		   jlo=jhi
		   inc=inc+inc
                   go to 1
	       end if
	else
	       jhi=jlo
   2           jlo=jhi-inc
	       if (jlo .lt. 1) then
	 	  jlo=1
	       else if (x .lt. xx(jlo) .eqv. ascnd) then
		  jhi=jlo
		  inc=inc+inc
                  go to 2
	       end if
	end if
  3     if (jhi-jlo .le. 1) then
   	    if (x .eq. xx(n)) jlo=n
	    if (x .eq. xx(1)) jlo=1
            return	
	else
   	    jm=(jhi+jlo)/2
	    if (x .ge. xx(jm) .eqv. ascnd) then
		jlo=jm
	    else
		jhi=jm
	    end if
	end if
        go to 3
        return
	END
C**************************************************************************
      SUBROUTINE READER_L2(MODFIL_MOD04,MODFIL_MOD05,MODFIL_L1B,Scan_No,
     &                     VAPOR,AEROSOL)

      IMPLICIT NONE
      Include 'PGS_MODIS_39500.f'
      Include 'mapi.inc'
      Include 'MOD05_CORR.inc'
C--------------------------------------------------------------------------
C
C !F77
C
C !DESCRIPTION:
C                       This subroutine reads MOD04 and MOD05 HDF
C                       arrays
C
C !INPUT PARAMETERS:
C
C         MODFIL_L1B    Value for reading L1B data
C         MODFIL_MOD04  Value for reading MOD04 data
C         MODFIL_MOD05  Value for reading MOD05 data
C         Scan_No       Scan number index (1- 200 typically)
C
C !OUTPUT PARAMETERS:
C
C         VAPOR         Total precipitable water (cm) retrieved from MOD05
C         AEROSOL       Aerosol Optical Thickness (0.67 micron) retrieved from MOD04
C
C !REVISION HISTORY:
C
C        WRITTEN BY
C        Dr. Allen Chu          11/25/97
C        Code 913
C        NASA Goddard Space Flight Center
C        Greenbelt, MD 20771
C
C !TEAM-UNIQUE HEADER:
C
C   This software was developed by the MODIS Atmosphere Science Team
C   for the National Aeronautics and Space Administration at
C   Goddard Space Flight Center.
C
C !END
C-----------------------------------------------------------------------

      CHARACTER*80 DATA_TYPE,ARRNM,GRPNM
      INTEGER RANK,NUMSCAN,RTN,Scan_No,IJ,IK,IS,L,I
      INTEGER MODFIL_L1B(MODFILLEN),MODFIL_MOD04(MODFILLEN),
     1        MODFIL_MOD05(MODFILLEN)
      INTEGER Dim_Size(3),START(3),EDGE_VAPOR(3),EDGE_AEROSOL(3)

      INTEGER*2 VAPOR_HDF(Max_Frames_Per_Line*No_Lines_Per_Scan)
      INTEGER*2 AEROSOL_HDF(Max_Grids_Per_Line*No_Grids_Per_Scan*3)

      REAL ADD_OFFSET,SCALE_FACTOR
      REAL VAPOR(Max_Frames_Per_Line,No_Lines_Per_Scan)
      REAL AEROSOL(Max_Grids_Per_Line,No_Grids_Per_Scan)

      EXTERNAL GMAR,GMARDM,GMFIN
      SAVE

C
C Call GMFIN to get 'Number of Scans' from L1B dataset (=200)
C

      Rank = 1
      data_type = 'INTEGER*4'
      arrnm= 'Number of Scans'

C      IF(Scan_No.EQ.1) THEN

      RTN = GMFIN(Modfil_L1B,arrnm,data_type,Rank,NUMSCAN)
      IF (rtn.ne.0) CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     &'MAPI function GMFIN failed',
     &'READER_L2 for L1B number of scans')

C      ENDIF

C
C Get 'Water_Vapor_Near_Infrared' from MOD05 HDF file
C

      rank=2
      data_type = 'INTEGER*4'
      arrnm='Water_Vapor_Near_Infrared'
      grpnm=' '

C
C Retrieve dimensions of sds data.
C

      If (GMARDM(MODFIL_MOD05,arrnm,grpnm,data_type,Rank,
     & Dim_Size).ne.MAPIOK)
     & CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,'GMARDM failed',
     & 'READER_L2')

      START(1)=0
      START(2)=(Scan_No-1)*10
      EDGE_VAPOR(1)=Dim_Size(1)
      EDGE_VAPOR(2)=Dim_Size(2)/NUMSCAN

      rtn=GMAR(MODFIL_MOD05,arrnm,grpnm,START,EDGE_VAPOR,VAPOR_HDF)

      IF(rtn.NE.MAPIOK)
     &  CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     &  'GMAR  failed for water','READER_L2')

C
C Save data in 2-D array
C

      add_offset=0.0
      scale_factor=1000.0

      L=0

      DO IK=1,EDGE_VAPOR(2)
      DO IS=1,EDGE_VAPOR(1)

        L=L+1
        VAPOR(IS,IK)=(REAL(VAPOR_HDF(L))-ADD_OFFSET)/SCALE_FACTOR

      ENDDO
      ENDDO

C
C Get aerosol optical thickness at ls micron from MOD04 HDF file,
C which is the 3rd dimension of the array 'Corrected_Optical_Depth_Land'
C

      rank=3
      data_type = 'INTEGER*4'
      arrnm='Corrected_Optical_Depth_Land'
      grpnm=' '

C
C RETRIEVE DIMENSIONS OF SDS DATA.
C


      If (GMARDM(MODFIL_MOD04,arrnm,grpnm,data_type,Rank,Dim_Size)
     & .ne.MAPIOK)
     & CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     &'GMARDM failed dim_info failed FOR aerosol optical thickness',
     &'READER_L2')

      START(1)=0
      START(2)=Scan_No-1
      START(3)=0
      EDGE_AEROSOL(1)=Dim_Size(1)
      EDGE_AEROSOL(2)=Dim_Size(2)/NUMSCAN
      EDGE_AEROSOL(3)=Dim_Size(3)

      rtn = GMAR(MODFIL_MOD04,arrnm,grpnm,START,EDGE_AEROSOL,
     &           AEROSOL_HDF)

      IF( rtn.NE.MAPIOK)
     & CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     & 'GMAR  failed for aerosol','READER_L2')

C
C Save data in 2-D array
C

      add_offset=0.0
      scale_factor=1000.0

      L=0

      DO IK=1,EDGE_AEROSOL(3)
      DO IJ=1,EDGE_AEROSOL(2)
      DO IS=1,EDGE_AEROSOL(1)
      
c  Pick up only for wavelength of 0.66 
        L=L+1
        IF(IK.EQ.EDGE_AEROSOL(3)) THEN
        AEROSOL(IS,IJ)=(REAL(AEROSOL_HDF(L))-ADD_OFFSET)/SCALE_FACTOR
        ENDIF

      ENDDO
      ENDDO
      ENDDO

      RETURN
      END


C*****************************************************************************
      SUBROUTINE MOD05_CORR_PutData(Modfil_MOD05,DataSize,Buf_Size1,
     &           Buf_Size2,Vapor_Corrected_Factor,Max_EV_Frame,
     &           LinesPerGranule)

      IMPLICIT NONE
      INCLUDE 'mapi.inc'
C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C
C   Parse and transfer a multiline data buffer as separate lines.
C
C !INPUT PARAMETERS:
C
C   INTEGER  Modfil_MOD05     Array containing HDF VS and SD ID numbers of
C                             MOD05 product file.
C   INTEGER  DataSize         Array specifying the dimension sizes of
C                             useful data within storage buffer
C   INTEGER  Buf_Size1/2      The dimension sizes of data buffer as
C                             allocated in calling routine
C   REAL     Vapor_Corrected_Factor  Buffer used for SDS "Water_Vapor_Correction_Factors"
C                             Buffer offset position to start of data line
C   INTEGER  Max_EV_Frame     Maximum value of earth view frame used in the
C			      level 1 B file
C   INTEGER  LinesPerGranule  Number of 1-km data lines in current granule
C
C !OUTPUT PARAMETERS:     NONE
C
C!Revision History:
C $Log: MOD05_IO.f,v $
c Revision 1.3  1997/11/25  13:29:50  DAC
c Baselined for Version 1.
c
c Revision 1.3  1996/11/01  19:33:43  jguu
c Baselined for Version 1.
c
c Revision 1.3  1996/11/01  19:33:43  jguu
c Remove some PRINT statements.
c
c Revision 1.2  1996/08/29  19:55:47  jguu
c The file name for the include file mod05_inc.f is changed to
c mod05.inc
c
c Revision 1.1  1996/07/08  18:32:08  vlin
c Initial revision
c
c Revision 1.2  1996/06/26  06:32:36  rhucek
c Added _FillValue attributes to MOD05 product file SDS arrays.
c
c Revision 1.1  1996/06/24  19:25:13  gao
c Initial revision
c
c Revision 1.2  1995/08/24  15:15:37  rhucek
c Updated for function return values
c
c Revision 1.1  1995/08/11  19:40:50  vlin
c Initial revision
C
C !TEAM-UNIQUE HEADER:
C
C      This software was developed by the MODIS Science Data Support
C      Team (SDST) for the National Aeronautics and Space Administration,
C      Goddard Space Flight Center, under contract NAS5-32373.
C
C !REFERENCES and CREDITS:
C
C   Written by Richard Hucek
C
C   Research and Data Systems Corporation
C   SAIC/GSC MODIS Science Data Support Office
C   7501 Forbes Blvd, Seabrook MD 20706
C
C   rhucek@ltpmail.gsfc.nasa.gov
C
C !DESIGN NOTES:
C
C Internals:
C
C   Functions and Subroutines:
C        MOD05_CORR_OUT
C
C   Variables:
C        iline       logical data line number
C        NPXLS       number of samples (pixels) in a data line
C        ioffset     buffer offset position at start of data line
C
C !END
C----------------------------------------------------------------------

      INTEGER Iline,NPXLS,Max_EV_Frame,LinesPerGranule,ioffset
      INTEGER Buf_Size1,Buf_Size2,DataSize(2),Modfil_MOD05(MODFILLEN)
      REAL    Vapor_Corrected_Factor(*)

      NPXLS = DataSize(1)

      Do Iline = 1, DataSize(2)

         ioffset = (iline-1)*Buf_Size1 + 1
         CALL MOD05_CORR_OUT(Modfil_MOD05,Vapor_Corrected_Factor(ioffset)
     2                       ,NPXLS,Max_EV_Frame,LinesPerGranule)

      End Do

      RETURN
      END


C**********************************************************************
      SUBROUTINE MOD05_CORR_OUT(MODFIL,DATA1,
     &                          NPXLS,Max_EV_Frame,NSCL)

      IMPLICIT NONE

      INCLUDE 'mapi.inc'
      INCLUDE 'MOD05_CORR.inc'
      INCLUDE 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C
C      Create and write Scientific Data Sets (SDS) for column water vapor
C      values above clear land surfaces and above clouds.  This sunroutine
C      will called every scan line but write data into HDF file at every
C      Lines_Per_Write.
C
C !INPUT PARAMETERS:
C
C   INTEGER  Max_EV_Frame    Maximum value of earth view frame used in the
C			     level 1 B file
C   INTEGER MODFIL(MODFILLEN)File handle structure for HDF files
C   REAL    DATA1*(NPXLS)    1-dimensional SDS "Column Water Vapor Correction Factor"
C   INTEGER NPXLS            Number of pixels.
C   INTEGER Max_EV_Frame     Maximum value of earth view frame used in the
C                            level 1 B file
C   INTEGER NSCL             Number of scanlines.
C
C !OUTPUT PARAMETERS:
C
C   NONE
C
C !REVISION HISTORY:
C Revision 2.0  1997/04/18  DAC
c
C        Take out the creation of HDF SDS structure, since the HDF file
C        has been created already create_mod05.
C
c Revision 1.2  1995/09/05  14:06:04  ding
c Updated the DimNames.
c
c Revision 1.1  1995/09/01  19:44:42  ding
c Initial revision
C
C !TEAM-UNIQUE HEADER:
C
C       This software is developed by the MODIS Science Data Support
C       Team for the National Aeronautics and Space Administration,
C       Goddard Space Flight Center, under contract NAS5-32373.
C
C !REFERENCES AND CREDITS
C
C       WRITTEN BY:
C       Xiao-Yang Ding                08/31/95
C       Research and Data systems Corporation
C       SAIC/GSC MODIS Science Data Support Office
C       7501 Forbes Blvd,  Seabrook MD 20706
C
C       ding@ltpmail.gsfc.nasa.gov
C
C !DESIGN NOTES:
C       This subroutine checks the return status of all MODIS Application
C       Program Interface (M-API) function calls.  A successful M-API
C       function call is indicated by a return value MAPIOK. If a M-API
C       function call is not successful, an error message is written
C       to the LogStatus file and the process aborts. This is achieved
C       by passing an fatal mnemonic error message (i.e. .._W_..)
C       to the function MODIS_SMF_SETDYNAMICMSG.
C
C  Externals:
C     Functions:
C        CRMAR        (libmapi.a)
C        PMAR         (libmapi.a)
C
C     Named Constants:
C        MODIS_W_GENERIC  (PGS_MODIS_39500.f)
C        MAPIOK           (mapi.inc)
C
C  Internals:
C     Functions and Subroutines:
C        DIMNAME      Write the dimension names of a SDS array.
C
C     Named Constants:
C        Glob_d*      Dimension name of glob metadata. (mod05_inc.f)
C        SDS*_N       Name of the SDS array. (mod05_inc.f)
C        Dim*_N       Dimension name for the SDS array. (mod05_inc.f)
C        Attr*_N      Attribute name for the SDS array. (mod05_inc.f)
C        DTYPE*       Data type for the SDS array. (mod05_inc.f)
C        Lines_Per_Write
C                     Scan line number for each PMAR call. (mod05_inc.f)
C
C     Variables:
C        START2(2)    Array specifying the starting location of each write.
C        Start(3)     Array specifying the starting location of EV data.
C        RANK2        The number of dimensions in an array.
C        Nmod         Remaindering of two integers.
C        Opmf_Ctr     Count of the times that this subroutine been called.
C        Buf(n),CLOUD(n),L_S_Flag(n)
C                     A temporary buffer to hold the subarray per write.
C                     n=1500*Lines_Per_Write
C
C !END
C-----------------------------------------------------------------------

C  Declaration

      INTEGER*2    Buf(1500*Lines_Per_Write)
      INTEGER      NPXLS,NSCL,RANK2,j,Opmf_Ctr,Max_EV_Frame,Nmod,l,
     &             MODFIL(MODFILLEN),DIMS(2),DIM2(2),START2(2)
      REAL DATA1(*)
      SAVE

      DATA RANK2/2/,Opmf_Ctr/0/,START2/2*0/

C  Initialization
C  DIMS specified the sizes of SDS array that CRMAR will create while
C  DIM2 specified the sizes of SDS array that PMAR will put in the
C  HDF file for each call.

      DIMS(1) = Max_EV_Frame
      DIMS(2) = NSCL
      DIM2(1) = NPXLS
      DIM2(2) = Lines_Per_Write
      Opmf_Ctr=Opmf_Ctr+1

C
C Calculate Buf for PMAR
C

      Nmod=Mod(Opmf_Ctr,Lines_Per_Write)
      if (Nmod.eq.0) Nmod=Lines_Per_Write
      do j=1,NPXLS
        l=(Nmod-1)*NPXLS+j
        Buf(l)=NINT(DATA1(j)*SCALE1+OFFSET1)
      enddo

CCCCCCCCCCCC    If Nmod equals to Lines_Per_Write   CCCCCCCCCCCCCCC

      IF (Nmod.eq.Lines_Per_Write) THEN

C PMAR places a multi-dimensional array of data into a HDF SDS array
C structure previously created using CRMAR.

          if (PMAR(MODFIL,SDS1_N,' ',START2,DIM2,Buf).ne.MAPIOK) call
     2    MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     3   'MAPI PMAR for SDS1 failed','MOD05_CORR_OUT')

          START2(2)=START2(2)+Lines_Per_Write

CCCCCCCCCCCC     Else if Opmf_Ctr equals to NSCL     CCCCCCCCCCCCCC

      ELSE IF (Opmf_Ctr.eq.NSCL) THEN

C Reset 2nd dimension for PMAR

          DIM2(2)=Nmod

C PMAR for the 1st SDS array

          if (PMAR(MODFIL,SDS1_N,' ',START2,DIM2,Buf).ne.MAPIOK) call
     2    MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     3   'MAPI PMAR for SDS1 failed','MOD05_CORR_OUT')

       ENDIF

      RETURN
      END
