      SUBROUTINE Process_Mod06CD(Error_Flag,idebug,SLOPE_MEAN_LAND,
     *                           SLOPE_MEAN_OCEAN)
C-----------------------------------------------------------------------
C !F77
C
C !Description: This program is to derive cirrus reflectance in the 0.4 - 1.0
C               micron region and contral index.
C
C !Input Parameters:
C
C   INTEGER     Refl_1         MODIS reflectance at channel 1
C   INTEGER     Refl_2         MODIS reflectance at channel 2
C   INTEGER     Refl_5         MODIS reflectance at channel 5
C   INTEGER     Refl_6         MODIS reflectance at channel 6
C   INTEGER     Refl_7         MODIS reflectance at channel 7
C   INTEGER     Refl_19        MODIS reflectance at channel 19
C   INTEGER     Refl_26        MODIS reflectance at channel 26
C   INTEGER     SolZen         Solar zenith angle
C   INTEGER     SatZen         Satellite view zenith angle
C   INTEGER     RelAZ          Relative azimuth angle between the Sun
C                              and the satellite
C   INTEGER     LandSea_Flag   index that indicates land when its value
C                              equals to 1 and ocean when 0
C   INTEGER     CldMsk         index that indicates clear pixel when its value
C                              equals to 1 and cloudy when 0.
C
C !Output Parameters:
C
C   INTEGER     Cirrus_Refl    Two dimensional (2-D) array for passing
C                              retrieved cirrus reflectance data
C   BYTE        I_QA_cirrus    Two dimensional (2-D) array for passing
C                              quality assurance data for cirrus reflectance
C                              and contrail index. Index 0 refers to bad data,
C                              index 2 refers to non-cirrus pixel, index 3 refers
C                              to cirrus pixel, and index 4 refers to contral pixel.
C  Logical      Error_Flag     Subroutine return status:  Error_Flag=.true. or
C               (not used)     .false. Currently assigned value .false.
C
C !Revision History:
C $Log: Process_mod06cd_V2.f,v $
C Modified by Liqun Ma    02/18/98
C Temporary output are moved out
C Prologs are updated
C Added flag for switching between testing binary files and MODIS data
C acoording to MODIS Science Team.
C
C !TEAM-UNIQUE HEADER:
C
C   This software was developed by the MODIS Science Team
C   (SDST) for the National Aeronautics and Space Administration,
C   Goddard Space Flight Center, under contract NAS5-31366.
C
C !REFERENCES and CREDITS:
C
C   Written by
C   Dr. Allen Chu                  November, 1997
C   SSAI
C   NASA GSFC Code 913
C   Greenbelt, MD 20771
C   email: achu@dustdevil.gsfc.nasa.gov
C   phone: 301-614-6237
C
C !DESIGN NOTES:
C
C   At present, this program is set up to process one granule of MODIS
C   sensor data. A granule consists of 200 MODIS scan swathes, each
C   containing 1354 1-km pixels in the scan direction and 10 pixels
C   along the spacecraft flight direction.
C
C Externals:
C   Named Constants:
C           MAPIOK                          (mapi.inc)
C           Max_Frames_Per_Line             (COMMONS_cirrus.inc)
C           No_Lines_Per_Scan               (COMMONS_cirrus.inc)
C
C Internals:
C
C   Functions and Subroutines:
C          FILEOC_MOD06CD
C          GetModisDat_Cirrus
C          READ_TABLE
C          MODIS_SMF_SETDYNAMICMSG
C          PROCESS_CIRRUS_REFLECTANCES
C
C !END
C-----------------------------------------------------------------------

      IMPLICIT NONE

      INCLUDE 'COMMONS_cirrus.inc'
      INCLUDE 'mapi.inc'
      INCLUDE 'PGS_MODIS_39500.f'
      INCLUDE 'PGS_PC.f'

C Liqun 02/18/98:  Added parameters FUNCNAME
C rhucek 01/04/99: Made variable PROCESS_FLAG a PARAMETER
C lzhou 08/18/00:  Changed FUNCNAME from 'Process_mod06cd_V2.f'
C                  to 'Process_Mod06CD'

      CHARACTER*(*) FUNCNAME
      PARAMETER    (FUNCNAME = 'Process_Mod06CD')

      CHARACTER*30  PROCESS_FLAG
      PARAMETER    (PROCESS_FLAG = 'MODIS')

      CHARACTER*5    ioflag
      CHARACTER*40   att_N,dtype
      CHARACTER*(PGSd_PC_FILE_PATH_MAX) FN_L1B

      INTEGER Max_EV_Frames

      CHARACTER*4  Scan_Type(Max_Scan_Per_Granule),ScanType
      CHARACTER*(PGSd_MET_GROUP_NAME_L-1) groups(PGSd_MET_NUM_OF_GROUPS)
      CHARACTER*(MAX_ECS_NAME_L-1) HDFAttNms(PGSd_MET_NUM_OF_GROUPS)

      INTEGER i, idebug, j, K, RTN, nms, Scans_Per_Granule,
     &        Total_No_Lines, BeginScan_No, EndScan_No, Scan_No,
     &        handle_QC,
     &        Buf_Size1, Buf_Size2,
     &        Modfil_Geo(MODFILLEN), Modfil_CldMsk(MODFILLEN),
     &        Modfil_L1B(MODFILLEN),
     &        handle_trans,Data_Size(2),
     &        Data_Size_For_Granule(2,Max_Scan_Per_Granule),
     &        No_Frames(Max_Scan_Per_Granule),
     &        bsize
      INTEGER GetModisDat_MOD06CD,IM,IN,IL

      INTEGER CldMsk      (Max_Frames_Per_Line, No_Lines_Per_Scan),
     2        LandSea_Flag(Max_Frames_Per_Line, No_Lines_Per_Scan),
     3        UFQ_Flag(Max_Frames_Per_Line, No_Lines_Per_Scan)

c Liqun Ma added logical data type
      LOGICAL Error_Flag

c-----Dr. Wei Han, han@han.nrl.navy.mil--------------------
      INTEGER  NumHandles
C
C What are the sample arrays ?
C
      BYTE    sample1(Max_Frames_Per_Line, No_Lines_Per_Scan),
     &        sample2(Max_Frames_Per_Line, No_Lines_Per_Scan),
     &        sample5(Max_Frames_Per_Line, No_Lines_Per_Scan),
     &        sample6(Max_Frames_Per_Line, No_Lines_Per_Scan),
     &        sample7(Max_Frames_Per_Line, No_Lines_Per_Scan)
      INTEGER Buf_cldmsk,Buf_cldmsk_QA
      PARAMETER(Buf_cldmsk=6,Buf_cldmsk_QA=10)
      BYTE Cloud(Buf_cldmsk,Max_Frames_Per_Line, No_Lines_Per_Scan),
     2 QA_Cloud(Buf_cldmsk_QA,Max_Frames_Per_Line, No_Lines_Per_Scan)


c-------end of han------------


C
C 1-km resolution data array (IR channels Rad_29, Rad_31, Rad_32,
C                             Height and Sfc_Temp are not used)
C

      REAL SatZen      (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     SolZen      (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     RelAz       (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Refl_1      (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Refl_2      (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Refl_5      (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Refl_6      (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Refl_7      (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Refl_19     (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Refl_26     (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Rad_29      (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Rad_31      (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Rad_32      (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Height      (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Sfc_Temp    (Max_Frames_Per_Line, No_Lines_Per_Scan)

C
C Uncertainty data
C

      BYTE UN_1        (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     UN_2        (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     UN_5        (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     UN_6        (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     UN_7        (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     UN_19       (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     UN_26       (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     UN_29       (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     UN_31       (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     UN_32       (Max_Frames_Per_Line, No_Lines_Per_Scan)

C
C Validity flag
C

      BYTE Vflag_1     (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Vflag_2     (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Vflag_5     (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Vflag_6     (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Vflag_7     (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Vflag_19    (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Vflag_26    (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Vflag_29    (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Vflag_31    (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Vflag_32    (Max_Frames_Per_Line, No_Lines_Per_Scan)

C
C ????? Uncertainty and saturation
C

      BYTE Buf_Un      (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Buf_Sa      (Max_Frames_Per_Line, No_Lines_Per_Scan)

      SAVE

C
C Error flag is not used right now and planned for future use
C

      Error_Flag   = .false.

C
C Initialize buffers to hold HDF VS and SD file Ids
C

      DO 20 i = 1, MODFILLEN

         modfil_L1B(i) = 0
         modfil_Geo(i) = 0
         modfil_CldMsk(i) = 0

   20 CONTINUE

C
C Initialize imaging arrays and parameters for image processing
C

      DO J = 1, Max_Lines_Per_Img
      DO I = 1, Max_Frames_Per_Line

         REFL_0P67_Img(I,J)    = 0.0
         REFL_0P86_Img(I,J)    = 0.0
         REFL_0P94_Img(I,J)    = 0.0
         REFL_1P24_Img(I,J)    = 0.0
         REFL_1P38_Img(I,J)    = 0.0
         REFL_1P64_Img(I,J)    = 0.0
         REFL_2P13_Img(I,J)    = 0.0
         ZENITH_SOL_Img(I,J)   = 0.0
         ZENITH_VIEW_Img(I,J)  = 0.0
         AZIMUTH_REL_Img(I,J)  = 0.0
         I_LAND_WATER_Img(I,J) = 0     !land/sea flag
         I_CLOUD_Img(I,J)      = 0     !cloud mask flag

      END DO
      END DO

C
C Files open for L1B, geolocation, and cloud mask (
C hanele_QC, and MOD06cd file should not be used)
C

      ioflag = 'OPEN'
      CALL FILEOC_MOD06CD(ioflag,handle_trans,handle_QC,FN_L1B,
     &            modfil_L1B,modfil_Geo,modfil_CldMsk,
     &            groups,HDFAttNms,NumHandles)

C
C Get metadata info (Scans_Per_Granule, Max_EV_Frames, No_Frames,
C                    and Scan_Type) from L1B data
C

      att_N = 'Number of Scans'
      dtype = 'INTEGER*4'
      nms = 1
      RTN = GMFIN(modfil_L1B,att_N,dtype,nms,Scans_Per_Granule)
      IF (rtn.ne.mapiok) CALL MODIS_SMF_SETDYNAMICMSG
     &(MODIS_W_GENERIC,'MAPI function 1st GMFIN failed', FUNCNAME)

      att_N = 'Max Earth View Frames'
      dtype = 'INTEGER*4'
      nms = 1
      RTN = GMFIN(modfil_L1B,att_N,dtype,nms,Max_EV_Frames)
      IF (rtn.ne.mapiok) CALL MODIS_SMF_SETDYNAMICMSG
     &(MODIS_W_GENERIC,'MAPI function 2nd GMFIN failed', FUNCNAME)

C
C No_Frames and Scan_Type are a function of scans
C

      bsize = Max_Scan_Per_Granule * 4
      RTN = GMTBL(modfil_L1B,'Level 1B Swath Metadata',' ',
     &      'EV_Frames',0,Scans_Per_Granule,bsize,No_Frames)
      IF (rtn.ne.mapiok) CALL MODIS_SMF_SETDYNAMICMSG
     &(MODIS_W_GENERIC,'MAPI function 1st GMTBL failed', FUNCNAME)

      bsize = Max_Scan_Per_Granule*4
      rtn = GMTBL(modfil_L1B,'Level 1B Swath Metadata',' ',
     &      'Scan Type',0,Scans_Per_Granule,bsize,Scan_Type)
      IF (rtn.ne.mapiok) CALL MODIS_SMF_SETDYNAMICMSG
     &(MODIS_W_GENERIC,'MAPI function 2nd GMTBL failed', FUNCNAME)

C
C We adopt the convention used in image processing with
C     x-dimension = No. of Samples, y-dimension = No. of Lines.
C Note: Because each scan in MODIS L1B data may have different
C       number of pixels in the x-dimension, the resulting images
C       from a MODIS data cube may be zig-zagged. We assume that
C       the x-dimension of images passed into our subroutines for cirrus
C       reflectance and contrail index derivations is the the
C       maximum of the allocated size. For the parts where MODIS has no
C       measured data, our assumed initial values will be used during
C       the derivations.
C

      Total_No_Lines = Scans_Per_Granule * 10

C
C Buf_Size1 and Buf_Size2 are used in getting MODIS data
C

      Buf_Size1 = Max_Frames_Per_Line
      Buf_Size2 = No_Lines_Per_Scan

C
C No_Of_Samples and No_Of_Lines are used to process the whole
C granule of data
C

      No_Of_Samples = Max_Frames_Per_Line
      No_Of_Lines   = Total_No_Lines

C
C L_Index_Img is an index used in the "Do 999" loop to build the
C whole image (typically 1354 x 2000 pixels for 1 km resolution)
C from individual MODIS scans (one MODIS scan has 10 lines in 1-km
C resolution).
C

      L_Index_Img = 0

C
C Begin reading MODIS data and ancillary data in HDF-EOS format
C

      BeginScan_No = 1
      EndScan_No   = Scans_Per_Granule

C
C Loop over the number of scan swathes.  For each swath, the MODIS
C inputs include sensor reflectance, geolocation, and cloud mask.
C The MODIS daytime measurements are only used.
C

      DO 999 Scan_No = BeginScan_No, EndScan_No

         ScanType=Scan_Type(Scan_No)

         IF (ScanType(1:1).eq.'D') THEN

C compared to IF(ScanType.eq.'D') THEN

            RTN=GetModisDat_MOD06CD(
     &          Modfil_L1B,Modfil_Geo,Modfil_CldMsk,FN_L1B,
     &          Scan_No,Buf_Size1,Buf_Size2,Data_Size,
     &          SatZen,SolZen,RelAz,Height,
     &          Refl_1, Refl_2, Refl_5, Refl_6, Refl_7, Refl_19, Refl_26, Rad_29,  Rad_31,  Rad_32,
     &          Un_1,   Un_2,   Un_5,   Un_6,   Un_7,   Un_19,   Un_26,   Un_29,   Un_31,   Un_32,
     &          Vflag_1,Vflag_2,Vflag_5,Vflag_6,Vflag_7,Vflag_19,Vflag_26,Vflag_29,Vflag_31,Vflag_32,
     &          Buf_Un,Buf_Sa,
     &          CldMsk,LandSea_Flag, Cloud,QA_Cloud, sample1,
     &          sample2, sample5, sample6, sample7,UFQ_Flag)

            IF (rtn.ne.0) CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     &                    'Call to GetModisDat_cirrus failed',FUNCNAME)

            IF (idebug .EQ. 1) THEN

C              CALL Write_QC(handle_QC,Scan_No,Refl_1,Refl_2,Refl_5,
C     &                      Refl_6,Refl_7,Refl_19,Refl_26,SatZen,
C     &                      SolZen,RelAZ,CldMsk,LandSea_Flag,Height)

            ENDIF

C
C What are the use of Data_Size array ???????
C

            Data_Size(1) = No_Frames(Scan_No)
            Data_Size(2) = No_Lines_Per_Scan

            DO J = 1, No_Lines_Per_Scan

               L_Index_Img = (Scan_No - 1)*No_Lines_Per_Scan + J

               DO I = 1, No_Frames(Scan_No)

C
C Temp Code [SHOULD BE REMOVED AFTER TESTS]
C---A flag for switching between using test binary images and using MODIS data
C---If PROCESS_FLAG = 'TEST', the testing binary image file is used as input.
C---Otherwise, MODIS data will be used.  W. Han, 02/17/98
C
               IF (PROCESS_FLAG .EQ. 'TEST') THEN

                  Refl_1(i,j)       =  0.9
                  Refl_5(I,J)       =  0.85
                  IF(J.LE.3)  REFL_26(I,J) =  0.2
                  IF(J.GT.3.and.J.LE.6)  REFL_26(I,J) =  0.000002
                  IF(J.GT.6)  REFL_26(I,J) =  BAD_FLAG
                  SolZen(i,j)       =  60.
                  SatZen(i,j)       =  60.
                  RelAZ(i,j)        =  60.
                  LandSea_Flag(i,j) =  1
                  CldMsk(i,j)       =  1

               END IF

C End of Temp Code ---

               REFL_0P67_Img(I,L_Index_Img)    = Refl_1(i,j)
               REFL_0P86_Img(I,L_Index_Img)    = Refl_2(i,j)
               REFL_1P24_Img(I,L_Index_Img)    = Refl_5(i,j)
               REFL_1P64_Img(I,L_Index_Img)    = Refl_6(i,j)
               REFL_2P13_Img(I,L_Index_Img)    = Refl_7(i,j)
               REFL_1P38_Img(I,L_Index_Img)    = Refl_26(i,j)
               REFL_0P94_Img(I,L_Index_Img)    = Refl_19(i,j)
               ZENITH_SOL_Img(I,L_Index_Img)   = SolZen(i,j)
               ZENITH_VIEW_Img(I,L_Index_Img)  = SatZen(i,j)
               AZIMUTH_REL_Img(I,L_Index_Img)  = RelAZ(i,j)
               I_LAND_WATER_Img(I,L_Index_Img) = LandSea_Flag(i,j)
               I_CLOUD_Img(I,L_Index_Img)      = CldMsk(i,j)

               END DO
            END DO

         ELSE

            Data_Size(1) = 10
            Data_Size(2) = No_Lines_Per_Scan

         END IF

            DO K=1,2
               Data_Size_For_Granule(K,Scan_No)=Data_Size(K)
            ENDDO

 999  CONTINUE

C
C Read in a pre-calculated lookup table. This table contains transmittance
C of the MODIS 1.375-micron channel. The transmittances are used in adjusting
C the transmittance factors due to the solar and viewing geometry.
C

        CALL READ_TABLE(handle_trans)

C
C Process MODIS data cube to get cirrus reflectance (path radiance)
C

        IF(I_PROCESS_MODE.EQ.1) CALL PR_CIRRUS_REFLECTANCES_MODE_1

        IF(I_PROCESS_MODE.EQ.2) THEN

C
C Multiply reflectance values (0-1.0) by 1000 to derive histogram
C

           DO J = 1, No_Of_Lines
              DO I = 1, No_Of_Samples
                 REFL_0P67_Img(I,J) = REFL_0P67_Img(I,J) * 1000.
                 REFL_0P86_Img(I,J) = REFL_0P86_Img(I,J) * 1000.
                 REFL_0P94_Img(I,J) = REFL_0P94_Img(I,J) * 1000.
                 REFL_1P24_Img(I,J) = REFL_1P24_Img(I,J) * 1000.
                 REFL_1P64_Img(I,J) = REFL_1P64_Img(I,J) * 1000.
                 REFL_2P13_Img(I,J) = REFL_2P13_Img(I,J) * 1000.
                 REFL_1P38_Img(I,J) = REFL_1P38_Img(I,J) * 1000.
              END DO
           END DO

           CALL PR_CIRRUS_REFLECTANCES_MODE_2(Scans_Per_Granule,
     *     No_Frames,SLOPE_IMG,INDEX_OCELAND)

C
C onvert back to reflectance units (0 - 1.0) by dividing 1000.
C

           DO J = 1, No_Of_Lines
              DO I = 1, No_Of_Samples
                 cirrus_refl(I,J) = cirrus_refl(I,J)/1000.
              END DO
           END DO

        END IF

C
C Calculate mean of slopes separating for ocean and land within the granule
C

      DO IL=1,3
        SLOPE_MEAN_LAND(IL)=0.0
        SLOPE_N_LAND(IL)=0.0
        SLOPE_MEAN_OCEAN(IL)=0.0
        SLOPE_N_OCEAN(IL)=0.0
        DO IM=1,7
          DO IN=1,5
            IF(SLOPE_IMG(IN,IM,IL).GT.0.0.AND.INDEX_OCELAND(IN,IM,IL).EQ.0) THEN
              SLOPE_N_LAND(IL)=SLOPE_N_LAND(IL)+1.0
              SLOPE_MEAN_LAND(IL)=SLOPE_MEAN_LAND(IL)+SLOPE_IMG(IN,IM,IL)
            ELSE IF(SLOPE_IMG(IN,IM,IL).GT.0.0.AND.INDEX_OCELAND(IN,IM,IL).EQ.1) THEN
              SLOPE_N_OCEAN(IL)=SLOPE_N_OCEAN(IL)+1.0
              SLOPE_MEAN_OCEAN(IL)=SLOPE_MEAN_OCEAN(IL)+SLOPE_IMG(IN,IM,IL)
            ENDIF
          ENDDO
        ENDDO
        IF(SLOPE_N_LAND(IL).GT.0.0) THEN
          SLOPE_MEAN_LAND(IL)=SLOPE_MEAN_LAND(IL)/SLOPE_N_LAND(IL)
          SLOPE_MEAN_LAND(IL)=1.0/SLOPE_MEAN_LAND(IL)
        ENDIF
        IF(SLOPE_N_OCEAN(IL).GT.0.0) THEN
          SLOPE_MEAN_OCEAN(IL)=SLOPE_MEAN_OCEAN(IL)/SLOPE_N_OCEAN(IL)
          SLOPE_MEAN_OCEAN(IL)=1.0/SLOPE_MEAN_OCEAN(IL)
        ENDIF
      ENDDO

C---End of temp coding


      RETURN
      END


C***********************************************************************
      SUBROUTINE READ_TABLE(handle_trans)

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION: This subroutine is to read geometrical and transmittance
C              factors for the MODIS 1.375-micron channel. These factors
C              are used in deriving cirrus reflectances by assuming
C              that the 1.375-micron channel transmittance for overhead
C              Sun and nadir looking geometry is known.
C
C              This subroutine also read in a pre-calculated
C              transmittance table for water vapor in six channels:
C              0.86, 0.905, 0.935, 0.94, 1.24, and 1.38 micron. These
C              quantities are used in deriving upper level water vapor
C              transmittances from imaging data themselves.
C
C !INPUT PARAMETERS:
C             None
C
C !OUTPUT PARAMETERS:
C   FLOAT     GH2O          Geometrical factors related to solar and view
C                           zenith angles.
C   FLOAT     T_1P38_0P50   Transmittance factors for MODIS' 1.375-micron
C                           channel correspond to a 2-way vertical
C                           transmittance of 0.50
C   FLOAT     T_1P38_0P55   Transmittance factors for MODIS' 1.375-micron
C                           channel correspond to a 2-way vertical
C                           transmittance of 0.55
C   FLOAT     T_1P38_0P60   Transmittance factors for MODIS' 1.375-micron
C                           channel correspond to a 2-way vertical
C                           transmittance of 0.60
C   FLOAT     T_1P38_0P65   Transmittance factors for MODIS' 1.375-micron
C                           channel correspond to a 2-way vertical
C                           transmittance of 0.65
C   FLOAT     T_1P38_0P70   Transmittance factors for MODIS' 1.375-micron
C                           channel correspond to a 2-way vertical
C                           transmittance of 0.70
C   FLOAT     T_1P38_0P75   Transmittance factors for MODIS' 1.375-micron
C                           channel correspond to a 2-way vertical
C                           transmittance of 0.75
C   FLOAT     T_1P38_0P80   Transmittance factors for MODIS' 1.375-micron
C                           channel correspond to a 2-way vertical
C                           transmittance of 0.80
C   FLOAT     VAPOR         A total of 220 Water vapor values used in
C                           building the transmittance table for 6 channels
C   FLOAT     TRANSM_0p86   Transmittances for MODIS' 0.865-micron channel
C                           correspond to the 220 water vapor values
C   FLOAT     TRANSM_0p905  Transmittances for MODIS' 0.905-micron channel
C                           correspond to the 220 water vapor values
C   FLOAT     TRANSM_0p935  Transmittances for MODIS' 0.935-micron channel
C                           correspond to the 220 water vapor values
C   FLOAT     TRANSM_0p94   Transmittances for MODIS' 0.94-micron channel
C                           correspond to the 220 water vapor values
C   FLOAT     TRANSM_1p24   Transmittances for MODIS' 1.24-micron channel
C                           correspond to the 220 water vapor values
C   FLOAT     TRANSM_1p38   Transmittances for MODIS' 1.38-micron channel
C                           correspond to the 220 water vapor values
C
C !REVISION HISTORY:
C $Log: Process_mod06cd_V2.f,v $
C
C !TEAM-UNIQUE HEADER: - NOT YET DEFINED
C
C !REFERENCES AND CREDITS
C    Written by Dr. Bo-Cai Gao    09/16/97
C    gao@neptune.nrl.navy.mil
C
C !DESIGN NOTES: data and variables are passed with common statements
C               contained in "COMMONS_cirrus.inc"
C
C !END
C----------------------------------------------------------------------

      IMPLICIT NONE
      SAVE

      INCLUDE 'COMMONS_cirrus.inc'
      INTEGER handle_trans,I


C---Read in pre-calculated geometrical and transmittance factors
C        for the MODIS 1.375-micron channel from a table, which was
C        obtained by concatenating six transmittance tables that
C        correspond to 2-way vertical transmittances of 0.50, 0.55, 0.60,
C        0.65, 0.70, 0.75, and 0.80.
C
c-----Dr. Wei Han, han@han.nrl.navy.mil--------------------

C----Read the table corresponding to 2-way vertical transmittance of 0.50
      READ(handle_trans,*)
      DO I = 1, N_GH2O
         READ(handle_trans,*) GH2O(I), T_1P38_0P50(I)
      END DO

C----Read the table corresponding to 2-way vertical transmittance of 0.55
      READ(handle_trans,*)
      DO I = 1, N_GH2O
         READ(handle_trans,*) GH2O(I), T_1P38_0P55(I)
      END DO

C----Read the table corresponding to 2-way vertical transmittance of 0.60
      READ(handle_trans,*)
      DO I = 1, N_GH2O
         READ(handle_trans,*) GH2O(I), T_1P38_0P60(I)
      END DO

C----Read the table corresponding to 2-way vertical transmittance of 0.65
      READ(handle_trans,*)
      DO I = 1, N_GH2O
         READ(handle_trans,*) GH2O(I), T_1P38_0P65(I)
      END DO

C----Read the table corresponding to 2-way vertical transmittance of 0.70
      READ(handle_trans,*)
      DO I = 1, N_GH2O
         READ(handle_trans,*) GH2O(I), T_1P38_0P70(I)
      END DO

C----Read the table corresponding to 2-way vertical transmittance of 0.75
      READ(handle_trans,*)
      DO I = 1, N_GH2O
         READ(handle_trans,*) GH2O(I), T_1P38_0P75(I)
      END DO

C----Read the table corresponding to 2-way vertical transmittance of 0.80
      READ(handle_trans,*)
      DO I = 1, N_GH2O
         READ(handle_trans,*) GH2O(I), T_1P38_0P80(I)
      END DO

C---Read the transmittance table for water vapor in six channels: 0.86,
C            0.905, 0.935, 0.94, 1.24, and 1.38 micron
      READ(handle_trans,*)
      DO I = 1, N_TRANSM_TABLE
         READ(handle_trans,*) VAPOR(I), TRANSM_0p86(I),
     $         TRANSM_0p905(I), TRANSM_0p935(I),
     $         TRANSM_0p94(I),  TRANSM_1p24(I),  TRANSM_1p38(I)
      END DO

      RETURN
      END


C**************************************************************************
      SUBROUTINE PR_CIRRUS_REFLECTANCES_MODE_1

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION: This subroutine derives cirrus reflectances in the
C              0.4 - 1.0 micron region using the special 1.375-micron
C              MODIS cirrus detection channel.
C
C !Input Parameters:
C
C   INTEGER     Refl_1         MODIS reflectance at channel 1
C   INTEGER     Refl_5         MODIS reflectance at channel 5
C   INTEGER     Refl_26        MODIS reflectance at channel 26
C   INTEGER     SolZen         Solar zenith angle
C   INTEGER     SatZen         Satellite view zenith angle
C   INTEGER     RelAZ          Relative azimuth angle between the Sun
C                              and the satellite
C
C !Output Parameters:
C
C   INTEGER     Cirrus_Refl    Two dimensional (2-D) array for passing
C                               retrieved cirrus reflectance data
C
C !REVISION HISTORY:
C $Log: Process_mod06cd_V2.f,v $
C
C !TEAM-UNIQUE HEADER: - NOT YET DEFINED
C
C !REFERENCES AND CREDITS
C    First version created by Dr. Bo-Cai Gao
C    (gao@neptune.nrl.navy.mil)
C                             1997 09 16
C
C !DESIGN NOTES: data and variables are passed with common statements
C               contained in "COMMONS_cirrus.inc"
C
C !Functions and Subroutines:
C      HUNT_MOD06CD
C      FIND_CIRRUS_REFLECTANCE
C
C
C---Note: The at launch version of the cirrus reflectance algorithm
C         is very simplified. The post-launch version of the algorithm
C         will be improved after the analysis of real MODIS data.
C
C !END
C-----------------------------------------------------------------------

      IMPLICIT NONE
      SAVE

      INCLUDE 'COMMONS_cirrus.inc'
      INTEGER No_Frames(Max_Scan_Per_Granule),Scans_Per_Granule

c     INTEGER I, J, K
c     REAL    PI, VAP, DEG_TO_RAD
c     DATA PI, DEG_TO_RAD, J_GH2O /3.1415926, 0.0174533, 1/

      INTEGER I,J,J_GH2O
      REAL DEG_TO_RAD,GEOM_FACTOR,T_1P38_PIXEL
      DATA DEG_TO_RAD, J_GH2O /0.0174533, 1/

      DO J = 1, No_Of_Lines
         DO I = 1, No_Of_Samples

C
C A geometrical factor related to the solar and view zenith angles
C
        GEOM_FACTOR = 1. / COS(ZENITH_SOL_Img(I,J)  * DEG_TO_RAD)
     &              + 1. / COS(ZENITH_VIEW_Img(I,J) * DEG_TO_RAD)
C

       CALL HUNT_MOD06CD(GH2O, N_GH2O, GEOM_FACTOR, J_GH2O)
       CALL FIND_CIRRUS_REFLECTANCE(GH2O, N_GH2O, GEOM_FACTOR,
     &           J_GH2O, T_1P38_0P70, T_1P38_PIXEL)

       CIRRUS_REFL(I,J) = REFL_1P38_Img(I,J) / T_1P38_PIXEL

       IF(CIRRUS_REFL(I,J).GT.REFL_1P24_Img(I,J))
     &                 CIRRUS_REFL(I,J)   =    REFL_1P24_Img(I,J)

C
C Taking care of BAD_FLAG for REFL_1P38_Img
C

       IF(IFIX(REFL_1P38_Img(I,J)).EQ.IFIX(BAD_FLAG))
     &                 CIRRUS_REFL(I,J)   =    0.0

        END DO
      END DO


      RETURN
      END


C***********************************************************************
      SUBROUTINE FIND_CIRRUS_REFLECTANCE(X, N, G, J, T_1P38, T)

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION: A linear interpolation routine to derive a transmittance
C          value. Specifically: giving an array X with N elements in
C          ascending or decending order, G a value satisfying the
C          relation X(J) < G < X(J+1), doing a linear interpolation using
C          two values in array T_1P38 (also with N elements) to find T.
C
C !Input Parameters:
C   INTEGER     N         Number of elements in array X.
C   FLOAT       X         An array with N elements in ascending or decending
C                            order.
C   INTEGER     J         An index.
C   FLOAT       G         A value satisfying the relation X(J) < G < X(J+1).
C   FLOAT       T_1P38    An array containing 1.38-micron channel transmitances
C
C !Output Parameters:
C   FLOAT       T         Transmittance factor for the 1.38-micron channel for
C                            one pixel.
C
C !REVISION HISTORY:
C $Log: Process_mod06cd_V2.f,v $
C     Explicit data type.
C
C !TEAM-UNIQUE HEADER: - NOT YET DEFINED
C
C !REFERENCES AND CREDITS
C     WRITTEN BY: Bo-Cai Gao
C
C !DESIGN NOTES: none
C
C !END
C-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER J, N
      REAL DLT, FJ, FJ_P1, G, T
      REAL X, T_1P38
      DIMENSION X(N), T_1P38(N)

      IF( (J.GT.0).AND.(J.LT.N) ) THEN

         DLT   = X(J+1) - X(J)
         FJ    = ( X(J+1) - G ) / DLT
         FJ_P1 = ( G - X(J) )   / DLT

         T     = FJ * T_1P38(J) + FJ_P1 * T_1P38(J+1)

      ELSE

         IF(J.LE.0) T = T_1P38(1)
         IF(J.GE.N) T = T_1P38(N)

      END IF


      RETURN
      END


C*******************************************************************************
      SUBROUTINE HUNT_MOD06CD(XX,N,X,JLO)

C-------------------------------------------------------------------------------
C !F77
C
C !DESCRIPTION: finds the element in array XX that is closest to value X.
C           Array XX must be monotonic, either increasing or decreasing.
C
C !INPUT PARAMETERS:
C   INTEGER     N         Number of elements in array XX.
C   FLOAT       XX        An array with N elements in ascending or decending
C                            order.
C   FLOAT       X         A number.
C
C !OUTPUT PARAMETERS:
C   INTEGER     JLO       An index of the closest matching element.
C
C !REVISION HISTORY:
C $Log: Process_mod06cd_V2.f,v $
C
C     Explicit data type.
C
C !TEAM-UNIQUE HEADER: - NOT YET DEFINED
C
C !REFERENCES AND CREDITS
C     This subroutine was copied from Numerical Recipes                *
C
C !DESIGN NOTES: None
C
C !END
C-------------------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER jlo,n
      REAL x,xx(n)
      INTEGER inc,jhi,jm
      LOGICAL ascnd
C---
      REAL ddx,sss
      data sss/1.0E-12/
      ddx = abs( xx(n) - xx(1) )
C---
      ascnd=xx(n).ge.xx(1)
      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          goto 2
        endif
      endif
3     if(jhi-jlo.eq.1)then
C---        if(x.eq.xx(n))jlo=n-1
      if(abs(x - xx(n)) .lt. ddx*sss ) jlo=n-1
C---        if(x.eq.xx(1))jlo=1
      if(abs(x - xx(1)) .lt. ddx*sss ) jlo=1
        return
      endif
      jm=(jhi+jlo)/2
      if(x.ge.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3

      END


C********************************************************************************
      SUBROUTINE Write_QC(handle_QC,Scan_No,Refl_1,Refl_2,Refl_5,
     &                    Refl_6,Refl_7,Refl_19,Refl_26,SatZen,
     &                    SolZen,RelAZ,CldMsk,LandSea_Flag,Height)

C-----------------------------------------------------------------------
C !F77
C
C !Description: This subroutine writes cirrus qc date.
C
C !Input Parameters:
C
C   INTEGER     Refl_1         MODIS reflectance at channel 1
C   INTEGER     Refl_5         MODIS reflectance at channel 5
C   INTEGER     Refl_26        MODIS reflectance at channel 26
C   INTEGER     SolZen         Solar zenith angle
C   INTEGER     SatZen         Satellite view zenith angle
C   INTEGER     RelAZ          Relative azimuth angle between the Sun
C                              and the satellite
C   INTEGER     LandSea_Flag   index that indicates land when its value
C                              equals to 1 and ocean when 0
C   INTEGER     CldMsk         index that indicates clear pixel when its value
C                              equals to 1 and cloudy when 0.
C
C !Output Parameters: None
C
C
C !Revision History:
C $Log: Process_mod06cd_V2.f,v $
C  Iniatial version
C
C !TEAM-UNIQUE HEADER:
C
C   This software was developed by the MODIS Science Team
C   (SDST) for the National Aeronautics and Space Administration,
C   Goddard Space Flight Center, under contract NAS5-31366.
C
C !REFERENCES and CREDITS:
C
C   Written by
C   Liqun Ma
C
C !DESIGN NOTES:
C
C
C Externals:
C
C Internals:
C
C
C !END
C-----------------------------------------------------------------------
C
      IMPLICIT NONE
      SAVE

      INCLUDE 'COMMONS_cirrus.inc'

      INTEGER Scan_No,handle_QC,i,j

C
C 1-km resolution data array
C

      INTEGER LandSea_Flag(Max_Frames_Per_Line, No_Lines_Per_Scan)
      INTEGER CldMsk(Max_Frames_Per_Line, No_Lines_Per_Scan)

      REAL SatZen      (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     SolZen      (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     RelAz       (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Refl_1      (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Refl_2      (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Refl_5      (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Refl_6      (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Refl_7      (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Refl_19     (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Refl_26     (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Rad_29      (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Rad_31      (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Rad_32      (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Height      (Max_Frames_Per_Line, No_Lines_Per_Scan)


            write(handle_QC,'(/,"Scan Number ",i3)') Scan_No


*/  Modified by JC Guu 07/25/96
*/  The delimiters for character strings are changed from
*/  quotation marks (") to apostrophes (') to comply to
*/  the ANSI FORTRAN 77 standard.  There are 26
*/  locations where this correction occurrs.


            write(handle_QC,*)
            write(handle_QC,'(A)') 'first 11 reflectances for ' //
     &                            'band 1 for 10 lines'
            write(handle_QC,'(11f11.5)') ((Refl_1(i,j),i=1,11),
     &                                        j=1,10)

            write(handle_QC,*)
            write(handle_QC,'(A)') 'first 11 reflectance for band 2 ' //
     &                                'for 10 lines'
            write(handle_QC,'(11f11.5)') ((Refl_2(i,j),i=1,11),j=1,10)


            write(handle_QC,*)
            write(handle_QC,'(A)') 'first 11 reflectance for band 5 ' //
     &                                'for 10 lines'
            write(handle_QC,'(11f11.5)') ((Refl_5(i,j),i=1,11),j=1,10)

            write(handle_QC,*)
            write(handle_QC,'(A)') 'first 11 reflectance for band 6 ' //
     &                                'for 10 lines'
            write(handle_QC,'(11f11.5)') ((Refl_6(i,j),i=1,11),j=1,10)

            write(handle_QC,*)
            write(handle_QC,'(A)') 'first 11 reflectance for band 7 ' //
     &                                'for 10 lines'
            write(handle_QC,'(11f11.5)') ((Refl_7(i,j),i=1,11),j=1,10)

            write(handle_QC,*)
            write(handle_QC,'(A)') 'first 11 reflectances for ' //
     &                             'band 19 for 10 lines'
            write(handle_QC,'(11f11.5)') ((Refl_19(i,j),i=1,11),j=1,10)

            write(handle_QC,*)
            write(handle_QC,'(A)') 'first 11 reflectances for ' //
     &                             'band 26 for 10 lines'
            write(handle_QC,'(11f11.5)') ((Refl_26(i,j),i=1,11),j=1,10)

            write(handle_QC,*)
            write(handle_QC,'(A)') 'first 11 sat. zenith angles ' //
     &                             'for 10 lines'
            write(handle_QC,'(11f11.5)') ((SatZen(i,j),i=1,11),j=1,10)

            write(handle_QC,*)
            write(handle_QC,'(A)') 'first 11 solar zenith angles ' //
     &                             'for 10 lines'
            write(handle_QC,'(11f11.5)') ((SolZen(i,j),i=1,11),j=1,10)

            write(handle_QC,*)
            write(handle_QC,'(A)') 'first 11 rel. azimuth angles ' //
     &                             'for 10 lines'
            write(handle_QC,'(11f11.5)') ((RelAz(i,j),i=1,11),j=1,10)

            write(handle_QC,*)
            write(handle_QC,'(A)') 'first 11 Cloud Mask Decisions ' //
     &                             'for 10 lines'
            write(handle_QC,'(10I4)') ((CldMsk(i,j),i=1,11),j=1,10)

            write(handle_QC,*)
            write(handle_QC,'(A)') 'first 11 Land/Sea Flags for ' //
     &                             '10 lines'
            write(handle_QC,'(11I4)') ((LandSea_Flag(i,j),i=1,11),
     &                                  j=1,10)

            write(handle_QC,*)
            write(handle_QC,'(A)') 'first 11 Heights for 10 lines'
            write(handle_QC,'(11F7.1)') ((Height(i,j),i=1,11),j=1,10)

            write(handle_QC,*)
C---            write(handle_QC,'(A)') 'first 11 Sfc_Temps for 10 lines'
C---            write(handle_QC,'(11F7.2)') ((sfc_temp(i,j),i=1,11),j=1,10)

      RETURN
      END



C***********************************************************************
       SUBROUTINE PR_CIRRUS_REFLECTANCES_MODE_2(Scans_Per_Granule,
     * No_Frames,SLOPE_IMG,INDEX_OCELAND)

C-----------------------------------------------------------------------
C !F77
C
C !Description: This program can be used to replace the subroutine
C             PR_CIRRUS_REFLECTANCES_MODE_1 in MOD_cirrus_V2.f,
C             which derives cirrus reflectance in the 0.4 - 1.0 micron
C             region for MODIS MOD06CD.
C
C !Input Parameters:
C      MODIS radiance data at 0.66, 0.86, 0.94, 1.24 and 1.38 micron in
C      reflectance unit, and water vapor transmittance table at 0.86,
C      0.94, and 1.38 micron.
C
C !Output Parameters:
C      Cirrus reflectance.
C
C !Revision History:
C $Log: Process_mod06cd_V2.f,v $
C  Modified by D. A. Chu 6/28/99
C  To separate land and ocean pixels using land/sea flag from cloud mask
C  and to handle the uneven length of boxes
C
C !TEAM-UNIQUE HEADER:
C   This software was developed by the MODIS Science Team
C   for the National Aeronautics and Space Administration,
C   Goddard Space Flight Center.
C
C !REFERENCES and CREDITS:
C      gao@rsd.nrl.navy.mil
C
C !DESIGN NOTES:
C      Data and variables are passed with common statements
C      contained in "COMMONS_INC".
C
C  Functions and Subroutines:
C     PR_CIRRUS_REFLECTANCES_MODE_2
C      GET_SLOPE
C         SORT2_MOD06CD
C         HUNT_MOD06CD
C         LOCATE
C         GET_SLOPE_OCEAN or GET_SLOPE_LAND
C            SORT2_MOD06CD
C         GET_TRANSM_0p94
C            SORT2_MOD06CD
C         GET_TRANSM_1p38
C            LOCATE_MOD06CD
C         CLEAR_OR_CIRRUS
C      GET_CIRRUS_REFL
C     END
C
C !END
C-----------------------------------------------------------------------

      IMPLICIT NONE
      SAVE

      INCLUDE 'COMMONS_cirrus.inc'

      INTEGER Scans_Per_Granule,No_Frames(Max_Scan_Per_Granule)
      INTEGER Min_No_Frames,Re_X_Frames,Re_Y_Lines,Index_Ocean
      INTEGER NJ,MJ,K,IK,I,N_Box_X,N_Box_Y,Jb,J_st,J_ed,Ib,I_st,I_ed
      INTEGER J
      REAL Percent_Ocean,Count_Ocean

      Min_No_Frames=Max_Frames_Per_Line

      DO I=1,Scans_Per_Granule

        IF(No_Frames(I).LE.Min_No_Frames) Min_No_Frames=No_Frames(I)

      ENDDO

      N_Box_X = Min_No_Frames/N_Sub_X
      N_Box_Y = Scans_Per_Granule*10/N_Sub_Y

      Re_X_Frames = Min_No_Frames - N_Box_X * N_Sub_X
      Re_Y_Lines = Scans_Per_Granule*10 - N_Box_Y * N_Sub_Y

      IF(Re_X_Frames.GT.0) N_Box_X=N_Box_X+1
      IF(Re_Y_Lines.GT.0) N_Box_Y=N_Box_Y+1

C
C Derive cirrus path radiance for given number of boxes
C
C Land_Sea_Flag to separate land and ocean and check the validity of
C the measurements (> 0. check). K should be saved for each box.
C
       DO IK=1,N

         DATA_FLOAT_066(IK) = 1000.0
         DATA_FLOAT_086(IK) = 1000.0
         DATA_FLOAT_094(IK) = 1000.0
         DATA_FLOAT_124(IK) = 1000.0
         DATA_FLOAT_138(IK) = 1000.0
         DATA_FLOAT_164(IK) = 1000.0
         DATA_FLOAT_213(IK) = 1000.0

       ENDDO

       DO Jb = 1, N_Box_Y

         J_st = (Jb - 1) * N_Sub_Y + 1
         J_ed =  Jb * N_Sub_Y
C         IF(Jb.EQ.N_Box_Y) J_ed = Jb * N_Sub_Y + Re_Y_Lines
         IF(Re_Y_Lines.GT.0.AND.Jb.EQ.N_Box_Y) J_ed = J_st-1+Re_Y_Lines

         DO Ib = 1, N_Box_X

           I_st = (Ib - 1) * N_Sub_X + 1
           I_ed =  Ib * N_Sub_X
C           IF(Ib.EQ.N_Box_X) I_ed = Ib * N_Sub_X + Re_X_Frames
           IF(Re_X_Frames.GT.0.AND.Ib.EQ.N_Box_X) I_ed = I_st-1+Re_X_Frames

           K = 0
           Count_Ocean=0.0
           Index_Ocean=0

           DO J = J_st, J_ed
             DO I = I_st, I_ed

               K = K + 1

                 IF(I_LAND_WATER_Img(I,J).EQ.0) THEN
                   Count_Ocean=Count_Ocean+1
                 ENDIF
                 DATA_FLOAT_066(K) = REFL_0P67_Img(I,J)
                 DATA_FLOAT_086(K) = REFL_0P86_Img(I,J)
                 DATA_FLOAT_094(K) = REFL_0P94_Img(I,J)
                 DATA_FLOAT_124(K) = REFL_1P24_Img(I,J)
                 DATA_FLOAT_138(K) = REFL_1P38_Img(I,J)
                 DATA_FLOAT_1380(K) = REFL_1P38_Img(I,J)
                 DATA_FLOAT_1381(K) = REFL_1P38_Img(I,J)
                 DATA_FLOAT_164(K) = REFL_1P64_Img(I,J)
                 DATA_FLOAT_213(K) = REFL_2P13_Img(I,J)

             END DO
           END DO

             Percent_Ocean=(Count_Ocean/FLOAT(K))*100.
C
C Assume index_ocean = 1 (ocean only)
C
             IF(Percent_Ocean.EQ.100.) Index_Ocean=1
C             IF(Percent_Ocean.GE.0.) Index_Ocean=1
C              Index_Ocean=1

             CALL GET_SLOPE(Index_Ocean,K)
C             CALL GET_CIRRUS_REFL(K)

C
C Modified by D. A. Chu
C For three different wavelengths over ocean
C

             IF(INDEX_OCEAN.EQ.1) MJ=3
             IF(INDEX_OCEAN.EQ.0) MJ=1
             DO NJ=1,MJ
               INDEX_OCELAND(Ib,Jb,NJ)=Index_Ocean
               SLOPE_IMG(Ib,Jb,NJ)=SLOPE(1,NJ)
             ENDDO

          END DO
       END DO

      RETURN
      END


C***********************************************************************
       SUBROUTINE GET_SLOPE(Index_Ocean,K)
C-----------------------------------------------------------------------
C!F77
C
C!DESCRIPTION: This subroutine is to derive a transmittance coefficient
C            due to water vapor at 1.38 micron.
C
C !Input Parameters:
C      MODIS radiance data at 0.66, 0.86, 0.94, 1.24, and 1.38 micron in
C      reflectance unit, water vapor transmittance table at 0.86, 0.94a
C      and 1.38 micron.
C
C !Output Parameters:
C      Transmittance coefficient due to water vapor at 1.38 micron.
C
C !Revision History:
C $Log: Process_mod06cd_V2.f,v $
c 10/20/1999 fhliang
c fixed prolog.
c
c!Team-Unique Header:
c
C !References and Credits:
C       gao@rsd.nrl.navy.mil
C
C  DESIGN NOTES: data and variables are passed with common statements
C              contained in 'COMMONS_cirrus.inc'.
C
C !Functions and Subroutines:
C            SORT2_MOD06CD
C            HUNT_MOD06CD
C            LOCATE
C            GET_SLOPE_OCEAN or GET_SLOPE_LAND
C              SORT2_MOD06CD
C            GET_TRANSM_0p94
C              SORT2_MOD06CD
C            GET_TRANSM_1p38
C              LOCATE
C            CLEAR_OR_CIRRUS
C
C !END
C-----------------------------------------------------------------------

      IMPLICIT NONE
      SAVE

      INCLUDE 'COMMONS_cirrus.inc'

      INTEGER Index_Ocean,I,J,K,MJ,N_INTERVALS,N_SUM_MAX
      INTEGER N_STEP_START,N_STEP_END
      REAL DELT_N_STEP,STEP

C- Step 1: Sort the image array DATA_FLOAT_138, then get DATA_FLOAT_138
C        in increasing order and re-arranged array DATA_FLOAT_124 or
C        DATA_FLOAT_066.
C
C      INDEX_OCEAN = 1
C
C--     INDEX_OCEAN = 1 for ocean and =0 for land

      IF (INDEX_OCEAN .EQ. 1) THEN
C         CALL SORT2 (N, DATA_FLOAT_138, DATA_FLOAT_124)
         CALL SORT2_MOD06CD (N, K, DATA_FLOAT_138, DATA_FLOAT_124)
         CALL SORT2_MOD06CD (N, K, DATA_FLOAT_1380, DATA_FLOAT_164)
         CALL SORT2_MOD06CD (N, K, DATA_FLOAT_1381, DATA_FLOAT_213)
      END IF

      IF (INDEX_OCEAN .EQ. 0) THEN
C         CALL SORT2 (N, DATA_FLOAT_138, DATA_FLOAT_066)
         CALL SORT2_MOD06CD (N, K, DATA_FLOAT_138, DATA_FLOAT_066)
      END IF

C-   Determine lower (REFL_TH_MIN_138) and upper (REFL_TH_MAX_138)
C      threshholds for image array DATA_FLOAT_138. The initial value
C      of REFL_TH_MIN_138 is given in file 'COMMONS_cirrus.inc'. In this
C      program, it was assumed that the input apparent reflectances were
C      multiplied by 1000.

C      REFL_TH_MAX_138 = DATA_FLOAT_138(IFIX(0.998*float(N)))
      REFL_TH_MAX_138 = DATA_FLOAT_138(IFIX(0.998*float(K)))
       IF(REFL_TH_MAX_138 .GT. 1000.) REFL_TH_MAX_138 = 1000.


C-  Determine the indices N_MIN_138 & N_MAX_138 for REFL_TH_MIN_138
C      and REFL_TH_MAX_138, respectively, in array DATA_FLOAT_138

C      CALL LOCATE(DATA_FLOAT_138, N, REFL_TH_MIN_138, N_MIN_138)
C      CALL LOCATE(DATA_FLOAT_138, N, REFL_TH_MAX_138, N_MAX_138)
      CALL LOCATE_MOD06CD(DATA_FLOAT_138, N, K, REFL_TH_MIN_138, N_MIN_138)
      CALL LOCATE_MOD06CD(DATA_FLOAT_138, N, K, REFL_TH_MAX_138, N_MAX_138)

C- Step 2: Calculate histogram of 1.38-micron image array DATA_FLOAT_138

      STEP = (REFL_TH_MAX_138 - REFL_TH_MIN_138) / FLOAT(N_STEP)

       DO I = 1, N_STEP + 1
          X_138(I) = REFL_TH_MIN_138 + float(I - 1) * STEP
      END DO

      DO I = 1, N_STEP
          N_SUM (I) = 0
      END DO

      J = 0
      DO I = N_MIN_138, N_MAX_138
          CALL HUNT_MOD06CD(X_138, N_STEP + 1, DATA_FLOAT_138(I), J)
          N_SUM (J)  = N_SUM (J) + 1
      END DO

C- Step 3: Find the peak of the histogram and related parameters

       I_HIST_PEAK = 0
      N_SUM_MAX   = 0

      DO I = 1, N_STEP
          IF (N_SUM(I) .GT. N_SUM_MAX) THEN
            N_SUM_MAX   = N_SUM(I)
            I_HIST_PEAK = I
          END IF
      END DO

C  In order to prevent an unlikely problem to occur ...
       IF(I_HIST_PEAK.GT.N_STEP - 10) I_HIST_PEAK = N_STEP - 10

C---       N_INTERVAL_TOT = N_STEP - I_HIST_PEAK
C---       DELT_N_STEP    = FLOAT( N_INTERVAL_TOT ) / 10.

C---       X_INTERVAL_TOT = X_138(N_STEP + 1) - X_138(I_HIST_PEAK)
C---       N_SLOPE = IFIX(X_INTERVAL_TOT / R_138_INTERVAL_SLOPE)

C---       IF (N_SLOPE .LT. 1) N_SLOPE = 1
C---       IF (N_SLOPE .GT. 4) N_SLOPE = 4

       N_INTERVALS = N_STEP - I_HIST_PEAK
       DELT_N_STEP = FLOAT( N_INTERVALS ) / 3.

       N_SLOPE = 1

C- Step 4: Get slopes

C
C Modified by D. A. Chu for deriving 3 slopes over ocean surface
C at 1.38, 1.64 and 2.1 micron (1 slope over land surface)
C
      IF(INDEX_OCEAN.EQ.1) MJ=3
      IF(INDEX_OCEAN.EQ.0) MJ=1

      IF(INDEX_OCEAN.EQ.1) THEN
        DO J=1,MJ
          DO I = 1, N_SLOPE + 1
C            IF(J.EQ.1) THEN
              N_STEP_START = IFIX(FLOAT(I)  * DELT_N_STEP) +I_HIST_PEAK +1
              N_STEP_END   = IFIX(FLOAT(I+1)* DELT_N_STEP) +I_HIST_PEAK
C            ENDIF

C---          N_STEP_START = IFIX( ((I - 1) * 8. / FLOAT(N_SLOPE) + 1. )
C---         $                      * DELT_N_STEP ) + I_HIST_PEAK
C---          N_STEP_END   = IFIX( ((I - 1) * 8. / FLOAT(N_SLOPE) + 2. )
C---         $                      * DELT_N_STEP ) + I_HIST_PEAK

            CALL GET_SLOPE_OCEAN(K,N_STEP_START,N_STEP_END,X1(I,J),Y1(I),J)
          END DO
        END DO
      ENDIF

      IF (INDEX_OCEAN.EQ.0) THEN
        DO J=1,MJ
          DO I = 1, N_SLOPE + 1
            N_STEP_START = IFIX(FLOAT(I)  * DELT_N_STEP) +I_HIST_PEAK +1
            N_STEP_END   = IFIX(FLOAT(I+1)* DELT_N_STEP) +I_HIST_PEAK

C---          N_STEP_START = IFIX( ((I - 1) * 8. / FLOAT(N_SLOPE) + 1. )
C---         $                      * DELT_N_STEP ) + I_HIST_PEAK
C---          N_STEP_END   = IFIX( ((I - 1) * 8. / FLOAT(N_SLOPE) + 2. )
C---         $                      * DELT_N_STEP ) + I_HIST_PEAK

            CALL GET_SLOPE_LAND(K,N_STEP_START,N_STEP_END,X1(I,J),Y1(I))
          END DO
        END DO
      ENDIF


      DO J=1,MJ
        DO I = 1, N_SLOPE
          IF (ABS(X1(I,J)-X1(I+1,J)).LE.ZEPS)THEN
            SLOPE(I,J)=0.0
          ELSE
            SLOPE(I,J)=(Y1(I+1)-Y1(I))/(X1(I+1,J)-X1(I,J))
          END IF

        END DO
      END DO


      CALL GET_TRANSM_0P94(K)
      CALL GET_TRANSM_1P38(K)
      CALL CLEAR_OR_CIRRUS(MJ)


      RETURN
      END


C***********************************************************************
       SUBROUTINE GET_SLOPE_LAND(N1,N_STEP_START, N_STEP_END, X, Y)

C-----------------------------------------------------------------------
C
C!F77
C
C!DESCRIPTION: This subroutine is to derive the water vapor transmittance
C            coefficient over ocean surface which is the slope of the
C            scatter plot of 1.38-micron reflectanc against 1.24-micron
C            reflectance.
C
C !Input Parameters:
C      MODIS radiance data in reflectance unit at 1.24 and 1.38 micron.
C
C !Output Parameters:
C      Slope of the scatter plot of 1.38-micron reflectanc against
C      0.66 micron reflectance.
C
C !Revision History:
C $Log: Process_mod06cd_V2.f,v $
c 10/20/1999 fhliang
c fixed prolog.
c
c!Team-Unique Header:
c
C !References and Credits:
C       gao@rsd.nrl.navy.mil
C
C  DESIGN NOTES: data and variables are passed with common statements
C              contained in 'COMMONS_cirrus.inc'.
C
C !Functions and Subroutines:
C              SORT2_MOD06CD
C
C !END
C-----------------------------------------------------------------------

       IMPLICIT NONE
       SAVE

       INCLUDE 'COMMONS_cirrus.inc'

       INTEGER I,K,N1,N_STEP_START,N_STEP_END,N_LOWER,N_UPPER
       INTEGER NSTART_L_138,NEND_L_138,NSTART_L_066,NEND_L_066
       INTEGER N_AVG_X_Y,N_SORT_066_L
       REAL X,Y,PER_LOWER1_138,PER_LOWER2_138

       N_LOWER = N_MIN_138
       DO K = 1, N_STEP_START
          N_LOWER = N_LOWER + N_SUM(K)
       END DO

C       PER_LOWER1_138 = FLOAT(N_LOWER)/FLOAT(N)
       PER_LOWER1_138 = FLOAT(N_LOWER)/FLOAT(N1)

       N_UPPER = N_MIN_138
       DO K = 1, N_STEP_END
          N_UPPER = N_UPPER + N_SUM(K)
       END DO

C       PER_LOWER2_138 = FLOAT(N_UPPER)/FLOAT(N)
       PER_LOWER2_138 = FLOAT(N_UPPER)/FLOAT(N1)

       NSTART_L_138 = N_LOWER
       NEND_L_138   = N_UPPER

C--    Creat the sub-image arrays BUF_L_066 & BUF_L_138 based on the
C--    indices NSTART_L_138 & NEND_L_138

       DO I = NSTART_L_138, NEND_L_138
          K = I - NSTART_L_138 + 1
          BUF_L_066(K) = DATA_FLOAT_066(I)
          BUF_L_138(K) = DATA_FLOAT_138(I)
       END DO

C--   Sort BUF_L_066 and get it in increasing order, and BUF_L_138
C--   rearranged in the lower part of the scatter plot

       N_SORT_066_L = NEND_L_138 - NSTART_L_138 + 1

              CALL SORT2_MOD06CD (N_SORT_066_L, N_SORT_066_L,
     *                    BUF_L_066, BUF_L_138)

C--    Creat the sub-sub-image arrays in 1.38 & 0.66 micron and calculate
C--    averaged X2,Y2 coordinates. The values of PER_LOWER1_066 and
C--    PER_LOWER2_066 are given in file 'COMMONS_cirrus.inc'.

       NSTART_L_066 = IFIX (FLOAT(N_SORT_066_L) * PER_LOWER_066)
       NEND_L_066   = IFIX (FLOAT(N_SORT_066_L) * PER_UPPER_066)

           N_AVG_X_Y = NEND_L_066 - NSTART_L_066 + 1

       IF (N_AVG_X_Y .GT. N_AVG_MAX) THEN
           NEND_L_066 = NSTART_L_066 + N_AVG_MAX - 1
           N_AVG_X_Y  = NEND_L_066 - NSTART_L_066 + 1
       END IF

       X = 0.0
       Y = 0.0

       DO I = NSTART_L_066, NEND_L_066
          X = BUF_L_066(I) + X
          Y = BUF_L_138(I) + Y
       END DO

       X = X / FLOAT(NEND_L_066 - NSTART_L_066 + 1)
       Y = Y / FLOAT(NEND_L_066 - NSTART_L_066 + 1)


       RETURN
       END


C***********************************************************************
       SUBROUTINE GET_SLOPE_OCEAN(N1,N_STEP_START, N_STEP_END, X, Y, J)

C-----------------------------------------------------------------------
C
C!F77
C
C!DESCRIPTION: This subroutine is to derive the water vapor transmittance
C            coefficient over ocean surface which is the slope of the
C            scatter plot of 1.38-micron reflectanc against 1.24-micron
C            reflectance.
C
C !Input Parameters:
C      MODIS radiance data in reflectance unit at 1.24 and 1.38 micron.
C
C !Output Parameters:
C      Slope of the scatter plot of 1.38-micron reflectance against
C      1.24-micron reflectance.
C
C !Revision History:
C $Log: Process_mod06cd_V2.f,v $
c 10/20/1999 fhliang
c fixed prolog.
c
c!Team-Unique Header:
c
C !References and Credits:
C       gao@rsd.nrl.navy.mil
C
C  DESIGN NOTES: data and variables are passed with common statements
C              contained in 'COMMONS_cirrus.inc'.
C
C !Functions and Subroutines:
C              SORT2_MOD06CD
C
C !END
C-----------------------------------------------------------------------

      IMPLICIT NONE
      SAVE

      INCLUDE 'COMMONS_cirrus.inc'

      INTEGER N_AVG_X_Y,N_STEP_START,N_STEP_END
      INTEGER N1,I,J,K,N_LOWER,N_UPPER,NSTART_L_138,NEND_L_138
      INTEGER N_SORT_124_L,NSTART_L_124,NEND_L_124
      REAL X,Y,PER_LOWER1_138,PER_LOWER2_138

      N_LOWER = N_MIN_138
      DO K = 1, N_STEP_START
          N_LOWER = N_LOWER + N_SUM(K)
      END DO

C      PER_LOWER1_138 = FLOAT(N_LOWER)/FLOAT(N)
      PER_LOWER1_138 = FLOAT(N_LOWER)/FLOAT(N1)

      N_UPPER = N_MIN_138
      DO K = 1, N_STEP_END
          N_UPPER = N_UPPER + N_SUM(K)
      END DO

C      PER_LOWER2_138 = FLOAT(N_UPPER)/FLOAT(N)
      PER_LOWER2_138 = FLOAT(N_UPPER)/FLOAT(N1)

       NSTART_L_138 = N_LOWER
       NEND_L_138   = N_UPPER

C--    Creat the sub-image arrays BUF_L_124 & BUF_L_138 based on the
C--    indices NSTART_L_138 & NEND_L_138

       DO I = NSTART_L_138, NEND_L_138
               K = I - NSTART_L_138 + 1
         BUF_L_124(K) = DATA_FLOAT_124(I)
         BUF_L_138(K) = DATA_FLOAT_138(I)
         IF(J.EQ.2) BUF_L_124(K) = DATA_FLOAT_164(I)
         IF(J.EQ.3) BUF_L_124(K) = DATA_FLOAT_213(I)
      END DO

C--   Sort BUF_L_124 and get it in increasing order, and BUF_L_138
C--   rearranged in the lower part of the scatter plot

      N_SORT_124_L = NEND_L_138 - NSTART_L_138 + 1

            CALL SORT2_MOD06CD (N_SORT_124_L, N_SORT_124_L,
     *                  BUF_L_124, BUF_L_138)

C--    Creat the sub-sub-image arrays in 1.38 & 1.24 micron and calculate
C--    averaged X2,Y2 coordinates. The values of PER_LOWER1_124 and
C--    PER_LOWER2_124 are given in file 'COMMONS_cirrus.inc'.

       NSTART_L_124 = IFIX (FLOAT(N_SORT_124_L) * PER_LOWER_124)
       NEND_L_124   = IFIX (FLOAT(N_SORT_124_L) * PER_UPPER_124)


           N_AVG_X_Y = NEND_L_124 - NSTART_L_124 + 1

       IF (N_AVG_X_Y .GT. N_AVG_MAX) THEN
           NEND_L_124 = NSTART_L_124 + N_AVG_MAX - 1
           N_AVG_X_Y  = NEND_L_124 - NSTART_L_124 + 1
       END IF

       X = 0.0
       Y = 0.0

       DO I = NSTART_L_124, NEND_L_124

C         IF(J.EQ.1) THEN
           X = BUF_L_124(I) + X
C         ELSE IF(J.EQ.2) THEN
C          X = BUF_L_164(I) + X
C         ELSE IF(J.EQ.3) THEN
C          X = BUF_L_213(I) + X
C         ENDIF

          Y = BUF_L_138(I) + Y
       END DO

       X = X / FLOAT(NEND_L_124 - NSTART_L_124 + 1)
       Y = Y / FLOAT(NEND_L_124 - NSTART_L_124 + 1)

      RETURN
      END


C***********************************************************************
      SUBROUTINE GET_TRANSM_0P94(N1)

C-----------------------------------------------------------------------
C!F77
C
C!DESCRIPTION: This subroutine derives water vapor transmittance at 0.94
C            micron for cloudy pixels, based on the fact that water vapor
C            transmittance for cloudy piexls is generally larger than
C            that for clear piexls.
C
C !Input Parameters:
C      MODIS radiance data in reflectance unit at 0.66, 0.86, 0.94 and
C      1.24 micron.
C
C !Output Parameters:
C      Water vapor transmittance at 0.94 micron for cloudy piexls.
C
C !Revision History:
C $Log: Process_mod06cd_V2.f,v $
c 10/20/1999 fhliang
c fixed prolog.
c
c!Team-Unique Header:
c
C !References and Credits:
C       gao@rsd.nrl.navy.mil
C
C DESIGN NOTES: data and variables are passed with common statements
C             contained in 'COMMONS_cirrus.inc'.
C
C !Functions and Subroutines:
C            SORT2_MOD06CD
C
C !END
C-----------------------------------------------------------------------

      IMPLICIT NONE
      SAVE

      INCLUDE 'COMMONS_cirrus.inc'

      INTEGER N1,I,NSTART,NEND

C      DO I = 1, N
      DO I = 1, N1
          TRANSM_R_0p94_MES(I) = DATA_FLOAT_094(I) /DATA_FLOAT_086(I)
      END DO

C       CALL SORT2 (N,N,TRANSM_R_0p94_MES,DATA_FLOAT_124)
       CALL SORT2_MOD06CD(N1,N1,TRANSM_R_0p94_MES,DATA_FLOAT_124)

C       NSTART = IFIX (FLOAT(N) * 0.99)
C       NEND   = IFIX (FLOAT(N) * 0.995)
       NSTART = IFIX (FLOAT(N1) * 0.99)
       NEND   = IFIX (FLOAT(N1) * 0.995)
      DO I = NSTART, NEND
          TRANSM_R_0p94_CLOUDY= TRANSM_R_0p94_CLOUDY +
     $                        TRANSM_R_0p94_MES(I)
      END DO

      TRANSM_R_0p94_CLOUDY = TRANSM_R_0p94_CLOUDY /
     $                      FLOAT(NEND - NSTART +1)

      RETURN
      END


C****************************************************************************
      SUBROUTINE GET_TRANSM_1P38(N1)

C-----------------------------------------------------------------------
C!F77
C
C!DESCRIPTION: This subroutine derives a threshold water vapor transmittance
C            coefficient at 1.38 micron based on the water vapor
C            transmittance at 0.94 micron for cloudy piexls derived in
C            GET_TRANSM_0p94, by a look-up table procedure.
C
C !Input Parameters:
C      MODIS radiance data in reflectance unit at 0.86, 0.94 and 1.38 micron.
C
C !Output Parameters:
C      Threshold water vapor transmittance at 1.38 micron.
C
C !Revision History:
C $Log: Process_mod06cd_V2.f,v $
c 10/20/1999 fhliang
c fixed prolog.
c
c!Team-Unique Header:
c
C !References and Credits:
C       gao@rsd.nrl.navy.mil
C
C DESIGN NOTES: data and variables are passed with common statements
C             contained in 'COMMONS_cirrus.inc'.
C
C !Functions and Subroutines:
C            LOCATE
C
C !END
C-----------------------------------------------------------------------

      IMPLICIT NONE
      SAVE

      INCLUDE 'COMMONS_cirrus.inc'

      INTEGER N1,I,N_TABLE_RESULT

      DO I = 1, N_TRANSM_TABLE
          TRANSM_R_0p94(I) = TRANSM_0p94(I) / TRANSM_0p86(I)
      END DO

      CALL LOCATE_MOD06CD(TRANSM_R_0p94,N_TRANSM_TABLE,N_TRANSM_TABLE,
     *     TRANSM_R_0p94_CLOUDY,N_TABLE_RESULT)

      TRANSM_1p38_THRESHOLD = TRANSM_1p38 (N_TABLE_RESULT)

      RETURN
      END


C***********************************************************************
      SUBROUTINE CLEAR_OR_CIRRUS(MJ)

C-----------------------------------------------------------------------
C!F77
C
C!DESCRIPTION: This subroutine distinguish cloudy or clear scene by
C            comparing the values of the slope derived in GET_SLOPE
C            and the threshold water vapor transmittance in 1.38 micron
C            derived in GET_TRANSM_1p38.
C
C !Input Parameters:
C             Threshold at 1.38 micron and the slope.
C
C !Output Parameters:
C             Water vapor transmittance coefficient in 1.38 miicron.
C
C !Revision History:
C $Log: Process_mod06cd_V2.f,v $
c 10/20/1999 fhliang
c fixed prolog.
c
c!Team-Unique Header:
c
C !References and Credits:
C       gao@rsd.nrl.navy.mil
C
C  DESIGN NOTES: data and variables are passed with common statements
C              contained in 'COMMONS_cirrus.inc'.
C
C !Functions and Subroutines: None.
C
C !END
C-----------------------------------------------------------------------

      IMPLICIT NONE
      SAVE

      INCLUDE 'COMMONS_cirrus.inc'

      INTEGER I,J,MJ

C--     In order to get rid of clear scenes, the value of SLOPE
C--     is compared with the value of TRANSM_1p38_THRESHOLD, based
C--     on the fact that the Slope is generally larger for cloudy
C--     scene than that of clear scene. The transmittance coefficient
C--     is set to be zero for clear scenes, which lead to zero of Cirrus
C--     reflectance.

      DO J=1,MJ
        DO I=1,N_SLOPE
               IF (SLOPE(I,J).GT.TRANSM_1p38_THRESHOLD) COEFF_FINAL(I,J)
     $                                      = 1./ SLOPE(I,J)
         IF (SLOPE(I,J) .LE. TRANSM_1p38_THRESHOLD) COEFF_FINAL(I,J) = 0.

        END DO
      END DO

      RETURN
      END


C***********************************************************************
      SUBROUTINE SORT2_MOD06CD(n,n1,arr,brr)

C-----------------------------------------------------------------------
C !F77
C
C!DESCRIPTION: This subroutine is to sort array of arr in an ascending
C              order and array brr will automatically reordered based
C              upon array arr.
C
C !Input Parameters:
C
C        N           Number of elements in array arr
C        N1          Number of elements in array arr
C        ARR         Input array ARR
C        BRR         Input array BRR
C
C !Output Parameters:
C
C       ARR          Array ARR in ascending order
C       BRR          Array BRR in ascending order
C
C !Revision History:
C $Log: Process_mod06cd_V2.f,v $
c 10/20/1999 fhliang
c fixed prolog.
c
c!Team-Unique Header:
c
C !References and Credits:
C
C  Numerical receipes
C
C !Functions and Subroutines: None.
C
C !END
C-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER n, M, NSTACK, n1, ir, i, j, k, l, jstack, im
      REAL arr(n), brr(n), a, b, temp
      PARAMETER (M = 7, NSTACK = 50)

c      Sorts an array arr(1:n) into ascending order using Quicksort, while
c      making the corresponding rearrangement of the array brr(1:n).

      INTEGER istack(NSTACK)
      jstack = 0
      l = 1
C      ir = n
      ir = n1
1      if (ir-l .lt. M) then
            do  j = l+1,ir
                  a = arr(j)
                  b = brr(j)
                  do  i = j-1, 1, -1
                        if (arr(i) .le. a) goto 2
                        arr(i+1) = arr(i)
                        brr(i+1) = brr(i)
                  end do
                  i       = 0
2                  arr(i+1) = a
                  brr(i+1) = b
            end do
            if (jstack .eq. 0) return
            ir     = istack(jstack)
            l      = istack(jstack-1)
            jstack = jstack-2
      else
            k = (l+ir)/2
            temp     = arr(k)
            arr(k)   = arr(l+1)
            arr(l+1) = temp
            temp     = brr(k)
            brr(k)   = brr(l+1)
            brr(l+1) = temp
            if (arr(l+1) .gt. arr(ir)) then
                  temp     = arr(l+1)
                  arr(l+1) = arr(ir)
                  arr(ir)  = temp
                  temp     = brr(l+1)
                  brr(l+1) = brr(ir)
                  brr(ir)  = temp
            endif
            if(arr(l) .gt. arr(ir)) then
                  temp    = arr(l)
                  arr(l)  = arr(ir)
                  arr(ir) = temp
                  temp    = brr(l)
                  brr(l)  = brr(ir)
                  brr(ir) = temp
            endif
            if (arr(l+1) .gt. arr(l)) then
                  temp     = arr(l+1)
                  arr(l+1) = arr(l)
                  arr(l)   = temp
                  temp     = brr(l+1)
                  brr(l+1) = brr(l)
                  brr(l)   = temp
            endif
            i = l+1
            j = ir
            a = arr(l)
            b = brr(l)
3            continue
                  i = i+1
            if (arr(i) .lt. a) goto 3
4            continue
                  j = j-1
            if (arr(j) .gt. a) goto 4
            if (j .lt. i) goto 5
            temp   = arr(i)
            arr(i) = arr(j)
            arr(j) = temp
            temp   = brr(i)
            brr(i) = brr(j)
            brr(j) = temp
            goto 3
5            arr(l) = arr(j)
            arr(j) = a
            brr(l) = brr(j)
            brr(j) = b
            jstack = jstack+2
C---            if(jstack .gt. NSTACK) pause 'NSTACK too small in sort2'
            if(jstack .gt. NSTACK) then
               stop
            end if
C---
            if(ir-i+1 .ge. j-1) then
                  istack(jstack)   = ir
                  istack(jstack-1) = i
                  ir             = j-1
            else
                  istack(jstack)   = j-1
                  istack(jstack-1) = l
                  l              = i
            endif
      endif
      goto 1

      end


C***********************************************************************
      SUBROUTINE LOCATE_MOD06CD(xx,n,n1,x,j)

C-----------------------------------------------------------------------
C !F77
C
C!DESCRIPTION: This subroutine is to locate value of x in an array xx
C
C !Input Parameters:
C
C        XX          Given an array
C        N           Number of elements in array xx
C        N1          Number of elements in array xx
C        X           Input value
C
C !Output Parameters:
C
C       J            Location of X in array XX
C
C !Revision History:
C $Log: Process_mod06cd_V2.f,v $
c 10/20/1999 fhliang
c fixed prolog.
c
c!Team-Unique Header:
c
C !References and Credits:
C
C  Numerical receipes
C
C !Functions and Subroutines: None.
C
C !END
C-----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER j,n,n1,jl,jm,ju
      REAL x,xx(n)

      jl=0
C      ju=n+1
      ju=n1+1
10    if(ju-jl.gt.1)then
       jm=(ju+jl)/2
       if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then
         jl=jm
       else
         ju=jm
       endif
      goto 10
      endif
      if(x.eq.xx(1))then
       j=1
      else if(x.eq.xx(n))then
       j=n-1
      else
       j=jl
      endif

      return
      END
