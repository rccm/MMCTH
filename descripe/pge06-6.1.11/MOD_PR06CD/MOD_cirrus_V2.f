      SUBROUTINE Process_Mod06CD(Error_Flag)

      IMPLICIT NONE
      INCLUDE 'COMMONS_cirrus.inc'
      INCLUDE 'mapi.inc'
      INCLUDE 'PGS_MODIS_39500.f'
      INCLUDE 'PGS_PC.f'
      INCLUDE 'PGS_IO.f'
      INCLUDE 'PGS_SMF.f'
      INCLUDE 'cirrus.inc'

C-----------------------------------------------------------------------
C !F77
C
C !Description: This driver program derives cirrus reflectance in the 0.4 - 1.0
C               micron region and contral index.
C
C !Input Parameters:
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
C !Output Parameters:
C   INTEGER     Cirrus_Refl    Two dimensional (2-D) array for passing
C                              retrieved cirrus reflectance data
C   INTEGER     I_QA_cirrus    Two dimensional (2-D) array for passing
C                              quality assurance data for cirrus reflectance
C                              and contrail index. Index 0 refers to bad data,
C                              index 1 refers to non-cirrus pixel, index 2 refers
C                              to cirrus pixel, and index 3 refers to contrail pixel.
C  Logical      Error_Flag     Subroutine return status:  Error_Flag=.true. or
C                              .false. Currently assigned value .false.
C
C
C
C
C
C
C
C !Revision History:
C
C  Modified by Ping Yang 11/2/2001
C
C
C $Log: MOD_PR05.f,v $
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
C   phone: 301-286-4080
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
C          FILEOC
C          GetModisDat_Cirrus
C          MODIS_SMF_SETDYNAMICMSG
C          MOD_cirrus_PutData
C          READ_TABLE
C          PROCESS_CIRRUS_REFLECTANCES
C          PROCESS_CONTRAIL_INDEX
C          OA_CIRRUS
C
C !END
C-----------------------------------------------------------------------

C Liqun 02/18/98:  Added parameters FUNCNAME
C rhucek 01/04/99: Made variable PROCESS_FLAG a PARAMETER

      CHARACTER*(*) FUNCNAME
      PARAMETER    (FUNCNAME = 'MOD_cirrus_V2.f')

      CHARACTER*30  PROCESS_FLAG
      PARAMETER    (PROCESS_FLAG = 'MODIS')


      CHARACTER*5    ioflag
      CHARACTER*40   att_N,dtype
      CHARACTER*1024 msgbuf
      CHARACTER*(PGSd_PC_FILE_PATH_MAX) FN_L1B
C
C ............ added by Ping Yang (11/2/2001)......................
C
      INTEGER N_PARTITION_X, N_PARTITION_Y
      PARAMETER(N_PARTITION_X=4, N_PARTITION_Y=4)
      INTEGER IPY_A, IPY_B, IPY_CENT, JPY_A, JPY_B, JPY_CENT

      REAL PY_SOLAR_ANG(N_PARTITION_X, N_PARTITION_Y), 
     &  PY_VIEW_ANG(N_PARTITION_X,N_PARTITION_Y)
 
C ............ end PY's addition  ......................

      INTEGER Max_EV_Frames

      CHARACTER*4  Scan_Type(Max_Scan_Per_Granule),ScanType

      INTEGER i, idebug, j, K, RTN, nms, Scans_Per_Granule,
     &        Total_No_Lines, BeginScan_No, EndScan_No, Scan_No,
     &        handle_QC,
     &        Buf_Size1, Buf_Size2,
     &        Modfil_Geo(MODFILLEN), Modfil_CldMsk(MODFILLEN),
     &        Modfil_MOD06cd(MODFILLEN),
     &        Modfil_L1B(MODFILLEN),
     &        handle_trans,Data_Size(2),
     &        Data_Size_For_Granule(2,Max_Scan_Per_Granule),
     &        No_Frames(Max_Scan_Per_Granule),
     &        bsize
      INTEGER GetModisDat_MOD06CD

C  Liqun Ma added 1 declaration line
      INTEGER handle_test,pgs_io_gen_openf

      INTEGER CldMsk      (Max_Frames_Per_Line, No_Lines_Per_Scan),
     2        LandSea_Flag(Max_Frames_Per_Line, No_Lines_Per_Scan),
     3        UFQ_Flag(Max_Frames_Per_Line, No_Lines_Per_Scan)

c Liqun Ma added logical data type
      LOGICAL Error_Flag

c-----Dr. Wei Han, han@han.nrl.navy.mil--------------------
      CHARACTER*1 work_buf
      INTEGER  Set_CoreMetadata_QC
      INTEGER  PGS_PC_GetConfigData
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
C 1-km resolution data array
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

      BYTE Buf_Un      (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Buf_Sa      (Max_Frames_Per_Line, No_Lines_Per_Scan)

      SAVE




C Initialize variables
      idebug       =  1
      Error_Flag   = .false.

c rhucek 01/04/99: made PROCESS_FLAG a PARAMETER (initialized in declarations)
c     PROCESS_FLAG = 'MODIS'


c-----Dr. Wei Han, han@han.nrl.navy.mil--------------------
c Check diagnostic parameter for debugging
c      (1:debugging: 0: no debugging)
       rtn=PGS_PC_GetConfigData(LRN_DIAG,work_buf)
       READ(work_buf,'(I1)')idebug
c-----------end of han---------------



C Initialize buffers to hold HDF VS and SD file Ids
      DO 20 i = 1, MODFILLEN
         modfil_L1B(i) = 0
         modfil_Geo(i) = 0
         modfil_CldMsk(i) = 0
         modfil_MOD06cd(i) = 0
   20 CONTINUE


C Initialize imaging arrays and parameters for image processing
C   ........ Ping Yang added ..........
C
      DO J = 1, Max_Lines_Per_Img
      DO I = 1, Max_Frames_Per_Line
C
         REFL_1P38_Img(I,J)    = 0.0
         REFL_0P67_Img(I,J)    = 0.0

         REFL_1P24_Img(I,J)    = 0.0

C .....end of PY modification ........

c         REFL_0P67_Img(I,J)    = 0.0
c         REFL_1P24_Img(I,J)    = 0.0
c         REFL_1P38_Img(I,J)    = 0.0

         ZENITH_SOL_Img(I,J)   = 0.0
         ZENITH_VIEW_Img(I,J)  = 0.0
         AZIMUTH_REL_Img(I,J)  = 0.0

         I_LAND_WATER_Img(I,J) = 0
         I_CLOUD_Img(I,J)      = 0
      END DO
      END DO



c rhucek 5/07/01:  added fileoc argument idebug
      ioflag = 'OPEN'
      CALL FILEOC(ioflag, idebug,
     &            handle_trans, handle_QC, FN_L1B,
     &            modfil_L1B, modfil_Geo,
     &            modfil_CldMsk, modfil_MOD06cd)



      att_N = 'Number of Scans'
      dtype = 'INTEGER*4'
      nms = 1
      RTN = GMFIN(modfil_L1B,att_N,dtype,nms,Scans_Per_Granule)
      IF (rtn.ne.mapiok) CALL MODIS_SMF_SETDYNAMICMSG
     &(MODIS_W_GENERIC,'MAPI function 1st GMFIN failed',
     &'FUNCNAME')


      att_N = 'Max Earth View Frames'
      dtype = 'INTEGER*4'
      nms = 1
      RTN = GMFIN(modfil_L1B,att_N,dtype,nms,Max_EV_Frames)
      IF (rtn.ne.mapiok) CALL MODIS_SMF_SETDYNAMICMSG
     &(MODIS_W_GENERIC,'MAPI function 2nd GMFIN failed',
     &'FUNCNAME')


C Call GMTBL to get L1B Scan Metadata, 'EV_Frames'
      bsize = Max_Scan_Per_Granule * 4
      RTN = GMTBL(modfil_L1B,'Level 1B Swath Metadata',' ',
     &      'EV_Frames',0,Scans_Per_Granule,bsize,No_Frames)
      IF (rtn.ne.mapiok) CALL MODIS_SMF_SETDYNAMICMSG
     &(MODIS_W_GENERIC,'MAPI function 1st GMTBL failed',
     &'FUNCNAME')


C Call GMTBL to get L1B Scan Metadata, 'Scan Type'
      bsize = Max_Scan_Per_Granule*4
      rtn = GMTBL(modfil_L1B,'Level 1B Swath Metadata',' ',
     &      'Scan Type',0,Scans_Per_Granule,bsize,Scan_Type)
      IF (rtn.ne.mapiok) CALL MODIS_SMF_SETDYNAMICMSG
     &(MODIS_W_GENERIC,'MAPI function 2nd GMTBL failed',
     &'FUNCNAME')




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

      Total_No_Lines = Scans_Per_Granule * 10
      Buf_Size1 = Max_Frames_Per_Line
      Buf_Size2 = No_Lines_Per_Scan

C ....         modified by Ping Yang April 11, 2001
C         the No_Of_Samples is fixed
C      No_Of_Samples = Max_Frames_Per_Line
      No_Of_Samples = 1354
C......... end the modification 

      No_Of_Lines   = Total_No_Lines


C L_Index_Img below is an index used in the "Do 999" loop to build images
C from individual MODIS scans (one MODIS scan has 10 lines).
      L_Index_Img = 0

C Begin reading MODIS data and ancillary data in HDF-EOS format
      BeginScan_No = 1
      EndScan_No   = Scans_Per_Granule


C Loop over the number of scan swathes.  For each swath, get MODIS
C inputs including sensor reflectances, geolocation, and cloud mask.

      DO 999 Scan_No = BeginScan_No, EndScan_No

         ScanType=Scan_Type(Scan_No)

C "IF - ELSE - ENDIF " added by vlin   06/21/96
         IF (ScanType(1:1).eq.'D'.or.ScanType(1:1).eq.'O') THEN

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


C Modified by Liqun Feb.1998
            IF (idebug .EQ. 1) call Write_QC(handle_QC,
     &                                       Scan_No,
     &                                       Refl_1,
     &                                       Refl_2,
     &                                       Refl_5,
     &                                       Refl_6,
     &                                       Refl_7,
     &                                       Refl_19,
     &                                       Refl_26,
     &                                       SatZen,
     &                                       SolZen,
     &                                       RelAZ,
     &                                       CldMsk,
     &                                       LandSea_Flag,
     &                                       Height)



c rhucek 01/04/99: commented 1 line
c           PROCESS_FLAG=' '

            Data_Size(1) = No_Frames(Scan_No)
            Data_Size(2) = No_Lines_Per_Scan

            DO J = 1, No_Lines_Per_Scan
               L_Index_Img = (Scan_No - 1)*No_Lines_Per_Scan + J

               DO I = 1, No_Frames(Scan_No)

               REFL_0P67_Img(I,L_Index_Img)    = Refl_1(i,j)
               REFL_1P24_Img(I,L_Index_Img)    = Refl_5(i,j)
               REFL_1P38_Img(I,L_Index_Img)    = Refl_26(i,j)

               ZENITH_SOL_Img(I,L_Index_Img)   = SolZen(i,j)
               ZENITH_VIEW_Img(I,L_Index_Img)  = SatZen(i,j)
               AZIMUTH_REL_Img(I,L_Index_Img)  = RelAZ(i,j)

C---Modified by Bo-Cai Gao in May 2004
C---               I_LAND_WATER_Img(I,L_Index_Img) = LandSea_Flag(i,j)
C---               I_CLOUD_Img(I,L_Index_Img)      = CldMsk(i,j)

               I_LAND_WATER_Img(I,L_Index_Img) = Vflag_1(i,j)
               I_CLOUD_Img(I,L_Index_Img)      = Vflag_26(i,j)

               END DO
            END DO

C rhucek 8/26/05:  removed else block
C        ELSE
C           Data_Size(1) = 10
C           Data_Size(2) = No_Lines_Per_Scan
C
         END IF

          DO K=1,2
             Data_Size_For_Granule(K,Scan_No)=Data_Size(K)
          ENDDO

 999  CONTINUE


C End of reading MODIS data and ancillary data in HDF-EOS format

C Read in a pre-calculated lookup table.
C       This table contains transmittances of the 1.375-micron MODIS
C       channel. The transmittances are used in adjusting the transmittance
C       factors for the solar and viewing geometry.
C        CALL READ_TABLE(handle_trans)


C Process MODIS data cube to get cirrus reflectance and contrail index images
C a. - Derive cirrus reflectances for a scene

C.....Added by Ping Yang (11/2/2001)

      DO I = 1, N_PARTITION_X
        IPY_A = IFIX(FLOAT(I-1) * (FLOAT(No_Of_Samples) 
     &     / FLOAT(N_PARTITION_X))) + 1
        IF(I .EQ. N_PARTITION_X ) THEN
         IPY_B =  No_Of_Samples
        ELSE
        IPY_B = IFIX(FLOAT(I) * (FLOAT(No_Of_Samples) 
     &    / FLOAT(N_PARTITION_X)))
        END IF
        IPY_CENT = ( IPY_A + IPY_B ) / 2 + 1 

      DO J = 1, N_PARTITION_Y
        JPY_A = IFIX(FLOAT(J-1) * (FLOAT(No_Of_Lines)
     &     / FLOAT(N_PARTITION_Y))) + 1
        IF(J .EQ. N_PARTITION_Y ) THEN
         JPY_B =  No_Of_Lines
        ELSE
        JPY_B = IFIX(FLOAT(J) * (FLOAT(No_Of_Lines)
     &     / FLOAT(N_PARTITION_Y))) 
        END IF
        JPY_CENT = ( JPY_A + JPY_B ) / 2 + 1 

        PY_SOLAR_ANG(I,J) = ZENITH_SOL_IMG(IPY_CENT,JPY_CENT)
        PY_VIEW_ANG(I,J) = ZENITH_VIEW_IMG(IPY_CENT,JPY_CENT)
        
        END DO
        END DO
         
        CALL CAL_CIRRUS(N_PARTITION_X, N_PARTITION_Y, 
     &    REFL_0P67_Img, REFL_1P38_Img,
     &    PY_SOLAR_ANG, PY_VIEW_ANG, 
     &    CIRRUS_REFL, No_Of_Samples, No_Of_Lines,
     &    Max_Frames_Per_Line,Max_Lines_Per_Img)

C.... end PY's add.....

C rhucek 8/26/05:  Commented call to subroutine PROCESS_CIRRUS_REFLECTANCES
C    by request of Bo-Cai Gao 
C       CALL PROCESS_CIRRUS_REFLECTANCES

C
C
C....................temporary output (for test purpose)....................
C
C---      OPEN(22,file='./zoutput/out_cirrus_refl.dat', 
C---     &        status='unknown',                        
C---     &    form='unformatted',access='direct' ,
C---     & recl=Max_Frames_Per_Line*Max_Lines_Per_Img)
C---          WRITE(22,REC=1) CIRRUS_REFL
C---      CLOSE(22)
C
C
C---      OPEN(23,file='./zoutput/app_refl_1p38.dat', 
C---     &        status='unknown',                        
C---     &    form='unformatted',access='direct',
C---     & recl=Max_Frames_Per_Line*Max_Lines_Per_Img)
C---          WRITE(23,REC=1) REFL_1P38_Img
C---      CLOSE(23)

C---     & recl=Max_Frames_Per_Line*Max_Lines_Per_Img*2)
C
C
C---      OPEN(24,file='./zoutput/app_refl_0p66.dat', 
C---     &        status='unknown',                        
C---     &    form='unformatted',access='direct',
C---     & recl=Max_Frames_Per_Line*Max_Lines_Per_Img)
C---          WRITE(24,REC=1) REFL_0P67_Img
C---      CLOSE(24)

C---      OPEN(44,file='./zoutput/out_ZENITH_SOL_Img.dat',
C---     &        status='unknown',
C---     &    form='unformatted',access='direct',
C---     & recl=Max_Frames_Per_Line*Max_Lines_Per_Img)
C---          WRITE(44,REC=1) ZENITH_SOL_Img
C---      CLOSE(44)

C

C b. - Derive contrail index values for a scene
C
C      Ping Yang's comment (11/2/2001):
C          After discussing with Bo-Cai Gao, we decided that the contrail
C      products should not be generated because the results are not 
C      reliable
C
C
C        CALL PROCESS_CONTRAIL_INDEX    ! commented out by PY (11/2/2001)

C c. - assign values to an array used in QA flags for our cirrus products
        CALL QA_CIRRUS

C---Modified by Bo-Cai Gao in May 2004 for additional assignment of QA
C---     parameters
   
        DO I = 1, No_Of_Samples
         DO J = 1, No_Of_Lines

C---IF VFlag_1(I,j) LE -1, assign I_QA_cirrus(I,J) to 0 ==> bad pixel
C---IF VFlag_26(I,j) LE -1,assign I_QA_cirrus(I,J) to 0 ==> bad pixel
           IF(I_LAND_WATER_Img(I,J) .LE. -1)  I_QA_cirrus(I,J) = 0
c          IF(I_LAND_WATER_Img(I,J) .LE. -1)  cirrus_refl(I,J) = -9.999
           IF(I_LAND_WATER_Img(I,J) .LE. -1)  cirrus_refl(I,J) = -1.9998
           IF(I_CLOUD_Img(I,J)      .LE. -1)  I_QA_cirrus(I,J) = 0
c          IF(I_CLOUD_Img(I,J)      .LE. -1)  cirrus_refl(I,J) = -9.999
          IF(I_CLOUD_Img(I,J)      .LE. -1)  cirrus_refl(I,J) = -1.9998

           IF(REFL_0P67_Img(I,J)    .LE. 0.001)  I_QA_cirrus(I,J) = 0
c          IF(REFL_0P67_Img(I,J)    .LT. 0.001)  cirrus_refl(I,J) = -9.999
           IF(REFL_0P67_Img(I,J)    .LT. 0.001)  cirrus_refl(I,J) = -1.9998
C---           IF(REFL_1P38_Img(I,J)    .LT. 0.000)  I_QA_cirrus(I,J) = 0
C---           IF(REFL_1P38_Img(I,J)    .LT. 0.000)  cirrus_refl(I,J) = -9.999

         END DO
        END DO

C---      OPEN(25,file='./zoutput/out_I_QA_cirrus.dat',
C---     &        status='unknown',
C---     &    form='unformatted',access='direct',
C---     & recl=Max_Frames_Per_Line*Max_Lines_Per_Img)
C---          WRITE(25,REC=1) I_QA_cirrus
C---      CLOSE(25)

C---      OPEN(26,file='./zoutput/out_Vflag_1.dat',
C---     &        status='unknown',
C---     &    form='unformatted',access='direct',
C---     & recl=Max_Frames_Per_Line*Max_Lines_Per_Img)
C---          WRITE(26,REC=1) I_LAND_WATER_Img
C---      CLOSE(26)

C---      OPEN(27,file='./zoutput/out_Vflag_26.dat',
C---     &        status='unknown',
C---     &    form='unformatted',access='direct',
C---     & recl=Max_Frames_Per_Line*Max_Lines_Per_Img)
C---          WRITE(27,REC=1) I_CLOUD_Img
C---      CLOSE(27)

C--------0---------0---------0---------0---------0---------0---------072
C Write the cirrus reflectance and contrail index images into an
C       output HDF file
C
            L_Index_Img = 0
      DO 1000 Scan_No = BeginScan_No, EndScan_No

C rhucek 8/26/05: Added 2 lines: The assignment of variable ScanType
C  and the IF test on ScanType
         ScanType = Scan_Type(Scan_No)
   
         IF (ScanType(1:1).eq.'D'.or.ScanType(1:1).eq.'O') THEN

            DO K=1,2
               Data_Size(K)=Data_Size_For_Granule(K,Scan_No)
            ENDDO
   
            DO J = 1, No_Lines_Per_Scan
               L_Index_Img = (Scan_No - 1)*No_Lines_Per_Scan + J
   
               DO I = 1, No_Frames(Scan_No)
                  cirrus_refl_sub(I,J) = cirrus_refl(I,L_Index_Img)

C the following one line commented out by PY (11/2/2001)
C                 I_CONTRAIL_sub(I,J)  =  I_CONTRAIL(I,L_Index_Img)
                  I_QA_cirrus_sub(I,J) = I_QA_cirrus(I,L_Index_Img)
C
C---Temp fixing of a bug in the code, so that the correct index
C    with values between 0 and 3 (0 = bad pixel; 1 = no cirrus;
C    2 = cirrus pixel; 3 = contrail pixel)

                  I_CONTRAIL_sub(I,J)  = I_QA_cirrus_sub(I,J)
C---End of Temp Fixing

               END DO
            END DO

            CALL MOD_cirrus_PutData(Scan_No,modfil_MOD06cd,Data_Size,
     2            Buf_Size1,Buf_Size2,cirrus_refl_sub,I_CONTRAIL_sub)

C rhucek 8/26/05: Added END IF
         END IF
1000  CONTINUE

C Finish the writing of cirrus reflectance and rail index images
C       into an output HDF file
C--------0---------0---------0---------0---------0---------0---------072

c-----Dr. Wei Han, han@han.nrl.navy.mil--------------------
c\* Adding metadata to mod06cd.qc file

      rtn = Set_CoreMetadata_QC(LRN_MCFQC,LRN_QCMET)
      if (rtn.ne.mapiok) then
         call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,
     2   'Create metadata file to mod06CD.qc failed' // char(10) //
     4   'Operator Action: Stage non-corrupt version' // char(10) //
     5   'or correct the PCF reference to QC file MCF.  Rerun PGE.',
     3   'MOD_PR06CD_V2')
      else
         call modis_smf_setdynamicmsg(MODIS_S_SUCCESS,
     2   'Create metadata file to mod06CD.qc without error',
     3   'MOD_PR06CD_V2')
      end if

c rhucek 05/07/01:  added debug condition on call to Set_CoreMetadata_QC
 
      IF (idebug .eq. 1) THEN
         rtn = Set_CoreMetadata_QC(LRN_MCFQC,LRN_QCMET)
 
         if (rtn.ne.mapiok) then
            call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,
     2      'Create metadata file to mod06CD.qc failed' // char(10) //
     4      'Operator Action: Stage non-corrupt version' // char(10) //
     5      'or correct the PCF reference to QC file MCF.  Rerun PGE.',
     3      'MOD_PR06CD_V2')
         else
            call modis_smf_setdynamicmsg(MODIS_S_SUCCESS,
     2      'Create metadata file to mod06CD.qc without error',
     3      'MOD_PR06CD_V2')
         end if
      END IF


c------------end of han-----------------


C Close the output file in HDF-EOS format
c rhucek 5/07/01:  added fileoc argument idebug
      ioflag = 'CLOSE'
      CALL FILEOC(ioflag, idebug,
     &            handle_trans, handle_QC, FN_L1B,
     &            modfil_L1B, modfil_Geo,
     &            modfil_CldMsk, modfil_MOD06cd)



      RETURN
      END



C-----------------------------------------------------------------------
C
      SUBROUTINE PROCESS_CIRRUS_REFLECTANCES
C
C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION: This subroutine derives cirrus reflectances in the
C              0.4 - 1.0 micron region using the special 1.375-micron
C              MODIS cirrus detection channel.
C
C !Input Parameters:
C   INTEGER     Refl_1         MODIS reflectance at channel 1
C   INTEGER     Refl_5         MODIS reflectance at channel 5
C   INTEGER     Refl_26        MODIS reflectance at channel 26
C   INTEGER     SolZen         Solar zenith angle
C   INTEGER     SatZen         Satellite view zenith angle
C   INTEGER     RelAZ          Relative azimuth angle between the Sun
C                              and the satellite
C
C !Output Parameters:
C   INTEGER     Cirrus_Refl    Two dimensional (2-D) array for passing
C                               retrieved cirrus reflectance data
C
C !REVISION HISTORY:
C
C !TEAM-UNIQUE HEADER: - NOT YET DEFINED
C
C !REFERENCES AND CREDITS
C    First version created by Dr. Bo-Cai Gao
C    (gao@rsd.nrl.navy.mil)
C                             1997 09 16
C
C !DESIGN NOTES: data and variables are passed with common statements
C               contained in "COMMONS_cirrus.inc"
C
C !Functions and Subroutines:
C      HUNT
C      FIND_CIRRUS_REFLECTANCE
C
C-----------------------------------------------------------------------
C
C---Note: The at launch version of the cirrus reflectance algorithm
C         is very simplified. The post-launch version of the algorithm
C         will be improved after the analysis of real MODIS data.
C
C !END-------------------------------------------------------------------

      INCLUDE 'COMMONS_cirrus.inc'
      SAVE

c     INTEGER I, J, K
c     REAL    PI, VAP, DEG_TO_RAD
c     DATA PI, DEG_TO_RAD, J_GH2O /3.1415926, 0.0174533, 1/

      INTEGER I
      REAL    DEG_TO_RAD
      DATA DEG_TO_RAD, J_GH2O /0.0174533, 1/

      DO J = 1, No_Of_Lines
         DO I = 1, No_Of_Samples

C ... Ping Yang (3/15/2001):
C     Here, the look_up process for H2O transmittance that was 
C      written by B.-C. Gao has been removed.
C     New algorithm is used for deriving cirrus clouds 
C     
C.....end PY's modification

C   ...... Ping Yang made a modification (4/18/2001)

       IF((CIRRUS_REFL(I,J).GT.REFL_1P24_Img(I,J)) .AND.
     &     (REFL_1P24_Img(I,J).GT.0.0) )  
     &                 CIRRUS_REFL(I,J)   =    REFL_1P24_Img(I,J)
C
C       IF(CIRRUS_REFL(I,J).GT.REFL_1P24_Img(I,J))
C     &                 CIRRUS_REFL(I,J)   =    REFL_1P24_Img(I,J)

C Taking care of BAD_FLAG for REFL_1P38_Img
C  Ping Yang comments out the following two lines (4/18/2001)
C       IF(IFIX(REFL_1P38_Img(I,J)).EQ.IFIX(BAD_FLAG))
C     &                 CIRRUS_REFL(I,J)   =    0.0


        END DO
      END DO


      RETURN
      END
C
C***********************************************************************
C
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
C     Explicit data type.
C
C !TEAM-UNIQUE HEADER: - NOT YET DEFINED
C
C !REFERENCES AND CREDITS
C     WRITTEN BY: Bo-Cai Gao
C
C !DESIGN NOTES: none
C
C !END------------------------------------------------------------------
C
      INTEGER J, N
      REAL DLT, FJ, FJ_P1, G, T
      REAL X, T_1P38
      DIMENSION X(N), T_1P38(N)
C
      IF( (J.GT.0).AND.(J.LT.N) ) THEN

         DLT   = X(J+1) - X(J)
         FJ    = ( X(J+1) - G ) / DLT
         FJ_P1 = ( G - X(J) )   / DLT

         T     = FJ * T_1P38(J) + FJ_P1 * T_1P38(J+1)

      ELSE

         IF(J.LE.0) T = T_1P38(1)
         IF(J.GE.N) T = T_1P38(N)

      END IF

C
      RETURN
      END
C
C*******************************************************************************

      SUBROUTINE HUNT(XX,N,X,JLO)

C-----------------------------------------------------------------------
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
C !END------------------------------------------------------------------
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

C-----------------------------------------------------------------------
C
      SUBROUTINE PROCESS_CONTRAIL_INDEX
C
C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION: This subroutine derives an index of aircraft
C              contrails from the image of the 1.375-micron MODIS channel.
C
C !Input Parameters:
C   INTEGER     Refl_26        MODIS reflectance at channel 26
C
C !OUTPUT PARAMETERS:
C   BYTE        I_contrail     Two dimensional (2-D) array for passing
C                              contrail index data. Index 0 refers to non-contrail
C                              pixel, and index 1 refers to contral pixel.
C
C !REVISION HISTORY:
C
C !TEAM-UNIQUE HEADER: - NOT YET DEFINED
C
C !REFERENCES AND CREDITS
C    First version created by Dr. Bill Ridgway
C    (ridgway@climate.gsfc.nasa.gov)
C                             1997 11 03
C
C !DESIGN NOTES: data and variables are passed with common statements
C               contained in "COMMONS_cirrus.inc"
C
C   Functions and Subroutines:
C     Find_Contrail_Pixels
C
C !END------------------------------------------------------------------

      INCLUDE 'COMMONS_cirrus.inc'
      INTEGER I, J

      REAL      REFL_1P38_Img_QTR
      DIMENSION REFL_1P38_Img_QTR(Max_Frames_Per_Line_HF,
     @                            Max_Lines_Per_Img_HF)

      INTEGER   I_CONTRAIL_QTR
      DIMENSION I_CONTRAIL_QTR   (Max_Frames_Per_Line_HF,
     @                            Max_Lines_Per_Img_HF)
      SAVE

C---Process 1st quarter of the 1.38-micron image for detecting contrails
C- Note: Because Bill Ridgway's contrail detection algorithm expects
C-       input images with digital numbers between about 0 and 255, while the
C-       MODIS Level_1B reflectances values are typically in the range of
C        0 - 1., we multiply 255 when passing images to Ridgway's contrail
C        algorithm.

      DO J = 1, Max_Lines_Per_Img_HF
         DO I = 1, Max_Frames_Per_Line_HF
            REFL_1P38_Img_QTR(I,J) = REFL_1P38_Img(I,J) * 255.
         END DO
      END DO

      CALL Find_Contrail_Pixels ( Max_Frames_Per_Line_HF,
     +                            Max_Lines_Per_Img_HF,
     +                            REFL_1P38_Img_QTR, I_CONTRAIL_QTR)

      DO J = 1, Max_Lines_Per_Img_HF
         DO I = 1, Max_Frames_Per_Line_HF
            I_CONTRAIL(I,J) = I_CONTRAIL_QTR(I,J)
         END DO
      END DO

C---Process 2nd quarter of the 1.38-micron image for detecting contrails

      DO J = 1, Max_Lines_Per_Img_HF
         DO I = Max_Frames_Per_Line_HF + 1, Max_Frames_Per_Line
            REFL_1P38_Img_QTR(I - Max_Frames_Per_Line_HF, J)
     +         =                             REFL_1P38_Img(I,J) * 255.
         END DO
      END DO

      CALL Find_Contrail_Pixels ( Max_Frames_Per_Line_HF,
     +                            Max_Lines_Per_Img_HF,
     +                            REFL_1P38_Img_QTR, I_CONTRAIL_QTR)

      DO J = 1, Max_Lines_Per_Img_HF
         DO I = Max_Frames_Per_Line_HF + 1, Max_Frames_Per_Line
            I_CONTRAIL(I,J)
     +          = I_CONTRAIL_QTR(I - Max_Frames_Per_Line_HF, J)
         END DO
      END DO

C---Process 3rd quarter of the 1.38-micron image for detecting contrails

      DO J = Max_Lines_Per_Img_HF + 1, Max_Lines_Per_Img
         DO I = 1, Max_Frames_Per_Line_HF
            REFL_1P38_Img_QTR(I,J - Max_Lines_Per_Img_HF)
     +         =                             REFL_1P38_Img(I,J) * 255.
         END DO
      END DO

      CALL Find_Contrail_Pixels ( Max_Frames_Per_Line_HF,
     +                            Max_Lines_Per_Img_HF,
     +                            REFL_1P38_Img_QTR, I_CONTRAIL_QTR)

      DO J = Max_Lines_Per_Img_HF + 1, Max_Lines_Per_Img
         DO I = 1, Max_Frames_Per_Line_HF
            I_CONTRAIL(I,J)
     +         =  I_CONTRAIL_QTR(I, J - Max_Lines_Per_Img_HF)
         END DO
      END DO

C---Process 4th quarter of the 1.38-micron image for detecting contrails

      DO J = Max_Lines_Per_Img_HF + 1, Max_Lines_Per_Img
         DO I = Max_Frames_Per_Line_HF + 1, Max_Frames_Per_Line
            REFL_1P38_Img_QTR(I - Max_Frames_Per_Line_HF,
     +                        J - Max_Lines_Per_Img_HF)
     +         =                             REFL_1P38_Img(I,J) * 255.
         END DO
      END DO

      CALL Find_Contrail_Pixels ( Max_Frames_Per_Line_HF,
     +                            Max_Lines_Per_Img_HF,
     +                            REFL_1P38_Img_QTR, I_CONTRAIL_QTR)

      DO J = Max_Lines_Per_Img_HF + 1, Max_Lines_Per_Img
         DO I = Max_Frames_Per_Line_HF + 1, Max_Frames_Per_Line
            I_CONTRAIL(I,J)
     +         =  I_CONTRAIL_QTR(I - Max_Frames_Per_Line_HF,
     +                           J - Max_Lines_Per_Img_HF)
         END DO
      END DO

C---End of processing all 4 quarters.

      RETURN
      END

C-----------------------------------------------------------------------
C
      SUBROUTINE QA_CIRRUS
C
C--------0---------0---------0---------0---------0---------0---------072
C !F77
C
C !DESCRIPTION: This subroutine is to fill an I*2 array to be used in
C              QA flags.
C
C !Input Parameters:
C   INTEGER     Refl_26        MODIS reflectance at channel 26
C   INTEGER     Cirrus_Refl    Two dimensional (2-D) array for passing
C                               retrieved cirrus reflectance data.
C   BYTE        I_contrail     Two dimensional (2-D) array for passing
C                              contrail index data. Index 0 refers to non-contrail
C                              pixel, and index 1 refers to contral pixel.
C
C !OUTPUT PARAMETERS:
C   INTEGER     I_QA_cirrus    Two dimensional (2-D) array for passing
C                              quality assurance data for cirrus reflectance
C                              and contrail index. Index 0 refers to bad data,
C                              index 1 refers to non-cirrus pixel, index 2 refers
C                              to cirrus pixel, and index 3 refers to contral pixel.
C
C !REVISION HISTORY:
C
C !TEAM-UNIQUE HEADER: - NOT YET DEFINED
C
C !REFERENCES AND CREDITS
C    First version created by Bo-Cai Gao
C    (gao@rsd.nrl.navy.mil)
C                             1997 11 05
C
C !DESIGN NOTES: data and variables are passed with common statements
C               contained in "COMMONS_cirrus.inc"
C
C !END---0---------0---------0---------0---------0---------0---------072

      INCLUDE 'COMMONS_cirrus.inc'
      SAVE

      INTEGER I, J

      DO J = 1, No_Of_Lines
        DO I = 1, No_Of_Samples

          I_QA_cirrus(I,J) = 0

         END DO
      END DO


      DO J = 1, No_Of_Lines
        DO I = 1, No_Of_Samples

          IF(cirrus_refl(I,J).GE.cirrus_refl_min) I_QA_cirrus(I,J) = 2
C the following one line is commented out by PY (11/2/2001)
C          IF(I_CONTRAIL(I,J).EQ.1)                I_QA_cirrus(I,J) = 3
          IF(cirrus_refl(I,J).LT.cirrus_refl_min) I_QA_cirrus(I,J) = 1

C... Ping Yang comments out the following two lines (4/18/2001)
C          IF(IFIX(REFL_1P38_Img(I,J)).EQ.IFIX(BAD_FLAG))
C     @                                            I_QA_cirrus(I,J) = 0

C--- Modified by Bo-Cai Gao in May 2004 to "LT  0.0"
C---          IF(REFL_1P38_Img(I,J).LT. -1.0)
C---     &                                            I_QA_cirrus(I,J) = 0
          IF(REFL_1P38_Img(I,J).LT.  0.0)
     &                                            I_QA_cirrus(I,J) = 0

         END DO
      END DO

      RETURN
      END
C
C--------0---------0---------0---------0---------0---------0---------072


C added by Liqun Feb.1998
      SUBROUTINE Write_QC(handle_QC,
     &                    Scan_No,
     &                    Refl_1,
     &                    Refl_2,
     &                    Refl_5,
     &                    Refl_6,
     &                    Refl_7,
     &                    Refl_19,
     &                    Refl_26,
     &                    SatZen,
     &                    SolZen,
     &                    RelAZ,
     &                    CldMsk,
     &                    LandSea_Flag,
     &                    Height)

      IMPLICIT NONE
      SAVE

      INCLUDE 'cirrus.inc'
      INCLUDE 'COMMONS_cirrus.inc'

C-----------------------------------------------------------------------
C !F77
C
C !Description: This subroutine writes cirrus qc date.
C
C !Input Parameters:
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
      INTEGER Scan_No,handle_QC,i,j
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

