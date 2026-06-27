      SUBROUTINE MOD_PR05_V2(modis_flag,idebug,MODFIL_MOD05,
     1                       MinSolarZenithAngle,RTN_NCEP)
C*************************************************************************
      IMPLICIT NONE

      INCLUDE 'COMMONS.inc'
      INCLUDE 'mapi.inc'
      INCLUDE 'L1B_Reader_V2.1.inc'
      INCLUDE 'PGS_PC.f'
      INCLUDE 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !Description: This program is to derive column water vapor values
C                above clear land surfaces and above clouds from MODIS
C                channels near 1 micron.
C
C !Input Parameters:
C      MODIS radiance data, solar and view zenith angles, cloud mask,
C      surface elevation, surface type (land or water).
C
C !Output Parameters:
C      Vertical column water vapor from ground to space or from cloud to
C      space.
C
C !Revision History:
C $Log: MOD_PR05.f,v $
c Revision 1.7  1996/12/31  16:18:01  jguu
c The size for the character string which holds the filename
c for the L1B granule is changed to PGSd_PC_FILE_PATH_MAX.
c
c Revision 1.6  1996/11/05  22:25:52  jguu
c The interface to ecsMeta_data is changed.
c
c Revision 1.5  1996/11/01  19:19:39  jguu
c Remove all the GO TO construct.
c Code meets MODIS coding standard.
c The code now can handle granule as large as 200 scans.
c
c Revision 1.4  1996/10/04  16:39:14  jguu
c New delivery form ST after the implementation of the
c ECS Core Metadata.
c
c Revision 1.3  1996/08/29  19:58:41  jguu
c The file name for the include file COMMONS_INC.f is changed
c to COMMONS.inc.
c
c Revision 1.2  1996/07/26  20:03:55  jguu
c The delimiters for character strings are changed from
c quotation marks to apostrophes.
c
c Revision 1.1  1996/07/08  18:32:08  vlin
c Initial revision
c
c Revision 1.1  1996/06/24  19:24:25  gao
c Initial revision
c
c Revision 1.12  1995/11/21  14:38:21  vlin
c Fix a bug around DO 999 loop
c
c Revision 1.11  1995/11/20  23:52:51  rhucek
c 1 - Identified code changes/additions made by rhucek
c 2 - Moved calls to FILEOC, SOLAR_IRRADIANCE, TRANSM_H2O_RATIOS, and
c     WEIGHTING_TABLE inside MODIS processing block.
c 3 - Updated design notes under PROGRAM VAPOR
c
c Revision 1.10  1995/11/02  18:57:25  vlin
c Get rid of tabs
c
c Revision 1.9  1995/11/02  18:04:33  vlin
c added declaration J on line #978
c
c Revision 1.8  1995/11/01  15:24:38  vlin
c Before code walkthrough
c
c Revision 1.7  1995/10/18  19:30:11  vlin
c Move all calls to subroutine "MODIS_IO_GEN_OPEN" & "MODIS_IO_GEN_CLOSE"
c to "mod05_fileoc.f"
c
c Revision 1.6  1995/10/18  17:16:43  vlin
c Wrote 2 more SDSs to the HDF file
c
c Revision 1.4  1995/10/06  14:26:47  vlin
c added two calls to subroutine mod05_fileoc
c
c Revision 1.3  1995/09/05  19:31:31  ding
c Updated with more comment lines.
c
c Revision 1.2  1995/09/01  19:38:51  ding
c added HDF file output
c
c Revision 1.1  1995/08/03  13:51:55  modiscm
c Initial revision
c
C !TEAM-UNIQUE HEADER:
C
C   This software was developed by the MODIS Science Team
C   for the National Aeronautics and Space Administration,
C   Goddard Space Flight Center.
C
C !REFERENCES and CREDITS:
C
C   Written by
C   Dr. Bo-Cai Gao       10/28/94
C   gao@imagecube.gsfc.nasa.gov
C
C !DESIGN NOTES:
C
C   At present, this program is set up to process one granule of MODIS
C   sensor data. A granule consists of 100 MODIS scan swathes, each
C   containing 1354 1-km pixels in the scan direction and 10 pixels
C   along the spacecraft flight direction.  This program also requires
C   6 atmospheric water vapor transmittance files (TRANSM_H2O.MDL_1,
C   TRANSM_H2O.MDL_2, TRANSM_H2O.MDL_3, TRANSM_H2O.MDL_4,
C   TRANSM_H2O.MDL_5, TRANSM_H2O.MDL_6) as input.
C
C Externals:
C   Named Constants:
C           MAPIOK                          (mapi.inc)
C           Max_Frames_Per_Line             (COMMONS.inc)
C           No_Lines_Per_Scan              (COMMONS.inc)
C
C Internals:
C
C   Functions and Subroutines:
C          TRANSM_H2O_RATIOS
C          SOLAR_IRRADIANCE
C          WEIGHTING_TABLE
C          FILEOC
C          GetModisDat_MOD05
C          MODIS_SMF_SETDYNAMICMSG
C          PROCESS_MODIS_DATA
C          MOD05_PutData
C
C !END
C-----------------------------------------------------------------------

c rhucek 01/19/98:  added parameter NO_QA_BYTES and FUNCNAME
      CHARACTER*(*) FUNCNAME
      PARAMETER    (FUNCNAME = 'MOD_PR05_V2')

      CHARACTER*5   ioflag
      CHARACTER*40  att_N,dtype
      CHARACTER*(PGSd_PC_FILE_PATH_MAX) FN_L1B

      INTEGER    NO_QA_BYTES
      PARAMETER (NO_QA_BYTES = 1)

      INTEGER Max_EV_Frames, Max_Scan_Per_Granule,No_PSA,RTN_NCEP
      PARAMETER (Max_Scan_Per_Granule = 550, No_PSA=15)

      CHARACTER*4  Scan_Type(Max_Scan_Per_Granule)
      CHARACTER*(PGSd_MET_GROUP_NAME_L) groups(PGSd_MET_NUM_OF_GROUPS)
      CHARACTER*(MAX_ECS_NAME_L-1) HDFAttNms(PGSd_MET_NUM_OF_GROUPS)

CDAC: for create_mod05
      character*512 usrlog
      INTEGER mrtn,create_mod05
CDAC: end
c
c  Parameters for cloud array
      INTEGER   Buf_cldmsk,  Buf_cldmsk_QA
      PARAMETER(Buf_cldmsk=6,Buf_cldmsk_QA=10)

c rhucek 01/15/98:  Added BYTE type declaration of QA_TEMP
        BYTE  QA_TEMP,
     1        Cloud(Buf_cldmsk,Max_Frames_Per_Line, No_Lines_Per_Scan),
     2        QA_Cloud(Buf_cldmsk_QA,Max_Frames_Per_Line,
     3        No_Lines_Per_Scan)

      INTEGER i, idebug, j, Rank, RTN, nms, Scans_Per_Granule,
     &        Total_No_Lines, BeginScan_No, EndScan_No, Scan_No,
     &        No_of_Scans, handle_weight, handle_QC,
     &        Buf_Size_QA, Buf_Size1, Buf_Size2,
     &        Modfil_Geo(MODFILLEN),Modfil_CldMsk(MODFILLEN),
     &        Modfil_MOD05(MODFILLEN),MODFIL_ANC(MODFILLEN),
     &        Modfil_L1B(MODFILLEN),
     &        handle_trans(6), Dim_Size(2), Data_Size(2),
     &        Data_Size_Anc(2),
     &        No_Frames(Max_Scan_Per_Granule),bsize
      INTEGER pgs_pc_getreference,GetModisDat_MOD05,ReadAnc_MOD05

c rhucek 01/15/98:  removed integer type declaration of variable
c                   QA_TEMP and replaced with byte type
      INTEGER CldMsk      (Max_Frames_Per_Line, No_Lines_Per_Scan),
     2        LandSea_Flag(Max_Frames_Per_Line, No_Lines_Per_Scan),
     3        CldMsk_QA   (Max_Frames_Per_Line, No_Lines_Per_Scan),
     4        QA_NIR      (Max_Frames_Per_Line, No_Lines_Per_Scan),
     5        SumQual_Flag,Qual_Flag,Ret_method,Suf_Flag,
     6        I_ret,I_total

*/  Modified by JC Guu  10/04/96
*/  The declaration for NumHandles is added.

      INTEGER  NumHandles

      LOGICAL modis_flag,done_HDF_EOS

      REAL Lat     (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Lon     (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     SatZen  (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     SolZen  (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     RelAz   (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Refl_2  (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Refl_5  (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Refl_17 (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Refl_18 (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Refl_19 (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Height  (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Sfc_Temp(Max_Frames_Per_Line, No_Lines_Per_Scan)
C
      BYTE UN_2    (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     UN_5    (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     UN_17   (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     UN_18   (Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     UN_19   (Max_Frames_Per_Line, No_Lines_Per_Scan)

C
      BYTE Vflag_2(Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Vflag_5(Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Vflag_17(Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Vflag_18(Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Vflag_19(Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Buf_Un(Max_Frames_Per_Line, No_Lines_Per_Scan),
     &     Buf_Sa(Max_Frames_Per_Line, No_Lines_Per_Scan)

      REAL QA_Metadata_MOD05(No_PSA)
      REAL MinSolarZenithAngle,SolarZenithAngleZEPS
      real Lat_center,Lon_center,sfctmp,ugrd,vgrd,sfctemp_clim,
     &     pwat,ozone

      PARAMETER(SolarZenithAngleZEPS=72.000001)

      SAVE

c rhucek 01/09/98:  added "Begin" status message
      usrlog = char(10)
     1// '------------------------------------------------------------'
     2// char(10)
     3// 'Begin Water Vapor Retrieval Algorithm w/o Aerosol Correction'
     4// char(10)
     5// '------------------------------------------------------------'

      Call MODIS_SMF_SETDYNAMICMSG(MODIS_M_GENERIC,usrlog,FUNCNAME)

C
C Initialize variables
C     idebug = 1

      sfctemp_clim=288.0
      I_ret=0
      I_total=0
      SumQual_Flag=0
      Qual_Flag=0
      Ret_Method=0
      Suf_Flag=0

C-----------------------------------------------------------------------
C If modis_flag = .true., process MODIS Data; otherwise, open
C "non-MODIS" data files.  Only MODIS option is presently available.
C-----------------------------------------------------------------------

C      modis_flag = .true.
      done_HDF_EOS= .false.
      IF (modis_flag) THEN

C      CALL MOD_PRO5_cre(done_HDF_EOS)
C      IF(done_HDF_EOS) write(*,*)'Finished creation'

C Initialize buffers to hold HDF VS and SD file Ids

         DO 20 i = 1, MODFILLEN
            modfil_L1B(i) = 0
            modfil_Geo(i) = 0
            modfil_CldMsk(i) = 0
            modfil_Anc(i) = 0
            modfil_MOD05(i) = 0
   20    CONTINUE

CDAC: creation of mod05
*/  Define the output swath file
      rtn=create_mod05()
      usrlog="Output product swath file defined - swath_mod5"

c rhucek 01/09/98:  replaced one line
c     call ckstatus_f_mod05(mrtn,usrlog,funcname)
      call ckstatus_f_mod05(rtn,usrlog,FUNCNAME)
CDAC: end


c.....process as night
      IF (MinSolarZenithAngle.GE.SolarZenithAngleZEPS) THEN
         Call MOD05_Process_Night()

c..... process as day
       ELSE
         ioflag = 'OPEN'

         CALL FILEOC(ioflag,modis_flag,handle_trans,
     &               handle_weight,handle_QC,FN_L1B,
     &               modfil_L1B,modfil_Geo,modfil_CldMsk,
     &               modfil_Anc,modfil_MOD05,groups,
     &               HDFAttNms,NumHandles)
c        WRITE(*,*)'Pass FILEOC'

         CALL TRANSM_H2O_RATIOS(handle_trans)
c        WRITE(*,*)'Pass TRANSM_H2O_RATIOS'

         CALL SOLAR_IRRADIANCE
c        WRITE(*,*)'Pass SOLAR_IRRADIANCE'

         CALL WEIGHTING_TABLE(handle_weight)
c        WRITE(*,*)'Pass WEIGHTING_TABLE'

C---       OPEN(61,file='QA_REFL_0p66_binary',status='unknown',
C---     &    form='unformatted',access='direct',recl=1354)

C
C Put behind opening of files
C modified by DAC, 9/19/96, to include ECS Core Metadata
C

C Call GMFIN to get  L1B file metadata, 'Number of Scans'
         att_N = 'Number of Scans'
         dtype = 'INTEGER*4'
         nms = 1
         RTN = GMFIN(modfil_L1B,att_N,dtype,nms,Scans_Per_Granule)

         IF (rtn.ne.mapiok) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'MAPI function 1st GMFIN failed','vapor.f')
c        WRITE(*,*)'Pass Scans_Per_Granule',Scans_Per_Granule

C Add 24 lines by vlin 06/24/96
C Call GMFIN to get L1B file metadata, 'Max Earth View Frames'
         att_N = 'Max Earth View Frames'
         dtype = 'INTEGER*4'
         nms = 1
         RTN = GMFIN(modfil_L1B,att_N,dtype,nms,Max_EV_Frames)

         IF (rtn.ne.mapiok) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'MAPI function 2nd GMFIN failed','vapor.f')
C         WRITE(*,*)'Pass Max_EV_Frames',Max_EV_Frames
CDAC
C Call GMFIN to get L1B file metadata, 'Dead MODIS Detectors'
C         att_N = 'Dead MODIS Detectors'
C         dtype = 'INTEGER*4'
C         nms = 490
C         RTN = GMFIN(modfil_L1B,att_N,dtype,nms,Dead_Detectors)
C
C         IF (rtn.ne.mapiok) CALL MODIS_SMF_SETDYNAMICMSG
C     &   (MODIS_W_GENERIC,'MAPI function 3rd GMFIN failed','vapor.f')
C
*/  Modified by JC Guu  11/01/96
*/  The variable Scans_Per_Granule in the argument list
*/  replaces the value 100.

C Call GMTBL to get L1B Scan Metadata, 'EV_Frames'
         bsize = Max_Scan_Per_Granule * 4
         RTN = GMTBL(modfil_L1B,'Level 1B Swath Metadata',' ',
     &         'EV_Frames',0,Scans_Per_Granule,bsize,No_Frames)

         IF (rtn.ne.mapiok) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'MAPI function 1st GMTBL failed','vapor.f')
C         WRITE(*,*)'Pass No_Frames'

C  Call GMTBL to get L1B Scan Metadata, 'Scan Type'
         bsize = Max_Scan_Per_Granule*4
         rtn = GMTBL(modfil_L1B,'Level 1B Swath Metadata',' ',
     &         'Scan Type',0,Scans_Per_Granule,bsize,Scan_Type)

         IF (rtn.ne.mapiok) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'MAPI function 2nd GMTBL failed','vapor.f')
C         WRITE(*,*)'Pass Scan_Type',Scan_Type

         Total_No_Lines = Scans_Per_Granule * 10

         Buf_Size1 = Max_Frames_Per_Line
         Buf_Size2 = No_Lines_Per_Scan

         BeginScan_No = 1
         EndScan_No   = Scans_Per_Granule

         IF (idebug .EQ. 1) THEN
            write(handle_QC,'(" ")')
            write(handle_QC,'("List of EV_Frames")')
            write(handle_QC,'(10i5)') No_Frames
         END IF

C Loop over the number of scan swathes.  For each swath, get MODIS
C inputs including sensor reflectances, geolocation, and cloud mask.


         DO 999 Scan_No = BeginScan_No, EndScan_No

c         WRITE(*,*) Scan_No
C
C "IF - ELSE - ENDIF " added by vlin   06/21/96
C
      IF (Scan_Type(Scan_No).eq.'D'.or.Scan_Type(Scan_No).eq.'O') THEN

            RTN=GetModisDat_MOD05(
     &        Modfil_L1B,Modfil_Geo,Modfil_CldMsk,FN_L1B,
     &        Scan_No,Buf_Size1,Buf_Size2,Data_Size,
     &        Lat,Lon,SatZen,SolZen,RelAz,Height,
     &        Refl_2,Refl_5,Refl_17,Refl_18,Refl_19,
     &        Un_2,Un_5,Un_17,Un_18,Un_19,Vflag_2,Vflag_5,
     &        Vflag_17,Vflag_18,Vflag_19,Buf_Un,Buf_Sa,
     &        CldMsk,LandSea_Flag,Cloud,QA_Cloud,CldMsk_QA)

            IF (rtn.ne.0) CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     &         'Call to GetModisDat_MOD05 failed','vapor.f')

C             RTN=ReadAnc_MOD05(MODFIL_ANC,Buf_Size1,Buf_Size2,
C     &          Scan_No,Data_Size_Anc,sfc_temp)
C
C          IF (RTN.NE.0) CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,
C     &         'ReadAnc_MOD05 Failed','vapor.f')

C Added 1 line by vlin   06/21/96

            Data_Size(1) = No_Frames(Scan_No)
            Data_Size(2) = No_Lines_Per_Scan

            DO J = 1, No_Lines_Per_Scan
            DO I = 1, No_Frames(Scan_No)
              Lat_center=Lat(I,J)
              Lon_center=Lon(I,J)
              IF(RTN_NCEP.EQ.0) THEN
C rhucek 11/22/02: Replaced call to read_anc_data with GetAncData_PGE04
                CALL GetAncData_PGE04(Lat_center,Lon_center,sfctmp,
     &                                ugrd,vgrd,pwat,ozone)
                sfc_temp(i,j)=sfctmp
              ELSE
                sfc_temp(i,j)=sfctemp_clim
              ENDIF
            END DO
            END DO

            IF (idebug .EQ. 1) THEN

               write(handle_QC,'(/,"Scan Number ",i3)') Scan_No
               write(handle_QC,*)


*/  Modified by JC Guu 07/25/96
*/  The delimiters for character strings are changed from
*/  quotation marks (") to apostrophes (') to comply to
*/  the ANSI FORTRAN 77 standard.  There are 26
*/  locations where this correction occurrs.


            write(handle_QC,'(A)') 'first 11 reflectances for ' //
     &                             'band 2 for 10 lines'
            write(handle_QC,'(11f11.5)') ((Refl_2(i,j),i=1,11),
     &                                        j=1,10)

            write(handle_QC,*)
            write(handle_QC,'(A)') 'first 11 reflectances for ' //
     &                             'band 5 for 10 lines'
            write(handle_QC,'(11f11.5)') ((Refl_5(i,j),i=1,11),j=1,10)

            write(handle_QC,*)
            write(handle_QC,'(A)') 'first 11 reflectances for ' //
     &                             'band 17 for 10 lines'
            write(handle_QC,'(11f11.5)') ((Refl_17(i,j),i=1,11),j=1,10)

            write(handle_QC,*)
            write(handle_QC,'(A)') 'first 11 reflectances for ' //
     &                             'band 18 for 10 lines'
            write(handle_QC,'(11f11.5)') ((Refl_18(i,j),i=1,11),j=1,10)

            write(handle_QC,*)
            write(handle_QC,'(A)') 'first 11 reflectances for ' //
     &                             'band 19 for 10 lines'
            write(handle_QC,'(11f11.5)') ((Refl_19(i,j),i=1,11),j=1,10)

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
            write(handle_QC,'(A)') 'first 11 Latitude for 10 lines'
            write(handle_QC,'(11f11.5)') ((Lat(i,j),i=1,11),j=1,10)

            write(handle_QC,*)
            write(handle_QC,'(A)') 'first 11 Longitude for 10 lines'
            write(handle_QC,'(11f11.5)') ((Lon(i,j),i=1,11),j=1,10)

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
            write(handle_QC,'(A)') 'first 11 Sfc_Temps for 10 lines'
            write(handle_QC,'(11F7.2)') ((sfc_temp(i,j),i=1,11),j=1,10)
            END IF

C--- Establish interface between variables used for reading in
C MODIS data (these variables are passed as subroutine arguments)
C and variables used locally in near-IR water vapor science
C algorithm (these variables are passed via COMMON statements).
C Interface established by Bo-Cai Gao

            DO J = 1, No_Lines_Per_Scan
            DO I = 1, No_Frames(Scan_No)

               RAD_0P86(I,J)     = Refl_2(i,j)
               RAD_0P905(I,J)    = Refl_17(i,j)
               RAD_0P935(I,J)    = Refl_18(i,j)
               RAD_0P94(I,J)     = Refl_19(i,j)
               RAD_1P24(I,J)     = Refl_5(i,j)

               ZENITH_SOL(I,J)   = SolZen(i,j)
               ZENITH_VIEW(I,J)  = SatZen(i,j)
               AZIMUTH_REL(I,J)  = RelAZ(i,j)

               I_LAND_WATER(I,J) = LandSea_Flag(i,j)
               I_CLOUD(I,J)      = CldMsk(i,j)
               I_GLINT(I,J)      = ibits( CldMsk_QA(I,J), 4, 1)

               SURF_ELEV(I,J)    = Height(i,j)
               SURF_T(I,J)       = sfc_temp(i,j)

            END DO
            END DO

            CALL PROCESS_MODIS_DATA(Data_Size)

            DO J = 1, No_Lines_Per_Scan
            DO I = 1, No_Frames(Scan_No)
                QA_Temp = 0

               IF(VAPVRT(I,J).GT.VAPVRT_MIN) THEN
                I_ret=I_ret+1
                                             SumQual_Flag = 1
                  
                IF(REFL_0P86(I,J).LT.0.1)       Qual_Flag = 1
                IF( (REFL_0P86(I,J).GE.0.1).AND.
     @              (REFL_0P86(I,J).LE.0.15) )  Qual_Flag = 2
                IF(REFL_0P86(I,J).GT.0.15)      Qual_Flag = 3


                IF(I_Land_Water(I,J).EQ.0)      Suf_Flag  = 1
                IF(I_GLINT(I,J).EQ.0)           Suf_Flag  = 3               
                IF(I_Land_Water(I,J).GE.1)      Suf_Flag  = 0
                IF(I_CLOUD(I,J).EQ.0)           Suf_Flag  = 2

                                               Ret_Method = 1
                IF(Suf_Flag.GE.1)              Ret_Method = 0


                CALL BYTE_SET(SumQual_Flag,0,QA_Temp)
                CALL BYTE_SET(Qual_Flag,1,QA_Temp)
                CALL BYTE_SET(Ret_Method,4,QA_Temp)
                CALL BYTE_SET(Suf_Flag,6,QA_Temp)

               ELSE
                                             SumQual_Flag = 0
                                                Qual_Flag = 0
                                                Suf_Flag  = 1
                                               Ret_Method = 2
                                         VAPVRT(I,J) = -9.999
                CALL BYTE_SET(SumQual_Flag,0,QA_Temp)
                CALL BYTE_SET(Qual_Flag,1,QA_Temp)
                CALL BYTE_SET(Ret_Method,4,QA_Temp)
                CALL BYTE_SET(Suf_Flag,6,QA_Temp)

               ENDIF

C---Modified by Bo-Cai Gao on 4/28/2004 to avoid problems with missing lines
C       and subsequent un-realistic cloud mask values, LandSea Flags, & SunGlint
C       flags in MOD35 idata and also to overpass cloud mask problems caused by 
C       the cloud mask reader (reflected in the re-assignment of values 
C       for the variable CldMsk_QA(p1,l1) in the subroutine "CldMsk_Info", which
C       is within the subroutine "GetModisData_MOD05_V2.f"
  
           IF(ibits(CldMsk_QA(I,J), 0,1) .EQ. 0) QA_Temp = 0

C---End of modification

C---Modified by Bo-Cai Gao on 4/29/2004 to avoid inconsistency in previous 
C       fill value and bad pixel flag assignment
           IF(Vflag_2(I,J)  .LE. -1)          QA_Temp = 0
           IF(Vflag_2(I,J)  .LE. -1)          VAPVRT(I,J) = -9.999

           IF(Vflag_5(I,J)  .LE. -1)          QA_Temp = 0
           IF(Vflag_5(I,J)  .LE. -1)          VAPVRT(I,J) = -9.999

           IF(Vflag_17(I,J) .LE. -1)          QA_Temp = 0
           IF(Vflag_17(I,J) .LE. -1)          VAPVRT(I,J) = -9.999

           IF(Vflag_18(I,J) .LE. -1)          QA_Temp = 0
           IF(Vflag_18(I,J) .LE. -1)          VAPVRT(I,J) = -9.999

           IF(Vflag_19(I,J) .LE. -1)          QA_Temp = 0
           IF(Vflag_19(I,J) .LE. -1)          VAPVRT(I,J) = -9.999

C---Note: The following 2 lines get rid of problems for pixels with
C---      solar zenith angles greater than "SolZen_Threshold", where
C---      all Vflags are assigned as zero (not "-1"), Refl_2, Refl_5,
C---      Refl_17, Refl_18 & Refl_19 were assigned to FV_L1, and incorrect
C---      near-IR water vapor retrievals were performed.
C--- 
           IF(REFL_0P86(I,J).LE.0.001)        QA_Temp = 0
           IF(REFL_0P86(I,J).LE.0.001)        VAPVRT(I,J) = -9.999

C---End of modification

C---Temp code
C---                         REFL_0p86(i,j)  = float(SumQual_Flag)
C---                         REFL_0p86(i,j)  = float( Vflag_2(i,j) )

               QA_NIR(I,J)=QA_Temp

            END DO
C---             write(61) (REFL_0P86(I,J),  i = 1, 1354)
            END DO


            I_total = I_total + No_Lines_Per_Scan*No_Frames(Scan_No)

            CALL MOD05_PutScanData(Scan_No,modfil_MOD05,
     1                 No_QA_Bytes,Data_Size,
     2                 Buf_Size_QA,Buf_Size1,Buf_Size2,
     3                 VAPVRT,CldMsk_QA,QA_NIR)

         END IF        !End check on Day/Night flag

 999     CONTINUE      !End loop on scans


*/  Modified by JC Guu  11/05/96
*/  Call MOD05_ecsMeta_data to handle the ECS metadata.

       DO I=1,No_PSA
         QA_Metadata_MOD05(I)=0.0
       ENDDO

C
C Set Successful Retrieval Rate at near infrared  to
C QA_Metadata_MOD05 array (1)
C

       QA_Metadata_MOD05(1)=100.0*float(I_ret)/float(I_total)

       CALL MOD05_ecsMeta_data(QA_Metadata_MOD05,RTN_NCEP,groups,
     * HDFAttNms,NumHandles)

C
C      CALL ecsMeta_data(Modfil_L1B,groups,HDFAttNms,NumHandles,
C     &                  Total_No_Lines,Data_Size(1))
C

C
C Close HDF file
C

      ioflag = 'CLOSE'

*/  Modified by JC Guu  10/11/96
*/  The variables groups, HDFAttNms, and NumHandles
*/  are added to the argument list.

      CALL FILEOC(ioflag,modis_flag,handle_trans,
     &            handle_weight,handle_QC,FN_L1B,
     &            modfil_L1B,modfil_Geo,modfil_CldMsk,
     &            modfil_Anc,modfil_MOD05,groups,
     &            HDFAttNms,NumHandles)

C \* Endif for MinSolarZenithAngle

      ENDIF

C \* Endif for modis_flag

      ENDIF

      RETURN
      END


C----------------------------------------------------------------------
      SUBROUTINE TRANSM_H2O_RATIOS(handle_trans)
C----------------------------------------------------------------------
C!F77
C
C!Description:This subroutine is to read in six water vapor
C              transmittance tables.
C
C!Input Parameters:six water vapor transmittance tables.
C
C!Output Parameters:transmittance ratios obtained by ratioing
C                    transmittances of absorption channels against
C                    those from nearby "window" chanels.
C
C!Revision History:
c 01/28/98 fhliang
c fixed prolog.
c
C   Added SDP Toolkit 4 file I/O calls, and explicit data types.
C
c!TEAM-UNIQUE HEADER:
c
c Developed by MODIS Aerosol/Solar Water Vapor Retrieval Team
c GSFC, Greenbelt, MD
c
C
C!REFERENCES AND CREDITS
C    Written by Dr. Bo-Cai Gao    10/28/94
C    gao@imagecube.gsfc.nasa.gov
C
C!DESIGN NOTES: data and variables are passed with common statements
C               contained in "COMMONS.inc"
C
C!END
C----------------------------------------------------------------------

*/  Modified by JC Guu  10/18/96
*/  The statement "IMPLICIT REAL(A-H, O-Z)" is
*/  replaced by "IMPLICIT NONE".

      IMPLICIT NONE

*/  Modified by JC Guu  08/23/96
*/  The file name for COMMONS.inc was COMMONS_INC.f.

      INCLUDE 'COMMONS.inc'
      SAVE

      INTEGER K,handle_trans(6)
C
C---Calculating weighting factors for 3-channel ratios

      WT_0P905_L = (WV_CNTR_1P24  - WV_CNTR_0P905) /
     &             (WV_CNTR_1P24  - WV_CNTR_0P86)

      WT_0P905_R = (WV_CNTR_0P905 - WV_CNTR_0P86)  /
     &             (WV_CNTR_1P24  - WV_CNTR_0P86)
C
C-      PRINT*, 'WT_0P905_L = ',WT_0P905_L, 'WT_0P905_R = ',WT_0P905_R
C
C
C
      WT_0P935_L = (WV_CNTR_1P24  - WV_CNTR_0P935) /
     &             (WV_CNTR_1P24  - WV_CNTR_0P86)

      WT_0P935_R = (WV_CNTR_0P935 - WV_CNTR_0P86)  /
     &             (WV_CNTR_1P24  - WV_CNTR_0P86)
C
C-      PRINT*, 'WT_0P935_L = ',WT_0P935_L, 'WT_0P935_R = ',WT_0P935_R
C
C
C
      WT_0P94_L  = (WV_CNTR_1P24  - WV_CNTR_0P94)  /
     &             (WV_CNTR_1P24  - WV_CNTR_0P86)

      WT_0P94_R  = (WV_CNTR_0P94  - WV_CNTR_0P86)  /
     &             (WV_CNTR_1P24  - WV_CNTR_0P86)
C
C-      PRINT*, 'WT_0P94_L  = ', WT_0P94_L, 'WT_0P94_R  = ', WT_0P94_R

C
C
C---Read in pre-calculated transmittance tables for six LOWTRAN7 standard
C        model atmospheres (Tropical, Mid-summer, Mid_winter,
C        Sub-arctic summer, Sub-arctic winter, US STD 1976).
C

      READ(handle_trans(1),*)
      DO K = 1, N_GH2O
         READ(handle_trans(1),*) VAPTOT_1(K),T_0P86_1(K),T_0P905_1(K),
     &                   T_0P935_1(K), T_0P94_1(K), T_1P24_1(K)
      END DO
C
C---Reset PCF_INDEX to second atmospheric transmission model
C
      READ(handle_trans(2),*)
      DO K = 1, N_GH2O
         READ(handle_trans(2),*) VAPTOT_2(K),T_0P86_2(K),T_0P905_2(K),
     &                 T_0P935_2(K), T_0P94_2(K), T_1P24_2(K)
      END DO
C
C---Reset PCF_INDEX to third atmospheric transmission model
C
      READ(handle_trans(3),*)
      DO K = 1, N_GH2O
         READ(handle_trans(3),*) VAPTOT_3(K),T_0P86_3(K),T_0P905_3(K),
     &                 T_0P935_3(K), T_0P94_3(K), T_1P24_3(K)
      END DO
C
C---Reset PCF_INDEX to fourth atmospheric transmission model
C

      READ(handle_trans(4),*)
      DO K = 1, N_GH2O
         READ(handle_trans(4),*) VAPTOT_4(K),T_0P86_4(K),T_0P905_4(K),
     &                 T_0P935_4(K), T_0P94_4(K), T_1P24_4(K)
      END DO
C
C---Reset PCF_INDEX to fifth atmospheric transmission model
C
      READ(handle_trans(5),*)
      DO K = 1, N_GH2O
         READ(handle_trans(5),*) VAPTOT_5(K),T_0P86_5(K),T_0P905_5(K),
     &                 T_0P935_5(K), T_0P94_5(K), T_1P24_5(K)
      END DO
C
C---Reset PCF_INDEX to sixth atmospheric transmission model
C
      READ(handle_trans(6),*)
      DO K = 1, N_GH2O
         READ(handle_trans(6),*) VAPTOT_6(K),T_0P86_6(K),T_0P905_6(K),
     &              T_0P935_6(K), T_0P94_6(K), T_1P24_6(K)
      END DO
C
C
C---Calculate 3-channel ratios (e.g., T(0.94) / {WT_0P94_L * T(0.86) +
C      WT_0P94_R * T(1.24)} ) and 2-channel ratios (e.g., T(0.94)/T(0.86) )
C      using transmittance values in the tables.
C
C
      DO K = 1, N_GH2O
         R_3CH_0P905_MDL1(K) = T_0P905_1(K) / (WT_0P905_L * T_0P86_1(K)
     &                                     +  WT_0P905_R * T_1P24_1(K))
         R_2CH_0P905_MDL1(K) = T_0P905_1(K) /  T_0P86_1(K)


         R_3CH_0P935_MDL1(K) = T_0P935_1(K) / (WT_0P935_L * T_0P86_1(K)
     &                                     +  WT_0P935_R * T_1P24_1(K))
         R_2CH_0P935_MDL1(K) = T_0P935_1(K) /  T_0P86_1(K)


         R_3CH_0P94_MDL1(K)  = T_0P94_1(K)  / (WT_0P94_L  * T_0P86_1(K)
     &                                     +  WT_0P94_R  * T_1P24_1(K))
         R_2CH_0P94_MDL1(K)  = T_0P94_1(K)  /  T_0P86_1(K)
      END DO

C
      DO K = 1, N_GH2O
         R_3CH_0P905_MDL2(K) = T_0P905_2(K) / (WT_0P905_L * T_0P86_2(K)
     &                                     +  WT_0P905_R * T_1P24_2(K))
         R_2CH_0P905_MDL2(K) = T_0P905_2(K) /  T_0P86_2(K)


         R_3CH_0P935_MDL2(K) = T_0P935_2(K) / (WT_0P935_L * T_0P86_2(K)
     &                                     +  WT_0P935_R * T_1P24_2(K))
         R_2CH_0P935_MDL2(K) = T_0P935_2(K) /  T_0P86_2(K)


         R_3CH_0P94_MDL2(K)  = T_0P94_2(K)  / (WT_0P94_L  * T_0P86_2(K)
     &                                     +  WT_0P94_R  * T_1P24_2(K))
         R_2CH_0P94_MDL2(K)  = T_0P94_2(K)  /  T_0P86_2(K)
      END DO


C
      DO K = 1, N_GH2O
         R_3CH_0P905_MDL3(K) = T_0P905_3(K) / (WT_0P905_L * T_0P86_3(K)
     &                                     +  WT_0P905_R * T_1P24_3(K))
         R_2CH_0P905_MDL3(K) = T_0P905_3(K) /  T_0P86_3(K)


         R_3CH_0P935_MDL3(K) = T_0P935_3(K) / (WT_0P935_L * T_0P86_3(K)
     &                                     +  WT_0P935_R * T_1P24_3(K))
         R_2CH_0P935_MDL3(K) = T_0P935_3(K) /  T_0P86_3(K)


         R_3CH_0P94_MDL3(K)  = T_0P94_3(K)  / (WT_0P94_L  * T_0P86_3(K)
     &                                     +  WT_0P94_R  * T_1P24_3(K))
         R_2CH_0P94_MDL3(K)  = T_0P94_3(K)  /  T_0P86_3(K)
      END DO

C
      DO K = 1, N_GH2O
         R_3CH_0P905_MDL4(K) = T_0P905_4(K) / (WT_0P905_L * T_0P86_4(K)
     &                                     +  WT_0P905_R * T_1P24_4(K))
         R_2CH_0P905_MDL4(K) = T_0P905_4(K) /  T_0P86_4(K)


         R_3CH_0P935_MDL4(K) = T_0P935_4(K) / (WT_0P935_L * T_0P86_4(K)
     &                                     +  WT_0P935_R * T_1P24_4(K))
         R_2CH_0P935_MDL4(K) = T_0P935_4(K) /  T_0P86_4(K)


         R_3CH_0P94_MDL4(K)  = T_0P94_4(K)  / (WT_0P94_L  * T_0P86_4(K)
     &                                     +  WT_0P94_R  * T_1P24_4(K))
         R_2CH_0P94_MDL4(K)  = T_0P94_4(K)  /  T_0P86_4(K)
      END DO

C
      DO K = 1, N_GH2O
         R_3CH_0P905_MDL5(K) = T_0P905_5(K) / (WT_0P905_L * T_0P86_5(K)
     &                                     +  WT_0P905_R * T_1P24_5(K))
         R_2CH_0P905_MDL5(K) = T_0P905_5(K) /  T_0P86_5(K)


         R_3CH_0P935_MDL5(K) = T_0P935_5(K) / (WT_0P935_L * T_0P86_5(K)
     &                                     +  WT_0P935_R * T_1P24_5(K))
         R_2CH_0P935_MDL5(K) = T_0P935_5(K) /  T_0P86_5(K)


         R_3CH_0P94_MDL5(K)  = T_0P94_5(K)  / (WT_0P94_L  * T_0P86_5(K)
     &                                     +  WT_0P94_R  * T_1P24_5(K))
         R_2CH_0P94_MDL5(K)  = T_0P94_5(K)  /  T_0P86_5(K)
      END DO

C
      DO K = 1, N_GH2O
         R_3CH_0P905_MDL6(K) = T_0P905_6(K) / (WT_0P905_L * T_0P86_6(K)
     &                                     +  WT_0P905_R * T_1P24_6(K))
         R_2CH_0P905_MDL6(K) = T_0P905_6(K) /  T_0P86_6(K)


         R_3CH_0P935_MDL6(K) = T_0P935_6(K) / (WT_0P935_L * T_0P86_6(K)
     &                                     +  WT_0P935_R * T_1P24_6(K))
         R_2CH_0P935_MDL6(K) = T_0P935_6(K) /  T_0P86_6(K)


         R_3CH_0P94_MDL6(K)  = T_0P94_6(K)  / (WT_0P94_L  * T_0P86_6(K)
     &                                     +  WT_0P94_R  * T_1P24_6(K))
         R_2CH_0P94_MDL6(K)  = T_0P94_6(K)  /  T_0P86_6(K)
      END DO



      RETURN
      END


C----------------------------------------------------------------------
      SUBROUTINE SOLAR_IRRADIANCE
C----------------------------------------------------------------------
C!F77
C
C!Description:This subroutine is to assign solar irradiance values for
C               five MODIS channels centered near 0.86, 0.905, 0.935m 0.94,
C               and 1.24 micron.
C
C!Input Parameters:solar irradiances
C
C!Output Parameters:solar irradiances
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
C!REFERENCES AND CREDITS
C    First version created by Dr. Bo-Cai Gao
C    (gao@imagecube.gsfc.nasa.gov)
C                             1994 10 28
C
C!DESIGN NOTES: data and variables are passed with common statements
C               contained in "COMMONS.inc"
C
C!END
C----------------------------------------------------------------------

*/  Modified by JC Guu  10/18/96
*/  IMPLICIT NONE is added.

      IMPLICIT NONE

*/  Modified by JC Guu  08/23/96
*/  The file name for COMMONS.inc was COMMONS_INC.f.

      INCLUDE 'COMMONS.inc'

      SAVE
C
C *** Solar IRR from Neckel & Labs (1984) in unit of watt/(m^2 um)
C     (weighted averaging using MODIS filter functions):
C
C*** Temporarily assign solar irradiance values to be 1 for the 5 MODIS
C     near-IR channels used in this water vapor algorithm.
C
C---      SOL_IRR_0P86  = 987.567
C---      SOL_IRR_0P905 = 916.687
C---      SOL_IRR_0P935 = 856.765
C---      SOL_IRR_0P94  = 857.691
C---      SOL_IRR_1P24  = 470.600

      SOL_IRR_0P86  = 1.0
      SOL_IRR_0P905 = 1.0
      SOL_IRR_0P935 = 1.0
      SOL_IRR_0P94  = 1.0
      SOL_IRR_1P24  = 1.0
C
      RETURN
      END
C
C----------------------------------------------------------------------
C
      SUBROUTINE WEIGHTING_TABLE(handle_weight)

C----------------------------------------------------------------------
C!F77
C
C!Description:This subroutine is to read in a table for performing
C              weighted average of column water vapor values derived
C              from the 0.905, 0.935, and 0.94 micron MODIS channels.
C
C!Input Parameters:a weighting table.
C
C!Output Parameters:arrays containing weighting factors.
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
C!REFERENCES AND CREDITS
C    First version created by Dr. Bo-Cai Gao
C    (gao@imagecube.gsfc.nasa.gov)
C                             1994 10 28
C
C!DESIGN NOTES: data and variables are passed with common statements
C               contained in "COMMONS.inc"
C
C!END
C----------------------------------------------------------------------

      IMPLICIT NONE

*/  Modified by JC Guu  08/23/96
*/  The file name for COMMONS.inc was COMMONS_INC.f.

      INCLUDE 'COMMONS.inc'
      SAVE
      INTEGER handle_weight,K

C *** The weighting factors are obtained in the following way: 1st the
C     transmittance table 'TRANSM_H2O.MDL_1' is read in; 2-channel
C     ratios, T(0.905)/T(0.86), T(0.935)/T(0.86), and T(0.94)/T(0.86),
C     are calculated, and plots of 2-channel ratios as a function of
C     total amount of water vapor are made; the absolute values of the
C     slopes of 2-channel ratios are computed; for a given water vapor
C     value, the weighting factor for one channel, e.g. Ch_0P905, is
C     proportional to the slope of the channel (S_0P905), and calculated
C     as:
C          WT_VAP_0P905 = S_0P905 / (S_0P905 + S_0P935 + S_0P94)

      READ(handle_weight,*)

      DO K = 1, N_GH2O
         READ(handle_weight,*) WT_VAP_0P905(K), WT_VAP_0P935(K),
     &       WT_VAP_0P94(K)
      END DO

      RETURN
      END

C-----------------------------------------------------------------------
C
      SUBROUTINE PROCESS_MODIS_DATA(Data_Size)
C
C-----------------------------------------------------------------------
C!F77
C
C!Description:This subroutine derives vertical column water vapor
C              values on a pixel-by-pixel basis using a look up table
C              procedure.  Data from one MODIS swath with 10 lines are
C              processed.
C
C!Input Parameters:channel ratios calculated from MODIS data and from
C                   theoretical calculations, cloud mask, land and water
C                   masks, solar and view zenith angles.
C
C!Output Parameters:column water vapor image.
C
C!Revision History:
c 01/28/98 fhliang
c fixed prolog.
c
C   Added SDP Toolkit 4 File I/O calls, explicit data types,
C   and repositioned DATA statement.
C
c!TEAM-UNIQUE HEADER:
c
c Developed by MODIS Aerosol/Solar Water Vapor Retrieval Team
c GSFC, Greenbelt, MD
c
C!REFERENCES AND CREDITS
C    First version created by Dr. Bo-Cai Gao
C    (gao@imagecube.gsfc.nasa.gov)
C                             1994 10 28
C
C!DESIGN NOTES: data and variables are passed with common statements
C               contained in "COMMONS.inc"
C
C   Internal subroutines called:
C               VAPOR_DRV_2CH
C               VAPOR_DRV_3CH
C
C!END
C-----------------------------------------------------------------------

*/  Modified by JC Guu  10/18/96
*/  IMPLICIT NONE is added.

      IMPLICIT NONE

*/  Modified by JC Guu  08/23/96
*/  The file name for COMMONS.inc was COMMONS_INC.f.

      INCLUDE 'COMMONS.inc'

      SAVE
C
      INTEGER I, J, K, Data_Size(2)
      REAL    PI, VAP, DEG_TO_RAD
      DATA PI, DEG_TO_RAD /3.1415926, 0.0174533/
C
C Convert MODIS radiances to apparent reflectances (( pi * Rad /
C         (cos(ZENITH_SOL) * SOL_IRR) )
C Temp code: because MODIS data in reflectance units (not in radiance units)
C         are currently read in from a standard MODIS data cube, the
C         following lines are temporarily commented out:
C-      DO J = 1, NL
C-         DO I = 1, NS
C-
C-            REFL_0P86(I,J)  = PI * RAD_0P86(I,J) /
C-     &           (COS(ZENITH_SOL(I,J)*DEG_TO_RAD) * SOL_IRR_0P86)
C-            REFL_0P905(I,J) = PI * RAD_0P905(I,J) /
C-     &           (COS(ZENITH_SOL(I,J)*DEG_TO_RAD) * SOL_IRR_0P905)
C-            REFL_0P935(I,J) = PI * RAD_0P935(I,J) /
C-     &           (COS(ZENITH_SOL(I,J)*DEG_TO_RAD) * SOL_IRR_0P935)
C-            REFL_0P94(I,J)  = PI * RAD_0P94(I,J) /
C-     &           (COS(ZENITH_SOL(I,J)*DEG_TO_RAD) * SOL_IRR_0P94)
C-            REFL_1P24(I,J)  = PI * RAD_1P24(I,J) /
C-     &           (COS(ZENITH_SOL(I,J)*DEG_TO_RAD) * SOL_IRR_1P24)
C-
C-         END DO
C-      END DO

      DO J = 1, Data_Size(2)
         DO I = 1, Data_Size(1)
            REFL_0P86(I,J)  =  RAD_0P86(I,J)
            REFL_0P905(I,J) =  RAD_0P905(I,J)
            REFL_0P935(I,J) =  RAD_0P935(I,J)
            REFL_0P94(I,J)  =  RAD_0P94(I,J)
            REFL_1P24(I,J)  =  RAD_1P24(I,J)
         END DO
      END DO

C
C Calculate 3-channel or 2-channel ratios from MODIS apparent reflectances and
C      derive column water vapor values.
C
C*!!!For pixels over water surfaces or over clouds, water vapor values should
C      be derived using 2-channel ratio methods, such as 0.905um/0.86um,
C      0.935um/0.86um, and 0.94um/0.86um. In the future, it is not necessary
C      to switch to 3-channel ratio mehtods or other 2-channel ratio methods.
C
C    For pixels over land surfaces, we are presently using 3-channel ratio
C      methods. However, in the future improvements to the present algorithm
C      for applications to land surfaces, other channel ratios, such as
C      0.935um/0.905um & 0.935um/0.94um, should also be considered.
C
C

      DO J = 1, Data_Size(2)
        DO K = 1, Data_Size(1)

C-- Note: I_cloud = 1 means clear in the current CldMsk algorithm
C--       I_cloud = 0 means cloudy

       IF ((I_LAND_WATER(K,J).GE.1).AND.(I_CLOUD(K,J).EQ.1)) THEN

         R_3CH_0P905_MODIS(K,J) = REFL_0P905(K,J) /
     &    (WT_0P905_L * REFL_0P86(K,J) + WT_0P905_R * REFL_1P24(K,J))

         R_3CH_0P935_MODIS(K,J) = REFL_0P935(K,J) /
     &    (WT_0P935_L * REFL_0P86(K,J) + WT_0P935_R * REFL_1P24(K,J))

         R_3CH_0P94_MODIS(K,J)  = REFL_0P94(K,J) /
     &    (WT_0P94_L  * REFL_0P86(K,J) + WT_0P94_R  * REFL_1P24(K,J))

         CALL VAPOR_DRV_3CH(K, J, VAP)
C
         VAPVRT(K, J) = VAP
C
       ELSE

         R_2CH_0P905_MODIS(K,J) = REFL_0P905(K,J) /  REFL_0P86(K,J)

         R_2CH_0P935_MODIS(K,J) = REFL_0P935(K,J) /  REFL_0P86(K,J)

         R_2CH_0P94_MODIS(K,J)  = REFL_0P94(K,J)  /  REFL_0P86(K,J)

         CALL VAPOR_DRV_2CH(K, J, VAP)
C
         VAPVRT(K, J) = VAP
C
       END IF
C-------------------
C Temp code:
C***      WRITE(IHANDLE,300) K, J, VAPVRT(K,J)
C***  300 FORMAT(2x,'K = ', I4, 2X, 'J = ', I4, 2X, 'VAPVRT(K) = ',f10.5)
C-------------------
        END DO
      END DO

      RETURN
      END


C----------------------------------------------------------------------
      SUBROUTINE VAPOR_DRV_3CH(K, J, VAP)
C-----------------------------------------------------------------------
C!F77
C
C!DESCRIPTION: To derive column water vapor value for one pixel using
C              a 3-channel ratio technique and look-up table procedure
C
C!INPUT PARAMETERS: channel ratios calculated from MODIS data and from
C                   theoretical calculations, solar and view zenith angles.
C
C
C!OUTPUT PARAMETERS: column water vapor value for one pixel
C
C!REVISION HISTORY:
c 01/28/98 fhliang
c fixed prolog.
c
C   Explicit data types, and repositioned
C   DATA statement.
C
c!TEAM-UNIQUE HEADER:
c
c Developed by MODIS Aerosol/Solar Water Vapor Retrieval Team
c GSFC, Greenbelt, MD
c
C REFERENCES AND CREDITS
C    First version created by Dr. Bo-Cai Gao  (gao@imagecube.gsfc.nasa.gov)
C                             1994 10 28
C
C DESIGN NOTES: data and variables are passed with common statements
C               contained in "COMMONS.inc"
C
C  Internal subroutines called:
C               HUNT
C               FIND_VAP
C
c!END-------------------------------------------------------------------

      INCLUDE 'COMMONS.inc'

      INTEGER K, J, N
      INTEGER J0P905, J0P935, J0P94
      REAL  DEG_TO_RAD, PI
      REAL  VAP, VAP_0P905, VAP_0P935, VAP_0P94

      SAVE
C
      DATA PI, DEG_TO_RAD /3.1415926, 0.0174533/
C
C Note: Based on solar zenith angles to select which model atmosphere to use.
C       Need to be changed in future improved algorithm.
C
C
C Use Model # 1 - Tropical Model
      J0P905=1
      J0P935=1
      J0P94=1
      IF(ZENITH_SOL(K,J).LE.20.) THEN
       CALL HUNT(R_3CH_0P905_MDL1, N_GH2O, R_3CH_0P905_MODIS(K,J),
     &           J0P905)
       CALL FIND_VAP(R_3CH_0P905_MDL1, N_GH2O, R_3CH_0P905_MODIS(K,J),
     &           J0P905, VAPTOT_1, VAP_0P905)

       CALL HUNT(R_3CH_0P935_MDL1, N_GH2O, R_3CH_0P935_MODIS(K,J),
     &           J0P935)
       CALL FIND_VAP(R_3CH_0P935_MDL1, N_GH2O, R_3CH_0P935_MODIS(K,J),
     &           J0P935, VAPTOT_1, VAP_0P935)

       CALL HUNT(R_3CH_0P94_MDL1,  N_GH2O, R_3CH_0P94_MODIS(K,J),
     &           J0P94)
       CALL FIND_VAP(R_3CH_0P94_MDL1,  N_GH2O, R_3CH_0P94_MODIS(K,J),
     &           J0P94,  VAPTOT_1, VAP_0P94)
      END IF
C
C Use Model # 2 - Mid-latitude Summer Model
      IF((ZENITH_SOL(K,J).GT.20.).AND.(ZENITH_SOL(K,J).LE.45.)) THEN
       CALL HUNT(R_3CH_0P905_MDL2, N_GH2O, R_3CH_0P905_MODIS(K,J),
     &           J0P905)
       CALL FIND_VAP(R_3CH_0P905_MDL2, N_GH2O, R_3CH_0P905_MODIS(K,J),
     &           J0P905, VAPTOT_2, VAP_0P905)

       CALL HUNT(R_3CH_0P935_MDL2, N_GH2O, R_3CH_0P935_MODIS(K,J),
     &             J0P935)
       CALL FIND_VAP(R_3CH_0P935_MDL2, N_GH2O, R_3CH_0P935_MODIS(K,J),
     &           J0P935, VAPTOT_2, VAP_0P935)

       CALL HUNT(R_3CH_0P94_MDL2,  N_GH2O, R_3CH_0P94_MODIS(K,J),
     &           J0P94)
       CALL FIND_VAP(R_3CH_0P94_MDL2,  N_GH2O, R_3CH_0P94_MODIS(K,J),
     &           J0P94,  VAPTOT_2, VAP_0P94)
      END IF
C
C Use Model # 3 - Mid-latitude Winter Model
      IF((ZENITH_SOL(K,J).GT.45.).AND.(ZENITH_SOL(K,J).LE.60.)) THEN
       CALL HUNT(R_3CH_0P905_MDL3, N_GH2O, R_3CH_0P905_MODIS(K,J),
     &           J0P905)
       CALL FIND_VAP(R_3CH_0P905_MDL3, N_GH2O, R_3CH_0P905_MODIS(K,J),
     &           J0P905, VAPTOT_3, VAP_0P905)

       CALL HUNT(R_3CH_0P935_MDL3, N_GH2O, R_3CH_0P935_MODIS(K,J),
     &           J0P935)
       CALL FIND_VAP(R_3CH_0P935_MDL3, N_GH2O, R_3CH_0P935_MODIS(K,J),
     &           J0P935, VAPTOT_3, VAP_0P935)

       CALL HUNT(R_3CH_0P94_MDL3,  N_GH2O, R_3CH_0P94_MODIS(K,J),
     &           J0P94)
       CALL FIND_VAP(R_3CH_0P94_MDL3,  N_GH2O, R_3CH_0P94_MODIS(K,J),
     &           J0P94,  VAPTOT_3, VAP_0P94)
      END IF
C
C Use Model # 4 - Sub-arctic Summer Model
      IF((ZENITH_SOL(K,J).GT.60.).AND.(ZENITH_SOL(K,J).LE.70.)) THEN
       CALL HUNT(R_3CH_0P905_MDL4, N_GH2O, R_3CH_0P905_MODIS(K,J),
     &           J0P905)
       CALL FIND_VAP(R_3CH_0P905_MDL4, N_GH2O, R_3CH_0P905_MODIS(K,J),
     &           J0P905, VAPTOT_4, VAP_0P905)

       CALL HUNT(R_3CH_0P935_MDL4, N_GH2O, R_3CH_0P935_MODIS(K,J),
     &           J0P935)
       CALL FIND_VAP(R_3CH_0P935_MDL4, N_GH2O, R_3CH_0P935_MODIS(K,J),
     &           J0P935, VAPTOT_4, VAP_0P935)

       CALL HUNT(R_3CH_0P94_MDL4,  N_GH2O, R_3CH_0P94_MODIS(K,J),
     &           J0P94)
       CALL FIND_VAP(R_3CH_0P94_MDL4,  N_GH2O, R_3CH_0P94_MODIS(K,J),
     &           J0P94,  VAPTOT_4, VAP_0P94)
      END IF
C
C Use Model # 5 - Sub-arctic Winter Model
      IF((ZENITH_SOL(K,J).GT.70.).AND.(SURF_T(K,J).LE.T_INV)) THEN
       CALL HUNT(R_3CH_0P905_MDL5, N_GH2O, R_3CH_0P905_MODIS(K,J),
     &           J0P905)
       CALL FIND_VAP(R_3CH_0P905_MDL5, N_GH2O, R_3CH_0P905_MODIS(K,J),
     &           J0P905, VAPTOT_5, VAP_0P905)

       CALL HUNT(R_3CH_0P935_MDL5, N_GH2O, R_3CH_0P935_MODIS(K,J),
     &           J0P935)
       CALL FIND_VAP(R_3CH_0P935_MDL5, N_GH2O, R_3CH_0P935_MODIS(K,J),
     &           J0P935, VAPTOT_5, VAP_0P935)

       CALL HUNT(R_3CH_0P94_MDL5,  N_GH2O, R_3CH_0P94_MODIS(K,J),
     &           J0P94)
       CALL FIND_VAP(R_3CH_0P94_MDL5,  N_GH2O, R_3CH_0P94_MODIS(K,J),
     &           J0P94,  VAPTOT_5, VAP_0P94)
      END IF
C
C Use Model # 6 - U.S. STD (1976)
      IF((ZENITH_SOL(K,J).GT.70.).AND.(SURF_T(K,J).GT.T_INV)) THEN
       CALL HUNT(R_3CH_0P905_MDL6, N_GH2O, R_3CH_0P905_MODIS(K,J),
     &           J0P905)
       CALL FIND_VAP(R_3CH_0P905_MDL6, N_GH2O, R_3CH_0P905_MODIS(K,J),
     &           J0P905, VAPTOT_6, VAP_0P905)

       CALL HUNT(R_3CH_0P935_MDL6, N_GH2O, R_3CH_0P935_MODIS(K,J),
     &           J0P935)
       CALL FIND_VAP(R_3CH_0P935_MDL6, N_GH2O, R_3CH_0P935_MODIS(K,J),
     &           J0P935, VAPTOT_6, VAP_0P935)

       CALL HUNT(R_3CH_0P94_MDL6,  N_GH2O, R_3CH_0P94_MODIS(K,J),
     &           J0P94)
       CALL FIND_VAP(R_3CH_0P94_MDL6,  N_GH2O, R_3CH_0P94_MODIS(K,J),
     &           J0P94,  VAPTOT_6, VAP_0P94)
      END IF
C
C A geometrical factor related to the solar and view zenith angles
C
      GEOM_FACTOR  =   1. / COS(ZENITH_SOL(K,J)  * DEG_TO_RAD)
     &               + 1. / COS(ZENITH_VIEW(K,J) * DEG_TO_RAD)
C
C A pressure scaling factor for water vapor amount based on the LOWTRAN
C     scaling algorithm ( w(z) P(z)/P0 ) with the scaling exponent of 0.9
C
      SCALE_FACTOR =  ( EXP( -SURF_ELEV(K,J)/SCALE_HEIGHT) )
     &                   **SCALE_EXPONENT

C\* Set SCALE_FACTOR to 1.0

C---      SCALE_FACTOR = 1.0
C
C Convert and average to get total column water vapor values.
C Note: for bad detectors in 1-km channels, the input radiances or
C     reflectance  will be negative. For the 250 m channel (0.865 micron)
C     and 500 m channel (1.24 micron), the aggregated radiances
C     at 1 km scale will likely be positive. Therefore, for pixels with
C     bad detectors, the corresponding 2-channel and 3-channel ratios
C     are mostly negative. The negative ratios will give very large column
C     water vapor values, and these values will not be included in our
C     final averaging process to get the output column water vapor
C     amounts for these pixels.
C
      VAP_0P905 = VAP_0P905 / (GEOM_FACTOR * SCALE_FACTOR )
      VAP_0P935 = VAP_0P935 / (GEOM_FACTOR * SCALE_FACTOR )
      VAP_0P94  = VAP_0P94  / (GEOM_FACTOR * SCALE_FACTOR )
C
      L_VAP_0P905 = (VAP_0P905.GT.VAPVRT_MIN).AND.
     &              (VAP_0P905.LT.VAPVRT_MAX)
      L_VAP_0P935 = (VAP_0P935.GT.VAPVRT_MIN).AND.
     &              (VAP_0P935.LT.VAPVRT_MAX)
      L_VAP_0P94  = (VAP_0P94.GT.VAPVRT_MIN).AND.
     &              (VAP_0P94.LT.VAPVRT_MAX)
C
C Several cases for averaging...
C
C   Case a:  VAP_0P905, VAP_0P935, and VAP_0P94 are all within
C            the bounds between VAPVRT_MIN and VAPVRT_MAX
C
      IF( (L_VAP_0P905 .and. L_VAP_0P935) .and. L_VAP_0P94 ) THEN
         WT_TOTAL = WT_VAP_0P905(J0P905) + WT_VAP_0P935(J0P935)
     &                                + WT_VAP_0P94(J0P94)
         VAP = (   VAP_0P905 * WT_VAP_0P905(J0P905)
     &           + VAP_0P935 * WT_VAP_0P935(J0P935)
     &           + VAP_0P94  * WT_VAP_0P94(J0P94) ) / WT_TOTAL
         GO TO 9999
      END IF

C
C   Case b:  VAP_0P905, VAP_0P935, and VAP_0P94 are all outside
C            the bounds between VAPVRT_MIN and VAPVRT_MAX
C
      IF( ((.not.L_VAP_0P905) .and. (.not.L_VAP_0P935)) .and.
     &      (.not.L_VAP_0P94) ) THEN
         VAP = VAP_FILL_FLOAT
         GO TO 9999
      END IF

C
C   Case c:  among VAP_0P905, VAP_0P935, and VAP_0P94, one is outside
C            the bounds between VAPVRT_MIN and VAPVRT_MAX, while
C            the other two are within the bounds.
C
      IF( ((.not.L_VAP_0P905) .and. L_VAP_0P935) .and. L_VAP_0P94) THEN
         WT_TOTAL =     WT_VAP_0P935(J0P935) + WT_VAP_0P94(J0P94)
         VAP      =  (  VAP_0P935 * WT_VAP_0P935(J0P935)
     &                + VAP_0P94  * WT_VAP_0P94(J0P94)    ) / WT_TOTAL
         GO TO 9999
      END IF

C
      IF( (L_VAP_0P905 .and. (.not.L_VAP_0P935)) .and. L_VAP_0P94) THEN
         WT_TOTAL =     WT_VAP_0P905(J0P905) + WT_VAP_0P94(J0P94)
         VAP      =  (  VAP_0P905 * WT_VAP_0P905(J0P905)
     &                + VAP_0P94  * WT_VAP_0P94(J0P94)    ) / WT_TOTAL
         GO TO 9999
      END IF

C
      IF( (L_VAP_0P905 .and. L_VAP_0P935) .and. (.not.L_VAP_0P94)) THEN
         WT_TOTAL = WT_VAP_0P905(J0P905) + WT_VAP_0P935(J0P935)
         VAP      =  (  VAP_0P905 * WT_VAP_0P905(J0P905)
     &                + VAP_0P935 * WT_VAP_0P935(J0P935)  ) / WT_TOTAL
         GO TO 9999
      END IF

C
C   Case d:  among VAP_0P905, VAP_0P935, and VAP_0P94, two are outside
C            the bounds between VAPVRT_MIN and VAPVRT_MAX, while
C            the other one is within the bounds.
C
      IF((L_VAP_0P905 .and. (.NOT.L_VAP_0P935)) .and. (.NOT.L_VAP_0P94))
     &   VAP      =     VAP_0P905

      IF(((.not.L_VAP_0P905) .and. L_VAP_0P935) .and. (.NOT.L_VAP_0P94))
     &   VAP      =     VAP_0P935

      IF(((.not.L_VAP_0P905) .and. (.not.L_VAP_0P935)) .and. L_VAP_0P94)
     &   VAP      =     VAP_0P94
C
C
 9999 CONTINUE
C
C  Make sure the water vapor fill value is assigned for all pixels having
C       column water vapor values outside the specified upper and lower bounds.
C
C-- Note: The mis-registrations among MODIS channels at 0.86, 0.905,
C         0.935, 0.94, and 1.24 micron (particularly over boundary
C         areas with contrasting brightnesses),  the bad, and missing data
C         in these channels will result in 3-channel ratio and
C         2-channels ratio values > 1., or < 0.  The current version of
C         vapor algorithm uses FIND_VAP, and HUNT routines for table searching.
C         For pixels with channel ratios > 1.0, the algorithm gives 0.0 cm
C         of column water vapor. For pixels with very small channel ratio
C         values or negative channel ratio values, the algorithm gives
C         column water vapor > 190. cm. Therefore, the following
C         statements: IF( VAP .GT. 30 )      VAP = -0.099
C         &           IF( VAP .LT. 0.00001 ) VAP = -0.099
C         will assign fill value ( - 0.99 * 1000 = - 99, where 1000 is
C         the scaling factor used in generating output HDF file ) to pixels
C         with channel ratios > 1, or < 0.
C
         IF( VAP .GT. VAPVRT_MAX) VAP = VAP_FILL_FLOAT
         IF( VAP .LT. VAPVRT_MIN) VAP = VAP_FILL_FLOAT
C
      RETURN
      END


C-------------------------------------------------------------------------------
      SUBROUTINE VAPOR_DRV_2CH(K, J, VAP)

C-----------------------------------------------------------------------
C!F77
C
C!DESCRIPTION: To derive column water vapor value for one pixel using
C              a 2-channel ratio technique and look-up table procedure
C
C!INPUT PARAMETERS: channel ratios calculated from MODIS data and from
C                   theoretical calculations, solar and view zenith angles.
C
C
C!OUTPUT PARAMETERS: column water vapor value for one pixel
C
C!REVISION HISTORY:
c 01/28/98 fhliang
c fixed prolog.
c
C   Explicit data types, and repositioned
C   DATA statement.
C
c!TEAM-UNIQUE HEADER:
c
c Developed by MODIS Aerosol/Solar Water Vapor Retrieval Team
c GSFC, Greenbelt, MD
c
C REFERENCES AND CREDITS
C    First version created by Dr. Bo-Cai Gao  (gao@imagecube.gsfc.nasa.gov)
C                             1994 10 28
C
C DESIGN NOTES: data and variables are passed with common statements
C               contained in "COMMONS.inc"
C
C  Internal subroutines called:
C               HUNT
C               FIND_VAP
C
c!END-------------------------------------------------------------------

      IMPLICIT NONE

      INCLUDE 'COMMONS.inc'

      INTEGER K,J
      INTEGER J0P905, J0P935, J0P94
      REAL  DEG_TO_RAD, PI
      REAL  VAP, VAP_0P905, VAP_0P935, VAP_0P94


      SAVE
C
      DATA PI, DEG_TO_RAD /3.1415926, 0.0174533/
C
C Note: Based on solar zenith angles to select which model atmosphere to use.
C       Need to be changed in future improved algorithm.
C
C
C Use Model # 1 - Tropical Model
      J0P905=1
      J0P935=1
      J0P94=1

      IF(ZENITH_SOL(K,J).LE.20.) THEN
         CALL HUNT(R_2CH_0P905_MDL1, N_GH2O, R_2CH_0P905_MODIS(K,J),
     &             J0P905)
         CALL FIND_VAP(R_2CH_0P905_MDL1, N_GH2O, R_2CH_0P905_MODIS(K,J),
     &             J0P905, VAPTOT_1, VAP_0P905)

         CALL HUNT(R_2CH_0P935_MDL1, N_GH2O, R_2CH_0P935_MODIS(K,J),
     &             J0P935)
         CALL FIND_VAP(R_2CH_0P935_MDL1, N_GH2O, R_2CH_0P935_MODIS(K,J),
     &             J0P935, VAPTOT_1, VAP_0P935)

         CALL HUNT(R_2CH_0P94_MDL1,  N_GH2O, R_2CH_0P94_MODIS(K,J),
     &             J0P94)
         CALL FIND_VAP(R_2CH_0P94_MDL1,  N_GH2O, R_2CH_0P94_MODIS(K,J),
     &             J0P94,  VAPTOT_1, VAP_0P94)
      END IF
C
C Use Model # 2 - Mid-latitude Summer Model
      IF((ZENITH_SOL(K,J).GT.20.).AND.(ZENITH_SOL(K,J).LE.45.)) THEN
         CALL HUNT(R_2CH_0P905_MDL2, N_GH2O, R_2CH_0P905_MODIS(K,J),
     &             J0P905)
         CALL FIND_VAP(R_2CH_0P905_MDL2, N_GH2O, R_2CH_0P905_MODIS(K,J),
     &             J0P905, VAPTOT_2, VAP_0P905)

         CALL HUNT(R_2CH_0P935_MDL2, N_GH2O, R_2CH_0P935_MODIS(K,J),
     &             J0P935)
         CALL FIND_VAP(R_2CH_0P935_MDL2, N_GH2O, R_2CH_0P935_MODIS(K,J),
     &             J0P935, VAPTOT_2, VAP_0P935)

         CALL HUNT(R_2CH_0P94_MDL2,  N_GH2O, R_2CH_0P94_MODIS(K,J),
     &             J0P94)
         CALL FIND_VAP(R_2CH_0P94_MDL2,  N_GH2O, R_2CH_0P94_MODIS(K,J),
     &             J0P94,  VAPTOT_2, VAP_0P94)
      END IF
C
C Use Model # 3 - Mid-latitude Winter Model
      IF((ZENITH_SOL(K,J).GT.45.).AND.(ZENITH_SOL(K,J).LE.60.)) THEN
         CALL HUNT(R_2CH_0P905_MDL3, N_GH2O, R_2CH_0P905_MODIS(K,J),
     &             J0P905)
         CALL FIND_VAP(R_2CH_0P905_MDL3, N_GH2O, R_2CH_0P905_MODIS(K,J),
     &             J0P905, VAPTOT_3, VAP_0P905)

         CALL HUNT(R_2CH_0P935_MDL3, N_GH2O, R_2CH_0P935_MODIS(K,J),
     &             J0P935)
         CALL FIND_VAP(R_2CH_0P935_MDL3, N_GH2O, R_2CH_0P935_MODIS(K,J),
     &             J0P935, VAPTOT_3, VAP_0P935)

         CALL HUNT(R_2CH_0P94_MDL3,  N_GH2O, R_2CH_0P94_MODIS(K,J),
     &             J0P94)
         CALL FIND_VAP(R_2CH_0P94_MDL3,  N_GH2O, R_2CH_0P94_MODIS(K,J),
     &             J0P94,  VAPTOT_3, VAP_0P94)
      END IF
C
C Use Model # 4 - Sub-arctic Summer Model
      IF((ZENITH_SOL(K,J).GT.60.).AND.(ZENITH_SOL(K,J).LE.70.)) THEN
         CALL HUNT(R_2CH_0P905_MDL4, N_GH2O, R_2CH_0P905_MODIS(K,J),
     &             J0P905)
         CALL FIND_VAP(R_2CH_0P905_MDL4, N_GH2O, R_2CH_0P905_MODIS(K,J),
     &             J0P905, VAPTOT_4, VAP_0P905)

         CALL HUNT(R_2CH_0P935_MDL4, N_GH2O, R_2CH_0P935_MODIS(K,J),
     &             J0P935)
         CALL FIND_VAP(R_2CH_0P935_MDL4, N_GH2O, R_2CH_0P935_MODIS(K,J),
     &             J0P935, VAPTOT_4, VAP_0P935)

         CALL HUNT(R_2CH_0P94_MDL4,  N_GH2O, R_2CH_0P94_MODIS(K,J),
     &             J0P94)
         CALL FIND_VAP(R_2CH_0P94_MDL4,  N_GH2O, R_2CH_0P94_MODIS(K,J),
     &             J0P94,  VAPTOT_4, VAP_0P94)
      END IF
C
C Use Model # 5 - Sub-arctic Winter Model
      IF((ZENITH_SOL(K,J).GT.70.).AND.(SURF_T(K,J).LE.T_INV)) THEN
         CALL HUNT(R_2CH_0P905_MDL5, N_GH2O, R_2CH_0P905_MODIS(K,J),
     &             J0P905)
         CALL FIND_VAP(R_2CH_0P905_MDL5, N_GH2O, R_2CH_0P905_MODIS(K,J),
     &             J0P905, VAPTOT_5, VAP_0P905)

         CALL HUNT(R_2CH_0P935_MDL5, N_GH2O, R_2CH_0P935_MODIS(K,J),
     &             J0P935)
         CALL FIND_VAP(R_2CH_0P935_MDL5, N_GH2O, R_2CH_0P935_MODIS(K,J),
     &             J0P935, VAPTOT_5, VAP_0P935)

         CALL HUNT(R_2CH_0P94_MDL5,  N_GH2O, R_2CH_0P94_MODIS(K,J),
     &             J0P94)
         CALL FIND_VAP(R_2CH_0P94_MDL5,  N_GH2O, R_2CH_0P94_MODIS(K,J),
     &             J0P94,  VAPTOT_5, VAP_0P94)
      END IF
C
C Use Model # 6 - U.S. STD Atmosphere (1976)
      IF((ZENITH_SOL(K,J).GT.70.).AND.(SURF_T(K,J).GT.T_INV)) THEN
         CALL HUNT(R_2CH_0P905_MDL6, N_GH2O, R_2CH_0P905_MODIS(K,J),
     &             J0P905)
         CALL FIND_VAP(R_2CH_0P905_MDL6, N_GH2O, R_2CH_0P905_MODIS(K,J),
     &             J0P905, VAPTOT_6, VAP_0P905)

         CALL HUNT(R_2CH_0P935_MDL6, N_GH2O, R_2CH_0P935_MODIS(K,J),
     &             J0P935)
         CALL FIND_VAP(R_2CH_0P935_MDL6, N_GH2O, R_2CH_0P935_MODIS(K,J),
     &             J0P935, VAPTOT_6, VAP_0P935)

         CALL HUNT(R_2CH_0P94_MDL6,  N_GH2O, R_2CH_0P94_MODIS(K,J),
     &             J0P94)
         CALL FIND_VAP(R_2CH_0P94_MDL6,  N_GH2O, R_2CH_0P94_MODIS(K,J),
     &             J0P94,  VAPTOT_6, VAP_0P94)
      END IF
C
C A geometrical factor related to the solar and view zenith angles
C
      GEOM_FACTOR  =   1. / COS(ZENITH_SOL(K,J)  * DEG_TO_RAD)
     &               + 1. / COS(ZENITH_VIEW(K,J) * DEG_TO_RAD)
C
C A pressure scaling factor for water vapor amount based on the LOWTRAN
C     scaling algorithm ( w(z) P(z)/P0 ) with the scaling exponent of 0.9
C
      SCALE_FACTOR =  ( EXP( -SURF_ELEV(K,J)/SCALE_HEIGHT) )
     &                   **SCALE_EXPONENT
C
C Convert and average to get total column water vapor values.
C Note: for bad detectors in 1-km channels, the input radiances or
C     reflectance  will be negative. For the 250 m channel (0.865 micron)
C     and 500 m channel (1.24 micron), the aggregated radiances
C     at 1 km scale will likely be positive. Therefore, for pixels with
C     bad detectors, the corresponding 2-channel and 3-channel ratios
C     are mostly negative. The negative ratios will give very large column
C     water vapor values, and these values will not be included in our
C     final averaging process to get the output column water vapor
C     amounts for these pixels.
C
      VAP_0P905 = VAP_0P905 / (GEOM_FACTOR * SCALE_FACTOR )
      VAP_0P935 = VAP_0P935 / (GEOM_FACTOR * SCALE_FACTOR )
      VAP_0P94  = VAP_0P94  / (GEOM_FACTOR * SCALE_FACTOR )
C
      L_VAP_0P905 = (VAP_0P905.GT.VAPVRT_MIN).AND.
     &              (VAP_0P905.LT.VAPVRT_MAX)
      L_VAP_0P935 = (VAP_0P935.GT.VAPVRT_MIN).AND.
     &              (VAP_0P935.LT.VAPVRT_MAX)
      L_VAP_0P94  = (VAP_0P94.GT.VAPVRT_MIN).AND.
     &              (VAP_0P94.LT.VAPVRT_MAX)
C
C Several cases for averaging...
C
C   Case a:  VAP_0P905, VAP_0P935, and VAP_0P94 are all within
C            the bounds between VAPVRT_MIN and VAPVRT_MAX
C
      IF( (L_VAP_0P905 .and. L_VAP_0P935) .and. L_VAP_0P94 ) THEN
         WT_TOTAL = WT_VAP_0P905(J0P905) + WT_VAP_0P935(J0P935)
     &                                + WT_VAP_0P94(J0P94)
         VAP = (   VAP_0P905 * WT_VAP_0P905(J0P905)
     &           + VAP_0P935 * WT_VAP_0P935(J0P935)
     &           + VAP_0P94  * WT_VAP_0P94(J0P94) ) / WT_TOTAL
         GO TO 9999
      END IF

C
C   Case b:  VAP_0P905, VAP_0P935, and VAP_0P94 are all outside
C            the bounds between VAPVRT_MIN and VAPVRT_MAX
C
      IF( ((.not.L_VAP_0P905) .and. (.not.L_VAP_0P935)) .and.
     &      (.not.L_VAP_0P94) ) THEN
         VAP = VAP_FILL_FLOAT
         GO TO 9999
      END IF

C
C   Case c:  among VAP_0P905, VAP_0P935, and VAP_0P94, one is outside
C            the bounds between VAPVRT_MIN and VAPVRT_MAX, while
C            the other two are within the bounds.
C
      IF( ((.not.L_VAP_0P905) .and. L_VAP_0P935) .and. L_VAP_0P94) THEN
         WT_TOTAL =     WT_VAP_0P935(J0P935) + WT_VAP_0P94(J0P94)
         VAP      =  (  VAP_0P935 * WT_VAP_0P935(J0P935)
     &                + VAP_0P94  * WT_VAP_0P94(J0P94)    ) / WT_TOTAL
         GO TO 9999
      END IF

C
      IF( (L_VAP_0P905 .and. (.not.L_VAP_0P935)) .and. L_VAP_0P94) THEN
         WT_TOTAL =     WT_VAP_0P905(J0P905) + WT_VAP_0P94(J0P94)
         VAP      =  (  VAP_0P905 * WT_VAP_0P905(J0P905)
     &                + VAP_0P94  * WT_VAP_0P94(J0P94)    ) / WT_TOTAL
         GO TO 9999
      END IF

C
      IF( (L_VAP_0P905 .and. L_VAP_0P935) .and. (.not.L_VAP_0P94)) THEN
         WT_TOTAL = WT_VAP_0P905(J0P905) + WT_VAP_0P935(J0P935)
         VAP      =  (  VAP_0P905 * WT_VAP_0P905(J0P905)
     &                + VAP_0P935 * WT_VAP_0P935(J0P935)  ) / WT_TOTAL
         GO TO 9999
      END IF

C
C   Case d:  among VAP_0P905, VAP_0P935, and VAP_0P94, two are outside
C            the bounds between VAPVRT_MIN and VAPVRT_MAX, while
C            the other one is within the bounds.
C
      IF((L_VAP_0P905 .and. (.NOT.L_VAP_0P935)) .and. (.NOT.L_VAP_0P94))
     &   VAP      =     VAP_0P905

      IF(((.not.L_VAP_0P905) .and. L_VAP_0P935) .and. (.NOT.L_VAP_0P94))
     &   VAP      =     VAP_0P935

      IF(((.not.L_VAP_0P905) .and. (.not.L_VAP_0P935)) .and. L_VAP_0P94)
     &   VAP      =     VAP_0P94
C
 9999 CONTINUE
C
C  Make sure the water vapor fill value is assigned for all pixels having
C       column water vapor values outside the specified upper and lower bounds.
C
C-- Note: The mis-registrations among MODIS channels at 0.86, 0.905,
C         0.935, 0.94, and 1.24 micron (particularly over boundary
C         areas with contrasting brightnesses),  the bad, and missing data
C         in these channels will result in 3-channel ratio and
C         2-channels ratio values > 1., or < 0.  The current version of
C         vapor algorithm uses FIND_VAP, and HUNT routines for table searching.
C         For pixels with channel ratios > 1.0, the algorithm gives 0.0 cm
C         of column water vapor. For pixels with very small channel ratio
C         values or negative channel ratio values, the algorithm gives
C         column water vapor > 190. cm. Therefore, the following
C         statements: IF( VAP .GT. 30 )      VAP = -0.099
C         &           IF( VAP .LT. 0.00001 ) VAP = -0.099
C         will assign fill value ( - 0.99 * 1000 = - 99, where 1000 is
C         the scaling factor used in generating output HDF file ) to pixels
C         with channel ratios > 1, or < 0.
C
         IF( VAP .GT. VAPVRT_MAX) VAP = VAP_FILL_FLOAT
         IF( VAP .LT. VAPVRT_MIN) VAP = VAP_FILL_FLOAT
C
      RETURN
      END


C***********************************************************************
      SUBROUTINE FIND_VAP(X, N, R, J, VAPTOT, VAP)

C-----------------------------------------------------------------------
C!F77
C
C!Description:A linear interpolation routine to derive a column water vapor
C          value. Specifically: giving an array X with N elements in
C          ascending or decending order, R a value satisfying the
C          relation X(J) < R < X(J+1), doing a linear interpolation using
C          two values in array VAPTOT (also with N elements) to find VAP.
C
C!Input Parameters:
C                  An array X with N elements in ascending or decending
C                  order, R a value satisfying the relation X(J) < R < X(J+1).
C
C!Output Parameters:find a value VAP based on a linear interpolation
C                    procedure.
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
C!REFERENCES AND CREDITS
C     WRITTEN BY:
C
C!DESIGN NOTES: none
C
C!END
C-----------------------------------------------------------------------
C

*/  Modified by JC Guu 10/18/96
*/  IMPLICIT NONE is added.

      IMPLICIT NONE
      INTEGER J, N
      REAL DLT, FJ, FJ_P1, R, VAP

*/  Modified by JC Guu  10/18/96
*/  A "DIMENSION" statement is deleted.

      REAL X(N), VAPTOT(N)
C
      IF( (J.GT.0).AND.(J.LT.N) ) THEN

         DLT   = X(J+1) - X(J)
         FJ    = ( X(J+1) - R ) / DLT
         FJ_P1 = ( R - X(J) )   / DLT

         VAP   = FJ * VAPTOT(J) + FJ_P1 * VAPTOT(J+1)

      ELSE

         IF(J.LE.0) VAP = VAPTOT(1)
         IF(J.GE.N) VAP = VAPTOT(N)

      END IF
C
C
C
      RETURN
      END


C*******************************************************************************

      SUBROUTINE HUNT(XX,N,X,JLO)

C-----------------------------------------------------------------------
C!F77
C
C!DESCRIPTION: finds the element in array XX that is closest to value X.
C           Array XX must be monotonic, either increasing or decreasing.
C
C!INPUT PARAMETERS: array XX with N elements and an element to
C                   search for closest match.
C
C!OUTPUT PARAMETERS: JLO - index of the closest matching element
C
C!REVISION HISTORY:
c 01/28/98 fhliang
c fixed prolog.
c
C     Explicit data type.
C
c!TEAM-UNIQUE HEADER:
c
c Developed by MODIS Aerosol/Solar Water Vapor Retrieval Team
c GSFC, Greenbelt, MD
c
C REFERENCES AND CREDITS
C     This subroutine was copied from Numerical Recipes                *
C
C DESIGN NOTES: None
C
c!END-------------------------------------------------------------------

      SAVE

      INTEGER jlo,n
      REAL x,xx(n),eps
      INTEGER inc,jhi,jm
      LOGICAL ascnd
      DATA eps/0.000001/
      ascnd=xx(n).gt.xx(1)

      if(jlo.le.0.or.jlo.gt.n)then
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      x=x+eps
      if((x-xx(jlo)).gt.eps.eqv.ascnd)then
C      if(x.ge.xx(jlo).eqv.ascnd)then
1       jhi=jlo+inc
        if(jhi.gt.n)then
          jhi=n+1
        else if((x-xx(jhi)).gt.eps.eqv.ascnd)then
C        else if(x.ge.xx(jhi).eqv.ascnd)then
          jlo=jhi
          inc=inc+inc
          goto 1
        endif
      else
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then
          jlo=0
        else if(x.lt.(xx(jlo)+eps).eqv.ascnd)then
          jhi=jlo
          inc=inc+inc
          goto 2
        endif
      endif
3     if(jhi-jlo.eq.1) return
C3     if(jhi-jlo.eq.1) then
C        if(x.eq.xx(n))jlo=n-1
C        if(x.eq.xx(1))jlo=1
C        if(abs(x-xx(n)).lt.eps)jlo=n-1
C        if(abs(x-xx(1)).lt.eps)jlo=1
C        return
C      endif
      jm=(jhi+jlo)/2
      if((x-xx(jm)).gt.eps.eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3
      END
