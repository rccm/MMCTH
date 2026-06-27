      INTEGER FUNCTION GetModisDat_MOD05_CORR(
     &        Modfil_L1B,Modfil_Geo,Modfil_CldMsk,FN_L1B,
     &        Scan_No,Buf_Size1,Buf_Size2,Data_Size,
     &        SatZen,SolZen,RelAz,Height,
     &        Refl_2,Refl_19,Ratio_Refl_19_2,
     &        Un_2,Un_19,Vflag_2,Vflag_19,
     &        Buf_Un,Buf_Sa,CldMsk,LandSea_Flag,Cloud,QA_Cloud)


      IMPLICIT NONE

      INCLUDE 'mapi.inc'
      INCLUDE 'MOD05_CORR.inc'
      INCLUDE 'PGS_MODIS_39500.f'
      INCLUDE 'L1B_Reader_V2.1.inc'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION: Retrieves one scan cube of MODIS Level 1B (L1B) sensor
C               data (bands 2 and 19), Geolocation data and
C               Cloud Mask data
C
C !INPUT PARAMETERS:
C    INTEGER    Modfil_*         File handle structure for HDF files
C    INTEGER    Modfil_L1B       File handle structure for L1B_1km files
C    CHARACTER  FN_L1B           L1B filename in 1 km resolution
C    INTEGER    Scan_No          1-based scan number (a single-sided mirror
C                                swath of earth) for the granule
C    INTEGER    Buf_Size1        Size of across swath index of L1B 1-km data
C                                the calling program
C    INTEGER    Buf_Size2        Size of along swath index of L1B 1km data
C                                the calling program
C
C !OUTPUT PARAMETERS:
C    INTEGER    Data_Size(2) Array specifying actual size of 1-km data
C                            block with within output buffer.
C                            Index 1 for columns; Index 2 for rows
C    REAL       SatZen(x,y)  Buffer to hold satellite zenith angle
C                            x=Buf_Size1, y=Buf_Size2
C    REAL       SolZen(x,y)  Buffer to hold solar zenith angle
C    REAL       Height(x,y)  Height above ellipsoid in meters
C    REAL       RelAz(x,y)   Buffer to hold relative azimuth angle
C                            between solar vector in the forward
C                            scattering direction and satellite vector.
C    REAL       Refl_*(x,y)  Buffer for band * reflectances
C    BYTE       sample*(x,y) The number of L1B (band number *) samples in
C                            aggregate data
C    REAL       Un_*(x,y)    Buffer for band * uncertainties
C    INTEGER    CldMsk(x,y)  Storage buffer to hold Cloud/No Cloud
C                            Flag: Cloud (0)/No Cloud (1)
C    INTEGER    LandSea_Flag(x,y) Storage buffer for Land/Sea flag:
C                                 Water(0),Coastal(1),Wetland(2),Land(3)
C
C !REVISION HISTORY:
C $Log: GetModisDat_MOD05_V1.f,v $
c 01/29/98 fhliang
c added NCSA acknowledgement;
c fixed prolog.
c
C Revisoin 2.0 for version 2 delivery
C New subroutine Read_L1B_V2
C
c Revision 1.6  1996/12/31  20:41:48  jguu
c A bug fix for the Read_L1B routine.
c
c Revision 1.4  1996/11/01  19:28:08  jguu
c A trivial value used to represent 0 in an equality
c test is normalized to reflect the proper magnitude.
c
c Revision 1.3  1996/10/04  16:37:04  jguu
c New delivery form ST after the implementation of the
c ECS Core Metadata.
c
c Revision 1.2  1996/08/29  19:52:43  jguu
c The file names for the include files are changed to COMMONS.inc
c and mod05.inc.
c
c Revision 1.1  1996/07/08  18:32:08  vlin
c Initial revision
c
c Revision 1.5  1996/06/26  06:29:19  rhucek
c Added functionality to identify dead detectors (in function
c GetModisDat_MOD05) and replace L1B values with -999.0.
c
c Revision 1.3  1996/04/29  19:26:02  vlin
c Read ancillary data set & added 2 arguments to output number of
c samples for band 2 and 5
c
c Revision 1.2  1996/04/26  16:50:35  vlin
c - 1. Return height from GetGeo_MOD05
c - 2. Check Fill values
c - 3. Use Aggregate_L1B function
c - 4. avoid multiple returns
c
C !TEAM-UNIQUE HEADER:
C
C    This software is developed by the MODIS Science Data Support
C    Team for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
c HDF portions developed at the National Center for Supercomputing
c Applications at the University of Illinois at Urbana-Champaign.
c
C !REFERENCES AND CREDITS
C
C    WRITTEN BY:
C
C    Richard Hucek            April 1996
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd, Seabrook MD 20706
C
C    rhucek@ltpmail.gsfc.nasa.gov
C
C !DESIGN NOTES:
C
C    GetModisDat_MOD05 checks the return status of all MODIS Application
C    Program Interface (M-API) function calls.  A successful M-API
C    function call is indicated by a return value of MAPIOK (0). If
C    unsuccessful, a warning error message (i.e. type .._W_..) is
C    written to the LogStatus file, and control reverts back to the
C    calling routine.  Subroutine MODIS_SMF_SETDYNAMICMSG is used for
C    message passing to the LogStatus file.
C
C   Externals:
C
C        MODIS_W_GENERIC            (MODIS_39500.f)
C
C   Internals:
C
C      Functions:
C        GetGeo_MOD05
C        ReadCldMsk_MOD05
C        Aggregate_L1B
C        Read_L1B
C
C      Variables:
C        Band_No     MODIS band number (1-36)
C        Scale       Spatial aggregation/subsampling factor
C        cossza      Cosine of Solar Zenith Angle
C        SolZen_Threshold  Solar Zenith Angle Processing Threshold
C
C !END
C-----------------------------------------------------------------------

      CHARACTER*(*) FN_L1B
      CHARACTER*4   Vtype
CDAC 9/4/97
      CHARACTER*1   Gain
      CHARACTER*6   Cal_Type
      INTEGER Band_No,Resol
CDAC 9/4/97
      INTEGER GetGeo_MOD05_CORR,ReadCldMsk_MOD05,
     &        Buf_Size1,Buf_Size2,Buf_Dim1,Buf_Dim2,Scan_No,rtn,
CDAC - change 3 to MODFILLEN
     &        Buf_Xtrack,Buf_Track,
     &        ipixel,iline,ip,il,Scale,Modfil_L1B(MODFILLEN),
     &        Modfil_Geo(MODFILLEN),Modfil_CldMsk(MODFILLEN),
     &        Data_Size(2),
CDAC
     &        data_size_out(2),CldMsk(Buf_Size1,Buf_Size2),
     &        LandSea_Flag(Buf_Size1,Buf_Size2),
     &        ii,ij,ik,jj,i_detector
      LOGICAL error_flag
      REAL pi,dtr,SolZen_Threshold,cossza,x_cossza,Rel_Equality_EPS,
CDAC - change 4*1500 to 1500 and 4*10 to 10
     &     Work_Buf_EV(1500,10),Work_Buf_Un(1500,10),
     &     Buf_EV(1500,10),Zero_Eps,
CDAC
     &     SatZen(Buf_Size1,Buf_Size2),SolZen(Buf_Size1,Buf_Size2),
     &     RelAz(Buf_Size1,Buf_Size2),Height(Buf_Size1,Buf_Size2),
     &     Refl_2(Buf_Size1,Buf_Size2),Refl_19(Buf_Size1,Buf_Size2),
     &     Ratio_Refl_19_2(Buf_Size1,Buf_Size2)
CDAC
      BYTE Un_2(Buf_Size1,Buf_Size2),Un_19(Buf_Size1,Buf_Size2)
      BYTE Vflag_2(Buf_Size1,Buf_Size2),Vflag_19(Buf_Size1,Buf_Size2)
      BYTE Buf_Un(Buf_Size1,Buf_Size2),Buf_Sa(Buf_Size1,Buf_Size2)
CDAC
      INTEGER DS_Dim1_CM,
     2        DS_Dim1_QA,DS_Dim2,DS_Dim3
      BYTE Cloud(Buf_cldmsk,Buf_Size1,Buf_Size2),
     2     QA_Cloud(Buf_cldmsk_QA,Buf_Size1,Buf_Size2)
C      LOGICAL Error_Flag

      PARAMETER (pi = 3.1415927, dtr = pi/180.0)
      PARAMETER (SolZen_Threshold = 80, Rel_Equality_Eps=0.000001,
     &           Zero_Eps = 0.000001)

C Set initial values
       SAVE
      GetModisDat_MOD05_CORR = -1
      error_flag = .false.

C----------------------------------------------------------------------
C Read Geolocation
C----------------------------------------------------------------------

      rtn = GetGeo_MOD05_CORR(Modfil_Geo,Scan_No,Buf_Size1,Buf_Size2,
     2      Data_Size,SolZen,SatZen,RelAz,Height)

      IF (rtn .ne. 0) THEN
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     2  'GetGeo_MOD05_CORR failed','GetModisDat_MOD05')
         error_flag = .true.
      END IF
C      WRITE(*,*) 'Pass geolocation'

C----------------------------------------------------------------------
C Read Cloud Mask data
C----------------------------------------------------------------------

      CALL Read_CldMsk(Modfil_CldMsk,Scan_No,Buf_cldmsk,
     1                 Buf_cldmsk_QA,Buf_Size1,Buf_Size2,
     2                 DS_Dim1_CM,DS_Dim1_QA,DS_Dim2,DS_Dim3,
     3                 Cloud,QA_Cloud,Error_Flag)


      IF(Error_Flag) THEN
        CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     1       'Read_CldMsk failed','GetModisDat_MOD04.f')
      ENDIF

C      WRITE(*,*) 'Pass cloud mask new'
C
C Decode cloud mask or generate our own cloud mask from L1B data
C
      CALL CldMsk_Info_CORR(Modfil_CldMsk,Buf_cldmsk,Buf_cldmsk_QA,
     1    Buf_Size1,Buf_Size2,Cloud,QA_Cloud,CldMsk,LandSea_Flag)

C
C      WRITE(*,*) 'Pass Cloud Mask info'

C--------------------------------------------------------------------
C Read L1B sensor data
CDAC
C--------------------------------------------------------------------
C
C Read MODIS band 2
C -----------------
      Band_No = 2
      Resol=Resol_Is_1km
      Cal_Type=Cal_Is_Refl
      Gain = ' '
      Buf_Xtrack=Buf_Size1
      Buf_Track=Buf_Size2

      Call Read_L1B(Modfil_L1B, FN_L1B, Scan_No, Band_No,
     &              Gain, Resol, Cal_type, Buf_Xtrack, Buf_Track,
     &              Refl_2, Un_2, Buf_Sa, Vflag_2, Data_size,
     &              error_flag)

C
C Read MODIS band 19
C ------------------
      Band_No = 19
      Resol=Resol_Is_1km
      Cal_Type=Cal_Is_Refl
      Gain = ' '
      Buf_Xtrack=Buf_Size1
      Buf_track=Buf_Size2

      Call Read_L1B(Modfil_L1B, FN_L1B, Scan_No, Band_No,
     &              Gain, Resol, Cal_type, Buf_Xtrack, Buf_Track,
     &              Refl_19, Un_19, Buf_Sa, Vflag_19, Data_size,
     &              error_flag)


      IF(Error_Flag) THEN
        CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     1       'Read_L1B failed','GetModisDat_MOD04.f')
      ENDIF

C----------------------------------------------------------------------
C Normalize sensor radiance data to reflectance units
C----------------------------------------------------------------------

      Do 40 il = 1, Data_Size(2)
      Do 45 ip = 1, Data_Size(1)
         cossza = cos(dtr*SolZen(ip,il))

         If (abs(cossza) .LT. Zero_Eps) Then
            x_cossza = 0.0
         Else
            x_cossza = 1.0/cossza
         End If

         If (abs((SolZen(ip,il)-FV_GEO)/FV_GEO) .LT.
     &        Rel_Equality_Eps) Then

            Refl_2(ip,il)  = FV_L1B
            Refl_19(ip,il) = FV_L1B
            Ratio_Refl_19_2(ip,il) = FV_L1B

         Else If (abs(cossza) .LT. Zero_Eps) Then

            Refl_2(ip,il)  = FV_L1B
            Refl_19(ip,il) = FV_L1B
            Ratio_Refl_19_2(ip,il) = FV_L1B

         Else If (SolZen(ip,il) .LT. SolZen_Threshold) Then

            If (abs((Refl_2(ip,il)-FV_L1B)/FV_L1B)
     &         .GT. Rel_Equality_Eps) then
               Refl_2(ip,il) = x_cossza*Refl_2(ip,il)
            Else
               Refl_2(ip,il) = FV_L1B
            End If

            If (abs((Refl_19(ip,il)-FV_L1B)/FV_L1B)
     &         .GT. Rel_Equality_Eps) then
               Refl_19(ip,il) = x_cossza*Refl_19(ip,il)
            Else
               Refl_19(ip,il) = FV_L1B
            End If
         Else

            Refl_2(ip,il)  = FV_L1B
            Refl_19(ip,il) = FV_L1B
            Ratio_Refl_19_2(ip,il) = FV_L1B

         End If

         If (Refl_2(ip,il) .GT. Zero_Eps) Then
           If(Refl_19(ip,il) .GT. Zero_Eps) Then
             Ratio_Refl_19_2(ip,il) = Refl_19(ip,il)/
     &                                 Refl_2(ip,il)
           Else
             Ratio_Refl_19_2(ip,il)= FV_L1B
           ENd If
         Else
           Ratio_Refl_19_2(ip,il) = FV_L1B
         End If

   45 continue
   40 continue

      if (.not.error_flag) GetModisDat_MOD05_CORR = 0

C      WRITE(*,*) 'Pass GetModisDat_MOD05_CORR'

      Return
      END


C***********************************************************************
      INTEGER FUNCTION GetGeo_MOD05_CORR(Modfil,Scan_No,Buf_Size1,
     &                 Buf_Size2,Data_Size,Sol_Zen,Sen_Zen,Rel_Az,
     &                 Height)

      IMPLICIT NONE
      SAVE

      INCLUDE 'mapi.inc'
      INCLUDE 'hdf.inc'
      INCLUDE 'MOD05_CORR.inc'
      INCLUDE 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION: Retrieves MODIS Geolocation data from an HDF target
C              array of typically 100 scan cubes (a granule).  This
C              routine is customized for MOD04 land algorithm.
C
C !INPUT PARAMETERS:
C    INTEGER Modfil(MODFILLEN)  SD_ID, File_ID and File access type.
C    INTEGER Scan_No     The scancube number to be read.
C    INTEGER Buf_Size1/2 Size of dimension 1/2 of geolocation output
C                        buffer as dimensioned in calling program
C
C !OUTPUT PARAMETERS:
C    INTEGER Data_Size(2) Size of geolocation data block within
C                         output buffer
C            In the definitions below x=Buf_Size1; y=Buf_Size2
C    REAL  Sol_Zen(x,y)   Buffer for solar zenith data
C    REAL  Sen_Zen(x,y)   Buffer for sensor zenith data.
C    REAL  Rel_Az(x,y)    Buffer for relative azimuth data
C    REAL  Height(x,y)    Buffer for height above ellipsoid data
C
C !REVISION HISTORY:
c 01/29/98 fhliang
c added NCSA acknowledgement.
c
c Revision 1.2  1995/10/18  22:41:24  rhucek
c revised GetGeo_MOD05 function to reflect data format of standard
c MOD03 product file.  The format of the synthetic MODIS data had
c been assumed prior to this update.
c
c Revision 1.1  1995/10/13  13:36:01  rhucek
c Initial revision
C
C !TEAM-UNIQUE HEADER:
C
C    This software is developed by the MODIS Science Data Support
C    Team for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
c HDF portions developed at the National Center for Supercomputing
c Applications at the University of Illinois at Urbana-Champaign.
c
C !REFERENCES AND CREDITS
C
C    WRITTEN BY:
C
C    Vicky Lin                  09/27/95
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd, Seabrook MD 20706
C
C    vlin@ltpmail.gsfc.nasa.gov
C
C !DESIGN NOTES:
C
C    GetGeo_MOD05 checks the return status of all MODIS Application
C    Program Interface (M-API) function calls.  A successful M-API
C    call is indicated by a return value of MAPIOK (0). If unsuccessful,
C    a warning error message (i.e. type .._W_..) is written to the
C    LogStatus file, and control reverts back to the calling routine.
C    Subroutine MODIS_SMF_SETDYNAMICMSG is used for message passing to
C    the LogStatus file.
C
C  Externals:
C
C    Function:
C      GMAR                       (libmapi.a)
C      GMARDM                     (libmapi.a)
C
C    Subroutines:
C      MODIS_SMF_SETDYNAMICMSG
C      CONCATENATE
C
C    Named Constant:
C      DFACC_READ                 (hdf.inc)
C      MAPIOK                     (mapi.inc)
C      MODIS_W_GENERIC            (MODIS_39500.f)
C      FV_*                       (mod06_inc.f)
C
C  Internal Variables:
C       arrnam      Name of the SDS array.
C       grpnm       Name of the data group
C       dtype       The data type of the array.
C       ret_dim(2)  Buffer holding the size of a 2-dimensional HDF
C                   SDS array.
C       Start(2)    The start position of data.
C       Fmax        Maximum frame number per scan line.
C       Lmax        Maximum line number per scan cube.
C       Rank        The number of dimensions in an array
C       MaxScan_No  Total Swath Number.
C       count       Temporary buffer for data of the target array.
C       LinesPerScan The number of lines per scan cube
C
C  Return Values:
C       GetGeo_MOD05_CORR     0 for success / -1 if failed
C
C !END
C-----------------------------------------------------------------------

      CHARACTER*13 grpnm,dtype
      CHARACTER*30 arrnm
      CHARACTER*80 msgbuf,msgbuf1,msgbuf2
      INTEGER   I,J,K,L,LinesPerScan,Rank,Fmax,Lmax,Scan_No,Buf_Size1,
     &          Buf_Size2,MaxScan_No,rtn,Rel_Equality_Eps,Start(2),
     &          ret_dim(2),Modfil(MODFILLEN),Data_Size(2)
      PARAMETER (Fmax=1500,Lmax=10,LinesPerScan=10)
      INTEGER*2 count(Fmax*Lmax)
      LOGICAL error_flag
      REAL pi
      PARAMETER (Rel_Equality_Eps = 0.000001, pi = 3.14159265)
      REAL Sol_Zen(Buf_Size1,Buf_Size2),Sen_Zen(Buf_Size1,Buf_Size2),
     &     Rel_Az(Buf_Size1,Buf_Size2),Height(Buf_Size1,Buf_Size2),
     &     Sol_Az(Fmax,Lmax),Sen_Az(Fmax,Lmax)

C Initialization

      GetGeo_MOD05_CORR=-1
      grpnm = ' '
      error_flag = .false.
      Start(1)=0
      Start(2)=(Scan_No-1)*LinesPerScan
      Data_Size(2)=LinesPerScan

      do 10 J = 1, Buf_Size2
      do 15 I = 1, Buf_Size1
         Sol_Zen(I,J) = 0.0
         Sen_Zen(I,J) = 0.0
         Sol_Az(I,J) = 0.0
         Sen_Az(I,J) = 0.0
         Rel_Az(I,J) = 0.0
         Height(I,J) = 0.0
   15 continue
   10 continue

C Checking for valid input
      IF (Modfil(1).le.0.or.Modfil(3).ne.DFACC_READ) THEN
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     2  'Invalid SD ID/file access type','GetGeo_MOD05_CORR')
         error_flag = .true.
      END IF

      do 999 I = 1, 5
         If (I.eq.1) then
            arrnm='SolarZenith'
         else If (I.eq.2) then
            arrnm='SensorZenith'
         else If (I.eq.3) then
            arrnm='SolarAzimuth'
         else if (I.eq.4) then
            arrnm='SensorAzimuth'
         else
            arrnm='Height'
         End If

C Retrieving the rank, dimensions and data type of SDS data.
         Rank = 2

         rtn = GMARDM(Modfil,arrnm,grpnm,dtype,Rank,ret_dim)
         IF (rtn.ne.MAPIOK) THEN
            CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     2     'GMARDM failed','GetGeo_MOD05_CORR')
            error_flag = .true.
         END IF

         Data_Size(1) = ret_dim(1)

C Additional input check of Scan_No.
         MaxScan_No=ret_dim(2)/LinesPerScan

         IF (Scan_No.lt.1 .or. Scan_No.gt.MaxScan_No) THEN
            write(msgbuf,'(i4)') MaxScan_No
            call Concatenate('Scan_No out of bounds : 1 - ',msgbuf,
     &      msgbuf1)
            CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,msgbuf1,
     &      'GetGeo_MOD05_CORR')
            error_flag = .true.
         END IF

C Read HDF target array into 'count' buffer
         rtn = GMAR(Modfil,arrnm,grpnm,Start,Data_Size,count)

         IF (rtn.ne. MAPIOK) THEN
            CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     &      'GMAR failed','GetGeo_MOD05_CORR')
            error_flag = .true.
         END IF

         IF (.not.error_flag) THEN

C Convert HDF 2-byte integer data to 4-byte real value

C Retrieve Solar Zenith Angle

            IF (I.eq.1) THEN
               L = 0

               do 20 K = 1, Data_Size(2)
               do 25 J = 1, Data_Size(1)
                  L = L + 1
                  Sol_Zen(J,K) = FV_GEO

                  if (count(L) .ne. FV_SolZen)
     &               Sol_Zen(J,K) = count(L)*0.01
C     &               Sol_Zen(J,K) = count(L)*0.018/pi
   25          continue
   20          continue

C Retrieve Satellite Zenith Angle

            ELSE IF (I.eq.2) THEN
               L = 0

               do 30 K = 1, Data_Size(2)
               do 35 J = 1, Data_Size(1)
                  L = L + 1
                  Sen_Zen(J,K) = FV_GEO

                  if (count(L) .ne. FV_SenZen)
     &               Sen_Zen(J,K) = count(L)*0.01
C     &               Sen_Zen(J,K) = count(L)*0.018/pi
   35          continue
   30          continue

C Retrieve Solar Azimuth

            ELSE IF (I.eq.3) THEN
               L = 0

               Do 40 K = 1, Data_Size(2)
               Do 45 J = 1, Data_Size(1)
                  L = L + 1
                  Sol_Az(J,K) = FV_GEO

                  if (count(L) .ne. FV_SolAz)
     &               Sol_Az(J,K) = count(L)*0.01
C     &               Sol_Az(J,K) = 180.0 + count(L)*0.01
C     &               Sol_Az(J,K) = 180.0 + count(L)*0.018/pi
   45          continue
   40          continue

C Retrieve Satellite Azimuth and compute Relative Azimuth

            ELSE IF (I.eq.4) THEN
               L = 0

               Do 50 K = 1, Data_Size(2)
               Do 55 J = 1, Data_Size(1)
                  L = L + 1
                  Sen_Az(J,K) = FV_GEO
                  Rel_Az(J,K) = FV_GEO

                  if (count(L) .ne. FV_SenAz)
     &               Sen_Az(J,K) = count(L)*0.01
C     &               Sen_Az(J,K) = count(L)*0.018/pi

                  If (abs((Sol_Az(J,K)-FV_GEO)/FV_GEO) .gt.
     &                Rel_Equality_Eps .AND.
     &                abs((Sen_Az(J,K)-FV_GEO)/FV_GEO) .gt.
     &                Rel_Equality_Eps) Then

                      Rel_Az(J,K) = abs(Sen_Az(J,K) - Sol_Az(J,K))
                      Rel_Az(J,K) = amod(Rel_Az(J,K),360.0)

                      if (Rel_Az(J,K) .gt. 180.0)
     &                   Rel_Az(J,K) = 360.0 - Rel_Az(J,K)
                  end if

   55          Continue
   50          Continue

C Retrieve height above ellipsoid

            ELSE
               L = 0

               Do 60 K = 1, Data_Size(2)
               Do 65 J = 1, Data_Size(1)
                  L = L + 1
                  Height(J,K) = FV_GEO
                  if (count(L) .ne. FV_Height) Height(J,K) = count(L)
   65          continue
   60          continue
C
            END IF
         END IF
  999 continue
      if (.not.error_flag) GetGeo_MOD05_CORR = 0
      Return
      END



C*********************************************************************
      SUBROUTINE CldMsk_Info_CORR(Modfil,No_Byte,No_Byte_QA,Fmax,Lmax,
     1    Count,QA_Count,CldMsk_1km,LandSea_Flag)

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C
C    Extract cloud mask flags values from MODIS Cloud Mask
C    HDF of typical 200 scan cubes (a granule).
C    Cloud mask QA flags are not used at this moment
C
C !INPUT PARAMETERS:
C
C    No_Byte    Number of byte of cloud mask storage (=6)
C    No_Byte_QA Number of byte of cloud mask storage (=10)
C    Fmax       Maximum frame number per scan line.
C    Lmax       Maximum line number per scan cube.
C    Count      Cloud mask flags
C    QA_Count   Cloud mask QA flags
C
C !OUTPUT PARAMETERS: (All integers)
C
C   C CldMsk(4*x,4*y)      Cloud Mask (250m resolution from 1km resolution)
C   C LandSea_FLag         Land/Sea Flag
C    DayNight_Flag        Day/Night Flag
C    SunGlint_Flag        SunGlint Flag
C    SnowIce_Flag         Snow/Ice Flag
C    LandSea_Flag         Land/Sea Flag
C    Non_CloudOb_Flag     Non-Cloud obstruction (dust) Flag
C    Thin_CirNIR_Flag     Thin Cirrus detected Flag (Solar)
C    Shadow_Flag          Shadow Flag
C    Thin_CirIR_Flag      Thin Cirrus detected Flag (IR)
C    Cloud_SimpIR_Flag    Cloud Flag-Simple Threshold Test
C    High_Cloud_Flag      High Cloud Flag - 1.38 Micron Test
C    Cloud_IRTemp_Flag    Cloud Flag - IR Temperature Difference
C    Cloud_3p75_11_Flag   Cloud Flag - 3.75-11 Micron Test
C    Cloud_VisRat_Flag    Cloud Flag - Visible Ratio Test
C    Cloud_SpatVar_Flag   Cloud Flag - Spatial Variability
C
c!REVISION HISTORY:
c 01/29/98 fhliang
c fixed prolog.
c
c!TEAM-UNIQUE HEADER:
c
c Developed by MODIS Aerosol/Solar Water Vapor Retrieval Team
c GSFC, Greenbelt, MD
c
C !DESIGN NOTES:
C
C   Internals:
C
C      Variables:
C        Dim_Size(3)     Array specifying the size of hdf SDS data array.
C        l,l1,l4,l2      Line and pixel counters p,p1,p4,p2
C        l4_offset       250-m line offset to 1-km cell position
C        p4_offset       250-m pixel offset to 1-km cell position
C        l2_offset       500-m line offset to 1-km cell position
C        p2_offset       500-m pixel offset to 1-km cell position
C        ibit            0-based bit postion (in INTEGER word)
C        count(15000)    A temporary buffer for data of the target array.
C        QA_count(15000) A temporary buffer for data of the target array.
C
C ! WRITTEN By
C
C      Dr. Allen Chu
C      Code 913/SSAI
C      NASA Goddard Space Flight Center
C      Greenbelt, MD 20771

C !END
C-----------------------------------------------------------------------

      IMPLICIT NONE

      INCLUDE 'mapi.inc'

      CHARACTER*80 arrnm,grpnm
      CHARACTER*13 data_type
      INTEGER  I,j,k
      INTEGER  l,l1,l4,l2,l2_offset,l4_offset,p,p1,p2,p2_offset,p4,
     1         p4_offset,ibit
      INTEGER  No_Byte,No_Byte_QA,Fmax,Lmax,index,Ibyte,Dim_Size(3)
     2         ,Modfil(MODFILLEN),LinesPerScancube,Rank
C
C      PARAMETER    (No_Byte=6,No_Byte_QA=10,Fmax=1500,Lmax=10)
C
      BYTE    Count(No_byte*Fmax*Lmax),
     2        QA_Count(No_Byte_QA*Fmax*Lmax)

      INTEGER CldMsk_1km(Fmax,Lmax),LandSea_Flag(Fmax,Lmax)
      LOGICAL error_flag

C      Dim_Size(1)=No_Byte
      Dim_Size(2)=Fmax
      Dim_Size(3)=Lmax
C
C Get 1-km Cloud MASK data
C
C-----------------------------------------------------------------------
C Derive cloud mask for 250-m pixels based on the results of the
C visible cloud test.  Also extract Land/Sea flag at 1-km resolution.
C Begin looping over 1-km pixels.
C-----------------------------------------------------------------------

       I = -5

       Do 80 l1 = 1, Dim_Size(3)
       Do 80 p1 = 1, Dim_Size(2)
C
C The Cloud Mask consists of 6 separate 1-byte words.  Increment
C memory buffer index by 6 to cycle through 1-km frames in a line
C

          I = I + 6

C-----------------------------------------------------------------------
C Examine first byte of cloud mask to determine whether cloud mask
C was even determined.
C (Zero-based bit 0 is 1 for determined, and 0 for not determined)
C
C If cloud mask not determined, set 250-m CldMsk(p4,l4) to -1.
C Also set 1-km LandSea_Flag to -1 (????).
C
C For Version 2, LandSea_Flag may take 4 values:
C [0 (water); 1 (coastal), 2 (desert), 3 (land)]
C-----------------------------------------------------------------------
C

          LandSea_Flag(p1,l1) = ibits(count(I),6,2)

          If (ibits(count(I),0,1) .EQ. 0) Then

             CldMsk_1km(p1,l1) = 0

          Else
C-----------------------------------------------------------------------
C Cloud Mask determined.  Extract value of 1-km LandSea_Flag.
C Also get results (a simple yes/no switch) of 250-m visible cloud test.
C-----------------------------------------------------------------------

             IF(ibits(count(I),1,2) .NE. 0 ) Then

               CldMsk_1km(p1,l1) = 1

             ELSE

               CldMsk_1km(p1,l1) = 0

             ENDIF

          End If

   80  continue

       Return
       END
