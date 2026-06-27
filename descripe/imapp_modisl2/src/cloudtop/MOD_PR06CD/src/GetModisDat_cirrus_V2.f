      INTEGER FUNCTION GetModisDat_MOD06CD(
     &        Modfil_L1B,Modfil_Geo,Modfil_CldMsk,FN_L1B,
     &        Scan_No,Buf_Size1,Buf_Size2,Data_Size,
     &        SatZen,SolZen,RelAz,Height,
     &        Refl_1,Refl_2,Refl_5,Refl_6,Refl_7,Refl_19,
     &        Refl_26,Rad_29,Rad_31,Rad_32,
     &        Un_1,Un_2,Un_5,Un_6,Un_7,Un_19,Un_26,Un_29,
     &        Un_31,Un_32,
     &        Vflag_1,Vflag_2,Vflag_5,Vflag_6,Vflag_7,Vflag_19,
     &        Vflag_26,Vflag_29,Vflag_31,Vflag_32,Buf_Un,Buf_Sa,
     &        CldMsk,LandSea_Flag, Cloud,QA_Cloud, sample1,
     &        sample2, sample5, sample6, sample7,UFQ_Flag)

      IMPLICIT NONE
      INCLUDE 'mapi.inc'
      INCLUDE 'cirrus.inc'
      INCLUDE 'PGS_MODIS_39500.f'
      INCLUDE 'L1B_Reader_V2.1.inc'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION: Retrieves one scan cube of MODIS Level 1B (L1B) sensor
C               data (bands 2, 5, 17, 18, 19), Geolocation data and
C               Cloud Mask data
C
C !INPUT PARAMETERS:
C    INTEGER    Modfil_*         File handle structure for HDF files
C   *INTEGER    Modfil_L1B_1km   File handle structure for L1B_1km files
C   *CHARACTER  L1B_fnamekm      L1B filename in 1 km resolution
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
C Modified by Liqun Ma   02/20/98
C Some unreferenced varaibles and unused includes are moved out
C Prologs are updated
C
C !TEAM-UNIQUE HEADER:
C
C    This software is developed by the MODIS Science Team
C    for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C !REFERENCES AND CREDITS
C
C    MODIFIED BY:
C
C    Dr. Allen Chu            April 1997
C    Code 913
C    NASA Goddard Space Flight Center
C    Greenbelt, MD 20771
C
C    achu@dustdevil.gsfc.nasa.gov
C
C !DESIGN NOTES:
C
C    GetModisDat_MOD06CD checks the return status of all MODIS Application
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
C        GetGeo_MOD06CD
C        ReadCldMsk_MOD06CD
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
      CHARACTER*1   Gain
      CHARACTER*6   Cal_Type
      INTEGER Band_No,Resol

      INTEGER GetGeo_MOD06CD,
     &        Buf_Size1,Buf_Size2,Scan_No,rtn
CDAC - change 3 to MODFILLEN
      INTEGER Buf_Xtrack,Buf_Track,
     &        ip,il,Modfil_L1B(MODFILLEN),
     &        Modfil_Geo(MODFILLEN),Modfil_CldMsk(MODFILLEN),
     &        Data_Size(2),
     &        CldMsk(Buf_Size1,Buf_Size2),
     &        UFQ_Flag(Buf_Size1,Buf_Size2),
     &        LandSea_Flag(Buf_Size1,Buf_Size2)
CDAC

      BYTE    sample1(Buf_Size1,Buf_Size2),sample2(Buf_Size1,Buf_Size2),
     &        sample5(Buf_Size1,Buf_Size2),sample6(Buf_Size1,Buf_Size2),
     &        sample7(Buf_Size1,Buf_Size2)
      LOGICAL error_flag
      REAL pi,dtr,SolZen_Threshold,cossza,x_cossza,Rel_Equality_EPS,
     &     Zero_Eps,
CDAC
     &     SatZen(Buf_Size1,Buf_Size2),SolZen(Buf_Size1,Buf_Size2),
     &     RelAz(Buf_Size1,Buf_Size2),Height(Buf_Size1,Buf_Size2)

      REAL Refl_1(Buf_Size1,Buf_Size2),Refl_2(Buf_Size1,Buf_Size2),
     &     Refl_5(Buf_Size1,Buf_Size2),Refl_6(Buf_Size1,Buf_Size2),
     &     Refl_7(Buf_Size1,Buf_Size2),Refl_26(Buf_Size1,Buf_Size2),
     &     Refl_19(Buf_Size1,Buf_Size2),Rad_29(Buf_Size1,Buf_Size2),
     &     Rad_31(Buf_Size1,Buf_Size2),Rad_32(Buf_Size1,Buf_Size2)

      BYTE Un_1(Buf_Size1,Buf_Size2),Un_2(Buf_Size1,Buf_Size2),
     &     Un_5(Buf_Size1,Buf_Size2),Un_6(Buf_Size1,Buf_Size2),
     &     Un_7(Buf_Size1,Buf_Size2),Un_26(Buf_Size1,Buf_Size2),
     &     Un_19(Buf_Size1,Buf_Size2),Un_29(Buf_Size1,Buf_Size2),
     &     Un_31(Buf_Size1,Buf_Size2),Un_32(Buf_Size1,Buf_Size2)

      BYTE Vflag_1(Buf_Size1,Buf_Size2),Vflag_2(Buf_Size1,Buf_Size2),
     &     Vflag_5(Buf_Size1,Buf_Size2),Vflag_6(Buf_Size1,Buf_Size2),
     &     Vflag_7(Buf_Size1,Buf_Size2),Vflag_26(Buf_Size1,Buf_Size2),
     &     Vflag_19(Buf_Size1,Buf_Size2),Vflag_29(Buf_Size1,Buf_Size2),
     &     Vflag_31(Buf_Size1,Buf_Size2),Vflag_32(Buf_Size1,Buf_Size2)

      BYTE Buf_Un(Buf_Size1,Buf_Size2),Buf_Sa(Buf_Size1,Buf_Size2)
CDAC
      INTEGER Buf_cldmsk,Buf_cldmsk_QA,DS_Dim1_CM,
     2        DS_Dim1_QA,DS_Dim2,DS_Dim3
      PARAMETER(Buf_cldmsk=6,Buf_cldmsk_QA=10)
      BYTE Cloud(Buf_cldmsk,Buf_Size1,Buf_Size2),
     2     QA_Cloud(Buf_cldmsk_QA,Buf_Size1,Buf_Size2)

      PARAMETER (pi = 3.1415927, dtr = pi/180.0)
      PARAMETER (SolZen_Threshold = 85, Rel_Equality_Eps=0.000001,
     &           Zero_Eps = 0.000001)

C Set initial values
      GetModisDat_MOD06CD = -1
      error_flag = .false.

C----------------------------------------------------------------------
C Read Geolocation
C----------------------------------------------------------------------

      rtn = GetGeo_MOD06CD(Modfil_Geo,Scan_No,Buf_Size1,Buf_Size2,
     2      Data_Size,SolZen,SatZen,RelAz,Height)

      IF (rtn .ne. 0) THEN
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     2  'GetGeo_MOD06CD failed','GetModisDat_MOD06CD')
         error_flag = .true.
      END IF

C----------------------------------------------------------------------
C Read Cloud Mask data
C----------------------------------------------------------------------

      CALL Read_CldMsk(Modfil_CldMsk,Scan_No,Buf_cldmsk,
     1                 Buf_cldmsk_QA,Buf_Size1,Buf_Size2,
     2                 DS_Dim1_CM,DS_Dim1_QA,DS_Dim2,DS_Dim3,
     3                 Cloud,QA_Cloud,Error_Flag)


      IF(error_flag) THEN
        CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     1       'Read_CldMsk failed','GetModisDat_MOD04.f')

      ELSE
        CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,
     1       'Read_CldMsk succeeded','GetModisDat_MOD04.f')

      ENDIF


C
C Decode cloud mask or generate our own cloud mask from L1B data
C
      CALL CldMsk_Info_MOD06CD(Modfil_CldMsk,Buf_cldmsk,
     1    Buf_cldmsk_QA,DS_Dim2,DS_Dim3,Cloud,QA_Cloud,CldMsk,
     2    LandSea_Flag,UFQ_Flag)


C--------------------------------------------------------------------
C Read L1B sensor data
C--------------------------------------------------------------------
C
C Read MODIS band 1
C -----------------
      Band_No = 1
      Resol=Resol_Is_1km
      Cal_Type=Cal_Is_Refl
      Gain = ' '
      Buf_Xtrack=Buf_Size1
      Buf_Track=Buf_size2

      Call Read_L1B(Modfil_L1B, FN_L1B, Scan_No, Band_No,
     &              Gain, Resol, Cal_type, Buf_Xtrack, Buf_Track,
     &              Refl_1, Un_1, Buf_Sa, Vflag_1, Data_size,
     &              error_flag)

C
C Read MODIS band 2
C -----------------
      Band_No = 2
      Resol=Resol_Is_1km
      Cal_Type=Cal_Is_Refl
      Gain = ' '
      Buf_Xtrack=Buf_Size1
      Buf_Track=Buf_size2

      Call Read_L1B(Modfil_L1B, FN_L1B, Scan_No, Band_No,
     &              Gain, Resol, Cal_type, Buf_Xtrack, Buf_Track,
     &              Refl_2, Un_2, Buf_Sa, Vflag_2, Data_size,
     &              error_flag)

C
C Read MODIS band 5
C -----------------
      Band_No = 5
      Resol=Resol_Is_1km
      Cal_Type=Cal_Is_Refl
      Gain = ' '
      Buf_Xtrack=Buf_Size1
      Buf_track=Buf_Size2

      Call Read_L1B(Modfil_L1B, FN_L1B, Scan_No, Band_No,
     &              Gain, Resol, Cal_type, Buf_Xtrack, Buf_Track,
     &              Refl_5, Un_5, Buf_Sa, Vflag_5, Data_size,
     &              error_flag)

C
C Read MODIS band 6
C -----------------
      Band_No = 6
      Resol=Resol_Is_1km
      Cal_Type=Cal_Is_Refl
      Gain = ' '
      Buf_Xtrack=Buf_Size1
      Buf_track=Buf_Size2

      Call Read_L1B(Modfil_L1B, FN_L1B, Scan_No, Band_No,
     &              Gain, Resol, Cal_type, Buf_Xtrack, Buf_Track,
     &              Refl_6, Un_6, Buf_Sa, Vflag_6, Data_size,
     &              error_flag)

C
C Read MODIS band 7
C -----------------
      Band_No = 7
      Resol=Resol_Is_1km
      Cal_Type=Cal_Is_Refl
      Gain = ' '
      Buf_Xtrack=Buf_Size1
      Buf_track=Buf_Size2

      Call Read_L1B(Modfil_L1B, FN_L1B, Scan_No, Band_No,
     &              Gain, Resol, Cal_type, Buf_Xtrack, Buf_Track,
     &              Refl_7, Un_7, Buf_Sa, Vflag_7, Data_size,
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

C
C Read MODIS band 26
C ------------------
      Band_No = 26
      Resol=Resol_Is_1km
      Cal_Type=Cal_Is_Refl
      Gain = ' '
      Buf_Xtrack=Buf_Size1
      Buf_track=Buf_Size2

      Call Read_L1B(Modfil_L1B, FN_L1B, Scan_No, Band_No,
     &              Gain, Resol, Cal_type, Buf_Xtrack, Buf_Track,
     &              Refl_26, Un_26, Buf_Sa, Vflag_26, Data_size,
     &              error_flag)

C
C Read MODIS band 29
C ------------------
      Band_No = 29
      Resol=Resol_Is_1km
      Cal_Type=Cal_Is_Rad
      Gain = ' '
      Buf_Xtrack=Buf_Size1
      Buf_track=Buf_Size2

      Call Read_L1B(Modfil_L1B, FN_L1B, Scan_No, Band_No,
     &              Gain, Resol, Cal_type, Buf_Xtrack, Buf_Track,
     &              Rad_29, Un_29, Buf_Sa, Vflag_29, Data_size,
     &              error_flag)

C
C Read MODIS band 31
C ------------------
      Band_No = 31
      Resol=Resol_Is_1km
      Cal_Type=Cal_Is_Rad
      Gain = ' '
      Buf_Xtrack=Buf_Size1
      Buf_track=Buf_Size2

      Call Read_L1B(Modfil_L1B, FN_L1B, Scan_No, Band_No,
     &              Gain, Resol, Cal_type, Buf_Xtrack, Buf_Track,
     &              Rad_31, Un_31, Buf_Sa, Vflag_31, Data_size,
     &              error_flag)

C
C Read MODIS band 32
C ------------------
      Band_No = 32
      Resol=Resol_Is_1km
      Cal_Type=Cal_Is_Rad
      Gain = ' '
      Buf_Xtrack=Buf_Size1
      Buf_track=Buf_Size2

      Call Read_L1B(Modfil_L1B, FN_L1B, Scan_No, Band_No,
     &              Gain, Resol, Cal_type, Buf_Xtrack, Buf_Track,
     &              Rad_32, Un_32, Buf_Sa, Vflag_32, Data_size,
     &              error_flag)

CDAC
C---
C---Modified by Bo-Cai Gao in May 2004 based on Rich Hucek's suggestions
C---
       Do il = 1, Data_Size(2)
         Do ip = 1, Data_Size(1)
          If (Vflag_1(Ip,il)  .lt. 0) Refl_1(Ip,il)  = FV_L1B
          If (Vflag_2(Ip,il)  .lt. 0) Refl_2(Ip,il)  = FV_L1B
          If (Vflag_5(Ip,il)  .lt. 0) Refl_5(Ip,il)  = FV_L1B
          If (Vflag_6(Ip,il)  .lt. 0) Refl_6(Ip,il)  = FV_L1B
          If (Vflag_7(Ip,il)  .lt. 0) Refl_7(Ip,il)  = FV_L1B
          If (Vflag_26(Ip,il) .lt. 0) Refl_26(Ip,il) = FV_L1B
         End do
       End do
C---End of modification

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

*/  Modified by JC Guu  10/18/96
*/  The value compared with a practical 0
*/  is normalized.

         If ( abs((SolZen(ip,il)-FV_GEO)/FV_GEO) .LT.
     &        Rel_Equality_Eps) Then
            Refl_1(ip,il)  = FV_L1B
            Refl_2(ip,il)  = FV_L1B
            Refl_5(ip,il)  = FV_L1B
            Refl_6(ip,il)  = FV_L1B
            Refl_7(ip,il)  = FV_L1B
            Refl_19(ip,il) = FV_L1B
            Refl_26(ip,il) = FV_L1B

         Else If (abs(cossza) .LT. Zero_Eps) Then
            Refl_1(ip,il)  = FV_L1B
            Refl_2(ip,il)  = FV_L1B
            Refl_5(ip,il)  = FV_L1B
            Refl_6(ip,il)  = FV_L1B
            Refl_7(ip,il)  = FV_L1B
            Refl_19(ip,il) = FV_L1B
            Refl_26(ip,il) = FV_L1B

         Else If (SolZen(ip,il) .LT. SolZen_Threshold) Then

            If ( abs((Refl_1(ip,il)-FV_L1B)/FV_L1B)
     &         .GT. Rel_Equality_Eps) then
               Refl_1(ip,il) = x_cossza*Refl_1(ip,il)
C---            Else
C---               Refl_1(ip,il) = FV_L1B
            End If

            If ( abs((Refl_2(ip,il)-FV_L1B)/FV_L1B)
     &         .GT. Rel_Equality_Eps) then
               Refl_2(ip,il) = x_cossza*Refl_2(ip,il)
C---            Else
C---               Refl_2(ip,il) = FV_L1B
            End If

            If ( abs((Refl_5(ip,il)-FV_L1B)/FV_L1B)
     &         .GT. Rel_Equality_Eps) then
               Refl_5(ip,il) = x_cossza*Refl_5(ip,il)
C---            Else
C---               Refl_5(ip,il) = FV_L1B
            End If

            If ( abs((Refl_6(ip,il)-FV_L1B)/FV_L1B)
     &         .GT. Rel_Equality_Eps) then
               Refl_6(ip,il) = x_cossza*Refl_6(ip,il)
C---            Else
C---               Refl_6(ip,il) = FV_L1B
            End If

            If ( abs((Refl_7(ip,il)-FV_L1B)/FV_L1B)
     &         .GT. Rel_Equality_Eps) then
               Refl_7(ip,il) = x_cossza*Refl_7(ip,il)
C---            Else
C---               Refl_7(ip,il) = FV_L1B
            End If

            If (abs((Refl_19(ip,il)-FV_L1B)/FV_L1B)
     &         .GT. Rel_Equality_Eps) then
               Refl_19(ip,il) = x_cossza*Refl_19(ip,il)
C---            Else
C---               Refl_19(ip,il) = FV_L1B
            End If

            If (abs((Refl_26(ip,il)-FV_L1B)/FV_L1B)
     &         .GT. Rel_Equality_Eps) then
               Refl_26(ip,il) = x_cossza*Refl_26(ip,il)
C---            Else
C---               Refl_26(ip,il) = FV_L1B
            End If

         Else
            Refl_1(ip,il)  = FV_L1B
            Refl_2(ip,il)  = FV_L1B
            Refl_5(ip,il)  = FV_L1B
            Refl_6(ip,il)  = FV_L1B
            Refl_7(ip,il)  = FV_L1B
            Refl_19(ip,il) = FV_L1B
            Refl_26(ip,il) = FV_L1B

         End If

   45 continue
   40 continue

      Do il = 1, Data_Size(2)
        Do ip = 1, Data_Size(1)

         If (SolZen(ip,il) .GE. SolZen_Threshold) Then
            Refl_1(ip,il)  = FV_L1B
            Refl_2(ip,il)  = FV_L1B
            Refl_5(ip,il)  = FV_L1B
            Refl_6(ip,il)  = FV_L1B
            Refl_7(ip,il)  = FV_L1B
            Refl_19(ip,il) = FV_L1B
            Refl_26(ip,il) = FV_L1B
         End if

        End do
      End do

      if (.not.error_flag) GetModisDat_MOD06CD = 0

      Return
      END

C***********************************************************************
      INTEGER FUNCTION GetGeo_MOD06CD(Modfil,Scan_No,Buf_Size1,
     &         Buf_Size2,Data_Size,Sol_Zen,Sen_Zen,Rel_Az,Height)

      IMPLICIT NONE
      SAVE
      INCLUDE 'mapi.inc'
      INCLUDE 'hdf.inc'
      INCLUDE 'cirrus.inc'
      INCLUDE 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION: Retrieves MODIS Geolocation data from an HDF target
C              array of typically 100 scan cubes (a granule).  This
C              routine is customized for MOD04 land algorithm.
C
C !INPUT PARAMETERS:
C    INTEGER Modfil(3)   SD_ID, File_ID and File access type.
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
c Revision 1.2  1995/10/18  22:41:24  rhucek
c revised GetGeo_MOD06CD function to reflect data format of standard
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
C    GetGeo_MOD06CD checks the return status of all MODIS Application
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
C      FV_*                       (*****_inc.f)
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
C       GetGeo_MOD06CD     0 for success / -1 if failed
C
C !END
C-----------------------------------------------------------------------

      CHARACTER*13 grpnm,dtype
      CHARACTER*30 arrnm
      CHARACTER*80 msgbuf,msgbuf1
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

      GetGeo_MOD06CD=-1
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
     2  'Invalid SD ID/file access type','GetGeo_MOD06CD')
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
     2     'GMARDM failed','GetGeo_MOD06CD')
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
     &      'GetGeo_MOD06CD')
            error_flag = .true.
         END IF

C Read HDF target array into 'count' buffer
         rtn = GMAR(Modfil,arrnm,grpnm,Start,Data_Size,count)

         IF (rtn.ne. MAPIOK) THEN
            CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     &      'GMAR failed','GetGeo_MOD06CD')
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
C---Modified by Bo-Cai Gao in May 2004 to correct the conversion factor error
C---     &               Sol_Zen(J,K) = count(L)*0.018/pi
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
C---Modified by Bo-Cai Gao in May 2004 to correct the conversion factor error
C---     &               Sen_Zen(J,K) = count(L)*0.018/pi
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
     &               Sol_Az(J,K) = 180.0 + count(L)*0.01
C---Modified by Bo-Cai Gao in May 2004 to correct the conversion factor error
C---     &               Sol_Az(J,K) = 180.0 + count(L)*0.018/pi
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
C---Modified by Bo-Cai Gao in May 2004 to correct the conversion factor error
C---     &               Sen_Az(J,K) = count(L)*0.018/pi

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

            END IF
         END IF
  999 continue
      if (.not.error_flag) GetGeo_MOD06CD = 0
      Return
      END

C*********************************************************************

      INTEGER FUNCTION ReadAnc_MOD06CD(MODFIL_ANC,Buf_Size1,
     &        Buf_Size2,scan,Data_Size,sfc_temp)

       implicit none
       include 'mapi.inc'
       include 'hdf.inc'
       include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION: This subroutine reads the ancillary data NMC/WISCONSIN
C               in HDF format
C
C !INPUT PARAMETERS:
C  integer MODFIL_ANC   File handle structure for HDF files.
C  integer Buf_Size1/2  Size of first/second index of surface
C                       temperature data buffer as dimensioned in the
C                       calling program
C  integer scan         Scan number
C
C !OUTPUT PARAMETERS:
C  integer Data_Size    Acture size of surface temperature data
C  real    sfc_temp(x,y) Buffer used to hold surface temperature data
C                        x = Buf_Size1, y = Buf_Size2
C
C !TEAM-UNIQUE HEADER:
C    This software is developed by the MODIS Science Data Support
C    Team for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C
C !REVISION HISTORY:
C
C !REFERENCES AND CREDITS:
C    Written by Vicky Lin            April 1996
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd, Seabrook MD 20706
C
C    vlin@ltpmail.gsfc.nasa.gov
C
C !DESIGN NOTES:
C
C !END
C-----------------------------------------------------------------------

      character*40 DATA_TYPE,ARRNM,GRPNM,msgbuf1
      character*80 msgbuf
      integer L,rank,ip,il,rtn,scan,MaxScan_No,Buf_Size1,Buf_Size2,
     &        Data_Size(2),start(2),count(2),MODFIL_ANC(MODFILLEN)
      logical error_flag
      real sfcT(15000),sfc_temp(Buf_Size1,Buf_Size2)

C Initialization

      ReadAnc_MOD06CD = -1
      error_flag = .false.
      arrnm='surface_temperature'
      grpnm=' '
      rank=2
      start(1) = 0
      start(2) = 10 * (scan-1)
      do 10 il = 1, Buf_Size2
      do 15 ip = 1, Buf_Size1
         sfc_temp(ip,il) = 0.0
   15 continue
   10 continue

C Check for valid input
      IF (Modfil_ANC(1).le.0.or.Modfil_ANC(3).ne.DFACC_READ) THEN
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     2  'Invalid Modfil_ANC','ReadAnc_MOD06CD')
         error_flag = .true.
      END IF

C Retrieve the rank, dimension and data type of surface temperature
      rtn = GMARDM(MODFIL_ANC,arrnm,grpnm,data_type,rank,Data_Size)
      if (rtn.ne.MAPIOK) THEN
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     2   'GMARDM failed','ReadAnc_MOD06CD')
         error_flag = .true.
      end if

C Additional input check
      if (Buf_Size1.lt.Data_Size(1)) then
         call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     &   'Buf_Size1 too small','ReadAnc_MOD06CD')
          error_flag = .true.
      end if

      MaxScan_No = Data_Size(2)/Buf_Size2
      if (scan.lt.1 .OR. scan.gt.MaxScan_No) then
         write(msgbuf1,'(i4)') MaxScan_No
         call Concatenate('Scan_No out of bounds : 1 - ',msgbuf1,
     &   msgbuf)
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,msgbuf,
     &   'ReadAnc_MOD06CD')
         error_flag = .true.
      end if

C Read HDF target array into 'count' buffer
      count(1) = Data_Size(1)
      count(2) = Buf_Size2
      rtn = GMAR(MODFIL_ANC,arrnm,grpnm,start,count,sfcT)
      if (rtn.ne. MAPIOK) then
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     &   'GMAR failed','ReadAnc_MOD06CD')
         error_flag = .true.
      end if

C Convert HDF 2-byte integer data to 4-byte real value

      if (.not.error_flag) then

         L = 0
         do 30 il = 1, Buf_Size2
         do 35 ip = 1, Buf_Size1
            L = L + 1
            sfc_temp(ip,il) = sfcT(L)
   35    continue
   30    continue

C ReadAnc_MOD06CD completed successfully !
         ReadAnc_MOD06CD = 0

      end if

      return
      end
C***********************************************************************
       SUBROUTINE CldMsk_Info_MOD06CD(Modfil,No_Byte,No_Byte_QA,Fmax,
     1    Lmax,Count,QA_Count,CldMsk_1km,LandSea_Flag,UFQ_Flag)
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
C !OUTPUT PARAMETERS:(All integers)
C
C   C CldMsk(4*x,4*y)      Cloud Mask (250m resolution from 1km resolution)
C   C LandSea_FLag         Land/Sea Flag
C    UFQ_Flag             Unobstructed FOV Quality Flag
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
C
C !Revision History:
C
C
C
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
c!Team-unique Header:
c
C   This software was developed by the MODIS Science Data Support Team
C   (SDST) for the National Aeronautics and Space Administration,
C   Goddard Space Flight Center, under contract NAS5-32373.
C
C
C !REFERENCES and CREDITS:
C   WRITTEN By
C
C      Dr. Allen Chu
C      Code 913/SSAI
C      NASA Goddard Space Flight Center
C      Greenbelt, MD 20771
C
C !END
C-----------------------------------------------------------------------
C
      IMPLICIT NONE
      INCLUDE 'mapi.inc'
      INTEGER  I
      INTEGER  l1,p1
      INTEGER  No_Byte,No_Byte_QA,Fmax,Lmax,Dim_Size(3)
     2         ,Modfil(MODFILLEN)

      BYTE    Count(No_byte*Fmax*Lmax),
     2        QA_Count(No_Byte_QA*Fmax*Lmax)

      INTEGER CldMsk_1km(Fmax,Lmax),LandSea_Flag(Fmax,Lmax),
     1        UFQ_Flag(Fmax,Lmax)

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
C
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
C If not determined, determined by Lorraine's scheme (Added later)
C
          If (ibits(count(I),0,1) .EQ. 0) Then

             LandSea_Flag(p1,l1) = -1

          Else
C-----------------------------------------------------------------------
C Cloud Mask determined.  Extract value of 1-km LandSea_Flag.
C Also get results (a simple yes/no switch) of 250-m visible cloud test.
C-----------------------------------------------------------------------

             UFQ_Flag(p1,l1) = ibits(count(I),1,2)
             LandSea_Flag(p1,l1) = ibits(count(I),6,2)

             IF(UFQ_Flag(p1,l1).NE.1.OR.UFQ_Flag(p1,l1).NE.0) THEN

               CldMsk_1km(p1,l1) = 1

             ENDIF

          End If

   80  continue

       Return
       END
