      INTEGER FUNCTION GetModisDat_MOD04(NUMSQ,
     &     Modfil_L1B_1km,Modfil_L1B_500,Modfil_L1B_250,
     &     Modfil_Geo,Modfil_CldMsk,
     &     FN_L1B_1km,FN_L1B_500,FN_L1B_250,
     &     ScanCube_No,Buf_Size1,Buf_Size2,
     &     Data_Size,Lat,Lon,SatZen,SatAz,SolZen,SolAz,RelAz,Height,
     &     Refl_1,Refl_2,Refl_3,Refl_4,Refl_5,Refl_6,Refl_7,Refl_8,
     &     Refl_9, Refl_12,Refl_13,Refl_16,Refl_26,Rad_20,Rad_31,Rad_32,
     &     Refl_15,Unc_1,Unc_2,Unc_3,Unc_4,Unc_5,Unc_6,Unc_7,Unc_8,Unc_9,
     &     Unc_12,Unc_13,Unc_16,Unc_26,Unc_20,Unc_31,Unc_32,Unc_15,
     &     Vflag_1,Vflag_2,Vflag_3,Vflag_4,Vflag_5,Vflag_6,Vflag_7,Vflag_8,
     &     Vflag_9,Vflag_12,Vflag_13,Vflag_16,Vflag_26,Vflag_20,
     &     Vflag_31,Vflag_32,Vflag_15,Buf_Un,Buf_Sa_1km,Buf_Sa_500,
     &     Buf_Sa_250,CldMsk_250,CldMsk_500,CldMsk_1km,DET_Flag,
     &     UFQ_Flag,DayNight_Flag,SunGlint_Flag,SnowIce_Flag,
     &     LandSea_Flag,Non_CloudOb_Flag,Thin_CirNIR_Flag,
     &     Shadow_Flag,Thin_CirIR_Flag,Cloud_SimpIR_Flag,
     &     High_Cloud_Flag,Cloud_IRTemp_Flag,Cloud_3p75_11_Flag,
     &     Cloud_VisRat_Flag,Cloud_SpatVar_Flag,Cloud,QA_CLOUD,Iscan)

      IMPLICIT NONE

      INCLUDE 'mapi.inc'
      INCLUDE 'mod04.inc'
      INCLUDE 'PGS_MODIS_39500.f'
      INCLUDE 'L1B_Reader_V2.1.inc'

C-----------------------------------------------------------------------
C!F77
C
C !DESCRIPTION:
C
C    Retrieves one scan cube of MODIS Level 1B (L1B) sensor data
C    (1-36 bands), Geolocation data and Cloud Mask data from an
C    HDF target array of typically 100 scan cubes (a granule).
C    In definitions below, x=Buf_Size1 and  y=Buf_Size2
C
C !INPUT PARAMETERS:
C
C    INTEGER Modfil_L1B_1km(MODFILLEN)
C                          Array containing File ID and access type
C                          for L1B data of 1-km resolution.
C    INTEGER Modfil_L1B_500(MODFILLEN)
C                          Array containing File ID and access type
C                          for L1B data of 500 m resolution.
C    INTEGER Modfil_L1B_250(MODFILLEN)
C                          Array containing File ID and access type
C                          for L1B data of 250 m resolution.
C    INTEGER Modfil_Geo(MODFILLEN)
C                          Array containing File ID and access type
C                          Geolocation data.
C    INTEGER Modfil_CldMsk(MODFILLEN)
C                          Array containing File ID and access type
C                          Cloud Mask data.
C    INTEGER ScanCube_No  1-based Scancube number (defined as 10 1-km
C                         scan lines) for the granule
C    INTEGER Buf_Size1    The buffer size allocated for the columns of
C                         1-km spatial data.
C    INTEGER Buf_Size2    The buffer size allocated for the rows of
C                         1-km spatial data.
C
C !OUTPUT PARAMETERS:
C    INTEGER  Data_Size(2) Size of retrieved 1-km data arrays.
C    INTEGER  Cloud(4x,4y) Storage buffer to hold Cloud/No Cloud Flag
C                          Index 1 for columns; Index 2 for rows.
C    INTEGER  LandSea_Flag Array containing Land/Sea flag retrieved
C             (x,y)        from Cloud Mask product.  0 water: 1 coastal;
C                          2 wetland; 3 land; -1 invalid.
C    REAL     Lat(x,y)     Storage buffer to hold pixel latitude.
C    REAL     Lon(x,y)     Storage buffer to hold pixel longitude.
C    REAL     SatZen(x,y)  Storage buffer to hold satellite zenith angle
C    REAL     SolZen(x,y)  Storage buffer to hold solar zenith angle.
C    REAL     RelAz(x,y)   Storage buffer to hold relative azimuth
C                          angle between solar vector in the forward
C                          scattering direction and satellite vector.
C    REAL     Refl_1(4x,4y) Storage buffer for band 1 reflectances.
C    REAL     Refl_2(4x,4y) Storage buffer for band 2 reflectances.
C    REAL     Refl_3(2x,2y) Storage buffer for band 3 reflectances.
C    REAL     Refl_4(2x,2y) Storage buffer for band 4 reflectances.
C    REAL     Refl_5(2x,2y) Storage buffer for band 5 reflectances.
C    REAL     Refl_6(2x,2y) Storage buffer for band 6 reflectances.
C    REAL     Refl_7(2x,2y) Storage buffer for band 7 reflectances.
C    REAL     Refl_9(x,y)   Storage buffer for band 9 reflectances.
C    REAL     Refl_12(x,y)  Storage buffer for band 12 reflectances.
C    REAL     Refl_13(x,y)  Storage buffer for band 13 reflectances.
C    REAL     Refl_16(x,y)  Storage buffer for band 16 reflectances.
C    REAL     Refl_26(x,y)  Storage buffer for band 26 reflectances.
C    REAL     Rad_20(x,y)   Storage buffer for band 20 radiances.
C    REAL     Rad_31(x,y)   Storage buffer for band 31 radiances.
C    REAL     Rad_32(x,y)   Storage buffer for band 32 radiances.
C    BYTE     Unc_1(4x,4y)  Storage buffer for band 1 uncertainty.
C    BYTE     Unc_2(4x,4y)  Storage buffer for band 2 uncertainty.
C    BYTE     Unc_3(2x,2y)  Storage buffer for band 3 uncertainty.
C    BYTE     Unc_4(2x,2y)  Storage buffer for band 4 uncertainty.
C    BYTE     Unc_5(2x,2y)  Storage buffer for band 5 uncertainty.
C    BYTE     Unc_6(2x,2y)  Storage buffer for band 6 uncertainty.
C    BYTE     Unc_7(2x,2y)  Storage buffer for band 7 uncertainty.
C    BYTE     Unc_9(x,y)    Storage buffer for band 9 uncertainty.
C    BYTE     Unc_12(x,y)   Storage buffer for band 12 uncertainty.
C    BYTE     Unc_13(x,y)   Storage buffer for band 13 uncertainty.
C    BYTE     Unc_16(x,y)   Storage buffer for band 16 uncertainty.
C    BYTE     Unc_26(x,y)   Storage buffer for band 26 uncertainty.
C    BYTE     Unc_20(x,y)   Storage buffer for band 20 uncertainty.
C    BYTE     Unc_31(x,y)   Storage buffer for band 31 uncertainty.
C    BYTE     Unc_32(x,y)   Storage buffer for band 32 uncertainty.
C    BYTE     Vflag_1(4x,4y)Storage buffer for band 1 valid flags.
C    BYTE     Vflag_2(4x,4y)Storage buffer for band 2 valid flags.
C    BYTE     Vflag_3(2x,2y)Storage buffer for band 3 valid flags.
C    BYTE     Vflag_4(2x,2y)Storage buffer for band 4 valid flags.
C    BYTE     Vflag_5(x,y)  Storage buffer for band 5 valid flags.
C    BYTE     Vflag_6(x,y)  Storage buffer for band 6 valid flags.
C    BYTE     Vflag_7(x,y)  Storage buffer for band 7 valid flags.
C    BYTE     Vflag_9(x,y)  Storage buffer for band 8 valid flags.
C    BYTE     Vflag_12(x,y) Storage buffer for band 9 valid flags.
C    BYTE     Vflag_13(x,y) Storage buffer for band 12 valid flags.
C    BYTE     Vflag_16(x,y) Storage buffer for band 16 valid flags.
C    BYTE     Vflag_26(x,y) Storage buffer for band 26 valid flags.
C    BYTE     Vflag_31(x,y) Storage buffer for band 31 valid flags.
C    BYTE     Vflag_32(x,y) Storage buffer for band 32 valid flags.
C
C!REVISION HISTORY:
C $Log: GetModisDat_MOD04_V1.f,v $
c 01/28/98 fhliang
c fixed prolog.
c
c Revision 1.3  1996/07/26  12:02:57  rhucek
c Retrieved 1-km Land/Sea flag from MODIS Cloud Mask product and
c return it as function argument.
c Replace L1B reader with new version that identifies missing packets
c data.  Values of -999.0 are returned for each sensor observation that
c is designated as missing in the L1B product.
c Revised the computation of reflectance for all sensor bands in FUNCTION
c GetModisDat_MOD04.  Before normalization by the cosine of the solar
c zenith angle is performed, a more careful check on the presence of "FILL"
c values in both "solar zenith angle (obtained from MOD03)" and L1B
c radiances is made.
c
c Revision 1.2  1996/03/29  17:42:01  vlin
c Checked "FILL values" in function GetGeo_MOD04
c
c Revision 1.1  1996/03/27  16:31:12  vlin
c Initial revision
C
C!TEAM-UNIQUE HEADER:
C
C    This software is developed by the MODIS Science Data Support
C    Team for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C REFERENCES AND CREDITS
C
C    Modified by:
C    Dr. Allen Chu     May 1997
C    Code 913/SSAI
C    NASA Goddard Flight Center
C    Greenbelt, MD 20771
C
C    Written by:
C    Vicky Lin         March 1996
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd, Seabrook, MD 20706
C    vlin@ltpmail.gsfc.nasa.gov
C
C DESIGN NOTES:
C
C    This subroutine checks the return status of all MODIS Application
C    Program Interface (M-API) function calls.  A successful M-API
C    function call is indicated by a return value mapiok.  If a M-API
C    function call is not successful, a warning message is written
C    to the LogStatus file by passing a mnemonic warning message
C    (i.e. .._W_..) to the function MODIS_SMF_SETDYNAMICMSG.
C
C   Externals:
C
C    MODIS_F_GENERIC            (MODIS_39500.f)
C    mapiok                     (mapi.inc)
C
C   Internals:
C
C    Subroutine:
C      Read_L1B
C    Function:
C      GetGeo_MOD04
C      ReadCldMsk_MOD04
C
C    Local Variables:
C
C      SDS_DIM(2)  Array containing spatial dimension sizes (number
C                  of frames/columns and lines/rows in an HDF SDS
C                  array
C      RTN         Return value
C      Band_No     MODIS band number (1-36)
C      ip1,ip2,isp2,iep2,ip4,isp4,iep4,il1,il2,isl2,iel2,il4,isl4,iel4
C                  Reference indices into 250m, 500m & 1 km data arrays
C      idebug      Debug flag: if 1, print debug data
C      SolIrr_21   Solar Spectral Irradiance for band 21
C      cossza      Cosine of Solar Zenith Angle
C      x_cossza    Inverse of cossza
C      SolZen_Threshold  Solar Zenith Angle Processing Threshold
C      MAX_EV_Frames Variable containing the value of L1B product
C                    metadata attribute "Max Earth View Frames"
C      dtype       Character variable containing data type of
C                  MAX_EV_Frames
C      att_N       Character variable to hold HDF attribute names
C      nms         Number of values contained in HDF attribute
C
C !END
C----------------------------------------------------------------------

      character *(*) FN_L1B_1km,FN_L1B_500,FN_L1B_250
      CHARACTER*100 att_N,dtype
CDAC 10/3/97
      CHARACTER*1   Gain
      CHARACTER*6   Cal_Type
      INTEGER Band_No,Resol
CDAC 10/3/97
      INTEGER GetGeo_MOD04,ReadCldMsk_MOD04
      INTEGER Modfil_L1B_1km(MODFILLEN),Modfil_L1B_500(MODFILLEN),
     1        Modfil_L1B_250(MODFILLEN),Modfil_Geo(MODFILLEN),
     2        Modfil_CldMsk(MODFILLEN),Buf_Size1,Buf_Size2,
     3        Buf_Xtrack,Buf_Track,I_flag,J_flag,
     4        Data_Size(2),SDS_DIM(2),Buf_Dims1,Buf_Dims2,
     5        ScanCube_No,RTN,NUMSQ
      INTEGER ip1,ip2,isp2,iep2,ip4,isp4,iep4,il1,il2,isl2,iel2,
     2        il4,isl4,iel4,i,j,nms,MAX_EV_Frames,iscan 
C
C Cloud mask arrays
C Old
C      INTEGER Cloud(4*Buf_Size1,4*Buf_Size2)
C      INTEGER LandSea_Flag(Buf_Size1,Buf_Size2)
C New
      INTEGER CldMsk_250(4*Buf_Size1,4*Buf_Size2),
     1    CldMsk_500(2*Buf_Size1,2*Buf_Size2),
     2    CldMsk_1km(Buf_Size1,Buf_Size2),
     2    DET_Flag(Buf_Size1,Buf_Size2),
     3    UFQ_Flag(Buf_Size1,Buf_Size2),
     4    DayNight_Flag(Buf_Size1,Buf_Size2),
     5    SunGlint_Flag(Buf_Size1,Buf_Size2),
     6    SnowIce_Flag(Buf_Size1,Buf_Size2),
     7    LandSea_Flag(Buf_Size1,Buf_Size2),
     8    Non_CloudOb_Flag(Buf_Size1,Buf_Size2),
     9    Thin_CirNIR_Flag(Buf_Size1,Buf_Size2),
     1    Shadow_Flag(Buf_Size1,Buf_Size2),
     2    Thin_CirIR_Flag(Buf_Size1,Buf_Size2),
     3    Cloud_SimpIR_Flag(Buf_Size1,Buf_Size2),
     4    High_Cloud_Flag(Buf_Size1,Buf_Size2),
     5    Cloud_IRTemp_Flag(Buf_Size1,Buf_Size2),
     6    Cloud_3p75_11_Flag(Buf_Size1,Buf_Size2),
     7    Cloud_VisRat_Flag(Buf_Size1,Buf_Size2),
     8    Cloud_SpatVar_Flag(Buf_Size1,Buf_Size2),k

      INTEGER DS_Dim1_CM,DS_Dim1_QA,DS_Dim2,DS_Dim3

      BYTE Cloud(Buf_cldmsk,Buf_Size1,Buf_Size2),
     2     QA_Cloud(Buf_cldmsk_QA,Buf_Size1,Buf_Size2)
C
C ?????
C
      INTEGER i_detector,ip,nnline
C
C L1B data arrays
C
C      REAL SolIrr_21,pi,dtr,SolZen_Threshold,cossza,x_cossza,
      REAL SolIrr_21,SolZen_Threshold,cossza,x_cossza,
     &     Zero_Eps,Rel_Equality_Eps,
     &     Lat(Buf_Size1,Buf_Size2),Lon(Buf_Size1,Buf_Size2),
     &     SatZen(Buf_Size1,Buf_Size2),SolZen(Buf_Size1,Buf_Size2),
     &     SatAz(Buf_Size1,Buf_Size2),SolAz(Buf_Size1,Buf_Size2),
     &     RelAz(Buf_Size1,Buf_Size2),Height(Buf_Size1,Buf_Size2),
     &     Refl_1(4*Buf_Size1,4*Buf_Size2),
     &     Refl_2(4*Buf_Size1,4*Buf_Size2),
     &     Refl_3(2*Buf_Size1,2*Buf_Size2),
     &     Refl_4(2*Buf_Size1,2*Buf_Size2),
     &     Refl_5(2*Buf_Size1,2*Buf_Size2),
     &     Refl_6(2*Buf_Size1,2*Buf_Size2),
     &     Refl_7(2*Buf_Size1,2*Buf_Size2),
     &     Refl_8(Buf_Size1,Buf_Size2),
     &     Refl_9(Buf_Size1,Buf_Size2),
     &     Refl_12(Buf_Size1,Buf_Size2),
     &     Refl_13(Buf_Size1,Buf_Size2),
     &     Refl_15(Buf_Size1,Buf_Size2),
     &     Refl_16(Buf_Size1,Buf_Size2),
     &     Refl_26(Buf_Size1,Buf_Size2),
     &     Rad_20(Buf_Size1,Buf_Size2),
     &     Rad_31(Buf_Size1,Buf_Size2),
     &     Rad_32(Buf_Size1,Buf_Size2)
          

      BYTE Unc_1(4*Buf_Size1,4*Buf_Size2),
     &     Unc_2(4*Buf_Size1,4*Buf_Size2),
     &     Unc_3(2*Buf_Size1,2*Buf_Size2),
     &     Unc_4(2*Buf_Size1,2*Buf_Size2),
     &     Unc_5(2*Buf_Size1,2*Buf_Size2),
     &     Unc_6(2*Buf_Size1,2*Buf_Size2),
     &     Unc_7(2*Buf_Size1,2*Buf_Size2),
     &     Unc_8(Buf_Size1,Buf_Size2),
     &     Unc_9(Buf_Size1,Buf_Size2),
     &     Unc_12(Buf_Size1,Buf_Size2),
     &     Unc_13(Buf_Size1,Buf_Size2),
     &     Unc_16(Buf_Size1,Buf_Size2),
     &     Unc_26(Buf_Size1,Buf_Size2),
     &     Unc_20(Buf_Size1,Buf_Size2),
     &     Unc_31(Buf_Size1,Buf_Size2),
     &     Unc_32(Buf_Size1,Buf_Size2),
     &     Unc_15(Buf_Size1,Buf_Size2) 

      BYTE Vflag_1(4*Buf_Size1,4*Buf_Size2),
     &     Vflag_2(4*Buf_Size1,4*Buf_Size2),
     &     Vflag_3(2*Buf_Size1,2*Buf_Size2),
     &     Vflag_4(2*Buf_Size1,2*Buf_Size2),
     &     Vflag_5(2*Buf_Size1,2*Buf_Size2),
     &     Vflag_6(2*Buf_Size1,2*Buf_Size2),
     &     Vflag_7(2*Buf_Size1,2*Buf_Size2),
     &     Vflag_8(Buf_Size1,Buf_Size2),
     &     Vflag_9(Buf_Size1,Buf_Size2),
     &     Vflag_12(Buf_Size1,Buf_Size2),
     &     Vflag_13(Buf_Size1,Buf_Size2),
     &     Vflag_16(Buf_Size1,Buf_Size2),
     &     Vflag_26(Buf_Size1,Buf_Size2),
     &     Vflag_20(Buf_Size1,Buf_Size2),
     &     Vflag_31(Buf_Size1,Buf_Size2),
     &     Vflag_32(Buf_Size1,Buf_Size2),
     &     Vflag_15(Buf_Size1,Buf_Size2)

      BYTE Buf_Un(Buf_Size1,Buf_Size2),
     &     Buf_Sa_1km(Buf_Size1,Buf_Size2),
     &     Buf_Sa_500(2*Buf_Size1,2*Buf_Size2),
     &     Buf_Sa_250(4*Buf_Size1,4*Buf_Size2)
            integer numnum1,numnum2,numnum3

      LOGICAL Error_Flag

C      PARAMETER (pi=3.1415927, dtr=pi/180.0, SolZen_Threshold=80)
      PARAMETER (SolZen_Threshold=80)
      PARAMETER (Rel_Equality_Eps=0.000001, Zero_Eps=0.000001)

C Set values of the extraterrestrial radiance at an earth-sun distance
C of 1 astronomical unit, and at the nominal midpoint frequencies of
C MODIS band 21.
      PARAMETER (SolIrr_21 = 9.82)
       SAVE 


C
C----------------------------------------------------------------------
C Read Geolocation
C----------------------------------------------------------------------

       if (GetGeo_MOD04(Modfil_Geo,ScanCube_No,
     2           Buf_Size1,Buf_Size2,SDS_DIM,
     3           Lat,Lon,SolZen,SolAz,SatZen,SatAz,RelAz,Height,
     4           LandSea_Flag).ne.0)
     4   CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,
     5  'GetGeo_MOD04 failed','GetModisDat_MOD04')


C----------------------------------------------------------------------
C Read L1B sensor data
C----------------------------------------------------------------------
C Read MODIS band 1
C -----------------
      Band_No=1
      Resol=Resol_Is_250m
      Cal_Type=Cal_Is_Refl
      Gain= ' '
      Buf_Xtrack = 4*Buf_Size1
      Buf_Track = 4*Buf_Size2

      Call Read_L1B(Modfil_L1B_250, FN_L1B_250, ScanCube_No, Band_No,
     &              Gain, Resol, Cal_type, Buf_Xtrack, Buf_Track,
     &              Refl_1, Unc_1, Buf_Sa_250, Vflag_1, Data_size,
     &              Error_Flag)

      if (error_flag) CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,
     &   'Read_L1B failed reading band 1','GetModisDat_MOD04')
C
C Identify invalid value (Vflag=-1) in band 1 by replacing the L1B data
C data with standard fill value (FV_L1B = -999.0)
C
      Do  J_flag = 1, Data_size(2)
        Do  I_flag= 1, Data_Size(1)
         If (Vflag_1(I_flag,J_flag).lt.0) Then
             Refl_1(I_flag,J_flag) = FV_L1B
         Endif
        Enddo
      Enddo

C -----------------
C Read MODIS band 2
C -----------------
      Band_No=2
      Resol=Resol_Is_250m
      Cal_Type=Cal_Is_Refl
      Gain= ' '
      Buf_Xtrack = 4*Buf_Size1
      Buf_Track = 4*Buf_Size2

      Call Read_L1B(Modfil_L1B_250, FN_L1B_250, ScanCube_No, Band_No,
     &              Gain, Resol, Cal_type, Buf_Xtrack, Buf_Track,
     &              Refl_2, Unc_2, Buf_Sa_250, Vflag_2, Data_size,
     &              error_flag)

      if (error_flag) CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,
     &   'Read_L1B failed reading band 2','GetModisDat_MOD04')
C
C Identify invalid value (Vflag=-1) in band 2 by replacing the L1B data
C data with standard fill value (FV_L1B = -999.0)
C
      Do  J_flag = 1, Data_size(2)
        Do  I_flag= 1, Data_Size(1)
         If (Vflag_2(I_flag,J_flag).lt.0) Then
             Refl_2(I_flag,J_flag) = FV_L1B
         Endif
        Enddo
      Enddo

C -----------------
C Read MODIS band 3
C -----------------
      Band_No=3
      Resol=Resol_Is_500m
      Cal_Type=Cal_Is_Refl
      Gain= ' '
      Buf_Xtrack = 2*Buf_Size1
      Buf_Track = 2*Buf_Size2

      Call Read_L1B(Modfil_L1B_500, FN_L1B_500, ScanCube_No, Band_No,
     &              Gain, Resol, Cal_type, Buf_Xtrack, Buf_Track,
     &              Refl_3, Unc_3, Buf_Sa_500, Vflag_3, Data_size,
     &              error_flag)

      if (error_flag) CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,
     &   'Read_L1B failed reading band 3','GetModisDat_MOD04')
C
C Identify invalid value (Vflag=-1) in band 3 by replacing the L1B data
C data with standard fill value (FV_L1B = -999.0)
C
      Do  J_flag = 1, Data_size(2)
        Do  I_flag= 1, Data_Size(1)
         If (Vflag_3(I_flag,J_flag).lt.0) Then
           numnum3=numnum3+1
             Refl_3(I_flag,J_flag) = FV_L1B
         Endif
        Enddo
      Enddo

C -----------------
C Read MODIS band 4
C -----------------
      Band_No=4
      Resol=Resol_Is_500m
      Cal_Type=Cal_Is_Refl
      Gain= ' '
      Buf_Xtrack = 2*Buf_Size1
      Buf_Track = 2*Buf_Size2

      Call Read_L1B(Modfil_L1B_500, FN_L1B_500, ScanCube_No, Band_No,
     &              Gain, Resol, Cal_type, Buf_Xtrack, Buf_Track,
     &              Refl_4, Unc_4, Buf_Sa_500, Vflag_4, Data_size,
     &              error_flag)

      if (error_flag) CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,
     &   'Read_L1B failed reading band 4','GetModisDat_MOD04')
C
C Identify invalid value (Vflag=-1) in band 4 by replacing the L1B data
C data with standard fill value (FV_L1B = -999.0)
C
      Do  J_flag = 1, Data_size(2)
        Do  I_flag= 1, Data_Size(1)
         If (Vflag_4(I_flag,J_flag).lt.0) Then
             Refl_4(I_flag,J_flag) = FV_L1B
         Endif
        Enddo
      Enddo
c              if( iscan .eq.15)then
c              write(43,*) Data_Size(1),Data_size(2)
c              do J_flag=1, Data_size(2)
c             write(43,*)'in GetModisDat_MOD04',
c     *       J_flag,(Refl_4(I_flag,j_flag),I_flag=1,Data_Size(1))
c              enddo
c              endif

C -----------------
C Read MODIS band 5
C -----------------
      Band_No=5
      Resol=Resol_Is_500m
      Cal_Type=Cal_Is_Refl
      Gain= ' '
      Buf_Xtrack = 2*Buf_Size1
      Buf_Track = 2*Buf_Size2

      Call Read_L1B(Modfil_L1B_500, FN_L1B_500, ScanCube_No, Band_No,
     &              Gain, Resol, Cal_type, Buf_Xtrack, Buf_Track,
     &              Refl_5, Unc_5, Buf_Sa_500, Vflag_5, Data_size,
     &              error_flag)

      if (error_flag) CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,
     &   'Read_L1B failed reading band 5','GetModisDat_MOD04')
C
C Identify invalid value (Vflag=-1) in band 5 by replacing the L1B data
C data with standard fill value (FV_L1B = -999.0)
C
      Do  J_flag = 1, Data_size(2)
        Do  I_flag= 1, Data_Size(1)
         If (Vflag_5(I_flag,J_flag).lt.0) Then
             Refl_5(I_flag,J_flag) = FV_L1B
         Endif
        Enddo
      Enddo

C -----------------
C Read MODIS band 6
C -----------------
      Band_No=6
      Resol=Resol_Is_500m
      Cal_Type=Cal_Is_Refl
      Gain= ' '
      Buf_Xtrack = 2*Buf_Size1
      Buf_Track = 2*Buf_Size2

      Call Read_L1B(Modfil_L1B_500, FN_L1B_500, ScanCube_No, Band_No,
     &              Gain, Resol, Cal_type, Buf_Xtrack, Buf_Track,
     &              Refl_6, Unc_6, Buf_Sa_500, Vflag_6, Data_size,
     &              error_flag)

      if (error_flag) CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,
     &   'Read_L1B failed reading band 6','GetModisDat_MOD04')
C
C Identify invalid value (Vflag=-1) in band 6 by replacing the L1B data
C data with standard fill value (FV_L1B = -999.0)
C
      Do  J_flag = 1, Data_size(2)
        Do  I_flag= 1, Data_Size(1)
         If (Vflag_6(I_flag,J_flag).lt.0) Then
            numnum1=numnum1+1
             Refl_6(I_flag,J_flag) = FV_L1B
         Endif
        Enddo
      Enddo

C -----------------
C Read MODIS band 7
C -----------------
      Band_No=7
      Resol=Resol_Is_500m
      Cal_Type=Cal_Is_Refl
      Gain= ' '
      Buf_Xtrack = 2*Buf_Size1
      Buf_Track = 2*Buf_Size2

      Call Read_L1B(Modfil_L1B_500, FN_L1B_500, ScanCube_No, Band_No,
     &              Gain, Resol, Cal_type, Buf_Xtrack, Buf_Track,
     &              Refl_7, Unc_7, Buf_Sa_500, Vflag_7, Data_size,
     &              error_flag)

      if (error_flag) CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,
     &   'Read_L1B failed reading band 7','GetModisDat_MOD04')
C
C Identify invalid value (Vflag=-1) in band 7 by replacing the L1B data
C data with standard fill value (FV_L1B = -999.0)
C

      Do  J_flag = 1, Data_size(2)
        Do  I_flag= 1, Data_Size(1)
         If (Vflag_7(I_flag,J_flag).lt.0) Then
             Refl_7(I_flag,J_flag) = FV_L1B
         Endif
        Enddo
      Enddo
C C -----------------
C Read MODIS band 8
C -----------------
      Band_No=8
      Resol=Resol_Is_1km
      Cal_Type=Cal_Is_Refl
      Gain= ' '
      Buf_Xtrack = Buf_Size1
      Buf_Track = Buf_Size2

      Call Read_L1B(Modfil_L1B_1km, FN_L1B_1km, ScanCube_No, Band_No,
     &              Gain, Resol, Cal_type, Buf_Xtrack, Buf_Track,
     &              Refl_8, Unc_8, Buf_Sa_1km, Vflag_8, Data_size,
     &              error_flag)

      if (error_flag) CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,
     &   'Read_L1B failed reading band 8','GetModisDat_MOD04')
C
C Identify invalid value (Vflag=-1) in band 9 by replacing the L1B data
C data with standard fill value (FV_L1B = -999.0)
C
      Do  J_flag = 1, Data_size(2)
        Do  I_flag= 1, Data_Size(1)
         If (Vflag_8(I_flag,J_flag).lt.0) Then
             Refl_8(I_flag,J_flag) = FV_L1B
         Endif
        Enddo
      Enddo


C -----------------
C Read MODIS band 9
C -----------------
      Band_No=9
      Resol=Resol_Is_1km
      Cal_Type=Cal_Is_Refl
      Gain= ' '
      Buf_Xtrack = Buf_Size1
      Buf_Track = Buf_Size2

      Call Read_L1B(Modfil_L1B_1km, FN_L1B_1km, ScanCube_No, Band_No,
     &              Gain, Resol, Cal_type, Buf_Xtrack, Buf_Track,
     &              Refl_9, Unc_9, Buf_Sa_1km, Vflag_9, Data_size,
     &              error_flag)

      if (error_flag) CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,
     &   'Read_L1B failed reading band 9','GetModisDat_MOD04')
C
C Identify invalid value (Vflag=-1) in band 9 by replacing the L1B data
C data with standard fill value (FV_L1B = -999.0)
C
      Do  J_flag = 1, Data_size(2)
        Do  I_flag= 1, Data_Size(1)
         If (Vflag_9(I_flag,J_flag).lt.0) Then
             Refl_9(I_flag,J_flag) = FV_L1B
         Endif
        Enddo
      Enddo

C ------------------
C Read MODIS band 12
C ------------------
      Band_No=12
      Resol=Resol_Is_1km
      Cal_Type=Cal_Is_Refl
      Gain= ' '
      Buf_Xtrack = Buf_Size1
      Buf_Track = Buf_Size2

      Call Read_L1B(Modfil_L1B_1km, FN_L1B_1km, ScanCube_No, Band_No,
     &              Gain, Resol, Cal_type, Buf_Xtrack, Buf_Track,
     &              Refl_12, Unc_12, Buf_Sa_1km, Vflag_12, Data_size,
     &              error_flag)

      if (error_flag) CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,
     &   'Read_L1B failed reading band 12','GetModisDat_MOD04')
C
C Identify invalid value (Vflag=-1) in band 12 by replacing the L1B data
C data with standard fill value (FV_L1B = -999.0)
C
      Do  J_flag = 1, Data_size(2)
        Do  I_flag= 1, Data_Size(1)
         If (Vflag_12(I_flag,J_flag).lt.0) Then
             Refl_12(I_flag,J_flag) = FV_L1B
         Endif
        Enddo
      Enddo

C ---------------------------------------------------------------------
C Read MODIS band 13 with low gain (saturation at high radiation level)
C ---------------------------------------------------------------------
      Band_No=13
      Resol=Resol_Is_1km
      Cal_Type=Cal_Is_Refl
      Gain= 'L'
      Buf_Xtrack = Buf_Size1
      Buf_Track = Buf_Size2

      Call Read_L1B(Modfil_L1B_1km, FN_L1B_1km, ScanCube_No, Band_No,
     &              Gain, Resol, Cal_type, Buf_Xtrack, Buf_Track,
     &              Refl_13, Unc_13, Buf_Sa_1km, Vflag_13, Data_size,
     &              error_flag)

      if (error_flag) CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,
     &   'Read_L1B failed reading band 13','GetModisDat_MOD04')
C
C Identify invalid value (Vflag=-1) in band 13 by replacing the L1B data
C data with standard fill value (FV_L1B = -999.0)
C
      Do  J_flag = 1, Data_size(2)
        Do  I_flag= 1, Data_Size(1)
         If (Vflag_13(I_flag,J_flag).lt.0) Then
             Refl_13(I_flag,J_flag) = FV_L1B
         Endif
        Enddo
      Enddo
      
      
CC ------------------
C Read MODIS band 15
C ------------------
      Band_No=15
      Resol=Resol_Is_1km
      Cal_Type=Cal_Is_Refl
      Gain= ' '
      Buf_Xtrack = Buf_Size1
      Buf_Track = Buf_Size2

      Call Read_L1B(Modfil_L1B_1km, FN_L1B_1km, ScanCube_No, Band_No,
     &              Gain, Resol, Cal_type, Buf_Xtrack, Buf_Track,
     &              Refl_15, Unc_15, Buf_Sa_1km, Vflag_15, Data_size,
     &              error_flag)

      if (error_flag) CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,
     &   'Read_L1B failed reading band 16','GetModisDat_MOD04')
C
C Identify invalid value (Vflag=-1) in band 16 by replacing the L1B data
C data with standard fill value (FV_L1B = -999.0)
C
      Do  J_flag = 1, Data_size(2)
        Do  I_flag= 1, Data_Size(1)
         If (Vflag_15(I_flag,J_flag).lt.0) Then
             Refl_15(I_flag,J_flag) = FV_L1B
         Endif
        Enddo
      Enddo

 
C ------------------
C Read MODIS band 16
C ------------------
      Band_No=16
      Resol=Resol_Is_1km
      Cal_Type=Cal_Is_Refl
      Gain= ' '
      Buf_Xtrack = Buf_Size1
      Buf_Track = Buf_Size2

      Call Read_L1B(Modfil_L1B_1km, FN_L1B_1km, ScanCube_No, Band_No,
     &              Gain, Resol, Cal_type, Buf_Xtrack, Buf_Track,
     &              Refl_16, Unc_16, Buf_Sa_1km, Vflag_16, Data_size,
     &              error_flag)

      if (error_flag) CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,
     &   'Read_L1B failed reading band 16','GetModisDat_MOD04')
C
C Identify invalid value (Vflag=-1) in band 16 by replacing the L1B data
C data with standard fill value (FV_L1B = -999.0)
C
      Do  J_flag = 1, Data_size(2)
        Do  I_flag= 1, Data_Size(1)
         If (Vflag_16(I_flag,J_flag).lt.0) Then
             Refl_16(I_flag,J_flag) = FV_L1B
         Endif
        Enddo
      Enddo

C ------------------
C Read MODIS band 26
C ------------------
      Band_No=26
      Resol=Resol_Is_1km
      Cal_Type=Cal_Is_Refl
      Gain= ' '
      Buf_Xtrack = Buf_Size1
      Buf_Track = Buf_Size2

      Call Read_L1B(Modfil_L1B_1km, FN_L1B_1km, ScanCube_No, Band_No,
     &              Gain, Resol, Cal_type, Buf_Xtrack, Buf_Track,
     &              Refl_26, Unc_26, Buf_Sa_1km, Vflag_26, Data_size,
     &              error_flag)

      if (error_flag) CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,
     &   'Read_L1B failed reading band 26','GetModisDat_MOD04')
C
C Identify invalid value (Vflag=-1) in band 26 by replacing the L1B data
C data with standard fill value (FV_L1B = -999.0)
C
      Do  J_flag = 1, Data_size(2)
        Do  I_flag= 1, Data_Size(1)
         If (Vflag_26(I_flag,J_flag).lt.0) Then
             Refl_26(I_flag,J_flag) = FV_L1B
         Endif
        Enddo
      Enddo

C ------------------
C Read MODIS band 20
C ------------------
      Band_No=20
      Resol=Resol_Is_1km
      Cal_Type=Cal_Is_Rad
      Gain= ' '
      Buf_Xtrack = Buf_Size1
      Buf_Track = Buf_Size2

      Call Read_L1B(Modfil_L1B_1km, FN_L1B_1km, ScanCube_No, Band_No,
     &              Gain, Resol, Cal_type, Buf_Xtrack, Buf_Track,
     &              Rad_20, Unc_20, Buf_Sa_1km, Vflag_20, Data_size,
     &              error_flag)

      if (error_flag) CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,
     &   'Read_L1B failed reading band 20','GetModisDat_MOD04')
C
C Identify invalid value (Vflag=-1) in band 20 by replacing the L1B data
C data with standard fill value (FV_L1B = -999.0)
C
      Do  J_flag = 1, Data_size(2)
        Do  I_flag= 1, Data_Size(1)
         If (Vflag_20(I_flag,J_flag).lt.0) Then
             Rad_20(I_flag,J_flag) = FV_L1B
         Endif
        Enddo
      Enddo

C ------------------
C Read MODIS band 31
C ------------------
      Band_No=31
      Resol=Resol_Is_1km
      Cal_Type=Cal_Is_Rad
      Gain= ' '
      Buf_Xtrack = Buf_Size1
      Buf_Track = Buf_Size2

      Call Read_L1B(Modfil_L1B_1km, FN_L1B_1km, ScanCube_No, Band_No,
     &              Gain, Resol, Cal_type, Buf_Xtrack, Buf_Track,
     &              Rad_31, Unc_31, Buf_Sa_1km, Vflag_31, Data_size,
     &              error_flag)

      if (error_flag) CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,
     &   'Read_L1B failed reading band 31','GetModisDat_MOD04')
C
C Identify invalid value (Vflag=-1) in band 31 by replacing the L1B data
C data with standard fill value (FV_L1B = -999.0)
C
      Do  J_flag = 1, Data_size(2)
        Do  I_flag= 1, Data_Size(1)
         If (Vflag_31(I_flag,J_flag).lt.0) Then
             Rad_31(I_flag,J_flag) = FV_L1B
         Endif
        Enddo
      Enddo

C ------------------
C Read MODIS band 32
C ------------------
      Band_No=32
      Resol=Resol_Is_1km
      Cal_Type=Cal_Is_Rad
      Gain= ' '
      Buf_Xtrack = Buf_Size1
      Buf_Track = Buf_Size2

      Call Read_L1B(Modfil_L1B_1km, FN_L1B_1km, ScanCube_No, Band_No,
     &              Gain, Resol, Cal_type, Buf_Xtrack, Buf_Track,
     &              Rad_32, Unc_32, Buf_Sa_1km, Vflag_32, Data_size,
     &              error_flag)

      If (error_flag) THEN
        CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,
     &     'Read_L1B failed reading band 32','GetModisDat_MOD04')
      Endif

C
C Identify invalid value (Vflag=-1) in band 32 by replacing the L1B data
C data with standard fill value (FV_L1B = -999.0)
C
      Do  J_flag = 1, Data_size(2)
        Do  I_flag= 1, Data_Size(1)
         If (Vflag_32(I_flag,J_flag).lt.0) Then
             Rad_32(I_flag,J_flag) = FV_L1B
         Endif
        Enddo
      Enddo
      NUMSQ=(Data_Size(1)/IGRIDX)*(Data_Size(2)/IGRIDY)
C
C----------------------------------------------------------------------
C Normalize sensor radiance data to reflectance units
C----------------------------------------------------------------------
C
C Loop over 1-km lines and pixels in scan cube
C
      Do 30 il1 = 1, Data_Size(2)
      Do 30 ip1 = 1, Data_Size(1)

         isl4 = (il1-1)*4 + 1
         iel4 = il1*4
         isp4 = (ip1-1)*4 + 1
         iep4 = ip1*4

         isl2 = (il1-1)*2 + 1
         iel2 = il1*2
         isp2 = (ip1-1)*2 + 1
         iep2 = ip1*2

         cossza = cos(dtr*SolZen(ip1,il1))

         If (abs(cossza) .LT. Zero_Eps) Then
            x_cossza = 0.0
         Else
            x_cossza = 1.0/cossza
         End If
C
C Normalize 1-km bands to reflectance units
C
         if ( abs((SolZen(ip1,il1)-FV_GEO)/FV_GEO)
     &       .LT. Rel_Equality_Eps ) then
            Refl_8(ip1,il1) = FV_L1B
            Refl_9(ip1,il1) = FV_L1B
            Refl_12(ip1,il1) = FV_L1B
            Refl_13(ip1,il1) = FV_L1B
            Refl_16(ip1,il1) = FV_L1B
            Refl_26(ip1,il1) = FV_L1B
            Refl_15(ip1,il1) = FV_L1B

         else if (abs(cossza) .LT. Zero_Eps) then
            Refl_8(ip1,il1) = FV_L1B
            Refl_9(ip1,il1) = FV_L1B
            Refl_12(ip1,il1) = FV_L1B
            Refl_13(ip1,il1) = FV_L1B
            Refl_16(ip1,il1) = FV_L1B
            Refl_26(ip1,il1) = FV_L1B
            Refl_15(ip1,il1) = FV_L1B

         else if (SolZen(ip1,il1) .LT. SolZen_Threshold) then
         
               if (abs((Refl_8(ip1,il1)-FV_L1B)/FV_L1B)
     &         .GT. Rel_Equality_Eps) then
               Refl_8(ip1,il1) = x_cossza*Refl_8(ip1,il1)
            end if

            if (abs((Refl_9(ip1,il1)-FV_L1B)/FV_L1B)
     &         .GT. Rel_Equality_Eps) then
               Refl_9(ip1,il1) = x_cossza*Refl_9(ip1,il1)
            end if
            if (abs((Refl_12(ip1,il1)-FV_L1B)/FV_L1B)
     &         .GT. Rel_Equality_Eps) then
               Refl_12(ip1,il1) = x_cossza*Refl_12(ip1,il1)
            end if
            if (abs((Refl_13(ip1,il1)-FV_L1B)/FV_L1B)
     &         .GT. Rel_Equality_Eps) then
               Refl_13(ip1,il1) = x_cossza*Refl_13(ip1,il1)
            end if
            if (abs((Refl_16(ip1,il1)-FV_L1B)/FV_L1B)
     &         .GT. Rel_Equality_Eps) then
               Refl_16(ip1,il1) = x_cossza*Refl_16(ip1,il1)
            end if
             if (abs((Refl_15(ip1,il1)-FV_L1B)/FV_L1B)
     &         .GT. Rel_Equality_Eps) then
               Refl_15(ip1,il1) = x_cossza*Refl_15(ip1,il1)
            end if
            if (abs((Refl_26(ip1,il1)-FV_L1B)/FV_L1B)
     &         .GT. Rel_Equality_Eps) then
              Refl_26(ip1,il1) = x_cossza*Refl_26(ip1,il1) 
            end if

         else
            Refl_8(ip1,il1) = FV_L1B
            Refl_9(ip1,il1) = FV_L1B
            Refl_12(ip1,il1) = FV_L1B
            Refl_13(ip1,il1) = FV_L1B
            Refl_16(ip1,il1) = FV_L1B
            Refl_15(ip1,il1) = FV_L1B
            Refl_26(ip1,il1) = FV_L1B
         endif
C
C Loop over 500-m lines and pixels between 1-km footprints.
C
         Do 40 il2 = isl2, iel2
         Do 40 ip2 = isp2, iep2

            if ( abs((SolZen(ip1,il1)-FV_GEO)/FV_GEO)
     &          .LT. Rel_Equality_Eps ) then
               Refl_3(ip2,il2) = FV_L1B
               Refl_4(ip2,il2) = FV_L1B
               Refl_5(ip2,il2) = FV_L1B
               Refl_6(ip2,il2) = FV_L1B
               Refl_7(ip2,il2) = FV_L1B

            else if (abs(cossza) .LT. Zero_Eps) then
               Refl_3(ip2,il2) = FV_L1B
               Refl_4(ip2,il2) = FV_L1B
               Refl_5(ip2,il2) = FV_L1B
               Refl_6(ip2,il2) = FV_L1B
               Refl_7(ip2,il2) = FV_L1B

            else if (SolZen(ip1,il1) .LT. SolZen_Threshold) then

               if (abs((Refl_3(ip2,il2)-FV_L1B)/FV_L1B)
     &            .GT. Rel_Equality_Eps) then
                  Refl_3(ip2,il2) = x_cossza*Refl_3(ip2,il2)
               end if
               if (abs((Refl_4(ip2,il2)-FV_L1B)/FV_L1B)
     &            .GT. Rel_Equality_Eps) then
                  Refl_4(ip2,il2) = x_cossza*Refl_4(ip2,il2)
               end if
               if (abs((Refl_5(ip2,il2)-FV_L1B)/FV_L1B)
     &            .GT. Rel_Equality_Eps) then
                  Refl_5(ip2,il2) = x_cossza*Refl_5(ip2,il2)
               end if
               if (abs((Refl_6(ip2,il2)-FV_L1B)/FV_L1B)
     &            .GT. Rel_Equality_Eps) then
                  Refl_6(ip2,il2) = x_cossza*Refl_6(ip2,il2)
               end if
               if (abs((Refl_7(ip2,il2)-FV_L1B)/FV_L1B)
     &            .GT. Rel_Equality_Eps) then
                  Refl_7(ip2,il2) = x_cossza*Refl_7(ip2,il2)
               end if
            else
               Refl_3(ip2,il2) = FV_L1B
               Refl_4(ip2,il2) = FV_L1B
               Refl_5(ip2,il2) = FV_L1B
               Refl_6(ip2,il2) = FV_L1B
               Refl_7(ip2,il2) = FV_L1B
            endif
   40    continue
C
C Loop over 250-m lines and pixels between 1-km footprints.
C
         Do 50 il4 = isl4, iel4
         Do 50 ip4 = isp4, iep4
C
C            if ((ip1.eq.1 .and. il1.eq.1) .and. (idebug .eq. 1)) then
C               print *, 'GetModisData before normalization'
C               print *, 'refl_1,refl_2',ip4,il4,refl_1(ip4,il4),
C     2                   refl_2(ip4,il4)
C            end if
C
            if ( abs((SolZen(ip1,il1)-FV_GEO)/FV_GEO)
     &          .LT. Rel_Equality_Eps ) then
               Refl_1(ip4,il4) = FV_L1B
               Refl_2(ip4,il4) = FV_L1B

            else if (abs(cossza) .LT. Zero_Eps) then
               Refl_1(ip4,il4) = FV_L1B
               Refl_2(ip4,il4) = FV_L1B

            else if (SolZen(ip1,il1) .LT. SolZen_Threshold) then

               if (abs((Refl_1(ip4,il4)-FV_L1B)/FV_L1B)
     &            .GT. Rel_Equality_Eps) then
                  Refl_1(ip4,il4) = x_cossza*Refl_1(ip4,il4)
               end if

               if (abs((Refl_2(ip4,il4)-FV_L1B)/FV_L1B)
     &            .GT. Rel_Equality_Eps) then
                  Refl_2(ip4,il4) = x_cossza*Refl_2(ip4,il4)
               end if

            else
               Refl_1(ip4,il4) = FV_L1B
               Refl_2(ip4,il4) = FV_L1B
            end if

   50    continue

   30 continue
 
        CALL Read_CldMsk(Modfil_CldMsk,ScanCube_No,Buf_cldmsk,
     1                 Buf_cldmsk_QA,Buf_Size1,Buf_Size2,
     2                 DS_Dim1_CM,DS_Dim1_QA,DS_Dim2,DS_Dim3,
     3                 Cloud,QA_Cloud,Error_Flag)

C
C Decode cloud mask or generate our own cloud mask from L1B data
C
      CALL CldMsk_Info_MOD04(Modfil_CldMsk,Buf_cldmsk,
     1    Buf_cldmsk_QA,Buf_Size1,Buf_Size2,
     2    Cloud,QA_Cloud,CldMsk_250,CldMsk_500,CldMsk_1km,DET_Flag,
     3    UFQ_Flag,DayNight_Flag,SunGlint_Flag,SnowIce_Flag,
     4   Non_CloudOb_Flag,Thin_CirNIR_Flag,Shadow_Flag,
c     4    LandSea_Flag,Non_CloudOb_Flag,Thin_CirNIR_Flag,Shadow_Flag,
     5    Thin_CirIR_Flag,Cloud_SimpIR_Flag,High_Cloud_Flag,
     6    Cloud_IRTemp_Flag,Cloud_3p75_11_Flag,Cloud_VisRat_Flag,
     7    Cloud_SpatVar_Flag)
C
C
      GetModisDat_MOD04 = 0
      
      
      
  
       

      return
      end



c***********************************************************************
      INTEGER FUNCTION GetGeo_MOD04(modfil,ScanCube_No,Buf_Size1,
     &        Buf_Size2,DataSize,Lat,Lon,Sol_Zen,Sol_Az,Sat_Zen,
     &        Sat_Az,Rel_Az,Height,LandSea_Flag)

      IMPLICIT NONE
      SAVE
      INCLUDE 'mapi.inc'
      INCLUDE 'hdf.inc'
      INCLUDE 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C!F77
C
C !DESCRIPTION: Retrieves MODIS Geolocation data from an HDF target
C              array of typically 100 scan cubes (a granule).  This
C              routine is customized for MOD04 land algorithm.
C              In the definitions below x=BufSize1; y=BufSize2
C
C !INPUT PARAMETERS:
C
C    INTEGER modfil(MODFILLEN)
C                          SD_ID, File_ID and File access type.
C    INTEGER ScanCube_No   The scancube number to be read.
C    INTEGER Buf_Size1/2   The size of the memory buffer to hold the
C                          geolocation data
C
C !OUTPUT PARAMETERS:
C
C    INTEGER DataSize(2)  Actual size of geolocation data
C    REAL  Lat(x,y)       Buffer for Latitude data.
C    REAL  Lon(x,y)       Buffer for Longitude data.
C    REAL  Sol_Zen(x,y)   Buffer for solar zenith data.
C    REAL  Sat_Zen(x,y)   Buffer for sensor zenith data.
C    REAL  Rel_Az(x,y)    Buffer to hold relative azimuth data.
C    REAL  Height(x,y)    Buffer for Height (m)
C
c!REVISION HISTORY:
c 01/29/98 fhliang
c added NCSA acknowledgement;
c fixed prolog.
c
c!TEAM-UNIQUE HEADER:
c
c Developed by MODIS Aerosol/Solar Water Vapor Retrieval Team
c GSFC, Greenbelt, MD
c
c HDF portions developed at the National Center for Supercomputing
c Applications at the University of Illinois at Urbana-Champaign.
c
C DESIGN NOTES:
C
C    This subroutine checks the return status of all MODIS Application
C    Program Interface (M-API) function calls.  A successful M-API
C    function call is indicated by a return value 0. If a M-API
C    function call is not successful, an error message is written
C    to the LogStatus file and the process aborts. This is achieved
C    by passing an fatal mnemonic error message (i.e. .._F_..)
C    to the function MODIS_SMF_SETDYNAMICMSG.
C
C Externals:
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
C      MODIS_F_GENERIC            (MODIS_39500.f)
C
C Internals:
C
C    Local Variables:
C       arrnam      Name of the SDS array.
C       grpnm       Name of the data group
C       dtype       The data type of the array.
C       msgbuf*     Temporary  Message buffer.
C       ret_dim(2)  Buffer holding the size of a 2-dimensional HDF
C                   SDS array.
C       Start(2)    The start position of data.
C       Fmax        Maximum frame number per scan line.
C       Lmax        Maximum line number per scan cube.
C       Rank        The number of dimensions in an array
C       I,j,k,L     Temporary variables for holding integer values.
C       TSCN        Total ScanCube Number.
C       count       Temporary buffer for data of the target array.
C       NofLine_PerCube The number of lines per scan cube
C
C    Return Values:
C       GetGeo_MOD04     0 for success / -1 if failed
C
C !END
C-----------------------------------------------------------------------

C Declarations
      CHARACTER*13 grpnm,dtype
      CHARACTER*30 arrnm
      CHARACTER*80 msgbuf,msgbuf1,msgbuf2
      INTEGER   I,J,K,L,ret_dim(2),NofLine_PerCube,Rank,Fmax,Lmax,TSCN,
     &          RTN,ScanCube_No,Buf_Size1,Buf_Size2,FILL_INT,Start(2),
     &          DataSize(2),modfil(MODFILLEN),
     &          LandSea_Flag(Buf_Size1,Buf_Size2)
      PARAMETER (Fmax=1500,Lmax=10,NofLine_PerCube=10,FILL_INT=-1000)
      INTEGER*2 count(Fmax*Lmax)
       byte     count_Byte(Fmax*Lmax)  
      REAL pi,FILL_FPT,Rel_equality_EPS,r4count(Fmax*Lmax),
     &     Lon(Buf_Size1,Buf_Size2),Lat(Buf_Size1,Buf_Size2),
     &     Sol_Zen(Buf_Size1,Buf_Size2),Sat_Zen(Buf_Size1,Buf_Size2),
     &     Rel_Az(Buf_Size1,Buf_Size2),Sol_Az(Buf_Size1,Buf_Size2),
     &     Sat_Az(Buf_Size1,Buf_Size2),Height(Buf_Size1,Buf_Size2)
      PARAMETER (FILL_FPT=-1000.0,Rel_equality_EPS=0.000001,
     &           pi = 3.14159265)

C Initialization

      GetGeo_MOD04=-1
      grpnm = ' '
      Start(1)=0
      Start(2)=(ScanCube_No-1)*NofLine_PerCube


      Do J = 1, Buf_Size2
      Do I = 1, Buf_Size1
         Lon(I,J) = 0.0
         Lat(I,J) = 0.0
         Sol_Zen(I,J) = 0.0
         Sat_Zen(I,J) = 0.0
         Sol_Az(I,J) = 0.0
         Sat_Az(I,J) = 0.0
         Rel_Az(I,J) = 0.0
         Height(I,J) = 0.0
         LandSea_Flag(I,J) =0
      End Do
      End Do

C Checking for valid input
      If (modfil(1).le.0.or.modfil(3).ne.DFACC_READ) CALL
     &    MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,
     &   'Invalid SD ID/file access type','GetGeo_MOD04')

       Do 100 I = 1, 8
          If (I.eq.1) Then
            arrnm='Longitude'
         Else If (I.eq.2) Then
            arrnm='Latitude'
         Else If (I.eq.3) Then
            arrnm='SolarZenith'
         Else If (I.eq.4) Then
            arrnm='SensorZenith'
         Else If (I.eq.5) Then
            arrnm='SolarAzimuth'
         Else if (I.eq.6) Then
            arrnm='SensorAzimuth'
         Else if (I.eq.7) Then
         arrnm='Height' 
           Else if (I.eq.8) Then
           arrnm='Land/SeaMask'  
          End If
          
         

C Retrieving the rank, dimensions and data type of SDS data.
         Rank = 2

         If (GMARDM(modfil,arrnm,grpnm,dtype,Rank,ret_dim).ne.MAPIOK)
     &   CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,'GMARDM failed',
     &   'GetGeo_MOD04')
         DataSize(1) = ret_dim(1)
         DataSize(2)=NofLine_PerCube


C Additional input check of ScanCube_No.
         TSCN=ret_dim(2)/NofLine_PerCube

         If (ScanCube_No.lt.1.or.ScanCube_No.gt.TSCN) Then
            write(msgbuf,'(i4)') TSCN
            call Concatenate('ScanCube_No out of bounds : 1 - ',msgbuf,
     &      msgbuf1)
            CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,msgbuf1,
     &      'GetGeo_MOD04')
         End If

C Read HDF target array into 'count' or 'r4count' buffer
         If (I.eq.1 .or. I.eq.2) Then
            RTN = GMAR(modfil,arrnm,grpnm,Start,DataSize,r4count)
        endif
        If (I.gt.2 .or. I.lt.8) Then
            RTN = GMAR(modfil,arrnm,grpnm,Start,DataSize,count)
         end if
         If (I.eq. 8) Then
           RTN = GMAR(modfil,arrnm,grpnm,Start,DataSize,count_Byte)
         end if
         If (RTN.ne.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_F_GENERIC,'GMAR failed','GetGeo_MOD04')
         DataSize(1) = ret_dim(1)
         DataSize(2)=NofLine_PerCube


C Convert HDF 2-byte integer data to 4-byte real value.

C Retrieve Longitude
         If (I.eq.1) Then
            L = 0

            Do 10 K = 1, DataSize(2)
            Do 15 J = 1, DataSize(1)
               L = L + 1
               if (abs((r4count(L)-FILL_FPT)/r4count(L)).lt.
     &             Rel_equality_EPS) then
                   Lon(J,K) = FILL_FPT
               else
                   Lon(J,K) = r4count(L)
               end if
   15       continue
   10       continue

C Retrieve latitude
         Else If (I.eq.2) Then
            L = 0

            Do 20 K = 1, DataSize(2)
            Do 25 J = 1, DataSize(1)
               L = L + 1
               if (abs((r4count(L)-FILL_FPT)/r4count(L)).lt.
     &            Rel_equality_EPS) then
                  Lat(J,K) = FILL_FPT
               else
                  Lat(J,K) = r4count(L)
               end if
   25       continue
   20       continue

C Retrieve Solar Zenith Angle
         Else If (I.eq.3) Then
            L = 0

            Do 30 K = 1, DataSize(2)
            Do 35 J = 1, DataSize(1)
               L = L + 1
               if (count(L).eq.FILL_INT) then
                  Sol_Zen(J,K) = FILL_INT
               else
C                  Sol_Zen(J,K) = count(L)*0.018/pi
                  Sol_Zen(J,K) = count(L)*0.01
               end if
   35       continue
   30       continue

C Retrieve Satellite Zenith Angle
         Else If (I.eq.4) Then
            L = 0

            Do 40 K = 1, DataSize(2)
            Do 45 J = 1, DataSize(1)
               L = L + 1
               if (count(L).eq.FILL_INT) then
                  Sat_Zen(J,K) = FILL_INT
               else
C                  Sat_Zen(J,K) = count(L)*0.018/pi
                  Sat_Zen(J,K) = count(L)*0.01
               end if
   45       continue
   40       continue

C Retrieve Solar Azimuth
         Else If (I.eq.5) Then
            L = 0

            Do 50 K = 1, DataSize(2)
            Do 55 J = 1, DataSize(1)
               L = L + 1
               if (count(L).eq.FILL_INT) then
                  Sol_Az(J,K) = FILL_INT
               else
C                  Sol_Az(J,K)=180.0 + count(L)*0.018/pi
C                  Sol_Az(J,K)=180.0 + count(L)*0.01
                  Sol_Az(J,K)=count(L)*0.01
               end if
   55       continue
   50       continue

C Retrieve Satellite Azimuth and compute Relative Azimuth
         Else If (I.eq.6) Then
            L = 0

            Do 60 K = 1, DataSize(2)
            Do 65 J = 1, DataSize(1)
               L = L + 1
               if (count(L).eq.FILL_INT) then 
                  Sol_Az(J,K) = FILL_INT
                  Rel_Az(J,K) = FILL_INT
               else
                   Sat_Az(J,K) = count(L)*0.01 
                   Rel_Az(J,K) = abs(Sat_Az(J,K) - Sol_Az(J,K))
                  Rel_Az(J,K) = amod(Sat_Az(J,K),360.0)
                  If (Rel_Az(J,K).gt.180.0)
     &            Rel_Az(J,K) = 360.0 - Rel_Az(J,K)
               end if
   65       continue
   60       continue

C Retrieve Height above ellipsoid
             Else If (I.eq.7) Then
               L = 0 
               Do 70 K = 1, DataSize(2)
               Do 75 J = 1, DataSize(1)
                  L = L + 1
                  Height(J,K) = FILL_INT
                  if (count(L) .ne. FILL_INT) Height(J,K) = count(L)
   75          continue
   70          continue
        
             Else If (I.eq.8) Then
                L = 0 
               Do 160 K = 1, DataSize(2)
               Do 165 J = 1, DataSize(1) 
                   L = L + 1  
              LandSea_Flag(J,K)= count_Byte(L) 
               
  165       continue
  160       continue
                EndIf
  100        Continue

          GetGeo_MOD04=0

      Return
      END



C*********************************************************************
       INTEGER FUNCTION ReadCldMsk_MOD04(modfil,ScanCube_No,
     2                  Buf_Size1,Buf_Size2,Data_Size,
     3                  Count,LandSea_Flag)

       IMPLICIT NONE
       INCLUDE 'mapi.inc'
       INCLUDE 'hdf.inc'
       INCLUDE 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C
C    Retrieves one scan cube of MODIS Cloud Mask data from
C    an HDF target array of typically 100 scan cubes (a granule).
C    In definitions below, x = Buf_Size1 and y = Buf_Size2
C
C !INPUT PARAMETERS:
C
C    INTEGER    modfil(MODFILLEN)
C                            Array containing SD_ID, File_ID and
C                            File access type, etc.
C    INTEGER    ScanCube_No  The ScanCube Number.
C    INTEGER    Buf_Size1/2  The sizes of 250-m CldMsk buffer.
C
C !OUTPUT PARAMETERS:
C
C    INTEGER    Data_Size(2) Array specifying the size of the data
C                            array placed in CldMsk buffer.
C    INTEGER    CldMsk(x,y)  Buffer storing Cloud Mask (0/cloudy,
C                            1/clear, -1 not determined.)
C    INTEGER    LandSea_FLag(x/4,y/4)  Buffer containing land/sea
C                                      flag: 0 water; 1 coastal;
C                                      2 wetland; 3 land; -1 invalid.
C
c!REVISION HISTORY:
c 01/29/98 fhliang
c added NCSA acknowledgement;
c fixed prolog.
c
c!TEAM-UNIQUE HEADER:
c
c Developed by MODIS Aerosol/Solar Water Vapor Retrieval Team
c GSFC, Greenbelt, MD
c
c HDF portions developed at the National Center for Supercomputing
c Applications at the University of Illinois at Urbana-Champaign.
c
C !DESIGN NOTES:
C
C    This subroutine checks the return status of all MODIS Application
C    Program Interface (M-API) function calls.  A successful M-API
C    function call is indicated by a return value 0. If a M-API
C    function call is not successful, a warning error message is written
C    to the LogStatus file (via function MODIS_SMF_SETDYNAMICMSG), but
C    the process continues and returns with a value of -1.
C
C   Externals:
C
C      Function:
C        GMAR                       (libmapi.a)
C        GMARDM                     (libmapi.a)
C
C      Subroutines:
C        MODIS_SMF_SETDYNAMICMSG
C        CONCATENATE
C
C      Named Constant:
C        DFACC_READ                 (hdf.inc)
C        MAPIOK                     (mapi.inc)
C        MODIS_W_GENERIC            (MODIS_39500.f)
C
C   Internals:
C
C      Variables:
C        arrnam      Name of the SDS array.
C        grpnm       Name of the data group containing the target
C        data_type   String describing the data type of the array.
C        msgbuf      Message buffer.
C        msgbuf1     Message buffer one.
C        msgbuf2     Message buffer two.
C        Dim_Size(2)  Array specifying the size of hdf SDS data array.
C        Edge(2)      Array specifying the number of data value to read.
C        Start(2)     Array specifying the starting location of data.
C        Fmax         Maximum frame number per scan line.
C        Lmax         Maximum line number per scan cube.
C        Rank         The number of dimensions in an array
C        l,l1,l4,     Line and pixel counters p,p1,p4
C        l4_offset    250-m line offset to 1-km cell position
C        p4_offset    250-m pixel offset to 1-km cell position
C        ibit         0-based bit postion (in INTEGER word)
C        idebug       Debug (1)/ No debug (0) flag
C        MaxScanCube_No  Total ScanCube Number.
C        count(15000) A temporary buffer for data of the target array.
C        LinesPerScancube The number of 1-km lines per scan cube
C
C !END
C-----------------------------------------------------------------------

C Declarations
      CHARACTER*80 arrnm,grpnm,data_type,msgbuf,msgbuf1,msgbuf2
      INTEGER      ScanCube_No,LinesPerScancube,Rank,I,j,k,idebug,
     2             MaxScanCube_No,Dim_Size(3)
      INTEGER      l,l1,l4,l4_offset,p,p1,p4,p4_offset,ibit
      INTEGER      No_byte,Fmax,Lmax,index,Ibyte
      PARAMETER    (No_byte=6,Fmax=1500,Lmax=10)
      INTEGER      count(No_byte*Fmax*Lmax)
      INTEGER      Temp,Start(3),Edge(3),Data_Size(2),Buf_Size1,
     2             Buf_Size2,modfil(MODFILLEN)
C     3             CldMsk(Buf_Size1,Buf_Size2),
C     4             LandSea_Flag(Buf_Size1/4,Buf_Size2/4)

      INTEGER CldMsk_250(4*Fmax,4*Lmax),CldMsk_500(2*Fmax,2*Lmax),
     1    CldMsk_1km(Fmax,Lmax),DET_Flag(Fmax,Lmax),UFQ_Flag(Fmax,Lmax),
     2    DayNight_Flag(Fmax,Lmax),SunGlint_Flag(Fmax,Lmax),
     3    SnowIce_Flag(Fmax,Lmax),LandSea_Flag(Fmax,Lmax),
     4    Non_CloudOb_Flag(Fmax,Lmax),Thin_CirNIR_Flag(Fmax,Lmax),
     5    Shadow_Flag(Fmax,Lmax),Thin_CirIR_Flag(Fmax,Lmax),
     6    Cloud_SimpIR_Flag(Fmax,Lmax),High_Cloud_Flag(Fmax,Lmax),
     7    Cloud_IRTemp_Flag(Fmax,Lmax),Cloud_3p75_11_Flag(Fmax,Lmax),
     8    Cloud_VisRat_Flag(Fmax,Lmax),Cloud_SpatVar_Flag(Fmax,Lmax)

      LOGICAL      error_flag

C Initialization
      grpnm = ' '
      arrnm = 'Cloud_Mask'
      ReadCldMsk_MOD04 = -1
      LinesPerScancube = 10
      Rank  =  3
      Start(1) = 0
      Start(2) = 0
      Start(3) = (ScanCube_No-1)*LinesPerScancube
      idebug = 0
      error_flag = .FALSE.

C Check for valid file and band numbers
      IF (modfil(1).le.0.OR.modfil(3).ne.DFACC_READ) THEN
        call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     & 'Invalid SD_ID or invalid file access type','ReadCldMsk_MOD04')
        error_flag = .TRUE.
      END IF


C Retrieve the rank, dimensions and data type of SDS data.
      IF (GMARDM(modfil, arrnm, grpnm, data_type, Rank, Dim_Size)
     &   .ne.MAPIOK) THEN
         call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     &   'GMARDM get dim_info failed','ReadCldMsk_MOD04')
         error_flag = .TRUE.
      END IF
C

C Set size of 250-m output cloud mask.
      Data_Size(1) = 4*Dim_Size(2)
      Data_Size(2) = 4*LinesPerScancube

C  Additional input check of ScanCube_No and buffer size
      MaxScanCube_No = Dim_Size(3)/LinesPerScancube

      IF (ScanCube_No.lt.1 .or. ScanCube_No.gt.MaxScanCube_No) THEN
         write(msgbuf,'(i4)') MaxScanCube_No
         CALL CONCATENATE('ScanCube_No out of bounds; range 1 -',
     &   msgbuf, msgbuf1)
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     &   msgbuf1,'ReadCldMsk_MOD04')
         error_flag = .TRUE.
      END IF

      IF (Buf_Size1.lt.Data_Size(1)) THEN
         CALL MODIS_SMF_SETDYNAMICMSG
     2   (MODIS_W_GENERIC,'Buffer size too small','ReadCldMsk_MOD04')
         error_flag = .TRUE.
      END IF

C Get 1-km Cloud MASK data
      Edge(1) = Dim_Size(1)
      Edge(2) = Dim_Size(2)
      Edge(3) = LinesPerScancube

C Read HDF target array into 'count' buffer
      if (GMAR(modfil, arrnm, grpnm, Start, Edge, count)
     &   .ne.MAPIOK) Then
         call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     &  'GMAR read arrdata failed','ReadCldMsk_MOD04')
         error_flag = .TRUE.
      End If

C-----------------------------------------------------------------------
C Derive cloud mask for 250-m pixels based on the results of the
C visible cloud test.  Also extract Land/Sea flag at 1-km resolution.
C Begin looping over 1-km pixels.
C-----------------------------------------------------------------------

       I = -5

       Do 80 l1 = 1, Edge(3)
       Do 80 p1 = 1, Edge(2)

C The Cloud Mask consists of 6 separate 1-byte words.  Increment
C memory buffer index by 6 to cycle through 1-km frames in a line
          I = I + 6

C         Set 250-m line and pixel offsets to corner of current 1-km
C         cell.  Then, loop over the 16 250-m pixels within the 1-km.
C         0-based bit positions 16-31 of the 1-km INTEGER cloud mask
C         value are associated with 250-m sublines (1-4) and samples
C         (1-4) with line varying most rapidly.

          l4_offset = (l1-1)*4
          p4_offset = (p1-1)*4


C-----------------------------------------------------------------------
C Examine first byte of cloud mask to determine whether cloud mask
C was even determined.  Zero-based bit 0 is 1 for determined, and 0 for
C not determined.  If cloud mask not determined, set 250-m
C CldMsk(p4,l4) to -1.
C
C Also set 1-km LandSea_Flag to -1.
C In Version 1, LandSea_Flag may take 5 values:  0 (water),
C 1 (coastal), 2 (wetland), 3 (land), and -1 (invalid data marker).
C-----------------------------------------------------------------------

          If (ibits(count(I),0,1) .EQ. 0) Then

             LandSea_Flag(p1,l1) = -1

             Do 90 l = 1,4
             Do 90 p = 1,4
                l4 = l4_offset + l
                p4 = p4_offset + p
                CldMsk_250(p4,l4) = -1
   90        continue

C-----------------------------------------------------------------------
C Cloud Mask determined.  Extract value of 1-km LandSea_Flag.
C Also get results (a simple yes/no switch) of 250-m visible cloud test.
C-----------------------------------------------------------------------
          Else
C
C go to bits 6 and 7 (of byte 1) to set land/sea flag:  0 for water;
C 1 coastal; 2 wetland, and 3 land.
C
             DET_Flag(p1,l1) = ibits(count(I),0,1)
             UFQ_Flag(p1,l1) = ibits(count(I),1,2)
             DayNight_Flag(p1,l1) = ibits(count(I),3,1)
             SunGlint_Flag(p1,l1) = ibits(count(I),4,1)
             SnowIce_Flag(p1,l1) = ibits(count(I),5,1)
!             LandSea_Flag(p1,l1) = ibits(count(I),6,2)
             Non_CloudOb_Flag(p1,l1) = ibits(count(I),8,1)
             Thin_CirNIR_Flag(p1,l1) = ibits(count(I),9,1)
             Shadow_Flag(p1,l1) = ibits(count(I),10,1)
             Thin_CirIR_Flag(p1,l1) = ibits(count(I),11,1)
             Cloud_SimpIR_Flag(p1,l1) = ibits(count(I),13,1)
             High_Cloud_Flag(p1,l1) = ibits(count(I),16,1)
             Cloud_IRTemp_Flag(p1,l1) = ibits(count(I),18,1)
             Cloud_3p75_11_Flag(p1,l1) = ibits(count(I),19,1)
             Cloud_VisRat_Flag(p1,l1) = ibits(count(I),21,1)
             Cloud_SpatVar_Flag(p1,l1) = ibits(count(I),25,1)

C             if (l1.eq.1 .and. (p1.eq.1 .or. p1.eq.1250)) then
C            print *, 'special Land/Sea write; I,l1,p1,Land/Sea flag',
C     &                    I,l1, p1, LandSea_Flag(p1,l1)
C             end if

             Do 96 l = 1,4
             Do 96 p = 1,4
                Ibyte = 5
                if (l.gt.2) Ibyte = 6
                l4 = l4_offset + l
                p4 = p4_offset + p
                ibit = mod(l+1,2) * 4 + p -1
                index = I + Ibyte - 1
                CldMsk_250(p4,l4) = ibits(count(index),ibit,1)

C                If (p1.eq.414.and.l1.eq.1.and.idebug.eq.1) THEN
C                   print *
C                   print *, '250-m cloud mask analysis at ',
C     &                      '1-km line 1, pixel 1'
C                   print *, 'p1,l1,p4,l4,ibit,index,CldMsk(p4,l4)'
C                   print *,  p1,l1,p4,l4,ibit,index,CldMsk(p4,l4)
C                End if

   96        continue

          End If

C          If (p1 .eq. 2 .AND. idebug .eq. 1) THEN
C             print *
C             print *, 'analysis at pixel #2 and line ',l1
C
C             Do 98 l = 1,4
C             Do 98 p = 1,4
C                l4 = l4_offset + l
C                p4 = p4_offset + p
C                print *, 'l1,p1,l4,p4,cldMsk(p4,l4)'
C                print *,  l1,p1,l4,p4,cldMsk(p4,l4)
C   98        continue
C          End If

   80  continue

C       If (idebug.eq.1) print *, (p,CldMsk(p,1),p=1652,1654)

       IF (.NOT. error_flag) ReadCldMsk_MOD04 = 0

       Return
       END


************************************************************************
      INTEGER FUNCTION ReadAnc_MOD04(MODFIL_ANC,Buf_Size1,Buf_Size2,
     &        scan,Data_Size,sfc_temp)

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
c!REVISION HISTORY:
c 01/29/98 fhliang
c added NCSA acknowledgement;
c fixed prolog.
c
C !TEAM-UNIQUE HEADER:
C    This software is developed by the MODIS Science Data Support
C    Team for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
c HDF portions developed at the National Center for Supercomputing
c Applications at the University of Illinois at Urbana-Champaign.
c
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

      ReadAnc_MOD04 = -1
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
     2  'Invalid Modfil_ANC','ReadAnc_MOD04')
         error_flag = .true.
      END IF

C Retrieve the rank, dimension and data type of surface temperature
      rtn = GMARDM(MODFIL_ANC,arrnm,grpnm,data_type,rank,Data_Size)
      if (rtn.ne.MAPIOK) THEN
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     2   'GMARDM failed','ReadAnc_MOD04')
         error_flag = .true.
      end if

C Additional input check
      if (Buf_Size1.lt.Data_Size(1)) then
         call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     &   'Buf_Size1 too small','ReadAnc_MOD04')
          error_flag = .true.
      end if

      MaxScan_No = Data_Size(2)/Buf_Size2
      if (scan.lt.1 .OR. scan.gt.MaxScan_No) then
         write(msgbuf1,'(i4)') MaxScan_No
         call Concatenate('Scan_No out of bounds : 1 - ',msgbuf1,
     &   msgbuf)
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,msgbuf,
     &   'ReadAnc_MOD04')
         error_flag = .true.
      end if

C Read HDF target array into 'count' buffer
      count(1) = Data_Size(1)
      count(2) = Buf_Size2
      rtn = GMAR(MODFIL_ANC,arrnm,grpnm,start,count,sfcT)
      if (rtn.ne. MAPIOK) then
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     &   'GMAR failed','ReadAnc_MOD04')
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

C ReadAnc_MOD04 completed successfully !
         ReadAnc_MOD04 = 0

      end if

      return
      end



C***********************************************************************
       SUBROUTINE CldMsk_Info_MOD04(Modfil,No_Byte,No_Byte_QA,Fmax,Lmax,
     1   Count,QA_Count,CldMsk_250,CldMsk_500,CldMsk_1km,DET_Flag,
c     2   UFQ_Flag,DayNight_Flag,SunGlint_Flag,SnowIce_Flag,LandSea_Flag,
     2   UFQ_Flag,DayNight_Flag,SunGlint_Flag,SnowIce_Flag,
     3   Non_CloudOb_Flag,Thin_CirNIR_Flag,Shadow_Flag,Thin_CirIR_Flag,
     4   Cloud_SimpIR_Flag,High_Cloud_Flag,Cloud_IRTemp_Flag,
     5   Cloud_3p75_11_Flag,Cloud_VisRat_Flag,Cloud_SpatVar_Flag)

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
C    No_Byte_QA Number of byte of cloud mask storage (=18)
C    Fmax       Maximum frame number per scan line.
C    Lmax       Maximum line number per scan cube.
C    Count      Cloud mask flags
C    QA_Count   Cloud mask QA flags
C
C !OUTPUT PARAMETERS: (All integers)
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
C
      IMPLICIT NONE

      INCLUDE 'mapi.inc'
      INCLUDE 'PGS_MODIS_39500.f'

      CHARACTER*80 arrnm,grpnm,data_type
      INTEGER  I,j,k
      INTEGER  l,l1,l4,l2,l2_offset,l4_offset,p,p1,p2,p2_offset,p4,
     1         p4_offset,ibit,p11,l11
      INTEGER  No_Byte,No_Byte_QA,Fmax,Lmax,index,Ibyte,Dim_Size(3)
     2         ,Modfil(MODFILLEN),LinesPerScancube,Rank,IR_cloud
C
C      PARAMETER    (No_Byte=6,No_Byte_QA=18,Fmax=1500,Lmax=10)
C
      BYTE    Count(No_byte*Fmax*Lmax),
     2        QA_Count(No_Byte_QA*Fmax*Lmax)

      INTEGER CldMsk_250(4*Fmax,4*Lmax),CldMsk_500(2*Fmax,2*Lmax),
     1    CldMsk_1km(Fmax,Lmax),DET_Flag(Fmax,Lmax),UFQ_Flag(Fmax,Lmax),
     2    DayNight_Flag(Fmax,Lmax),SunGlint_Flag(Fmax,Lmax),
     3    SnowIce_Flag(Fmax,Lmax), 
     4    Non_CloudOb_Flag(Fmax,Lmax),Thin_CirNIR_Flag(Fmax,Lmax),
     5    Shadow_Flag(Fmax,Lmax),Thin_CirIR_Flag(Fmax,Lmax),
     6    Cloud_SimpIR_Flag(Fmax,Lmax),High_Cloud_Flag(Fmax,Lmax),
     7    Cloud_IRTemp_Flag(Fmax,Lmax),Cloud_3p75_11_Flag(Fmax,Lmax),
     8    Cloud_VisRat_Flag(Fmax,Lmax),Cloud_SpatVar_Flag(Fmax,Lmax)
C     9    Cloud_CO2_Flag(Fmax,Lmax),Cloud_6p7_Flag(Fmax,Lmax)
      LOGICAL error_flag
C
C Initialization
C
      grpnm = ' '
      arrnm = 'Cloud_Mask'
      LinesPerScancube = 10
      Rank  =  3
C
C Retrieve the rank, dimensions and data type of SDS data.
C
      IF (GMARDM(modfil, arrnm, grpnm, data_type, Rank, Dim_Size)
     &   .ne.MAPIOK) THEN
         call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     &   'GMARDM get dim_info failed','ReadCldMsk_MOD04')
         error_flag = .TRUE.
      END IF

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
C       Do 80 l1 = 1, 10
C       Do 80 p1 = 1, 1354
C
C The Cloud Mask consists of 6 separate 1-byte words.  Increment
C memory buffer index by 6 to cycle through 1-km frames in a line
C

          I = I + 6
C
C         Set 250-m line and pixel offsets to corner of current 1-km
C         cell.  Then, loop over the 16 250-m pixels within the 1-km.
C         0-based bit positions 16-31 of the 1-km INTEGER cloud mask
C         value are associated with 250-m sublines (1-4) and samples
C         (1-4) with line varying most rapidly.
C
          l4_offset = (l1-1)*4
          p4_offset = (p1-1)*4
          l2_offset = (l1-1)*2
          p2_offset = (p1-1)*2

          DET_Flag(p1,l1) = ibits(count(I),0,1)
          UFQ_Flag(p1,l1) = ibits(count(I),1,2)
          DayNight_Flag(p1,l1) = ibits(count(I),3,1)
          SunGlint_Flag(p1,l1) = ibits(count(I),4,1)
          SnowIce_Flag(p1,l1) = ibits(count(I),5,1)
!          LandSea_Flag(p1,l1) = ibits(count(I),6,2)
          Non_CloudOb_Flag(p1,l1) = 1
          Thin_CirNIR_Flag(p1,l1) = ibits(count(I+1),1,1)
          Shadow_Flag(p1,l1) = 1
c    following four tests are used for ocean to compute High_cloud flag
c     Thin_CirIR_Flag                   Thin cirrus test
c     Cloud_SimpIR_Flag )            co2 threshold test
c     Cloud_3p75_11_Flag           6.7 um test
c     Cloud_IRTemp_Flag             IR temperature difference test
c
          Thin_CirIR_Flag(p1,l1) = ibits(count(I+1),3,1)
          Cloud_SimpIR_Flag(p1,l1) = ibits(count(I+1),6,1)
          Cloud_3p75_11_Flag(p1,l1) = ibits(count(I+1),7,1)
          Cloud_IRTemp_Flag(p1,l1) = ibits(count(I+2),2,1)
c
          Cloud_VisRat_Flag(p1,l1) = ibits(count(I+2),5,1)
          Cloud_SpatVar_Flag(p1,l1) = ibits(count(I+3),4,1)

C          WRITE(*,*) l1,p1,UFQ_Flag(p1,l1),LandSea_Flag(p1,l1),count(I),I
C          WRITE(*,*) 'Thin_CirIR,Cloud_IRTemp,Cloud_6p7,Cloud_CO2',
C     1  Thin_CirIR_Flag(p1,l1),Cloud_IRTemp_Flag(p1,l1),
C     2  Cloud_SimpIR_Flag(p1,l1),Cloud_3p75_11_Flag(p1,l1)
C
C-----------------------------------------------------------------------
C Examine first byte of cloud mask to determine whether cloud mask
C was even determined.
C (Zero-based bit 0 is 1 for determined, and 0 for not determined)
C
C If cloud mask not determined, set 250-m CldMsk(p4,l4) to 0
C Also set 1-km LandSea_Flag to 0, not processed.
C
C For Version 2, LandSea_Flag may take 4 values:
C [0 (water); 1 (coastal), 2 (desert), 3 (land)]
C-----------------------------------------------------------------------
C
C Default set to clear (in favor of clear pixels)
C 

             CldMsk_1km(p1,l1) = 1
             Cloud_VisRat_Flag(p1,l1)= 1

             Do 90 l = 1,4
             Do 90 p = 1,4
               l4 = l4_offset + l
               p4 = p4_offset + p
               CldMsk_250(p4,l4) = 1
   90        continue

             Do 91 l = 1,2
             Do 91 p = 1,2
               l2 = l2_offset + l
               p2 = p2_offset + p
               CldMsk_500(p2,l2) = 1
   91        continue

C
C High cloud flag for ocean; default set to cloudy if cloud mask is not determined 
C (0: cloudy; 1: clear)
C
c1. identifying the IR cirrus (mainly sensitive to thin cirrus clouds), 
c2.homogeneous thin cirrus and high clouds are the CO2 test
c(using 13.9 m band, sensitive mainly to high clouds in cold regions of the atmosphere), 
c3.the 6.7m test (sensitive to high clouds in cold regions of the atmosphere), 
c 4.and the Delta-IR test, all described in details by Ackerman et al. (1998).

          High_Cloud_Flag(p1,l1)=0
          IF(ibits(count(I),0,1).EQ.1) THEN
            IR_cloud=0
            High_Cloud_Flag(p1,l1)=1
            IF(Thin_CirIR_Flag(p1,l1).EQ.0)IR_cloud=IR_cloud+1
c            IF(Cloud_SimpIR_Flag(p1,l1).EQ.0)IR_cloud=IR_cloud+1
             IF(Cloud_3p75_11_Flag(p1,l1).EQ.0)IR_cloud=IR_cloud+1
            IF(Cloud_IRTemp_Flag(p1,l1).EQ.0)IR_cloud=IR_cloud+1
            IF(IR_cloud.GE.1) High_Cloud_Flag(p1,l1)=0
          ENDIF

C
C Using Wisonsin University Cloud Mask over land 
C

          If (ibits(count(I),0,1) .EQ. 0) Then

C If not determined, set cloud mask to be cloudy (=0) --> NOT processed
C Need to set to cloudy

!             LandSea_Flag(p1,l1) = -1

             CldMsk_1km(p1,l1) = 0
             Cloud_VisRat_Flag(p1,l1)=0

             Do 88 l = 1,4
             Do 88 p = 1,4
               l4 = l4_offset + l
               p4 = p4_offset + p
               CldMsk_250(p4,l4) = 0
   88        continue

             Do 89 l = 1,2
             Do 89 p = 1,2
               l2 = l2_offset + l
               p2 = p2_offset + p
               CldMsk_500(p2,l2) = 0
   89        continue

          Else
C-----------------------------------------------------------------------
C Cloud Mask determined.  Extract value of 1-km LandSea_Flag.
C Also get results (a simple yes/no switch) of 250-m visible cloud test.
C-----------------------------------------------------------------------
C
C Extract bits from 0 to 25 for use in land and ocean aerosol algorithm
C
C
C IF UFQ_Flag(p1,l1) NE 0 and 1 ---> indicating clear pixel
C Also for both 250m and 500m resolution pixels (see cloud mask ATBD)
C

            IR_cloud=0

            IF(UFQ_Flag(p1,l1).EQ.0.OR.UFQ_Flag(p1,l1).EQ.1) THEN
C     &         .OR.UFQ_Flag(p1,l1).EQ.2) THEN

                CldMsk_1km(p1,l1) = 0
                Cloud_VisRat_Flag(p1,l1)=0

               Do 95 l = 1,2
               Do 95 p = 1,2
                 l2 = l2_offset + l
                 p2 = p2_offset + p
                 CldMsk_500(p2,l2) = 0
   95          continue

               Do 96 l = 1,4
               Do 96 p = 1,4
                 l4 = l4_offset + l
                 p4 = p4_offset + p
                 CldMsk_250(p4,l4) = 0
   96          continue
              
            ELSE
C
C IF UFQ_Flag(p1,l1) EQ 0 and 1 ---> indicating cloudy pixel
C IR test will not overwirte the decision if more than 2
C tests indicate high-cloud free. 
C

              IF(Thin_CirIR_Flag(p1,l1).EQ.0)IR_cloud=IR_cloud+1
              IF(Cloud_SimpIR_Flag(p1,l1).EQ.0)IR_cloud=IR_cloud+1
              IF(Cloud_3p75_11_Flag(p1,l1).EQ.0)IR_cloud=IR_cloud+1
              IF(Cloud_IRTemp_Flag(p1,l1).EQ.0)IR_cloud=IR_cloud+1

              IF(IR_cloud.GE.2) THEN

                CldMsk_1km(p1,l1) = 0

                Do 97 l = 1,2
                Do 97 p = 1,2
                  l2 = l2_offset + l
                  p2 = p2_offset + p
                  CldMsk_500(p2,l2) = 0
   97           continue

                Do 98 l = 1,4
                Do 98 p = 1,4
                  l4 = l4_offset + l
                  p4 = p4_offset + p
                  CldMsk_250(p4,l4) = 0
   98           continue

              ENDIF

            ENDIF

C
C Non_cloud obstruction bit will overwrite the cloudy pixel (disabled for now)
C Non_CloudOb_Flag(p1,l1) is set to 1 (see above)
C

               IF(High_Cloud_Flag(p1,l1).EQ.1) THEN

 !               IF(LandSea_Flag(p1,l1).NE.0.AND.Non_CloudOb_Flag(p1,l1).EQ.0) THEN

                 CldMsk_1km(p1,l1) = 1

                 Do 92 l = 1,2
                 Do 92 p = 1,2
                   l2 = l2_offset + l
                   p2 = p2_offset + p
                   CldMsk_500(p2,l2) = 1
   92            continue


                 Do 93 p = 1,4
                   l4 = l4_offset + l
                   p4 = p4_offset + p
                   CldMsk_250(p4,l4) = 1
   93            continue

!              ENDIF
     
               ENDIF

          End If

   80  continue

C
C 3 x 3 windhowing: Cloud_SpatVar_Flag(p1,l1).EQ.0 is now replaced with dust bit (bit 28)
C

       Do 888 l1 = 1, Dim_Size(3)
       Do 888 p1 = 1, Dim_Size(2)

               IF(Cloud_VisRat_Flag(p1,l1).EQ.0) THEN
                 IF(l1.GT.1.AND.l1.LT.Dim_Size(3)) THEN
                 IF(p1.GT.1.AND.P1.LT.Dim_Size(2)) THEN
                   DO l11=l1-1,l1+1
                   DO p11=p1-1,p1+1
                     CldMsk_1km(p11,l11) = 0
                     CldMsk_1km(p1,l1) = 0
                   ENDDO
                   ENDDO
                 ENDIF 
                 ENDIF
               ENDIF 

  888 continue

       Do 85 l1 = 1, Dim_Size(3)
       Do 85 p1 = 1, Dim_Size(2)

          l4_offset = (l1-1)*4
          p4_offset = (p1-1)*4
          l2_offset = (l1-1)*2
          p2_offset = (p1-1)*2

          IF(CldMsk_1km(p1,l1).EQ.0) THEN

          Do 86 l = 1,2
          Do 86 p = 1,2
            l2 = l2_offset + l
            p2 = p2_offset + p
            CldMsk_500(p2,l2) = 0
   86     continue

          Do 87 l = 1,4
          Do 87 p = 1,4
            l4 = l4_offset + l
            p4 = p4_offset + p
            CldMsk_250(p4,l4) = 0
   87     continue

          ENDIF  

   85 continue

       Do 81 l1 = 1, Dim_Size(3)
       Do 81 p1 = 1, Dim_Size(2)

               IF(Cloud_SpatVar_Flag(p1,l1).EQ.0) THEN
                 IF(l1.GT.1.AND.l1.LT.Dim_Size(3)) THEN
                 IF(p1.GT.1.AND.P1.LT.Dim_Size(2)) THEN
                   DO l11=l1-1,l1+1
                   DO p11=p1-1,p1+1
                     CldMsk_1km(p11,l11) = 1
                   ENDDO
                   ENDDO
                 ENDIF 
                 ENDIF
               ENDIF 

   81 continue

       Do 82 l1 = 1, Dim_Size(3)
       Do 82 p1 = 1, Dim_Size(2)

          l4_offset = (l1-1)*4
          p4_offset = (p1-1)*4
          l2_offset = (l1-1)*2
          p2_offset = (p1-1)*2

          IF(CldMsk_1km(p1,l1).EQ.1) THEN

          Do 83 l = 1,2
          Do 83 p = 1,2
            l2 = l2_offset + l
            p2 = p2_offset + p
            CldMsk_500(p2,l2) = 1
   83     continue

          Do 84 l = 1,4
          Do 84 p = 1,4
            l4 = l4_offset + l
            p4 = p4_offset + p
            CldMsk_250(p4,l4) = 1
   84     continue

          ENDIF  

   82 continue

C            write(*,*)
C            write(*,*)
C            write(*,'(A)') 'first 11 UFQ Flags for ' //
C     &                             '10 lines'
C            write(*,'(11I4)') ((UFQ_Flag(i,j),i=1,11),
C     &                                  j=1,10)
C            write(98,'(30I4)') ((CldMsk_1km(i,j),i=1,30),
C     &                                  j=1,10)

       Return
       END
