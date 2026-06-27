      subroutine Get_Design_Indices(Resol,Band_No,Gain,Cal_Type,
     &           R_Index,A_Index,Ch_Index,Cal_Index,LinesPerScan)
      implicit none
      include 'L1B_Reader.inc'
  
C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION: 
C
C   Retrieve the Read_L1B software internal design indices used for 
C   referencing MODIS L1B sensor data by spatial resolution, 
C   calibration type, and (HDF SDS) array name and channel.
C   Subroutine Get_Design_Indices converts the character string and 
C   integer arguments passed from the application program to the 
C   Read_L1B software into integer array indices.  These indices 
C   represent classifications of the data in terms of MODIS channel 
C   number, calibration type (reflectance, radiance or counts), 
C   and spatial resolution.
C
C !INPUT PARAMETERS:
C
C  integer Resol      Spatial resolution factor relative to 1 km.  
C                     Possible values are 250 M (16), 500 M (4) and 
C                     1 KM (1).
C  integer Band_No    MODIS spectral band number from 1 to 36
C  character Gain     Radiometric gain specification, high ("H") or 
C                     low ("L").  It is ignored for all bands except 
C                     13 and 14.
C  character Cal_type Calibration format in one of three types:  
C                     reflectance ("Refl"), radiances ("Rad"), or
C                     corrected counts ("Counts").
C
C !OUTPUT PARAMETERS:
C
C  integer  R_index   Spatial resolution index of MODIS sensor data, 
C                     including aggregate data: either 1 (250 m), 
C                     2 (500 m) or 3 (1 km).
C  integer  A_index   L1B array identifier from 1 to 4 representing the
C                     250 m, 500 m, 1 km solar reflectance, and 1 km 
C                     thermal emissive band data arrays, respectively.
C  integer  Ch_index  Dimension index identifying the plane within the 
C                     multidimensional L1B array containing the 
C                     specified channel data.  Ch_index may range from 
C                     1 to 16.
C  integer  Cal_index Calibration type representing 1 (reflectances), 
C                     2 (radiances), and 3 (corrected counts).
C  integer  LinesPerScan    
C                     The number of detector lines in a scan of 
C                     aggregated MODIS data.  May have values of 
C                     40 (250 m), 20 (500 m) or 10 (1 km)
C
C !REVISION HISTORY:
C Revision 1.2  1997/08/12  13:01:27  rhucek
C Added explicitly type declarations for dummy arguments
C Gain and Cal_type [character*(*)], and Band_No and Resol
C (integer).
C
C Revision 1.1  1997/04/09  19:49:29  vlin
C Initial revision
C
C !TEAM-UNIQUE HEADER:
C
C    This software was developed by the MODIS Science Data Support Team
C    for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C !REFERENCES AND CREDITS:
C
C    Written by    
C    S. Vicky Lin                  March 1997 
C    *** Under guidance of Richard Hucek  ***
C
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd, Seabrook MD 20706
C    rhucek@ltpmail.gsfc.nasa.gov
C    vlin@modis1.gsfc.nasa.gov
C
C !DESIGN NOTES:
C
C  Externals:
C    Named Constants:
C        Band_Is_*         (L1B_Reader.inc)
C        Cal_Is_*          (L1B_Reader.inc)
C        First_*_Band      (L1B_Reader.inc)
C        High_Gain         (L1B_Reader.inc)
C        Last_*_Band       (L1B_Reader.inc)
C        LinesPerScan_*    (L1B_Reader.inc)
C        Resol_Is_*        (L1B_Reader.inc)
C
C !END
C---------------------------------------------------------------------

      character*(*) Cal_Type, Gain

      integer Band_No, Resol
      integer R_Index,A_Index,Ch_Index,Cal_Index

c ... Determine design indices R_index, A_index, Ch_index, and the 
c     number of aggregate data lines in a MODIS scan LinesPerScan

      if (Resol.eq.Resol_Is_250m) then

         R_Index = 1
         A_Index = 1
         Ch_Index = Band_No
         LinesPerScan = LinesPerScan_250m

      else if (Resol.eq.Resol_Is_500m) then

         R_Index = 2
         LinesPerScan = LinesPerScan_500m
         if (Band_No.le.Last_250m_Band) then
             A_Index = 1
             Ch_Index = Band_No
         else
             A_Index = 2
             Ch_Index = Band_No - Last_250m_Band
         end if

      else 

         R_Index = 3
         LinesPerScan = LinesPerScan_1km

         if (Band_No.le.Last_250m_Band) then
             A_Index = 1
             Ch_Index = Band_No
         else if (Band_No.le.Last_500m_Band) then
             A_Index = 2
             Ch_Index = Band_No - Last_250m_Band
         else if (Band_No.lt.First_Thermal_Band .OR.
     &            Band_No.eq.Last_Solar_Band) then
             A_Index = 3
             if (Band_No.lt.Band_Is_13) then
                 Ch_Index = Band_No - 7
             else if (Band_No.eq.Band_Is_13) then
                 Ch_Index = 6
                 if (Gain.eq.High_Gain) Ch_Index = 7
             else if (Band_No.eq.Band_Is_14) then
                 Ch_Index = 8
                 if (Gain.eq.High_Gain) Ch_Index = 9
             else
                 Ch_Index = Band_No - 5
                 if (Band_No.eq.Last_Solar_Band) Ch_Index = 15
             end if
         else
             A_Index = 4
             Ch_Index = Band_No - 19
             if (Band_No.gt.Last_Solar_Band) 
     &       Ch_Index = Band_No - First_Thermal_Band
         end if

      end if

c ... Assign calibration format design index, Cal_index

      if (Cal_Type.eq.Cal_Is_Refl) then
          Cal_Index = 1
      else if (Cal_Type.eq.Cal_Is_Rad) then
          Cal_Index = 2
      else 
          Cal_Index = 3
      end if

      return
      end
