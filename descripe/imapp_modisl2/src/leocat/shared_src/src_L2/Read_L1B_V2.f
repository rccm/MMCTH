      subroutine Read_L1B(L1B_MAPI_fhandle, L1B_fname,
     &           Scan_No, Band_No, Gain, Resol, Cal_type,
     &           BufSize_Xtrack, BufSize_Track, Buf_EV, Buf_Un, Buf_Sa,
     &           V_flag, Data_size, error_flag)

      implicit none

      include 'PGS_PC.f'
      include 'PGS_MODIS_39500.f'
      include 'L1B_Reader.inc'
      include 'mapi.inc'

C----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C
C    Retrieve one scan of one spectral channel of MODIS Level 1B sensor 
C    data.  Calibrated data are returned in a choice of radiance, 
C    reflectance or corrected counts units.
C    The user identifies a unique MODIS data channel by passing a 
C    spectral band number between 1 and 36, and the channel gain spec-
C    ification.  The gain is high (H) or low (L) for bands 13 and 14,
C    and for all other bands it is ignored.  Additional inputs are MODIS
C    scan number for the granule, and one of the following valid cali-
C    bration types: reflectance (refl), radiance (rad), or corrected 
C    counts (counts) units.  Reflectance units are an invalid type for 
C    the MODIS thermal bands.  To allow processing on more than one L1B 
C    granule, the user also supplies the name of the L1B HDF file, and 
C    its HDF SD and VS interface Ids and access mode through an array 
C    structure.  The user also supplies four 2-dimensional buffers and 
C    their dimension sizes to hold the retrieved L1B sensor data, an 
C    uncertainty risk index, aggregate data sample size, and a data 
C    validity flag.  The validity flag is 0 (valid) or -1 (invalid). 
C    The aggregate data sample size contains the number of high spatial 
C    resolution L1B observations integrated up to form the aggregate 
C    pixel values.  On requests for non-aggregated bands (that is, for 
C    native sensor data), the number-of-samples buffer is still required
C    and contains the value 1 for all pixels except for missing data.
C    The dimension sizes are used internally to compare the buffer size 
C    against the size of the retrieved L1B arrays.  If they are too 
C    small to hold the L1B data, an error message is reported to the SDP
C    LogStatus file, and control is returned back to the calling routine.
C
C !INPUT PARAMETERS:
C
C    integer L1B_MAPI_fhandle(*)   
C                         Array containing HDF SD (element 1) and 
C                         VS (element 2) interface IDs, and the L1B file 
C                         access mode (read; read/write, or create in
C                         element 3).  Additional elements carrying 
C                         pointers to C language structures may also be
C                         passed, and are used internally by the MODIS
C                         Application Program Interface (M-API) library 
C                         to maximize CPU efficiency.
C
C    character*(*) L1B_fname  
C                         Name of L1B product input file
C
C    integer Scan_No      1-based MODIS scan number relative to start 
C                         of granule
C
C    integer Band_No      MODIS spectral band number from 1 to 36
C
C    character*(*) Gain   Radiometric Gain specification, high ("H") or
C                         low ("L").  It is ignored for all bands 
C                         except 13 and 14.
C
C    integer Resol        Spatial resolution factor relative to 1 km.
C                         Possible values are 250 m (16), 500 m (4), 
C                         and 1 km (1).
C
C    character*(*) Cal_type   
C                         Calibration format passed as one of three 
C                         types:  reflectances ("refl"), radiances 
C                         ("rad"), or corrected counts ("counts").
C
C    integer BufSize_Xtrack 
C                         Size of first dimension of output buffer as 
C                         declared in calling program.  All three output
C                         buffers (Buf_EV, Buf_Un and Buf_Sa) are 
C                         assumed to have the same dimension sizes.
C
C    integer BufSize_Track  
C                         Size of second dimension of output buffer as 
C                         declared in calling program.  All three output 
C                         buffers (Buf_EV, Buf_Un and Buf_Sa) are 
C                         assumed to have the same dimension sizes.
C
C !OUTPUT PARAMETERS: 
C
C    real      Buf_EV     Two dimensional (2-D) buffer for passing L1B 
C                         Earth View (EV) data.  Buffer 
C                         index 1 represents frame number within a scan,
C                         index 2 is line number within scan.
C    byte      Buf_Un     2-D buffer for passing L1B uncertainty index.
C                         Buffer index 1 represents frame number within 
C                         a scan.  Index 2 is line number within scan.
C    byte      Buf_Sa     2-D buffer for passing sample sizes used to
C                         compute L1B aggregate data.  Buffer 
C                         index 1 represents frame number within a scan.
C                         Index 2 is line number within scan.
C    byte      V_flag     2-D buffer identifying locations of valid (0)
C                         and invalid (-1) L1B data (see L1B file spec.
C                         for definition of non-valid data).
C    integer   Data_size  Array containing size of data segment within 
C                         output buffers (Buf_EV, Buf_Un and Buf_Sa).
C                         Element 1 contains X-track dimension size and
C                         element 2 the along track dimension size.
C    logical  error_flag  Error indicator
C
C !REVISION HISTORY:
C Revision 1.3  1998/06/18  13:55:44  rhucek
C Added "Operator Action" strings to error messages.  Optimized code
C according to Fred Gunther's suggestions, i.e. eliminated assignment
C statements TempEV within do loop 90.  Fixed problem with NewFileUpdate
C error flag as recommended by Liam Gumley.
C
C Revision 1.2  1998/04/23 15:38:37  gumley
C Fixed problem with NewFileUpdate error flag not being set.
C
C Revision 1.1.1.1  1998/04/22 13:54:33  gumley
C V2 baselined code from Rich Hucek
C
C Revision 1.2  1997/08/12  16:29:49  rhucek
C Added explicit type declarations for dummy arguments
C Gain and Cal_type [character*(*)], and Band_No and Resol
C (integer).
C
C Revision 1.1  1997/04/15  15:36:41  vlin
C Initial revision
C
C !TEAM-UNIQUE HEADER:
C
C    This software is developed by the MODIS Science Data Support Team
C    for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C !REFERENCES AND CREDITS:
C
C    Written by Vicky Lin          April 1997
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
C  Externals
C  ---------
C    GMAR                         (mapi.inc; libmapi.a)
C    MAPIOK                       (mapi.inc)
C
C  Partial Data Dictionary
C  -----------------------
C  Sa_data        Buffer large enough to hold number of samples in 
C                 aggregate data for a scan of any band of MODIS data.
C  Un_data        Buffer large enough to hold L1B uncertainty index for 
C                 a scan of any band of MODIS data
C  EV_name        Array containing names of L1B EV arrays for 4 arrays 
C                 and 3 spatial resolutions.
C  Last_L1B_fname Array containing names of last L1B file accessed for 
C                 each spatial resolution.
C  Sa_array_name  L1B number of samples array name
C  Sa_Name        Array containing names of L1B sample size arrays for 
C                 4 arrays and 3 spatial resolutions
C  Un_array_name  L1B Uncertainty array name
C  Un_name        Array containing name of L1B uncertainty arrays for 
C                 4 arrays and 3 spatial resolutions
C  EV_data        Buffer large enough to hold a scan of L1B Earth View 
C                 data for any band
C  L1B_Missing    Missing L1B data marker, PARAMETER = -1
C  A_index        Array type index:  1 (250m solar), 2 (500m solar), 
C                                    3 (1 km solar), 4 (1 km thermal).
C  Cal_index      Calibration index:  1 reflectance, 2 radiance, 
C                                     3 corrected counts.
C  Ch_index       Array channel index: 
C                               1-2 (250m solar), 1-5 (500m solar), 
C                               1-15 (1km solar) and 1-16 (1km thermal).
C  Dim_size       Dimension sizes of the L1B arrays
C  cal_size       The slab size of L1B data to extract
C  LinesPerScan   Number of lines in scan of MODIS data:  40 (250m), 
C                 20 (500m) and 10 (1km)
C  No_Scans       Number of scans in L1B granule
C  R_index        Spatial resolution index: 1(250m), 2(500m) and 3(1km)
C  Sa_flag        Flag indicating whether sample size array exists for 
C                 band of L1B data
C  EV_offset      Array of L1B offsets for 16 channels, 3 calibrations, 
C                 4 array types, and 3 spatial resolutions
C  EV_slope       Array of L1B slopes for 16 channels, 3 calibrations, 
C                 4 array types, and 3 spatial resolutions
C
C !END
C---------------------------------------------------------------------

      SAVE Last_L1B_fname,Dim_size,EV_offset,EV_slope
      integer i,j,element,BufSize_Xtrack,BufSize_Track

c-----PARAMETER declarations
      character*(*)  FUNCNAME
      parameter     (FUNCNAME = 'Read_L1B')

      byte       L1B_INVALID,      L1B_VALID
      parameter (L1B_INVALID = -1, L1B_VALID = 0)

      integer*2  L1B_MISSING 
      parameter (L1B_MISSING = -1)

      integer    SUCCEED
      parameter (SUCCEED = 0)

      integer    CH_DIM_IDX,     TRACK_DIM_IDX,     XTRACK_DIM_IDX
      parameter (CH_DIM_IDX = 3, TRACK_DIM_IDX = 2, XTRACK_DIM_IDX = 1)

      byte Sa_data(1500*10*16),
     &     UN_data(1500*10*16),Buf_Un(BufSize_Xtrack,BufSize_Track),
     &     Buf_Sa(BufSize_Xtrack,BufSize_Track),
     &     V_flag(BufSize_Xtrack,BufSize_Track)

      character*(*) Cal_type,Gain,L1B_fname
      character*6   msg6
      character*10  msg10,grpnm
      character*2048 msgbuf
      character*(PGSd_PC_FILE_PATH_MAX) Last_L1B_fname(3)

      integer*2 EV_data(1500*10*16)
      integer A_Index,Cal_Index,Ch_index,No_Scans,R_Index,
     &        Scan_No,rtn,cal_size(3),
     &        Dim_size(3,4,3),Start(3),Data_Size(*),
     &        Band_No,L1B_MAPI_fhandle(*),Resol
      integer Initial_Input_Checks

      logical Sa_flag,error_flag,file_ErrFlag

      real scale,add_offset,Buf_EV(BufSize_Xtrack,BufSize_Track),
     &     EV_slope(16,3,4,3),EV_offset(16,3,4,3)

      data Last_L1B_fname / 3*' ' /, Dim_size / 36*0 /, scale / 0.0 /,
     &     add_offset / 0.0 /

c ... Initialization
      grpnm = ' '
      Sa_flag = .false.
      error_flag = .false.

c ... Conduct quality assurance tests on Read_L1B input arguments.  
c     If the return value is not zero, return to calling routine
 
      rtn = Initial_Input_Checks(L1B_MAPI_fhandle,L1B_fname,
     1                           Band_No,Gain,Cal_Type,Resol)

      if (rtn .ne. SUCCEED) then
          error_flag = .true.

          msgbuf = 'Invalid input data detected by Initial_Input_Checks'
     1     // char(10) // 'Operator Action:  Refer to prior low level LogStatus '
     2     // char(10) // 'messages originating within routine Initial_Input_Checks. ' 

          call Modis_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
          return
      end if
 

c ... Retrieve Read_L1B software internal design indices for 
c     referencing elements of the L1B product.

      Call Get_Design_Indices(Resol,Band_No,Gain,Cal_Type,
     &                        R_index,A_Index,Ch_index,Cal_Index,
     &                        LinesPerScan)

      EV_array_name = EV_name(A_index, R_index)
      Un_array_name = Un_name(A_index, R_index)
      Sa_array_name = Sa_name(A_index, R_index)

c ... If there has been a change of L1B granules, then 
c     update the arrays containing SDS dimension sizes and 
c     the integer to real data conversions factors

      if (L1B_fname.ne.Last_L1B_fname(R_index)) then
         file_ErrFlag = .false.

         call NewFileUpdate(L1B_MAPI_fhandle,R_index,
     1        EV_offset,EV_slope,Dim_size,file_ErrFlag)

         if (file_ErrFlag) then
            error_flag = .true.

            msgbuf = 'Read_L1B detected error from routine NewFileUpdate. '
     1       // char(10) // 'Operator Action:  Refer to prior low level LogStatus '
     2       // char(10) // 'messages originating within routine NewFileUpdate. ' 

            call Modis_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
            return
         end if
      endif


c ... Perform range checking on input scan number and on the 
c     data buffer size provided by calling routine.

      No_Scans = Dim_size(TRACK_DIM_IDX,A_index,R_index)/LinesPerScan

      if (Scan_No.lt.1 .OR. Scan_No.gt.No_Scans) then
          error_flag = .true.
          write(msg6,'(i6)') Scan_No

          msgbuf = 'Scan number out of bounds: = ' // msg6
     1     // char(10) // 'Operator Action:  Suspect system problem.  If a fault is '
     2     // char(10) // 'identified, correct and rerun PGE.  Otherwise, notify SDST.' 

          call Modis_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      end if


      Data_size(XTRACK_DIM_IDX) = Dim_size(XTRACK_DIM_IDX,A_index,R_index)
      Data_size(TRACK_DIM_IDX)  = LinesPerScan

      If (Data_size(XTRACK_DIM_IDX).gt.BufSize_Xtrack) then
          error_flag = .true.
          write(msg10,'(i10)') BufSize_Xtrack
          msgbuf = 'Across track data size greater than buffer size: = ' // msg10
     1     // char(10) // 'Operator Action:  Suspect system problem.  If a fault is '
     2     // char(10) // 'identified, correct and rerun PGE.  Otherwise, notify SDST.' 

          call Modis_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      end if


      if (Data_size(TRACK_DIM_IDX).gt.BufSize_Track) then
          error_flag = .true.
          write(msg10,'(i10)') BufSize_Track
          msgbuf = 'Along track data size greater than buffer size: = ' // msg10
     1     // char(10) // 'Operator Action:  Suspect system problem.  If a fault is '
     2     // char(10) // 'identified, correct and rerun PGE.  Otherwise, notify SDST.' 

          call Modis_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      end if


c ... Return if additional QA checks reveal an error.
      if (error_flag) return

c ... Set start point and edge dimensions of desired L1B data slab.
      Start(XTRACK_DIM_IDX)    = 0
      Start(TRACK_DIM_IDX)     = (Scan_No-1) * LinesPerScan
      Start(CH_DIM_IDX)        = Ch_index - 1
      cal_size(XTRACK_DIM_IDX) = Data_size(XTRACK_DIM_IDX)
      cal_size(TRACK_DIM_IDX)  = Data_size(TRACK_DIM_IDX)
      cal_size(CH_DIM_IDX)     = 1


c ... Get Earth View (EV) and Uncertainty data slabs.  
c     If data are aggregates, also get Sample sizes.

c ... Get Earth View data ... 
      rtn = gmar(L1B_MAPI_fhandle,EV_array_name,grpnm,start,cal_size,EV_data)

      if (rtn.ne.MAPIOK) then
          error_flag = .true.
          write(msg6,'(i6)') Band_No
          msgbuf = 'M-API function GMAR detected error while retrieving L1B EV data '
     1     // 'for band ' // msg6
     2     // char(10) // 'Operator Action:  Suspect system problem.  If a fault is '
     3     // char(10) // 'identified, correct and rerun PGE.  Otherwise, notify SDST.' 

          call Modis_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      end if


c ... Get Uncertainty data ...
      rtn = gmar(L1B_MAPI_fhandle,UN_array_name,grpnm,start,cal_size,UN_data)

      if (rtn.ne.MAPIOK) then
          error_flag = .true.
          write(msg6,'(i6)') Band_No
          msgbuf = 'M-API function GMAR detected error while retrieving L1B uncertainty '
     1     // 'data for band ' // msg6
     2     // char(10) // 'Operator Action:  Suspect system problem.  If a fault is '
     3     // char(10) // 'identified, correct and rerun PGE.  Otherwise, notify SDST.' 

          call Modis_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      end if


c ... Get aggregate data sample sizes ...

      if (Sa_array_name.ne.' ') then
          Sa_flag = .true.
          rtn = gmar(L1B_MAPI_fhandle,Sa_array_name,grpnm,start,cal_size,Sa_data)

          if (rtn.ne.MAPIOK) then
             error_flag = .true.
             write(msg6,'(i6)') Band_No
             msgbuf = 'M-API function GMAR detected error while retrieving L1B '
     1        // char(10) // 'aggregate data Sample Sizes for band ' // msg6
     2        // char(10) // 'Operator Action:  Suspect system problem.  If a fault is '
     3        // char(10) // 'identified, correct and rerun PGE.  Otherwise, notify SDST.' 

             call Modis_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
          end if
      end if


      if (error_flag) return

c ... Set EV data calibration slope and offset.

      scale      = EV_slope(Ch_index, Cal_Index, A_index, R_index)
      add_offset = EV_offset(Ch_index, Cal_Index, A_index, R_index)

c ... Loop over lines and pixels; check for missing or invalid data.

      do 90 j = 1, Data_size(TRACK_DIM_IDX)
      do 90 i = 1, Data_size(XTRACK_DIM_IDX)

         element = (j-1) * Data_size(XTRACK_DIM_IDX) + i
         if (EV_data(element).lt.0) then

c        ... Data are missing or are other invalid type ...
             V_flag(i,j) = L1B_INVALID

             if (EV_data(element).eq.L1B_MISSING) then

c            ... Data are missing ...
                 Buf_EV(i,j) = -9999.0
                 Buf_Un(i,j) = -99
                 Buf_Sa(i,j) = -99

             else

c            ... Data are other invalid types; calibrate anyway ...
                 Buf_EV(i,j) = scale * ibits(EV_data(element),0,15) + add_offset
                 Buf_Un(i,j) = ibits(Un_data(element),0,4)
                 Buf_Sa(i,j) = 1
                 if (Sa_flag) Buf_Sa(i,j) = Sa_data(element)

             end if

         else

c        ... Data are valid; set validity flag to 0. ...
             V_Flag(i,j) = L1B_VALID
             Buf_EV(i,j) = scale * ibits(EV_data(element),0,15) + add_offset
             Buf_Un(i,j) = ibits(Un_data(element),0,4)
             Buf_Sa(i,j) = 1
             if (Sa_flag) Buf_Sa(i,j) = Sa_data(element)

         end if

   90 continue

      Last_L1B_fname(R_index) = L1B_fname

      return
      end
