      Subroutine Read_L1B( L1B_MAPI_fhandle, L1B_fname,
     1                     Scan_NUM, Band_NUM, Gain, Resol, Cal_Type,
     2                     BufSize_Xtrack, BufSize_Track,
     3                     Buf_EV, Buf_Un, Buf_Sa,
     4                     V_flag, Data_size, error_flag )

      implicit none

      include 'PGS_PC.f'
      include 'PGS_MODIS_39500.f'
      include 'Read_L1B_Aerosol.inc'
      include 'extra.inc'
      include 'mapi.inc'

C----------------------------------------------------------------------
C !F77
C
C
C !DESCRIPTION:
C
C    Retrieve one scan of one spectral channel of MODIS Level 1B sensor
C    data.  Calibrated Earth View (EV) data are returned in a choice of
C    reflectance, radiance, or corrected counts units for solar 
C    reflective bands (bands 1-19 and 26), and radiance units only for 
C    emissive bands (bands 20-36, but not 26).  Uncertainty flags are 
C    also returned (See Design Notes).
C
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
C    integer Scan_NUM     1-based MODIS scan number relative to start
C                         of granule
C
C    integer Band_NUM     MODIS spectral band number from 1 to 36
C
C    character*(*) Gain   Radiometric Gain specification, high ("H") or
C                         low ("L").  It is ignored for all bands
C                         except 13 and 14.
C
C    integer Resol        Spatial resolution factor relative to 1 km.
C                         Possible values are 16 (250 m), 4 (500 m),
C                         and 1 (1 km).
C
C    character*(*) Cal_Type
C                         Calibration format passed as one of three
C                         types:  reflectance ("refl"), radiance
C                         ("rad"), or corrected counts ("counts").
C                         Only "rad" is valid for emissive band data
C                         (20-36, but not 26). 
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
C
C !REVISION HISTORY:
C
C 05/19/04 rhucek:  
C   Corrected an error in which SDS scale and offset terms were applied 
C   to "unusable" Earth View stored integers.  Unusable Earth View data 
C   are represented by reserved integer values identifying the reason 
C   for non-usability.  No detector signal is included in stored integer
C   value. 
C
C
C !TEAM-UNIQUE HEADER:
C
C    This software is developed by the MODIS Science Data Support Team
C    for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C
C !REFERENCES AND CREDITS:
C
C
C !DESIGN NOTES:
C
C    The user identifies a unique MODIS data channel by passing a
C    spectral band number between 1 and 36, and the channel gain spec-
C    ification.  The gain is high (H) or low (L) for bands 13 and 14,
C    and for all other bands it is ignored.  Additional inputs are MODIS
C    scan number for the granule, and one of the following valid cali-
C    bration types: reflectance ("refl"), radiance ("rad"), or corrected
C    counts ("counts") units.  Reflectance and corrected counts units are 
C    invalid types for the MODIS thermal bands.  To allow processing on more 
C    than one L1B granule, the user also supplies the name of the L1B HDF 
C    file, and its HDF SD and VS interface Ids and access mode through an 
C    array structure.  The user also supplies four 2-dimensional buffers 
C    and their dimension sizes to hold the retrieved L1B sensor data, an
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
C  Externals
C  ---------
C    GMAR                         (mapi.inc; libmapi.a)
C    MAPIOK                       (mapi.inc)
C
C    Named Constants:
C     BAND_IS_*                   (L1B_Reader_V2.1.inc)
C     DIM_NUM_CH                  (L1B_Reader_V2.1.inc)
C     DIM_NUM_TRACK               (L1B_Reader_V2.1.inc)
C     DIM_NUM_XTRACK              (L1B_Reader_V2.1.inc)
C     FAIL                        (L1B_Reader_V2.1.inc)
C     MAX_L1B_NAME_SIZE           (L1B_Reader_V2.1.inc)
C     LINES_PER_SCAN_1KM          (L1B_Reader_V2.1.inc)
C     PIXELS_PER_FRAME_250M       (L1B_Reader_V2.1.inc)
C     MAX_FRAMES_PER_SCAN         (L1B_Reader_V2.1.inc)
C     MAX_NUM_ARRAYS              (L1B_Reader_V2.1.inc)
C     MAX_NUM_CHANNELS            (L1B_Reader_V2.1.inc)
C     MAX_NUM_CAL_TYPES           (L1B_Reader_V2.1.inc)
C     NUM_L1B_FILES               (L1B_Reader_V2.1.inc)
C     SUCCEED                     (L1B_Reader_V2.1.inc)
C
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
C  Dim_Size       Dimension sizes of the L1B arrays
C  Slab_Size      The slab size of L1B data to extract
C  LinesPerScan   Number of lines in scan of MODIS data:  40 (250m),
C                 20 (500m) and 10 (1km)
C  NUM_Scans      Number of scans in L1B granule
C  F_Index        L1B file index: 1 (250m), 2 (500m) and 3 (1km)
C  Sa_flag        Flag indicating whether sample size array exists for
C                 band of L1B data
C  EV_offset      Array of L1B offsets for MAX_NUM_CHANNELS,
C                 MAX_NUM_CAL_TYPES, MAX_NUM_ARRAYS, and NUM_L1B_FILES.
C  EV_slope       Array of L1B slopes for MAX_NUM_CHANNELS,
C                 MAX_NUM_CAL_TYPES, MAX_NUM_ARRAYS, and NUM_L1B_FILES.
C
C !END
C---------------------------------------------------------------------

      SAVE Last_L1B_fname, Dim_Size, EV_offset, EV_slope

c-----PARAMETER declarations
      character*(*)  FUNCNAME
      parameter     (FUNCNAME = 'Read_L1B')


c-----Function argument declarations
      integer BufSize_Xtrack,BufSize_Track

      byte Buf_Un(BufSize_Xtrack,BufSize_Track),
     1     Buf_Sa(BufSize_Xtrack,BufSize_Track),
     2     V_flag(BufSize_Xtrack,BufSize_Track)

      character*(*)  Cal_Type, Gain, L1B_fname

      integer L1B_MAPI_fhandle(*), Data_Size(*)

      integer Band_NUM, Resol, Scan_NUM

      logical error_flag

      real Buf_EV(BufSize_Xtrack,BufSize_Track)


c-----Local variable and array declarations
      byte Sa_data( MAX_FRAMES_PER_SCAN * LINES_PER_SCAN_1KM * PIXELS_PER_FRAME_250M ),
     1     UN_data( MAX_FRAMES_PER_SCAN * LINES_PER_SCAN_1KM * PIXELS_PER_FRAME_250M )

      character*6    msg6
      character*10   msg10
      character*2048 msgbuf
      character*(MAX_L1B_NAME_SIZE)  EV_array_name, Sa_array_name, UN_array_name, grpnm
      character*(PGSd_PC_FILE_PATH_MAX) Last_L1B_fname(NUM_L1B_FILES)

      integer A_Index,Cal_Index,Ch_index,LinesPerScan,NUM_Scans,F_Index,
     1        element,i,j,rtn,
     2        Start(MAX_NUM_DIMS),
     3        Slab_Size(MAX_NUM_DIMS),
     4        Dim_Size(MAX_NUM_DIMS,MAX_NUM_ARRAYS,NUM_L1B_FILES)

      integer Initial_Input_Checks

      integer*2 EV_data(MAX_FRAMES_PER_SCAN * LINES_PER_SCAN_1KM * PIXELS_PER_FRAME_250M)

      logical Sa_flag,file_ErrFlag

      real scale,add_offset,
     1     EV_slope (MAX_NUM_CHANNELS, MAX_NUM_CAL_TYPES, MAX_NUM_ARRAYS, NUM_L1B_FILES),
     2     EV_offset(MAX_NUM_CHANNELS, MAX_NUM_CAL_TYPES, MAX_NUM_ARRAYS, NUM_L1B_FILES)


      DATA Last_L1B_fname / 3*' ' /,
     1     Dim_Size       / 36*0 /,
     2     scale          / 0.0 /,
     3     add_offset     / 0.0 /


c-------------------------------------------------------------------------------
c Initialization
c-------------------------------------------------------------------------------

      grpnm      = ' '
      Sa_flag    = .false.
      error_flag = .false.
        

c-------------------------------------------------------------------------------
c Test for valid input arguments
c-------------------------------------------------------------------------------

      rtn = Initial_Input_Checks( L1B_MAPI_fhandle, L1B_fname,
     1                            Band_NUM, Gain, Cal_Type, Resol )

      if (rtn .ne. SUCCEED) then
          error_flag = .true.

          msgbuf =
     1    'Invalid input data detected by Initial_Input_Checks'
     2    // char(10) // 'Operator Action:  Refer to prior low level LogStatus '
     3    // char(10) // 'messages originating within routine Initial_Input_Checks.'

          call Modis_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
          return
      end if
       
     
c-------------------------------------------------------------------------------
c Get L1B reader internal design indices
c-------------------------------------------------------------------------------

      Call Get_Design_Indices( Resol, Band_NUM, Gain, Cal_Type,
     1                         F_Index, A_Index, Ch_index, Cal_Index,
     2                         LinesPerScan )

      EV_array_name = EV_name(A_index, F_Index)
      Un_array_name = Un_name(A_index, F_Index)
      Sa_array_name = Sa_name(A_index, F_Index)

       
c-------------------------------------------------------------------------------
c Renew slope and offset factors, and dimension sizes when changing granules
c-------------------------------------------------------------------------------

      if (L1B_fname .ne. Last_L1B_fname(F_Index)) then
         file_ErrFlag = .false.

         Call NewFileUpdate( L1B_MAPI_fhandle, F_Index,
     1                       EV_offset, EV_slope,
     2                       Dim_Size, file_ErrFlag )

         if (file_ErrFlag) then
            error_flag = .true.

            msgbuf = 'Error detected by NewFileUpdate.'
     1       // char(10) // 'Operator Action:  Refer to prior low level LogStatus '
     2       // char(10) // 'messages originating within routine NewFileUpdate.'

            call Modis_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
            return
         end if
      endif
       

c-------------------------------------------------------------------------------
c Perform range checking on input scan number and on data buffer size supplied
c calling code
c-------------------------------------------------------------------------------

      NUM_Scans = Dim_Size(DIM_NUM_TRACK,A_index,F_Index)/LinesPerScan
       
      if (Scan_NUM.lt.1 .OR. Scan_NUM.gt.NUM_Scans) then
          error_flag = .true.
          write(msg6,'(i6)') Scan_NUM

          msgbuf = 'Scan number out of bounds: = ' // msg6
     1     // char(10) // 'Operator Action: Notify SDST.'

          call Modis_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      end if


       Data_size(DIM_NUM_XTRACK) = Dim_Size(DIM_NUM_XTRACK,A_index,F_Index)
       Data_size(DIM_NUM_TRACK)  = LinesPerScan
        
      If (Data_size(DIM_NUM_XTRACK).gt.BufSize_Xtrack) then
          error_flag = .true.
          write(msg10,'(i10)') BufSize_Xtrack

          msgbuf = 'Across track data size greater than buffer size: = ' // msg10
     1     // char(10) // 'Operator Action: Notify SDST.'

          call Modis_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      end if


      if (Data_size(DIM_NUM_TRACK).gt.BufSize_Track) then
          error_flag = .true.
          write(msg10,'(i10)') BufSize_Track

          msgbuf = 'Along track data size greater than buffer size: = ' // msg10
     1     // char(10) // 'Operator Action: Notify SDST.'

          call Modis_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      end if


c ... Return if additional QA checks reveal an error.
      if (error_flag) return


c-------------------------------------------------------------------------------
c retrieve L1B Earth View (EV), Sample Sizes,  and Uncertainty data
c-------------------------------------------------------------------------------

c.....Set start point and edge dimensions of desired L1B data slab.
      Start(DIM_NUM_XTRACK)     = 0
      Start(DIM_NUM_TRACK)      = (Scan_NUM-1) * LinesPerScan
      Start(DIM_NUM_CH)         = Ch_index - 1
      Slab_Size(DIM_NUM_XTRACK) = Data_size(DIM_NUM_XTRACK)
      Slab_Size(DIM_NUM_TRACK)  = Data_size(DIM_NUM_TRACK)
      Slab_Size(DIM_NUM_CH)     = 1


c.....Get Earth View (EV) and Uncertainty data slabs.  If aggregates, also
c.....get sample sizes.

c.....Earth View data ...
      rtn = gmar(L1B_MAPI_fhandle,EV_array_name,grpnm,start,Slab_Size,EV_data)

      if (rtn.ne.MAPIOK) then
          error_flag = .true.
          write(msg6,'(i6)') Band_NUM

          msgbuf =
     1    'M-API function GMAR detected error while retrieving L1B EV data '
     2    // 'for band ' // msg6
     3    // char(10) // 'Operator Action: Notify SDST.'

          call Modis_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      end if


c.....Uncertainty data ...
      rtn = gmar(L1B_MAPI_fhandle,UN_array_name,grpnm,start,Slab_Size,UN_data)

      if (rtn.ne.MAPIOK) then
          error_flag = .true.
          write(msg6,'(i6)') Band_NUM

          msgbuf =
     1    'M-API function GMAR detected error while retrieving L1B uncertainty '
     2    // 'data for band ' // msg6
     3    // char(10) // 'Operator Action: Notify SDST.'

          call Modis_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      end if


c.....Sample sizes
      if (Sa_array_name.ne.' ') then
          Sa_flag = .true.
          rtn = gmar(L1B_MAPI_fhandle,Sa_array_name,grpnm,start,Slab_Size,Sa_data)

          if (rtn.ne.MAPIOK) then
             error_flag = .true.
             write(msg6,'(i6)') Band_NUM

             msgbuf = 'M-API function GMAR detected error while retrieving L1B '
     1        // char(10) // 'aggregate data Sample Sizes for band ' // msg6
     2     // char(10) // 'Operator Action: Notify SDST.'

             call Modis_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
          end if
      end if


      if (error_flag) return


c-------------------------------------------------------------------------------
c Calibrate data
c-------------------------------------------------------------------------------

c.....Set EV data calibration slope and offset.
      scale      = EV_slope(Ch_index, Cal_Index, A_index, F_Index)
      add_offset = EV_offset(Ch_index, Cal_Index, A_index, F_Index)

c.....Loop over lines and pixels; check for missing or invalid data.
      do 90 j = 1, Data_size(DIM_NUM_TRACK)
      do 90 i = 1, Data_size(DIM_NUM_XTRACK)
         element = (j-1) * Data_size(DIM_NUM_XTRACK) + i

c........Data are missing or are other invalid types. 
         if (EV_data(element).lt.0) then
            V_flag(i,j) = L1B_INVALID

c           Return stored integer without application of SDS scale and offset
            Buf_EV(i,j) = EV_data(element) 
  
c           Return stored uncertainty index without application of SDS scale and offset.  
            Buf_Un(i,j) = ibits(Un_data(element),0,4)
            if (Un_data(element) .eq. L1B_MISSING) Buf_Un(i,j) = L1B_MISSING 

c           Return sample size of -99 if unusable data
            Buf_Sa(i,j) = -99

c........Data are valid; set validity flag to 0. Apply scale and offset to
c        Earth View (EV) data only. 
         else
             V_Flag(i,j) = L1B_VALID
             Buf_EV(i,j) = scale * (ibits(EV_data(element),0,15) - add_offset)
             Buf_Un(i,j) = ibits(Un_data(element),0,4)
             Buf_Sa(i,j) = 1
             if (Sa_flag) Buf_Sa(i,j) = Sa_data(element)
         end if

90    continue

      Last_L1B_fname(F_Index) = L1B_fname
      
      return
      end

c***********************************************************************

      Function Initial_Input_Checks(L1B_MAPI_fhandle,L1B_fname,
     1                              Band_NUM,Gain,Cal_Type,Resol)


      implicit none
      include 'mapi.inc'
      include 'PGS_MODIS_39500.f'
      include 'Read_L1B_Aerosol.inc'
      include 'extra.inc'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C
C  Perform quality assurance testing on the input arguments passed
C  to function Initial_Input_Checks.  These tests come mainly in the
C  form of out of bounds and valid range checking.
C
C  Function Initial_Input_Checks separately compares the input
C  arguments L1B_MAPI_fhandle, Gain, Cal_type and Resol against a
C  discrete set of permissible values.  If for each case no match is
C  found, an error message is reported and an error flag is set.
C  Out of bounds checking is performed on arguments L1B_fname and
C  Band_NUM.  If either satisfies the out of bounds criteria, an error
C  message is reported to the SDP toolkit LogStatus file.
C
C !INPUT PARAMETERS:
C
C  integer L1B_MAPI_fhandle  Array containing HDF SD (element 1) and
C                            VS (element 2) interface IDs, and
C                            the L1B file access mode (element 3).
C                            Additional elements carrying pointers to
C                            C language structures may also be passed,
C                            and are used internally by the MAPI
C                            library to maximize CPU efficiency.
C  character*(*) Cal_type    Calibration type applied to stored integers.
C                            Retrieve data as reflectance ("refl"),
C                            radiance ("rad") or
C                            corrected counts ("counts").
C  character*(*) Gain        Band radiometric gain factor,
C                            high ("H") or low ("L")
C  character*(*) L1B_fname   Variable containing full path name of
C                            L1B data file
C  integer   Band_NUM        MODIS spectral band number from 1 to 36
C  integer   Resol           Data spatial resolution factor relative to
C                            1 km.  Possible values are 16 (250 m),
C                            4 (500 m) and 1 (1 km).
C
C !OUTPUT PARAMETERS:      NONE
C
C !REVISION HISTORY:
C
C !TEAM-UNIQUE HEADER:
C
C    This software is developed by the MODIS Science Data Support
C    Team for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C !REFERENCES AND CREDITS:
C
C
C !DESIGN NOTES:
C
C  Return Value:     0 for success; -1 for FAIL
C
C  Externals:
C    Named Constants:
C        BAND_IS_*         (L1B_Reader_V2.1.inc)
C        CAL_IS_*          (L1B_Reader_V2.1.inc)
C        FAIL              (L1B_Reader_V2.1.inc)
C        FIRST_*_BAND      (L1B_Reader_V2.1.inc)
C        HIGH_GAIN         (L1B_Reader_V2.1.inc)
C        LAST_*_BAND       (L1B_Reader_V2.1.inc)
C        LOW_GAIN          (L1B_Reader_V2.1.inc)
C        P_SDID            (mapi.inc)
C        RESOL_IS_*        (L1B_Reader_V2.1.inc)
C        SUCCEED           (L1B_Reader_V2.1.inc)
C
C  Internals:
C    Subroutines:
C        MODIS_SMF_SetDynamicMsg
C
C !END
C---------------------------------------------------------------------

c PARAMETER declarations
      character*(*) FUNCNAME
      PARAMETER    (FUNCNAME = 'Initial_Input_Checks')


c Input argument declarations
      character*(*) Cal_Type,Gain,L1B_fname

      integer       Band_NUM,L1B_MAPI_fhandle(*),Resol


c Local variable and array declarations
      character*6   msg6
      character*135 msgbuf

      integer string_length,fbyte,lbyte
      integer Initial_Input_Checks


c-------------------------------------------------------------------------------
c Initialization
c-------------------------------------------------------------------------------

      string_length        = -1
      Initial_Input_Checks = SUCCEED


c-------------------------------------------------------------------------------
c Check for valid L1B HDF file ID (SD interface)
c-------------------------------------------------------------------------------

      If (L1B_MAPI_fhandle(P_SDID) .le. 0) Then
          Initial_Input_Checks = FAIL

          msgbuf = 'Invalid (< or = to 0) L1B HDF file ID (SD interface). '
     1             // char(10) // 'Operator Action:  Notify SDST.'

          call Modis_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf


c-------------------------------------------------------------------------------
c Check for file name out of bounds
c-------------------------------------------------------------------------------

      call CUT_NAME_L1B(L1B_fname,fbyte,lbyte)

      string_length = lbyte - fbyte

      If (string_length .le. 0) Then
          Initial_Input_Checks = FAIL
          msgbuf = 'L1B file name is a blank string, READ_L1B unable to continue.'
     1             // char(10) // 'Operator Action:  Notify SDST.'

          call Modis_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf


c-------------------------------------------------------------------------------
c check for valid, consistent "Resol" and "Band_NUM" input argument values.
c-------------------------------------------------------------------------------

      If (Resol .eq. RESOL_IS_250M) Then

         If ( (Band_NUM .lt. FIRST_250M_BAND) .or. (Band_NUM .gt. LAST_250M_BAND) ) Then
             Initial_Input_Checks = FAIL
             write(msg6,'(i6)') Band_NUM

             msgbuf = 'MODIS L1B 250M band set does not include band' // msg6
     1                // char(10) // 'Operator Action:  Notify SDST.'

             call Modis_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         EndIf

      ElseIf (Resol .eq. RESOL_IS_500M) Then

         If ( (Band_NUM .lt. FIRST_250M_BAND) .or. (Band_NUM .gt. LAST_500M_BAND) ) Then
            Initial_Input_Checks = FAIL
            write(msg6,'(i6)') Band_NUM

            msgbuf = 'MODIS L1B 500M band set does not include band' // msg6
     1               // char(10) // 'Operator Action:  Notify SDST.'

            call Modis_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         EndIf

      ElseIf (Resol .eq. RESOL_IS_1KM) Then

         If ( (Band_NUM .lt. FIRST_MODIS_BAND) .or. (Band_NUM .gt. LAST_MODIS_BAND) ) Then
            Initial_Input_Checks = FAIL
            write(msg6,'(i6)') Band_NUM

            msgbuf = 'MODIS L1B 1KM band set does not include band' // msg6
     1               // char(10) // 'Operator Action:  Notify SDST.'

            call Modis_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

         ElseIf ( (Band_NUM .eq. BAND_IS_13) .or. (Band_NUM .eq. BAND_IS_14) ) Then

            If ( (Gain .ne. HIGH_GAIN) .and. (Gain.ne.LOW_GAIN) ) Then
               Initial_Input_Checks = FAIL
               write(msg6,'(i6)') Band_NUM

               msgbuf = 'Gain is not one of "H" or "L" for MODIS L1B band' // msg6
     1                  // char(10) // 'Operator Action:  Notify SDST.'

               call Modis_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

            EndIf
         EndIf

      Else
         Initial_Input_Checks = FAIL
         write(msg6,'(i6)') Resol

         msgbuf = 'The resolution factor is not set to 16, 4, or 1: = ' // msg6
     1            // char(10) // 'Operator Action:  Notify SDST.'

         call Modis_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf


c-------------------------------------------------------------------------------
c Check for acceptable calibration type
c-------------------------------------------------------------------------------

      If (Cal_Type.eq.CAL_IS_REFL .or. Cal_Type.eq.CAL_IS_COUNT) Then

          If ( (Band_NUM .ge. First_Thermal_Band)   .and.
     1         (Band_NUM .le. Last_Thermal_Band )   .and.
     2         (Band_NUM .ne. Last_Solar_Band   ) ) Then

             Initial_Input_Checks = FAIL
             write(msg6,'(i6)') Band_NUM

             If (Cal_Type.eq.CAL_IS_REFL) Then
                msgbuf = 
     1          '"refl" is invalid calibration type for emissive band ' // msg6
     2          // char(10) // 'Operator Action:  Notify SDST.'
             Else
                msgbuf = 
     1          '"counts" is invalid calibration type for emissive band ' // msg6
     2          // char(10) // 'Operator Action:  Notify SDST.'
             EndIf

             Call Modis_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
          EndIf

      ElseIf ( (Cal_Type .ne. CAL_IS_RAD) .and. (Cal_Type .ne. CAL_IS_COUNT) ) Then
         Initial_Input_Checks = FAIL

         msgbuf = 'Calibration type is not one of "refl", "rad" or "counts"'
     1            // char(10) // 'Operator Action:  Notify SDST.'

         call Modis_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
      EndIf

      return
      end


c***********************************************************************

      SUBROUTINE CUT_NAME_L1B(string, fbyte, lbyte)
      IMPLICIT NONE

C----------------------------------------------------------------------
C !F77
C
C !Description:
C
C    Find the position (first and last bytes) of the file name within a
C    string buffer containing both path (file location) and file name.
C
C !INPUT PARAMETERS:
C
C    character string  A string variable of arbitrary length which
C                      shall consist of a unix path and file name.
C
C !OUTPUT PARAMETERS:
C
C    integer fbyte     The byte location of the first nonblank character
C                      of the input string.
C    integer lbyte     The byte location of the last nonblank character
C                      of the input string.
C
C !REVISION HISTORY:
c Revision 1.2  1997/02/25  18:55:14  vlin
c Changed from function to a subroutine
C
C !TEAM-UNIQUE HEADER:
C
C    This software was developed by the MODIS Science Data Support Team
C    (SDST) for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C !REFERENCES AND CREDITS
C
C    Written by Richard Hucek
C
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd, Seabrook MD 20706
C
C    rhucek@ltpmail.gsfc.nasa.gov
C
C !DESIGN NOTES:
C
C    There is no check to identify the unpredictable value of an
C    undefined string.  Consequently, users must take care to initialize
C    all string variables before passing them to CUT_NAME.
C
C !END
C-----------------------------------------------------------------------

      CHARACTER*(*) string
      INTEGER fbyte, lbyte, string_len

C Initialize variables
      string_len=len(string)
      fbyte=1
      lbyte = string_len

C Determine byte position of last non-blank, non-slash character.
C This is last character in file name
      DO WHILE ( (string(lbyte:lbyte).eq.' ').and.(lbyte.ge.1) )
         lbyte=lbyte-1
      END DO

C Determine byte position of last slash (/) in character string
      fbyte = lbyte - 1

      DO WHILE ( (fbyte.ge.1) .and. (string(fbyte:fbyte).ne.'/') )
         fbyte=fbyte-1
      END DO

      fbyte = fbyte + 1

      RETURN
      END

c***********************************************************************

      subroutine Get_Design_Indices( Resol,Band_NUM,Gain,Cal_Type,
     1                               F_Index,A_Index,Ch_Index,Cal_Index,
     2                               LinesPerScan )

      implicit none
      include 'Read_L1B_Aerosol.inc'
      include 'extra.inc'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C
C   Retrieve the L1B reader software internal design indices which are
C   used for referencing MODIS L1B sensor data by spatial resolution,
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
C  integer Band_NUM   MODIS spectral band number from 1 to 36
C  character Gain     Radiometric gain specification, high ("H") or
C                     low ("L").  It is ignored for all bands except
C                     13 and 14.
C  character Cal_Type Calibration format in one of three types:
C                     reflectance ("refl"), radiance ("rad"), or
C                     corrected counts ("counts").  Only "rad" is
C                     available for emissive band data (bands 20-36,
C                     but not 26). 
C
C
C !OUTPUT PARAMETERS:
C
C  integer  F_Index   L1B file index: either 1 (250 m), 2 (500 m) or
C                     3 (1 km).
C  integer  A_Index   L1B array identifier from 1 to 4 representing the
C                     250 m, 500 m, 1 km solar reflectance, and 1 km
C                     thermal emissive band data arrays, respectively.
C  integer  Ch_Index  Dimension index identifying the plane within the
C                     multidimensional L1B array containing the
C                     specified channel data.  Ch_index may have values
C                     from 1 to 16.
C  integer  Cal_Index Calibration type representing 1 (reflectance),
C                     2 (radiance), and 3 (corrected counts).
C  integer  LinesPerScan
C                     The number of detector lines in a scan of
C                     aggregated MODIS data.  May have values of
C                     40 (250 m), 20 (500 m) or 10 (1 km)
C
C !REVISION HISTORY:
C
C !TEAM-UNIQUE HEADER:
C
C    This software was developed by the MODIS Science Data Support Team
C    for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C !REFERENCES AND CREDITS:
C
C
C !DESIGN NOTES:
C
C  Externals:
C    Named Constants:
C        BAND_IS_*         (L1B_Reader_V2.1.inc)
C        CAL_IS_*          (L1B_Reader_V2.1.inc)
C        FIRST_*_BAND      (L1B_Reader_V2.1.inc)
C        HIGH_GAIN         (L1B_Reader_V2.1.inc)
C        LAST_*_BAND       (L1B_Reader_V2.1.inc)
C        LINES_PER_SCAN_*  (L1B_Reader_V2.1.inc)
C        RESOL_IS_*        (L1B_Reader_V2.1.inc)
C
C !END
C---------------------------------------------------------------------

C Declare function arguments
      character*(*) Cal_Type, Gain

      integer Band_NUM, LinesPerScan, Resol,
     1        F_Index, A_Index, Ch_Index, Cal_Index

C-------------------------------------------------------------------------------
C Determine design indices F_Index, A_index, Ch_index, and LinesPerScan
C-------------------------------------------------------------------------------

      If (Resol .eq. RESOL_IS_250M) Then
         F_Index      = 1
         A_Index      = 1
         Ch_Index     = Band_NUM
         LinesPerScan = LINES_PER_SCAN_250M

      ElseIf (Resol .eq. RESOL_IS_500M) Then
         F_Index      = 2
         LinesPerScan = LINES_PER_SCAN_500M

         If (Band_NUM .le. LAST_250M_BAND) Then
             A_Index  = 1
             Ch_Index = Band_NUM
         Else
             A_Index  = 2
             Ch_Index = Band_NUM - LAST_250M_BAND
         EndIf

      Else
         F_Index      = 3
         LinesPerScan = LINES_PER_SCAN_1KM

         If (Band_NUM .le. LAST_250M_BAND) Then
             A_Index  = 1
             Ch_Index = Band_NUM

         Else if (Band_NUM .le. LAST_500M_BAND) Then
             A_Index  = 2
             Ch_Index = Band_NUM - LAST_250M_BAND

         ElseIf (Band_NUM.lt.FIRST_THERMAL_BAND .OR. Band_NUM.eq.LAST_SOLAR_BAND) Then
             A_Index = 3

             If (Band_NUM .lt. BAND_IS_13) Then
                 Ch_Index = Band_NUM - 7
             ElseIf (Band_NUM .eq. BAND_IS_13) Then
                 Ch_Index = 6
                 If (Gain .eq. HIGH_GAIN) Ch_Index = 7
             ElseIf (Band_NUM .eq. BAND_IS_14) Then
                 Ch_Index = 8
                 If (Gain .eq. HIGH_GAIN) Ch_Index = 9
             Else
                 Ch_Index = Band_NUM - 5
                 If (Band_NUM .eq. LAST_SOLAR_BAND) Ch_Index = 15
             EndIf

         Else
             A_Index  = 4
             Ch_Index = Band_NUM - 19

             If (Band_NUM .gt. LAST_SOLAR_BAND) Ch_Index = Band_NUM - FIRST_THERMAL_BAND
         EndIf

      EndIf

c ... Assign calibration format design index, Cal_index

      If (Cal_Type .eq. CAL_IS_REFL) Then
          Cal_Index = 1
      ElseIf (Cal_Type .eq. CAL_IS_RAD) Then
          Cal_Index = 2
      Else
          Cal_Index = 3
      EndIf
      
      Return
      End

c***********************************************************************

      subroutine NewFileUpdate( L1B_MAPI_fhandle,
     1                          F_Index,
     2                          EV_offset,
     3                          EV_slope,
     4                          Dim_Size,
     5                          error_flag )

      implicit none

      include 'mapi.inc'
      include 'PGS_MODIS_39500.f'
      include 'Read_L1B_Aerosol.inc'
      include 'extra.inc'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C
C  Read and return the dimension sizes of all Earth View (EV) arrays in
C  a single L1B file.  In addition, the scale and offset factors needed
C  to transform the stored EV data to calibrated floating point numbers
C  are also provided (See Design Notes).
C
C !INPUT PARAMETERS:
C
C  integer  L1B_MAPI_fhandle  Array containing HDF SD (element 1) and VS
C                     (element 2) interface IDs, and the L1B file access
C                     mode (read, read/write, or create - element 3).
C                     Additional elements carrying pointers to C
C                     language structures are also passed, but they are
C                     used internally only by the MODIS Application
C                     Program Interface (M-API) library for maximizing
C                     CPU efficiency.
C  integer  F_Index   Spatial resolution index of L1B file.  Valid
C                     values are 1 (250 m), 2 (500 m) or 3 (1 km).
C
C !OUTPUT PARAMETERS:
C
C  integer Dim_Size   Array containing dimension sizes of L1B SDS
C                     arrays as a function of dimension index (the
C                     L1B arrays are 3 dimensional), L1B SDS array
C                     (up to 4), and spatial resolution (3).
C  real    EV_offset  Array containing Earth View data INTEGER to REAL
C                     conversion offset as a function of channel index
C                     (up to 16), calibration type (3), L1B SDS array
C                     (up to 4) and spatial resolution (3).
C  real    EV_slope   Array containing Earth View data INTEGER to REAL
C                     conversion slope as a function of channel index
C                     (up to 16), calibration type (3), L1B SDS array
C                     (up to 4) and spatial resolution (3).
C  logical error_flag Flag indicating error (TRUE) or no error (FALSE)
C                     detected.
C
C !REVISION HISTORY:
C Revision 1.4  1998/06/08  20:39:09  fhliang
C Update error messages with "Operator Action" strings;
C Modified error messages and description of input parameters.
C
C Revision 1.3  1997/11/26  18:19:13  rhucek
C Variable error_flag initialized to .false.
C
C Revision 1.2  1997/08/12  13:08:45  rhucek
C Changed type declaration of array L1B_MAPI_fhandle from fixed
C to assumed size.
C
C Changed range of DO 30 loop Cal_Index from
C Start_Cal_Type, No_Cal_Types to Start_Cal_Type, 2.
C This eliminates reading the calibration slope and
C offset factors for "corrected counts"
C
C Revision 1.1  1997/04/07  12:26:48  vlin
C Initial revision
C
C !TEAM-UNIQUE HEADER:
C
C    This software was developed by the MODIS Science Data Support
C    Team for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C !REFERENCES AND CREDITS:
C
C    Written by Vicky Lin          March 1997
C
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd, Seabrook MD 20706
C    rhucek@ltpmail.gsfc.nasa.gov
C    vlin@modis1.gsfc.nasa.gov
C
C !DESIGN NOTES:
C
C  A granule of a MODIS L1B product consists of 3 separate EV files
C  containing: 250 m spatial resolution data, 500 m data, and 1 km data,
C  respectively.  The EV observations in all files are stored as 2-byte
C  integers in three dimensional HDF arrays.  The most rapidly varying
C  dimension of these arrays is frame number, followed by scan line,
C  and then channel index.  NewFileUpdate returns the sizes of these
C  three dimensions as functions of array type (up to 4) and spatial
C  resolution (3).  The slope and offset values are returned in arrays
C  of dimensions 16 channels (max), 3 calibration modes, 4 L1B SDS
C  arrays (max) and 3 file types (250m, 500m, and 1km).
C
C  Externals
C  ---------
C   GMARDM, GMARIN               (mapi.inc, libmapi.a)
C   MAPIOK                       (mapi.inc)
C
C    Named Constants:
C     DIM_NUM_CH                 (L1B_Reader_V2.1.inc)
C     MAX_L1B_NAME_SIZE          (L1B_Reader_V2.1.inc)
C     MAX_NUM_ARRAYS             (L1B_Reader_V2.1.inc)
C     MAX_NUM_CHANNELS           (L1B_Reader_V2.1.inc)
C     MAX_NUM_DIMS               (L1B_Reader_V2.1.inc)
C     MAX_NUM_CAL_TYPES          (L1B_Reader_V2.1.inc)
C     NUM_L1B_FILES              (L1B_Reader_V2.1.inc)
C
C    Variables and Arrays
C     EV_Name                    (L1B_Reader_V2.1.inc)
C     Name_Of_Offset             (L1B_Reader_V2.1.inc)
C     Name_Of_Scale              (L1B_Reader_V2.1.inc)
C
C
C  Local Variables
C  ---------------
C   A_Index        Array type index:  1 (250m solar), 2 (500m solar),
C                  3 (1km solar), and 4 (1 km thermal).
C   A_Name         SDS array name
C   Attr_Name      A dimension attribute name
C   Cal_Index      Calibration type index: 1 reflectance, 2 radiance,
C                  and 3 corrected counts.
C   DIM_NUM_CH     channel dimension number
C   Dim_Name       SDS's dimension name
C   No_of_Arrays   Array containing number of EV arrays in L1B file as a
C                  function of spatial resolution: 250 m (1), 500 m (2),
C   End_Cal_Type   variable containing index of last calibration type
C                  for a specified L1B EV array.
C   Start_Cal_Type variable containing index of first calibration type
C                  for a specified L1B EV array.
C
C !END
C---------------------------------------------------------------------

C.....PARAMETER declarations
      character*(*)  FUNCNAME
      PARAMETER     (FUNCNAME = 'NewFileUpdate')


C.....Function argument declarations
      integer L1B_MAPI_fhandle(*),
     1        F_Index,
     2        Dim_Size(MAX_NUM_DIMS, MAX_NUM_ARRAYS, NUM_L1B_FILES)

      logical error_flag

      real EV_slope (MAX_NUM_CHANNELS, MAX_NUM_CAL_TYPES, MAX_NUM_ARRAYS, NUM_L1B_FILES),
     1     EV_offset(MAX_NUM_CHANNELS, MAX_NUM_CAL_TYPES, MAX_NUM_ARRAYS, NUM_L1B_FILES)


C.....Local variable and array declarations
      character*6    msg6,msg62,msg63
      character*35   data_type
      character*2048 msgbuf
      character*(MAX_L1B_NAME_SIZE) A_Name,Attr_Name,grpnm

      integer A_Index,Cal_Index,cal_rank,fbyte,lbyte,rtn,
     1        No_Channels,End_Cal_Type,Start_Cal_Type,
     2        No_of_Arrays(NUM_L1B_FILES)

      integer string_loc


c.....DATA statements
      data No_of_Arrays / 1, 2, 4 /


C-------------------------------------------------------------------------------
C Initialization
C-------------------------------------------------------------------------------

      grpnm          = ' '
      cal_rank       = 3
      error_flag     = .false.


C-------------------------------------------------------------------------------
C Retrieve dimension sizes and data conversion slope and offsets for all EV
C arrays in file.
C
C Note: there is 1 EV array in 250m file, 2 in 500m file and 4 in 1km file
C-------------------------------------------------------------------------------

      Do 30 A_Index = 1, No_of_Arrays(F_Index)
         A_Name = EV_Name(A_Index,F_Index)
         rtn    = string_loc(A_Name,fbyte,lbyte)

c........Get EV array dimension sizes
         rtn = GMARDM( L1B_MAPI_fhandle, A_Name, grpnm,
     1                 data_type, cal_rank, Dim_Size(1,A_Index,F_Index))

         If (rtn.ne.MAPIOK) Then
            write(msg6,'(i6)')  A_Index
            write(msg62,'(i6)') F_Index

            msgbuf = 'GMARDM detected error retrieving dimension sizes of L1B EV '
     1      // 'data array, ' // A_Name(fbyte:lbyte) // '.'
     2      // char(10) // 'Array index is ' // msg6 // '; spatial resolution index is ' // msg62
     3      // char(10) // 'Operator Action: Notify SDST.'

            call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

            error_flag = .true.
         EndIf


C-------------------------------------------------------------------------------
C Retrieve EV array slope and offsets for all calibration types applicable to
C array
C
C Note:  there are 3 calibration types (reflectance, radiance and counts) for
C        all EV arrays except the 1km thermal emissive data (A_Index=4) for
C        which there is just one, radiance. 
C-------------------------------------------------------------------------------

         Start_Cal_Type = 1
         End_Cal_Type = MAX_NUM_CAL_TYPES

c........Only the "radiance" calibration type is available for the thermal emissive array 
         If (A_Index .eq. 4) Start_Cal_Type = 2
         If (A_Index .eq. 4) End_Cal_Type = 2

         Do 40 Cal_Index = Start_Cal_Type, End_Cal_Type
            Attr_Name   = Name_Of_Offset(Cal_Index)
            data_type   = 'REAL*4'
            No_Channels = Dim_Size(DIM_NUM_CH,A_Index,F_Index)

c...........Get data conversion offsets
            rtn = GMARIN( L1B_MAPI_fhandle, A_Name, grpnm,
     1                    Attr_Name, data_type, No_Channels,
     2                    EV_offset(1,Cal_Index,A_Index,F_Index) )

            If (rtn.ne.MAPIOK) Then
               write(msg6,'(i6)')  A_Index
               write(msg62,'(i6)') F_Index
               write(msg63,'(i6)') Cal_Index

               msgbuf =
     1         'GMARIN detected error retrieving integer-to-real conversion '
     2         // 'offsets in array ' // A_Name(fbyte:lbyte)
     3         // char(10) // 'Array index is '// msg6 // '; spatial resolution index is '
     4         // msg62 // '; calibration index is ' // msg63
     5         // char(10) // 'Operator Action: Notify SDST.'

               call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

               error_flag = .true.
            EndIf


            Attr_Name = Name_Of_Scale(Cal_Index)

c.......... Get data conversion slopes
            rtn = GMARIN( L1B_MAPI_fhandle, A_Name, grpnm,
     1                    Attr_Name, data_type, No_Channels,
     2                    EV_slope(1,Cal_Index,A_Index,F_Index) )

            If (rtn.ne.MAPIOK) Then
               write(msg6,'(i6)')  A_Index
               write(msg62,'(i6)') F_Index
               write(msg63,'(i6)') Cal_Index

               msgbuf =
     1         'GMARIN detected error retrieving integer-to-real conversion '
     2         // char(10) // 'scale factors in array ' // A_Name(fbyte:lbyte) // char(10)
     3         // 'Array index is '// msg6 // '; spatial resolution index is ' // msg62
     4         // '; calibration index is ' // msg63
     5         // char(10) // 'Operator Action: Notify SDST.'

               call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

               error_flag = .true.
               return
            EndIf

40       Continue
30    Continue

      return
      end
