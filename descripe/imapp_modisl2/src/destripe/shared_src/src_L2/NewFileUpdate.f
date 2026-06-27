      subroutine NewFileUpdate(L1B_MAPI_fhandle,R_index,
     &                         EV_offset,EV_scale,Dim_size,error_flag)

      implicit none

      include 'mapi.inc'
      include 'PGS_MODIS_39500.f'
      include 'L1B_Reader.inc'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION: 
C
C  Read and return the dimension sizes of all Earth View (EV) arrays in
C  a single L1B file.  In addition, the scale and offset factors needed
C  to transform the stored EV data to calibrated floating point numbers
C  are also provided.
C  A granule of a MODIS L1B product is stored in 3 separate files 
C  containing: 250 m spatial resolution data, 500 m data, and 1 km data,
C  respectively.  The EV observations in all files are stored as 2-byte 
C  integers in three dimensional HDF arrays.  The most rapidly varying 
C  dimension of these arrays is frame number, followed by scan line, 
C  and then channel index.  NewFileUpdate returns the sizes of these 
C  three dimensions as functions of array type (up to 4) and spatial 
C  resolution (3).  The slope and offset values are returned in arrays 
C  of dimensions 16 channels (max), 3 calibration modes, 4 L1B SDS 
C  arrays and 3 spatial resolutions.
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
C  integer  R_index   Spatial resolution index of L1B sensor data 
C                     including aggregate data.  The set of valid values
C                     is 1 (250 m), 2 (500 m) or 3 (1 km).
C
C !OUTPUT PARAMETERS:
C
C  integer Dim_size   Array containing dimension sizes of L1B SDS 
C                     arrays as a function of dimension index (the 
C                     L1B arrays are 3 dimensional), L1B SDS array
C                     (up to 4), and spatial resolution (3).
C  real    EV_offset  Array containing Earth View data INTEGER to REAL 
C                     conversion offset as a function of channel index
C                     (up to 16), calibration type (3), L1B SDS array 
C                     (up to 4) and spatial resolution (3).
C  real    EV_scale   Array containing Earth View data INTEGER to REAL 
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
C Start_Cal_type, No_Cal_Types to Start_Cal_type, 2.
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
C   GMARDM, GMDMIN               (mapi.inc, libmapi.a)
C   MAPIOK                       (mapi.inc)
C
C  Data Dictionary
C  ---------------
C   A_index        Array type index:  1 (250m solar), 2 (500m solar),
C                  3 (1km solar), and 4 (1 km thermal).
C   A_name         SDS array name
C   Attr_name      A dimension attribute name
C   Cal_index      Calibration type index: 1 reflectance, 2 radiance,
C                  and 3 corrected counts.
C   CH_DIM_IDX     Channel dimension array index
C   CH_DIM_NO      HDF L1B EV array channel dimension index
C   Dim_Name       SDS's dimension name 
C   EV_name(4,3)   Names of L1B EV arrays for 4 array types and
C                  3 spatial resolutions
C   Name_of_Offset Name of dimension offset attribute for 3 different 
C                  calibration types: 1 reflectance, 2 radiance and 
C                  3 corrected counts.
C   Name_of_Scale  Name of dimension scale factor attribute for 3 
C                  different calibration types: 1 reflectance, 
C                  2 radiance and 3 corrected counts.
C   No_Channels    Size of array channel dimension
C   No_of_Arrays   Array containing number of EV arrays in L1B file as a
C                  function of spatial resolution: 250 m (1), 500 m (2),
C                                                  and 1 km (3).
C   Start_Cal_type Local index variable representing first allowable 
C                  calibration type for a specified L1B array.
C
C !END
C---------------------------------------------------------------------

C.....PARAMETER declarations
      character*(*)  FUNCNAME
      PARAMETER     (FUNCNAME = 'NewFileUpdate')

      integer    CH_DIM_IDX,     CH_DIM_NO
      PARAMETER (CH_DIM_IDX = 3, CH_DIM_NO = 0)


C.....Variable declarations
      character*6    msg6,msg62,msg63
      character*30   data_type,grpnm,Attr_Name,
     1               Name_of_Scale(3)  /'reflectance_scales',
     2                    'radiance_scales', 'corrected_counts_scales'/,
     3               Name_of_Offset(3) /'reflectance_offsets',
     4                    'radiance_offsets', 'corrected_counts_offsets'/
      character*70   A_Name
      character*2048 msgbuf

      integer cal_rank,fbyte,lbyte,rtn,
     1        A_index,R_index,Cal_Index,No_Channels, Start_Cal_type,
     2        No_of_Arrays(3),Dim_size(3,4,3),L1B_MAPI_fhandle(*)
      integer string_loc

      logical error_flag

      real EV_offset(16,3,4,3),EV_scale(16,3,4,3)


c.....DATA statements 
      data No_of_Arrays / 1, 2, 4 /

c ... Initialization
      grpnm = ' '
      cal_rank = 3      
      Start_Cal_type = 1
      error_flag = .false.
      
c ... Loop over the number of (EV) arrays in the L1B file.  Remember 
c     that the number of arrays in a L1B file varies with the spatial 
c     resolution of the data.  there is 1 EV array in the 250 m files,
c     2 in the 500 m files and 4 in the 1 km files.

      do 30 A_index = 1, No_of_Arrays(R_Index)
         A_Name = EV_Name(A_index,R_Index)
         rtn = string_loc(A_Name,fbyte,lbyte)

c ...... Get EV array dimension sizes
         rtn = GMARDM(L1B_MAPI_fhandle,A_Name,grpnm,data_type,cal_rank,
     1                Dim_size(1,A_index,R_Index))

         if (rtn.ne.MAPIOK) then
            write(msg6,'(i6)') A_index
            write(msg62,'(i6)') R_index

            msgbuf = 'GMARDM detected error retrieving dimension sizes of L1B EV ' 
     1      // 'data array, ' // A_Name(fbyte:lbyte) // '.'
     2      // char(10) // 'Array index is ' // msg6 // '; spatial resolution index is ' // msg62
     3      // char(10) // 'Operator Action:  Suspect system problem.  If a fault is '
     4      // char(10) // 'identified, correct and rerun PGE.  Otherwise, notify SDST.'

            call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

            error_flag = .true.
         end if


c ...... Loop over the possible calibration types of reflectance, radiance, 
c        and counts.  For the 1km thermal emissive data (A_index = 4) 
c        though, only two calibrations are possible, radiance and counts.

         A_Name = Dim_Name(A_index)
         rtn = string_loc(A_Name,fbyte,lbyte)

         if (A_index.eq.4) Start_Cal_type = 2

         do 40 Cal_Index = Start_Cal_type, 2 
            Attr_Name = Name_of_Offset(Cal_Index)
            data_type = 'REAL*4'
            No_Channels = Dim_Size(CH_DIM_IDX,A_index,R_index)

c...........Get integer to real data conversion offset
            rtn = GMDMIN(L1B_MAPI_fhandle,A_Name,grpnm,CH_DIM_NO,Attr_Name,data_type,
     &                   No_Channels,EV_offset(1,Cal_Index,A_index,R_Index))

            if (rtn.ne.MAPIOK) then
               write(msg6,'(i6)') A_index
               write(msg62,'(i6)') R_index
               write(msg63,'(i6)') Cal_index

               msgbuf = 'GMDMIN detected error retrieving integer-to-real conversion '
     1         // char(10) // 'offsets in array ' // A_Name(fbyte:lbyte) // char(10)
     2         // 'Array index is '// msg6 // '; spatial resolution index is ' // msg62 
     3         // '; calibration index is ' // msg63
     4         // char(10) // 'Operator Action:  Suspect system problem.  If a fault is '
     5         // char(10) // 'identified, correct and rerun PGE.  Otherwise, notify SDST.'

               call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

               error_flag = .true.
            end if


            Attr_Name = Name_of_Scale(Cal_Index)
               
c.......... Get integer to real data conversion scale factor
            rtn = GMDMIN(L1B_MAPI_fhandle,A_Name,grpnm,CH_DIM_NO,Attr_Name,data_type,
     1                   No_Channels,EV_scale(1,Cal_Index,A_index,R_Index))

            if (rtn.ne.MAPIOK) then
               write(msg6,'(i6)') A_Index
               write(msg62,'(i6)') R_index
               write(msg63,'(i6)') Cal_index

               msgbuf = 'GMDMIN detected error retrieving integer-to-real conversion '
     1         // char(10) // 'scale factors in array ' // A_Name(fbyte:lbyte) // char(10)
     2         // 'Array index is '// msg6 // '; spatial resolution index is ' // msg62 
     3         // '; calibration index is ' // msg63
     4         // char(10) // 'Operator Action:  Suspect system problem.  If a fault is '
     5         // char(10) // 'identified, correct and rerun PGE.  Otherwise, notify SDST.'

               call MODIS_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

               error_flag = .true.
               return
            end if

   40    continue
   30 continue

      return
      end
