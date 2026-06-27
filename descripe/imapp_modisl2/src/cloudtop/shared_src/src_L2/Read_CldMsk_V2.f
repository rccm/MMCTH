      SUBROUTINE Read_CldMsk(Modfil,Scan_No,
     &                       Dim1_CM,Dim1_QA,Dim2,Dim3,
     &                       DS_Dim1_CM,DS_Dim1_QA,DS_Dim2,DS_Dim3,
     &                       CM,QA,Error_Flag)

C----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:     
C
C   Read_CldMsk retrieves one scan of MODIS Cloud Mask and Cloud Mask
C   QA data from the MOD35 HDF product file.  It is assumed that the
C   spatial dimensions (number of frames and lines) of the HDF Cloud 
C   Mask and QA SDSs are equal (See Design Notes below). 
C
C !INPUT PARAMETERS:
C   INTEGER  Modfil      M-API file handle structure for HDF files
C
C   INTEGER  Scan_No     1-based instrument scan number 
C
C   INTEGER  Dim1_CM     Size of dimension 1 of Cloud Mask buffer as 
C                        dimensioned in calling program 
C
C   INTEGER  Dim1_QA     Size of dimension 1 of Cloud Mask QA buffer 
C                        as dimensioned in calling program 
C
C   INTEGER  Dim2/Dim3   Size of dimensions 2 and 3 of Cloud Mask and 
C                        Cloud Mask QA buffers as dimensioned in the  
C                        calling program.   
C
C !OUTPUT PARAMETERS:
C   BYTE     CM          Three dimensional (3-D) array for passing 
C                        cloud mask data.  Index 1 is byte number, 
C                        index 2 is (1-km) frame number, and index 3 is 
C                        relative (1-km) line number within scan. 
C 
C   BYTE     QA          Three dimensional (3-D) array for passing cloud 
C                        mask quality assurance data.  Index 1 is QA byte 
C                        number, index 2 is (1-km) frame number, and 
C                        index 3 is relative (1-km) line number within
C                        scan. 
C
C   INTEGER  DS_Dim1_CM  Size of retrieved Cloud Mask data block along
C                        dimension 1 of output buffer.  It is as large
C                        as the byte dimension of the HDF Cloud Mask 
C                        SDS
C
C   INTEGER  DS_Dim1_QA  Size of retrieved Cloud Mask QA data block 
C                        along dimension 1 of output buffer.  It is 
C                        as large as the byte dimension of the HDF Cloud
C                        Mask SDS 
C
C   INTEGER  DS_Dim2     Size of retrieved Cloud Mask and QA data 
C                        blocks along dimension 2 of output buffers.  
C                        It is as large as the frame (across track) 
C                        dimension of the HDF Cloud Mask and QA SDS 
C                        data arrays, which are assumed equal. 
C
C   INTEGER  DS_Dim3     Size of retrieved Cloud Mask and QA data 
C                        blocks along dimension 3 of output buffers.  
C                        It is equal to 10, the number of 1-km
C                        detector lines in a MODIS instrument scan. 
C
C   LOGICAL  Error_Flag  variable that is set to .TRUE. if an error is 
C                        detected. It is set to .FALSE. if no errors 
C                        are identified. 
C
C !REVISION HISTORY:
c Revision 1.4  1998/06/18  16:01:18  fhliang
c Update error messages with "Operator Action" strings.
c
c Revision 1.3  1997/12/15  17:04:45  fhliang
c changed 'mapic.inc' to 'mapi.inc', and added 'hdf.inc'.
c
c Revision 1.2  1997/06/13  19:37:09  rhucek
c Updated to return the data sizes of cloud mask and QA
c arrays retrieved from the MOD35 product file.
c
C !TEAM-UNIQUE HEADER:
C
C    This software is developed by the MODIS Science Data Support
C    Team for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C !REFERENCES AND CREDITS
C
C    Written by Vicky Lin         May 1997
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd, Seabrook MD 20706
C
C    vlin@ltpmail.gsfc.nasa.gov
C
C !DESIGN NOTES:
C
C    Subroutine Read_CldMsk checks the return status of all internal 
C    function calls.  If any call returns a fail indicator, Read_CldMsk 
C    reports an error message to the LogStatus file, and sets the output
C    argument Error_Flag to .TRUE..  Additional checks on Scan_No
C    and comparision of the dimension sizes of the 'Cloud_Mask' and 
C    'Quality_Assurance' arrays are made.  If incompatibilities are
C    found, Error_Flag = .TRUE. will be returned 
C
C    If all function calls are successful and no other discrepancies 
C    in the input parameters and dimensions size are found,  
C    Read_CldMsk runs to completion and returns Error_Flag = .FALSE.. 
C     
C
C  Externals:
C   Function:
C     GMAR                       (libmapi.a)
C     GMARDM                     (libmapi.a)
C
C   Named Constant:
C     P_SDID, P_ACCESS           (mapi.inc)
C     DFACC_READ                 (hdf.inc)
C     MAPIOK                     (mapi.inc)
C     MODIS_W_GENERIC            (MODIS_39500.f)
C
C  Internals:
C   Variables:
C    arrnam       SDS array name
C    grpnm        SDS group name
C    Edge(3)      Array specifying the number of data value to read.
C    Start(3)     Array specifying the starting location of data.
C    Max_QA_Bytes Maximum number of Cloud Mask QA bytes
C    Max_Frames   Maximum number of frames per scan line.
C    Max_Lines    Maximum number of 1-km lines per scan cube.
C    Rank         Number of dimensions in an array
C    MaxScan_No   Total Swath Number
C    count        A temporary buffer for data of the target array.
C    LinesPerScan Number of lines per scan cube
C    fbyte        Byte location of 1st nonblank character of the input
C                 string.
C    lbyte        Byte location of the last nonblank character of the
C                 input string.
C
C   Subroutines:
C     MODIS_SMF_SETDYNAMICMSG
C     STRING_LOC
C
C !END
C-----------------------------------------------------------------------

      IMPLICIT NONE

      INCLUDE 'mapi.inc'
      INCLUDE 'hdf.inc'
      INCLUDE 'PGS_MODIS_39500.f'

C Function argument declarations
      INTEGER Modfil(*),Scan_No, Dim1_CM, Dim1_QA, Dim2, Dim3,
     *        DS_Dim1_CM, DS_Dim1_QA, DS_Dim2, DS_Dim3 
 
      BYTE    CM(Dim1_CM,Dim2,Dim3), QA(Dim1_QA,Dim2,Dim3)

      LOGICAL Error_Flag

C Local variable declarations
      CHARACTER*4  msg4
      CHARACTER*13 data_type 
      CHARACTER*25 msg25
      CHARACTER*80 arrnm, grpnm
      CHARACTER*255 msgbuf
      CHARACTER*255 lmsgbuf

      CHARACTER*(*) NAME_CM_SDS, NAME_QA_SDS
      PARAMETER ( NAME_CM_SDS='Cloud_Mask', 
     *            NAME_QA_SDS='Quality_Assurance' )
      
      INTEGER LinesPerScan, Max_Frames, Max_Lines, Max_QA_Bytes 
      PARAMETER (LinesPerScan=10, Max_QA_Bytes=10, Max_Frames=1500, 
     *           Max_Lines=10)

      BYTE count(Max_Frames*Max_Lines*Max_QA_Bytes) 

      INTEGER Dim_Size_CM(3), Dim_Size_QA(3), Edge(3), fbyte, i, 
     *        indx, j, k, lbyte, MaxScan_No, Rank, rtn, Start(3)
      INTEGER STRING_LOC


C Initialization
      Error_Flag = .FALSE. 
      grpnm = ' '
      Rank = 3


C Check for valid file and access mode
      IF (Modfil(P_SDID).le.0 .or. Modfil(P_ACCESS).ne.DFACC_READ) THEN
         lmsgbuf = 'Invalid SD_ID or file access type detected'
     +   // char(10) // 'Operator Action:  Suspect system problem.  If a fault is '
     +   // char(10) // 'identified, correct and rerun PGE.  Otherwise, notify SDST.'     

         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,lmsgbuf,'Read_CldMsk')

         Error_Flag = .TRUE. 

         RETURN
      End If


C Retrieve dimensions of SDS array "Cloud_Mask" 
      arrnm = NAME_CM_SDS

      rtn = GMARDM(Modfil, arrnm, grpnm, data_type, Rank, Dim_Size_CM)

      IF (rtn .NE. MAPIOK) then
         lmsgbuf = 'MAPI routine GMARDM failed during access to Cloud_Mask array'
     +   // char(10) // 'Operator Action:  Suspect system problem.  If a fault is '
     +   // char(10) // 'identified, correct and rerun PGE.  Otherwise, notify SDST.'     

         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,lmsgbuf,'Read_CldMsk')

         Error_Flag = .TRUE. 
      ENDIF



C Retrieve dimensions of SDS array "Quality_Assurance"
      arrnm = NAME_QA_SDS 

      rtn = GMARDM(Modfil, arrnm, grpnm, data_type, Rank, Dim_Size_QA)

      IF (rtn .NE. MAPIOK) THEN
         lmsgbuf = 'MAPI routine GMARDM for Quality_Assurance failed'
     +   // char(10) // 'Operator Action:  Suspect system problem.  If a fault is '
     +   // char(10) // 'identified, correct and rerun PGE.  Otherwise, notify SDST.'     

         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,lmsgbuf,'Read_CldMsk')

         Error_Flag = .TRUE. 
      ENDIF

      If (Error_Flag) Return


      DS_Dim1_CM = Dim_Size_CM(3)   
      DS_Dim1_QA = Dim_Size_QA(1) 
      DS_Dim2 = Dim_Size_CM(1)
      DS_Dim3 = LinesPerScan 
  
C Compare line and frame dimension sizes of "Quality_Assurance"
C and "Cloud_Mask" arrays.  First compare frames, then lines.

      IF (Dim_Size_CM(1) .NE. Dim_Size_QA(2)) THEN 
         WRITE(msg25, '(2(2x, I6))') Dim_Size_CM(1), Dim_Size_QA(2)
         rtn = STRING_LOC(msg25,fbyte,lbyte) 
         msgbuf = 'Cloud_Mask and Quality_Assurance '
     +   // 'frame dimension sizes do not match: ' // msg25(fbyte:lbyte)
     +   // char(10) // 'Operator Action:  Suspect system problem.  If a fault is '
     +   // char(10) // 'identified, correct and rerun PGE.  Otherwise, notify SDST.'     

         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,'Read_CldMsk')

         Error_Flag = .TRUE.
      End If
      
      IF ( Dim_Size_CM(2) .NE. Dim_Size_QA(3) ) THEN 
         WRITE(msg25, '(2(2x, I6))') Dim_Size_CM(2), Dim_Size_QA(3)
         rtn = STRING_LOC(msg25,fbyte,lbyte) 
         msgbuf = 'Cloud_Mask and Quality_Assurance '
     +   // 'line dimension sizes do not match: ' // msg25(fbyte:lbyte)
     +   // char(10) // 'Operator Action:  Suspect system problem.  If a fault is '
     +   // char(10) // 'identified, correct and rerun PGE.  Otherwise, notify SDST.'     

         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,'Read_CldMsk')

         Error_Flag = .TRUE.
      End If


      IF (Error_Flag) Return


      
C Check for valid input value for variable "Scan_No" 
      MaxScan_No=Dim_Size_CM(2)/LinesPerScan

      IF (Scan_No .LT. 1 .OR. Scan_No .GT. MaxScan_No) THEN
         WRITE(msg4,'(i4)') MaxScan_No
         WRITE(msg25,'(i15)') Scan_No
         rtn = STRING_LOC(msg25,fbyte,lbyte) 
         msgbuf = 'Scan_No out of bounds.  It should be in range 1 -' // msg4
     +   // char(10) // 'Scan_No = ' // msg25(fbyte:lbyte)
     +   // char(10) // 'Operator Action:  Suspect system problem.  If a fault is '
     +   // char(10) // 'identified, correct and rerun PGE.  Otherwise, notify SDST.'     

         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,'Read_CldMsk')

         Error_Flag = .TRUE. 
      ENDIF



C Check for adequate output buffer size to store a scan of 
C Cloud_Mask data. 

      IF ( Dim1_CM .LT. Dim_Size_CM(3)) THEN
         WRITE(msg25, '(i15)') Dim1_CM
         rtn = STRING_LOC(msg25,fbyte,lbyte) 
         msgbuf = '1st dimension of output buffer too small'
     +   // ' to hold Cloud_Mask array'
     +   // char(10) // 'Dim1_CM = ' // msg25(fbyte:lbyte)
     +   // char(10) // 'Operator Action:  Suspect system problem.  If a fault is '
     +   // char(10) // 'identified, correct and rerun PGE.  Otherwise, notify SDST.'     

         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,'Read_CldMsk')

         Error_Flag = .TRUE. 
      END IF

      IF (Dim2 .LT. Dim_Size_CM(1)) THEN
         WRITE(msg25,'(i15)') Dim2
         rtn = STRING_LOC(msg25,fbyte,lbyte) 
         msgbuf = '2nd dimension of output buffer too small'
     +   // ' to hold Cloud_Mask array'
     +   // char(10) // 'Dim2 = ' // msg25(fbyte:lbyte)
     +   // char(10) // 'Operator Action:  Suspect system problem.  If a fault is '
     +   // char(10) // 'identified, correct and rerun PGE.  Otherwise, notify SDST.'     

         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,'Read_CldMsk')

         Error_Flag = .TRUE. 
      END IF

      IF (Dim3 .LT. LinesPerScan) THEN
         WRITE(msg25,'(i15)') Dim3
         rtn = STRING_LOC(msg25,fbyte,lbyte) 
         msgbuf = '3rd dimension of output buffer too small to hold Cloud_Mask array'
     +   // char(10) // 'Dim3 = ' // msg25(fbyte:lbyte)
     +   // char(10) // 'Operator Action:  Suspect system problem.  If a fault is '
     +   // char(10) // 'identified, correct and rerun PGE.  Otherwise, notify SDST.'     

         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,'Read_CldMsk')

         Error_Flag = .TRUE. 
      END IF

      IF (Error_Flag) RETURN 
     
      

C Retrieve SDS "Cloud_Mask" data
      arrnm = NAME_CM_SDS 
      Start(1) = 0
      Start(2) = (Scan_No-1)*LinesPerScan
      Start(3) = 0
      Edge(1) = Dim_Size_CM(1)
      Edge(2) = LinesPerScan
      Edge(3) = Dim_Size_CM(3)

      rtn = GMAR(Modfil, arrnm, grpnm, Start, Edge, count)

      IF (rtn .NE. MAPIOK) THEN
         write(msg25,'(3(2x,I6))') Start
         rtn = STRING_LOC(msg25,fbyte,lbyte) 
         msgbuf = 'MAPI routine GMAR failed during access to Cloud_Mask array'
     +   // char(10) // 'Read Dimension Offsets = ' // msg25(fbyte:lbyte)
     +   // char(10) // 'Operator Action:  Suspect system problem.  If a fault is '
     +   // char(10) // 'identified, correct and rerun PGE.  Otherwise, notify SDST.'     

         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,'Read_CldMsk')

         Error_Flag = .TRUE. 

         RETURN
      ENDIF

C Rebuffer 3-dimension cloud mask data with byte dimension
C varying most rapidly.

       Do 30 k=1,Edge(3)
       Do 30 j=1,Edge(2)
       Do 30 i=1,Edge(1)
          indx = (k-1)*Edge(1)*Edge(2) + (j-1)*Edge(1) + i
          CM(k,i,j) = count(indx)
   30 continue



C Check for adequate output buffer size to store a scan of  
C Quality_Assurance data. 

      IF (Dim1_QA .LT. Dim_Size_QA(1)) THEN
         WRITE(msg25,'(i15)') Dim1_QA
         rtn = STRING_LOC(msg25,fbyte,lbyte) 
         msgbuf = '1st dimension of output buffer too small'
     +   // ' to hold Quality_Assurance array'
     +   // char(10) // 'Dim1_QA = ' // msg25(fbyte:lbyte)
     +   // char(10) // 'Operator Action:  Suspect system problem.  If a fault is '
     +   // char(10) // 'identified, correct and rerun PGE.  Otherwise, notify SDST.'     

         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,'Read_CldMsk')

         Error_Flag = .TRUE. 
      END IF


      IF (Error_Flag) RETURN 



C Retrieve SDS "Quality_Assurance" data
      arrnm = NAME_QA_SDS
      Start(1) = 0
      Start(2) = 0
      Start(3) = (Scan_No-1)*LinesPerScan
      Edge(1) = Dim_Size_QA(1)
      Edge(2) = Dim_Size_QA(2)
      Edge(3) = LinesPerScan

      rtn = GMAR(Modfil, arrnm, grpnm, Start, Edge, count)

      IF (rtn .NE. MAPIOK) THEN
         write(msg25,'(3(2x,I6))') Start
         msgbuf = 'MAPI routine GMAR failed during access to Cloud Mask QA array'
     +   // char(10) // 'Read Dimension Offsets = ' // msg25
     +   // char(10) // 'Operator Action:  Suspect system problem.  If a fault is '
     +   // char(10) // 'identified, correct and rerun PGE.  Otherwise, notify SDST.'     

         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,'Read_CldMsk')

         Error_Flag = .TRUE. 

         RETURN
      ENDIF


C Move scan of QA data from work to output buffer. 
       Do 40 k = 1, Edge(3)
       Do 40 j = 1, Edge(2)
       Do 40 i = 1, Edge(1)
          indx = (k-1)*Edge(1)*Edge(2) + (j-1)*Edge(1) + i
          QA(i,j,k) = count(indx)
   40 Continue


      RETURN 
      END
