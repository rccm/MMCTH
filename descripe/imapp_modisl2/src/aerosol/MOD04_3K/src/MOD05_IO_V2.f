      SUBROUTINE MOD05_PutScanData(Scan_No,
     1                             Modfil_MOD05,
     2                             No_QA_Bytes,
     3                             DataSize_XY,
     4                             BufSize_QA,
     5                             BufSize_XTrack,
     6                             BufSize_Track,
     7                             VAPVRT,
     8                             CldMsk_QA,
     9                             QA_NIR)

      IMPLICIT NONE

      INCLUDE 'mod05.inc'
      INCLUDE 'mapi.inc'
      INCLUDE 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION: 
C
C   Write a scan of NIR water vapor parameters to MOD05 product file. 
C
C !INPUT PARAMETERS:   
C
C   INTEGER  Scan_NO       Instrument scan number relative to start of
C                          of granule. 
C   INTEGER  Modfil_MOD05  Array containing HDF VS and SD ID numbers of
C                          MOD05 product file.
C   INTEGER  DataSize_XY   Array specifying (spatial) dimension sizes 
C                          of data block within storage buffers.   
C   INTEGER  BufSize_XTRack
C                          Across track dimension size of data 
C                          buffers as allocated in calling routine.
C   INTEGER  BufSize_TRack Along track dimension size of data buffers 
C                          as allocated in calling routine.
C   REAL     VAPVRT        NIR total column WV data 
C   INTEGER  Cloud         Cloud Mask QA data 
C   INTEGER  QA_NIR        NIR Quality Assurance data  
C
C !OUTPUT PARAMETERS:     NONE
C
C!Revision History:
C
C!TEAM-UNIQUE HEADER: 
C
C   This software was developed by the MODIS Science Data Support Team
C   (SDST) for the National Aeronautics and Space Administration, 
C   Goddard Space Flight Center, under contract NAS5-32373.
C
C !REFERENCES and CREDITS:
C
C   Written by Richard Hucek
C
C   Research and Data Systems Corporation
C   SAIC/GSC MODIS Science Data Support Office
C   7501 Forbes Blvd, Seabrook MD 20706
C
C   rhucek@ltpmail.gsfc.nasa.gov
C
C !DESIGN NOTES:
C
C !END
C----------------------------------------------------------------------

c-----PARAMETER definitions
      CHARACTER*(*) FUNCNAME
      PARAMETER    (FUNCNAME = 'MOD05_IO_V2.f') 

      INTEGER     MAX_QA_BYTES
      PARAMETER ( MAX_QA_BYTES=10 )

      INTEGER    MAX_XTRACK,      MAX_TRACK
      PARAMETER (MAX_XTRACK=1500, MAX_TRACK=10)

      INTEGER    No_Lines_Per_Scan
      PARAMETER (No_Lines_Per_Scan = 10)

c-----Argument declarations
      INTEGER Scan_NO,Modfil_MOD05(MODFILLEN),
     1        No_QA_Bytes,DataSize_XY(2),
     2        BufSize_QA,BufSize_XTrack,BufSize_Track,
     3        CldMsk_QA(BufSize_XTrack,BufSize_Track),
     4        QA_NIR(No_QA_Bytes,BufSize_XTrack,BufSize_Track)

      REAL VAPVRT(BufSize_XTrack,BufSize_Track)


c-----Local declarations
      BYTE Buf_Byte(MAX_QA_BYTES*MAX_XTRACK*MAX_TRACK)
      
      CHARACTER*25   msg25
      CHARACTER*1024 msgbuf

      INTEGER*2 Buf_Int2(MAX_XTRACK*MAX_TRACK)

      INTEGER fbyte,lbyte,i,ibyte,iframe,iline,rtn,
     1        DIM2(2),DIM3(3),START2(2),START3(3),
     2        String_Loc

C-----------------------------------------------------------------------
C Buffer, scale and insert a scan of NIR water vapor parameters into
C MOD05 product file. 
C-----------------------------------------------------------------------

      write(msg25,'(i25)') Scan_No
      rtn = String_Loc(msg25,fbyte,lbyte)

      Do iline =  1, DataSize_XY(2)
      Do iframe = 1, DataSize_XY(1)
         i = (iline-1)*DataSize_XY(1) + iframe
         Buf_Int2(i) = NINT(VAPVRT(iframe,iline)*SCALE1+OFFSET1)
         Buf_Byte(i) = CldMsk_QA(iframe,iline)*SCALE2+OFFSET2
      End Do
      End Do 

      START2(1) = 0
      START2(2) = (Scan_No-1)*No_Lines_Per_Scan
      DIM2(1) = DataSize_XY(1)
      DIM2(2) = DataSize_XY(2)

C-----Insert scan of Water Vapor (1st SDS) into HDF array structure 
      If (PMAR(MODFIL_MOD05,SDS1_N,' ',START2,DIM2,Buf_Int2) .NE.
     1    MAPIOK) Then 
         msgbuf = 'MAPI PMAR returned error writing Water Vapor '
     1           // 'array (SDS1).'
     2           // char(10) // 'Scan No is ' // msg25(fbyte:lbyte)
     3   // char(10) // 'Operator Action:  Check system resources/environment, '
     4   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     5   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
         Call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,msgbuf,FUNCNAME)
      End If


C-----Insert Cloud Mask QA data (2ond SDS) into HDF array structure 
      If (PMAR(MODFIL_MOD05,SDS2_N,' ',START2,DIM2,Buf_Byte) .NE.
     1    MAPIOK) Then 
         msgbuf = 'MAPI PMAR returned error writing Cloud Mask '
     1           // 'QA data (SDS2).'
     2           // char(10) // 'Scan No is ' // msg25(fbyte:lbyte)
     3   // char(10) // 'Operator Action:  Check system resources/environment, '
     4   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     5   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
         Call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,msgbuf,FUNCNAME)
      End If

 
c-----Reuse Buf_Byte data buffer for NIR Water Vapor QA data
      Do iline =  1, DataSize_XY(2)
      Do iframe = 1, DataSize_XY(1)
      Do ibyte =  1, No_QA_Bytes 
         i = (iline-1)*DataSize_XY(1) + (iframe-1)*No_QA_Bytes + ibyte 
         Buf_Byte(i) = QA_NIR(ibyte,iframe,iline)*SCALE3 + OFFSET3
      End Do
      End Do 
      END Do

      START3(1) = 0
      START3(2) = 0 
      START3(3) = (Scan_NO-1)*No_Lines_Per_Scan
      DIM3(1) = No_QA_Bytes 
      DIM3(2) = DataSize_XY(1)
      DIM3(3) = DataSize_XY(2)

C-----Insert NIR QA data (3rd SDS) into HDF array structure 

      If (PMAR(MODFIL_MOD05,SDS3_N,' ',START3,DIM3,Buf_Byte) .NE. 
     1    MAPIOK) Then
         msgbuf = 'MAPI PMAR returned error writing NIR '
     1           // 'QA data (SDS3).'
     2           // char(10) // 'Scan No is ' // msg25(fbyte:lbyte)
     3   // char(10) // 'Operator Action:  Check system resources/environment, '
     4   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     5   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
         Call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,msgbuf,FUNCNAME)
      End If

      RETURN

      END
