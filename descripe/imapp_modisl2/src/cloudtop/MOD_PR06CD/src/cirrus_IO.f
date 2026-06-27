      SUBROUTINE MOD_cirrus_PutData(Scan_No,
     1                             Modfil_MOD06cd,
     3                             DataSize_XY,
     5                             BufSize_XTrack,
     6                             BufSize_Track,
     7                             cirrus_refl,
     8                             I_CONTRAIL)

      IMPLICIT NONE

      INCLUDE 'cirrus.inc'
      INCLUDE 'mapi.inc'
      INCLUDE 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION: 
C
C   Write a scan of MOD06 Cirrus parameters to MOD06 product file. 
C
C !INPUT PARAMETERS:   
C
C   INTEGER  Scan_NO        Instrument scan number relative to start of
C                           of granule. 
C   INTEGER  Modfil_MOD06cd Array containing HDF VS and SD ID numbers of
C                           MOD06CD product file.
C   INTEGER  DataSize_XY    Array specifying (spatial) dimension sizes 
C                           of data block within storage buffers.   
C   INTEGER  BufSize_XTRack
C                           Across track dimension size of data 
C                           buffers as allocated in calling routine.
C   INTEGER  BufSize_TRack  Along track dimension size of data buffers 
C                           as allocated in calling routine.
C   REAL     cirrus_refl    cirrus reflectance data 
C   INTEGER  I_CONTRAIL     cirrus reflectance flag
C
C !OUTPUT PARAMETERS:     NONE
C
C!Revision History:
C Modified by Liqun Ma    02/18/98
C Modified from MOD05_PutScanData.
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
      PARAMETER    (FUNCNAME = 'MOD_cirrus_PutData') 

      INTEGER    MAX_XTRACK,      MAX_TRACK
      PARAMETER (MAX_XTRACK=1500, MAX_TRACK=10)

      INTEGER    No_Lines_Per_Scan
      PARAMETER (No_Lines_Per_Scan = 10)

c-----Argument declarations
      INTEGER Scan_NO,Modfil_MOD06CD(MODFILLEN),
     1        DataSize_XY(2),
     2        BufSize_XTrack,BufSize_Track,
     3        I_CONTRAIL(BufSize_XTrack,BufSize_Track)

      REAL cirrus_refl(BufSize_XTrack,BufSize_Track)


c-----Local declarations
      BYTE Buf_Byte(MAX_XTRACK*MAX_TRACK)
      
      CHARACTER*25   msg25
      CHARACTER*1024 msgbuf

      INTEGER*2 Buf_Int2(MAX_XTRACK*MAX_TRACK)

      INTEGER fbyte,lbyte,i,ibyte,iframe,iline,rtn,
     1        DIM2(2),START2(2),
     2        String_Loc

C-----------------------------------------------------------------------
C Buffer, scale and insert a scan of NIR water vapor parameters into
C MOD06CD product file. 
C-----------------------------------------------------------------------

      write(msg25,'(i25)') Scan_No
      rtn = String_Loc(msg25,fbyte,lbyte)

      Do iline =  1, DataSize_XY(2)
      Do iframe = 1, DataSize_XY(1)
         i = (iline-1)*DataSize_XY(1) + iframe
         Buf_Int2(i) = NINT(cirrus_refl(iframe,iline)*SCALE1+OFFSET1)
         Buf_Byte(i) = I_CONTRAIL(iframe,iline)*SCALE2+OFFSET2
      End Do
      End Do 

      START2(1) = 0
      START2(2) = (Scan_No-1)*No_Lines_Per_Scan
      DIM2(1) = DataSize_XY(1)
      DIM2(2) = DataSize_XY(2)

C-----Insert scan of Cirrus Reflectance(1st SDS) into HDF array structure 
      If (PMAR(Modfil_MOD06cd,SDS1_N,' ',START2,DIM2,Buf_Int2) .NE.
     1    MAPIOK) Then 
         msgbuf = 'MAPI PMAR returned error writing Cirrus Reflectance'
     1       // 'array (SDS1).'
     2       // char(10) // 'Scan No is ' // msg25(fbyte:lbyte)
     3       // char(10) // 'Operator Action: check system disk' 
     4       // char(10) // 'resources; If adequate, contact SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,msgbuf,FUNCNAME)
      End If


C-----Insert  Contrail_Index_Flag data (2ond SDS) into HDF array structure 
      If (PMAR(MODFIL_MOD06cd,SDS2_N,' ',START2,DIM2,Buf_Byte) .NE.
     1    MAPIOK) Then 
      msgbuf='MAPI PMAR returned error writing Cirrus Reflectance Flag'
     1           // 'array (SDS2).'
     2           // char(10) // 'Scan No is ' // msg25(fbyte:lbyte)
     3           // char(10) // 'Operator Action: check system disk'
     4           // char(10) // 'resources; If adequate, contact SDST.'

         Call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,msgbuf,FUNCNAME)
      End If

 

      RETURN

      END
