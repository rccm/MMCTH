C*********************************************************************
       INTEGER FUNCTION GetQAFlags(QAFlags,QAFlags_O,DTaot,dims)

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
C    INTEGER    Buf_Size1/2  The sizes of 250-m CldMsk buffer.
C
C !OUTPUT PARAMETERS:
C
C    INTEGER    CldMsk(x,y)  Buffer storing Cloud Mask (0/cloudy,
C                            1/clear, -1 not determined.)
C    INTEGER    QAFlags(x/4,y/4)  Buffer containing land/sea
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
C        Dim_Size(2)  Array specifying the size of hdf SDS data array.
C        Edge(2)      Array specifying the number of data value to read.
C        Start(2)     Array specifying the starting location of data.
C        Fmax         Maximum frame number per scan line.
C        Lmax         Maximum line number per scan cube.
C        Rank         The number of dimensions in an array
C        l,l1,l4,     Line and pixel counters p,p1,p4
C        count(15000) A temporary buffer for data of the target array.
C
C !END
C-----------------------------------------------------------------------

C Declarations
      CHARACTER*80 arrnm,grpnm,data_type
      CHARACTER*256 fname
      CHARACTER*25 istr
      CHARACTER*512 funcname, usrlog
      parameter (funcname='GetQAFlags')
      INTEGER      Rank,I,Dim_Size(3), fbyte1, lbyte1, dims(3)
      INTEGER      rtn, string_loc
      INTEGER      p1, l1, file_version, LUN_sw, LOGFLAG, m1
      INTEGER      prtn, pgs_pc_getreference
      INTEGER      No_byte,Fmax,Lmax, No_byte_O
      PARAMETER    (No_byte=6,Fmax=150,Lmax=550,No_byte_O=5)
      PARAMETER    (LUN_sw=405000)
      BYTE         QAFlags(No_byte,Fmax,Lmax), QAFlags_O(No_byte_O,Fmax,Lmax)
      BYTE         count(No_byte*Fmax*Lmax)
      INTEGER      Dim_Size2d(2), Start2d(2), Edge2d(2)
      INTEGER*2    DTaot(Fmax,Lmax)
      INTEGER*2    count2d(Fmax*Lmax)
      INTEGER      Start(3),Edge(3),
     2             modfil(MODFILLEN), ierr

      LOGICAL      error_flag

      external OPMFIL

C Initialization
      grpnm = ' '
      arrnm = 'Quality_Assurance_Land'
      GetQAFlags = -1
      Rank  =  3
      Start(1) = 0
      Start(2) = 0
      Start(3) = 0
      error_flag = .FALSE.

C Open file

*/  Retrieve the filename of the output swath file.
      file_version=1
      prtn=pgs_pc_getreference(LUN_sw,file_version,fname)

      WRITE (istr, '(I25)') LUN_sw
      rtn = string_loc(istr, fbyte1, lbyte1)
      usrlog = "Retrieving filename for mod04, LUN "
     +          //istr(fbyte1:lbyte1)// " - pgs_pc_getreference"
      CALL ckstatus_s(prtn,usrlog,funcname,LOGFLAG)

      IF (OPMFIL(fname, 'r', MODFIL) .eq. MAPIOK) THEN
      ELSE
        call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,
     &               'Failed opening MOD04 file','GetQAFlags')
      ENDIF

C Check for valid file and band numbers
      IF (modfil(1).le.0.OR.modfil(3).ne.DFACC_READ) THEN
        call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     & 'Invalid SD_ID or invalid file access type','GetQAFlags')
        error_flag = .TRUE.
        return
      END IF


C Retrieve the rank, dimensions and data type of SDS data.
      IF (GMARDM(modfil, arrnm, grpnm, data_type, Rank, Dim_Size)
     &   .ne.MAPIOK) THEN
         call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     &   'GMARDM get dim_info failed','GetQAFlags')
         error_flag = .TRUE.
      END IF
C

      dims(1) = Dim_Size(2)
      dims(2) = Dim_Size(3)

C Get 1-km Cloud MASK data
      Edge(1) = Dim_Size(1)
      Edge(2) = Dim_Size(2)
      Edge(3) = Dim_Size(3)

C Read HDF target array into 'QAFlags' buffer
      if (GMAR(modfil, arrnm, grpnm, Start, Edge, count)
     &   .ne.MAPIOK) Then
         call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     &  'GMAR read arrdata failed','GetQAFlags')
         error_flag = .TRUE.
      End If

C Transfer flags from count array to QAflags

       I = 0
       Do 80 l1 = 1, Edge(3)
       Do 80 p1 = 1, Edge(2)
       Do 80 m1 = 1, Edge(1)

          I = I + 1
          QAflags(m1,p1,l1) = count(I)
          if (m1 .eq. No_byte) QAflags(m1,p1,l1) = 0

   80  continue
      
      arrnm = 'Quality_Assurance_Ocean'
C Retrieve the rank, dimensions and data type of SDS data.
      IF (GMARDM(modfil, arrnm, grpnm, data_type, Rank, Dim_Size)
     &   .ne.MAPIOK) THEN
         call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     &   'GMARDM get dim_info failed','GetQAFlags')
         error_flag = .TRUE.
      END IF
			
			Edge(1) = Dim_Size(1)
      Edge(2) = Dim_Size(2)
      Edge(3) = Dim_Size(3)
      
      if (GMAR(modfil, arrnm, grpnm, Start, Edge, count)
     &   .ne.MAPIOK) Then
         call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     &  'GMAR read arrdata failed','GetQAFlags')
         error_flag = .TRUE.
      End If
C Transfer flags from count array to QAflags

       I = 0
       Do 90 l1 = 1, Edge(3)
       Do 90 p1 = 1, Edge(2)
       Do 90 m1 = 1, Edge(1)

          I = I + 1
          QAflags_O(m1,p1,l1) = count(I)
          if (m1 .eq. No_byte_O) QAflags_O(m1,p1,l1) = 0

   90  continue

C [2] Dark Target AOT
C Re-initialization

      grpnm = ' '
      arrnm = 'Optical_Depth_Land_And_Ocean'

      GetQAFlags = -1
      Rank  =  2
      Start2d(1) = 0
      Start2d(2) = 0
      error_flag = .FALSE.

C Retrieve the rank, dimensions and data type of SDS data.
      IF (GMARDM(modfil, arrnm, grpnm, data_type, Rank, Dim_Size2d)
     &   .ne.MAPIOK) THEN
         call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     &   'GMARDM get dim_info failed','GetQAFlags')
         error_flag = .TRUE.
      END IF

C Get Dark Target AOT
      Edge2d(1) = Dim_Size2d(1)
      Edge2d(2) = Dim_Size2d(2)
      
C      print *, 'Start2d,Edge2d',Start2d,Edge2d
      
C Read HDF target array into 'DTaot' buffer
      if (GMAR(modfil, arrnm, grpnm, Start2d, Edge2d, count2d)
     &   .ne.MAPIOK) Then
         call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     &  'GMAR read arrdata failed','GetQAFlags')
         error_flag = .TRUE.
      End If
      IF (.NOT. error_flag) GetQAFlags = 0

C Transfer flags from count array to DTaot - have to do it this way because 
C don't know exact dimensions beforehand
C        print *, 'Edge(3),Edge(2),Edge(1)',Edge(3),Edge(2),Edge(1)
C        print *, 'Fmax,Lmax,Edge2d(2),Edge2d(1)',Fmax,Lmax,Edge2d(2),Edge2d(1)
       
       I = 0
       Do 85 l1 = 1, Edge2d(2)
       Do 85 p1 = 1, Edge2d(1)
C          print *, 'p1,l1,count2d(I)',p1,l1,count2d(I)
          I = I + 1
          DTaot(p1,l1) = count2d(I)

   85  continue
C
C Close up shop and go home

       IF (.NOT. error_flag) GetQAFlags = 0
       ierr = CLMFIL(MODFIL)

       Return
       END
