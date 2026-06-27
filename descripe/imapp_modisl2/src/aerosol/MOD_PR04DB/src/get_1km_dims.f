C*********************************************************************
       INTEGER FUNCTION get_1km_dims(dims, is_subset, p_offsets)

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
      CHARACTER*256 fname, attrname
      CHARACTER*25 istr
      CHARACTER*512 funcname, usrlog
      parameter (funcname='get_1km_dims')
      INTEGER      Rank,I,Dim_Size(3), fbyte1, lbyte1, dims(3)
      INTEGER      rtn, string_loc
      INTEGER      p1, l1, file_version, LRN_L1B_1km, LRN_L1B_1km_RA, LOGFLAG, m1
      INTEGER      prtn, pgs_pc_getreference
      INTEGER      No_byte,Fmax,Lmax
      PARAMETER    (No_byte=5,Fmax=150,Lmax=550)
      PARAMETER    (LRN_L1B_1km=700002, LRN_L1B_1km_RA=430001)
      BYTE         QAFlags(No_byte,Fmax,Lmax)
      BYTE         count(No_byte*Fmax*Lmax)
      INTEGER      Start(3),Edge(3),
     2             modfil(MODFILLEN), ierr

      LOGICAL      error_flag
      LOGICAL      is_subset
      INTEGER      p_offsets(2)
      
      external OPMFIL

      integer  num_args
      integer  FlagRA
      character FlagBuff*10
      integer  iargc


      num_args = iargc ( )
      
      if(num_args .eq. 1) then
         call getarg ( 1, FlagBuff )
         read (FlagBuff,*) FlagRA
      else
         !     This is the default value
         FlagRA = 0
      endif

C Initialization
      grpnm = ' '
      arrnm = 'EV_1KM_RefSB'
      get_1km_dims = -1
      Rank  =  3
      Start(1) = 0
      Start(2) = 0
      Start(3) = 0
      error_flag = .FALSE.

C Open file

*/  Retrieve the filename of the output swath file.
      file_version=1
      if( FlagRA .eq. 1) then
         prtn=pgs_pc_getreference(LRN_L1B_1km_RA,file_version,fname)
         WRITE (istr, '(I25)') LRN_L1B_1km_RA
      else
         prtn=pgs_pc_getreference(LRN_L1B_1km,file_version,fname)
         WRITE (istr, '(I25)') LRN_L1B_1km
      endif

      rtn = string_loc(istr, fbyte1, lbyte1)
      usrlog = "Retrieving filename for mod021km, LUN "
     +          //istr(fbyte1:lbyte1)// " - pgs_pc_getreference"
      CALL ckstatus_s(prtn,usrlog,funcname,LOGFLAG)

      IF (OPMFIL(fname, 'r', MODFIL) .eq. MAPIOK) THEN
      ELSE
        call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,
     &               'Failed opening MOD021KM file','get_1km_dims')
      ENDIF

C Check for valid file and band numbers
      IF (modfil(1).le.0.OR.modfil(3).ne.DFACC_READ) THEN
        call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     & 'Invalid SD_ID or invalid file access type','get_1km_dims')
        error_flag = .TRUE.
        return
      END IF


C Retrieve the rank, dimensions and data type of SDS data.
      IF (GMARDM(modfil, arrnm, grpnm, data_type, Rank, dims)
     &   .ne.MAPIOK) THEN
         call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     &   'GMARDM get dim_info failed','get_1km_dims')
         error_flag = .TRUE.
      END IF

C

C      dims(1) = Dim_Size(2)
C      dims(2) = Dim_Size(3)

C Get 1-km Cloud MASK data
C      Edge(1) = Dim_Size(1)
C      Edge(2) = Dim_Size(2)
C      Edge(3) = Dim_Size(3)

C Read HDF target array into 'QAFlags' buffer
C      if (GMAR(modfil, arrnm, grpnm, Start, Edge, count)
C    &   .ne.MAPIOK) Then
C         call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
C     &  'GMAR read arrdata failed','GetQAFlags')
C         error_flag = .TRUE.
C      End If

C Transfer flags from count array to QAflags

C       I = 0
C       Do 80 l1 = 1, Edge(3)
C       Do 80 p1 = 1, Edge(2)
C       Do 80 m1 = 1, Edge(1)

C          I = I + 1
C          QAflags(m1,p1,l1) = count(I)
C          if (m1 .eq. No_byte) QAflags(m1,p1,l1) = 0

C   80  continue

C
C Close up shop and go home

       IF (.NOT. error_flag) get_1km_dims = 0
       ierr = CLMFIL(MODFIL)
       Return
       END
