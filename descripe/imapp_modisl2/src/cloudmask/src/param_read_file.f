      INTEGER FUNCTION PARAM_READ_FILE( PCF_NUM, PARAM_MAX,
     &  PARAM_NUM, PARAM_LIST )

C-------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C     Read a parameter file. A parameter file is an ASCII text file
C     containing 1 or more name/value pairs of the form
C
C     NAME : VALUE
C
C     A valid name/value pair must contain
C     - a name containing at least one character,
C     - a colon,
C     - at least one value. More than one value
C     may be defined by using commas to separate values, e.g.
C
C     ANGLES : 0.0, 15.0, 30.0, 45.0, 60.0, 75.0, 90.0
C
C     Comments are identified by the '!' character, which may occur
C     at the beginning of a line, or after a name/value pair, thus
C
C     ! This is a comment
C     PI : 3.1415    ! This is also a comment
C
C     are both valid comments. Blank lines are ignored. 
C
C !INPUT PARAMETERS:
C     PCF_NUM       PCF number for parameter file
C     PARAM_MAX     Maximum number of parameters
C                   (dimension of output array PARAM_LIST)
C
C !OUTPUT PARAMETERS:
C     PARAM_NUM     Number of parameters read from FILE
C     PARAM_LIST    Array of parameter strings read from FILE
C
C !REVISION HISTORY:
C
C !TEAM-UNIQUE HEADER:
C     Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
C
C !DESIGN NOTES:
C     Original version by Liam.Gumley@ssec.wisc.edu
C
C !END
C--------------------------------------------------------------------
      
      IMPLICIT NONE

c ... Input arguments

      INTEGER pcf_num
      INTEGER param_max
      
c ... Output arguments

      INTEGER param_num
      CHARACTER*(*) param_list( param_max )

c ... Local variables
      
      INTEGER status, lun
      INTEGER param_len
      CHARACTER*255 string
      INTEGER count
      integer record_length
      integer file_version

c ... Include file required for toolkit I/O

      INCLUDE 'PGS_IO.f'
      INCLUDE 'PGS_SMF.f'

c ... External functions

      integer pgs_io_gen_openf, pgs_io_gen_closef
      external pgs_io_gen_openf, pgs_io_gen_closef

c ... Set number of parameters found
      
      param_num = 0

c ... Open input file

      record_length = 1
      file_version = 1
      lun = -1
      status = pgs_io_gen_openf( pcf_num, PGSd_IO_Gen_RSeqFrm,
     &  record_length, lun, file_version )
      if ( status .ne. PGS_S_SUCCESS ) then
        param_read_file = -1
        return
      endif

c ... Get string length of parameter list

      param_len = len( param_list( 1 ) )

c ... Check that string length of parameter list does not exceed
c ... internal string length

      if ( param_len .gt. len( string ) ) then
        param_read_file = -2
        return
      endif

c ... Read all lines from the input file, checking that maximum
c ... parameter element number is not exceeded

      count = 0
20    continue
        read( lun, '(a)', end = 40 ) string
        count = count + 1
        if ( count .gt. param_max ) then
          param_read_file = -3
          return
        endif
        param_list( count ) = string( 1 : param_len )
      goto 20
40    continue

c ... Close input file

      status = pgs_io_gen_closef( lun )
      if ( status .ne. PGS_S_SUCCESS ) then
        param_read_file = -4
        return
      endif

c ... Set return values

      param_num = count
      param_read_file = 0
      
      END
