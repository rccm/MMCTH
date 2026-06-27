      SUBROUTINE PARAM_REAL( PCF_NUM, NAME, NUMBER, ARRAY )

C-------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C     Extract the real value(s) corresponding to a parameter
C     name in a parameter file. See function PARAM_READ_FILE for a
C     description of the file format.
C
C !INPUT PARAMETERS:
C     PCF_NUM    PCF number for parameter file
C     NAME       Name of parameter to extract
C     NUMBER     Number of values to extract for this parameter
C                (Maximum number of values in this version is 20).
C
C !OUTPUT PARAMETERS:
C     ARRAY      Array containing NUMBER extracted values
C
C !REVISION HISTORY:
C
C !TEAM-UNIQUE HEADER:
C     Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
C
C !DESIGN NOTES:
C     Original version by Liam.Gumley@ssec.wisc.edu
C
C Example:
C
C      implicit none
C
C      integer pcf_num, number, i
C      real array(4)
C      
C      pcf_num = 600300
C      number = 4
C      call param_real( pcf_num, 'pds4_11', number, array )
C      write(*,*) ( array( i ), i = 1, number )
C      
C      end
C
C !END
C--------------------------------------------------------------------

      IMPLICIT NONE
      
c ... Arguments

      integer pcf_num
      character*(*) name
      integer number
      real array( number )

c ... Local variables

      integer status
      integer param_max
      parameter( param_max = 10000 )
      character*255 param_list( param_max )
      integer param_num
      character*255 param_value
      integer ios
      integer i                  
      character*12 pcf_text
      character*40 string
      
c ... External variables
      
      integer param_read_file, param_get_value
      external param_read_file, param_get_value

c ... Read parameter file      

      status = param_read_file( pcf_num, param_max, 
     &  param_num, param_list )
      if ( status .ne. 0 ) then
        write( pcf_text, '(i12)' ) pcf_num
        call message( 'param_real',
     &    'Error reading parameter file PCF#' // pcf_text //
     &    ' [OPERATOR ACTION: Contact SDST]', status, 2 )
      endif

c ... Get parameter value

      status = param_get_value( param_num, param_list, name,
     &  param_value )
      if ( status .ne. 0 ) then
        write( string, '(a)' ) name( 1 : len( name ) )
        call message( 'param_real',
     &    'Error getting parameter value ' // string //
     &    ' [OPERATOR ACTION: Contact SDST]', status, 2 )
      endif

c ... Check that requested number of values does not exceed the
c ... number allowed by the FORMAT in the next READ statement

      if ( number .gt. 20 ) then
        write( string, '(a)' ) name( 1 : len( name ) )
        call message( 'param_real',
     &    'Too many values requested for parameter value ' // string //
     &    ' [OPERATOR ACTION: Contact SDST]', ios, 2 )
      endif
      
c ... Read parameter values into array

      read( param_value, '(20f20.10)', iostat=ios )
     &  ( array( i ), i = 1, number )
      if ( ios .ne. 0 ) then
        write( string, '(a)' ) name( 1 : len( name ) )
        call message( 'param_real',
     &    'Error reading data for parameter value ' // string //
     &    ' [OPERATOR ACTION: Contact SDST]', ios, 2 )
      endif
      
      END
