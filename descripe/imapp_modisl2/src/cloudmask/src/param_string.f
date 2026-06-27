      SUBROUTINE PARAM_STRING( PCF_NUM, NAME, STRING )

C-------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C     Extract the string corresponding to a parameter
C     name in a parameter file. See function PARAM_READ_FILE for a
C     description of the file format.
C
C !INPUT PARAMETERS:
C     PCF_NUM    PCF number for parameter file
C     NAME       Name of parameter to extract
C
C !OUTPUT PARAMETERS:
C     STRING     String containing the extracted value
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
C      integer pcf_num
C      character*80 string
C      
C      pcf_num = 600300
C      call param_string( pcf_num, 'rcs_id', string )
C      write(*,'(a)') string
C      
C      end
C
C !END
C--------------------------------------------------------------------

      IMPLICIT NONE
      
c ... Arguments

      integer pcf_num
      character*(*) name
      character*(*) string

c ... Local variables

      integer status
      integer param_max
      parameter( param_max = 10000 )
      character*255 param_list( param_max )
      integer param_num
      character*255 param_value
      integer ios
      character*12 pcf_text
      character*40 temp_string
            
c ... External variables
      
      integer param_read_file, param_get_value
      external param_read_file, param_get_value

c ... Read parameter file      

      status = param_read_file( pcf_num, param_max, 
     &  param_num, param_list )
      if ( status .ne. 0 ) then
        write( pcf_text, '(i12)' ) pcf_num
        call message( 'param_string',
     &    'Error reading parameter file PCF#' // pcf_text //
     &    ' [OPERATOR ACTION: Contact SDST]', status, 2 )
      endif

c ... Get parameter value

      status = param_get_value( param_num, param_list, name,
     &  param_value )
      if ( status .ne. 0 ) then
        write( temp_string, '(a)' ) name( 1 : len( name ) )
        call message( 'param_string',
     &    'Error getting parameter value ' // temp_string //
     &    ' [OPERATOR ACTION: Contact SDST]', status, 2 )
      endif
      
c ... Read parameter values into array

      read( param_value, '(a)', iostat = ios ) string
      if ( status .ne. 0 ) then
        write( temp_string, '(a)' ) name( 1 : len( name ) )
        call message( 'param_string',
     &    'Error reading string for parameter value ' // temp_string //
     &    ' [OPERATOR ACTION: Contact SDST]', status, 2 )
      endif
     
      END
