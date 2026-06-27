      INTEGER FUNCTION PARAM_GET_VALUE( PARAM_NUM, PARAM_LIST,
     &  PARAM_NAME, PARAM_VALUE )
           
C-------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C     Get the value string for a named parameter contained in an
C     array of parameter strings read by PARAM_READ_FILE.
C
C !INPUT PARAMETERS:
C     PARAM_NUM     Number of parameters in PARAM_LIST
C     PARAM_LIST    Array of parameter strings
C     PARAM_NAME    Name of parameter to extract
C
C !OUTPUT PARAMETERS:
C     PARAM_VALUE   String containing value for parameter PARAM_NAME
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

      INTEGER param_num
      CHARACTER*(*) param_list( param_num )
      CHARACTER*(*) param_name
      
c ... Output arguments

      CHARACTER*(*) param_value
      
c ... Local variables

      CHARACTER*255 name_string, curr_string, string
      INTEGER name_length, curr_length
      INTEGER i
      INTEGER start_pos, end_pos, exc_pos

c ... External functions

      INTEGER strpos, strlen
      EXTERNAL strpos, strlen
            
c ... Set return values

      param_get_value = -1
      param_value = ' '

c ... If parameter name is empty, return

      if ( len( param_name ) .eq. 0 ) return

c ... Get lowercase compressed version of parameter name

      name_string( 1 : len( name_string ) ) = ' '
      name_string( 1 : len( param_name ) ) = param_name
      call strlower( name_string )
      call strcompress( name_string, .TRUE., name_length )

c ... Loop through parameter list until parameter name is found

      do i = 1, param_num

c ...   Get lowercase compressed version of current parameter      
      
        curr_string( 1 : len( curr_string ) ) = ' '
        curr_string( 1 : len( param_list( i ) ) ) = param_list( i )
        call strlower( curr_string )
        call strcompress( curr_string, .TRUE., curr_length )

c ...   Check that current parameter is valid

        if ( strpos( curr_string, ':' ) .ge. 1 .and.
     &       curr_string( 1 : 1 ) .ne. ':' .and.           
     &       curr_string( 1 : 1 ) .ne. '!' ) then

c ...     Check if current parameter matches parameter name

          if ( curr_string( 1 : name_length ) .eq.
     &         name_string( 1 : name_length ) ) then

c ...       Get start and end positions of parameter value

            start_pos = strpos( param_list( i ), ':' ) + 1
            end_pos = strlen( param_list( i ) )
            exc_pos = strpos( param_list( i ), '!' )
            if ( exc_pos .ge. 1 ) end_pos = exc_pos - 1

c ...       Extract parameter value and return

            string( 1 : len( string ) ) = ' '
            string = param_list( i )
            param_value( 1 : len( param_value ) ) = ' '
            param_value = string( start_pos : end_pos )
            param_get_value = 0
            return

          endif
          
        endif

      end do
      
      END             
