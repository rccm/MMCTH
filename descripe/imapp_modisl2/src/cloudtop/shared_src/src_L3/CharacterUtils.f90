module CharacterUtils
  ! Module Revision History:
  !   $Log: CharacterUtils.f90,v $
  ! Revision 1.3  1998/02/13  16:25:44  pincus
  ! Added function CharToReal, which returns a real value from a charater string.
  !
  ! Revision 1.2  1997/11/20  21:56:28  pincus
  ! Added function IntToChar.
  !
  ! Revision 1.1  1997/07/14  23:46:31  pincus
  ! Initial revision
  
  implicit none
  private
  public :: CharToInt, IntToChar, CharToReal
contains
  function CharToInt(inputString)
    character(len = *), intent( in) :: inputString
    integer                         :: CharToInt
    ! !F90
    !
    ! !Description:
    !   Reads a character string into an integer variable 
    !  
    ! !Input Parameters:
    !    inputString: character string. Should contain the character representation
    !      of a single value. Traling blanks are allowed. Example: " 30   "
    !
    ! !Output Parameters:
    !    CharToInt returns a single integer value. 
    ! 
    ! !Revision History:
    !   See module revision history in data section of module. 
    !
    ! !Team-Unique Header:
    !   Cloud Retrieval Group, NASA Goddard Space Flight Center
    !
    ! !References and Credits:
    !   Written by
    !    Robert Pincus
    !    Climate and Radiation Branch, Code 913
    !    NASA/GSFC
    !    Greenbelt MD 20771
    !    Robert.Pincus@gsfc.nasa.gov
    !
    ! !Design Notes:
    !    No internal error checking is done - if an incorrect string
    !      is passed in the entire program will fail with a Fortran runtime error. 
    !
    ! !END
    
    ! Local variables
    character(len = len_trim(inputString)) tempString

    ! Remove trailing blanks, which might be read as 0s on some processors
    ! Use list-directed read. 
    tempString = trim(inputString)
    read(tempString, *) CharToInt
  end function CharToInt
  ! --------------------------------------------------------------------------------
  function IntToChar(integerValue)
    integer, intent( in) :: integerValue
    character(len = 25)  :: IntToChar
    ! !F90
    !
    ! !Description:
    !   Reads the character representation of an integer.  
    !  
    ! !Input Parameters:
    !    integerValue: an single integer value. 
    !
    ! !Output Parameters:
    !    CharToInt returns a 25 character string containing the input argument 
    !       with trailing blanks.  
    ! 
    ! !Revision History:
    !   See module revision history in data section of module. 
    !
    ! !Team-Unique Header:
    !   Cloud Retrieval Group, NASA Goddard Space Flight Center
    !
    ! !References and Credits:
    !   Written by
    !    Robert Pincus
    !    Climate and Radiation Branch, Code 913
    !    NASA/GSFC
    !    Greenbelt MD 20771
    !    Robert.Pincus@gsfc.nasa.gov
    !
    ! !Design Notes:
    !    No internal error checking is done, although the string should be 
    !      big enough to contain a 64 bit integer.
    !
    ! !END

    write(IntToChar, *) integerValue
    IntToChar = AdjustL(IntToChar)
  end function IntToChar
  ! --------------------------------------------------------------------------------
  function CharToReal(inputString)
    character(len = *), intent( in) :: inputString
    real                            :: CharToReal
    ! !F90
    !
    ! !Description:
    !   Reads a character string into an integer variable 
    !  
    ! !Input Parameters:
    !    inputString: character string. Should contain the character representation
    !      of a single value. Traling blanks are allowed. Example: " 30   "
    !
    ! !Output Parameters:
    !    CharToInt returns a single integer value. 
    ! 
    ! !Revision History:
    !   See module revision history in data section of module. 
    !
    ! !Team-Unique Header:
    !   Cloud Retrieval Group, NASA Goddard Space Flight Center
    !
    ! !References and Credits:
    !   Written by
    !    Robert Pincus
    !    Climate and Radiation Branch, Code 913
    !    NASA/GSFC
    !    Greenbelt MD 20771
    !    Robert.Pincus@gsfc.nasa.gov
    !
    ! !Design Notes:
    !    No internal error checking is done - if an incorrect string
    !      is passed in the entire program will fail with a Fortran runtime error. 
    !
    ! !END
    
    ! Local variables
    character(len = len_trim(inputString)) tempString

    ! Remove trailing blanks, which might be read as 0s on some processors
    ! Use list-directed read. 
    tempString = trim(AdjustL(inputString))
    read(tempString, *) CharToReal
  end function CharToReal
  
end module CharacterUtils
