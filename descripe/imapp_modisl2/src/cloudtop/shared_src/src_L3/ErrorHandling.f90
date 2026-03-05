module ErrorHandling
  ! Module Revision History
  !   $Log: ErrorHandling.f90,v $
  ! Revision 1.7  1998/02/13  18:43:07  pincus
  ! Added prologs to hidden generic procedures underlying public generic
  ! procedures to comply with ECS requirements.
  !
  ! Revision 1.6  1997/11/20  21:55:11  pincus
  ! Changed WarningMessage to write error code MODIS_W_GENERIC, added
  ! ErrorMessage to write MODIS_E_GENERIC.
  !
  ! Revision 1.5  1997/11/11  18:08:52  pincus
  ! Added trailing space to messages - easier to read LogStatus.
  !
  ! Revision 1.4  1997/11/07  22:06:41  pincus
  ! Commented out unneeded include files.
  !
  ! Revision 1.3  1997/11/04  16:48:45  pincus
  ! Removed tabs to comply with ANSI standard.
  !
  ! Revision 1.2  1997/10/30  21:58:27  pincus
  ! Added internal Fortran 90 implementation of Modis_SMF_SetDynamicMsg.
  ! Changed comments to be more instructive.
  !
  ! Revision 1.1  1997/07/14  23:15:08  pincus
  ! Initial revision
  !
  implicit none
  private

  ! Toolkit include file - defines MODIS message codes
  include "PGS_MODIS_39500.f"

  ! This interface allows each of ModisSuccessMessageNone, ModisSuccessMessageReal, and
  !   ModisSuccessMessageInt to be called using the name SuccessMessage. The specific names 
  !   (i.e. are ModisSuccessMessageInt) are hidden via the PRIVATE statement; the procedure
  !   can ONLY be called using SuccessMessage. 
  interface SuccessMessage
    module procedure ModisSuccessMessageNone, ModisSuccessMessageReal, ModisSuccessMessageInt
  end interface
  ! !F90
  !
  ! !Description:
  !    This function is a simple wrapper for MODIS_SMF_SETDYNAMICMSG. It 
  !      exists to isolate the science code from the toolkit calls. 
  !    SuccessMessage is invoked to write a MODIS_S_SUCCESS message to LogStatus.
  !  
  ! !Input Parameters:
  !   procedure, message: Character strings containing the name of the invoking 
  !     function and the error message, respectively. 
  !   realValue, intValue: optional arguments, only one of which may be supplied. 
  !     realValue is of type real, intValue of type integer. 
  !     The values are appended to the error message string before printing
  ! !Output Parameters:
  !   none
  ! 
  ! !Revision History:
  !   See module revision history at top of file. 
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
  !
  ! !END

  interface WarningMessage
    module procedure ModisWarningMessageNone, ModisWarningMessageReal, ModisWarningMessageInt
  end interface
  ! !F90
  !
  ! !Description:
  !    This function is a simple wrapper for MODIS_SMF_SETDYNAMICMSG. It 
  !      exists to isolate the science code from the toolkit calls. 
  !    WarningMessage is invoked to write a MODIS_W_GENERIC message to LogStatus.
  !  
  ! !Input Parameters:
  !   procedure, message: Character strings containing the name of the invoking 
  !     function and the error message, respectively. 
  !   realValue, intValue: optional arguments, only one of which may be supplied. 
  !     realValue is of type real, intValue of type integer. 
  !     The values are appended to the error message string before printing
  ! !Output Parameters:
  !   none
  ! 
  ! !Revision History:
  !   See module revision history at top of file. 
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
  !
  ! !END

  interface ErrorMessage
    module procedure ModisErrorMessageNone, ModisErrorMessageReal, ModisErrorMessageInt
  end interface
  ! !F90
  !
  ! !Description:
  !    This function is a simple wrapper for MODIS_SMF_SETDYNAMICMSG. It 
  !      exists to isolate the science code from the toolkit calls. 
  !    ErrorMessage is invoked to write a MODIS_E_GENERIC message to LogStatus.
  !  
  ! !Input Parameters:
  !   procedure, message: Character strings containing the name of the invoking 
  !     function and the error message, respectively. 
  !   realValue, intValue: optional arguments, only one of which may be supplied. 
  !     realValue is of type real, intValue of type integer. 
  !     The values are appended to the error message string before printing
  ! !Output Parameters:
  !   none
  ! 
  ! !Revision History:
  !   See module revision history at top of file. 
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
  !
  ! !END

  interface FatalErrorMessage
    module procedure ModisFatalErrorMessageNone, ModisFatalErrorMessageReal, ModisFatalErrorMessageInt
  end interface
  ! !F90
  !
  ! !Description:
  !    This function is a simple wrapper for MODIS_SMF_SETDYNAMICMSG. It 
  !      exists to isolate the science code from the toolkit calls. 
  !    FatalErrorMessage is invoked to write a MODIS_F_GENERIC message to LogStatus.
  !      MODIS_SMF_SETDYNAMICMSG should immediately stop execution. 
  !  
  ! !Input Parameters:
  !   procedure, message: Character strings containing the name of the invoking 
  !     function and the error message, respectively. 
  !   realValue, intValue: optional arguments, only one of which may be supplied. 
  !     realValue is of type real, intValue of type integer. 
  !     The values are appended to the error message string before printing
  ! !Output Parameters:
  !   none
  ! 
  ! !Revision History:
  !   See module revision history at top of file. 
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
  !
  ! !END

  public :: SuccessMessage, WarningMessage, ErrorMessage, FatalErrorMessage
  
contains
  subroutine Modis_SMF_SetDynamicMsg(code, messageString, procedureName)
    integer,            intent( in) :: code
    character(len = *), intent( in) :: messageString, procedureName
    ! !F90
    !
    ! !Description:
    !    Use toolkit calls to write messages to LogStatus file. Test for fatal errors
    !      and abort execution if necessary. 
    !  
    ! !Input Parameters:
    !   code: Modis message code, as defined in message seed file "PGS_MODIS_39500.f"
    !   procedure, message: Character strings containing the name of the invoking 
    !     function and the error message, respectively. 
    !
    ! !Output Parameters:
    !   none
    ! 
    ! !Revision History:
    !   See module revision history at top of file. 
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
    !   This code is a Fortran 90 implementation of a function by the same name
    !     written in Fortran 77 by Vicky Lin. 
    !
    ! !Design Notes:
    !
    ! !END
    
    ! Toolkit include files
    INCLUDE 'PGS_SMF.f'
    ! Documentation claims I need these include files but nothing in them get used...
!     INCLUDE 'PGS_PC.f'
!     INCLUDE 'PGS_PC_9.f'
!     INCLUDE 'PGS_IO.f'
!     INCLUDE 'PGS_IO_1.f'
    
    ! Local variables
    character(len = 128) :: ModisMessage
    integer              :: status
    
    ! Adding an explicit interface (copied from the SDP toolkit user's guide)
    !   lets the compiler catch mistakes in calls to toolkit routines. 
    interface 
      integer function PGS_SMF_GetMsgByCode(code, messageString)
        integer,            intent( in) :: code
        character(len = *), intent(out) :: messageString
     end function PGS_SMF_GetMsgByCode
      
      integer function PGS_SMF_SetDynamicMsg(code, messageString, procedureName)
        integer,            intent( in) :: code
        character(len = *), intent( in) :: messageString, procedureName
      end function PGS_SMF_SetDynamicMsg
      
      integer function PGS_SMF_TestFatalLevel(code)
        integer, intent( in) :: code
      end function PGS_SMF_TestFatalLevel
    end interface
    
    ! Execution starts here.
    
    ! Get predefined message from seed file. Log error and abort if we can't get the message.
    if( PGS_SMF_GetMsgByCode(code, ModisMessage) /= PGS_S_Success) then
      status = PGS_SMF_SetDynamicMsg(MODIS_F_GENERIC,                                           &
                                     "Unable to invoke SDP TK function 'PGS_SMF_GetMsgByCode'", &
                                     "Modis_SMF_SetDynamicMessage")
      call exit(1)
    end if
    
    ! Write user and Modis messages to LogStatus 
    if( PGS_SMF_SetDynamicMsg(code,                                             & 
                              trim(ModisMessage) // " " // trim(messageString), &
                              trim(procedureName)) /= PGS_S_Success ) call exit(1)
    
    ! Test for fatal error.                           
    if( PGS_SMF_TestFatalLevel(code) == PGS_True ) call exit(1)
    
  end subroutine Modis_SMF_SetDynamicMsg
  
  ! -----------------------------------------------------------------
  !  Specific functions for SuccessMessage
  ! -----------------------------------------------------------------
  subroutine ModisSuccessMessageNone(procedureName, message)
    character(len = *), intent( in) :: procedureName, message
    ! !F90
    !
    ! !Description:
    !    SuccessMessage is invoked to write a MODIS_S_SUCCESS message to LogStatus.
    !  
    ! !Input Parameters:
    !   procedure, message: Character strings containing the name of the invoking 
    !     function and the error message, respectively. 
    ! !Output Parameters:
    !   none
    ! 
    ! !Revision History:
    !   See module revision history at top of file. 
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
    !
    ! !END
    
    call MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,  &
                                 trim(message),    &
                                 trim(procedureName))
  end subroutine ModisSuccessMessageNone
  ! ------------------------------------
  subroutine ModisSuccessMessageReal(procedureName, message, realValue)
    character(len = *), intent( in) :: procedureName, message
    real,               intent( in) :: realValue
    ! !F90
    !
    ! !Description:
    !    SuccessMessage is invoked to write a MODIS_S_SUCCESS message to LogStatus.
    !  
    ! !Input Parameters:
    !   procedure, message: Character strings containing the name of the invoking 
    !     function and the error message, respectively. 
    !   realValue: an realValue appended to the message string
    !
    ! !Output Parameters:
    !   none
    ! 
    ! !Revision History:
    !   See module revision history at top of file. 
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
    !
    ! !END

    ! Local variable
    character(len = 13) :: buf

    write(buf, '(f13.5)') realValue
    call MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,                   &
                                 trim(message) // " " // trim(buf), &
                                 trim(procedureName))
  end subroutine ModisSuccessMessageReal
  ! ------------------------------------
  subroutine ModisSuccessMessageInt(procedureName, message, intValue)
    character(len = *), intent( in) :: message, procedureName
    integer,            intent( in) :: intValue
    ! !F90
    !
    ! !Description:
    !    SuccessMessage is invoked to write a MODIS_S_SUCCESS message to LogStatus.
    !  
    ! !Input Parameters:
    !   procedure, message: Character strings containing the name of the invoking 
    !     function and the error message, respectively. 
    !   intValue: an integer appended to the message string
    !
    ! !Output Parameters:
    !   none
    ! 
    ! !Revision History:
    !   See module revision history at top of file. 
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
    !
    ! !END

    ! Local variable
    character(len = 9) :: buf

    write(buf, '(i9)') intValue
    call MODIS_SMF_SETDYNAMICMSG(MODIS_S_SUCCESS,                   &
                                 trim(message) // " " // trim(buf), &
                                 trim(procedureName))
  end subroutine ModisSuccessMessageInt
  
  ! -----------------------------------------------------------------
  !  Specific functions for WarningMessage
  ! -----------------------------------------------------------------
  subroutine ModisWarningMessageNone(procedureName, message)
    character(len = *), intent( in) :: procedureName, message
    ! !F90
    !
    ! !Description:
    !    WarningMessage is invoked to write a MODIS_W_GENERIC message to LogStatus.
    !  
    ! !Input Parameters:
    !   procedure, message: Character strings containing the name of the invoking 
    !     function and the error message, respectively. 
    !   
    ! !Output Parameters:
    !   none
    ! 
    ! !Revision History:
    !   See module revision history at top of file. 
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
    !
    ! !END

    call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,  &
                                 trim(message),    &
                                 trim(procedureName))
  end subroutine ModisWarningMessageNone
  ! ------------------------------------
  subroutine ModisWarningMessageReal(procedureName, message, realValue)
    character(len = *), intent( in) :: procedureName, message
    real,               intent( in) :: realValue
    ! !F90
    !
    ! !Description:
    !    WarningMessage is invoked to write a MODIS_W_GENERIC message to LogStatus.
    !  
    ! !Input Parameters:
    !   procedure, message: Character strings containing the name of the invoking 
    !     function and the error message, respectively. 
    !   realValue: an realValue appended to the message string
    !   
    ! !Output Parameters:
    !   none
    ! 
    ! !Revision History:
    !   See module revision history at top of file. 
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
    !
    ! !END

    ! Local variable
    character(len = 13) :: buf

    write(buf, '(f13.5)') realValue
    call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,                   &
                                 trim(message) // " " // trim(buf), &
                                 trim(procedureName))
  end subroutine ModisWarningMessageReal
  ! ------------------------------------
  subroutine ModisWarningMessageInt(procedureName, message, intValue)
    character(len = *), intent( in) :: message, procedureName
    integer,            intent( in) :: intValue
    ! !F90
    !
    ! !Description:
    !    WarningMessage is invoked to write a MODIS_W_GENERIC message to LogStatus.
    !  
    ! !Input Parameters:
    !   procedure, message: Character strings containing the name of the invoking 
    !     function and the error message, respectively. 
    !   intValue: an integer appended to the message string
    !   
    ! !Output Parameters:
    !   none
    ! 
    ! !Revision History:
    !   See module revision history at top of file. 
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
    !
    ! !END

    ! Local variable
    character(len = 9) :: buf

    write(buf, '(i9)') intValue
    call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,                   &
                                 trim(message) // " " // trim(buf), &
                                 trim(procedureName))
  end subroutine ModisWarningMessageInt
  
  ! -----------------------------------------------------------------
  !  Specific functions for ErrorMessage
  ! -----------------------------------------------------------------
  subroutine ModisErrorMessageNone(procedureName, message)
    character(len = *), intent( in) :: procedureName, message
    ! !F90
    !
    ! !Description:
    !    ErrorMessage is invoked to write a MODIS_E_GENERIC message to LogStatus.
    !  
    ! !Input Parameters:
    !   procedure, message: Character strings containing the name of the invoking 
    !     function and the error message, respectively. 
    ! !Output Parameters:
    !   none
    ! 
    ! !Revision History:
    !   See module revision history at top of file. 
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
    !
    ! !END

    call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,  &
                                 trim(message),    &
                                 trim(procedureName))
  end subroutine ModisErrorMessageNone
  ! ------------------------------------
  subroutine ModisErrorMessageReal(procedureName, message, realValue)
    character(len = *), intent( in) :: procedureName, message
    real,               intent( in) :: realValue
    ! !F90
    !
    ! !Description:
    !    ErrorMessage is invoked to write a MODIS_E_GENERIC message to LogStatus.
    !  
    ! !Input Parameters:
    !   procedure, message: Character strings containing the name of the invoking 
    !     function and the error message, respectively. 
    !   realValue: an realValue appended to the message string
    !
    ! !Output Parameters:
    !   none
    ! 
    ! !Revision History:
    !   See module revision history at top of file. 
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
    !
    ! !END

    ! Local variable
    character(len = 13) :: buf

    write(buf, '(f13.5)') realValue
    call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,                   &
                                 trim(message) // " " // trim(buf), &
                                 trim(procedureName))
  end subroutine ModisErrorMessageReal
  ! ------------------------------------
  subroutine ModisErrorMessageInt(procedureName, message, intValue)
    character(len = *), intent( in) :: message, procedureName
    integer,            intent( in) :: intValue
    ! !F90
    !
    ! !Description:
    !    ErrorMessage is invoked to write a MODIS_E_GENERIC message to LogStatus.
    !  
    ! !Input Parameters:
    !   procedure, message: Character strings containing the name of the invoking 
    !     function and the error message, respectively. 
    !   intValue: an integer appended to the message string
    !
    ! !Output Parameters:
    !   none
    ! 
    ! !Revision History:
    !   See module revision history at top of file. 
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
    !
    ! !END

    ! Local variable
    character(len = 9) :: buf

    write(buf, '(i9)') intValue
    call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,                   &
                                 trim(message) // " " // trim(buf), &
                                 trim(procedureName))
  end subroutine ModisErrorMessageInt
  
  ! -----------------------------------------------------------------
  !  Specific functions for FatalErrorMessage
  ! -----------------------------------------------------------------
  subroutine ModisFatalErrorMessageNone(procedureName, message)
    character(len = *), intent( in) :: procedureName, message
    ! !F90
    !
    ! !Description:
    !    FatalErrorMessage is invoked to write a MODIS_F_GENERIC message to LogStatus.
    !      MODIS_SMF_SETDYNAMICMSG should immediately stop execution. 
    !  
    ! !Input Parameters:
    !   procedure, message: Character strings containing the name of the invoking 
    !     function and the error message, respectively. 
    !
    ! !Output Parameters:
    !   none
    ! 
    ! !Revision History:
    !   See module revision history at top of file. 
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
    !
    ! !END

    call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,  &
                                 trim(message),    &
                                 trim(procedureName))
  end subroutine ModisFatalErrorMessageNone
  ! ------------------------------------
  subroutine ModisFatalErrorMessageReal(procedureName, message, realValue)
    character(len = *), intent( in) :: procedureName, message
    real,               intent( in) :: realValue
    ! !F90
    !
    ! !Description:
    !    FatalErrorMessage is invoked to write a MODIS_F_GENERIC message to LogStatus.
    !      MODIS_SMF_SETDYNAMICMSG should immediately stop execution. 
    !  
    ! !Input Parameters:
    !   procedure, message: Character strings containing the name of the invoking 
    !     function and the error message, respectively. 
    !   realValue: an realValue appended to the message string
    !
    ! !Output Parameters:
    !   none
    ! 
    ! !Revision History:
    !   See module revision history at top of file. 
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
    !
    ! !END

    ! Local variable
    character(len = 13) :: buf

    write(buf, '(f13.5)') realValue
    call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,                   &
                                 trim(message) // " " // trim(buf), &
                                 trim(procedureName))
  end subroutine ModisFatalErrorMessageReal
  ! ------------------------------------
  subroutine ModisFatalErrorMessageInt(procedureName, message, intValue)
    character(len = *), intent( in) :: message, procedureName
    integer,            intent( in) :: intValue
    ! !F90
    !
    ! !Description:
    !    FatalErrorMessage is invoked to write a MODIS_F_GENERIC message to LogStatus.
    !      MODIS_SMF_SETDYNAMICMSG should immediately stop execution. 
    !  
    ! !Input Parameters:
    !   procedure, message: Character strings containing the name of the invoking 
    !     function and the error message, respectively. 
    !   intValue: an integer appended to the message string
    !
    ! !Output Parameters:
    !   none
    ! 
    ! !Revision History:
    !   See module revision history at top of file. 
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
    !
    ! !END

    ! Local variable
    character(len = 9) :: buf

    write(buf, '(i9)') intValue
    call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,                   &
                                 trim(message) // " " // trim(buf), &
                                 trim(procedureName))
  end subroutine ModisFatalErrorMessageInt
  ! -----------------------------------------------------------------
  
end module ErrorHandling
