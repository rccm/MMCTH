module TimeHandling
  ! Module revision history
  !   $Log: TimeHandling.f90,v $
  !
  ! Revision 1.3  2000/02/29  Rich Hucek
  ! Changed implementation of GMT offset in function getExecutionTime 
  ! from "+" to "-".
  !
  ! Revision 1.2  1998/02/13  19:37:08  pincus
  ! Added function to report current time during execution, convert from toolkit
  ! internal time to CCSDS time format A; to convert between CCSDS time formats
  ! A and B; to get start and stop times of files and compare them to the
  ! time window for the current PCF.
  ! Allowed for leap second ambiguity in StringToTaiTime.
  ! Added toolkit-defined LUNs for PCF collection start and stop times as
  ! public parameters.
  !
  ! Revision 1.1  1997/10/29  22:58:23  pincus
  ! Initial revision

  ! General shared code
  use typeSizes,      only: EightByteRealKind
  use CharacterUtils, only: CharToInt
  ! Toolkit shared code
  use ErrorHandling,  only: WarningMessage, ErrorMessage, FatalErrorMessage
  use FileHandling,   only: getNumberOfFiles, readRuntimeParameter
  use Metadata,       only: getCoreMetadata  
  implicit none
  private
  
  ! Parameters defining the lengths of the CCSDS time formats. See, for example, page
  !   241 of the SDP User's Guide.
  integer, public, parameter :: FormatAStringLength = 27,    &
                                FormatBStringLength = 26,    &
                                ! These logical unit numbers are defined by the toolkit. 
                                StartTimeLUN        = 10258, &
                                StopTimeLUN         = 10259
  
  ! Public functions
  public :: StringToTaiTime, TaiTimeToString, TimeFormatAtoB,               &
            getPGETime, getExecutionTime, getFileBeginTime, getFileEndTime, &
            compareFileTimesToPGE

  ! Toolkit include files.
  include "PGS_PC.f"
  INCLUDE 'PGS_SMF.f'
  INCLUDE "PGS_TD_3.f"
  ! Interface definitions copied from the SDP toolkit user's guide, version 
  !   5.2. Putting the interface in the module makes it harder to make 
  !   mistakes when calling toolkit functions. 
  interface 
    integer function PGS_TD_UtcToTAI(AsciiUtcTime, secondsAfter1993)
      use typeSizes,     only: EightByteRealKind
      character (len = *),             intent( in) :: AsciiUtcTime ! Time code in first 27 characters.
      real (kind = EightByteRealKind), intent(out) :: secondsAfter1993
    end function PGS_TD_UtcToTAI
    
    integer function PGS_TD_TAIToUtc(secondsAfter1993, AsciiUtcTime)
      use typeSizes,     only: EightByteRealKind
      real (kind = EightByteRealKind), intent( in) :: secondsAfter1993
      character (len = *),             intent(out) :: AsciiUtcTime ! Time code in first 27 characters.
    end function PGS_TD_TAIToUtc
    
    integer function PGS_TD_AsciiTime_AtoB(AsciiUtcTimeA, AsciiUTCTimeB)
      ! This function changes from CCSDS Format A to B.
      character (len = *),             intent( in) :: AsciiUtcTimeA ! Time code in first 27 characters.
      character (len = *),             intent(out) :: AsciiUtcTimeB ! Time code in first 26 characters.
    end function PGS_TD_AsciiTime_AtoB
  end interface
  
contains
  ! -------------------------------------------------
  function StringToTaiTime(AsciiUtcTime)
    use typeSizes,     only: EightByteRealKind
    character (len = *), intent( in) :: AsciiUtcTime ! Time code in first 27 characters.
    real (kind = EightByteRealKind)  :: StringToTaiTime
    ! !F90
    !
    ! !Description:
    !    Converts a string containing UTC time in CCSDS Ascii Time Code A or B 
    !      to toolkit internal time.  Error trapping is included.
    !  
    ! !Input Parameters:
    !   AsciiUtcTime: The first 27 characters should contain UTC time in CCSDS Ascii Time Code 
    !     in A or B format.
    !
    ! !Output Parameters:
    !    Function value is a eight byte real value containing toolkit internal time (the 
    !      number of seconds since midnight on Jan 1, 1993). 
    ! 
    ! !Revision History:
    !    See module revision history at top of file. 
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
    !   This function is simply a wrapper for PGS_TD_UtcToTAI().
    !
    ! !END

    
    ! Local variables
    integer                          :: status
    character (len = 2048)           :: msgbuffer

    status = PGS_TD_UtcToTAI(trim(AsciiUtcTime) // "Z", StringToTaiTime)
    ! The Toolkit function can have several non-success outcomes that are still
    !   acceptable. 
    select case (status)
      case (PGS_S_SUCCESS)
        status = PGS_S_SUCCESS ! No operation
      case (PGSTD_W_PRED_LEAPS)
        call WarningMessage   ("StringToTaiTime", "Warning from toolkit function.")
      case (PGSTD_E_NO_LEAP_SECS)
        call ErrorMessage     ("StringToTaiTime", "Non-fatal error from toolkit function.")
      case default
        msgbuffer = "Conversion of string " // &
                    trim(AsciiUtcTime) // "Z" // " to time failed." &
                    // char(10) // "Operator Action:  Notify SDST."
        call FatalErrorMessage("StringToTaiTime", msgbuffer)
    end select
                           
  end function StringToTaiTime
  ! -------------------------------------------------
  function TaiTimeToString(UtcTimeString)
    use typeSizes,     only: EightByteRealKind
    real (kind = EightByteRealKind), intent( in)  :: UtcTimeString
    character (len = FormatAStringLength)         :: TaiTimeToString
    ! !F90
    !
    ! !Description:
    !    Converts toolkit internal time to CCSDS Ascii Time Code A format.  
    !      Error trapping is included.
    !  
    ! !Input Parameters:
    !    UtcTimeString: toolkit internal time (the number of seconds since midnight on Jan 1, 1993). 
    !
    ! !Output Parameters:
    !    None. Function result is the UTC time in CCSDS Ascii Time Code format A.
    ! 
    ! !Revision History:
    !    See module revision history at top of file. 
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
    !   This function is simply a wrapper for PGS_TD_UtcToTAI().
    !
    ! !END

    
    ! Local variables
    integer                          :: status
    character (len = 2048)           :: msgbuffer

    status = PGS_TD_TaiToUTC(UtcTimeString, TaiTimeToString)
    ! The Toolkit function can have several non-success outcomes that are still
    !   acceptable. 
    select case (status)
      case (PGS_S_SUCCESS)
        status = PGS_S_SUCCESS ! No operation
      case (PGSTD_W_PRED_LEAPS)
        call WarningMessage   ("TaiTimeToString", "Warning from toolkit function.")
      case (PGSTD_E_NO_LEAP_SECS)
        call ErrorMessage     ("TaiTimeToString", "Non-fatal error from toolkit function.")
      case default
        msgbuffer = "Conversion of time to string failed." &
                    // char(10) // "Operator Action:  Notify SDST."
        call FatalErrorMessage("TaiTimeToString", msgbuffer)
    end select
                           
  end function TaiTimeToString
  ! -------------------------------------------------
  function TimeFormatAtoB(FormatATimeString)

  use CharacterUtils, only: IntToChar

    character (len = *),      intent( in) :: FormatATimeString
    character (len = FormatBStringLength) :: TimeFormatAtoB

    ! !F90
    !
    ! !Description:
    !  Convert from CCSDS time string format A to format B.
    !
    ! !Input Parameters:
    !   FormatATimeString: string in CCSDS time string format A
    !
    ! !Output Parameters:
    !   Function result is a character string with the same time represented 
    !     in CCSDS format B.
    ! 
    ! !Revision History:
    !    See module revision history at top of file. 
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
    integer                              :: status
    character (len = 2048)               :: msgbuffer

    if(len_trim(FormatATimeString) > FormatAStringLength) then
      msgbuffer = "Non-fatal error: Input string " // &
                   trim(FormatATimeString) // &
                  " exceeds allowed length of " // &
                   trim(IntToChar(FormatAStringLength)) &
                 // char(10) // "Operator Action:  Notify SDST."
      call ErrorMessage("TimeFormatAtoB", msgbuffer)
      TimeFormatAtoB = ""
    else
      status = PGS_TD_AsciiTime_AtoB(FormatATimeString, TimeFormatAtoB)
    end if
    if (status /= 0) then
      msgbuffer = "Non-fatal error: Can't convert input string " // &
                   trim(FormatATimeString) // & 
                  " to CCSDS Format B." &
                   // char(10) // "Operator Action:  Notify SDST."
      call ErrorMessage("TimeFormatAtoB", msgbuffer)
      TimeFormatAtoB = ""
    end if
  end function 
  ! -------------------------------------------------
  function getPGETime(TimeLUN) result(secsAfter1993)
    integer,           intent( in) :: TimeLUN
    real(kind = EightByteRealKind) :: secsAfter1993
    ! !F90
    !
    ! !Description:
    !  Find the time in seconds after 1993 (TAI) associated with 
    !    a given LUN in the PCF.
    !
    ! !Input Parameters:
    !   TimeLUN: Logical Unit Number in the Process Control File containing the time. 
    !
    ! !Output Parameters:
    !   secondsAfter1993: Time in TAI (seconds after 00:00, Jan 1, 1993). 
    ! 
    ! !Revision History:
    !    See module revision history at top of file. 
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
    
    ! Local variables
    character (len = PGSd_PC_VALUE_LENGTH_MAX) :: stringValue
    integer                                    :: status
    
    call readRuntimeParameter(TimeLUN, stringValue)
    secsAfter1993 = StringToTaiTime(trim(stringValue) // "Z")
    
  end function getPGETime
  ! -------------------------------------------------
  function getExecutionTime()
    real(kind = EightByteRealKind) :: getExecutionTime
    ! !F90
    !
    ! !Description:
    !  Returns the time in seconds after 1993 (TAI time) at the time the function is called. 
    !
    ! !Input Parameters:
    !   none
    !
    ! !Output Parameters:
    !   Function result is the time in TAI (seconds after 00:00, Jan 1, 1993). 
    ! 
    ! !Revision History:
    !    See module revision history at top of file. 
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
    
    ! Local variables
    character (len = 10)  :: date, time
    integer, dimension(8) :: intValues
    
    call Date_And_Time(date = date, time = time, values = intValues)
    ! Construct CCSDS time format on the fly from results of system call.
    getExecutionTime = StringToTaiTime(                                                &
                       date(1:4) // "-" // date(5:6) // "-" // date(7:8 ) // "T" //    &
                       time(1:2) // ":" // time(3:4) // ":" // time(5:10) ) -          &
                       intValues(4) * 60 ! GMT Offset in minutes from Date_And_Time intrinsic
  end function getExecutionTime
  ! -----------------------------------------------------
  function getFileBeginTime(LogicalUnitNumber, version)
    integer, intent( in) :: LogicalUnitNumber, version
    real(kind = EightByteRealKind) :: getFileBeginTime
    ! !F90
    !
    ! !Description:
    !    Gets the starting time of a file in seconds after 1993 from file core metadata. The
    !      eight byte real value is computed from the RangeBeginningDate and RangeBeginningTime 
    !      fields.
    !  
    ! !Input Parameters:
    !   LogicalUnitNumber: the Logical Unit Number of the file in the Process Control File (PCF)
    !   version: the instance of the LUN to read
    !
    ! !Output Parameters:
    !   Function value is the starting time of the file in seconds after 00:00:00.0, Jan 1, 1993
    !   (Toolkit internal time)
    ! 
    ! !Revision History:
    !    See module revision history at top of file. 
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
    
    character(len = FormatAStringLength) :: BeginDateString, BeginTimeString
    
    call getCoreMetadata(fileID = LogicalUnitNumber, version = version, &
                         ECSparamName = "RangeBeginningDate", paramValue = BeginDateString)
    call getCoreMetadata(fileID = LogicalUnitNumber, version = version, &
                         ECSparamName = "RangeBeginningTime", paramValue = BeginTimeString)
    getFileBeginTime = StringToTaiTime(trim(BeginDateString) // "T" // trim(BeginTimeString) // "Z")
    
  end function getFileBeginTime
  ! -----------------------------------------------------
  function getFileEndTime(LogicalUnitNumber, version)
    integer, intent( in) :: LogicalUnitNumber, version
    real(kind = EightByteRealKind) :: getFileEndTime
    ! !F90
    !
    ! !Description:
    !    Gets the ending time of a file in seconds after 1993 from file core metadata. The
    !      eight byte real value is computed from the RangeEndingDate and RangeEndingTime 
    !      fields.
    !  
    ! !Input Parameters:
    !   LogicalUnitNumber: the Logical Unit Number of the file in the Process Control File (PCF)
    !   version: the instance of the LUN to read
    !
    ! !Output Parameters:
    !   Function value is the ending time of the file in seconds after 00:00:00.0, Jan 1, 1993
    !   (Toolkit internal time)
    ! 
    ! !Revision History:
    !    See module revision history at top of file. 
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
    
    ! Local variables
    character(len = FormatAStringLength) :: EndDateString,   EndTimeString

    call getCoreMetadata(fileID = LogicalUnitNumber, version = version, &
                         ECSparamName = "RangeEndingDate", paramValue = EndDateString)
    call getCoreMetadata(fileID = LogicalUnitNumber, version = version, &
                         ECSparamName = "RangeEndingTime", paramValue = EndTimeString)
    getFileEndTime = StringToTaiTime(trim(EndDateString) // "T" // trim(EndTimeString) // "Z")
    
  end function getFileEndTime
  ! -------------------------------------------------
  function compareFileTimesToPGE(InputProductLUN, numVersions)
    integer,                         intent( in) :: InputProductLUN, numVersions
    logical, dimension(numVersions)              :: compareFileTimesToPGE
    
    ! Local variables
    real (kind = EightByteRealKind)                      :: PGEStartTime, PGEStopTime
    real (kind = EightByteRealKind), &
                               dimension(numVersions)    :: FileBeginTimes, FileEndTimes
    integer                                              :: file

    character (len = 2048)                               :: msgbuffer

    ! !F90
    !
    ! !Description:
    !    Compares the range beginning and ending time for all the files on a LUN  
    !      with the collection start and end time of the PGE. 
    !  
    ! !Input Parameters:
    !   LogicalUnitNumber: the Logical Unit Number of the file in the Process Control File (PCF)
    !   numVersions: the number of instances on the LUN 
    !
    ! !Output Parameters:
    !   None. Function result is a logical vector of length numVersions. The vector is .true.
    !     if either the begin or end time of the file falls within the time window defined by the
    !     PGE; it is .false. otherwise. 
    ! 
    ! !Revision History:
    !    See module revision history at top of file. 
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

    ! Have we been passed the correct number of versions?

    if(getNumberOfFiles(inputProductLUN) /= numVersions) then
      msgbuffer = "Number of versions in PCF must agree with numVersions." &
     // char(10) // "Operator Action:  Check for valid PCF file.  If wrong or " &
     // char(10) // "corrupted, stage correct PCF and rerun PGE.  Otherwise, " &
     // char(10) // "notify SDST. "
      call FatalErrorMessage("compareFileTimesToPGE", msgbuffer)
    end if

    PGEStartTime = getPGETime(StartTimeLUN)
    PGEStopTime  = getPGETime( StopTimeLUN)
    
    ! Now check the starting and ending times for each file. 
    checkFileTimes: do file = 1, numVersions
      ! Get the times...
      FileBeginTimes(file) = getFileBeginTime(InputProductLUN, file)
      FileEndTimes  (file) = getFileEndTime  (InputProductLUN, file)
    end do checkFileTimes
    
    ! We accept any Level 2 granule which overlaps the processing time window at all. 
    compareFileTimesToPGE(:) = (FileBeginTimes(:) >= PGEStartTime .and.      &
                                FileBeginTimes(:) <  PGEStopTime)       .or. &
                               (FileEndTimes  (:) >= PGEStartTime .and.      &
                                FileEndTimes  (:) <  PGEStopTime)
  end function compareFileTimesToPGE
end module TimeHandling
