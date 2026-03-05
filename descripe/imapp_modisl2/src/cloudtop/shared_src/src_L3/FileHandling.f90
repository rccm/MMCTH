module FileHandling
  ! Module Revision History:
  !   $Log: FileHandling.f90,v $
  ! Revision 1.4  1998/02/13  18:24:18  pincus
  ! Modified getNumberOfFiles to return 0 if there are no files on a LUN.
  ! Added constants MaxFileNameLength, MaxUniversalRefNameLength; added
  ! prologs to specific procedures underlying generic procedures to comply
  ! with ECS requirements.
  !
  ! Revision 1.3  1997/10/31  23:54:28  pincus
  ! Changed readAttributeValue to generic subroutine; added
  ! optional versionNumber argument to getSingleHdfFileName; made
  ! public definition of ToolkitMaxStringLength.
  !
  ! Revision 1.2  1997/07/29  22:37:38  pincus
  ! It turns out that the version argument to PGS_PC_GetReference is
  ! used for both input AND output. Changed the calling routines to reflect
  ! this; most importantly I stopped passing this variable a literal.
  !
  ! Revision 1.1  1997/07/14  23:24:37  pincus
  ! Initial revision
  !
  use ErrorHandling,  only: WarningMessage, FatalErrorMessage
  use CharacterUtils, only: CharToInt
  implicit none
  private
  INCLUDE 'PGS_SMF.f'
  INCLUDE 'PGS_PC.f'
  INCLUDE 'PGS_PC_9.f'
  
  integer, parameter, public :: MaxFileNameLength         = PGSd_PC_FILE_PATH_MAX,   &
                                MaxUniversalRefNameLength = PGSd_PC_UREF_LENGTH_MAX, &
                                ToolkitMaxStringLength    = PGSd_PC_VALUE_LENGTH_MAX

  ! Publically available functions. 
  public :: readRuntimeParameter, getNumberOfFiles, getHdfFileNames, getSingleHdfFileName
            
  interface readRuntimeParameter
    module procedure readRuntimeParameterInt, readRuntimeParameterChar
  end interface 
  ! !F90
  !
  ! !Description:
  !    Reads a single parameter associated with a given LUN in a PCF file. 
  !    Generic function - versions exist for integer, character.
  !  
  ! !Input Parameters:
  !    LogicalUnitNumber: The LUN used in the PCF file. 
  !
  ! !Output Parameters:
  !    Returns the value of the parameter associated with the LUN. May be integer or 
  !      character string of length PGSd_PC_VALUE_LENGTH_MAX
  ! 
  ! !Revision History:
  !   See module revision history in data section. 
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
  !   Named Constants:
  !     PGS_S_SUCCESS            (PGS_SMF.f)
  !     PGSd_PC_VALUE_LENGTH_MAX (PGS_PC.f)
  !
  ! !END

  ! Interface definitions copied from the SDP toolkit user's guide, version 
  !   5.2. Putting the interface in the module makes it harder to make 
  !   mistakes when calling toolkit functions. 
  interface 
    function PGS_PC_GetNumberOfFiles(LogicalUnitNumber, numberOfFiles)
      integer, intent( in) :: LogicalUnitNumber
      integer, intent(out) :: numberOfFiles
      integer              :: PGS_PC_GetNumberOfFiles
    end function PGS_PC_GetNumberOfFiles

    function PGS_PC_GetReference(HDF_INFILE, version, PhysicalFileName)
      integer, intent( in)               :: HDF_INFILE, version
      character (len = *), intent(inout) :: PhysicalFileName
      integer                            :: PGS_PC_GetReference
    end function PGS_PC_GetReference

    function PGS_PC_GetConfigData(LogicalUnitNumber, runtimeString)
      integer,             intent( in) :: LogicalUnitNumber
      character (len = *), intent(out) :: runtimeString
      integer                          :: PGS_PC_GetConfigData
    end function PGS_PC_GetConfigData
  end interface 

contains
! ---------------
  subroutine readRuntimeParameterInt(LogicalUnitNumber, intValue)

  use CharacterUtils, only: IntToChar

    integer, intent( in) :: LogicalUnitNumber
    integer, intent(out) :: intValue

    ! !F90
    !
    ! !Description:
    !    Reads a single parameter associated with a given LUN in a PCF file. 
    !  
    ! !Input Parameters:
    !    LogicalUnitNumber: The LUN used in the PCF file. 
    !
    ! !Output Parameters:
    !    Returns the integer value of the parameter associated with the LUN.
    ! 
    ! !Revision History:
    !   See module revision history in data section. 
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
    !   Named Constants:
    !     PGS_S_SUCCESS            (PGS_SMF.f)
    !     PGSd_PC_VALUE_LENGTH_MAX (PGS_PC.f)
    !
    ! !END

    ! Local variables
    character (len = PGSd_PC_VALUE_LENGTH_MAX) :: runtimeString
    integer                                    :: status
    character (len = 2048)                     :: msgbuffer
  

    status = PGS_PC_GetConfigData(LogicalUnitNumber, runtimeString)

    if(status /= PGS_S_Success) then
         msgbuffer = "Error reading runtime parameter from LUN" // &
                      IntToChar(LogicalUnitNumber) &
  // char(10) // "Operator Action:  Check for valid PCF file.  If wrong or " &
  // char(10) // "corrupted, stage correct PCF and rerun PGE.  Otherwise, " &
  // char(10) // "notify SDST. "
         call FatalErrorMessage("readRuntimeParameterInt", msgbuffer)
    end if

    intValue = CharToInt(trim(runtimeString))
  end subroutine readRuntimeParameterInt
! ---------------
  subroutine readRuntimeParameterChar(LogicalUnitNumber, charValue)

    use CharacterUtils, only: IntToChar

    integer,             intent( in) :: LogicalUnitNumber
    character (len = *), intent(out) :: charValue
    ! !F90
    !
    ! !Description:
    !    Reads a single parameter associated with a given LUN in a PCF file. 
    !  
    ! !Input Parameters:
    !    LogicalUnitNumber: The LUN used in the PCF file. 
    !
    ! !Output Parameters:
    !    Returns the value of the parameter associated with the LUN as a 
    !      character string of length PGSd_PC_VALUE_LENGTH_MAX
    ! 
    ! !Revision History:
    !   See module revision history in data section. 
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
    !   Named Constants:
    !     PGS_S_SUCCESS            (PGS_SMF.f)
    !     PGSd_PC_VALUE_LENGTH_MAX (PGS_PC.f)
    !
    ! !END

    ! Local variables
    character (len = PGSd_PC_VALUE_LENGTH_MAX) :: runtimeString
    integer                                    :: status
    character (len = 2048)                     :: msgbuffer

    status = PGS_PC_GetConfigData(LogicalUnitNumber, runtimeString)

    if(status /= PGS_S_Success) then
         msgbuffer = "Error reading runtime parameter from LUN" // &
                      IntToChar(LogicalUnitNumber) &
  // char(10) // "Operator Action:  Check for valid PCF file.  If wrong or " &
  // char(10) // "corrupted, stage correct PCF and rerun PGE.  Otherwise, " &
  // char(10) // "notify SDST. "
         call FatalErrorMessage("readRuntimeParameterChar", msgbuffer)
    end if

    charValue = ""

    if(len_trim(runtimeString) > len(charValue)) then
      msgbuffer = "Not enough space supplied for character runtime parameter." &
// char(10) // "Operator Action:  Suspect system problem.  If a fault is " &
// char(10) // "identified, correct and rerun PGE.  Otherwise, notify SDST." 
      call FatalErrorMessage("readRuntimeParameterChar", msgbuffer)
    end if
    
    charValue = runtimeString
  end subroutine readRuntimeParameterChar
! ---------------
  function getNumberOfFiles(LogicalUnitNumber)
    integer, intent( in) :: LogicalUnitNumber
    integer              :: getNumberOfFiles
    ! !F90
    !
    ! !Description:
    !  
    ! !Input Parameters:
    ! LogicalUnitNumber: Number associated with input product from PCF file
    !
    ! !Output Parameters:
    ! getNumberOfFiles returns an integer descirbing how many physical files are assocated
    !   with the given input product.
    ! 
    ! !Revision History:
    !   See module revision history in data section. 
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
    !   Named Constants:
    !     PGS_S_SUCCESS   (PGS_SMF.f)
    !
    ! !END
  
   ! Local variables 
    integer                            :: status
    character (len = 2048)             :: msgbuffer

    status = PGS_PC_GetNumberOfFiles(LogicalUnitNumber, getNumberOfFiles)
    select case(status)
      case(PGS_S_Success)
        status = PGS_S_Success ! No operation
      case(PgsPC_W_No_Files_For_ID)
        call WarningMessage("getNumberOfFiles", "No files available on LUN ", LogicalUnitNumber)
        getNumberOfFiles = 0
      case default
        msgbuffer = "Error return from PGS_PC_GetNumberOfFiles." &
        // char(10) // "Operator Action:  Notify SDST."
        call FatalErrorMessage("getNumberOfFiles", msgbuffer)
    end select    
  end function getNumberOfFiles
  ! ---------------
  subroutine getHdfFileNames(LogicalUnitNumber, fileNames)

  use CharacterUtils, only: IntToChar

    integer,                          intent( in) :: LogicalUnitNumber
    character(len = *), dimension(:), intent(out) :: fileNames

    ! !F90
    !
    ! !Description:
    !  
    ! !Input Parameters:
    !   LogicalUnitNumber: Number associated with input product from PCF file
    !
    ! !Output Parameters:
    !   fileNames: a one-dimensional array of character stings. Each element is filled with the
    !     physical file name associated with the corresponding version of the input file. 
    ! 
    ! !Revision History:
    !   See module revision history in data section. 
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
    !   Named Constants:
    !     PGS_S_SUCCESS   (PGS_SMF.f)
    !
    ! !END

  ! Local variables 
    integer                       :: status, fileNum, version, numFiles
    character (len = 2048)        :: msgbuffer
  
    numFiles = getNumberOfFiles(LogicalUnitNumber)

    if(numFiles > size(fileNames)) then
      msgbuffer = "Haven't allowed enough space for filenames." &
                   // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("getHdfFileNames", msgbuffer)
    end if

    if(len(fileNames(1)) < MaxFileNameLength) then
      msgbuffer = "String length for file names may be too short." &
                  // char(10) // "Operator Action:  Notify SDST."
      call WarningMessage("getHdfFileNames", msgbuffer)
    end if

    do fileNum = 1, numFiles
      version = fileNum
      status = PGS_PC_GetReference(LogicalUnitNumber, version, fileNames(fileNum))

      if (status /= PGS_S_Success) then
        msgbuffer = "Error return from PGS_PC_GetReference on file number" // &
                     IntToChar(fileNum) &
                  // char(10) // "Operator Action:  Notify SDST."
        call FatalErrorMessage("getHdfFileNames", msgbuffer)
      end if

    end do
    fileNames(numFiles + 1:) = ""
  end subroutine getHdfFileNames
! ---------------
  subroutine getSingleHdfFileName(LogicalUnitNumber, VersionNumber, fileName)
    integer,             intent( in) :: LogicalUnitNumber
    integer,   optional, intent( in) :: versionNumber
    character(len = *),  intent(out) :: fileName
    ! !F90
    !
    ! !Description:
    !  
    ! !Input Parameters:
    !   LogicalUnitNumber: Number associated with input product from PCF file
    !   VersionNumber: Optionally, the version number of the file to read. 
    !
    ! !Output Parameters:
    !   fileName: a single character string with the physical file name of the first instance
    !     of the file. 
    ! 
    ! !Revision History:
    !   See module revision history in data section. 
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
    !   Named Constants:
    !     PGS_S_SUCCESS   (PGS_SMF.f)
    !
    ! !END

  ! Local variables 
    integer                          :: status, version
    character (len = 2048)           :: msgbuffer
    version = 1
    if(present(VersionNumber)) version = VersionNumber

    status = PGS_PC_GetReference(LogicalUnitNumber, version, fileName)

    if (status /= PGS_S_Success) then
      msgbuffer = "Error return from PGS_PC_GetReference." &
          // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("getSingleHdfFileName",  msgbuffer)
    end if

  end subroutine getSingleHdfFileName
! ---------------
end module FileHandling
