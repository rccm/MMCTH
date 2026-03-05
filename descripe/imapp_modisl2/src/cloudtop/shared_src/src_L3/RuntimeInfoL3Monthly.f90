module RuntimeInfoL3Monthly
  use typeSizes,     only: EightByteRealKind, FourByteIntKind
  use hdf,           only: SFsnatt, DFNT_INT32, SUCCEED, MAX_NC_NAME, sfscatt, DFNT_CHAR8
  use ErrorHandling, only: FatalErrorMessage, ErrorMessage
! MODIFICATION 2/3/99 P.HUBANKS - Added MaxFileNameLength
  use FileHandling,  only: getNumberOfFiles, getHdfFileNames,          & 
                           getSingleHdfFileName, readRuntimeParameter, &
                           ToolkitMaxStringLength, MaxUniversalRefNameLength, &
                           MaxFileNameLength
  use TimeHandling,  only: getPGETime, getExecutionTime,                              &
                           TaiTimeToString, StringToTaiTime, TimeFormatAtoB,          &
                           FormatAStringLength, FormatBStringLength
  use Metadata,      only: getArchiveMetadata, getCoreMetadata,                       &
                           getMetadataAttribute, getUniversalReference,               &
                           initializeMetadata, writeMetadata,                         &
                           setMetadataAttribute,                                      &
                           ! Constants
                           ToolkitEndDataString, ECSParameterNameLength,              &
                           metadataGroupNameLength, metadataNumberOfGroups,           &
                           ! Metadata names
                           MCORE_EAST_BOUND,  MCORE_NORTH_BOUND,                      &
                           MCORE_SOUTH_BOUND, MCORE_WEST_BOUND,                       &
                           MCORE_RANGE_BEG_DATE,    MCORE_RANGE_BEG_TIME,             &
                           MCORE_RANGE_ENDING_DATE, MCORE_RANGE_ENDING_TIME,          &
                           MCORE_DAYNIGHTFLAG, MCORE_INPUT_POINTER, MCORE_SHORT_NAME, &
                           MCORE_LOCALINPUTGRANULEID,                                 &
                           ArchAttrStringLengthMax, ArchiveMetadataHDFAttrName,       &
                           InventoryMetadataHDFAttrName, MAX_ESDT_Name_LEN

  implicit none
  private

  ! The Logical Unit Numbers that we'll look for in the PCF are defined here. 
  ! fhliang 09/18/98: changed MonthlyMetadataLUN from 461000 to 460001
  ! rhucek  06/06/02: renamed RuntimeParametersLUN to RuntimeParametersLUN_I and added
  !                   parameter RuntimeParametersLUN_A (archive RP start LUN)
  ! rhucek  05/24/03: Added new parameter, TemplateFileLUN

  integer, parameter :: DailyFileLUN            = 450000, &
                        MonthlyFileLUN          = 460000, &
                        MonthlyMetadataLUN      = 460001, &
                        TemplateFileLUN         = 460002, &
                        RuntimeParametersLUN_I  = 462000, &
                        RuntimeParametersLUN_A  = 462500, &
                        LocalVersionLUN         = 463000, &
                        StartTimeLUN            = 10258,  &
                        StopTimeLUN             = 10259
                                                                               
   ! Metadata parameters are defined in the SDP User's Guide on page 6-52.
   integer, parameter :: InventoryMetadata = 2, ArchiveMetadata = 3, &
                         LocalVersionStringMaxLen = 3, & ! per code from Rich Hucek
                         LocalGranuleIDLength     = MAX_ESDT_Name_LEN + 45, &
                         ! From the tile process MCF; better make darn any changes there get reflected here.
                         MaxL2FilesPerTile  = 250
   
   ! Which procedures will be publically available? 

  public :: getRuntimeInfo, setOutputMetadata
  
contains
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
    !   $Log: RuntimeInfoL3Monthly.f90,v $
    ! 
    !# Revision 1.5  1998/02/19  15:10:49  xliang
    !# This is the first delivery version.
    !#
    !# Revision 1.4  1998/01/23  16:30:38  xliang
    !# This version is the SAME as previous version except that it has
    !# some of the unnecessary parts (used for testing) taken out.
    !#
    !# Revision 1.3  1998/01/23  15:16:40  xliang
    !# The write Metadata part is added, but couldn't be tested due to the
    !# problems of PGS toolkit.
    !#
    !# Revision 1.2  1997/12/31  14:39:38  xliang
    !# This version has the limits enlarged and so it can run the entire SDS's list.
    !#
    !# Revision 1.1  1997/12/30  21:27:49  xliang
    !# Initial revision
    !#
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
    
    character(len = ToolkitMaxStringLength) :: BeginDateString, BeginTimeString
    
    call getCoreMetadata(fileID = LogicalUnitNumber, version = version, &
                         ECSparamName = "RangeBeginningDate", paramValue = BeginDateString)
    call getCoreMetadata(fileID = LogicalUnitNumber, version = version, &
                         ECSparamName = "RangeBeginningTime", paramValue = BeginTimeString)
    getFileBeginTime = StringToTaiTime(trim(BeginDateString) // "T" // trim(BeginTimeString) // "Z")
    !getFileBeginTime = StringToTaiTime("1996-08-01"  // "T" // "00:00:00" // "Z") !testing without CoreMetadata
    
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
    !   $Log: RuntimeInfoL3Monthly.f90,v $
    !# Revision 1.5  1998/02/19  15:10:49  xliang
    !# This is the first delivery version.
    !#
    !# Revision 1.4  1998/01/23  16:30:38  xliang
    !# This version is the SAME as previous version except that it has
    !# some of the unnecessary parts (used for testing) taken out.
    !#
    !# Revision 1.3  1998/01/23  15:16:40  xliang
    !# The write Metadata part is added, but couldn't be tested due to the
    !# problems of PGS toolkit.
    !#
    !# Revision 1.2  1997/12/31  14:39:38  xliang
    !# This version has the limits enlarged and so it can run the entire SDS's list.
    !#
    !# Revision 1.1  1997/12/30  21:27:49  xliang
    !# Initial revision
    !#
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
    character(len = ToolkitMaxStringLength) :: EndDateString,   EndTimeString

    call getCoreMetadata(fileID = LogicalUnitNumber, version = version, &
                         ECSparamName = "RangeEndingDate", paramValue = EndDateString)
    call getCoreMetadata(fileID = LogicalUnitNumber, version = version, &
                         ECSparamName = "RangeEndingTime", paramValue = EndTimeString)
    getFileEndTime = StringToTaiTime(trim(EndDateString) // "T" // trim(EndTimeString) // "Z")
    !getFileEndTime = StringToTaiTime("1996-08-03" // "T" // "23:59:59" // "Z") !testing without CoreMetadata
    
  end function getFileEndTime
  ! -------------------------------------------------
  function compareFileTimesToPGE(InputProductLUN, numVersions)
    integer,                         intent( in) :: InputProductLUN, numVersions
    logical, dimension(numVersions)              :: compareFileTimesToPGE
    ! !F90
    !
    ! !Description:
    !    Determines which files on a given LUN fall within the time boundaries specified in the PGE. 
    !  
    ! !Input Parameters:
    !  InputProductLUN: The Logical Unit Number from which to get the file names.
    !  numVersions: the number of files expected on the LUN. This is passed as an input
    !    argument so we can return a vector as the function result. 
    !
    ! !Output Parameters:
    !    None. Function result is a rank one logical array the same length as the number of
    !      versions of the input file. 
    ! 
    ! !Revision History:
    !   $Log: RuntimeInfoL3Monthly.f90,v $
    !# Revision 1.5  1998/02/19  15:10:49  xliang
    !# This is the first delivery version.
    !#
    !# Revision 1.4  1998/01/23  16:30:38  xliang
    !# This version is the SAME as previous version except that it has
    !# some of the unnecessary parts (used for testing) taken out.
    !#
    !# Revision 1.3  1998/01/23  15:16:40  xliang
    !# The write Metadata part is added, but couldn't be tested due to the
    !# problems of PGS toolkit.
    !#
    !# Revision 1.2  1997/12/31  14:39:38  xliang
    !# This version has the limits enlarged and so it can run the entire SDS's list.
    !#
    !# Revision 1.1  1997/12/30  21:27:49  xliang
    !# Initial revision
    !#
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
    !    getRuntimeInfo calls function in module FileHandling. The LUNs for various 
    !      parameters are set in the data section of this module (RuntimeInfo)
    !
    ! !END
    
    ! Local variables
    character (len = 2048)                               :: msgbuffer
    real (kind = EightByteRealKind)                      :: PGEStartTime, PGEStopTime
    real (kind = EightByteRealKind), &
                               dimension(numVersions)    :: FileBeginTimes, FileEndTimes
    integer                                              :: file

    !
    ! MODIFICATION 2/3/99 P.HUBANKS replaced next line with the line after
    ! character(len = MAX_NC_NAME), dimension(numVersions) :: allFileNames
    character(len = MaxFileNameLength), dimension(numVersions) :: allFileNames
    !

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
    call getHdfFileNames     (InputProductLUN, allFileNames)
    
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
! --------------------------------------------------------------------------------
  subroutine getRuntimeInfo(numL3Dailies, L3DailyFileNames, L3MonthlyFileName, reportWarnings)
    integer,                           optional, intent(out) :: numL3Dailies
    character (len = *), dimension(:), optional, intent(out) :: L3DailyFileNames
    character (len = *),               optional, intent(out) :: L3MonthlyFileName
    logical,                           optional, intent( in) :: reportWarnings
    ! !F90
    !
    ! !Description:
    !    Returns information needed at runtime, obtained from PCF file. 
    !  
    ! !Input Parameters:
    !   reportErrors: if present and value is true, Level 3 daily files in the PCF which 
    !     fall outside the PCF processing time window are reported. 
    !
    ! !Output Parameters:
    !   All argurments are optional
    !   numL3Dailies: the number of daily files in the PCF which fall within the PGE processing time
    !     window.
    !   L3DailyFileNames:  the names of the daily files in the PCF which fall within the PGE processing time
    !     window.
    !   L3MonthlyFileName: the name of the Level 3 monthly output file. 
    ! 
    ! !Revision History:
    !   $Log: RuntimeInfoL3Monthly.f90,v $
    !# Revision 1.5  1998/02/19  15:10:49  xliang
    !# This is the first delivery version.
    !#
    !# Revision 1.4  1998/01/23  16:30:38  xliang
    !# This version is the SAME as previous version except that it has
    !# some of the unnecessary parts (used for testing) taken out.
    !#
    !# Revision 1.3  1998/01/23  15:16:40  xliang
    !# The write Metadata part is added, but couldn't be tested due to the
    !# problems of PGS toolkit.
    !#
    !# Revision 1.2  1997/12/31  14:39:38  xliang
    !# This version has the limits enlarged and so it can run the entire SDS's list.
    !#
    !# Revision 1.1  1997/12/30  21:27:49  xliang
    !# Initial revision
    !#
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
    !    getRuntimeInfo calls function in module FileHandling. The LUNs for various 
    !      parameters are set in the data section of this module (RuntimeInfo)
    !
    ! !END

    ! Local variables
    character (len = 2048)                  :: msgbuffer
    integer                                 :: numFiles, file
    integer,      dimension(2)              :: status = 0
    logical,      dimension(:), allocatable :: withinTimeBounds

    !
    ! MODIFICATION 2/3/99 P.HUBANKS replaced next line with the line after
    ! character(len = MAX_NC_NAME), &
    !              dimension(:), allocatable :: allFileNames, badFileNames
    character(len = MaxFileNameLength), &
                  dimension(:), allocatable :: allFileNames, badFileNames
                  
    if(present(L3MonthlyFileName))  &
      call getSingleHdfFileName(MonthlyFileLUN, VersionNumber = 1, fileName = L3MonthlyFileName)

    if(present(numL3Dailies) .or. present(L3DailyFileNames) ) then
      ! We return the number and the names of only those Level 2 granules which fall within
      !   the processing time window for this PCF. 
      numFiles = getNumberOfFiles(DailyFileLUN)
      allocate(withinTimeBounds(numFiles), stat = status(1))
      allocate(allFileNames    (numFiles), stat = status(2))

      if(any(status /= 0)) then
        msgbuffer = "Can't allocate memory." &
  // char(10) // "Operator Action:  Suspect system problem.  If a fault is " &
  // char(10) // "identified, correct and rerun PGE.  Otherwise, notify SDST." 
        call FatalErrorMessage("getRuntimeInfo", msgbuffer)
      end if
        
      withinTimeBounds(:) = compareFileTimesToPGE(InputProductLUN = DailyFileLUN, &
                                                  numVersions     = numFiles)
      call getHdfFileNames(DailyFileLUN, allFileNames)

      ! We report the names of those file which don't belong in the PGE if we're asked to.
      reportOutOfBoundFiles: if(present(reportWarnings)) then
        if(reportWarnings .and. .not. all(withinTimeBounds) ) then
          allocate(badFileNames(numFiles - count(withinTimeBounds)), stat = status(1))

          if (status(1) /= 0) then
            msgbuffer = "Can't allocate memory." &
  // char(10) // "Operator Action:  Suspect system problem.  If a fault is " &
  // char(10) // "identified, correct and rerun PGE.  Otherwise, notify SDST." 
            call FatalErrorMessage("getRuntimeInfo", msgbuffer)
          end if

           badFileNames = pack(array = allFileNames, mask = .not. withinTimeBounds)

           do file = 1, size(badFileNames)

             msgbuffer = "File " // trim(badFileNames(file)) // &
                         " does not fall within PCF processing time window." &
                         // char(10) // "Operator Action:  Notify SDST."
             call ErrorMessage("getRuntimeInfo", msgbuffer)

           end do

           deallocate(badFileNames, stat = status(1))

           if (status(1) /= 0) then
               msgbuffer = "Can't deallocate memory." &
  // char(10) // "Operator Action:  Suspect system problem.  If a fault is " &
  // char(10) // "identified, correct and rerun PGE.  Otherwise, notify SDST." 
               call FatalErrorMessage("getRuntimeInfo", msgbuffer)
           end if

        end if
      end if reportOutOfBoundFiles

      ! Both the number (here) and the names (below) of the input files reflect only those
      !   that fall within the processing time window for the PGE. 
      if(present(numL3Dailies))      &
        numL3Dailies = count(withinTimeBounds(:))
        
      returnDailyFileNames: if(present(L3DailyFileNames)) then

        if(size(L3DailyFileNames) < count(withinTimeBounds(:))) then
          msgbuffer = "Not enough space supplied for file names." &
  // char(10) // "Operator Action:  Suspect system problem.  If a fault is " &
  // char(10) // "identified, correct and rerun PGE.  Otherwise, notify SDST." 
          call FatalErrorMessage("getRuntimeInfo", msgbuffer)
        end if

        L3DailyFileNames(:) = pack(array = allFileNames(:), mask = withinTimeBounds(:))
        L3DailyFileNames(count(withinTimeBounds(:)) + 1: ) = ""
      end if returnDailyFileNames
      
      deallocate(withinTimeBounds, stat = status(1))
      deallocate(allFileNames,     stat = status(2))

      if(any(status /= 0)) then
        msgbuffer = "Can't deallocate memory." &
  // char(10) // "Operator Action:  Suspect system problem.  If a fault is " &
  // char(10) // "identified, correct and rerun PGE.  Otherwise, notify SDST." 
        call FatalErrorMessage("getRuntimeInfo", msgbuffer)
      end if

    end if
  end subroutine getRuntimeInfo
   ! -------------------------------------------------
  function GenerateModisFilename(EsdtShortName, LocalVersionString) 
    character(len = *),         intent( in) :: EsdtShortName, LocalVersionString
    character(len = LocalGranuleIDLength)   :: GenerateModisFilename
    ! !F90
    !
    ! !Description:
    !    Generate a MODIS file name based on the Level 2 convention. 
    !  
    ! !Input Parameters:
    !   EsdtShortName: Short name of the output product
    !   LocalVersionString: Version of the file being produced, read as a string
    !
    ! !Output Parameters:
    !   Function result is the MODIS file name. 
    ! 
    ! !Revision History:
    !   $Log: RuntimeInfoL3Monthly.f90,v $
    !# Revision 1.5  1998/02/19  15:10:49  xliang
    !# This is the first delivery version.
    !#
    !# Revision 1.4  1998/01/23  16:30:38  xliang
    !# This version is the SAME as previous version except that it has
    !# some of the unnecessary parts (used for testing) taken out.
    !#
    !# Revision 1.3  1998/01/23  15:16:40  xliang
    !# The write Metadata part is added, but couldn't be tested due to the
    !# problems of PGS toolkit.
    !#
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
    character (len = len("BadFileName")), &
                                parameter      :: ErrorString = "BadFileName"
    integer                                    :: counter
    character (len = FormatBStringLength)      :: PGEStartString, RunTimeStartString
    character (len = LocalVersionStringMaxLen) :: versionString
    character (len = 2048)                     :: msgbuffer

    GenerateModisFilename = ""                           
    ! Check validity and length of input arguments
    if(len_trim(EsdtShortName) > MAX_ESDT_Name_LEN) then
       GenerateModisFilename = ErrorString
           msgbuffer = "EsdtShortName " // trim(EsdtShortName) // " exceeds allowed length of " // &
                       "MAX_ESDT_Name_LEN" &
  // char(10) // "Operator Action:  Check for valid PCF file.  If wrong or " &
  // char(10) // "corrupted, stage correct PCF and rerun PGE.  Otherwise, " &
  // char(10) // "notify SDST. "
           call ErrorMessage("GenerateModisFilename", msgbuffer)
    end if
    if(EsdtShortName == "") then
       GenerateModisFilename = ErrorString
       msgbuffer = "EsdtShortName is blank." &
  // char(10) // "Operator Action:  Check for valid PCF file.  If wrong or " &
  // char(10) // "corrupted, stage correct PCF and rerun PGE.  Otherwise, " &
  // char(10) // "notify SDST. "
       call ErrorMessage("GenerateModisFilename", msgbuffer)
    end if
    if(len_trim(LocalVersionString) > LocalVersionStringMaxLen) then
       GenerateModisFilename = ErrorString
          msgbuffer = "LocalVersionString " // trim(LocalVersionString) // &
                      " exceeds allowed length of " // & 
                      "LocalVersionStringMaxLen" &
  // char(10) // "Operator Action:  Check for valid PCF file.  If wrong or " &
  // char(10) // "corrupted, stage correct PCF and rerun PGE.  Otherwise, " &
  // char(10) // "notify SDST. "
          call ErrorMessage("GenerateModisFilename", msgbuffer)
    end if
    if(LocalVersionString == "") then
        GenerateModisFilename = ErrorString
        msgbuffer = "LocalVersionString is blank." &
  // char(10) // "Operator Action:  Check for valid PCF file.  If wrong or " &
  // char(10) // "corrupted, stage correct PCF and rerun PGE.  Otherwise, " &
  // char(10) // "notify SDST. "
        call ErrorMessage("GenerateModisFilename", msgbuffer)
     end if
     do counter = 1, len_trim(LocalVersionString)
       if(llt(LocalVersionString(counter:counter), '0') .or. &
          lgt(LocalVersionString(counter:counter), '9') ) then
         GenerateModisFilename = ErrorString
           msgbuffer = "LocalVersionString " // trim(LocalVersionString) // &
                       " contains non-numeric characters." &
  // char(10) // "Operator Action:  Check for valid PCF file.  If wrong or " &
  // char(10) // "corrupted, stage correct PCF and rerun PGE.  Otherwise, " &
  // char(10) // "notify SDST. "
           call ErrorMessage("GenerateModisFilename", msgbuffer)

       end if
     end do
     
     if(GenerateModisFilename /= ErrorString) then
       select case(len_trim(LocalVersionString))
         case(1)
           versionString = "00" // trim(LocalVersionString)
         case(2)
           versionString =  "0" // trim(LocalVersionString)
         case(3)
           versionString =         trim(LocalVersionString)
         ! Default not necessary - we checked above that the length was 1, 2, or 3.
       end select
       ! Get PGE Collection Start Time and execution time (from the Fortran 90 Date_And_Time 
       !   intrinsic) in CCSDS Format B.
       PGEStartString     = TimeFormatAtoB(TaiTimeToString(getPGETime(StartTimeLUN)))
       RunTimeStartString = TimeFormatAtoB(TaiTimeToString(getExecutionTime(      )))
     end if
     GenerateModisFilename = trim(EsdtShortName) // ".A" // &
                              PGEStartString    (1 :4)  // PGEStartString    (6 :8 ) // "." // &
                              versionString             //                              "." // &
! fhliang 07/31/98: removed '.':
!                             RunTimeStartString(1 :4)  // RunTimeStartString(6 :8 ) // "." // &
                              RunTimeStartString(1 :4)  // RunTimeStartString(6 :8 ) //        &
                              RunTimeStartString(10:11) // RunTimeStartString(13:14) //        &
                              RunTimeStartString(16:17) // ".hdf"
     if(len_trim(GenerateModisFilename) < len(GenerateModisFilename)) &
       GenerateModisFilename(len_trim(GenerateModisFilename) + 1:len(GenerateModisFilename)) = char(0)
    
  end function GenerateModisFilename
  ! -------------------------------------------------
  subroutine setOutputMetadata(outputFileID)
    integer (kind = FourByteIntKind) :: outputFileID
    ! !F90
    !
    ! !Description:
    !    Writes required metadata to the output file. 
    !  
    ! !Input Parameters:
    !   outputFileID: the HDF file ID as used by the SD interface.
    !
    ! !Output Parameters:
    !   None. Core metadata is set in output file.
    ! 
    ! !Revision History:
    !   $Log: RuntimeInfoL3Monthly.f90,v $
    !
    !  Revision 1.7  2003/05/24 rhucek
    !  The product template file is entered into the PCF as a separate input
    !  file LUN (460002).  A writable copy of the template is made and renamed 
    !  by the system to serve as product output file.  The template file is 
    !  identified in the ECS metadata (under the InputPointer field) providing 
    !  an unambiguous record of the product version used. 
    !  
    !  Revision 1.7  2002/06/06 rhucek
    !  Added code to set Archive metadata attributes stored as runtime
    !  parameters.  The new code is very nearly a copy of original Pincus 
    !  statements used for setting Inventory metadata runtime parameters.
    !
    !  Revision 1.6  2000/03/05 rhucek   
    !  New code was added to identify the product file template as an additional 
    !  element in the ECS INPUTPOINTER metadata array. 
    !
    !# Revision 1.5  1998/02/19  15:10:49  xliang
    !# This is the first delivery version.
    !#
    !# Revision 1.4  1998/01/23  16:30:38  xliang
    !# This version is the SAME as previous version except that it has
    !# some of the unnecessary parts (used for testing) taken out.
    !#
    !# Revision 1.3  1998/01/23  15:16:40  xliang
    !# The write Metadata part is added, but couldn't be tested due to the
    !# problems of PGS toolkit.
    !#
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
    !    The individual ECS inventory metadata fields written by 
    !    setOutputMetadata, the source of the field values, and the  
    !    "Data Location" parameter of Metadata Configuration File (MCF) are 
    !    listed below.  Data sources are identified as follows:
    ! 
    !    MCF      - Metadata Configuration File
    !    TK       - set by SDP toolkit
    !    PGE      - Passed in from calling routine
    !    PCF      - Process Control File (RP denotes Runtime Parameters)
    !    internal - computed internally within function setOutputMetadata
    !    static   - predefined value set within function setOutputMetadata
    !
    !                         ECS Inventory Metadata Fields
    !                     -----------------------------
    !                                           Source of    Data Location
    !           Metadata Object Name              Value      Value in MCF 
    !    -------------------------------        ---------    -------------
    !       1    ShortName                         MCF           MCF
    !       2    VersionID                         MCF           MCF
    !       3    ProductionDateTime                TK            TK
    !       4    InputPointer                      PGE           PGE
    !       5    LocalGranuleID                  internal        PGE
    !       6    DayNightFlag                     static         PGE
    !       7    LocalVersionID                    PCF (RP)      PCF
    !       8    PGEVersion                        PCF (RP)      PCF
    !       9    RangeBeginningTime                PCF           PGE
    !      10    RangeEndingTime                   PCF           PGE
    !      11    RangeBeginningDate                PCF           PGE
    !      12    RangeEndingDate                   PCF           PGE
    !      13    EastBoundingCoordinate           static         PGE
    !      14    WestBoundingCoordinate           static         PGE 
    !      15    NorthBoundingCoordinate          static         PGE
    !      16    SouthBoundingCoordinate          static         PGE 
    !      36    AssociatedPlatformShortName.1     PCF (RP)      PCF
    !      37    AssociatedInstrumentShortName.1   PCF (RP)      PCF
    !      38    AssociatedSensorShortName.1       PCF (RP)      PCF
    ! 
    ! 
    !                     ECS Archive Metadata Fields
    !                     -----------------------------
    !                                           Source of    Data Location
    !        Metadata Object Name                 Value      Value in MCF
    !    -------------------------------        ---------    -------------
    !       1    LocalInputGranuleID               PGE           PGE 
    !       2    ProcessingEnvironment             PGE           PGE 
    !
    ! !END
      
    ! Local variables
    character (len = 2048)                    :: msgbuffer
    integer                                   :: numTileFiles
    integer, dimension(:), allocatable        :: versions
    logical, dimension(:), allocatable        :: withinTimeBounds
    character(len = LocalGranuleIDLength), &
             dimension(:), allocatable        :: inputGranuleIDs
    character(len = MaxUniversalRefNameLength), &
             dimension(:), allocatable        :: universalRefs
    integer, dimension(2)                     :: status = 0
    integer                                   :: numGoodFiles,                        &
                                                 numRuntimeMetadataElements, counter, &
                                                 currentPosition, endOfDataPosition
    !character(len = ToolkitMaxStringLength), &
    !         dimension(:), allocatable       :: universalRefs, inputGranuleIDs
    !integer, dimension(2)                    :: status = 0
    !integer                                  :: numGoodFiles, numRuntimeMetadataElements, counter
    
    character (len = ToolkitMaxStringLength)  :: CollectionStartString, CollectionEndString, &
                                                 EsdtShortName, localVersionString,          &
                                                 RuntimeMetadataFieldName, RuntimeMetadataValue

    ! rhucek 12/04/03: declared scalars InputProdHistory and ProductionHistory
    character (len = ArchAttrStringLengthMax)  :: InputProdHistory  = "", &
                                                  ProductionHistory = "" 

    character (len = metadataGroupNameLength), &
            dimension(metadataNumberOfGroups) :: metadataHandles
  
    integer pgs_pc_getconfigdata
    character*4 :: pcf_satid
    character(len = 40), dimension(40) :: doi_char
    character(len = 30), dimension(30) :: doiAuth
    integer :: nRet
    character*40 AttrN

    CHARACTER*(*) MCORE_ADD_ATTRIBUTENAME
    PARAMETER ( MCORE_ADD_ATTRIBUTENAME = 'ADDITIONALATTRIBUTENAME')

    CHARACTER*(*) MCORE_PARAMETERVALUE
    PARAMETER (MCORE_PARAMETERVALUE = 'PARAMETERVALUE')

    ! Find out how many files there are.
    numTileFiles = getNumberOfFiles(DailyFileLUN)
    
    allocate(versions        (numTileFiles), stat = status(1))
    allocate(withinTimeBounds(numTileFiles), stat = status(2))

    if(any(status /= 0)) then
        msgbuffer = "Can't allocate memory." &
  // char(10) // "Operator Action:  Suspect system problem.  If a fault is " &
  // char(10) // "identified, correct and rerun PGE.  Otherwise, notify SDST." 
        call FatalErrorMessage("setOutputMetadata", msgbuffer)
    end if

    ! Determine which versions of which product fall within the processing time window 
    !   of the PGE. Discard those that don't.                                            
    versions        (:) = (/ ( counter, counter = 1, numTileFiles) /) 
    withinTimeBounds(:) = compareFileTimesToPGE(InputProductLUN = DailyFileLUN, &
                                                numVersions     = numTileFiles)
        
    numGoodFiles            = count( withinTimeBounds(:) )
    versions(:numGoodFiles) = pack(array = versions, mask = withinTimeBounds(:))
    versions(numGoodFiles + 1:) = 0
    
    ! Create arrays for the Universal References and Local Granule IDs and populate them. 
    !   We include an end-of-data marker from the toolkit. 
    !allocate(universalRefs  (numGoodFiles * MaxL2FilesPerTile + 1), stat = status(1) ) 
    !allocate(inputGranuleIDs(numGoodFiles * MaxL2FilesPerTile + 1), stat = status(2) )

    ! rhucek 03/21/00:  increased allocation of URs to numGoodFiles + 2 
    allocate(universalRefs  (numGoodFiles + 2), stat = status(1) ) 
    allocate(inputGranuleIDs(numGoodFiles + 1), stat = status(2) )

    if(any(status /= 0)) then
        msgbuffer = "Can't allocate memory." &
  // char(10) // "Operator Action:  Suspect system problem.  If a fault is " &
  // char(10) // "identified, correct and rerun PGE.  Otherwise, notify SDST." 
        call FatalErrorMessage("setOutputMetadata", msgbuffer)
    end if

    ! rhucek 03/21/00:  set first UR to identify product template file
    ! First get Universal Reference for product template file  
    ! call getUniversalReference( MonthlyFileLUN, 1, universalRefs(1) )

    ! rhucek 05/24/03: reference the product template file on LUN TemplateFileLUN, 
    !                  not the product output file on MonthlyFileLUN.  Also added
    !                  template file to inputgranuleIDs
    call getUniversalReference( TemplateFileLUN, 1, universalRefs(1) )

    do counter = 1, numGoodFiles
      ! rhucek 03/21/00:  reset array index for URs to counter+1
      call getUniversalReference(DailyFileLUN, versions(counter), universalRefs(counter+1) )
      call getCoreMetadata      (DailyFileLUN, versions(counter), "LOCALGRANULEID", & 
                                                                  inputGranuleIDs(counter) )
      ! rhucek 12/04/03: Retrieve "Archive attribute PRODUCTIONHISTORY" from first 
      ! available MOD08_D3 daily input file 
      if (InputProdHistory == "")  call getArchiveMetadata                        & 
         (DailyFileLUN, versions(counter), "PRODUCTIONHISTORY", InputProdHistory)

    end do

    ! rhucek 03/21/00:  reset UR array index to numGoodFiles+2 
    universalRefs  (numGoodFiles + 2) = ToolkitEndDataString
    inputGranuleIDs(numGoodFiles + 1) = ToolkitEndDataString
    
    ! Free up unneeded arrays...
    deallocate(versions,         stat = status(1))
    deallocate(withinTimeBounds, stat = status(2))

    if(any(status /= 0)) then
         msgbuffer = "Can't allocate memory." &
  // char(10) // "Operator Action:  Suspect system problem.  If a fault is " &
  // char(10) // "identified, correct and rerun PGE.  Otherwise, notify SDST." 
         call ErrorMessage("setOutputMetadata", msgbuffer)
    end if
                                
    ! ---- Begin setting metadata ---
    ! Initialize metadata tool defining array Handles, and set 
    ! inventory attribute fields whose values are provided in the MCF.
    call initializeMetadata(MonthlyMetadataLUN, metadataHandles)
    
    ! -- Inventory metadata --
    ! Set ECS Bounding Coordinates. These fields include:
    ! "EastBoundingCoordinate", "NorthBoundingCoordinate" 
    ! "SouthBoundingCoordinate", "WestBoundingCoordinate"
    call setMetadataAttribute(metadataHandles(InventoryMetadata), MCORE_EAST_BOUND,  dble( 180))
    call setMetadataAttribute(metadataHandles(InventoryMetadata), MCORE_WEST_BOUND,  dble(-180))
    call setMetadataAttribute(metadataHandles(InventoryMetadata), MCORE_NORTH_BOUND, dble(  90))
    call setMetadataAttribute(metadataHandles(InventoryMetadata), MCORE_SOUTH_BOUND, dble(- 90))
    
    ! Read the collection start "datetime" from PCF and set 
    !   "RangeBeginningDate" and "RangeBeginningTime"
    call readRuntimeParameter(StartTimeLUN, CollectionStartString)

    if(len_trim(CollectionStartString) == 0) then
        msgbuffer = "Can't read range begin time from LUN " // &
                    "StartTimeLUN" &
  // char(10) // "Operator Action:  Check for valid PCF file.  If wrong or " &
  // char(10) // "corrupted, stage correct PCF and rerun PGE.  Otherwise, " &
  // char(10) // "notify SDST. "
        call FatalErrorMessage("setOutputMetata", msgbuffer)
    end if

    CollectionStartString = AdjustL(CollectionStartString)
    call setMetadataAttribute(metadataHandles(InventoryMetadata), MCORE_RANGE_BEG_DATE,  &
                              CollectionStartString(1:10))
    call setMetadataAttribute(metadataHandles(InventoryMetadata), MCORE_RANGE_BEG_TIME,  &
                              CollectionStartString(12:len_trim(CollectionStartString)))
                              
    ! Read the collection stop "datetime" from PCF and set 
    !   "RangeEndingDate" and "RangeEndingTime"
    call readRuntimeParameter(StopTimeLUN, CollectionEndString)

    if(len_trim(CollectionEndString) == 0) then
        msgbuffer = "Can't read range end time from LUN: " // &
                    "StopTimeLUN" &
  // char(10) // "Operator Action:  Check for valid PCF file.  If wrong or " &
  // char(10) // "corrupted, stage correct PCF and rerun PGE.  Otherwise, " &
  // char(10) // "notify SDST. "
        call FatalErrorMessage("setOutputMetata", msgbuffer)
    end if

    CollectionEndString = AdjustL(CollectionEndString)
    call setMetadataAttribute(metadataHandles(InventoryMetadata), MCORE_RANGE_ENDING_DATE,  &
                              CollectionEndString(1:10))
    call setMetadataAttribute(metadataHandles(InventoryMetadata), MCORE_RANGE_ENDING_TIME,  &
                              CollectionEndString(12:len_trim(CollectionEndString)))
                              
    ! Set Day/Night flag
    call setMetadataAttribute(metadataHandles(InventoryMetadata), MCORE_DAYNIGHTFLAG, "Both")
    
    ! Construct and write out LocalGranuleID. The version ID is also written to 
    !  the inventory metadata.
    call readRuntimeParameter(LocalVersionLUN, localVersionString)
    call setMetadataAttribute(metadataHandles(InventoryMetadata), "LOCALVERSIONID",  &
                              trim(AdjustL(localVersionString)))
    call getMetadataAttribute(metadataHandles(InventoryMetadata), MCORE_SHORT_NAME, EsdtShortName)

    call setMetadataAttribute(metadataHandles(InventoryMetadata), "LOCALGRANULEID",  &
                              GenerateModisFilename(EsdtShortName, LocalVersionString))
    
    ! Write out optional metadata, stored as runtime parameters on alternating LUNs
    !  We expect a fieldname of "None" when we are finished. (Otherwise the read will fail.)
    counter = RuntimeParametersLUN_I
    runtimeMetadata: do
      call readRuntimeParameter(counter,     RuntimeMetadataFieldName)
      if(trim(RuntimeMetadataFieldName) == "None" .or. trim(RuntimeMetadataFieldName) == "NONE") &
        exit runtimeMetadata
      call readRuntimeParameter(counter + 1, RuntimeMetadataValue)
      call setMetadataAttribute(metadataHandles(InventoryMetadata), &
                                RuntimeMetadataFieldName,  &
                                RuntimeMetadataValue)
      counter = counter + 2
    end do runtimeMetadata
    
    ! Set array InputPointers which we determined above. 
    call setMetadataAttribute(metadataHandles(InventoryMetadata), MCORE_INPUT_POINTER,  &
                              universalRefs)

    ! -- Archive metadata --
    ! rhucek 06/06/02:  Added code to set archive metadata stored as runtime parameters.
    !                   This is very nearly a copy of the original source implemented by
    !                   Pincus for inventory metadata.
    !
    ! Write out "Archive" metadata stored as runtime parameters on alternating LUNs
    ! We expect a fieldname of "None" when we are finished. (Otherwise the read will fail.)
    counter = RuntimeParametersLUN_A
    runtimeMetadata_A: do
      call readRuntimeParameter(counter,     RuntimeMetadataFieldName)
      if(trim(RuntimeMetadataFieldName) == "None" .or. trim(RuntimeMetadataFieldName) == "NONE") &
        exit runtimeMetadata_A
      call readRuntimeParameter(counter + 1, RuntimeMetadataValue)

      ! rhucek 12/04/03: Read PROCESSINGPGE and append to prior PRODUCTIONHISTORY 
      if (trim(RuntimeMetadataFieldName) == "PROCESSINGPGE" .or.   &
          trim(RuntimeMetadataFieldName) == "ProcessingPGE" .or.   &
          trim(RuntimeMetadataFieldName) == "processingPGE" .or.   &
          trim(RuntimeMetadataFieldName) == "processingpge" ) then 
          ProductionHistory = trim(InputProdHistory) // ';' // trim(RuntimeMetadataValue)  

          call setMetadataAttribute(metadataHandles(ArchiveMetadata),  &
                                    "PRODUCTIONHISTORY",               & 
                                    ProductionHistory)
      else
          call setMetadataAttribute(metadataHandles(ArchiveMetadata),  &
                                    RuntimeMetadataFieldName,          &
                                    RuntimeMetadataValue)
      endif

      counter = counter + 2
    end do runtimeMetadata_A

    ! Set array InputPointers which we determined above. 
    call setMetadataAttribute(metadataHandles(ArchiveMetadata), MCORE_LOCALINPUTGRANULEID,  &
                              inputGranuleIDs)
    
    doi_char =  "10.5067/MODIS/" // trim(EsdtShortName) // ".006"
    doiAuth = 'http://dx.doi.org'

    AttrN = MCORE_ADD_ATTRIBUTENAME// '.' // '1'
    call setMetadataAttribute(metadataHandles(InventoryMetadata), AttrN, &
         'identifier_product_doi')
    AttrN = MCORE_PARAMETERVALUE// '.' // '1'
    call setMetadataAttribute(metadataHandles(InventoryMetadata), AttrN, &
         doi_char)

    AttrN = MCORE_ADD_ATTRIBUTENAME// '.' // '2'
    call setMetadataAttribute(metadataHandles(InventoryMetadata), AttrN, &
         'identifier_product_doi_authority')
    AttrN = MCORE_PARAMETERVALUE// '.' // '2'
    call setMetadataAttribute(metadataHandles(InventoryMetadata), AttrN, &
         doiAuth)

    ! Write 'em.
    call writeMetadata(handle           = metadataHandles(InventoryMetadata), &
                       hdfAttributeName = InventoryMetadataHDFAttrName,       &
                       hdfFileID        = outputFileID)
    call writeMetadata(handle           = metadataHandles(  ArchiveMetadata), &
                       hdfAttributeName =   ArchiveMetadataHDFAttrName,       &
                       hdfFileID        = outputFileID)


    nRet = sfscatt(outputFileID, 'identifier_product_doi', DFNT_CHAR8, &
         26, doi_char) 
    if (nRet .eq. -1) then
       msgbuffer = 'Error in sfscatt Writing the global Attribute' // &
            ' identifier_product_doi in the output file.'
       call ErrorMessage("setOutputMetadata", msgbuffer)
    endif
    
    nRet = sfscatt(outputFileID, 'identifier_product_doi_authority', &
         DFNT_CHAR8, 17, doiAuth)
    if (nRet .eq. -1) then
       msgbuffer = 'Error in sfscatt Writing the global Attribute' // &
            ' identifier_product_doi_authority in the output file.'
       call ErrorMessage("setOutputMetadata", msgbuffer)
    endif
    ! ---- End setting metadata ---

    ! Free up remaining memory.
    deallocate(universalRefs,   stat = status(1))
    deallocate(inputGranuleIDs, stat = status(2))

    if(any(status /= 0)) then
       msgbuffer = "Can't allocate memory." &
  // char(10) // "Operator Action:  Suspect system problem.  If a fault is " &
  // char(10) // "identified, correct and rerun PGE.  Otherwise, notify SDST." 
       call ErrorMessage("setOutputMetadata", msgbuffer)
    end if

  end subroutine setOutputMetadata
! --------------------------------------------------------------------------------
end module RuntimeInfoL3Monthly
