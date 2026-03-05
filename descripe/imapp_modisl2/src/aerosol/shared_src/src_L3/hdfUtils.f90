module hdfUtils
  ! Module Revision History
  !   $Log: hdfUtils.f90,v $
  ! Revision 1.7  1998/02/13  16:40:08  pincus
  ! Separated functionality of subroutine getDataSetList into function
  ! getNumDataSets (which gets the number of data sets) and subroutine
  ! getDataSetList (which returns the names of the datasets).
  ! Added prologs to hidden specific procedures underlying generic procedures
  ! to conform with ECS requirements.
  !
  ! Revision 1.6  1997/12/30  15:36:00  pincus
  ! netCDF constants are now taken from the hdf module; removed reference
  ! to module netcdf. Names of the constants changed.
  !
  ! Revision 1.5  1997/11/20  21:53:48  pincus
  ! Removed explicit conversion of integer types to placate FORCHECK;
  ! removed unneeded variables.
  !
  ! Revision 1.4  1997/11/01  20:29:42  pincus
  ! Specified the entities used from module hdf.
  !
  ! Revision 1.3  1997/10/23  21:37:04  pincus
  ! Corrected check of argument sizes in vector versions of readAttribute.
  !
  ! Revision 1.2  1997/10/23  19:04:43  pincus
  ! Added function getAttributeLength, generic procedure readAttribute;
  ! changed endDataSet to subroutine.
  !
  ! Revision 1.1  1997/07/15  17:22:13  pincus
  ! Initial revision
  use typeSizes,     only: OneByteIntKind, TwoByteIntKind, FourByteIntKind, FourByteRealKind
  use ErrorHandling, only: FatalErrorMessage
  use hdf,           only: SFstart, SFend, SFselect, SFendacc, SFfattr, Sffinfo, &
                           SFn2index, SFginfo, SFgainfo, SFrcatt, SFrnatt,       &
                           DFACC_READ, DFACC_WRITE,                              &
                           DFNT_CHAR, DFNT_INT8, DFNT_INT16, DFNT_INT32, DFNT_FLOAT32, & 
                           FAIL,                                                 &
                           MAX_NC_NAME, MAX_VAR_DIMS
  implicit none
  private
  
  interface readAttribute
    module procedure readCharacterAttribute
    module procedure readByteAttributeScalar,     readByteAttributeVector
    module procedure readShortIntAttributeScalar, readShortIntAttributeVector
    module procedure readIntAttributeScalar,      readIntAttributeVector
    module procedure readRealAttributeScalar,     readRealAttributeVector
    ! !F90
    !
    ! !Description:
    !   This generic procedure has the following interface:
    !     subroutine readAttribute(dataSetID, attributeName, attributeValue)
    !        integer(kind = FourByteIntKind), intent( in) :: dataSetID
    !        character(len = *)               intent( in) :: attributeName
    !        *various*, *various*,            intent(out) :: attributeValue
    !     end subroutine readAttribute
    !
    !    Read the value of an attribute from an HDF file. Error trapping
    !      is internal to the subroutine. 
    !  
    ! !Input Parameters:
    !   dataSetID: the HDF ID of the file, array, or dimension from which to read the attribute
    !   attributeName: a character string holding the name of the attribute
    !
    ! !Output Parameters:
    !   attributeValue: the value of the attribute. This may be a character string, 
    !     a scalar or vector of one, two, or four byte integers or of four byte reals.
    !     If four byte integer(s) are supplied, one and two byte integer attributes are
    !     converted before output. 
    ! 
    ! !Revision History:
    !    See module revision history in data section of module. 
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
  end interface ! readAttribute
  
  public :: attributeExists, readAttribute
  public :: getDataSetID, endDataSet, &
            getNumDataSets, getDataSetList, getDataSetDimSizes, getDataSetType
  public :: getFillValue, getAttributeLength
  
contains
  ! ----------------------------------------------------------------------------------
  function getDataSetID(fileID, dataSetName)

  use CharacterUtils, only: IntToChar

    integer (kind = FourByteIntKind), intent( in) :: fileID
    character (len = *),              intent( in) :: dataSetName
    integer (kind = FourByteIntKind)              :: getDataSetID
    character (len = 2048)                        :: msgbuffer


    ! !F90
    !
    ! !Description:
    !   Get HDF dataset ID given HDF file ID and name of data set. Serves mostly 
    !     as an error-checking wrapper to HDF call SFselect. 
    !  
    ! !Input Parameters:
    !   fileID: HDF file ID, from a call to SFstart.
    !   dataSetName: character string containing name of data set to get ID. 
    !
    ! !Output Parameters:
    !    getDataSetID returns the HDF dataset ID of the dataset. If the ID can't be found 
    !      the function stops execution. 
    ! 
    ! !Revision History:
    !    See module revision history in data section of module. 
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
    !   Named constants:
    !     FAIL       (hdf.f90)
    !
    ! !END
    getDataSetID = SFselect(fileID, SFn2index(fileID, trim(dataSetName)))

    if(getDataSetID == FAIL) then
      msgbuffer = "Can't get dataSetID for dataSet " // &
                   trim(dataSetName) // " in file ID" // &
                   IntToChar(fileID) &
                // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("getDataSetID", msgbuffer)
    end if

  end function getDataSetID
  ! ----------------------------------------------------------------------------------
  subroutine endDataSet(dataSetID)
    integer (kind = FourByteIntKind), intent( in) :: dataSetID
    ! !F90
    !
    ! !Description:
    !   End access to an HDF SDS. This function is an error-checking wrapper to 
    !     HDF call SFselect. 
    !  
    ! !Input Parameters:
    !    dataSetID: HDF data set ID.
    !
    ! !Output Parameters:
    !    None. 
    ! 
    ! !Revision History:
    !    See module revision history in data section of module. 
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
    integer                       :: status
    character (len = 2048)        :: msgbuffer

    status = SFendacc(dataSetID)

    if(status == FAIL) then
      msgbuffer = "Error ending access to data set." &
           // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("endDataSet", msgbuffer)
    end if

  end subroutine endDataSet

  ! ----------------------------------------------------------------------------------
  function getNumDataSets(fileID)
    integer (kind = FourByteIntKind), intent( in) :: fileID
    integer                                       :: getNumDataSets
    ! !F90
    !
    ! !Description:
    !   Returns the number of SDSs in a file.
    !  
    ! !Input Parameters:
    !    fileID: the file handle as returned by the SD interface to the HDF library.
    !
    ! !Output Parameters:
    !    Function returns the number of SDSs defined through the SD interface.
    ! 
    ! !Revision History:
    !    See module revision history in data section of module. 
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
    integer (kind = FourByteIntKind) :: numAttributes
    integer                          :: status
    character (len = 2048)           :: msgbuffer

    status = SFfinfo(fileID, getNumDataSets, numAttributes)

    if (status == FAIL) then
      msgbuffer = "Can't get info for file. " &
        // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("getNumDataSets", msgbuffer)
    end if

  end function getNumDataSets
  ! ----------------------------------------------------------------------------------
  subroutine getDataSetList(fileID, dataSetNames)
    implicit none
    integer (kind = FourByteIntKind),  intent( in) :: fileID
    character (len = *), dimension(:), intent(out) :: dataSetNames
    ! !F90
    !
    ! !Description:
    !    Gets the number of parameters and the corresponding names from an HDF file. 
    !      The file is opened and closed within the subroutine. All error trapping
    !      is internal to the subroutine.
    !  
    ! !Input Parameters:
    !    fileID: the file handle as returned by the SD interface to the HDF library.
    !
    ! !Output Parameters:
    !    dataSetNames: a vector of character strings with the names of the parameters 
    !      in the file. Only the first numDataSets elements are set. 
    ! 
    ! !Revision History:
    !    See module revision history in data section of module. 
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
    !   Named constants:
    !     DFACC_READ (hdf.f90)
    !     FAIL       (hdf.f90)
    ! !END
    
    ! Local variables 
    integer (kind = FourByteIntKind), &
                 dimension(MAX_VAR_DIMS) :: dimensionIDs
    integer (kind = FourByteIntKind)     :: dataSet, dataSetID, &
                                            numDimensions, dataSetType, numAttributes
    integer                              :: status 
    character (len = 2048)               :: msgbuffer

    if (getNumDataSets(fileID) > size(dataSetNames)) then
      msgbuffer = "More datasets in file than room in array dataSetNames." &
                   // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("getDataSetList", msgbuffer)
    end if

    ! Loop over variables, get name of each.
    getDataSetInfo: do dataSet = 1, getNumDataSets(fileID) 
      dataSetID = SFselect(fileID, dataSet - 1)

    ! rhucek 02/05/08: Replaced undefined variable "start" below with 
    ! "dataSetID".  dataSetID is either -1 (FAIL) or a positive integer
    !
    ! if(status == FAIL) then 

      if(dataSetID == FAIL) then
        msgbuffer = "Can't get expected data set ID." &
              // char(10) // "Operator Action:  Notify SDST."
        call FatalErrorMessage("getDataSetList", msgbuffer)
      end if

      status = SFginfo(dataSetID, dataSetNames(dataSet), &
                       numDimensions, dimensionIDs, dataSetType, numAttributes)

      if(status == FAIL) then
        msgbuffer = "Failure getting data set info." &
            // char(10) // "Operator Action:  Notify SDST."
        call FatalErrorMessage("getDataSetList", msgbuffer)
      end if

      call endDataSet(dataSetID)
    end do getDataSetInfo

  end subroutine getDataSetList

! ----------------------------------------------------------------------------------  

  function getDataSetDimSizes(dataSetID, fileName, dataSetName)
    integer (kind = FourByteIntKind), intent( in) :: dataSetID
    character (len = *), optional,    intent( in) :: fileName, dataSetName
    integer (kind = FourByteIntKind), dimension(MAX_VAR_DIMS) :: getDataSetDimSizes
    ! !F90
    !
    ! !Description:
    !   Gets the dimension sizes (number of elements along each dimension) of an HDF SDS. 
    !  
    ! !Input Parameters:
    !   dataSetID: HDF data set ID, as returned from SFselect of getDataSetID. 
    !   fileName, dataSetName: optional character strings with the name of the file and
    !     the name of the data set. Used in error reporting.
    !
    ! !Output Parameters:
    !    getDataSetDimSizes returns a 4 byte integer array of length MAX_VAR_DIMS. The first elements
    !      contain the number of elements along the corresponding dimension. 
    ! 
    ! !Revision History:
    !    See module revision history in data section of module. 
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
    !    Named constants
    !      MAX_VAR_DIMS (netcdf.f90)
    !
    ! !END
    
    ! Local variables
    integer                                       :: status
    integer (kind = FourByteIntKind)              :: thisNumDims, thisDataSetType, &
                                                     numAttributes
    character(len = MAX_NC_NAME)                  :: thisDataSetName
    character (len = 2048)                        :: msgbuffer

    status =  SFginfo(dataSetID, thisDataSetName, thisNumDims, &
                      getDataSetDimSizes, thisDataSetType, numAttributes)
    if (status == FAIL) then
      if(present(fileName) .and. present(dataSetName)) then
        msgbuffer = "Can't get information from variable " // &
             trim(dataSetName) // " in file " // trim(fileName) &
             // char(10) // "Operator Action:  Notify SDST."
        call FatalErrorMessage("getDataSetDimSizes", msgbuffer)
      else
        msgbuffer = "Failed to get needed information from variable." &
                    // char(10) // "Operator Action:  Notify SDST."
        call FatalErrorMessage("getDataSetDimSizes", msgbuffer)
      end if
    end if
    getDataSetDimSizes(thisNumDims + 1 : ) = 0
        
  end function getDataSetDimSizes

! ----------------------------------------------------------------------------------
  function getDataSetType(dataSetID, fileName, dataSetName)
    integer (kind = FourByteIntKind), intent( in) :: dataSetID
    character (len = *), optional,    intent( in) :: fileName, dataSetName
    integer (kind = FourByteIntKind)              :: getDataSetType
    ! !F90
    !
    ! !Description:
    !   Returns the HDF data set type (i.e. DFNT_FLOAT32) of a given HDF data set. This function
    !     is an error checking wrapper for SFginfo. 
    !  
    ! !Input Parameters:
    !   dataSetID: HDF data set ID, as returned from SFselect of getDataSetID. 
    !   fileName, dataSetName: optional character strings with the name of the file and
    !     the name of the data set. Used in error reporting.
    !
    ! !Output Parameters:
    !    getDataSetType returns a scalar 4 byte integer with the data set type. The definitions 
    !      of the types are in the hdf module. 
    ! 
    ! !Revision History:
    !    See module revision history in data section of module. 
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
    !    Named constants
    !      MAX_NC_NAME, MAX_VAR_DIMS (hdf.f90)
    !
    ! !END
    
    ! Local variables
    integer                                       :: status
    integer (kind = FourByteIntKind)              :: thisNumDims, numAttributes
    character(len = MAX_NC_NAME)                  :: thisDataSetName
    integer (kind = FourByteIntKind), &
                              dimension(MAX_VAR_DIMS) :: thisDimensionSizes
    character (len = 2048)                        :: msgbuffer

    status =  SFginfo(dataSetID, thisDataSetName, thisNumDims, &
                      thisDimensionSizes, getDataSetType, numAttributes)
    if (status == FAIL) then
      if(present(fileName) .and. present(dataSetName)) then
        msgbuffer = "Failed to get needed information from variable " // &
               trim(dataSetName) // " in file " // trim(fileName) &
                 // char(10) // "Operator Action:  Notify SDST."
        call FatalErrorMessage("getDataSetType", msgbuffer)
      else
        msgbuffer = "Failed to get needed information from variable." &
                 // char(10) // "Operator Action:  Notify SDST."
        call FatalErrorMessage("getDataSetType", msgbuffer)
      end if
    end if
        
  end function getDataSetType
! ----------------------------------------------------------------------------------
  function getAttributeLength(dataSetID, attributeName)
    integer (kind = FourByteIntKind), intent( in) :: dataSetID
    character (len = *),              intent( in) :: attributeName
    integer                                       :: getAttributeLength
    ! !F90
    !
    ! !Description:
    !   Gets from an HDF attribute the number of value contained in the attribute. 
    !     Error checking included in the function.
    !
    ! !Input Parameters:
    !   dataSetID: the HDF ID of the file, array, or dimension from which to read the attribute
    !   attributeName: a character string holding the name of the attribute
    !
    ! !Output Parameters:
    !   Function value is the number of values in the attribute
    !
    ! !Revision History:
    !    See module revision history in data section of module. 
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
    integer                            :: status
    character(len = MAX_NC_NAME)       :: name
    integer(kind = FourByteIntKind)    :: attributeType, numValues
    character (len = 2048)             :: msgbuffer

    
    status    = SFgainfo(dataSetID, SFfattr(dataSetID, trim(attributeName)), &
                         name, attributeType, numValues)
    if(status == FAIL) then
      msgbuffer = "Can't get information for attribute " // trim(attributeName) &
                   // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("getAttributeLength", msgbuffer)
    end if
    getAttributeLength = numValues
    
  end function getAttributeLength
! ----------------------------------------------------------------------------------
  function attributeExists(dataSetID, attributeName)
    integer(kind = FourByteIntKind), intent( in) :: dataSetID
    character(len = *),              intent( in) :: attributeName
    logical                                      :: attributeExists
    ! !F90
    !
    ! !Description:
    !    
    !  
    ! !Input Parameters:
    !   dataSetID: the HDF ID of the file, array, or dimension from which to read the attribute
    !   attributeName: a character string holding the name of the attribute
    !
    ! !Output Parameters:
    !   The function returns .true. if the index to the attribute is found (i.e. 
    !     if the attribute exists) and .false. otherwise
    ! 
    ! !Revision History:
    !    See module revision history in data section of module. 
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
    
      
    !attributeExists = (SFfattr(dataSetID, trim(attributeName) == FAIL)
! richhucek:  modified below line to use '/=' rather than '=='
!    attributeExists = (SFfattr(dataSetID, trim(attributeName)) == FAIL)
     attributeExists = (SFfattr(dataSetID, trim(attributeName)) /= FAIL)

  end function attributeExists
! ----------------------------------------------------------------------------------
  subroutine readCharacterAttribute(dataSetID, attributeName, attributeValue)
    integer(kind = FourByteIntKind), intent( in) :: dataSetID
    character(len = *),              intent( in) :: attributeName
    character(len = *),              intent(out) :: attributeValue
    ! !F90
    !
    ! !Description:
    !    Read the value of an attribute from an HDF file. Error trapping
    !      is internal to the subroutine. 
    !  
    ! !Input Parameters:
    !   dataSetID: the HDF ID of the file, array, or dimension from which to read the attribute
    !   attributeName: a character string holding the name of the attribute
    !
    ! !Output Parameters:
    !   attributeValue: the value of the attribute; a character string.
    ! 
    ! !Revision History:
    !    See module revision history in data section of module. 
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
    integer                            :: status
    character(len = MAX_NC_NAME)       :: name
    integer(kind = FourByteIntKind)    :: attributeType, numValues
    character (len = 2048)             :: msgbuffer
    
    
    attributeValue = ""
    status    = SFgainfo(dataSetID, SFfattr(dataSetID, trim(attributeName)), &
                         name, attributeType, numValues)

    if(status == FAIL) then
      msgbuffer = "Can't read information about attribute " // trim(attributeName) &
                  // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("readCharacterAttribute", msgbuffer)
    end if

    if(numValues > len(attributeValue)) then
      msgbuffer = "Not enough space supplied to read character attribute " // &
                   trim(attributeName) &
                   // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("readCharacterAttribute", msgbuffer)
    end if

    if(attributeType /= DFNT_CHAR ) then
      msgbuffer = "Attribute " // trim(attributeName) // &
                  " is not of type character." &
                 // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("readCharacterAttribute", msgbuffer)
    end if

    status    = SFrcatt (dataSetID, SFfattr(dataSetID, trim(attributeName)), attributeValue)

    if(status == FAIL) then
      msgbuffer = "Can't read value of attribute " // trim(attributeName) &
                   // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("readCharacterAttribute", msgbuffer)
    end if

  end subroutine readCharacterAttribute
  ! ---------------------------------------------------------------------------------- 
  subroutine readByteAttributeScalar(dataSetID, attributeName, attributeValue)
    integer(kind = FourByteIntKind), intent( in) :: dataSetID
    character(len = *),              intent( in) :: attributeName
    integer(kind = OneByteIntKind),  intent(out) :: attributeValue
    ! !F90
    !
    ! !Description:
    !    Read the value of an attribute from an HDF file. Error trapping
    !      is internal to the subroutine. 
    !  
    ! !Input Parameters:
    !   dataSetID: the HDF ID of the file, array, or dimension from which to read the attribute
    !   attributeName: a character string holding the name of the attribute
    !
    ! !Output Parameters:
    !   attributeValue: the value of the attribute; a scalar one byte integer.
    ! 
    ! !Revision History:
    !    See module revision history in data section of module. 
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
    integer                            :: status
    character(len = MAX_NC_NAME)       :: name
    integer(kind = FourByteIntKind)    :: attributeType, numValues
    character (len = 2048)             :: msgbuffer


    
    status    = SFgainfo(dataSetID, SFfattr(dataSetID, trim(attributeName)), &
                         name, attributeType, numValues)

    if(status == FAIL) then
      msgbuffer = "Can't read information about attribute " // trim(attributeName) &
                   // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("readByteAttributeScalar", msgbuffer)
    end if

    if(numValues /= 1) then
      msgbuffer = "Attribute " //  trim(attributeName) //  &
                  "has more than one value." &
                  // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("readByteAttributeScalar", msgbuffer)
    end if

    if(attributeType /= DFNT_INT8) then
      msgbuffer = "Attribute " // trim(attributeName) //  &
                  " is not of type byte." &
                  // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("readByteAttributeScalar", msgbuffer)
    end if

    status    = SFrnatt (dataSetID, SFfattr(dataSetID, trim(attributeName)), attributeValue)

    if(status == FAIL) then
      msgbuffer = "Can't read value of attribute " // trim(attributeName) &
                  // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("readByteAttributeScalar", msgbuffer)
    end if

  end subroutine readByteAttributeScalar
  ! ----------------------------------------------------------------------------------
  subroutine readByteAttributeVector(dataSetID, attributeName, attributeValue)
    integer(kind = FourByteIntKind),               intent( in) :: dataSetID
    character(len = *),                            intent( in) :: attributeName
    integer(kind = OneByteIntKind),  dimension(:), intent(out) :: attributeValue
    ! !F90
    !
    ! !Description:
    !    Read the value of an attribute from an HDF file. Error trapping
    !      is internal to the subroutine. 
    !  
    ! !Input Parameters:
    !   dataSetID: the HDF ID of the file, array, or dimension from which to read the attribute
    !   attributeName: a character string holding the name of the attribute
    !
    ! !Output Parameters:
    !   attributeValue: the value of the attribute; a vector of one byte integers.
    ! 
    ! !Revision History:
    !    See module revision history in data section of module. 
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
    integer                            :: status
    character(len = MAX_NC_NAME)       :: name
    integer(kind = FourByteIntKind)    :: attributeType, numValues
    character (len = 2048)             :: msgbuffer

    
    status    = SFgainfo(dataSetID, SFfattr(dataSetID, trim(attributeName)), &
                         name, attributeType, numValues)
 
   if(numValues > size(attributeValue)) then
      msgbuffer = "Not enough space supplied for attribute " // &
                   trim(attributeName) &
                   // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("readByteAttributeVector", msgbuffer)
   end if

    if(attributeType /= DFNT_INT8) then
      msgbuffer = "Attribute " // trim(attributeName) // &
                  " is not of type byte. " & 
                   // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("readByteAttributeVector", msgbuffer)
    end if

    status    = SFrnatt (dataSetID, SFfattr(dataSetID, trim(attributeName)), attributeValue)

    if(status == FAIL) then
      msgbuffer = "Can't read value of attribute " // trim(attributeName) &
                   // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("readByteAttributeVector", msgbuffer)
    end if

  end subroutine readByteAttributeVector
  ! ---------------------------------------------------------------------------------- 
  subroutine readShortIntAttributeScalar(dataSetID, attributeName, attributeValue)
    integer(kind = FourByteIntKind), intent( in) :: dataSetID
    character(len = *),              intent( in) :: attributeName
    integer(kind = TwoByteIntKind),  intent(out) :: attributeValue
    ! !F90
    !
    ! !Description:
    !    Read the value of an attribute from an HDF file. Error trapping
    !      is internal to the subroutine. 
    !  
    ! !Input Parameters:
    !   dataSetID: the HDF ID of the file, array, or dimension from which to read the attribute
    !   attributeName: a character string holding the name of the attribute
    !
    ! !Output Parameters:
    !   attributeValue: the value of the attribute; a scalar two byte integer.
    ! 
    ! !Revision History:
    !    See module revision history in data section of module. 
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
    integer                            :: status
    character(len = MAX_NC_NAME)       :: name
    integer(kind = FourByteIntKind)    :: attributeType, numValues
    character (len = 2048)             :: msgbuffer

    
    status    = SFgainfo(dataSetID, SFfattr(dataSetID, trim(attributeName)), &
                         name, attributeType, numValues)

    if(status == FAIL) then
      msgbuffer = "Can't read information about attribute " // trim(attributeName) &
                   // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("readShortIntAttributeScalar", msgbuffer)
    end if

    if(numValues /= 1) then
      msgbuffer = "Attribute " //  trim(attributeName) //  &
                  " has more than one value." &
                  // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("readShortIntAttributeScalar", msgbuffer)
    end if

    if(attributeType /= DFNT_INT16) then
      msgbuffer = "Attribute " // trim(attributeName) //  &
                   " is not of type integer (two byte)." &
                   // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("readShortIntAttributeScalar", msgbuffer)
    end if

    status    = SFrnatt (dataSetID, SFfattr(dataSetID, trim(attributeName)), attributeValue)

    if(status == FAIL) then
      msgbuffer = "Can't read value of attribute " // trim(attributeName) &
                   // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("readShortIntAttributeScalar", msgbuffer)
    end if

  end subroutine readShortIntAttributeScalar
! ----------------------------------------------------------------------------------
  subroutine readShortIntAttributeVector(dataSetID, attributeName, attributeValue)
    integer(kind = FourByteIntKind),               intent( in) :: dataSetID
    character(len = *),                            intent( in) :: attributeName
    integer(kind = TwoByteIntKind),  dimension(:), intent(out) :: attributeValue
    ! !F90
    !
    ! !Description:
    !    Read the value of an attribute from an HDF file. Error trapping
    !      is internal to the subroutine. 
    !  
    ! !Input Parameters:
    !   dataSetID: the HDF ID of the file, array, or dimension from which to read the attribute
    !   attributeName: a character string holding the name of the attribute
    !
    ! !Output Parameters:
    !   attributeValue: the value of the attribute; a vector of two byte integers.
    ! 
    ! !Revision History:
    !    See module revision history in data section of module. 
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
    integer                            :: status
    character(len = MAX_NC_NAME)       :: name
    integer(kind = FourByteIntKind)    :: attributeType, numValues
    character (len = 2048)             :: msgbuffer

    
    status    = SFgainfo(dataSetID, SFfattr(dataSetID, trim(attributeName)), &
                         name, attributeType, numValues)

    if(status == FAIL) then
      msgbuffer = "Can't read information about attribute " // trim(attributeName) &
                   // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("readShortIntAttributeVector", msgbuffer)
    end if

    if(numValues > size(attributeValue)) then
      msgbuffer = "Not enough space supplied for attribute " // &
                   trim(attributeName) &
                   // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("readShortIntAttributeVector", msgbuffer)
    end if

    if(attributeType /= DFNT_INT16) then
      msgbuffer = "Attribute " // trim(attributeName) // &
                  " is not of type integer (two byte)." &
                  // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("readShortIntAttributeVector", msgbuffer)
    end if

    status    = SFrnatt (dataSetID, SFfattr(dataSetID, trim(attributeName)), attributeValue)

    if(status == FAIL) then
      msgbuffer = "Can't read value of attribute " // trim(attributeName) &
                  // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("readShortIntAttributeVector", msgbuffer)
    end if

  end subroutine readShortIntAttributeVector
! ----------------------------------------------------------------------------------
  subroutine readIntAttributeScalar(dataSetID, attributeName, attributeValue)
    integer(kind = FourByteIntKind), intent( in) :: dataSetID
    character(len = *),              intent( in) :: attributeName
    integer(kind = FourByteIntKind), intent(out) :: attributeValue
    ! !F90
    !
    ! !Description:
    !    Read the value of an attribute from an HDF file. Error trapping
    !      is internal to the subroutine. 
    !  
    ! !Input Parameters:
    !   dataSetID: the HDF ID of the file, array, or dimension from which to read the attribute
    !   attributeName: a character string holding the name of the attribute
    !
    ! !Output Parameters:
    !   attributeValue: the value of the attribute; a scalar four byte integer. Attributes of other
    !     integer types are converted to four byte integers.
    ! 
    ! !Revision History:
    !    See module revision history in data section of module. 
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
    integer                         :: status
    character(len = MAX_NC_NAME)       :: name
    integer(kind = FourByteIntKind) :: attributeType, numValues
    integer(kind =  OneByteIntKind) :: byteValue
    integer(kind =  TwoByteIntKind) :: shortValue
    character (len = 2048)          :: msgbuffer

    
    status    = SFgainfo(dataSetID, SFfattr(dataSetID, trim(attributeName)), &
                         name, attributeType, numValues)

    if(status == FAIL) then
      msgbuffer = "Can't read information about attribute " // trim(attributeName) &
                   // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("readIntAttributeScalar", msgbuffer)
    end if

    if(numValues /= 1) then
      msgbuffer = "Attribute " //  trim(attributeName) // &
                  " has more than one value." &
                  // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("readIntAttributeScalar", msgbuffer)
    end if

    select case (attributeType)
      case(DFNT_INT8 ) 
       status    = SFrnatt (dataSetID, SFfattr(dataSetID, trim(attributeName)), byteValue     )
       attributeValue = byteValue ! Type conversion
      case(DFNT_INT16) 
       status    = SFrnatt (dataSetID, SFfattr(dataSetID, trim(attributeName)), shortValue    )
       attributeValue = shortValue ! Type conversion
      case(DFNT_INT32)
       status    = SFrnatt (dataSetID, SFfattr(dataSetID, trim(attributeName)), attributeValue)
      case default

        msgbuffer = "Attribute " // trim(attributeName) // &
                    " is not of type integer." & 
                  // char(10) // "Operator Action:  Notify SDST."
        call FatalErrorMessage("readIntAttributeScalar", msgbuffer)

    end select 

    if(status == FAIL) then
      msgbuffer = "Can't read value of attribute " // trim(attributeName) &
                  // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("readIntAttributeScalar", msgbuffer)
    end if

  end subroutine readIntAttributeScalar
! ----------------------------------------------------------------------------------
  subroutine readIntAttributeVector(dataSetID, attributeName, attributeValue)
    integer(kind = FourByteIntKind),               intent( in) :: dataSetID
    character(len = *),                            intent( in) :: attributeName
    integer(kind = FourByteIntKind), dimension(:), intent(out) :: attributeValue
    ! !F90
    !
    ! !Description:
    !    Read the value of an attribute from an HDF file. Error trapping
    !      is internal to the subroutine. 
    !  
    ! !Input Parameters:
    !   dataSetID: the HDF ID of the file, array, or dimension from which to read the attribute
    !   attributeName: a character string holding the name of the attribute
    !
    ! !Output Parameters:
    !   attributeValue: the value of the attribute; a vector of four byte integers. Attributes of
    !     other integer types are converted to four byte integers.
    ! 
    ! !Revision History:
    !    See module revision history in data section of module. 
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
    integer                                                          :: status
    character(len = MAX_NC_NAME)                                     :: name
    integer(kind = FourByteIntKind)                                  :: attributeType, numValues
    integer(kind =  OneByteIntKind), dimension(size(attributeValue)) :: byteValues
    integer(kind =  TwoByteIntKind), dimension(size(attributeValue)) :: shortValues
    character (len = 2048)                                           :: msgbuffer

    
    status    = SFgainfo(dataSetID, SFfattr(dataSetID, trim(attributeName)), &
                         name, attributeType, numValues)

    if(status == FAIL) then
      msgbuffer = "Can't read information about attribute " // trim(attributeName) &
                  // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("readIntAttributeVector", msgbuffer)
    end if

    if(numValues > size(attributeValue)) then
      msgbuffer = "Not enough space supplied for attribute " // &
                   trim(attributeName) &
                   // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("readIntAttributeVector", msgbuffer)
    end if

    select case (attributeType)
      case(DFNT_INT8 ) 
       status    = SFrnatt (dataSetID, SFfattr(dataSetID, trim(attributeName)), byteValues    )
       attributeValue = byteValues(:) ! Type conversion
      case(DFNT_INT16) 
       status    = SFrnatt (dataSetID, SFfattr(dataSetID, trim(attributeName)), shortValues   )
       attributeValue = shortValues(:) ! Type conversion
      case(DFNT_INT32)
       status    = SFrnatt (dataSetID, SFfattr(dataSetID, trim(attributeName)), attributeValue)
      case default
        msgbuffer = "Attribute " // trim(attributeName) //  &
                    "is not of type integer" &
                    // char(10) // "Operator Action:  Notify SDST."
        call FatalErrorMessage("readIntAttributeVector", msgbuffer)
    end select 

    if(status == FAIL) then
      msgbuffer = "Can't read value of attribute " // trim(attributeName) &
                   // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("readIntAttributeVector", msgbuffer)
    end if

  end subroutine readIntAttributeVector
! ---------------------------------------------------------------------------------- 
  subroutine readRealAttributeScalar(dataSetID, attributeName, attributeValue)
    integer(kind = FourByteIntKind), intent( in) :: dataSetID
    character(len = *),              intent( in) :: attributeName
    real  (kind = FourByteRealKind), intent(out) :: attributeValue
    ! !F90
    !
    ! !Description:
    !    Read the value of an attribute from an HDF file. Error trapping
    !      is internal to the subroutine. 
    !  
    ! !Input Parameters:
    !   dataSetID: the HDF ID of the file, array, or dimension from which to read the attribute
    !   attributeName: a character string holding the name of the attribute
    !
    ! !Output Parameters:
    !   attributeValue: the value of the attribute; a scalar four byte real value.
    ! 
    ! !Revision History:
    !    See module revision history in data section of module. 
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
    integer                            :: status
    character(len = MAX_NC_NAME)       :: name
    integer(kind = FourByteIntKind)    :: attributeType, numValues
    character (len = 2048)             :: msgbuffer

    status    = SFgainfo(dataSetID, SFfattr(dataSetID, trim(attributeName)), &
                         name, attributeType, numValues)

    if(status == FAIL) then
      msgbuffer = "Can't read information about attribute " // trim(attributeName) &
                   // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("readRealAttributeScalar", msgbuffer)
    end if

    if(numValues /= 1) then
      msgbuffer = "Attribute " // trim(attributeName) // &
                  " has more than one value." &
                  // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("readRealAttributeScalar", msgbuffer)
    end if

    if(attributeType /= DFNT_FLOAT32) then
      msgbuffer = "Attribute " // trim(attributeName) // &
                  " is not of type real." & 
                  // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("readRealAttributeScalar", msgbuffer)
    end if

    status    = SFrnatt (dataSetID, SFfattr(dataSetID, trim(attributeName)), attributeValue)

    if(status == FAIL) then
      msgbuffer = "Can't read value of attribute " // trim(attributeName) &
                  // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("readRealAttributeScalar", msgbuffer)
    end if

  end subroutine readRealAttributeScalar
! ----------------------------------------------------------------------------------
  subroutine readRealAttributeVector(dataSetID, attributeName, attributeValue)
    integer(kind = FourByteIntKind),               intent( in) :: dataSetID
    character(len = *),                            intent( in) :: attributeName
    real  (kind = FourByteRealKind), dimension(:), intent(out) :: attributeValue
    ! !F90
    !
    ! !Description:
    !    Read the value of an attribute from an HDF file. Error trapping
    !      is internal to the subroutine. 
    !  
    ! !Input Parameters:
    !   dataSetID: the HDF ID of the file, array, or dimension from which to read the attribute
    !   attributeName: a character string holding the name of the attribute
    !
    ! !Output Parameters:
    !   attributeValue: the value of the attribute; a vector of four byte real values.
    ! 
    ! !Revision History:
    !    See module revision history in data section of module. 
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
    integer                            :: status
    character(len = MAX_NC_NAME)       :: name
    integer(kind = FourByteIntKind)    :: attributeType, numValues
    character (len = 2048)             :: msgbuffer

    
    status    = SFgainfo(dataSetID, SFfattr(dataSetID, trim(attributeName)), &
                         name, attributeType, numValues)

    if(status == FAIL) then
      msgbuffer = "Can't read information about attribute " // trim(attributeName) &
                   // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("readRealAttributeVector", msgbuffer)
    end if

    if(numValues > size(attributeValue)) then
      msgbuffer = "Not enough space supplied for attribute " // &
                   trim(attributeName) &
                   // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("readRealAttributeVector", msgbuffer)
    end if

    if(attributeType /= DFNT_FLOAT32) then
      msgbuffer = "Attribute " // trim(attributeName) //  &
                  "is not of type real." &
                   // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("readRealAttributeVector", msgbuffer)
    end if

    status    = SFrnatt (dataSetID, SFfattr(dataSetID, trim(attributeName)), attributeValue)

    if(status == FAIL) then
      msgbuffer = "Can't read value of attribute " // trim(attributeName) &
                   // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("readRealAttributeVector", msgbuffer)
    end if

  end subroutine readRealAttributeVector
! ---------------------------------------------------------------------------------- 
  subroutine getFillValue(dataSetID,  oneByteIntFillValue, &
                                      twoByteIntFillValue, &
                                     fourByteIntFillValue, &
                                     fourByteRealFillValue )
    integer (kind = FourByteIntKind ),           intent( in) :: dataSetID
    integer (kind =  OneByteIntKind ), optional, intent(out) ::  oneByteIntFillValue
    integer (kind =  TwoByteIntKind ), optional, intent(out) ::  twoByteIntFillValue
    integer (kind = FourByteIntKind ), optional, intent(out) :: fourByteIntFillValue
    real    (kind = FourByteRealKind), optional, intent(out) :: fourByteRealFillValue
    character (len = 2048)                                   :: msgbuffer

    ! !F90
    !
    ! !Description:
    !    Return the _FillValue attribute of an HDF data set. The data set type of the SDS 
    !      must be known before hand, and the correct input argument (only one of 
    !      oneByteIntFillValue, twoByteIntFillValue, fourByteIntFillValue, fourByteRealFillValue)
    !      supplied. 
    !  
    ! !Input Parameters:
    !   dataSetID: HDF data set ID, as returned from SFselect of getDataSetID. 
    !
    ! !Output Parameters:
    !   Only one of oneByteIntFillValue, twoByteIntFillValue, fourByteIntFillValue, 
    !     fourByteRealFillValue. The _FillValue attribute should have the same data type
    !     as the data set to which it applies. 
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
    !    Named constants
    !      MAX_VAR_DIMS (netcdf.f90)
    !
    ! !END

    if(count((/ present( oneByteIntFillValue), present( twoByteIntFillValue),          &
                present(fourByteIntFillValue), present(fourByteRealFillValue)/)) /= 1) then
      msgbuffer = "Must request exactly one fill value." &
                   // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("getFillValue", msgbuffer)
    end if

    if(present(oneByteIntFillValue)) &
      call readAttribute(dataSetID, "_FillValue", oneByteIntFillValue)

    if(present(twoByteIntFillValue)) &
      call readAttribute(dataSetID, "_FillValue", twoByteIntFillValue)

    if(present(fourByteIntFillValue)) &
      call readAttribute(dataSetID, "_FillValue", fourByteIntFillValue)

    if(present(fourByteRealFillValue)) &
      call readAttribute(dataSetID, "_FillValue", fourByteRealFillValue)

  end subroutine getFillValue
! ----------------------------------------------------------------------------------
end module hdfUtils
