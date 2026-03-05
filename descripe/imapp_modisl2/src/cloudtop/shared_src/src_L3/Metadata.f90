module Metadata
  ! Module Revision History
  !   $Log: Metadata.f90,v $
  !
  ! Revision 2.4 2005/11/07  R. Hucek
  !  Made the following 4 MAPI parameter names public:
  !  MCORE_PARAMETERNAME, MCORE_PERCENT_MISSING, MCORE_AUTO_QUALITY,
  !  and MCORE_AQUAL_FLG.
  ! 
  ! Revision 2.3 2000/02/22  Paul Hubanks
  !  expanded error message for core metadata read
  !
  ! Revision 2.2 1998/05/18  Liqun Ma 
  !  mapi.inc was replaced by mapi.f90.inc
  !
  ! Revision 2.1  1998/02/13  19:01:12  pincus
  !  Added functions for reading values from archive metadata and metadata that
  !  has been previously set (for example, in the MCF).
  !  Replaced explicit references to CoreMetadata.o and ArchiveMetadata.o with
  !  symbolic names.
  !  Added prologs to private specific functions underlying prublic generic functions
  !  to comply with ECS requirements.
  !
  ! Revision 2.0  1997/11/21  19:24:51  pincus
  !  No changes - incremented version number to indicate code
  !  as delivered to SDST.
  !
  ! Revision 1.3  1997/11/21  18:47:05  pincus
  !  Added support for writing out archive and inventory metadata fields.
  !
  ! Revision 1.1  1997/10/29  22:59:12  pincus
  !  Initial revision
  !
  use ErrorHandling,  only: ErrorMessage, FatalErrorMessage
  use CharacterUtils, only: IntToChar
  implicit none
  private
  ! Include files 
  ! include "PGS_MET.f"    ! Included in mapi.inc
  include "PGS_MET_13.f"
  include "PGS_SMF.f"
  include "mapi.f90.inc"
    
  ! Public constants
  character (len = len(PGSd_MET_STR_END)),  parameter, public :: &
                               ToolkitEndDataString = PGSd_MET_STR_END

  integer, parameter, public :: metadataGroupNameLength = PGSd_MET_GROUP_NAME_L,   &
                                metadataNumberOfGroups  = PGSd_MET_NUM_OF_GROUPS,  &
                                ECSParameterNameLength  = MAX_ECS_NAME_L,          &

  ! rhucek 12/04/03: added parameter ArchAttrStringLengthMax 
                                ArchAttrStringLengthMax = 255, & ! ECS B.0 Data Model 
                                MAX_ESDT_Name_LEN       = 8      ! From code by Rich Hucek.
  ! The following constants are included in mapi.inc; these lines make them available 
  ! outside this module. 
  public :: MCORE_EAST_BOUND,             MCORE_NORTH_BOUND,            &
            MCORE_SOUTH_BOUND,            MCORE_WEST_BOUND,             &
            MCORE_RANGE_BEG_DATE,         MCORE_RANGE_BEG_TIME,         &
            MCORE_RANGE_ENDING_DATE,      MCORE_RANGE_ENDING_TIME,      &
            MCORE_PARAMETERNAME,          MCORE_PERCENT_MISSING,        &
            MCORE_AUTO_QUALITY,           MCORE_AQUAL_FLG,              & 
            MCORE_DAYNIGHTFLAG,           MCORE_INPUT_POINTER,          &
            MCORE_LOCALINPUTGRANULEID,    MCORE_PGEVERSION,             &
            MCORE_SHORT_NAME
 
  character (len = ECSParameterNameLength), parameter, public :: &
            InventoryMetadataHDFAttrName = MECS_CORE,            &
            ArchiveMetadataHDFAttrName   = MECS_ARCHIVE

  ! Public functions
  public     :: getCoreMetadata, getArchiveMetadata, getUniversalReference, &
                initializeMetadata, setMetadataAttribute, getMetadataAttribute, writeMetadata
  
  interface getCoreMetadata
    module procedure getCoreMetadata_s, getCoreMetadata_s_v, getCoreMetadata_i, getCoreMetadata_d
  end interface ! getCoreMetadata
  ! Generic procedure 
  !   subroutine getCoreMetadata(LUN, version, ECSparamName, paramValue)
  !     integer,            intent( in) :: fileID
  !     integer,            intent( in) :: version
  !     character(len = *), intent( in) :: hdfAttrName, ECSParamName
  !     *various*,          intent(out) :: paramValue
  !  end subroutine getCoreMetadata
  ! !F90
  !
  ! !Description:
  !    Read a value from CoreMetadata.0. Error trapping included. 
  !  
  ! !Input Parameters:
  !   LUN: the Logical Unit Number of the file in the Process Control File (PCF)
  !   version: the instance of the LUN to read
  !   ECSparamName: the field within the HDf attribute "CoreMetadata.0" which is to be read.
  !
  ! !Output Parameters:
  !   paramValue: the value of the field within CoreMetadata.0. May be string, integer, or real
  !     valued. 
  ! 
  ! !Revision History:
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
  
  interface getArchiveMetadata
    module procedure getArchiveMetadata_s, getArchiveMetadata_s_v, &
                     getArchiveMetadata_i, getArchiveMetadata_d
  end interface ! getArchiveMetadata
  ! Generic procedure 
  !   subroutine getArchiveMetadata(LUN, version, ECSparamName, paramValue)
  !     integer,            intent( in) :: fileID
  !     integer,            intent( in) :: version
  !     character(len = *), intent( in) :: ECSParamName
  !     *various*,          intent(out) :: paramValue
  !  end subroutine getCoreMetadata
  ! !F90
  !
  ! !Description:
  !    Read a value from ArchiveMetadata.0. Error trapping included. 
  !  
  ! !Input Parameters:
  !   LUN: the Logical Unit Number of the file in the Process Control File (PCF)
  !   version: the instance of the LUN to read
  !   ECSparamName: the field within the HDf attribute "ArchiveMetadata.0" which is to be read.
  !
  ! !Output Parameters:
  !   paramValue: the value of the field within ArchiveMetadata.0. May be string, integer, or real
  !     valued. 
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
  
  interface setMetadataAttribute
    module procedure setMetadataAttribute_s,        setMetadataAttribute_i,        &
                     setMetadataAttribute_d
    module procedure setMetadataAttribute_s_scalar, setMetadataAttribute_i_scalar, &
                     setMetadataAttribute_d_scalar
  end interface ! setMetadataAttribute
  ! Generic procedure 
  !     subroutine setMetadataAttribute(handle, attributeName, attributeValue)
  !       character (len = *),           intent ( in) :: handle, attributeName
  !       *user defined*,      *varies*, intent ( in) :: attributeValue
  !  end subroutine setMetadataAttribute
  ! !F90
  !
  ! !Description:
  !    Set a metadata attribute after the MCF has been initialized. Error trapping included. 
  !  
  ! !Input Parameters:
  !   handle: the metadata handle, as retruned by PGS_MET_Init.
  !   attributeName: the ECS name of the attribute.
  !   attributeValue: the value to be assigned to the attribute. The underlying PGS 
  !     toolkit function takes one or multiple values based on entries in the PCF file. 
  !     There are separate toolkit routines for strings (_s), reals/doubles (_d), and 
  !     integers (_i). 
  !
  ! !Output Parameters:
  !   none.
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
  
  interface getMetadataAttribute
    module procedure getMetadataAttribute_s,        getMetadataAttribute_i,        &
                     getMetadataAttribute_d
    module procedure getMetadataAttribute_s_scalar, getMetadataAttribute_i_scalar, &
                     getMetadataAttribute_d_scalar
  end interface ! getMetadataAttribute
  ! Generic procedure 
  !     subroutine getMetadataAttribute(handle, attributeName, attributeValue)
  !       character (len = *),           intent ( in) :: handle, attributeName
  !       *user defined*,      *varies*, intent ( in) :: attributeValue
  !  end subroutine getMetadataAttribute
  ! !F90
  !
  ! !Description:
  !    Get the value of a metadata attribute after the MCF has been initialized. 
  !    Error trapping included. 
  !  
  ! !Input Parameters:
  !   handle: the metadata handle, as retruned by PGS_MET_Init.
  !   attributeName: the ECS name of the attribute.
  !   attributeValue: the value to be assigned to the attribute. The underlying PGS 
  !     toolkit function takes one or multiple values based on entries in the PCF file. 
  !     There are separate toolkit routines for strings (_s), reals/doubles (_d), and 
  !     integers (_i). 
  !
  ! !Output Parameters:
  !   none.
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
  
  ! Interface definitions copied from the SDP toolkit user's guide, version 
  !   5.2. Putting the interface in the module makes it harder to make 
  !   mistakes when calling toolkit functions. 
  interface 
    ! Runtime attribute functions
    integer function PGS_MET_GetPCAttr_s(productID, versionNumber, hdfAttrName, &
                                         ECSParamName, paramValue)
      integer,                          intent( in) :: productID
      integer,                          intent( in) :: versionNumber
      character(len = *),               intent( in) :: hdfAttrName, ECSParamName
      character(len = *), dimension(*), intent(out) :: paramValue
    end function PGS_MET_GetPCAttr_s

    integer function PGS_MET_GetPCAttr_i(productID, versionNumber, hdfAttrName, &
                                         ECSParamName, paramValue)
      integer,            intent( in) :: productID
      integer,            intent( in) :: versionNumber
      character(len = *), intent( in) :: hdfAttrName, ECSParamName
      integer,            intent(out) :: paramValue
    end function PGS_MET_GetPCAttr_i

    integer function PGS_MET_GetPCAttr_d(productID, versionNumber, hdfAttrName, &
                                         ECSParamName, paramValue)
      integer,            intent( in) :: productID
      integer,            intent( in) :: versionNumber
      character(len = *), intent( in) :: hdfAttrName, ECSParamName
      real,               intent(out) :: paramValue
    end function PGS_MET_GetPCAttr_d
    
    ! Universal reference function
    integer function PGS_PC_GetUniversalRef(productID, versionNumber, &
                                            universalReference)
      integer,            intent(   in) :: productID
      integer,            intent(inout) :: versionNumber
      character(len = *), intent(  out) :: universalReference
    end function PGS_PC_GetUniversalRef
    
    ! Metadata functions
    integer function PGS_MET_Init(fileID, handles)
      include "PGS_MET.f"
      integer,                                  intent( in) :: fileID
      character (len = PGSd_MET_GROUP_NAME_L), &
            dimension(PGSd_MET_NUM_OF_GROUPS),  intent(out) :: handles
    end function PGS_MET_Init
    
    ! Functions for setting metadata attributes
    integer function PGS_MET_SetAttr_s(handle, attributeName, attributeValue)
      character (len = *),               intent ( in) :: handle, attributeName
      character (len = *), dimension(*), intent ( in) :: attributeValue
    end function PGS_MET_SetAttr_s
    
    integer function PGS_MET_SetAttr_i(handle, attributeName, attributeValue)
      character (len = *),               intent ( in) :: handle, attributeName
      integer,             dimension(*), intent ( in) :: attributeValue
    end function PGS_MET_SetAttr_i
    
    integer function PGS_MET_SetAttr_d(handle, attributeName, attributeValue)
      character (len = *),               intent ( in) :: handle, attributeName
      double precision,    dimension(*), intent ( in) :: attributeValue
    end function PGS_MET_SetAttr_d
    
    ! Functions for reading (getting) metadata attributes
    integer function PGS_MET_GetSetAttr_s(handle, attributeName, attributeValue)
      character (len = *),               intent ( in) :: handle, attributeName
      character (len = *), dimension(*), intent ( in) :: attributeValue
    end function PGS_MET_GetSetAttr_s
    
    integer function PGS_MET_GetSetAttr_i(handle, attributeName, attributeValue)
      character (len = *),               intent ( in) :: handle, attributeName
      integer,             dimension(*), intent ( in) :: attributeValue
    end function PGS_MET_GetSetAttr_i
    
    integer function PGS_MET_GetSetAttr_d(handle, attributeName, attributeValue)
      character (len = *),               intent ( in) :: handle, attributeName
      double precision,    dimension(*), intent ( in) :: attributeValue
    end function PGS_MET_GetSetAttr_d
    
   integer function PGS_MET_Write(handle, hdfAttributeName, hdfFileID)
      character (len = *), intent( in) :: handle, hdfAttributeName
      integer,             intent( in) :: hdfFileID
    end function PGS_MET_Write
  end interface 
  
contains
  !  --------------------------------------------------------------------
  subroutine getCoreMetadata_s(fileID, version, ECSparamName, paramValue)
    integer,            intent( in) :: fileID
    integer,            intent( in) :: version
    character(len = *), intent( in) :: ECSParamName
    character(len = *), intent(out) :: paramValue
    ! !F90
    !
    ! !Description:
    !    Read a value from CoreMetadata.0. Error trapping included. 
    !  
    ! !Input Parameters:
    !   LUN: the Logical Unit Number of the file in the Process Control File (PCF)
    !   version: the instance of the LUN to read
    !   ECSparamName: the field within the HDf attribute "CoreMetadata.0" which is to be read.
    !
    ! !Output Parameters:
    !   paramValue: the value of the field within CoreMetadata.0, reported as a character string
    ! 
    ! !Revision History:
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
    integer                                         :: status
    character (len = len(paramValue)), dimension(1) :: localString
    character (len = 2048)                          :: msgbuffer

    status = PGS_MET_GetPCAttr_s(fileID, version, InventoryMetadataHDFAttrName, &
                                 trim(ECSparamName), localString)

    if (status /= PGS_S_SUCCESS) then
      msgbuffer = "Error reading string core metadata " // trim(ECSparamName) &
                  // char(10) // "from FileID = " // trim(IntToChar(fileID)) &
                  // " FileVersion = " // trim(IntToChar(version)) &
                  // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("getCoreMetadata_s", msgbuffer)
    end if

    paramValue = localString(1)                         
  end subroutine getCoreMetadata_s
  !  --------------------------------------------------------------------
  subroutine getCoreMetadata_s_v(fileID, version, ECSparamName, paramValue)
    integer,                          intent( in) :: fileID
    integer,                          intent( in) :: version
    character(len = *),               intent( in) :: ECSParamName
    character(len = *), dimension(:), intent(out) :: paramValue
    ! !F90
    !
    ! !Description:
    !    Read a value from CoreMetadata.0. Error trapping included. 
    !  
    ! !Input Parameters:
    !   LUN: the Logical Unit Number of the file in the Process Control File (PCF)
    !   version: the instance of the LUN to read
    !   ECSparamName: the field within the HDf attribute "CoreMetadata.0" which is to be read.
    !
    ! !Output Parameters:
    !   paramValue: the value of the field within CoreMetadata.0, reported as a vector of 
    !     character strings
    ! 
    ! !Revision History:
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
    integer                           :: status
    character (len = 2048)            :: msgbuffer

    status = PGS_MET_GetPCAttr_s(fileID, version, InventoryMetadataHDFAttrName, &
                                 trim(ECSparamName), paramValue)

    if (status /= PGS_S_SUCCESS) then
      msgbuffer = "Error reading string core metadata " // trim(ECSparamName) &
                  // char(10) // "from FileID = " // trim(IntToChar(fileID)) &
                  // " FileVersion = " // trim(IntToChar(version)) &
                  // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("getCoreMetadata_s_v", msgbuffer)
    end if

  end subroutine getCoreMetadata_s_v
  ! --------------------------------------------------------------------
  subroutine getCoreMetadata_i(fileID, version, ECSparamName, paramValue)
    integer,            intent( in) :: fileID
    integer,            intent( in) :: version
    character(len = *), intent( in) :: ECSParamName
    integer,            intent(out) :: paramValue
    ! !F90
    !
    ! !Description:
    !    Read a value from CoreMetadata.0. Error trapping included. 
    !  
    ! !Input Parameters:
    !   LUN: the Logical Unit Number of the file in the Process Control File (PCF)
    !   version: the instance of the LUN to read
    !   ECSparamName: the field within the HDf attribute "CoreMetadata.0" which is to be read.
    !
    ! !Output Parameters:
    !   paramValue: the value of the field within CoreMetadata.0, reported as a (default) integer.
    ! 
    ! !Revision History:
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
    integer                          :: status
    character (len = 2048)           :: msgbuffer

    status = PGS_MET_GetPCAttr_i(fileID, version, InventoryMetadataHDFAttrName, &
                                 trim(ECSparamName), paramValue)

    if (status /= PGS_S_SUCCESS) then
      msgbuffer = "Error reading string core metadata " // trim(ECSparamName) &
                  // char(10) // "from FileID = " // trim(IntToChar(fileID)) &
                  // " FileVersion = " // trim(IntToChar(version)) &
                  // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("getCoreMetadata_i", msgbuffer)
    end if

  end subroutine getCoreMetadata_i
  ! --------------------------------------------------------------------
  subroutine getCoreMetadata_d(fileID, version, ECSparamName, paramValue)
    integer,            intent( in) :: fileID
    integer,            intent( in) :: version
    character(len = *), intent( in) :: ECSParamName
    real,               intent(out) :: paramValue
    ! !F90
    !
    ! !Description:
    !    Read a value from CoreMetadata.0. Error trapping included. 
    !  
    ! !Input Parameters:
    !   LUN: the Logical Unit Number of the file in the Process Control File (PCF)
    !   version: the instance of the LUN to read
    !   ECSparamName: the field within the HDf attribute "CoreMetadata.0" which is to be read.
    !
    ! !Output Parameters:
    !   paramValue: the value of the field within CoreMetadata.0, reported as a real value. 
    ! 
    ! !Revision History:
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
    integer                              :: status
    character (len = 2048)               :: msgbuffer

    status = PGS_MET_GetPCAttr_d(fileID, version, InventoryMetadataHDFAttrName, &
                                 trim(ECSparamName), paramValue)

    if (status /= PGS_S_SUCCESS) then
      msgbuffer = "Error reading string core metadata " // trim(ECSparamName) &
                  // char(10) // "from FileID = " // trim(IntToChar(fileID)) &
                  // " FileVersion = " // trim(IntToChar(version)) &
                  // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("getCoreMetadata_d", msgbuffer)
    end if

  end subroutine getCoreMetadata_d
  !  --------------------------------------------------------------------
   subroutine getArchiveMetadata_s(fileID, version, ECSparamName, paramValue)
    integer,            intent( in) :: fileID
    integer,            intent( in) :: version
    character(len = *), intent( in) :: ECSParamName
    character(len = *), intent(out) :: paramValue
    ! !F90
    !
    ! !Description:
    !    Read a value from ArchiveMetadata.0. Error trapping included. 
    !  
    ! !Input Parameters:
    !   LUN: the Logical Unit Number of the file in the Process Control File (PCF)
    !   version: the instance of the LUN to read
    !   ECSparamName: the field within the HDf attribute "ArchiveMetadata.0" which is to be read.
    !
    ! !Output Parameters:
    !   paramValue: the value of the field within ArchiveMetadata.0, reported as a character string.
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
  
    ! Local variables
    integer                                         :: status
    character (len = 2048)                          :: msgbuffer
    character (len = len(paramValue)), dimension(1) :: localString
    
    status = PGS_MET_GetPCAttr_s(fileID, version, ArchiveMetadataHDFAttrName, &
                                 trim(ECSparamName), localString)

    if (status /= PGS_S_SUCCESS) then
      msgbuffer = "Error reading string Archive metadata " // trim(ECSparamName) &
                   // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("getArchiveMetadata_s", msgbuffer)
    end if

    paramValue = localString(1)                         
  end subroutine getArchiveMetadata_s
 ! --------------------------------------------------------------------
  subroutine getArchiveMetadata_s_v(fileID, version, ECSparamName, paramValue)
    integer,                          intent( in) :: fileID
    integer,                          intent( in) :: version
    character(len = *),               intent( in) :: ECSParamName
    character(len = *), dimension(:), intent(out) :: paramValue
    ! !F90
    !
    ! !Description:
    !    Read a value from ArchiveMetadata.0. Error trapping included. 
    !  
    ! !Input Parameters:
    !   LUN: the Logical Unit Number of the file in the Process Control File (PCF)
    !   version: the instance of the LUN to read
    !   ECSparamName: the field within the HDf attribute "ArchiveMetadata.0" which is to be read.
    !
    ! !Output Parameters:
    !   paramValue: the value of the field within ArchiveMetadata.0, reported as a vector of
    !     character string.
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
  
    ! Local variables
    integer                              :: status
    character (len = 2048)               :: msgbuffer

    status = PGS_MET_GetPCAttr_s(fileID, version, ArchiveMetadataHDFAttrName, &
                                 trim(ECSparamName), paramValue)

    if (status /= PGS_S_SUCCESS) then
      msgbuffer = "Error reading string Archive metadata " // trim(ECSparamName) &
                   // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("getArchiveMetadata_s_v", msgbuffer)
    end if

  end subroutine getArchiveMetadata_s_v
  ! --------------------------------------------------------------------
  subroutine getArchiveMetadata_i(fileID, version, ECSparamName, paramValue)
    integer,            intent( in) :: fileID
    integer,            intent( in) :: version
    character(len = *), intent( in) :: ECSParamName
    integer,            intent(out) :: paramValue
    ! !F90
    !
    ! !Description:
    !    Read a value from ArchiveMetadata.0. Error trapping included. 
    !  
    ! !Input Parameters:
    !   LUN: the Logical Unit Number of the file in the Process Control File (PCF)
    !   version: the instance of the LUN to read
    !   ECSparamName: the field within the HDf attribute "ArchiveMetadata.0" which is to be read.
    !
    ! !Output Parameters:
    !   paramValue: the value of the field within ArchiveMetadata.0, reported as a default integer.
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
  
    ! Local variables
    integer                              :: status
    character (len = 2048)               :: msgbuffer

    status = PGS_MET_GetPCAttr_i(fileID, version, ArchiveMetadataHDFAttrName, &
                                 trim(ECSparamName), paramValue)
    if (status /= PGS_S_SUCCESS) then
      msgbuffer = "Error reading integer Archive metadata " // trim(ECSparamName) &
                   // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("getArchiveMetadata_i", msgbuffer)
    end if

  end subroutine getArchiveMetadata_i
  ! --------------------------------------------------------------------
  subroutine getArchiveMetadata_d(fileID, version, ECSparamName, paramValue)
    integer,            intent( in) :: fileID
    integer,            intent( in) :: version
    character(len = *), intent( in) :: ECSParamName
    real,               intent(out) :: paramValue
    ! !F90
    !
    ! !Description:
    !    Read a value from ArchiveMetadata.0. Error trapping included. 
    !  
    ! !Input Parameters:
    !   LUN: the Logical Unit Number of the file in the Process Control File (PCF)
    !   version: the instance of the LUN to read
    !   ECSparamName: the field within the HDf attribute "ArchiveMetadata.0" which is to be read.
    !
    ! !Output Parameters:
    !   paramValue: the value of the field within ArchiveMetadata.0, reported as a real value.
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
  
    ! Local variables
    integer                              :: status
    character (len = 2048)               :: msgbuffer

    status = PGS_MET_GetPCAttr_d(fileID, version, ArchiveMetadataHDFAttrName, &
                                 trim(ECSparamName), paramValue)
    if (status /= PGS_S_SUCCESS) then
      msgbuffer = "Error reading real Archive metadata " // trim(ECSparamName) &
                   // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("getArchiveMetadata_d", msgbuffer)
    end if

  end subroutine getArchiveMetadata_d
  ! --------------------------------------------------------------------
  subroutine getUniversalReference(productID, versionNumber, universalReference)
    integer,            intent( in) :: productID
    integer,            intent( in) :: versionNumber
    character(len = *), intent(out) :: universalReference
    ! !F90
    !
    ! !Description:
    !    Read a UR (Universal Reference) given a version and file number. 
    !  
    ! !Input Parameters:
    !   productID: the Logical Unit Number of the file in the Process Control File (PCF)
    !   version: the instance of the LUN to read
    !
    ! !Output Parameters:
    !   universalReference: string containing the UR of the file in question
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
    !   This is esentially an error-trapping wrapper for PGS_PC_GetUniversalReference. 
    !
    ! !END
    
    ! Local variable
    integer                              :: status, tempVersion
    character (len = 2048)               :: msgbuffer

    ! We pass a local copy of the version number to PGS_PC_GetUniversalRef because that 
    !   function changes the value on return. 
    tempVersion = versionNumber
    status = PGS_PC_GetUniversalRef(productID, tempVersion, universalReference)

    if( status /= PGS_S_SUCCESS) then
      msgbuffer = "Can't get universal reference for product ID " // &
                  trim(IntToChar(productID)) &
                 // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("getUniversalReference", msgbuffer)
    end if

  end subroutine getUniversalReference
  ! --------------------------------------------------------------------
  subroutine initializeMetadata(fileID, handles)
    integer,                            intent( in) :: fileID
    ! character (len = *), dimension(:),  intent(out) :: handles
    ! whuang 09/18/01: character len is redefined to PGSd_MET_GROUP_NAME_L instead
    !                  of *
    character (len = PGSd_MET_GROUP_NAME_L), dimension(:),  intent(out) :: handles
    ! !F90
    !
    ! !Description:
    !    Initialize the MCF data. Error trapping included.
    !  
    ! !Input Parameters:
    !   fileID: the HDF file ID as returned by the SD interface.
    !
    ! !Output Parameters:
    !   handles: the list of metasdata handles
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
    !   Input arguments are identical to those for PGS_MET_Init.
    !
    ! !END

    ! Local variables
    integer                              :: status
    character (len = 2048)               :: msgbuffer

    ! Begin with some error checking.

    if(size(handles, 1) < PGSd_MET_NUM_OF_GROUPS) then
      msgbuffer = "Must provide at least " // &
                   trim(IntToChar(PGSd_MET_NUM_OF_GROUPS)) // " handles." &
                   // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("initializeMetadata", msgbuffer)
    end if

    if(len(handles(1)) < PGSd_MET_GROUP_NAME_L) then
      msgbuffer = "Handles must be at least " // &
                   trim(IntToChar(PGSd_MET_GROUP_NAME_L)) // " characters long." &
                   // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("initializeMetadata", msgbuffer)
    end if
                             
    ! Initialize the metadata
    status = PGS_MET_Init(fileID, handles)

    if(status /= PGS_S_SUCCESS) then
      msgbuffer = "Can't initialize metedata. Return code " // &
                   trim(IntToChar(status)) &
                   // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("initializeMetadata", msgbuffer)
    end if

  end subroutine initializeMetadata
  ! --------------------------------------------------------------------
  
  ! --------------------------------------------------------------------
  ! Functions for getting a setting a metadata attribute.
  ! --------------------------------------------------------------------
  subroutine setMetadataAttribute_s(handle, attributeName, attributeValue)
    character (len = *),               intent ( in) :: handle, attributeName
    character (len = *), dimension(*), intent ( in) :: attributeValue
    character (len = 2048)                          :: msgbuffer

    ! !F90
    !
    ! !Description:
    !    Set a metadata attribute after the MCF has been initialized. Error trapping included. 
    !  
    ! !Input Parameters:
    !   handle: the metadata handle, as retruned by PGS_MET_Init.
    !   attributeName: the ECS name of the attribute.
    !   attributeValue: the value to be assigned to the attribute 
    !
    ! !Output Parameters:
    !   none.
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
    
    if(PGS_MET_SetAttr_s(handle, attributeName, attributeValue) /= PGS_S_SUCCESS) then
      msgbuffer = "Error setting metadata attribute " &
                  // trim(attributeName) &
                  // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("setMetadataAttribute_s", msgbuffer)
    end if

  end subroutine setMetadataAttribute_s
  ! --------------------------------------------------------------------
  subroutine setMetadataAttribute_i(handle, attributeName, attributeValue)
    character (len = *),               intent ( in) :: handle, attributeName
    integer,             dimension(*), intent ( in) :: attributeValue

    ! !F90
    !
    ! !Description:
    !    Set a metadata attribute after the MCF has been initialized. Error trapping included. 
    !  
    ! !Input Parameters:
    !   handle: the metadata handle, as retruned by PGS_MET_Init.
    !   attributeName: the ECS name of the attribute.
    !   attributeValue: the value to be assigned to the attribute 
    !
    ! !Output Parameters:
    !   none.
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
    
    character (len = 2048)               :: msgbuffer

    if(PGS_MET_SetAttr_i(handle, attributeName, attributeValue) /= PGS_S_SUCCESS) then
      msgbuffer = "Error setting metadata attribute " &
                  // trim(attributeName) &
                  // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("setMetadataAttribute_i", msgbuffer)
    end if

  end subroutine setMetadataAttribute_i
  ! --------------------------------------------------------------------
  subroutine setMetadataAttribute_d(handle, attributeName, attributeValue)
    character (len = *),               intent ( in) :: handle, attributeName
    double precision,    dimension(*), intent ( in) :: attributeValue
    character (len = 2048)                          :: msgbuffer

    ! !F90
    !
    ! !Description:
    !    Set a metadata attribute after the MCF has been initialized. Error trapping included. 
    !  
    ! !Input Parameters:
    !   handle: the metadata handle, as retruned by PGS_MET_Init.
    !   attributeName: the ECS name of the attribute.
    !   attributeValue: the value to be assigned to the attribute 
    !
    ! !Output Parameters:
    !   none.
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
    
    if(PGS_MET_SetAttr_d(handle, attributeName, attributeValue) /= PGS_S_SUCCESS) then
      msgbuffer = "Error setting metadata attribute " &
                   // trim(attributeName) &
                   // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("setMetadataAttribute_d", msgbuffer)
    end if

  end subroutine setMetadataAttribute_d
  ! --------------------------------------------------------------------
  subroutine setMetadataAttribute_s_scalar(handle, attributeName, attributeValue)
    character (len = *), intent ( in) :: handle, attributeName
    character (len = *), intent ( in) :: attributeValue

    ! !F90
    !
    ! !Description:
    !    Set a metadata attribute after the MCF has been initialized. Error trapping included. 
    !  
    ! !Input Parameters:
    !   handle: the metadata handle, as retruned by PGS_MET_Init.
    !   attributeName: the ECS name of the attribute.
    !   attributeValue: the value to be assigned to the attribute 
    !
    ! !Output Parameters:
    !   none.
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
    
    call setMetadataAttribute_s(handle, attributeName, (/ attributeValue /) )
  end subroutine setMetadataAttribute_s_scalar
  ! --------------------------------------------------------------------
  subroutine setMetadataAttribute_i_scalar(handle, attributeName, attributeValue)
    character (len = *), intent ( in) :: handle, attributeName
    integer,             intent ( in) :: attributeValue
    ! !F90
    !
    ! !Description:
    !    Set a metadata attribute after the MCF has been initialized. Error trapping included. 
    !  
    ! !Input Parameters:
    !   handle: the metadata handle, as retruned by PGS_MET_Init.
    !   attributeName: the ECS name of the attribute.
    !   attributeValue: the value to be assigned to the attribute 
    !
    ! !Output Parameters:
    !   none.
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
    
    call setMetadataAttribute_i(handle, attributeName, (/ attributeValue /) )
  end subroutine setMetadataAttribute_i_scalar
  ! --------------------------------------------------------------------
  subroutine setMetadataAttribute_d_scalar(handle, attributeName, attributeValue)
    character (len = *), intent ( in) :: handle, attributeName
    double precision,    intent ( in) :: attributeValue
    ! !F90
    !
    ! !Description:
    !    Set a metadata attribute after the MCF has been initialized. Error trapping included. 
    !  
    ! !Input Parameters:
    !   handle: the metadata handle, as retruned by PGS_MET_Init.
    !   attributeName: the ECS name of the attribute.
    !   attributeValue: the value to be assigned to the attribute 
    !
    ! !Output Parameters:
    !   none.
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
    
    call setMetadataAttribute_d(handle, attributeName, (/ attributeValue /) )
  end subroutine setMetadataAttribute_d_scalar
  ! --------------------------------------------------------------------
  
  ! --------------------------------------------------------------------
  ! Functions for getting a previously set metadata attribute.
  ! --------------------------------------------------------------------
  subroutine getMetadataAttribute_s(handle, attributeName, attributeValue)
    character (len = *),               intent ( in) :: handle, attributeName
    character (len = *), dimension(*), intent (out) :: attributeValue
    character (len = 2048)                          :: msgbuffer

    ! !F90
    !
    ! !Description:
    !    Get the value of a metadata attribute after the MCF has been initialized. 
    !    Error trapping included. 
    !  
    ! !Input Parameters:
    !   handle: the metadata handle, as retruned by PGS_MET_Init.
    !   attributeName: the ECS name of the attribute.
    !   attributeValue: the value(s) to be assigned to the attribute. 
    !
    ! !Output Parameters:
    !   none.
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
    
    if(PGS_MET_GetSetAttr_s(handle, attributeName, attributeValue) /= PGS_S_SUCCESS) then
      msgbuffer = "Error setting metadata attribute " &
                   // trim(attributeName) &
                   // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("getMetadataAttribute_s", msgbuffer)
    end if

  end subroutine getMetadataAttribute_s
  ! --------------------------------------------------------------------
  subroutine getMetadataAttribute_i(handle, attributeName, attributeValue)
    character (len = *),               intent ( in) :: handle, attributeName
    integer,             dimension(*), intent (out) :: attributeValue
    character (len = 2048)                          :: msgbuffer

    ! !F90
    !
    ! !Description:
    !    Get the value of a metadata attribute after the MCF has been initialized. 
    !    Error trapping included. 
    !  
    ! !Input Parameters:
    !   handle: the metadata handle, as retruned by PGS_MET_Init.
    !   attributeName: the ECS name of the attribute.
    !   attributeValue: the value(s) to be assigned to the attribute. 
    !
    ! !Output Parameters:
    !   none.
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
    
    if(PGS_MET_GetSetAttr_i(handle, attributeName, attributeValue) /= PGS_S_SUCCESS) then
      msgbuffer = "Error setting metadata attribute " &
                   // trim(attributeName) &
                   // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("getMetadataAttribute_i", msgbuffer)
    end if

  end subroutine getMetadataAttribute_i
  ! --------------------------------------------------------------------
  subroutine getMetadataAttribute_d(handle, attributeName, attributeValue)
    character (len = *),               intent ( in) :: handle, attributeName
    double precision,    dimension(*), intent (out) :: attributeValue
    character (len = 2048)                          :: msgbuffer

    ! !F90
    !
    ! !Description:
    !    Get the value of a metadata attribute after the MCF has been initialized. 
    !    Error trapping included. 
    !  
    ! !Input Parameters:
    !   handle: the metadata handle, as retruned by PGS_MET_Init.
    !   attributeName: the ECS name of the attribute.
    !   attributeValue: the value(s) to be assigned to the attribute. 
    !
    ! !Output Parameters:
    !   none.
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
    
    if(PGS_MET_GetSetAttr_d(handle, attributeName, attributeValue) /= PGS_S_SUCCESS) then
      msgbuffer = "Error setting metadata attribute " &
                   // trim(attributeName) &
                   // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("getMetadataAttribute_d", msgbuffer)
    end if

  end subroutine getMetadataAttribute_d
  ! --------------------------------------------------------------------
  subroutine getMetadataAttribute_s_scalar(handle, attributeName, attributeValue)
    character (len = *), intent ( in) :: handle, attributeName
    character (len = *), intent (out) :: attributeValue
    ! !F90
    !
    ! !Description:
    !    Get the value of a metadata attribute after the MCF has been initialized. 
    !    Error trapping included. 
    !  
    ! !Input Parameters:
    !   handle: the metadata handle, as retruned by PGS_MET_Init.
    !   attributeName: the ECS name of the attribute.
    !   attributeValue: the value(s) to be assigned to the attribute. 
    !
    ! !Output Parameters:
    !   none.
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
    character(len = len(attributeValue)), dimension(1) :: localValue
    
    call getMetadataAttribute_s(handle, attributeName, localValue)
    attributeValue = localValue(1)
  end subroutine getMetadataAttribute_s_scalar
  ! --------------------------------------------------------------------
  subroutine getMetadataAttribute_i_scalar(handle, attributeName, attributeValue)
    character (len = *), intent ( in) :: handle, attributeName
    integer,             intent (out) :: attributeValue
    ! !F90
    !
    ! !Description:
    !    Get the value of a metadata attribute after the MCF has been initialized. 
    !    Error trapping included. 
    !  
    ! !Input Parameters:
    !   handle: the metadata handle, as retruned by PGS_MET_Init.
    !   attributeName: the ECS name of the attribute.
    !   attributeValue: the value(s) to be assigned to the attribute. 
    !
    ! !Output Parameters:
    !   none.
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
    integer, dimension(1) :: localValue
    
    call getMetadataAttribute_i(handle, attributeName, localValue)
    attributeValue = localValue(1)
  end subroutine getMetadataAttribute_i_scalar
  ! --------------------------------------------------------------------
  subroutine getMetadataAttribute_d_scalar(handle, attributeName, attributeValue)
    character (len = *), intent ( in) :: handle, attributeName
    double precision,    intent (out) :: attributeValue
    ! !F90
    !
    ! !Description:
    !    Get the value of a metadata attribute after the MCF has been initialized. 
    !    Error trapping included. 
    !  
    ! !Input Parameters:
    !   handle: the metadata handle, as retruned by PGS_MET_Init.
    !   attributeName: the ECS name of the attribute.
    !   attributeValue: the value(s) to be assigned to the attribute. 
    !
    ! !Output Parameters:
    !   none.
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
    double precision, dimension(1) :: localValue
    
    call getMetadataAttribute_d(handle, attributeName, localValue)
    attributeValue = localValue(1)
  end subroutine getMetadataAttribute_d_scalar
  ! --------------------------------------------------------------------
  
  ! --------------------------------------------------------------------
  subroutine writeMetadata(handle, hdfAttributeName, hdfFileID)
    character (len = *), intent( in) :: handle, hdfAttributeName
    integer,             intent( in) :: hdfFileID
    character (len = 2048)           :: msgbuffer

    ! !F90
    !
    ! !Description:
    !    Write metadata to HDF file attribute. Error trapping included.
    !  
    ! !Input Parameters:
    !   handle: metadata group in MCF
    !   hdfAttributeName: HDF attribute name to contain metadata.
    !   hdfFileID: HDF file ID as returned by the SD interface. 
    !
    ! !Output Parameters:
    !   None.
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
    !   Input arguments are identical to those for PGS_MET_Write.
    !
    ! !END

    ! Local variable
    integer :: status
    
    status = PGS_MET_Write(handle, hdfAttributeName, hdfFileID)
    if(status /= PGS_S_SUCCESS) then
      msgbuffer = "Error writing " // &
                   trim(hdfAttributeName) // &
                   "Return code is " // &
                   trim(IntToChar(status)) &
                   // char(10) // "Operator Action:  Notify SDST."
      call ErrorMessage("writeMetadata", msgbuffer)
    end if  
  
  end subroutine writeMetadata
end module Metadata
