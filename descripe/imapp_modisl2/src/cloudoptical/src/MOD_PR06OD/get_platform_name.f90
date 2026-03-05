      function get_platform_name(filename, NAME)
      
!-----------------------------------------------------------------------
! !F90
!
! !DESCRIPTION:
!     Given the filename of a MODIS MOD021KM or MYD021KM product file,
!     return the name of the satellite platform (Terra or Aqua).
!    This is a rewrite of a toolkit UWisc Madison function that does the same thing
!    only this one is without toolkits
!
! !INPUT PARAMETERS:
!     filename     name of a MOD021KM or MYD021KM file
!           
! !OUTPUT PARAMETERS:
!     GET_PLATFORM_NAME  Return status
!                        (0=success, -1=failure)
!     NAME               Platform name
!                        ('Terra' or 'Aqua' if successful)
!                        (' ' otherwise)
!
! NOTES:
!     (1) For a description of STATUS values, see the SDP toolkit
!         documentation for 'PGS_MET_GetPCAttr'. The status values
!         are defined in the include file PGS_MET_13.f.
!
! !REVISION HISTORY:
!
! !TEAM-UNIQUE HEADER:
!		Written by Gala Wind
!	    Cloud Retrieval Group
!		SSAI Inc, / NASA GSFC
!		Greenbelt MD 20771
!		wind@climate.gsfc.nasa.gov
!
! !END
!-----------------------------------------------------------------------
	implicit none

	include "hdf.f90"
	include "dffunc.f90"
	

 integer :: get_platform_name
 integer :: file_id, err_code, attr_id
 integer :: attr_size, attr_type
 
 character(len=*), intent(in) :: filename
 character(len=*), intent(inout) :: NAME
 character(len=20000) :: metadata
 character(len=200) :: attr_name, temp 
 
 integer :: idx, ln, idx2
  
 file_id = sfstart(filename, DFACC_READ)
 
 attr_id = sffattr(file_id, "CoreMetadata.0")
 
 err_code = sfgainfo(file_id, attr_id, attr_name, attr_type, attr_size)
 
 
 err_code = sfrattr(file_id, attr_id, metadata)
  
 err_code = sfend(file_id)
 
 idx = index(metadata, "ASSOCIATEDPLATFORMSHORTNAME")
 ln = len(metadata)

 temp = metadata(idx:idx+190)
 
 idx2 = index(temp, "Terra")
 
 if (idx2 > 0) then
	NAME =  "Terra"
	get_platform_name = 0
	return
 endif
 
 if (idx2 == 0) then 
	idx2 = index(temp, "Aqua")
	if (idx2 == 0) then 
		NAME = " "
		get_platform_name = -1
		return 
	endif
	NAME = "Aqua"
	get_platform_name = 0
	return
 endif


      END function get_platform_name
