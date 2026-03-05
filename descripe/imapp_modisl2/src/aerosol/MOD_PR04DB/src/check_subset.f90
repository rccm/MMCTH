integer function check_subset (geo_file, l1b_file, is_subset, offsets) result(status)

   use hdf

   implicit none

!   include 'PGS_MODIS_39500.f'
    
   character(len=*), intent(in)                 ::  geo_file
   character(len=*), intent(in)                 ::  l1b_file
   logical, intent(inout)                       ::  is_subset
   integer, dimension(2), intent(inout)         ::  offsets
   
   integer, dimension(2)                        ::  g_offsets
   integer, dimension(2)                        ::  l_offsets
   
   integer                                      ::  attr_id
   character(len=256)                           ::  attr_name
	 integer                                      ::  hdfid
   
   offsets = (/-999,-999/)
   is_subset = .false.
   
   print *, 'geolocation name: ', geo_file
   hdfid = SFstart(geo_file,DFACC_READ)
   
   attr_name = "originalPixelOffsets"
   attr_id = sffattr(hdfid, attr_name)
   if (attr_id == FAIL) then
    print *, "WARNING: Failed to get id of attribute for geolocation file."// &
             " Assuming full granule." 
    status = 0
    return
   end if
   
   status = SFrnatt(hdfid, attr_id, g_offsets)
   if (status == FAIL) then
    print *, "ERROR: Failed to read offset attribute in geolocation file: ", status
    return
   end if
   status = sfend(hdfid)
   if (status == FAIL) then
    print *, "ERROR: Failed to close geolocation file: ", status
    return
   end if
   
   print *, 'l1b name: ', l1b_file
   hdfid = SFstart(geo_file,DFACC_READ)
   attr_name = "originalPixelOffsets"
   attr_id = sffattr(hdfid, attr_name)
   if (attr_id == FAIL) then
    print *, "ERROR: Geolocation has originalPixelOffset attribute, but L1B does not: ", status
    status = attr_id
    return
   endif
   status = SFrnatt(hdfid, attr_id, l_offsets)
   if (status == FAIL) then
    print *, "ERROR: Failed to read offset attribute in l1b file: ", status
    return
   endif
   status = sfend(hdfid)
   if (status == FAIL) then
    print *, "ERROR: Failed to close l1b file: ", status
    return
   end if
   
! check that the offsets are the same
   if (g_offsets(1) /= l_offsets(1) .OR. g_offsets(2) /= l_offsets(2)) then
    print *, "ERROR: Got different offsets from geolocation and l1b files."
    status = -1
    return
   end if
   
   is_subset = .true.
   offsets = g_offsets
   status = 0
   return
   
end function check_subset

