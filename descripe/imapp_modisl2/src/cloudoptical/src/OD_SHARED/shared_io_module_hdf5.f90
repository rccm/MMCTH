 module shared_io_module_hdf5

  use hdf5
 
 implicit none
 
 include "hdf.f90"
 include "dffunc.f90"
 
  
 contains
 
 subroutine read_byte_array(filedata, sdsname, start, stride, edge, dataholder, status) 
 
	use nonscience_parameters
 
	integer, dimension(:), intent(in) :: filedata
	character(*), intent(in) :: sdsname
	integer, dimension(:), intent(in) :: start, stride, edge
	integer, intent(inout) :: status
	integer*1, dimension(:,:), intent(inout) :: dataholder
 
 
	integer :: file_id, var_id, err_code, attr_id
	integer*1 :: fill_val
	
	file_id = filedata(1)
	var_id = sfselect(file_id, sfn2index(file_id, sdsname))
			
	fill_val = 0
	
	status = sfrdata(var_id, start, stride, edge, dataholder)	
	err_code = sfendacc(var_id)
		
 end subroutine read_byte_array	 

 subroutine read_float_array_h4(filedata, sdsname, start, stride, edge, dataholder, status) 
 
	use nonscience_parameters
 
	integer, dimension(:), intent(in) :: filedata
	character(*), intent(in) :: sdsname
	integer, dimension(:), intent(in) :: start, stride, edge
	integer, intent(inout) :: status
	real, dimension(:,:), intent(inout) :: dataholder
 
 
	integer :: file_id, var_id, err_code, attr_id
	integer*1 :: fill_val
	
	file_id = filedata(1)
	var_id = sfselect(file_id, sfn2index(file_id, sdsname))
			
	fill_val = 0
	
	status = sfrdata(var_id, start, stride, edge, dataholder)	
	err_code = sfendacc(var_id)
		
 end subroutine read_float_array_h4	 


 
   subroutine read_byte_scaled_array(filedata, sdsname, start, stride, edge, use_offset, dataholder, status) 
 
	use nonscience_parameters
	use names
 
	integer, dimension(:), intent(in) :: filedata
	logical, intent(in) :: use_offset
	character(*), intent(in) :: sdsname
	integer, dimension(:), intent(in) :: start, stride, edge
	integer, intent(inout) :: status
	real, dimension(:,:), intent(inout) :: dataholder
 
 
	integer :: file_id, var_id, err_code, attr_id
	real :: scale_factor, offset
	integer*1 :: fill_val
	integer*1, dimension(:,:), allocatable :: temp
	integer:: i, j
	
	allocate(temp(edge(1), edge(2)))
	
	file_id = filedata(1)
	var_id = sfselect(file_id, sfn2index(file_id, trim(sdsname)))
	
	attr_id = sffattr(var_id, "SCALE_FACTOR")
	if (attr_id == -1) attr_id = sffattr(var_id, "scale_factor")
	err_code = sfrattr(var_id, attr_id, scale_factor)

	if (use_offset) then 
		attr_id = sffattr(var_id, "ADD_OFFSET")
		if (attr_id == -1) attr_id = sffattr(var_id, "add_offset")
		err_code = sfrattr(var_id, attr_id, offset)
	else
		offset = 0.
	endif
		
	attr_id = sffattr(var_id, "_FILLVALUE")
	if (attr_id == -1) attr_id = sffattr(var_id, "_FillValue")
	err_code = sfrattr(var_id, attr_id, fill_val)
	
	status = sfrdata(var_id, start, stride, edge, temp)

	err_code = sfendacc(var_id)

	do i=1, edge(1)
	  do j=1, edge(2)

		if (temp(i,j) /= fill_val) then 
			dataholder(i,j) = temp(i,j)*scale_factor + offset
		else
			dataholder(i,j) = fillvalue_real
  		endif

	 end do
       end do
	
	deallocate(temp)
	
 end subroutine read_byte_scaled_array	 

   subroutine read_int_scaled_array(filedata, sdsname, start, stride, edge, use_offset, dataholder, status) 
 
	use nonscience_parameters
	use names
 
	integer, dimension(:), intent(in) :: filedata
	logical, intent(in) :: use_offset
	character(*), intent(in) :: sdsname
	integer, dimension(:), intent(in) :: start, stride, edge
	integer, intent(inout) :: status
	real, dimension(:,:), intent(inout) :: dataholder
 
 
	integer :: file_id, var_id, err_code, attr_id
	real :: scale_factor, offset
	integer*2 :: fill_val
	integer*2, dimension(:,:), allocatable :: temp
	integer:: i, j
	
	allocate(temp(edge(1), edge(2)))
	
	file_id = filedata(1)
	var_id = sfselect(file_id, sfn2index(file_id, trim(sdsname)))
	
	attr_id = sffattr(var_id, "SCALE_FACTOR")
	if (attr_id == -1) attr_id = sffattr(var_id, "scale_factor")
	err_code = sfrattr(var_id, attr_id, scale_factor)

	if (use_offset) then 
		attr_id = sffattr(var_id, "ADD_OFFSET")
		if (attr_id == -1) attr_id = sffattr(var_id, "add_offset")
		err_code = sfrattr(var_id, attr_id, offset)
	else
		offset = 0.
	endif
		
	attr_id = sffattr(var_id, "_FILLVALUE")
	if (attr_id == -1) attr_id = sffattr(var_id, "_FillValue")
	err_code = sfrattr(var_id, attr_id, fill_val)
	
	status = sfrdata(var_id, start, stride, edge, temp)

	err_code = sfendacc(var_id)

	do i=1, edge(1)
	  do j=1, edge(2)

		if (temp(i,j) /= fill_val) then 
			dataholder(i,j) = temp(i,j)*scale_factor + offset
		else
			dataholder(i,j) = fillvalue_real
  		endif

	 end do
       end do
	
	deallocate(temp)
	
 end subroutine read_int_scaled_array	 


 
  subroutine read_int_array(filedata, sdsname, start, stride, edge, use_offset, dataholder, status) 
 
	use nonscience_parameters
 
	integer, dimension(:), intent(in) :: filedata
	logical, intent(in) :: use_offset
	character(*), intent(in) :: sdsname
	integer, dimension(:), intent(in) :: start, stride, edge
	integer, intent(inout) :: status
	real, dimension(:,:), intent(inout) :: dataholder
 
 
	integer :: file_id, var_id, err_code, attr_id
	real*8 :: scale_factor, offset
	integer*2 :: fill_val
	integer*2, dimension(:,:), allocatable :: temp
	integer:: i, j
	
	allocate(temp(edge(1), edge(2)))
	
	file_id = filedata(1)
	var_id = sfselect(file_id, sfn2index(file_id, sdsname))
	
	attr_id = sffattr(var_id, "scale_factor")
	err_code = sfrattr(var_id, attr_id, scale_factor)

	if (use_offset) then 
		attr_id = sffattr(var_id, "add_offset")
		err_code = sfrattr(var_id, attr_id, offset)
	else
		offset = 0.
	endif
		
	attr_id = sffattr(var_id, "_FillValue")
	err_code = sfrattr(var_id, attr_id, fill_val)
	
	err_code = sfrdata(var_id, start, stride, edge, temp)

	
	err_code = sfendacc(var_id)

	do i=1, edge(1)
	  do j=1, edge(2)

		if (temp(i,j) /= fill_val) then 
			dataholder(i,j) = (temp(i,j) - offset)*scale_factor
		else
			dataholder(i,j) = fillvalue_real
  		endif

	 end do
       end do
	
	deallocate(temp)
	
	status = 0

 end subroutine read_int_array	 

 subroutine read_float_array(filedata, sdsname, start, stride, edge, dataholder, status) 
 
	integer, dimension(:), intent(in) :: filedata
	character(*), intent(in) :: sdsname
	integer, dimension(:), intent(in) :: start, stride, edge
	integer, intent(inout) :: status
	real, dimension(:,:), intent(inout) :: dataholder
 
 
	integer :: file_id
	
	integer :: dataspace_id, memspace_id
	integer(HSIZE_T) :: edge1d(1)
	integer(HID_T) :: var_id, err_code

	real :: temp_fac(2)
	integer(HSIZE_T) :: offset(2), count(2), localstride(2), offset_mem(2)
	
	character(len=2) :: band_tag

	
	file_id = filedata(1)


	call h5dopen_f(file_id, trim(sdsname), var_id, err_code)
	call h5dget_space_f(var_id, dataspace_id, err_code)

	offset = start
	offset_mem = (/0,0/)
	count = edge
	localstride = stride

	call h5sselect_hyperslab_f(dataspace_id, H5S_SELECT_SET_F, offset, count, err_code, localstride)
	call h5screate_simple_f(2, count, memspace_id, err_code)
	call h5sselect_hyperslab_f(memspace_id, H5S_SELECT_SET_F, offset_mem, count, err_code, localstride)

	call h5dread_f(var_id, H5T_NATIVE_REAL, dataholder, count, err_code, memspace_id, dataspace_id)
		
	call h5sclose_f(dataspace_id, err_code)
	call h5sclose_f(memspace_id, err_code)
		
	call h5dclose_f(var_id, err_code)
		
	status = 0

 end subroutine read_float_array	 


 subroutine set_processing_attributes(filedata)
	
		use mod06_run_settings
		use specific_other
		implicit none
	
	integer, dimension(:), intent(in) :: filedata
	integer :: file_id, var_id, err_code, attr_id
	character :: yes_str(3), no_str(2)
	
	file_id = filedata(1)

	yes_str = "y"
	no_str = "n"

	if (DO_CSR) then 
		err_code = sfsattr(file_id, "Clear_Sky_Restoral_Status", DFNT_CHAR, 1, yes_str)
	else
		err_code = sfsattr(file_id, "Clear_Sky_Restoral_Status", DFNT_CHAR, 1, no_str)
	endif	

	if (DO_C4PhaseTest) then 
		err_code = sfsattr(file_id, "Collection_4_Phase_Used", DFNT_CHAR, 1, yes_str)
	else
		err_code = sfsattr(file_id, "Collection_4_Phase_Used", DFNT_CHAR, 1, no_str)
	endif	
	
	if (FORCE_ICE) then 
		err_code = sfsattr(file_id, "Ice_Phase_Forced", DFNT_CHAR, 1, yes_str)
	else
		err_code = sfsattr(file_id, "Ice_Phase_Forced", DFNT_CHAR, 1, no_str)
	endif	

	if (FORCE_WATER) then 
		err_code = sfsattr(file_id, "Water_Phase_Forced", DFNT_CHAR, 1, yes_str)
	else
		err_code = sfsattr(file_id, "Water_Phase_Forced", DFNT_CHAR, 1, no_str)
	endif	
	
	call set_process_time(file_id)
	
	call set_processing_extra(file_id)
	
 
 end subroutine set_processing_attributes

 
 
  subroutine openclose_files (l1bname,             &
                            cloudmask_name,      &
                            geolocation_name,    &
                            mod06_name,          &
                            directive,           &
                            l1b_filedata,        &
                            cloudmask_filedata,  &
                            geolocation_filedata,&
                            mod06_filedata,      &
                            status)

  use nonscience_parameters
  use mod06_run_settings
  use hdf5
  
  implicit none

	include "hdf.f90"
	include "dffunc.f90"

  character(*), intent(in) :: l1bname(:), cloudmask_name, geolocation_name, directive,mod06_name
  integer, intent(inout)    :: l1b_filedata(:), cloudmask_filedata(:),  &
                               geolocation_filedata(:),mod06_filedata(:)
  integer, intent(out)      :: status

  integer :: err_code, i

  status = success 

  if (directive == 'open') then

	if (MODIS_MODE .or. IFF_yes) then 

		l1b_filedata(1) = sfstart(l1bname(1), DFACC_READ)
	    if(l1b_filedata(1) == -1) then
	       CALL MODIS_SMF_SETDYNAMICMSG(2,'MAPI open failed for L1B file, Check PCF.','openclose_files')
	       status = failure
	    else
	       call MODIS_SMF_SETDYNAMICMSG(0,'MAPI open for L1B file OK','openclose_files')
	    endif    
		l1b_filedata(2) = -1
	
	else
	
		call h5open_f(err_code)

		do i=1, NUM_L1B
			call h5fopen_f(l1bname(i), H5F_ACC_RDONLY_F, l1b_filedata(i), err_code)
			if(l1b_filedata(1) == -1) then
				CALL MODIS_SMF_SETDYNAMICMSG(2,'MAPI open failed for L1B file, Check PCF.','openclose_files')
				status = failure
			else    
				call MODIS_SMF_SETDYNAMICMSG(0,'MAPI open for L1B file OK','openclose_files')
			endif
		end do
	endif

	cloudmask_filedata(1) = sfstart(cloudmask_name, DFACC_READ)
    if(cloudmask_filedata(1) == -1) then
       CALL MODIS_SMF_SETDYNAMICMSG(2,'MAPI open failed for cloudmask file, Check PCF.','openclose_files')
       status = failure
    else    
       call MODIS_SMF_SETDYNAMICMSG(0,'MAPI open for cloudmask file OK','openclose_files')
    endif

	if (IFF_yes) then 
		geolocation_filedata(1) = sfstart(geolocation_name, DFACC_READ)
	    if(geolocation_filedata(1) == -1) then
	       CALL MODIS_SMF_SETDYNAMICMSG(2,'MAPI open failed for geolocation file, Check PCF.','openclose_files')
	       status = failure
	    else    
	       call MODIS_SMF_SETDYNAMICMSG(0,'MAPI open for geolocation file OK','openclose_files')
	    endif
	
	else
		geolocation_filedata(1) = cloudmask_filedata(1)
	endif
	
	mod06_filedata(1) = sfstart(mod06_name, DFACC_WRITE)
    if(mod06_filedata(1) == -1) then
       CALL MODIS_SMF_SETDYNAMICMSG(2,'MAPI open failed for MOD06_L2, Check PCF.','openclose_files')
       status = failure
    else
       call MODIS_SMF_SETDYNAMICMSG(0,'MAPI open for MOD06_L2 file OK','openclose_files')
    endif
  else

	if (MODIS_MODE .or. IFF_yes) then 
		err_code = sfend(l1b_filedata(1))
	else
		do i=1, NUM_L1B
			call h5fclose_f(l1b_filedata(i), err_code)
		end do
		
		call h5close_f(err_code)
		
	endif
	
	err_code = sfend(cloudmask_filedata(1))
    if (err_code /= 0) then
      CALL MODIS_SMF_SETDYNAMICMSG(2,'MAPI close failed for cloudmask file, Check PCF.','openclose_files')
    else
      call MODIS_SMF_SETDYNAMICMSG(0,'MAPI close for cloudmask file OK','openclose_files')
    endif

	if (IFF_yes) then 

		err_code = sfend(geolocation_filedata(1))
	    if (err_code /= 0) then
	      CALL MODIS_SMF_SETDYNAMICMSG(2,'MAPI close failed for geolocation file, Check PCF.','openclose_files')
	    else
	      call MODIS_SMF_SETDYNAMICMSG(0,'MAPI close for geolocation file OK','openclose_files')
	    endif


	endif


	err_code = sfend(mod06_filedata(1))
    if (err_code /= 0) then
      CALL MODIS_SMF_SETDYNAMICMSG(2,'MAPI close failed for MOD06_L2, Check PCF.','openclose_files')
    else
      call MODIS_SMF_SETDYNAMICMSG(0,'MAPI close for MOD06_L2 file OK','openclose_files')
    endif 
	
	
  endif



 end subroutine openclose_files

 subroutine set_quality_data(xsize, ysize)

   use core_arrays
   use libraryarrays, only: water_radii, ice_radii, number_waterradii, number_iceradii

   integer , intent(in) :: xsize, ysize

   integer :: i, j

   integer, parameter :: marginal = 1, good = 2, very_good = 3, no_confidence = 0
   
   integer*1, parameter :: no_cloud_mask = 0, no_cloud = 1, liquid = 2, &
   							ice = 3, unknown = 4, good_retrieval = 8

   do i = 1, xsize
     do j = 1, ysize

       if (cloudsummary(i,j)%cloudmask_determined) then

            processing_information(i,j)%path_and_outcome = no_cloud
            processing_information(i,j)%path_and_outcome_16 = no_cloud
            processing_information(i,j)%path_and_outcome_37 = no_cloud
            processing_information(i,j)%path_and_outcome_PCL = no_cloud
            processing_information(i,j)%path_and_outcome_16_PCL = no_cloud
            processing_information(i,j)%path_and_outcome_37_PCL = no_cloud
       
          if(cloudsummary(i,j)%watercloud) then
            processing_information(i,j)%path_and_outcome = liquid
            processing_information(i,j)%path_and_outcome_16 = liquid
            processing_information(i,j)%path_and_outcome_37 = liquid
            processing_information(i,j)%path_and_outcome_PCL = liquid
            processing_information(i,j)%path_and_outcome_16_PCL = liquid
            processing_information(i,j)%path_and_outcome_37_PCL = liquid
          elseif(cloudsummary(i,j)%icecloud) then
            processing_information(i,j)%path_and_outcome = ice
            processing_information(i,j)%path_and_outcome_16 = ice
            processing_information(i,j)%path_and_outcome_37 = ice
            processing_information(i,j)%path_and_outcome_PCL = ice
            processing_information(i,j)%path_and_outcome_16_PCL = ice
            processing_information(i,j)%path_and_outcome_37_PCL = ice
          elseif(cloudsummary(i,j)%unknowncloud) then
            processing_information(i,j)%path_and_outcome = unknown
            processing_information(i,j)%path_and_outcome_16 = unknown
            processing_information(i,j)%path_and_outcome_37 = unknown
            processing_information(i,j)%path_and_outcome_PCL = unknown
            processing_information(i,j)%path_and_outcome_16_PCL = unknown
            processing_information(i,j)%path_and_outcome_37_PCL = unknown
          endif
       else 
          processing_information(i,j)%path_and_outcome = no_cloud_mask
          processing_information(i,j)%path_and_outcome_16 = no_cloud_mask
          processing_information(i,j)%path_and_outcome_37 = no_cloud_mask
          processing_information(i,j)%path_and_outcome_PCL = no_cloud_mask
          processing_information(i,j)%path_and_outcome_16_PCL = no_cloud_mask
          processing_information(i,j)%path_and_outcome_37_PCL = no_cloud_mask

       endif 

       if (cloudsummary(i,j)%cloudmask_determined .and. &
            ( cloudsummary(i,j)%ocean_surface .or.      &
              cloudsummary(i,j)%snowice_surface)      )  then

            processing_information(i,j)%path_and_outcome_1621 = no_cloud
            processing_information(i,j)%path_and_outcome_1621_PCL = no_cloud
			
          if(cloudsummary(i,j)%watercloud) then
            processing_information(i,j)%path_and_outcome_1621 = liquid
            processing_information(i,j)%path_and_outcome_1621_PCL = liquid
          elseif(cloudsummary(i,j)%icecloud) then
            processing_information(i,j)%path_and_outcome_1621 = ice
            processing_information(i,j)%path_and_outcome_1621_PCL = ice
          elseif(cloudsummary(i,j)%unknowncloud) then
            processing_information(i,j)%path_and_outcome_1621 = unknown
            processing_information(i,j)%path_and_outcome_1621_PCL = unknown
          endif 
		  
       else
          processing_information(i,j)%path_and_outcome_1621 = no_cloud_mask
          processing_information(i,j)%path_and_outcome_1621_PCL = no_cloud_mask
       endif

     
		!  confidence = 11, usefulness = 1
       	if (optical_thickness_final(i,j) > 0. .and. effective_radius_21_final(i,j) > 0.) then 
       		processing_information(i,j)%optical_thickness_GC = 7
       		processing_information(i,j)%effective_radius_GC = 7
       		processing_information(i,j)%water_path_GC = 7
       	
       	    processing_information(i,j)%path_and_outcome = &
         				processing_information(i,j)%path_and_outcome + good_retrieval  
	
       	else 
       		processing_information(i,j)%optical_thickness_GC = 0
       		processing_information(i,j)%effective_radius_GC = 0
       		processing_information(i,j)%water_path_GC = 0
		endif       		

		!  confidence = 11, usefulness = 1
       	if (optical_thickness_final_PCL(i,j) > 0. .and. effective_radius_21_final_PCL(i,j) > 0.) then 
       		processing_information(i,j)%optical_thickness_GC = 7
       		processing_information(i,j)%effective_radius_GC = 7
       		processing_information(i,j)%water_path_GC = 7

         	processing_information(i,j)%path_and_outcome_PCL = &
         				processing_information(i,j)%path_and_outcome_PCL + good_retrieval  
       	endif 
     		

		!  confidence = 11, usefulness = 1
       	if (optical_thickness_1621_final(i,j) > 0. .and. effective_radius_1621_final(i,j) > 0.)  then
      		processing_information(i,j)%optical_thickness_1621_GC = 7
       		processing_information(i,j)%effective_radius_1621_GC = 7
       		processing_information(i,j)%water_path_1621_GC = 7

	        processing_information(i,j)%path_and_outcome_1621 = &
     					processing_information(i,j)%path_and_outcome_1621 + good_retrieval  
		else 
      		processing_information(i,j)%optical_thickness_1621_GC = 0
       		processing_information(i,j)%effective_radius_1621_GC = 0
       		processing_information(i,j)%water_path_1621_GC = 0
		endif

		!  confidence = 11, usefulness = 1
       	if (optical_thickness_1621_final_PCL(i,j) > 0. .and. effective_radius_1621_final_PCL(i,j) > 0.)  then
      		processing_information(i,j)%optical_thickness_1621_GC = 7
       		processing_information(i,j)%effective_radius_1621_GC = 7
       		processing_information(i,j)%water_path_1621_GC = 7
         	
         	processing_information(i,j)%path_and_outcome_1621_PCL = &
         				processing_information(i,j)%path_and_outcome_1621_PCL + good_retrieval  
		endif

	   if (optical_thickness_16_final_PCL(i,j) > 0. .and. effective_radius_16_final_PCL(i,j) > 0.) & 
         processing_information(i,j)%path_and_outcome_16_PCL = &
         					processing_information(i,j)%path_and_outcome_16_PCL + good_retrieval  
	   if (optical_thickness_37_final_PCL(i,j) > 0. .and. effective_radius_37_final_PCL(i,j) > 0.) & 
         processing_information(i,j)%path_and_outcome_37_PCL = &
         					processing_information(i,j)%path_and_outcome_37_PCL + good_retrieval  
         					
	   if (optical_thickness_16_final(i,j) > 0. .and. effective_radius_16_final(i,j) > 0.) & 
         processing_information(i,j)%path_and_outcome_16 = &
         					processing_information(i,j)%path_and_outcome_16 + good_retrieval  
	   if (optical_thickness_37_final(i,j) > 0. .and. effective_radius_37_final(i,j) > 0.) & 
         processing_information(i,j)%path_and_outcome_37 = &
         					processing_information(i,j)%path_and_outcome_37 + good_retrieval  
		   
		   
        
       if (processing_information(i,j)%optical_thickness_GC /= 0 .and.  &
           processing_information(i,j)%band_used_for_optical_thickness ==1 )then
         processing_information(i,j)%rayleigh_correction = 1 
       else
         processing_information(i,j)%rayleigh_correction = 0
       endif 


! initialize the variable first of all. -- GW 4.6.05
       processing_information(i,j)%multi_layer_cloud = 0

       if (cloudsummary(i,j)%cloudmask_determined) then

         if (processing_information(i,j)%path_and_outcome == 1)then ! dec.tree stop
            processing_information(i,j)%multi_layer_cloud = 1

         elseif (processing_information(i,j)%path_and_outcome == 2 .or. &
         		 processing_information(i,j)%path_and_outcome == 10 ) then ! water cloud
            if (cloud_layer_flag(i,j) < 2 .or. ml_test_flag(i,j) == 16) then 
               processing_information(i,j)%multi_layer_cloud = 2 ! SL water
            else 
               processing_information(i,j)%multi_layer_cloud = 3 ! ML water
            endif

         elseif(processing_information(i,j)%path_and_outcome == 3 .or. &
         		processing_information(i,j)%path_and_outcome == 11 ) then !ice cloud
            if (cloud_layer_flag(i,j) < 2 .or. ml_test_flag(i,j) == 16) then 
               processing_information(i,j)%multi_layer_cloud = 4 ! SL ice
            else 
               processing_information(i,j)%multi_layer_cloud = 5 ! ML ice
            endif

         elseif (processing_information(i,j)%path_and_outcome == 4 .or. &
         		 processing_information(i,j)%path_and_outcome == 12 ) then ! unknown cloud
            if (cloud_layer_flag(i,j) < 2 .or. ml_test_flag(i,j) == 16) then 
               processing_information(i,j)%multi_layer_cloud = 6 ! SL unknown
            else 
               processing_information(i,j)%multi_layer_cloud = 7 ! ML unknown
            endif
         endif

       else
          processing_information(i,j)%multi_layer_cloud = 0 
       endif

! set the information for CSR QA -- GW. 4.7.05
       processing_information(i,j)%CSR_flag = 0
       processing_information(i,j)%CSR_flag = CSR_flag_array(i,j)

! set the information for the ML Test QA -- GW. 5.13.09
	   processing_information(i,j)%ml_test_mark = 0
	   processing_information(i,j)%ml_test_mark = ml_test_flag(i,j)

     enddo
   enddo
 end subroutine set_quality_data

 subroutine convert_binary_qa( quality_assurance_1km, &
                             status)
 use core_arrays, only: processing_information
 use nonscience_parameters

 implicit none
 
 integer*1 , intent(out) ::  quality_assurance_1km(:,:,:)
 integer, intent(inout) :: status

 integer   :: i,j, blah, cm_wid, cm_ht
 
 quality_assurance_1km = 0

 cm_wid = size(processing_information, 1)
 cm_ht = size(processing_information, 2)


  do j= 1, cm_ht  
	do i = 1, cm_wid


!    Quality Assurance 1KM, byte 1 ----------------------------------------------------------------------------
		quality_assurance_1km(1,i,j) = processing_information(i,j)%optical_thickness_GC

     if(processing_information(i,j)%optical_thickness_outofbounds== 1) then
           quality_assurance_1km(1,i,j) =  ibset(quality_assurance_1km(1,i,j),3)
     elseif(processing_information(i,j)%optical_thickness_outofbounds== 2) then
           quality_assurance_1km(1,i,j) =  ibset(quality_assurance_1km(1,i,j),4)
     elseif(processing_information(i,j)%optical_thickness_outofbounds == 3) then 
           quality_assurance_1km(1,i,j) =  ibset(quality_assurance_1km(1,i,j),3)
           quality_assurance_1km(1,i,j) =  ibset(quality_assurance_1km(1,i,j),4)
     endif

	 if (processing_information(i,j)%effective_radius_GC /= 0) then 
	 	quality_assurance_1km(1,i,j) =  ibset(quality_assurance_1km(1,i,j),5)
	 	quality_assurance_1km(1,i,j) =  ibset(quality_assurance_1km(1,i,j),6)
	 	quality_assurance_1km(1,i,j) =  ibset(quality_assurance_1km(1,i,j),7)
	 endif



!    Quality Assurance 1KM, byte 2 ----------------------------------------------------------------------------

	 quality_assurance_1km(2,i,j) = ishft(processing_information(i,j)%path_and_outcome_1621, 3)

	 if (processing_information(i,j)%water_path_GC /= 0) then 
           quality_assurance_1km(2,i,j) =  ibset(quality_assurance_1km(2,i,j),0)
           quality_assurance_1km(2,i,j) =  ibset(quality_assurance_1km(2,i,j),1)
           quality_assurance_1km(2,i,j) =  ibset(quality_assurance_1km(2,i,j),2)
	 endif	 	

!    Quality Assurance 1KM, byte 3 ----------------------------------------------------------------------------

	  quality_assurance_1km(3,i,j) = processing_information(i,j)%path_and_outcome

     if(processing_information(i,j)%rayleigh_correction == 1) &
          quality_assurance_1km(3,i,j) = ibset(quality_assurance_1km(3,i,j),4) 

! atmospheric correction is always done. 
     quality_assurance_1km(3,i,j) = ibset(quality_assurance_1km(3,i,j),5) 

     if(processing_information(i,j)%band_used_for_optical_thickness ==1 ) then
          quality_assurance_1km(3,i,j) = ibset(quality_assurance_1km(3,i,j),6) 
     elseif(processing_information(i,j)%band_used_for_optical_thickness == 2 ) then
          quality_assurance_1km(3,i,j) = ibset(quality_assurance_1km(3,i,j),7)
     elseif(processing_information(i,j)%band_used_for_optical_thickness == 3 ) then
          quality_assurance_1km(3,i,j) = ibset(quality_assurance_1km(3,i,j),6)
          quality_assurance_1km(3,i,j) = ibset(quality_assurance_1km(3,i,j),7) 
     endif

!    Quality Assurance 1KM, byte 4 ----------------------------------------------------------------------------

	 quality_assurance_1km(4,i,j) = processing_information(i,j)%optical_thickness_1621_GC

	 if (processing_information(i,j)%effective_radius_1621_GC /= 0) then 
          quality_assurance_1km(4,i,j) = ibset(quality_assurance_1km(4,i,j),3)
          quality_assurance_1km(4,i,j) = ibset(quality_assurance_1km(4,i,j),4)
          quality_assurance_1km(4,i,j) = ibset(quality_assurance_1km(4,i,j),5)
	 endif
	 
     ! CSR QA added by G.Wind 4.7.05
     if (processing_information(i,j)%CSR_flag == 1) &
          quality_assurance_1km(4,i,j) = ibset(quality_assurance_1km(4,i,j),6)
     if (processing_information(i,j)%CSR_flag == 2) &
          quality_assurance_1km(4,i,j) = ibset(quality_assurance_1km(4,i,j),7)
     if (processing_information(i,j)%CSR_flag == 3) then
          quality_assurance_1km(4,i,j) = ibset(quality_assurance_1km(4,i,j),6)
          quality_assurance_1km(4,i,j) = ibset(quality_assurance_1km(4,i,j),7)
     endif
	

!    Quality Assurance 1KM, byte 5 ----------------------------------------------------------------------------

	 quality_assurance_1km(5,i,j) = processing_information(i,j)%water_path_1621_GC

     if (processing_information(i,j)%multi_layer_cloud == 1) then
       quality_assurance_1km(5,i,j) = ibset(quality_assurance_1km(5,i,j),3)
     elseif(processing_information(i,j)%multi_layer_cloud == 2) then
       quality_assurance_1km(5,i,j) = ibset(quality_assurance_1km(5,i,j),4)
     elseif(processing_information(i,j)%multi_layer_cloud == 3) then
       quality_assurance_1km(5,i,j) = ibset(quality_assurance_1km(5,i,j),3)
       quality_assurance_1km(5,i,j) = ibset(quality_assurance_1km(5,i,j),4)
     elseif(processing_information(i,j)%multi_layer_cloud == 4) then
       quality_assurance_1km(5,i,j) = ibset(quality_assurance_1km(5,i,j),5)
     elseif(processing_information(i,j)%multi_layer_cloud == 5) then
       quality_assurance_1km(5,i,j) = ibset(quality_assurance_1km(5,i,j),3)
       quality_assurance_1km(5,i,j) = ibset(quality_assurance_1km(5,i,j),5)
     elseif(processing_information(i,j)%multi_layer_cloud == 6) then
       quality_assurance_1km(5,i,j) = ibset(quality_assurance_1km(5,i,j),4)
       quality_assurance_1km(5,i,j) = ibset(quality_assurance_1km(5,i,j),5)
     elseif(processing_information(i,j)%multi_layer_cloud == 7) then
       quality_assurance_1km(5,i,j) = ibset(quality_assurance_1km(5,i,j),3)
       quality_assurance_1km(5,i,j) = ibset(quality_assurance_1km(5,i,j),4)
       quality_assurance_1km(5,i,j) = ibset(quality_assurance_1km(5,i,j),5)
     endif 

     !This is a COPY of the retrieval_outcome bit for Level 3 compatibility
     if(processing_information(i,j)%path_and_outcome > 4)  &
          quality_assurance_1km(5,i,j) = ibset(quality_assurance_1km(5,i,j),6)

!    Quality Assurance 1KM, byte 6 ----------------------------------------------------------------------------

		if (processing_information(i,j)%path_and_outcome_16 > 8) &
			quality_assurance_1km(6,i,j) = ibset(quality_assurance_1km(6,i,j), 0)			
		if (processing_information(i,j)%path_and_outcome_16_PCL > 8) &
			quality_assurance_1km(6,i,j) = ibset(quality_assurance_1km(6,i,j), 1)			

		if (processing_information(i,j)%path_and_outcome_37 > 8) &
			quality_assurance_1km(6,i,j) = ibset(quality_assurance_1km(6,i,j), 2)			
		if (processing_information(i,j)%path_and_outcome_37_PCL > 8) &
			quality_assurance_1km(6,i,j) = ibset(quality_assurance_1km(6,i,j), 3)			

		if (processing_information(i,j)%path_and_outcome_1621 > 8) &
			quality_assurance_1km(6,i,j) = ibset(quality_assurance_1km(6,i,j), 4)			
		if (processing_information(i,j)%path_and_outcome_PCL > 8) &
			quality_assurance_1km(6,i,j) = ibset(quality_assurance_1km(6,i,j), 5)			


  enddo
 enddo   
  
 end subroutine convert_binary_qa
	

 
 end module shared_io_module_hdf5
