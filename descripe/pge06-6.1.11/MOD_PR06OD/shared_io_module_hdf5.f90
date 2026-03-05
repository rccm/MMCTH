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
	
	err_code = sfrdata(var_id, start, stride, edge, dataholder)
	
	err_code = sfendacc(var_id)
		
	
	status = 0

 end subroutine read_byte_array	 

 
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
	err_code = sfrattr(var_id, attr_id, scale_factor)

	if (use_offset) then 
		attr_id = sffattr(var_id, "ADD_OFFSET")
		err_code = sfrattr(var_id, attr_id, offset)
	else
		offset = 0.
	endif
		
	attr_id = sffattr(var_id, "_FILLVALUE")
	err_code = sfrattr(var_id, attr_id, fill_val)
	
	err_code = sfrdata(var_id, start, stride, edge, temp)

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
	
	status = 0

 end subroutine read_byte_scaled_array	 

 
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

	cloudmask_filedata(1) = sfstart(cloudmask_name, DFACC_READ)
    if(cloudmask_filedata(1) == -1) then
       CALL MODIS_SMF_SETDYNAMICMSG(2,'MAPI open failed for cloudmask file, Check PCF.','openclose_files')
       status = failure
    else    
       call MODIS_SMF_SETDYNAMICMSG(0,'MAPI open for cloudmask file OK','openclose_files')
    endif

	call h5fopen_f(geolocation_name, H5F_ACC_RDONLY_F, geolocation_filedata(1), err_code)	
	if(geolocation_filedata(1) == -1) then
       CALL MODIS_SMF_SETDYNAMICMSG(2,'MAPI open failed for geolocation, Check PCF.','openclose_files')
       status = failure 
    else    
       call MODIS_SMF_SETDYNAMICMSG(0,'MAPI open for geolocation file OK','openclose_files')
    endif
	
	mod06_filedata(1) = sfstart(mod06_name, DFACC_WRITE)
    if(mod06_filedata(1) == -1) then
       CALL MODIS_SMF_SETDYNAMICMSG(2,'MAPI open failed for MOD06_L2, Check PCF.','openclose_files')
       status = failure
    else
       call MODIS_SMF_SETDYNAMICMSG(0,'MAPI open for MOD06_L2 file OK','openclose_files')
    endif
  else

	do i=1, NUM_L1B
		
		call h5fclose_f(l1b_filedata(i), err_code)

	end do

	err_code = sfend(cloudmask_filedata(1))
    if (err_code /= 0) then
      CALL MODIS_SMF_SETDYNAMICMSG(2,'MAPI close failed for cloudmask file, Check PCF.','openclose_files')
    else
      call MODIS_SMF_SETDYNAMICMSG(0,'MAPI close for cloudmaskfile OK','openclose_files')
    endif


	call h5fclose_f(geolocation_filedata(1), err_code)
  
	err_code = sfend(mod06_filedata(1))
    if (err_code /= 0) then
      CALL MODIS_SMF_SETDYNAMICMSG(2,'MAPI close failed for MOD06_L2, Check PCF.','openclose_files')
    else
      call MODIS_SMF_SETDYNAMICMSG(0,'MAPI close for MOD06_L2 file OK','openclose_files')
    endif 
	
	
	call h5close_f(err_code)
	
  endif



 end subroutine openclose_files

 subroutine set_quality_data(xsize, ysize)

   use core_arrays
   use libraryarrays, only: water_radii, ice_radii, number_waterradii, number_iceradii

   integer , intent(in) :: xsize, ysize

   integer :: i, j

   integer, parameter :: marginal = 1, good = 2, very_good = 3, no_confidence = 0

   do i = 1, xsize
     do j = 1, ysize
       if (optical_thickness_final(i,j) > 0. .and. effective_radius_21_final(i,j) > 0. .and. CSR_flag_array(i,j) /= 2) then
         processing_information(i,j)%optical_thickness_general = 1
       else
         processing_information(i,j)%optical_thickness_general = 0
       endif
      
       if (effective_radius_21_final(i,j) > 4. .and. CSR_flag_array(i,j) /= 2 ) then
         processing_information(i,j)%effective_radius_general = 1
       else
         processing_information(i,j)%effective_radius_general = 0
       endif 


!    Determine QA confidence level for tau, re, and water path - several changes
!        made in way logic constructed as well as updates to QA matrix for collection 5
!         - GTA 2/17/05

!     Set default to no confidence answer
      processing_information(i,j)%optical_thickness_confidence = no_confidence
      processing_information(i,j)%effective_radius_confidence = no_confidence
      processing_information(i,j)%water_path_confidence = no_confidence

!     Process according to QA matrix (12/30/04) for water cloud

       if( (cloudsummary(i,j)%watercloud .or. cloudsummary(i,j)%unknowncloud) .and. CSR_flag_array(i,j) == 0) then
       
			if ( effective_radius_21_final(i,j) < water_radii(number_waterradii) .and. &
				 effective_radius_21_final(i,j) > water_radii(1) ) then
               processing_information(i,j)%optical_thickness_confidence = very_good
               processing_information(i,j)%effective_radius_confidence = very_good
               processing_information(i,j)%water_path_confidence = very_good
             endif
             
       endif

!     Process according to QA matrix (12/30/04) for ice cloud

       if(cloudsummary(i,j)%icecloud .and. CSR_flag_array(i,j) == 0) then
	   
          if (effective_radius_21_final(i,j) < ice_radii(number_iceradii) .and. &
              effective_radius_21_final(i,j) > ice_radii(1) ) then
               processing_information(i,j)%optical_thickness_confidence = very_good
               processing_information(i,j)%effective_radius_confidence = very_good
               processing_information(i,j)%water_path_confidence = very_good
          endif
       
	   endif


       if (optical_thickness_1621_final(i,j) > 0. .and. effective_radius_1621_final(i,j) > 0. .and. CSR_flag_array(i,j) /= 2)  then
         processing_information(i,j)%optical_thickness_1621_gen = 1
       else
         processing_information(i,j)%optical_thickness_1621_gen = 0
       endif
  
       if (effective_radius_1621_final(i,j) > 4.0 .and. CSR_flag_array(i,j) /= 2) then
         processing_information(i,j)%effective_radius_1621_gen = 1
       else
         processing_information(i,j)%effective_radius_1621_gen = 0
       endif


!     Process confidence levels for 16-21 retrieval

!     Set default to no confidence answer
      processing_information(i,j)%optical_thickness_1621_conf = no_confidence
      processing_information(i,j)%effective_radius_1621_conf = no_confidence
      processing_information(i,j)%water_path_1621_conf = no_confidence


      if( (cloudsummary(i,j)%watercloud .or. cloudsummary(i,j)%unknowncloud) .and. CSR_flag_array(i,j) == 0) then

			if ( effective_radius_1621_final(i,j) < water_radii(number_waterradii) .and. &
				 effective_radius_1621_final(i,j) > water_radii(1) ) then
               processing_information(i,j)%optical_thickness_1621_conf= very_good
               processing_information(i,j)%effective_radius_1621_conf= very_good
               processing_information(i,j)%water_path_1621_conf= very_good
             endif

      endif

       if(cloudsummary(i,j)%icecloud .and. CSR_flag_array(i,j) == 0) then
          if (effective_radius_1621_final(i,j) < ice_radii(number_iceradii) .and. &
              effective_radius_1621_final(i,j) > ice_radii(1) ) then
             processing_information(i,j)%optical_thickness_1621_conf= very_good
             processing_information(i,j)%effective_radius_1621_conf= very_good
             processing_information(i,j)%water_path_1621_conf= very_good
         endif
       endif


! Added by Riedi Mon Apr 11 20:21:52 CEST 2005 
!     Process according to cloudsuspicious flag for all clouds and all retrievals
! We don't want to set any QA whatsoever for failed retrievals. G.Wind: 11.17.2011
       if(cloudsummary(i,j)%cloudsuspicious .and. CSR_flag_array(i,j) == 0) then
	   
			if (optical_thickness_final(i,j) > 0. .and. effective_radius_21_final(i,j) > 0.) then	
				processing_information(i,j)%optical_thickness_confidence = marginal
				processing_information(i,j)%effective_radius_confidence = marginal
				processing_information(i,j)%water_path_confidence = marginal
			endif

			if (optical_thickness_1621_final(i,j) > 0. .and. effective_radius_1621_final(i,j) > 0.) then	
				processing_information(i,j)%optical_thickness_1621_conf= marginal
				processing_information(i,j)%effective_radius_1621_conf= marginal
				processing_information(i,j)%water_path_1621_conf= marginal
			endif
       endif 
! End of Riedi modifs


       if (liquid_water_path(i,j) > 0. .and. CSR_flag_array(i,j) /= 2) then
         processing_information(i,j)%water_path_general = 1
       else
         processing_information(i,j)%water_path_general = 0
       endif

       if (liquid_water_path_1621(i,j) > 0. .and. CSR_flag_array(i,j) /= 2) then
         processing_information(i,j)%water_path_1621_gen = 1
       else
          processing_information(i,j)%water_path_1621_gen = 0
       endif


       if (cloudsummary(i,j)%cloudmask_determined) then

            processing_information(i,j)%cloud_retrieval_processing_path = 1
       
          if(cloudsummary(i,j)%watercloud) then
            processing_information(i,j)%cloud_retrieval_processing_path = 2
          elseif(cloudsummary(i,j)%icecloud) then
            processing_information(i,j)%cloud_retrieval_processing_path = 3
          elseif(cloudsummary(i,j)%unknowncloud) then
            processing_information(i,j)%cloud_retrieval_processing_path = 4
          endif
       else 
          processing_information(i,j)%cloud_retrieval_processing_path = 0
       endif 

       if (processing_information(i,j)%optical_thickness_general /= 0 .and. &
           processing_information(i,j)%effective_radius_general /= 0 ) then
         processing_information(i,j)%retrieval_outcome = 1 
       else
         processing_information(i,j)%retrieval_outcome = 0
       endif
        
       processing_information(i,j)%atmospheric_correction = 1
       if (processing_information(i,j)%optical_thickness_general /= 0 .and.  &
           processing_information(i,j)%band_used_for_optical_thickness ==1 )then
         processing_information(i,j)%rayleigh_correction = 1 
       else
         processing_information(i,j)%rayleigh_correction = 0
       endif 

       if (cloudsummary(i,j)%cloudmask_determined .and. &
            ( cloudsummary(i,j)%ocean_surface .or.      &
              cloudsummary(i,j)%snowice_surface)      )  then

            processing_information(i,j)%cloud_retrieval_proc_path_1621= 1
			
          if(cloudsummary(i,j)%watercloud) then
            processing_information(i,j)%cloud_retrieval_proc_path_1621= 2
          elseif(cloudsummary(i,j)%icecloud) then
            processing_information(i,j)%cloud_retrieval_proc_path_1621= 3
          elseif(cloudsummary(i,j)%unknowncloud) then
            processing_information(i,j)%cloud_retrieval_proc_path_1621= 4
          endif 
		  
       else
         processing_information(i,j)%cloud_retrieval_proc_path_1621= 0    
       endif

       if (processing_information(i,j)%optical_thickness_1621_gen /= 0 .and. &
           processing_information(i,j)%effective_radius_1621_gen /= 0) then
           processing_information(i,j)%retrieval_outcome_1621 = 1
       else
           processing_information(i,j)%retrieval_outcome_1621 = 0 
       endif
 


! initialize the variable first of all. -- GW 4.6.05
       processing_information(i,j)%multi_layer_cloud = 0

       if (cloudsummary(i,j)%cloudmask_determined) then

         if (processing_information(i,j)%cloud_retrieval_processing_path == 1)then ! dec.tree stop
            processing_information(i,j)%multi_layer_cloud = 1

         elseif (processing_information(i,j)%cloud_retrieval_processing_path == 2) then ! water cloud
            if (cloud_layer_flag(i,j) < 2 .or. ml_test_flag(i,j) == 16) then 
               processing_information(i,j)%multi_layer_cloud = 2 ! SL water
            else 
               processing_information(i,j)%multi_layer_cloud = 3 ! ML water
            endif

         elseif(processing_information(i,j)%cloud_retrieval_processing_path == 3) then !ice cloud
            if (cloud_layer_flag(i,j) < 2 .or. ml_test_flag(i,j) == 16) then 
               processing_information(i,j)%multi_layer_cloud = 4 ! SL ice
            else 
               processing_information(i,j)%multi_layer_cloud = 5 ! ML ice
            endif

         elseif (processing_information(i,j)%cloud_retrieval_processing_path == 4) then ! unknown cloud
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


 subroutine convert_binary_qa(cloud_mask_1km,         &
                             quality_assurance_1km, &
                             status)
 use core_arrays, only: processing_information
 use nonscience_parameters

 implicit none
 
 integer*1 , intent(out) :: cloud_mask_1km(:,:,:), quality_assurance_1km(:,:,:)
 integer, intent(inout) :: status

 integer   :: i,j, blah, cm_wid, cm_ht
 
 cloud_mask_1km = 0
 quality_assurance_1km = 0

 cm_wid = size(cloud_mask_1km, 2)
 cm_ht = size(cloud_mask_1km, 3)

  do j= 1, cm_ht  
	do i = 1, cm_wid

!    Cloud Mask 1KM, Byte 1 -----------------------------------------------------------------------------------

     if(processing_information(i,j)%cloud_mask_usefulness == 1) cloud_mask_1km(1,i,j) = ibset(cloud_mask_1km(1,i,j),0)
    
     if (processing_information(i,j)%cloud_mask_unobs_fov ==1) then
         cloud_mask_1km(1,i,j) = ibset(cloud_mask_1km(1,i,j),1)
     elseif (processing_information(i,j)%cloud_mask_unobs_fov == 2) then
         cloud_mask_1km(1,i,j) = ibset(cloud_mask_1km(1,i,j),2)
     elseif  (processing_information(i,j)%cloud_mask_unobs_fov == 3) then
         cloud_mask_1km(1,i,j) = ibset(cloud_mask_1km(1,i,j),1)
         cloud_mask_1km(1,i,j) = ibset(cloud_mask_1km(1,i,j),2)
     endif

     if(processing_information(i,j)%day_night == 1) cloud_mask_1km(1,i,j) = ibset(cloud_mask_1km(1,i,j),3)
     if(processing_information(i,j)%sun_glint == 1) cloud_mask_1km(1,i,j) = ibset(cloud_mask_1km(1,i,j),4)
     if(processing_information(i,j)%snow_ice == 1)  cloud_mask_1km(1,i,j) = ibset(cloud_mask_1km(1,i,j),5)

     if(processing_information(i,j)%surface_type == 1) then
        cloud_mask_1km(1,i,j) = ibset(cloud_mask_1km(1,i,j),6)
     elseif(processing_information(i,j)%surface_type == 2) then
        cloud_mask_1km(1,i,j) = ibset(cloud_mask_1km(1,i,j),7)
     elseif(processing_information(i,j)%surface_type == 3) then
        cloud_mask_1km(1,i,j) = ibset(cloud_mask_1km(1,i,j),6)
        cloud_mask_1km(1,i,j) = ibset(cloud_mask_1km(1,i,j),7)
     endif

!    Cloud Mask 1KM, Byte 2 -----------------------------------------------------------------------------------

     if(processing_information(i,j)%heavy_aerosol == 1) cloud_mask_1km(2,i,j) = ibset(cloud_mask_1km(2,i,j),0)
     if(processing_information(i,j)%thin_cirrus == 1)  cloud_mask_1km(2,i,j) = ibset(cloud_mask_1km(2,i,j),1)
     if(processing_information(i,j)%shadow == 1) cloud_mask_1km(2,i,j) = ibset(cloud_mask_1km(2,i,j),2)

!    Quality Assurance 1KM, byte 1 ----------------------------------------------------------------------------

     if(processing_information(i,j)%optical_thickness_general == 1)  &
           quality_assurance_1km(1,i,j) =  ibset(quality_assurance_1km(1,i,j),0)
     if(processing_information(i,j)%optical_thickness_confidence == 1) then
           quality_assurance_1km(1,i,j) =  ibset(quality_assurance_1km(1,i,j),1)
     elseif(processing_information(i,j)%optical_thickness_confidence == 2) then
           quality_assurance_1km(1,i,j) =  ibset(quality_assurance_1km(1,i,j),2)
     elseif(processing_information(i,j)%optical_thickness_confidence == 3) then
           quality_assurance_1km(1,i,j) =  ibset(quality_assurance_1km(1,i,j),1)
           quality_assurance_1km(1,i,j) =  ibset(quality_assurance_1km(1,i,j),2)
     endif

     if(processing_information(i,j)%optical_thickness_outofbounds== 1) then
           quality_assurance_1km(1,i,j) =  ibset(quality_assurance_1km(1,i,j),3)
     elseif(processing_information(i,j)%optical_thickness_outofbounds== 2) then
           quality_assurance_1km(1,i,j) =  ibset(quality_assurance_1km(1,i,j),4)
     elseif(processing_information(i,j)%optical_thickness_outofbounds == 3) then 
           quality_assurance_1km(1,i,j) =  ibset(quality_assurance_1km(1,i,j),3)
           quality_assurance_1km(1,i,j) =  ibset(quality_assurance_1km(1,i,j),4)
     endif

     if(processing_information(i,j)%effective_radius_general == 1) &
           quality_assurance_1km(1,i,j) =  ibset(quality_assurance_1km(1,i,j),5)
 
     if(processing_information(i,j)%effective_radius_confidence== 1) then
           quality_assurance_1km(1,i,j) =  ibset(quality_assurance_1km(1,i,j),6)
     elseif(processing_information(i,j)%effective_radius_confidence== 2) then
           quality_assurance_1km(1,i,j) =  ibset(quality_assurance_1km(1,i,j),7)
     elseif(processing_information(i,j)%effective_radius_confidence==3) then
           quality_assurance_1km(1,i,j) =  ibset(quality_assurance_1km(1,i,j),6)
           quality_assurance_1km(1,i,j) =  ibset(quality_assurance_1km(1,i,j),7)
     endif

!    Quality Assurance 1KM, byte 2 ----------------------------------------------------------------------------

     if(processing_information(i,j)%water_path_general== 1) &
           quality_assurance_1km(2,i,j) =  ibset(quality_assurance_1km(2,i,j),0)

     if(processing_information(i,j)%water_path_confidence== 1) then
           quality_assurance_1km(2,i,j) =  ibset(quality_assurance_1km(2,i,j),1)
     elseif(processing_information(i,j)%water_path_confidence== 2) then
           quality_assurance_1km(2,i,j) =  ibset(quality_assurance_1km(2,i,j),2)
     elseif(processing_information(i,j)%water_path_confidence==3) then
           quality_assurance_1km(2,i,j) =  ibset(quality_assurance_1km(2,i,j),1)
           quality_assurance_1km(2,i,j) =  ibset(quality_assurance_1km(2,i,j),2)
     endif

     if(     processing_information(i,j)%cloud_retrieval_proc_path_1621 == 1 ) then
          quality_assurance_1km(2,i,j) = ibset(quality_assurance_1km(2,i,j),3)
     elseif (processing_information(i,j)%cloud_retrieval_proc_path_1621 == 2 ) then
          quality_assurance_1km(2,i,j) = ibset(quality_assurance_1km(2,i,j),4)
     elseif (processing_information(i,j)%cloud_retrieval_proc_path_1621 == 3 ) then
          quality_assurance_1km(2,i,j) = ibset(quality_assurance_1km(2,i,j),3)
          quality_assurance_1km(2,i,j) = ibset(quality_assurance_1km(2,i,j),4)
     elseif (processing_information(i,j)%cloud_retrieval_proc_path_1621 == 4 ) then
          quality_assurance_1km(2,i,j) = ibset(quality_assurance_1km(2,i,j),5)
     endif

     if(processing_information(i,j)%retrieval_outcome_1621 == 1)  &
          quality_assurance_1km(2,i,j) = ibset(quality_assurance_1km(2,i,j),6)

!    Quality Assurance 1KM, byte 3 ----------------------------------------------------------------------------

     if(processing_information(i,j)%cloud_retrieval_processing_path ==1 ) then
          quality_assurance_1km(3,i,j) = ibset(quality_assurance_1km(3,i,j),0)  
     elseif (processing_information(i,j)%cloud_retrieval_processing_path == 2) then
          quality_assurance_1km(3,i,j) = ibset(quality_assurance_1km(3,i,j),1)  
     elseif (processing_information(i,j)%cloud_retrieval_processing_path == 3) then
          quality_assurance_1km(3,i,j) = ibset(quality_assurance_1km(3,i,j),0)  
          quality_assurance_1km(3,i,j) = ibset(quality_assurance_1km(3,i,j),1)   
     elseif (processing_information(i,j)%cloud_retrieval_processing_path == 4) then
          quality_assurance_1km(3,i,j) = ibset(quality_assurance_1km(3,i,j),2)    
     endif

     if(processing_information(i,j)%retrieval_outcome == 1)  &
          quality_assurance_1km(3,i,j) = ibset(quality_assurance_1km(3,i,j),3) 

     if(processing_information(i,j)%rayleigh_correction == 1) &
          quality_assurance_1km(3,i,j) = ibset(quality_assurance_1km(3,i,j),4) 

     if(processing_information(i,j)%atmospheric_correction == 1)  &
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

     if (processing_information(i,j)%optical_thickness_1621_gen == 1) &
          quality_assurance_1km(4,i,j) = ibset(quality_assurance_1km(4,i,j),0)

     if (processing_information(i,j)%optical_thickness_1621_conf ==1)  then
          quality_assurance_1km(4,i,j) = ibset(quality_assurance_1km(4,i,j),1)
     elseif (processing_information(i,j)%optical_thickness_1621_conf ==2)  then
          quality_assurance_1km(4,i,j) = ibset(quality_assurance_1km(4,i,j),2)
     elseif (processing_information(i,j)%optical_thickness_1621_conf ==3)  then
          quality_assurance_1km(4,i,j) = ibset(quality_assurance_1km(4,i,j),1)
          quality_assurance_1km(4,i,j) = ibset(quality_assurance_1km(4,i,j),2)
     endif

     if (processing_information(i,j)%effective_radius_1621_gen== 1) &
          quality_assurance_1km(4,i,j) = ibset(quality_assurance_1km(4,i,j),3)

     if (processing_information(i,j)%effective_radius_1621_conf==1)  then
          quality_assurance_1km(4,i,j) = ibset(quality_assurance_1km(4,i,j),4)
     elseif (processing_information(i,j)%effective_radius_1621_conf==2)  then
          quality_assurance_1km(4,i,j) = ibset(quality_assurance_1km(4,i,j),5)
     elseif (processing_information(i,j)%effective_radius_1621_conf==3)  then
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

     if (processing_information(i,j)%water_path_1621_gen == 1) &
          quality_assurance_1km(5,i,j) = ibset(quality_assurance_1km(5,i,j),0)

     if (processing_information(i,j)%water_path_1621_conf==1)  then
          quality_assurance_1km(5,i,j) = ibset(quality_assurance_1km(5,i,j),1)
     elseif (processing_information(i,j)%water_path_1621_conf==2)  then
          quality_assurance_1km(5,i,j) = ibset(quality_assurance_1km(5,i,j),2)
     elseif (processing_information(i,j)%water_path_1621_conf==3)  then
          quality_assurance_1km(5,i,j) = ibset(quality_assurance_1km(5,i,j),1)
          quality_assurance_1km(5,i,j) = ibset(quality_assurance_1km(5,i,j),2)
     endif

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
     if(processing_information(i,j)%retrieval_outcome == 1)  &
          quality_assurance_1km(5,i,j) = ibset(quality_assurance_1km(5,i,j),6)

!    Quality Assurance 1KM, byte 6 ----------------------------------------------------------------------------
		quality_assurance_1km(6,i,j) = processing_information(i,j)%ml_test_mark
	

  enddo
 enddo   
  
 end subroutine convert_binary_qa
 

 
 end module shared_io_module_hdf5
