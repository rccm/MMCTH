 module modis_cloudhandler

 private

 public :: modis_cloudprocessing

 integer*1, allocatable, dimension(:,:,:) ::  cloud_mask_1km

 contains

 subroutine modis_cloudprocessing(cloudmask_filedata, mod06_filedata, current_scan,start, stride, edge, debug, status)

   use GeneralAuxType, only: integer_fourbyte, integer_onebyte
   use core_arrays
   use modis_io_module, only: writeqaarray_toolkit, write_int2_array
    
   implicit none
   
   integer, dimension(:), intent(in) :: cloudmask_filedata, mod06_filedata
   integer, intent(in)                                  :: current_scan
   integer(integer_fourbyte), dimension(2), intent(in)  :: start, stride, edge
   logical, intent(in)                                  :: debug
   integer(integer_fourbyte), intent(out)               :: status


   integer(integer_onebyte), allocatable, dimension(:,:,:) :: packedcloudmask, packedcloudmaskQA
   integer(integer_fourbyte) :: xdimension, ydimension, cloudmask_bits, checkvariable,i,j,cloudmaskQA_bits

	
   cloudmask_bits = 48
   cloudmaskQA_bits = 80

   xdimension =  edge(1)  
   ydimension =  edge(2) 
 
   allocate (packedcloudmask(xdimension,ydimension,cloudmask_bits), stat = checkvariable)
   if ( checkvariable /= 0 ) then
      status=1
   endif
   allocate (packedcloudmaskQA(xdimension,ydimension,cloudmaskQA_bits), stat = checkvariable)
   if ( checkvariable /= 0 ) then
      status=1
   endif

   allocate(cloud_mask_1km(2, xdimension, ydimension))
   

   ! read libraries in preparation for later interpolation
   call modis_readcloudmask( cloudmask_filedata,    &
                             current_scan,      &
                             start,             &
                             edge,              &
                             stride,            &
                             packedcloudmask,   &
                             packedcloudmaskQA, &
                             debug, status)

   call modis_unpackcloudmask ( packedcloudmask, packedcloudmaskQA, debug, status)

   deallocate(packedcloudmask, stat = checkvariable)
   if ( checkvariable /= 0 ) then
      status=1
   endif
   deallocate(packedcloudmaskQA, stat = checkvariable)
   if ( checkvariable /= 0 ) then
      status=1
   endif
   
! write it to disk here and forget about it for the rest of the run  
   call writeqaarray_toolkit(mod06_filedata,  &
                           current_scan, 1,         &
                           cloud_mask_1km,   &
                           'Cloud_Mask_1km', &
                           status )
   
   deallocate(cloud_mask_1km)

! same goes for cloud_mask_SPI

    call write_int2_array(mod06_filedata, "Cloud_Mask_SPI", &
    						(/ 0, start(1), start(2) /), &
    						(/ 1, stride(1),  stride(2) /), &
    						(/ 2, edge(1),    edge(2) /), &
    						cloud_mask_SPI, status) 
	print*, "written: ", status  
 
 end subroutine modis_cloudprocessing
!f90 modis_readcloudmask
!
!Description:
!
!  Read MAS 48 bit cloud mask from a hdf file previously opened for reading.
!  Convert individual tests from bits to short integers and pass short integer
!  array out.
!
!input parameters:
!
!output parameters:
!
!revision history:
!
! v.1 September 2001 Initial work mag gray@climate.gsfc.nasa.gov
!
!team-unique header:
!
! Cloud Retrieval Group, NASA GSFC, Greenbelt, Maryland, USA
!
!references and credits:
!
! Mark Gray
! gray@climate.gsfc.nasa.gov
! EmergentIT
! Code 913, NASA GSFC
! Greenbelt, MD 20771
!
!----------------------------------------------------------------------

 subroutine modis_readcloudmask(filedata, current_scan,start, edge, stride, mydata, QAdata, debug,status)
   use nonscience_parameters
   use GeneralAuxType,only: integer_onebyte,  integer_fourbyte, single
   use core_arrays, only: cloud_mask_SPI
   use hdf_mod 
   implicit none

   integer, intent(in) :: current_scan
   integer(integer_fourbyte),  dimension (2),   intent(in)     :: start,edge,stride
   integer, dimension(:),                   intent(in)                  :: filedata
   integer(integer_onebyte) ,  dimension(:,:,:),intent(out)    :: mydata, QAdata
   integer(integer_fourbyte),                   intent(out)    :: status
   logical, intent(in)                                  :: debug

   integer(integer_fourbyte)                       :: sds_id,sds_index,attr_index, hdfid
   integer(integer_fourbyte),dimension (3)         :: dimsizes, localstart,localedge,localstride
   integer(integer_fourbyte)                       :: rank,  number_type, nattrs
   character(MAX_SDS_NAME_LEN)                                   :: local_sds_name,sds_name
   integer                                         :: xdimension, ydimension , checkvariable, &
                                                      numberofbands, index,cloudmask_bytes
   integer                                         :: bit, byte, i, j, output_bits, cloudmaskQA_bytes
   integer(integer_onebyte), allocatable, dimension(:,:,:)     :: rawbytedata , rawQAbytedata, tempbytedata
!   integer*2, dimension(:,:,:), allocatable :: tempSPI
   
   integer ::  DS_Dim1_Cldmsk, DS_Dim1_Cldmsk_QA,DS_Dim2,DS_Dim3
   logical :: Error_Flag
   
   integer :: gfile_id, gvar_id, gerr_code, iii
   integer*1, parameter :: seven = 7
   

!   hdfid = sfstart(filename,DFACC_READ)

!  read the cloudmask sds

   xdimension =   edge(1)
   ydimension =   edge(2)


   cloudmask_bytes=6
   localstart  = (/ start(1), start(2), 0/)
   localedge   = (/ edge(1),    edge(2), cloudmask_bytes/)
   localstride = (/ stride(1),  stride(2), 1/)
  
   allocate (rawbytedata(cloudmask_bytes,  xdimension,ydimension ),stat=checkvariable)
   allocate(tempbytedata(xdimension, ydimension, cloudmask_bytes))
   
!	allocate(tempSPI(2, xdimension, ydimension))

!  read the cloudmask QA sds

   cloudmaskQA_bytes = 10

   allocate (rawQAbytedata(cloudmaskQA_bytes, xdimension,ydimension ),stat=checkvariable)

   gfile_id = filedata(1)
   
   gvar_id = sfselect(gfile_id, sfn2index(gfile_id, "Cloud_Mask"))
   gerr_code = sfrdata(gvar_id, localstart, localstride, localedge, tempbytedata)
   if(gerr_code == -1) then
	gerr_code = sfendacc(gvar_id)
	gerr_code = sfend(gfile_id)
    status = failure
    call local_message_handler('Failure detected in Cloud Mask read. Check previous message/s', status,'modis_readcloudmask')
   endif 
   gerr_code = sfendacc(gvar_id)
   
   do iii=1, cloudmask_bytes
	rawbytedata(iii, :, :) = tempbytedata(:, :, iii)
   end do
   
   deallocate(tempbytedata)
   
! save the first two bytes of cloud mask for output to disk
! but zero out the top five bits of byte 2
   cloud_mask_1km(1, :,:) = rawbytedata(1, :,:)   
   cloud_mask_1km(2, :,:) = iand(rawbytedata(2, :,:), seven)
      
   localstart  = (/ 0, start(1), start(2)/)
   localedge   = (/ cloudmaskQA_bytes, edge(1),    edge(2)/)
   localstride = (/ 1, stride(1),  stride(2)/)
   
   
    gvar_id = sfselect(gfile_id, sfn2index(gfile_id, "Quality_Assurance"))
   gerr_code = sfrdata(gvar_id, localstart, localstride, localedge, rawQAbytedata)
   if(gerr_code == -1) then
	gerr_code = sfendacc(gvar_id)
	gerr_code = sfend(gfile_id)
    status = failure
    call local_message_handler('Failure detected in Cloud Mask QA read. Check previous message/s', status,'modis_readcloudmask')
   endif 
   gerr_code = sfendacc(gvar_id)
  
   localstart  = (/ 0, start(1), start(2) /)
   localedge   = (/ 2, edge(1),    edge(2) /)
   localstride = (/ 1, stride(1),  stride(2) /)

	
    gvar_id = sfselect(gfile_id, sfn2index(gfile_id, "Cloud_Mask_SPI"))
   gerr_code = sfrdata(gvar_id, localstart, localstride, localedge, cloud_mask_SPI)
   
    print*, gerr_code
   
   if(gerr_code == -1) then
	gerr_code = sfendacc(gvar_id)
	gerr_code = sfend(gfile_id)
    status = failure
    call local_message_handler('Failure detected in SPI read. Check previous message/s', status,'modis_readcloudmask')
   endif 
   gerr_code = sfendacc(gvar_id)


!	cloud_mask_SPI = 0.
!	
!	do j=1, ydimension
!		do i=1, xdimension
!		
!			cloud_mask_SPI(1,i,j) = tempSPI(2, i,j) / tempSPI(1,i,j) * 100.
!			cloud_mask_SPI(2,i,j) = tempSPI(4, i,j) / tempSPI(3,i,j) * 100.
!		
!		end do
!	end do


!	deallocate(tempSPI)


!  convert n x m x 6 byte(48 bit) data to n x m 48 integer array

   mydata(:,:,:) = 0

   do j = 1, ydimension
	do i = 1, xdimension
       output_bits = 0
       do byte = 1 , cloudmask_bytes
         do bit = 0, 7
            output_bits = output_bits + 1
            if ( btest(rawbytedata(byte, i, j ), bit) ) mydata(i, j, output_bits) = 1
         enddo   
       enddo
     enddo
   enddo

   QAdata(:,:,:) = 0
   do j = 1, ydimension
	  do i = 1, xdimension
       output_bits = 0
       do byte = 1 , cloudmaskQA_bytes
         do bit = 0, 7
            output_bits = output_bits + 1
            if ( btest(rawQAbytedata(byte, i, j ), bit) ) QAdata(i, j, output_bits) = 1
         enddo
       enddo
     enddo
   enddo
   deallocate(rawbytedata, stat = checkvariable)
   deallocate(rawQAbytedata, stat = checkvariable)

 end subroutine modis_readcloudmask
 
 
 subroutine modis_unpackcloudmask(packedcloudmask,packedcloudmaskQA,debug, status)
  use GeneralAuxType, only: integer_onebyte, integer_fourbyte
  use core_arrays, only: cloudmask, processing_information

  implicit none

  integer(integer_fourbyte), intent(out)                  :: status
  integer(integer_onebyte), intent(in), dimension(:,:,:)  :: packedcloudmask, packedcloudmaskQA
  logical, intent(in)                                     :: debug
  integer(integer_fourbyte)                               :: xdimension, ydimension, bytedimension,i,j, k, cm_count
  integer*1 :: pixel250(16), pixel250qa(16)
  status=0

  xdimension    = size(packedcloudmask,1)
  ydimension    = size(packedcloudmask,2)
  bytedimension = size(packedcloudmask,3)

! top level flags
  cloudmask(:,:)%cloudmask_determined = .false.  
  where(packedcloudmask(:,:,1) == 1) cloudmask(:,:)%cloudmask_determined = .true.

  cloudmask(:,:)%confident_cloudy = .false.
  cloudmask(:,:)%probablyclear_66 = .false.
  cloudmask(:,:)%probablyclear_95 = .false.
  cloudmask(:,:)%probablyclear_99 = .false.
  where(packedcloudmask(:,:,2) == 0 .and. packedcloudmask(:,:,3) ==0) &
           cloudmask(:,:)%confident_cloudy = .true.
  where(packedcloudmask(:,:,2) == 1 .and. packedcloudmask(:,:,3) ==0) &
           cloudmask(:,:)%probablyclear_66 = .true.
  where(packedcloudmask(:,:,2) == 0 .and. packedcloudmask(:,:,3) ==1) &
           cloudmask(:,:)%probablyclear_95 = .true.
  where(packedcloudmask(:,:,2) == 1 .and. packedcloudmask(:,:,3) ==1) &
           cloudmask(:,:)%probablyclear_99 = .true.

! processing path flags
 
  cloudmask(:,:)%night = 0
  where(packedcloudmask(:,:,4) == 0) cloudmask(:,:)%night = 1

  cloudmask(:,:)%sunglint = 0
  where(packedcloudmask(:,:,5) == 0) cloudmask(:,:)%sunglint = 1


  cloudmask(:,:)%snowice_surface = .false.
  where(packedcloudmask(:,:,6) == 0) cloudmask(:,:)%snowice_surface = .true.

  cloudmask(:,:)%water_surface = .false.
  cloudmask(:,:)%coastal_surface = .false.
  cloudmask(:,:)%desert_surface = .false.
  cloudmask(:,:)%land_surface =  .false.

  where(packedcloudmask(:,:,7) == 0 .and. packedcloudmask(:,:,8) == 0) &
           cloudmask(:,:)%water_surface= .true.
  where(packedcloudmask(:,:,7) == 1 .and. packedcloudmask(:,:,8) == 0) &
           cloudmask(:,:)%coastal_surface = .true.
  where(packedcloudmask(:,:,7) == 0 .and. packedcloudmask(:,:,8) == 1) &
           cloudmask(:,:)%desert_surface = .true.
  where(packedcloudmask(:,:,7) == 1 .and. packedcloudmask(:,:,8) == 1) &
           cloudmask(:,:)%land_surface = .true.

! **** FOR SCREWY CLOUD MASK IN TC4 ONLY
!  where(cloudmask(:,:)%sunglint .and. cloudmask(:,:)%probablyclear_95 &
!		.and. cloudmask(:,:)%water_surface) cloudmask(:,:)%confident_cloudy = .true.



! byte 1

  cloudmask(:,:)%test_high_138 =  0
  where(packedcloudmask(:,:,17) == 0) cloudmask(:,:)%test_high_138  = 1

  cloudmask(:,:)%test_visiblereflectance = 0
  where(packedcloudmask(:,:,21) == 0) cloudmask(:,:)%test_visiblereflectance= 1

  cloudmask(:,:)%test_visnirratio =  0
  where(packedcloudmask(:,:,22) == 0) cloudmask(:,:)%test_visnirratio  = 1

  !250m cloud flags
   do j = 1, ydimension
	do i = 1, xdimension
		
		cm_count = 16
		pixel250(:) = packedcloudmask(i,j,33:48)
		pixel250qa(:) = packedcloudmaskQA(i,j,33:48)
		
		do k = 1, 16 
			if (pixel250(k) == 1 .and. pixel250qa(k) == 1 ) &
				cm_count = cm_count - 1
		end do
	
		cloudmask(i,j)%visible_cloudtest_250m = cm_count
	
	end do
  end do


! QA bit fields

  cloudmask(:,:)%applied_highcloud138 = 0
  where(packedcloudmaskQA(:,:,17) == 1) cloudmask(:,:)%applied_highcloud138= 1

  cloudmask(:,:)%applied_visiblereflectance = 0
  where(packedcloudmaskQA(:,:,21) == 1) cloudmask(:,:)%applied_visiblereflectance= 1

  cloudmask(:,:)%applied_visnirratio = 0
  where(packedcloudmaskQA(:,:,22) == 1) cloudmask(:,:)%applied_visnirratio= 1

 end subroutine modis_unpackcloudmask

 end module modis_cloudhandler
