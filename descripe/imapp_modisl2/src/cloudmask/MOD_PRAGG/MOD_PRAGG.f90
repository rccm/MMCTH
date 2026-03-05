
!*************************************************************************** 
!*************************************************************************** 
!*************************************************************************** 
MODULE module_mod02
  
	IMPLICIT NONE
	SAVE

	include 'hdf.f90'
	include 'dffunc.f90'
        include 'modagg.inc'
        include 'PGS_MODIS_39500.f'
        include 'PGS_SMF.f'

	INTEGER, PARAMETER  :: rows = 4*1354
	! changed this to 2200 to deal with potentially larger granule sizes occuring in Terra MOD06
	! RB 11/15/08
        ! changed this to 5500 to deal with potentially larger granule sizes occuring in direct broadcast 
        ! RMC 1/28/14
	INTEGER, PARAMETER  :: cols = 4*5500

	!*************************************************************************** 
	!*************************************************************************** 
	!*************************************************************************** 
	TYPE modis_data_field
	   INTEGER*2, DIMENSION(rows,cols,2)  :: values
	   INTEGER  				              :: size_x
	   INTEGER  				              :: size_y
	END TYPE modis_data_field

	CONTAINS

	!*************************************************************************** 
	!*************************************************************************** 
	!*************************************************************************** 
	SUBROUTINE read_field(sd_id,fieldname,datafield) 

	   CHARACTER(len=*), INTENT(IN) 			:: fieldname
	   INTEGER(4)      , INTENT(IN) 			:: sd_id
	   TYPE(modis_data_field), INTENT(out)   	:: datafield

	   INTEGER			:: 	istat,sds_index,sds_id,RANK,err,num_attrs,data_type,count,attr_index
	   integer, DIMENSION(3)				        ::  dims, start, edges, stride,dimsizes
	   REAL*8								        :: scale_Factor,add_offset
	   character(len=256)					        ::  sds_name,attr_name
	   INTEGER*2,DIMENSION(:,:,:),ALLOCATABLE       :: i_dummy_data_field
           
           character(len=256) :: msg
           character(len=256) :: location
    
           location = 'read_field, MOD_PRAGG.f90'

	   sds_index   = sfn2index(sd_id,fieldname)
	   sds_id	   = sfselect(sd_id, sds_index)

	   err = sfginfo(sds_id, sds_name, rank,dimsizes,data_type,num_attrs)
	   	   
	   start(1)    = 0;
	   start(2)    = 0;
	   start(3)    = 0;
	   stride(1)   = 1;
	   stride(2)   = 1;
	   stride(3)   = 1;
	   edges(1)    = dimsizes(1)
	   edges(2)    = dimsizes(2)
	   edges(3)    = dimsizes(3)

	   write(*,*) ' Reading         : ', fieldname
	   write(*,*) '    Type         : ', data_type
	   write(*,*) '    Rank         : ', rank
	   write(*,*) '    Dimsize      : ', dimsizes

	   !c  read the $#@$-HDF-4 Data
	   ALLOCATE(i_dummy_data_field(dimsizes(1),dimsizes(2),2))

	   istat = sfrdata (sds_id, start, stride, edges, i_dummy_data_field)

	! Now write out
!	OPEN(1,FILE='csd.f77.dat',FORM='UNFORMATTED')
!	WRITE(1)i_dummy_data_field
!	CLOSE(1) 

	   datafield%values(1:dimsizes(1),1:dimsizes(2),1:2) = i_dummy_data_field
	   datafield%size_x=dimsizes(1)
	   datafield%size_y=dimsizes(2)

	   DEALLOCATE(i_dummy_data_field)
 
	   !c  End the access to the open variable
	   istat = sfendacc(sds_id)
	   if (istat .eq. FAIL) then
                   write(msg, "(A)") "Error: Error Running sfendacc in read_field soubroutine"
                   call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC, msg, location)
    		   istat=heprnt(0)
	   end if  

	END SUBROUTINE read_field

	!*************************************************************************** 
	!*************************************************************************** 
	!*************************************************************************** 
	SUBROUTINE read_data(filename,sds_name,hr_data) 

	   CHARACTER(len=*), INTENT(IN) 	    ::  filename,sds_name
	   TYPE(modis_data_field),  INTENT(out)	::	hr_data		
	   INTEGER	                            ::	sd_id,istat,sds_index,sds_id,RANK
           character(len=256) :: msg
           character(len=256) :: location
           
           location = 'read_data, MOD_PRAGG.f90'


		! open HDF File	
    	sd_id = sfstart(filename, DFACC_RDONLY)
		if (sd_id .eq. FAIL) then
                           write(msg, "(A,A,A)") "Error: Could not open input file:", filename, " for reading"
                           call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC, msg, location)
			   istat=heprnt(0)
		 else
				print*,'... sd.hdf opened with READ access'
    	end if 

		CALL read_field(sd_id,sds_name,hr_data)	

		!  Terminate access to the SD interface and close the file.
    	istat = sfend(sd_id)
    	if (istat .eq. FAIL) then
                         write(msg, "(A,A)") "Error: Could not terminate access for file:", filename
                         call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC, msg, location)
    	   		 istat=heprnt(0)
		 else
				print*,'... sd.hdf closed'
    	end if	

	END SUBROUTINE read_data

	!*************************************************************************** 
	!*************************************************************************** 
	!*************************************************************************** 
	 SUBROUTINE convolve_hr_lr(hr_data,lr_data) 

	   TYPE(modis_data_field),  INTENT(IN)   ::  hr_data	 
	   TYPE(modis_data_field),  INTENT(OUT)  ::  lr_data 
	   	
	   INTEGER                           	 ::  i,j,k,x_ind,y_ind
	   REAL                                  ::  value(8,8),kernel(8,8)

           ! Ridgway 6/24/2013 updated Wolfe Kernel based on later MCST data

	   REAL ::  ref_kernel(8,8) = RESHAPE(   &
            (/ 0.0000,       0.0000,       0.0000,       0.0000,       0.0000,       0.0000,       0.0000,       0.0000,         &
            0.0000,       0.0000,       0.0000,       0.0000,       0.0000,       0.0000,       0.0000,       0.0000,         &
            0.0048,       0.0104,       0.0162,       0.0212,       0.0179,       0.0122,       0.0065,       0.0009,         &
            0.0132,       0.0289,       0.0449,       0.0588,       0.0497,       0.0340,       0.0180,       0.0026,         &
            0.0132,       0.0289,       0.0449,       0.0588,       0.0497,       0.0340,       0.0180,       0.0026,         &
            0.0132,       0.0289,       0.0449,       0.0588,       0.0497,       0.0340,       0.0180,       0.0026,         &
            0.0085,       0.0185,       0.0287,       0.0376,       0.0318,       0.0217,       0.0115,       0.0017,         &
            0.0000,       0.0000,       0.0000,       0.0000,       0.0000,       0.0000,       0.0000,       0.0000      /), &
            (/ 8, 8/))

	   ! RB 2/19/2013 added modified kernel for use with every 10th scan line

           ! Ridgway 6/24/2013 updated Wolfe Kernel based on later MCST data

    REAL ::  ref_kernel_10(8,8) = RESHAPE(   &
         (/ 0.0000,       0.0000,       0.0000,       0.0000,       0.0000,       0.0000,       0.0000,       0.0000,         &
         0.0000,       0.0000,       0.0000,       0.0000,       0.0000,       0.0000,       0.0000,       0.0000,         &
         0.0048,       0.0104,       0.0162,       0.0212,       0.0179,       0.0122,       0.0065,       0.0009,         &
         0.0132,       0.0289,       0.0449,       0.0588,       0.0497,       0.0340,       0.0180,       0.0026,         &
         0.0132,       0.0289,       0.0449,       0.0588,       0.0497,       0.0340,       0.0180,       0.0026,         &
         0.0,          0.0,          0.0,          0.0,          0.0,          0.0,          0.0,          0.0           , &
         0.0,          0.0,          0.0,          0.0,          0.0,          0.0,          0.0,          0.0           , &
         0.0,          0.0,          0.0,          0.0,          0.0,          0.0,          0.0,          0.0         /), &
         (/ 8, 8/))
    
	   ref_kernel = ref_kernel / SUM(ref_kernel)

	   ! RB 2/19/2013 need to normmalize modified kernel
	   ref_kernel_10 = ref_kernel_10 / SUM(ref_kernel_10)
	   	   
           ! RB 2/19/2013 CHANGED THIS TO START AT j=2 to avoid negative HiRes indices for first scan position 
           ! Ridgway 3/12/2013 updated to [DO j=2,lr_data%size_y-1] to keep y_ind-4 and y_ind+3 in bounds
           ! Ridgway first, last scan and 1st pixel in each scan default to previous EV_250_Aggr1km_RefSB values 

	   DO k=1,2
		 DO j=2,lr_data%size_y-1
	           y_ind = j * 4
	           DO i=2,lr_data%size_x 

	   ! RB Select correct kernel - use modified kernel for every 10th scan line corresponding to last lo_res detector
			IF (MOD(j,10).EQ.0) THEN 
			      kernel = ref_kernel_10 
			ELSE 
			      kernel = ref_kernel 
			ENDIF

 	              x_ind = (i-1) * 4 + 1
			      value = REAL(RESHAPE(hr_data%Values(x_ind-4:x_ind+3,y_ind-4:y_ind+3,k),(/8,8/)))
				  WHERE(value .LE. 0.0) kernel = 0.0
				  IF (SUM(kernel) .GT. 0.0) lr_data%values(i,j,k) = INT((SUM(value*kernel) / SUM(kernel)) + 0.5)
	    	ENDDO
		 ENDDO
	   ENDDO


	END SUBROUTINE convolve_hr_lr
!*************************************************************************** 
 SUBROUTINE write_lr(outfile,sds_name,lr_data)
!*************************************************************************** 

	   character(len=256), 	  INTENT(in)		::  OUTFILE
	   character(len=*), 	  INTENT(in)		::  sds_name
	   TYPE(modis_data_field),INTENT(in)		::	lr_data

 
 	   INTEGER*2,DIMENSION(:,:,:),ALLOCATABLE       :: idum
	   INTEGER			:: 	sd_id,istat,sds_index,sds_id,RANK,err,num_attrs,data_type,count,attr_index
	   integer, DIMENSION(3)				        ::  dims, start, edges, stride,dimsizes

           character(len=256) :: msg
           character(len=256) :: location
    
           location = 'write_lr, MOD_PRAGG.f90'
    
	   ALLOCATE(idum(lr_data%size_x,lr_data%size_y,2))
      idum=lr_data%values(1:lr_data%size_x,1:lr_data%size_y,1:2)

	   sd_id	 = sfstart  (outfile,DFACC_WRITE)
           if(sd_id .eq. FAIL) then
              write(msg, "(A,A,A)") "Error: Could not open output file:", outfile, " for writting"
              call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC, msg, location)
           end if
	   sds_index = sfn2index(sd_id,sds_name)
	   sds_id    = sfselect (sd_id, sds_index)

      ! Define the dimensions of the array to be created (N)
	   start(1)    = 0;
	   start(2)    = 0;
	   start(3)    = 0;
	   stride(1)   = 1;
	   stride(2)   = 1;
	   stride(3)   = 1;
	   edges(1)    = lr_data%size_x
	   edges(2)    = lr_data%size_y
	   edges(3)    = 2

	   write(*,*) ' Writing         : ', sds_name
	   write(*,*) '    Dimsize      : ', dimsizes

     
	   istat  = sfwdata(sds_id, start,stride,edges,idum)

!      End access to the hdf file
       istat = sfend(sd_id)
	   DEALLOCATE(idum)

     END SUBROUTINE write_lr

END MODULE module_mod02 


!*************************************************************************** 
!*************************************************************************** 
!*************************************************************************** 

	PROGRAM SAMPLE_NEW
	!
	!
	! Ralf Bennartz, 2/10/2009
	!

    USE module_mod02

!	CHARACTER(len=256)		 ::	MYD02_QKM_FILE,MYD02_1KM_FILE
	TYPE(modis_data_field)   :: hr_data
	TYPE(modis_data_field)   :: lr_data


        INTEGER fileid, version, rtn, status, pgs_pc_getreference
        character(len=256) :: MYD02_1KM_FILE, MYD02_QKM_FILE
        character(len=256) :: msg
        character(len=256) :: location
        
        location = 'module_mod02, MOD_PRAGG.f90'

        fileid = LRN_L1B_1km
        version = 1
        MYD02_1KM_FILE = ' '

        rtn = pgs_pc_getreference( fileid, version, MYD02_1KM_FILE)
        if(rtn .ne. PGS_S_SUCCESS) then
           write(msg, "(A,I6)") "Error from PGS_PC_GetReference() obtaining the filename for LUN", fileid
           call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC, msg, location)
        end if
        
        fileid = LRN_L1B_250
        version = 1
        MYD02_QKM_FILE = ' '

        rtn = pgs_pc_getreference(fileid, version, MYD02_QKM_FILE)
        if(rtn .ne. PGS_S_SUCCESS) then
           write(msg, "(A,I6)") "Error from PGS_PC_GetReference() obtaining the filename for LUN", fileid
           call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC, msg, location)
        end if

!	! Get filenames	
!	CALL getarg(1,MYD02_QKM_FILE)
!	CALL getarg(2,MYD02_1KM_FILE)
		
	! Read input data
	CALL read_data(MYD02_QKM_FILE,"EV_250_RefSB"        ,hr_data)	
	CALL read_data(MYD02_1KM_FILE,"EV_250_Aggr1km_RefSB",lr_data)	
	
	! convolve HR --> LR
	CALL convolve_hr_lr(hr_data,lr_data)
	
	CALL write_lr(MYD02_1KM_FILE,"EV_250_Aggr1km_RefSB",lr_data)	

    END
	
