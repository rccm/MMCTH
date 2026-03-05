      program Driver_MOD_PR06OD

!-------------------------------------------------------------------
!
!!F77
!
!!Description:
! Driver program for MODIS cloud property algorithm.
!
!!Input Parameters:
! None
!
!!Output Parameters:
! None
!
!!Revision History:
! Revision 3.0  2004/08/24 mag
! the new driver
!
!!Team-unique Header:
! Cloud Retrieval Group, NASA/GSFC, Code 913, Greenbelt, Maryland, USA
!
!!Credit and Reference:
! Programmed by:
! Mark Gray (mag)
! L-3 GSI
! NASA/GSFC, code 913
! Greenbelt, Maryland, USA
!
!!End
!---------------------------------------------------------------------
	use names
	use mod06_run_settings

      implicit none
 
!      include 'PGS_SMF.f'
!      include 'PGS_PC.f'
!      include 'PGS_MODIS_39500.f'
!      include 'mapi.inc'

!.....Parameter Declarations


!*****************************************************
!     Variable Declarations
!*****************************************************


!.....ECS inventory metadata variables

	  character(len = 200) :: par_file, dummy
	character(len=300) :: dirname, dayname, yearname, timename, outdirname
	character(len=500) :: temp
	real :: statistics(7)
	integer :: ErrorFlag
	integer :: mylun, i, idx

!--------------------------------------------------------------------------------
! Begin executable code by writing "Begin" message to LogStatus file 
!--------------------------------------------------------------------------------
	
	mylun = 8


	call getarg(1, par_file)
	print*, par_file

	open (unit=mylun, file=trim(par_file), status='old', form='formatted')
	read(mylun, *) dummy
	read(mylun, *) yearname, dayname, timename
	
	read(yearname,*) MYYEAR
	read(dayname, *) MYDAY
	read(timename, *) MYTIME
	dirname = './inputs-' // trim(yearname) // trim(dayname) // '.' // trim(timename) // '/'
	outdirname = './outputs-' // trim(yearname) // trim(dayname) // '.' // trim(timename) // '/OD/'
	print*, trim(dirname) 
	
	read(mylun, *) dummy
	read(mylun, *) temp
	Alevel1b_name(1) = trim(dirname) // trim(temp)
	print*, trim(Alevel1b_name(1))

	read(mylun, *) dummy
	read(mylun, *) temp
	Acloudmask_name = trim(dirname) // trim(temp)
	print*, trim(Acloudmask_name)
	
	read(mylun, *) dummy
	read(mylun, *) temp
	Ageolocation_name = trim(dirname) // trim(temp)
	print*, trim(Ageolocation_name)
	
	
	read(mylun, *) dummy
	read(mylun, *) temp
	Amod06_name = trim(outdirname) // trim(temp)
	print*, trim(Amod06_name)
!	MY_TEXT_FILE = trim(Amod06_name) // ".txt"
	
	read(mylun, *) dummy
	read(mylun, *) temp
	Agdas_name = trim(dirname) // trim(temp)
	read(mylun, *) temp
	Agdas_name2 = trim(dirname) // trim(temp)
	print*, trim(Agdas_name)
	print*, trim(Agdas_name2)

#ifdef USE_TOAST
	read(mylun, *) dummy
	read(mylun, *) temp
	Aozone_name = trim(dirname) // trim(temp)
	print*, trim(Aozone_name)
#else
!	Aozone_name = "none"
#endif	
	
	read(mylun, *) dummy
	read(mylun, *) temp
	Ancepice_name = trim(dirname) // trim(temp)
	print*, trim(Ancepice_name)
	
	read(mylun, *) dummy
	read(mylun, *) temp
	Anise_name = trim(dirname) // trim(temp)
	print*, trim(Anise_name)

! skip the CSR bias file name if present
	
	read(mylun, *) dummy
	print*, trim(dummy)
	idx = index (dummy, "CSR")
	if (idx /= 0) then 
		read(mylun,*) dummy
		read(mylun,*) dummy
	endif


	read(mylun, '(a)') temp
	Aice_library =  trim(temp)
	read(mylun, '(a)') temp
	Awater_library = trim(temp)
	
	do i=1, 3

		read(mylun, '(a)') temp
		Alibnames_ice(i) =  trim(temp)
		read(mylun, '(a)') temp
		Alibnames_ice_sdev(i) =  trim(temp)
	
	end do
	do i=1, 3

		read(mylun, '(a)') temp
		Alibnames_water(i) =  trim(temp)
		read(mylun, '(a)') temp
		Alibnames_water_sdev(i) =  trim(temp)
	
	end do


	read(mylun, '(a)') temp
	Aphase_library = trim(temp)
    print*, "ICE: ", trim(Aice_library)
	print*, trim(Awater_library)
	do i=1, 3
		print*, trim(Alibnames_ice(i))
		print*, trim(Alibnames_ice_sdev(i))
		print*, trim(Alibnames_water(i))
		print*, trim(Alibnames_water_sdev(i))
	end do

	print*, trim(Aphase_library)
	
	read(mylun, *) dummy
	read(mylun, '(a)') temp
	Atransmittance_library =  trim(temp)
	print*, trim(Atransmittance_library)
	
	read(mylun, *) dummy
	read(mylun, '(a)') temp
	Aecosystem_data_name =  trim(temp)
	print*, trim(Aecosystem_data_name)

	read(mylun, *) dummy
	read(mylun, '(a)') temp
	Asnowicealbedo_data_name =  trim(temp)
	print*, trim(Asnowicealbedo_data_name)

	read(mylun, *) dummy
	read(mylun, '(a)') temp
	Asurfacealbedo_lib_659 =  trim(temp)
	read(mylun, '(a)') temp
	Asurfacealbedo_lib_858 =  trim(temp)
	read(mylun, '(a)') temp
	Asurfacealbedo_lib_124 =  trim(temp)
	read(mylun, '(a)') temp
	Asurfacealbedo_lib_164 =  trim(temp)
	read(mylun, '(a)') temp
	Asurfacealbedo_lib_21 =  trim(temp)

	read(mylun, *) dummy
	read(mylun, '(a)') temp
	ACT_lib_path = trim(temp)
	read(mylun, *) dummy
	read(mylun, '(a)') temp
	Aemissivity_name = trim(temp)

	
	print*, trim(Asurfacealbedo_lib_659)
	print*, trim(Asurfacealbedo_lib_858)
	print*, trim(Asurfacealbedo_lib_124)
	print*, trim(Asurfacealbedo_lib_164)
	print*, trim(Asurfacealbedo_lib_21)
	
	
	close(mylun)


!...........Perform cloud optical properties retrieval


               call MOD_PR06OD(statistics,ErrorFlag)



      End program Driver_MOD_PR06OD
