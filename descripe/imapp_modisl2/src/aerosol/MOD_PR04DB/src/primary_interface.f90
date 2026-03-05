module primary_interface

! We don't really need this, but it's a good practise to get into. It should make compiling and
! debugging the interfaces easier
! interface is good for compilation purpose, easier to detect bugs for compiler
! provide the subroutine with uncertain arrays
!
	
       interface
           subroutine allocate_arrays ( edge, status )
              use GeneralAuxtype
              use core_arrays
              use hdf
              implicit none
              integer, dimension(:), intent(in)     :: edge
              integer,               intent(out)    :: status
           end subroutine allocate_arrays
        end interface


        interface
           subroutine get_modis_data_cube ( level1b_name,geolocation_name, start, edge, stride, spbands, status)
                  use GeneralAuxType
                  use core_arrays
                  use hdf
		  implicit none
                  integer, dimension (2), intent (in)       :: start, edge, stride
                  integer, dimension (:), intent (in)       :: spbands
                  character(*),           intent (in)       :: level1b_name, geolocation_name
                  integer,                intent (out)      :: status
           end subroutine get_modis_data_cube 
        end interface

		 
end module primary_interface

