! In this module important paramters that address program control
! and messaging are set

module nonscience_parameters

  implicit none
 
  integer :: success, warning, error, failure  


  parameter (success = 0, warning = 1, error = 2, failure = 3)

  real,parameter :: fillvalue_real = -999.
  integer*2, parameter :: fillvalue_int2 = -9999

  real,parameter ::  MAX_TAU_RTRIEVED = 150.0  
  real,parameter ::  MAX_COST_FUNCTION = 2.0  !not in percent

  integer, parameter :: Max_Scan_Per_Granule = 250

  integer, parameter :: MAX_SDS_NAME_LEN      = 64, &
                        MAX_FILE_NAME_LEN     = 1000 

end module nonscience_parameters
