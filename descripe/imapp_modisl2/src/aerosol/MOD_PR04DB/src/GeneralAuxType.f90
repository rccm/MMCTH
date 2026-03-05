  module GeneralAuxType
!-----------------------------------------------------------------------
! !F90
!
! !DESCRIPTION:
!
!    GeneralAuxType provides data type declarations and several
!    parameter definitions.
!
! !INPUT PARAMETERS:  none
!
! !OUTPUT PARAMETERS:  none
!
! !REVISION HISTORY:
!
!    Initial Version by Jennifer Wei, July 8, 2004
!    Jeremy Warner added prolog, December 1, 2006
!
! !TEAM-UNIQUE HEADER:
!
!    This software is developed by the Deep Blue Science Team
!    for the National Aeronautics and Space Administration,
!    Goddard Space Flight Center, under contract NAS5-02041.
!
! !REFERENCES AND CREDITS
!
! !DESIGN NOTES:
!
! !END
!-----------------------------------------------------------------------

!... symbolic names for data types:

 integer, parameter :: integer_fourbyte= selected_int_kind(9)
 integer, parameter :: integer_twobyte = selected_int_kind(4)
 integer, parameter :: integer_onebyte = selected_int_kind(2)
 integer, parameter :: single = kind(1.0)
 integer, parameter :: double = kind(1.0d0)
 integer, parameter :: singlecomplex = kind((1.0,1.0))
 integer, parameter :: doublecomplex = kind((1.0d0,1.0d0))
 integer, parameter :: logical = kind(.true.)


!... useful mathematical constants

 real(single), parameter :: pi=3.141592653589793238462643383279502884197_single
 real(single), parameter :: pio2=1.57079632679489661923132169163975144209858_single
 real(single), parameter :: twopi=6.283185307179586476925286766559005768394_single
 real(single), parameter :: d2r = 0.017453292519943295_single

 real(double), parameter :: pi_d=3.141592653589793238462643383279502884197_double
 real(double), parameter :: pio2_d=1.57079632679489661923132169163975144209858_double
 real(double), parameter :: twopi_d=6.283185307179586476925286766559005768394_double
 real(double), parameter :: d2r_d = 0.017453292519943295_double

 real(single),       parameter :: missingvalue = -1.0

end module GeneralAuxType

