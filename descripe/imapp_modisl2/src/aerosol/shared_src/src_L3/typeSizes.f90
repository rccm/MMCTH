! !F90
!
! !Description:
!   Provide named kind parameters for use in declarations of real and integer 
!    variables with specific byte sizes (i.e. one, two, four, and eight byte
!    integers; four and eight byte reals). The parameters can then be used
!    in (KIND = XX) modifiers in declarations. Specific sizes are needed for
!    interface with the HDF libraries. 
!   A single function (byteSizesOK()) is provided to ensure that the selected 
!    kind parameters are correct.
!  
! !Input Parameters:
!   None.
!
! !Output Parameters:
!   Public parameters, fixed at compile time:
!     OneByteIntKind, TwoByteIntKind, FourByteIntKind, EightByteIntKind
!                                     FourByteRealKind, EightByteRadlKind
! 
! !Revision History:
!   $Log: typeSizes.f90,v $
! Revision 1.1  1997/06/03  16:36:21  pincus
! Initial revision
!
!
! !Team-Unique Header:
!
! !References and Credits:
!   Written by
!    Robert Pincus
!    Climate and Radiation Branch, Code 913
!    NASA/GSFC
!    Greenbelt MD 20771
!    Robert.Pincus@gsfc.nasa.gov
!
!
! !Design Notes:
!   Fortran 90 doesn't allow one to check the number of bytes in a real variable;
!     we check only that four byte and eight byte reals have different kind parameters. 
!
! !END
module typeSizes
  implicit none
  public
  integer, parameter ::   OneByteIntKind = selected_int_kind(2), &
                          TwoByteIntKind = selected_int_kind(4), &
                         FourByteIntKind = selected_int_kind(9), &
                        EightByteIntKind = selected_int_kind(18)

  integer, parameter ::                                          &
                        FourByteRealKind = selected_real_kind(P =  6, R =  37), &
                       EightByteRealKind = selected_real_kind(P = 15, R = 307)
contains
  logical function byteSizesOK()

! Local variables
    integer (kind =  OneByteIntKind) :: One
    integer (kind =  TwoByteIntKind) :: Two
    integer (kind = FourByteIntKind) :: Four

    if (bit_size( One) == 8  .and. bit_size( Two) == 16 .and.  &
        bit_size(Four) == 32 .and.                             &
        FourByteRealKind > 0 .and. EightByteRealKind > 0 .and. &
        FourByteRealKind /= EightByteRealKind) then
      byteSizesOK = .true.
    else
      byteSizesOK = .false.
    end if
  end function byteSizesOK
end module typeSizes
