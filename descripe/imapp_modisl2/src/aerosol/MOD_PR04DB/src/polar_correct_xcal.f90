! ----------------------------------------------------------------------
! modis_polcor() - Polarization correction for MODIS.
!
!
! Algorithm by:
!
!	Howard Gordon and Menghua Wang
!
! Modified for MSl12 by:
!
!       B. Franz, SAIC, January 2003.
!
! Notes:
!
!   This code reads and applies the standard MODIS polarization correction
!   tables developed by Howard Gordon, Univ. of Miami.  The routine was 
!   originally adapted from new_modis_pol_corr_sub.f (revision 1.12) of the 
!   modcol software developed at RSMAS, Univ. of Miami.
!
!   Collection #4 MODIS processing assumes that only the Rayleigh part of fred@seawifs.gsfc.nasa.gov
!   the total radiance is polarized, i.e., the aerosol and the water parts 
!   have no polarization associated with them.  This is reproduced by 
!   setting L_x to the TOA radiances (whitecap and ozone corrected). 
!
!   If L_x is set to Lr_i, this would give the polarization correction 
!   assuming all components are polarized, but in the same way as the Rayleigh.
!
! History:
!
!   Changed take polarization file as input arg, easing redirection to 
!   alternate files.  Eliminated Miami option to switch files between north
!   and south. BAF, 16 March 2004.
!
!
! Modified by JCW
! version 1: keep same way as original file
!
! Modified by Myeong-Jae Jeong (MJ), August 2008
!   - RVS and polarization sensitivity derived by the SeaWiFS team 
!    (Ewa Kwiatkowska, Bryan Franz, ..) are used here instead of pre-launch
!     polarization correction parameters. The name of subroutine was changed
!     to clarify this change.
!   - The definition of polarization correction factor conforms the one
!     used in Meister et al.(2005; Appl. Optics, Vol44, No26) and Gordon
!     et al.(1997; Appl. Optics, 1997, Vol36, No27). This definition and
!     the definition of "alpha" are consistent with those used in Seadas5.2.
!   - The definition of "polcor" was changed into the inverse of the old
!     definition to be consistent with the one used in seadas 5.2; i.e.,
!     "polcor"_new = 1.0 / "polcor"_old
! ----------------------------------------------------------------------

!... L_x:  TOA radiance corrected for ozone absorption (solar and view
!...       paths) -- Bryan Franz
!...       L_x is measured MODIS reflectance (after vicarious
!...       calibration applied) -- Gerhard Meister
!... ********************************************
!... Final reflectance (L_unpolarized) to be used for Deep Blue
!... aerosol retrieval can be calculated as follows:
!... "L_unpolarized = L_measured / RVS / polcor",
!... where L_measured is Terra/MODIS measured reflectance,
!... RVS is vicarious calibration coeff.(M11), and
!... polcor is polarization correction factor.
!... ********************************************

	subroutine polar_correct_xcal( lat, phi0, xam12, xam13, xrvs,  &
                                  alpha, Lr_q, Lr_u, L_x, polcor )

	use GeneralAuxType
	use lut_arrays

        implicit none

!       inputs
        real(single), intent(in)  ::  	lat	               ! latitude
        real(single), intent(in)  ::  	phi0	               ! solar azimuth
        real(double), intent(in)  ::  	alpha	               ! polarization frame rotation angle
        real(single), dimension(nbands),intent(in)  ::  Lr_q   ! Rayleigh radiance, q-component 
        real(single), dimension(nbands),intent(in)  ::  Lr_u   ! Rayleigh radiance, u-component 
	real(single), dimension(nbands),intent(in)  ::  L_x    ! normalization radiance
        real(double), dimension(nbands),intent(in)  :: 	xam12, xam13, xrvs  ! coeff m12, m13, and M11 @chg@

!       ouputs
        real(single), dimension(nbands),intent(out) ::  polcor ! polarization correction

!       working vars
	integer	:: ib, iib, iblut
        real(single) ::  Lr_qp, Lr_up

!	
!	Compute the polarization correction for each band.
!	 
 	 do iib = 1, nbands
!	    ib=1     ! Band 08(412nm) 
!	    ib=2     ! Band 03(466nm) 
!	    ib=3     ! Band 01(650nm) 
            ib=iib
	    if (L_x(ib) .gt. 0.0) then
               Lr_qp = Lr_q(ib)*cos(2*alpha*d2r) + Lr_u(ib)*sin(2*alpha*d2r)  
               Lr_up = Lr_u(ib)*cos(2*alpha*d2r) - Lr_q(ib)*sin(2*alpha*d2r)  
!               if (ib == 1) print *, "Lr_qp: ", xam12(ib)*Lr_qp/L_x(ib)
               polcor(ib) = 1.0 / (1.0 - xam12(ib)*Lr_qp/L_x(ib) & 
                                       - xam13(ib)*Lr_up/L_x(ib))  
!... ***********************************************************
!... ***********************************************************
!... ***********************************************************
!...  Polarization correction applies to B08(412nm) & B03*(466nm)
!...   --> *parameters for B03 are set to those for B10 (for sensitivity test)
               if(ib.ge.3) polcor(ib) = 1.0    ! @chg@
!... ***********************************************************
!... ***********************************************************
!... ***********************************************************
	    else
               polcor(ib) = 1.0
	    endif
         enddo   

	end subroutine polar_correct_xcal    

