       subroutine put_irphase_out_arrays(npixels, icode, indat_bt, var)

C!F77------------------------------------------------------------------
c!Description: Place the processing variables into the
c    output variables.  This requires going from one and two
c    dimension processing arrays into a two or 3 dimensioned
c    output arrays.
c
c!Input Parameters:
c     npixels                    Number of pixels in scan
c
c           Processing arrays:
c    icode                       holds 1D array of phase product
c    indat_bt                    holds bt values for processing
c    var                         holds variance values for processing
c
c!Output Parameters:   None
c
c!Revision History:
c $Id: put_irphase_out_arrays.f,v 1.5 1999/04/16 22:58:02 kis Exp $
c
c!Team-unique Header:
c
c!References and Credits:
c   K. Strabala (kathys@ssec.wisc.edu)
c   See Cloud Mask ATBD-MOD-04 and HDF specs.
c
c!DESIGN NOTES:
c
c    The following arrays in COMMON /MOD06_DATA/ are filled:
c    radiance_variance           Variance of radiances over box
c                                 for all three ir phase bands
c    cloud_phase                 Tri-spectral irphase result

c!END-------------------------------------------------------------------
      implicit none
      save

      include 'mod06uw_data.inc'
      include 'mod06uw_debug.inc'
 
c ... scalar arguments
      integer npixels
c     array arguments
      integer indat_bt( max_samp_line*max_samp_pixel , irphase_band+5 ),
     +        icode( max_samp_line*max_samp_pixel )
      real var ( max_samp_line*max_samp_pixel , irphase_band )


c ... local variables
      integer line, pixel, kount, kount1, kount2, bnd

c ... initalize
      kount = 0
      kount1 = 0
      kount2 = 0

c ... first place actual phase code result and ecosystme into the 
c ...  output product array
c ...
      do line = 1, ( max_line / isamp ) 
         do pixel = 1, ( npixels / isamp ) 
            kount = kount + 1
            if (icode(kount) .ne. 100) 
     +      cloud_phase( pixel, line ) = real(icode(kount))
         end do 
       end do

c ... Next the radiance variance 
      do line = 1, ( max_line / isamp ) 
         do pixel = 1, ( npixels / isamp ) 
            kount1 = kount1 + 1 
            do bnd = 1 , irphase_band
               radiance_variance(pixel,line,bnd) =  var(kount1,bnd)
            end do
         end do
      end do

      return
      end
