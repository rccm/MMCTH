      subroutine get_nl_thresholds(btdiff,locut,hicut,midpt,power)

      implicit none
      save

c---------------------------------------------------------------------
C!F77
c
c!Description:
c     Routine for setting cloud test thresholds for 11 minus 3.9 um
c     test in polar night conditions.
c
c!Input parameters:
c btdiff        11-12 um BTD for current pixel
c
c!Output Parameters:
c locut         Zero confidence value of clear sky confidence interval
c hicut         1.0 confidence value of clear sky confidence interval
c midpt         0.5 confidence value of clear sky confidence interval
c power         Power of clear sky confidence curve
c
c!Revision History:
c Original version - 03/07/2002     
c Rich Frey              MODIS Science Group at UW-Madison
c
c!Team-unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-35.
c
c!END
c---------------------------------------------------------------------

c     Include files.
      include 'LandNite_thr.inc'

c     Scalar arguments.
      real btdiff,locut,hicut,midpt,power

c     Local scalars.
      real lo_val,hi_val,lo_val_thr,hi_val_thr,a,conf_range

c     Compute clear sky confidence interval boundaries and mid-point.

      if(btdiff .gt. bt_diff_bounds(1)) then

        hicut = nl_11_4l(3)
        locut = nl_11_4l(1)
        midpt = nl_11_4l(2)
        power = nl_11_4l(4)

      else if(btdiff .lt. bt_diff_bounds(2)) then

        hicut = nl_11_4h(3)
        locut = nl_11_4h(1)
        midpt = nl_11_4h(2)
        power = nl_11_4h(4)

      else

        lo_val = bt_diff_bounds(1)
        hi_val = bt_diff_bounds(2)
        lo_val_thr = nl_11_4m(1)
        hi_val_thr = nl_11_4m(2)
        conf_range = nl_11_4m(3)
        power = nl_11_4m(4)

        a = (btdiff - lo_val) / (hi_val - lo_val)
        midpt = lo_val_thr + a*(hi_val_thr - lo_val_thr)
        hicut = midpt - conf_range
        locut = midpt + conf_range
       
      end if

c     write(*,'(10f7.2)') btdiff,bt_diff_bounds,lo_val_thr,hi_val_thr,
c    *     conf_range,a,midpt,hicut,locut

      return
      end
