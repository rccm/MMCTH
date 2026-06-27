      subroutine get_sg_thresholds(refang,locut,hicut,midpt,power)

      implicit none
      save

c---------------------------------------------------------------------
C!F77
c
c!Description:
c     Routine for setting cloud test thresholds for band 2 in sun-glint
c     conditions.
c
c!Input parameters:
c refang        Reflectance angle
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
      include 'snglntr_thr.inc'

c     Scalar arguments.
      real refang,locut,hicut,midpt,power

c     Local scalars.
      integer debug, h_output
      real lo_ang,hi_ang,lo_ang_val,hi_ang_val,a,conf_range

c     Common statement for debug purposes
      common / bug / debug, h_output

c---------------------------------------------------------------------

c     Compute clear sky confidence interval boundaries and mid-point.

      if(refang .le. snglnt_bounds(2)) then

        hicut = snglnt0(3)
        locut = snglnt0(1)
        midpt = snglnt0(2)
        power = snglnt0(4)

      else

        if(refang .gt. snglnt_bounds(2) .and. refang .le. 
     +                 snglnt_bounds(3)) then

          lo_ang = snglnt_bounds(2)
          hi_ang = snglnt_bounds(3)
          lo_ang_val = snglnt10(1)
          hi_ang_val = snglnt10(2)
          power = snglnt10(4)
          conf_range = snglnt10(3)

        else if(refang .gt. snglnt_bounds(3) .and. refang .le. 
     +                      snglnt_bounds(4)) then

          lo_ang = snglnt_bounds(3)
          hi_ang = snglnt_bounds(4)
          lo_ang_val = snglnt20(1)
          hi_ang_val = snglnt20(2)
          power = snglnt20(4)
          conf_range = snglnt20(3)

        end if

        a = (refang - lo_ang) / (hi_ang - lo_ang)
        midpt = lo_ang_val + a*(hi_ang_val - lo_ang_val)
        hicut = midpt - conf_range
        locut = midpt + conf_range
       
      end if

c---------------------------------------------------------------------

c     Debug statement
      if(debug .gt. 1) then
        write(h_output,'(5f10.3)') snglnt_bounds,refang
        write(h_output,'(5f10.3)') a,midpt,hicut,locut,power
      end if

c---------------------------------------------------------------------

      return
      end
