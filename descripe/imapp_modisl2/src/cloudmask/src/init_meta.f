      subroutine init_meta(East_Lim,West_Lim,North_Lim,South_Lim,
     +                     Begin_G_Date,Begin_G_Time,End_G_Date,
     +                     End_G_Time,TAI_start,TAI_end,nscans,
     +                     max_pixel,lines_in_edge,pixels_in_edge)

      implicit none
      save

c----------------------------------------------------------------------
C!F77
c
c!Description:
C     Routine for initializing variables needed in metadata routines.
C
c!Input parameters:
c
c East_Lim      East longitude limit of given granule
c West_Lim      West longitude limit of given granule
c North_Lim     North latitude limit of given granule
c South_Lim     South latitude limit of given granule
c Begin_G_Date  Beginning date of granule
c Begin_G_Time  Beginning time of granule
c End_G_Date    Ending date of granule
c End_G_Time    Ending time of granule
c TAI_start     Beginning TAI time of granule
c TAI_end       Ending TAI time of granule
c nscans        Number of scans for this granule
C max_pixel     Maximum number of pixels per scan in granule
c lines_in_edge Number of lines outside of processing region
c pixels_in_edge Number of elements outside of processing region
c
c!Output Parameters:
c None
c  
c!Revision History:
c
c Revision 1.1  1997/03/07  20:10:56  wolf
c Initial revision
c
c!Team-unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c
c!END
c----------------------------------------------------------------------
c
      INCLUDE 'global.inc'

c     scalar arguments
      character*70 Begin_G_Time,End_G_Time,Begin_G_Date,End_G_Date
      integer nscans,lines_in_edge,pixels_in_edge,max_pixel
      double precision East_Lim,West_Lim,North_Lim,South_Lim,
     +                 TAI_start,TAI_end

c     Initialize character time/data holders
      Begin_G_Time = '  '
      End_G_Time = '  '
      Begin_G_Date = '  '
      End_G_Date = '  '

      max_pixel = 0
      nscans = 0
      lines_in_edge = (nlcntx-1)/2
      pixels_in_edge = (necntx-1)/2

c     Initialize granule lat/lon bound variables
      East_Lim = 0.0d0
      West_Lim = 0.0d0
      North_Lim = 0.0d0
      South_Lim = 0.0d0
      TAI_start = 0.0d0
      TAI_end = 0.0d0

      return
      end



