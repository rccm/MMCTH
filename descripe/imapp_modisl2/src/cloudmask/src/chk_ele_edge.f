      subroutine chk_ele_edge(nc,s_pixels,pixels_in_edge,klin,nlin,
     +                        lines_in_edge,ele_edge)

      implicit none
      save

c---------------------------------------------------------------------
C!F77
c
c!Description:
C     Routine for determining if the current element you are processing
c     is a border pixel.
C
c!Input parameters:
c nc            Counter for current pixel being processed
c s_pixels      Context of pixels in current scan
c pixels_in_edge Number of elements outside of processing region
c klin          Counts number of lines output to bit file
c nlin          Total number of lines to process
c lines_in_edge Number of lines outside of processing region
c
c!Output Parameters:
c ele_edge     Logical variable - true if processing border pixel
c
c!Revision History:
c
c!Team-Unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c
c!END
c---------------------------------------------------------------------
      include 'global.inc'

c     scalar arguments
      integer nc,pixels_in_edge,nlin,klin,lines_in_edge,i
      logical ele_edge

c     array arguments
      integer s_pixels(nlcntx)

c     Initialize ele_edge to false
      ele_edge = .false.

c     always take the middle line
      if (klin .le. lines_in_edge) then
         i = klin
      else if (klin .gt. nlin-lines_in_edge) then
         i = nlcntx
      else
         i = ((nlcntx-1) / 2) + 1
      endif

c     Compare current elements to border pixels
      if ((nc .le. pixels_in_edge) .or. 
     +   (s_pixels(i) .le. 2*pixels_in_edge) .or. 
     +   (nc .gt. s_pixels(i)-pixels_in_edge)) then
          ele_edge = .true.
      end if
  
      return
      end
