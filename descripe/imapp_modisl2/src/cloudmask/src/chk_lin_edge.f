       subroutine chk_lin_edge(klin,nlin,lines_in_edge,line_edge)

      implicit none
      save

c--------------------------------------------------------------------
c!F77
c
c!Description:
c     Routine for determining if the current line you are processing
c     is a border line.
c
c!Input Parameters:
c klin          Counts number of lines output to bit file
c nlin          Total number of lines to process
c lines_in_edge Number of lines outside of processing region
c
c!Output Parameters:
c line_edge     Logical variable - true if processing border line
c
c!Revision History:
c
c!Team-Unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c
c!END
c--------------------------------------------------------------------

c     scalar arguments
      integer klin,lines_in_edge,nlin
      logical line_edge

c     Initialize line_edge to false
      line_edge = .false.

c     Compare current line to border values
      if ((klin .lt. lines_in_edge) .or. 
     +   (nlin .le. 2*lines_in_edge) .or. 
     +   (klin .ge. nlin-lines_in_edge)) then
          line_edge = .true.
      end if
  
      return
      end
