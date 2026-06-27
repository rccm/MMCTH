      subroutine store_meta(mlin,scan_type,scan_pixels,nadir,
     +                      s_type,s_pixels,n_nadir)

      implicit none
      save

c ... Common statement for debug purposes
      common / bug / debug, h_output

c---------------------------------------------------------------------
c!F77 
c
c!Description:
c     Routine which takes scan cube's worth of a certain variable
c     and puts it into a context of lines (usually 3 lines worth).
c     Use routine bright to convert radiances to brightness
c     temperatures.
C
c!Input Parameters:
c mlin         Context line number
c scan_pixels  Number of pixels in this scan
c scan_type    Scan type - Day or night
c nadir        Nadir pixel number
c 
c!Output Parameters:
c s_pixel      Number of pixels in this scan for lines in a context
c s_type       Scan type - Day or night for lines in a context
c n_nadir      Nadir pixel number for lines in a context
c
c!Revision History:
c
c!Team-Unique Header:
c
c    This software is developed by the MODIS Science Data Support Team
c    for the National Aeronautics and Space Administration,
c    Goddard Space Flight Center, under contract NAS5-32373.
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c
c!END
c---------------------------------------------------------------------

      include 'global.inc'

c     scalar arguments
      character*1 scan_type
      integer nadir,scan_pixels,mlin

c     array arguments
      character*1 s_type(nlcntx)
      integer n_nadir(nlcntx),s_pixels(nlcntx)

c     local scalars
      integer i,h_output,debug

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(10x/,''Subroutine store_meta  '',/)')
      endif
c ................................................................

c     Store nlcntx number of metadata values in holders.  This will
c     make sure you don't screw up your processing by crossing a 
c     scan
      do 190 i  = 1 , nlcntx

            s_type(mlin) = scan_type
            s_pixels(mlin) = scan_pixels
            n_nadir(mlin) = nadir

 190  continue
  

c ... debug statement ............................................
      if (debug .gt. 3) then
        write(h_output,'(15x,'' META FOR 3 CUBES '',I4)') mlin
        write(h_output,'(2x,'' s_type '',(3(7x,A1)/))') 
     +       s_type(1),s_type(2),s_type(3)
        write(h_output,'(2x,'' s_pixels '',(3I8/))') 
     +       s_pixels(1),s_pixels(2),s_pixels(3)
        write(h_output,'(2x,'' n_nadir '',(3I8/))') 
     +       n_nadir(1),n_nadir(2),n_nadir(3)
      endif
c ................................................................

      return 
      end
