       subroutine get_irphase(npts,indat_bt,icode)
C!F77-------------------------------------------------------------------
c!Description: This routine is the meat of the tri_spectral cloud phase
c technique.  It performs a cluster analysis on points not labelled
c clear by the cloud mask product, which aids in opaque 
c cloud identification.  It continues by calling the maximum likelihood
c estimator which works well in indentifying single phase cloud.
c Finally, points not identified by the opaque and single phase
c cloud algorithms are phase identified single pixel by single
c pixel.  
c 
c!Input Parameters:
c npts     Number of points in this particular set
c indat_bt Array which holds the brightness temperature information
c          and the spatial variance info. for the current scan
c 
c!Output Parameters: 
c icode    Output cloud phase identifier
c
c!Revision History:
c $Id: get_irphase.f,v 1.8.2.3 2001/05/02 14:04:50 raf Exp $
c
c!Team-unique Header:
c
c!References and Credits:
c
c See Cloud Mask ATBD-MOD-04 and 
c Strabala, K.I., S.A. Ackerman and W.P.Menzel, 1994: Cloud
c properties inferred from 8-12 micron data.  J. App. Meteor,
c Vol. 33, No. 2, 212-229.
c
c!Design Notes:
c!End------------------------------------------------------------------ 
      implicit none
      save
 
      include 'mod06uw_data.inc'
      include 'mod06uw_debug.inc'

c ... scalar arguments
      integer npts

c ... array arguments
      integer icode( max_samp_line*max_samp_pixel ),
     +        indat_bt( max_samp_line*max_samp_pixel , irphase_band+5 )

c ... local scalars
      integer maxc,loop,knt
      real bt8,bt11,bt12

c ... local arrays ..
      real btd8m11
 

c ...  debug statement ............................................
       if (debug .gt. 0) then
         write( h_output, '(72(''-''))' )
         write(h_output,'(2x/,''Subroutine get_irphase  '',/)')
       endif

c ...   Initialize output

        maxc = npts

        ! initialize the icode vector 
        do knt=1,maxc
                icode(knt) = -1
                ! flag boxes as being clear sky
                if (indat_bt(knt,1) .eq. clr_data) then
                  icode(knt) = 0                   ! value of clr_data is -999
                endif
                ! flag scenes where the input data were bad
                if (indat_bt(knt,1) .eq. out_misg) then
                        icode(knt) =  clr_misg     ! value of clr_misg is 100
                endif
        enddo

! ...
! ...   ************* point by point identification *********

        do loop = 1,maxc

! ...     set up brightness temperatures and brightness temperature differences
          bt8  = indat_bt(loop,1)/100.
          bt11 = indat_bt(loop,2)/100.
          bt12 = indat_bt(loop,3)/100.
  
          btd8m11 = bt8 - bt11
  
!************************************
!   Output phase identification codes
!************************************
!   Value  Meaning
!   -----  ------- 
!     0    cloud free 
!     1    water cloud 
!     2    ice cloud
!     3    mixed phase cloud
!     4    also ice cloud (redundant)
!     5    also water cloud (redundant)
!     6    undecided
!     100  means missing data (called clr_misg)

! 23 February -
! attempting to code the FORTRAN equivalent of SN's matlab
! simphase.m code

          if (icode(loop).eq.-1) then
 
              if (bt11 .le. 238.0 .or. btd8m11 .ge. 0.5) then
                  icode(loop) = 2  ! flag as ice cloud

              else if (bt11 .gt. 238.0  .and.
     +                 bt11 .lt. 268.0  .and.
     +  	       btd8m11 .ge. -0.25 .and.
     +                 btd8m11 .lt. 0.5) then
!                 icode(loop) = 3  ! flag as mixed-phase cloud
!
! For Collection 6, eliminate separate "mixed phase" category and 
! include with "undecided" (6).

                  icode(loop) = 6  ! flag as mixed-phase cloud

              else if (bt11 .gt. 238.0  .and.
     +                 bt11 .lt. 268.0  .and.
     + 	      	       btd8m11 .lt. -0.25  .and. 
     +  	       btd8m11 .ge. -1.0) then
                  icode(loop) = 6  ! flag as uncertain phase 

              else if (bt11 .gt. 238.0  .and.
     +  	       btd8m11 .le. -1.0 ) then                  
     		  icode(loop) = 1  ! flag as water cloud
            	            
              else if (bt11 .ge. 285.0   .and.
     +  	       btd8m11 .le. -0.5 ) then                  
     		  icode(loop) = 1  ! flag as water cloud	  

              else 
                  icode(loop) = 6  ! flag what's left as uncertain

              endif
     
      
          endif  ! end classification
  
c ....... debug statement ..........................................
          if (debug .eq. 4) then
            write( h_output, '(72(''-''))' )
            write(h_output,'(2x/,''IR phase output '',/)')
            write(h_output,'(''bt11,bt8,bt8m11,icode: '')')
            write(h_output,'(3f10.3,i10)') bt11,bt8,btd8m11,
     *        icode(loop)
          endif

  	enddo
     
      return
      end
