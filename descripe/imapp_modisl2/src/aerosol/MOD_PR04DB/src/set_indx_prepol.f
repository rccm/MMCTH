!c...  -------------------------------------------------------------------
!c...  -------------------------------------------------------------------
      subroutine set_indx_prepol(aoiref,zaoi,iaoi1,iaoi2)
!c... 
!c... written by Myeong-Jae Jeong (MJ) GEST/UMBC & NASA GSFC
!c...    October 2010  
!c...    Finds indices for pre-launch polarization parameters
!c...    to be interpolated. 
!c...   
!c...    input: 
!c...        zaoi   : angle of incidence (AOI)
!c...    output:
!c...        aoiref : AOI node points in the LUT
!c...        iaoi1  : 1st index
!c...        iaoi2  : 2nd index
!c...   
      real*4 aoiref(5)
      real*4 zaoi
      integer iaoi, iaoi1, iaoi2

!    data aoiref/15.500,26.750,38.000,49.250,60.500/
      aoiref(1)=15.500
      aoiref(2)=26.750
      aoiref(3)=38.000
      aoiref(4)=49.250
      aoiref(5)=60.500
!c...  
!c... Initialization
      iaoi1=0
      iaoi2=0
!c...
      if(zaoi.lt.aoiref(1)) then 
         iaoi1=1
         iaoi2=2
      elseif(zaoi.ge.aoiref(1).and.zaoi.lt.aoiref(5)) then
         do iaoi=1,4
            if(zaoi.ge.aoiref(iaoi).and.zaoi.lt.aoiref(iaoi+1)) then
              iaoi1=iaoi
              iaoi2=iaoi1+1
            endif
         enddo
      else
         iaoi1=4
         iaoi2=5
      endif
!c... 
!c...  -------------------------------------------------------------------
      return
      end
!c...  -------------------------------------------------------------------
!c...  -------------------------------------------------------------------


