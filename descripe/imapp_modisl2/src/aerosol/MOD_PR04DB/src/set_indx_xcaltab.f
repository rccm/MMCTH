      subroutine set_indx_xcaltab(sectab, sec_inp, it1, it2)

! subroutine to interpolate m12, m13 coeff. LUT (from the SeaWiFS team)
! Written by MJ, Aug, 2008
!c... modified Oct 2010 to change array size (Terra/MODIS)

      implicit none

      real*8 sec_inp
!     real*8 sectab(200)
      real*8 sectab(127) ! added 20101012
      integer it1, it2, itt

!     Initialization
      it1=0
      it2=0

      if(sec_inp.le.sectab(1)) then
         it1=1
         it2=2
!c    elseif (sec_inp.ge.sectab(92)) then
      elseif (sec_inp.ge.sectab(127)) then
!c       it1=91
!c       it2=92
         it1=126
         it2=127
!...     note: this situation could be dangerous if applied.
      else
!c       do itt=2,91
         do itt=2,127
            if(sec_inp.ge.sectab(itt-1).and.sec_inp.lt.sectab(itt)) then
               it1=itt-1
               it2=itt
!              print *, itt-1, 'th ',sec_inp, sectab(itt-1), sectab(itt)
            endif
         enddo
      endif

      return
      end

