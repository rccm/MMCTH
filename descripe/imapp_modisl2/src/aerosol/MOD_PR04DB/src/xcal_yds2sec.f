!===========================================================
      subroutine xcal_yds2sec(iyear,ijday,xhr,xcalsec)
!     subroutine xcal_yds2sec(itid,iyear,ijday,xhr,xcalsec)
!-----------------------------------------------------------
!  Purpose: To convert year and DoY to seconds since
!           1993/01/01 00:00:00
!  Inputs: 
!     iyear: Year from xcal data file (xcal_modist.hdf) 
!     ijday: Day of year from xcal data file 
!     xhr  : hour (min, sec in decimal)
! 
!  Output:
!     xcalsec: seconds since 1993/01/01 
! 
!  ver. 1.0  Written by MJ (Aug 2008) 
! 
      parameter (iyrbeg=1993)  ! Starting year Do not change(for MODIS)
      integer iyear,ijday,iyr,naccdays
!     integer iyear,ijday,iyr,naccdays,itid
      integer idleap, onedleap, ndays
      real*8  xcalsec
      real    xhr
!
      xcalsec=0.0D0
      naccdays=0
      if(iyear.ge.1993) then
         do iyr=iyrbeg,iyear
            idleap = mod(iyr,4)
            if(idleap.eq.0) then 
               onedleap=1
            else
               onedleap=0
            endif 

            ndays = 365 + onedleap
            naccdays = naccdays + ndays
         enddo
         naccdays = naccdays - ndays + ijday - 1
!  
!...     xcalsec = (86400.0*naccdays + 3600.0*xhr)*1.0D0
         xcalsec = 86400.0D0*naccdays + 3600.0D0*xhr
!
      else
         print *, 'Error in xcal_yds2sec.f '
         print *, iyear, '=iyear < 1993 - stop!'
         stop
      endif
!
!-----------------------------------------------------------
      return
      end 
!===========================================================

