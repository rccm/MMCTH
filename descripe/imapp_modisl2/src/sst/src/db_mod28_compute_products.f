C--------------------------------------------------------------------
C  Copyright (C) 2002,  Space Science and Engineering Center, 
C  University C  of Wisconsin-Madison, Madison WI.
C      
C  This program is free software; you can redistribute it 
C  and/or modify it under the terms of the GNU General 
C  Public License as published by the Free Software Foundation; 
C  either version 2 of the License, or (at your option) any 
C  later version.
C
C  This program is distributed in the hope that it will be 
C  useful, but WITHOUT ANY WARRANTY; without even the implied 
C  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
C  See the  GNU General Public License for more details.
C
C  You should have received a copy of the GNU General Public 
C  License along with this program; if not, write to the Free 
C  Software Foundation, Inc., 59 Temple Place, Suite 330, 
C  Boston, MA  02111-1307 USA
C--------------------------------------------------------------------
C
C
      subroutine db_mod28_compute_products(line,pixel,qual_flag)

c-----------------------------------------------------------------------
c
c!F77
c
c!Description:
c    Compute MOD28 retrieval products.
c
c!Input Parameters:
c    cfgname          Configuration file name
c    LINE             Line number within a swath (1-10)
c    PIXEL            Pixel number within a 1km scan (1-1500)
c
c    The following in COMMON /MOD28_DATA/ are used:
c    RADIANCE1        Radiances for IR bands
c
c!Output Parameters:
c    QUAL_FLAG        SST quality flag 
c                       (0=unusable,1=usable)
c    CONF_FLAG        SST confidence flag 
c                       (0-3, low to high)
c
c    The following are output to COMMON /MOD28_DATA/:
c    PROCESSING_FLAG  Flag indicating which technique was used to 
c                      retrieve SST
c    SST1             the 11/12 micron SST
c    SST4             the 3.9/4.1 micron SST
c
c!Revision History:
c
c calls:
c function modis_bright    Compute brightness temperature from radiances
c
c!Team-unique Header:
c
c!End
c
c-----------------------------------------------------------------------
      
      implicit none
      save
      
      include 'db_mod28uw_data.inc'
      include 'db_mod28uw_debug.inc'
      include 'platform_name.inc'
      
c ... arguments

      integer line, pixel, qual_flag

c ... local variables

      integer k, proc_flag, units, isat, ip, il, in, i, j
      real rad,temp,d2r,mu
      real b20,b22,b23,b31,b31_b32,bsst
      real sstday,sstnight,w1,w2,reyn_sst,tenv 
      real k1(3,2),k2(3,2),k3(3,2),k4(3,2)
      logical day

c ... external functions

      real modis_bright
      external modis_bright

c ... data (move to file later)
      data k1/1.052,1.886,-0.065,1.152,2.133,0.987/
      data k2/0.984,0.938,1.034,0.960,0.926,1.031/
      data k3/0.130,0.128,0.723,0.151,0.125,0.349/
      data k4/1.860,1.094,0.972,2.021,1.198,1.766/
      data d2r/0.0174533/


      qual_flag = 0

c ... Determine the coefficient set (Terra or Aqua)
      if (platform_name(1:5) .eq. 'Terra' .or.
     &    platform_name(1:5) .eq. 'terra' .or.
     &    platform_name(1:5) .eq. 'TERRA') then
        isat = 1
      else if (platform_name(1:4) .eq. 'Aqua' .or.
     &         platform_name(1:4) .eq. 'aqua' .or.
     &         platform_name(1:4) .eq. 'AQUA') then
        isat = 2
      else
        isat = 0
        call message('db_mod28_compute_products.f',
     &    'Platform name not recognized ' ,
     &     0, 2)
      endif

c     indices (1,2,3,4,5) = (20,21,23,31,32) bands
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      mu = cos(VIEW(pixel,line)*d2r)

c     convert radiance to brightness temp
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~          
      units = 1
      do k = 1 , n_output
        rad = radiance1( pixel, line, sstbandptr(k))
        raw_radiance(pixel,line,k) = rad 
	if (rad.lt.0.0) then 
	  qual_flag = qual_flag - 1
          temp =  bad_value
	else
          temp = modis_bright(rad,sstband(k),units)
	end if
        brightness_temp(pixel,line,k) = temp
      end do

c     first attempt algorithm - no bias correction
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (isat.eq.0) goto 999 
      b20 = brightness_temp(pixel,line,1) - 273.15
      b22 = brightness_temp(pixel,line,2) - 273.15
      b23 = brightness_temp(pixel,line,3) - 273.15
      b31 = brightness_temp(pixel,line,4) - 273.15
      b31_b32 = 0.0
      in = 0
      do ip = MAX0(pixel-1,1),MIN0(pixel+1,MAX_PIXEL)
        do il = MAX0(line-1,1),MIN0(line+1,MAX_LINE)
          rad = radiance1(ip,il,sstbandptr(4))
	  if (rad.gt.0.0) then
            b31_b32 = b31_b32 + modis_bright(rad,sstband(4),units)
            rad = radiance1(ip,il,sstbandptr(5))
            if (rad.gt.0.0) then
              b31_b32 = b31_b32 - modis_bright(rad,sstband(5),units)
              in = in + 1
	    else
	      b31_b32 = b31_b32 - modis_bright(rad,sstband(4),units)
	    end if
	  end if
        end do
      end do
      if (in.ge.1) then
        b31_b32 = b31_b32/FLOAT(in)
      else
        qual_flag = qual_flag - 1
      end if

      if (qual_flag.lt.0) goto 900

c ... Get day night flag from solar zenith angle
      day = .false.
      if (nint(solz(pixel,line)) .ne. nint(bad_value)) then
          if (solz(pixel,line).le.90.0 .and. solz(pixel,line).gt.0.0) then
              day = .true.
          endif
      endif

c     (1)  SST4 algorithm (3.9/4.1 um)
c     ................................

      sst4(pixel,line) = k1(3,isat) + k2(3,isat)*b22 
     &  + k3(3,isat)*(b22-b23) + k4(3,isat)*(1./mu-1.)

c     (2)  SST algorithm (11/12 um)
c     .............................

      i = NINT(LON1(pixel,line) + 180.5)
      j = NINT(LAT1(pixel,line) + 90.5)
      if ((i.ge.1).and.(i.le.360)) then
        if ((j.ge.1).and.(j.le.180)) then
          if (oisst(i,j).gt.0.0) then
            reyn_sst = oisst(i,j)
          end if    
        end if    
      end if

c ... Added logic to use oisst in sun glint regions
      if (day) then
          bsst = reyn_sst
          tenv = reyn_sst
      else
          bsst = b20
          if (reyn_sst .gt. 0.0) bsst = reyn_sst
          tenv = sst4(pixel,line)
      endif

      sstday = k1(1,isat) + k2(1,isat)*b31 + k3(1,isat)*b31_b32*bsst 
     &  + k4(1,isat)*b31_b32*(1./mu-1.)
      sstnight = k1(2,isat) + k2(2,isat)*b31 
     &  + k3(2,isat)*b31_b32*tenv
     &  + k4(2,isat)*b31_b32*(1./mu-1.)

      if (b31_b32 .lt. 0.5) then
        sst1(pixel,line) = sstday
      else if (b31_b32 .gt. 0.9) then
        sst1(pixel,line) = sstnight
      else
        w1 = (0.9-b31_b32)/0.4
        w2 = 1.- w1
        sst1(pixel,line) = w1*sstday + w2*sstnight
      end if      
 
  900 proc_flag = mod28_proc_opt

c-----------------------------------------------------------------------

c ...   Write debug information about input data.
  999 if(debug .gt. 0) then
        if(line .eq. 6 .and. pixel .eq. 676) then
          write(h_output,'(1x)')
          write(h_output,'(''Routine mod28_compute_products'')')
          write(h_output,'(''Line, pixel are '',2i5)') line,pixel
          write(h_output,'(''Process option is '',i3)')
     &       mod28_proc_opt
          write(h_output,'(''Input radiances (radiance1)'')')
          do k = 1,n_output
            write(h_output,'(''Ch '',i5,'' Pos '',i5)')
     &        sstband(1),sstbandptr(1)
            write(h_output,'(1x,5f10.3)')
     &        radiance1(pixel,line,sstbandptr(1))
          end do
          write( h_output,
     &    ' (''sst1 4='',2f9.3,'', bsst ='',f10.3,'', tenv ='',f10.3)')
     &    sst1( pixel, line ), sst4(pixel,line), bsst, tenv

        end if
      end if

c-----------------------------------------------------------------------

      return
      end
