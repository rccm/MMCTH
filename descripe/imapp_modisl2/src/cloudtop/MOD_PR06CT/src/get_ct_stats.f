      subroutine get_ct_stats( line, pixel, qual_flag, co2_flag,
     &  bias_flag, pct, eca, ssctpr,
     &  slcd, smcd, shcd, sthncd, sthkcd, sopcd, scicd,sco2,sbias )
      
c-----------------------------------------------------------------------
c!F77
c
c!Description:
c This routine sums the number of occurrances of particular cloud types
c determined by the CO2-slicing method over an entire granule of MODIS
c data.
c
c!Input Parameters:
c line                   Line of current scan (1-10)
c pixel                  Pixel in current scan (1-1500)
c qual_flag              Indicates quality of CO2-slicing cloud height
c                          (0=not usable, 1=usable)
c pct                    cloud top pressure (mb)
c eca                    cloud top efective emissivity (0-1.0)
c
c!Output Parameters:
c ssctpr                 number of usable cloud top pressure retrievals
c slcd                   number of low cloud heights found 
c smcd                   number of mid-level cloud heights found
c shcd                   number of high cloud heights found
c sthncd                 number of "thin" clouds found
c sthkcd                 number of "thick" clouds found
c sopcd                  number of opaque clouds found
c scicd                  number of cirrus clouds found
c
c!Revision History:
c $Id: get_ct_stats.f,v 1.1.1.1 2005/02/22 17:15:54 gumley Exp $
c
c!Team-unique Header:
c
c!End
c-----------------------------------------------------------------------

      implicit none
      save 

      include 'mod06uw_debug.inc'

c ... scalar arguments
      integer qual_flag,ssctpr,slcd,smcd,shcd,sthncd,sthkcd,sopcd,scicd,
     &  line,pixel, co2_flag,sco2,bias_flag,sbias
      real pct,eca

c-----------------------------------------------------------------------

c ... Get sums of the various cloud categories

      if(qual_flag .eq. 1) ssctpr = ssctpr + 1
      if(co2_flag .eq. 1)  sco2 = sco2 + 1
      if(bias_flag .eq. 1) sbias = sbias + 1
      if(pct .gt. 680.0) slcd = slcd + 1
      if(pct .le. 680.0 .and. pct .ge. 440.0) smcd = smcd + 1
      if(pct .lt. 440.0) shcd = shcd + 1
      if(eca .le. 0.50) sthncd = sthncd + 1
      if(eca .gt. 0.50 .and. eca .le. 0.95) sthkcd = sthkcd + 1
      if(eca .gt. 0.95) sopcd = sopcd + 1
      if(eca .le. 0.95) scicd = scicd + 1

c-----------------------------------------------------------------------

c ... Write debug information
      if(debug .gt. 0) then
       if(line .eq. 6 .and. pixel .eq. 676) then
        write(h_output,'(1x)')
        write(h_output,'(''Subroutine get_ct_stats'')')
        write(h_output,'(''qual_flag, ssctpr, slcd, '',
     &                 ''smcd, shcd, sthncd, sthkcd, sopcd, scicd'')')
        write(h_output,'(9i7)') qual_flag,ssctpr,slcd,
     &      smcd,shcd,sthncd,sthkcd,sopcd,scicd
       end if
      end if

c-----------------------------------------------------------------------

      return
      end
