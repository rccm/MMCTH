      subroutine height(p,t,w,zs,nl, z)
c * Compute geopotential height ... profiles stored from top down
c.. hong.. Feb.27,2002   heights in meters*0.1
c   Rich Frey Sept. 25, 2008  heights in km

      real rog,fac
      parameter (rog=29.2898)
      parameter (fac=rog*0.5)
      real p(*),t(*),w(*),z(*)
      real vl,algpl,h,zs,vu,algpu
      integer  l,nl,m
c	  Input:
c		p ... profile of pressure (mb)
c		t ... profile of temperature (degK)
c		w ... profile of h2o mixing ratio (g/kg)
c		zs .. surface height -- input 0. if not known
c		nl .. number of levels
c	   Output:
c		z ... profile of geopotential height (m)

      save

      vl=t(nl)*(1.+.00061*w(nl))
      algpl=alog(p(nl))
c * start with surface height
	h=zs
c * do integration
      do l=2,nl
         m=nl+1-l
         vu=t(m)
c ..... virtual temperature adjustment
         if(p(m).ge.300.) vu=vu*(1.+.00061*w(m))
         algpu=alog(p(m))
         h=h+fac*(vu+vl)*(algpl-algpu)
         algpl=algpu
         vl=vu
c **** store heights like temp from top down to surface (level 'nl')
c         z(m)=h
c ... R. Frey:    change to km
         z(m)=h * 0.001
	enddo
      return
      end
