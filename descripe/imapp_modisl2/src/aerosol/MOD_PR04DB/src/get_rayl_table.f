      subroutine get_rayl_table(nvalc_in, x1, x2, x3, xq, xu, xi)
c
c			include 'newaottbl.inc'

      real nvalc_in(3,10,46,31), yy(3)
      real theta0(10), theta(46), phi(31)
      character*6 xname412(20)
c     data pi      /3.14159/
      data pi      /3.141592653589793/
      data theta0  /0.0,8.0,16.0,24.0,32.0,40.0,48.0,56.0,64.0,72.0/
			
			
c-----------------------------------------------------------------------
c Load Tables  --> Now, the tables are read from the main program.
c-----------------------------------------------------------------------
!...  open (1,file='./Sub/table_rayl_412nm_Ref_07', status='old')
!...  read (1,*) nval
!...  close (1)
c-----------------------------------------------------------------------
c     End Loading Tables
c-----------------------------------------------------------------------

!...  x1 = 31.35  ! SZA 
!...  x2 = 23.95  ! VZA 
!...  x3 = 177.67 ! RAA 

      do i = 1, 46
      theta(i) = 2.*float(i-1)
      enddo

      do i = 1, 31
      phi(i) = 0. + 6.*float(i-1)
      enddo

      mm = 10     ! solar zenith
      nn = 46     ! satellite zenith
      ll = 31     ! rel. azimuth

c---     interpolating View Geometry tables
c...  Note: This code goes wrong when x3 (RAA) = 180.  
      if(x3.ge.179.9999) x3=179.9999

      do 100 ia = 1, 3 
      call intep(theta0, theta, phi, nvalc_in, mm, nn, ll, ia,
     1          x1,x2,x3,y,dy)
      yy(ia) = y/pi
100   continue

      XI = yy(1)      
      XQ = yy(2)
      XU = yy(3)

c      print *,'I, Q, U =', XI,XQ,XU
!      print 12345,'I, Q, U =', XI,XQ,XU,100.0*XQ/XI, 100.0*XU/XI  
!12345 format(a9,3(2x,f11.8),2(2x,f10.6))

      return
      end
c
c-------------------------------------------------------------
c
      subroutine intep(x1a,x2a,x3a,ya,m,n,l,ia,x1,x2,x3,y,dy)
      parameter (nmax=46,mmax=10,lmax=31)
      dimension x1a(m),x2a(n),x3a(l),ya(3,m,n,l)
      dimension xx2a(4), xx1a(4)
      dimension yntmp(4),ymtmp(4),yltmp(lmax)

      integer i, j, k, o

      call search3(x3,x3a,l,ii)
      frac = (x3-x3a(ii))/(x3a(ii+1)-x3a(ii))

      nsm=1
      dif=x1-x1a(1)
      do i=1,m
        dift=x1-x1a(i)
        if (dift.gt.0. .and. dift.lt.dif) then
          nsm=i
          dif=dift
        endif
      enddo
      mbeg = nsm - 2
      if (mbeg.le.0) then
        mbeg = 0
      else if (mbeg.gt.m-4) then
        mbeg = m-4
      endif

      nsn=1
      dif=x2-x2a(1)
      do i=1,n
        dift=x2-x2a(i)
        if (dift.gt.0. .and. dift.lt.dif) then
          nsn=i
          dif=dift
        endif
      enddo
      nbeg = nsn - 2
      if (nbeg.le.0) then
        nbeg = 0
      else if (nbeg.gt.n-4) then
        nbeg = n-4
      endif

      do 12 j=1,4
        do 11 k=1,4
          do 10 o=1,l
           yltmp(o)=ya(ia,j+mbeg,k+nbeg,o)
10        continue
          yntmp(k) = yltmp(ii)*(1.-frac) + yltmp(ii+1)*frac
          xx2a(k) = x2a(k+nbeg)
11      continue
        call polint2(xx2a,yntmp,4,x2,ymtmp(j),dy)
        xx1a(j) = x1a(j+mbeg)
12    continue
      call polint2(xx1a,ymtmp,4,x1,y,dy)
      return
      end
c
c     --------------
c
      subroutine polint2(xa,ya,n,x,y,dy)
      parameter (nmax=50) 
      dimension xa(n),ya(n),c(nmax),d(nmax)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n 
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
c          if(den.eq.0.) print *,'den=',den
c          if(den.eq.0.)pause
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      end
c
c     --------------
c
      subroutine  search3(xbar,x,n,i)
c
c     purpose
c       locate position in table of point at which interpolation is
c       required
c
c     usage
c       call  search3 (xbar, x, n, i)
c
c     description of parameters
c       xbar   - point at which interpolation is required
c       x      - array containing independent variable
c       n      - number of points in x array
c       i      - index specifying segment containing xbar
c
      dimension x(n)
      data b/.69314718/
      if (n.lt.2) go to 15
      if(x(1).gt.x(2)) go to 17
      m = int((log(float(n)))/b)
      i=2**m
      k=i
   10 k=k/2
      if (xbar.ge.x(i).and.xbar.lt.x(i+1))return
      if (xbar.gt.x(i)) go to 12
      i = i-k
      go to 10
   12 i = i+k
      if (i.le.n) go to 10
      i=n
      go to 10
   15 print *, "Search n is less than 2."
      return
   17 print *, "Search table is not in increasing order."
      return
  100 format(5x,6a8)
      end 
