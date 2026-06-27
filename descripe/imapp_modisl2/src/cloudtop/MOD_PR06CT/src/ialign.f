      subroutine ialign(b,a,ni,m)

      implicit none
      save

c-----------------------------------------------------------------------
c!F77
c
c!Description:
c  This routine orders the data in the array 'b' from smallest to
c     largest.  The re-ordered data is returned in the array 'a'.
c     Upon completion of the ordering, the array 'ni' contains the
c     original positional subscripts of 'b' in the new order.  This
c     index array can then be used to order other data associated
c     with the data in 'b'.
c     'm' is the number of elements to be ordered.
c
c!Input Parameters:
c b     integer array with data set you wish to order
c
c!Output Parameters:
c a     integer array with data set after you have ordered it
c ni    original positional subscripts of b 
c m     number of elements to be ordered
c
c!Revision History:
c $Id: ialign.f,v 1.5 1999/04/16 22:45:06 kis Exp $
c
c!Team-unique Header:
c
c!End
c-----------------------------------------------------------------------
 
c ... scalar arguments
      integer m

c ... array arguments
      integer a(*),b(*),ni(*)

c ... local scalars
      double precision d
      real f,fn
      integer i,ns,s,mm,l,n1,n2,k,n,j

c ... intrinsic functions
      intrinsic int

c ... initialization
      do 100 i=1,m
        ni(i)=i
        a(i)=b(i)
 100  continue

c ... order a(1) and a(2)
      if(a(2) .le. a(1)) then 
        ns=ni(1)
        ni(1)=ni(2)
        ni(2)=ns
        s=a(1)
        a(1)=a(2)
        a(2)=s
      endif

      if (m .gt. 2) then
        mm=m-1
c ...   bisect a(1)---a(n) repeatedly ....
        l=0
        n1=1
        n2 = 0
        do 300 while (n2 .lt. mm)
          l=l+1
          n1=n1*2
          n2=n1*2-1
          if(n2.gt.mm) n2=mm
          do 260 n=n1,n2
            ns=ni(n+1)
            s=a(n+1)
            fn=n
            f=1.0+fn/2.0
            k=int(f)
            d=2.0d0
            if (s-a(k) .ne. 0) then
              do 210 i=1,l
                d=d*2.0d0
                if (s-a(k) .gt. 0) then
                  f=sngl(dble(f)+dble(fn)/d)
                else
                  f=sngl(dble(f)-dble(fn)/d)
                endif
                k=int(f)
  210         continue
            endif
c ...       insert a(n+1)  ....
            if(s.gt.a(k)) k=k+1
            do 250 i=k,n
              j=n+k-i
              ni(j+1)=ni(j)
              a(j+1)=a(j)
  250       continue
            ni(k)=ns
            a(k)=s
  260     continue
  300   continue
      endif

      return
      end
