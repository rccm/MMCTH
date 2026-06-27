      subroutine Dstar(btd8,btd11,Dstar1)
c***********************************************************************
c
c     purpose
c        calculate D* dust parameter 
c
c     calling sequence
c        call Dstar 
c
c     variables
c       name      type  i/o  description
c       ----      ----  ---  -----------
c      btd8       r*4   i   brightness temperature difference (8-11)
c      btd11      r*4   i   brightness temperature difference(11-12)
c      Dstar1     r*4   o   D* parameter
c
c***********************************************************************
c
c     -- input parameters
      real btd8, btd11

c     -- output parameters
      real Dstar1
      
      real ratio
      
      Dstar1 = -999

c     -- intialize thermal constants 
      A = -0.05
      B = 10.0

c     --calculate D*. If ratio is high, D* will overflow. Set to -999 to 
c     --signal craziness.
      ratio = (btd11 - A)/ (btd8 - B)
      if (ratio .GE. 5.0) then
        Dstar1 = -999.0       
      else
        Dstar1 = exp(ratio)
      end if 

      return
      
c         write (6,1000) btd8,btd11 
c         write (6,1100) Dstar1
c1000  format (/'Subroutine Dstar'/'Input:   btd8 = ',f8.4,
c     1 ' btd11 = ',f8.4)
c1100  format ('Output:  Dstar1  = ',f8.4) 

      end
