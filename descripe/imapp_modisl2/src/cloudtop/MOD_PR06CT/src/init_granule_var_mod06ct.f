      subroutine init_granule_var_mod06ct( ssctpr, slcd, smcd, 
     *  shcd, sthncd, sthkcd, sopcd, scicd,sco2,sbias )
      
c-----------------------------------------------------------------------
c
c!F77
c
c!Description:
c    Initialize counters for metadata.
c
c!Input Parameters:
c
c    Counters for metadata:
c    ssctpr                      number of successful retrievals
c    slcd                        number of low clouds
c    smcd                        number of mid-level clouds
c    shcd                        number of high clouds
c    sthncd                      number of "thin clouds"
c    sthkcd                      number of "thick clouds"
c    sopcd                       number of opaque clouds
c    scicd                       number of cirrus clouds
c    sco2                        number of clouds from co2 slicing
c    scbias                      number of clouds with bias adjustment
c
c!Output Parameters:
c
c    Counters for metadata:
c    ssctpr                      number of successful retrievals
c    slcd                        number of low clouds
c    smcd                        number of mid-level clouds
c    shcd                        number of high clouds
c    sthncd                      number of "thin clouds"
c    sthkcd                      number of "thick clouds"
c    sopcd                       number of opaque clouds
c    scicd                       number of cirrus clouds
c    sco2                        number of clouds from co2 slicing
c    scbias                      number of clouds with bias adjustment
c
c!Revision History:
c $Id: init_granule_var_mod06ct.f,v 1.1.1.1 2005/02/22 17:15:54 gumley Exp $
c
c!Team-unique Header:
c
c!End
c
c-----------------------------------------------------------------------

      implicit none
      save
          
c ... arguments

      integer ssctpr,slcd,smcd,shcd,sthncd,sthkcd,sopcd,scicd,sco2,sbias


c     Initialize counters
      ssctpr = 0
      slcd = 0
      smcd = 0
      shcd = 0
      sthncd = 0
      sthkcd = 0
      sopcd = 0
      scicd = 0
      sco2 = 0
      sbias = 0

      return
      end
