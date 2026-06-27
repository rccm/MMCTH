      subroutine init_granule_var_mod06ir( ni, nw, nx, nu, ng_irp )

c-----------------------------------------------------------------------
c
c!F77
c
c!Description:
c    Initialize granule based counter variables for MOD_PR06IR
c
c!Input Parameters:
c
c    Counters for metadata:
c    ni                          number of ice clouds
c    nw                          number of water clouds
c    nx                          number of mixed phase clouds
c    nu                          number of uncertain phase clouds
c    ng_irp                      number if good irphase retrievals
c
c!Output Parameters:
c
c    Initialized counters for metadata:
c    ni                          number of ice clouds
c    nw                          number of water clouds
c    nx                          number of mixed phase clouds
c    nu                          number of uncertain phase clouds
c    ng_irp                      number if good irphase retrievals
c
c!Revision History:
c $Id: init_granule_var_mod06ir.f,v 1.4 1999/04/16 22:45:06 kis Exp $
c
c!Team-unique Header:
c
c!End
c
c-----------------------------------------------------------------------
c
c ... arguments
      integer ni, nw, nx, nu, ng_irp

c ... Initialize counters
      ni = 0
      nw = 0
      nx = 0
      nu = 0
      ng_irp = 0

      return
      end

