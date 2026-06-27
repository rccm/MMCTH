      subroutine gphite( p, t, w, z_sfc, n_levels, i_dir, z)
c-----------------------------------------------------------------------
c
c!F77
c
c!Description:
c  Routine to compute geopotential height given the atmospheric state.
c    Includes virtual temperature adjustment.
c    Dimension of height array may not not be the same as that of the
c    input profile data.
c
c!Input Parameters:
c       p     - REAL*4 pressure array (mb)
c       t     - REAL*4 temperature profile array (K)
c       w     - REAL*4 water vapour profile array (g/kg)
c     z_sfc   - REAL*4 surface height (m).  0.0 if not known.
c    n_levels - INT*4 number of elements used in passed arrays
c     i_dir   - INT*4 direction of increasing layer number
c                 i_dir = +1, Level(1) == p(top)         } satellite/AC
c                             Level(n_levels) == p(sfc)  }    case
c
c                 i_dir = -1, Level(1) == p(sfc)         } ground-based
c                             Level(n_levels) == p(top)  }    case
c
c!Output Parameters:
c       z     - REAL*4 pressure level height array (m)
c
c!Revision History:
c    06-DEC-1997 Liam Gumley, CIMSS/SSEC
c                Modified version from Hal Woolf to meet ECS standards
c                and use PGS toolkit open/close routines.
c
c!Team-unique Header:
c
c!End
c
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c              -- Prevent implicit typing of variables --
c-----------------------------------------------------------------------

      implicit none

c-----------------------------------------------------------------------
c                           -- Arguments --
c-----------------------------------------------------------------------

c -- Arrays

      real*4    p(*), t(*), w(*), 
     +          z(*)

c -- Scalars

      integer*4 n_levels, i_dir

      real*4    z_sfc

c-----------------------------------------------------------------------
c                         -- Local variables --
c-----------------------------------------------------------------------

c -- Parameters

      real*4    rog, fac
      parameter ( rog = 29.2898, 
     +            fac = 0.5 * rog )

c -- Scalars

      integer*4 i_start, i_end, l

      real*4    v_lower, v_upper, algp_lower, algp_upper, hgt

      save
      
c***********************************************************************
c                         ** Executable code **
c***********************************************************************

c-----------------------------------------------------------------------
c  -- Calculate virtual temperature adjustment and exponential       --
c  -- pressure height for level above surface.  Also set integration --
c  -- loop bounds                                                    --
c-----------------------------------------------------------------------

      if( i_dir .gt. 0 )then

c       --------------------
c       Data stored top down
c       --------------------

        v_lower = t(n_levels) * ( 1.0 + ( 0.00061 * w(n_levels) ) )

        algp_lower = alog( p(n_levels) )

        i_start = n_levels-1
        i_end   = 1

      else

c       ---------------------
c       Data stored bottom up
c       ---------------------

        v_lower = t(1) * ( 1.0 + ( 0.00061 * w(1) ) )

        algp_lower = alog( p(1) )

        i_start = 2
        i_end   = n_levels

      end if

c-----------------------------------------------------------------------
c                     -- Assign surface height --
c-----------------------------------------------------------------------

      hgt = z_sfc

      if(i_dir .gt. 0) then
        z(n_levels) = z_sfc
      else
        z(1) = z_sfc
      end if

c-----------------------------------------------------------------------
c             -- Loop over layers always from sfc -> top --
c-----------------------------------------------------------------------

      do l = i_start, i_end, -1*i_dir

c       ----------------------------------------------------
c       Apply virtual temperature adjustment for upper level
c       ----------------------------------------------------

        v_upper = t(l)
        if( p(l) .ge. 300.0 )
     +    v_upper = v_upper * ( 1.0 + ( 0.00061 * w(l) ) )

c       ----------------------------------------------------- 
c       Calculate exponential pressure height for upper layer
c       ----------------------------------------------------- 

        algp_upper = alog( p(l) )

c       ----------------
c       Calculate height
c       ----------------

        hgt = hgt + ( fac*(v_upper+v_lower)*(algp_lower-algp_upper) )

c       -------------------------------
c       Overwrite values for next layer
c       -------------------------------

        v_lower = v_upper
        algp_lower = algp_upper

c       ---------------------------------------------
c       Store heights in same direction as other data
c       ---------------------------------------------

        z(l) = hgt

      end do

      return
      end
