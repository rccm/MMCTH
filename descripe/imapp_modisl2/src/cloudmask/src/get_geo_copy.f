      SUBROUTINE GET_GEO_COPY( FNAME_GEO, FNAME_MOD35, MAX5KM_LIN,
     &  MAX5KM_ELE )


c#######################################################################
C!F77
C
C!DESCRIPTION:
c    This is the driver program which passes the information
c    to the geolocation copy subroutine for atmospheric
c    MODIS products.
C
C!INPUT PARAMETERS:
C    FNAME_GEO      Name of MODIS L1B geolocation file
C    FNAME_MOD35    Name of existing UW MODIS L2 product file
c    MAX5KM_LIN     Maximum number of 5x5 km pixels in a line
c    MAX5KM_ELE     Maximum number of 5x5 km elements in a line
C
C!OUTPUT PARAMETERS:
c    None
C
C!REVISION HISTORY:
C
C!TEAM-UNIQUE HEADER:
C
C!END
C#######################################################################

      implicit none

c ... Scalar arguements

      integer max5km_lin,max5km_ele

c ... Array arguments

      character*(*)  fname_geo,fname_mod35

c ... Local scalars

      integer rtn

c ... Local arrays

      character*100 error_text
      integer start( 2 ), stride( 2 ), edge( 2 )

c ... External functions

      integer copy_info
      external copy_info

c ... Initialize

      error_text = ' '

c ... Set input parameters for UW MODIS L2 5x5 products
c ... (NOTE! start array is zero based)

      start( 1 )  = 2
      start( 2 )  = 2
      stride( 1 ) = 5
      stride( 2 ) = 5
      edge( 1 )   = max5km_ele
      edge( 2 )   = max5km_lin

c ... Call copy routine

      rtn = copy_info( fname_geo, fname_mod35,
     &  start, stride, edge, error_text )
      if ( rtn .ne. 0 )
     &  call message( 'get_geo_copy', error_text, rtn, 2 )

      end
