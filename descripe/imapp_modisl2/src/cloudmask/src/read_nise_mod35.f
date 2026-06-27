      INTEGER FUNCTION READ_NISE( FILENAME, GRIDSIZE,
     &  NISE_NORTH, NISE_SOUTH )
      
c-----------------------------------------------------------------------
c !F77
c
c !DESCRIPTION:
c    To read a "Near Real-Time SSM/I EASE-Grid Daily Global Ice
c    Concentration and Snow Extent (NISE)" data file from NSIDC.
c    Information on these files is available at
c    http://www-nsidc.colorado.edu/NSIDC/CATALOG/catalog_index.html
c
c    This function reads the northern and southern hemisphere azimuthal
c    equal area grids which are stored at 25 km resolution in HDF-EOS.
c
c    To obtain NISE files (which are updated daily), contact NSIDC at
c    nsidc@kryos.colorado.edu
c
c !INPUT PARAMETERS:
c    FILENAME      Name of the NISE file
c    GRIDSIZE      Dimension for output arrays NISE_N and NISE_S
c                  (GRIDSIZE=721 for 25 km azimuthal grid)
c
c !OUTPUT PARAMETERS:
c    READ_NISE     Success flag
c                   0 => Success
c                  -1 => Error opening FILENAME
c                  -2 => Error reading northern hemisphere grid
c                  -3 => Error reading southern hemisphere grid
c    NISE_NORTH    Northern hemisphere data grid
c    NISE_SOUTH    Southern hemisphere data grid
c
c    Data grid value     Meaning
c             0          snow-free land
c         1-100          sea ice concentration percentage
c           101          permanent ice (Greenland, Antarctica)
c           102          not used
c           103          dry snow
c           104          wet snow
c       105-251          not used
c           252          mixed pixels at coastlines
c           253          suspect ice value
c           254          corners (undefined)
c           255          ocean
c           
c !REVISION HISTORY:
c
c !TEAM-UNIQUE HEADER:
c     Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
c
c !END
c-----------------------------------------------------------------------

      IMPLICIT NONE

c ... Input arguments

      CHARACTER filename*(*)
      INTEGER gridsize

c ... Output arguments

      CHARACTER*1 nise_north( gridsize, gridsize )
      CHARACTER*1 nise_south( gridsize, gridsize )

c ... Local variables
      
      INTEGER file_id, grid_id, status
      
      INTEGER xsize, ysize
c rhucek 052500: resized upleft and lowrgt to 2 elements
c     DOUBLE PRECISION upleft, lowrgt
      DOUBLE PRECISION upleft(2), lowrgt(2)
      
      INTEGER start( 2 ), stride( 2 ), edge( 2 )
      
c ... HDF-EOS functions
      
      INTEGER gdopen, gdattach, gddetach, gdclose, gdgridinfo, gdrdfld
      EXTERNAL gdopen, gdattach, gddetach, gdclose, gdgridinfo, gdrdfld

      INTEGER DFACC_READ
      PARAMETER ( DFACC_READ = 1 )
      
c-----------------------------------------------------------------------
c     OPEN FILE
c-----------------------------------------------------------------------

c ... Open file

      file_id = gdopen( filename, DFACC_READ )
      if ( file_id .eq. -1 ) then
        read_nise = -1
        return
      endif
      
c-----------------------------------------------------------------------
c     READ NORTHERN HEMISPHERE GRID
c-----------------------------------------------------------------------

c ... Open northern hemisphere grid

      grid_id = gdattach( file_id, 'Northern Hemisphere' )
      
c ... Get grid information

      status = gdgridinfo( grid_id, xsize, ysize, upleft, lowrgt )

c ... Read grid data

      start( 1 ) = 0
      start( 2 ) = 0
      stride( 1 ) = 1
      stride( 2 ) = 1
      edge( 1 ) = xsize
      edge( 2 ) = ysize
      status = gdrdfld( grid_id, 'Extent', start, stride, edge,
     &  nise_north )
      if ( status .lt. 0 ) then
        read_nise = -2
        return
      endif
      
c     Call routine to fill in snow/ice info for coastal regions as 
c     much as possible (northern hemisphere).
      call massage_snowice(nise_north,xsize,ysize)

c ... Close northern hemisphere grid

      status = gddetach( grid_id )

c-----------------------------------------------------------------------
c     READ SOUTHERN HEMISPHERE GRID
c-----------------------------------------------------------------------

c ... Open southern hemisphere grid

      grid_id = gdattach( file_id, 'Southern Hemisphere' )

c ... Get grid information

      status = gdgridinfo( grid_id, xsize, ysize, upleft, lowrgt )

c ... Read grid data

      start( 1 ) = 0
      start( 2 ) = 0
      stride( 1 ) = 1
      stride( 2 ) = 1
      edge( 1 ) = xsize
      edge( 2 ) = ysize
      status = gdrdfld( grid_id, 'Extent', start, stride, edge,
     &  nise_south )
      if ( status .lt. 0 ) then
        read_nise = -3
        return
      endif
      
c     Call routine to fill in snow/ice info for coastal regions as 
c     much as possible (southern hemisphere).
      call massage_snowice(nise_south,xsize,ysize)

c ... Close southern hemisphere grid

      status = gddetach( grid_id )

c-----------------------------------------------------------------------
c     CLOSE FILE
c-----------------------------------------------------------------------

c ... Close file and return success flag

      status = gdclose( file_id )
      read_nise = 0

      END



      subroutine massage_snowice(map_nise,xsize,ysize)

c---------------------------------------------------------------------
c!F77
c
c!Description:
c     Subroutine for filling in snow and ice information in coastal
c     regions.  For all coast pixels in the NISE data, search both
c     latitudinally and longitudinally for non-coast pixels which
c     indicate snow or ice.  Stop search when non-snow/ice, non-coast
c     pixels are found.  For those coast pixels which have snow/ice
c     pixels indicated on their boundaries, set flag to indicate 
c     snow/ice.
c
c!Input Parameters:
c map_nise    Input/output snow/ice data grid
c xsize       Grid size in the horizontal direction
c ysize       Grid size in the vertical direction
c
c!Output Parameters:
c map_nise    Input/output snow/ice data grid
c
c!Revision History:
c Original subroutine:
c     03-09-04  R. Frey
c
c!Team-Unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c
c!End
c----------------------------------------------------------------------

      implicit none
      save

c     Arguments.
      integer xsize, ysize
      character*1 map_nise(xsize,ysize)

c     Local arrays.
      integer k(4)

c     Local scalars.
      integer temp_nise, i, j, m, n

c     Loop through map data.

      do j = 1, xsize

        do i = 1, ysize

          temp_nise = 0

c         Check for coast pixel (value = 252).

          if( ichar( map_nise(j,i) ) .eq. 252 ) then

c           Check if on map edges.  Cannot search on all four sides.
c           No data in corners.

            if(j .eq. 1) then
              k(1) = ichar(map_nise(j+1,i))
              k(2) = ichar(map_nise(j,i-1))
              k(3) = ichar(map_nise(j,i+1))
              k(4) = 0
         
            else if(i .eq. 1) then
              k(1) = ichar(map_nise(j,i+1))
              k(2) = ichar(map_nise(j+1,i))
              k(3) = ichar(map_nise(j-1,i))
              k(4) = 0

            else if(j .eq. xsize) then
              k(1) = ichar(map_nise(j-1,i))
              k(2) = ichar(map_nise(j,i-1))
              k(3) = ichar(map_nise(j,i+1))
              k(4) = 0
         
            else if(i .eq. ysize) then
              k(1) = ichar(map_nise(j,i-1))
              k(2) = ichar(map_nise(j+1,i))
              k(3) = ichar(map_nise(j-1,i))
              k(4) = 0

            else

c             Get all four adjacent values.
              k(1) = ichar(map_nise(j,i-1))
              k(2) = ichar(map_nise(j,i+1))
              k(3) = ichar(map_nise(j-1,i))
              k(4) = ichar(map_nise(j+1,i))

            end if

c           Fill in missing snow/ice values (=200).
            do m = 1,4

              n = k(m)
              if(n .ne. 0) then

                if((n .eq. 103 .or. n .eq. 104 .or. n .eq. 200) .or. 
     *                             (n .gt. 25 .and. n .lt. 102)) then
                  temp_nise = 200
                end if

              end if
            enddo

            if (temp_nise .eq. 200) map_nise(j,i) = char(temp_nise)

          end if

        enddo

      enddo

      return
      end
