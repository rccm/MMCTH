      SUBROUTINE READ_GOODE( NX, NY, LAT, LON, LUN, RESULT )

c-----------------------------------------------------------------------
c !F77
c
c !DESCRIPTION:
c     Given a set of lat/lon points, read an Interrupted Goode
c     Homolosine dataset from USGS of the type described at
c     http://edcwww.cr.usgs.gov/landdaac/glcc/globe_int.html and
c     http://edcwww.cr.usgs.gov/landdaac/1KM/goodesarticle.html
c
c     The data dimensions of the Interrupted Goode Homolosine
c     projection for the global land cover characteristics data set are
c     17347 lines (rows) and 40031 samples (columns).
c
c     Files of this type are exactly 17347 x 40031 = 694417757 bytes
c     in size, and files in use at UW include goge1_2_img, lst1km.v3.
c
c !INPUT PARAMETERS:
c     NX      First dimension of arrays LAT, LON, RESULT
c     NY      Second dimension of arrays LAT, LON, RESULT
c     LAT     Array of latitude values (degrees)
c     LON     Array of longitude values (degrees, -180W to +180E)
c     LUN     Logical unit number for input file
c             (must be opened by the caller for direct access with
c              a record length of 512 bytes. Caller must close the file).
c           
c !OUTPUT PARAMETERS:
c     RESULT  Array of land type values
c
c !REVISION HISTORY:
c
c !TEAM-UNIQUE HEADER:
c     Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
c
c !END
c-----------------------------------------------------------------------
      
      IMPLICIT NONE

c ... Arguments
      
      INTEGER nx, ny
      REAL lat( nx, ny ), lon( nx,ny )
      INTEGER lun
      BYTE result( nx, ny )

c ... Parameters

      INTEGER max_row, max_col
      PARAMETER ( max_row = 17347, max_col = 40031 )

      INTEGER recsize
      PARAMETER ( recsize = 512 )

      INTEGER max_count
      PARAMETER ( max_count = 3000000 )
      
      INTEGER max_rec
      PARAMETER ( max_rec = ( max_row * max_col ) / recsize )

c ... Local variables
      
      INTEGER max_pixel, max_line, pixel, line, row, col

      DOUBLE PRECISION dlat, dlon, drow, dcol

      INTEGER loc_value( max_count ), loc_index( max_count )
            
      INTEGER loc, count, kflag, i, recnum, elenum, oldrec, ios, rectest

      BYTE recdata( recsize )

      LOGICAL opentest
      
c-----------------------------------------------------------------------
c     CHECK ARGUMENTS
c-----------------------------------------------------------------------

c ... Check that input arrays do not exceed local work array size
                     
      if ( nx * ny .gt. max_count ) then
        call message( 'read_goode.f',
     &    'NX*NY greater than MAX_COUNT ' //
     &    '[OPERATOR ACTION: Contact SDST]', 0, 2 )
      endif

c ... Check that input LUN is open

      inquire( lun, opened=opentest )
      if ( .not. opentest ) then
        call message( 'read_goode.f',
     &    'Input file LUN is not open ' //
     &    '[OPERATOR ACTION: Contact SDST]', 0, 2 )
      endif
      
c ... Check that input LUN was opened with the correct record size

      inquire( lun, recl=rectest )
      if ( rectest .ne. recsize ) then
        call message( 'read_goode.f',
     &    'Record size of input file does not match RECSIZE ' //
     &    '[OPERATOR ACTION: Contact SDST]', 0, 2 )
      endif
      
c-----------------------------------------------------------------------
c     COMPUTE BYTE LOCATIONS FOR ALL LINES/PIXELS
c-----------------------------------------------------------------------

c ... Begin loop over all lines and pixels

      max_pixel = nx
      max_line  = ny
      count = 0
      
      do line = 1, max_line
      
        do pixel = 1, max_pixel

c ...     Increment pixel count

          count = count + 1

c ...     Compute row/column using USGS routine (original version from
c ...     ftp://edcftp.cr.usgs.gov/pub/software/misc/gihll2ls.c)

          dlat = dble( lat( pixel, line ) )
          dlon = dble( lon( pixel, line ) )
          call getcoord( dlat, dlon, drow, dcol )
          row = nint( drow )
          col = nint( dcol )
          
c ...     Constrain row/column to allowed range
          
          row = max( min( row, max_row ), 1 )
          col = max( min( col, max_col ), 1 )

c ...     Compute byte location

          loc = ( row - 1 ) * max_col + col

c ...     Save location value and index

          loc_value( count ) = loc
          loc_index( count ) = count

        end do
        
      end do

c ... Sort the location value array in increasing order,
c ... and re-order the index array the same way

      kflag = 2
      call isort( loc_value, loc_index, count, kflag )

c-----------------------------------------------------------------------
c     READ DATA AT SORTED BYTE LOCATIONS
c-----------------------------------------------------------------------

c ... Loop through sorted byte locations

      oldrec = 0
      
      do i = 1, count
      
c ...   Compute record number and element number within record,
c ...   assuming the first record in the file is record one,
c ...   and the first element in each record is element one.
      
        recnum = ( loc_value( i ) - 1 ) / recsize + 1
        elenum = mod( loc_value( i ) - 1, recsize ) + 1

c ...   Constrain record number to allowed range
c ...   (note that this only has an effect on the bottom right hand
c ...    side of the projection, where all data values are 'ocean'.
c ...    See http://edcwww.cr.usgs.gov/landdaac/1KM/FIG4.gif).
          
        row = max( min( recnum, max_rec ), 1 )

c ...   If record number is new, read the data

        if ( recnum .ne. oldrec ) then

          read( lun, rec=recnum, iostat=ios ) recdata
          if ( ios .ne. 0 ) then
            call message( 'read_goode.f',
     &        'Error reading input file ' //
     &        '[OPERATOR ACTION: Contact SDST]', 0, 2 )
          endif
          
          oldrec = recnum

        endif

c ...   Compute the line and pixel number

        line  = ( loc_index( i ) - 1 ) / max_pixel + 1
        pixel = mod( loc_index( i ) - 1, max_pixel ) + 1

c ...   Store the data for this pixel/line

        result( pixel, line ) = recdata( elenum )
                    
      end do

      END
