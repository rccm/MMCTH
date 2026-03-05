module MOD06AlbedoEcoModule

! HI
!****************************************************************************
  ! !F90
  !
  ! !Description:
  !    This module contains the subroutines used to provide albedo and ecosystem
  !     maps for processing a MOD06 granule.
  !    There is only one callable routine, getAlbedoEco, which returns albedo 
  !     values and IGBP ecosystem classification for each pixel of the MOD06 granule,
  !     and for the wavelengths specified.  Also returned are values of snow albedos
  !     by ecosystem class for the wavelengths specified.
  !  
  ! !Callable routines:
  !    getAlbedoEco()
  ! 
  ! !Revision History:
  !
  ! Revision 1.0  2003/12/18  12:43:43  EGMoody
  ! Initial revision.
  !
  ! !Team-Unique Header:
  !   Cloud Retrieval Group, NASA Goddard Space Flight Center
  !
  ! !References and Credits:
  !   Written by
  !    Eric Moody
  !    Climate and Radiation Branch, Code 913
  !    NASA/GSFC
  !    Greenbelt MD 20771
  !    Eric.Moody@gsfc.nasa.gov
  !
  ! !Design Notes:
  !
  ! !END
  !*****************************************************************************

  !Dependencies:
!  use hdf_mod
  !use hdf,           only :  MAX_VAR_DIMS, SFstart, SFend,                    &
  !                           SFn2index, SFselect, SFrdata,                    &
  !                           SFdimid, SFendacc,                               &
  !                           SFfinfo, SFginfo, SFgdinfo, SFgainfo,            & 
  !                           SFfattr, SFrcatt, SFscatt, SFrnatt, SFsnatt,     &
  !                           DFACC_CREATE, DFACC_WRITE, DFACC_READ, FAIL,     &
  !                           DFNT_CHAR, DFNT_FLOAT, DFNT_DOUBLE,              &
  !                           DFNT_INT8, DFNT_UINT8,                           &
  !                           DFNT_INT16, DFNT_UINT16,                         &
  !                           DFNT_INT32, DFNT_UINT32
!  use GeneralAuxType
  use nonscience_parameters
  implicit none

   include "hdf.f90"
   include "dffunc.f90"
 
  
  private

  public :: GetNISEType, ReadSnowAlbStats, init_NISE_processing


  !local variables:
  !counters:
  integer                          :: i,j,k,l,m,n

  !HDF variables:
  integer                           :: HDFstatus
  integer, dimension(10)                         :: hdfStart, hdfStride, hdfEdge
  integer                          :: sds_id, sds_index
  !  HDFstatus: Used for error checking.
  !  hdfstart : Follows HDF conventions. The array element from which to begin reading the
  !              HDF array. Zero based.
  !  hdfstride: Follows HDF conventions. The frequency with which data should be read in
  !              each dimension. Stride 1 means read every element.
  !  hdfedge  : Follows HDF conventions. The number of elements to read in each dimension. 
  !              If start(:) == 0 and stride(:) == 1, setting edge equal to the shape of 
  !              the underlying array reads all elements. 


  INTEGER, PARAMETER                              :: gridsize = 721  
  integer*2, DIMENSION(gridsize, gridsize) :: niseNorth, niseSouth


contains

	subroutine init_NISE_processing(nise_file)
	
		character(len=*), intent(in) :: nise_file
		integer :: iret
	
		iret = Read_NISE(nise_file, gridsize, niseNorth, niseSouth)


	
	end subroutine init_NISE_processing



  ! ----------------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------------

  SUBROUTINE GetNISEType(lon, lat, snowIceType, errorLevel)
  !...............................................................................
  !
  !...............................................................................
  ! !F90
  !
  !   !DESCRIPTION: 
  !        return snow ice type using "Near Real-Time SSM/I EASE-Grid Daily 
  !        Global Ice Concentration and Snow Extent (NISE)" data file from NSIDC.
  !        If NISE is not available, leave the snowIceType array untouched. 
  ! 
  !
  !   !INPUT PARAMETERS:
  !        
  !        lat: latitudes  (in degrees:  -90 to  90)
  !        lon: longitudes (in degrees: -180 to 180)
  !
  !   !OUTPUT PARAMETERS:
  !        snowIceType: 
  !                  0          snow-free land
  !              1-100          sea ice concentration percentage
  !                101          permanent ice (Greenland, Antarctica)
  !                102          not used
  !                103          dry snow
  !                104          Not Used / 2000-March 2002 data Wet Snow
  !            105-251          not used
  !                252          mixed pixels at coastlines
  !                253          suspect ice value
  !                254          corners (undefined)
  !                255          ocean
  ! 
  !        errorLevel: return error code, 0 (success) or 1 (fail)
  !...............................................................................
  !   !Revision History:
  !
  !...............................................................................
  !   !Team-unique Header:
  !    Cloud Retrieval Group, NASA/GSFC
  !
  !   !PROGRAMMER:
  !            Jason Li (SM&A)
  !            Climate and Radiation Branch
  !            NASA Goddard Space Flight Center
  !            Greenbelt, Maryland 20771, U.S.A.
  !
  !-------------------------------------------------------------------------------
  !   !END

	implicit none

 
  REAL,        DIMENSION(:,:) :: lat, lon
  INTEGER,    DIMENSION(:,:) :: snowIceType
  INTEGER,    INTENT(OUT)    :: errorLevel


  INTEGER :: version, xdim, ydim, iret, nise
  REAL    :: x, y, pixel_lon, pixel_lat


  ! Liam Gumley's NISE reader:
  INTEGER, external :: EZLH_Convert


  !**********

  !Set error level
  errorLevel = 0

  !Determine x, y sizes:
  xdim = SIZE(snowIceType, dim=1)
  ydim = SIZE(snowIceType, dim=2)


  ! Get cell coordinates for southern or northern hemisphere.
  DO m = 1, xdim
  DO n = 1, ydim
  
        pixel_lon = lon(m,n)
        pixel_lat = lat(m,n)
  
        x = min( max( pixel_lon, -179.99 ),  179.99 )
        y = min( max( pixel_lat,  -89.99 ),   89.99 )
  
        IF ( y < 0.0 ) THEN
             iret = EZLH_Convert( 'Sl', y, x, i, j )
        ELSE
             iret = EZLH_Convert( 'Nl', y, x, i, j )
        ENDIF
 
        !  Save output NISE data for southern or northern hemisphere
        ! NISE data grid value definitions:
        !
        !             Data Grid Value     Meaning
        !                  0          snow-free land
        !              1-100          sea ice concentration percentage
        !                101          permanent ice (Greenland, Antarctica)
        !                102          not used
        !                103          dry snow
        !                104          Not Used / 2000-March 2002 data Wet Snow
        !            105-251          not used
        !                252          mixed pixels at coastlines
        !                253          suspect ice value
        !                254          corners (undefined)
        !                255          ocean
  
        IF ( iret /= 0 ) THEN
             errorLevel = 1
        ELSE
             IF( y < 0.0 ) THEN
                 nise = niseSouth( i, j ) 
             ELSE
                 nise = niseNorth( i, j ) 
             ENDIF
             
			 
             !Set the snow/ice type and flag:             
             SELECT CASE (nise)
                    CASE(1:100)
                         !sea ice:
                         snowIceType(m,n) = nise
                    CASE(101)
                         !perm ice:
                         snowIceType(m,n) = 101
                    CASE(103)
                         !dry snow:
                         snowIceType(m,n) = 103
                    CASE(104)
                         !wet snow:
                         !For the 2000-March 2002 data, the snow is further divided by
                         ! wet and dry snow.  Per the NISE documentation, the wet class
                         ! was inaccurate and thereby no longer used after March 2002.
                         ! For reprocesing of 2000-March 2002 data, set this value to 
                         ! dry snow so that the collection has consistent snow albedo
                         ! values.  This is the change for collection 5.
                         snowIceType(m,n) = 103
                    CASE(252)
                         !mixed scene pixels at coasts
                         snowIceType(m,n) = 252
                    CASE(200)
                         !Fill scene pixels:
                         snowIceType(m,n) = 200
                    CASE DEFAULT                  
                         !no snow/ice, so set this type to 0:
                         snowIceType(m,n) = 0
             END SELECT
             
        ENDIF
  END DO
  END DO


  END SUBROUTINE GetNISEType

  ! ----------------------------------------------------------------------------------
  ! ----------------------------------------------------------------------------------

  subroutine ReadSnowAlbStats  ( StatsFN, NumSnowTypes, NumAlbBnds, numEco, &
                                 AlbedoMean,  errorLevel                    )
                                 
    character (len = *),                 intent( in)  :: StatsFN
    real      ,  &
              dimension(:,0:,:),         intent(inout):: AlbedoMean
    integer   , intent( in)  :: NumAlbBnds,       &
                                                         numEco,           &
                                                         NumSnowTypes
    INTEGER,            INTENT(OUT) :: errorLevel
                                                         
    ! !Description:
    !   This routine will output only a single band of albedo, for the specified amount
    !    of data, to the specified position within the hdf file.
    !  
    ! !Input Parameters:
    !    StatsFN : The name, only, of the new HDF file.
    !    Albedo          : Contains the single band of albedo to be stored.
    !    NumSnowTypes    : 2, Dry or wet, via NISE classes.
    !    NumAlbBnds      : Number of Albedo wavelengths.
    !    numEco          : Number of Ecosystem Classes.
    !
    ! !Output Parameters:
    !    AlbedoMean      : Mean statistic
    !
    ! !Revision History:
    !   See Module revision history at the beginning of the file.
    !
    ! !Team-Unique Header:
    !   Cloud Retrieval Group, NASA Goddard Space Flight Center
    !
    ! !References and Credits:
    !   Written by
    !    Eric Moody
    !    Climate and Radiation Branch, Code 913
    !    NASA/GSFC
    !    Greenbelt MD 20771
    !    Eric.Moody@gsfc.nasa.gov
    !
    ! !Design Notes:
    !
    ! !END
    
    !local variables: 
    !HDF variables:
    integer               :: status
    integer , dimension(10) :: hdfStart, hdfStride, hdfEdge
    integer               :: sds_id, newHDFID, sds_index
    character(len =  200)                :: SDSName
    real      ,  &
              dimension(1:2,1:NumAlbBnds,0:numEco,1:1,1:2)   :: DummyAlb
    !  status    : Used for error checking.
    !  hdfstart  : Follows HDF conventions. The array element from which to begin reading the
    !               HDF array. Zero based.
    !  hdfstride : Follows HDF conventions. The frequency with which data should be read in
    !               each dimension. Stride 1 means read every element.
    !  hdfedge   : Follows HDF conventions. The number of elements to read in each dimension. 
    !               If start(:) == 0 and stride(:) == 1, setting edge equal to the shape of 
    !               the underlying array reads all elements.
    !  sds_id,
    !   newHDFID : HDF SDS ID.
    !  sds_index : Index of the SDS in the HDF file.
    ! SDSName    : Stores the name of the SDS being procesed.


  !Set error level
  errorLevel = success 

    !************************************************************************************
    ! Open the input file:
    !************************************************************************************
    !Open the albedo HDF file:
    newHDFID = SFstart( trim(StatsFN), DFACC_READ )
    if (newHDFID == Fail) then
       Print *, "Can't Reopen the NEW HDF File", trim(StatsFN)
       errorLevel = fail
       return
    end if


    !************************************************************************************
    ! Input the Albedo mean statistic:
    !************************************************************************************   
    !Set up the output hdf variables, in this case we are writing a Belt:
    ! hdfStart, note that we are storing the entire Longitude, Albedo Bands and Ecosystems
    ! so these start at 0, however we are storing a box of Latitude, so this starts at the
    ! starty counter.
    hdfStart ( : ) = 0
    hdfStart ( 5 ) = 0  !Latitude, 0 = NH, 1 = SH
    hdfStart ( 4 ) = 0  !Longitude, only 1, so 0
    hdfStart ( 3 ) = 0  !Eco
    hdfStart ( 2 ) = 0  !AlbBands
    hdfStart ( 1 ) = 0  !NumSnowTypes

    !Strides = 1, since we are not skiping points
    hdfStride( : ) = 1
    
    !Edge is the total number of points being read:
    ! It is the number of Albedo bands, the number of Ecosystems, the Number of Longitude
    ! boxes, and 1 lat box:
    hdfEdge  ( : ) = 1
    hdfEdge  ( 5 ) = 2
    hdfEdge  ( 4 ) = 1
    hdfEdge  ( 3 ) = numEco
    hdfEdge  ( 2 ) = NumAlbBnds
    hdfEdge  ( 1 ) = NumSnowTypes
    
    !Read the Mean data:
    SDSName = 'Snow_Albedo_Year_Mean'
    !determine the sds_index for the SDS:
    sds_index = SFn2index( newHDFID, trim(SDSName) )
    !get access to this SDS:
    sds_id = SFselect(newHDFID,sds_index)
        
    !Read the data:
    status = SFrdata(sds_id, hdfStart, hdfStride, hdfEdge, DummyAlb)

    !Store the dummy data in the final array:
    do i = 1, NumAlbBnds
      do j = 0, numEco-1
        do k = 1, NumSnowTypes
           AlbedoMean(i,j,k) = DummyAlb(k,i,j,1,1)
        end do
      end do
    end do

    !Error Checking:
    if (sds_index == fail .or. &
        sds_id    == fail .or. &
        status    == fail)     then
        errorLevel = fail
        return
     end if
    
    !End Access to the SDS
    status = sfendacc(sds_id)
    if (status  == fail) then
       errorLevel = fail
       return
    end if
    !************************************************************************************
    ! Close the file:
    !************************************************************************************
    !Close the  HDF file:
    status = SFend(newHDFID)
    if (status == fail) then
       Print *, "Can't End access to the NEW HDF File", trim(StatsFN)
       errorLevel = fail
       return
    end if


  end subroutine ReadSnowAlbStats



      INTEGER FUNCTION READ_NISE( FILENAME, GRIDSIZE, NISE_NORTH, NISE_SOUTH )
      
!-----------------------------------------------------------------------
! !F77
!
! !DESCRIPTION:
!    To read a "Near Real-Time SSM/I EASE-Grid Daily Global Ice
!    Concentration and Snow Extent (NISE)" data file from NSIDC.
!    Information on these files is available at
!    http://www-nsidc.colorado.edu/NSIDC/CATALOG/catalog_index.html
!
!    This function reads the northern and southern hemisphere azimuthal
!    equal area grids which are stored at 25 km resolution in HDF-EOS.
!
!    To obtain NISE files (which are updated daily), contact NSIDC at
!    nsidc@kryos.colorado.edu
!
! !INPUT PARAMETERS:
!    FILENAME      Name of the NISE file
!    GRIDSIZE      Dimension for output arrays NISE_N and NISE_S
!                  (GRIDSIZE=721 for 25 km azimuthal grid)
!
! !OUTPUT PARAMETERS:
!    READ_NISE     Success flag
!                   0 => Success
!                  -1 => Error opening FILENAME
!                  -2 => Error reading northern hemisphere grid
!                  -3 => Error reading southern hemisphere grid
!    NISE_NORTH    Northern hemisphere data grid
!    NISE_SOUTH    Southern hemisphere data grid
!
!    Data grid value     Meaning
!             0          snow-free land
!         1-100          sea ice concentration percentage
!           101          permanent ice (Greenland, Antarctica)
!           102          not used
!           103          dry snow
!           104          wet snow
!       105-251          not used
!           252          mixed pixels at coastlines
!           253          suspect ice value
!           254          corners (undefined)
!           255          ocean
!           
! !REVISION HISTORY:
!
! !TEAM-UNIQUE HEADER:
!     Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
!
! !END
!-----------------------------------------------------------------------

      IMPLICIT NONE

	include "hdf.f90"
	include "dffunc.f90"

! ... Input arguments

      CHARACTER filename*(*)
      INTEGER gridsize

! ... Output arguments

      integer*2 nise_north( gridsize, gridsize )
      integer*2 nise_south( gridsize, gridsize )

! ... Local variables
      
      INTEGER file_id, grid_id, status
      
      INTEGER xsize, ysize
! rhucek 052500: resized upleft and lowrgt to 2 elements
!     DOUBLE PRECISION upleft, lowrgt
      DOUBLE PRECISION upleft(2), lowrgt(2)
      
      INTEGER start( 2 ), stride( 2 ), edge( 2 )
      
	  integer*1 temp(gridsize, gridsize)
	  
	  character*80 dummy_name
	  integer rank, tp, nattr, i, j
	  
! ... HDF-EOS functions
      
!      INTEGER gdopen, gdattach, gddetach, gdclose, gdgridinfo, gdrdfld
!      EXTERNAL gdopen, gdattach, gddetach, gdclose, gdgridinfo, gdrdfld

      
!-----------------------------------------------------------------------
!     OPEN FILE
!-----------------------------------------------------------------------

! ... Open file

!      file_id = gdopen( filename, DFACC_READ )
	  file_id = sfstart(filename, DFACC_READ)
      if ( file_id .eq. -1 ) then
        read_nise = -1
        return
      endif
      
!-----------------------------------------------------------------------
!     READ NORTHERN HEMISPHERE GRID
!-----------------------------------------------------------------------

! ... Open northern hemisphere grid

!      grid_id = gdattach( file_id, 'Northern Hemisphere' )
	  grid_id = sfselect(file_id, 0)
      
! ... Get grid information

!      status = gdgridinfo( grid_id, xsize, ysize, upleft, lowrgt )
	  status = sfginfo(grid_id, dummy_name, rank, edge, tp, nattr)

! ... Read grid data

      start( 1 ) = 0
      start( 2 ) = 0
      stride( 1 ) = 1
      stride( 2 ) = 1
!      edge( 1 ) = xsize
!      edge( 2 ) = ysize
!      status = gdrdfld( grid_id, 'Extent', start, stride, edge, nise_north )
	  status = sfrdata(grid_id, start, stride, edge, temp)
      if ( status .lt. 0 ) then
        read_nise = -2
        return
      endif
		
		do i=1, edge(1)
			do j=1, edge(2)
				nise_north(i,j) = temp(i,j)
				if (temp(i,j) == -4) nise_north(i,j) = 252
				if (temp(i,j) == -3) nise_north(i,j) = 253
				if (temp(i,j) == -2) nise_north(i,j) = 254
				if (temp(i,j) == -1) nise_north(i,j) = 255
			end do
		end do
	

	  
!     Call routine to fill in snow/ice info for coastal regions as 
!     much as possible (northern hemisphere).
      call massage_snowice(nise_north,edge(1),edge(2))

! ... Close northern hemisphere grid

!      status = gddetach( grid_id )
	status = sfendacc(grid_id)

!-----------------------------------------------------------------------
!     READ SOUTHERN HEMISPHERE GRID
!-----------------------------------------------------------------------

! ... Open southern hemisphere grid

!      grid_id = gdattach( file_id, 'Southern Hemisphere' )
	grid_id = sfselect(file_id, 2)

! ... Get grid information

!      status = gdgridinfo( grid_id, xsize, ysize, upleft, lowrgt )
	  status = sfginfo(grid_id, dummy_name, rank, edge, tp, nattr)

! ... Read grid data

      start( 1 ) = 0
      start( 2 ) = 0
      stride( 1 ) = 1
      stride( 2 ) = 1
!      edge( 1 ) = xsize
!      edge( 2 ) = ysize
!      status = gdrdfld( grid_id, 'Extent', start, stride, edge, nise_south )
	  status = sfrdata(grid_id, start, stride, edge, temp)
      if ( status .lt. 0 ) then
        read_nise = -3
        return
      endif
      
		do i=1, edge(1)
			do j=1, edge(2)
				nise_south(i,j) = temp(i,j)
				if (temp(i,j) == -4) nise_south(i,j) = 252
				if (temp(i,j) == -3) nise_south(i,j) = 253
				if (temp(i,j) == -2) nise_south(i,j) = 254
				if (temp(i,j) == -1) nise_south(i,j) = 255
			end do
		end do
	  
!     Call routine to fill in snow/ice info for coastal regions as 
!     much as possible (southern hemisphere).
      call massage_snowice(nise_south,edge(1),edge(2))

! ... Close southern hemisphere grid

!      status = gddetach( grid_id )
	status = sfendacc(grid_id)

!-----------------------------------------------------------------------
!     CLOSE FILE
!-----------------------------------------------------------------------

! ... Close file and return success flag

!      status = gdclose( file_id )
	status = sfend(file_id)
      read_nise = 0



      END FUNCTION READ_NISE

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------


      subroutine massage_snowice(map_nise,xsize,ysize)

!---------------------------------------------------------------------
!!F77
!
!!Description:
!     Subroutine for filling in snow and ice information in coastal
!     regions.  For all coast pixels in the NISE data, search both
!     latitudinally and longitudinally for non-coast pixels which
!     indicate snow or ice.  Stop search when non-snow/ice, non-coast
!     pixels are found.  For those coast pixels which have snow/ice
!     pixels indicated on their boundaries, set flag to indicate 
!     snow/ice.
!
!!Input Parameters:
! map_nise    Input/output snow/ice data grid
! xsize       Grid size in the horizontal direction
! ysize       Grid size in the vertical direction
!
!!Output Parameters:
! map_nise    Input/output snow/ice data grid
!
!!Revision History:
! Original subroutine:
!     03-09-04  R. Frey
!
!!Team-Unique Header:
!
!!References and Credits:
! See Cloud Mask ATBD-MOD-06.
!
!!End
!----------------------------------------------------------------------

      implicit none
      save

      !Arguments.
      integer xsize, ysize
      integer*2 map_nise(xsize,ysize)

      !Local arrays.
      integer k(4)

      !Local scalars.
      integer temp_nise, i, j, m, n

      !Loop through map data.
 
      do j = 1, xsize

        do i = 1, ysize

          temp_nise = 0

         !Check for coast pixel (value = 252).


          if( map_nise(j,i) .eq. 252 ) then

           !Check if on map edges.  Cannot search on all four sides.
           !No data in corners.

            if(j .eq. 1) then
              k(1) = map_nise(j+1,i)
              k(2) = map_nise(j,i-1)
              k(3) = map_nise(j,i+1)
              k(4) = 0
         
            else if(i .eq. 1) then
              k(1) = map_nise(j,i+1)
              k(2) = map_nise(j+1,i)
              k(3) = map_nise(j-1,i)
              k(4) = 0

            else if(j .eq. xsize) then
              k(1) = map_nise(j-1,i)
              k(2) = map_nise(j,i-1)
              k(3) = map_nise(j,i+1)
              k(4) = 0
         
            else if(i .eq. ysize) then
              k(1) = map_nise(j,i-1)
              k(2) = map_nise(j+1,i)
              k(3) = map_nise(j-1,i)
              k(4) = 0

            else

             !Get all four adjacent values.
              k(1) = map_nise(j,i-1)
              k(2) = map_nise(j,i+1)
              k(3) = map_nise(j-1,i)
              k(4) = map_nise(j+1,i)

            end if

           !Fill in missing snow/ice values (=200).
            do m = 1,4

              n = k(m)
              if(n .ne. 0) then

                if((n .eq. 103 .or. n .eq. 104 .or. n .eq. 200) .or. &
                                  (n .gt. 25 .and. n .lt. 102)) then
                  temp_nise = 200
                end if

              end if
            enddo


            if (temp_nise .eq. 200) map_nise(j,i) = temp_nise

          end if

		  
        enddo

      enddo

      return
      end subroutine massage_snowice

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------


end module MOD06AlbedoEcoModule
