module Map_to_CMGrid
  ! Module revision history
  ! $Log: Map_to_CMG.f90,v $
  ! Revision 1.3  1998/06/30  11:56:38  rhucek
  ! parameter statements reordered by SCF
  !
! Revision 1.8  1998/06/29  20:42:25  pincus
! Reorded parameter statements.
!
! Revision 1.7  1998/06/23  18:08:22  pincus
! Added prologs to individual subroutines to meet ECS requirements. Some cosmetic
! cleanup as well.
!
  ! Revision 1.6  1997/11/11  22:47:47  pincus
  ! Removed warning messages for out of bounds geolocation pairs
  ! as per Rich Hucek's suggestion. We'll report those errors higher
  ! up the food chain.
  !
  ! Revision 1.5  1997/11/10  20:36:01  pincus
  ! Removed tabs for ANSI compliance.
  !
  ! Revision 1.4  1997/07/08  17:28:06  pincus
  ! Added FatalErrorMessage to the list of entities from ErrorHandling that
  ! are USEd by Map_to_CMGrid.
  !
  ! Revision 1.3  1997/07/03  22:31:48  pincus
  ! Changed calls to ErrorMessage to reflect changes in the ErrorHandling module.
  ! Calls are now made to either FatalErrorMessage or WarningMessage.
  ! Also changed a few typos in the messages.
  !
  ! Revision 1.2  1997/06/24  20:01:31  pincus
  ! Changed naming convention to refer to latBin and lonBin
  ! instead of row and column. The old convention was causing some
  ! confusion because Modis arrays are set up with longitude varying
  ! fastest.
  !
  ! Revision 1.1  1997/06/10  18:05:47  pincus
  ! Initial revision
  !
  use typeSizes
  use ErrorHandling, only: FatalErrorMessage
  IMPLICIT NONE
  private

  ! A tile as used here is defined as a latitudinal band, or ring, 
  !   around the globe.  Tiles are numbered consecutively beginning 
  !   with 1 at the North Pole. 
  
  ! Parameters
  ! DEG_PER_TILE specifies tile width in degrees of latitude. 
  integer, parameter, public :: DEG_PER_TILE = 5, IN_BOUNDS = 0, OUT_OF_BOUNDS = -1
  ! GRID_RESOL, sets the CMG grid projection resolution.  
  !  Possible GRID_RESOL values for the AM1 program are 0.25, 0.5 and 
  !  1.0 degree latitude and longitude.
  real,    parameter, public :: GRID_RESOL     = 1.0
  real,    parameter         ::  XGRID_RESOL = 1.0/GRID_RESOL,            &
                                   MaxLonBin = 360.0 * XGRID_RESOL + 0.1, &
                                   MaxLatBin = 180.0 * XGRID_RESOL + 0.1, &
                                  LON_IS_360 = 360.0, COLAT_IS_180   = 180.0
  integer, parameter, public :: LatBinsPerTile = DEG_PER_TILE * XGRID_RESOL + 0.1, &
                                LonBinsPerTile = 360.0        * XGRID_RESOL + 0.1
                                
  ! Procedures
  interface Map_to_CMG
    module procedure Map_to_CMG_scalar, Map_to_CMG_vector, Map_to_CMG_matrix
  end interface
  public :: Map_to_CMG

CONTAINS
  !----------------------------------------------------------------------
  ! Three versions of the code - scalar, and for rank 1 and rank two arrays. 
  !----------------------------------------------------------------------
  SUBROUTINE Map_to_CMG_scalar (lat, lon, TileNum, latBin, lonBin, status)
    IMPLICIT NONE
    real,                            intent( in) :: lat, lon
    integer (kind = TwoByteIntKind), intent(out) :: TileNum, latBin, lonBin, status
    ! !F90
    !
    ! !DESCRIPTION:  
    !  Identify the tile number containing the geographic location of a latitude 
    !    and longitude coordinate pair, and the tile row (latitude) and column (longitude) numbers 
    !    of the observation on the Earth Observing System equal-angle Climate Modeling Grid 
    !    (CMG) projections.  
    !
    ! !INPUT PARAMETERS:
    !    lat: Latitude in degrees (+90 to -90)
    !    lon: Longitude in degrees (-180 to +180)
    !
    ! !OUTPUT PARAMETERS:
    !    status: Indicates whether input "lat" and "lon"  represent a valid 
    !      geographic location. Possible values are 0 (in bounds) and -1 (out of bounds) 
    !    TileNum: Tile number containing lat/lon observation. 
    !      If either "lat" or "lon" is invalid, TileNum set to -1. 
    !    latBin: Relative tile row (latitude) number of observation. 
    !       If either "lat" or "lon" is invalid, latBin is set to -1. 
    !    lonBin: Relative tile column (longitude) number of observation.
    !       If either "lat" or "lon" is invalid, lonBin is set to -1. 
    !
    ! !REVISION HISTORY:
    !    See module revision history at top of file. 
    !
    ! !TEAM-UNIQUE HEADER:
    !    Cloud Retrieval Group, NASA Goddard Space Flight Center
    !
    ! !REFERENCES AND CREDITS:
    !    Written by 
    !      Robert Pincus
    !      Climate and Radiation Branch, Code 913
    !      NASA Goddard Space Flight Center
    !      Greenbelt MD 20771
    !      Robert.Pincus@gsfc.nasa.gov
    !
    ! !DESIGN NOTES
    !   Derived from code written by Richard Hucek.
    !
    ! !END
 
    !...  Declarations
    REAL :: colat, lonP180

    !...  compute colatitude, longitude-plus-180
    colat = 90.0 - lat
    lonP180 = 180.0 + lon

    If (colat   < -tiny(colat)                         .or. &
        colat   > COLAT_IS_180 + spacing(COLAT_IS_180) .or. &
        lonP180 < -tiny(lonP180)                       .or. &
        lonP180 > LON_IS_360 + spacing(LON_IS_360) ) Then 
      TileNum = OUT_OF_BOUNDS
      latBin  = OUT_OF_BOUNDS
      lonBin  = OUT_OF_BOUNDS
      status  = OUT_OF_BOUNDS 
    else
      !...  colat and longitude-plus-180 are within epsilon of the ranges 
      !     0-180 and 0-360, respectively.  Safe to compute TileNum, and 
      !     tile latBin and column numbers on CMG equal-angle grids 
      !     (0.25, 0.5 and 1 degree). 
      latBin = int(colat * XGRID_RESOL) + 1
      If (latBin < 1)         latBin = 1 
      If (latBin > MaxLatBin) latBin = MaxLatBin 

      TileNum = (latBin-1)/LatBinsPerTile + 1
      latBin  = latBin - (TileNum - 1)*LatBinsPerTile

      lonBin  = int(lonP180 * XGRID_RESOL) + 1
      If ( (lonBin < 1) .or. (lonBin > MaxLonBin) ) lonBin = 1 
  
      status = IN_BOUNDS 
    end if
  End subroutine Map_to_CMG_scalar
  ! ----------------------------------------------------------------------
  SUBROUTINE Map_to_CMG_vector (lat, lon, TileNum, latBin, lonBin, status)
    IMPLICIT NONE
    real,    dimension(:), intent( in) :: lat, lon
    integer (kind = TwoByteIntKind), &
             dimension(:), intent(out) :: TileNum, latBin, lonBin, status
    ! !F90
    !
    ! !DESCRIPTION:  
    !  Identify the tile number containing the geographic location of a latitude 
    !    and longitude coordinate pair, and the tile row (latitude) and column (longitude) numbers 
    !    of the observation on the Earth Observing System equal-angle Climate Modeling Grid 
    !    (CMG) projections.  
    !
    !  Note that all input and output vectors passed to this subroutine must be the same size. 
    !
    ! !INPUT PARAMETERS:
    !    lat: Latitude in degrees (+90 to -90)
    !    lon: Longitude in degrees (-180 to +180)
    !
    ! !OUTPUT PARAMETERS:
    !    status: Indicates whether input "lat" and "lon"  represent a valid 
    !      geographic location. Possible values are 0 (in bounds) and -1 (out of bounds) 
    !    TileNum: Tile number containing lat/lon observation. 
    !      If either "lat" or "lon" is invalid, TileNum set to -1. 
    !    latBin: Relative tile row (latitude) number of observation. 
    !       If either "lat" or "lon" is invalid, latBin is set to -1. 
    !    lonBin: Relative tile column (longitude) number of observation.
    !       If either "lat" or "lon" is invalid, lonBin is set to -1. 
    !
    ! !REVISION HISTORY:
    !    See module revision history at top of file. 
    !
    ! !TEAM-UNIQUE HEADER:
    !    Cloud Retrieval Group, NASA Goddard Space Flight Center
    !
    ! !REFERENCES AND CREDITS:
    !    Written by 
    !      Robert Pincus
    !      Climate and Radiation Branch, Code 913
    !      NASA Goddard Space Flight Center
    !      Greenbelt MD 20771
    !      Robert.Pincus@gsfc.nasa.gov
    !
    ! !DESIGN NOTES
    !   Derived from code written by Richard Hucek.
    !
    ! !END
    
    !...  Local variables
    integer                      :: i           ! Counter
    character (len = 2048)       :: msgbuffer

    ! Automatic arrays - allocated on entry, deallocated on exit.
    REAL,    dimension(size(lat, 1)) :: colat, lonP180
    logical, dimension(size(lat, 1)) :: latInBounds, lonInBounds

    !...  Check that all the vectors have the same length

    if(.not. all( (/ size(lat, 1), size(lon, 1), size(TileNum, 1), &
                     size(latBin, 1), size(lonBin, 1) /) == size(status, 1) )) then
      msgbuffer = "Argument vectors are not all the same length." &
                   // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("Map_to_CMG_vector", msgbuffer)
    end if

    !...  compute colatitude,  longitude-plus-180, mark out-of-bounds elements
    colat(:) = 90.0 - lat(:)
    latInBounds(:) = (colat(:) >= 0) .and. (colat(:) <= COLAT_IS_180)

    lonP180(:) = 180.0 + lon(:)
    lonInBounds(:) = (lonP180(:) >= 0) .and. (lonP180(:) <= LON_IS_360)

    !...  colat and longitude-plus-180 are within epsilon of the ranges 
    !     0-180 and 0-360, respectively.  Safe to compute TileNum, and 
    !     tile latBin and column numbers on CMG equal-angle grids 
    !     (0.25, 0.5 and 1 degree). 

    latBin(:) = int(colat(:) * XGRID_RESOL) + 1
    where (latBin < 1)         latBin = 1 
    where (latBin > MaxLatBin) latBin = MaxLatBin 

    TileNum(:) = (latBin(:) - 1)/LatBinsPerTile + 1
    latBin(:)  =  latBin(:) - (TileNum - 1) * LatBinsPerTile

    lonBin(:)  = int(lonP180 * XGRID_RESOL) + 1
    where ( (lonBin < 1) .or. (lonBin > MaxLonBin) ) lonBin = 1 

! Insert bad value flags where appropriate
    where ( (.not. lonInBounds(:)) .or. (.not. latInBounds(:)) )
      TileNum(:) = OUT_OF_BOUNDS
      latBin(:)  = OUT_OF_BOUNDS
      lonBin(:)  = OUT_OF_BOUNDS
      status(:)  = OUT_OF_BOUNDS 
    elsewhere 
      status(:)  = IN_BOUNDS 
    end where
  End subroutine Map_to_CMG_vector
  ! ----------------------------------------------------------------------
  SUBROUTINE Map_to_CMG_matrix (lat, lon, TileNum, latBin, lonBin, status)
    IMPLICIT NONE
    real,    dimension(:,:), intent( in) :: lat, lon
    integer (kind = TwoByteIntKind), &
             dimension(:,:), intent(out) :: TileNum, latBin,lonBin,status
    ! !F90
    !
    ! !DESCRIPTION:  
    !  Identify the tile number containing the geographic location of a latitude 
    !    and longitude coordinate pair, and the tile row (latitude) and column (longitude) numbers 
    !    of the observation on the Earth Observing System equal-angle Climate Modeling Grid 
    !    (CMG) projections.  
    !
    !  Note that all input and output arrays passed to this subroutine must be the same size. 
    !
    ! !INPUT PARAMETERS:
    !    lat: Latitude in degrees (+90 to -90)
    !    lon: Longitude in degrees (-180 to +180)
    !
    ! !OUTPUT PARAMETERS:
    !    status: Indicates whether input "lat" and "lon"  represent a valid 
    !      geographic location. Possible values are 0 (in bounds) and -1 (out of bounds) 
    !    TileNum: Tile number containing lat/lon observation. 
    !      If either "lat" or "lon" is invalid, TileNum set to -1. 
    !    latBin: Relative tile row (latitude) number of observation. 
    !       If either "lat" or "lon" is invalid, latBin is set to -1. 
    !    lonBin: Relative tile column (longitude) number of observation.
    !       If either "lat" or "lon" is invalid, lonBin is set to -1. 
    !
    ! !REVISION HISTORY:
    !    See module revision history at top of file. 
    !
    ! !TEAM-UNIQUE HEADER:
    !    Cloud Retrieval Group, NASA Goddard Space Flight Center
    !
    ! !REFERENCES AND CREDITS:
    !    Written by 
    !      Robert Pincus
    !      Climate and Radiation Branch, Code 913
    !      NASA Goddard Space Flight Center
    !      Greenbelt MD 20771
    !      Robert.Pincus@gsfc.nasa.gov
    !
    ! !DESIGN NOTES
    !   Derived from code written by Richard Hucek.
    !
    ! !END

    !...  Local variables
    integer              :: i     ! Counter

    ! Automatic arrays - allocated on entry, deallocated on exit.

    real,    dimension(size(lat, 1), size(lat, 2)) :: colat, lonP180
    real,    dimension(size(lat))                  :: badValues
    logical, dimension(size(lat, 1), size(lat, 2)) :: latInBounds, lonInBounds
    character (len = 2048)                         :: msgbuffer

    !...  Check that all the matrices have the same size

    If(.not. all( (/ size(lat, 1),    size(lon, 1), size(TileNum, 1),          &
                     size(latBin, 1), size(lonBin, 1) /) == size(status, 1) )  &
       .or.                                                                    &
       .not. all( (/ size(lat, 2),    size(lon, 2), size(TileNum, 2),          &
                     size(latBin, 2), size(lonBin, 2) /) == size(status, 2) ) ) then
      msgbuffer = "Argument matrices are not all the same length." &
                  // char(10) // "Operator Action:  Notify SDST."
      call FatalErrorMessage("Map_to_CMG_matrix", msgbuffer)
    End If

    !...  compute colatitude,  longitude-plus-180, mark out-of-bounds elements
    colat(:, :) = 90.0 - lat(:, :)
    latInBounds(:, :) = (colat(:, :) >= 0) .and. (colat(:, :) <= COLAT_IS_180)

    lonP180(:, :) = 180.0 + lon(:, :)
    lonInBounds(:, :) = (lonP180(:, :) >= 0) .and. (lonP180(:, :) <= LON_IS_360)

    !...  colat and longitude-plus-180 are within epsilon of the ranges 
    !     0-180 and 0-360, respectively.  Safe to compute TileNum, and 
    !     tile latBin and column numbers on CMG equal-angle grids 
    !     (0.25, 0.5 and 1 degree). 

    latBin(:, :) = int(colat(:, :) * XGRID_RESOL) + 1
    where (latBin < 1)         latBin = 1 
    where (latBin > MaxLatBin) latBin = MaxLatBin 

    TileNum(:, :) = (latBin(:, :) - 1)/LatBinsPerTile + 1
    latBin(:, :)  =  latBin(:, :) - (TileNum - 1) * LatBinsPerTile

    lonBin(:, :)  = int(lonP180 * XGRID_RESOL) + 1
    where ( (lonBin < 1) .or. (lonBin > MaxLonBin) ) lonBin = 1 

    ! Insert bad value flags where appropriate
    where ( (.not. lonInBounds(:, :)) .or. (.not. latInBounds(:, :)) )
      TileNum(:, :) = OUT_OF_BOUNDS
      latBin(:, :)  = OUT_OF_BOUNDS
      lonBin(:, :)  = OUT_OF_BOUNDS
      status(:, :)  = OUT_OF_BOUNDS 
    elsewhere 
      status(:, :)  = IN_BOUNDS 
    end where
  End subroutine Map_to_CMG_matrix
  ! ----------------------------------------------------------------------
end module Map_to_CMGrid
