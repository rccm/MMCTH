      SUBROUTINE GET_CORE_METADATA( MODFIL_GEO, MODFIL_LUN, MCF_LUN,
     &  EAST_LIM, WEST_LIM, NORTH_LIM, SOUTH_LIM,
     &  BEGIN_DATE, BEGIN_TIME, END_DATE, END_TIME ) 

C-----------------------------------------------------------------------
C !F77
C
C !PURPOSE:
C    Extract Core Metadata fields from the MODIS Geolocation product.
C
C !DESCRIPTION:
C    The extracts information from both the the Geolocation file and
C    the the MCF file. The output fields represent the areal and time
C    bounds of the input granule.
C
C !INPUT PARAMETERS:
C    MODFIL_GEO      Geolocation file information array containing
C                    HDF SD_ID (index 1)
C                    VS file identifier (index 2)
C                    File access mode (index 3)
C    MODFIL_LUN      Variable containing the PCF logical
C                    reference number (LRN) to the Geolocation file.
C    MCF_LUN         Variable containing the PCF logical
C                    reference number (LRN) to the metadata
C                    configuration file (MCF)
C
C !OUTPUT PARAMETERS:
C    EAST_LIM        East longitude limit of given granule
C    WEST_LIM        West longitude limit of given granule
C    NORTH_LIM       North latitude limit of given granule
C    SOUTH_LIM       South latitude limit of given granule
C    BEGIN_DATE      Beginning date of granule
C    BEGIN_TIME      Beginning time of granule
C    END_DATE        Ending date of granule
C    END_TIME        Ending time of granule
C
C !REVISION HISTORY
C
C !TEAM-UNIQUE HEADER:
C    Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
C
C !END
C----------------------------------------------------------------------

      IMPLICIT NONE

      INCLUDE 'debug.inc'
      INCLUDE 'PGS_SMF.f'
      INCLUDE 'mapi.inc'
      INCLUDE 'hdf.inc'
      
c ... Arguments 

      INTEGER modfil_geo(modfillen), modfil_lun, mcf_lun
      DOUBLE PRECISION east_lim, west_lim, north_lim, south_lim
      CHARACTER*(*) begin_date, begin_time, end_date, end_time

c ... Local variables

      CHARACTER*36 attr, field_name(4)
      CHARACTER*49 master_grp_hndls(pgsd_met_num_of_groups)
      CHARACTER*70 buf_char(4)

      INTEGER i, rtn, version

      DOUBLE PRECISION buf_dbl(4)

c ... External functions

      INTEGER pgs_met_init, pgs_met_getpcattr_d, pgs_met_getpcattr_s
      EXTERNAL pgs_met_init, pgs_met_getpcattr_d, pgs_met_getpcattr_s
      
c ... Initialize holding arrays

      do i = 1 , 4
         buf_char(i) = '  '
         buf_dbl(i) = 0.0D0
      end do

      do i = 1 , PGSd_MET_NUM_OF_GROUPS
         Master_Grp_Hndls(i) = '  '
      end do

c ... Check inputs

      if ( modfil_Geo(P_SDID) .le. 0 ) then
        call message( 'Get_Core_Metadata',
     &  'Failed to get Geo Metadata [OPERATOR ACTION:' //
     &  ' Check input Geolocation file]', 0, 2 )
      else if ( modfil_Geo(P_ACCESS) .ne. DFACC_READ ) then
        call message( 'Get_Core_Metadata',
     &  'Failed to read Geo Metadata file [OPERATOR ACTION:' //
     &  ' Check Geolocation file read permissions]', 0, 2 )
      endif

c ... Initialize metadata tool defining Master_Grp_Hndls, and set 
c ... CoreMetadata.0 attribute fields that have values provided in MCF.

      rtn = pgs_met_init( MCF_lun, Master_Grp_Hndls )
      if ( rtn .ne. PGS_S_SUCCESS ) then
        call message( 'Get_Core_Metadata',
     &  'Failed to initialize Geo Metadata file [OPERATOR ACTION:' //
     &  ' Check MCF file]', 0, 2 )
      endif

c-----------------------------------------------------------------------
c     First call accesses the matadata parameters 
c     "EastBoundingCoordinate","NorthBoundingCoordinate", 
c     "SouthBoundingCoordinate" and "WestBoundingCoordinate"
c     in "mapi.inc" and sets them to the values retrieved 
c     from Geo product.
c     Note that the use of MAPI parameters are defined in "mapi.inc".
c-----------------------------------------------------------------------

      Field_Name(1) = MCORE_EAST_BOUND
      Field_Name(2) = MCORE_WEST_BOUND
      Field_Name(3) = MCORE_NORTH_BOUND
      Field_Name(4) = MCORE_SOUTH_BOUND
      
      do i = 1, 4

         attr = Field_Name(i)
         version = 1
         rtn = pgs_met_getpcattr_d( Modfil_lun, version,
     &     'ArchiveMetadata.0', attr, buf_dbl(i) )
     
         if ( rtn .ne. PGS_S_SUCCESS ) then
           call message( 'Get_Core_Metadata',
     &      'Failed to extract metadata from Geolocation file' //
     &      ' [OPERATOR ACTION: Make sure correct hdf file is loaded' //
     &      ' and correct MCF file is being use.]', 0, 2 )
         endif

      end do

c ... Place extracted values into correct variable names

      East_Lim  = buf_dbl(1)
      West_Lim  = buf_dbl(2)
      North_Lim = buf_dbl(3)
      South_Lim = buf_dbl(4)

c-----------------------------------------------------------------------
c     Second call accesses "RangeBeginningDateTime", 
c     "RangeEndingDateTime" in "mapi.inc" and
c     sets them to the values retrieved from Geo product. 
c-----------------------------------------------------------------------

      Field_Name(1) = 'RANGEBEGINNINGDATE'
      Field_Name(2) = 'RANGEBEGINNINGTIME'
      Field_Name(3) = 'RANGEENDINGDATE'
      Field_Name(4) = 'RANGEENDINGTIME'

      do i = 1, 4

         attr = Field_Name(i)
         version = 1
         rtn = pgs_met_getpcattr_s( Modfil_lun, version,
     &     MECS_CORE, attr, buf_char(i) )

         if ( rtn .ne. PGS_S_SUCCESS ) then
           call message( 'Get_Core_Metadata',
     &      'Failed to extract time bounds from L1B Metadata file' //
     &      ' [OPERATOR ACTION: Check input and MCF files]', 0, 2 )
         endif

      end do
     
c ... Place extracted values into correct variable names

      Begin_Date = buf_char(1)
      Begin_Time = buf_char(2)
      End_Date   = buf_char(3)
      End_Time   = buf_char(4)

c ... Write debug output if required

      if (debug .gt. 0) then

        write( h_output, '(10x,'' Areal Extent of Granule '')' )
        write( h_output,
     &    '(10x,''   East        West        North        South'','//
     &    '/,8x,4f12.6,/)')
     &    East_Lim, West_Lim, North_Lim, South_Lim
     
        write( h_output, '(10x,'' Temporal Bounds of Granule '')' )
        write( h_output,
     &     '(10x,'' Beginning Date & Time       Ending Date & Time'','//
     &     '/,5x,2A15,5x,2A15/)')
     &     Begin_Date, Begin_Time, End_Date, End_Time

      endif

      END
