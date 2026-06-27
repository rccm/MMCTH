      Subroutine mod06_Get_Core_Metadata(Modfil_Geo,Modfil_lun,MCF_lun,
     +                             East_Lim,West_Lim,North_Lim,
     +                             South_Lim,Begin_Date,Begin_Time,
     +                             End_Date,End_Time) 
      implicit none
      save

      include 'mapi.inc'
      include 'hdf.inc'
      include 'PGS_PC.f'
      include 'PGS_MODIS_39500.f'
      include 'PGS_SMF.f'
      include 'mod06uw_debug.inc'
      
c-----------------------------------------------------------------------
c!F77
c
c!Purpose:
c   To strip out the needed EOSDIS Core System inventory metada fields
c   needed by the UW science production software packages.  
c
c!Description:
c
c   The program uses information found in the MCF file as well as
c   the Geolocation file to extract 6 fields from ECS core metadata.
c   The fields represent the areal and time bounds of the given
c   granule. This program originally used the L1B file to strip
c   out values, but between May and July, a new L1B file was
c   created with changes that made this revision necessary.
c   
c
c!Input Parameters:
c  integer  Modfil_Geo      Array containing HDF SD (index 1) and VS 
c                           (index 2) file identifiers and L1B file 
c                           access mode (index 3) for Geolocation file.
c  integer  Modfil_lun      Variable containing the PCF logical 
c                           reference number (LRN) to the Geolocation
c                           file.
c  integer  MCF_lun         Variable containing the PCF logical 
c                           reference number (LRN) to the metadata 
c                           configuration file (MCF)
c!Output Parameters:
c
c  double   East_Lim        East longitude limit of given granule
c   precision
c
c  double   West_Lim        West longitude limit of given granule
c   precision
c
c  double   North_Lim       North latitude limit of given granule
c   precision
c
c  double   South_Lim       South latitude limit of given granule
c   precision
c
c  Character Begin_Date     Beginning date of granule 
c     *70
c 
c  Character Begin_Time     Beginning time of granule 
c     *70
c
c  Character End_Date       Ending date of granule 
c     *70
c 
c  Character End_Time       Ending time of granule 
c     *70
c
c !REVISION HISTORY:
c $Id: mod06_Get_Core_Metadata.f,v 1.7 1999/05/21 21:59:16 kis Exp $
c
c!Team-Unique Header:
c
c    This software is developed by the MODIS Science Data Support Team
c    for the National Aeronautics and Space Administration,
c    Goddard Space Flight Center, under contract NAS5-32373.
c
c HDF portions developed at the National Center for Supercomputing
c Applications at the University of Illinois at Urbana-Champaign.
c
c !REFERENCES AND CREDITS:
c
c!Design Notes:
c  
c  Externals:
c
c    Functions and Subroutines:
c        pgs_met_init          (libPGSTK.a)
c        pgs_met_setattr_*     (libPGSTK.a)
c        message  writes error message to LogStatus file
c
c    Named Constant:
c        DFACC_READ            (hdf.inc)
c        MAPIOK                (mapi.inc)
c        MCORE_*               (mapi.inc)
c        MECS_CORE             (mapi.inc)
c        MFAIL                 (mapi.inc)
c        MODFILLEN             (mapi.inc)
c        MODIS_E_GENERIC       (PGS_MODIS_39500.f)
c        MODIS_W_GENERIC       (PGS_MODIS_39500.f)
c        PGSd_MET_NUM_OF_GROUPS(PGS_MET.f)
c        PGSd_PC_FILE_PATH_MAX (PGS_PC.f)
c        PGS_S_SUCCESS         (PGS_SMF.f)
c        P_SDID                (mapi.inc)
c        P_ACCESS              (mapi.inc)
c    
c!END
c----------------------------------------------------------------------

c     Scalar Arguments 
      character*70 Begin_Date,End_Date,Begin_Time,End_Time
      integer Modfil_lun,MCF_lun
      double precision East_Lim,West_Lim,North_Lim,South_Lim

c     Local variables
      character*36  attr,Field_Name(4)
      character*49  Master_Grp_Hndls(PGSd_MET_NUM_OF_GROUPS)
      character*70  buf_char(4)
      integer    i,rtn,version,modfil_Geo(MODFILLEN)
      integer    pgs_met_init,pgs_met_getpcattr_d,pgs_met_getpcattr_s
      double precision buf_dbl(4)

c ... set program name for error messaging
      character*32 routine_name
      parameter ( routine_name = 'mod06_Get_Core_Metadata' )

C     Initialize holding arrays
      do 5 i = 1 , 4
         buf_char(i) = '  '
   5  continue

      do 10 i = 1, 4
         buf_dbl(i) = 0.0D0
   10 continue

      do 11 i = 1 , PGSd_MET_NUM_OF_GROUPS
         Master_Grp_Hndls(i) = '  '
   11 continue

c     Check inputs
      if (modfil_Geo(P_SDID).le.0) then
        call message( routine_name,
     &  'Failed to get Geo Metadata. [OPERATOR ACTION:' //
     &  ' Check input Geolocation file]',
     &  0, 2 )

      else if (modfil_Geo(P_ACCESS).ne.DFACC_READ) then
        call message( routine_name,
     &  'Failed to read Geo Metadata file. [OPERATOR ACTION:' //
     &  ' Check Geolocation file read permissions]',
     &  0, 2 )
      endif

C     Initialize metadata tool defining Master_Grp_Hndls, and set 
C     CoreMetadata.0 attribute fields that have values provided in MCF.
      rtn = pgs_met_init(MCF_lun, Master_Grp_Hndls)
      if (rtn.ne.PGS_S_SUCCESS) then
        call message( routine_name,
     &  'Failed to initialize Geo Metadata file [OPERATOR ACTION:' //
     &  ' Check MCF file]',
     &  0, 2 )
      endif

c----------------------------------------------------------------
c     Set CoreMetadata.0 attribute fields that have values retrieved 
c     from Geo product.  Note that the use of M-API parameters are 
c     defined in "mapi.inc".
c     First call accesses the matadata parameters 
c     "EastBoundingCoordinate","NorthBoundingCoordinate", 
c     "SouthBoundingCoordinate" and "WestBoundingCoordinate"
c     in "mapi.inc" and sets them to the values retrieved 
c     from Geo product.
c----------------------------------------------------------------

      Field_Name(1) = MCORE_EAST_BOUND
      Field_Name(2) = MCORE_WEST_BOUND
      Field_Name(3) = MCORE_NORTH_BOUND
      Field_Name(4) = MCORE_SOUTH_BOUND
      do 20 i = 1, 4
         attr = Field_Name(i)
         version = 1
         rtn = pgs_met_getpcattr_d(Modfil_lun,version,
     +                             'ArchiveMetadata.0',attr,buf_dbl(i))
         if (rtn.ne.PGS_S_SUCCESS) then
           call message( routine_name,
     &      'Failed to grab lat/lon bounds from Geo Metadata file' //
     &      ' [OPERATOR ACTION: check input and MCF files]',
     &      0, 2 )
         endif
   20 continue

c     Place extracted values into correct variable names
      East_Lim = buf_dbl(1)
      West_Lim = buf_dbl(2)
      North_Lim = buf_dbl(3)
      South_Lim = buf_dbl(4)

c----------------------------------------------------------------

c     The second call access  "RangeBeginningDateTime", 
c     "RangeEndingDateTime" in "mapi.inc" and
c     sets them to the values retrieved from Geo product. 
      Field_Name(1) = 'RANGEBEGINNINGDATE'
      Field_Name(2) = 'RANGEBEGINNINGTIME'
      Field_Name(3) = 'RANGEENDINGDATE'
      Field_Name(4) = 'RANGEENDINGTIME'
      do 30 i = 1, 4
         attr = Field_Name(i)
         version = 1
         rtn = pgs_met_getpcattr_s(Modfil_lun,version,
     +                             MECS_CORE,attr,buf_char(i))
         if (rtn.ne.PGS_S_SUCCESS) then
           call message( routine_name,
     &      'Failed to grab time bounds from L1B Metadata file' //
     &      ' [OPERATOR ACTION: check input and MCF files]',
     &      0, 2 )
         endif
   30 continue

      Begin_Date = buf_char(1)
      Begin_Time = buf_char(2)
      End_Date = buf_char(3)
      End_Time = buf_char(4)

c----------------------------------------------------------------
c ....................................................................
      if (debug .gt. 0) then
       write( h_output, '(72(''-''))' )
       write( h_output, '(a,/)' ) 'mod06_Get_Core_Metadata debug'
       write(h_output,'(10x,'' Areal Extent of Granule '')')
       write(h_output,'(10x,''   East        West        North        So
     +uth'',
     +      /,8x,4f12.6,/)') East_Lim, West_Lim, North_Lim, South_Lim
       write(h_output,'(10x,'' Temporal Bounds of Granule '')')
       write(h_output,'(10x,'' Beginning Date & Time        Ending Date 
     +& Time'',
     +          /,5x,2A15,5x,2A15/)') Begin_Date,Begin_Time,
     +                    End_Date,End_Time
      endif
c ....................................................................


      return
      end
