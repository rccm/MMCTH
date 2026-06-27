             PROGRAM MOD_PR04CR

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION: This program builds an HDF file from a file
C               specification document written in the tiny language
C               known as CDL (network Common data form Description
C               Language). The HDF file produced contains the
C               definitions of all data structures listed in the product
C               specification document (except for ECS metadata), but
C               none of its arrays or tables objects are populated at
C               this time.
C
C
C !INPUT PARAMETERS:    N/A
C
C
C
C !OUTPUT PARAMETERS:   N/A
C
C
C
C !REVISION HISTORY:
C  Initial Version: rhucek 2006/12/12
C
C !TEAM-UNIQUE HEADER:
C
C    This software is developed by the MODIS Science Data Support Team
C    for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C    HDF portions developed at the National Center for Supercomputing
C    Applications at the University of Illinois at Urbana-Champaign.
C
C !REFERENCES AND CREDITS:
C
C    Written by         Richard Hucek
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd,  Seabrook MD 20706
C
C !DESIGN NOTES:
C
C  Externals:
C
C    Functions and Subroutines:
C        mod_create_hdf            (science code)
C        ehidinfo                  (hdfeos lib)
C        get_spec_info             (mod04_get_spec_info.f)
C        mod04_GetL1B_Info         (mod04_GetL1B_Info.f)
C        pgs_pc_getreference       (sdptk lib)
C        message                   (src_UW)
C        SetGlobalAttr             (mod04_SetGlobalAttr.f)
C        swattach                  (hdfeos lib)
C        swclose                   (hdfeos lib)
C        swcreate                  (hdfeos lib)
C        swdefdim                  (hdfeos lib)
C        swdetach                  (hdfeos lib)
C        swopen                    (hdfeos lib)
C        swwrfld                   (hdfeos lib)
C
C
C    Named Constant:
C        PGSd_PC_FILE_PATH_MAX     (PGS_PC.f)
C        PGS_S_SUCCESS             (PGS_SMF.f)
C
C !END
C----------------------------------------------------------------------

      implicit none
      save

      include 'mapi.inc'
      include 'hdf.inc'
      include 'dffunc.inc'
      include 'PGS_PC.f'
      include 'PGS_SMF.f'
      include 'mod04_create2.inc'
      include 'mod04_pcfnum2.inc'

c   Symbolic constants for HDF-EOS
      integer        HDFE_NOMERGE
      parameter    ( HDFE_NOMERGE = 0 )

      character*32   ROUTINE_NAME
      parameter    ( ROUTINE_NAME = 'MOD_PR04CR' )

      integer        LIMIT
      parameter    ( LIMIT = 1000)


c   Character data
      character*80   ProgramName, FieldNames(100), DimNames(100)
      character*80   DataNames(100)
      character*256  attrname, title
      character*256  swathname
      character*2048 descrip
      character*( PGSd_PC_FILE_PATH_MAX ) productfilename

c   Functions defined in this module
      integer  SetGlobalAttr, create_hdf, create_vdata,
     &         pgs_pc_getreference
      external SetGlobalAttr, create_hdf, create_vdata,
     &         pgs_pc_getreference

c   Variables for HDF-EOS routine
      integer swopen,swcreate,swdefdim
      integer swdetach,swattach,swclose,swwrfld
      integer ehidinfo
      integer fid,swathid, srtn, merge, sdid

      integer nscans,    npixels
      integer lines1km,  elements1km
      integer lines10km, elements10km

      integer i, j, fileversion, rtn
      integer FieldStart, FieldCount, NumOfData
      integer DimCount, DataCount(LIMIT), DimNumbers(LIMIT)
      integer LongValues(LIMIT)

      real    DataValues(LIMIT, LIMIT)
      real    RealValues(LIMIT)

      integer*2 IntValues(LIMIT)

c   Variables for writing vdata to the HDF file
      integer vfid, vdata_id
      integer VDataStart, VDataCount

c   Variables for locally developed code.
      integer count, Data(20), numbertype, sl, strlen

      integer DataType(20) /7*2, 13*0/

c   Variable DataType maps an integer value to a Fortran data type
c   The association is:
c
c      Integer Value     Fortran Data Type
c      -------------     -----------------
c           1                Integer*2
c           2                Integer
c           3                Real

c     initialize variables
      sdid        = 0
      ProgramName = "MOD04_L2"


c=======================================================================
c
c The procedures for the software development of this code follow the
c method prescribed for writing HDF-EOS swath products as presented in
c the technical paper entitled:
c
c                 Writing HDF-EOS Swath Products for
c                     Optimum Subsetting Services
c
c                           170-TP-009-001
c
c The following development steps are observed.
c
c  1 - open the file using function swopen and create a swath using
c      swcreate
c
c  2 - Set global metadata
c
c  3 - define all dimensions needed for the data and geolocation
c      variables using swdefdim
c
c  4 - define each geolocation field in the swath using swdefgfld.
c      To enable subsetting by region, two geolocation fields must exist
c      in the swath and be named "Longitude" and either "Latitude" or
c      "Colatitude".
c
c  5 - define each data field in the swath using swdefdfld.
c
c  6 - create dimension map for the swath fields using swdefmap. If the
c      geolocation and data fields shared dimensions, no mapping is
c      needed.
c
c  7 - detach from the swath using swdetach.  This fixes the dimension
c      maps and data fields before writing data
c
c  8 - re-attach to the swath using swattach
c
c  9 - write the science and geolocation data into the proper arrays
c      using swwrfld
c
c 10 - detach from the swath using swdetach and close the file using
c      swclose
c
c=======================================================================
c ... Get output filename from process control file
c
c
c=======================================================================
      fileversion    = 1
      productfilename = ' '
      rtn             = PGS_S_SUCCESS - 1
      rtn             = pgs_pc_getreference( out_pcfnum, fileversion,
     &                                       productfilename )

c=======================================================================
c Item 1
c
c Open product file and create a swath
c
c=======================================================================
c ... Open the output swath file.
      fid=swopen(productfilename,DFACC_CREATE)

      if (fid .eq. -1) then
         call message( 'MOD_PR04CR',
     $   'Error creating L2 HDF product file' //
     $   '[OPERATOR ACTION: Check that file creation is possible]',0,2)
      endif

c ... Create a swath.
      swathname = "mod04"
      swathid   = swcreate(fid,swathname)

c ... Get the SDS interface ID and HDF file ID for the created swath.
      srtn=ehidinfo(fid,vfid,sdid)

      if (srtn .ne. 0) then
         call message ('MOD_PR04CR',
     $   'Error getting output file id-s using ehidinfo ' //
     $   '[OPERATOR ACTION: Notify SDST.]',0,2)
      endif



c=======================================================================
c Item 2
c
c Get the number of 1km MODIS scans and pixels from the L1B file. Then
c set global attributes
c
c=======================================================================
      call mod04_GetL1B_Info( nscans, npixels )
       
      lines1km     = (nscans *10) / lines_per_scan_1km
      lines10km    = lines1km / isamp
      elements1km  = npixels
      elements10km = elements1km / isamp
      
c ... Set global attributes
c      srtn = SetGlobalAttr(sdid, lines10km, elements10km)


c=======================================================================
c Items 3 and 4
c
c Define dimensions of data and geolocation fields.  Call get_spec_data
c to retrieve all dimension names listed in MOD04 file specification.
c
c=======================================================================
      call get_spec_data(ProgramName, FieldNames, DimNames,
     &                   DimNumbers, DimCount, DataNames, DataCount,
     &                   DataValues, NumOfData, LIMIT, FieldCount)


      do i = 1,DimCount
          lines1km     = (nscans *10) / lines_per_scan_1km
c ...... dimensions "Cell_Across_Swath" and "Cell_Along_Swath" are
c        variable quantities that are set dynamically from information
c        retrieved from the L1B product - not from the filespec.

         if (DimNames(i) .eq. "Cell_Across_Swath") Then
            DimNumbers(i) = elements10km
         else If (DimNames(i) .eq. "Cell_Along_Swath") Then
            DimNumbers(i) = lines1km
         end If
         
         if (DimNames(i) .eq. "Cell_Across_Swath_500") Then
            DimNumbers(i) = elements10km*20
            
         else If (DimNames(i) .eq. "Cell_Along_Swath_500") Then
            DimNumbers(i) =lines10km*20
             
         end If
         
    
         if (DimNumbers(i) .eq. 0) Then
            srtn=swdefdim(swathid,DimNames(i),SD_UNLIMITED)
         else
            srtn=swdefdim(swathid,DimNames(i),DimNumbers(i))
         end if

         if (srtn .ne. 0) then
            call message ('MOD_PR04CR',
     &      'Error defining dimensions for output file: swdefdim ' //
     &      '[OPERATOR ACTION: Notify SDST.]',0,2)
         endif

      end do


c=======================================================================
c Items 5 & 6 (no dimension mapping need)
c
c                Define Science Data and Geolocation Fields
c
c Define data fields one at a time, both geolocation and science data.
c Geolocation fields use function swdefgfld while data fields use
c swdefdfld.  Local attributes are also assigned to fields in this step.
c
c The previous call to get_spec_data provides the scalar variable
c FieldCount and the FieldNames array for use in the do loop.
c
c=======================================================================
      FieldStart = 1
       do i = FieldStart, FieldCount
      srtn = SetGlobalAttr(sdid, lines10km, elements10km)
      lines1km     = (nscans * 10)/lines_per_scan_1km
      lines10km    = lines1km / isamp
      elements1km  = npixels
      elements10km = elements1km / isamp   
       
      if( FieldNames(i) .eq.  "Aerosol_Cldmask_Land_Ocean" .or.
     &   FieldNames(i) .eq.   "Cloud_Distance_Land_Ocean") then 
         lines1km=lines1km *2
         lines10km=lines10km*20
         elements10km=  elements10km*20
         srtn = SetGlobalAttr(sdid, lines10km, elements10km)
         endif
      
         merge=HDFE_NOMERGE
       
         srtn=create_hdf(vfid,sdid,swathid,lines1km,
     &                   ProgramName, FieldNames(i),merge)
          
         if (srtn .ne. 0) then
           call message ('MOD_PR04CR',
     &     'Error creating output hdf file: create_hdf ' //
     &     '[OPERATOR ACTION: Check system resources.]',0,2)
         endif
      end do


c=======================================================================
c Items 7 & 8
c
c "Fix" the data fields (geolocation and science) by detaching and
c re-attaching to the same swath prior to writing any data.
c
c Note that the MOD04 application shares dimension fields between the
c geolocation and data fields and this obviates the need to provide
c explicit dimension mapping.
c
c=======================================================================
c ... Detach from the swath.
      srtn=swdetach(swathid)

      if (srtn .ne. 0) then
         call message ('MOD_PR04CR',
     &   'Error detaching from swath: swdetach ' //
     &   '[OPERATOR ACTION: Notify SDST.]',0,2)
      endif

c ... Attach to swath
      swathid=swattach(fid,swathname)


c=======================================================================
c Item 9
c
c Write integer data to each 1-dimensional SDS (vdata) field.  These
c data identify planes of 3-dimensional arrays and permit subsetting by
c level.
c
c=======================================================================
      do i=1,NumOfData

         do j = 1, LIMIT
            IntValues(j)  = 0
            LongValues(j) = 0
            RealValues(j) = 0
         end do

         if (DataType(i) .eq. 1) Then
            do j = 1, DataCount(i)
               IntValues(j) = INT(DataValues(i,j))
            end do

            vdata_id=swwrfld(swathid, DataNames(i), 0, 1, DataCount(i), IntValues)

            if (vdata_id .lt. 0) Then
               call message ('MOD_PR04CR',
     &           'Error writing datanames to file:vhfsd1 ' //
     &           '[OPERATOR ACTION: Notify SDST.]',0,2)
            endif

         else if (DataType(i) .eq. 2) Then

            do j = 1, DataCount(i)
               LongValues(j) = INT(DataValues(i,j))
            end do

            vdata_id=swwrfld(swathid, DataNames(i), 0, 1, DataCount(i), LongValues)

            if (vdata_id .lt. 0) Then
               call message ('MOD_PR04CR',
     &           'Error writing datanames to file:vhfsd2 ' //
     &           '[OPERATOR ACTION: Notify SDST.]',0,2)
            endif

         else if (DataType(i) .eq. 3) Then
            do j = 1, DataCount(i)
               RealValues(j) = Real(DataValues(i,j))
            end do

            vdata_id=swwrfld(swathid, DataNames(i), 0, 1, DataCount(i), RealValues)

            if (vdata_id .lt. 0) Then
               call message ('MOD_PR04CR',
     &           'Error writing datanames to file:vhfsd3 ' //
     &           '[OPERATOR ACTION: Notify SDST.]',0,2)
            endif

         end If
      end do


c=======================================================================
c Item 10
c
c Detach from the swath and close the MOD04 product file
c
c=======================================================================

c ... Detach from the swath
      srtn=swdetach(swathid)

      if (srtn .ne. 0) then
         call message ('MOD_PR04CR',
     &   'Error ending vdata writes to output file: vfend ' //
     &   '[OPERATOR ACTION: Notify SDST.]',0,2)
      endif


c ... Close the file
      srtn=swclose(fid)

      if (srtn .ne. 0) then
         call message ('MOD_PR04CR',
     &   'Error close output file: swclose ' //
     &   '[OPERATOR ACTION: Check System Resources.]',0,2)
      endif


      END
