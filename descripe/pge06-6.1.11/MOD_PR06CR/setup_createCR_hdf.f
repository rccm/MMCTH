      SUBROUTINE SETUP_CREATE_HDF (PROG_NUMBER, PROG_NAME,
     $                                   LINES1KM, ELEMENTS1KM,
     $                                   LINES5KM, ELEMENTS5KM,
     $                                   MODFIL)

c---------------------------------------------------------------------
C!F77
c
c!Description:
c     Routine for creating a UW MODIS HDF output file for the
c     given swath based upon the L1B file and the product
c     (MOD35,MOD07,etc.) file specification.  At the end of the
c     create routine, which uses mostly EOS hdf calls, the
c     output file is opened using MAPI commands.
c
c!Input parameters:
c Prog_Number   EOS MODIS program number.
c Prog_Name     Output file name - extracted from pcf file
c lines1km      Maximum number of lines in file (from L1B file)
c elements1km   Maximum number of elements in file (from L1B file)
c lines5km      Maximum number of lines in 5km file (lines1km/5)
c elements5km   Maximum number of elements 5km fil (elements1km/5)
c
c!Output Parameters:
c modfil        Output array carrying HDF and sds info
c
c!Revision History:
c
c!Team-unique Header:
c    Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
c
c HDF portions developed at the National Center for Supercomputing
c Applications at the University of Illinois at Urbana-Champaign.
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c   External subroutines of importance:  create_hdf.f
c                                        get_spec_data.f
c                                        message.f
c
c!END
c---------------------------------------------------------------------

      implicit none

      save

      include 'mapi.inc'
      include 'hdf.inc'
      include 'dffunc.inc'


*/  Functions defined in this module

      character*(*) Prog_Number, Prog_Name
      integer create_hdf, create_vdata
      integer Lines1km, Elements1km
      integer Lines5km, Elements5km
      integer modfil(MODFILLEN)


*/  Symbolic constants for HDF-EOS

      integer HDFE_NOMERGE
      parameter (HDFE_NOMERGE=0)


*/  Variables for HDF-EOS routine

      integer swopen,swcreate,swdefdim
      integer swdetach,swattach,swclose,swwrfld
      integer ehidinfo, swdefmap
      integer fid,swathid, srtn, merge, sdid
      character*256 swathname

*/  Variables for the creation of the HDF output file.
*/  A number of these variables come from the get_spec_info routine

      character*80 ProgramName, FieldNames(150), DimNames(150)
      character*80 DataNames(150)

      integer i, j, rtn
      integer Limit, FieldStart, FieldCount, NumOfData
      parameter (Limit=1000)

      integer DimCount, DataCount(Limit), DimNumbers(Limit)
      real DataValues(Limit, Limit)
      integer*2 IntValues(Limit)
      integer LongValues(Limit)
      real    RealValues(Limit)

*/  Variables for writing vdata to the HDF file

      integer vfid, vdata_id
      integer VDataStart, VDataCount

*/  Variables for locally develpoed code.

      integer Data(5)

*-------------------------------------------------------------------------*

*/    initialize variables

      sdid = 0

*/  Open the output swath file.

      fid=swopen(Prog_Name,DFACC_CREATE)
      if( fid .eq. -1 ) then
         call message( 'setup_create_hdf',
     $   'Error creating L2 HDF product file' //
     $   ' [OPERATOR ACTION: Check that file creation is possible]',
     $   0, 2 )
      endif

*/  Create a swath.

      If (Prog_Number .eq. "35_CLRRAD") Then
         swathname="mod35_clrrad"
      else if (Prog_Number .eq. "_PRCSRD") Then
         swathname="mod_prcsrd"
      else if (Prog_Number .eq. "_PRCSR8") Then
         swathname="mod_prcsr8"
      else if (Prog_Number .eq. "_PRCSRB") Then
         swathname="mod_prcsrb"
      Else
         swathname="mod" // Prog_Number
      End If
      swathid=swcreate(fid,swathname)

*/  Get the SDS interface ID and HDF file ID for the created swath.

      srtn=ehidinfo(fid,vfid,sdid)
      if( srtn .ne. 0) then
         call message( 'setup_create_hdf',
     $   'Error getting output file id-s using ehidinfo ' //
     $   ' [OPERATOR ACTION: Notify SDST.]',
     $   0,2)
      endif

*/  Get the field names, variable names, and dimension names from
*/  the spec files


      ProgramName="MOD" // Prog_Number
      call get_spec_data(ProgramName, FieldNames, DimNames,
     $        DimNumbers, DimCount, DataNames, DataCount,
     $        DataValues, NumOfData, Limit, FieldCount)

*/--  Define dimensions for data and geo. fields.  --/*

      do i = 1,DimCount
         if (DimNames(i) .eq. "Cell_Across_Swath_5km") Then
            DimNumbers(i) = Elements5km
         else if (DimNames(i) .eq. "Cell_Along_Swath_5km") Then
            DimNumbers(i) = Lines5km
         else if (DimNames(i) .eq. "Cell_Across_Swath_1km") Then
            DimNumbers(i) = Elements1km
         else if (DimNames(i) .eq. "Cell_Along_Swath_1km") Then
            DimNumbers(i) = Lines1km
         else if (DimNames(i) .eq. "Cell_Across_Swath") Then
            DimNumbers(i) = Elements5km
         else if (DimNames(i) .eq. "Cell_Along_Swath") Then
            DimNumbers(i) = Lines5km
         else if (DimNames(i) .eq. "Number_Cells_Day") Then
            DimNumbers(i) = Lines1km
         else if (DimNames(i) .eq. "Number_Cells_Night") Then
            DimNumbers(i) = Elements1km
         End if
         if (DimNumbers(i) .eq. 0) Then
            srtn=swdefdim(swathid,DimNames(i),SD_UNLIMITED)
         else
            srtn=swdefdim(swathid,DimNames(i),DimNumbers(i))
         end if
         if( srtn .ne. 0) then
            call message( 'setup_create_hdf',
     &      'Error defining dimensions for output file: swdefdim ' //
     &      ' [OPERATOR ACTION: Notify SDST.]',
     &      0,2)
         endif
      end do

*/  Define the global variables

      If (ProgramName .eq. "MOD35_CLRRAD") Then
         srtn=sfsattr(sdid,"Granule_Start_TAI",
     &                DFNT_FLOAT64,1,DBLE(Lines5km))
         if( srtn .ne. 0) then
            call message( 'setup_create_hdf',
     &      'Error defining global var. for output file: swdefdim ',0,2)
         endif
         srtn=sfsattr(sdid,"Granule_End_TAI",
     &                DFNT_FLOAT64,1,DBLE(Elements5km))
         if( srtn .ne. 0) then
            call message( 'setup_create_hdf',
     &      'Error defining global var. for output file: swdefdim ',0,2)
         endif
      Else If (ProgramName .eq. "MOD05" .OR.
     $    ProgramName .eq. "MOD06" .OR.
     $    ProgramName .eq. "MOD35") Then

         srtn=sfsattr(sdid,"Number_of_Instrument_Scans",
     &                DFNT_INT32,1,Lines1km)
         if( srtn .ne. 0) then
            call message( 'setup_create_hdf',
     &         'Error defining dimensions for output file: swdefdim ' //
     &         ' [OPERATOR ACTION: Notify SDST.]',
     &         0,2)
         endif
         srtn=sfsattr(sdid,"Maximum_Number_of_1km_Frames",
     &                DFNT_INT32,1,Elements1km)
         if( srtn .ne. 0) then
            call message( 'setup_create_hdf',
     &         'Error defining dimensions for output file: swdefdim ' //
     &         ' [OPERATOR ACTION: Notify SDST.]',
     &         0,2)
         endif
      End If


*/  Define geolocation fields and set the attributes.
*/--------  Define data fields and set the attributes.  --------/*


*/ Data arrays

      FieldStart = 1

      do i = FieldStart, FieldCount
         merge=HDFE_NOMERGE

c rhucek 09/15/03: added argument LINES1KM to function call create_hdf
         srtn=create_hdf(LINES1KM,sdid,swathid,ProgramName,
     &                   FieldNames(i),merge)

         if( srtn .ne. 0) then
           call message( 'setup_create_hdf',
     &     'Error creating output hdf file: create_hdf ' //
     &     ' [OPERATOR ACTION: Check system resources.]',
     &     0,2)
         endif
      end do

*/  Write the 1-D SDS (vdata) attributes to the output HDF file

      VDataStart = 1
      VDataCount = 0
      if (ProgramName .eq. "MOD06") Then
         VDataStart = 1
         VDataCount = 2
      else if (ProgramName .eq. "MOD07") Then
         VDataStart = 1
         VDataCount = 2
      else if (ProgramName .eq. "MOD35") Then
         VDataStart = 1
         VDataCount = 0
      else if (ProgramName .eq. "MOD35_CLRRAD") Then
         VDataStart = 1
         VDataCount = 1
      else if (ProgramName .eq. "MOD_PRCSRD") Then
         VDataStart = 1
         VDataCount = 1
      else if (ProgramName .eq. "MOD_PRCSR8") Then
         VDataStart = 1
         VDataCount = 1
      else if (ProgramName .eq. "MOD_PRCSRB") Then
         VDataStart = 1
         VDataCount = 1
      End If

*/  Write the attributes for each field name

      srtn = vfstart(vfid)

      do i = VDataStart, VDataCount
         srtn = create_vdata(vfid, ProgramName, FieldNames(i))
         if( srtn .ne. 0) then
            call message( 'setup_create_hdf',
     &      'Error writing vdata attributes: create_vdata ' //
     &      ' [OPERATOR ACTION: Notify SDST.]',
     &      0,2)
         endif
      end do

      srtn = vfend(vfid)
      if( srtn .ne. 0) then
         call message( 'setup_create_hdf',
     &   'Error closing output file: vfend ' //
     &   ' [OPERATOR ACTION: Notify SDST.]',
     &   0,2)
      endif

*/  Detach from the swath.

      srtn=swdetach(swathid)
      if( srtn .ne. 0) then
         call message( 'setup_create_hdf',
     &   'Error detaching from swath: swdetach ' //
     &   ' [OPERATOR ACTION: Notify SDST.]',
     &   0,2)
      endif


*/  Attach to the swath

      swathid=swattach(fid,swathname)

*/  Write the 1-D data to the SDS (vdata) for each data type

      if (ProgramName .eq. "MOD06") Then
         Data(1) = 2
      else if (ProgramName .eq. "MOD07") Then
         Data(1) = 2
         Data(2) = 3
      else if (ProgramName .eq. "MOD35") Then
         Data(1) = 1
      else if (ProgramName .eq. "MOD35_CLRRAD") Then
         Data(1) = 2
      else if (ProgramName .eq. "MOD_PRCSRD") Then
         Data(1) = 2
      else if (ProgramName .eq. "MOD_PRCSR8") Then
         Data(1) = 2
      else if (ProgramName .eq. "MOD_PRCSRB") Then
         Data(1) = 2
      End If

      do i=1,NumOfData
         do j = 1, Limit
            IntValues(j) = 0
            LongValues(j) = 0
            RealValues(j) = 0
         end do
         if (Data(i) .eq. 1) Then
            do j = 1, DataCount(i)
               IntValues(j) = INT(DataValues(i,j))
            end do
            vdata_id=swwrfld(swathid, DataNames(i), 0, 1, DataCount(i),
     $                       IntValues)
            if (vdata_id .lt. 0) Then
               call message( 'setup_create_hdf',
     &           'Error writing datanames to file:vhfsd1 ' //
     &           ' [OPERATOR ACTION: Notify SDST.]',
     &           0,2)
            endif
         else if (Data(i) .eq. 2) Then
            do j = 1, DataCount(i)
               LongValues(j) = INT(DataValues(i,j))
            end do
            vdata_id=swwrfld(swathid, DataNames(i), 0, 1, DataCount(i),
     $                       LongValues)
            if (vdata_id .lt. 0) Then
               call message( 'setup_create_hdf',
     &           'Error writing datanames to file:vhfsd2 ' //
     &           ' [OPERATOR ACTION: Notify SDST.]',
     &           0,2)
            endif
         else if (Data(i) .eq. 3) Then
            do j = 1, DataCount(i)
               RealValues(j) = Real(DataValues(i,j))
            end do
            vdata_id=swwrfld(swathid, DataNames(i), 0, 1, DataCount(i),
     $                       RealValues)
            if (vdata_id .lt. 0) Then
               call message( 'setup_create_hdf',
     &           'Error writing datanames to file:vhfsd3 ' //
     &           ' [OPERATOR ACTION: Notify SDST.]',
     &           0,2)
            endif
         End If
      end do

*/  Detach from the swath

      srtn=swdetach(swathid)
      if( srtn .ne. 0) then
         call message( 'setup_create_hdf',
     &   'Error ending vdata writes to output file: vfend ' //
     &   ' [OPERATOR ACTION: Notify SDST.]',
     &   0,2)
      endif

C
C The next section sets up the Dimension Map information in the MetaData
C It defines mapping between Geolocation and Data Dimensions
C  9/20/97  WWW
C

*/  Attach to the swath

      swathid=swattach(fid,swathname)

      If (ProgramName .eq. "MOD05" .OR.
     $    ProgramName .eq. "MOD06" .OR.
     $    ProgramName .eq. "MOD35") Then
         srtn=swdefmap(swathid, "Cell_Across_Swath_5km",
     $                 "Cell_Across_Swath_1km", 2, 5)
         if( srtn .ne. 0) then
          call message( 'setup_create_hdf',
     &     'Error defing dimension map to output file: swdefmap1 ' //
     &     ' [OPERATOR ACTION: Notify SDST.]',
     &     0,2)
         endif
         srtn=swdefmap(swathid, "Cell_Along_Swath_5km",
     $                 "Cell_Along_Swath_1km", 2, 5)
         if( srtn .ne. 0) then
          call message( 'setup_create_hdf',
     &     'Error defing dimension map to output file: swdefmap2 ' //
     &     ' [OPERATOR ACTION: Notify SDST.]',
     &     0,2)
         endif
      End If

*/  Detach from the swath.

      srtn=swdetach(swathid)
      if( srtn .ne. 0) then
         call message( 'setup_create_hdf',
     &   'Error detaching from swath: swdetach ' //
     &   ' [OPERATOR ACTION: Notify SDST.]',
     &   0,2)
      endif

*/  Close the file

      srtn=swclose(fid)
      if( srtn .ne. 0) then
         call message( 'setup_create_hdf',
     &   'Error close output file: swclose ' //
     &   ' [OPERATOR ACTION: Check System Resources.]',
     &   0,2)
      endif

*/    Finally, use mapi commands to open output hdf file again

c     ... Open the output data file

      rtn = opmfil( Prog_Name, 'a',  modfil)
      if( rtn .eq. 0 ) then
        call message( 'setup_create_hdf',
     &  'Success opening output product file ', 0, 3 )
      else
        call message( 'setup_create_hdf',
     &  'Error opening output product file '  //
     &   ' [OPERATOR ACTION: Check that file was created.]',
     &   0,2)
      endif

      END
