      INTEGER FUNCTION CREATE_HDF(VFID,SDID,SWATHID,LINES1KM,
     &                            PROGRAMNAME,FIELDNAME,MERGE)

c---------------------------------------------------------------------
C!F77
c
c!Description:
c  This function reads the attribute information for a particular
c  variable (using get_spec_info) and creates an output HDF formatted
c  file.
c
c!Input parameters:
c vfid          HDF vdata interface id 
c sdid          HDF SDS interface id number
c swathid       HDF swath id number
c lines1km      number of 1km lines in MODIS granule
c programname   MODIS product file name
c fieldname     HDF SDS field name
c merge         EOSHDF merge option
c
c!Output Parameters:
c create_hdf returns 0 if exited correctly
c
c!Revision History:
c
c Rev 2: Dec 29, 2004 gfireman
c   Added and set variable messagetext to avoid illegal concatenation of 
c   variable-length string in function call.
c
c Rev 1: Sept 15, 2003 rhucek
c   Added integer "LINES1KM" to function argument list. It is used to 
c   define array element Cell_Along_Swath_Sampling[2], the effective 
c   granule position (in 1km lines) of the last line sample in an SDS. 
c
c!Team-unique Header:
c    Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
c
c HDF portions developed at the National Center for Supercomputing
c Applications at the University of Illinois at Urbana-Champaign.
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c   External subroutines of importance:  get_spec_info.f
c                                        set_array_attr.f
c                                        message.f
c
c!END
c---------------------------------------------------------------------

      implicit none

      save

      include 'dffunc.inc'
      include 'hdf.inc'
      include 'atmos_hdf_eos.inc'

c rhucek 
      integer sl, count, strlen

      character*(*) ProgramName, FieldName
      integer swathid, numbertype, merge
      integer lines1km, sdid
      byte charrange(2), charfillval
      integer*2 shortrange(2), shortfillval
      integer longrange(2), longfillval
      real floatrange(2), floatfillval
      double precision doublerange(2), scale_factor,
     &                 add_offset, doublefillval
      double precision valid_range(2), fillvalue
      character*256 attr
      integer swdefdfld, set_array_attr, swdefgfld, swsetfill,
     &        srtn, hrtn, rtn, nms, type
      integer sdsidx,sdsid
      character*85 title,history, Geolocation_Pointer
      character*80 dimlist, unit, Param_Type
      character*360 longname, messagetext, val_char
      character*85 Description(320)
      integer Cell_Across_Swath_Sampling(3)
      integer Cell_Along_Swath_Sampling(3)
      integer DescCount, TitleCount, HistoryCount
      integer i, start, stride 

      integer create_vdata, strpos, vfid

*----------------------------------------------------------------*

c   Initial all the variables to be retrieved

      valid_range(1) = 0
      valid_range(2) = 0
      longname = 'none'
      unit = 'none'
      scale_factor = 1.0
      add_offset = 0.0
      fillvalue = 0.0
      DescCount = 0
      TitleCount = 0
      HistoryCount = 0
      dimlist = 'none'
      Param_Type = 'none'
      Geolocation_Pointer = 'none'
      Cell_Across_Swath_Sampling(1) = -1
      Cell_Across_Swath_Sampling(2) = -1
      Cell_Across_Swath_Sampling(3) = -1
      Cell_Along_Swath_Sampling(1) = -1
      Cell_Along_Swath_Sampling(2) = -1
      Cell_Along_Swath_Sampling(3) = -1


c   Get the variables from the specification files

      DescCount = 0
      TitleCount = 0
      HistoryCount = 0

      call get_spec_info (ProgramName, fieldname,
     $        longname, dimlist, unit, valid_range, fillvalue,
     $        NumberType, scale_factor, add_offset, Param_Type,
     $        Cell_Across_Swath_Sampling, Cell_Along_Swath_Sampling,
     $        Description, DescCount, title, TitleCount,
     $        history, HistoryCount, Geolocation_Pointer)

c   Adjust array element Cell_Along_Swath_Sampling[2] to reflect actual
c   sampling end line 

      if (Cell_Along_Swath_Sampling(1) .gt. -1) then 
        start  = Cell_Along_Swath_Sampling(1)
        stride = Cell_Along_Swath_Sampling(3)
        Cell_Along_Swath_Sampling(2) = 
     &  ( (LINES1KM - start)/stride )*stride + start 
      end if


c=======================================================================
c                Define Science Data and Geolocation Fields  
c
c Define data fields one at a time, both geolocation and science data. 
c Geolocation fields use function swdefgfld, data fields use swdefdfld. 
c 
c The previous call to get_spec_data provides the FieldCount for the do
c loop and for the array of FieldNames.
c
c=======================================================================

c ... latitude/longitude fields
      if ( fieldname(1:9) .eq. "Longitude" .or. 
     &     fieldname(1:8) .eq. "Latitude" ) then

         srtn=swdefgfld(swathid,fieldname,dimlist,numbertype,merge)

         if ( srtn .ne. 0) then
            call message( 'create_hdf',
     &      'Error defining lat/lon fieldname for output file.' //
     &      ' [OPERATOR ACTION: Notify SDST.]', 0, 2)
         endif
      else

c ... science data fields
         srtn=swdefdfld(swathid,fieldname,dimlist,numbertype,merge)
         
         if ( srtn .ne. 0) then
            messagetext = 'Error defining fieldname for output file.' //
     &      fieldname // ' [OPERATOR ACTION: Notify SDST.]'
            call message('create_hdf',messagetext, 0, 2)
         endif
      endif


c=======================================================================
c               Set local SDS and vdata attributes  
c
c One dimensional HDF data arrays are implemented as vdata objects (or 
c tables) by the HDF library. In contrast, multi-dimensional arrays are
c written as SDS 
c objects.  Due to this difference,  
c  
c=======================================================================

      if ( strpos(dimlist, ',' ) .eq. -1) then
         srtn = vfstart(vfid)
         srtn = create_vdata(vfid, 'MOD04', fieldname)
           
         if (srtn .ne. 0) then
            call message ('setup_create_hdf',
     &      'Error writing vdata attributes: create_vdata ' //
     &      '[OPERATOR ACTION: Notify SDST.]', 0, 2)
         endif

         srtn = vfend(vfid)

         if (srtn .ne. 0) then
            call message ('setup_create_hdf',
     &      'Error closing output file: vfend ' //
     &      ' [OPERATOR ACTION: Notify SDST.]', 0, 2)
         endif
      
       else

c   Retrieve the SDS index for the data field.
      sdsidx=sfn2index(sdid,fieldname)
      sdsid=sfselect(sdid,sdsidx)
      
c   Set the valid range and fill value of the data field.
      attr="valid_range"
      type=numbertype
      nms=2

      if (valid_range(1) .ne. 0.0 .or. valid_range(2) .ne. 0.0) Then

         if (numbertype .eq. DFNT_INT8) Then
            do i = 1, 2
               charrange(i) = valid_range(i)
            end do

            hrtn=sfsattr(sdsid,attr,type,nms,charrange)

            if (hrtn .ne. 0) then
               call message( 'create_hdf',
     &         'Error setting valid range for output sds: sfsattr1.' //
     &         ' [OPERATOR ACTION: Notify SDST.]', 0, 2)
            endif

            charfillval = fillvalue
            hrtn=swsetfill(swathid,fieldname,charfillval)

            if (hrtn .ne. 0) then
               call message( 'create_hdf',
     &         'Error setting fill value for output sds: swsetfill1.' //
     &         ' [OPERATOR ACTION: Notify SDST.]', 0, 2)
            endif

         else if (numbertype .eq. DFNT_INT16) Then
            do i = 1, 2
               shortrange(i) = INT(valid_range(i))
            end do

            hrtn=sfsattr(sdsid,attr,type,nms,shortrange)

            if (hrtn .ne. 0) then
               call message( 'create_hdf',
     &         'Error setting valid range for output sds: sfsattr2.' //
     &         ' [OPERATOR ACTION: Notify SDST.]', 0, 2)
            endif

            shortfillval = fillvalue
            hrtn=swsetfill(swathid,fieldname,shortfillval)

            if (hrtn .ne. 0) then
               call message( 'create_hdf',
     &         'Error setting fill value for output sds: swsetfill2.' //
     &         ' [OPERATOR ACTION: Notify SDST.]', 0, 2)
            endif

         else if (numbertype .eq. DFNT_INT32) Then
            do i = 1, 2
               longrange(i) = INT(valid_range(i))
            end do

            hrtn=sfsattr(sdsid,attr,type,nms,longrange)

            if (hrtn .ne. 0) then
               call message( 'create_hdf',
     &         'Error setting valid range for output sds: sfsattr3.' //
     &         ' [OPERATOR ACTION: Notify SDST.]', 0, 2)
            endif

            longfillval = fillvalue
            hrtn=swsetfill(swathid,fieldname,longfillval)

            if (hrtn .ne. 0) then
               call message( 'create_hdf',
     &         'Error setting fill value for output sds: swsetfill3.' //
     &         ' [OPERATOR ACTION: Notify SDST.]', 0, 2)
            endif

         else if (numbertype .eq. DFNT_FLOAT32) Then
            do i = 1, 2
               floatrange(i) = REAL(valid_range(i))
            end do

            hrtn=sfsattr(sdsid,attr,type,nms,floatrange)

            if (hrtn .ne. 0) then
               call message( 'create_hdf',
     &         'Error setting valid range for output sds: sfsattr4.' //
     &         ' [OPERATOR ACTION: Notify SDST.]', 0, 2)
            endif

            floatfillval = REAL(fillvalue)
            hrtn=swsetfill(swathid,fieldname,floatfillval)

            if (hrtn .ne. 0) then
               call message( 'create_hdf',
     &         'Error setting fill value for output sds: swsetfill4.' //
     &         ' [OPERATOR ACTION: Notify SDST.]', 0, 2)
            endif

         else if (numbertype .eq. DFNT_FLOAT64) Then
            do i = 1, 2
               doublerange(i) = valid_range(i)
            end do

            hrtn=sfsattr(sdsid,attr,type,nms,doublerange)

            if (hrtn .ne. 0) then
               call message( 'create_hdf',
     &         'Error setting valid range for output sds: sfsattr5.' //
     &         ' [OPERATOR ACTION: Notify SDST.]', 0, 2)
            endif

            doublefillval = REAL(fillvalue)
            hrtn=swsetfill(swathid,fieldname,doublefillval)

            if (hrtn .ne. 0) then
               call message( 'create_hdf',
     &         'Error setting fill value for output sds: swsetfill5.' //
     &         ' [OPERATOR ACTION: Notify SDST.]', 0, 2)
            endif

         end if
      endif


c ... Set the attributes with same types throughout the swath.
      if (sdsid .ne. -1) then
                                         
         rtn=set_array_attr(sdsid,longname,unit,scale_factor,add_offset,
     &                      Param_Type, Cell_Across_Swath_Sampling,
     &                      Cell_Along_Swath_Sampling, Geolocation_Pointer,
     &                      Description, DescCount)

         if (rtn .ne. 0) then
            call message( 'create_hdf',
     &      'Error using set_array_attr to set attribs of output file.' //
     &      ' [OPERATOR ACTION: Notify SDST.]', 0, 2)
         endif
      endif
      endif


c ... Successfully returns.
      create_hdf=SUCCESS

      END

c ------------------------------------------------------------------/*

      INTEGER FUNCTION SET_ARRAY_ATTR(SDSID,LONGNAME,UNIT,SFACTOR,
     &                                AOFFSET,PTYPE,ACSWATH,ALSWATH,
     &                                GEOLOCATION_POINTER,
     &                                DESCRIP, DESCCOUNT)

c---------------------------------------------------------------------
C!F77
c
c!Description:
c  This function sets up some of the HDF variable attributes.
c
c!Input parameters:
c sdsid         HDF SDS id number
c
c!Output Parameters:
c longname      Longname of sds
c unit          SDS unit
c sfactor       SDS scale factor
c aoffset       SDS add/offset factor
c ptype         SDS parameter type
c acswath       SDS across swath sampling
c alswath       SDS along swath sampling
c Geolocation_Pointer  SDS Geolocation pointer
c descrip       SDS description
c DescCount     SDS description counter
c set_array_attr returns 0 if exited correctly
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
c   External subroutines of importance:  strlen.f
c                                        message.f
c
c!END
c---------------------------------------------------------------------

      implicit none

      save

      include 'hdf.inc'
      include 'dffunc.inc'
      include 'atmos_hdf_eos.inc'

      character*(*) longname,unit,ptype, Geolocation_Pointer
      character*256 attr
      character*(*) descrip(*)
      integer sdsid,acswath(3),alswath(3)
      double precision sfactor,aoffset
      integer strlen,hrtn,nms,type,DescCount,i

*--------------------------------------------------------------------*

      attr="long_name"
      type=DFNT_CHAR8
      nms=strlen(longname)
      hrtn=sfsattr(sdsid,attr,type,nms,longname)
      if( hrtn .ne. 0) then
        call message( 'set_array_attr',
     &  'Error setting long name into ouput file using sfsattr.' //
     &  ' [OPERATOR ACTION: Notify SDST.]',
     &   0,2)
      endif

      attr="units"
      type=DFNT_CHAR8
      nms=strlen(unit)
      hrtn=sfsattr(sdsid,attr,type,nms,unit)
      if( hrtn .ne. 0) then
        call message( 'set_array_attr',
     &  'Error setting units into ouput file using sfsattr.' //
     &  ' [OPERATOR ACTION: Notify SDST.]',
     &   0,2)
      endif

      attr="scale_factor"
      type=DFNT_FLOAT64
      nms=1
      hrtn=sfsattr(sdsid,attr,type,nms,sfactor)
      if( hrtn .ne. 0) then
        call message( 'set_array_attr',
     &  'Error setting scale_factor into ouput file using sfsattr.' //
     &  ' [OPERATOR ACTION: Notify SDST.]',
     &   0,2)
      endif

      attr="add_offset"
      type=DFNT_FLOAT64
      nms=1
      hrtn=sfsattr(sdsid,attr,type,nms,aoffset)
      if( hrtn .ne. 0) then
        call message( 'set_array_attr',
     &  'Error setting add_offset into ouput file using sfsattr.' //
     &  ' [OPERATOR ACTION: Notify SDST.]',
     &   0,2)
      endif

      attr="Parameter_Type"
      type=DFNT_CHAR8
      nms=strlen(ptype)
      hrtn=sfsattr(sdsid,attr,type,nms,ptype)
      if( hrtn .ne. 0) then
        call message( 'set_array_attr',
     &  'Error setting Parameter_Type into ouput file: sfsattr.' //
     &  ' [OPERATOR ACTION: Notify SDST.]',
     &   0,2)
      endif

      if(alswath(1).gt.0) then
         attr="Cell_Along_Swath_Sampling"
         type=DFNT_INT32
         nms=3
         hrtn=sfsattr(sdsid,attr,type,nms,alswath)
         if( hrtn .ne. 0) then
           call message( 'set_array_attr',
     &     'Error setting Cell Sampling into ouput file: sfsattr.' //
     &     ' [OPERATOR ACTION: Notify SDST.]',
     &      0,2)
         endif
      endif

      if(acswath(1).gt.0) then
          attr="Cell_Across_Swath_Sampling"
          type=DFNT_INT32
          nms=3
          hrtn=sfsattr(sdsid,attr,type,nms,acswath)
          if( hrtn .ne. 0) then
            call message( 'set_array_attr',
     &     'Error setting Cell Sampling into ouput file: sfsattr.' //
     &     ' [OPERATOR ACTION: Notify SDST.]',
     &     0,2)
          endif
      endif

      attr="Geolocation_Pointer"
      type=DFNT_CHAR8
      nms=strlen(Geolocation_Pointer)
      hrtn=sfsattr(sdsid,attr,type,nms,Geolocation_Pointer)
      if( hrtn .ne. 0) then
        call message( 'set_array_attr',
     &  'Error setting Geolocation Pointer to ouput file: sfsattr.' //
     &  ' [OPERATOR ACTION: Notify SDST.]',
     &  0,2)
      endif

      if (DescCount .gt. 0) then
         attr="description"
         type=DFNT_CHAR8
         do i = 1, DescCount
            descrip(i)(85:85) = char(10)
         end do
         hrtn=sfsattr(sdsid,attr,type,DescCount*85,descrip)
         if( hrtn .ne. 0) then
           call message( 'set_array_attr',
     &     'Error setting Descip.Counter to ouput file: sfsattr.' //
     &     ' [OPERATOR ACTION: Notify SDST.]',
     &     0,2)
         endif
         do i = 1, DescCount
            descrip(i)(1:85) = ' '
         end do
      endif

      set_array_attr=SUCCESS

      END

c ------------------------------------------------------------------/*

      INTEGER FUNCTION CREATE_VDATA(VFID,PROGRAMNAME,FIELDNAME)

c---------------------------------------------------------------------
C!F77
c
c!Description:
c  This function reads the attribute information for a particular
c  variable (using get_spec_info) and writes VData attributes to
C  the HDF file.
c
c!Input parameters:
c vfid          HDF VData interface id 
c programname   MODIS product file name
c fieldname     HDF SDS field name
c
c!Output Parameters:
c create_vdata returns 0 if exited correctly
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
c   External subroutines of importance:  get_spec_info.f
c                                        message.f
c
c!END
c---------------------------------------------------------------------

      implicit none

      save

      include 'hdf.inc'
      include 'dffunc.inc'
      include 'atmos_hdf_eos.inc'

      character*(*) ProgramName, FieldName
      integer i,vfid
      double precision scale_factor, valid_range(2)
      double precision fillvalue, add_offset
      integer srtn, nms, NumberType
      integer vdata_id, vdata_ref
      character*85 title, history, Geolocation_Pointer
      character*85 dimlist, longname, unit, Parameter_Type
      character*85 Description(320)
      integer Cell_Across_Swath_Sampling(3)
      integer Cell_Along_Swath_Sampling(3)
      integer DescCount, TitleCount, HistoryCount
      integer strlen

*---------------------------------------------------------------------------*

c   Get the variables from the specification files

      DescCount = 0
      TitleCount = 0
      HistoryCount = 0
      fillvalue = 0.0


c ... Get the VData attributes
      call get_spec_info (ProgramName, FieldName,
     $        longname, dimlist, unit, valid_range, fillvalue,
     $        NumberType, scale_factor, add_offset, Parameter_Type,
     $        Cell_Across_Swath_Sampling, Cell_Along_Swath_Sampling,
     $        description, DescCount, title, TitleCount,
     $        history, HistoryCount, Geolocation_Pointer)


c ... Get the VData ID
      vdata_ref = vsffnd(vfid, FieldName)
      vdata_id  = vsfatch(vfid, vdata_ref, "w")


c-----------------------------------------------------------------------
c
c ... Write the Vdata attributes to the output file
c
c-----------------------------------------------------------------------

c     Long Name:
      nms = strlen(longname)

      srtn = vsfscat(vdata_id, -1, 'long_name', DFNT_CHAR8, nms, 
     &               longname)

      if( srtn .ne. 0) then
         call message( 'create_vdata',
     &     'Error setting VData long name attrib. to file: vsfscat.' //
     &     ' [OPERATOR ACTION: Notify SDST.]',
     &     0,2)
      endif


c ... Units:
      nms = strlen(unit)

      srtn = vsfscat(vdata_id, -1, 'units', DFNT_CHAR8, nms, unit)

      if( srtn .ne. 0) then
         call message( 'create_vdata',
     &     'Error setting VData units attribute to file: vsfscat.' //
     &     ' [OPERATOR ACTION: Notify SDST.]',
     &     0,2)
      endif


c ... Geolocation Pointer
      if (ProgramName .eq. "MOD06") Then
         nms = strlen(Geolocation_Pointer)
         srtn = vsfscat(vdata_id, -1, 'Geolocation_Pointer',
     $                  DFNT_CHAR8, nms, Geolocation_Pointer)

         if( srtn .ne. 0) then
           call message( 'create_vdata',
     &    'Error setting VData Geo. Pointer attr to file: vsfscat.' //
     &     ' [OPERATOR ACTION: Notify SDST.]',
     &     0,2)
         endif
      end if


c  rhucek 09/02/99:  added description attribute section immediately below
c    Description 

      if (DescCount .gt. 0) then
         do i = 1, DescCount
            description(i)(85:85) = char(10)
         end do

         srtn = vsfscat(vdata_id, -1, 'description', DFNT_CHAR8, DescCount*85, description)

         if( srtn .ne. 0) then
           call message( 'create_vdata',
     &     'Error setting vdata "description" attribute: vsfscat.' //
     &     ' [OPERATOR ACTION: Notify SDST.]',
     &     0,2)
         endif

         do i = 1, DescCount
            description(i)(1:85) = ' '
         end do
      endif



c   Detach from the file
      srtn = vsfdtch (vdata_id)

      create_vdata = SUCCESS

      END
