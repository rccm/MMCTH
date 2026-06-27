      integer function create_mod04()

C !F77------------------------------------------------------------------
C
C !Description:
C   This routine define the swath components for MOD_PR04. This routine
C   creates a HDF file which contains the swath elements such as the swath
C   dimensions, the SDSs defined in terms of the swath dimensions, the
C   Vdata objects, the global attributes of the swath and the local
C   attributes of the SDSs and Vdata objects.
C
C
C !Input Parameters: None
C
C !Output Parameters: None
C
C !Revision History:
c 01/29/98 fhliang
c added NCSA acknowledgement.
C
C
C !Team-Unique Header:
C
C   This software was developed by the MODIS Science Data Support Team
C   (SDST) for the National Aeronautics and Space Administration,
C   Goddard Space Flight Center, under contract NAS5-32373.
C
C   HDF portions developed at the National Center for Supercomputing
C   Applications at the University of Illinois at Urbana-Champaign.
C
C
C
C !References and Credits:
C
C   Developed by JC Guu 03/10/97 jguu@ltpmail.gsfc.nasa.gov
C   Modified by Fay Liang & Rich Hucek September 1997
C   fhliang@ltpmail.gsfc.nasa.gov & rhucek@ltpmail.gsfc.nasa.gov
C
C !Externals:
C   Functions:
C    pgs_pc_getreference     (PGS_PC_GetReference.c)
C    swopen                  (SWapi.c)
C    swcreate                (SWapi.c)
C    swdefdim                (SWapi.c)
C    swdefdfld               (SWapi.c)
C    swdetach                (SWapi.c)
C    swclose                 (SWapi.c)
C    ehidinfo                (SWapi.c)
C    swsetfill               (SWapi.c)
C    sfdimid                 (dffunc.inc)
C    sfscatt                 (dffunc.inc)
C    sfsdmname               (dffunc.inc)
C    sfsnatt                 (dffunc.inc)
C
C  Named Constants:
C    PGSd_PC_FILE_PATH_MAX   (PGS_PC.f)
C    MODFILLEN               (mapi.inc)
C
C !Internals:
C   Variables:
C    acswath(3)       Across swath sampling
C    alswath(3)       Along swath sampling
C    aoffset          Added offset
C    attr             Attribute name
C    band_land(7)     Data array
C    band_ocean(7)    Data array
C    crtn             Return value from one of create_array_... function
C    descrip          Description
C    dimid            Dimension ID
C    dimn             Dimemsion value
C    dimname          Dimension name
C    dimlist          Dimension list
C    dtype            Data type of the array
C    fid              HDF-EOS file ID
C    fieldname        SDS array name
C    file_version     Version number for the file
C    fillval_f32      Fill value for 32-bit floating point type
C    fillval_f64      Fill value for 64-bit floating point type
C    fillval_i16      Fill value for 16-bit signed integer type
C    fillval_i32      Fill value for 32-bit signed integer type
C    fillval_i8       Fill value for 8-bit signed integer type
C    fname_l1b        File name of the input L1B file
C    fname_sw         File name of the output swath file
C    funcname         Character string for function names
C    hdfid            HDF File ID
C    hrtn             Return value from an HDF function
C    istr             String obtained from int2str()
C    longname         Long name of SDS
C    merge            Merge code
C    modfil_l1b       Array used to reference the L1B file in all
C                     other M-API routines
C    modfil_sw        Array used to reference the swath file in
C                     all other M-API routines
C    mrtn             Return value from a M-API function
C    nframe           Number of frames possible
C    nms              Number of values for an attribute
C    numbertype       Data type code
C    numscan          Number of scan in granule
C    prtn             Return value from an SDP Toolkit function
C    ptype            Parameter type
C    rtn              Return code from a function
C    sdid             SDS interface ID
C    sdsid            SDS ID
C    sds_unit         Data units
C    sfactor          Scaling factor
C    sl               Actual string length
C    srtn             Return code from Swath API
C    swathid          Swath structure ID
C    swathname        Swath structure name
C    vrange_f64(2)    Valid range for 64-bit floating point type
C    vrange_f32(2)    Valid range for 32-bit floating point type
C    vrange_i32(2)    Valid range for 32-bit signed integer type
C    vrange_i16(2)    Valid range for 16-bit signed integer type
C    vrange_i8(2)     Valid range for 8-bit signed integer type
C    usrlog           Character string for log messages
C
C   Named Constants:
C    HDFE_NOMERGE     Code for 'no merge' action
C    HDFE_AUTOMERGE   Code for 'merge' action
C
C   Functions and Subroutines:
C    create_array_dbl    Creates a 64-bit-floating-point array
C    create_array_flt    Creates a 32-bit-floating-point array
C    create_array_lng    Creates a 32-bit-signed-integer array
C    create_array_sht    Creates a 16-bit-signed-integer array
C    create_array_byt    Creates an 8-bit-signed-integer array
C    string_loc          Finds the position (first and last bytes) of the
C                        nonblank characters within a string buffer
C    strlen              Returns actual length of a string in a
C                        string buffer
C !END--------------------------------------------------------------


      implicit none

      include 'PGS_PC.f'
      include 'PGS_MODIS_39500.f'
      include 'mapi.inc'
      include 'hdf.inc'
      include 'dffunc.inc'
      include 'mapi_hdfeos.inc'
      include 'ckstatus.inc'
      include 'extra.inc'
c  Number of models
       integer num_index_models
       parameter(num_index_models=9)

*/  Functions
      integer pgs_pc_getreference
      integer swopen,swcreate,swdefdim,swdefdfld,swwrdmeta,swdetach,
     &        swsetfill,swwrfld,swclose,ehidinfo

c rhucek 11/11/04: removed declarations of sfscatt and sfsnatt.  They 
c are redundant with the declaration present in include file dffunc.inc. 
c             sfscatt,sfsnatt     


c     integer sfdimid,sfsdmname,
c    &        vfstart,vfend,vsfatch,vsfdtch,vsffnd,vsfgid,vsfscat,
c    &        vfend
      integer create_array_dbl,create_array_flt,create_array_lng,
     &        create_array_sht,create_array_byt,
     &        addtext
      integer strlen, string_loc


*/  Variables for SDP toolkit routines
      character*(PGSd_PC_FILE_PATH_MAX) fname_l1b,fname_sw
      integer file_version,prtn


*/  Variables for M-API routines
      character*256 dtype,attr,sdsattr
      integer modfil_l1b(MODFILLEN),modfil_sw(MODFILLEN)
      integer numscan,nframe,nms
      integer mrtn


*/  Symbolic constants for HDF-EOS
      integer HDFE_NOMERGE
      parameter (HDFE_NOMERGE=0)
      integer HDFE_AUTOMERGE
      parameter (HDFE_AUTOMERGE=1)


*/  Variables for HDF-EOS routines
      integer fid,swathid,hdfid,sdid,numbertype,merge
      integer srtn,dimn,val_i32,count
      character*256 swathname,attrname,fieldname,dimname,dimlist,
     &              longname,sds_unit,ptype,geoptr
      character*1024 text
      character*2048 val_char,descrip
      double precision fillval_f64,vrange_f64(2),sfactor,aoffset
      real fillval_f32,vrange_f32(2)
      integer fillval_i32,vrange_i32(2),acswath(3),alswath(3)
      integer*2 land_1(2), land_2(3), ocean_data(2), band_land(7),
     &  land_4(4),band_ocean(7),land_3(3),ocean_data_index(9)
      integer*2 fillval_i16,vrange_i16(2)
      byte fillval_i8,vrange_i8(2)


*/  Variables for native HDF routines
      integer sdsid,hrtn,dimid,vsid,vsref,vrtn


*/  Variables for locally develpoed code.
      character*512 funcname,usrlog
      parameter (funcname = "create_mod04")
      character*25 istr
      integer fbyte1, lbyte1
      integer fbyte2, lbyte2

      character*20 hlog
      integer rtn, crtn, sl


*/  Array initilization
      data hlog /"Defining data field "/
      data land_1/470,660/
      data land_2/470,550,660/
      data land_3/470,660,2130/
      data land_4/470,555,660,2130/
      data ocean_data/1,2/
      data ocean_data_index/1,2,3,4,5,6,7,8,9/
      data band_land/470,555,659,865,1240,1640,2130/
      data band_ocean/470,555,659,865,1240,1640,2130/
      data acswath/1,1354,3/
      data alswath/1,2023,3/


*----------------------------------------------------------------------*



*/  Get the filenames of the input L1B file.
      file_version=1
      prtn=pgs_pc_getreference(LUN_l1b,file_version,fname_l1b)

      WRITE (istr, '(I25)') LUN_l1b
      rtn = string_loc(istr, fbyte1, lbyte1)
      usrlog="Retrieving filename for LUN "//istr(fbyte1:lbyte1)//
     +       " - pgs_pc_getreference"
      CALL ckstatus_s(prtn,usrlog,funcname,LOGFLAG)


*/  Opening the input file
      mrtn=opmfil(fname_l1b,"r",modfil_l1b)

      rtn = string_loc(fname_l1b, fbyte2, lbyte2)

      IF (mrtn .EQ. MAPIOK) THEN
        usrlog = "Opening input file "//fname_l1b(fbyte2:lbyte2)//
     +         " - opmfil"
        CALL modis_smf_setdynamicmsg(modis_s_generic, usrlog, funcname)
      ELSE
        usrlog="Error opening input file "//fname_l1b(fbyte2:lbyte2)//
     1         " - opmfil"
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
        CALL modis_smf_setdynamicmsg(modis_e_generic, usrlog, funcname)
        create_mod04 = FAIL
        RETURN
      ENDIF


*/  Get attribute "Number of Scans" of the granule.
      attr="Number of Scans"
      dtype="INTEGER*4"
      nms=1
      mrtn=gmfin(modfil_l1b,attr,dtype,nms,numscan)

      rtn = string_loc(attr, fbyte1, lbyte1)
      usrlog="Retrieving the value of attribute " //attr(fbyte1:lbyte1)
     +       // " from "  //fname_l1b(fbyte2:lbyte2)// " - gmfin"
      CALL ckstatus_f(mrtn,usrlog,funcname,LOGFLAG)


*/  Get the attribute of "Max Earth View Frames" of the granule.
      attr="Max Earth View Frames"
      mrtn=gmfin(modfil_l1b,attr,dtype,nms,nframe)

      rtn = string_loc(attr, fbyte1, lbyte1)
      usrlog="Retrieving the value of attribute " //attr(fbyte1:lbyte1)
     +       // " from "  //fname_l1b(fbyte2:lbyte2)// " - gmfin"
      CALL ckstatus_f(mrtn,usrlog,funcname,LOGFLAG)


*/  Close the L1B file.
      mrtn=clmfil(modfil_l1b)

      IF (mrtn .EQ. MAPIOK) THEN
        usrlog="Closing input file "//fname_l1b(fbyte2:lbyte2)//
     +         " - clmfil"
        CALL modis_smf_setdynamicmsg(modis_s_generic, usrlog, funcname)
      ELSE
        usrlog="Error closing input file "//fname_l1b(fbyte2:lbyte2)//
     +         " - clmfil"
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
        CALL modis_smf_setdynamicmsg(modis_e_generic, usrlog, funcname)
        create_mod04 = FAIL
        RETURN
      ENDIF


*/  Retrieve the filename of the output swath file.
      file_version=1
      prtn=pgs_pc_getreference(LUN_sw,file_version,fname_sw)

      WRITE (istr, '(I25)') LUN_sw
      rtn = string_loc(istr, fbyte1, lbyte1)
      usrlog = "Retrieving filename for mod04, LUN "
     +          //istr(fbyte1:lbyte1)// " - pgs_pc_getreference"
      CALL ckstatus_s(prtn,usrlog,funcname,LOGFLAG)


*/  Open the output swath file.
      fid=swopen(fname_sw,DFACC_CREATE)

      rtn = string_loc(fname_sw, fbyte1, lbyte1)
      usrlog="Opening output file "//char(10)//
     &       fname_sw(fbyte1:lbyte1)//" with HDF-EOS"//" - swopen"
      CALL ckstatus_f(fid,usrlog,funcname,LOGFLAG)


*/  Create a swath.
      swathname="mod04"
      swathid=swcreate(fid,swathname)

      rtn = string_loc(swathname, fbyte1, lbyte1)
      usrlog="Creating a swath named "//swathname(fbyte1:lbyte1)//
     &       " - swcreate"
      CALL ckstatus_f(swathid,usrlog,funcname,LOGFLAG)


*/  Get the SDS interface ID and HDF file ID for the created sawth.
      srtn=ehidinfo(fid,hdfid,sdid)

      usrlog="Obtaining the HDF ID and SDS interface ID - ehidinfo"
      CALL ckstatus_f(srtn,usrlog,funcname,LOGFLAG)


*/---------  Define dimensions for data and geo. fields.  ------------/*

      dimname="Cell_Along_Swath"
      dimn=(numscan*10)/Iline-1
      srtn=swdefdim(swathid,dimname,dimn)

      usrlog="Defining dimension Cell_Along_Swath - swdefdim"
      CALL ckstatus_f(srtn,usrlog,funcname,LOGFLAG)
        
      dimname="Cell_Across_Swath"
      dimn= nframe/iline
      srtn=swdefdim(swathid,dimname,dimn)

      usrlog="Defining dimension Cell_Across_Swath - swdefdim"
      CALL ckstatus_f(srtn,usrlog,funcname,LOGFLAG)
        

      dimname="Solution_1_Land"
      dimn=2
      srtn=swdefdim(swathid,dimname,dimn)

      usrlog="Defining dimension Solution_1_Land - swdefdim"
      CALL ckstatus_f(srtn,usrlog,funcname,LOGFLAG)


      dimname="Solution_2_Land"
      dimn=3
      srtn=swdefdim(swathid,dimname,dimn)

      usrlog="Defining dimension Solution_2_Land - swdefdim"
      CALL ckstatus_f(srtn,usrlog,funcname,LOGFLAG)

      dimname="Solution_3_Land"
      dimn=3
      srtn=swdefdim(swathid,dimname,dimn)

      usrlog="Defining dimension Solution_3_Land - swdefdim"
      CALL ckstatus_f(srtn,usrlog,funcname,LOGFLAG)
      
      dimname="Solution_4_Land"
      dimn=4
      srtn=swdefdim(swathid,dimname,dimn)

      usrlog="Defining dimension Solution_4_Land - swdefdim"
      CALL ckstatus_f(srtn,usrlog,funcname,LOGFLAG)


      dimname="Solution_Ocean"
      dimn=2
      srtn=swdefdim(swathid,dimname,dimn)

      usrlog="Defining dimension Solution_Ocean - swdefdim"
      CALL ckstatus_f(srtn,usrlog,funcname,LOGFLAG)


      dimname="Solution_Index"
      dimn=9
      srtn=swdefdim(swathid,dimname,dimn)

      usrlog="Defining dimension Solution_Index - swdefdim"
      CALL ckstatus_f(srtn,usrlog,funcname,LOGFLAG)


      dimname="MODIS_Band_Land"
      dimn=7
      srtn=swdefdim(swathid,dimname,dimn)

      usrlog="Defining dimension MODIS_Band_Land - swdefdim"
      CALL ckstatus_f(srtn,usrlog,funcname,LOGFLAG)


      dimname="MODIS_Band_Ocean"
      dimn=7
      srtn=swdefdim(swathid,dimname,dimn)

      usrlog="Defining dimension MODIS_Band_Ocean - swdefdim"
      CALL ckstatus_f(srtn,usrlog,funcname,LOGFLAG)
      
       dimname="Num_By_Products"
       dimn=7
      srtn=swdefdim(swathid,dimname,dimn)

      usrlog="Defining dimension MODIS_Band_Ocean - swdefdim"
      CALL ckstatus_f(srtn,usrlog,funcname,LOGFLAG)


       

      dimname="QA_Byte_Land"
      dimn=5
      srtn=swdefdim(swathid,dimname,dimn)

      usrlog="Defining dimension QA_Byte_Land - swdefdim"
      CALL ckstatus_f(srtn,usrlog,funcname,LOGFLAG)


      dimname="QA_Byte_Ocean"
      dimn=5
      srtn=swdefdim(swathid,dimname,dimn)

      usrlog="Defining dimension QA_Byte_Ocean - swdefdim"
      CALL ckstatus_f(srtn,usrlog,funcname,LOGFLAG)



*/---------  Define global attributes of the swath.  -----------------/*

      attrname="Number_of_Instrument_Scans"
      numbertype=DFNT_INT32
      count=1
      val_i32=numscan

      srtn=sfsnatt(sdid,attrname,numbertype,count,val_i32)

      usrlog="Defining global attribute Number_of_Instrument_Scans"//
     &       " - sfsnatt"
      CALL ckstatus_f(srtn,usrlog,funcname,LOGFLAG)


      attrname="Maximum_Number_of_1km_Frames"
      numbertype=DFNT_INT32
      count=1
      val_i32=nframe

      srtn=sfsnatt(sdid,attrname,numbertype,count,val_i32)

      usrlog="Defining global attribute Maximum_Number_of_1km_Frames"//
     &       " - sfsnatt"
      CALL ckstatus_f(srtn,usrlog,funcname,LOGFLAG)


      attrname="title"
      numbertype=DFNT_CHAR
      val_char=char(10)//
     &         " MODIS HDF File Specification MOD04_L2: MODIS Level"//
     &         " 2 Aerosol over Land and  "//char(10)//
     &         " Ocean Product                                     "//
     &         "                          "//char(10)
      sl=strlen(val_char)
      count=sl

      srtn=sfscatt(sdid,attrname,numbertype,count,val_char)

      usrlog="Defining global attribute title - sfscatt"
      CALL ckstatus_f(srtn,usrlog,funcname,LOGFLAG)


c rhucek 02/05/98:  Commented "history" attribute
c     attrname="history"
c     numbertype=DFNT_CHAR
c     val_char=char(10)//
c    &         " 20-3-1997   Allen Chu/GSFC                         "//
c    &         "                         "//char(10)//
c    &         "             Shana Mattoo/GSFC                      "//
c    &         "                         "//char(10)
c     sl=strlen(val_char)
c     count=sl

c     srtn=sfscatt(sdid,attrname,numbertype,count,val_char)

c     usrlog="Defining global attribute history - sfscatt"
c     CALL ckstatus_f(srtn,usrlog,funcname,LOGFLAG)


c rhucek 02/05/98:  Added global attribute 'Slope_and_Offset_Usage'
      attrname = 'Slope_and_Offset_Usage'
      numbertype = DFNT_CHAR
      descrip =
     1 char(10)
     2 // 'The local SDS scale_factor and add_offset attributes are used for the      ' // char(10)
     3 // 'conversion of stored integer data to geophysical floating point numbers.   ' // char(10)
     4 // 'The implementation follows conventional HDF usage (See HDF Users Guide):   ' // char(10)
     5 // '                                                                           ' // char(10)
     6 // '      float value = scale_factor*(stored integer - add_offset)             ' // char(10)
     7 // '                                                                           ' // char(10)
     8 // 'The unit of the derived floating point value is indicated in the ''units''   ' // char(10)
     9 // 'local attribute which is also provided.                                    ' // char(10)

       sl = strlen(descrip)
       count = sl
       srtn = sfscatt(sdid, attrname, numbertype, count, descrip)



*/  Two geolocation fields and 52 data fields will be defined
*/  as SDSs. For each SDS, it's attributes such as dimension
*/  list, data type, merge flag, long name, data unit, valid
*/  range, fill value, scaling factor, offset, parameter
*   type, description, and geolocation pointer are set by
*/  variables and these variables are passed to an appropriate
*/  creation routines according to the data type of the SDS.
*/  4 creation routines are used here: create_array_flt,
*/  create_array_dbl,create_array_sht,create_array_byt.
*/  After callng an appropriate creation routine, a routine
*/  for status report, ckstatus_f, will be called to handle
*/  the status log.

*/-----  Define geolocation fields and set the attributes.  ----------/*

      fieldname="Longitude"
      dimlist="Cell_Across_Swath,Cell_Along_Swath"
      numbertype=DFNT_FLOAT32
      merge=HDFE_NOMERGE
      longname="Geodetic Longitude"
      sds_unit="Degrees_east"
      vrange_f32(1)= -180.0
      vrange_f32(2)=  180.0
      fillval_f32= -999.0
      sfactor=1.0
      aoffset=0.0
      ptype="MODIS Input"
      descrip=" "
      geoptr="Geolocation data not applicable"

      crtn=create_array_flt(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_f32,
     &                    fillval_f32,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      usrlog="Defining geolocation field Longitude"//
     &       " - create_array_flt"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)


      fieldname="Latitude"
      dimlist="Cell_Across_Swath,Cell_Along_Swath"
      numbertype=DFNT_FLOAT32
      merge=HDFE_NOMERGE
      longname="Geodetic Latitude"
      sds_unit="Degrees_north"
      vrange_f32(1)= -90.0
      vrange_f32(2)=  90.0
      fillval_f32= -999.0
      sfactor=1.0
      aoffset=0.0
      ptype="MODIS Input"
      descrip=" "
      geoptr="Geolocation data not applicable"

      crtn=create_array_flt(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_f32,
     &                    fillval_f32,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      usrlog="Defining geolocation field Latitude"//
     &       " - create_array_flt"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)



*/---------  Define data fields and set the attributes.  -------------/*

*/ Data arrays with the rank of 2 or higher. There are SDSs.

      fieldname="Scan_Start_Time"
      dimlist="Cell_Across_Swath,Cell_Along_Swath"
      numbertype=DFNT_FLOAT64
      merge=HDFE_NOMERGE
      longname="TAI Time at Start of Scan replicated across the "//
     &         "swath"
      sds_unit="Seconds since 1993-1-1 00:00:00.0 0"
      vrange_f64(1)= 0.0
      vrange_f64(2)= 3.1558e9
      fillval_f64= -999.0
      sfactor=1.0
      aoffset=0.0
      ptype="MODIS Input"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_dbl(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_f64,
     &                    fillval_f64,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_dbl"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)


      fieldname="Solar_Zenith"
      dimlist="Cell_Across_Swath,Cell_Along_Swath"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Solar Zenith Angle, Cell to Sun"
      sds_unit="Degrees"
      vrange_i16(1)= 0
      vrange_i16(2)= 18000
      fillval_i16= -9999
      sfactor=0.01
      aoffset=0.0
      ptype="MODIS Input"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)


      fieldname="Solar_Azimuth"
      dimlist="Cell_Across_Swath,Cell_Along_Swath"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Solar_Azimuth Angle, Cell to Sun"
      sds_unit="Degrees"
      vrange_i16(1)= -18000
      vrange_i16(2)=  18000
      fillval_i16= -9999
      sfactor=0.01
      aoffset=0.0
      ptype="MODIS Input"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)


      fieldname="Sensor_Zenith"
      dimlist="Cell_Across_Swath,Cell_Along_Swath"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Sensor_Zenith Angle, Cell to Sensor"
      sds_unit="Degrees"
      vrange_i16(1)= 0
      vrange_i16(2)= 18000
      fillval_i16= -9999
      sfactor=0.01
      aoffset=0.0
      ptype="MODIS Input"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)


      fieldname="Sensor_Azimuth"
      dimlist="Cell_Across_Swath,Cell_Along_Swath"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Sensor_Azimuth Angle, Cell to Sun"
      sds_unit="Degrees"
      vrange_i16(1)= -18000
      vrange_i16(2)=  18000
      fillval_i16= -9999
      sfactor=0.01
      aoffset=0.0
      ptype="MODIS Input"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)


      fieldname="Cloud_Mask_QA"
      dimlist="Cell_Across_Swath,Cell_Along_Swath"
      numbertype=DFNT_INT8
      merge=HDFE_NOMERGE
      longname="Cloud Mask info on 10x10 km resolution"
      sds_unit="None"
      vrange_i8(1)= 0
      vrange_i8(2)= -1
      fillval_i8= 0
      sfactor=1.0
      aoffset=0.0
      ptype="MODIS Input"


*/ Since the maximum number of continuation lines is 19
*/ for FORTRAN 77, the long description is written to
*/ the buffer in a few times with aceptable length.
*/ Function addtext is used for this purpose.

      descrip(1:1)=char(10)

      text=char(10)//
     &"Cloud_mask_QA flags:                                          "
     &//char(10)//
     &"                                                              "
     &//char(10)//
     &"                                                              "
     &//char(10)//
     &"QA Flag Name     Number of    Bit Value   Description         "
     &//char(10)//
     &"                 Bits                                         "
     &//char(10)//
     &"--------------------------------------------------------------"
     &//char(10)//
     &"Cloud Mask           1            0       Undetermined        "
     &//char(10)//
     &"                                  1       Determined          "
     &//char(10)//
     &"                                                              "
     &//char(0)

      rtn=addtext(text,descrip)

      text=char(10)//
     &"Cloud Mask           2            0       0-25% Cloudy pixels "
     &//char(10)//
     &"Quality Flag                      1       25-50% cloudy pixels"
     &//char(10)//
     &"                                  2       50-75% cloudy pixels"
     &//char(10)//
     &"                                  3       75-100%cloudy pixels"
     &//char(10)//
     &"                                                              "
     &//char(10)//
     &"Day/Night            1            0       Night               "
     &//char(0)

      rtn=addtext(text,descrip)

      text=char(10)//
     &"flag                              1       Day                 "
     &//char(10)//
     &"                                                              "
     &//char(10)//
     &"Sun glint            1            0       Yes                 "
     &//char(10)//
     &"flag                              1       No                  "
     &//char(10)//
     &"                                                              "
     &//char(10)//
     &"Snow/Ice flag        1            0       Yes                 "
     &//char(10)//
     &"                                  1       No                  "
     &//char(10)//
     &"                                                              "
     &//char(10)//
     &"Land/Water           2            0       Water (ocean)       "
     &//char(0)

      rtn=addtext(text,descrip)

      text=char(10)//
     &"flag                              1       Coastal             "
     &//char(10)//
     &"                                  2       Desert              "
     &//char(10)//
     &"                                  3       Land                "
     &//char(10)//
     &"---------------------- 1 byte total --------------------------"
     &//char(10)//char(0)

      rtn=addtext(text,descrip)


      geoptr="Internal geolocation arrays"

      crtn=create_array_byt(sdid,swathid,fieldname,dimlist,
     &                     numbertype,merge,longname,sds_unit,vrange_i8,
     &                     fillval_i8,sfactor,aoffset,ptype,acswath,
     &                     alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_byt"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)

Change
      fieldname="Scattering_Angle"
      dimlist="Cell_Across_Swath,Cell_Along_Swath"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Scattering Angle"
      sds_unit="Degrees"
      vrange_i16(1)= 0
      vrange_i16(2)= 18000
      fillval_i16= -9999
      sfactor=0.01
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)
Cendchange


      fieldname="Optical_Depth_Land_And_Ocean"
      dimlist="Cell_Across_Swath,Cell_Along_Swath"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="AOT at 0.55 micron for both ocean (best) and"//
     &      " land (corrected) with best quality data(Quality flag=3)"
      sds_unit="None"
      vrange_i16(1)= -100
      vrange_i16(2)= 5000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)
      
c   New SDS 2/1/2006       
      
      fieldname="Image_Optical_Depth_Land_And_Ocean"
      dimlist="Cell_Across_Swath,Cell_Along_Swath"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="AOT at 0.55 micron for both ocean (best) and"//
     &    " land (corrected) with all quality data (Quality flag=1,2,3)"
      sds_unit="None"
      vrange_i16(1)= -100
      vrange_i16(2)= 5000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)
c   New SDS 1/2002  
      fieldname="Optical_Depth_Ratio_Small_Land_And_Ocean"
      dimlist="Cell_Across_Swath,Cell_Along_Swath"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Ratio of small mode optical depth at 0.55 micron"
      sds_unit="None"
      vrange_i16(1)= 0
      vrange_i16(2)= 1000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"
       crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)
c
      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG) 

 



      fieldname="Aerosol_Type_Land"
      dimlist="Cell_Across_Swath,Cell_Along_Swath"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Aerosol Type"
      sds_unit="None"
      vrange_i16(1)= 0
      vrange_i16(2)= 4
      fillval_i16= -9999
      sfactor=1.0
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)
c
       fieldname="Fitting_Error_Land"
      dimlist="Cell_Across_Swath,Cell_Along_Swath"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Fitting error for Land"
      sds_unit="None"
      vrange_i16(1)= 0
      vrange_i16(2)= 1000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG) 
      
      

      fieldname="Surface_Reflectance_Land"
      dimlist="Cell_Across_Swath,Cell_Along_Swath,Solution_2_Land"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Surface Reflectance at 0.47,0.66 and 2.13micron"
      sds_unit="None"
      vrange_i16(1)= 0
      vrange_i16(2)= 5000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)
      
      
  
      


      fieldname="Corrected_Optical_Depth_Land"
      dimlist="Cell_Across_Swath,Cell_Along_Swath,Solution_3_Land"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Corrected optical thickness at 0.47, 0.55,0.66  "//
     &         " micron"
      sds_unit="None"
      vrange_i16(1)= -100
      vrange_i16(2)= 5000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)


      fieldname="Corrected_Optical_Depth_Land_wav2p1"
      dimlist="Cell_Across_Swath,Cell_Along_Swath"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Corrected optical thickness at 2.13  "//
     &         " micron"
      sds_unit="None"
      vrange_i16(1)= -100
      vrange_i16(2)= 5000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)

C New Sds 1/2/2006  

      fieldname="Optical_Depth_Small_Land"
      dimlist="Cell_Across_Swath,Cell_Along_Swath,Solution_4_Land"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname=" optical thickness for small Mode at 0.47, 0.55,0.66 & 2.13 "//
     &         " micron"
      sds_unit="None"
      vrange_i16(1)= -100
      vrange_i16(2)= 5000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"


      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)



      fieldname="Optical_Depth_Ratio_Small_Land"
      dimlist="Cell_Across_Swath,Cell_Along_Swath"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Small mode aerosol fraction"
      sds_unit="None"
      vrange_i16(1)= 0
      vrange_i16(2)= 1000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)


      fieldname="Number_Pixels_Used_Land"
      dimlist="Cell_Across_Swath,Cell_Along_Swath,Solution_1_Land"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Number of good pixels used for retrieval at 470 ,659 nm" 
      sds_unit="None"
      vrange_i16(1)= 1
      vrange_i16(2)= 400
      fillval_i16= -9999
      sfactor=1.0
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)


      fieldname="Mean_Reflectance_Land"
      dimlist="Cell_Across_Swath,Cell_Along_Swath,MODIS_Band_Land"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Mean reflectance at 0.47,0.55,0.66,0.87,1.24 ,1.64 and 2.13 Microns"
      sds_unit="None"
      vrange_i16(1)= 0
      vrange_i16(2)= 10000
      fillval_i16= -9999
      sfactor=0.0001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)


      fieldname="STD_Reflectance_Land"
      dimlist="Cell_Across_Swath,Cell_Along_Swath,MODIS_Band_Land"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Standard deviation of reflectance at 0.47,0.55,0.66,0.87,1.64 and 2.13 Microns"
      sds_unit="None"
      vrange_i16(1)= 0
      vrange_i16(2)= 20000
      fillval_i16= -9999
      sfactor=0.0001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)
      
      
      fieldname="Mass_Concentration_Land"
      dimlist="Cell_Across_Swath,Cell_Along_Swath"
      numbertype=DFNT_FLOAT32
      merge=HDFE_NOMERGE
      longname="Mass concentration"
      sds_unit="1.0e-6g/cm^2"
      vrange_f32(1)= 0.0
      vrange_f32(2)= 1000.0
      fillval_f32= -999.0
      sfactor=1.0
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_flt(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_f32,
     &                    fillval_f32,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_flt"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)


      fieldname="Angstrom_Exponent_Land"
      dimlist="Cell_Across_Swath,Cell_Along_Swath"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Angstrom exponent for 0.47 and 0.67 micron"
      sds_unit="None"
      vrange_i16(1)= -1000
      vrange_i16(2)= 5000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)




      fieldname="Cloud_Fraction_Land"
      dimlist="Cell_Across_Swath,Cell_Along_Swath"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Cloud fraction from Land aerosol cloud mask from retrieved "//
     &"and overcast pixels not including cirrus mask"
      sds_unit="None"
      vrange_i16(1)= 0
      vrange_i16(2)= 1000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)



      fieldname="Quality_Assurance_Land"
      dimlist="QA_Byte_Land,Cell_Across_Swath,Cell_Along_Swath"
      numbertype=DFNT_INT8
      merge=HDFE_NOMERGE
      longname="Runtime QA flags"
      sds_unit="None"
      vrange_i8(1)= 0
      vrange_i8(2)= -1
      fillval_i8= 0
      sfactor=1.0
      aoffset=0.0
      ptype="Output"
      descrip="see MODIS atmosphere QA plan for details"
      geoptr="Internal geolocation arrays"

      crtn=create_array_byt(sdid,swathid,fieldname,dimlist,
     &                     numbertype,merge,longname,sds_unit,vrange_i8,
     &                     fillval_i8,sfactor,aoffset,ptype,acswath,
     &                     alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_byt"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)

      
      
      
 


      fieldname="Solution_Index_Ocean_Small"
      dimlist="Cell_Across_Swath,Cell_Along_Swath,Solution_Ocean"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Solution Number index small particles"
      sds_unit="None"
      vrange_i16(1)= 1
      vrange_i16(2)= 4
      fillval_i16= -9999
      sfactor=1.0
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)


      fieldname="Solution_Index_Ocean_Large"
      dimlist="Cell_Across_Swath,Cell_Along_Swath,Solution_Ocean"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Solution Number index large particles"
      sds_unit="None"
      vrange_i16(1)= 5
      vrange_i16(2)= 9
      fillval_i16= -9999
      sfactor=1.0
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)



      fieldname="Effective_Optical_Depth_Best_Ocean"
      dimlist="Cell_Across_Swath,Cell_Along_Swath,MODIS_Band_Ocean"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="AOT at 7 bands for best solution for 0.47, 0.55,0.66,0.86,1.24,1.63,2.13 um"
      sds_unit="None"
      vrange_i16(1)= -100
      vrange_i16(2)= 5000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)


      fieldname="Effective_Optical_Depth_Average_Ocean"
      dimlist="Cell_Across_Swath,Cell_Along_Swath,MODIS_Band_Ocean"
       numbertype=DFNT_INT16
       merge=HDFE_NOMERGE
      longname="AOT  at 7 bands for Average solution for 0.47, 0.55,0.66,0.86,1.24,1.63,2.13 um"
      sds_unit="None"
      vrange_i16(1)= -100
      vrange_i16(2)= 5000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)


      fieldname="Optical_Depth_Small_Best_Ocean"
      dimlist="Cell_Across_Swath,Cell_Along_Swath,MODIS_Band_Ocean"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="AOT at 7 bands for small mode of best solution for 0.47, 0.55,0.66,0.86,1.24,1.63,2.13 um"
      sds_unit="None"
      vrange_i16(1)= -100
      vrange_i16(2)= 5000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)


      fieldname="Optical_Depth_Small_Average_Ocean"
      dimlist="Cell_Across_Swath,Cell_Along_Swath,MODIS_Band_Ocean"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="AOT at 7 bands for small mode of average solution for 0.47, 0.55,0.66,0.86,1.24,1.63,2.13 um"
      sds_unit="None"
      vrange_i16(1)= -100
      vrange_i16(2)= 5000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)


      fieldname="Optical_Depth_Large_Best_Ocean"
      dimlist="Cell_Across_Swath,Cell_Along_Swath,MODIS_Band_Ocean"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="AOT at 7 bands for large mode of best solution for 0.47, 0.55,0.66,0.86,1.24,1.63,2.13 um"
      sds_unit="None"
      vrange_i16(1)= -100
      vrange_i16(2)= 5000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)


      fieldname="Optical_Depth_Large_Average_Ocean"
      dimlist="Cell_Across_Swath,Cell_Along_Swath,MODIS_Band_Ocean"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="AOT at 7 bands for large mode of average solution for 0.47, 0.55,0.66,0.86,1.24,1.63,2.13 um"
      sds_unit="None"
      vrange_i16(1)= -100
      vrange_i16(2)= 5000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)


      fieldname="Mass_Concentration_Ocean"
      dimlist="Cell_Across_Swath,Cell_Along_Swath,Solution_Ocean"
      numbertype=DFNT_FLOAT32
      merge=HDFE_NOMERGE
      longname="Mass_Concentration for best and average solutions"
      sds_unit="1.0e-6g/cm^2"
      vrange_f32(1)= 0.0
      vrange_f32(2)= 1000.0
      fillval_f32= -999.0
      sfactor=1.0
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_flt(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_f32,
     &                    fillval_f32,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_flt"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)


      fieldname="Effective_Radius_Ocean"
      dimlist="Cell_Across_Swath,Cell_Along_Swath,Solution_Ocean"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Effective_Radius at 0.55 micron of both solutions"
      sds_unit="micron"
      vrange_i16(1)= 0
      vrange_i16(2)= 5000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)


      fieldname="Cloud_Condensation_Nuclei_Ocean"
      dimlist="Cell_Across_Swath,Cell_Along_Swath,Solution_Ocean"
      numbertype=DFNT_FLOAT32
      merge=HDFE_NOMERGE
      longname="Column number of CCN at 0.55 micron of both solutions"
      sds_unit="CCN/cm^2"
      vrange_f32(1)= 0.0
      vrange_f32(2)= 10.0e10
      fillval_f32= -999.0
      sfactor=1.0
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_flt(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_f32,
     &                    fillval_f32,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_flt"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)


      fieldname="Asymmetry_Factor_Best_Ocean"
      dimlist="Cell_Across_Swath,Cell_Along_Swath,MODIS_Band_Ocean"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Asymmetry_Factor at 7 bands for best solution for 0.47, 0.55,0.66,0.86,1.24,1.63,2.13 um"
      sds_unit="None"
      vrange_i16(1)= 0
      vrange_i16(2)= 3000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)


      fieldname="Asymmetry_Factor_Average_Ocean"
      dimlist="Cell_Across_Swath,Cell_Along_Swath,MODIS_Band_Ocean"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Asymmetry_Factor at 7 bands for average solution for 0.47, 0.55,0.66,0.86,1.24,1.63,2.13 um"
      sds_unit="None"
      vrange_i16(1)= 0
      vrange_i16(2)= 3000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)


      fieldname="Backscattering_Ratio_Best_Ocean"
      dimlist="Cell_Across_Swath,Cell_Along_Swath,MODIS_Band_Ocean"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Backscattering ratio at 7 bands for best solution for 0.47, 0.55,0.66,0.86,1.24,1.63,2.13 um"
      sds_unit="None"
      vrange_i16(1)= 0
      vrange_i16(2)= 3000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)


      fieldname="Backscattering_Ratio_Average_Ocean"
      dimlist="Cell_Across_Swath,Cell_Along_Swath,MODIS_Band_Ocean"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Backscattering ratio at 7 bands for average solution for 0.47, 0.55,0.66,0.86,1.24,1.63,2.13 um"
      sds_unit="None"
      vrange_i16(1)= 0
      vrange_i16(2)= 3000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)


      fieldname="Angstrom_Exponent_1_Ocean"
      dimlist="Cell_Across_Swath,Cell_Along_Swath,Solution_Ocean"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Angstrom Exponent for 0.55 and 0.86 micron"
      sds_unit="None"
      vrange_i16(1)= -1000
      vrange_i16(2)= 5000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)


      fieldname="Angstrom_Exponent_2_Ocean"
      dimlist="Cell_Across_Swath,Cell_Along_Swath,Solution_Ocean"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Angstrom exponent for 0.865 and 2.130 micron"
      sds_unit="None"
      vrange_i16(1)= -1000
      vrange_i16(2)= 5000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)

 

   

      fieldname="Least_Squares_Error_Ocean"
      dimlist="Cell_Across_Swath,Cell_Along_Swath,Solution_Ocean"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Least square error estimated"
      sds_unit="None"
      vrange_i16(1)= 0
      vrange_i16(2)= 1000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)


      fieldname="Optical_Depth_Ratio_Small_Ocean_0.55micron"
      dimlist="Cell_Across_Swath,Cell_Along_Swath,Solution_Ocean"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Ratio of small mode optical depth at 0.55 microns for best (1) and average (2) solutions"
      sds_unit="None"
      vrange_i16(1)= 0
      vrange_i16(2)= 1000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)


      fieldname="Optical_Depth_by_models_ocean"
      dimlist="Cell_Across_Swath,Cell_Along_Swath,Solution_Index"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="optical depth for small and large modes"//
     &         " placed in model index"
      sds_unit="None"
      vrange_i16(1)= -100
      vrange_i16(2)= 5000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)


      fieldname="Cloud_Fraction_Ocean"
      dimlist="Cell_Across_Swath,Cell_Along_Swath"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Cloud fraction from Ocean  aerosol cloud mask from retrieved "//
     &   "and overcast pixels not including cirrus mask"
      sds_unit="None"
      vrange_i16(1)= 0
      vrange_i16(2)= 1000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)


      fieldname="Number_Pixels_Used_Ocean"
      dimlist="Cell_Across_Swath,Cell_Along_Swath"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Number of good pixels used for retrieval at at 865 nm"
      sds_unit="None"
      vrange_i16(1)= 1
      vrange_i16(2)= 400
      fillval_i16= -9999
      sfactor=1.0
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)


      fieldname="Mean_Reflectance_Ocean"
      dimlist="Cell_Across_Swath,Cell_Along_Swath,MODIS_Band_Ocean"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Mean reflectances at 7 bands for 0.47, 0.55,0.66,0.86,1.24,1.63,2.13 micron"
      sds_unit="None"
      vrange_i16(1)= 0
      vrange_i16(2)= 10000
      fillval_i16= -9999
      sfactor=0.0001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)


      fieldname="STD_Reflectance_Ocean"
      dimlist="Cell_Across_Swath,Cell_Along_Swath,MODIS_Band_Ocean"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Standard devaition of reflectances at 7 bands for 0.47, 0.55,0.66,0.86,1.24,1.63,2.13 um"
      sds_unit="None"
      vrange_i16(1)= 0
      vrange_i16(2)= 20000
      fillval_i16= -9999
      sfactor=0.0001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)


      fieldname="Quality_Assurance_Ocean"
      dimlist="QA_Byte_Ocean,Cell_Across_Swath,Cell_Along_Swath"
      numbertype=DFNT_INT8
      merge=HDFE_NOMERGE
      longname="Run time QA flags"
      sds_unit="None"
      vrange_i8(1)= 0
      vrange_i8(2)= -1
      fillval_i8= 0
      sfactor=1.0
      aoffset=0.0
      ptype="Output"
      descrip="(see MODIS atmosphere QA plan for details)"
      geoptr="Internal geolocation arrays"

      crtn=create_array_byt(sdid,swathid,fieldname,dimlist,
     &                     numbertype,merge,longname,sds_unit,vrange_i8,
     &                     fillval_i8,sfactor,aoffset,ptype,acswath,
     &                     alswath,descrip,geoptr)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_byt"
      CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)
      
      

  


*/  One-dimensional data fields are implemented as Vdata objetcs
*/  by the HDF-EOS library. They are used to to construct 3D or
*/  higher-dimensional arrays from 2D arrays. The HDF-EOS library
*/  routine swdefdfld is used to define a data field.

      fieldname="Solution_1_Land"
      dimlist="Solution_1_Land"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE

      srtn=swdefdfld(swathid,fieldname,dimlist,numbertype,merge)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - swdefdfld"
      CALL ckstatus_f(srtn,usrlog,funcname,LOGFLAG)



      fieldname="Solution_2_Land"
      dimlist="Solution_2_Land"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE

      srtn=swdefdfld(swathid,fieldname,dimlist,numbertype,merge)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - swdefdfld"
      CALL ckstatus_f(srtn,usrlog,funcname,LOGFLAG)



c Added SDS Solution_3_Land 
      fieldname="Solution_3_Land"
      dimlist="Solution_3_Land"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE

      srtn=swdefdfld(swathid,fieldname,dimlist,numbertype,merge)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - swdefdfld"
      CALL ckstatus_f(srtn,usrlog,funcname,LOGFLAG)
      
c Added SDS Solution_4_Land 
      fieldname="Solution_4_Land"
      dimlist="Solution_4_Land"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE

      srtn=swdefdfld(swathid,fieldname,dimlist,numbertype,merge)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - swdefdfld"
      CALL ckstatus_f(srtn,usrlog,funcname,LOGFLAG)
      
      



      fieldname="Solution_Ocean"
      dimlist="Solution_Ocean"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE

       srtn=swdefdfld(swathid,fieldname,dimlist,numbertype,merge)

       rtn = string_loc(fieldname, fbyte1, lbyte1)
       usrlog=hlog//fieldname(fbyte1:lbyte1)//" - swdefdfld"
       CALL  ckstatus_f(srtn,usrlog,funcname,LOGFLAG)



c Added SDS Solution_Index
      fieldname="Solution_Index"
      dimlist="Solution_Index"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE

       srtn=swdefdfld(swathid,fieldname,dimlist,numbertype,merge)

       rtn = string_loc(fieldname, fbyte1, lbyte1)
       usrlog=hlog//fieldname(fbyte1:lbyte1)//" - swdefdfld"
       CALL  ckstatus_f(srtn,usrlog,funcname,LOGFLAG)

     

      fieldname="MODIS_Band_Land"
      dimlist="MODIS_Band_Land"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE

      srtn=swdefdfld(swathid,fieldname,dimlist,numbertype,merge)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - swdefdfld"
      CALL ckstatus_f(srtn,usrlog,funcname,LOGFLAG)



      fieldname="MODIS_Band_Ocean"
      dimlist="MODIS_Band_Ocean"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE

      srtn=swdefdfld(swathid,fieldname,dimlist,numbertype,merge)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog=hlog//fieldname(fbyte1:lbyte1)//" - swdefdfld"
      CALL ckstatus_f(srtn,usrlog,funcname,LOGFLAG)




*/----- Set the attributes for 1-D arrays (Vdata objects) ------------/*

*/  Native HDF routines are used to set the attributes.
      vrtn=vfstart(hdfid)

      usrlog="Starting vdata interface - vfstart"
      CALL ckstatus_f(vrtn,usrlog,funcname,LOGFLAG)



*/ Vdata MODIS_Band_Land
      fieldname="MODIS_Band_Land"

c Getting reference number for vdata
      vsref=vsffnd(hdfid,fieldname)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog="Getting reference number for vdata "//
     &        fieldname(fbyte1:lbyte1)//" - vsffnd"
      CALL ckstatus_f(vsref,usrlog,funcname,LOGFLAG)


c Attaching vdata
      vsid=vsfatch(hdfid,vsref,"w")

      usrlog="Attaching vdata "//
     &        fieldname(fbyte1:lbyte1)//" - vsfatch"
      CALL ckstatus_f(vsid,usrlog,funcname,LOGFLAG)


      attr="long_name"
      longname="Center Wavelengths of MODIS Bands Used in Land Retrieval Algorithms"//char(0)
      nms=strlen(longname)
      numbertype=DFNT_CHAR

      vrtn=vsfscat(vsid,-1,attr,numbertype,nms,longname)

      usrlog="Setting long name for "//
     &        fieldname(fbyte1:lbyte1)//" - vsfscat"
      CALL ckstatus_f(vrtn,usrlog,funcname,LOGFLAG)


      attr="units"
      sds_unit="Nanometers"//char(0)
      nms=strlen(sds_unit)
      numbertype=DFNT_CHAR

      vrtn=vsfscat(vsid,-1,attr,numbertype,nms,sds_unit)

      usrlog="Setting units for "//
     &        fieldname(fbyte1:lbyte1)//" - vsfscat"
      CALL ckstatus_f(vrtn,usrlog,funcname,LOGFLAG)


      attr="Geolocation_Pointer"
      geoptr="Geolocation data not applicable"//char(0)
      nms=strlen(geoptr)
      numbertype=DFNT_CHAR

      vrtn=vsfscat(vsid,-1,attr,numbertype,nms,geoptr)

      usrlog="Setting geolocation pointer for "//
     &        fieldname(fbyte1:lbyte1)//" - vsfscat"
      CALL ckstatus_f(vrtn,usrlog,funcname,LOGFLAG)


c Detaching vdata
      vrtn=vsfdtch(vsid)

      usrlog="Detaching vdata"//fieldname(fbyte1:lbyte1)//" - vsfdtch"
      CALL ckstatus_f(vrtn,usrlog,funcname,LOGFLAG)



*/ Vdata MODIS_Band_Ocean
      fieldname="MODIS_Band_Ocean"

c Getting reference number for vdata
      vsref=vsffnd(hdfid,fieldname)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog="Getting reference number for vdata "//
     &        fieldname(fbyte1:lbyte1)//" - vsffnd"
      CALL ckstatus_f(vsref,usrlog,funcname,LOGFLAG)


c Attaching vdata
      vsid=vsfatch(hdfid,vsref,"w")

      usrlog="Attaching vdata "//
     &        fieldname(fbyte1:lbyte1)//" - vsfatch"
      CALL ckstatus_f(vsid,usrlog,funcname,LOGFLAG)


      attr="long_name"
      longname="Center Wavelengths of MODIS Bands Used in Ocean Retrieval Algorithms"//char(0)
      nms=strlen(longname)
      numbertype=DFNT_CHAR

      vrtn=vsfscat(vsid,-1,attr,numbertype,nms,longname)

      usrlog="Setting long name for "//
     &        fieldname(fbyte1:lbyte1)//" - vsfscat"
      CALL ckstatus_f(vrtn,usrlog,funcname,LOGFLAG)


      attr="units"
      sds_unit="Nanometers"//char(0)
      nms=strlen(sds_unit)
      numbertype=DFNT_CHAR

      vrtn=vsfscat(vsid,-1,attr,numbertype,nms,sds_unit)

      usrlog="Setting units for "//
     &        fieldname(fbyte1:lbyte1)//" - vsfscat"
      CALL ckstatus_f(vrtn,usrlog,funcname,LOGFLAG)


      attr="Geolocation_Pointer"
      geoptr="Geolocation data not applicable"//char(0)
      nms=strlen(geoptr)
      numbertype=DFNT_CHAR

      vrtn=vsfscat(vsid,-1,attr,numbertype,nms,geoptr)

      usrlog="Setting geolocation pointer for "//
     &        fieldname(fbyte1:lbyte1)//" - vsfscat"
      CALL ckstatus_f(vrtn,usrlog,funcname,LOGFLAG)


c Detaching vdata
      vrtn=vsfdtch(vsid)

      usrlog="Detaching vdata"//fieldname(fbyte1:lbyte1)//" - vsfdtch"
      CALL ckstatus_f(vrtn,usrlog,funcname,LOGFLAG)



*/ Vdata Solution_1_Land
      fieldname="Solution_1_Land"

c Getting reference number for vdata
      vsref=vsffnd(hdfid,fieldname)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog="Getting reference number for vdata "//
     &        fieldname(fbyte1:lbyte1)//" - vsffnd"
      CALL ckstatus_f(vsref,usrlog,funcname,LOGFLAG)


c Attaching vdata
      vsid=vsfatch(hdfid,vsref,"w")

      attr="long_name"
      longname="Central Wavelength of MODIS Bands Used in "
     +//"Continental Model Retrieval"//char(0)
      nms=strlen(longname)
      numbertype=DFNT_CHAR

      vrtn=vsfscat(vsid,-1,attr,numbertype,nms,longname)

      usrlog="Setting long name for "//
     &        fieldname(fbyte1:lbyte1)//" - vsfscat"
      CALL ckstatus_f(vrtn,usrlog,funcname,LOGFLAG)


      attr="units"
      sds_unit="Nanometers"//char(0)
      nms=strlen(sds_unit)
      numbertype=DFNT_CHAR

      vrtn=vsfscat(vsid,-1,attr,numbertype,nms,sds_unit)

      usrlog="Setting units for "//
     &        fieldname(fbyte1:lbyte1)//" - vsfscat"
      CALL ckstatus_f(vrtn,usrlog,funcname,LOGFLAG)


      attr="Geolocation_Pointer"
      geoptr="Geolocation data not applicable"//char(0)
      nms=strlen(geoptr)
      numbertype=DFNT_CHAR

      vrtn=vsfscat(vsid,-1,attr,numbertype,nms,geoptr)

      usrlog="Setting geolocation pointer for "//
     &        fieldname(fbyte1:lbyte1)//" - vsfscat"
      CALL ckstatus_f(vrtn,usrlog,funcname,LOGFLAG)


c Detaching vdata Solution_1_Land
      vrtn=vsfdtch(vsid)

      usrlog="Detaching vdata"//fieldname(fbyte1:lbyte1)//" - vsfdtch"
      CALL ckstatus_f(vrtn,usrlog,funcname,LOGFLAG)


*/ Vdata Solution_2_Land
      fieldname="Solution_2_Land"

c Getting reference number for vdata
      vsref=vsffnd(hdfid,fieldname)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog="Getting reference number for vdata "//
     &        fieldname(fbyte1:lbyte1)//" - vsffnd"
      CALL ckstatus_f(vsref,usrlog,funcname,LOGFLAG)


c Attaching vdata
      vsid=vsfatch(hdfid,vsref,"w")

      usrlog="Attaching vdata "//
     &        fieldname(fbyte1:lbyte1)//" - vsfatch"
      CALL ckstatus_f(vsid,usrlog,funcname,LOGFLAG)


      attr="long_name"
      longname="Central Wavelength of MODIS Bands Used in "
     +//"Corrected Model Retrieval"//char(0)
      nms=strlen(longname)
      numbertype=DFNT_CHAR

      vrtn=vsfscat(vsid,-1,attr,numbertype,nms,longname)

      usrlog="Setting long name for "//
     &        fieldname(fbyte1:lbyte1)//" - vsfscat"
      CALL ckstatus_f(vrtn,usrlog,funcname,LOGFLAG)


      attr="units"
      sds_unit="Nanometers"//char(0)
      nms=strlen(sds_unit)
      numbertype=DFNT_CHAR

      vrtn=vsfscat(vsid,-1,attr,numbertype,nms,sds_unit)

      usrlog="Setting units for "//
     &        fieldname(fbyte1:lbyte1)//" - vsfscat"
      CALL ckstatus_f(vrtn,usrlog,funcname,LOGFLAG)


      attr="Geolocation_Pointer"
      geoptr="Geolocation data not applicable"//char(0)
      nms=strlen(geoptr)
      numbertype=DFNT_CHAR

      vrtn=vsfscat(vsid,-1,attr,numbertype,nms,geoptr)

      usrlog="Setting geolocation pointer for "//
     &        fieldname(fbyte1:lbyte1)//" - vsfscat"
      CALL ckstatus_f(vrtn,usrlog,funcname,LOGFLAG)


c Detaching vdata Solution_2_Land
      vrtn=vsfdtch(vsid)

      usrlog="Detaching vdata"//fieldname(fbyte1:lbyte1)//" - vsfdtch"
      CALL ckstatus_f(vrtn,usrlog,funcname,LOGFLAG)


c rhucek 05/10/01:  Added attributes for "Solution_3_Land"
*/ Vdata Solution_3_Land
      fieldname="Solution_3_Land"

c Getting reference number for vdata
      vsref=vsffnd(hdfid,fieldname)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog="Getting reference number for vdata "//
     &        fieldname(fbyte1:lbyte1)//" - vsffnd"
      CALL ckstatus_f(vsref,usrlog,funcname,LOGFLAG)


c Attaching vdata
      vsid=vsfatch(hdfid,vsref,"w")

      usrlog="Attaching vdata "//
     &        fieldname(fbyte1:lbyte1)//" - vsfatch"
      CALL ckstatus_f(vsid,usrlog,funcname,LOGFLAG)


      attr="long_name"
      longname="Central Wavelength of MODIS Bands Used in "
     +//"Scatter Plot Solution"//char(0)
      nms=strlen(longname)
      numbertype=DFNT_CHAR

      vrtn=vsfscat(vsid,-1,attr,numbertype,nms,longname)

      usrlog="Setting long name for "//
     &        fieldname(fbyte1:lbyte1)//" - vsfscat"
      CALL ckstatus_f(vrtn,usrlog,funcname,LOGFLAG)


      attr="units"
      sds_unit="Nanometers"//char(0)
      nms=strlen(sds_unit)
      numbertype=DFNT_CHAR

      vrtn=vsfscat(vsid,-1,attr,numbertype,nms,sds_unit)

      usrlog="Setting units for "//
     &        fieldname(fbyte1:lbyte1)//" - vsfscat"
      CALL ckstatus_f(vrtn,usrlog,funcname,LOGFLAG)


      attr="Geolocation_Pointer"
      geoptr="Geolocation data not applicable"//char(0)
      nms=strlen(geoptr)
      numbertype=DFNT_CHAR

      vrtn=vsfscat(vsid,-1,attr,numbertype,nms,geoptr)

      usrlog="Setting geolocation pointer for "//
     &        fieldname(fbyte1:lbyte1)//" - vsfscat"
      CALL ckstatus_f(vrtn,usrlog,funcname,LOGFLAG)


c Detaching vdata Solution_3_Land
      vrtn=vsfdtch(vsid)

      usrlog="Detaching vdata"//fieldname(fbyte1:lbyte1)//" - vsfdtch"
      CALL ckstatus_f(vrtn,usrlog,funcname,LOGFLAG)


*/ Vdata Solution_Ocean
      fieldname="Solution_Ocean"

c Getting reference number for vdata
      vsref=vsffnd(hdfid,fieldname)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog="Getting reference number for vdata "//
     &        fieldname(fbyte1:lbyte1)//" - vsffnd"
      CALL ckstatus_f(vsref,usrlog,funcname,LOGFLAG)


c Attaching vdata
      vsid=vsfatch(hdfid,vsref,"w")

      usrlog="Attaching vdata "//
     &        fieldname(fbyte1:lbyte1)//" - vsfatch"
      CALL ckstatus_f(vsid,usrlog,funcname,LOGFLAG)


      attr="long_name"
      longname="Index of Ocean Solution Types: "
     +//"1 - Best; 2 - Average"//char(0)
      nms=strlen(longname)
      numbertype=DFNT_CHAR

      vrtn=vsfscat(vsid,-1,attr,numbertype,nms,longname)

      usrlog="Setting long name for "//
     &        fieldname(fbyte1:lbyte1)//" - vsfscat"
      CALL ckstatus_f(vrtn,usrlog,funcname,LOGFLAG)


      attr="units"
      sds_unit="None"//char(0)
      nms=strlen(sds_unit)
      numbertype=DFNT_CHAR

      vrtn=vsfscat(vsid,-1,attr,numbertype,nms,sds_unit)

      usrlog="Setting units for "//
     &        fieldname(fbyte1:lbyte1)//" - vsfscat"
      CALL ckstatus_f(vrtn,usrlog,funcname,LOGFLAG)


      attr="Geolocation_Pointer"
      geoptr="Geolocation data not applicable"//char(0)
      nms=strlen(geoptr)
      numbertype=DFNT_CHAR

      vrtn=vsfscat(vsid,-1,attr,numbertype,nms,geoptr)

      usrlog="Setting geolocation pointer for "//
     &        fieldname(fbyte1:lbyte1)//" - vsfscat"
      CALL ckstatus_f(vrtn,usrlog,funcname,LOGFLAG)


c Detaching vdata Solution_Ocean
      vrtn=vsfdtch(vsid)

      usrlog="Detaching vdata "//fieldname(fbyte1:lbyte1)//" - vsfdtch"
      CALL ckstatus_f(vrtn,usrlog,funcname,LOGFLAG)

c rhucek 05/10/01:  Added attributes for "Solution_Index"
*/ Vdata Solution_Index
      fieldname="Solution_Index"

c Getting reference number for vdata
      vsref=vsffnd(hdfid,fieldname)

      rtn = string_loc(fieldname, fbyte1, lbyte1)
      usrlog="Getting reference number for vdata "//
     &        fieldname(fbyte1:lbyte1)//" - vsffnd"
      CALL ckstatus_f(vsref,usrlog,funcname,LOGFLAG)


c Attaching vdata
      vsid=vsfatch(hdfid,vsref,"w")

      usrlog="Attaching vdata "//
     &        fieldname(fbyte1:lbyte1)//" - vsfatch"
      CALL ckstatus_f(vsid,usrlog,funcname,LOGFLAG)


      attr="long_name"
      longname="Solution Index for Ocean Best Small (1-4) "
     +//"and Large (5-9) Modes"//char(0)
      nms=strlen(longname)
      numbertype=DFNT_CHAR

      vrtn=vsfscat(vsid,-1,attr,numbertype,nms,longname)

      usrlog="Setting long name for "//
     &        fieldname(fbyte1:lbyte1)//" - vsfscat"
      CALL ckstatus_f(vrtn,usrlog,funcname,LOGFLAG)


      attr="units"
      sds_unit="None"//char(0)
      nms=strlen(sds_unit)
      numbertype=DFNT_CHAR

      vrtn=vsfscat(vsid,-1,attr,numbertype,nms,sds_unit)

      usrlog="Setting units for "//
     &        fieldname(fbyte1:lbyte1)//" - vsfscat"
      CALL ckstatus_f(vrtn,usrlog,funcname,LOGFLAG)


      attr="Geolocation_Pointer"
      geoptr="Geolocation data not applicable"//char(0)
      nms=strlen(geoptr)
      numbertype=DFNT_CHAR

      vrtn=vsfscat(vsid,-1,attr,numbertype,nms,geoptr)

      usrlog="Setting geolocation pointer for "//
     &        fieldname(fbyte1:lbyte1)//" - vsfscat"
      CALL ckstatus_f(vrtn,usrlog,funcname,LOGFLAG)


c Detaching vdata Solution_Index
      vrtn=vsfdtch(vsid)

      usrlog="Detaching vdata "//fieldname(fbyte1:lbyte1)//" - vsfdtch"
      CALL ckstatus_f(vrtn,usrlog,funcname,LOGFLAG)



*/ Releasing resources for vdata interface
      vrtn=vfend(hdfid)

      usrlog="Releasing resources for vdata interface - vfend"
      CALL ckstatus_f(vrtn,usrlog,funcname,LOGFLAG)



*/ Write the vdata

      CALL write_vdata(fid,swathname,swathid,
     &                       "Solution_1_Land",2,land_1)
      CALL write_vdata(fid,swathname,swathid,
     &                       "Solution_2_Land",3,land_2)
      CALL write_vdata(fid,swathname,swathid,
     &                       "Solution_3_Land",3,land_3)
       CALL write_vdata(fid,swathname,swathid,
     &                       "Solution_4_Land",4,land_4)
      CALL write_vdata(fid,swathname,swathid,
     &                       "Solution_Ocean",2,ocean_data)
       CALL write_vdata(fid,swathname,swathid,
     &                       "Solution_Index",9,ocean_data_index)
      CALL write_vdata(fid,swathname,swathid,
     &                       "MODIS_Band_Land",7,band_land)
      CALL write_vdata(fid,swathname,swathid,
     &                       "MODIS_Band_Ocean",7,band_ocean)



*/ Detach from the swath
      srtn=swdetach(swathid)

      usrlog="Detaching from the swath - swdetach"
      CALL ckstatus_f(srtn,usrlog,funcname,LOGFLAG)



*/  Close the output file.
      srtn=swclose(fid)

      rtn = string_loc(fname_sw, fbyte1, lbyte1)
      usrlog="Closing the output file "//char(10)//
     &       fname_sw(fbyte1:lbyte1)//" - swclose"
      CALL ckstatus_f(srtn,usrlog,funcname,LOGFLAG)



*/  Return to the calling routine.
      create_mod04=SUCCEED

      return
      end





***/----------------------------------------------------------------/***
      subroutine write_vdata(fid,swathname,swathid,fieldname,num,data)

C !F77------------------------------------------------------------------
C
C !Description: This routine writes Vdata objects values for every
C               fieldname passed in.
C
C !Input Parameters:
C
C    integer         fid            File ID of the swath file
C    character*(*)   swathname      Swath name in character string
C    integer         swathid        Swath ID
C    character*(*)   fieldname      Name for the Vdata in character
C                                   string
C    integer         num            Number of values for the Vdata
C    integer(*)      data           Buffer for the values for the Vdata
C    character*512   funcname       Character string for function names
C    character*512   usrlog         Character string for log messages
C
C
C !Output Parameters: None
C
C !Revision History:
C
C !Team-Unique Header:
C
C   This software was developed by the MODIS Science Data Support Team
C   (SDST) for the National Aeronautics and Space Administration,
C   Goddard Space Flight Center, under contract NAS5-32373.
C
C
C !References and Credits:
C
C   Developed by JC Guu 03/10/97 jguu@ltpmail.gsfc.nasa.gov
C   Modified by Fay Liang & Rich Hucek September 1997
C   fhliang@ltpmail.gsfc.nasa.gov & rhucek@ltpmail.gsfc.nasa.gov
C
C !Externals:
C
C   Functions:
C      swdetach           SWapi.c
C      swattach           SWapi.c
C      swwrfld            SWapi.c
C
C !Internals:
C      integer       srtn          return code from library routines
C
C !END------------------------------------------------------------------


      implicit none

      include 'ckstatus.inc'


      character*(*) swathname,fieldname
      integer fid,swathid,num
      integer swdetach,swattach,swwrfld,srtn
      integer*2 data(*)

      character*512 funcname,usrlog
      parameter (funcname = "write_vdata")

      srtn=swdetach(swathid)

      usrlog="Detaching from the swath - swdetach"
      CALL ckstatus_f(srtn,usrlog,funcname,LOGFLAG)


      swathid=swattach(fid,swathname)

      usrlog="Attaching to the swath - swattach"
      CALL ckstatus_f(swathid,usrlog,funcname,LOGFLAG)


      srtn=swwrfld(swathid,fieldname,0,1,num,data)

      usrlog="Writing data to swath field - swwrfld"
      CALL ckstatus_f(srtn,usrlog,funcname,LOGFLAG)



      end





