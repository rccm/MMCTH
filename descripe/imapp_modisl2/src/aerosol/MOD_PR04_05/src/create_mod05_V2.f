      integer function create_mod05()

C!F77-------------------------------------------------------------------
C
C!Description: This routine defines the swath components for MOD_PR05.
C
C!Input Parameters: None
C
C!Output Parameters: None
C
C!Revision History:
C $Log: create_mod05.f,v $
C
C!Team-Unique Header:
C This software was developed by the MODIS Science Data Support Team
C (SDST) for the National Aeronautics and Space Administration,
C Goddard Space Flight Center, under contract NAS5-32373.
C
C HDF functions were developed at the National Center for Supercomputing
C Applications at the Univ. of Illinois at Urbana-Champaign.
C
C!References and Credits:
C   Developed by JC Guu 03/10/97 jguu@ltpmail.gsfc.nasa.gov
C   Modified by Fay Liang fhliang@ltpmail.gsfc.nasa.gov
C           and Rich Hucek rhucek@ltpmail.gsfc.nasa.gov
C
C!Externals (non-SDST code):
C   Functions:
C    pgs_pc_getreference     (PGS_PC_GetReference.c)
C    ehidinfo                (SWapi.c)
C    swattach                (SWapi.c)
C    swclose                 (SWapi.c)
C    swcreate                (SWapi.c)
C    swdefdfld               (SWapi.c)
C    swdefdim                (SWapi.c)
C    swdefmap                (SWapi.c)
C    swdetach                (SWapi.c)
C    swdiminfo               (SWapi.c)
C    swopen                  (SWapi.c)
C    sfdimid                 (dffunc.inc)
C    sfscattr
C    sfsnattr
C    sfsdmname               (dffunc.inc)
C
C  Named Constants:
C    MAPIOK                  (mapi.inc)
C    MODFILLEN               (mapi.inc)
C    MODIS_E_GENERIC         (PGS_MODIS_39500.f)
C    MODIS_S_GENERIC         (PGS_MODIS_39500.f)
C    PGSd_PC_FILE_PATH_MAX   (PGS_PC.f)
C    FAIL                    (hdf.inc)
C    SUCCEED                 (hdf.inc)
C
C!Internals (SDST code):
C   Variables:
C    acswath_1km      Across swath sampling for 1 km resolution
C    acswath_5km      Across swath sampling for 5 km resolution
C    alswath_1km      Along swath sampling for 1 km resolution
C    alswath_5km      Along swath sampling for 5 km resolution
C    aoffset          Added offset
C    attrname         Attribute name
C    descrip          Description
C    dimsize          Dimemsion value
C    dimid            Dimension ID
C    dimlist_1km      Dimension list for 1 km resolution, 2 dimensions
C                     would be listed in it
C    dimlist_5km      Dimension list for 5 km resolution, 2 dimensions
C                     would be listed in it
C    dimlist_NIR      Dimension list for 1 km resolution, 3 dimensions
C                     would be listed in it, including 'QA_Bytes_NIR'
C    dimlist_IR       Dimension list for 5 km resolution, 3 dimensions
C                     would be listed in it, including 'QA_Bytes_IR'
C    dimname          Dimension name
C    dtype            Data type of the array
C    fbyte            1st byte of nonblank text of string
C    fid_05           HDF-EOS file ID for mod05
C    fid_07           HDF-EOS file ID for mod07
C    file_version     Version number for the file
C    fillval_f32      Fill value for 32-bit floating point type
C    fillval_f64      Fill value for 64-bit floating point type
C    fillval_i16      Fill value for 16-bit signed integer type
C    fillval_i8       Fill value for 8-bit signed integer type
C    fname_geo        File name of the input Geolocation file
C    fname_sw         File name of the output swath file
C    funcname         Character string for function names
C    gpointer         Geolocation pointer
C    hdfid_05         HDF File ID for mod05
C    hdfid_07         HDF File ID for mod07
C    hrtn             Return value from an HDF function
C    istr             String used in WRITE statement & usrlog
C    lbyte            Last byte of nonblank text of string
C    longname         Long name of SDS
C    merge            Merge code
C    modfil_geo       Array used to reference the geolocation file 
C                     iin M-API routines
C    modfil_sw        Array used to reference the swath file in
C                     all other M-API routines
C    mrtn             Return value from a M-API function
C    nframe           Number of frames possible
C    numbertype       Data type code
C    numscan          Number of scan in granule
C    num_val          Number of values for an attribute
C    prtn             Return value from an SDP Toolkit function
C    ptype            Parameter type
C    rank             Number of dimensions in an SDS
C    rtn              Return code from a function
C    sd_id_05         SDS interface ID for mod05
C    sd_id_07         SDS interface ID for mod07
C    sds_id_05        SDS ID in mod05
C    sds_id_07        SDS ID in mod07
C    sds_name_05      SDS array name of mod05
C    sds_name_07      SDS array name of mod07
C    sds_unit         Data unit
C    sfactor          Scaling factor
C    sl               Actual string length
C    srtn             Return code from Swath API
C    swathid_05       Swath structure ID for mod05
C    swathname_05     Swath structure name for mod05
C    swathname_07     Swath structure name for mod07
C    usrlog           Character string for log messages
C    val_i32          Buffer to hold global attribute valus
C    vrange_f32       Valid range for 32-bit floating point type
C    vrange_f64       Valid range for 64-bit floating point type
C    vrange_i16       Valid range for 16-bit signed integer type
C    vrange_i32       Valid range for 32-bit signed integer type
C    vrange_i8        Valid range for 8-bit signed integer type
C
C   Named Constants:
C    HDFE_NOMERGE     Code for 'no merge' action
C    HDFE_AUTOMERGE   Code for 'merge' action
C
C   Functions and Subroutines:
C    addtext                  (atmos shared code)
C    Copy_Swath_SDS           (atmos shared code)
C    create_array_sht_mod05   (Utility_V2.f)
C    modis_smf_setdynamicmsg  (atmos shared code)
C    set_array_attr_mod05     (Utility_V2.f)
C    string_loc               (atmos shared code)
C    strlen                   (atmos shared code)
C
C!END-------------------------------------------------------------------

      implicit none

      include 'mapi.inc'
      include 'hdf.inc'
      include 'PGS_SMF.f'
      include 'PGS_PC.f'
      include 'PGS_MODIS_39500.f'
      include 'mod05_hdfeos.inc'


*/ external functions
      integer pgs_pc_getreference
      integer ehidinfo, swattach, swclose, swcreate, swdefdfld
      integer swdefdim, swdefmap, swdetach, swdiminfo, swopen
      integer sfsnatt, sfscatt


*/ internal functions
      integer Copy_Swath_SDS 
      integer create_array_sht_mod05, set_array_attr_mod05
      integer create_array_byte
      integer string_loc, strlen, addtext


*/ variables for SDP toolkit routines

c rhucek 01/22/98: replaced fname_l1b with fname_geo
c     character*(PGSd_PC_FILE_PATH_MAX) fname_l1b,fname_sw,
      character*(PGSd_PC_FILE_PATH_MAX) fname_geo,fname_sw,
     +                                  fname_mod07
      integer file_version, prtn


*/ variables for M-API routines
      character*80 dtype

c rhucek 01/22/98:  replaced modfil_l1b with modfil_geo
c     integer modfil_l1b(MODFILLEN), modfil_sw(MODFILLEN),
      integer modfil_geo(MODFILLEN), modfil_sw(MODFILLEN),
     +        modfil_mod07(MODFILLEN)
      integer mrtn


*/ variables for HDF-EOS routines
      integer fid_05, swathid_05, hdfid_05, sd_id_05
      integer fid_07, swathid_07, hdfid_07, sd_id_07
      integer numbertype, merge, numscan, nframe
      integer srtn
      integer rank

      character*256 sds_name_05, sds_name_07
      character*256 swathname_05, swathname_07
      character*256 dimlist_1km, dimlist_5km, dimlist_IR, dimlist_NIR
      character*256 attrname, dimname, dimlist, longname, ptype,
     +              gpointer, sds_unit
      character*2048 descrip 
      character*1024 text 
      integer val_i32

      double precision fillval_f64,vrange_f64(2),sfactor,aoffset
      real fillval_f32,vrange_f32(2)
      integer vrange_i32(2)
      integer acswath_1km(3),alswath_1km(3),acswath_5km(3),
     +        alswath_5km(3)
      integer*2 fillval_i16,vrange_i16(2)
      byte fillval_i8,vrange_i8(2)


*/ variables for M-API routines & HDF-EOS routines
      integer num_val


*/ variables for native HDF routines
      integer sds_id_05, sds_index_05, hrtn, dimid

*/ variables for locally develpoed code.
      character*512 funcname, usrlog
      parameter (funcname = "create_mod05")
      character*10 istr
      integer dimsize
      integer rtn, sl
      integer fbyte, lbyte

*/ data array initilization (Rechcek the increment ????)
      data acswath_1km/1,1354,1/
      data alswath_1km/1,2030,1/
      data acswath_5km/1,1350,5/
      data alswath_5km/1,2000,5/


      create_mod05 = SUCCEED

*/ get the filename of the input geolocation file.
      file_version = 1
      prtn = pgs_pc_getreference(LUN_geo, file_version, fname_geo)

c rhucek 01/22/98: replaced LUN_l1b with LUN_geo
c     WRITE (istr, '(I6)') LUN_l1b
      WRITE (istr, '(I6)') LUN_geo

c rhucek 12/18/97:  commented 5 lines and added new "IF" line
c      IF (prtn .EQ. pgs_s_success) THEN
c        usrlog="Retrieving filename for LUN "//istr//
c     +       " - pgs_pc_getreference"
c        CALL modis_smf_setdynamicmsg(modis_s_generic, usrlog, funcname)
c      ELSE

      IF (prtn .NE. pgs_s_success) THEN
        usrlog="Error retrieving filename for LUN "//istr//
     +       " - pgs_pc_getreference"
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
        CALL modis_smf_setdynamicmsg(modis_e_generic, usrlog, funcname)
        create_mod05 = FAIL
        RETURN
      ENDIF


*/  Opening Geolocation  file
c rhucek 01/22/98: replaced fname_l1b and modfil_l1b with 
c fname_geo and modfil_geo
c     mrtn = opmfil(fname_l1b, "r", modfil_l1b)
      mrtn = opmfil(fname_geo, "r", modfil_geo)

c     rtn = string_loc(fname_l1b, fbyte, lbyte)
      rtn = string_loc(fname_geo, fbyte, lbyte)

c rhucek 12/18/97:  commented 5 lines and added new "IF" line
c      IF (mrtn .EQ. MAPIOK) THEN
c        usrlog = "Opening input file "//fname_geo(fbyte:lbyte)//
c     +         " - opmfil"
c        CALL modis_smf_setdynamicmsg(modis_s_generic, usrlog, funcname)
c      ELSE

      IF (mrtn .NE. MAPIOK) THEN
        usrlog="Error opening input file "//fname_geo(fbyte:lbyte)//
     +         " - opmfil"
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
        CALL modis_smf_setdynamicmsg(modis_e_generic, usrlog, funcname)
        create_mod05 = FAIL
        RETURN
      ENDIF


*/  Get attribute "Number of Scans" in granule.
      attrname ="Number of Scans"
      dtype = "INTEGER*4"
      num_val = 1
      mrtn = gmfin(modfil_geo, attrname, dtype, num_val, numscan)

c rhucek 12/18/97:  commented 6 lines and added new "IF" line
c      IF (mrtn .EQ. MAPIOK) THEN
c        WRITE (istr, '(I6)') numscan
c        usrlog="Getting numscan " //istr// "from "
c     +          //fname_geo(fbyte:lbyte)// " - gmfin"
c        CALL modis_smf_setdynamicmsg(modis_s_generic, usrlog, funcname)
c      ELSE

      IF (mrtn .NE. MAPIOK) THEN
        usrlog="Error getting numscan from "
     +          //fname_geo(fbyte:lbyte)// " - gmfin"
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
        CALL modis_smf_setdynamicmsg(modis_e_generic, usrlog, funcname)
        create_mod05 = FAIL
        RETURN
      ENDIF
      

c rhucek 01/23/98: changed "Max Earth View Frames" to "Max Earth Frames"
*/  Get the attribute of "Max Earth Frames" of the granule.
      attrname ="Max Earth Frames"
      mrtn=gmfin(modfil_geo, attrname, dtype, num_val, nframe)

c rhucek 12/18/97:  commented 6 lines and added new "IF" line
c      IF (mrtn .EQ. MAPIOK) THEN
c        WRITE (istr, '(I6)') nframe
c        usrlog="Getting nframe " //istr// "from "
c     +          //fname_geo(fbyte:lbyte)// " - gmfin"
c        CALL modis_smf_setdynamicmsg(modis_s_generic, usrlog, funcname)
c      ELSE

      IF (mrtn .NE. MAPIOK) THEN
        usrlog="Error getting nframe from "
     +          //fname_geo(fbyte:lbyte)// " - gmfin"
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
        CALL modis_smf_setdynamicmsg(modis_e_generic, usrlog, funcname)
        create_mod05 = FAIL
        RETURN
      ENDIF


*/  Close the Geolocation file.
      mrtn=clmfil(modfil_geo)

c rhucek 12/18/97:  commented 5 lines and added new "IF" line
c      IF (mrtn .EQ. MAPIOK) THEN
c        usrlog="Closing input file "//fname_geo(fbyte:lbyte)//
c     +         " - clmfil"
c        CALL modis_smf_setdynamicmsg(modis_s_generic, usrlog, funcname)
c      ELSE

      IF (mrtn .NE. MAPIOK) THEN
        usrlog="Error closing input file "//fname_geo(fbyte:lbyte)//
     +         " - clmfil"
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
        CALL modis_smf_setdynamicmsg(modis_e_generic, usrlog, funcname)
        create_mod05 = FAIL
        RETURN
      ENDIF



*/ retrieve file name of output swath file mod05
      file_version = 1
      prtn = pgs_pc_getreference(LUN_sw, file_version, fname_sw)

      WRITE (istr, '(I6)') LUN_sw

c rhucek 12/18/97:  commented 5 lines and added new "IF" line
c      IF (prtn .EQ. pgs_s_success) THEN
c        usrlog="Retrieving filename formod05, LUN "//istr//
c     +       " - pgs_pc_getreference"
c        CALL modis_smf_setdynamicmsg(modis_s_generic, usrlog, funcname)
c      ELSE

      IF (prtn .NE. pgs_s_success) THEN
        usrlog="Error retrieving filename for mod05, LUN "//istr//
     +       " - pgs_pc_getreference"
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
        CALL modis_smf_setdynamicmsg(modis_e_generic, usrlog, funcname)
        create_mod05 = FAIL
      ENDIF


*/ open the output swath file
      fid_05 = swopen(fname_sw, DFACC_CREATE)

      If (fid_05 .gt. 0) Then
         usrlog='Opened mod05 swath product for create - swopen.'
         CALL modis_smf_setdynamicmsg(modis_s_generic, usrlog, funcname)
      Else
         usrlog='Failed to open MOD05 product file.'
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
         CALL modis_smf_setdynamicmsg(modis_e_generic, usrlog,
     1        funcname)
         create_mod05 = FAIL
      End If

*/ If unable to open MOD05 product file, return FAIL
      If (create_mod05 .eq. FAIL) RETURN

*/ create a swath named 'mod05'
      swathname_05 ="mod05"
      swathid_05 = swcreate(fid_05, swathname_05)


*/ get SDS interface ID and HDF file ID for the created swath mod05
      srtn = ehidinfo(fid_05, hdfid_05, sd_id_05)


*/ retrieve file name of input swath file mod07
      file_version=1
      prtn=pgs_pc_getreference(LUN_mod07,file_version,fname_mod07)

      WRITE (istr, '(I6)') LUN_mod07

c rhucek 12/18/97:  commented 5 lines and added new "IF" line
c      IF (prtn .EQ. pgs_s_success) THEN
c        usrlog="Retrieving filename for mod07, LUN "//istr//
c     +       " - pgs_pc_getreference"
c        CALL modis_smf_setdynamicmsg(modis_s_generic, usrlog, funcname)
c      ELSE

      IF (prtn .NE. pgs_s_success) THEN
        usrlog="Error retrieving filename for mod07, LUN "//istr//
     +       " - pgs_pc_getreference"
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
        CALL modis_smf_setdynamicmsg(modis_e_generic, usrlog, funcname)
        create_mod05 = FAIL
        RETURN
      ENDIF


*/ open swath mod07 for READ
      fid_07 = swopen(fname_mod07, DFACC_READ)

      IF (fid_07 .EQ. FAIL) THEN
        rtn = string_loc(fname_mod07, fbyte, lbyte)
        usrlog="Error opening swath file "//fname_mod07(fbyte:lbyte)//
     +         " - swopen"
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
        CALL modis_smf_setdynamicmsg(modis_e_generic, usrlog, funcname)
        create_mod05 = FAIL
        RETURN
      ENDIF


*/ get mod07's swathid
      swathid_07 = swattach(fid_07, 'mod07')

*/ Get the SDS interface ID and HDF file ID for sawth mod07
      srtn = ehidinfo(fid_07, hdfid_07, sd_id_07)


*/ Define dimensions for data and geo. fields ----------------
      dimname="Cell_Along_Swath_1km"
      dimsize=numscan*10
      srtn=swdefdim(swathid_05, dimname, dimsize)


      dimname="Cell_Across_Swath_1km"
      dimsize=nframe
      srtn=swdefdim(swathid_05, dimname, dimsize)


c following 3 dimensions from/of mod07 to be defined in swath_id_05
c retrieve size of dimensions from mod07
      dimname="Cell_Along_Swath_5km"
      dimsize = swdiminfo(swathid_07, 'Cell_Along_Swath')

      IF (dimsize .EQ. FAIL) THEN
        rtn = string_loc(dimname, fbyte, lbyte)
        usrlog="Error retrieving dimsize for "//dimname(fbyte:lbyte)//
     +         " - swdiminfo"
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
        CALL modis_smf_setdynamicmsg(modis_e_generic, usrlog, funcname)
        create_mod05 = FAIL
        RETURN
      ELSE
        srtn=swdefdim(swathid_05, dimname, dimsize)
      ENDIF


      dimname="Cell_Across_Swath_5km"
      dimsize = swdiminfo(swathid_07, 'Cell_Across_Swath')

      IF (dimsize .EQ. FAIL) THEN
        rtn = string_loc(dimname, fbyte, lbyte)
        usrlog="Error retrieving dimsize for "//dimname(fbyte:lbyte)//
     +         " - swdiminfo"
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
        CALL modis_smf_setdynamicmsg(modis_e_generic, usrlog, funcname)
        create_mod05 = FAIL
        RETURN
      ELSE
        srtn=swdefdim(swathid_05, dimname, dimsize)
      ENDIF


      dimname="QA_Bytes_IR"
      dimsize = swdiminfo(swathid_07, 'Water_Vapor_QA_Bytes')

      IF (dimsize .EQ. FAIL) THEN
        rtn = string_loc(dimname, fbyte, lbyte)
        usrlog="Error retrieving dimsize for "//dimname(fbyte:lbyte)//
     +         " - swdiminfo"
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
        CALL modis_smf_setdynamicmsg(modis_e_generic, usrlog, funcname)
        create_mod05 = FAIL
        RETURN
      ELSE
        srtn=swdefdim(swathid_05, dimname, dimsize)
      ENDIF


      dimname="QA_Bytes_NIR"
      dimsize=1
      srtn=swdefdim(swathid_05, dimname, dimsize)



*/ create dimension maps relating 1 km data to 5 km geolocation
      srtn = swdefmap(swathid_05, 'Cell_Across_Swath_5km',
     +                'Cell_Across_Swath_1km', 2, 5)

      srtn = swdefmap(swathid_05, 'Cell_Along_Swath_5km',
     +                'Cell_Along_Swath_1km', 2, 5)


*/ define dimlists to be used in Copy_Swath_SDS()
      dimlist_1km = 'Cell_Across_Swath_1km,Cell_Along_Swath_1km'
      dimlist_5km = 'Cell_Across_Swath_5km,Cell_Along_Swath_5km'
      dimlist_NIR = 'QA_Bytes_NIR,' //
     +             'Cell_Across_Swath_1km,Cell_Along_Swath_1km'
      dimlist_IR = 'QA_Bytes_IR,' //
     +             'Cell_Across_Swath_5km,Cell_Along_Swath_5km'



*/ define global attributes of swath mod05 using native HDF routines
      attrname = 'Maximum_Number_of_1km_Frames'
      numbertype = DFNT_INT32
      num_val = 1
      val_i32 = nframe 

      hrtn=sfsnatt(sd_id_05, attrname, numbertype, num_val, val_i32)


      attrname = 'Number_of_Instrument_Scans'
      numbertype = DFNT_INT32
      num_val = 1
      val_i32 = numscan 

      hrtn = sfsnatt(sd_id_05, attrname, numbertype, num_val, val_i32)


      attrname = "title"
      numbertype = DFNT_CHAR
      descrip = char(10)//
     +         " MODIS HDF File Specification MOD05_L2: MODIS Level"//
     +         " 2 NIR and IR Total Precipitable Water "//char(10)
      sl = strlen(descrip)
      num_val = sl

      srtn=sfscatt(sd_id_05, attrname, numbertype, num_val, descrip)


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
       num_val = sl
       srtn = sfscatt(sd_id_05, attrname, numbertype, num_val, descrip)

c rhucek 01/22/98:  Commented "history" file attribute in accordance with
c                   CCR 376 
c     attrname = "history"
c     numbertype = DFNT_CHAR
c     descrip = char(10)//
c    +         " 4-15-1997   Allen Chu/GSFC                         "//
c    +         "                         "//char(10)//
c    +         "             Kathy Strabale/Univ. of Wisconsin      "//
c    +         "                         "//char(10)
c     sl = strlen(descrip)
c     num_val = sl

c     srtn = sfscatt(sd_id_05, attrname, numbertype, num_val, descrip)



*/ define dimension QA_Parameter_NIR

C      sds_name_05="QA_Parameter_NIR"
C      numbertype=DFNT_INT16
C      sds_id_05=sfcreate(sd_id_05,sds_name_05,numbertype,1,1)

C      dimid=sfdimid(sds_id_05,0)
C      hrtn=sfsdmname(dimid, sds_name_05)



*/ define dimension QA_Parameter_IR

C      sds_name_05="QA_Parameter_IR"
C      numbertype=DFNT_INT16
C      sds_id_05=sfcreate(sd_id_05,sds_name_05,numbertype,1,4)

C      dimid=sfdimid(sds_id_05,0)
C      hrtn=sfsdmname(dimid, sds_name_05)



*/ copy data & geolocation fields & set attributes from mod07

*/ copying 'Scan_Start_Time'
      sds_name_07 = 'Scan_Start_Time'
      sds_name_05 = 'Scan_Start_Time'
      dimlist = dimlist_5km
      rtn = Copy_Swath_SDS(sd_id_07, sd_id_05, sds_name_07, sds_name_05,
     +                     fid_05, swathid_05, swathname_05, dimlist)

*/ copying 'Latitude'
      sds_name_07 = 'Latitude'
      sds_name_05 = 'Latitude'
      dimlist = dimlist_5km
      rtn = Copy_Swath_SDS(sd_id_07, sd_id_05, sds_name_07, sds_name_05,
     +                     fid_05, swathid_05, swathname_05, dimlist)

*/ copying 'Longitude'
      sds_name_07 = 'Longitude'
      sds_name_05 = 'Longitude'
      dimlist = dimlist_5km
      rtn = Copy_Swath_SDS(sd_id_07, sd_id_05, sds_name_07, sds_name_05,
     +                     fid_05, swathid_05, swathname_05, dimlist)

*/ copying 'Solar_Zenith'
      sds_name_07 = 'Solar_Zenith'
      sds_name_05 = 'Solar_Zenith'
      dimlist = dimlist_5km
      rtn = Copy_Swath_SDS(sd_id_07, sd_id_05, sds_name_07, sds_name_05,
     +                     fid_05, swathid_05, swathname_05, dimlist)

*/ copying 'Solar_Azimuth'
      sds_name_07 = 'Solar_Azimuth'
      sds_name_05 = 'Solar_Azimuth'
      dimlist = dimlist_5km
      rtn = Copy_Swath_SDS(sd_id_07, sd_id_05, sds_name_07, sds_name_05,
     +                     fid_05, swathid_05, swathname_05, dimlist)

*/ copying 'Sensor_Zenith'
      sds_name_07 = 'Sensor_Zenith'
      sds_name_05 = 'Sensor_Zenith'
      dimlist = dimlist_5km
      rtn = Copy_Swath_SDS(sd_id_07, sd_id_05, sds_name_07, sds_name_05,
     +                     fid_05, swathid_05, swathname_05, dimlist)

*/ copying 'Sensor_Azimuth'
      sds_name_07 = 'Sensor_Azimuth'
      sds_name_05 = 'Sensor_Azimuth'
      dimlist = dimlist_5km
      rtn = Copy_Swath_SDS(sd_id_07, sd_id_05, sds_name_07, sds_name_05,
     +                     fid_05, swathid_05, swathname_05, dimlist)


*/ create 2-D array 'Cloud_Mask_QA'
      sds_name_05="Cloud_Mask_QA"
      dimlist = dimlist_1km
      numbertype=DFNT_INT8
      merge=HDFE_NOMERGE
      longname="MODIS Cloud Mask, First Byte"
      sds_unit="none"
      vrange_i8(1)= 0
      vrange_i8(2)= -1
      fillval_i8= 0
      sfactor=1.0
      aoffset=0.0
      ptype="MODIS Input"
      gpointer = "External MODIS geolocation product"

      descrip = char(10) //
     1"                                                                "
     2// " " // char(10) // 
     3" bit field   Description                     Key                "
     4// " " // char(10) //
     5" ---------   -----------                     ---                "
     6// " " // char(10) //
     7"                                                                "
     8// " " //char(10) //
     9" 0           Cloud Mask Flag                 0 = not determined "
     +// " " //char(10) //
     1"                                             1 = determined     "
     2// " " //char(10) //
     3"                                                                "
     4// " " //char(10) //
     5" 1-2         Unobstructed FOV Quality Flag   00 = cloud         "
     6// " " //char(10) //
     7"                                             01 = 66% prob.     "
     8// " " //char(10)

      text = 
     1"                                                    clear       "
     2// " " //char(10) //
     3"                                             10 = 95% prob.     "
     4// " " //char(10) //
     5"                                                    clear       "
     6// " " //char(10) //
     7"                                             11 = 99% prob.     "
     8// " " //char(10) //
     9"                                                    clear       "
     +// " " //char(10)

c-----Append text to descrip
      rtn = addtext(text, descrip)

      text =
     1"             PROCESSING PATH FLAGS                              "
     2// " " //char(10) //
     3"             ---------------------                              "
     4// " " //char(10) //
     5" 3           Day / Night Flag                0 = Night / 1 = Day"
     6// " " //char(10) //
     7" 4           Sunglint Flag                   0 = Yes / 1 = No   "
     8// " " //char(10) //
     9" 5           Snow / Ice Background Flag      0 = Yes / 1 = No   "
     +// " " //char(10) //
     1" 6-7         Land / Water Flag               00=Water/01=Coastal"
     2// " " //char(10) //
     3"                                             10=Desert/11=Land  "
     4// " " //char(10)

c-----Append text to descrip
      rtn = addtext(text, descrip)

      rtn=create_array_byte(sd_id_05,swathid_05, sds_name_05, dimlist,
     +                     numbertype, merge, longname, sds_unit,
     +                     vrange_i8, fillval_i8, sfactor, aoffset,
     +                     ptype, gpointer, descrip,
     +                     acswath_1km, alswath_1km)


*/ create 2-D array 'Water_Vapor_Near_Infrared'
      sds_name_05="Water_Vapor_Near_Infrared"
      dimlist = dimlist_1km
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Total Column Precipitable Water Vapor - Near Infrared"//
     +         " Retrieval"
      sds_unit="cm"
      vrange_i16(1)= 0
      vrange_i16(2)= 20000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      gpointer = "Internal geolocation arrays"
      descrip=" "

      rtn=create_array_sht_mod05(sd_id_05,swathid_05, sds_name_05, 
     +                     dimlist,
     +                     numbertype, merge, longname, sds_unit,
     +                     vrange_i16, fillval_i16, sfactor, aoffset,
     +                     ptype, gpointer, descrip,
     +                     acswath_1km, alswath_1km)

*/ create 2-D array 'Water_Vapor_CORRECTION_FACTORS'
      sds_name_05="Water_Vapor_Correction_Factors "
      dimlist = dimlist_1km
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Aerosol Correction Factor for Water Vapor - Near Infrared"//
     +         " Retrieval"
      sds_unit="none"
      vrange_i16(1)= 0
      vrange_i16(2)= 20000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      gpointer = "Internal geolocation arrays"
      descrip=" "

      rtn=create_array_sht_mod05(sd_id_05,swathid_05, sds_name_05,
     +                     dimlist,
     +                     numbertype, merge, longname, sds_unit,
     +                     vrange_i16, fillval_i16, sfactor, aoffset,
     +                     ptype, gpointer, descrip,
     +                     acswath_1km, alswath_1km)


*/ copying 'Water_Vapor_Infrared'
      sds_name_07 = 'Water_Vapor'
      sds_name_05 = 'Water_Vapor_Infrared'
      dimlist = dimlist_5km
      rtn = Copy_Swath_SDS(sd_id_07, sd_id_05, sds_name_07, sds_name_05,
     +                     fid_05, swathid_05, swathname_05, dimlist)


*/ create 3-D array 'Quality_Assurance_Near_Infrared'
      sds_name_05="Quality_Assurance_Near_Infrared"
      dimlist = dimlist_NIR
      numbertype=DFNT_INT8
      merge=HDFE_NOMERGE
      longname="Run time QA flags"
      sds_unit="none"
      vrange_i8(1)= 0
      vrange_i8(2)= -1
      fillval_i8= 0
      sfactor=1.0
      aoffset=0.0
      ptype="Output"
      gpointer = "External MODIS geolocation product"
      descrip="See MODIS atmosphere QA plan for details"

      rtn=create_array_byte(sd_id_05,swathid_05, sds_name_05, dimlist,
     +                     numbertype, merge, longname, sds_unit,
     +                     vrange_i8, fillval_i8, sfactor, aoffset,
     +                     ptype, gpointer, descrip,
     +                     acswath_1km, alswath_1km)


*/ copying 'Quality_Assurance_Infrared'
      sds_name_07 = 'Quality_Assurance_Infrared'
      sds_name_05 = 'Quality_Assurance_Infrared'
      dimlist = dimlist_IR
      rtn = Copy_Swath_SDS(sd_id_07, sd_id_05, sds_name_07, sds_name_05,
     +                     fid_05, swathid_05, swathname_05, dimlist)



*/ detach from the swath mod05
      srtn = swdetach(swathid_05)


*/ close the swath file mod05.
      srtn = swclose(fid_05)

      IF (srtn .EQ. SUCCEED) THEN
        usrlog="Closed mod05 swath product after create - swclose."
        CALL modis_smf_setdynamicmsg(modis_s_generic, usrlog,funcname)
      ELSE
        usrlog="Error closing swath mod05 - swclose"
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
        CALL modis_smf_setdynamicmsg(modis_e_generic, usrlog,funcname)
        create_mod05 = FAIL
        RETURN
      ENDIF



*/ detach from the mod07 swath
      srtn = swdetach(swathid_07)


*/ close the mod07 swath
      srtn = swclose(fid_07)

      IF (srtn .EQ. SUCCEED) THEN
        usrlog="Closing mod07 swath product - swclose"
        CALL modis_smf_setdynamicmsg(modis_s_generic, usrlog,funcname)
      ELSE
        usrlog="Error closing swath mod07 - swclose"
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
        CALL modis_smf_setdynamicmsg(modis_e_generic, usrlog,funcname)
        create_mod05 = FAIL
        RETURN
      ENDIF



*/ return to the calling routine
      create_mod05 = SUCCEED

      RETURN

      END



******/----------------------------------------------------------/******

      integer function create_array_byte(sd_id_05,swathid_05,
     +        sds_name_05,
     +        dimlist, numbertype, merge, longname, sds_unit, vrange,
     +        fillval, sfactor, aoffset, ptype, gpointer, descrip,
     +        acswath, alswath)

C!F77-------------------------------------------------------------------
C
C!Description:
C   This routine creates an 8-bit-signed-integer (byte) array using
C   HDF-EOS Swath API and calls set_array_attr() and HDF native routines
C   to set the attributes of the SDS.
C
C!Input Parameters:
C   acswath      Across swath sampling
C   alswath      Along swath sampling 
C   aoffset      Added offset
C   descrip      Description
C   dimlist      Dimension list
C   fillval      Fill value
C   gpointer     Geolocation pointer
C   longname     Long name of SDS
C   merge        Merge code
C   numbertype   Data type code
C   ptype        Parameter type
C   sds_name_05  SDS array name
C   sds_unit     Data unit
C   sd_id_05     SDS interface ID
C   sfactor      Scaling factor
C   swathid_05   Swath ID
C   vrange       Valid range
C
C
C!Output Parameters: None
C
C!Revision History:
C $Log: create_mod05_V2.f,v $
C
C!Team-Unique Header:
C   This software was developed by the MODIS Science Data Support Team
C   (SDST) for the National Aeronautics and Space Administration,
C   Goddard Space Flight Center, under contract NAS5-32373.
C
C HDF functions were developed at the National Center for Supercomputing
C Applications at the Univ. of Illinois at Urbana-Champaign.
C
C!References and Credits:
C   Developed by JC Guu jguu@ltpmail.gsfc.nasa.gov 03/10/97
C   Modified by Fay Liang fhliang@ltpmail.gsfc.nasa.gov
C           and Rich Hucek rhucek@ltpmail.gsfc.nasa.gov
C
C!Externals (non-SDST code):
C  Functions:
C    swdefdfld               (SWapi.c)
C    swsetfill               (SWapi.c)
C    sfscattr                (dffunc.inc)
C    sfsnattr                (dffunc.inc)
C
C  Named Constants:
C   MODIS_E_GENERIC         (PGS_MODIS_39500.f)
C   FAIL                    (hdf.inc)
C   SUCCEED                 (hdf.inc)
C
C!Internals (SDST code):
C  Variables:
C   attr       Character string for attribute names
C   fbyte      1st byte of nonblank text of string
C   funcname   Character string for fnction names
C   hrtn       Return code from HDF function
C   lbyte      last byte of nonblank text of string
C   num_val    Number of values for an attribute
C   rtn        Return code from a function
C   sdsidx     SDS index
C   sdsid      SDS ID
C   srtn       Return code from Swath API
C   usrlog     Character string for log messages
C
C  Functions & Subroutines:
C   set_array_attr  (Utility_V2.f)
C   string_loc      (atmos shared code)
C   strlen          (atmos shared code)
C
C!END-------------------------------------------------------------------

      IMPLICIT NONE

      include 'PGS_MODIS_39500.f'
      include 'hdf.inc'
      include 'dffunc.inc'


*/ passed-in arguments
      character*(*) sds_name_05, dimlist, longname, sds_unit
      character*(*) ptype, gpointer, descrip
      integer swathid_05, numbertype, merge
      integer sd_id_05, acswath(*), alswath(*)
      byte vrange(*), fillval
      double precision sfactor, aoffset

*/ functions
      integer swdefdfld, swsetfill
      integer set_array_attr_mod05, strlen, string_loc

*/ variables
      character*256 attr
      character*512 funcname, usrlog
      integer srtn, hrtn, rtn, num_val
      integer fbyte, lbyte
      integer sdsidx, sdsid

      data sdsidx, sdsid/-1,-1/


      funcname = "create_array_byte"

      rtn = string_loc(sds_name_05, fbyte, lbyte)

*/ define the data field
      srtn = swdefdfld(swathid_05, sds_name_05(fbyte:lbyte),
     +                 dimlist, numbertype, merge)

      IF (srtn .EQ. FAIL) THEN
        usrlog="Error defining data field "//sds_name_05(fbyte:lbyte)//
     +         " - swdefdfld"
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
        CALL modis_smf_setdynamicmsg(modis_e_generic, usrlog,funcname)
        create_array_byte = FAIL
        RETURN
      ENDIF


*/ retrieve the SDS index for the data field
      sdsidx = sfn2index(sd_id_05, sds_name_05)

      IF (sdsidx .EQ. FAIL) THEN
        usrlog = 
     1  'MOD05 HDF procedure sfn2index returned FAIL for array ' 
     2  // sds_name_05(fbyte:lbyte)
     3  // char(10) // 'Operator Action:  Check system resources/environment, '
     4  // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     5  // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
        CALL modis_smf_setdynamicmsg(MODIS_E_GENERIC, usrlog, funcname)
        create_array_byte = FAIL
        RETURN
      ENDIF


*/ retrieve the SDS ID for the data field
      sdsid = sfselect(sd_id_05, sdsidx)


*/ set the attributes common to all arrays throughout the swath
      rtn = set_array_attr_mod05(sdsid, longname, sds_unit, sfactor,
     +      aoffset, ptype, gpointer, descrip, acswath, alswath)


*/ set the valid range of the data field
      attr="valid_range"
      num_val=2

      hrtn=sfsnatt(sdsid,attr,numbertype,num_val,vrange)


*/ set the fillvalue of the data field
      srtn=swsetfill(swathid_05, sds_name_05(fbyte:lbyte), fillval)


*/ successfully return to the calling routine
      create_array_byte=SUCCEED

      RETURN

      END
