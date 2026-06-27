      integer function add_db_mod04()

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
      
c  Number of models
       integer num_index_models,LRN_L1B_1km,LRN_L1B_1km_RA
       parameter(num_index_models=9)
       PARAMETER (LRN_L1B_1km=700002,LRN_L1B_1km_RA=430001)

*/  Functions
      integer pgs_pc_getreference
      integer swopen,swcreate,swdefdim,swdefdfld,swwrdmeta,swdetach,
     &        swsetfill,swwrfld,swclose,ehidinfo,swattach

      integer  num_args
      integer  FlagRA
      character FlagBuff*10
      integer  iargc

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
C			integer sds_id, sds_index
      integer sdsidx

*/  Variables for locally develpoed code.
      character*512 funcname,usrlog
      parameter (funcname = "add_db_mod04")
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
      data acswath/5,1345,10/
      data alswath/5,2025,10/


*----------------------------------------------------------------------*

*/  Get the filenames of the input L1B file.
      file_version=1
      if( FlagRA .eq. 1) then
         prtn=pgs_pc_getreference(LRN_L1B_1km_RA,file_version,fname_l1b)
      else
         prtn=pgs_pc_getreference(LRN_L1B_1km,file_version,fname_l1b)
      endif

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
     2   // char(10) // "Operator Action:  Check system "//
     2         "resources/environment, "
     3   // char(10) // "PCF, and SDPTK configuration.  If a fault "//
     2         "is identified, "
     4   // char(10) // "correct and rerun PGE.  Otherwise, "//
     2         "notify SDST."
        CALL modis_smf_setdynamicmsg(modis_e_generic, usrlog, funcname)
        add_db_mod04 = FAIL
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
c       CALL modis_smf_setdynamicmsg(modis_s_generic, usrlog, funcname)
      ELSE
        usrlog="Error closing input file "//fname_l1b(fbyte2:lbyte2)//
     +         " - clmfil"
     2   // char(10) // "Operator Action:  Check system "//
     2         "resources/environment, "
     3   // char(10) // "PCF, and SDPTK configuration.  If a fault "//
     2         "is identified, "
     4   // char(10) // "correct and rerun PGE.  Otherwise, notify "//
     2         "SDST."
        CALL modis_smf_setdynamicmsg(modis_e_generic, usrlog, funcname)
        add_db_mod04 = FAIL
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
      fid=swopen(fname_sw,DFACC_RDWR)

      rtn = string_loc(fname_sw, fbyte1, lbyte1)
      usrlog="Opening output file "//char(10)//
     &       fname_sw(fbyte1:lbyte1)//" with HDF-EOS"//" - swopen"
      CALL ckstatus_f(fid,usrlog,funcname,LOGFLAG)

*/  Attach to swath.
      swathname="mod04"
      swathid=swattach(fid,swathname)

      rtn = string_loc(swathname, fbyte1, lbyte1)
      usrlog="Attaching to swath named "//swathname(fbyte1:lbyte1)//
     &       " - swattach"
      CALL ckstatus_f(swathid,usrlog,funcname,LOGFLAG)

*/  Get the SDS interface ID and HDF file ID for the created sawth.
      srtn=ehidinfo(fid,hdfid,sdid)

      usrlog="Obtaining the HDF ID and SDS interface ID - ehidinfo"
      CALL ckstatus_f(srtn,usrlog,funcname,LOGFLAG)

      dimname="Num_DeepBlue_Wavelengths"
      dimn=3
      srtn=swdefdim(swathid,dimname,dimn)

      usrlog="Defining dimension Num_DeepBlue_Wavelengths - swdefdim"
      CALL ckstatus_f(srtn,usrlog,funcname,LOGFLAG)

*/---------  Define data fields and set the attributes.  -------------/*

      fieldname="Deep_Blue_Aerosol_Optical_Depth_550_Land"
      dimlist="Cell_Across_Swath,Cell_Along_Swath"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="AOT at 0.55 micron for "//
     &      "land with all quality data (QA=1,2,3)"
      sds_unit="None"
      vrange_i16(1)= 0
      vrange_i16(2)= 5000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"
	 
      sdsidx = sfn2index(sdid , fieldname)
		  IF (sdsidx .lt. 0) THEN
				print *, 'Adding sds: ',fieldname
      	crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)
    
				rtn = string_loc(fieldname, fbyte1, lbyte1)
				usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
				CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)
			ENDIF
			
			fieldname="Deep_Blue_Aerosol_Optical_Depth_550_Land_Best_Estimate"
      dimlist="Cell_Across_Swath,Cell_Along_Swath"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="AOT at 0.55 micron for "//
     &      "land with QA=3 only"
      sds_unit="None"
      vrange_i16(1)= 0
      vrange_i16(2)= 5000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"
	 
      sdsidx = sfn2index(sdid , fieldname)
		  IF (sdsidx .lt. 0) THEN
				print *, 'Adding sds: ',fieldname
      	crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)
				rtn = string_loc(fieldname, fbyte1, lbyte1)
				usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
				CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)
			ENDIF
			
      fieldname="Deep_Blue_Spectral_Aerosol_Optical_Depth_Land"
      dimlist="Cell_Across_Swath,Cell_Along_Swath,Num_DeepBlue_Wavelengths"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="AOT at 0.412, 0.47, and 0.66 micron for "//
     &      "land with all quality data (QA=1,2,3)"
      sds_unit="None"
      vrange_i16(1)= 0
      vrange_i16(2)= 5000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      sdsidx = sfn2index(sdid , fieldname)
		  IF (sdsidx .lt. 0) THEN
				print *, 'Adding sds: ',fieldname
      	crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)
    
				rtn = string_loc(fieldname, fbyte1, lbyte1)
				usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
				CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)
			ENDIF
      
      
c   New SDS 2/1/2006       
      
      fieldname="Deep_Blue_Angstrom_Exponent_Land"
      dimlist="Cell_Across_Swath,Cell_Along_Swath"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Deep Blue Angstrom Exponent for "//
     &    "land (0.412-0.47 micron) with all quality data (QA=1,2,3)"
      sds_unit="None"
      vrange_i16(1)= -500
      vrange_i16(2)= 5000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      sdsidx = sfn2index(sdid , fieldname)
		  IF (sdsidx .lt. 0) THEN
				print *, 'Adding sds: ',fieldname
      	crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)
    
				rtn = string_loc(fieldname, fbyte1, lbyte1)
				usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
				CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)
			ENDIF

      fieldname="Deep_Blue_Spectral_Single_Scattering_Albedo_Land"
      dimlist="Cell_Across_Swath,Cell_Along_Swath,Num_DeepBlue_Wavelengths"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Deep Blue Single Scattering Albedo at 0.412, 0.47, and 0.66 micron for "//
     &    "land with all quality data (QA=1,2,3)"
      sds_unit="None"
      vrange_i16(1)= 700
      vrange_i16(2)= 1000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      sdsidx = sfn2index(sdid , fieldname)
		  IF (sdsidx .lt. 0) THEN
				print *, 'Adding sds: ',fieldname
      	crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)
    
				rtn = string_loc(fieldname, fbyte1, lbyte1)
				usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
				CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)
			ENDIF

      fieldname="Deep_Blue_Spectral_Surface_Reflectance_Land"
      dimlist="Cell_Across_Swath,Cell_Along_Swath,Num_DeepBlue_Wavelengths"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Deep Blue Surface Reflectance at 0.412, 0.47, and 0.66 micron for land "//
     &    "with all quality data (QA=1,2,3)"
      sds_unit="None"
      vrange_i16(1)= 0
      vrange_i16(2)= 1000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      sdsidx = sfn2index(sdid , fieldname)
		  IF (sdsidx .lt. 0) THEN
				print *, 'Adding sds: ',fieldname
      	crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)
    
				rtn = string_loc(fieldname, fbyte1, lbyte1)
				usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
				CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)
			ENDIF

c   New SDS 3/11/2008      

      fieldname="Deep_Blue_Spectral_TOA_Reflectance_Land"
      dimlist="Cell_Across_Swath,Cell_Along_Swath,Num_DeepBlue_Wavelengths"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Average measured top-of-atmosphere reflectance after cloud screening at "//
     &    "0.412, 0.47, and 0.66 micron for land"
      sds_unit="None"
      vrange_i16(1)= 0
      vrange_i16(2)= 10000        
      fillval_i16= -9999
      sfactor=0.0001              ! different that the other SDS's but correct!
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      sdsidx = sfn2index(sdid , fieldname)
		  IF (sdsidx .lt. 0) THEN
				print *, 'Adding sds: ',fieldname
      	crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)
    
				rtn = string_loc(fieldname, fbyte1, lbyte1)
				usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
				CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)
			ENDIF
 
      fieldname="Deep_Blue_Number_Pixels_Used_550_Land"
      dimlist="Cell_Across_Swath,Cell_Along_Swath"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Number of pixels used for AOT retrievals at 0.55 micron for land"
      sds_unit="None"
      vrange_i16(1)= 0
      vrange_i16(2)= 100
      fillval_i16= -999
      sfactor=1.0
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      sdsidx = sfn2index(sdid , fieldname)
		  IF (sdsidx .lt. 0) THEN
				print *, 'Adding sds: ',fieldname
      	crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)
    
				rtn = string_loc(fieldname, fbyte1, lbyte1)
				usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
				CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)
			ENDIF
     
      fieldname="Deep_Blue_Aerosol_Optical_Depth_550_Land_STD"
      dimlist="Cell_Across_Swath,Cell_Along_Swath"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Standard deviation of AOT at 0.55 micron for "//
     &      "land with all quality data (QA=1,2,3)"
      sds_unit="None"
      vrange_i16(1)= 0
      vrange_i16(2)= 10000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      sdsidx = sfn2index(sdid , fieldname)
		  IF (sdsidx .lt. 0) THEN
				print *, 'Adding sds: ',fieldname
      	crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)
    
				rtn = string_loc(fieldname, fbyte1, lbyte1)
				usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
				CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)
			ENDIF
      
c      fieldname="Deep_Blue_Aerosol_Optical_Depth_Land_STD"
c      dimlist="Cell_Across_Swath,Cell_Along_Swath,Num_DeepBlue_Wavelengths"
c      numbertype=DFNT_INT16
c      merge=HDFE_NOMERGE
c      longname="Standard deviation of AOT at 0.412, 0.47, and 0.66 micron for "//
c     &      "land with all quality data (Quality flag=1,2,3)"
c      sds_unit="None"
c      vrange_i16(1)= 0
c      vrange_i16(2)= 10000
c      fillval_i16= -9999
c      sfactor=0.001
c      aoffset=0.0
c      ptype="Output"
c      descrip=" "
c      geoptr="Internal geolocation arrays"
c
c      sdsidx = sfn2index(sdid , fieldname)
c		  IF (sdsidx .lt. 0) THEN
c				print *, 'Adding sds: ',fieldname
c      	crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
c     &                    numbertype,merge,longname,sds_unit,vrange_i16,
c     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
c     &                    alswath,descrip,geoptr)
c    
c				rtn = string_loc(fieldname, fbyte1, lbyte1)
c				usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
c				CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)
c			ENDIF

      fieldname="Deep_Blue_Algorithm_Flag_Land"
      dimlist="Cell_Across_Swath,Cell_Along_Swath"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Algorithm used to perform retrievals (0=DeepBlue,1=Vegetation,2=Both)"
      sds_unit="None"
      vrange_i16(1)= 0
      vrange_i16(2)= 2
      fillval_i16= -999
      sfactor=1.0
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"
	 
      sdsidx = sfn2index(sdid , fieldname)
		  IF (sdsidx .lt. 0) THEN
				print *, 'Adding sds: ',fieldname
      	crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)
    
				rtn = string_loc(fieldname, fbyte1, lbyte1)
				usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
				CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)
			ENDIF
			
			fieldname="Deep_Blue_Aerosol_Optical_Depth_550_Land_QA_Flag"
      dimlist="Cell_Across_Swath,Cell_Along_Swath"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Quality assurance (QA) flag for Deep Blue aerosol optical thickness value at 0.55 micron "//
     &           "(0=No confidence/Fill value, 1=Marginal, 2=Good, 3=Very good)"
      sds_unit="None"
      vrange_i16(1)= 0
      vrange_i16(2)= 3
      fillval_i16= -999
      sfactor=1.0
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"
	 
      sdsidx = sfn2index(sdid , fieldname)
		  IF (sdsidx .lt. 0) THEN
				print *, 'Adding sds: ',fieldname
      	crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)
    
				rtn = string_loc(fieldname, fbyte1, lbyte1)
				usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
				CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)
			ENDIF
			
			fieldname="Deep_Blue_Aerosol_Optical_Depth_550_Land_Estimated_Uncertainty"
      dimlist="Cell_Across_Swath,Cell_Along_Swath"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Estimated uncertainty per cell for Deep Blue aerosol optical depth" //
     &  " over land for all QA."
      sds_unit="None"
      vrange_i16(1)= 0
      vrange_i16(2)= 500
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"
	 
      sdsidx = sfn2index(sdid , fieldname)
		  IF (sdsidx .lt. 0) THEN
				print *, 'Adding sds: ',fieldname
      	crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)
    
				rtn = string_loc(fieldname, fbyte1, lbyte1)
				usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
				CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)
			ENDIF
			
			fieldname="Deep_Blue_Cloud_Fraction_Land"
      dimlist="Cell_Across_Swath,Cell_Along_Swath"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Deep Blue cloud fraction over land"
      sds_unit="None"
      vrange_i16(1)= 0
      vrange_i16(2)= 1000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"

      sdsidx = sfn2index(sdid , fieldname)
		  IF (sdsidx .lt. 0) THEN
				print *, 'Adding sds: ',fieldname
      	crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)
    
				rtn = string_loc(fieldname, fbyte1, lbyte1)
				usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
				CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)
			ENDIF
			
! New for C006 - Combined Product

      fieldname="Dark_Target_Deep_Blue_Optical_Depth_550_Land_And_Ocean"
      fieldname="AOD_550_Dark_Target_Deep_Blue_Combined"
      dimlist="Cell_Across_Swath,Cell_Along_Swath"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Aerosol optical depth at 0.55 micron from Dark Target " // 
     &         "over ocean using QA=1,2,3 and combined Dark Target, " // 
     &         "Deep Blue over land with best QA (DB_QA>=2, DT_QA=3)"
      sds_unit="None"
      vrange_i16(1)= 0
      vrange_i16(2)= 5000
      fillval_i16= -9999
      sfactor=0.001
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"
	 
      sdsidx = sfn2index(sdid , fieldname)
		  IF (sdsidx .lt. 0) THEN
				print *, 'Adding sds: ',fieldname
      	crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)
    
				rtn = string_loc(fieldname, fbyte1, lbyte1)
				usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
				CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)
			ENDIF

      fieldname="Dark_Target_Deep_Blue_Optical_Depth_550_Land_And_Ocean_QA"
      fieldname="AOD_550_Dark_Target_Deep_Blue_Combined_QA_Flag"
      dimlist="Cell_Across_Swath,Cell_Along_Swath"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Quality assurance (QA) flag of aerosol optical depth at 0.55 micron from Dark Target " //
     &          "over ocean and combined Dark Target, Deep Blue over land " //
     &          "(0=No confidence, 1=Marginal, 2=Good, 3=Very good)"
      sds_unit="None"
      vrange_i16(1)= 0
      vrange_i16(2)= 3
      fillval_i16= -999
      sfactor=1.0
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"
	 
      sdsidx = sfn2index(sdid , fieldname)
		  IF (sdsidx .lt. 0) THEN
				print *, 'Adding sds: ',fieldname
      	crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)
    
				rtn = string_loc(fieldname, fbyte1, lbyte1)
				usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
				CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)
			ENDIF
			
			fieldname="Dark_Target_Deep_Blue_Optical_Depth_550_Land_And_Ocean_AlgFlag"
			fieldname="AOD_550_Dark_Target_Deep_Blue_Combined_Algorithm_Flag"
      dimlist="Cell_Across_Swath,Cell_Along_Swath"
      numbertype=DFNT_INT16
      merge=HDFE_NOMERGE
      longname="Algorithm flag indicating source of Dark_Target_Deep_Blue_Optical_Depth_550_Land_And_Ocean " //
     &          "data (0=Dark Target, 1=Deep Blue, 2=Mixed)"
      sds_unit="None"
      vrange_i16(1)= 0
      vrange_i16(2)= 2
      fillval_i16= -999
      sfactor=1.0
      aoffset=0.0
      ptype="Output"
      descrip=" "
      geoptr="Internal geolocation arrays"
	 
      sdsidx = sfn2index(sdid , fieldname)
		  IF (sdsidx .lt. 0) THEN
				print *, 'Adding sds: ',fieldname
      	crtn=create_array_sht(sdid,swathid,fieldname,dimlist,
     &                    numbertype,merge,longname,sds_unit,vrange_i16,
     &                    fillval_i16,sfactor,aoffset,ptype,acswath,
     &                    alswath,descrip,geoptr)
    
				rtn = string_loc(fieldname, fbyte1, lbyte1)
				usrlog=hlog//fieldname(fbyte1:lbyte1)//" - create_array_sht"
				CALL ckstatus_f(crtn,usrlog,funcname,LOGFLAG)
			ENDIF
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
      add_db_mod04=SUCCEED

      return
      end
