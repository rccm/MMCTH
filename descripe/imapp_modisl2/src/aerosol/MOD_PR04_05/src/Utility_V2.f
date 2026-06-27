      integer function set_array_attr_mod05(sdsid, longname, sds_unit,
     +        sfactor, aoffset, ptype, gpointer, descrip, acswath,
     +        alswath)

C!F77-------------------------------------------------------------------
C
C!Description: This routine sets the attributes common to all SDS.
C
C!Input Parameters:
C   acswath      Across swath sampling
C   alswath      Along swath sampling
C   aoffset      Added offset
C   descrip      Description
C   gpointer     Geolocation pointer
C   longname     Long name of SDS
C   ptype        Parameter type
C   sds_unit     Data unit
C   sdsid        SDS ID
C   sfactor      Scaling factor
C
C
C!Output Parameters: None
C
C!Revision History:
C $Log: Utility_V2.f,v $
C
C!Team-Unique Header:
C
C   This software was developed by the MODIS Science Data Support Team
C   (SDST) for the National Aeronautics and Space Administration,
C   Goddard Space Flight Center, under contract NAS5-32373.
C
C HDF functions were developed at the National Center for Supercomputing
C Applications at the Univ. of Illinois at Urbana-Champaign.
C
C!References and Credits:
C
C   Developed by JC Guu 03/10/97 jguu@ltpmail.gsfc.nasa.gov
C   Modified by Fay Liang fhliang@ltpmail.gsfc.nasa.gov
C           and Rich Hucek rhucek@ltpmail.gsfc.nasa.gov
C
C!Externals (non-SDST code):
C   Functions:
C    sfscattr
C    sfsnattr
C
C  Named Constants:
C
C!Internals: (SDST code)
C   Variables:
C   attr             Attribute name
C   hrtn             Return value from an HDF function
C   numbertype       Data type code
C   num_val          Number of values for an attribute
C
C   Named Constants: None
C
C   Functions and Subroutines:
C    strlen          (atmos shared code)
C
C!END-------------------------------------------------------------------

      implicit none

      include 'hdf.inc'


*/ Passed-in arguments

      character*(*) longname, sds_unit
      character*(*) ptype, gpointer, descrip
      integer sdsid, acswath(*), alswath(*)
      double precision sfactor, aoffset

*/ Functions

      integer strlen
      integer sfscatt, sfsnatt

*/ Variables

      character*256 attr
      integer hrtn, num_val, numbertype


*/ Attributes "long_name","sds_unit","scale_factor","add_offset",
*/ "Parameter_Type","Cell_Along_Swath_Sampling",
*/ "Cell_Across_Swath_Sampling" and "description" are set.
*/ For some data fields, "Cell_Across_Swath_Sampling" and/or "description"
*/ are not applicable, therefore an IF construct is provided to see if
*/ these attributes are needed by the data field in question.

      attr="long_name"
      numbertype=DFNT_CHAR8
      num_val=strlen(longname)
      hrtn=sfscatt(sdsid,attr,numbertype,num_val,longname)


      attr="unit"
      numbertype=DFNT_CHAR8
      num_val=strlen(sds_unit)
      hrtn=sfscatt(sdsid,attr,numbertype,num_val,sds_unit)


      attr="scale_factor"
      numbertype=DFNT_FLOAT64
      num_val=1
      hrtn=sfsnatt(sdsid,attr,numbertype,num_val,sfactor)


      attr="add_offset"
      numbertype=DFNT_FLOAT64
      num_val=1
      hrtn=sfsnatt(sdsid,attr,numbertype,num_val,aoffset)


      attr="Parameter_Type"
      numbertype=DFNT_CHAR8
      num_val=strlen(ptype)
      hrtn=sfscatt(sdsid,attr,numbertype,num_val,ptype)


      attr="Cell_Along_Swath_Sampling"
      numbertype=DFNT_INT32
      num_val=3
      hrtn=sfsnatt(sdsid,attr,numbertype,num_val,alswath)


      IF (acswath(1) .GT. 0) THEN
          attr="Cell_Across_Swath_Sampling"
          numbertype=DFNT_INT32
          num_val=3
          hrtn=sfsnatt(sdsid,attr,numbertype,num_val,acswath)
      ENDIF


      attr = "Geolocation_Pointer"
      numbertype = DFNT_CHAR
      num_val = strlen(gpointer)

      hrtn = sfscatt(sdsid, attr, numbertype, num_val, gpointer)


      IF (descrip(1:1) .NE. " ") THEN
          attr="description"
          numbertype=DFNT_CHAR8
          num_val=strlen(descrip)
          hrtn=sfscatt(sdsid,attr,numbertype,num_val,descrip)
      ENDIF


*/ return to the calling routine
      set_array_attr_mod05=SUCCEED

      RETURN

      END

******/---------------------------------------------------------------/******

      integer function create_array_sht_mod05(sd_id_05,swathid_05,
     +        sds_name_05,
     +        dimlist, numbertype, merge, longname, sds_unit, vrange,
     +        fillval, sfactor, aoffset, ptype, gpointer, descrip,
     +        acswath, alswath)

C!F77-------------------------------------------------------------------
C
C!Description:
C   This routine creates a 16-bit-signed-integer (short) array using
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
C $Log: Utility_V2.f,v $
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
C    sfscattr
C    sfsnattr
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
      integer*2 vrange(*), fillval
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


      funcname = "create_array_sht_mod05"

      rtn = string_loc(sds_name_05, fbyte, lbyte)

*/ define the data field
      srtn = swdefdfld(swathid_05, sds_name_05(fbyte:lbyte),
     +                 dimlist, numbertype, merge)

      IF (srtn .EQ. FAIL) THEN
        usrlog="Error defining data field "//sds_name_05(fbyte:lbyte)//
     1         " - swdefdfld"
     2  // char(10) // 'Operator Action:  Check system resources/environment, '
     3  // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4  // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
        CALL modis_smf_setdynamicmsg(MODIS_E_GENERIC, usrlog, funcname)
        create_array_sht_mod05 = FAIL
        RETURN
      ENDIF


*/ retrieve the SDS index for the data field
      sdsidx = sfn2index(sd_id_05, sds_name_05)

      IF (sdsidx .EQ. FAIL) THEN
        usrlog = 'HDF procedure sfn2index returned FAIL for array ' //
     1           sds_name_05(fbyte:lbyte)
     2  // char(10) // 'Operator Action:  Check system resources/environment, '
     3  // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4  // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
        CALL modis_smf_setdynamicmsg(MODIS_E_GENERIC, usrlog, funcname)
        create_array_sht_mod05 = FAIL
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
      create_array_sht_mod05=SUCCEED

      RETURN

      END

******/---------------------------------------------------------------/******

      integer function create_array_dbl(sdid,swathid,fieldname,
     &                                  dimlist,numbertype,merge,
     &                                  longname,datau,vrange,fillval,
     &                                  sfactor,aoffset,ptype,acswath,
     &                                  alswath,descrip,geoptr)

C----------------------------------------------------------------------
C !F77
C
C !Description:
C
C   This routine creates a double-precision array using HDF-EOS Swath
C   API and calls set_array_attr() and HDF native routines to set the
C   attributes of the SDS.
C
C !Input Parameters:
C
C   Type               Name         Description
C   ----               ----         -----------
C   integer            sdid         SDS interface ID
C   integer            swathid      Swath ID
C   character*(*)      fieldname    SDS array name
C   character*(*)      dimlist      Dimension list
C   integer            numbertype   Data type code
C   integer            merge        Merge code
C   character*(*)      longname     Long name of SDS
C   character*(*)      datau        Units for the SDS data
C   double precision   vrange(*)    Valid range
C   double precision   fillval      Fill value
C   double precision   sfactor      Scaling factor
C   double precision   aoffset      Added offset
C   character*(*)      ptype        Parameter type
C   integer            acswath(*)   Across swath sampling
C   integer            alswath(*)   Along swath sampling
C   character*(*)      descrip      Description
C
C
C !Output Parameters: None
C
C !Revision History:
C $Log: Utility_V2.f,v $
c 01/29/98 fhliang
c added NCSA acknowledgement.
C
C !Team-Unique Header:
C
C   This software was developed by the MODIS Science Data Support Team
C   (SDST) for the National Aeronautics and Space Administration,
C   Goddard Space Flight Center, under contract NAS5-32373.
C
c HDF portions developed at the National Center for Supercomputing
c Applications at the University of Illinois at Urbana-Champaign.
c
C !References and Credits:
C
C   Developer: JC Guu 03/10/97 jguu@ltpmail.gsfc.nasa.gov
C
C !Externals:
C
C   Type            Name           Source               Description
C   ----            ----           ------               -----------
C   integer         SUCCEED        hdf.inc              Return code for
C                                                       succeed.
C
C !Internals:
C
C   Type              Name             Description
C   ----              ----             -----------
C   character*256     attr             Character string for attribute names
C   character*512     funcname         Character string for fnction names
C   character*512     usrlog           Character string for log messages
C   integer           srtn             Hold return code from Swath API
C   integer           hrtn             Hold return code from HDF.
C   integer           rtn              Hold generic return code.
C   integer           nms              Number of values for an attribute
C   integer           type             Type code
C   integer           sl               True string length
C   integer           sdsidx           SDS index
C   integer           sdsid            SDS ID
C
C
C !END
C----------------------------------------------------------------------

      implicit none

      include 'hdf.inc'
      include 'dffunc.inc'
      include 'ckstatus.inc'


*/  Passed-in arguments

      character*(*) fieldname,dimlist,longname,datau,ptype,descrip,
     &              geoptr
      integer swathid,numbertype,merge
      integer sdid,acswath(*),alswath(*)
      double precision vrange(*),fillval,sfactor,aoffset

*/  Functions

      integer swdefdfld,swsetfill,set_array_attr,strlen

*/  Variables

      character*256 attr
      character*512 funcname,usrlog
      integer srtn,hrtn,rtn,nms,type,sl
      integer sdsidx,sdsid

      data sdsidx,sdsid/-1,-1/

*---------------------------------------------------------------------*


*/  Define the data field.

      srtn=swdefdfld(swathid,fieldname,dimlist,numbertype,merge)

      funcname="create_array_dbl"
      sl=strlen(fieldname)
      usrlog="Creating array "//fieldname(1:sl)//
     &       " - create_array_dbl"
      call ckstatus_f(srtn,usrlog,funcname,LOGFLAG)


*/  Retrieve the SDS index for the data field.

      sdsidx=sfn2index(sdid,fieldname)

      usrlog="Retrieving SDS index for "//fieldname(1:sl)//
     &       " - sfn2index"
      call ckstatus_f(sdsidx,usrlog,funcname,LOGFLAG)


*/  Retrieve the SDS ID for the data field.

      sdsid=sfselect(sdid,sdsidx)

      usrlog="Retrieving SDS ID for "//fieldname(1:sl)//
     &       " - sfselect"
      call ckstatus_f(sdsid,usrlog,funcname,LOGFLAG)


*/  Set the attributes with same types throughout the swath.

      rtn=set_array_attr(sdsid,fieldname,longname,datau,sfactor,aoffset,
     &                   ptype,acswath,alswath,descrip,geoptr)

      usrlog="Setting attributes for "//fieldname(1:sl)//
     &       " - set_array_attr"
      call ckstatus_f(rtn,usrlog,funcname,LOGFLAG)


*/  The following code slice for setting the fill value is commented
*/  out because currently the FORTRAN HDF_EOS routine for assigning
*/  fill values, swsetfill, is not working properly in the way that
*/  a fill value can be assigned to a SDS but this fill value is
*/  not accessible to a data reader.


*/  Set the fill value of the data field.

      srtn=swsetfill(swathid,fieldname(1:sl),fillval)

      usrlog="Setting fillvalue for "//fieldname(1:sl)//
     &       " - swsetfill"
      call ckstatus_f(srtn,usrlog,funcname,LOGFLAG)


*/  Set the valid range of the data field.

      attr="valid_range"
      type=numbertype
      nms=2
      hrtn=sfsattr(sdsid,attr,type,nms,vrange)

      usrlog="Setting valid range for "//fieldname(1:sl)//
     &       " - sfsattr"
      call ckstatus_f(hrtn,usrlog,funcname,LOGFLAG)


*/  Successfully returns.

      create_array_dbl=SUCCEED

      return
      end


******/---------------------------------------------------------------/******

      integer function create_array_flt(sdid,swathid,fieldname,
     &                                  dimlist,numbertype,merge,
     &                                  longname,datau,vrange,fillval,
     &                                  sfactor,aoffset,ptype,acswath,
     &                                  alswath,descrip,geoptr)

C----------------------------------------------------------------------
C !F77
C
C !Description:
C
C   This routine creates a floating-point array using HDF-EOS Swath
C   API and calls set_array_attr() and HDF native routines to set the
C   attributes of the SDS.
C
C !Input Parameters:
C
C   Type               Name         Description
C   ----               ----         -----------
C   integer            sdid         SDS interface ID
C   integer            swathid      Swath ID
C   character*(*)      fieldname    SDS array name
C   character*(*)      dimlist      Dimension list
C   integer            numbertype   Data type code
C   integer            merge        Merge code
C   character*(*)      longname     Long name of SDS
C   character*(*)      datau        Units for the SDS data
C   real               vrange(*)    Valid range
C   real               fillval      Fill value
C   double precision   sfactor      Scaling factor
C   double precision   aoffset      Added offset
C   character*(*)      ptype        Parameter type
C   integer            acswath(*)   Across swath sampling
C   integer            alswath(*)   Along swath sampling
C   character*(*)      descrip      Description
C
C
C !Output Parameters: None
C
C !Revision History:
C $Log: Utility_V2.f,v $
c 01/29/98 fhliang
c added NCSA acknowledgement.
C
C !Team-Unique Header:
C
C   This software was developed by the MODIS Science Data Support Team
C   (SDST) for the National Aeronautics and Space Administration,
C   Goddard Space Flight Center, under contract NAS5-32373.
C
c HDF portions developed at the National Center for Supercomputing
c Applications at the University of Illinois at Urbana-Champaign.
c
C !References and Credits:
C
C   Developer: JC Guu 03/10/97 jguu@ltpmail.gsfc.nasa.gov
C
C !Externals:
C
C   Type            Name           Source               Description
C   ----            ----           ------               -----------
C   integer         SUCCEED        hdf.inc              Return code for
C                                                       succeed.
C
C !Internals:
C
C   Type              Name             Description
C   ----              ----             -----------
C   character*256     attr             Character string for attribute names
C   character*512     funcname         Character string for fnction names
C   character*512     usrlog           Character string for log messages
C   integer           srtn             Hold return code from Swath API
C   integer           hrtn             Hold return code from HDF.
C   integer           rtn              Hold generic return code.
C   integer           nms              Number of values for an attribute
C   integer           type             Type code
C   integer           sl               True string length
C   integer           sdsidx           SDS index
C   integer           sdsid            SDS ID
C
C
C !END
C----------------------------------------------------------------------
      implicit none

      include 'hdf.inc'
      include 'dffunc.inc'
      include 'ckstatus.inc'


*/  Passed-in arguments

      character*(*) fieldname,dimlist,longname,datau,ptype,descrip,
     &              geoptr
      integer swathid,numbertype,merge
      integer sdid,acswath(*),alswath(*)
      real vrange(*),fillval
      double precision sfactor,aoffset

*/  Functions

      integer swdefdfld,swdefgfld,swsetfill,set_array_attr,strlen,
     &        swgetfill

*/  Variables

      character*256 attr
      character*512 funcname,usrlog
      integer srtn,hrtn,rtn,nms,type,sl
      integer sdsidx,sdsid

      data sdsidx,sdsid/-1,-1/

*---------------------------------------------------------------------*


*/  Define the data or geo. field.

      if(fieldname(1:9).eq."Longitude" .or.
     &   fieldname(1:8).eq."Latitude") then
          srtn=swdefgfld(swathid,fieldname,dimlist,numbertype,merge)
      else
          srtn=swdefdfld(swathid,fieldname,dimlist,numbertype,merge)
      endif

      funcname="create_array_flt"
      sl=strlen(fieldname)
      usrlog="Creating array "//fieldname(1:sl)//
     &       " - swdef[dg]fld"
      call ckstatus_f(srtn,usrlog,funcname,LOGFLAG)


*/  Retrieve the SDS index for the data or geo. field.

      sdsidx=sfn2index(sdid,fieldname)

      usrlog="Retrieving SDS index for "//fieldname(1:sl)//
     &       " - sfn2index"
      call ckstatus_f(sdsidx,usrlog,funcname,LOGFLAG)


*/  Retrieve the SDS ID for the data or geo. field.

      sdsid=sfselect(sdid,sdsidx)

      usrlog="Retrieving SDS ID for "//fieldname(1:sl)//
     &       " - sfselect"
      call ckstatus_f(sdsid,usrlog,funcname,LOGFLAG)


*/  Set the attributes with same types throughout the swath.

      rtn=set_array_attr(sdsid,fieldname,longname,datau,sfactor,aoffset,
     &                   ptype,acswath,alswath,descrip,geoptr)

      usrlog="Setting attributes for "//fieldname(1:sl)//
     &       " - set_array_attr"
      call ckstatus_f(rtn,usrlog,funcname,LOGFLAG)


*/  The following code slice for setting the fill value is commented
*/  out because currently the FORTRAN HDF_EOS routine for assigning
*/  fill values, swsetfill, is not working properly in the way that
*/  a fill value can be assigned to a SDS but this fill value is
*/  not accessible to a data reader.


*/  Set the fillvalue of the data or geo. field.

      srtn=swsetfill(swathid,fieldname(1:sl),fillval)

      usrlog="Setting fillvalue for "//fieldname(1:sl)//
     &       " - swsetfill"
      call ckstatus_f(srtn,usrlog,funcname,LOGFLAG)


*/  Set the valid range of the data field.

      attr="valid_range"
      type=numbertype
      nms=2
      hrtn=sfsattr(sdsid,attr,type,nms,vrange)

      usrlog="Setting valid range for "//fieldname(1:sl)//
     &       " - sfsattr"
      call ckstatus_f(hrtn,usrlog,funcname,LOGFLAG)


*/  Successfully returns.

      create_array_flt=SUCCEED

      return
      end


******/---------------------------------------------------------------/******

      integer function create_array_lng(sdid,swathid,fieldname,
     &                                  dimlist,numbertype,merge,
     &                                  longname,datau,vrange,fillval,
     &                                  sfactor,aoffset,ptype,acswath,
     &                                  alswath,descrip,geoptr)

C----------------------------------------------------------------------
C !F77
C
C !Description:
C
C   This routine creates a long integer array using HDF-EOS Swath
C   API and calls set_array_attr() and HDF native routines to set the
C   attributes of the SDS.
C
C !Input Parameters:
C
C   Type               Name         Description
C   ----               ----         -----------
C   integer            sdid         SDS interface ID
C   integer            swathid      Swath ID
C   character*(*)      fieldname    SDS array name
C   character*(*)      dimlist      Dimension list
C   integer            numbertype   Data type code
C   integer            merge        Merge code
C   character*(*)      longname     Long name of SDS
C   character*(*)      datau         Data datau
C   integer            vrange(*)    Valid range
C   integer            fillval      Fill value
C   double precision   sfactor      Scaling factor
C   double precision   aoffset      Added offset
C   character*(*)      ptype        Parameter type
C   integer            acswath(*)   Across swath sampling
C   integer            alswath(*)   Along swath sampling
C   character*(*)      descrip      Description
C
C
C !Output Parameters: None
C
C !Revision History:
C $Log: Utility_V2.f,v $
c 01/29/98 fhliang
c added NCSA acknowledgement.
C
C !Team-Unique Header:
C
C   This software was developed by the MODIS Science Data Support Team
C   (SDST) for the National Aeronautics and Space Administration,
C   Goddard Space Flight Center, under contract NAS5-32373.
C
c HDF portions developed at the National Center for Supercomputing
c Applications at the University of Illinois at Urbana-Champaign.
c
C !References and Credits:
C
C   Developer: JC Guu 03/10/97 jguu@ltpmail.gsfc.nasa.gov
C
C !Externals:
C
C   Type            Name           Source               Description
C   ----            ----           ------               -----------
C   integer         SUCCEED        hdf.inc              Return code for
C                                                       succeed.
C
C !Internals:
C
C   Type              Name             Description
C   ----              ----             -----------
C   character*256     attr             Character string for attribute names
C   character*512     funcname         Character string for fnction names
C   character*512     usrlog           Character string for log messages
C   integer           srtn             Hold return code from Swath API
C   integer           hrtn             Hold return code from HDF.
C   integer           rtn              Hold generic return code.
C   integer           nms              Number of values for an attribute
C   integer           type             Type code
C   integer           sl               True string length
C   integer           sdsidx           SDS index
C   integer           sdsid            SDS ID
C
C
C !END
C----------------------------------------------------------------------
      implicit none

      include 'hdf.inc'
      include 'dffunc.inc'
      include 'ckstatus.inc'


*/ Passed-in arguments

      character*(*) fieldname,dimlist,longname,datau,ptype,descrip,
     &              geoptr
      integer swathid,numbertype,merge
      integer sdid,acswath(*),alswath(*)
      integer vrange(*),fillval
      double precision sfactor,aoffset

*/  Functions

      integer swdefdfld,swsetfill,set_array_attr,strlen


*/ Variables

      character*256 attr
      character*512 funcname,usrlog
      integer srtn,hrtn,rtn,nms,type,sl
      integer sdsidx,sdsid

      data sdsidx,sdsid/-1,-1/

*---------------------------------------------------------------------*


*/  Define the data field.

      srtn=swdefdfld(swathid,fieldname,dimlist,numbertype,merge)

      funcname="create_array_lng"
      sl=strlen(fieldname)
      usrlog="Creating array "//fieldname(1:sl)//
     &       " - swdefdfld"
      call ckstatus_f(srtn,usrlog,funcname,LOGFLAG)


*/  Retrieve the SDS index for the data field.

      sdsidx=sfn2index(sdid,fieldname)

      usrlog="Retrieving SDS index for "//fieldname(1:sl)//
     &       " - sfn2index"
      call ckstatus_f(sdsidx,usrlog,funcname,LOGFLAG)


*/  Retrieve the SDS ID for the data field.

      sdsid=sfselect(sdid,sdsidx)

      usrlog="Retrieving SDS ID for "//fieldname(1:sl)//
     &       " - sfselect"
      call ckstatus_f(sdsid,usrlog,funcname,LOGFLAG)


*/  Set the attributes with same types throughout the swath.

      rtn=set_array_attr(sdsid,fieldname,longname,datau,sfactor,aoffset,
     &                   ptype,acswath,alswath,descrip,geoptr)

      usrlog="Setting attributes for "//fieldname(1:sl)//
     &       " - set_array_attr"
      call ckstatus_f(rtn,usrlog,funcname,LOGFLAG)


*/  The following code slice for setting the fill value is commented
*/  out because currently the FORTRAN HDF_EOS routine for assigning
*/  fill values, swsetfill, is not working properly in the way that
*/  a fill value can be assigned to a SDS but this fill value is
*/  not accessible to a data reader.


*/  Set the fillvalue of the data field.

      srtn=swsetfill(swathid,fieldname(1:sl),fillval)

      usrlog="Setting fillvalue for "//fieldname(1:sl)//
     &       " - swsetfill"
      call ckstatus_f(srtn,usrlog,funcname,LOGFLAG)


*/  Set the valid range of the data field.

      attr="valid_range"
      type=numbertype
      nms=2
      hrtn=sfsattr(sdsid,attr,type,nms,vrange)

      usrlog="Setting valid range for "//fieldname(1:sl)//
     &       " - sfsattr"
      call ckstatus_f(hrtn,usrlog,funcname,LOGFLAG)


*/  Successfully returns.

      create_array_lng=SUCCEED

      return
      end


******/---------------------------------------------------------------/******

      integer function create_array_sht(sdid,swathid,fieldname,
     &                                  dimlist,numbertype,merge,
     &                                  longname,datau,vrange,fillval,
     &                                  sfactor,aoffset,ptype,acswath,
     &                                  alswath,descrip,geoptr)

C----------------------------------------------------------------------
C !F77
C
C !Description:
C
C   This routine creates a short integer array using HDF-EOS Swath
C   API and calls set_array_attr() and HDF native routines to set the
C   attributes of the SDS.
C
C !Input Parameters:
C
C   Type               Name         Description
C   ----               ----         -----------
C   integer            sdid         SDS interface ID
C   integer            swathid      Swath ID
C   character*(*)      fieldname    SDS array name
C   character*(*)      dimlist      Dimension list
C   integer            numbertype   Data type code
C   integer            merge        Merge code
C   character*(*)      longname     Long name of SDS
C   character*(*)      datau         Data datau
C   integer*2          vrange(*)    Valid range
C   integer*2          fillval      Fill value
C   double precision   sfactor      Scaling factor
C   double precision   aoffset      Added offset
C   character*(*)      ptype        Parameter type
C   integer            acswath(*)   Across swath sampling
C   integer            alswath(*)   Along swath sampling
C   character*(*)      descrip      Description
C
C
C !Output Parameters: None
C
C !Revision History:
C $Log: Utility_V2.f,v $
c 01/29/98 fhliang
c added NCSA acknowledgement.
C
C !Team-Unique Header:
C
C   This software was developed by the MODIS Science Data Support Team
C   (SDST) for the National Aeronautics and Space Administration,
C   Goddard Space Flight Center, under contract NAS5-32373.
C
c HDF portions developed at the National Center for Supercomputing
c Applications at the University of Illinois at Urbana-Champaign.
c
C !References and Credits:
C
C   Developer: JC Guu 03/10/97 jguu@ltpmail.gsfc.nasa.gov
C
C !Externals:
C
C   Type            Name           Source               Description
C   ----            ----           ------               -----------
C   integer         SUCCEED        hdf.inc              Return code for
C                                                       succeed.
C
C !Internals:
C
C   Type              Name             Description
C   ----              ----             -----------
C   character*256     attr             Character string for attribute names
C   character*512     funcname         Character string for fnction names
C   character*512     usrlog           Character string for log messages
C   integer           srtn             Hold return code from Swath API
C   integer           hrtn             Hold return code from HDF.
C   integer           rtn              Hold generic return code.
C   integer           nms              Number of values for an attribute
C   integer           type             Type code
C   integer           sl               True string length
C   integer           sdsidx           SDS index
C   integer           sdsid            SDS ID
C
C
C !END
C----------------------------------------------------------------------
      implicit none

      include 'hdf.inc'
      include 'dffunc.inc'
      include 'ckstatus.inc'


*/  Passed-in arguments

      character*(*) fieldname,dimlist,longname,datau,ptype,descrip,
     &              geoptr
      integer swathid,numbertype,merge
      integer sdid,acswath(*),alswath(*)
      integer*2 vrange(*),fillval
      double precision sfactor,aoffset

*/  Functions

      integer swdefdfld,swsetfill,set_array_attr,strlen

*/  Variables

      character*256 attr
      character*512 funcname,usrlog
      integer srtn,hrtn,rtn,nms,type,sl
      integer sdsidx,sdsid

      data sdsidx,sdsid/-1,-1/

*---------------------------------------------------------------------*


*/  Define the data field.

      srtn=swdefdfld(swathid,fieldname,dimlist,numbertype,merge)

      funcname="create_array_sht"
      sl=strlen(fieldname)
      usrlog="Creating array "//fieldname(1:sl)//
     &       " - swdefdfld"
      call ckstatus_f(srtn,usrlog,funcname,LOGFLAG)


*/  Retrieve the SDS index for the data field.

      sdsidx=sfn2index(sdid,fieldname)

      usrlog="Retrieving SDS index for "//fieldname(1:sl)//
     &       " - sfn2index"
      call ckstatus_f(sdsidx,usrlog,funcname,LOGFLAG)


*/  Retrieve the SDS ID for the data field.

      sdsid=sfselect(sdid,sdsidx)

      usrlog="Retrieving SDS ID for "//fieldname(1:sl)//
     &       " - sfselect"
      call ckstatus_f(sdsid,usrlog,funcname,LOGFLAG)


*/  Set the attributes with same types throughout the swath.

      rtn=set_array_attr(sdsid,fieldname,longname,datau,sfactor,aoffset,
     &                   ptype,acswath,alswath,descrip,geoptr)

      usrlog="Setting attributes for "//fieldname(1:sl)//
     &       " - set_array_attr"
      call ckstatus_f(rtn,usrlog,funcname,LOGFLAG)


*/  The following code slice for setting the fill value is commented
*/  out because currently the FORTRAN HDF_EOS routine for assigning
*/  fill values, swsetfill, is not working properly in the way that
*/  a fill value can be assigned to a SDS but this fill value is
*/  not accessible to a data reader.


*/  Set the fillvalue of the data field.

      srtn=swsetfill(swathid,fieldname(1:sl),fillval)

      usrlog="Setting fillvalue for "//fieldname(1:sl)//
     &       " - swsetfill"
      call ckstatus_f(srtn,usrlog,funcname,LOGFLAG)


*/  Set the valid range of the data field.

      attr="valid_range"
      type=numbertype
      nms=2
      hrtn=sfsattr(sdsid,attr,type,nms,vrange)

      usrlog="Setting valid range for "//fieldname(1:sl)//
     &       " - sfsattr"
      call ckstatus_f(hrtn,usrlog,funcname,LOGFLAG)


*/  Successfully returns.

      create_array_sht=SUCCEED

      return
      end


*/-------------------------------------------------------------------/*

      integer function create_array_byt(sdid,swathid,fieldname,
     &                                  dimlist,numbertype,merge,
     &                                  longname,datau,vrange,fillval,
     &                                  sfactor,aoffset,ptype,acswath,
     &                                  alswath,descrip,geoptr)

C----------------------------------------------------------------------
C !F77
C
C !Description:
C
C   This routine creates a array with type byte using HDF-EOS Swath
C   API and calls set_array_attr() and HDF native routines to set the
C   attributes of the SDS.
C
C !Input Parameters:
C
C   Type               Name         Description
C   ----               ----         -----------
C   integer            sdid         SDS interface ID
C   integer            swathid      Swath ID
C   character*(*)      fieldname    SDS array name
C   character*(*)      dimlist      Dimension list
C   integer            numbertype   Data type code
C   integer            merge        Merge code
C   character*(*)      longname     Long name of SDS
C   character*(*)      datau         Data datau
C   integer*1          vrange(*)    Valid range
C   integer*1          fillval      Fill value
C   double precision   sfactor      Scaling factor
C   double precision   aoffset      Added offset
C   character*(*)      ptype        Parameter type
C   integer            acswath(*)   Across swath sampling
C   integer            alswath(*)   Along swath sampling
C   character*(*)      descrip      Description
C
C
C !Output Parameters: None
C
C !Revision History:
C $Log: Utility_V2.f,v $
c 01/29/98 fhliang
c added NCSA acknowledgement.
C
C !Team-Unique Header:
C
C   This software was developed by the MODIS Science Data Support Team
C   (SDST) for the National Aeronautics and Space Administration,
C   Goddard Space Flight Center, under contract NAS5-32373.
C
c HDF portions developed at the National Center for Supercomputing
c Applications at the University of Illinois at Urbana-Champaign.
c
C !References and Credits:
C
C   Developer: JC Guu 03/10/97 jguu@ltpmail.gsfc.nasa.gov
C
C !Externals:
C
C   Type            Name           Source               Description
C   ----            ----           ------               -----------
C   integer         SUCCEED        hdf.inc              Return code for
C                                                       succeed.
C
C !Internals:
C
C   Type              Name             Description
C   ----              ----             -----------
C   character*256     attr             Character string for attribute names
C   character*512     funcname         Character string for fnction names
C   character*512     usrlog           Character string for log messages
C   integer           srtn             Hold return code from Swath API
C   integer           hrtn             Hold return code from HDF.
C   integer           rtn              Hold generic return code.
C   integer           nms              Number of values for an attribute
C   integer           type             Type code
C   integer           sl               True string length
C   integer           sdsidx           SDS index
C   integer           sdsid            SDS ID
C
C
C !END
C----------------------------------------------------------------------
      implicit none

      include 'hdf.inc'
      include 'dffunc.inc'
      include 'ckstatus.inc'


*/  Passed-in arguments

      character*(*) fieldname,dimlist,longname,datau,ptype,descrip,
     &              geoptr
      integer swathid,numbertype,merge
      integer sdid,acswath(*),alswath(*)
      byte vrange(*),fillval
      double precision sfactor,aoffset


*/  Functions

      integer swdefdfld,swdefgfld,swsetfill,set_array_attr,strlen


*/  Variables

      character*256 attr
      character*512 funcname,usrlog
      integer srtn,hrtn,rtn,nms,type,sl
      integer sdsidx,sdsid

      data sdsidx,sdsid/-1,-1/

*---------------------------------------------------------------------*


*/  Define the data field.

      srtn=swdefdfld(swathid,fieldname,dimlist,numbertype,merge)

      funcname="create_array_byt"
      sl=strlen(fieldname)
      usrlog="Creating array "//fieldname(1:sl)//
     &       " - swdefdfld"
      call ckstatus_f(srtn,usrlog,funcname,LOGFLAG)


*/  Retrieve the SDS index for the data field.

      sdsidx=sfn2index(sdid,fieldname)

      usrlog="Retrieving SDS index for "//fieldname(1:sl)//
     &       " - sfn2index"
      call ckstatus_f(sdsidx,usrlog,funcname,LOGFLAG)


*/  Retrieve the SDS ID for the data field.

      sdsid=sfselect(sdid,sdsidx)

      usrlog="Retrieving SDS ID for "//fieldname(1:sl)//
     &       " - sfselect"
      call ckstatus_f(sdsid,usrlog,funcname,LOGFLAG)


*/  Set the attributes with same types throughout the swath.

      rtn=set_array_attr(sdsid,fieldname,longname,datau,sfactor,aoffset,
     &                   ptype,acswath,alswath,descrip,geoptr)

      usrlog="Setting attributes for "//fieldname(1:sl)//
     &       " - set_array_attr"
      call ckstatus_f(rtn,usrlog,funcname,LOGFLAG)


*/  The following code slice for setting the fill value is commented
*/  out because currently the FORTRAN HDF_EOS routine for assigning
*/  fill values, swsetfill, is not working properly in the way that
*/  a fill value can be assigned to a SDS but this fill value is
*/  not accessible to a data reader.


*/  Set the fillvalue of the data field.

      srtn=swsetfill(swathid,fieldname(1:sl),fillval)

      usrlog="Setting fillvalue for "//fieldname(1:sl)//
     &       " - swsetfill"
      call ckstatus_f(srtn,usrlog,funcname,LOGFLAG)


*/ Set the valid range of the data field.

      attr="valid_range"
      type=numbertype
      nms=2
      hrtn=sfsattr(sdsid,attr,type,nms,vrange)

      usrlog="Setting valid range for "//fieldname(1:sl)//
     &       " - sfsattr"
      call ckstatus_f(hrtn,usrlog,funcname,LOGFLAG)


*/  Successfully returns.

      create_array_byt=SUCCEED

      return
      end

C****************************************************************************
      subroutine ckstatus_f_mod05(rtn,usrlog,funcname)

C----------------------------------------------------------------------------
C !F77
C
C !Description:
C
C   This routine is used to check the return code of functions. If
C   a function only returns SUCCEED(0) or DAIL(-1), then this routine
C   should be used instead of ckstatus_s().
C
C !Input Parameters:
C
C   integer         rtn        The return code from a function
C   character*(*)   usrlog     The message to be written to the LogStatus
C   character*(*)   funcname   The name of the function which returns.
C
C !Output Parameters: None
C
C !Revision History:
C $Log: Utility_V2.f,v $
c 11/01/1999 rhucek
c resolved forchk problem 'referenced character elements assigned'.
c
c 01/29/98 fhliang
c added NCSA acknowledgement.
C
C !Team-Unique Header:
C
C   This software was developed by the MODIS Science Data Support Team
C   (SDST) for the National Aeronautics and Space Administration,
C   Goddard Space Flight Center, under contract NAS5-32373.
C
c HDF portions developed at the National Center for Supercomputing
c Applications at the University of Illinois at Urbana-Champaign.
c
C !References and Credits:
C
C   Developer: JC Guu 03/10/97 jguu@ltpmail.gsfc.nasa.gov
C
C !Externals:
C   FAIL                       (hdf.inc)
C   MODIS_E_GENERIC            (PGS_MODIS_39500.f)
C   MODIS_S_GENERIC            (PGS_MODIS_39500.f)
C
C !Internals:
C
C !END
C----------------------------------------------------------------------------

      implicit none

      include 'PGS_MODIS_39500.f'
      include 'hdf.inc'

      integer    MAXMSGLEN,         MAXUSRLOGLEN
      parameter (MAXMSGLEN = 10000, MAXUSRLOGLEN = 9000)


      integer rtn
      character*(*)         funcname,usrlog
      character*(MAXMSGLEN) message

*---------------------------------------------------------------------------*

c.....message too long to handle
      if ( LEN(usrlog) .gt. MAXUSRLOGLEN ) then
         message = 'Input message from routine too long; ckstatus_f_mod05 unable to print. '
     2   // char(10) // 'Operator Action:  Notify SDST.'

         call modis_smf_setdynamicmsg(MODIS_E_GENERIC,message,funcname)

c.....fail scenario
      else if (rtn.eq.FAIL) then
         message = usrlog
     1   // char(10) // 'Operator Action:  Check system resources/environment, '
     2   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     3   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'

         call modis_smf_setdynamicmsg(MODIS_E_GENERIC,message,funcname)

c.....success
      else
         call modis_smf_setdynamicmsg(MODIS_S_GENERIC,usrlog,funcname)
      end if

      end



C***************************************************************************

      integer function set_array_attr(sdsid,fieldname,longname,datau,
     &                                sfactor,aoffset,ptype,acswath,
     &                                alswath,descrip,geoptr)

C!F77-------------------------------------------------------------------
C
C!Description: This routine sets the attributes common to all SDS.
C
C!Input Parameters:
C   acswath      Across swath sampling
C   alswath      Along swath sampling
C   aoffset      Added offset
C   datau        Data unit
C   descrip      Description
C   fieldname    SDS array name
C   geoptr       Geolocation pointer
C   longname     Long name of SDS
C   ptype        Parameter type
C   sdsid        SDS ID
C   sfactor      Scaling factor
C
C
C!Output Parameters: None
C
C!Revision History:
C $Log: Utility_V2.f,v $
c 01/29/98 fhliang
c added prolog.
C
C!Team-Unique Header:
C
C   This software was developed by the MODIS Science Data Support Team
C   (SDST) for the National Aeronautics and Space Administration,
C   Goddard Space Flight Center, under contract NAS5-32373.
C
C HDF functions were developed at the National Center for Supercomputing
C Applications at the Univ. of Illinois at Urbana-Champaign.
C
C!References and Credits:
C
C   Developed by JC Guu 03/10/97 jguu@ltpmail.gsfc.nasa.gov
C   Modified by Fay Liang fhliang@ltpmail.gsfc.nasa.gov
C           and Rich Hucek rhucek@ltpmail.gsfc.nasa.gov
C
C!Externals (non-SDST code):
C   Functions:
C    sfscattr
C    sfsnattr
C
C  Named Constants:
C
C!Internals: (SDST code)
C   Variables:
C   attr             Attribute name
C   funcname         Character string for function names
C   hrtn             Return value from an HDF function
C   nms              Number of values for an attribute
C   sla              Length of attr
C   slf              Length of fieldname
C   type             Data type code
C   usrlog           Character string for log messages
C
C   Named Constants: None
C
C   Functions and Subroutines:
C    strlen          (atmos shared code)
C
C!END-------------------------------------------------------------------

      implicit none

      include 'hdf.inc'
      include 'dffunc.inc'
      include 'ckstatus.inc'


*/  Passed-in arguments

      character*(*) fieldname,longname,datau,ptype,descrip,geoptr
      integer sdsid,acswath(*),alswath(*)
      double precision sfactor,aoffset

*/  Functions

      integer strlen,swwratt

*/  Variables

      character*256 attr
      character*512 funcname,usrlog
      integer hrtn,nms,type,slf,sla

*---------------------------------------------------------------------------*


*/  Attributes "long_name","units","scale_factor","add_offset",
*/  "Parameter_Type","Cell_Along_Swath_Sampling",
*/  "Cell_Across_Swath_Sampling","description" a and
*/  "Geolocation_Pointer" are set.
*/  For some data fields, "Cell_Across_Swath_Sampling" and/or
*/  "description" are not applicable, therefore a IF construct
*/  is provided to see if these attributes are needed by the
*/  data field in question.



      attr="long_name"
      type=DFNT_CHAR8
      nms=strlen(longname)
      hrtn=sfsattr(sdsid,attr,type,nms,longname)

      slf=strlen(fieldname)
      funcname="set_array_attr"
      sla=strlen(attr)
      usrlog="Setting attribure "//attr(1:sla)//" for "//
     &       fieldname(1:slf)//" - sfsattr"
      call ckstatus_f(hrtn,usrlog,funcname,LOGFLAG)



      attr="units"
      type=DFNT_CHAR8
      nms=strlen(datau)
      hrtn=sfsattr(sdsid,attr,type,nms,datau)

      sla=strlen(attr)
      usrlog="Setting attribure "//attr(1:sla)//" for "//
     &       fieldname(1:slf)//" - sfsattr"
      call ckstatus_f(hrtn,usrlog,funcname,LOGFLAG)



      attr="scale_factor"
      type=DFNT_FLOAT64
      nms=1
      hrtn=sfsattr(sdsid,attr,type,nms,sfactor)

      sla=strlen(attr)
      usrlog="Setting attribure "//attr(1:sla)//" for "//
     &       fieldname(1:slf)//" - sfsattr"
      call ckstatus_f(hrtn,usrlog,funcname,LOGFLAG)



      attr="add_offset"
      type=DFNT_FLOAT64
      nms=1
      hrtn=sfsattr(sdsid,attr,type,nms,aoffset)

      sla=strlen(attr)
      usrlog="Setting attribure "//attr(1:sla)//" for "//
     &       fieldname(1:slf)//" - sfsattr"
      call ckstatus_f(hrtn,usrlog,funcname,LOGFLAG)



      attr="Parameter_Type"
      type=DFNT_CHAR8
      nms=strlen(ptype)
      hrtn=sfsattr(sdsid,attr,type,nms,ptype)

      sla=strlen(attr)
      usrlog="Setting attribure "//attr(1:sla)//" for "//
     &       fieldname(1:slf)//" - sfsattr"
      call ckstatus_f(hrtn,usrlog,funcname,LOGFLAG)



      if(acswath(1).gt.0) then

          attr="Cell_Across_Swath_Sampling"
          type=DFNT_INT32
          nms=3
          hrtn=sfsattr(sdsid,attr,type,nms,acswath)

          sla=strlen(attr)
          usrlog="Setting attribure "//attr(1:sla)//" for "//
     &           fieldname(1:slf)//" - sfsattr"
          call ckstatus_f(hrtn,usrlog,funcname,LOGFLAG)

      endif



      attr="Cell_Along_Swath_Sampling"
      type=DFNT_INT32
      nms=3
      hrtn=sfsattr(sdsid,attr,type,nms,alswath)

      sla=strlen(attr)
      usrlog="Setting attribure "//attr(1:sla)//" for "//
     &       fieldname(1:slf)//" - sfsattr"
      call ckstatus_f(hrtn,usrlog,funcname,LOGFLAG)



      if(strlen(descrip).gt.0) then

          attr="description"
          type=DFNT_CHAR
          nms=strlen(descrip)
          hrtn=sfsattr(sdsid,attr,type,nms,descrip)

          sla=strlen(attr)
          usrlog="Setting attribure "//attr(1:sla)//" for "//
     &           fieldname(1:slf)//" - sfsattr"
          call ckstatus_f(hrtn,usrlog,funcname,LOGFLAG)

      endif



      attr="Geolocation_Pointer"
      type=DFNT_CHAR8
      nms=strlen(geoptr)
      hrtn=sfsattr(sdsid,attr,type,nms,geoptr)

      sla=strlen(attr)
      usrlog="Setting attribure "//attr(1:sla)//" for "//
     &       fieldname(1:slf)//" - sfsattr"
      call ckstatus_f(hrtn,usrlog,funcname,LOGFLAG)



      set_array_attr=SUCCEED

      return
      end
******/---------------------------------------------------------------/******

      subroutine ckstatus_f(rtn,usrlog,funcname,runmode)

      implicit none

      include 'ckstatus.inc'
      include 'PGS_MODIS_39500.f'
      include 'hdf.inc'

C----------------------------------------------------------------------------
C !F77
C
C !Description:
C
C   This routine is used to check the return code of functions. If
C   a function only returns SUCCEED(0) or FAIL(-1), then this routine
C   should be used instead of ckstatus_s().
C
C !Input Parameters:
C
C   integer         rtn        The return code form a function
C   character*(*)   usrlog     The message to be written to the LogStatus
C   character*(*)   funcname   The name of the function which returns.
C   integer         runmode    The flag for turning on or off the writing
C                              of SUCCESS messages to LogStatus
C
C !Output Parameters: None
C
C !Revision History:
C $Log: Utility_V2.f,v $
c 11/01/1999 rhucek
c resolved forchk problem 'referenced character elements assigned'.
c
c 01/29/98 fhliang
c added NCSA acknowledgement.
c
c Revision 1.3  1997/08/26  15:57:45  rhucek
c removed functions int2str and colseq
c
c Revision 1.2  1997/06/18  13:45:07  vlin
c used to run V2 Metadata writer
c
C
C !Team-Unique Header:
C
C   This software was developed by the MODIS Science Data Support Team
C   (SDST) for the National Aeronautics and Space Administration,
C   Goddard Space Flight Center, under contract NAS5-32373.
C
c HDF portions developed at the National Center for Supercomputing
c Applications at the University of Illinois at Urbana-Champaign.
c
C !References and Credits:
C
C   Developer: JC Guu 03/10/97 jguu@ltpmail.gsfc.nasa.gov
C
C !Externals:
C
C   FAIL             (hdf.inc)
C   MODIS_E_GENERIC  (PGS_MODIS_39500.f)
C   MODIS_S_GENERIC  (PGS_MODIS_39500.f)
C   DEBUG            (ckstatus.inc)
C
C !Internals:
C
C !END
C----------------------------------------------------------------------------

c.....PARAMETER definitions
      integer    MAXMSGLEN,         MAXUSRLOGLEN
      parameter (MAXMSGLEN = 10000, MAXUSRLOGLEN = 9000)

c.....function arguments
      character*(*) funcname,usrlog
      integer rtn,runmode

c.....local variables
      character*(MAXMSGLEN) message


c.....message too long to handle
      if ( LEN(usrlog) .gt. MAXUSRLOGLEN ) then
         message = 'Input message from routine too long; ckstatus_f unable to print. '
     2   // char(10) // 'Operator Action:  Notify SDST.'

         call modis_smf_setdynamicmsg(MODIS_E_GENERIC,message,funcname)

c.....fail scenario
      elseif (rtn.eq.FAIL) then
         message = usrlog
     1   // char(10) // 'Operator Action:  Check system resources/environment, '
     2   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     3   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'

         call modis_smf_setdynamicmsg(MODIS_E_GENERIC,message,funcname)

      elseif(runmode.eq.DEBUG) then
         call modis_smf_setdynamicmsg(MODIS_S_GENERIC,usrlog,funcname)
      end if


      end
******/---------------------------------------------------------------/******

      subroutine ckstatus_s(rtn,usrlog,funcname,runmode)

      implicit none

      include 'PGS_SMF.f'
      include 'PGS_MODIS_39500.f'
      include 'ckstatus.inc'

C----------------------------------------------------------------------------
C !F77
C
C !Description:
C
C   This routine is used to check the return code of functions. If
C   a function returns PGS_S_SUCCESS for success, this routine should
C   be used instead of ckstatus_f().
C
C !Input Parameters:
C
C   integer         rtn        The return code form a function
C   character*(*)   usrlog     The message to be written to the LogStatus
C   character*(*)   funcname   The name of the function which returns.
C   integer         runmode    The flag for turning on or off the writing
C                              of SUCCESS messages to LogStatus
C
C !Output Parameters: None
C
C !Revision History:
C $Log: Utility_V2.f,v $
c 11/01/1999 rhucek
c resolved forchk problem 'referenced character elements assigned'.
c
c 01/29/98 fhliang
c added NCSA acknowledgement.
c
c Revision 1.3  1997/08/26  15:57:45  rhucek
c removed functions int2str and colseq
c
c Revision 1.2  1997/06/18  13:45:07  vlin
c used to run V2 Metadata writer
c
C
C !Team-Unique Header:
C
C   This software was developed by the MODIS Science Data Support Team
C   (SDST) for the National Aeronautics and Space Administration,
C   Goddard Space Flight Center, under contract NAS5-32373.
C
c HDF portions developed at the National Center for Supercomputing
c Applications at the University of Illinois at Urbana-Champaign.
c
C !References and Credits:
C
C   Developer: JC Guu 03/10/97 jguu@ltpmail.gsfc.nasa.gov
C
C !Externals:
C
C   FAIL                       (hdf.inc)
C   MODIS_E_GENERIC            (PGS_MODIS_39500.f)
C   MODIS_S_GENERIC            (PGS_MODIS_39500.f)
C   DEBUG                      (ckstatus.inc)
C
C !Internals:
C
C !END
C----------------------------------------------------------------------------

c.....PARAMETER definitions
      integer    MAXMSGLEN,         MAXUSRLOGLEN
      parameter (MAXMSGLEN = 10000, MAXUSRLOGLEN = 9000)

c.....function arguments
      character*(*) funcname,usrlog
      integer       rtn,runmode

c.....local variables
      character*(MAXMSGLEN) message


c.....input usrlog too long
      if ( LEN(usrlog) .gt. MAXUSRLOGLEN ) then
         message = 'Input message from routine too long; ckstatus_f unable to print. '
     2   // char(10) // 'Operator Action:  Notify SDST.'

         call modis_smf_setdynamicmsg(MODIS_E_GENERIC,message,funcname)

c.....success scenario
      elseif(rtn.eq.PGS_S_SUCCESS) then

          if(runmode.eq.DEBUG) then
              call modis_smf_setdynamicmsg(MODIS_S_GENERIC,usrlog,funcname)
          end if

c.....fail
      else
        message = usrlog
     1  // char(10) // 'Operator Action:  Check system resources/environment, '
     2  // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     3  // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'

        call modis_smf_setdynamicmsg(MODIS_E_GENERIC,message,funcname)
      end if


      end
