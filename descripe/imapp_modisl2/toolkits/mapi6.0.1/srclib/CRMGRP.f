      FUNCTION CRMGRP(mfile, grpnm, clsnm)
C
C--------------------------------------------------------------------------
C!F77
C
C!Description: Function CRMGRP is part of a large software system called the 
C              MODIS Application Programming Interface(API) utility, 
C              abbreviated M-API.  The M-API Utility consists of subroutines 
C              which allow MODIS Science Team-supplied software to read in 
C              Level 1B radiance bands and write out output products and 
C              metadata to HDF files.  The functionality of the M-API is 
C              defined in the MODIS API User's Guide, Version 1, dated 4/3/95.
C
C              CRMGRP creates an HDF Vgroup structure in a MODIS HDF file to
C              store table and array structures.  It must be called before any
C              of the data objects to be aggregated in it are created.  The use
C              of data groups is optional, but data objects stored in them are
C              documented in the MODIS Product File Definitions in Appendix C.
C              A data gruop with the name grpnm must be unique in a file.  This
C              prevents confusion that is cuased by multiple data groups with 
C              the same name.
C
C!Input Parameters:
C mfile        Array that is used to reference the MODIS-HDF file.
C grpnm        ASCII string that will be the name of the data group.
C clsnm        (Optional) additional ASCII string that will be the class name
C              of the group.  Set to a blank string " " for default.
C
C!Output Parameters: None.
C 
C Returns:     MAPIOK if successful, MFAIL if an error occurs.
C
C External references:
C              nccrmgrp           (mapic.h)
C
C!Revision History:
C
C
C !Team-unique Header:
C             This software is developed by the MODIS Science Data Support
C             Team for the National Aeronautics and Space Administration,
C             Goddard Space Flight Center, under contract NAS5-32373.
C
C !References and Credits:
C     WRITTEN BY:    Xia W. Chang,  xia@ltpmail.gsfc.nasa.gov
C                    Version 2.0 design From Prototype by Liping Di
C                    MODIS Science Data Support Team
C                    Maryland Corporate Center
C                    7501 Forbes Blvd - Suite 103
C                    Seabrook, MD  20706
C
C----------------------------------------------------------------------
C!END
C
      IMPLICIT NONE
      INCLUDE  'mapi.inc'
      INTEGER        mfile(MODFILLEN)
      CHARACTER*(*)  grpnm
      CHARACTER*(*)  clsnm

C     Declare local variable
C
      INTEGER        ret
C
      call ccrmgrp(mfile,grpnm,len(grpnm),clsnm,len(clsnm),ret)
      
      CRMGRP = ret

      RETURN
      END
