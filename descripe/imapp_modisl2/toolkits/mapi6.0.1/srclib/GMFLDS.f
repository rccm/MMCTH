      FUNCTION GMFLDS(mfile, tbname, group, strln, recno,
     *                fldno, fldnm, dtype, clss)
C
C-------------------------------------------------------------------------
C!F77
C
C!Description: Function GMFLDS is part of a large software system called the 
C              MODIS Applications Programming Interface (API) utility, 
C              abbreviated M-API.  The M-API Utility consists of subroutines 
C              which allow MODIS Science Team-supplied software to read and 
C              write data and metadata from/to HDF files.  The functionality 
C              of the M-API is defined in the MODIS API specification.
C
C              GMFLDS retrieves the essential characteristics of an HDF Vdata
C              table structure contained in a MODIS-HDF file. This provides 
C              the information needed for properly reading data from the table
C              structure using GMTBL or to write it using PMTBL. An error 
C              (MFAIL) will be returned if 1) The output strings are not
C              long enough to contain the data type or field name strings for
C              all the Vdata's fields, 2) an unknown (e.g. not supported by 
C              the MODIS API) number type is encountered or 3) an HDF routine
C              FAILs.  The data type string will be returned truncated to the
C              point where the fault occurred.  The strln will return with the
C              actual string length required to hold the larger of the two
C              output strings.  If the latter two errors occur, however, 
C              strln is set to 0.
C
C              The group string provides the facility to select a table 
C              structure existing in a particular HDF 'Vgroup' data group.
C              Alternatively, the entire file will be searched for a table
C              structure named tbname if group=' ' in FORTRAN.
C
C!Input Parameters:
C mfile        M-API FORTRAN file array that is used to reference the 
C              MODIS-HDF file containing the target table structure.
C tbname       ASCII string name of the target table structure.
C group        ASCII string name of the data group containing the target 
C              table structure.  If set to ' ', the entire file will be 
C              searched for the table structure named tbname.
C
C!Output Parameters:
C strln        Returns the minimum string length actually required to hold 
C              the longer of the two strings. It is set to 0 if a functional 
C              error occurs.
C recno        Number of records (row) present in the table structure.
C fldno        Number of fields (columns) present in the table structure.
C fldmn        ASCII string of comma-delimited ASCII string table headers.
C dtype        ASCII string of comma-delimited data types for each table 
C              fields.  Permitted FORTRAN data types:
C                  'CHARACTER*(*)'
C                  'INTEGER*1'
C                  'UINTEGER*1'
C                  'INTEGER*2'
C                  'UINTEGER*2'
C                  'INTEGER*4'
C                  'UINTEGER*4'
C                  'INTEGER*8'
C                  'UINTEGER*8'
C                  'REAL*4'
C                  'REAL*8'
C clss        ASCII string for the class name of the table structure. 
C             Provided array should be at least VSNAMELENMAX(64) bytes long.
C
C Returns:    MAPIOK if successful, MFAIL if unable to retrieve all the
C             requested outputs or unable to attach to the Vdata at all.
C
C External references:
C              ncgmflds           (mapic.h)
C
C!Revision History:
C $Log: GMFLDS.f,v $
C Revision 1.1  1999/04/15 19:26:51  jayshree
C Initial revision
C
C
C !Team-unique Header:
C             This software is developed by the MODIS Science Data Support
C             Team for the National Aeronautics and Space Administration,
C             Goddard Space Flight Center, under contract NAS5-32373.
C
C !References and Credits:
C     WRITTEN BY:    Xia W. Chang,  xia@ltpmail.gsfc.nasa.gov
C                    Version 1.4 original development From Prototype by 
C                    Liping Di, Dec 12, 1995
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
      INTEGER        strln, recno, fldno
      CHARACTER*(*)  tbname, group, fldnm, dtype, clss

C     Declare local variable
      INTEGER  ret

      call cgmflds(mfile, tbname, len(tbname), group, len(group),
     *              strln, recno, fldno, fldnm, len(fldnm), 
     *              dtype, len(dtype), clss, len(clss), ret)
 
      GMFLDS = ret
   
      RETURN
      END
