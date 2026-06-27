      FUNCTION CRMTBL(mfile,tbname,class,grpnm,field,dtype)

C--------------------------------------------------------------------------
C!F77
C
C!Purpose:     Initializes a table structure in a MODIS-HDF file to store a multi-
C              column table of data.
C
C!Description: Function GMTBL is part of a large software system called the 
C              MODIS Application Programming Interface(API) utility, 
C              abbreviated M-API.  The M-API Utility consists of subroutines 
C              which allow MODIS Science Team-supplied software to read in 
C              Level 1B radiance bands and write out output products and 
C              metadata to HDF files.  The functionality of the M-API is 
C              defined in the MODIS API User's Guide, Version 1, dated 4/3/95.
C

C              CRMTBL (CRMTBL) creates an HDF Vdata structure in a MODIS HDF 
C              file to store a new data table.  It must be called before the data 
C              may be written to the file using putMODIStable (PMTBL).  The text 
C              headers for each field (column) and the data type stored in each 
C              field must be provided.
C              
C              The grpnm  string provides the facility to place the new table in 
C              an HDF  ‘Vgroup’ data group.  If a Vgroup with the name grpnm 
C              does not exist, the table structure will not be created.  The table 
C              may be placed in the file outside of any Vgroup by setting grpnm 
C              = "" or " " in FORTRAN.
C              
C              An empty Vdata may not be created, so a dummy record is inserted 
C              into it.  This dummy record should be overwritten with the first 
C              call from PMTBL.  
C              
C              CRMTBL calls C createMODIStable using the wrapping method. 
C              createMODIStable makes extensive use of dynamically allocated 
C              memory.  All of this memory is freed, however, prior to exiting 
C              the routine.
C              
C              If a table with the name tbname  is created outside of a Vgroup, it 
C              must not already exist in the file.  This is to prevent the confusion 
C              caused by multiple data objects with the same name.  Tables with 
C              the same name may be stored in the same file, however, if they 
C              are placed in separate Vgroups.
C              
C FORTRAN:     INTEGER FUNCTION CRMTBL(mfile, tbname, class, grpnm,  field, 
C              dtype)
C              INTEGER mfile(3)
C              CHARACTER*(*) tbname, class, group, field, dtype
C              
C!Input Parameters:
C              mfile	IN: 	Array that is used to reference the MODIS-HDF 
C                               file receiving the new table.
C              tbname	IN:	ASCII string that will be the name of the 
C                               table.
C              class	IN:	ASCII string that will be the class name of the 
C                               table.  If set to an empty string, the table will 
C                               have no class.
C              group	IN:	ASCII string name of the data group to place 
C                               the new table in.  If set to a empty string, the 
C                               table will not be placed in a data group.
C              field	IN:	comma-delimited ASCII string table headers.  
C                               The headers should be in the same order that 
C                               the data for each table row will subsequently 
C                               be written in.
C              dtype	IN:	String of comma-delimited data types for each 
C                               table field.  The data type strings should be in 
C                               the same order that the data for each table 
C                               row will subsequently be written in.
C              
C                       Permitted FORTRAN data types:
C                               'CHARACTER*(*)'
C				'INTEGER*1'
C				'UINTEGER*1'
C				'INTEGER*2'
C				'UINTEGER*2'
C				'INTEGER*4'
C				'UINTEGER*4'
C				'REAL*4'
C				'REAL*8'
C
C!Output Parameters:	NONE.
C
C Returns:	MAPIOK if successful, MFAIL if an error occurs.
C
C External references:
C                      ccrmtbl           (mapic.h)
C		       MODFILLEN	 (mapi.inc)
C
C!Revision History:
C $Log: CRMTBL.f,v $
C Revision 1.1  1999/04/15 19:26:51  jayshree
C Initial revision
C
C
C !Team-unique Header:
C             This software is developed by the MODIS Science Data Support
C             Team for the National Aeronautics and Space Administration,
C             Goddard Space Flight Center, under contract NAS5-32373.
C
C             Portions developed at the National Center for Supercomputing
C             Applications at the Univ. of Illinois at Urbana-Champaign.
C
C !References and Credits:
C
C----------------------------------------------------------------------
C!END
C
      IMPLICIT NONE
      INCLUDE  'mapic.inc'
      INTEGER  mfile(MODFILLEN) 
      CHARACTER*(*) tbname, class, grpnm, field, dtype
      
      INTEGER  ret

      CALL ccrmtbl(mfile,tbname,len(tbname),class,len(class),grpnm,
     1             len(grpnm),field,len(field),dtype,len(dtype),ret)
      CRMTBL = ret
      RETURN
      END
