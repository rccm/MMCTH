       FUNCTION EMOBJ(modfil, nname, group, type)
C
C------------------------------------------------------------------------------
C!F77
C
C!Purpose: ends the access to an individual or a group of opened objects
C            (close the object and release resources).
C
C!Description: Function EMOBJ is part of a larger software system called
C                the MODIS Applications Programming Interface (API) Utility,
C                abbreviated M-API.  The M-API Utility consists of subroutines
C                which allow MODIS Science Team-supplied software to read and
C                write data and metadata from/to HDF files. The functionality
C                of M-API is defined in the MODIS Application Program Interface
C                (API) Specification.
C
C                A Vdata or SDS is opened and remains opened when application
C                programs call any of M-API routines to access the object. The
C                M-API access routines automatically call M-API internal routine
C                addid to insert the HDF object ID and related information into
C                the ring super structure in MODFILE. In the subsequent calls to
C                access the same object, M-API routines will use the id store in
C                the ring super structure for fast data access.   
C
C                EMOBJ ends the access to an individual or a group of opened
C                objects by deleting objects' DATAID structure from the ring super
C                structure, releasing memory, and detaching (for Vdata) or
C                end-accessing (for SDS) the objects. All opened objects will be
C                closed automatically once an application program calls either CLMFIL
C                or CPMFIL. Therefore,  As long as an application program calls either
C                CLMFIL or CPMFIL, the program does not need to worry  when or how
C                to use this routine to close an object or a group of object. However,
C                if an application program determines the an object will no longer
C                be accessed and wish to end  the access to the object or releasing
C                computer resources, the application program can call this routines
C                any time before the application program closes the opened HDF file.
C
C !Input parameters:
C
C         modfil  IN: Array that is used to refer a MODIS HDF file containing
C                     objects to be closed (Vdata or SDS).
C
C         nname    IN: The name of the object. If the name is set to ' ', all
C                     objects matched with group and object type will be closed.
C
C         group   IN: The name of the group to which objects belongs. If group is
C                     set to ' ', all lone objects (objects belonging to no group)
C                     matched with  name and type will be closed. If both name and
C                     group are ' ', all objects matched with type will be closed.
C         type    IN: The object type: MODIS_ARRAY(for SDS, the numerical value is
C                     720, the same value as DFTAG_NDG),  MODIS_TABLE(for Vdata, the
C                     numerical value is 1962, the same value as DFTAG_VH), or
C                     ALL_MODIS_TYPES (the nuemrical value is  0). If type is 720,
C                     only SDS objects will be closed.  If type is 1962, only Vdata
C                     will be closed. If type is  0, both Vdata and SDS objects
C                     specified by name and group will be closed. Therefore, to close
C                     all opened objects, set both name and group to NULL and set type
C                     to 0.
C
C!Output Parameters:  None
C
C  Returns:  Number of objects closed, or MFAIL if errors.
C
C  External References: 
C
C !Revision History:
C $Log: EMOBJ.f,v $
C Revision 1.1  1999/04/15 19:26:51  jayshree
C Initial revision
C
C
C !Team-unique header:
C  Portions developed at the National Center for Supercomputing
C  Applications at the Univ. of Illinois at Urbana-Champaign.
C
C !References and Credits:
C  This software is developed by the MODIS Science Data Support
C  Team for the National Aeronautics and Space Administration,
C  Goddard Space Flight Center, under contract NAS5-32373.
C
C !!Design Notes:
C
C!END-------------------------------------------------------------------
C

      IMPLICIT NONE
      INCLUDE  'mapic.inc'

      INTEGER   modfil(MODFILLEN)
      CHARACTER*(*) nname, group 
      INTEGER   type, ret

C     Call cemobj  to close MODIS object. 

      CALL cemobj( modfil, nname, LEN(nname),
     &             group, LEN(group), type, ret) 

      EMOBJ = ret

      RETURN
      END

