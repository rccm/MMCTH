       FUNCTION DF2DT(dfnt, dtype) 
C
C-----------------------------------------------------------------------
C!F77
C 
C !Description:Subroutine DF2DT is part of a larger software system 
C               called the MODIS Applications Programming Interface 
C               (API) Utility, abbreviated M-API.  The M-API Utility con-
C               sists of subroutines which allow MODIS Science Team-sup-
C               plied software to read in Level 1B radiance bands and 
C               write out output products and metadata to HDF files. 
C 
C               DF2DT takes an HDF number type (DFNT_*) used by HDF
C               routines to identify number types and returns one of the
C               following M-API data type strings, if possible. If no
C               match can be made, a null string is returned. An eight
C               or more element character string must be provided for
C               datatype.
C
C !Input Parameters:INTEGER dfnt       :  HDF number type constant
C
C !Output Parameters:CHARACTER*(*) dtype:  M-API data type string.
C
C !Revision History:
C
C $Id: DF2DT.f,v 1.1 1999/04/15 19:26:51 jayshree Exp $
C
C $Log: DF2DT.f,v $
C Revision 1.1  1999/04/15 19:26:51  jayshree
C Initial revision
C
c
C
C !Team-unique Header:
C
C               This software is developed by the MODIS Science Data Support
C               Team for the National Aeronautics and Space Administration,
C               Goddard Space Flight Center, under contract NAS5-32373.
C
C               Portions developed at the National Center for Supercomputing
C               Applications at the Univ. of Illinois at Urbana-Champaign.
C
C !References and Credits:
C
C               WRITTEN BY:
C               Qi Huang      06/14/95
C               Research and Data Systems Corporation
C               SAIC/GSC MODIS Science DATA Support Office
C               7501 Forbes Blvd
C               Seabrook, MD 20706
C
C               Mitchell Weiss      07/26/95
C               Modification / Testing
C
C
C !Design Notes:
C
C       Externals:
C
C          Functions:
C               DF2DT             (mapic.inc)
C
C          Variables:
C               DFNT_INT8         (hdf.inc included in "mapic.inc")
C               DFNT_INT16        (hdf.inc included in "mapic.inc")
C               DFNT_INT32        (hdf.inc included in "mapic.inc")
C               DFNT_INT64        (hdf.inc included in "mapic.inc")
C               DFNT_UINT8        (hdf.inc included in "mapic.inc")
C               DFNT_UINT16       (hdf.inc included in "mapic.inc")
C               DFNT_UINT32       (hdf.inc included in "mapic.inc")
C               DFNT_FLOAT32      (hdf.inc included in "mapic.inc")
C               DFNT_FLOAT64      (hdf.inc included in "mapic.inc")
C               DFNT_CHAR         (hdf.inc included in "mapic.inc")
C               I8                (mapi.inc included in "mapic.inc")
C               I16               (mapi.inc included in "mapic.inc")
C               I32               (mapi.inc included in "mapic.inc")
C               I64               (mapi.inc included in "mapic.inc")
C               UI8               (mapi.inc included in "mapic.inc")
C               UI16              (mapi.inc included in "mapic.inc")
C               UI32              (mapi.inc included in "mapic.inc")
C               R32               (mapi.inc included in "mapic.inc")
C               R64               (mapi.inc included in "mapic.inc")
C               TXT               (mapi.inc included in "mapic.inc")
C               MFAIL             (mapi.inc included in "mapic.inc")
C               MAPIOK            (mapi.inc included in "mapic.inc")
C
C       Internals:
C
C          Variables:
C
C              DF2DT Return Values:
C              returns: MFAIL on error and MAPIOK if successful
C
C-----------------------------------------------------------------------
C!END
C
       IMPLICIT NONE
       INCLUDE 'mapic.inc'
 
            INTEGER            dfnt 
            CHARACTER*(*)      dtype

C      SET DF2DT TO MAPIOK, SUCCESSFUL MATCH
 
       DF2DT = MAPIOK
       dtype = ' '

C---------------------------------------------------------------------
C
C      CASE OF HDF NUMBER TYPE          RETURNS M-API data type
C
C           DFNT_INT8          set datatype to I8
C           DFNT_INT16         set datatype to I16
C           DFNT_INT32         set datatype to I32
C           DFNT_INT64         set datatype to I64
C           DFNT_FLOAT32       set datatype to R32
C           DFNT_FLOAT64       set datatype to R64
C           DFNT_UINT8         set datatype to UI8
C           DFNT_UINT16        set datatype to UI16
C           DFNT_UINT32        set datatype to UI32
C           DFNT_CHAR          set datatype to TXT
C           if no match        set datatype to an empty string and
C                              DF2DT to MFAIL, no match.
C
C---------------------------------------------------------------------
   
       IF (dfnt .EQ. DFNT_INT8) THEN
         dtype = I8
       ELSE IF (dfnt .EQ. DFNT_INT16) THEN
         dtype = I16
       ELSE IF (dfnt .EQ. DFNT_INT32) THEN
         dtype = I32
       ELSE IF (dfnt .EQ. DFNT_INT64) THEN
         dtype = I64
       ELSE IF (dfnt .EQ. DFNT_UINT8) THEN
         dtype = UI8
       ELSE IF (dfnt .EQ. DFNT_UINT16) THEN
         dtype = UI16
       ELSE IF (dfnt .EQ. DFNT_UINT32) THEN
         dtype = UI32
       ELSE IF (dfnt .EQ. DFNT_FLOAT32) THEN
         dtype = R32
       ELSE IF (dfnt .EQ. DFNT_FLOAT64) THEN
         dtype = R64
       ELSE IF (dfnt .EQ. DFNT_CHAR) THEN
         dtype = TXT
       ELSE
         dtype = ' '
         DF2DT = MFAIL
       END IF
  
C      RETURN DF2DT

       RETURN
       END
