      FUNCTION DT2DF(dtype)
C
C-------------------------------------------------------------------------------
C!F77
C
C!Description:  Subroutine DT2DF is part of a larger software system called the
C               MODIS Applications Programming Interface (API) Utility,
C              	abbreviated M-API.  The M-API Utility consists of subroutines 
C		which allow MODIS Science Team-supplied software to read in 
C		Level 1B radiance bands and write out output products and 
C		metadata to HDF files. 
C
C	        DT2DF takes one of the following M-API data type strings and 
C		returns the HDF number type integer used by the HDF routines 
C		to identify number types.
C
C !Input Parameters:CHARACTER dtype  :  M-API data type string
C
C
C				         Valid data types :
C						`INTEGER*1`
C                                       	`INTEGER*2`
C                                       	`INTEGER*4`
C "COMMENTED OUT: NOT SUPPORTED BY HDF 3.3r4"   `INTEGER*8`
C                                             * 'UINTEGER*1`
C                                      l      * 'UINTEGER*2`
C                                             * 'UINTEGER*4`
C "COMMENTED OUT: NOT SUPPORTED BY HDF 3.3r4" * 'UINTEGER*8`
C                                       	`REAL*4`
C                                       	`REAL*8`
C                                       	`CHARACTER*(*)`
C * (unsigned integer)
C
C !Output Parameters:NONE
C
C !Revision History:
C
C $Id: DT2DF.f,v 1.1 1999/04/15 19:26:51 jayshree Exp $
C
C $Log: DT2DF.f,v $
C Revision 1.1  1999/04/15 19:26:51  jayshree
C Initial revision
C
C
C !Team-unique Header:
C
C               This software is developed by the MODIS Science Data Support
C            	Team for the National Aeronautics and Space Administration,
C             	Goddard Space Flight Center, under contract NAS5-32373.
C
C             	Portions developed at the National Center for Supercomputing
C               Applications at the Univ. of Illinois at Urbana-Champaign.
C
C !References and Credits:
C
C     		WRITTEN BY:   
C
C		Dave Lorenzi  05/16/95
C               Original Development / Testing
C
C		Vicky Lin          06/08/95
C		Modification / Testing
C		vlin@modis-xl.gsfc.nasa.gov
C
C               Mitchell Weiss     07/26/95
C               Modification / Testing
C
C !Design Notes:
C
C	Externals:
C
C          Functions:
C               DT2DF           (mapic.inc)
C
C	   Variables:
C		DFNT_INT8	(hdf.inc: included in "mapic.inc")	
C		DFNT_INT16	(hdf.inc: included in "mapic.inc")
C		DFNT_INT32	(hdf.inc: included in "mapic.inc")
C NOT USED      DFNT_INT64	(hdf.inc: included in "mapic.inc")
C               DFNT_UINT8      (hdf.inc: included in "mapic.inc")
C               DFNT_UINT16     (hdf.inc: included in "mapic.inc")
C               DFNT_UINT32     (hdf.inc: included in "mapic.inc")
C NOT USED      DFNT_UINT64     (hdf.inc: included in "mapic.inc")
C		DFNT_FLOAT32	(hdf.inc: included in "mapic.inc")
C		DFNT_FLOAT64	(hdf.inc: included in "mapic.inc")
C		DFNT_CHAR	(hdf.inc: included in "mapic.inc")
C               I8              (mapi.inc included in "mapic.inc")
C               I16             (mapi.inc included in "mapic.inc")
C               I32             (mapi.inc included in "mapic.inc")
C               I64             (mapi.inc included in "mapic.inc")
C               UI8             (mapi.inc included in "mapic.inc")
C               UI16            (mapi.inc included in "mapic.inc")
C               UI32            (mapi.inc included in "mapic.inc")
C               UI64            (mapi.inc included in "mapic.inc")
C               R32             (mapi.inc included in "mapic.inc")
C               R64             (mapi.inc included in "mapic.inc")
C               TXT             (mapi.inc included in "mapic.inc")
C               MFAIL           (mapi.inc: included in "mapic.inc")
C
C	Internals:
C
C       Variables:
C
C	        DT2DF Return Values:
C               returns: MFAIL on error, HDF number type if successful
C
C-----------------------------------------------------------------------
C!END
C
      IMPLICIT NONE
      INCLUDE 'mapic.inc'
      CHARACTER*(*) dtype

C-----------------------------------------------------------------------
C
C       CASE OF M-API DATA TYPE STRING    RETURNS HDF # type (integer)
C
C               I8	set to HDF NUMBER TYPE DFNT_INT8
C               I16     set to HDF NUMBER TYPE DFNT_INT16
C               I32     set to HDF NUMBER TYPE DFNT_INT32
C               UI8     set to HDF NUMBER TYPE DFNT_UINT8
C              UI16     set to HDF NUMBER TYPE DFNT_UINT16
C              UI32     set to HDF NUMBER TYPE DFNT_UINT32
C               R32     set to HDF NUMBER TYPE DFNT_FLOAT32
C               R64     set to HDF NUMBER TYPE DFNT_FLOAT64
C               TXT     set to HDF NUMBER TYPE DFNT_CHAR
C
C-----------------------------------------------------------------------
C  initialize DT2DF to MFAIL
       DT2DF = MFAIL    

       IF (dtype.eq.I8) THEN
          DT2DF = DFNT_INT8

       ELSE IF (dtype.eq.I16) THEN
          DT2DF= DFNT_INT16

       ELSE IF (dtype.eq.I32) THEN
          DT2DF = DFNT_INT32

CCC    ELSE IF (dtype.eq.I64) THEN 
CCC       DT2DF = DFNT_INT64

       ELSE IF (dtype.eq.UI8) THEN
          DT2DF= DFNT_UINT8

       ELSE IF (dtype.eq.UI16) THEN
          DT2DF = DFNT_UINT16

       ELSE IF (dtype.eq.UI32) THEN
          DT2DF= DFNT_UINT32

CCC    ELSE IF (dtype.eq.UI64) THEN
CCC       DT2DF = DFNT_UINT64

       ELSE IF (dtype.eq.R32) THEN
          DT2DF = DFNT_FLOAT32

       ELSE IF (dtype.eq.R64) THEN
          DT2DF = DFNT_FLOAT64

       ELSE IF (dtype.eq.TXT) THEN
          DT2DF = DFNT_CHAR

       ELSE
          DT2DF = MFAIL

       ENDIF
 
       RETURN 
       END
