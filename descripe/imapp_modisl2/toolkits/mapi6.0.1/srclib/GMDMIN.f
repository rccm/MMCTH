      FUNCTION GMDMIN(modfil, arrnm, grpnm, moddim, attr, dtype,
     &     nms, value)
      IMPLICIT NONE  
      INCLUDE 'mapic.inc'

C--------------------------------------------------------------------------
C !F77
C 
C !Purpose: Reads the value(s) of a local attribute attached to a specific dimension
C           of an SDS (array) from a MODIS HDF file.
C
C !Description:	Function GMDMIN is part of a larger software system called the MODIS
C               Applications Programming Interface (API) Utility, abbreviated M-API. 
C               The M-API Utility consists of subroutines which allow MODIS Science
C               Team-supplied software to read  and write data and metadata from/to
C               HDF files. The functionality of the M-API is defined in the MODIS
C               Application Program Interface (API) Specification.
C
C               GMDMIN retrieves the value(s) associated with an attribute = values
C               pair for a dimension attribute by giving the array's name and
C               dimension number and attribute name.  If the attribute cannot be
C               found, the routine will return MFAIL and the passed variable unchanged.
C
C               The routine will also fail if the provided dtype is found to be different
C               from the attribute's data type or the number of elements (nms) is found
C               to be too small to contain the attribute's value.  GMDMIN replaces this
C               input information with the actual data type and number of elements
C               contained in the attribute value (in the case of character data, it is
C               the length of the string).  These attribute's attribute may be used to
C               properly retrieve the attribute value with a second call to the routine.
C               If a function failure occurs, or specified array or dimension does not
C               exist,  nms will be set to zero.
C
C               A variable of the proper data type should be passed for the value
C               parameter.  The data type information required to properly use either
C               routine may be found in Appendix A, M-API-Supplied Constants, and
C               Appendix C, MODIS Data Product File Definitions of M-API User's Guide.
C               Appendix A has a listing for each M-API provided attribute that
C               includes the data type, the format, and/or specific values associated
C               with it.
C
C !Input parameters:
C     modfil  IN:  array  that is used to reference a MODIS HDF file containing the
C                  attribute.
C     arrnm   IN:  ASCII string name of the array.  Provided macros for accepted MODIS
C                  HDF file array names are listed in Appendix A, M-API-Supplied
C                  Constants.
C     grpnm   IN:  ASCII string name of the data group containing the array structure
C                  to which the attribute is attached.  If set to NULL the entire file
C                  will be searched for the array structure named arrnm.
C     moddim     IN:  The dimension number from which the attribute values will be
C                  retrieved (0-based).
C     attr    IN:  ASCII string name of the attribute.  Provided macros for accepted
C                  MODIS HDF file attribute names are listed in Appendix A, M-API
C                  Supplied Constants.
C     dtype IN/OUT: ASCII string of data type of the value output.  Output replaces
C                   with the data type of the retrieved attribute. The memory size
C                   of dtype should be at least 13 characters long. 
C     nms   IN/OUT: The number of elements available  in the value array.  Output 
C                   replaces with the number of elements required to contain the
C                   attribute. If a function failure occurs, the value will be set
C                   to zero.
C
C !Output parameters:
C     value OUT:  Values associated with the attribute.
C
C !Returns: MAPIOK if successful, MFAIL if value cannot contain the retrieved
C           attribute value, the data type is different, the attr cannot be found,
C           or an error occurs.
C
C !Revision History:
C
C !Team-unique header:
C     This software is developed by the MODIS Science Data Support Team for the
C     National Aeronautics and Space Administration, Goddard Space Flight Center,
C     under contract NAS5-32373.
C--------------------------------------------------------------------------------
C !END
C
      INTEGER modfil(MODFILLEN), nms, moddim
      CHARACTER*(*) grpnm, arrnm, dtype, attr
      BYTE value(*)
      INTEGER ret

      CALL cgmdmin(modfil, arrnm, len(arrnm), grpnm, len(grpnm), moddim,
     &       attr, len(attr), dtype, len(dtype), nms, value, ret)

      GMDMIN = ret 
      return
      end

