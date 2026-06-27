      fUNCTION GMTBL(mfile, tbname, group, field, 
     *               start, recno,  bsize, ddata)

C--------------------------------------------------------------------------
C!F77
C
C!Description: Function GMTBL is part of a large software system called the 
C              MODIS Application Programming Interface(API) utility, 
C              abbreviated M-API.  The M-API Utility consists of subroutines 
C              which allow MODIS Science Team-supplied software to read in 
C              Level 1B radiance bands and write out output products and 
C              metadata to HDF files.  The functionality of the M-API is 
C              defined in the MODIS API User's Guide, Version 1, dated 4/3/95.
C
C              GMTBL retrieves one or more fields of data from one or more 
C              records stored in an HDF Vdata table structure contained in a 
C              MODIS-HDF file.  The data are placed in the data buffer in
C              consecutive records and in the order that the input fields are
C              listed.  The length of this buffer must be able to contain
C              all the fields requested times the number of records requested.
C              If the bsize input indicates that it is too small to contain
C              the actual quantity of data requested, GMTBL will fail, but
C              it will return the actual bsize required.  The output data
C              buffer must be at least this size.
C
C              The group string provides the facility to select a table
C              structure placed in a particular HDF 'Vgroup' data group.
C              Alternatively, the entire file will be search for a table
C              structure named tbname if group is set to ' '.
C
C!Input Parameters:
C mfile        Array that is used to reference the MODIS-HDF file.
C tbname       ASCII string name of the source table structure.
C group        ASCII string name of data group containing the source table
C              structure.  If set to ' ', the entire file will be search for
C              the table structure named tbname.
C field        Array of comma-delimited ASCII string table headers.  The data
C              from each field will appear in the same order as the headers.
C start        Zero-based record location to begin reading the data from the
C              table structure.
C recno        Number of records to retrieve from the table structure.
C
C bsize        data buffer size out input, in bytes.  The buffer must be
C              at least this size. bsize will normally return the number
C              of bytes of data successfully retrieved.  If the buffer is 
C              too small, however, the routine returns MFAIL and bsize will
C              contain the size a buffer must be to contain the output data.
C              It is set to 0 if a functional error occurs and makes this
C              output size determination unreliable.
C 
C !Output Parameters: 
C ddata         the data buffer.
C
C
C Returns:
C              MOPIOK if successful, MFAIL if an error occurs.
C
C External references:
C              ncmtbl           (mapic.h)
C
C!Revision History:
C$Log: GMTBL.f,v $
CRevision 1.1  1999/04/15 19:26:51  jayshree
CInitial revision
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
C     WRITTEN BY:    Xia W. Chang,  xia@ltpmail.gsfc.nasa.gov
C                    From Prototype by Liping Di
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
      CHARACTER*(*)  tbname
      CHARACTER*(*)  group
      CHARACTER*(*)  field
      INTEGER        start
      INTEGER        recno
      INTEGER        bsize
      BYTE           ddata(*)
C     Declare local variable
C
      INTEGER        ret
C
      call cgmtbl(mfile, tbname, len(tbname), group, len(group),
     *           field, len(field), start, recno, bsize, ddata, ret)
      
      GMTBL = ret

      RETURN
      END
