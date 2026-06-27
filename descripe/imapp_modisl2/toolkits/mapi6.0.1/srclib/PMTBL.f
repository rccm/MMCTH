      FUNCTION PMTBL(mfile, tbname, group, start, recno, ddata)
C
C--------------------------------------------------------------------------
C!F77
C
C!Description: Function PMTBL is part of a large software system called the 
C              MODIS Application Programming Interface(API) utility, 
C              abbreviated M-API.  The M-API Utility consists of subroutines 
C              which allow MODIS science Team-supplied software to read and 
C              write data and metadata to HDF files.  The functionality of 
C              the M-API is defined in the MODIS API specification.
C
C              PMTBL places one or more data record into an HDF Vdata table 
C              structure previously created using CRMTBL. The data to be 
C              inserted into a buffer array.  The length of this stored in the
C              Vdata byte array must be an integral number of the table 
C              structure's record length.  The various data that make up a 
C              record should be inserted into the buffer in the same order as 
C              the field headers were ordered in the CRMTBL call.  See Chapter
C              7 of the User's Guide, "Accessing Tables" for additional
C              information.  This routine may be called multiple times
C              to fill the table structure.  Data previously in the 
C              table structure may be overwritten.
C
C              An empty Vdata may not be created, so a dummy record is
C              inserted into it.  This dummy record should be overwritten
C              with the first call from PMTBL.  If this initial write is
C              performed with a single record an ambiguity is resolved 
C              by the creation of a semaphore metadata to indicate that
C              the dummy record has been overwritten.
C
C              The group string provides the facility to select a table
C              structure placed in a particular HDF 'Vgroup' data group.
C              The entire file will be searched for a table structure
C              named tbname if group = ' '
C
C!Input Parameters:
C mfile        Array that is used to reference the MODIS-HDF file.
C tbname       ASCII string name of the table structure.
C group        ASCII string name of the data containing the target structure.
C              If set to ' ', the entire file will be searched for the table
C              structure named tbname.
C start        Zero-based record location to begin putting the data into
C              the table structure.
C recno        Number of records to put into the table structure.
C ddata         Address of the data buffer.
C
C!Output Parameters:None.
C
C Returns:     MAPIOK if successful, MFAIL if an error occurs.
C
C External reference:
C              ncpmtbl           (mapic.h)
C
C!Revision History:
C   
C
C
C !Team-unique Header:
C             This software is developed by the MODIS Science Data Support
C             Team for the National Aeronautics and Space Administration,
C             Goddard Space Flight Center, under contract NAS5-32373.
C             Portions developed at the National Center for Supercomputing
C             Applications at the Univ. of Illinois at Urbana-Champaign.
C
C !References and Credits:
C     WRITTEN BY:    Xia W. Chang,  xia@ltpmail.gsfc.nasa.gov
C                    from prototype by Liping Di.
C                    MODIS Science Data Support Team
C                    Maryland Corporate Center
C                    7501 Forbes Blvd - Suite 103
C                    Seabrook, MD  20706
C
C----------------------------------------------------------------------------
C!END
C
      IMPLICIT NONE
      INCLUDE  'mapi.inc'
      INTEGER        mfile(MODFILLEN)
      CHARACTER*(*)  tbname
      CHARACTER*(*)  group
      INTEGER        start
      INTEGER        recno
      BYTE           ddata(*)

      INTEGER    ret
C
      call cpmtbl(mfile, tbname, len(tbname), group, len(group),
     *            start, recno, ddata, ret)
      PMTBL = ret
C
      RETURN
      END
