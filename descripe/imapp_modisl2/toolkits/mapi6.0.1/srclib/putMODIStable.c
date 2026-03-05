#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mapic.h"

int putMODIStable(MODFILE *file, char const *tablename, char const *groupname,
long int start, long int recno, unsigned char const *data)

/*
!C**********************************************************************
* 
*!Purpose:     Inserts data records into a MODIS-HDF file table structure. 
* 
*!Description: Function putMODIStable is part of a larger software sys-
*              tem called the MODIS Applications Programming Interface
*              (API) Utility, abbreviated M-API.  The M-API Utility con-
*              sists of subroutines which allow MODIS Science Team-sup-
*              plied software to read and write data and metadata from/to
*              HDF files.  The functionality of the M-API is defined in 
*              the MODIS API Specification.
*
*              putMODIStable (PMTBL) places one or more data records into
*              an HDF Vdata table structure previously created using  
*              createMODIStable (CRMTBL).  The data to be inserted into 
*              the table must be inserted into a buffer string.  The length
*              of this string must be an integral number of the table stru-
*              ture's record length.  The various data that make up a re-
*              cord should be inserted into the buffer in the same order 
*              as the field headers were ordered in the createMODIStable 
*              call.  See Chapter 7 of the User's Guide, "Accessing Tables"
*              for additional information.  This routine may be called mul-
*              tiple times to fill the table structure.  Data previously in
*              the table structure may be overwritten.
* 
*              An empty Vdata may not be created, so a dummy record is in-
*              serted into it.  This dummy record should be overwritten 
*              with the first call from putMODIStable.  If this initial 
*              write is performed with a single record an ambiguity would 
*              exist whether this record was actual data or the dummy.  The
*              ambiguity is resolved by the creation of a semaphore meta-
*              data to indicate that the dummy record has been overwritten.
*
*              The groupname string provides the facility to select a table
*              structure placed in a particular HDF 'Vgroup' data group.
*              The entire file will be searched for a table structure named
*              tablename if grouopname =NULL in C and grpnm ='' in FORTRAN.
*
* !Input Parameters:
*            file       IN:   Address of MODFILE structure that is used to
*                       reference the MODIS-HDF file containing the target
*                       table structure.
*            tablename  IN:   ASCII string name of the target table
*                       structure.
*            groupname  IN:   ASCII string name of the data group contai-   
*                       ning the target table structure.  If set to NULL
*                       the entire file will be searched for the table 
*                       structure named tablename.
*            start      IN:  Zero-based record location to begin placing
*                       the data into the table structure.  If start = -1
*                       data records will be appended to the end of the   
*                       table structure.
*            recno      IN:   Number of records being inserted into the
*                       table structure.  The product of recno ane the
*                       table structure's record length must have the same
*                       length as the buffer addressed by data.
*            data       IN:   Address of the data buffer.
*
* !Output Parameters:NONE.
*         
* Returns:                 MAPIOK if successful, MFAIL if an error 
*                          occurs.
*
* External references:
*			   MODFILE			(mapi.h)
*			   PGS_SMF_MAX_MSGBUF_SIZE	(PGS_SMF.h)
*			   MAPIERR			(mapic.h)
*			   NULLMODFIL			(mapic.h)
*			   VSNAMELENMAX			(hdf.h)
*                          MFAIL               		(mapi.h)
*                          MAPIOK              		(mapi.h)
*			   getMODIStableid		(mapic.h)
*			   APPEND_VDATA			(mapic.h)
*			   FULL_INTERLACE		(hdf.h)
*			   VSwrite			(vproto.h)
*                          NULLstr             		(mapic.h)
*                          DFACC_READ          		(hdf.h)
*                          VSinquire           		(vproto.h)
*			   VSQueryref			(vg.h)
*                          emptyVdata          		(mapic.h)
*                          VSread              		(vproto.h)
*                          VSseek              		(vproto.h)
*                          set_Vhasdata        		(mapic.h)
*			   MAPIWARN			(mapic.h)
*
* !Revision History:
*		Qi Huang	1996/09/26
*		Version 2.2
*		Implementing ring super structure for the Vdata.
*
*		Qi Huang  1996/02/06
*		M-API Version 2.0
*
* $Log: putMODIStable.c,v $
* Revision 5.1  2005/04/04 18:49:11  vlin
* constant safe for pointer arguments.
*
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
 * Revision 1.5  1996/10/28  18:49:10  qhuang
 * Version 2.2
 *
 * Revision 1.4  1996/02/06  21:25:34  qhuang
 * Version 2.0. Added group handling, new messages and input checking.
 *
 * Revision 1.4 1995/12/29   Xia W. Chang, xia@modis-xl.gsfc.nasa.gov
 * Replace the HDF routine VSfind to mapi routine MVSfind
 *
 * Revision 1.3  1995/11/09  20:58:08  qhuang
 * minor changes
 *
 * Revision 1.2  1995/11/01  18:48:23  qhuang
 * Added capability to pass status messages to log files.
 *
 * Revision 1.1  1995/11/01  18:47:27  qhuang
 * Initial revision
 *
*
*              Revision 01.00 1995/9/26
*              Qi Huang
*              Version 1.3 original development/testing.
*
*              Liping Di, Hughes STX  Oct., 17, 1996
*              Using a local int32 variable nrecs32 intead of long int nrecs32 in the call
*              to VSinquire for the portability.
*
* !Team-unique Header:
*       This software is developed by the MODIS Science Data Support
*       Team for the National Aeronautics and Space Administration,
*       Goddard Space Flight Center, under contract NAS5-32373.
* 
*      HDF portions developed at the National Center for Supercomputing
*      Applications at the University of Illinois at Urbana-Champaign.
*
* !References and Credits:
*
* !Design Notes:
*
!END*********************************************************************
*/

{
  char  buff[PGS_SMF_MAX_MSGBUF_SIZE];  /* buffer to hold the error/warning
                                           message */
  char *funcname="putMODIStable";       /* name of this module */

  DATAID	*did;
  int32		nrecs;             /* local variable to hold the number of 
                                                 records in a Vdata */
  int32		recsize;           /* The size in bytes of a record 
                                     in the Vdata for the local machine.*/
  int		status;            /* variable its value will be returned
                                     to the calling routine */
  int32         ref_id;            /* Vdata reference number */
  uint8 	*databuf;          /* temporal buffer */

  status = MFAIL;                   /* error */

  /* Input checks: */
  if ( NULLstr(tablename) )
  {
    sprintf(buff,"ERROR: putMODIStable unable to write to a table\n"
			"\t without a table name input.\n");
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( NULLMODFIL(file) )
  {
    sprintf(buff,"ERROR: putMODIStable unable to write to the %.*s\n"
			"\t table with a NULL file MODFILE structure.\n",
			VSNAMELENMAX,tablename);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( data == NULL )
  {
    sprintf(buff,"ERROR: putMODIStable unable to write to the %.*s\n"
			"\t table without a data buffer.\n",
			VSNAMELENMAX,tablename);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( file->access == DFACC_READ )
  {
    sprintf(buff,"ERROR: putMODIStable unable to write to the %.*s\n"
			"\t table in file opened for read only.\n",
			VSNAMELENMAX,tablename);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  if ( (did = getMODIStableid(file,tablename,groupname,"w")) == NULL )
  {
    sprintf(buff,"ERROR: putMODIStable detected FAIL from procedure\n"
			"\t getMODIStableid attempting to access the\n"
			"\t %.*s table.\n",VSNAMELENMAX,tablename);
    MAPIERR(buff,funcname);
    return(MFAIL);
  }

  /* Check to determine that write will append or overwrite the existing
     Vdata records: */

  /* Get the number of records(nrecs) in the Vdata and the record size(recsize)
   using VSinquire.  */
  if ( VSinquire((int32)did->id,&nrecs,NULL,NULL,&recsize,NULL) == FAIL ) 
  {
    sprintf(buff,"ERROR: putMODIStable detected FAIL from HDF procedure\n"
                          "\t VSinquire while attempting to write to the\n"
			  "\t %.*s table.\n",VSNAMELENMAX,tablename);
    MAPIERR(buff,funcname);
  }
  
  else  /* VSinquire == SUCCESS */
  {
    ref_id = VSQueryref((int32)did->id);
    if ( (nrecs == 1) && (emptyVdata(file,ref_id) == 1) )
      nrecs = 0;   /* Set the number of records currently in
                             the Vdata to 0. */

    if (start == APPEND_VDATA)
      start = nrecs;

    /* If start, the requested record position to begin writing data 
      to the Vdata, is greater than the number of records currently
      in the Vdata OR less than 0. */
    if ( (start > nrecs) || (start < 0) )
    {
      sprintf(buff,"ERROR: putMODIStable unable to write data to the\n"
                            "\t %.*s table to invalid table structure\n"
                            "\t record %ld\n",VSNAMELENMAX,tablename,start);
      MAPIERR(buff,funcname);
    }
 
    else  /* start <= nrecs  && start >= 0 */
    { 
      if ( (start == 1) && (nrecs == 1) )
      { /* (NOTE:Perform special Vdata append procedure required to avoid
                 reported VSseek bug that prevents it from pointing to
                 appending position when there is only one record in the
                 Vdata. It is not necessary to set the field when reading
                 whole record) */
        if ( (databuf = (uint8 *)malloc((size_t)recsize)) == NULL )
        {
          sprintf(buff,"ERROR: putMODIStable memory allocation\n"
                                "\t failure while attempting to write\n"
				"\t data to the %.*s table.\n",
				VSNAMELENMAX,tablename);
          MAPIERR(buff,funcname);
        }
     
  	else if (VSseek((int32)did->id,0) == FAIL)
	{
          free(databuf);
	  sprintf(buff,"ERROR: putMODIStable detected FAIL from\n"
                             "\t HDF procedure VSseek while\n"
                             "\t attempting to write to the %.*s\n"
                             "\t table.\n",VSNAMELENMAX,tablename);
          MAPIERR(buff,funcname);
        }
 
        else if ( VSread((int32)did->id,databuf,1,FULL_INTERLACE) == FAIL )
        {
          free(databuf);	
          sprintf(buff,"ERROR: putMODIStable detected FAIL from\n"
                             "\t HDF procedure VSread while\n"
			     "\t attempting to write to the %.*s\n"
			     "\t table.\n",VSNAMELENMAX,tablename);
          MAPIERR(buff,funcname);
        }

        else
        {
          free(databuf);
          status = MAPIOK;
        }
      } /* end of if; if ( (start == 1) && (nrecs == 1) ) */

      else /* start != 1 || nrecs != 1 */
      {
        if ( VSseek((int32)did->id,(int32)start) == FAIL )
        {
          sprintf(buff,"ERROR: putMODIStable detected FAIL from\n"
                                "\t HDF procedure VSseek attempting to\n"
                                "\t write  to the %.*s table.\n",
				VSNAMELENMAX,tablename);
          MAPIERR(buff,funcname);
        }

        else
          status = MAPIOK;
      }

      if ( status == MAPIOK )
      {
        if ( VSwrite((int32)did->id,data,(int32)recno,FULL_INTERLACE) == FAIL )
        {
          sprintf(buff,"ERROR: putMODIStable detected FAIL from\n"
                                "\t HDF procedure VSwrite while\n"
				"\t attempting to write to the %.*s\n"
				"\t table.\n",VSNAMELENMAX,tablename);
          MAPIERR(buff,funcname);
          status = MFAIL;
        }

        else 
	  if ( (nrecs == 0) && (recno == 1) ) /* there were no real records in the Vdata
						 AND only a single record was written */
            if ( set_Vhasdata(file,ref_id) == MFAIL )
          {
            sprintf(buff,"WARNING: putMODIStable detected MFAIL\n"
                                  "\t from M-API procedure\n"
                                  "\t set_Vhasdata.\n");
            MAPIWARN(buff,funcname);
          }

      } /* end of if; if ( status == MAPIOK ) */

    }/* end of else; start <= nrecs && start >= 0 */

  }/* end of else; VSinquire == SUCCESS */

  return(status);
}  
