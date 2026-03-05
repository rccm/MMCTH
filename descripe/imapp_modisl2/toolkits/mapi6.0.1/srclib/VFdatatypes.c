#include <stdio.h>
#include <stdlib.h>
#include "mapic.h"

int VFdatatypes(int32 vdata_id, long int *stringlen, char *data_type) 
/*
!C**********************************************************************
* 
*!Purpose:     Retrieves a comma delimited string of the MODIS datatypes
*               of a Vdata's fields.           
* 
*!Description: Function VFdatatypes is part of a larger software sys-
*              tem called the MODIS Applications Programming Interface
*              (API) Utility, abbreviated M-API.  The M-API Utility con-
*              sists of subroutines which allow MODIS Science Team-sup-
*              plied software to read in Level 1B radiance bands and 
*              write out output products and metadata to HDF files.  The
*              functionality of the M-API is defined in the MODIS API 
*              User's Guide, Version 1, dated 4/3/95.
*
*              VFdatatypes retrieves the data type for each field of a 
*              Vdata.  The HDF number type of each field is converted 
*              into its corresponding M-API data type and this is appen-
*              ded to the output string.  An error (MFAIL) will be re-
*              turned if 1) The output string is not long enough to con
*              tain the data type strings for all the Vdata's fields,
*              2) an unknown (e.g. not supported by the MODIS API) num-
*              ber type is encountered or 3) an HDF routine FAILs.  The
*              data type string will be returned truncated to the point
*              where the fault occurred.  stringlen, the address of the
*              length of the data_type output string, is normally re- 
*              vised to contain the actual array length required to 
*              hold the output string.  If the latter two errors occur,
*              however, it is set to 0.
*
* !Input Parameters:
*       
*            vdata_id      Vdata's access identifier returned from 
*                          Vsattach.
*
*            stringlen     IN/OUT:   Address of the minimum length of the 
*                          data_type array.  Returns the minimum array
*			   length actually required to hold the data_type
*			   string.  It is set to 0 if a functional error 
*			   occurs.
*
* !Output Parameters:
*            data_type     OUT:   Array of comma-delimited data types
*                          for each table field.
*
*                          Possible C data types:
*                                 "int8"
*                                 "uint8"
*                                 "int16"
*                                 "uint16"
*                                 "int32"
*                                 "uint32"
*                                 "int64"
*                                 "float32"
*                                 "float64"
*
* Returns:                 MAPIOK if successful, MFAIL if an error 
*                          occurs.
*
* External references:
*                          DATATYPELENMAX      (mapic.h)
*                          MAPIOK              (mapi.h)
*                          MFAIL               (mapi.h)
*                          FAIL                (hdf.h)
*                          VFnfields           (vproto.h)
*                          VFfieldtype         (vproto.h)
*                          DFNT_to_datatype    (mapic.h)
*
* !Revision History:
*
* $Log: VFdatatypes.c,v $
* Revision 6.1  2010/07/13 19:29:18  kuyper
* Corrected format specifier and argument type to match each other.
*
* Revision 5.1  2005/04/04 16:50:23  vlin
* cast "data_typeMAXLEN" to long.
*
* Revision 1.1  1998/02/06 22:26:06  fshaw
* Initial revision
*
 * Revision 1.4  1996/02/26  22:01:46  qhuang
 * Changed the definition of stringlen as an in/out parameter in prolog.
 *
 * Revision 1.3  1995/11/07  18:27:06  qhuang
 * minor changes
 *
 * Revision 1.2  1995/10/31  14:56:23  qhuang
 * Added capability to pass status messages to log files.
 *
*
*              Revision 01.03 1995/08/30
*              Qi Huang
*              Add corrections per code walkthru sheet.
*
*              Revision 01.02 1995/08/28
*              Qi Huang
*              Version 1.3 original development. 
*
*              Revision 01.01 1995/05/03
*              Qi Huang
*              Add corrections per code walkthru sheet.
*
*              revision 01.00 1995/05/01
*              Qi Huang
*              Original develpment.
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
!END********************************************************************
*/

{
  char  buff[PGS_SMF_MAX_MSGBUF_SIZE];   /* buffer to hold the error/warning message */
  char *funcname="VFdatatypes";          /* name of this routine */
  long int  data_typeMAXLEN;             /* local variable to hold the
                                           size of array data_type */
  int  n;                                /* local variable to hold the
                                           the number of fields in the
                                           Vdata */
  int  i;                                /* variable to control FOR
                                             loop */
  int  accumulator=0;                      /* local variable to sum up
                                           the minimum length for an
                                           array capable of holding
                                           the whole data type string*/         
  char  field_datatype[DATATYPELENMAX]; /* array to contain a field
                                           data type in the Vdata */
  char  *comma=",";            /* character string which has a comma */

  int   output_status;                   /* variable its value will be
                                           returned to the routine */
  int32  VFnumtype;                   /* field number type in a Vdata*/

  if ( stringlen == NULL )
    return(MFAIL);

  data_typeMAXLEN= *stringlen;       /* hold the size of data_type */

  if ( data_type ==NULL )
  {
    *stringlen=0;
    return(MFAIL);
  }

  if ( *stringlen > 0)
    data_type[0]='\0';

  output_status= MAPIOK;                /* successful retrieval */

  if ( (n=VFnfields(vdata_id)) == FAIL) 
  { 
    sprintf(buff,"ERROR: VFdatatypes detected FAIL from HDF \n"
                          "\t routine VFnfields\n");
    MAPIERR(buff,funcname);
    output_status= MFAIL;               /* Fail */
  } 
  else
  {
    VFnumtype= VFfieldtype(vdata_id,0); /* Get the first field's number                                                  type */

    /* If the field's number type conversion fails, ... */
    if (DFNT_to_datatype (VFnumtype,field_datatype) == FAIL )
    {
      sprintf(buff,"ERROR: VFdatatypes detected unrecognized HDF\n"
                            "\t number type %ld\n", (long)VFnumtype);
      MAPIERR(buff,funcname);
      output_status=MFAIL;                /* Fail */
    }
    else
    {
      accumulator = (int)strlen(field_datatype) + 1; /* Set the accumulator
                                                to the length required 
                                                to hold the first 
                                                field's data type */
      if ( accumulator <= data_typeMAXLEN )
        strcpy(data_type,field_datatype);    /* Copy the data type into                                                 the data_type output string */
    }

    /* For each remaining field in the Vdata AND while the status code
        is still MAPIOK */
    for ( i=1; (i<n) && (output_status == MAPIOK); i++)
    {
      VFnumtype= VFfieldtype(vdata_id,i); /* Get the ith  field's 
                                            number type */
      /* If the field's number type conversion fails, ... */
      if (DFNT_to_datatype (VFnumtype,field_datatype) == FAIL )
      {
        sprintf(buff,"ERROR: VFdatatypes detected unrecognized HDF\n"
                                "\t number type %ld\n", (long)VFnumtype);
        MAPIERR(buff,funcname);
        output_status=MFAIL;                /* Fail */
      }
      else
      {
        /* Update the accumulator with the additional length required
           to hold the field's data type + a comma. */
        accumulator += (int)strlen(field_datatype) + 1;
        if ( accumulator <= data_typeMAXLEN)
        {
          /* Copy the data type into the data_type output string, 
            preceeded by a comma. */
          strcat(data_type,comma);
          strcat(data_type,field_datatype);
        }
      }
    }
  }

  if ( output_status == MFAIL )
    *stringlen = 0;
  else	
  {
    if ( data_typeMAXLEN < accumulator ) /* the whole data type string
                                         will NOT fit in the data_type                                                  output string */
    { 
      sprintf(buff,"ERROR: VFdatatypes unable to fit a \n"
                     "\t%d byte data type string into\n"
              "\t a %ld byte array\n",accumulator, (long)data_typeMAXLEN);
      MAPIERR(buff,funcname);
      output_status=MFAIL;               /* Fail */
    }

    *stringlen=accumulator;              /* Set *stringlen to the
                                         quantity in the accumulator */
  }

  return(output_status);
}
