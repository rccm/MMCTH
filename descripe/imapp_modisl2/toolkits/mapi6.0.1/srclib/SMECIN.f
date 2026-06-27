        FUNCTION SMECIN(cvalue,nelmnt,nstrs,substr)
C-----------------------------------------------------------------------
C!F77
C!Purpose:  Decomposes a multiple character string ECS metadata (retrieved by 
C GMECIN) into individual strings.
C 
C!Description:  ECS metadata values may be integer, floating point, or character 
C string values or arrays of values.  Some may be multiple strings.  The routine  
C GMECIN  retrieves such strings into a one-dimension character array with the 
C individual strings separated by nulls ('\0').   SMECIN  breaks this 'packed' 
C character array into its constituent substrings.  SMECIN copies these 
C substrings into separate rows of a FORTRAN character string array.
C 
C!Input parameters:
C  		cvalue	IN: 	Character string containing the 'packed' 	
C 				multiple substrings of ECS metadata retrieved with  
C 				GMECIN.
C  		nelmnt	IN: The composite output dimensions, from GMECIN,
C 				containing (in the case of character string metadata
C  				the total length (in bytes) of the string in cvalue in
C 				its lower two bytes and the number of substrings 	
C 				packed into  cvalue less one in the upper two bytes.  	
C 				The calculations
C 				int   n_strings  =  n_elements/65536  + 1
C 					n_bytes  =  n_elements%65536
C 				provide the number of substrings and the total 	
C 				length, respectively, of the data in  cvalue.  When 	
C 				there is only one string in  cvalue, nelmnt  will be 	
C 				less than 65536 and there is no need to use  SMECIN.
C 		nstrs		IN/OUT:   Number of  elements available  in the 	
C 				substr array.  The substr will not be set to the 		
C 				substrings in  cvalue  unless there are sufficient  	
C 				elements available in the  substr array.   SMECIN 	
C 				replaces this input with the number of substrings  
C 				already set   in the  cvalue  array.   nstrs will be set 	
C 				to 0 if a function error occurs.   .
C!Output parameters:
C 		substr		OUT:	Array of  substrings obtained from thecvalue.
C 				array.	
C  
C Returns:			MAPIOK (0) if each substring in cvalue  is in the 	
C 				substr  array,  MFAIL (-1) if substr  contains 		
C 				insufficient  elements for all of the substrings in 	
C 				cvalue  or an error occurs. 
C 
C External references:
C		MAPI_SMF_MAX_MSGBUF_SIZE	(mapic.inc)
C		pgs_smf_setdynamicmsg		(SDP Toolkit function)
C		MFAIL				(mapi.inc)
C		MAPIOK				(mapi.inc)
C
C!Revision History:
C 	Qi Huang, RDC		May 16, 1996
C	Version 2.0
C	Original development and testing
C 
C!Team-unique header:
C This software is developed by the MODIS Science Data SupportTeam for the 
C National Aeronautics and Space Administration, Goddard Space Flight Center, 
C under contract NAS5-32373.
C
C!REFERENCES AND CREDITS
C
C!DESIGN NOTES:
C
C-----------------------------------------------------------------------
C!END
C
      IMPLICIT NONE
      INCLUDE  'mapic.inc'
      INTEGER   nelmnt, nstrs
      CHARACTER*(*)     cvalue, substr(*)

C
C     Local variable declaration
C
C     lnstrs		local variable for number of substrings
C     i			outer WHILE loop counter
C     j			subscript of substr array
C     first		byte location of first non_blank character
C			of cvalue array
C     last		byte location of last non_blank characte
C			of cvalue array
C     n_bytes		total length (in bytes) of the string in
C			cvalue
C     bbyte		byte location of first character of the 
C			substring in cvalue array
C     len_substr_in     length of substr array (number of columns)
C     len_substr	actual length of the substring in cvalue array
C     returnstatus      return value of pgs_smf_setdynamicmsg
C     pgs_smf_setdynamicmsg
C                   SMF function that writes error messages to log
C     int_cvalue	return value of ICHAR
C     ICHAR		FORTRAN 77 standard library function to convert
C			a character to an integer

      INTEGER           lnstrs
      INTEGER           i,j
      INTEGER           n_bytes
      INTEGER           bbyte
      INTEGER           len_substr_in
      INTEGER           len_substr
      CHARACTER*(MAPI_SMF_MAX_MSGBUF_SIZE) BUF
      CHARACTER*(*)     FUNC
      PARAMETER (FUNC = 'SMECIN')
      INTEGER           returnstatus
      INTEGER           pgs_smf_setdynamicmsg

      SMECIN = MFAIL
C
C     Input checks
C
      if ( cvalue .eq. ' ') then
        write(BUF,'(A,A)') 'ERROR: SMECIN unable to continue without',
     +                     'cvalue input.'
        returnstatus = pgs_smf_setdynamicmsg(MAPI_E_ERR,BUF,FUNC)
        return
      endif

      if ( nstrs .lt. 1 ) then
        write(BUF,'(A,A)') 'ERROR: SMECIN unable to continue without',
     +                     ' nstrs input.'
        returnstatus = pgs_smf_setdynamicmsg(MAPI_E_ERR,BUF,FUNC)
        return
      endif

      if ( nelmnt .lt. 1 )  then
        write(BUF,'(A,A,I10,A)') 'ERROR: SMECIN unable to continue ',
     +                          'with invalid ', nelmnt,' nelmnt.'
        returnstatus = pgs_smf_setdynamicmsg(MAPI_E_ERR,BUF,FUNC)
        return
      endif

      lnstrs = nelmnt/65536 + 1
C     Calculate n_bytes.
      n_bytes = MOD(nelmnt,65536)

      if ( nstrs .lt. lnstrs ) then
        write(BUF,'(A,I10,A,I10,A)')'ERROR: SMECIN unable to fit ',
     +          lnstrs,' substrings into ',nstrs,
     +          ' element substring array.'
        returnstatus = pgs_smf_setdynamicmsg(MAPI_E_ERR,BUF,FUNC)
        nstrs = lnstrs
        return
      endif

C     Get the dimension size (number of columns) of substr array.
      len_substr_in = len(substr(1))

      bbyte = 1
      j = 1
C     Loop to break multistrings in cvalue array into substrings and 
C     copy the first (n_strs-1) substrings to  substr array.
      do i=1,n_bytes

        if ( ICHAR(cvalue(i:i)) .eq. 0 ) then

          len_substr = len(cvalue(bbyte:i-1))
          if (len_substr_in .lt. len_substr) then
            write(BUF,'(A,I10,A,I10,A)')'ERROR: SMECIN unable to fit ',
     +          len_substr,' byte string into ',len_substr_in,
     +          ' byte string array.'
            returnstatus = pgs_smf_setdynamicmsg(MAPI_E_ERR,BUF,FUNC)
            return
          else
            substr(j) = cvalue(bbyte:i-1)
            bbyte = i + 1
            j = j + 1
          end if

        end if

      end do
 
C     Copy the last substring to substr array.
      substr(j) = cvalue(bbyte:n_bytes)

      nstrs = lnstrs
      SMECIN = MAPIOK

      RETURN
      END      
