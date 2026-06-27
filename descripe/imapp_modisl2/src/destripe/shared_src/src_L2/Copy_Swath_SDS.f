      integer function Copy_Swath_SDS( sd_id_in, 
     1                                 sd_id_out, 
     2                                 sds_name_in,
     3                                 sds_name_out, 
     4                                 fid_out, 
     5                                 swathid_out, 
     6                                 swathname,
     7                                 dimlist )


C-------------------------------------------------------------------------------
C!F77
C
C!Description: Copy one SDS array and all its attributes from one HDF-EOS file 
C              (input) to another (output).  Only two and three dimensional 
C              numeric data arrays are supported.  
C   
C!Input Parameters: 
C
C integer       fid_out          HDF-EOS file ID to output file 
C integer       sd_id_out        HDF SDS interface ID to output file 
C integer       sd_id_in         HDF SDS interface ID to input file 
C character*(*) sds_name_out     HDF SDS name in output file 
C character*(*) sds_name_in      HDF SDS array name in input file 
C integer       swathid_out      HDF-EOS Swath ID in output file 
C character*(*) swathname        HDF-EOS Swath name in output file 
C character*(*) dimlist          Comma delimited list of array dimension 
C                                names
C
C!Output Parameters: None
C
C!Revision History:
c Revision 1.1  1999/02/08  17:00:05  vlin
c Initial revision
c
C
C!Team-Unique Header:
C This software was developed by the MODIS Science Data Support Team
C (SDST) for the National Aeronautics and Space Administration,
C Goddard Space Flight Center, under contract NAS5-32373.
C
C HDF functions were developed at the National Center for Supercomputing
C Applications at the Univ. of Illinois at Urbana-Champaign.
C
C!References and Credits:
C   Developed in May 1997 by Fay Liang fhliang@ltpmail.gsfc.nasa.gov
C                        and Rich Hucek rhucek@ltpmail.gsfc.nasa.gov
C
C!Externals 
C   Functions:
C    swattach                (SWapi.c)
C    swdefdfld               (SWapi.c)
C    swdefgfld               (SWapi.c)
C    swdetach                (SWapi.c)
C    sfendacc                (dffunc.inc)
C    sfginfo                 (dffunc.inc)
C    sfn2index               (dffunc.inc)
C    sfselect                (dffunc.inc)
C    Copy_SDS                (science code)
C    string_loc              (science code)
C
C  Named Constants:
C    FAIL                    (hdf.inc)
C    MODIS_E_GENERIC         (PGS_MODIS_39500.f)
C    SUCCEED                 (hdf.inc)
C
C!Internals: (SDST code)
C  Variables:
C  data_type          Data type for data/attribute in an SDS
C  dimsizes           Array to store size of each dimension in an SDS
C                     be read along each dimension
C  fbyte_dim          1st byte of nonblank text of dimlist
C  fbyte_sds_out      1st byte of nonblank text of sds_name_out
C  fbyte_sds_in       1st byte of nonblank text of sds_name_in
C  lbyte_dim          Last byte of nonblank text of dimlist
C  lbyte_sds_out      Last byte of nonblank text of sds_name_out
C  lbyte_sds_in       Last byte of nonblank text of sds_name_in
C                     data_buf
C  num_attr           Number of attributes in an SDS
C  rank               Number of dimensions in an SDS
C  rtn                Return code from a function
C                     each dimension
C  sds_id_out         HDF SDS ID in output file 
C  sds_id_in          HDF SDS ID in input file 
C  sds_index_out      SDS index in output file 
C  sds_index_in       SDS index in input file
C  msgbuf             Character string for log messages
C  word_size          Size, in byte, of data_type
C
C  Named Constants:
C  HDFE_NOMERGE       Code for 'no merge' action
C  HDFE_AUTOMERGE     Code for 'merge' action
C  FUNCNAME           Character string for function names
C    
C
C  Functions and Subroutines:
C   string_loc       Finds positions of both 1st & last nonblank
C                    character of a character string
C
C!END-------------------------------------------------------------------

      IMPLICIT NONE

      include 'hdf.inc'
      include 'PGS_MODIS_39500.f'

C-----FORTRAN PARAMETER definitions
      character*(*)  BLANK
      parameter     (BLANK = ' ')

      character*(*)  FUNCNAME
      parameter     (FUNCNAME = 'Copy_Swath_SDS')

      integer        HDFE_NOMERGE        ! Symbolic constants for HDF-EOS
      parameter     (HDFE_NOMERGE = 0)

      integer        HDFE_AUTOMERGE
      parameter     (HDFE_AUTOMERGE = 1)


C-----Function Argument Declarations
      character*(*) dimlist, sds_name_in, sds_name_out, swathname
      integer fid_out, sd_id_in, sd_id_out, swathid_out


C-----Local Variable Declarations
      character*9    msg9
      character*2048 msgbuf

      integer sfn2index, sfendacc, sfselect, sfginfo, swdefdfld, swdefgfld,
     1        swdetach,  swattach, get_sds_id, copy_sds, string_loc

      integer fbyte_sds_out, fbyte_sds_in, fbyte_dim,
     1        lbyte_sds_out, lbyte_sds_in, lbyte_dim,
     2        sds_index_in,  sds_id_in,
     3        sds_index_out, sds_id_out

      integer data_type, num_attr, rank, rtn, dimsizes(3) 

      logical error_flag


C-------------------------------------------------------------------------------
C Initialization
C-------------------------------------------------------------------------------

      Copy_Swath_SDS =  FAIL
      error_flag     = .FALSE.


C-------------------------------------------------------------------------------
C Input argument checks 
C-------------------------------------------------------------------------------

c-----input file SDS name check
      IF (sds_name_in .EQ. BLANK) THEN
         error_flag = .TRUE.

         msgbuf = 
     1   'Input argument string "sds_name_in" is blank. '
     2   // char(10) // 'Operator Action:  Notify SDST. '

         CALL modis_smf_setdynamicmsg(MODIS_E_GENERIC, msgbuf, FUNCNAME)
      ENDIF

c-----output file SDS name check
      IF (sds_name_out .EQ. BLANK) THEN
         error_flag = .TRUE.

         msgbuf = 
     1   'Input argument string "sds_name_out" is blank. '
     2   // char(10) // 'Operator Action:  Notify SDST. '

         CALL modis_smf_setdynamicmsg(MODIS_E_GENERIC, msgbuf, FUNCNAME)
      ENDIF

c-----output file swath name check 
      IF (swathname .EQ. BLANK) THEN
         error_flag = .TRUE.

         msgbuf = 
     1   'Input argument string "swathname" is blank. '
     2   // char(10) // 'Operator Action:  Notify SDST. '

         CALL modis_smf_setdynamicmsg(MODIS_E_GENERIC, msgbuf, FUNCNAME)
      ENDIF

c-----array dimension names check 
      IF (dimlist .EQ. BLANK) THEN
         error_flag = .TRUE.

         msgbuf = 
     1   'Input argument string "dimlist" is blank. '
     2   // char(10) // 'Operator Action:  Notify SDST. '

         CALL modis_smf_setdynamicmsg(MODIS_E_GENERIC, msgbuf, FUNCNAME)
      ENDIF

c-----return if error
      IF (error_flag) RETURN


c-----determine positions of leading & trailing blanks in character variables 
      rtn = string_loc(sds_name_out, fbyte_sds_out, lbyte_sds_out)
      rtn = string_loc(sds_name_in,  fbyte_sds_in,  lbyte_sds_in)
      rtn = string_loc(dimlist,      fbyte_dim,     lbyte_dim)


C-------------------------------------------------------------------------------
C Get input file HDF array sds_id 
C-------------------------------------------------------------------------------

      sds_id_in = Get_SDS_ID( sd_id_in, sds_name_in )

      IF (sds_id_in .EQ. FAIL) THEN
         error_flag = .TRUE.

         msgbuf = 
     1   'Get_SDS_ID detected error retrieving input file array "'
     2   // sds_name_in(fbyte_sds_in:lbyte_sds_in) // '" sds_id.' 
     3   // char(10) // 'Operator Action:  Notify SDST. '

         CALL modis_smf_setdynamicmsg(MODIS_E_GENERIC, msgbuf, FUNCNAME)
      ENDIF


C-------------------------------------------------------------------------------
C Get input SDS information: name, rank, data type, dimension sizes and 
C number of attributes 
C-------------------------------------------------------------------------------

      rtn = sfginfo(sds_id_in,sds_name_in,rank,dimsizes,data_type,num_attr)

      IF (rtn .EQ. FAIL) Then
         error_flag = .TRUE.

         msgbuf = 
     1   'sfginfo returned FAIL retrieving information for input file SDS ' 
     2   // sds_name_in(fbyte_sds_in:lbyte_sds_in) // '.' 
     3   // char(10) // 'Operator Action:  Notify SDST.'

         CALL modis_smf_setdynamicmsg(MODIS_E_GENERIC, msgbuf, FUNCNAME)
      ENDIF   

      IF (error_flag) RETURN


C-------------------------------------------------------------------------------
C Define output SDS data field 
C-------------------------------------------------------------------------------

c-----data is geolocation field
      IF (sds_name_in(fbyte_sds_in:lbyte_sds_in) .EQ. 'Longitude' .OR.
     1    sds_name_in(fbyte_sds_in:lbyte_sds_in) .EQ. 'Latitude') THEN

         msg9 = 'swdefgfld'

         rtn = swdefgfld( swathid_out,
     1                    sds_name_out,
     2                    dimlist(fbyte_dim:lbyte_dim), 
     3                    data_type, 
     4                    HDFE_NOMERGE )


c-----data is science parameter field 
      ELSE
         msg9 = 'swdefdfld'

         rtn = swdefdfld( swathid_out,
     1                    sds_name_out,
     2                    dimlist(fbyte_dim:lbyte_dim), 
     3                    data_type, 
     4                    HDFE_NOMERGE )
      ENDIF


      IF (rtn .EQ. FAIL) THEN
         error_flag = .TRUE.

         msgbuf = 
     1   msg9 // ' detected error defining output SDS ' 
     2   // sds_name_out(fbyte_sds_out:lbyte_sds_out) 
     3   // char(10) // 'Operator Action:  Notify SDST.'
         CALL modis_smf_setdynamicmsg(MODIS_E_GENERIC, msgbuf, FUNCNAME)
      ENDIF

C-------------------------------------------------------------------------------
C detach and reattach to output file swath before copying input file SDS data 
C-------------------------------------------------------------------------------

      rtn = swdetach(swathid_out)

      IF (rtn .EQ. FAIL) THEN
         error_flag = .TRUE.
         
         msgbuf = 
     1   'swdetach failed to detach from output file swath "mod05".'
     2   //char(10)// 'Operator Action:  Notify SDST.' 

         CALL modis_smf_setdynamicmsg(MODIS_E_GENERIC, msgbuf, FUNCNAME)

c-----successfully detach from output file swath; now reattach
      ELSE
         swathid_out = swattach(fid_out, swathname)

         IF (swathid_out .EQ. FAIL) THEN
            error_flag = .TRUE.
         
            msgbuf = 
     1      'swattach failed to reattach to output file swath "mod05".'
     2      // char(10) // 'Operator Action:  Notify SDST.' 

            CALL modis_smf_setdynamicmsg(MODIS_E_GENERIC, msgbuf, FUNCNAME)
         ENDIF
      ENDIF

      IF (error_flag) RETURN


C-------------------------------------------------------------------------------
C Get output file HDF array sds_id 
C-------------------------------------------------------------------------------

      sds_id_out = Get_SDS_ID( sd_id_out, sds_name_out )

      IF (sds_id_out .EQ. FAIL) THEN
         error_flag = .TRUE.

         msgbuf = 
     1   'Get_SDS_ID detected error retrieving output file array "'
     2   // sds_name_out(fbyte_sds_out:lbyte_sds_out) // '" sds_id.' 
     3   // char(10) // 'Operator Action:  Notify SDST. '

         CALL modis_smf_setdynamicmsg(MODIS_E_GENERIC, msgbuf, FUNCNAME)
      ENDIF


C-------------------------------------------------------------------------------
C Copy input file SDS and all it's attributes to output file SDS  
C-------------------------------------------------------------------------------

      rtn = copy_sds( sds_id_in, sds_id_out )

      IF (rtn .eq. FAIL) THEN
         error_flag = .TRUE.

         msgbuf =
     1   'copy_sds detected error copying input file array ' 
     2   // sds_name_in(fbyte_sds_in:lbyte_sds_in) // ' to ' 
     3   // char(10) // 'output file array ' // sds_name_out(fbyte_sds_out:lbyte_sds_out) // '.'
     4   // char(10) // 'Operator Action:  Refer to prior low level LogStatus error' 
     5   // char(10) // 'messages produced by routine Copy_SDS.' 

         CALL modis_smf_setdynamicmsg(MODIS_E_GENERIC, msgbuf, FUNCNAME)
      ENDIF


C-------------------------------------------------------------------------------
C Terminate access to input and output file sds arrays 
C-------------------------------------------------------------------------------

c-----end access to input file SDS
      rtn = sfendacc(sds_id_in)

      IF (rtn .eq. FAIL) THEN
         error_flag = .TRUE.

         msgbuf = 
     1   'sfendacc detected error ending access to input file array ' 
     2   // sds_name_in(fbyte_sds_in:lbyte_sds_in)
     3   // char(10) // 'Operator Action:  Notify SDST.'

         CALL modis_smf_setdynamicmsg(MODIS_E_GENERIC, msgbuf, FUNCNAME)
      ENDIF

c-----end access to output file array
      rtn = sfendacc(sds_id_out)

      IF (rtn .eq. FAIL) THEN
         error_flag = .TRUE.

         msgbuf = 
     1   'sfendacc detected error ending access to output file array ' 
     2   // sds_name_out(fbyte_sds_out:lbyte_sds_out)
     3   // char(10) // 'Operator Action:  Notify SDST.'

         CALL modis_smf_setdynamicmsg(MODIS_E_GENERIC, msgbuf, FUNCNAME)
      ENDIF


c-----return successfully if no errors
      If (.NOT. error_flag) Copy_Swath_SDS = SUCCEED

      RETURN

      END

C**********************************************************************

      integer function Copy_SDS( sds_id_in, sds_id_out )


C-------------------------------------------------------------------------------
C!F77
C
C!Description: Copy an HDF SDS array and all its attributes from one HDF file 
C              to another.  Only two and three dimensional numeric data arrays 
C              are supported.  
C   
C!Input Parameters: 
C integer  sds_id_in        HDF sds_id of input array (copied from)
C integer  sds_id_out       HDF sds_id of output array (copied to)
C
C!Output Parameters: None
C
C!Revision History:
c Revision 1.1  1999/02/08  17:00:05  vlin
c Initial revision
C
C
C!Team-Unique Header:
C This software was developed by the MODIS Science Data Support Team
C (SDST) for the National Aeronautics and Space Administration,
C Goddard Space Flight Center, under contract NAS5-32373.
C
C HDF functions were developed at the National Center for Supercomputing
C Applications at the Univ. of Illinois at Urbana-Champaign.
C
C!References and Credits:
C   Developed in May 1997 by Fay Liang fhliang@ltpmail.gsfc.nasa.gov
C                        and Rich Hucek rhucek@ltpmail.gsfc.nasa.gov
C
C!Externals 
C   Functions:
C    sfgainfo                (dffunc.inc)
C    sfginfo                 (dffunc.inc)
C    sfrcatt                 
C    sfrdata                 (dffunc.inc)
C    sfrnatt
C    sfscatt
C    sfsnatt
C    sfwdata                 (dffunc.inc)
C    string_loc              (science code)
C
C  Named Constants:
C    FAIL                    (hdf.inc)
C    MODIS_E_GENERIC         (PGS_MODIS_39500.f)
C    SUCCEED                 (hdf.inc)
C
C!Internals: (SDST code)
C  Variables:
C  attr_buf_num       Buffer for values of attribute with numerical type
C  attr_buf_char      Buffer for values of attribute with character type
C  attr_index         Attribute index of an sds_id_in
C  buf_size           Size of a buffer
C  data_buf           Buffer the data to be read into or written to
C  data_type          Data type for data/attribute in an SDS
c
C  dimsizes           Array to store size of each dimension in an SDS
C  edges              Array containing number of data element that will
C                     be read along each dimension
C  FUNCNAME           Character string for function names
C  i                  Index in DO-loop
C  nblock             Number of block to be read/written
C  no_lines_per_buf   Number of lines could be read into data_buf
C  no_planes_per_buf  Number of planes in an SDS could be read into
C                     data_buf
C  no_reads_per_array Number of reads needed for each SDS
C  no_reads_per_plane Number of reads in a plane
C  num_of_attr        Number of attributes in an SDS
C  num_val            Number of values in an attribute of an SDS
C  plane_no           Plane number
C  plane_read_no      Plane being read
C  rank               Number of dimensions in an SDS
C  rtn                Return code from a function
C  start              Array containing position the read will start for
C                     each dimension
C  stride             Array containing interval between each read for
C                     each dimension
C  msgbuf             Character string for log messages
C  word_size          Size, in byte, of data_type
C
C  Named Constants:
C  HDFE_NOMERGE       Code for 'no merge' action
C  HDFE_AUTOMERGE     Code for 'merge' action
C  LINES_PER_READ     Number of lines to be read
C  MAX_FRAMES 
C  MAX_BYTES_PER_WORD 
C    
C
C!END-------------------------------------------------------------------

      IMPLICIT NONE

      include 'hdf.inc'
      include 'PGS_MODIS_39500.f'


C-----FORTRAN PARAMETER definitions
      character *1   BLANK
      parameter     (BLANK = ' ')

      character*(*)  FUNCNAME
      parameter     (FUNCNAME = 'Copy_SDS')

      integer        LINES_PER_READ,   MAX_FRAMES,     MAX_BYTES_PER_WORD
      parameter     (LINES_PER_READ=50, MAX_FRAMES=275, MAX_BYTES_PER_WORD=8) 

      integer        BUF_SIZE
      parameter     (BUF_SIZE = MAX_FRAMES*LINES_PER_READ*MAX_BYTES_PER_WORD)


C-----Function Argument Declarations
      integer sds_id_in, sds_id_out 


C-----Local Variable Declarations
      byte attr_buf_num(100), data_buf(buf_size), word_size

      character*25     msg25a, msg25b
      character*255    attr_name, sds_name_in
      character*1024   msgbuf
      character*100000 attr_buf_char

      integer sfginfo, sfgainfo, sfrcatt,  sfscatt, sfrnatt,  
     1        sfsnatt, sfrdata,  sfwdata,  string_loc

      integer fbyte_a, fbyte_b, fbyte_name, fbyte_aname, 
     1        lbyte_a, lbyte_b, lbyte_name, lbyte_aname

      integer attr_index, attr_data_type, num_of_attr, num_val, rank, 
     1        sds_data_type

      integer i, nblock, rtn, rtn_loc

      integer no_lines_per_buf, no_reads_per_array, no_reads_per_plane,
     1        no_planes_per_buf, plane_no, plane_read_no

      integer dimsizes(3), start(3), edges(3), stride(3)

      logical error_flag, loop_exit


C-------------------------------------------------------------------------------
C Initialization
C-------------------------------------------------------------------------------

      Copy_SDS     =  FAIL
      error_flag   = .FALSE.
      sds_name_in  =  BLANK
      loop_exit    = .false.

C-------------------------------------------------------------------------------
C Get input SDS name, rank, data type, dimsizes, and number of attributes 
C-------------------------------------------------------------------------------

c-----get SDS information
      rtn     = sfginfo(sds_id_in,sds_name_in,rank,dimsizes,sds_data_type,num_of_attr)
      rtn_loc = string_loc(sds_name_in,fbyte_name,lbyte_name)

      IF (rtn .EQ. FAIL) Then
         error_flag = .TRUE.
 
         msgbuf =
     1   'sfginfo returned FAIL retrieving input SDS information: name, rank, etc. '    
     2   // char(10) // 'Operator Action:  Notify SDST.'

         CALL modis_smf_setdynamicmsg(MODIS_E_GENERIC, msgbuf, FUNCNAME)
      ENDIF

      IF (error_flag) RETURN


c-----set variable "word_size" based on input SDS data type 
      IF (sds_data_type .EQ. DFNT_CHAR) THEN
         word_size = 1
      ELSE IF (sds_data_type .EQ. DFNT_INT8) THEN
         word_size = 1
      ELSE IF (sds_data_type .EQ. DFNT_INT16) THEN
         word_size = 2
      ELSE IF (sds_data_type .EQ. DFNT_INT32) THEN
         word_size = 4
      ELSE IF (sds_data_type .EQ. DFNT_INT64) THEN
         word_size = 8
      ELSE IF (sds_data_type .EQ. DFNT_FLOAT32) THEN
         word_size = 4
      ELSE IF (sds_data_type .EQ. DFNT_FLOAT64) THEN
         word_size = 8
      ENDIF


C-------------------------------------------------------------------------------
C Copy input SDS attributes to output SDS 
C-------------------------------------------------------------------------------

c-----get mod07 attribute data type and number of values in attribute
      DO 100 attr_index = 0, num_of_attr - 1

         IF (.NOT. error_flag) THEN
            rtn = sfgainfo(sds_id_in, attr_index, attr_name, attr_data_type, num_val)

            write(msg25a, '(I25)') attr_index 
            write(msg25b, '(I25)') attr_data_type

            rtn_loc = string_loc(msg25a,fbyte_a,lbyte_a)
            rtn_loc = string_loc(msg25b,fbyte_b,lbyte_b)

            IF (rtn .eq. FAIL) THEN
               error_flag = .TRUE.

               msgbuf = 'sfgainfo failed to read info on attribute ' // msg25a(fbyte_a:lbyte_a) 
     1         // char(10) // 'of SDS ' // sds_name_in(fbyte_name:lbyte_name) // '.' 
     2         // char(10) // 'Operator Action:  Notify SDST. ' 

               CALL modis_smf_setdynamicmsg(MODIS_E_GENERIC, msgbuf, FUNCNAME)


            ELSE
               rtn_loc = string_loc(attr_name,fbyte_aname,lbyte_aname)

c--------------data type is character
C--------------lipo testing. attr_data_type = 101        
c              attr_data_type = 101

               IF (attr_data_type .EQ. DFNT_CHAR) THEN
                  rtn = sfrcatt(sds_id_in, attr_index, attr_buf_char)

                  IF (rtn .eq. FAIL) THEN
                     error_flag = .TRUE.

		     msgbuf = 
     1               'sfrcatt failed to read attribute ' // attr_name(fbyte_aname:lbyte_aname) 
     2               // char(10) // 'on input SDS ' // sds_name_in(fbyte_name:lbyte_name) 
     3    	     // char(10) // 'Operator Action:  Notify SDST. '

                     CALL modis_smf_setdynamicmsg(MODIS_E_GENERIC, msgbuf, FUNCNAME)
                  ELSE

                     rtn = sfscatt(sds_id_out, attr_name, attr_data_type, num_val,attr_buf_char)

                     IF (rtn .eq. FAIL) THEN
                        error_flag = .TRUE.

		        msgbuf = 
     1                  'sfscatt failed to set attribute ' // attr_name(fbyte_aname:lbyte_aname)
     2                  // char(10) // 'to output SDS ' // sds_name_in(fbyte_name:lbyte_name) 
     3                  // char(10) // 'Operator Action:  Notify SDST. '

                        CALL modis_smf_setdynamicmsg(MODIS_E_GENERIC, msgbuf, FUNCNAME)
                     ENDIF  ! test sfscatt rtn

                  ENDIF   ! test sfrcatt rtn


c--------------data type is numeric
               ELSE IF ((attr_data_type .EQ. DFNT_INT8) .OR.
     +                  (attr_data_type .EQ. DFNT_INT16) .OR.
     +                  (attr_data_type .EQ. DFNT_INT32) .OR.
     +                  (attr_data_type .EQ. DFNT_INT64) .OR.
     +                  (attr_data_type .EQ. DFNT_FLOAT32) .OR.
     +                  (attr_data_type .EQ. DFNT_FLOAT64)) THEN

                  rtn = sfrnatt(sds_id_in, attr_index, attr_buf_num)

                  IF (rtn .eq. FAIL) THEN
                     error_flag = .TRUE.

                      msgbuf =
     1                'sfrnatt failed to read attribute ' // attr_name(fbyte_aname:lbyte_aname)
     2                // char(10) // 'on input SDS ' // sds_name_in(fbyte_name:lbyte_name)
     3                // char(10) // 'Operator Action:  Notify SDST. '

                      CALL modis_smf_setdynamicmsg(MODIS_E_GENERIC, msgbuf, FUNCNAME)
                  ELSE

                      rtn = sfsnatt(sds_id_out, attr_name, attr_data_type, num_val,attr_buf_num)
                      IF (rtn .eq. FAIL) THEN
                         error_flag = .TRUE.

                         msgbuf =
     1                   'sfsnatt failed to set attribute ' // attr_name(fbyte_aname:lbyte_aname)
     2                   // char(10) // 'to output SDS ' // sds_name_in(fbyte_name:lbyte_name)
     3                   // char(10) // 'Operator Action:  Notify SDST. '

                         CALL modis_smf_setdynamicmsg(MODIS_E_GENERIC, msgbuf, FUNCNAME)
                      ENDIF   !  test sfrnatt rtn

                  ENDIF  ! test sfrnatt rtn

c--------------data type not supported by routine 
               ELSE

                  error_flag = .TRUE.

                  msgbuf = 
     1            'Copy_SDS.f does not support attribute data type ('
     2            // msg25b(fbyte_b:lbyte_b) // ') returned by sfgainfo '
     3            // char(10) // 'Operator Action:  Notify SDST. '
c     1            'attribute data type (=' // msg25b(fbyte_b:lbyte_b) // ') returned by '
c     2            // 'sfgainfo not on ' // FUNCNAME // ' list.'
c     3            // char(10) // 'Operator Action:  Notify SDST. ' 

                  CALL modis_smf_setdynamicmsg(MODIS_E_GENERIC, msgbuf, FUNCNAME)

               ENDIF    ! test data type 

            ENDIF   ! test sfgainfo rtn

         ENDIF   ! test error_flag

100   CONTINUE


C-------------------------------------------------------------------------------
C populate output file (mod05) array sds_id_out.  Allow for 2D and 3D arrays only  
C-------------------------------------------------------------------------------

c-----2D case
      IF (rank .EQ. 2) THEN
         nblock = (dimsizes(2) - 1)/LINES_PER_READ + 1

         DO 200 i = 1, nblock

            start(1) = 0
            start(2) = (i - 1) * LINES_PER_READ
            stride(1) = 1
            stride(2) = 1
            edges(1) = dimsizes(1)
            edges(2) = LINES_PER_READ

            IF (i .EQ. nblock) THEN
              edges(2) = MOD((dimsizes(2)-1), LINES_PER_READ) + 1
            ENDIF

            rtn = sfrdata(sds_id_in, start, stride, edges, data_buf)
            IF (rtn .eq. FAIL) THEN
               error_flag = .TRUE.

               msgbuf =
     1         'sfrdata failed to read attribute ' // attr_name(fbyte_aname:lbyte_aname)
     2         // char(10) // 'on input SDS ' // sds_name_in(fbyte_name:lbyte_name)
     3         // char(10) // 'Operator Action:  Notify SDST. '

               CALL modis_smf_setdynamicmsg(MODIS_E_GENERIC, msgbuf, FUNCNAME)
            ENDIF   !  test sfrdata rtn

            If (rtn .ne. SUCCEED) loop_exit = .true.

            If (.not.loop_exit) Then
              rtn = sfwdata(sds_id_out, start, stride, edges, data_buf)
              IF (rtn .eq. FAIL) THEN
                 error_flag = .TRUE.

                 msgbuf =
     1           'sfwdata failed to write attribute ' // attr_name(fbyte_aname:lbyte_aname)
     2           // char(10) // 'to output SDS ' // sds_name_in(fbyte_name:lbyte_name)
     3           // char(10) // 'Operator Action:  Notify SDST. '

                 CALL modis_smf_setdynamicmsg(MODIS_E_GENERIC, msgbuf, FUNCNAME)
              ENDIF  !  test sfwdata rtn

              If (rtn .ne. SUCCEED) loop_exit = .true.
            EndIf

            If (loop_exit) return

200      CONTINUE

c-----3D case 
      ELSE IF (rank .EQ. 3) THEN
         stride(1) = 1
         stride(2) = 1
         stride(3) = 1

         no_lines_per_buf = buf_size / (dimsizes(1) * word_size)

c--------case1: whole array fits into data_buf
         IF (no_lines_per_buf .GE. (dimsizes(2) * dimsizes(3))) THEN
            no_reads_per_array = 1
            start(1) = 0
            start(2) = 0
            start(3) = 0
            edges(1) = dimsizes(1)
            edges(2) = dimsizes(2)
            edges(3) = dimsizes(3)

            rtn = sfrdata(sds_id_in, start, stride, edges, data_buf)
            IF (rtn .eq. FAIL) THEN
               error_flag = .TRUE.

               msgbuf =
     1         'sfrdata failed to read attribute ' // attr_name(fbyte_aname:lbyte_aname)
     2         // char(10) // 'on input SDS ' // sds_name_in(fbyte_name:lbyte_name)
     3         // char(10) // 'Operator Action:  Notify SDST. '

               CALL modis_smf_setdynamicmsg(MODIS_E_GENERIC, msgbuf, FUNCNAME)
            ENDIF   !  test sfrdata rtn

            If (rtn .ne. SUCCEED) loop_exit = .true.

            If (.not.loop_exit) Then
              rtn = sfwdata(sds_id_out, start, stride, edges, data_buf)
              IF (rtn .eq. FAIL) THEN
                 error_flag = .TRUE.

                 msgbuf =
     1           'sfwdata failed to write attribute ' // attr_name(fbyte_aname:lbyte_aname)
     2           // char(10) // 'to output SDS ' // sds_name_in(fbyte_name:lbyte_name)
     3           // char(10) // 'Operator Action:  Notify SDST. '

                 CALL modis_smf_setdynamicmsg(MODIS_E_GENERIC, msgbuf, FUNCNAME)
              ENDIF  !  test sfwdata rtn

              If (rtn .ne. SUCCEED) loop_exit = .true.
            EndIf

            If (loop_exit) return

c-------case2: only # of planes of array could be fit into data_buf
         ELSE IF (no_lines_per_buf .GE. dimsizes(2)) THEN
            no_planes_per_buf = INT(no_lines_per_buf/dimsizes(2))
            no_reads_per_array=INT((dimsizes(3) -1) / no_planes_per_buf)+1

            DO 300 i = 1, no_reads_per_array
               start(1) = 0
               start(2) = 0
               start(3) = (i-1) * no_planes_per_buf
               edges(1) = dimsizes(1)
               edges(2) = dimsizes(2)

               IF (i .EQ. no_reads_per_array) THEN
                  edges(3) = dimsizes(3) - start(3)
               ELSE
                  edges(3) = no_planes_per_buf
               ENDIF

               rtn = sfrdata(sds_id_in, start, stride, edges, data_buf)
               IF (rtn .eq. FAIL) THEN
                  error_flag = .TRUE.

                  msgbuf =
     1            'sfrdata failed to read attribute ' // attr_name(fbyte_aname:lbyte_aname)
     2            // char(10) // 'on input SDS ' // sds_name_in(fbyte_name:lbyte_name)
     3            // char(10) // 'Operator Action:  Notify SDST. '

                  CALL modis_smf_setdynamicmsg(MODIS_E_GENERIC, msgbuf, FUNCNAME)
               ENDIF   !  test sfrdata rtn


               If (rtn .ne. SUCCEED) loop_exit = .true.

               If (.not.loop_exit) Then
                  rtn = sfwdata(sds_id_out, start, stride, edges, data_buf)

                  IF (rtn .eq. FAIL) THEN
                     error_flag = .TRUE.

                     msgbuf =
     1               'sfwdata failed to write attribute ' // attr_name(fbyte_aname:lbyte_aname)
     2               // char(10) // 'to output SDS ' // sds_name_in(fbyte_name:lbyte_name)
     3               // char(10) // 'Operator Action:  Notify SDST. '

                     CALL modis_smf_setdynamicmsg(MODIS_E_GENERIC, msgbuf, FUNCNAME)
                  ENDIF  !  test sfwdata rtn

                  If (rtn .ne. SUCCEED) loop_exit = .true.
               EndIf

               If (loop_exit) return

300         CONTINUE


c-------case3: only # of lines in a plate of array could be fit into data_buf
         ELSE
            no_reads_per_plane=INT((dimsizes(2) -1) / no_lines_per_buf) +1
            no_reads_per_array = no_reads_per_plane * dimsizes(3)

            DO 400 i = 1, no_reads_per_array
               plane_no = (i-1) / no_reads_per_plane + 1
               plane_read_no = MOD(i -1, no_reads_per_plane) + 1

               start(1) = 0
               start(2) = (plane_read_no - 1) * no_lines_per_buf
               start(3) = plane_no - 1
               edges(1) = dimsizes(1)

               IF (plane_read_no .EQ. no_reads_per_plane) THEN
                  edges(2) = dimsizes(2) - start(2)
               ELSE
                  edges(2) = no_lines_per_buf
               ENDIF

               edges(3) = 1

               rtn = sfrdata(sds_id_in, start, stride, edges, data_buf)
               IF (rtn .eq. FAIL) THEN
                  error_flag = .TRUE.

                  msgbuf =
     1            'sfrdata failed to read attribute ' // attr_name(fbyte_aname:lbyte_aname)
     2            // char(10) // 'on input SDS ' // sds_name_in(fbyte_name:lbyte_name)
     3            // char(10) // 'Operator Action:  Notify SDST. '

                  CALL modis_smf_setdynamicmsg(MODIS_E_GENERIC, msgbuf, FUNCNAME)
               ENDIF   !  test sfrdata rtn


               If (rtn .ne. SUCCEED) loop_exit = .true.

               If (.not.loop_exit) Then
                 rtn = sfwdata(sds_id_out, start, stride, edges, data_buf)

                 IF (rtn .eq. FAIL) THEN
                    error_flag = .TRUE.

                    msgbuf =
     1              'sfwdata failed to write attribute ' // attr_name(fbyte_aname:lbyte_aname)
     2              // char(10) // 'to output SDS ' // sds_name_in(fbyte_name:lbyte_name)
     3              // char(10) // 'Operator Action:  Notify SDST. '

                    CALL modis_smf_setdynamicmsg(MODIS_E_GENERIC, msgbuf, FUNCNAME)
                 ENDIF  !  test sfwdata rtn

                 If (rtn .ne. SUCCEED) loop_exit = .true.
               EndIf

               If (loop_exit) return

400         CONTINUE

         ENDIF  !end if (no_lines_per_buf .GE. (dimsizes(2) * dimsizes(3)))
      ENDIF !end if (rank .EQ. 2)


      If (.NOT. error_flag) Copy_SDS = SUCCEED

      RETURN

      END

C***********************************************************************

      integer function Get_SDS_ID( sd_id, sds_name ) 

C----------------------------------------------------------------------
C!F77
C
C!Description: Return HDF SDS ID given file ID and SDS Name. 
C   
C!Input Parameters: 
C
C integer       sd_id        HDF SD interface file ID
C character*(*) sds_name     HDF SDS name 
C
C!Output Parameters: 
C  
C integer   Get_SDS_ID       HDF SDS ID 
C
C!Revision History:
c Revision 1.1  1999/02/08  17:00:05  vlin
c Initial revision
c
C
C!Team-Unique Header:
C This software was developed by the MODIS Science Data Support Team
C (SDST) for the National Aeronautics and Space Administration,
C Goddard Space Flight Center, under contract NAS5-32373.
C
C HDF functions were developed at the National Center for Supercomputing
C Applications at the Univ. of Illinois at Urbana-Champaign.
C
C!References and Credits:
C   Developed in December, 1998 by Vicky Lin,  vlin@ltpmail.gsfc.nasa.gov
C   under guidance of Rich Hucek, rhucek@ltpmail.gsfc.nasa.gov
C
C!Externals 
C   Functions:
C    sfn2index               (dffunc.inc)
C    sfselect                (dffunc.inc)
C
C  Named Constants:
C    FAIL                    (hdf.inc)
C    MODIS_E_GENERIC         (PGS_MODIS_39500.f)
C    SUCCEED                 (hdf.inc)
C
C!Internals: (SDST code)
C  Variables:
C  data_type          Data type for data/attribute in an SDS
C  dimsizes           Array to store size of each dimension in an SDS
C                     be read along each dimension
C  fbyte          1st byte of nonblank text of sds_name
C  lbyte          Last byte of nonblank text of sds_name
C  rtn                Return code from a function
C                     each dimension
C  sds_id_out         HDF SDS ID in output file 
C  sds_index_out      SDS index in output file 
C  sds_index_in       SDS index in input file
C  msgbuf             Character string for log messages
C
C  Named Constants:
C  FUNCNAME           Character string for function names
C    
C!END-------------------------------------------------------------------

      IMPLICIT NONE

      include 'hdf.inc'
      include 'PGS_MODIS_39500.f'

C-----FORTRAN PARAMETER definitions
      character*(*)  BLANK
      parameter     (BLANK = ' ')

      character*(*)  FUNCNAME
      parameter     (FUNCNAME = 'Get_SDS_ID')


C-----Function Argument Declarations
      character*(*) sds_name
      integer sd_id


C-----Local Variable Declarations
      character*25   msg25, msg25a
      character*2048 msgbuf

      integer fbyte,     lbyte, 
     1        fbytea,    lbytea,
     2        sds_index, sds_id,
     3        sfn2index, sfselect, string_loc,
     4        rtn_loc

      logical error_flag


C-------------------------------------------------------------------------------
C Initialization
C-------------------------------------------------------------------------------

      Get_SDS_ID =  FAIL
      error_flag = .FALSE.
      rtn_loc    =  string_loc(sds_name,fbyte,lbyte)

c-----lipo 03/03/99 set sd_id negative for testing
c      sd_id = -2

      write(msg25a,'(I25)') sd_id
      rtn_loc    =  string_loc(msg25a,fbytea,lbytea)


C-------------------------------------------------------------------------------
C Check for non-positive sd_id and BLANK sds_name
C-------------------------------------------------------------------------------

      IF (sds_name .EQ. BLANK) THEN
         error_flag = .TRUE.

         msgbuf = 
     1   'Input argument "sds_name" is blank. '
     2   // char(10) // 'Operator Action:  Notify SDST. '

         CALL modis_smf_setdynamicmsg(MODIS_E_GENERIC, msgbuf, FUNCNAME)
      ENDIF

      IF (sd_id .LE. 0) THEN
         error_flag = .TRUE.

         msgbuf = 
     1   'Input argument "sd_id" (' // msg25a(fbytea:lbytea) // ') non-positive. '
     2   // char(10) // 'Operator Action:  Notify SDST. '

         CALL modis_smf_setdynamicmsg(MODIS_E_GENERIC, msgbuf, FUNCNAME)
      ENDIF

      IF (error_flag) RETURN


C-------------------------------------------------------------------------------
C Retrieve SDS ID
C-------------------------------------------------------------------------------

c-----first get HDF SDS Index 
      sds_index = sfn2index(sd_id, sds_name)

      IF (sds_index .EQ. FAIL) THEN
         error_flag = .TRUE.

         msgbuf = 
     1   'sfn2index unable to get sds_index of array ' // sds_name(fbyte:lbyte) 
     3   // char(10) // 'Operator Action:  Notify SDST.' 

        CALL modis_smf_setdynamicmsg(MODIS_E_GENERIC, msgbuf, FUNCNAME)

      ELSE

c--------get HDF SDS ID 
         sds_id = sfselect(sd_id,sds_index)

         IF (sds_id .EQ. FAIL) THEN
            error_flag = .TRUE.

            msgbuf = 
     1      'sfselect unable to get sds_id of array ' // sds_name(fbyte:lbyte) 
     3      // char(10) // 'Operator Action:  Notify SDST.' 

            CALL modis_smf_setdynamicmsg(MODIS_E_GENERIC, msgbuf, FUNCNAME)
         ENDIF

      ENDIF

c-----override default return value
      If (.NOT. error_flag) Get_SDS_ID = sds_id

      RETURN

      END
