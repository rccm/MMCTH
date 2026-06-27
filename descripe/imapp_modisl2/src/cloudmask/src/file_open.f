      subroutine file_open(qa_bits,h_eco1,h_eco2,modfil_L1B_250,
     +                     cubefile_1,modfil_L1B_1km,cubefile_3,
     +                     modfil_Geo,geofile,outfile,no_250)
      implicit none
      save

      include 'mapi.inc'
      include 'mod35.inc'

      include 'PGS_PC.f'
      include 'PGS_IO.f'
      include 'PGS_SMF.f'

c     Common block used for writing debug output statements
      common / bug / debug, h_output

C-----------------------------------------------------------------------
C!F77
C
C!DESCRIPTION:
c Program which will open files and return unit numbers and
c hdf file id's for all files used in the main program.  It
c uses mapi and sdp toolkit functions as the open calls.
c This file contains a description of where the L3 snow mask
c and clear radiance file would be accessed.
C
C
C!INPUT PARAMETERS:    None.
C
C!OUTPUT PARAMETERS:
C
c   qa_bits           Byte array containing qa bit results
C   integer h_*       File handles used to manipulate files opened
C                     for READ/WRITE
c   cubefile_1        L1B 250m file name
c   cubefile_3        L1B 1km file name
c   geofile           Geolocation product name
c   outfile           MOD35 product file name
C   integer modfil_*  HDF file id's
c   no_250            Logical flag set to true if 250 data is used
C
C!REVISION HISTORY:
c
C!TEAM-UNIQUE HEADER:
C
C!REFERENCES AND CREDITS
C
C!DESIGN NOTES:
C
C       TO USE THIS MODULE, THE ENVIRONMENT VARIABLE PGS_PC_INFO_FILE
C       MUST BE SET TO YOUR PROCESS CONTROL FILE (PCF).  IN ADDITION,
C       A LINK TO THE SDP TOOLKIT INCLUDE FILES MUST BE PRESENT IN THE
C       DIRECTORY USED FOR EXECUTION.
C
C EXTERNALS:
C
C       NAMED CONSTANTS:
C               PGSd_IO_Gen_RDirUnf             (PGS_IO.f)
C               PGSd_IO_Gen_WSeqFrm             (PGS_IO.f)
C               PGSd_IO_Gen_WDirUnf             (PGS_IO.f)
C               PGS_S_SUCCESS                   (PGS_SMF.f)
C
C       FUNCTIONS AND SUBROUTINES:
C               pgs_pc_getreference             (libPGSTK.a)
c               pgs_io_gen_open                 (libPGSTK.a)
c               opmfil                          (mapi.a)
c               message.f
c               set_qa_bit.f
C
C!END-------------------------------------------------------------------

C     Variables
      character*(PGSd_PC_FILE_PATH_MAX) geofile,
     *           outfile, cubefile_1,cubefile_3
c rhucek 4/16/02: ddded character*80 declarations 
      character*80  hdfattrname, parmname, parmvalue
      integer modfil_Geo(MODFILLEN),
     *        modfil_L1B_250(MODFILLEN),modfil_L1B_1km(MODFILLEN),
     *        file_version,file_access,record_length,
     *        h_eco1,h_eco2,h_output,
     *        rtn,ier,debug,file_id
      logical open_flag, no_250

      byte qa_bits(10)

c ... set program name for error messaging
      character*32 routine_name
      parameter ( routine_name = 'file_open' )

C     Functions
c rhucek 04/16/02:  added function pgs_met_getpcattr_s
      integer pgs_pc_getreference,pgs_io_gen_openf,pgs_met_getpcattr_s

      integer num_args
      integer FlagRA
      character FlagBuff*10
      num_args = command_argument_count()
      
      if(num_args == 1) then
         call get_command_argument(1,FlagBuff)
         read (FlagBuff,*) FlagRA
      else
      !This is the default value
         FlagRA = 0
      endif

C     Initialize file handles
      h_eco1 = -5555
      h_eco2 = -5555
      h_output= -5555

c ... Initialize flag that will tell us if 250 file is to be used
      no_250 = .false.

c ... initialize successful open flag
      open_flag = .false.

C ... First, open the generic files
c ... Open ecosystem files

c ... First open the Global 1 km ecosystem file
      file_version = 1
      file_access = PGSd_IO_Gen_RDirUnf
      record_length = 512

      ier = pgs_io_gen_openf(LRN_eco1, file_access,
     &    record_length, h_eco1, file_version)
      if( ier .eq. PGS_S_SUCCESS ) then
        call message( routine_name,
     &  'Successfully opened global 1km ecosystem file.',
     &  0, 3 )
        open_flag = .false.
        inquire(h_eco1,opened=open_flag)
        if (.not. open_flag) then
          call message( routine_name,
     &    'Could not open open global 1km ecosystem file. ' //
     &    '  Will attempt to use the 10 minute global file. ' //
     &    ' [OPERATOR ACTION: Check input file and PCF file]',
     &    0, 3 )
          h_eco1 = -5555
         endif
      else
        call message( routine_name,
     &  'Could not open global 1km ecosystem file.' //
     &  '  Will attempt to use the 10 minute global file. ' //
     &  ' [OPERATOR ACTION: Check input file and PCF file]',
     &  0, 3 )
        h_eco1 = -5555
      endif

c ... Open global 10 minute file if 1 km file is not available
      if (h_eco1 .eq. -5555) then
        file_version = 1
        file_access = PGSd_IO_Gen_RDirUnf
        record_length = 2

        ier = pgs_io_gen_openf(LRN_eco2, file_access,
     &      record_length, h_eco2, file_version)
        if( ier .eq. PGS_S_SUCCESS ) then
          call message( routine_name,
     &   'Successfully opened Global 10 minute ecosystem file. ',
     &    0, 3 )
        else
          call message( routine_name,
     &    'Could not open global 10 minute ecosystem file.' //
     &    ' [OPERATOR ACTION: Check input file and PCF file]',
     &    0, 2 )
        endif
        open_flag = .false.
        inquire(h_eco2,opened=open_flag)
        if (.not. open_flag) then
          call message( routine_name,
     &    'Could not open global 10 minute ecosystem file.' //
     &    ' [OPERATOR ACTION: Check input file and PCF file]',
     &    0, 2 )
         endif
      endif


c ... Open runtime output debug file
      file_version = 1
      file_access = PGSd_IO_Gen_WSeqFrm
      record_length = 1

      ier = pgs_io_gen_openf(LRN_out, file_access,
     &    record_length, h_output, file_version)
      if( ier .eq. PGS_S_SUCCESS ) then
        call message( routine_name,
     &  'Successfully opened output debug file. ',
     &  0, 3 )
      else
        call message( routine_name,
     &  'Could not open output debug file.' //
     &  ' [OPERATOR ACTION: Remove existing debug file]', 0, 2 )
      endif
      open_flag = .false.
      inquire(h_output,opened=open_flag)
      if ( .not. open_flag) then
        call message( routine_name,
     &  'Could not open output debug file. ' //
     &  ' [OPERATOR ACTION: Check input file and PCF file]',
     &  0, 2 )
       endif


c ... Now for the hdf files

c ... Get geolocation file name from process control file

      file_id = LRN_Geo
      file_version = 1
      geofile = ' '
      ier = pgs_pc_getreference( file_id, file_version, geofile )
      if( ier .eq. PGS_S_SUCCESS ) then
        call message( routine_name,
     &  'Success getting geoloc. filename from process control file',
     &  0, 3 )
      else
        call message( routine_name,
     &  'Error getting MODIS geolocation filename from PCF' //
     &  ' [OPERATOR ACTION: Check PCF entry]', 0, 2 )
      endif

c ... Open the geolocation file

      rtn = opmfil( geofile, 'r',  modfil_Geo)
      if( rtn .eq. 0 ) then
        call message( routine_name,
     &  'Success opening geolocation data file ', 0, 3 )
      else
        call message( routine_name,
     &  'Error opening MODIS geolocation input file' //
     &  ' [OPERATOR ACTION: Check that file exists]', 0, 2 )

      endif

c ... Retrieve DayNightFlag from geolocation product or exit if unsuccessful 
      file_id      = LRN_Geo
      file_version = 1  !prior call to pgs_pc_getreference returned file_version=0;
      hdfattrname  = 'CoreMetadata.0'
      parmname     = 'DAYNIGHTFLAG'
      parmvalue    = ' '

      ier = pgs_met_getpcattr_s( file_id, file_version, 
     &                           hdfattrname, parmname, parmvalue)

      if ( ier .ne. PGS_S_SUCCESS ) then
         call message( routine_name,  
     &   'Unable to retrieve DayNightFlag from MODIS geolocation file.'
     &   // char(10) // ' [OPERATOR ACTION: Check PCF entry]', 0, 2 ) 
      endif


c ... Get Level 1B 250m data if swath not sampled entirely in MODIS 
c ... instrument "Night" mode 

c ... MODIS instrument at least partially in "Day" mode  
      if ( (parmvalue .eq. 'Day') .or. (parmvalue .eq. 'Both') ) then
         file_id = LRN_L1B_250
         file_version = 1
         cubefile_1 = ' '

c ...... Get pathname of L1B 250m product
         ier = pgs_pc_getreference( file_id, file_version, cubefile_1 )

         if( ier .ne. PGS_S_SUCCESS) then
            no_250 = .true.

            call message( routine_name,
     &      'Error getting MODIS L1B 250m filename from PCF: ' //
     &      'L1B 250m data will not be used. ' // char(10) //
     &      '[OPERATOR ACTION: Check PCF entry.  If error ' // 
     &      'persists, contact SDST]', 0, 0 )
         else
            call message( routine_name,
     &      'Success getting cube L1B 250m filename from process ' //
     &      'control file', 0, 3 )

c ........  Open the L1B 250m product file
            rtn = opmfil( cubefile_1, 'r',  modfil_L1B_250)

            if ( rtn .eq. 0 ) then
               call message( routine_name,
     &         'Success opening L1B 250m scan cube data file ', 0, 3 )
            else
               no_250 = .true.

               call message( routine_name,
     &         'Error opening MODIS L1B 250m input file: L1B 250m ' //
     &         'data will not be used. ' // char(10) // 
     &         '[OPERATOR ACTION: In future, stage correct 250m L1B ' //
     &         ' granule.' // ' If error persists, contact SDST]', 0, 0 )
            endif ! open L1B 250m product

         endif ! retrieve L1B 250m product name from PCF 

      else 
         no_250 = .true.

         call message( routine_name,
     &   'MODIS in "Night" mode: L1B 250m data not used.', 0, 3 )
      endif ! check for illuminated swath 


c ... Get Level 1B 1km scan cube file name from process control file

      if( FlagRA == 1) then
         file_id = LRN_DSL1B_1km_RA
      else
         file_id = LRN_DSL1B_1km
      endif
      file_version = 1
      cubefile_3 = ' '
      ier = pgs_pc_getreference( file_id, file_version, cubefile_3 )
      if( ier .eq. PGS_S_SUCCESS) then
        call message( routine_name,
     &  'Success getting cube 1km filename from process control file',
     &  0, 3 )
      else
        call message( routine_name,
     &  'Error getting MODIS L1B 1km filename from PCF' //
     &  ' [OPERATOR ACTION: Check PCF entry]', 0, 2 )
      endif

c ... Open the L1B file

      rtn = opmfil( cubefile_3, 'r',  modfil_L1B_1km)
      if( rtn .eq. 0 ) then
        call message( routine_name,
     &  'Success opening 1km scan cube data file ', 0, 3 )
      else
        call message( routine_name,
     &  'Error opening MODIS L1B 1km input file' //
     &  ' [OPERATOR ACTION: Check that file exists]', 0, 2 )
      endif


c ... get output file name from process control file
      file_id = LRN_mod35
      file_version = 1
      outfile = ' '
      ier = pgs_pc_getreference( file_id, file_version, outfile )

      if( ier .eq. PGS_S_SUCCESS ) then
        call message( routine_name,
     &  'Success getting output filename from process control file',
     &  0, 3 )
      else
        call message( routine_name,
     &  'Error getting MODIS output filename from PCF' //
     &  ' [OPERATOR ACTION: Check PCF entry]', 0, 2 )
      endif

c ..................................................................
c     Add qa_bit initializations for files expected to be used
c     in future cloud mask versions
c     Clear radiance files
      call set_qa_bit(qa_bits,56)
      call set_qa_bit(qa_bits,57)
c     Digital Elevation Model
      call set_qa_bit(qa_bits,72)
c .................................................................

c .....debug statement.............................................
c ... write output file names to debug file if desired
      write(h_output,'(15x,''BEGIN RUNTIME OUTPUT DEBUG '',/)')
      write(h_output,'(15x,''  Files Opened '',/)')
      write(h_output,'(2x,'' 1km  L1b file: '',A70/)') cubefile_3
      write(h_output,'(2x,'' 250m L1b file: '',A70/)') cubefile_1
      write(h_output,'(2x,'' Geolocation file: '',A70/)') geofile
c .................................................................

      return
      end
