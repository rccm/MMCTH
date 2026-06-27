      PROGRAM MOD_PR07

c-----------------------------------------------------------------------
c !F77
c
c !DESCRIPTION:
c     Driver program for UW MODIS Atmospheric Profiles.
c
c !INPUT PARAMETERS:
c     See the accompanying README
c
c !OUTPUT PARAMETERS:
c     See the accompanying README
c
c !REVISION HISTORY:
c EvaB, October 23, 2009 QA initialization was missing (respons for Paul Hubanks's note 
c Kevin Bagett solved it.      

c !TEAM-UNIQUE HEADER:
c     Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
c
c !END
c-----------------------------------------------------------------------

      IMPLICIT NONE

c ... PGS toolkit include files

      include 'PGS_PC.f'
      include 'mapi.inc'

c ... UW include files

      include 'mod07_pcfnum.inc'
      include 'mod07_data.inc'
      include 'mod07_meta.inc'
      include 'debug.inc'
      include 'platform_name.inc'
      
c rhucek 12/17/02:  added SDST include file 
c ... SDST include files

      include 'Atmos_AncData.inc'

c ... Local variables

      character*128 cube_file, geo_file, cld_file, out_file, bug_file
      integer cube_handle( MODFILLEN ), geo_handle( MODFILLEN ),
     &  cld_handle( MODFILLEN ), out_handle( MODFILLEN )

      double precision east_lim, west_lim, north_lim, south_lim
      character*70 begin_date, begin_time, end_date, end_time

      integer nscans, npixels, scan, line, pixel, ngood, nmissing,
     &  ncloudy, nland, flag, good_total, month
      real landfrac

      character*2 out_number
      integer out_lines_1km, out_elements_1km, out_lines_5km,
     &  out_elements_5km

      integer start( 2 ), stride( 2 ), edge( 2 ), rtn
      character*100 text

c rhucek 12/17/02:  define sdptk temp file luns 
c rhucek 08/06/04:  removed SEA_ICE from temp file lun list 
      integer NumTempFiles /1/, 
     &        TempFile_LUNList(1) /LUN_TEMP_GDAS_0ZF/

      character*(PGSd_MET_GROUP_NAME_L)
     &  MET_Groups( PGSd_MET_NUM_OF_GROUPS )
      real pcn4s, pcnnd, pcnnt, pcnng, pcnni, pcnnl, pcbad,
     &  pcnnw, pcnns, pcnnv, pcnnr, pcnnc, max_sol, min_sol
c gfireman 12/27/04: declared FileNameSuffix 
      character*3 FileNameSuffix


      integer DFACC_WRITE, DFNT_CHAR8
      parameter   (DFACC_WRITE = 2, 
     +             DFNT_CHAR8  = 4)

      integer sfstart, sfscatt 
      integer sfend, pgs_pc_getconfigdata
      
      integer sd_id
      character*4 pcf_satid
      character*255 doi_char
      integer doiLen
      
      integer LUN_Sat_Instrument
      parameter (LUN_Sat_Instrument = 800510)


c ... Declarations for function 'get_platform_name'

      integer fileid, version, status
      integer get_platform_name
      external get_platform_name

c ... External functions

      integer copy_info, pgs_io_gen_closef
      external clmfil, copy_info, pgs_io_gen_closef

      integer delete_commontempfiles_atmos
      external delete_commontempfiles_atmos

      integer Write_ECSMetadata_nonHDF
      external Write_ECSMetadata_nonHDF

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

C-----------------------------------------------------------------------
C     INITIALIZATION
C-----------------------------------------------------------------------

c ... Set MODIS IR band numbers used by MOD07
c ... (Can't use a DATA statement because BANDS is in a COMMON block)

      bands( 1) = 24
      bands( 2) = 25
      bands( 3) = 27
      bands( 4) = 28
      bands( 5) = 29
      bands( 6) = 30
      bands( 7) = 31
      bands( 8) = 32
      bands( 9) = 33
      bands(10) = 34
      bands(11) = 35
      bands(12) = 36

c ... Get platform name from MOD021KM file ('Terra' or 'Aqua')

      if( FlagRA == 1) then
         fileid = cube_pcfnum_RA
      else
         fileid = cube_pcfnum
      endif
      version = 1
      rtn = get_platform_name(fileid, version, status, platform_name)

c ... Check return value from get_platform_name

      if (rtn .ne. 0) then
        call message( 'MOD_PR07', 'Error getting platform name' //
     &    ' [OPERATOR ACTION: Contact SDST]', 0, 2)
      endif

c ... Check the platform name (it should be 'Terra' or 'Aqua')

      if (platform_name(1:5) .ne. 'Terra' .and.
     &    platform_name(1:4) .ne. 'Aqua') then
        call message( 'MOD_PR07', 'Platform name not recognized' //
     &    ' [OPERATOR ACTION: Contact SDST]', 0, 2)
      endif

c ... Get debug value from PCF

      call get_debug( debug )

c ... Get input and output filenames,
c ... Check file times against PCF collection time

      call mod07_open_files( cube_file, geo_file, cld_file, out_file,
     &  bug_file, cube_handle, geo_handle, cld_handle, h_output )

c ... Read ECS Core Metadata

      call get_core_metadata( geo_handle, geo_pcfnum, mcf_pcfnum,
     &  east_lim, west_lim, north_lim, south_lim,
     &  begin_date, begin_time, end_date, end_time )

c ... Get the month (begin_date format is YYYY-MM-DD)

      read( begin_date, '(5x, i2)' ) month

c ... Get number of scans and pixels for 1km bands

      call get_global_metadata( cube_handle, nscans, npixels )

c ... Create the MOD07 output HDF file, and open it using MAPI

      out_number = '07'
      out_lines_1km = nscans * 10
      out_elements_1km = npixels
      out_lines_5km = out_lines_1km / isamp
      out_elements_5km = out_elements_1km / isamp

      call setup_create_hdf( out_number, out_file, out_lines_1km,
     &  out_elements_1km, out_lines_5km, out_elements_5km, out_handle )

      if ( debug .eq. 1 ) then
        write( h_output, '(/,72(''-''),/)' )
        write( h_output, '(''setup_create_hdf debug info'',/)' )
        write( h_output, '(''Output file was created successfully'')' )
      endif

C-----------------------------------------------------------------------
C     PROCESSING
C-----------------------------------------------------------------------

      good_total = 0

c ... Begin scan loop

      do scan = 1, nscans

c ...   Write debug info

        if ( debug .eq. 1 ) then
          write( h_output, '(/,72(''-''),/)' )
          write( h_output, '(''Processing scan '',i4)' ) scan
        endif

c ...   Extract radiance, geolocation, mask, and ancillary data
        call mod07_get_data( cube_handle, cube_file, geo_handle,
     &    cld_handle, scan, begin_date, begin_time )

c ...   Initialize output product arrays for this scan
        call mod07_initialize_output()

c ...   Start loop over pixel boxes (line and pixel is at corner of box)
        do line = 1, ( max_line / isamp ) * isamp, isamp

          do pixel = 1, ( npixels / isamp ) * isamp, isamp

c ...     Check this box for good data
          call mod07_check_box( line, pixel, ngood, nmissing, ncloudy,
     &      nland )

c EvaB, October 23, 2009 QA initialization was missing (respons for Paul Hubanks's note 
c Kevin Bagett solved it.      
          flag=-1    

c ...     If this box has enough good pixels, compute products
          if ( ngood .ge. min_ngood ) then
    
c ...       Compute land fraction
            landfrac = real(nland) / real(ngood)

c ...       Compute retrieval 
            call mod07_compute_products( line, pixel, month, landfrac,
     &        flag )

c ...       Increment retrieval counter            
            if ( flag .eq. 0 ) good_total = good_total + 1
            
          endif

c ...     Set the QA information for this box
          call mod07_set_qa( line, pixel, ngood, nmissing, ncloudy,
     &      flag )

          end do

        end do

c ...   Write products for this scan

        call mod07_write_products( out_handle, scan, npixels / isamp )

c ... End of loop over all scans

      end do

C-----------------------------------------------------------------------
C     TERMINATION
C-----------------------------------------------------------------------

c ... Get metadata from the cloud mask input file

      call copy_cldmask_meta( cld_pcfnum, pcn4s, pcnnd,
     &  pcnnt, pcnng, pcnni, pcnnl, pcnnw, pcnns,
     &  pcnnv, pcnnr, pcnnc, max_sol, min_sol )

c ... Get metadata for MOD07

      pcbad = 100.0 * ( 1.0 - real( good_total ) /
     &  real( out_lines_5km * out_elements_5km ) )
      call mod07_prepare_ecs_metadata( pcbad, pcn4s, pcnnd, pcnnt,
     &  pcnng, pcnni, pcnnl, pcnnw, pcnns, pcnnv, pcnnr, pcnnc,
     &  max_sol, min_sol, met_groups )

c ... Set the runtime debug metadata

c gfireman 12/27/04: changed to call Write_ECSMetadata_nonHDF
      FileNameSuffix = 'txt'
      rtn = Write_ECSMetadata_nonHDF( mcf_qcnum,
     1                                bugRPref_pcfnum,
     2                                LRN_OF_NUM_INV_RP_PAIRS,
     3                                0,
     4                                FileNameSuffix )
      if ( rtn .ne. 0 ) then
         call message( 'MOD_PR07',
     &   'Error writing runtime debug metadata' //
     &   ' [OPERATOR ACTION: Notify SDST]', 0, 2 )
      endif

c ... Close input files, and close/complete the output file

      call mod07_close_files( cube_handle, geo_handle, cld_handle,
     &  out_handle, met_groups )

c ... Copy geolocation information to output file

      start( 1 ) = isamp / 2
      start( 2 ) = isamp / 2
      stride( 1 ) = isamp
      stride( 2 ) = isamp
      edge( 1 ) = out_elements_5km
      edge( 2 ) = out_lines_5km
      text = ' '

      rtn = -1
      rtn = copy_info( geo_file, out_file, start, stride, edge, text )
      if ( debug .eq. 1 ) then
        write( h_output, '(/,72(''-''),/)' )
        write( h_output, '(''copy_info debug info'',/)' )
        if ( rtn .ne. 0 ) then
          write( h_output, '(''Geolocation copy failed:'',/,a)' ) text
        else
          write( h_output, '(''Geolocation copy successful'')' )
        endif
      endif
      if ( rtn .ne. 0 ) call message( 'MOD_PR07',
     &  'Geolocation copy failed:' // text //
     &  ' [OPERATOR ACTION: Check available disk space.' //
     &  ' If problem persists, contact SDST]', 0, 2 )

c ... Write debug info and then close debug file

      if ( debug .eq. 1 ) then
        write( h_output, '(/,72(''-''),/)' )
        write( h_output, '(''MOD_PR07 terminated normally'',/)' )
      endif

      rtn = -1
      rtn = pgs_io_gen_closef( h_output )
      if ( rtn .ne. 0 ) call message( 'MOD_PR07',
     &  'Close failed on runtime debug file ' // bug_file //
     &  ' [OPERATOR ACTION: Check available disk space.' //
     &  ' If problem persists, contact SDST]', 0, 2 )

c ... Delete temporary file entries in PCF

c rhucek 12/12/02: referencing new version of temp file cleanup function
      rtn = -1
      rtn = Delete_CommonTempFiles_Atmos(TempFile_LUNList, NumTempFiles)

      if ( rtn .ne. 0 ) call message( 'MOD_PR07',
     &  'Failed to delete temporary files' // bug_file //
     &  ' [OPERATOR ACTION: Contact SDST]', 0, 1 )


c           Get satellite instrument name.
      rtn = pgs_pc_getconfigdata(LUN_Sat_Instrument,pcf_satid)
      if (rtn .ne. 0) then
        call message( 'MOD_PR07',
     &  'Error reading instrument name from pcf LUN 800510.' //
     &  ' [OPERATOR ACTION: Check PCF file, fix if corrupt.]',
     &  0, 1 )
      endif

      if (pcf_satid .eq. 'AM1M') then
         doi_char = '10.5067/MODIS/MOD07_L2.006'
         doiLen = 26
      else
         doi_char = '10.5067/MODIS/MYD07_L2.006'
         doiLen = 26
      endif

      sd_id = sfstart(out_file, DFACC_WRITE)
      if(sd_id .eq. -1) then
         call message( 'MOD_PR07',
     &        'Problem openning the file', 0, 2 )
      endif

      rtn = sfscatt(sd_id, 'identifier_product_doi', DFNT_CHAR8, doiLen,  
     +                 doi_char) 
      if (rtn .eq. -1) then
         call message( 'MOD_PR07',
     &        'Problem writting the global attribute identifier_product_doi', 0, 2 )
      endif

      rtn = sfscatt(sd_id, 'identifier_product_doi_authority', DFNT_CHAR8, 17,  
     +                 'http://dx.doi.org') 
      if (rtn .eq. -1) then
         call message( 'MOD_PR07',
     &        'Problem writting the global attribute identifier_product_doi_authority', 0, 2 )
      endif

      rtn = sfend(sd_id)
      if (rtn .eq. -1) then
         call message( 'MOD_PR07',
     &        'Problem closing the file', 0, 2 )
      endif


c ... Return exit code 0 to shell

      call exit( 0 )

      END
