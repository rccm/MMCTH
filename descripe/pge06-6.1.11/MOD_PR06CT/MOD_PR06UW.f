      SUBROUTINE MOD_PR06UW( SUCCESS_CO2, SUCCESS_IRP, PLOW, PMID,
     &                       PHIGH, PTHIN, PTHICK, POPQ, PCIRRUS,
     &                       PICE, PWATER, PMIXED, PUNC, POCEAN,
     &                       PLAND, PSNOW, ERROR_FLAG )

c-----------------------------------------------------------------------
c!F77
c
c!Description:
c    This program is the main subroutine which will return all
c    the metadata information and an error flag back to the
c    calling program.
c
c!Input Parameters:
c
c None.
c
c!Output Parameters:
c success_co2   Percentage of pixels with successful CO2 retrievals
c success_irp   Percentage of pixels with successful irp retrievals
c plow          Percentage of pixels found to have low cloud
c pmid          Percentage of pixels found to have mid cloud
c phigh         Percentage of pixels found to have high cloud
c pthin         Percentage of pixels found to have thin cloud
c pthick        Percentage of pixels found to have thick cloud
c popq          Percentage of pixels found to have opaque cloud
c pcirrus       Percentage of pixels found to have cirrus cloud
c pice          Percentage of pixels found to have ice cloud
c pwater        Percentage of pixels found to have water cloud
c pmixed        Percentage of pixels found to have mixed phase cloud
c punc          Percentage of pixels found to have uncertain phase cloud
c pocean        Percentage of pixels with ocean background
c pland         Percentage of pixels with land background
c psnow         Percentage of pixels with snow background
c error_flag    Flag to show that subroutine completed nornally
c               or not (0 = normal, 1 = abnormal)
c
c!Revision History:
c
c!Team-unique Header:
c
c!End
c-----------------------------------------------------------------------

C#######################################################################
C     DECLARATIONS
C#######################################################################

      IMPLICIT NONE

c ... PGS toolkit include files

      include 'PGS_SMF.f'
      include 'PGS_PC.f'
      include 'PGS_IO.f'
      include 'mapi.inc'

c ... UW include files

      include 'mod06uw_data.inc'
      include 'mod06uw_debug.inc'
      include 'mod06uw_pcfnum.inc'
      include 'platform_name.inc'
      
c ... Scalar arguments

      real success_co2, success_irp, plow, pmid, phigh,
     &     pthin, pthick, popq, pcirrus, pice, pwater,
     &     pmixed, punc, pocean, pland, psnow
      integer error_flag

c ... Local arrays

      character*10 routine_name
      parameter ( routine_name = 'MOD_PR06UW' )
      character*(PGSd_PC_FILE_PATH_MAX)
     &  cube_file, geo_file, cld_file, out_file, bug_file
      integer cube_handle( MODFILLEN ), geo_handle( MODFILLEN ),
     &        cld_handle( MODFILLEN ), out_handle( MODFILLEN ),
     &        indat_bt( max_samp_line*max_samp_pixel , irphase_band+5 ),
     &        icode( max_samp_line*max_samp_pixel ), met_date(4)
      real var ( max_samp_line*max_samp_pixel , irphase_band )
      byte qual_flag_ir ( max_samp_line*max_samp_pixel ),
     &     conf_flag_ir ( max_samp_line*max_samp_pixel )

c ... Local scalars

      double precision east_lim, west_lim, north_lim, south_lim
      character*70 begin_date, begin_time, end_date, end_time
      character*4 pcf_satid
      real pcn4s, pcnnd, pcnnt, pcnng, pcnni, pcnnl, 
     &     pcnnw, pcnns, pcnnv, pcnnr, pcnnc, max_sol, min_sol,
     &     gfit, Chisq
      integer  h_eco1, nscans, npixels, scan,
     &         npix, ng_irp, ssctpr, slcd, smcd, shcd, sthncd,
     &         sthkcd, sopcd, scicd, ni, nw, nx, nu, qual_flag,
     &         conf_flag, ci_flag, hc_flag, out_lines_5km, 
     &         out_elements_5km,maxl, clus, nmaxl, line, pixel,
     &         ncloud, ncloud_co2, nmissing, npix_scan, ngood_msk, 
     &         beg_lin, nlins, beg_ele, neles, ibes, beg_scan,
     &         ngood_dayIRPhase, ngood_IRPhase, ngood_co2,
     &         sco2, sbias, os_top_flag, cldhgt_cat, nearnad_flag,
     &         cldhgtmet_flag
      integer start( 2 ), stride( 2 ), edge( 2 ), rtn
      character*100 text

c ... External functions

      integer copy_info_mod06
      external copy_info_mod06
      integer pgs_io_gen_closef
      external pgs_io_gen_closef
      
c ... Declarations for function 'get_platform_name'

      integer fileid, version, status
      integer get_platform_name
      external get_platform_name

      save /platform_name_common/

c ... whuang 05/09/01:
c ... Added save statement for common /mod06_data/ and /mod06_debug/
      save /mod06_data/, /mod06_debug/

      integer num_args
      integer FlagRA
      character FlagBuff*10
      integer iargc

      num_args = iargc ( )
      
      if(num_args .eq. 1) then
         call getarg ( 1, FlagBuff )
         read (FlagBuff,*) FlagRA
      else
C      This is the default value
         FlagRA = 0
      endif

C#######################################################################
C     INITIALIZATION
C#######################################################################

      error_flag = 0

c 03/10/98 fhliang added following block for data initialization:
c (bands, ir_mband, and ir_aband are defined in mod06uw_debug.inc)

      bands(1) = 20
      bands(2) = 22
      bands(3) = 23
      bands(4) = 24
      bands(5) = 25
      bands(6) = 27
      bands(7) = 28
      bands(8) = 29
      bands(9) = 30
      bands(10) = 31
      bands(11) = 32
      bands(12) = 33
      bands(13) = 34
      bands(14) = 35
      bands(15) = 36
      bands(16) = 1
      bands(17) = 6
      bands(18) = 26

! list bands required for IR cloud phase algorithm
      ir_mband(1) = 29
      ir_mband(2) = 31
      ir_mband(3) = 32

      ir_aband(1) = 8
      ir_aband(2) = 10
      ir_aband(3) = 11

! list bands needed for CO2 slicing      
      ir_co2band(1) = 36
      ir_co2band(2) = 35
      ir_co2band(3) = 34
      ir_co2band(4) = 33
      ir_co2band(5) = 31
      ir_co2band(6) = 29
      
      ir_co2bandptr(1) = 15
      ir_co2bandptr(2) = 14
      ir_co2bandptr(3) = 13
      ir_co2bandptr(4) = 12
      ir_co2bandptr(5) = 10
      ir_co2bandptr(6) = 8

! list bands for daytime cloud phase algorithm      
      day_band(1) = 29
      day_band(2) = 31
      day_band(3) = 32
      day_band(4) =  1
      day_band(5) =  6
      day_band(6) = 26

      day_bandptr(1) = 8
      day_bandptr(2) = 10
      day_bandptr(3) = 11
      day_bandptr(4) = 16
      day_bandptr(5) = 17
      day_bandptr(6) = 18

c ... Get platform name from MOD021KM file ('Terra' or 'Aqua')
      if( FlagRA .eq. 1) then
         fileid = cube_pcfnum_RA
      else
         fileid = cube_pcfnum
      endif
      version = 1
      rtn = get_platform_name(fileid, version, status, platform_name)
      write(*,'(''platform name:'', 3i10, a50)') fileid, version, status, platform_name

c ... Check return value from get_platform_name
      if (rtn .ne. 0) then
        call message( 'MOD_PR06CT', 'Error getting platform name' //
     &    ' [OPERATOR ACTION: Contact SDST]', 0, 2)
      endif

c ... Check the platform name (it should be 'Terra' or 'Aqua')
      if (platform_name(1:5) .ne. 'Terra' .and.
     &    platform_name(1:4) .ne. 'Aqua') then
        call message( 'MOD_PR06CT', 'Platform name not recognized' //
     &    ' [OPERATOR ACTION: Contact SDST]', 0, 2)
      endif
      
c ... Get debug value and processing interval information from PCF

      call mod06_get_rp(debug,beg_lin,nlins,beg_ele,neles,pcf_satid)

c ... Open required files for this processing section

      call mod06uw_open_files( cube_file, geo_file, cld_file,
     &  out_file, bug_file, cube_handle, geo_handle, cld_handle,
     &  out_handle, h_eco1, h_output, error_flag )
      

c ... If there was a problem opening required files,
c ... then print a warning message and return to the caller
c ... (this would allow another algorithm to execute after this one)

      if ( error_flag .ne. 0 ) then
        call message( 'MOD_PR06UW',
     &    'Failed to open a required input file. ' //
     &    ' [OPERATOR ACTION: See previous error messages.' //
     &    ' If problem persists, contact SDST]', 0, 1 )
        return
      endif

c ... Find bounds of granule from reading ECS Core Metadata

      call mod06_Get_Core_Metadata(geo_handle,geo_pcfnum,mcf_pcfnum,
     +                             East_Lim,West_Lim,North_Lim,
     +                             South_Lim,Begin_Date,Begin_Time,
     +                             End_Date,End_Time)

c ... Get number of scans and pixels for 1km bands

      call mod06_Get_Global_Metadata( nscans, npixels )


c ... Check User defined processing intervals versus scan boundaries

      call mod06_chk_input(nscans,npixels,beg_lin,nlins,beg_ele,
     +                     neles,pcf_satid,beg_scan,ibes,out_lines_5km,
     +                     out_elements_5km)

C#######################################################################
C     INITIALIZATION OF GRANULE BASED VARIABLES
C#######################################################################

      npix = 0

c ... Initialize metadata counters for both products

      call init_granule_var_mod06ir( ni, nw, nx, nu, ng_irp )

      call init_granule_var_mod06ct( ssctpr, slcd, smcd,
     &  shcd, sthncd, sthkcd, sopcd, scicd, sco2, sbias )

C#######################################################################
C     START SCIENCE PROCESSING
C#######################################################################

c ... Begin scan loop

      do scan = beg_scan , nscans

c ...   Write debug info

        if ( debug .gt. 0 ) then
          write( h_output, '(/,72(''-''),/)' )
          write( h_output, '(''Processing scan '',i4)' ) scan
        endif

c ...   Extract radiance, geolocation, mask, and ancillary data

        call mod06_get_data( cube_handle, cube_file, geo_handle,
     &    cld_handle, h_eco1, scan , begin_date, begin_time, met_date)

c ...   Initialize output product arrays for this scan

        call mod06CT_initialize_output()

        call mod06IR_initialize_output( qual_flag_ir, conf_flag_ir,
     &    maxl, clus, gfit, Chisq, nmaxl, icode, indat_bt, var)

        npix_scan = 0

c ...   Start loop over pixel boxes (line and pixel is at corner of box)

        do line = 1 , ( max_line / isamp ) * isamp, isamp

          do pixel = ibes , ( npixels / isamp ) * isamp, isamp

             npix = npix + 1
             npix_scan = npix_scan + 1

c ...        Initialize scan based qa variables

             call init_scan_var_mod06ct( line, pixel, qual_flag, conf_flag,
     &         ci_flag, hc_flag, os_top_flag, cldhgt_cat, nearnad_flag,
     &         cldhgtmet_flag )

c ...        Check this box for good data

             call mod06_check_box( line, pixel, nmissing, ncloud, ncloud_co2,
     &         ngood_msk, ngood_co2, ngood_IRphase, ngood_dayIRPhase )

c ...        If this box has enough good pixels, compute products
c ...        First CO2 slicing

             if( ngood_co2 .ge. min_ngoodct ) then
               call mod06CT_compute_products( line, pixel, ngood_msk,
     &           ncloud, ncloud_co2, ngood_co2, met_date, qual_flag, conf_flag, 
     &           hc_flag, ci_flag, os_top_flag, cldhgt_cat, nearnad_flag,
     &           cldhgtmet_flag, ssctpr, slcd, smcd, shcd, sthncd, sthkcd,
     &           sopcd, scicd, sco2, sbias)
             end if

c ...        For IR phase, fill processing variables

             if ( ngood_IRphase .ge. min_ngood )
     &         call put_proc_var_irphase( npix_scan, line, pixel,
     &           ncloud, nmissing,ngood_msk,indat_bt,var)

c ...        Set the QA information for this box

             call mod06ct_set_qa( line, pixel, nmissing,
     &         ncloud, qual_flag, conf_flag, hc_flag, ci_flag, os_top_flag,
     &         cldhgt_cat, nearnad_flag, cldhgtmet_flag )

c ...        Fill up the common output arrays

             call fill_common_arrays(line,pixel)

          end do

        end do

c ...   Compute the cloud phase product for this scan

        call get_irphase( npix_scan, indat_bt, icode)

c ...   Set QA information for cloud phase and sum metadada counters

        call mod06ir_set_qa( npixels,qual_flag_ir, conf_flag_ir,
     &    maxl, clus, gfit, Chisq, nmaxl, icode, ng_irp,
     &    ni, nw, nx, nu )

c ...   Place phase processing arrays in correct output formatted arrays

        call put_irphase_out_arrays(npixels, icode, indat_bt, var)

c ...   Write products for this scan

        call mod06_write_products( out_handle, scan, out_elements_5km )

c ... End of loop over all scans

      end do

C#######################################################################
C     END SCIENCE PROCESSING
C#######################################################################

c ... Get metadata from the cloud mask input file
c 03/13/98 fhliang commented out mcf_cmnum.
c     call copy_cldmask_M( cld_handle, cld_pcfnum, mcf_cmnum,


      call copy_cldmask_M( cld_handle, cld_pcfnum,
     &        pcn4s, pcnnd, pcnnt, pcnng, pcnni, pcnnl, pcnnw, pcnns,
     &        pcnnv, pcnnr, pcnnc, max_sol, min_sol )

c ... Get metadata for MOD06

      call mod06uw_stats_out( pcnnl, pcnnw, pcnni, ssctpr, ng_irp,
     &  npix, slcd, smcd, shcd, sthncd, sthkcd, sco2, sbias,
     &  sopcd, scicd, ni, nw, nx, nu, success_co2,
     &  success_irp, plow, pmid, phigh, pthin,
     &  pthick, popq, pcirrus, pice, pwater,
     &  pmixed, punc, pocean, pland, psnow )

c ... close files for this section of code

      call mod06_close_files( cube_handle, geo_handle, cld_handle,
     &  out_handle, h_eco1 )

c ... Set the runtime debug file metadata

      call mod06_set_runtime_meta()

c ... Copy geolocation information to output file

      start( 1 ) = isamp / 2
      start( 2 ) = isamp / 2
      stride( 1 ) = isamp
      stride( 2 ) = isamp
      edge( 1 ) = out_elements_5km
      edge( 2 ) = out_lines_5km
      text = ' '
      rtn = -1

      rtn = copy_info_mod06( geo_file, out_file, start, stride, edge, text )
      if ( debug .gt. 0 ) then
        write( h_output, '(/,72(''-''),/)' )
        write( h_output, '(''copy_info_mod06 debug info'',/)' )
        if ( rtn .ne. 0 ) then
          write( h_output, '(''Geolocation copy failed:'',/,a)' ) text
        else
          write( h_output, '(''Geolocation copy successful'')' )
        endif        
      endif
      
      if ( rtn .ne. 0 ) call message( routine_name,
     &  'Geolocation copy failed:' // text //
     &  ' [OPERATOR ACTION: Check available disk space.' //
     &  ' If problem persists, contact SDST]', 0, 2 )

c ... Write debug info and then close debug file

      if ( debug .gt. 0 ) then
        write( h_output, '(/,72(''-''),/)' )
        write( h_output, '(''MOD_PR06UW terminated normally'',/)' )
      endif

      rtn = -1
      rtn = pgs_io_gen_closef( h_output )
      if ( rtn .ne. 0 ) call message( 'MOD_PR06UW',
     &  'Close failed on runtime debug file ' // bug_file //
     &  ' [OPERATOR ACTION: Check available disk space.' //
     &  ' If problem persists, contact SDST]', 0, 2 )

      END
