      PROGRAM Write_Alg29Metadata

c-----------------------------------------------------------------------
c !F77
c
c !DESCRIPTION:
c     Driver program for writing LEOCAT Algorithm 29 Metadata
c
c !INPUT PARAMETERS:
c     See the accompanying README
c
c !OUTPUT PARAMETERS:
c     See the accompanying README
c
c !REVISION HISTORY:
c
c !TEAM-UNIQUE HEADER:
c     Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
c
c !END
c-----------------------------------------------------------------------

      IMPLICIT NONE

c ... PGS toolkit include files

      include 'PGS_PC.f'
      include 'PGS_SMF.f'
      INCLUDE 'PGS_PC_9.f'
      INCLUDE 'PGS_IO.f'
      INCLUDE 'PGS_IO_1.f'
      INCLUDE 'PGS_MODIS_39500.f'
      include 'mapi.inc'

c ... UW include files

      include 'mod35.inc'
      include 'global.inc'
      include 'debug.inc'
      include 'platform_name.inc'
      
c rhucek 12/17/02:  added SDST include file 
c ... SDST include files

      include 'Atmos_AncData.inc'

c ... Local variables

      character*128 leo_file, cld_file, bug_file
      character*255 debug_line
      integer cldrd_handle( MODFILLEN ), cldwrt_handle( MODFILLEN )

      character*(PGSd_MET_GROUP_NAME_L)
     &  met_groups( PGSd_MET_NUM_OF_GROUPS )
      real pcbad, pcgood, pcn1s, pcn2s, pcn3s, pcn4s, pcnnd, pcnnt,  
     &  pcnng, pcnni, pcnnl, pcnnw, pcnnv, pcnnr, pcnnc, pcnncd, pcnncr,
     &  max_sol, min_sol
      character*34 ndvifil
      character*23 thresfil
      logical no_250
      real pcnns
      integer rtn
      integer i

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

c ... External functions

      integer pgs_met_init, pgs_io_gen_closef
      external pgs_met_init, pgs_io_gen_closef
      external alg29_open_files, get_global_metadata_leocat,
     +        Prepare_ECS_Metadata, alg29_close_files, 
     +        test_for_250m_file

C-----------------------------------------------------------------------
C     INITIALIZATION
C-----------------------------------------------------------------------

c ... Get debug value from PCF

      call get_debug( debug )

c ... Initialize file handles
      do i = 1, MODFILLEN
         cldrd_handle(i) = 0
         cldwrt_handle(i) = 0
      end do

c ... Get input and output filenames,

      call alg29_open_files( leo_file, cld_file, 
     &  bug_file, cldrd_handle, cldwrt_handle, h_output )
      
      no_250 = .false.
      call test_for_250m_file(no_250)

c ... Get global metadata from LEOCAT

      call get_global_metadata_leocat( cldrd_handle, pcgood,
     &  pcn1s, pcn2s, pcn3s, pcn4s, pcnncd, pcnncr, pcnnd,
     &  pcnnt, pcnng, pcnni, pcnnl, pcnnw, pcnnv, pcnnr, pcnnc,
     &  max_sol, min_sol, ndvifil,thresfil)
 
      pcbad = 100.0 - pcgood
c     set these data to missing
      pcnns=Meta_misg
      write(debug_line,891) 'pcgood=',pcgood,'pcn1s=',pcn1s
891   format(a7,f9.2,1x,a7,f9.2)
      write( h_output, '(a)' ) debug_line
      write(debug_line,891) 'pcn2s=',pcn2s,'pcn3s=',pcn3s
      write( h_output, '(a)' ) debug_line
      write(debug_line,891) 'pcn4s=',pcn4s,'pcnncd=',pcnncd
      write( h_output, '(a)' ) debug_line
      write(debug_line,891) 'pcnncr=',pcnncr,'pcnnd=',pcnnd
      write( h_output, '(a)' ) debug_line
      write(debug_line,891) 'pcnnt=',pcnnt,'pcnng=',pcnng
      write( h_output, '(a)' ) debug_line
      write(debug_line,891) 'pcnni=',pcnni,'pcnnl=',pcnnl
      write( h_output, '(a)' ) debug_line
      write(debug_line,891) 'pcnnw=',pcnnw,'pcnnv=',pcnnv
      write( h_output, '(a)' ) debug_line
      write(debug_line,891) 'pcnnr=',pcnnr,'pcnnc=',pcnnc
      write( h_output, '(a)' ) debug_line
      write(debug_line,891) 'mxsol=',max_sol,'mnsol=',min_sol
      write( h_output, '(a)' ) debug_line
      write( h_output, '(a)' ) ndvifil
      write( h_output, '(a)' ) debug_line
      write( h_output, '(a)' ) thresfil


      call Prepare_ECS_Metadata(pcbad,pcn1s,pcn2s,pcn3s,pcn4s,
     +                          pcnncd,pcnncr,pcnnd,pcnnt,pcnng,
     +                          pcnni,pcnnl,pcnnw,pcnns,pcnnv,pcnnr,
     +                          pcnnc,max_sol,min_sol,no_250,met_groups,
     +                          ndvifil,thresfil)

c ... Set the debug runtime file metadata
      call Set_Runtime_Meta


c ... Close input files, and close/complete the output file

      call alg29_close_files( cldrd_handle, cldwrt_handle, met_groups )

c ... Write debug info and then close debug file

      if ( debug .eq. 1 ) then
        write( h_output, '(/,72(''-''),/)' )
        write( h_output, '(''Write_Alg29 terminated normally'',/)' )
      endif


c           Get satellite instrument name.
      rtn = pgs_pc_getconfigdata(LUN_Sat_Instrument,pcf_satid)
      if (rtn .ne. 0) then
        call message( 'Write_Alg29Metadata',
     &  'Error reading instrument name from pcf LUN 800510.' //
     &  ' [OPERATOR ACTION: Check PCF file, fix if corrupt.]',
     &  0, 1 )
      endif

      if (pcf_satid .eq. 'AM1M') then
         doi_char = '10.5067/MODIS/MOD35_L2.006'
         doiLen = 26
      else
         doi_char = '10.5067/MODIS/MYD35_L2.006'
         doiLen = 26
      endif

      sd_id = sfstart(cld_file, DFACC_WRITE)
      if(sd_id .eq. -1) then
         call message( 'Write_Alg29Metadata',
     &        'Problem openning the file', 0, 2 )
      endif

      rtn = sfscatt(sd_id, 'identifier_product_doi', DFNT_CHAR8, doiLen,  
     +                 doi_char) 
      if (rtn .eq. -1) then
         call message( 'Write_Alg29Metadata',
     &        'Problem writting the global attribute identifier_product_doi', 0, 2 )
      endif

      rtn = sfscatt(sd_id, 'identifier_product_doi_authority', DFNT_CHAR8, 17,  
     +                 'http://dx.doi.org') 
      if (rtn .eq. -1) then
         call message( 'Write_Alg29Metadata',
     &        'Problem writting the global attribute identifier_product_doi_authority', 0, 2 )
      endif

      rtn = sfend(sd_id)
      if (rtn .eq. -1) then
         call message( 'Write_Alg29Metadata',
     &        'Problem closing the file', 0, 2 )
      endif


c      rtn = -1
c      rtn = pgs_io_gen_closef( h_output )
c      if ( rtn .ne. 0 ) call message( 'Write_Alg29Metadata',
c     &  'Close failed on runtime debug file ' // bug_file //
c     &  ' [OPERATOR ACTION: Check available disk space.' //
c     &  ' If problem persists, contact SDST]', 0, 2 )
c
c      if ( rtn .ne. 0 ) call message( 'Write_Alg29Metadata',
c     &  'Failed to delete temporary files' // bug_file //
c     &  ' [OPERATOR ACTION: Contact SDST]', 0, 1 )

c ... Return exit code 0 to shell

      call exit( 0 )

      END
