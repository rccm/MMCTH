      program Driver_MOD_PR06UW

      implicit none
      
c-----------------------------------------------------------------------
c
c!F77
c
c!Description:
c    Driver program for UW MODIS Version 2 Cloud Top Properties.
c
c!Input Parameters:
c    See the accompanying README
c
c!Output Parameters:
c    See the accompanying README
c
c!Revision History:
c
c $Id: Driver_MOD_PR06UW.f,v 1.1.1.1 2005/02/22 17:15:54 gumley Exp $
c
c!Team-unique Header:
c
c!End
c
c
c-----------------------------------------------------------------------

C#######################################################################
C     DECLARATIONS
C#######################################################################

c ... Include files

      include 'mod06uw_pcfnum.inc'
      include 'mod06uw_meta.inc'
      include 'mapi.inc'
      include 'Atmos_AncData.inc'
      include 'PGS_PC.f'
      include 'PGS_SMF.f'
      
c ... Local variables

      character*32 routine_name
      parameter ( routine_name = 'Driver_MOD_PR06UW' )

      character*( PGSd_PC_FILE_PATH_MAX ) out_file, csbias_file

      integer file_version, csb_vrsn /1/
      integer out_handle( MODFILLEN )      

      integer error_flag
      
c ... Variables needed for core and archive metadata insertion
c ... (array sizes set in mod06uw_meta.inc include file)

      character*( PGSd_MET_GROUP_NAME_L )
     &  MET_Groups( PGSd_MET_NUM_OF_GROUPS )

      character*30 Auto_QA_Flag( No_ParaName )

      real PSAValue( No_PSA )

      integer LRN_MODIS_InputFiles( NUM_MODIS_InputFiles )
      integer Input_LRN(No_Input_LRN)
      integer QA_Miss(No_ParaName)
      integer rtn, i

c ... Product Specific Attributes (PSA's)

      real success_co2, success_irp, plow, pmid, phigh, pthin,
     &     pthick, popq, pcirrus, pice, pwater, pmixed, punc,
     &     pocean, pland, psnow

c ... Metadata missing value

      real        meta_misg
      parameter ( meta_misg = -99999.0 )

c     Metadata variables

      integer     INVENTORYMETADATA
      parameter ( INVENTORYMETADATA = 2 )

      integer     ARCHIVEMETADATA
      parameter ( ARCHIVEMETADATA = 3 )

      integer     NUM_of_HDFAttrNms
      parameter ( NUM_of_HDFAttrNms = 2 )

c ... rhucek 11/18/05 - metadata keywork that enables skipping an 
c     entry in the InputPointer list 
      integer     NO_SET_UR 
      parameter ( NO_SET_UR = -99999 )

      character*( MAX_ECS_NAME_L - 1 )
     &  HDFAttNms( PGSd_MET_NUM_OF_GROUPS )

c ... Variables for sdptk temporary file cleanup 

c ... rhucek 12/12/02:  added temp file cleanup variables 
c     rhucek 11/18/05: sea_ice not used by MOD_PR06CT, removed from 
c     temporary files list. 
c  
c     integer NumTempFiles /2/, 
c    &        TempFile_LUNList(2) /LUN_TEMP_GDAS_0ZF, LUN_TEMP_SEA_ICE/
      integer NumTempFiles /1/, TempFile_LUNList(1) /LUN_TEMP_GDAS_0ZF/

c ... External functions

c rhucek 12/12/02:  added Delete_CommonTempFiles_Atmos to external
c function list
      integer  Delete_CommonTempFiles_Atmos, Set_InvMet_MOD06, 
     &         Set_ArchMet_MOD06
      external Delete_CommonTempFiles_Atmos, Set_InvMet_MOD06, 
     &         Set_ArchMet_MOD06

      external clmfil, cpmfil
      
      integer  pgs_pc_getreference
      external pgs_pc_getreference

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

c ... Get output file name from process control file

      file_version = 1
      out_file = ' '
      rtn = PGS_S_SUCCESS - 1
      rtn = pgs_pc_getreference( out_pcfnum, file_version, out_file )
      if( rtn .ne. PGS_S_SUCCESS ) call message( routine_name,
     &  'Error getting output filename from PCF.' //
     &  ' [OPERATOR ACTION: Check PCF entry]', 0, 2 )
      
C#######################################################################
C     START SCIENCE PROCESSING
C#######################################################################

c ..  initialize the output metadata fields

      success_co2 = 0.0
      success_irp = 0.0
      plow    =  meta_misg
      pmid    =  meta_misg
      phigh   =  meta_misg
      pthin   =  meta_misg
      pthick  =  meta_misg
      popq    =  meta_misg
      pcirrus =  meta_misg
      pice    =  meta_misg
      pwater  =  meta_misg
      pmixed  =  meta_misg
      punc    =  meta_misg
      pocean  =  meta_misg
      pland   =  meta_misg
      psnow   =  meta_misg

c ... Run the IR Cloud Top Pressure and IR Cloud Phase algorithms

      call MOD_PR06UW( success_co2, success_irp, plow, pmid,
     &                 phigh, pthin, pthick, popq, pcirrus,
     &                 pice, pwater, pmixed, punc, pocean,
     &                 pland, psnow, error_flag )

C#######################################################################
C     END SCIENCE PROCESSING
C######################################################################

C######################################################################
C     START UW PRODUCT MOD06 CORE AND METADATA WRITING
C######################################################################

c ... Initialization of the output holder of the ODL metadata

      do i = 1, PGSd_MET_NUM_OF_GROUPS
         MET_Groups(i) = ' '
      end do

c ... Assign correct PSA's to inventory metadata

      PSAValue(1) = success_co2
      PSAValue(2) = success_irp
      PSAValue(3) = plow
      PSAValue(4) = pmid
      PSAValue(5) = phigh
      PSAValue(6) = pthin
      PSAValue(7) = pthick
      PSAValue(8) = popq
      PSAValue(9) = pcirrus
      PSAValue(10) = pice
      PSAValue(11) = pwater
      PSAValue(12) = pmixed
      PSAValue(13) = punc
      PSAValue(14) = pocean
      PSAValue(15) = pland
      PSAValue(16) = psnow
      
c ... Assign our 1 measured parameter percent missing

      do i = 1 , No_ParaName
        QA_Miss(i) = nint( 100.0 - success_co2 )
      enddo

c ... Now for the fun part.  Assign the Automatic QA flag
c ... Decide based upon bad data percentage
c ... rhucek 02/26/98:  Changed valid values of Auto_QA_Flag from
c ... "1"/"0"/ to "Passed"/"Failed" to conform with ECS data model. 

      do i = 1 , No_ParaName
        if ( success_co2 .lt. 10.0 ) then
          Auto_QA_Flag(i) =  'Failed'
        else
          Auto_QA_Flag(i) =  'Passed'
        endif
      enddo


c ... Now assign PCF LUNs and Versions to MODIS L1B, geolocation, MOD35 
c ... and CSBIAS files.  

      if( FlagRA .eq. 1) then
         LRN_MODIS_InputFiles(1) = cube_pcfnum_RA
      else
         LRN_MODIS_InputFiles(1) = cube_pcfnum
      endif
      LRN_MODIS_InputFiles(2) = geo_pcfnum
      LRN_MODIS_InputFiles(3) = cld_pcfnum
      LRN_MODIS_InputFiles(4) = csb_pcfnum

c ... Now assign the pcf number of input files to the correct
c ... variables 
      if( FlagRA .eq. 1) then
         Input_LRN(  1 ) = cube_pcfnum_RA
      else
         Input_LRN(  1 ) = cube_pcfnum
      endif
      Input_LRN(  2 ) = geo_pcfnum
      Input_LRN(  3 ) = cld_pcfnum
      Input_LRN(  4 ) = anc1_pcfnum
      Input_LRN(  5 ) = NO_SET_UR   ! reserved for anc2_pcfnum
      Input_LRN(  6 ) = anc3_pcfnum
      Input_LRN(  7 ) = NO_SET_UR   ! reserved for anc4_pcfnum
      Input_LRN(  8 ) = eco1_pcfnum
      Input_LRN(  9 ) = dry_pcfnum
      Input_LRN( 10 ) = ozo_pcfnum
      Input_LRN( 11 ) = wts_pcfnum
      Input_LRN( 12 ) = wtl_pcfnum
      Input_LRN( 13 ) = wco_pcfnum
      Input_LRN( 14 ) = csb_pcfnum

c ... Check to verify that optional CSBIAS product is present in PCF.
c ... If not, reset CSBIAS arrays elements to NO_SET_UR.

      if ( pgs_pc_getreference(csb_pcfnum, csb_vrsn, csbias_file) 
     &     .ne. PGS_S_SUCCESS ) then
         Input_LRN( 14 )         = NO_SET_UR
         LRN_MODIS_InputFiles(4) = NO_SET_UR
      endif

           
c ... Set Core Metadata

      rtn = MAPIOK - 1
      rtn = Set_InvMet_MOD06( ExtGeoPntr_Flag,
     &                        No_Input_LRN, Input_LRN, Input_Vrsn,
     &                        No_ParaName, Value_Para_Name, QA_Miss,
     &                        Auto_QA_Flag, Auto_QA_Expl,
     &                        No_PSA, PSAName, PSAValue,
     &                        MET_Groups )
      if ( rtn .ne. MAPIOK ) call message( routine_name,
     &  'Set_CoreMetadata failed. ' //
     &  ' [OPERATOR ACTION:  Contact the SDST]', 0, 1 )

c ... Set Archive Metadata

      rtn = MAPIOK - 1
      rtn = Set_ArchMet_MOD06( MET_Groups,
     &                         NUM_MODIS_InputFiles,
     &                         LRN_MODIS_InputFiles,
     &                         VRSN_MODIS_InputFiles,
     &                         NUM_ArchivePSA_SC,
     &                         Name_ArchivePSA_SC,
     &                         Value_ArchivePSA_SC )
      if ( rtn .ne. MAPIOK ) call message( 'routine_name',
     &  'Set_ArchiveMetadata failed.' //
     &  ' [OPERATOR ACTION:  Contact the SDST]', 0, 1 )

c ... Open output file

      rtn = -1
      rtn = opmfil( out_file, 'a', out_handle )
      if ( rtn .ne. MAPIOK ) call message( routine_name,
     &  'Error opening output MOD06_L2 file.' //
     &  ' [OPERATOR ACTION: Check system resources]', 0, 2 )

c ... Complete the output file with all metadata fields

      HDFAttNms( INVENTORYMETADATA ) = MECS_CORE
      HDFAttNms( ARCHIVEMETADATA   ) = MECS_ARCHIVE
      rtn = -1
      rtn = cpmfil( out_handle, MET_Groups, HDFAttNms,
     &  NUM_of_HDFAttrNms )
      if( rtn .ne. MAPIOK ) call message( routine_name,
     &  'Error completing output MOD06_L2 file.' //
     &  ' [OPERATOR ACTION: Check system resources]', 0, 2 )

C######################################################################
C     END UW PRODUCT MOD06 CORE AND METADATA WRITING
C######################################################################

c ... remove ancillary temporary files in PCF

c rhucek 12/12/02: referencing new version of sdptk temp file cleanup 
c function
      rtn = Delete_CommonTempFiles_Atmos(TempFile_LUNList, NumTempFiles)


c ... Return exit code 0 to shell

      call exit( 0 )

      end
