      subroutine Prepare_ECS_Metadata(pcbad,pcn1s,pcn2s,pcn3s,
     +                                pcn4s,pcnncd,pcnncr,pcnnd,
     +                                pcnnt,pcnng,pcnni,pcnnl,
     +                                pcnnw,pcnns,pcnnv,pcnnr,
     +                                pcnnc,max_sol,min_sol,
     +                                no_250,MET_groups)

      implicit none

      include 'mapi.inc'
      include 'mod35_meta.inc'
      include 'mod35.inc'
      include 'global.inc'



C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C
C    This program is designed to driver the Core, Archive and
c    Inventory Metadata subroutines which place the values
c    into storage until they are attached to the product file
c    using mapi CPMFIL.
c
c!Input Parameters:
c pcbad         Percentage of bad pixels
c pcn1s         Percentage of pixels in category 1 per granule
c pcn2s         Percentage of pixels in category 2 per granule
c pcn3s         Percentage of pixels in category 3 per granule
c pcn4s         Percentage of pixels in category 4 per granule
c pcnncd        Number of 250m pixels found to be cloudy per granule
c               (out of all possible pixels)
c pcnncr        Number of 250m pixels found to be clear per granule
c               (out of all possible pixels)
c pcnnd         Percentage of day processed pixels per granule
c pcnnt         Percentage of night processed pixels per granule
c pcnng         Percentage of sunglint found pixles per granule
c pcnni         Percentage of snow/ice processed pixels per granule
c pcnnl         Percentage of land processed pixels per granule
c pcnnw         Percentage of water processed pixels per granule
c pcnns         Percentage of shadow pixels found per granule
c pcnnv         Percentage of thin cirrus (vis) found pixels per granule
c pcnnr         Percentage of thin cirrus (ir) found pixels per granule
c pncnc         Percentage of non-cloud obstruction pixels per granule
c max_sol       Maximum solar zenith angle for this granule
c min_sol       Minimum solar zenith angle for this granule
c no_250        logical variable set to .true. if L1B 250 data unavailable
c
c!Output Parameters:
c MET_Groups    Array containing the names of the "MASTER" groups
c               as defined in the MCF.  Should be allocated in calling 
c               routine as: 
c               character*(PGSd_MET_GROUP_NAME_L) MET_Groups(PGSd_MET_NUM_OF_GROUPS)
c
C !REVISION HISTORY:
c
c !Team-unique Header:
c    This software is developed by the MODIS Science Data Support Team
c    for the National Aeronautics and Space Administration,
c    Goddard Space Flight Center, under contract NAS5-32373.
c
C !REFERENCES AND CREDITS:
C
c    Modified by Kathy Strabala from a driver program written by
C       Vicky Lin                    June 1997
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd,  Seabrook MD 20706
C    vlin@ltpmail.gsfc.nasa.gov
C
C !DESIGN NOTES:
C    Externals:
C          PGSd_MET_NUM_OF_GROUPS (PGS_MET.f:included in "mapi.inc")
c          A whole mess of variables are defined and set in
c            mod35_meta.inc
C
c
C !END
C---------------------------------------------------------------------

c ... scalar arguments
      logical no_250

      real pcbad,pcn1s,pcn2s,pcn3s,pcn4s,pcnncd,pcnncr,pcnnd,pcnnt,
     *     pcnng,pcnni,pcnnl,pcnnw,pcnns,pcnnv,pcnnr,pcnnc,max_sol,
     *     min_sol

c ..  array arguments
      character*(*) MET_Groups(*)

c ... parameter statements
      integer     MAX_NUM_INPUT
      parameter ( MAX_NUM_INPUT = 25 )

      integer     MAX_NUM_MODIS_Inputs
      parameter ( MAX_NUM_MODIS_Inputs = 30)

      integer     NO_SET_UR 
      parameter ( NO_SET_UR = -99999 )

c ... local scalars
      integer i,rtn,LRN_Prod_MCF,
     &        NUM_Input_LRN, 
     &        NUM_MODIS_InputFiles

c ... local arrays
      character*30 Auto_QA_Flag(No_ParaName),
     *             Value_Para_Name(No_ParaName)
      character*50 Auto_QA_Expl(No_ParaName)

      integer QA_Miss(No_ParaName),
     &        Input_LRN(MAX_NUM_INPUT),
     &        Input_Vrsn(MAX_NUM_INPUT),
     &        MODIS_InputFiles_LRN(MAX_NUM_MODIS_Inputs),
     &        VRSN_MODIS_InputFiles(MAX_NUM_MODIS_Inputs)

      real PSAValue(No_PSA),Value_ArchivePSA_SC(NUM_ArchivePSA_SC)


c ... external functions
      integer Set_ArchiveMetadata, Set_CoreMetadata
      external Set_ArchiveMetadata, Set_CoreMetadata

c ... external subroutines
      external message

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

c ... Initialization of the output holder of the ODL metadata
      Value_Para_Name(1) = 'Cloud_Mask'
      do i = 1, PGSd_MET_NUM_OF_GROUPS
         MET_Groups(i) = ' '
      end do

c ................... PLACE META VALUES IN CORRECT VARIABLES ..........
c ....................................................................
c ... Assign correct PSA's to inventory metadata
c ... Have already assigned the missing data value to the
c ...  PSA's if all bad data.
      PSAValue(1) = 100. - pcbad
      PSAValue(2) = pcn1s
      PSAValue(3) = pcn2s
      PSAValue(4) = pcn3s
      PSAValue(5) = pcn4s
      PSAValue(6) = pcnncd
      PSAValue(7) = pcnncr
      PSAValue(8) = pcnnd
      PSAValue(9) = pcnnt
      PSAValue(10) = pcnng
      PSAValue(11) = pcnni
      PSAValue(12) = pcnnl
      PSAValue(13) = pcnnw
      PSAValue(14) = pcnns
      PSAValue(15) = pcnnv
      PSAValue(16) = pcnnr
      PSAValue(17) = pcnnc
      if (max_sol .lt. -999.0) max_sol = Meta_misg
      if (min_sol .gt.  999.0) min_sol = Meta_misg
      PSAValue(18) = max_sol
      PSAValue(19) = min_sol


c ... Assign our 1 measured paramter which is percent missing
      QA_Miss(1) = nint(pcbad)

c ... Now for the fun part.  Assign the Automatic QA flag based
c     based on bad data percentage
      if (pcbad .gt. 90.) then
        Auto_QA_Flag(1) =  'Failed'
      else
        Auto_QA_Flag(1) =  'Passed'
      endif
      Auto_QA_Expl(1) =  'Passed if useable, Failed if not useable'

c ... Set our 1 archive PSA - code version number - set in include file
      Value_ArchivePSA_SC(1) = mod35_ver
c ....................................................................
c ....................................................................

C ... Assign LRN of MCF.
      LRN_Prod_MCF = LRN_MCF

c ... Assign PCF LUN and Version numbers of all input files 
c     (ECS InputPointer and LocalInputGranuleID)
      NUM_Input_LRN = 12 
      Input_LRN(1)  = LRN_Geo 
      if( FlagRA == 1) then
         Input_LRN(2)  = LRN_L1B_1km_RA
      else
         Input_LRN(2)  = LRN_L1B_1km
      endif
      Input_LRN(3)  = LRN_L1B_250 
      Input_LRN(4)  = LRN_eco1
      Input_LRN(5)  = LRN_eco2
      Input_LRN(6)  = LRN_THR_PAR 
      Input_LRN(7)  = LRN_ANC_GDAS 
c rhucek 06/22/04: retain 11 input file references, but do not set
c number 8 which is reserved for TOVS Ozone 
      Input_LRN(8)  = NO_SET_UR    !this LUN reserved for LRN_ANC_TOVS 
      Input_LRN(9)  = LRN_ANC_SST   
      Input_LRN(10) = LRN_ANC_ICE 
      Input_LRN(11) = LRN_ANC_NISE 
c rhucek 11/30/05: added LRN_DSCONF (L1B destriper config file)
      Input_LRN(12) = LRN_DSCONF 

      do i = 1, NUM_Input_LRN
         Input_Vrsn(i) = 1
      enddo

c ... exclude L1B 250m product from ECS metadata if file not staged 
      if (no_250) then
         Input_LRN(3) = NO_SET_UR 
      endif

C---------------------------------------------------------------------
C Set CoreMetadata
C---------------------------------------------------------------------

      rtn = Set_CoreMetadata( LRN_Prod_MCF,ExtGeoPntr_Flag,
     *                        NUM_Input_LRN,Input_LRN,Input_Vrsn,
     *                        No_ParaName,Value_Para_Name,
     *                        QA_Miss,Auto_QA_Flag,Auto_QA_Expl,
     *                        No_PSA,PSAName,PSAValue,
     *                        LUN_OF_NUM_INV_RP_PAIRS, 
     *                        MET_Groups )

      if (rtn.eq.mapiok) then
         call message('Prepare_ECS_Metadata',
     &   'Successfully completed Set_CoreMetadata without error ',
     &   0, 3 )
      else
         call message( 'Prepare_ECS_Metadata',
     &   'Set_CoreMetadata failed ' //
     &   ' [OPERATOR ACTION: Notify SDST]', 0, 1 )
      endif


C---------------------------------------------------------------------
C Set Archive Metadata
C---------------------------------------------------------------------

       rtn = Set_ArchiveMetadata( LRN_NUM_Of_ArchiveRP_Pairs,
     +                            NUM_Input_LRN, 
     +                            Input_LRN,
     +                            Input_Vrsn, 
     +                            NUM_ArchivePSA_SC,
     +                            Name_ArchivePSA_SC,
     +                            Value_ArchivePSA_SC,
     +                            MET_Groups )

      if (rtn.eq.mapiok) then
         call message('Prepare_ECS_Metadata',
     &   'Successfully completed Set_ArchiveMetadata without error ',
     &   0, 3 )
      else
         call message( 'Prepare_ECS_Metadata',
     &   'Set_ArchiveMetadata failed ' //
     &   ' [OPERATOR ACTION: Notify SDST]', 0, 1 )
      endif

      end
