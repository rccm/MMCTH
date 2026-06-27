       PROGRAM MOD_PR35

C----------------------------------------------------------------------
C!F77
C
C!DESCRIPTION:   Calculate the probability that a single MODIS pixel is 
C                unobstructed.  This particular software was developed 
C                using mas 50 channel data and tested upon the MODIS
C                V2 simulated data set. 
C                MODIS radiances and reflectances are passed into 
C                processing routines based upon time of day, latitude, 
C                and surface to determine a confidence level of 
C                onobstructed.It outputs the results into a 48 bit file 
C                containing the final confidence plus information on 
C                individual test results. 
C                This version uses a sliding 3x3 box 
C                approach to regional processing, to clarify scenes of 
C                lower confidence.
C
C!INPUT PARAMETERS:
c See file_open.f for static inputs.  MOD_PR35.pcf contains a list
c of all input files.
c 
c       The 19 MODIS channels used in this processing are:
c               4  -  .565 microns
c               1  -  .659
c               2  -  .865
c               17 -  .905
c               18 -  .936
c               19 -  .940
c               5  - 1.240
c               26 - 1.375
c               6  - 1.640
c               7  - 2.130
c               20 - 3.750
c               22 - 3.959
c               27 - 6.715
c               28 - 7.225
c               29 - 8.550
c               31 - 11.030
c               32 - 12.020
c               33 - 13.335
c               35 - 13.935
c
c!Output Parameters:
c MOD35.hdf    Output MODIS cloud mask in 48 bit structure
c runtime.output Output from write statements in code
c 
c!Revision History:
c
c!Team-unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-35.  
c 
c!Design Notes:
c
c   Externals:     
c     Subroutines init_input,file_open,init_meta,Get_Core_Metadata,
c        Get_Global_Metadata,chk_input,init_var,Get_Swath_Metadata,
c        get_L1b,get_angles,get_anc_data,store_context,store_meta,
c        chk_lin_edge,chk_ele_edge,pxinit,get_pxldat,
c        check_reg_uniformity,polar_day,polar_nite,
c        land_day,land_nite,water_day,water_nite,shadows,noncld_obs_chk,
c        thin_ci_chk_ir,proc_path,set_unused_bits,get_stats,set_confdnc,
c        set_quality_A,fill_bit_line,mod35_out,reorder_array,stats_out,
c        Prepare_ECS_Metadata,file_close,setup_create_hdf,get_rp,
c        Set_Runtime_Meta,get_geo_copy
c 
c!END
c----------------------------------------------------------------------

      IMPLICIT NONE

c     R. Frey added 12/09/04
c...  HDF include file
      INCLUDE 'hdf.inc'

c...  PGS toolkit include files
      INCLUDE 'PGS_SMF.f'
      INCLUDE 'PGS_PC.f'
      INCLUDE 'PGS_PC_9.f'
      INCLUDE 'PGS_IO.f'
      INCLUDE 'PGS_IO_1.f'
      INCLUDE 'PGS_MODIS_39500.f'
      INCLUDE 'mapi.inc'
c rhucek 12/17/02: added new include file
      INCLUDE 'Atmos_AncData.inc'


c ... UW include files      
      INCLUDE 'mod35.inc'
      INCLUDE 'global.inc'
      INCLUDE 'debug.inc'
      INCLUDE 'platform_name.inc'

c ... Local scalars

      character*(PGSd_PC_FILE_PATH_MAX) fname_1km, fname_250m,
     +           fname_mod35,fname_geo
      character*(PGSd_MET_GROUP_NAME_L)
     +           MET_Groups(PGSd_MET_NUM_OF_GROUPS)
      character*1 scan_type
      character*2 Prog_number
      character*70 Begin_G_Date,End_G_Date,Begin_G_Time,End_G_Time
      character*4 pcf_satid

      double precision Sstart_time,East_Lim,West_Lim,North_Lim,
     +                 South_Lim,TAI_start,TAI_end
      real confdnc,precip_water,vza,plat,plon,sfctmp,pmsl,u_wind,v_wind,
     +     pcbad,pcn1s,pcn2s,pcn3s,pcn4s,pcnnc,pcnns,pcnnd,pcnnt,pcnng,
     +     pcnni,pcnnl,pcnnw,pcnnr,pcnnv,pcnncr,pcnncd,max_sol,min_sol,
     +     pcnu,pcnn,pcnz,pcna,pcne,tbadj,refang
      integer outrec,cube,maxlin,maxele,beg_lin,nlins,beg_ele,neles,
     +        lines_in_edge,pixels_in_edge,max_count,ilin,nu,nn,nz,na,
     +        mlin,klin,I_ind,prev_cube,h_eco1,ne,
     +        h_eco2,line,line_count,nc,npix,n1sm,lsf,
     +        n2sm,n3sm,n4sm,num_bad,no,ns,nd,nt,ng,ni,nl,nw,
     +        nr,nv,n250,ncr,ncd,nmtests,nbands,nadir,scan_pixels,
     +        scan_number,nbad_1km,nbad_250,nscans,max_pixel,
     +        max5km_ele,max5km_lin,rtn

      integer Delete_CommonTempFiles_Atmos
      external Delete_CommonTempFiles_Atmos
      
      logical line_edge,ele_edge,polar,day,night,land,water,snglnt,
     +        visusd,snow,ice,uniform,bad_value,shadow,coast,vrused,
     +        desert,smoke,cirrus_ir,cirrus_vis,no_250,bad_geo,process,
     +        hi_elev,antarctic,sh_ocean,sg_bad_data,map_ice,map_snow,
     +        sh_lake
      byte eco_type 

c ... Local arrays

      character*1  s_type(nlcntx)
 
      real pxldat(inband),pxl250(vis_band,num250_per_1km),
     +     indat(npixel,nlcntx,inband),v250_dat(nx*4,nlcntx,vis_band),
     +     contx_lat(npixel,nlcntx),contx_lon(npixel,nlcntx),
     +     contx_elev(npixel,nlcntx),contx_sst(npixel,nlcntx),
     +     contx_vzen(npixel,nlcntx),contx_szen(npixel,nlcntx),
     +     contx_sfctmp(npixel,nlcntx),contx_msl(npixel,nlcntx),
     +     contx_ugrd(npixel,nlcntx),contx_vgrd(npixel,nlcntx),
     +     contx_rel(npixel,nlcntx),rad_band(nx,ny,inband),
     +     vis250_band(nx,ny,vis_band),rlat(npixel,scans_cube),
     +     rlon(npixel,scans_cube),view_zen(npixel,scans_cube),
     +     elev(npixel,scans_cube),
     +     solar_zen(npixel,scans_cube),rel_angle(npixel,scans_cube),
     +     cube_pw(npixel,scans_cube),cube_ice(npixel,scans_cube),
     +     cube_sfctmp(npixel,scans_cube),cube_msl(npixel,scans_cube),
     +     cube_ugrd(npixel,scans_cube),cube_vgrd(npixel,scans_cube),
     +     cube_sst(npixel,scans_cube),
     +     contx_pw(npixel,nlcntx),contx_ice(npixel,nlcntx),
     +     scan_confdnc(npixel)

      integer buf_size(2),buf_geo_size(2),buf_anc_size(2),
     +        modfil_L1B_1km(MODFILLEN),modfil_L1B_250m(MODFILLEN),
     +        modfil_mod35(MODFILLEN),modfil_Geo(MODFILLEN),
     +        cube_topog(npixel,scans_cube),
     +        contx_topog(npixel,nlcntx),s_pixels(nlcntx),geo_flag(5),
     +        n_nadir(nlcntx),
     +        cube_snow(npixel,scans_cube),contx_snow(npixel,nlcntx)

c rhucek 12/17/02:  define sdptk temp file luns
      integer NumTempFiles /2/,
     &        TempFile_LUNList(2) /LUN_TEMP_GDAS_0ZF, LUN_TEMP_SEA_ICE/

      logical scan_ice(npixel),scan_snow(npixel),scan_desert(npixel),
     +        scan_coast(npixel),scan_land(npixel),
     +        scan_snglnt(npixel),cm1km(num250_per_1km),
     +        ice_250(num250_per_1km),snow_250(num250_per_1km),
     +        desert_250(num250_per_1km),coast_250(num250_per_1km),
     +        land_250(num250_per_1km),
     +        snglnt_250(num250_per_1km),scan_visusd(npixel),
     +        scan_process(npixel),day_250(num250_per_1km),
     +        proc_250(num250_per_1km),visusd_250(num250_per_1km),
     +        qa_250(num250_per_1km),scan_day(npixel)

      byte bitarray(npixel,6),testbits(6),cube_eco(npixel,scans_cube),
     +     contx_eco(npixel,nlcntx),v250_band(nx,ny,vis_band),
     +     v1km_band(nx,ny,inband),qa_bits(10),qa_bitarray(npixel,10)

c      character*80 thresholds_version
      
c R. Frey 12/09/04 Added variables for reading/writing destriping flag.
      integer rcs_id,config_id,sffattr,sfrcatt,sfscatt,strlen
      character rcs_string*132,config_string*132
      external sffattr,sfrcatt,sfscatt,strlen

c ... External subroutines and functions

      external init_input,file_open,init_meta,Get_Core_Metadata,
     +        Get_Global_Metadata,chk_input,init_var,Get_Swath_Metadata,
     +        get_L1b,get_angles,get_anc_data,store_context,store_meta,
     +        chk_lin_edge,chk_ele_edge,pxinit,get_pxldat,
     +        check_reg_uniformity,polar_day,polar_nite,
     +        land_day,land_nite,water_day,water_nite,shadows,
     +        noncld_obs_chk,thin_ci_chk_ir,proc_path,set_unused_bits,
     +        get_stats,set_confdnc,set_quality_A,fill_bit_line,
     +        mod35_out,reorder_array,stats_out,Prepare_ECS_Metadata,
     +        file_close,setup_create_hdf,get_rp,Set_Runtime_Meta,
     +        get_geo_copy

c ... Declarations for function 'get_platform_name'

      integer fileid, version, status
      integer get_platform_name
      external get_platform_name
      
c ... Data statements

c ... Initialize output record number and output bitarray value
      data outrec/0/

c ... Setup up program number
      data Prog_number /'35'/

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
c ... End of declarations
C-----------------------------------------------------------------------

c ... Get platform name from MOD021KM file ('Terra' or 'Aqua')
      if( FlagRA == 1) then
         fileid = LRN_DSL1B_1km_RA
      else
         fileid = LRN_DSL1B_1km
      endif

      version = 1
      rtn = get_platform_name(fileid, version, status, platform_name)

c ... Check return value from get_platform_name
      if (rtn .ne. 0) then
        call message( 'MOD_PR35', 'Error getting platform name' //
     &    ' [OPERATOR ACTION: Contact SDST]', 0, 2)
      endif

c ... Check the platform name (it should be 'Terra' or 'Aqua')
      if (platform_name(1:5) .ne. 'Terra' .and.
     &    platform_name(1:4) .ne. 'Aqua') then
        call message( 'MOD_PR35', 'Platform name not recognized' //
     &    ' [OPERATOR ACTION: Contact SDST]', 0, 2)
      endif
      
c ... Extract debug value and processing bounds from pcf file
      call get_rp(debug,beg_lin,nlins,beg_ele,neles,pcf_satid)

c ... Intialize variables needed for inputs and outputs
      call init_input(modfil_L1B_1km,modfil_L1B_250m,
     +                modfil_mod35,modfil_Geo,
     +                buf_size,buf_geo_size,buf_anc_size,
     +                qa_bits,bitarray,qa_bitarray)

c ... Open all input data files and return unit numbers
      call file_open(qa_bits,h_eco1,h_eco2,modfil_L1B_250m,
     +               fname_250m,modfil_L1B_1km,fname_1km,
     +               modfil_Geo,fname_geo,fname_mod35,
     +               no_250)

c ... Initialize variables needed to perform metadata extraction
      call init_meta(East_Lim,West_Lim,North_Lim,South_Lim,
     +               Begin_G_Date,Begin_G_Time,End_G_Date,
     +               End_G_Time,TAI_start,TAI_end,nscans,
     +               max_pixel,lines_in_edge,pixels_in_edge)

c ... Find bounds of granule from reading ECS Core Metadata
      call Get_Core_Metadata(Modfil_Geo,LRN_Geo,LRN_MCF,East_Lim,
     +                       West_Lim,North_Lim,South_Lim,
     +                       Begin_G_Date,Begin_G_Time,End_G_Date,
     +                       End_G_Time)

c ... Get the element and scan numbers from L1B file
      call Get_Global_Metadata(modfil_L1B_1km,nscans,max_pixel)

c ... check consistency of input information.
      call chk_input(max_pixel,nscans,lines_in_edge,no_250,
     +               Begin_G_Date,Begin_G_Time,End_G_Date,End_G_Time,
     +               beg_lin,nlins,beg_ele,neles,pcf_satid,maxele,maxlin,
     +               max_count,max5km_lin,max5km_ele,TAI_start,TAI_end)

c ... Create the output HDF file and open the file using mapi routines
      call setup_create_hdf(Prog_Number,fname_mod35,maxlin,max_pixel,
     +                      max5km_lin,max5km_ele,modfil_mod35)

c ... Read thresholds from parameter file
c      call thresholds_read( thresholds_version )

c ... Save thresholds version string to output file
c      rtn = pmfin( modfil_mod35, 'RCS_Id', 'CHARACTER*(*)',
c     &  len( thresholds_version ), thresholds_version )
c      if ( rtn .ne. MAPIOK ) then
c        call message( 'MOD_PR35', 'Error saving thresholds version' //
c     &    ' [OPERATOR ACTION: Contact SDST]', 0, 2 )
c      endif
              
c ... Finally, initialize the granule based processing variables
      call init_var(cube,pxldat,pxl250,indat,v250_dat,
     +              contx_lat,contx_lon,contx_elev,contx_vzen,
     +              contx_szen,contx_rel,contx_pw,contx_eco,
     +              contx_topog,contx_ice,contx_snow,contx_sfctmp,
     +              contx_msl,contx_ugrd,contx_vgrd,contx_sst,s_type,
     +              s_pixels,n_nadir,nadir,scan_pixels,scan_number,
     +              scan_type,Sstart_time,max_sol,min_sol,ilin,mlin,klin,
     +              I_ind,prev_cube)

c ... ******************************************************************
c ... Begin input scan line loop.
c ... Monitor the number of lines ('nlcntx') needed to fill a context 
c ... region with the variable 'mlin'.  Determine cloud mask at the 
c ... center pixel of each context region. 
c ...
c ... For lines and pixels with incomplete context regions at the 
c ... boundaries of our processing region, determine mask, but do
c ... not use spatial variability routines.
c ... ******************************************************************

      do 100 line = beg_lin , max_count
         ilin = ilin + 1
         line_count = beg_lin - 1 + ilin
         I_ind = I_ind + 1
         cube = (line_count-1)/10 + 1
         mlin = mlin + 1

c-----------------------------------------------------------------------
c ... Top of if statment for reading a scan cube of data out of the
c ... granule - have to get new cube and geo info every 10 lines.
c-----------------------------------------------------------------------
         if ((cube .ne. prev_cube) .and. (ilin .le. maxlin)) then
            prev_cube = cube

c ...       Extract scan based metadata
            call Get_Swath_Metadata(modfil_L1B_1km,cube,scan_number,
     +                              scan_pixels,scan_type,
     +                              Sstart_time,nadir)

c ...       Extract radiance and reflectance information 
            call get_L1b(modfil_L1B_1km,fname_1km,modfil_L1B_250m,
     +                   fname_250m,no_250,cube,buf_size,scan_type,
     +                   v250_band,v1km_band,vis250_band,rad_band)

c ...       Get the geolocation file variables 
            call  get_angles(modfil_Geo,cube,scan_type,max_sol,
     +                       min_sol,buf_geo_size,rlat,rlon,elev,
     +                       view_zen,solar_zen,cube_topog,rel_angle)

c ...       Get the ancillary data
            call get_anc_data(h_eco1,h_eco2,cube,buf_anc_size,max_pixel,
     +                        rlat,rlon,TAI_end,cube_pw,cube_eco,
     +                        cube_ice,cube_snow,cube_sfctmp,cube_msl,
     +                        cube_ugrd,cube_vgrd,cube_sst,
     +                        qa_bits)

            I_ind = MOD(line_count-1,10) + 1 
         end if
   
C-----------------------------------------------------------------------
c ... read scan line radiance and reflectance data, and if necessary
c ... convert to brightness temperature for use in threshold testing
C-----------------------------------------------------------------------

         if (ilin .le. maxlin) then
            call store_context(mlin,I_ind,maxele,rlat,rlon,elev,
     +                        view_zen,solar_zen,rel_angle,v250_band,
     +                        vis250_band,rad_band,v1km_band,indat,
     +                        v250_dat,cube_pw,cube_eco,
     +                        cube_topog,cube_ice,cube_snow,cube_sfctmp,
     +                        cube_msl,cube_ugrd,cube_vgrd,cube_sst,
     +                        contx_lat,
     +                        contx_lon,contx_elev,contx_vzen,contx_szen,
     +                        contx_rel,contx_pw,contx_eco,contx_topog,
     +                        contx_ice,contx_snow,contx_sfctmp,
     +                        contx_msl,contx_ugrd,contx_vgrd,
     +                        contx_sst)

            call store_meta(mlin,scan_type,scan_pixels,nadir,
     +                      s_type,s_pixels,n_nadir)

         end if 
 
C-----------------------------------------------------------------------
C                   PROCESSING SECTION
c
c        The Cloud Mask is processed pixel by pixel using a sliding
c        box approach. Regional context data is stored for use
c        in uniformity tests for uncertain pixels.  The edge of 
c        the area in which you are processing, (outline of region)
c        will be processed but will not include uniformity tests
c        (exception is nighttime water scenes).
C----------------------------------------------------------------------

c        Check if current line is a part of the border
         call chk_lin_edge(klin,maxlin,lines_in_edge,line_edge)

         if (mlin .eq. nlcntx .or. line_edge) then

            klin = klin + 1

            go to 5000

c           Do 300 nc = beg_ele , maxele
   
c               Check if current element is a part of the border
c               call chk_ele_edge(nc,s_pixels,pixels_in_edge,klin,
c    +                            maxlin,lines_in_edge,ele_edge)

c               Initialize regional variables
c               call pxinit(testbits,qa_bits,geo_flag,precip_water,vza,
c    +                  sfctmp,pmsl,u_wind,v_wind,plat,plon,lsf,polar,
c    +                  day,night,land,water,coast,snglnt,visusd,vrused,
c    +                  snow,ice,desert,bad_value,bad_geo,uniform,
c    +                  shadow,smoke,cirrus_ir,cirrus_vis,nmtests,
c    +                  nbands,nbad_1km,nbad_250,hi_elev,antarctic,
c    +                  sh_ocean,sg_bad_data,map_ice,map_snow,sh_lake)

c               Get pixel values for processing variables
c               call get_pxldat(indat,lines_in_edge,contx_lat,
c    +                          contx_lon,contx_vzen,contx_szen,
c    +                          contx_elev,contx_rel,contx_pw,contx_eco,
c    +                          contx_topog,contx_ice,contx_snow,
c    +                          contx_sfctmp,contx_msl,contx_ugrd,
c    +                          contx_vgrd,contx_sst,
c    +                          s_type,klin,beg_lin,
c    +                          maxlin,nc,nbands,nbad_1km,pxldat,
c    +                          plat,plon,precip_water,vza,
c    +                          sfctmp,pmsl,u_wind,v_wind,eco_type,
c    +                          lsf,geo_flag,polar,day,night,land,water,
c    +                          coast,snglnt,visusd,ice,vrused,desert,
c    +                          snow,bad_value,bad_geo,hi_elev,tbadj,
c    +                          antarctic,sh_ocean,sg_bad_data,map_ice,
c    +                          map_snow,refang,sh_lake)

c               Check for consistency of variables within processing
c               region
c               call check_reg_uniformity(line_edge,ele_edge,contx_eco,
c    +                                   contx_topog,contx_snow,n_nadir,
c    +                                   nc,eco_type,day,land,
c    +                                   water,coast,snow,ice,uniform)

c               Decision to process or not to process this pixel.
c               if(sg_bad_data .and. water .and. snglnt) then
c                 process = .false.
c               else
c                 process = .true.
c               end if

c               if(process) then

c                 Decision tree for processing paths 
c                 if (polar) then

c                   First polar processing

c                   if (day) then
c                     Daytime processing
c                     call polar_day(pxldat,vza,snglnt,
c    +                               visusd,refang,vrused,cirrus_vis,
c    +                               land,ice,snow,desert,coast,
c    +                               eco_type,uniform,hi_elev,nc,indat,
c    +                               nmtests,testbits,tbadj,antarctic,
c    +                               sh_ocean,sfctmp,qa_bits,confdnc)
c                   else
c                     Nighttime processing
c                     call polar_nite(pxldat,vza,land,ice,snow,desert,
c    +                                hi_elev,sfctmp,eco_type,nmtests,testbits,
c    +                                uniform,indat,nc,antarctic,sh_ocean,
c    +                                qa_bits,confdnc)
c                   endif
                     
c                 
c                 else if (land) then

c                   Primarily land surface.

c                   if (day) then
c                     Daytime processing.
c                     call land_day(pxldat,vza,visusd,vrused,
c    +                              cirrus_vis,desert,coast,snow,
c    +                              ice,hi_elev,tbadj,eco_type,
c    +                              testbits,qa_bits,nmtests,confdnc)

c                   else
c                     Nighttime processing.
c                     call land_nite(pxldat,plat,vza,ice,snow,coast,tbadj,
c    +                               desert,hi_elev,sh_lake,sfctmp,eco_type,
c    +                               nmtests,testbits,qa_bits,confdnc)
c                   endif

c                 else if (water) then

c                   Primarily water surface.

c                   if (day) then
c                     Daytime processing.
c                     call water_day(pxldat,vza,snglnt,visusd,refang,
c    +                               cirrus_vis,sfctmp,hi_elev,uniform,
c    +                               ice,maxele,klin,line_edge,
c    +                               sh_ocean,indat,nc,nmtests,
c    +                               testbits,qa_bits,confdnc)
c                   else
c                     Nighttime processing.
c                     call water_nite(pxldat,vza,uniform,ice,indat,sfctmp,
c    +                                sh_ocean,nc,nmtests,testbits,qa_bits,
c    +                                confdnc)
c                   end if

c                 end if


c                 Test for shadows,if necessary. Set bit to indicate 
c                 no shadow was found.
c                 if(.not.water .and. .not.coast .and. day .and. 
c    +               .not.polar .and. confdnc.ge.0.66) then
c                   call shadows(pxldat,shadow,visusd,qa_bits)
c                 end if

c                 Test for possible non-cloud obstruction.
c                 if(land .and. day .and. .not. snow) then
c                   call noncld_obs_chk(indat,pxldat,confdnc,nc,
c    +                                  maxele,line_edge,klin,qa_bits,
c    +                                  testbits,smoke)
c                 end if

c                 Test of thin cirrus in the infrared
c                 if ( (.not. snow) .and. (.not. ice) ) then
c                   call thin_ci_chk_ir(pxldat,vza,cirrus_ir,
c    +                                  qa_bits,testbits)
c                 end if

c               endif

c               Set bits which indicate processing path through 
c               algorithm.
c               call proc_path(water,land,day,ice,snow,snglnt,coast,
c    +                         desert,smoke,shadow,testbits)


c               Set the bits which are not used in output array.
c               call set_unused_bits(testbits)


c               Get cloud mask statistics
c               call get_stats(day,land,water,snglnt,snow,coast,desert,
c    +                         ice,shadow,smoke,cirrus_ir,cirrus_vis,
c    +                         confdnc,nmtests,nbands,bad_value,bad_geo,
c    +                         geo_flag,no,ns,nd,nt,ng,
c    +                         ni,nl,nw,nr,nv,nu,nn,nz,na,
c    +                         ne,npix,n1sm,n2sm,n3sm,n4sm,num_bad)
 

c               Set cloud mask quality bit flags
c               call set_confdnc(confdnc,testbits)

c               Set remaining QA bits 
c               call set_quality_A(nmtests,nbands,lsf,h_eco1,qa_bits)

c               Store element in line array 
c               call fill_bit_line(nc,nmtests,nbands,bad_value,bad_geo,
c    +                             snglnt,desert,testbits,qa_bits,
c    +                             bitarray,qa_bitarray)

c               call store_cm_scan(nc,confdnc,day,visusd,snglnt,snow,ice,
c    +                             coast,desert,land,process,
c    +                             scan_confdnc,scan_visusd,scan_snglnt,
c    +                             scan_snow,scan_ice,scan_coast,scan_desert,
c    +                             scan_land,scan_process,scan_day)

c 300        continue

c            Create the 250-m cloud mask corresponding to the current
c            1-km scan line.

c            do 200 nc = beg_ele, maxele

c              call get_250m_data(maxlin,lines_in_edge,klin,nc,maxele,
c    +                            v250_dat,scan_confdnc,scan_ice,
c    +                            scan_snow,scan_desert,scan_coast,
c    +                            scan_land,pxl250,cm1km,
c    +                            ice_250,snow_250,desert_250,coast_250,
c    +                            land_250,scan_snglnt,
c    +                            snglnt_250,scan_day,scan_process,
c    +                            scan_visusd,qa_bitarray,day_250,
c    +                            proc_250,visusd_250,qa_250)

c              call create_250m_cm(nc,pxl250,cm1km,ice_250,snow_250,
c    +                             desert_250,coast_250,land_250,
c    +                             bitarray,qa_bitarray,n250,ncr,ncd,
c    +                             day_250,proc_250,visusd_250,qa_250,
c    +                             snglnt_250)

c 200        continue

 5000        continue

c ...        write scanline of data into output HDF file
c            call mod35_out(s_pixels,maxlin,klin,lines_in_edge,
c    +                      modfil_mod35,bitarray,qa_bitarray)

c ...        keep track of number of lines processed
             outrec = outrec + 1

c ...        Check to see if processing first lines in edge of
c ...        context.  Don't want to rebuffer or change counter.
c            if (klin .gt. lines_in_edge) then
c ...           if not processing in single-line mode, re-buffer data 
c ...           in order to add the next scan line.
c               call reorder_array(contx_lat,contx_lon,contx_vzen,
c    +                             contx_szen,contx_rel,contx_pw,
c    +                             contx_eco,contx_topog,contx_ice,
c    +                             contx_snow,contx_sfctmp,contx_msl,
c    +                             contx_ugrd,contx_vgrd,contx_elev,
c    +                             s_type,s_pixels,n_nadir,maxele,indat,
c    +                             v250_dat,contx_sst)

c ...           Set context counter to 2 so that when processing 
c ...           next line will be in right order. 
                mlin = 2
c            end if 

         end if 

  100 continue
c ... ******************************************************************

c     Get stats needed for Product Specific Attributes (Metadata)
c     call stats_out(npix,num_bad,n1sm,n2sm,n3sm,n4sm,no,ns,
c    +               nd,nt,ng,ni,nl,nw,nr,nv,n250,ncr,
c    +               ncd,nu,nn,nz,na,ne,pcbad,pcn1s,pcn2s,pcn3s,pcn4s,
c    +               pcnnc,pcnns,pcnnd,pcnnt,pcnng,pcnni,pcnnl,pcnnw,
c    +               pcnnr,pcnnv,pcnncr,pcnncd,pcnu,pcnn,pcnz,pcna,
c    +               pcne)

      call Prepare_ECS_Metadata(pcbad,pcn1s,pcn2s,pcn3s,pcn4s,
     +                          pcnncd,pcnncr,pcnnd,pcnnt,pcnng,
     +                          pcnni,pcnnl,pcnnw,pcnns,pcnnv,pcnnr,
     +                          pcnnc,max_sol,min_sol,no_250,MET_groups)

c ... Set the debug runtime file metadata
      call Set_Runtime_Meta

c R. Frey 12/09/04 Added destripe flag copy to output file.
c     Copy UW "destripe" flag from L1b file to output MOD35 file.
c     Check for destriping attribute in L1b file.
      rcs_id = sffattr(modfil_L1B_1km(1), 'UW_DESTRIPE_LWIR')
      config_id = sffattr(modfil_L1B_1km(1), 'UW_DESTRIPE_CONFIG')
c     If attribute was found, write it to MOD35 output file.
      if (rcs_id .ne. FAIL .and. config_id .ne. FAIL) then
        rtn = sfrcatt(modfil_L1B_1km(1), rcs_id, rcs_string)
        rtn = sfrcatt(modfil_L1B_1km(1), config_id, config_string)
        rtn = sfscatt(modfil_mod35(1), 'UW_DESTRIPE_LWIR', DFNT_CHAR8,
     +                strlen(rcs_string), rcs_string)
        rtn = sfscatt(modfil_mod35(1), 'UW_DESTRIPE_CONFIG', DFNT_CHAR8,
     +                strlen(config_string), config_string)
      end if

      write(h_output,fmt='(/,20x,''MODIS CLOUD MASK FINISHED'')')

c ... Close files.
      call file_close(h_eco1,h_eco2,no_250,modfil_Geo,modfil_mod35,
     +                modfil_L1B_250m,modfil_L1B_1km,MET_Groups)

c ... Copy the 6 geolocation parameters into the output file
      call get_geo_copy(fname_geo,fname_mod35,max5km_lin,
     +                   max5km_ele)

c ... Delete temporary files
c rhucek 12/12/02: referencing new version of temp file cleanup function
      rtn = Delete_CommonTempFiles_Atmos(TempFile_LUNList, NumTempFiles)

c ... Cloud Mask completed successfully, return exit code of 0 to shell
      call exit(0)

      END
