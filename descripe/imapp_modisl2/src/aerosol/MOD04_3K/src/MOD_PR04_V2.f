       SUBROUTINE MOD_PR04_V2(modis_flag,idebug,modfil_mod04,RTN_NCEP,
     *                       SLOPE_MEAN_LAND,SLOPE_MEAN_OCEAN)

C----------------------------------------------------------------------
C!F77
C
C!DESCRIPTION:
C
C  The MODIS ocean aerosol product consists of aerosol optical thickness
C  and size parameters estimates derived on a 10x10 (1-km) pixel spatial
C  array.  The measured radiances in a wide spectral range (0.47-2.13
C  microns) are inverted assuming a bi-lognormal size distribution.
C  The volume and the mean particle size for each log-normal mode are
C  determined.  When fully developed, the aerosol ocean algorithm with
C  use the seven MODIS bands: 1, 2, 3, 4, 5, 6, and 7.
C
C!INPUT PARAMETERS:
C  modis_flag
C  idebug
C  modfil_mod04
C
C!OUTPUT PARAMETERS:  NONE
C
C!REVISION HISTORY:
C $Log: MOD_PR04_V2.f,v $
c 10/26/1999 fhliang
c moved the original #26 argument to the current #28 position in the
c calls to subroutine FILLVALUE_Ocean (L.635, L.692, and L.728).
c
c Revision 2.1  1996/05/28  15:34:05  vlin
c Updated to write out one scancube at a time
c
C
C!TEAM-UNIQUE HEADER:  NOT YET DEFINED
C
C!REFERENCES AND CREDITS:
C
C  WRITTEN BY:
C  SHANA MATTOO                E-MAIL:mattoo@climate.gsfc.nasa.gov
C  APPLIED RESEARCH CORPORATION          PHONE:  (301) 286-1025
C  NASA/GODDARD SPACE FLIGHT CENTER      FAX:    (301) 286-1759
C  CODE 913                             OFFICE: BLDG.22  RM.C109E
C  GREENBELT, MD 20771
C
C!DESIGN NOTES:
C
C Internals Variables:
C
C      W470_SYN(2*IGRIDX,2*IGRIDY)       Reflectance for wav=0.47um
C      W550_SYN(2*IGRIDX,2*IGRIDY)       Reflectance for wav=0.55um
C      W660_SYN(4*IGRIDX,4*IGRIDY)       Reflectance for wav=0.66um
C      W860_SYN(4*IGRIDX,4*IGRIDY)       Reflectance for wav=0.86um
C      W124_SYN(2*IGRIDX,2*IGRIDY)       Reflectance for wav=1.24um
C      W164_SYN(2*IGRIDX,2*IGRIDY)       Reflectance for wav=1.64um
C      W213_SYN(2*IGRIDX,2*IGRIDY)       Reflectance for wav=2.13um
C      W375_SYN(IGRIDX,IGRIDY)           Reflectance for wav=3.75um
C      CLDMSK_250(4*IGRIDX,4*IGRIDY) Cloud mask of 250 m resolution
C                                    (0=cloudy,1=clear)
C      CLDMSK_500(2*IGRIDX,2*IGRIDY) Cloud mask of 500 m resolution
C                                    (0=cloudy,1=clear)
C      CLDMSK_1km(2*IGRIDX,2*IGRIDY) Cloud mask of 1000 m resolution
C                                    (0=cloudy,1=clear)
C      MTHET0                        Measured solar zenith angle in deg.
C      MTHET                         Measured viewangle from ground in deg.
C      MPHI                          Measured azimuth  in deg.
C      Numdata                       Number of input data sets.
C
C Land and Ocean SDS_ARRAY for both land and ocean
C
C      SDS_Tau_Land_Ocean   Optical thickness at 0.55 micron for both land and ocean
C
C Ocean SDS_ARRAYS..........FOR ocean ONLY
C
C      SDS_ref_STD          Standard deviation of reflectances at 7 bands
C      SDSTAU_best          Optical thickness for best solution
C      SDSTAUS_best         Optical thickness contribution small particles for best solution
C      SDSTAUB_best         Optical thickness contribution large particles for best solution
C      SDSTAU_average       Optical thickness for best solution
C      SDSTAUS_average      Optical thickness contribution small particles for best solution
C      SDSTAUB_average      Optical thickness contribution large particles for best solution
C      SDS_Least_error      Least square error estimated
C      SDS_small_weighting  Small mode weighting factor
C      SDS_sol_INDX_small   Solution Number index small particles
C      SDS_sol_INDX_large   Solution Number index large particles
C      SDSASSY_best         Asymmetry_Factor at 7 bands for best solution
C      SDSASSY_average      Asymmetry_Factor at 7 bands for average solution
C      SDSBACK_best         Backscattering ratio at 7 bands of best solution
C      SDSBACK_average      Backscattering ratio at 7 bands of average solution
C      SDS_effrad           Effective_Radius at 0.55 micron of both solutions
C      SDS_AOT_model   Ratio of optical depth of small mode vs effective optical depth at 550
C      SDS_RefF_best        Normalized reflected_flux at 7 bands of best solution
C      SDS_RefF_average     Normalized reflected_flux at 7 bands of average solution
C      SDS_TranF_best       Normalized Transmitted_flux at 7 bands of best solution
C      SDS_TranF_average    Normalized Transmitted_flux at 7 bands of average solution
C      SDS_SCAT_ANGLE_OCEAN Scattering angle ocean
C      SDS_QCONTROL         Quality control SDS array
C      SDS_NUMPIXELS        Number of Pixels used for 0.55 micron
C      SDS_ccn              Cloud_Fraction in percentage
C      SDS_mass_conc        Mass concentration
C      SDS_angs_coeff1      Angstrom Exponent for 0.550 and 0.865 miron
C      SDS_angs_coeff2      Angstrom exponent for 0.865 and 2.130 micron
C
C LAND SDS_ARRAYS..........FOR LAND ONLY
C
C      SDS_QCONTROL_land    Quality control SDS array
C      SDS_Aerosol_Type     Index of Aerosol type
C      SDS_SCAT_ANGLE_land  Scattering Angle
C      SDS_angs_coeff_land  Angstrom exponent for 0.47 and 0.67 micron
C      SDS_CLDFRC_land      Cloud fraction (%)
C      SDS_dust_weighting   Dust aerosol weighting factor
C      SDS_est_uncer        Uncertainty of optical thickness at 0.47 and 0.66 micron
C      SDS_RefF_land        Normalized reflected flux at 0.47 and 0.66 microns
C      SDS_TranF_land       Normalized Transmitted flux at 0.47 and 0.66 microns
C      SDS_NUMPIXELS_land   Number of pixels with desired percentile
C      SDSTAU_corrected     Corrected optical thickness at 0.47 0.55 and 0.66 micron
C      SDS_ref_land         Mean reflectance at five bands
C      SDS_ref_STD_land     Standard deviation of reflectance at five bands
C      SDS_mass_conc_land   Mass concentration of aerosol
C EXTRA LAND SDS_ARRAYS..........FOR Statistics LAND ONLY
C
C      SDS_Mean_Reflectance_Land_All     
C      SDS_SDev_Reflectance_Land_All      
C      SDS_Path_Radiance_Land   
C      SDS_Critical_Reflectance_Land  
C      SDS_Error_Crit_Reflectance_Land      
C      SDS_Error_Path_Radiance_Land  
C      SDS_QFlag_Critical_Ref_Land
C      SDS_QFlag_Path_Radiance_Land
C
C!END
C-----------------------------------------------------------------------
C
      IMPLICIT  NONE
      SAVE

      include 'PGS_PC.f'
      include 'PGS_SMF.f'
      include 'PGS_MODIS_39500.f'
      include 'mapi.inc'
      include 'hdf.inc'
      include 'mod04.inc' 
c rhucek 08/09/02:  defined FORTRAN PARAMETER FUNCNAME
      CHARACTER*(*)   FUNCNAME
      PARAMETER      (FUNCNAME = 'MOD_PR04_V2')

c
c PGStoolkit,mapi etc
c
      CHARACTER*(PGSd_PC_FILE_PATH_MAX) L1B_FILE
c rhucek 08/09/02: added declaration of MSGBUF
      CHARACTER*(PGS_SMF_MAX_MSGBUF_SIZE) MSGBUF
      CHARACTER*(PGSd_MET_GROUP_NAME_L) groups(PGSd_MET_NUM_OF_GROUPS)
      CHARACTER*(MAX_ECS_NAME_L-1) HDFAttNms(PGSd_MET_NUM_OF_GROUPS)
      CHARACTER*(PGSD_PC_FILE_PATH_MAX) FN_L1B_1km,FN_L1B_500,FN_L1B_250

      INTEGER NumHandles,NUMSCAN_new
C
C HDF array names..........COMMON FOR BOTH OCEAN AND LAND
C
      BYTE       Cloud(Buf_cldmsk,ISWATH,ILINE),
     2           QA_Cloud(Buf_cldmsk_QA,ISWATH,ILINE)
      INTEGER *2 SDS_MTHET0(NUMCELLS),SDS_MTHET(NUMCELLS),
     &           SDS_MPHI0(NUMCELLS),SDS_MPHI(NUMCELLS),
     &           SDS_Scattering_Angle(NUMCELLS),
     &           SDS_Tau_Land_Ocean(NUMCELLS),
     &           SDS_Tau_Land_Ocean_img(NUMCELLS),
     &           SDS_ratio_small_Land_Ocean(NUMCELLS) ,
     &           SDS_Reflected_flux_Land_Ocean(NUMCELLS),
     &           SDS_Sea_Land_Flag(ISWATH,Tot_scan),
     &           SDS_Sea_Sunglint_Flag(ISWATH,Tot_scan),
     &           SDS_NCEP_Wspeed(ISWATH,Tot_scan), 
     &           SDS_Topo_Altitude(ISWATH,Tot_scan)
      REAL       SDSLAT(NUMCELLS),SDSLON(NUMCELLS),
     &           SDS_mass_conc_land(NUMCELLS),
     &           SLOPE_MEAN_LAND(3),SLOPE_MEAN_OCEAN(3)
      REAL*8     SDS_SST(NUMCELLS),SST
      BYTE       SDS_CldMskQA(NUMCELLS),QA_Temp

C
C LAND SDS_ARRAYS..........FOR LAND ONLY
C
      BYTE         SDS_QCONTROL_land(QA_LAND,NUMCELLS)
      BYTE         SDS_QCONTROL_CritRef_land(QA_LAND,NUMCELLS)
      BYTE         SDS_QCONTROL_ocean(QA_Ocean,NUMCELLS)
      INTEGER *2   SDS_Aerosol_Type(NUMCELLS),
     &             SDS_SCAT_ANGLE_land(NUMCELLS),
     &             SDS_angs_coeff_land(NUMCELLS),
     &             SDS_CLDFRC_land(NUMCELLS),
     &             SDS_dust_weighting(NUMCELLS),
     &             SDS_NUMPIXELS_land(NUMCELLS,Land_Sol1),
     &             SDSTAU_corrected(NUMCELLS,Land_Sol3),
     &             SDS_ref_land(NUMCELLS,Band_land),
     &             SDS_ref_STD_land(NUMCELLS,Band_land),
C 9/2005 ( two new SDS's)......     
     &            SDS_Surface_Reflectance_Land(NUMCELLS,Land_Sol3),
     &            SDS_Fitting_Error_Land(NUMCELLS),
     &             SDSTAU_corrected_213(NUMCELLS),
     &             SDSTAU_small_land(NUMCELLS,Land_Sol4)

C
C EXTRA LAND SDS_ARRAYS..........FOR LAND Statistics ONLY
C
      INTEGER *2  SDS_Mean_Reflectance_Land_All(NUMCELLS,Land_Sol3), 
     &SDS_SDev_Reflectance_Land_All(NUMCELLS,Land_Sol3),     
     &SDS_Path_Radiance_Land(NUMCELLS,Land_Sol1),   
     &  SDS_Critical_Reflectance_Land(NUMCELLS,Land_Sol1),  
     &  SDS_Error_Crit_Reflectance_Land(NUMCELLS,Land_Sol1),      
     &  SDS_Error_Path_Radiance_Land(NUMCELLS,Land_Sol1),  
     &  SDS_QFlag_Critical_Ref_Land(NUMCELLS,Land_Sol1),
     &  SDS_QFlag_Path_Radiance_Land(NUMCELLS,Land_Sol1),
     &  SDS_Quality_Land_Ocean(NUMCELLS)

C
C Obsolete (02/2006) Land SDS Arrays
C
      INTEGER *2   
     &            SDS_est_uncer(NUMCELLS,Land_Sol1),
     &            SDS_RefF_land(NUMCELLS,Land_Sol2),
     &            SDS_TranF_land(NUMCELLS,Land_Sol1)

C
C OCEAN SDS_ARRAYS.........FOR OCEAN ONLY
C
      REAL         SDS_ccn(NUMCELLS,NUM_solutions),
     &             SDS_mass_conc(NUMCELLS,NUM_solutions)
      BYTE         SDS_QCONTROL(QA_OCEAN,NUMCELLS)
      INTEGER *2   SDS_ref(NUMCELLS,NWAV_S),
     &             SDS_ref_STD(NUMCELLS,NWAV_S),
     &             SDSTAU_best(NUMCELLS,NWAV_S),
     &             SDSTAU_average(NUMCELLS,NWAV_S),
     &             SDSTAUB_best(NUMCELLS,NWAV_S),
     &             SDSTAUB_average(NUMCELLS,NWAV_S),
     &             SDSTAUS_best(NUMCELLS,NWAV_S),
     &             SDSTAUS_average(NUMCELLS,NWAV_S),
     &             SDSASSY_best(NUMCELLS,NWAV_S),
     &             SDSASSY_average(NUMCELLS,NWAV_S),
     &             SDSBACK_best(NUMCELLS,NWAV_S),
     &             SDSBACK_average(NUMCELLS,NWAV_S),
     &             SDS_RefF_best(NUMCELLS,NWAV_S)
      INTEGER *2   SDS_RefF_average(NUMCELLS,NWAV_S),
     &             SDS_TranF_best(NUMCELLS,NWAV_S),
     &             SDS_TranF_average(NUMCELLS,NWAV_S),
     &             SDS_small_weighting(NUMCELLS,NUM_solutions),
     &             SDS_correc_small_weighting(NUMCELLS),
     &             SDS_Least_error(NUMCELLS,NUM_solutions),
     &             SDS_effrad(NUMCELLS,NUM_solutions),
     &             SDS_sol_INDX_small(NUMCELLS,NUM_solutions),
     &             SDS_sol_INDX_large(NUMCELLS,NUM_solutions),
     &             SDS_angs_coeff1(NUMCELLS,NUM_solutions),
     &             SDS_angs_coeff2(NUMCELLS,NUM_solutions),
     &             SDS_AOT_model(NUMCELLS,num_model_index),
     &             SDS_NUMPIXELS(NUMCELLS),SDS_CLDFRC(NUMCELLS),
     &             SDS_SCAT_ANGLE_OCEAN(NUMCELLS),
     &             SDS_CLDFRC_ocean(NUMCELLS)
C
C Input and Output file names for PGSTK etc.......
C
      INTEGER iik,RTN_MOD05,RTN_MOD07,RTN_NCEP,RTN_DAO,RTN_TOMS,
     1        RTN_TOVS,QCONTROL_land_wav1,QCONTROL_land_wav2
      INTEGER MODFIL_L1B_1KM(MODFILLEN),MODFIL_L1B_500(MODFILLEN),
     1        MODFIL_L1B_250(MODFILLEN),MODFIL_GEO(MODFILLEN),
     2        MODFIL_CLDMSK(MODFILLEN),MODFIL_MOD04(MODFILLEN),
     2        MODFIL_MOD05(MODFILLEN),MODFIL_MOD07(MODFILLEN),
     3        MODFIL_ANC(MODFILLEN),MODFIL_MOD04S(MODFILLEN)
C
C Position of image in the processing of 10x10 km box
C
      INTEGER START_500,END_500,START_250,END_250,START_1KM,END_1KM
      INTEGER Iscan,NUMSCAN,Buf_Size1,Buf_Size2,Data_Size(2),i,j
      INTEGER IMONTH(maxnum_scans),IDATA,Water,Land,cloud_num,
     1        cloud_num_land,Glint,Snowice,NUMSQ,ij,ik,Qcontrol_special,
     2         Land_CLDMSK_forfraction(ISWATH,ILINE),
     3         Qcontrol_special_land,Quality_to_pass,xsize,ysize
C
C
C Define arrays of Geoloation, cloud mask and solar and satellite angles
C
      REAL MTHET0,MTHET,MPHI0,MPHI,MDPHI,MSCATT,MHGHT,Lat_center,Lon_center
      REAL sfctmp,ugrd,vgrd,pwat,ozone
      REAL Lat(ISWATH,ILINE),Lon(ISWATH,ILINE),
     2     SatZen(ISWATH,ILINE),SolZen(ISWATH,ILINE),
     3     SatAz(ISWATH,ILINE),SolAz(ISWATH,ILINE),
     4     RelAz(ISWATH,ILINE),Height(ISWATH,ILINE),G_factor
C      INTEGER L_S_Flag(ISWATH,ILINE),Sunglint_Flag(ISWATH,ILINE)
C
C CLoud mask
C Old
C      INTEGER Cloud(4*ISWATH,4*ILINE)
C      INTEGER CLDMSK_syn_land(4*ISWATH,4*ILINE)
       INTEGER ierror, Quality_Flag_for_retr,New_QA_Flag_Land(19)
        REAL        WSPEED 
C New
       INTEGER CldMsk_250(4*ISWATH,4*ILINE),
     1         CldMsk_500(2*ISWATH,2*ILINE),
     2         CldMsk_1km(ISWATH,ILINE),
     2         DET_Flag(ISWATH,ILINE),
     3         UFQ_Flag(ISWATH,ILINE),
     4         DayNight_Flag(ISWATH,ILINE),
     5         SunGlint_Flag(ISWATH,ILINE),
     6         SnowIce_Flag(ISWATH,ILINE),
     6         SnowMsk_Ratio(ISWATH,ILINE),
     6         SnowMsk_500m(2*ISWATH,2*ILINE) 
          INTEGER       LandSea_Flag(ISWATH,ILINE),
     8         Non_CloudOb_Flag(ISWATH,ILINE),
     9         Thin_CirNIR_Flag(ISWATH,ILINE),
     1         Shadow_Flag(ISWATH,ILINE),
     2         Thin_CirIR_Flag(ISWATH,ILINE),
     3         Cloud_SimpIR_Flag(ISWATH,ILINE),
     4         High_Cloud_Flag(ISWATH,ILINE),
     5         Cloud_IRTemp_Flag(ISWATH,ILINE),
     6         Cloud_3p75_11_Flag(ISWATH,ILINE),
     7         Cloud_VisRat_Flag(ISWATH,ILINE),
     8         Cloud_SpatVar_Flag(ISWATH,ILINE),
     9         High_Cloud_Flag_500(2*ISWATH,2*ILINE),
     9         NROWS(ISWATH*2),index_wspeed

C
C L1B reflectance and radiance data (temporary)
C

      REAL Refl_1 (4*ISWATH,4*ILINE),Refl_2 (4*ISWATH,4*ILINE),
     *     Refl_3 (2*ISWATH,2*ILINE),Refl_4 (2*ISWATH,2*ILINE),
     *     Refl_5 (2*ISWATH,2*ILINE),Refl_6 (2*ISWATH,2*ILINE),
     *     Refl_7 (2*ISWATH,2*ILINE),
     *     Refl_9 (ISWATH,ILINE),Refl_12(ISWATH,ILINE),
     *     Refl_13(ISWATH,ILINE),Refl_16(ISWATH,ILINE),
     *     Refl_26(ISWATH,ILINE),
     *     Rad_20(ISWATH,ILINE),Rad_31(ISWATH,ILINE),
     *     Rad_32(ISWATH,ILINE)

      BYTE Unc_1 (4*ISWATH,4*ILINE),Unc_2 (4*ISWATH,4*ILINE),
     *     Unc_3 (2*ISWATH,2*ILINE),Unc_4 (2*ISWATH,2*ILINE),
     *     Unc_5 (2*ISWATH,2*ILINE),Unc_6 (2*ISWATH,2*ILINE),
     *     Unc_7 (2*ISWATH,2*ILINE),
     *     Unc_9 (ISWATH,ILINE),Unc_12(ISWATH,ILINE),
     *     Unc_13(ISWATH,ILINE),Unc_16(ISWATH,ILINE),
     *     Unc_26(ISWATH,ILINE),
     *     Unc_20(ISWATH,ILINE),Unc_31(ISWATH,ILINE),
     *     Unc_32(ISWATH,ILINE)

      BYTE Vflag_1(4*ISWATH,4*ILINE),
     &     Vflag_2(4*ISWATH,4*ILINE),
     &     Vflag_3(2*ISWATH,2*ILINE),
     &     Vflag_4(2*ISWATH,2*ILINE),
     &     Vflag_5(2*ISWATH,2*ILINE),
     &     Vflag_6(2*ISWATH,2*ILINE),
     &     Vflag_7(2*ISWATH,2*ILINE),
     &     Vflag_9(ISWATH,ILINE),
     &     Vflag_12(ISWATH,ILINE),
     &     Vflag_13(ISWATH,ILINE),
     &     Vflag_16(ISWATH,ILINE),
     &     Vflag_26(ISWATH,ILINE),
     &     Vflag_20(ISWATH,ILINE),
     &     Vflag_31(ISWATH,ILINE),
     &     Vflag_32(ISWATH,ILINE)

      BYTE Buf_Un(ISWATH,ILINE),
     &     Buf_Sa_1km(ISWATH,ILINE),
     &     Buf_Sa_500(2*ISWATH,2*ILINE),
     &     Buf_Sa_250(4*ISWATH,4*ILINE)

      REAL W659_SYN(4*ISWATH,4*ILINE),W865_SYN(4*ISWATH,4*ILINE),
     *     W470_SYN(2*ISWATH,2*ILINE),W550_SYN(2*ISWATH,2*ILINE),
     *     W124_SYN(2*ISWATH,2*ILINE),W164_SYN(2*ISWATH,2*ILINE),
     *     W213_SYN(2*ISWATH,2*ILINE),
     *     W443o_SYN(ISWATH,ILINE),W551o_SYN(ISWATH,ILINE),
     *     W667o_SYN(ISWATH,ILINE),W869o_SYN(ISWATH,ILINE),
     *     W395_SYN(ISWATH,ILINE),W138_SYN(ISWATH,ILINE),
     *     W1100_SYN(ISWATH,ILINE),W1100_SYN_500m(2*ISWATH,2*ILINE),
     *     W1200_SYN(ISWATH,ILINE),AA

C
C Ancillary Data of Total precipitable water and ozone 
C
      REAL Total_H2O(ISWATH,ILINE),Total_O3(ISWATH,ILINE)
C
C QA arrays (see MODIS Atmosphere QA Plan)
C
      INTEGER QA_Flag_Land(19),QA_Flag_Ocean(12),Success_Ret_Land,
     *        Fail_Ret_Land,NO_Ret_Land,Total_Grids,
     *        Success_Ret_Ocean,Fail_Ret_Ocean,NO_Ret_Ocean,
     *       quality_cirrus(ISWATH,ILINE),Ret_Quality_cirrus,
     *       Quality_flag_forJoint
C
C MetaData arrays
c
        INTEGER No_PSA,JX,JY
        PARAMETER(No_PSA=15)
        REAL QA_Metadata_MOD04(No_PSA)
C
C Arrays for miscellaneous use
C
      INTEGER Set_Counter_Land,Set_Counter_Ocean,GetModisDat_MOD04
      INTEGER maskoption,num_resol,rtn,Set_Counter_Ocean_cloud
      INTEGER create_mod04,idebug,Indx_wspeed1,Indx_wspeed2
      CHARACTER*4 choice_flag
      LOGICAL error_flag,modis_flag
      Real    Glint_angle
      integer CldMsk_500_Ocean(2*ISWATH,2*ILINE+1)
      integer index_wave,quality_land,ix,iy,coastal_Pix
      integer savecldmask(2*ISWATH,2*ILINE),Pure_Land
c    New cldmask
      REAL RMED(ISWATH*2 ),RMEDSQ(ISWATH*2 )
      REAL RMED_1km(ISWATH),RMEDSQ_1km(ISWATH),Wind(Lut_indx) 
      
c For scan Time      
       real*8  xx(maxnum_Scans),yy(maxnum_Scans),xbar,y1
       integer Exrea_scan,set_counter_for_Gread,snow_present_Flag
       real   Multi_factor(NWAV_S+3)


      integer sfstart, sfscatt 
      integer sfend, pgs_pc_getconfigdata
      
      integer sd_id
      character*4 pcf_satid
      character*26 doi_char
      
      integer LUN_Sat_Instrument
      parameter (LUN_Sat_Instrument = 800510)

      integer file_version

      character*(PGSd_PC_FILE_PATH_MAX) HDF_FILENAME

      integer pgs_io_gen_openf,pgs_pc_getreference
      character*1024 msg

C
C Valid choice_flag: 'MA', 'ML', 'MS', 'NS', and 'NL'
C (MA=Mode of aerosol over land and ocean)
C (ML=Mode of aerosol over land)
C (MS=Mode of aerosol over ocean)
C (NS= ??)
C (NL= ??)
C

C rhucek 08/09/02:  Added "Begin..." message 
      MSGBUF = char(10) 
     1// '------------------------------------------------------------'
     2// char(10)
     3// 'Begin Land and Ocean Aerosol Retrievals' 
     4// char(10)
     5// '------------------------------------------------------------'

      Call MODIS_SMF_SETDYNAMICMSG(MODIS_M_GENERIC,MSGBUF,FUNCNAME)


      choice_flag = 'MA'
      error_flag = .false.
C
C Create MOD_PR04 product file with HDF AND HDF-EOS libraries
C temporarily disabled
C
c      rtn = create_mod04() 
      RTN_TOMS=1
      RTN_TOVS=1
      RTN_DAO=1
C      RTN_NCEP=0
      RTN_MOD05=1
      RTN_MOD07=1
      Success_Ret_Land=0
      Fail_Ret_Land=0
      NO_Ret_Land=0
      Success_Ret_Ocean=0
      Fail_Ret_Ocean=0
      NO_Ret_Ocean=0
      set_counter_for_Gread =0
      
      IF(rtn.ne.SUCCEED) CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,
     &    'create_mod04 failed','MOD_PR04_V2.f')
C
C Open input and output files
C
      CALL FILE_OPEN(choice_flag,error_flag,
     &     MODFIL_L1B_1km,MODFIL_L1B_500,MODFIL_L1B_250,
     &     FN_L1B_1km,FN_L1B_500,FN_L1B_250,
     &     MODFIL_Geo,MODFIL_Cldmsk,MODFIL_MOD04,MODFIL_MOD05,
     &     MODFIL_MOD07,MODFIL_ANC,
     &     HANDLE_LUT466,HANDLE_LUT553,HANDLE_LUT644,HANDLE_LUT213,
     &     HANDLE_LUTMAP,HANDLE_S,HANDLE_L,
     &     handle_gascoeff,HANDLE_OUTSCI,HANDLE_QC)
C
C Call to Get_Metatable for Dead_detectors,Scan_Type, EV_Frames
C
       CALL Get_Metatable(MODFIL_L1B_1km,NUMSCAN,Scan_Type,
     &                    EV_Frames,Scan_Start_Time,IMONTH)

        Buf_Size1 = ISWATH
        Buf_Size2 = ILINE
C
C NUMSCAN gives number of scans of a granule
C
c                  NUMSCAN=1
c      Set_Counter_Ocean is set to read the table once as first entry
C      into ocean algorithm
       Set_Counter_Ocean=0
       Set_Counter_Land=0
! Extrapolationg for Scan time       
       
            NUMSCAN_new=(NUMSCAN*10)/Iline  
              DO   Iscan = 1,NUMSCAN
                XX(iscan)=Iscan
                yy(Iscan)=Scan_Start_Time(Iscan) 
              enddo
               DO  Exrea_scan = 1,NUMSCAN_new
                  Xbar= Exrea_scan
             CALL INTERPOlation_R_eight(NUMSCAN,Xbar,XX,YY,Y1,1)
              New_Scan_Start_Time(Exrea_scan)=Y1 
             Enddo
             
            
 
              DO 9999 Iscan = 1,NUMSCAN_new
                 
 
              

c       Set_Counter_Ocean_cloud is set to make cloud mask once for
C       each scan inside ocean algorithm

        Set_Counter_Ocean_cloud=0
c        SST=Scan_Start_Time(Iscan)
        SST=New_Scan_Start_Time(Iscan)
C  Check  NE Night is because we have sometimes both and we need to process those granues so leave it NE night
           IF(Scan_Type(Iscan).NE.'N') THEN    
C
C
C Replace with new subroutine to read L1B data
C
         RTN=GetModisDat_MOD04(NUMSQ,
     &      Modfil_L1B_1km,Modfil_L1B_500,Modfil_L1B_250,
     &      Modfil_Geo,Modfil_CldMsk,
     &      FN_L1B_1km,FN_L1B_500,FN_L1B_250,
     &      Iscan,Buf_Size1,Buf_Size2,Data_Size,
     &      Lat,Lon,SatZen,SatAz,SolZen,SolAz,RelAz,Height,
     &      Refl_1,Refl_2,Refl_3,Refl_4,Refl_5,Refl_6,Refl_7,
     &      Refl_9,Refl_12,Refl_13,Refl_16,Refl_26,
     &      Rad_20,Rad_31,Rad_32,
     &      Unc_1,Unc_2,Unc_3,Unc_4,Unc_5,Unc_6,UNc_7,Unc_9,
     &      Unc_12,Unc_13,Unc_16,Unc_26,Unc_20,Unc_31,Unc_32,
     &      Vflag_1,Vflag_2,Vflag_3,Vflag_4,Vflag_5,Vflag_6,Vflag_7,
     &      Vflag_9,Vflag_12,Vflag_13,Vflag_16,Vflag_26,Vflag_20,
     &      Vflag_31,Vflag_32,Buf_Un,Buf_Sa_1km,Buf_Sa_500,Buf_Sa_250,
     &      CldMsk_250,CldMsk_500,CldMsk_1km,DET_Flag,UFQ_Flag,
     &      DayNight_Flag,SunGlint_Flag,SnowIce_Flag,
     &      LandSea_Flag,Non_CloudOb_Flag,Thin_CirNIR_Flag,
     &      Shadow_Flag,Thin_CirIR_Flag,Cloud_SimpIR_Flag,
     &      High_Cloud_Flag,Cloud_IRTemp_Flag,Cloud_3p75_11_Flag,
     &      Cloud_VisRat_Flag,Cloud_SpatVar_Flag,Cloud,QA_CLOUD)
              
      
         IF (rtn.ne.0) CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     &                'Call to GetModisDat_MOD04 failed',
     &                'MODIS_AEROSOL_LAND_OCEAN.f')
  
C Vanderlie  Cloud mask...........
      CALL CldMsk_Land(Data_Size,ISWATH,ILINE,Refl_3,Refl_5,
     &  Refl_26,CldMsk_250,CldMsk_500,CldMsk_1km, RMED,RMEDSQ, 
     &RMED_1km,RMEDSQ_1km,iscan,nrows,quality_cirrus,
     &Land_CLDMSK_forfraction)
     
CC ...    Vanderlie  Cloud mask  OCEAN >>>>>>>>>>>>       
               Call cldMsk_Ocean(Refl_3 ,Refl_4,Refl_1,
     &         Refl_26,Refl_5,CldMsk_500_Ocean,High_Cloud_Flag,
     &         data_size,savecldmask) 
         
      DO 9000 IDATA = 1,NUMSQ

      SDS_Tau_Land_Ocean_img(IDATA)=-9999
      SDS_Tau_Land_Ocean(IDATA)=-9999
      SDS_Quality_Land_Ocean(IDATA)=-9999
      SDS_Sea_Sunglint_Flag(Idata,iscan)=-9999
      SDS_NCEP_Wspeed(Idata,iscan) =-9999
      SDS_Topo_Altitude(Idata,iscan) =-9999
      Quality_to_pass=-9999
      Quality_flag_forJoint=-9999
      SDS_Sea_Land_Flag(Idata,iscan)=-9999
C
C Call to SET_INDEX sets the indexing for Modis Channels
C

        CALL SET_INDEX(START_500,END_500,START_250,END_250,
     &         START_1KM,END_1KM,IDATA)
         
C
C Call to geoLoc_angle computes geolocation and geometry in
C       center of 10* 10 bin.
C

        CALL GEOLOC_ANGLE(LAT,LON,SatZen,SatAz,SolZen,SolAz,RelAz,
     &       Height,MTHET0,MTHET,MPHI0,MPHI,MDPHI,MSCATT,MHGHT,IDATA,
     &       START_500,END_500,START_250,END_250,START_1KM,END_1KM,
     &       Lat_center,Lon_center,iscan)

  

C
C Call to GetAncData_PGE04 reads NCEP data modified by UW group 
C tailored to PGE04 
C At this point it is only reading the data
C

       IF(RTN_NCEP.EQ.0) THEN
         CALL GetAncData_PGE04(Lat_center,Lon_center,sfctmp,ugrd,vgrd,pwat,ozone)
          
       ENDIF 
                 set_counter_for_Gread=set_counter_for_Gread+1
         
          call compute_Gascorrection(pwat,ozone,
     +   MTHET0,MTHET,set_counter_for_Gread,Multi_factor,RTN_NCEP,
     +   handle_gascoeff,QA_Flag_Land,QA_Flag_Ocean)  
       
C     Bilinearly interpolate NCEP wind speed to exact
C     Lat/Lon on Lat_center, lon_center - Leigh Munchak, 03/27/2013 
C
      CALL interpol_wspeed_spatial(Lat_center,Lon_center,ugrd,vgrd) 

      CALL indx_wspeed(ugrd,vgrd,Indx_wspeed1,Indx_wspeed2,WSPEED,Wind)
         

       CALL TRANS_2WAY(QA_Flag_Land,QA_Flag_Ocean,START_250,END_250,
     *       START_500,END_500,START_1KM,END_1KM,RTN_MOD05,
     *       RTN_MOD07,RTN_NCEP,RTN_DAO,RTN_TOMS,RTN_TOVS,
     *       ISWATH,ILINE,MTHET,MTHET0,Refl_1,Refl_2,Refl_3,Refl_4,
     *       Refl_5,Refl_6,Refl_7,Refl_9,Refl_12,Refl_13,Refl_16,
     *       Refl_26,Rad_20,Rad_31,Rad_32,LandSea_Flag,CldMsk_1km,
     *       SnowMsk_500m,SnowMsk_Ratio,
     *       pwat,ozone,SLOPE_MEAN_LAND,SLOPE_MEAN_OCEAN,W659_SYN,
     *       W865_SYN,W470_SYN,W550_SYN,W124_SYN,W164_SYN,W213_SYN,
     *       W395_SYN,W138_SYN,W1100_SYN,W1100_SYN_500m,
     *       W1200_SYN,G_factor,Multi_factor)
               

C
C Call to determine counts of land, water, cloudy, sun-glint and
C snow/ice pixels
C
          
  
       
       CALL LAND_WATER_CLOUD(DET_Flag,CldMsk_1km,CldMsk_500,CldMsk_250,
     &      High_Cloud_Flag,High_Cloud_Flag_500,LandSea_Flag,
     &      SunGlint_Flag,SnowIce_Flag,SnowMsk_Ratio,Shadow_Flag,cloud_num,
     &      cloud_num_land,Water,Land,Glint,Snowice,START_1KM,END_1KM,
     &      QA_Flag_Land,QA_Flag_Ocean,QA_Temp,snow_present_Flag, 
     *       Quality_cirrus, Ret_Quality_cirrus,Land_CLDMSK_forfraction,
     *       coastal_Pix,Pure_Land)
              
     
C
C Call to POPULATE_common_arrays sets the SDS for HDF write for
C        Geolocation and angles for 10*10 bin size.

       CALL POPULATE_COMMON_ARRAYS(IDATA,Lat_center,Lon_center,
     &         SST,MTHET0,MTHET,MPHI0,MPHI,MSCATT,QA_Temp,SDSLAT,
     &         SDSLON,SDS_SST,SDS_MTHET0,SDS_MTHET,SDS_MPHI0,SDS_MPHI,
     &         SDS_Scattering_Angle,SDS_CldMskQA)

C
C Restore cloud fraction in QA_Flag_Land to pass to Process_Land
C

  

      Call COMPUTE_GLINTANGLE(MTHET0,MTHET,MDPHI,
     *             GLINT_ANGLE,QA_Flag_Ocean)
      SDS_Sea_Sunglint_Flag(Idata,iscan) =(GLINT_ANGLE*SCALE2+OFFSET2)
 
C If all pixels of 10*10 bin are water, then processing ocean
C algorithm, else  processing land algorithm; if the number
C of cloudy pixels >90 then reject any processing
C
C
           Qcontrol_special_land=0
C 
   
         cloud_num=0
         Qcontrol=-1
         IF(water.ge.(Iline*iline) .or. land .gt. 0) THEN
           IF(WATER .GE. (Iline*iline)) THEN 
           SDS_Sea_Land_Flag(Idata,iscan)=0
             Set_Counter_Ocean=Set_Counter_Ocean+1
            Set_Counter_Ocean_cloud=Set_Counter_Ocean_cloud+1
            CALL PROCESS_ocean(HANDLE_S,HANDLE_L,
     &        ISCAN,IDATA,NUMSQ,MTHET0,MTHET,MDPHI,START_500,
     &        END_500,START_250,END_250,START_1KM,END_1KM,
     &        W659_SYN,W865_SYN,W470_SYN,W550_SYN,W124_SYN,
     &        W164_SYN,W213_SYN,W551o_SYN,W667o_SYN,W869o_SYN,
     &        Sunglint_Flag,CldMsk_500,Set_Counter_Ocean,
     &        SDSTAU_best,SDSTAUS_best,SDSTAUB_best,
     &        SDSTAU_average,SDSTAUS_average,SDSTAUB_average,
     &        SDS_Least_error,SDS_small_weighting,SDS_sol_INDX_small,
     &        SDS_sol_INDX_large,SDS_ref,SDS_ref_STD,SDSASSY_best,
     &        SDSASSY_average,SDSBACK_best,SDSBACK_average,SDS_effrad,
     &        SDS_RefF_best,SDS_ccn,SDS_mass_conc,
     &        SDS_RefF_average,SDS_TranF_best,SDS_TranF_average,
     &        SDS_QCONTROL_ocean,SDS_NUMPIXELS,SDS_SCAT_ANGLE_OCEAN,
     &        SDS_AOT_model,SDS_CLDFRC_ocean,
     &        Set_Counter_Ocean_cloud,QA_Flag_ocean,GLINT_ANGLE,
     &        SDS_angs_coeff1,SDS_angs_coeff2,SDS_Tau_Land_Ocean,Refl_4,
     &  Qcontrol_special,SDS_correc_small_weighting,SDS_Tau_Land_Ocean_img,
     &  High_Cloud_Flag_500,W138_SYN,Refl_3,Refl_1,data_size,CldMsk_500_Ocean,
     &  savecldmask,Quality_to_pass, Indx_wspeed1,Indx_wspeed2, WSPEED, Wind)
          
C

              if((SDSTAU_average(IDATA,2)/(SCALE3+OFFSET3)).GT.-0.01) then
                 SDS_Quality_Land_Ocean(IDATA)=Quality_to_pass
                SDS_NCEP_Wspeed(Idata,iscan) = WSPEED*SCALE2+OFFSET2
               Endif
C Call Subroutine to get total number of successful and failed retreivals
C
         CALL Total_retrieval_ocean(Success_Ret_Ocean,
     *        Fail_Ret_Ocean)
C
C1/2002   Fill small weighting for ocean from average value
          if((SDSTAU_average(IDATA,2)/(SCALE3+OFFSET3)).GT.-0.01) then
                 aa= ((SDSTAUS_average(IDATA,2)/(SCALE3+OFFSET3))
     *                       /(SDSTAU_average(IDATA,2)/SCALE3+OFFSET3))
                SDS_ratio_small_Land_Ocean(IDATA) = (aa*SCALE3+OFFSET3)
                else
               SDS_ratio_small_Land_Ocean(IDATA)=SDSTAU_average(IDATA,2)
               ENDIF
           SDS_Reflected_flux_Land_Ocean(IDATA)= SDS_RefF_average(IDATA,2)
          CALL POPULATE_TAU_LAND_OCEAN(IDATA,SDS_Tau_Land_Ocean_img,
     *   SDSTAU_average,NUMCELLS,NWAV_S)
     
         if((SDSTAUS_average(IDATA,2)/(SCALE3+OFFSET3)).GT.-0.01 .and. 
     *     Quality_to_pass .gt.0 .and. Quality_to_pass .le.3) 
     *    CALL POPULATE_TAU_LAND_OCEAN(IDATA,SDS_Tau_Land_Ocean,
     *    SDSTAU_average,NUMCELLS,NWAV_S)
 
     
C
C Filled with Fill_Value for land
C
             SDS_CLDFRC_land(IDATA)=-99
            CALL FILLVALUE_LAND(IDATA,SDS_Tau_Land_Ocean_img,
     *        SDS_Aerosol_Type,SDSTAU_corrected_213,
     *        SDS_SCAT_ANGLE_land,SDS_mass_conc_land,
     *        SDS_angs_coeff_land,SDS_CLDFRC_land,
     *        SDS_dust_weighting,
     *        SDS_est_uncer,SDS_RefF_land,SDS_TranF_land,
     *        SDS_NUMPIXELS_land,SDSTAU_corrected,SDS_ref_land,
     *        SDS_ref_STD_land,SDS_QCONTROL_land,SDSTAU_small_land, 
     *        SDS_Surface_Reflectance_Land,
     *        SDS_Fitting_Error_Land,Qcontrol_special_land)
                   NO_Ret_Land=NO_Ret_Land+1
              index_wave=3
            QCONTROL_land_wav1=0
            QCONTROL_land_wav2=0
      
   
C ELSE FOR LAND>>>>>>>>>>>>>>>>>

          ELSE
            if( Land .GT.0) then
 
C
C IF not all pixels detected as water pixels, proceed land algorithm
C
            quality_land=1
            
            Set_Counter_Land=Set_Counter_Land+1
!           Degrade the quality of Land Ret if water           
              if( water .ge.8) quality_land=0  
! For Sea land flag, set land +intermittent water(Pureland), 50% coastal as coastal and less than 50% back to be called land            
              if(Pure_land .gt.0)SDS_Sea_Land_Flag(Idata,iscan)=1
              if(coastal_Pix.le.5)SDS_Sea_Land_Flag(Idata,iscan)=1
              if(coastal_Pix.gt.5)SDS_Sea_Land_Flag(Idata,iscan)=2
              
            CALL PROCESS_Land(HANDLE_LUT466,HANDLE_LUT553,HANDLE_LUT644,
     &  HANDLE_LUT213,HANDLE_LUTMAP,IMONTH(ISCAN),ISCAN,IDATA,NUMSQ,MTHET0,
     &  MTHET,MDPHI,MHGHT,Lat_center,Lon_center,START_500,END_500, 
     &  START_250,END_250,START_1KM,END_1KM,W470_syn,W550_SYN,
     &  W659_syn,W865_syn,W124_SYN,W164_SYN,W213_syn,
     &  CldMsk_250,Set_Counter_Land,QA_Flag_Land,Success_Ret_Land,
     &  Fail_Ret_Land,SDSLAT,SDSLON,SDS_MTHET0,SDS_MTHET,SDS_MPHI,
     &  SDS_Tau_Land_Ocean,CLDMSK_500,SDS_Tau_Land_Ocean_img,
     &  SDS_Aerosol_Type,SDS_SCAT_ANGLE_land,SDS_mass_conc_land,
     &  SDS_angs_coeff_land,SDS_CLDFRC_land,SDS_dust_weighting,
     &  SDS_est_uncer,SDS_RefF_land,
     &  SDS_TranF_land,SDS_NUMPIXELS_land,SDSTAU_corrected,
     &  SDS_ref_land,SDS_ref_STD_land,SDS_QCONTROL_land,
     &  SDS_Mean_Reflectance_Land_All, SDS_SDev_Reflectance_Land_All,     
     &  SDS_Path_Radiance_Land, SDS_Critical_Reflectance_Land,  
     &  SDS_Error_Crit_Reflectance_Land,SDS_Error_Path_Radiance_Land,  
     &  SDS_QFlag_Critical_Ref_Land,SDS_QFlag_Path_Radiance_Land,SDS_QCONTROL_CritRef_land,
     &  G_factor,quality_land,Ret_Quality_cirrus,cloud_num_land,
     &  SDS_Surface_Reflectance_Land,SDS_Fitting_Error_Land,Qcontrol_special_land,
     &  SDSTAU_corrected_213,Quality_flag_forJoint,SDSTAU_small_land,New_QA_Flag_Land,
     &  snow_present_Flag,coastal_Pix)
          
C1/2002   Fill small weighting for land 

            SDS_ratio_small_Land_Ocean(IDATA) =SDS_dust_weighting(idata) 
            SDS_Reflected_flux_Land_Ocean(IDATA)= SDS_RefF_land(IDATA,2) 

            CALL POPULATE_TAU_LAND_OCEAN(IDATA,SDS_Tau_Land_Ocean_img,
     *     SDSTAU_corrected,NUMCELLS,Land_Sol3) 
     
C If retreival populate Quality flag                
        if((SDSTAU_corrected(Idata,2)/SCALE3+OFFSET3) .ge.-0.05) then
            SDS_Quality_Land_Ocean(IDATA)=Quality_flag_forJoint
            SDS_topo_Altitude(Idata,iscan) =MHGHT*SCALE2+OFFSET2
             Endif
        if((SDSTAU_corrected(Idata,2)/SCALE3+OFFSET3) .ge.-0.05 .and.
     *       Quality_flag_forJoint.eq.3) 
     *      CALL POPULATE_TAU_LAND_OCEAN(IDATA,
     *     SDS_Tau_Land_Ocean,SDSTAU_corrected,NUMCELLS,Land_Sol3) 
           
    
C Filled with Fill_Value for land
C
            Qcontrol=-7
            Call Fill_QAflag_ocean(QA_Flag_Ocean,SDS_QCONTROL_Ocean,
     *              Idata)
            SDS_CLDFRC_OCEAN(IDATA)=-99
            CALL FILLVALUE_Ocean(IDATA,SDS_ref,SDS_ref_STD,
     *        SDSTAU_best,SDSTAUS_best,SDSTAUB_best,SDSTAU_average,
     *        SDSTAUS_average,SDSTAUB_average,SDS_Least_error,
     *        SDS_small_weighting,SDS_sol_INDX_small,SDS_sol_INDX_large,
     *        SDSASSY_best,SDSASSY_average,SDS_ccn,sds_mass_conc,
     *        SDSBACK_best,SDSBACK_average,SDS_effrad,SDS_AOT_model,
     *        SDS_RefF_best,SDS_RefF_average,SDS_TranF_best,SDS_TranF_average,
     *        SDS_angs_coeff1,SDS_angs_coeff2,SDS_SCAT_ANGLE_OCEAN,
     *        SDS_QCONTROL_ocean,SDS_NUMPIXELS,SDS_CLDFRC_OCEAN,
     *        SDS_Tau_Land_Ocean_img,Qcontrol_special,
     *        SDS_correc_small_weighting)
C endif for Land pixels
          ENDIF
C endif for water at ocean
          ENDIF

C
c endif for clouds
          ELSE
           NO_Ret_Land=NO_Ret_Land+1
             SDS_CLDFRC_land(IDATA)=-99
            CALL FILLVALUE_LAND(IDATA,SDS_Tau_Land_Ocean_img,
     *        SDS_Aerosol_Type,SDSTAU_corrected_213,
     *        SDS_SCAT_ANGLE_land,SDS_mass_conc_land,
     *        SDS_angs_coeff_land,SDS_CLDFRC_land,
     *     SDS_dust_weighting,
     *     SDS_est_uncer,SDS_RefF_land,SDS_TranF_land,
     *     SDS_NUMPIXELS_land,SDSTAU_corrected,SDS_ref_land,
     *     SDS_ref_STD_land,SDS_QCONTROL_land,SDSTAU_small_land,
     *      SDS_Surface_Reflectance_Land,
     *     SDS_Fitting_Error_Land,Qcontrol_special_land)
              index_wave=3
              QCONTROL_land_wav1=0
              QCONTROL_land_wav2=0
    
  
C1/2002   Fill small weighting for land
           SDS_ratio_small_Land_Ocean(IDATA) =SDS_dust_weighting(idata)
           SDS_Reflected_flux_Land_Ocean(IDATA)= SDS_RefF_land(IDATA,2) 
             CALL POPULATE_TAU_LAND_OCEAN(IDATA,SDS_Tau_Land_Ocean_img,
     *     SDSTAU_corrected,NUMCELLS,Land_Sol3)
          if(Quality_flag_forJoint.eq.3)
     *    CALL POPULATE_TAU_LAND_OCEAN(IDATA,
     *    SDS_Tau_Land_Ocean,SDSTAU_corrected,NUMCELLS,Land_Sol3)
           
        
            QCONTROL=-2
            CALL Fill_QAflag_ocean(QA_Flag_Ocean,SDS_QCONTROL_Ocean,
     *              Idata)
            SDS_CLDFRC_OCEAN(IDATA)=-99
            CALL FILLVALUE_Ocean(IDATA,SDS_ref,SDS_ref_STD,
     *        SDSTAU_best,SDSTAUS_best,SDSTAUB_best,SDSTAU_average,
     *        SDSTAUS_average,SDSTAUB_average,SDS_Least_error,
     *        SDS_small_weighting,SDS_sol_INDX_small,SDS_sol_INDX_large,
     *        SDSASSY_best,SDSASSY_average,SDS_ccn,sds_mass_conc,
     *        SDSBACK_best,SDSBACK_average,SDS_effrad,SDS_AOT_model,
     *        SDS_RefF_best,SDS_RefF_average,SDS_TranF_best,
     *        SDS_TranF_average,
     *        SDS_angs_coeff1,SDS_angs_coeff2,SDS_SCAT_ANGLE_OCEAN,
     *        SDS_QCONTROL_ocean,SDS_NUMPIXELS,SDS_CLDFRC_OCEAN,
     *        SDS_Tau_Land_Ocean_img,Qcontrol_special,
     *        SDS_correc_small_weighting)
C1/2002   Fill small weighting for ocean from average value
                if(SDSTAU_average(IDATA,2).gt.0) then
                 aa= ((SDSTAUS_average(IDATA,2)/(SCALE3+OFFSET3))
     *                       /(SDSTAU_average(IDATA,2)/SCALE3+OFFSET3))
                SDS_ratio_small_Land_Ocean(IDATA) = (aa*SCALE3+OFFSET3)
                else
               SDS_ratio_small_Land_Ocean(IDATA)=SDSTAU_average(IDATA,2)
               ENDIF
            SDS_Reflected_flux_Land_Ocean(IDATA)= SDS_RefF_average(IDATA,2)
           CALL POPULATE_TAU_LAND_OCEAN(IDATA,SDS_Tau_Land_Ocean_img,
     *       SDSTAU_average,NUMCELLS,NWAV_S)
            if(Quality_to_pass .gt.0 .and. Quality_to_pass .le.3)
     *       CALL POPULATE_TAU_LAND_OCEAN(IDATA,
     *     SDS_Tau_Land_Ocean,SDSTAU_average,NUMCELLS,NWAV_S)  
         ENDIF
     
   
          
 9000  CONTINUE

            Else 
          NUMSQ=EV_Frames(iscan)/3 
        call FILLVALUE_missing(SDS_ref,SDS_ref_STD,SDSTAU_best,SDSTAUS_best,SDSTAUB_best,
     *SDSTAU_average,SDSTAUS_average,SDSTAUB_average,SDS_Least_error,SDS_small_weighting,
     *SDS_sol_INDX_small,SDS_sol_INDX_large,SDSASSY_best,SDSASSY_average,SDS_ccn,sds_mass_conc,
     *SDSBACK_best,SDSBACK_average,SDS_effrad,SDS_AOT_model,SDS_RefF_best,SDS_RefF_average,
     *SDS_TranF_best,SDS_TranF_average,SDS_angs_coeff1,SDS_angs_coeff2,SDS_SCAT_ANGLE_OCEAN,
     *SDS_QCONTROL_ocean,SDS_NUMPIXELS,SDS_CLDFRC_ocean,SDS_Tau_Land_Ocean_img,
     *Qcontrol_special,SDS_correc_small_weighting,SDS_Aerosol_Type,SDSTAU_corrected_213,
     *SDS_SCAT_ANGLE_land,SDS_mass_conc_land,SDS_angs_coeff_land,SDS_CLDFRC_land,
     *SDS_dust_weighting,SDS_est_uncer,SDS_RefF_land,SDS_TranF_land,SDS_NUMPIXELS_land,
     *SDSTAU_corrected,SDS_ref_land,SDS_ref_STD_land,SDS_QCONTROL_land,SDSTAU_small_land,
     *SDS_Surface_Reflectance_Land,SDS_Fitting_Error_Land,Qcontrol_special_land,
     *index_wave,SDS_Mean_Reflectance_Land_All,SDS_SDev_Reflectance_Land_All,     
     *SDS_Path_Radiance_Land,SDS_Critical_Reflectance_Land,SDS_Error_Crit_Reflectance_Land,
     *SDS_Error_Path_Radiance_Land,SDS_QFlag_Critical_Ref_Land,
     *SDS_QFlag_Path_Radiance_Land,SDSLAT,SDSLON,SDS_SST,SDS_MTHET0,SDS_MTHET,
     *SDS_MPHI0,SDS_MPHI,SDS_Scattering_Angle,SDS_CldMskQA,SDS_Tau_Land_Ocean,
     *SDS_ratio_small_Land_Ocean,NUMSQ,SDS_Sea_Land_Flag,iscan)
        Endif
C
C Call to MOD04_OUT writes the SDS HDF arrays and writes the output file
C
        
       CALL MOD04_OUT(MODFIL_MOD04,NUMSQ,NUMSCAN_new,SDSLAT,SDSLON,
     *SDS_Tau_Land_Ocean_img,SDS_Scattering_Angle,SDS_SST,SDS_MTHET0,SDS_MTHET,SDS_MPHI0,SDS_MPHI,
     *SDS_CldMskQA,SDS_QCONTROL_ocean,SDS_NUMPIXELS,SDS_SCAT_ANGLE_OCEAN,SDS_CLDFRC_ocean,
     *SDS_small_weighting,SDS_Least_error,SDS_effrad,SDS_sol_INDX_small,SDS_sol_INDX_large,SDS_ccn, 
     *SDS_mass_conc,SDS_angs_coeff1, SDS_angs_coeff2,SDS_AOT_model,SDS_ref,SDS_ref_STD, SDSTAU_best,
     *SDSTAU_average,SDSTAUS_best,SDSTAUS_average,SDSTAUB_best,SDSTAUB_average,SDSASSY_best, 
     *SDSASSY_average, SDSBACK_best,SDSBACK_average,SDS_RefF_best, SDS_RefF_average,SDS_TranF_best,
     *SDS_TranF_average,SDS_Aerosol_Type,SDS_SCAT_ANGLE_land,SDS_mass_conc_land,SDS_angs_coeff_land,
     *SDS_CLDFRC_land,SDS_dust_weighting,SDS_est_uncer,SDS_RefF_land,SDS_TranF_land,
     *SDS_NUMPIXELS_land,SDSTAU_corrected,SDS_ref_land,SDS_ref_STD_land,SDS_QCONTROL_land, 
     *SDS_Mean_Reflectance_Land_All, SDS_SDev_Reflectance_Land_All,SDS_Path_Radiance_Land,
     *SDS_Critical_Reflectance_Land,SDS_Error_Crit_Reflectance_Land,SDS_Error_Path_Radiance_Land,  
     *SDS_QFlag_Critical_Ref_Land,SDS_QFlag_Path_Radiance_Land,SDS_correc_small_weighting,
     *SDS_QCONTROL_CritRef_land,SDS_ratio_small_Land_Ocean,SDSTAU_corrected_213,SDSTAU_small_land,
     *SDS_Reflected_flux_Land_Ocean,SDS_Surface_Reflectance_Land,SDS_Fitting_Error_Land,SDS_Tau_Land_Ocean,
     *SDS_Quality_Land_Ocean)
C
C endif for scantype
c      ENDIF

9999  CONTINUE

       Xsize=NUMSQ
       Ysize=NUMSCAN_new
         
C Write Cloud mask and cloud distance to HDF file........      
        CALL MOD04_OUT_Extra(MODFIL_MOD04,Xsize,Ysize,SDS_Sea_Land_Flag,
     &   Fid_Land_sea_Flag)
       CALL MOD04_OUT_Extra(MODFIL_MOD04,Xsize,Ysize,SDS_Sea_Sunglint_Flag,
     &   Fid_Glint_Angle)
       CALL MOD04_OUT_Extra(MODFIL_MOD04,Xsize,Ysize,SDS_NCEP_Wspeed,
     &  Fid_Wind_Speed) 
       CALL MOD04_OUT_Extra(MODFIL_MOD04,Xsize,Ysize,SDS_topo_Altitude,
     &  Fid_Altitude_Land)
      DO I=1,No_PSA
        QA_Metadata_MOD04(I)=0.0
      ENDDO

      Total_Grids=NUMSCAN_new*NUMSQ 
      IF(NO_Ret_Land.LT.Total_Grids) THEN
        QA_Metadata_MOD04(1)=100.*float(Success_Ret_Land)/
     *                      fLoat(Total_Grids-NO_Ret_Land)
      ELSE
        QA_Metadata_MOD04(1)=0.0
      ENDIF 
        NO_Ret_Ocean=Total_Grids-(Success_Ret_Ocean+ Fail_Ret_Ocean)
C
      IF(NO_Ret_Ocean.LT.Total_Grids) THEN
        QA_Metadata_MOD04(2)=100.*float(Success_Ret_Ocean)/
     *                      fLoat(Total_Grids-NO_Ret_Ocean)
      ELSE
        QA_Metadata_MOD04(2)=0.0
      ENDIF
 
        CALL MOD04_Metadata(QA_Metadata_MOD04,Modfil_L1B_1km,
     *  RTN_NCEP,groups,HDFAttNms,NumHandles)
 
C Close all input/output files.
C
 
1111  CALL FILE_CLOSE(choice_flag,error_flag,
     &  MODFIL_L1B_1km,MODFIL_L1B_500,MODFIL_L1B_250,
     &  MODFIL_Geo,MODFIL_Cldmsk,MODFIL_MOD04,MODFIL_MOD05,
     &  MODFIL_MOD07,MODFIL_ANC,
     &  HANDLE_LUT466,HANDLE_LUT553,HANDLE_LUT644,HANDLE_LUT213,
     &  HANDLE_LUTMAP,HANDLE_S,HANDLE_L,
     &  handle_gascoeff,HANDLE_OUTSCI,HANDLE_QC,
     &  groups,HDFAttNms,NumHandles)

c     Get satellite instrument name.
         rtn = pgs_pc_getconfigdata(LUN_Sat_Instrument,pcf_satid)
         if (rtn .ne. 0) then
            call message( 'MOD_PR04_V2',
     &           'Error reading instrument name from pcf LUN 800510.' //
     &           ' [OPERATOR ACTION: Check PCF file, fix if corrupt.]',
     &           0, 1 )
         endif


      file_version = 1
      HDF_FILENAME = ' '
      rtn = -1
      rtn = pgs_pc_getreference(LRN_MOD04,FILE_VERSION,HDF_FILENAME)
      if (rtn.ne.PGS_S_SUCCESS) then
         error_flag = .true.
         msg=
     1        'pgs_pc_getreference for MOD_PR04 file'
     2        // char(10) // 'Operator Action:  Check system resources/environment, '
     3        // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4        // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msg,
     &        'MOD_PR04_V2.f')
      else
         
         if (pcf_satid .eq. 'AM1M') then
            doi_char = '10.5067/MODIS/MOD04_3K.006'
         else
            doi_char = '10.5067/MODIS/MYD04_3K.006'
         endif
         
         sd_id = sfstart(HDF_FILENAME, DFACC_WRITE)
         if(sd_id .eq. -1) then
            call message( 'MOD_PR04_V2',
     &           'Problem openning the file', 0, 2 )
         endif
         
         rtn = sfscatt(sd_id, 'identifier_product_doi', DFNT_CHAR8, 26,  
     +        doi_char) 
         if (rtn .eq. -1) then
            call message( 'MOD_PR04_V2',
     &           'Problem writting the global attribute identifier_product_doi', 0, 2 )
         endif
         
         rtn = sfscatt(sd_id, 'identifier_product_doi_authority', DFNT_CHAR8, 17,  
     +        'http://dx.doi.org') 
         if (rtn .eq. -1) then
            call message( 'MOD_PR04_V2',
     &           'Problem writting the global attribute identifier_product_doi_authority', 0, 2 )
         endif
         
         rtn = sfend(sd_id)
         if (rtn .eq. -1) then
            call message( 'MOD_PR04_V2',
     &           'Problem closing the file', 0, 2 )
         endif
      end if
 
      RETURN
      END

C*********************************************************************
      SUBROUTINE SET_INDEX(START_500,END_500,START_250,END_250,
     &                   START_1KM,END_1KM,IDATA)
C-----------------------------------------------------------------------
C!F77
C
C!DESCRIPTION:
C
C        This subroutine sets up the indexing for 10* 10 pixel boxes
C
C!INPUT PARAMETERS:
C
C        IDATA         10*10 array index on the swath
C
C!OUTPUT PARAMETERS:
C
C        START_500     Starting  Index for 500 meter resolution data
C        END_500       Ending    Index for 500 meter resolution data
C        START_250     Starting  Index for 250 meter resolution data
C        END_250       Ending    Index for 250 meter resolution data
C        START_1KM     Starting  Index for 1 km  resolution data
C        END_1KM       Ending    Index for 1 km  resolution data
C
C
C!REVISION HISTORY:
c 10/15/1999 fhliang
c fixed prolog.
C
C!TEAM-UNIQUE HEADER:
C
c!END
C----------------------------------------------------------------------

      IMPLICIT  NONE
      INCLUDE 'mod04.inc'
      INTEGER START_500,END_500,START_250,END_250,IDATA,START_1KM
      INTEGER END_1KM
      IF( IDATA .EQ.1) THEN
        START_1KM=1
        END_1KM=iline
        START_500=1
        END_500=iline*2
        START_250=1
        END_250=iline*4
      ELSE
        START_500=START_500+iline*2
        END_500=END_500+iline*2
        START_250=START_250+iline*4
        END_250=END_250+iline*4
        START_1KM=START_1KM+iline
        END_1KM=END_1KM+iline
      ENDIF

      RETURN
      END


C*********************************************************************
      SUBROUTINE LAND_WATER_CLOUD(DET_Flag,CldMsk_1km,CldMsk_500,
     *         CldMsk_250,High_Cloud_Flag,High_Cloud_Flag_500,
     *         LandSea_Flag,SunGlint_Flag,SnowIce_Flag,SnowMsk_Ratio,
     *         Shadow_Flag,cloud_num,cloud_num_land,Water,Land,Glint,
     *         Snowice,START_1KM,END_1KM,QA_Flag_Land,QA_Flag_Ocean,
     *         QA_Temp,snow_present_Flag,
     *  Quality_cirrus, Ret_Quality_cirrus,Land_CLDMSK_forfraction,
     *  coastal_Pix,Pure_Land)

C---------------------------------------------------------------------
C!F77
C
C!DESCRIPTION: This subroutine counts water, land, cloudy, glint and
C             snow/ice pixels in a 10x10 km area
C
C             (1) If 100% water pixels go to water process,
C                 otherwise goes to land.
C             (2) If 100% cloudy pixels, no aerosol retrieval is done.
C             (3) Glint and snow/ice pixels will be set to as cloudy pixel
C                 to prohibit any retrievals.
C
C!INPUT PARAMETERS:
C
C       CldMsk_1km     1km cloud mask
C       LandSea_Flag   Land_sea Flag from cloud mask
C       START_1KM      Starting Index for 1 km  resolution data
C       END_1KM        Ending Index for 1 km  resolution data
C
C!OUTPUT PARAMETERS:
C
C        Cloud_num     Number of cloudy pixels
C        Water         Number of pixels which identify water
C        Land          Number of pixels which identify land
C
c
C!REVISION HISTORY:
c 10/20/1999 fhliang
c fixed prolog.
C!TEAM-UNIQUE HEADER:
C
C This software was developed by the MODIS Atmosphere Science Team
C for the National Aeronautics and Space Administration at
C Goddard Space Flight Center.
C
C!END
C----------------------------------------------------------------------

      IMPLICIT NONE
      SAVE

      INCLUDE 'mod04.inc'

      INTEGER Water,Land,Desert,Glint,Snowice,cloud_num,cloud_num_land,
     &        IYY,IXX,I,J,START_1KM,END_1KM,Det_cldmsk,Y4_offset,X4_offset,
     &        Y2_offset,X2_offset,Y4,X4,Y2,X2,l,p,l11,p11,N,coast
      INTEGER LandSea_Flag(ISWATH,ILINE),CldMsk_1km(ISWATH,ILINE),
     &        CldMsk_250(4*ISWATH,4*ILINE),CldMsk_500(2*ISWATH,2*ILINE),
     &        SunGlint_Flag(ISWATH,ILINE),SnowIce_Flag(ISWATH,ILINE),
     &        SnowMsk_Ratio(ISWATH,ILINE),     
     &        DET_Flag(ISWATH,ILINE),Shadow_Flag(ISWATH,ILINE),
     &        High_Cloud_Flag(ISWATH,ILINE),
     &        High_Cloud_Flag_500(2*ISWATH,2*ILINE)
       Integer Land_CLDMSK_forfraction(ISWATH,ILINE)
      INTEGER QA_Flag_Land(19),QA_Flag_Ocean(12),for_percent
      integer   Ret_Quality_cirrus,Quality_cirrus(ISWATH,ILINE)
      integer snow_present_Flag,coastal_Pix,Pure_Land
      BYTE QA_Temp
      Ret_Quality_cirrus=1
      QA_Temp=0
      Det_cldmsk=0
      cloud_num=0
      cloud_num_land=0
      Water=0
      Land=0
      Desert=0
      Glint=0
      Snowice=0
      snow_present_Flag=0
      coastal_Pix =0
      Pure_Land =0
C cloud_num and cloud_num _land after cloud mask at 1km  before ice,snow flag 
         DO IYY = 1,IGRIDY
          DO IXX=START_1KM,END_1KM 
            IF(CldMsk_1km(IXX,IYY).EQ.0) THEN
              cloud_num=cloud_num+1
            ENDIF

            IF(Land_CLDMSK_forfraction(IXX,IYY).EQ.0)then
            cloud_num_land =cloud_num_land +1
            ENDIF
       ENDDO
          ENDDO
            

C Calculating number of land, water, cloudy, sun-glint and snow/ice

C pixels
C
          DO IYY = 1,IGRIDY
          DO IXX=START_1KM,END_1KM
c       IF(DET_Flag(IXX,IYY) .eq.1 .and.LandSea_Flag(IXX,IYY).EQ.0 )THEN
c       ELSE IF((DET_Flag(IXX,IYY) .eq.1.and. LandSea_Flag(IXX,IYY).EQ.3) 
c     *  .or.(DET_Flag(IXX,IYY) .eq.1 .And.LandSea_Flag(IXX,IYY).EQ. 2 ))THEN
C   Changed flag 11/26/2012
c		0:	Shallow Ocean (Ocean <5k from coast OR <50m deep).
c        	1:	Land (not anything else).
c        	2:	Ocean Coastlines and Lake Shorelines.
c		3:	Shallow Inland Water (Inland Water < 5km from shore
c				OR < 50m deep).
c		4:	Ephemeral (intermittent) Water.
c		5:	Deep Inland Water (Inland water > 5km from shoreline
c				AND > 50m deep).
c		6:	Moderate or Continental Ocean (Ocean > 5km from coast
c				AND > 50m deep AND < 500m deep).
c		7:	Deep Ocean (Ocean > 500m deep).

C Count Coastal pixels
        
 
           If(LandSea_Flag(IXX,IYY) . eq. 0 .or.
     &        LandSea_Flag(IXX,IYY) . eq. 3 .or.
     &        LandSea_Flag(IXX,IYY) . eq. 5 .or.
     &        LandSea_Flag(IXX,IYY) . eq. 6 .or.
     &        LandSea_Flag(IXX,IYY) . eq. 7) then
              water=water+1 
           Elseif(LandSea_Flag(IXX,IYY) . eq. 1 .or.
     &            LandSea_Flag(IXX,IYY) . eq. 2 .or.
     &            LandSea_Flag(IXX,IYY) . eq. 4)then
               Land=Land+1 
               if( LandSea_Flag(IXX,IYY) . eq. 2 )coastal_Pix=coastal_Pix+1
               if(LandSea_Flag(IXX,IYY) . eq. 1 .or.
     &          LandSea_Flag(IXX,IYY) . eq. 4)Pure_Land=Pure_Land+1 
             ENDIF
c  Enddo for 10 * 10   box.             
           ENDDO
         ENDDO 
C  Count snowy pixels based on SnowIce_Flag( MOD35) and  SnowMsk_Ratio (Rong_rong) when clear sky Shana Mattoo( 4/2013)         
         If(Land .gt.0) then         
          DO IYY = 1,IGRIDY
          DO IXX=START_1KM,END_1KM    
! save flag if cloud mask is clear and  snow is present           
       IF(CldMsk_1km(IXX,IYY) .gt.0 .and.
     *   (SnowIce_Flag(IXX,IYY).EQ.0 .OR.
     *    SnowMsk_Ratio(IXX,IYY).EQ.0))THEN
           snow_present_Flag=snow_present_Flag+1
         Endif
C Enddo for 10 * 10   box.
          ENDDO
          ENDDO 
C Exclude snow-ice and Shadow_Flag( read from MOD35) and 
C SnowMsk_Ratio computed in Trans_twoway
Change 2/21/2001 5x5 window neighboring pixel screening
        
          DO IYY = 1,IGRIDY
          DO IXX=START_1KM,END_1KM 
           IF(SnowIce_Flag(IXX,IYY).EQ.0.OR.
     &         SnowMsk_Ratio(IXX,IYY).EQ.0.OR.
     &         Shadow_Flag(IXX,IYY).EQ.0) THEN 
                CldMsk_1km(IXX,IYY)  =0 
             ENDIF 
c  Enddo for 10 * 10   box.
          ENDDO
          ENDDO  
         
C Setting 250 and 500 meter resolution  cloud masks created earlier with 0 after adding snow flags.
        DO IYY = 1,IGRIDY
          DO IXX=START_1KM,END_1KM

          IF(CldMsk_1km(IXX,IYY).EQ.0) THEN

             Y4_offset = (IYY-1)*4
             X4_offset = (IXX-1)*4

             DO l = 1,4
             DO p = 1,4
               Y4 = Y4_offset + l
               X4 = X4_offset + p
               CldMsk_250(X4,Y4) = 0
             ENDDO
             ENDDO

             Y2_offset = (IYY-1)*2
             X2_offset = (IXX-1)*2

             DO l = 1,2
             DO p = 1,2
               Y2 = Y2_offset + l
               X2 = X2_offset + p
               CldMsk_500(X2,Y2) = 0
             ENDDO
             ENDDO

          ENDIF

          ENDDO
        ENDDO
         
c   Setting quality Flag to retrive with lower quality for cirrus or smoke?
       DO IYY = 1,IGRIDY
          DO IXX=START_1KM,END_1KM
      
             Y4_offset = (IYY-1)*4
             X4_offset = (IXX-1)*4
             N=0
             DO l = 1,4
             DO p = 1,4
               Y4 = Y4_offset + l
               X4 = X4_offset + p
               IF(quality_cirrus(IXX,IYY).EQ.1 .and.
     *      CldMsk_250(X4,Y4) .eq.1) N=N+1
             ENDDO
             ENDDO
            IF(N .EQ. 16) Ret_Quality_cirrus=0
       ENDDO 
       ENDDO 
! Endif for Land only       
       Endif  
      RETURN
      END


C**********************************************************************
      SUBROUTINE Get_Metatable(Modfil_L1B_1km,NUMSCAN,Scan_Type,
     *EV_Frames,Scan_Start_Time,Imonth)

C----------------------------------------------------------------------
C!F77
C
C!DESCRIPTION:   This subroutine reads the Meta data table
C
C!INPUT PARAMETERS:
C
C      modfil_L1B_L1B
C
C!OUTPUT PARAMETERS:
C
C      NUMSCAN       Total number of scan ina granuale
C      Dead_detectors  Value to indicate bad detector
C      Scan_Type      day or night scan type
C      EV_Frames      Last number of pixel which is good.
C
c
C!REVISION HISTORY:
c 10/20/1999 fhliang
c fixed prolog.
C
C!Modifie by
C
C      Dr. Allen Chu            5/14/97
C      Code 913/SSAI
C      NASA Goddard Space Flight Center
C      Greenbelt, MD 20771
C
C!Team-Unique Header:
C
C!End
C-----------------------------------------------------------------------

       IMPLICIT NONE
       SAVE

       INCLUDE 'mod04.inc'
       INCLUDE 'mapi.inc'
       INCLUDE 'PGS_MODIS_39500.f'

       CHARACTER*20 att_n,dtype
       CHARACTER*27 asciiutc
       INTEGER pgs_td_taitoutc
       INTEGER nms,RTN,Modfil_L1B_1km(MODFILLEN),ik
       INTEGER NUMSCAN,IYear,IMonth(maxnum_scans),I,J,ITemp
       REAL*8 Day_Secs,Temp1,Days,Normal_Days(12),ANormal_Days(12)
       DATA Day_Secs/86400.0/
       DATA Normal_Days/31.,59.,90.,120.,151.,181.,212.,243,
     &                  273.,304.,334.,365./
       DATA ANormal_Days/31.,60.,91.,121.,152.,182.,213.,244,
     &                  274.,305.,335.,366./

C
C Call GMFIN to get 'Number of Scans' from L1B dataset
C

       att_N = 'Number of Scans'
       dtype = 'INTEGER*4'
       nms = 1
       RTN = GMFIN(Modfil_L1B_1km,att_N,dtype,nms,NUMSCAN)
       IF (rtn.ne.0) CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     &   'MAPI function GMFIN failed',
     &   'MOD_PR04_V2.f SUBROUTINE Get_Metatable for number of scans')

 
        bsize = maxnum_scans*4

        RTN = GMTBL(Modfil_L1B_1km,'Level 1B Swath Metadata',' ',
     &     'EV_Frames',0,NUMSCAN,bsize,EV_Frames)
C
C
      IF (RTN.ne.0) CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     &   'MAPI function GMTBL failed',
     &   'MOD_PR04_V2.f SUBROUTINE Get_Metatable for EV_Frames')
C
C Call GMTBL to get 'scan type' ( day, night or otherwise' from L1B cdataset
C
        bsize = maxnum_scans*4

        RTN = GMTBL(Modfil_L1B_1km,'Level 1B Swath Metadata',' ',
     &     'Scan Type',0,NUMSCAN,bsize,Scan_Type)
C
C
      IF (RTN.ne.0) CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     &   'MAPI function GMTBL failed',
     &   'MOD_PR04_V2.f SUBROUTINE Get_Metatable for Scan Type')

 
C
C Call GMTBL to get 'scan start time' (day, night, mixed or otherwise)
C
      bsize = maxnum_scans*8

      RTN = GMTBL(Modfil_L1B_1km,'Level 1B Swath Metadata',' ',
     &   'EV Sector Start Time',0,NUMSCAN,bsize,Scan_Start_Time)
C
C
      IF (RTN.ne.0) CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     & 'MAPI function GMTBL failed',
     & 'MOD_PR04_V2.f SUBROUTINE Get_Metatable for Scan_Start_Time')

 
C
C Derive MONTH from seconds since 1993-1-1
C
      DO J=1,NUMSCAN

      rtn=pgs_td_taitoutc(Scan_Start_Time(j),asciiutc)

      IF(rtn.eq.0) THEN
        Read(asciiutc(6:7),'(I2)') IMonth(J)
      ELSE
        call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,
     & 'Convertion of TAI to UTC time is incorrect',
     & 'MOD_PR04_V2.f')
      ENDIF

      ENDDO

C
C      DO J=1,NUMSCAN
C      Days=Scan_Start_Time(J)/Day_Secs
C      Temp1=Days/365.
C      IYear=Temp1
C
C      DO I=1,IYear
C        Days=Days-365
C      ENDDO
C
C Correction for leap-year
C
C      IF(IYEAR.GT.4)Days=Days-1
C      IF(IYEAR.GT.8)Days=Days-2
C      IF(IYEAR.GT.12)Days=Days-3
C
C Estimation of the month
C
C      ITemp=Days/30.
C      IF(IYEAR.EQ.3.OR.IYEAR.EQ.6.OR.IYEAR.EQ.9) THEN
C        Days=Days-ANormal_Days(ITemp)
C      ELSE
C        Days=Days-Normal_Days(ITemp)
C      ENDIF
C      IMonth(J)=ITEMP
C      IF(Days.GT.0.0) IMonth(J)=ITemp+1
C
C      ENDDO
C

      RETURN
      END


C*********************************************************************
      SUBROUTINE POPULATE_COMMON_ARRAYS(IDATA,Lat_center,Lon_center,
     &         SST,MTHET0,MTHET,MPHI0,MPHI,MSCATT,QA_Temp,SDSLAT,
     &         SDSLON,SDS_SST,SDS_MTHET0,SDS_MTHET,SDS_MPHI0,SDS_MPHI,
     &         SDS_Scattering_Angle,SDS_CldMskQA)
C
C-----------------------------------------------------------------------
c!F77
C
c!Description:
C This subroutine populates common arrays such as latitude,
C longitude, and the angles of sun, satellite and the relative
C azimuth
c!INPUT PARAMETERS:
C
c        Lat_center    Latitude  at center of 10*10 box
C        Lon_center    Londitude at center of 10*10 box
C              SST     Time
c           MTHET0     Solar Zenith Angle at center of 10*10 box
C            MTHET     View  Angle at center of 10*10 box
C            MPHI0     Satellite Azimuth  at center of 10*10 box
C             MPHI     Solar Azimuth  at center of 10*10 box
C            MDPHI     Diff Satellite Azimuth  and Solar Azimuth
C
c!OUTPUT PARAMETERS:
C        all SDS* for all above variables for HDF write.
C
c!Revision History:
C
C
c!Team-unique Header:
C
C
c!End
C----------------------------------------------------------------------

      IMPLICIT  NONE
      SAVE

      INCLUDE 'mod04.inc'

      INTEGER IDATA
      REAL Lat_center,Lon_center,MTHET0,MTHET,MPHI0,MPHI,MSCATT
      REAL SDSLAT(NUMCELLS),SDSLON(NUMCELLS)
      INTEGER*2 SDS_MTHET0(NUMCELLS),SDS_MTHET(NUMCELLS),
     &          SDS_MPHI0(NUMCELLS),SDS_MPHI(NUMCELLS),
     &          SDS_Scattering_Angle(NUMCELLS)
      REAL*8    SDS_SST(NUMCELLS),SST
      BYTE      SDS_CldMskQA(NUMCELLS),QA_Temp



      SDSLAT(IDATA) = LAT_center*SCALE1+OFFSET1
      SDSLON(IDATA) = Lon_center*SCALE1+OFFSET1
      SDS_SST(IDATA)= SST*SCALE1+OFFSET1
      SDS_MTHET0(IDATA) = MTHET0*SCALE2+OFFSET1
      SDS_MTHET(IDATA) =  MTHET*SCALE2+OFFSET1
      SDS_MPHI0(IDATA) = MPHI0*SCALE2+OFFSET1
      SDS_MPHI(IDATA) = MPHI*SCALE2+OFFSET1
      SDS_Scattering_Angle(IDATA)=MSCATT*SCALE2+OFFSET1
      SDS_CldMskQA(IDATA) = QA_Temp*SCALE1+OFFSET1

      RETURN
      END




C************************************************************************
      SUBROUTINE TRANS_2WAY(QA_Flag_Land,QA_Flag_Ocean,START_250,END_250,
     *   START_500,END_500,START_1KM,END_1KM,RTN_MOD05,
     *   RTN_MOD07,RTN_NCEP,RTN_DAO,RTN_TOMS,RTN_TOVS,
     *   ISWATH,ILINE,SatZen,SolZen,Refl_1,Refl_2,Refl_3,Refl_4,Refl_5,
     *   Refl_6,Refl_7,Refl_9,Refl_12,Refl_13,Refl_16,Refl_26,Rad_20,
     *   Rad_31,Rad_32,LandSea_Flag,CldMsk_1km,SnowMsk_500m,
     *   SnowMsk_Ratio,Total_H2O,Total_O3,
     *   SLOPE_MEAN_LAND,SLOPE_MEAN_OCEAN,W659_SYN,W865_SYN,W470_SYN,
     *   W550_SYN,W124_SYN,W164_SYN,W213_SYN,
     *   W395_SYN,W138_SYN,W1100_SYN,W1100_SYN_500m,
     *   W1200_SYN,G_factor,Multi_factor)

C-----------------------------------------------------------------------
C!F77
C
C!Description:
C
C         To correct reflectance by gases (H2O,O3,CO2) absorption and
C         cirrus clouds.
C
C!Input Parameters:
C
C         QA_Flag_Land     : QA flags over land
C         QA_Flag_Ocean    : QA flags over ocean
C         RTN_MOD05         : Return number of MOD05 product
C         RTN_MOD07         : Return number of MOD07 product
C         RTN_NCEP          : Return number of NCEP data
C         RTN_DAO           : Return number of DAO data
C         RTN_TOMS          : Return number of TOMS data
C         RTN_TOVS          : Return number of TOVS data
C         ISWATH            : Number of pixels across scan
C         ILINE             : NUmber of pixels along scan
C         SatZen            : Satellite zenith angle
C         SolZen            : Solar zenith angle
C         Refl_1            : Reflectance from MODIS channel 1
C         Refl_2            : Reflectance from MODIS channel 2
C         Refl_3            : Reflectance from MODIS channel 3
C         Refl_4            : Reflectance from MODIS channel 4
C         Refl_5            : Reflectance from MODIS channel 5
C         Refl_6            : Reflectance from MODIS channel 6
C         Refl_7            : Reflectance from MODIS channel 7
C         Refl_9            : Reflectance from MODIS channel 9
C         Refl_12           : Reflectance from MODIS channel 12
C         Refl_13           : Reflectance from MODIS channel 13
C         Refl_16           : Reflectance from MODIS channel 16
C         Refl_26           : Reflectance from MODIS channel 26
C         Rad_20            : Radiance from MODIS channel 20
C         Rad_31            : Radiance from MODIS channel 31
C         Refl_32           : Radiance from MODIS channel 32
C         Total_H2O         : Total precipitable water
C         Total_O3          : Total ozone
C
C!Output Parameters:
C
C         W659_SYN          : Corrected reflectance for MODIS channel 1
C         W865_SYN          : Corrected reflectance for MODIS channel 2
C         W470_SYN          : Corrected reflectance for MODIS channel 3
C         W550_SYN          : Corrected reflectance for MODIS channel 4
C         W124_SYN          : Corrected reflectance for MODIS channel 5
C         W164_SYN          : Corrected reflectance for MODIS channel 6
C         W213_SYN          : Corrected reflectance for MODIS channel 7
C         W443o_SYN         : Corrected reflectance for MODIS channel 9
C         W551o_SYN         : Corrected reflectance for MODIS channel 12
C
C         W869o_SYN         : Corrected reflectance for MODIS channel 16
C         W138_SYN          : Corrected reflectance for MODIS channel 26
C         W395_SYN          : Corrected radiance for MODIS channel 20
C         W1100_SYN         : Corrected radiance for MODIS channel 31
C         W1200_SYN         : Corrected radiance for MODIS channel 32
C
C!Revision History:
C
C         WRITTEN BY:
C         Dr. Allen Chu
C         Code 913/SSAI
C         NASA Goddard Space Flight Center
C         Greenbelt, MD 20771
C
C!Team-unique Header:
C
C
C!End
C-----------------------------------------------------------------------

      IMPLICIT NONE

C *\ Global variables

      INTEGER START_500,END_500,START_250,END_250,START_1KM,END_1KM
      INTEGER ISWATH,ILINE,RTN_MOD05,RTN_MOD07,RTN_NCEP,RTN_DAO,
     *        RTN_TOMS,RTN_TOVS,
     *        QA_Flag_Land(19),QA_Flag_Ocean(12),ISNOW
      INTEGER LandSea_Flag(ISWATH,ILINE),CldMsk_1km(ISWATH,ILINE),
     *        SnowMsk_500m(2*ISWATH,2*ILINE),SnowMsk_Ratio(ISWATH,ILINE)

       REAL SatZen,SolZen,Total_H2O,Total_O3,SLOPE_MEAN_LAND(3),
     *      SLOPE_MEAN_OCEAN(3)

      REAL Refl_1 (4*ISWATH,4*ILINE),Refl_2 (4*ISWATH,4*ILINE),
     *     Refl_3 (2*ISWATH,2*ILINE),Refl_4 (2*ISWATH,2*ILINE),
     *     Refl_5 (2*ISWATH,2*ILINE),Refl_6 (2*ISWATH,2*ILINE),
     *     Refl_7 (2*ISWATH,2*ILINE), 
     *     Refl_9 (ISWATH,ILINE),Refl_12(ISWATH,ILINE),
     *     Refl_13(ISWATH,ILINE),Refl_16(ISWATH,ILINE),
     *     Refl_26(ISWATH,ILINE),
     *     Rad_20(ISWATH,ILINE),Rad_31(ISWATH,ILINE),
     *     Rad_32(ISWATH,ILINE)

      REAL W659_SYN(4*ISWATH,4*ILINE),W865_SYN(4*ISWATH,4*ILINE),
     *     W470_SYN(2*ISWATH,2*ILINE),W550_SYN(2*ISWATH,2*ILINE),
     *     W124_SYN(2*ISWATH,2*ILINE),W164_SYN(2*ISWATH,2*ILINE),
     *     W213_SYN(2*ISWATH,2*ILINE), 
     *     W395_SYN(ISWATH,ILINE),W138_SYN(ISWATH,ILINE),
     *     W1100_SYN(ISWATH,ILINE),W1100_SYN_500m(2*ISWATH,2*ILINE),
     *     W1200_SYN(ISWATH,ILINE),Multi_factor(10)


C *\ Local variables

      INTEGER IL,IS,IL2,IS2,IL4,IS4,IMASK2,JMASK2,IMASK4,JMASK4,
     *        N_Red,N_0p86,I
      REAL G_factor,LOGCON,LOGCON2,EXPONENT
      REAL Total_Red, Corr_Cirrus(3),Corr_Cirrus_S,Total_0p86

      REAL Ratio
      INTEGER ij, Y2_offset, X2_offset,l,p,y2,x2,iscan,idata

c temperature conversion variables
      REAL  Planck_constant,Speed_light,Boltz_cons
      REAL   w_meter,wave
      REAL  wav1,wav2,wav3,c1,c2
      PARAMETER(Planck_constant=6.6260755e-34,
     &   Speed_light=2.9979246e+8, Boltz_cons=1.380658e-23,
     &   wav1=3.75,wav2=11.0,wav3=12.0)

      REAL DEGRAD
C Intialize the data    
 
          
              DO    IL=1,ILINE
              DO    IS=START_1KM,END_1KM 
           W395_SYN(IS,IL)=-999.0
           W138_SYN(IS,IL)=-999.0 
           W1100_SYN(IS,IL)=-999.0 
           W1200_SYN(IS,IL)=-999.0
           Enddo
           Enddo 
              DO    IL=1,ILINE*2
              DO    IS=START_1KM*2,END_1KM*2
              W470_SYN(IS,IL)=-999.0
              W550_SYN(IS,IL)=-999.0
              W124_SYN(IS,IL)=-999.0
              W164_SYN(IS,IL)=-999.0
              W213_SYN(IS,IL)=-999.0
               W1100_SYN_500m(IS,IL)=-999.0
              Enddo
              Enddo 
              
              DO    IL=1,ILINE*4
              DO    IS=START_1KM*4,END_1KM*4 
              W659_SYN(IS,IL)=-999.0
              W865_SYN(IS,IL)=-999.0
              Enddo
              Enddo
              
 
         DO 100 IL=1,ILINE
            DO 110 IS=START_1KM,END_1KM

        IMASK2=2*IS -1
        JMASK2=2*IL -1
        IMASK4=4*IS -3
        JMASK4=4*IL -3

C
C Correct radiance due to gaseous absorption and cirrus clouds contamination
C

          DO 200 IL4=JMASK4,4*IL
            DO 210 IS4=IMASK4,4*IS
              IF(Refl_1(IS4,IL4).GT.0.0) THEN
                W659_SYN(IS4,IL4)=Refl_1(IS4,IL4)*Multi_factor(3)
              ELSE
                W659_SYN(IS4,IL4)=Refl_1(IS4,IL4) 
              ENDIF
              IF(Refl_2(IS4,IL4).GT.0.0) THEN
                W865_SYN(IS4,IL4)=Refl_2(IS4,IL4)*Multi_factor(4) 
              ELSE
                W865_SYN(IS4,IL4)=Refl_2(IS4,IL4)
              ENDIF
210         CONTINUE
200       CONTINUE

          DO 300 IL2=JMASK2,2*IL
            DO 310 IS2=IMASK2,2*IS
              IF(Refl_3(IS2,IL2).GT.0.0) THEN
                W470_SYN(IS2,IL2)=Refl_3(IS2,IL2)*Multi_factor(1) 
              ELSE
                W470_SYN(IS2,IL2)=Refl_3(IS2,IL2)
              ENDIF
              IF(Refl_4(IS2,IL2).GT.0.0) THEN
                W550_SYN(IS2,IL2)=Refl_4(IS2,IL2)*Multi_factor(2) 
              ELSE
                W550_SYN(IS2,IL2)=Refl_4(IS2,IL2)
              ENDIF
              IF(Refl_5(IS2,IL2).GT.0.0) THEN
                W124_SYN(IS2,IL2)=Refl_5(IS2,IL2)*Multi_factor(5) 
              ELSE
                W124_SYN(IS2,IL2)=Refl_5(IS2,IL2)
              ENDIF
              IF(Refl_6(IS2,IL2).GT.0.0) THEN
                W164_SYN(IS2,IL2)=Refl_6(IS2,IL2)*Multi_factor(6) 
              ELSE
                W164_SYN(IS2,IL2)=Refl_6(IS2,IL2)
              ENDIF
              IF(Refl_7(IS2,IL2).GT.0.0) THEN
                W213_SYN(IS2,IL2)=Refl_7(IS2,IL2)*Multi_factor(7) 
              ELSE
                W213_SYN(IS2,IL2)=Refl_7(IS2,IL2)
              ENDIF
310         CONTINUE
300       CONTINUE 

C
C No change for the rest of 4 channels
C
           W138_SYN(IS,IL)=Refl_26(IS,IL)

 


C compute Temprature channels ( convert from radiance to temperature)
c Derive constants 
         c1=2.0*Planck_constant*(Speed_light*Speed_light)
         c2=(Planck_constant*Speed_light)/Boltz_cons
c convert wavelength to meters
           do ij=1,3
           if( ij .eq.1)wave=wav1
           if( ij .eq.2)wave=wav2
           if( ij .eq.3)wave=wav3
           w_meter=(1.0e-6*wave) 
          if( ij .eq.1)W395_SYN(IS,IL)=
     &   c2/(w_meter*alog(c1/(1.0e+6*Rad_20(IS,IL)*w_meter**5)+1.0)) 
          if( ij .eq.2)W1100_SYN(IS,IL)=
     &  c2/(w_meter*alog(c1/(1.0e+6*Rad_31(IS,IL)*w_meter**5)+1.0))
          if( ij .eq.3)W1200_SYN(IS,IL)=
     &  c2/(w_meter*alog(c1/(1.0e+6*Rad_32(IS,IL)*w_meter**5)+1.0))  
          enddo

110     CONTINUE
100   CONTINUE
c convert to 500 meter resolution
          DO IL = 1,ILINE
          DO IS=START_1KM,END_1KM
             Y2_offset =( IL-1)*2
             X2_offset = (IS-1)*2
        IF(W1100_SYN(IS,IL) .gt.0) then
                  DO l = 1,2
                  DO p = 1,2
                    Y2 = Y2_offset + l
                   X2 = X2_offset + p
                   W1100_SYN_500m(X2,Y2) =W1100_SYN(IS,IL)
               ENDDO
               ENDDO
          ENDIF
       ENDDO
       ENDDO
C Snow-masking scheme based upon Ref(1.24 micron)/Ref(0.86 micron) < 1 at 500 m resolution
C 2/26/2000

      DO 120 IL=1,ILINE*2
        DO 130 IS=START_500,END_500

          IMASK2=2*IS -1
          JMASK2=2*IL -1

          N_0p86=0
          Total_0p86=0.0
          DO 420 IL2=JMASK2,2*IL
            DO 430 IS2=IMASK2,2*IS
              IF(Refl_2(IS2,IL2).GT.0.0) THEN
                Total_0p86=Total_0p86+Refl_2(IS2,IL2)
                N_0p86=N_0p86+1
              ENDIF
430         CONTINUE
420       CONTINUE 
             Ratio=0
          IF(Total_0p86.GT.0.0.AND.Refl_5(IS,IL).GT.0.0) THEN
            IF(N_0p86.GT.0) THEN
c rong -rong def (ref86-ref1.24)/(ref86+ref1.24)

         Ratio=((Total_0p86/N_0p86)-Refl_5(IS,IL))/
     &         ((Total_0p86/N_0p86)+Refl_5(IS,IL))
            ELSE
              Total_0p86=0.0
                Ratio=0.0
            ENDIF 
          ELSE
            Total_0p86=0.0
             Ratio=0.0
          ENDIF

           SnowMsk_500m(IS,IL)=1
Cccc    IF(ratio.gt.0.20 .and.(Total_0p86/N_0p86).gt.0.08) THEN

        IF(ratio.gt.0.01 .and. W1100_SYN_500m(IS,IL) .lt. 285) THEN
            SnowMsk_500m(IS,IL)=0
          ENDIF 
         
130     CONTINUE
120   CONTINUE



   

      DO 140 IL=1,ILINE
        DO 150 IS=START_1KM,END_1KM

          IMASK2=2*IS -1
          JMASK2=2*IL -1
          ISNOW=0

          DO 440 IL2=JMASK2,2*IL
            DO 450 IS2=IMASK2,2*IS

              IF(SnowMsk_500m(IS2,IL2).EQ.0) THEN
                ISNOW=ISNOW+1
              ENDIF
      
450         CONTINUE
440       CONTINUE

          SnowMsk_Ratio(IS,IL)=1
          IF(ISNOW.GE.1) THEN
            SnowMsk_Ratio(IS,IL)=0
          ENDIF 
c           write(38,*)'sc',IL,IS,SnowMsk_Ratio(IS,IL)
150     CONTINUE
140   CONTINUE
C 2/26/2001
      RETURN
      END
C****************************************************************
      SUBROUTINE BYTE_SET(QA_V,Bit_SP,QA_Byte)
      IMPLICIT NONE
      SAVE
C-------------------------------------------------------------------
C
C !F77
C
C !DESCRIPTION:
C                 This subroutine is to set QA bit given nth byte
C                 Note: the QA byte setting starts from the rightmost
C                       (not the leftmost) bit
C
C        1st word       7th  6th  5th  4th  3rd  2nd  1st  0th  (bit)
C        2nd word       7th  6th  5th  4th  3rd  2nd  1st  0th  (bit)
C        3rd word       7th  6th  5th  4th  3rd  2nd  1st  0th  (bit)
C        4th word       7th  6th  5th  4th  3rd  2nd  1st  0th  (bit)
C        5th word       7th  6th  5th  4th  3rd  2nd  1st  0th  (bit)
C
C !INPUT PARAMETERS:
C
C        QA_V       QA parameter value
C        Bit_SP     Bit starting position of Ith QA parameter
C                    (see MODIS atmosphere QA plan)
C
C !OUTPUT PARAMETERS:
C
C        QA_Byte    Byte set for quality control
C
C !REVISION HISTORY:
C
C        WRITTEN BY
C        Dr. Allen Chu          11/25/97
C        Code 913
C        NASA Goddard Space Flight Center
C        Greenbelt, MD 20771
C
C !TEAM-UNIQUE HEADER:
C
C   This software was developed by the MODIS Atmosphere Science Team
C   for the National Aeronautics and Space Administration at
C   Goddard Space Flight Center.
C
C !END
C--------------------------------------------------------------------
C
      Intrinsic ibset,ibclr
      INTEGER QA_V,Bit_SP,Byte_Temp
      BYTE QA_Byte
      Byte_Temp=QA_Byte

      IF(QA_V.EQ.0) THEN

        Byte_Temp = ibclr(Byte_Temp,Bit_SP)

      ELSE IF(QA_V.EQ.1) THEN

        Byte_Temp = ibset(Byte_Temp,Bit_SP)

      ELSE IF(QA_V.EQ.2) THEN

        Byte_Temp = ibclr(Byte_Temp,Bit_SP)
        Byte_Temp = ibset(Byte_Temp,Bit_SP+1)

      ELSE IF(QA_V.EQ.3) THEN

        Byte_temp = ibset(Byte_Temp,Bit_SP)
        Byte_temp = ibset(Byte_Temp,Bit_SP+1)

      ELSE IF(QA_V.EQ.4) THEN

        Byte_Temp = ibclr(Byte_Temp,Bit_SP)
        Byte_Temp = ibclr(Byte_Temp,Bit_SP+1)
        Byte_temp = ibset(Byte_Temp,Bit_SP+2)

      ELSE IF(QA_V.EQ.5) THEN

        Byte_Temp = ibset(Byte_Temp,Bit_SP)
        Byte_Temp = ibclr(Byte_Temp,Bit_SP+1)
        Byte_temp = ibset(Byte_Temp,Bit_SP+2)

      ELSE IF(QA_V.EQ.6) THEN

        Byte_Temp = ibclr(Byte_Temp,Bit_SP)
        Byte_Temp = ibset(Byte_Temp,Bit_SP+1)
        Byte_temp = ibset(Byte_Temp,Bit_SP+2)

      ELSE IF(QA_V.EQ.7) THEN

        Byte_Temp = ibset(Byte_Temp,Bit_SP)
        Byte_Temp = ibset(Byte_Temp,Bit_SP+1)
        Byte_temp = ibset(Byte_Temp,Bit_SP+2)

      ELSE IF(QA_V.EQ.8) THEN

        Byte_Temp = ibclr(Byte_Temp,Bit_SP)
        Byte_Temp = ibclr(Byte_Temp,Bit_SP+1)
        Byte_temp = ibclr(Byte_Temp,Bit_SP+2)
        Byte_temp = ibset(Byte_Temp,Bit_SP+3)

      ELSE IF(QA_V.EQ.9) THEN

        Byte_Temp = ibset(Byte_Temp,Bit_SP)
        Byte_Temp = ibclr(Byte_Temp,Bit_SP+1)
        Byte_temp = ibclr(Byte_Temp,Bit_SP+2)
        Byte_temp = ibset(Byte_Temp,Bit_SP+3)

      ELSE IF(QA_V.EQ.10) THEN

        Byte_Temp = ibclr(Byte_Temp,Bit_SP)
        Byte_Temp = ibset(Byte_Temp,Bit_SP+1)
        Byte_temp = ibclr(Byte_Temp,Bit_SP+2)
        Byte_temp = ibset(Byte_Temp,Bit_SP+3)

      ELSE IF(QA_V.EQ.11) THEN

        Byte_Temp = ibset(Byte_Temp,Bit_SP)
        Byte_Temp = ibset(Byte_Temp,Bit_SP+1)
        Byte_temp = ibclr(Byte_Temp,Bit_SP+2)
        Byte_temp = ibset(Byte_Temp,Bit_SP+3)

      ELSE IF(QA_V.EQ.12) THEN

        Byte_Temp = ibclr(Byte_Temp,Bit_SP)
        Byte_Temp = ibclr(Byte_Temp,Bit_SP+1)
        Byte_temp = ibset(Byte_Temp,Bit_SP+2)
        Byte_temp = ibset(Byte_Temp,Bit_SP+3)

      ELSE IF(QA_V.EQ.13) THEN

        Byte_Temp = ibset(Byte_Temp,Bit_SP)
        Byte_Temp = ibclr(Byte_Temp,Bit_SP+1)
        Byte_temp = ibset(Byte_Temp,Bit_SP+2)
        Byte_temp = ibset(Byte_Temp,Bit_SP+3)

      ELSE IF(QA_V.EQ.14) THEN

        Byte_Temp = ibclr(Byte_Temp,Bit_SP)
        Byte_Temp = ibset(Byte_Temp,Bit_SP+1)
        Byte_temp = ibset(Byte_Temp,Bit_SP+2)
        Byte_temp = ibset(Byte_Temp,Bit_SP+3)

      ELSE IF(QA_V.EQ.15) THEN

        Byte_Temp = ibset(Byte_Temp,Bit_SP)
        Byte_Temp = ibset(Byte_Temp,Bit_SP+1)
        Byte_temp = ibset(Byte_Temp,Bit_SP+2)
        Byte_temp = ibset(Byte_Temp,Bit_SP+3)

      ENDIF

      QA_Byte=Byte_Temp

      RETURN
      END



C***********************************************************************
       SUBROUTINE  COMPUTE_GLINTANGLE(MTHET0,MTHET,MDPHI,
     *             GLINT_ANGLE,QA_Flag_Ocean)
C
C-----------------------------------------------------------------------
C!F77
C
C!DESCRIPTION: This subroutine computes Glint_angle
C
C
C!INPUT PARAMETERS:
C
C           MTHET0     Solar Zenith Angle at center of a 10*10 km^2 box
C            MTHET     View  Angle at center of a 10*10 km^2 box
C            MDPHI     Difference of Satellite Azimuth and Solar Azimuth
C                      Angle of a 10*10 km^2 box
C
C!OUTPUT PARAMETERS:
C
C      GLINT_ANGLE     Glint_angle
C    QA_Flag_Ocean     Runtime QA Flags for Ocean
C
C !REVISION HISTORY:
C
C        WRITTEN BY
C        Shana Mattoo
C        Code 913
C        NASA Goddard Space Flight Center
C        Greenbelt, MD 20771
C
C !TEAM-UNIQUE HEADER:
C
C   This software was developed by the MODIS Atmosphere Science Team
C   for the National Aeronautics and Space Administration at
C   Goddard Space Flight Center.
C
C !END
C
C----------------------------------------------------------------------

       IMPLICIT NONE
       SAVE

       INCLUDE 'mod04.inc'

       REAL GLINT_ANGLE, MTHET0,MTHET,MDPHI
       INTEGER QA_Flag_Ocean(12)

      GLINT_ANGLE = 0.0
        IF(MTHET0.GT.0.0.AND.MTHET.GT.0.0.AND.MDPHI.GT.0.0) THEN
         GLINT_ANGLE=(COS(MTHET0*DTR))*(COS(MTHET*DTR))
     *            +((SIN(MTHET0*DTR))*(SIN(MTHET*DTR))
     *            *( COS(MDPHI*DTR)))
         GLINT_ANGLE = (ACOS(GLINT_ANGLE))*RTD
        ENDIF

       IF(GLINT_ANGLE .GE.GLINT_THRESHOLD) then
         QA_Flag_Ocean(9)=1
       ELSE
         QA_Flag_Ocean(9)=0
       ENDIF

       RETURN
       END



C**************************************************************************
       SUBROUTINE POPULATE_TAU_LAND_OCEAN(IDATA,SDS_Tau_Land_Ocean_img,
     *     SDSTAU_array,NX,NY)
C-----------------------------------------------------------------------
C!F77
C
C!DESCRIPTION: This subroutine populates optical depth at 0.5 micron
C              for both land and ocean
C
C
C!INPUT PARAMETERS:
C
C              IDATA     index of grid box
C       SDSTAU_ARRAY     Array of optical depth at 0.5 micron for input
C
C!OUTPUT PARAMETERS:
C
C SDS_Tau_Land_Ocean    Common Array of optical depth at 0.5 micron for
C                       both land and ocean
C
C !REVISION HISTORY:
C
C        WRITTEN BY
C        Shana Mattoo
C        Code 913
C        NASA Goddard Space Flight Center
C        Greenbelt, MD 20771
C
C !TEAM-UNIQUE HEADER:
C
C   This software was developed by the MODIS Atmosphere Science Team
C   for the National Aeronautics and Space Administration at
C   Goddard Space Flight Center.
C
C !END
C
C----------------------------------------------------------------------

       IMPLICIT NONE
       SAVE

C       INCLUDE 'mod04.inc'

       INTEGER IDATA,NX,NY
       INTEGER*2 SDSTAU_array(NX,NY),SDS_Tau_Land_Ocean_img(NX)

       SDS_Tau_Land_Ocean_img(IDATA)= SDSTAU_array(IDATA,2)

       RETURN
       END
C**************************************************************************
       SUBROUTINE CldMsk_Land(Data_Size,ISWATH,ILINE,Refl_3,refl_5
     * ,Refl_26,CldMsk_250,CldMsk_500,CldMsk_1km, RMED,RMEDSQ, 
     &RMED_1km,RMEDSQ_1km,iscan,nrows,quality_cirrus,
     & Land_CLDMSK_forfraction)
C----------------------------------------------------------------------
C!F77
C
C!DESCRIPTION:   
C               Generate cloud mask over land using spatial variability
C               of 0.47 (>0.01) and 1.38 um (> 0.007) reflectance as well 
C               as absolute  value of 1.38 um > 0.1
C
C!INPUT PARAMETERS:
C
C             ISWATH      Number of pixels at 1 km resolution along scan
C              ILINE      Number of pixels at 1 km resolution against scan
C             Refl_3      Reflectance at 0.47 um
C            Refl_26      Reflectance at 1.38 um
C
C!OUTPUT PARAMETERS:
C
C         CLDMSK_1KM      Cloud mask at 1 km resolution
C         CldMsk_500      Cloud mask at 500m resolution
C         CldMsk_250      Cloud mask at 250m resolution
C
C !REVISION HISTORY:
C $Log: Process_ocean_V2.f,v $
c 10/18/1999 fhliang
c fixed prolog;
C
C !TEAM-UNIQUE HEADER:
C
C !REFERENCES AND CREDITS:
C
C !DESIGN NOTES:
C
C !END
C-----------------------------------------------------------------------
C
      IMPLICIT  NONE
      SAVE
      INTEGER ISTART,IEND,ISWATH,ILINE,Data_Size(2)
      INTEGER IX,IY,l,p,l2,p2,ll,pp,ll2,pp2,N,iscan,NROWS(ISWATH*2)
      INTEGER X2_offset,Y2_offset,X4_offset,Y4_offset
      Integer Land_CLDMSK_forfraction(ISWATH,ILINE)
      REAL Refl_3(ISWATH*2,ILINE*2),Refl_26(ISWATH,ILINE)
      real Refl_5(ISWATH*2,ILINE*2) 
      REAL THRHLD1380_1,THRHLD1380_2
      REAL cloud_threhold_land47_1,cloud_threhold_land47_2
      REAL RMED(ISWATH*2),RMEDSQ(ISWATH*2)
      REAL RMED_1km(ISWATH),RMEDSQ_1km(ISWATH)
      INTEGER CldMsk_250(ISWATH*4,ILINE*4),CldMsk_500(ISWATH*2,ILINE*2),
     *        CldMsk_1km(ISWATH,ILINE),quality_cirrus(ISWATH,ILINE)
c         real  Real_CldMsk_500(ISWATH*2,ILINE*2)
       DATA cloud_threhold_land47_1/0.0025/,cloud_threhold_land47_2/0.4/
c        DATA THRHLD1380_1/0.007/,THRHLD1380_2/0.01/
c       DATA cloud_threhold_land47_1/0.01/,cloud_threhold_land47_2/0.4/
        DATA THRHLD1380_1/0.003/,THRHLD1380_2/0.01/


C
C Initialize cloud mask arrays (cloudy)
C
      
      DO IY=1,Data_Size(2)*4
        DO  IX=1,Data_Size(1)*4
          CldMsk_250(IX,IY)=0
        ENDDO
      ENDDO

      DO IY=1,Data_Size(2)*2
        DO  IX=1,Data_Size(1)*2
          CldMsk_500(IX,IY)=0
        ENDDO
      ENDDO

      DO IY=1,Data_Size(2)
        DO  IX=1,Data_Size(1)
          CldMsk_1km(IX,IY)=0
           Quality_cirrus(IX,IY)=0
         Land_CLDMSK_forfraction(IX,IY)=0
        ENDDO
      ENDDO

C
C  All pixels clear(option 1 of cloud mask) ?????
C
C
C      RMED(1)=0
C      RMEDSQ(1)=0

C
C Cloud mask based upon spatial variability of 0.47 micron 500 m
C resolution reflectance data
C

 
         CALL CldMsk_3by3_HKM(ISWATH,ILINE,Refl_3,
     &cloud_threhold_land47_1,NROWS,RMED,RMEDSQ,CldMsk_500,data_size)
 
c reflactance test for 0.47 um
        DO IY=1,Data_Size(2)*2
         DO IX=1,Data_Size(1)*2
      if(Refl_3(ix,iy) .gt.cloud_threhold_land47_2)CldMsk_500(IX,IY)=0
         Enddo
       Enddo

      CALL CldMsk_3by3_1KM(Data_Size,ISWATH,ILINE,Refl_26,THRHLD1380_1,
     *THRHLD1380_2,CldMsk_1km,RMED_1km,RMEDSQ_1km,refl_5,quality_cirrus,
     * Land_CLDMSK_forfraction)
     
      DO IY=1,Data_Size(2)
       DO IX=1,Data_Size(1)
          Y2_offset = (IY-1)*2
          X2_offset = (IX-1)*2 
          DO l = 1,2
          DO p = 1,2 
            l2 = Y2_offset + l
            p2 = X2_offset + p
           IF(CldMsk_1km(IX,IY).EQ.1 .and.CldMsk_500(p2,l2).EQ.1)THEN
            CldMsk_500(p2,l2)=1
           ELSE
            CldMsk_500(p2,l2)=0
            ENDIF
          ENDDO
          ENDDO
          Y4_offset = (IY-1)*4
          X4_offset = (IX-1)*4
          N=0
          DO l = 1,2
          DO p = 1,2
            l2 = Y2_offset + l
            p2 = X2_offset + p
            IF(CldMsk_500(p2,l2).EQ.1) THEN
            N=N+1
            ENDIF
          ENDDO
          ENDDO

C
C  N=4 to overwrite 1.38 micro channel cloud screening results
C

          IF(N.EQ.4) THEN 

            DO ll = 1,4
            DO pp = 1,4

              ll2 = Y4_offset + ll
              pp2 = X4_offset + pp

              CldMsk_250(pp2,ll2) = 1

            ENDDO
            ENDDO
  
          ENDIF 

      ENDDO
      ENDDO
      
      RETURN
      END
C***************************************************************************
     
C***************************************************************************
        SUBROUTINE CldMsk_3by3_HKM( ISWATH,ILINE,REF1KM,THRHLD,
     *                                NROWS,RMED,RMEDSQ,CLDMSK,data_size )
C-----------------------------------------------------------------------
C!F77
C
C!DESCRIPTION:
C              This subroutine derive cloud mask using 3 x 3 spatial 
C              variability based upon MODIS 500 m resolution data
C
C!INPUT PARAMETERS:
C
C             Data_Size   Number of pixels in along and against scan
C                ISWATH   Number of pixels at 1 km resolution along scan
C                 ILINE   Number of pixels at 1 km resolution against scan
C                REF1KM   Reflectance at 500 meter resolution
C               THRHLD1   Threshold1 of 3x3 spatial variability
C               THRHLD2   Threshold2 of pixel reflectance value
C
C!OUTPUT PARAMETERS:
C
C                CLDMSK  Cloud mask
C                  RMED  Mean
C                RMEDSQ  Standard deviation
C
C!REVISION HISTORY:
C
C!TEAM-UNIQUE HEADER:
C
C Developed by MODIS Aerosol Team at NASA GSFC, Greenbelt, MD
C
C!DESIGN NOTES:
C Edited 21 Dec 2011 by Leigh Munchak.
C Implemented overwrite of cloud mask based on standard deviation
C threshold of .0075 - only called if cloud mask is considered
C cloudy. This keeps bright but spatially homogenous smoke plumes
C from being classified as cloudy.
C 
C!END
C-----------------------------------------------------------------------
C
        IMPLICIT NONE
        INTEGER IY,JY,JX,NY,ISWATH,ILINE,Data_Size(2),ISTART
        integer IEND
        INTEGER NROWS(ISWATH*2),CLDMSK(ISWATH*2,ILINE*2),add_numcld
        REAL REF1KM(ISWATH*2,ILINE*2),THRHLD
        REAL RMED(ISWATH*2),RMEDSQ(ISWATH*2),STD,STD2,VAR
        SAVE
        
C     
C     Initialize work arrays 
        DO JX = 1, data_size(1)*2
           NROWS(JX) = 0
           RMED(JX) = 0
           RMEDSQ(JX) = 0
        ENDDO
        
C     
C     Checking 500 meter resolution 
        
        NY=0
        DO IY=2,ILINE*2-1
           ISTART=IY-1
           IEND=ISTART+2
           
         DO JY=ISTART,IEND
            NY=NY+1 
            
            DO JX=2,data_size(1)*2-1
               
c     IF(REF1KM(JX,JY).GT.0.0.AND.REF1KM(JX+1,JY).GT.0.0.AND.REF1KM(JX-1,JY).GT.0.0) THEN
               NROWS(JX)=NROWS(JX) + 1
               RMED(JX)=RMED(JX)+(REF1KM(JX,JY)+REF1KM(JX+1,JY)
     1              +REF1KM(JX-1,JY))
               
               RMEDSQ(JX)=RMEDSQ(JX)+(REF1KM(JX,JY)*REF1KM(JX,JY)
     1              + REF1KM(JX+1,JY)*REF1KM(JX+1,JY)
     2              + REF1KM(JX-1,JY)*REF1KM(JX-1,JY))
c     ENDIF
               
C...........make clear determination where possible and re-initialize work array elements
               IF(NY.EQ.3) THEN
                  
C..............make clear/cloud determination only when all 9 pixels in 3x3 array are valid 
                  IF(NROWS(JX) .EQ. 3) THEN
                     VAR=9.0/8.0*(RMEDSQ(JX)/9.0-RMED(JX)*RMED(JX)/81.0)
                     
                     IF(VAR.GT.0.0) THEN
                        STD=SQRT(VAR)
c     New definitation.....
                        STD2=STD* (RMED(JX)/3.0)

C     New cloud masking in collection 6
                        IF(STD2 .LT. THRHLD) THEN                           
                           CLDMSK(JX,IY)=1                           
                        ELSE                         
                           IF(STD .LT. .0075) THEN
                              CLDMSK(JX,IY)=1
                           ELSE
                              CLDMSK(JX,IY)=0
                           ENDIF
                        ENDIF
                     ENDIF
                  ENDIF
                  
                  NROWS(JX)=0
                  RMED(JX)=0.0
                  RMEDSQ(JX)=0.0
               ENDIF
               
            ENDDO
         ENDDO
         
         NY=0
         
      ENDDO 
      
      DO JX=1,data_size(1)*2
         CLDMSK(JX,1)=CLDMSK(JX,2)
         CLDMSK(JX,ILINE*2)=CLDMSK(JX,ILINE*2-1)
      ENDDO
      
      DO JY=1,ILINE*2
         CLDMSK(1,JY)=CLDMSK(2,JY)
         CLDMSK(data_size(1)*2,JY)=CLDMSK(data_size(1)*2-1,JY)
      ENDDO
      RETURN
      END

C***************************************************************************
      SUBROUTINE CldMsk_3by3_1KM(Data_Size,ISWATH,ILINE,REF1KM,THRHLD1,
     *THRHLD2,CLDMSK,RMED,RMEDSQ,refl_5, Quality_cirrus, 
     *Land_CLDMSK_forfraction)
C-----------------------------------------------------------------------
C!F77
C
C!DESCRIPTION:
C              This subroutine derive cloud mask using 3 x 3 spatial 
C              variability based upon MODIS 500 m resolution data
C
C!INPUT PARAMETERS:
C
C             Data_Size   Number of pixels in along and against scan
C                ISWATH   Number of pixels at 1 km resolution along scan
C                 ILINE   Number of pixels at 1 km resolution against scan
C                REF1KM   Reflectance at 500 meter resolution
C               THRHLD1   Threshold1 of 3x3 spatial variability
C               THRHLD2   Threshold2 of pixel reflectance value
C
C!OUTPUT PARAMETERS:
C
C                CLDMSK  Cloud mask
C                  RMED  Mean
C                RMEDSQ  Standard deviation
C
C!REVISION HISTORY:
C
C!TEAM-UNIQUE HEADER:
C
C Developed by MODIS Aerosol Team at NASA GSFC, Greenbelt, MD
C
C!DESIGN NOTES:
C
C!END
C-----------------------------------------------------------------------
C
      IMPLICIT NONE
      INTEGER IY,JY,NY,JX,ISWATH,ILINE,Data_Size(2),ISTART,IEND
      INTEGER CLDMSK(ISWATH,ILINE), Quality_cirrus(ISWATH,ILINE)
      integer Land_CLDMSK_forfraction(ISWATH,ILINE)
      REAL REF1KM(ISWATH,ILINE),THRHLD1,THRHLD2
      REAL RMED(ISWATH),RMEDSQ(ISWATH),STD,VAR
      INTEGER Y2_offset,X2_offset,l2,p2,l,p
      Real Refl_5(ISWATH*2,ILINE*2)
      SAVE

C
C Checking for 1.38 micron spatial variability
C

      NY=0

      DO IY=2,Data_Size(2)-1
        ISTART=IY-1
        IEND=ISTART+2

        DO JY=ISTART,IEND
          NY=NY+1 

          DO JX=2,Data_Size(1)-1

c            IF(REF1KM(JX,JY).GT.0.0.AND.REF1KM(JX+1,JY).GT.0.0.AND.REF1KM(JX-1,JY).GT.0.0) THEN

              RMED(JX)=RMED(JX)+(REF1KM(JX,JY)+REF1KM(JX+1,JY)+REF1KM(JX-1,JY))
              RMEDSQ(JX)=RMEDSQ(JX)+(REF1KM(JX,JY)*REF1KM(JX,JY)
     1                             + REF1KM(JX+1,JY)*REF1KM(JX+1,JY)
     2                             + REF1KM(JX-1,JY)*REF1KM(JX-1,JY))

              IF(NY.EQ.3) THEN

                VAR=9.0/8.0*(RMEDSQ(JX)/9.0-RMED(JX)*RMED(JX)/81.0)
                IF(VAR.GT.0.0) THEN
                  STD=SQRT(VAR) 
                  
                     IF(STD.LT.THRHLD1) THEN
                    CLDMSK(JX,IY)=1
                  ENDIF
                ENDIF

                RMED(JX)=0.0
                RMEDSQ(JX)=0.0

               ENDIF

c            ENDIF

          ENDDO
        ENDDO

        NY=0

      ENDDO 

      DO JX=1,Data_Size(1)
        CLDMSK(JX,1)=CLDMSK(JX,2)
        CLDMSK(JX,Data_size(2))=CLDMSK(JX,Data_Size(2)-1)
      ENDDO
   
      DO JY=1,Data_Size(2)
        CLDMSK(1,JY)=CLDMSK(2,JY)
        CLDMSK(Data_Size(1),JY)=CLDMSK(Data_Size(1)-1,JY)
      ENDDO

C
C Checking for 1.38 micron reflectance
C
  
      DO JY=1,Data_Size(2)
        DO JX=1,Data_Size(1) 
       Land_CLDMSK_forfraction(JX,JY)=CLDMSK(JX,JY)
      IF(REF1KM(JX,JY).gt.0.025) then
             CLDMSK(JX,JY)=0
c save for cloud fraction
       Land_CLDMSK_forfraction(JX,JY)=CLDMSK(JX,JY)
      Else
c if Variabily cloud mask says cloud free and smoke or cirrus reduce quality
         IF(CLDMSK(JX,JY) .eq.1 .and.REF1KM(JX,JY).gt. 0.01 
     *    .and.  REF1KM(JX,JY) .le.0.025 ) then
c  cloud free  may be cirrus  or smoke , will be retrived but put a quality value  
              Quality_cirrus(JX,JY)=1  
      
       ENDIF
cENDIF   cloud free  may be cirrus  or smoke 
       ENDIF   
      
        ENDDO
      ENDDO

      RETURN
      END
C*********************************************************************
      SUBROUTINE GEOLOC_ANGLE(LAT,LON,SatZen,SatAz,SolZen,
     &SolAz,RelAz,Height,MTHET0,MTHET,MPHI0,MPHI,MDPHI,MSCATT,MHGHT,
     &IDATA,START_500,END_500,START_250,END_250,START_1KM,END_1KM,
     &Lat_center,Lon_center,iscan)

C---------------------------------------------------------------------
C!F77
C
C!DESCRIPTION:
C
C          This subroutine findes the averaged latitude, longitude, height
C          and geometrical angles at the center of 10x10 km^2 area
C
C!INPUT PARAMETERS:
C
C        LAT         Lat array for 10*10 box
C        LON         Lon array for 10*10 box
C        SatZen      Satellite zenith angle
C        SatAz       Satellite Azimuth
C        SolZen      Solar Zenith Angle
C        SolAz       Solar Azimuth
C        RelAz       Relative Azimuth
c        START_500   Starting  Index for 500 m resolution data
C        END_500     Ending    Index for 500 m resolution data
C        START_250   Starting  Index for 250 m resolution data
C        END_250     Ending    Index for 250 m resolution data
C        START_1KM   Starting  Index for 1 km  resolution data
C        END_1KM     Ending    Index for 1 km  resolution data
C
C!OUTPUT PARAMETERS:
C
C        Lat_center  Averaged latitude (in degree)
C        Lon_center  Averaged londitude (in degree)
C        MTHET0      Averaged solar zenith angle (in degree)
C        MTHET       Averaged satellite zenith angle (in degree)
C        MPHI0       Averaged solar azimuth angle (in degree)
C        MPHI        Averaged satellite azimuth angle (in degree)
C        MDPHI       Averaged relative azimuth angle (in degree)
C        MHGHT       Averaged topographic altitude (in km)
c
C!REVISION HISTORY:
c 10/15/1999 fhliang
c fixed prolog.
C
C!TEAM-UNIQUE HEADER:
C
c!END
C----------------------------------------------------------------------

      IMPLICIT NONE
      SAVE

      INCLUDE 'mod04.inc'
        REAL diff, sumdif
      INTEGER START_500,END_500,START_250,END_250,START_1KM,END_1KM
C     INTEGER IX,IY,NUMDATA,IDATA,WATER,JMASK,IMASK,CLDFREE,II,JJ
      INTEGER IX1,IX2,IY1,IY2,IDATA,I,iscan
      REAL Lat_center,Lon_center,Lon_Min,Lon_Max,Lon_4(4)
      REAL MTHET0,MTHET,MPHI0,MPHI,MDPHI,MSCATT,MHGHT,AVE1,AVE2
      REAL Lat(ISWATH,ILINE),Lon(ISWATH,ILINE),
     2     SatZen(ISWATH,ILINE),SolZen(ISWATH,ILINE),
     3     SatAz(ISWATH,ILINE),SolAz(ISWATH,ILINE),
     4     RelAz(ISWATH,ILINE),Height(ISWATH,ILINE)
       Real SolAz_4(4),SolAz_Min,SolAz_Max
       Real SatAz_4(4),SatAz_Min,SatAz_Max

! since it is 3* 3, to keep the software consistant taking  centeral point  
      IX1=START_1KM+1
      IX2=START_1KM+1
      IY1=IGRIDY/2+1
      IY2=IGRIDY/2+1
      
       
     
C
C Finding minimum and maximum of longitudes
C

      Lon_4(1)=LON(IX1,IY1)
      Lon_4(2)=LON(IX2,IY1)
      Lon_4(3)=LON(IX1,IY2)
      Lon_4(4)=LON(IX2,IY2)
      Lon_Min=Lon_4(1)
      Lon_Max=Lon_4(1)
      
      SolAz_4(1)=SolAz(IX1,IY1)
      SolAz_4(2)=SolAz(IX2,IY1) 
      SolAz_4(3)=SolAz(IX1,IY2)
      SolAz_4(4)=SolAz(IX2,IY2) 
      SolAz_Min=SolAz_4(1) 
      SolAz_Max=SolAz_4(1)
      
      SatAz_4(1)=SatAz(IX1,IY1)
      SatAz_4(2)=SatAz(IX2,IY1) 
      SatAz_4(3)=SatAz(IX1,IY2)
      SatAz_4(4)=SatAz(IX2,IY2) 
      SatAz_Min=SatAz_4(1) 
      SatAz_Max=SatAz_4(1)

      DO I=1,4
        IF(Lon_4(I).LE.Lon_Min) Lon_Min=Lon_4(I)
        IF(Lon_4(I).GE.Lon_Max) Lon_Max=Lon_4(I)
        IF(SolAz_4(I).LE.SolAz_Min) SolAz_Min=SolAz_4(I)
        IF(SolAz_4(I).GE.SolAz_Max) SolAz_Max=SolAz_4(I)
        IF(SatAz_4(I).LE.SatAz_Min) SatAz_Min=SatAz_4(I)
        IF(SatAz_4(I).GE.SatAz_Max) SatAz_Max=SatAz_4(I)
      ENDDO

C
C If Lon_Max <-180 means fill values are found and if Lon_Min <-180 means at least one'
C fill value is found, under these two condition fill value is set to longitude,
C latitude, and all the angles
C
      IF(Lon_Max.LT.-180.0.OR.Lon_Min.LT.-180.0) THEN

        Lon_center=-999.0
        Lat_center=-999.0
        MTHET0=-99.99
        MTHET=-99.99
        MPHI0=-99.99
        MPHI=-99.99
        MSCATT=-99.99

      ELSE

C
C Otherwise, check for other condisitons that may contain fill value
C
      sumdif = 0.

c.....take longitude difference of pixels 2 through 4 relative to pixel 1
      do i = 2,4
          diff = Lon_4(i) - Lon_4(1)

c..........we are working on a sphere, take shortest arc between points
          if ( diff  .gt.  180. ) diff = diff   - 360.
          if ( diff  .lt.  -180. ) diff = 360. + diff

           sumdif = sumdif + diff
      enddo

c.....mathematically equivalent to {Lon_4(1)+Lon_4(2)+Lon_4(3)+Lon_4(4)} / 4
c    using shortest arcs
      lon_center = Lon_4(1) + sumdif/4.

      if ( lon_center  .gt. +180. ) lon_center = lon_center  - 360.
      if ( lon_center  .lt.   -180. ) lon_center = lon_center + 360.

       
             
         
           

      Lat_center=(LAT(IX1,IY1)+LAT(IX2,IY1)+
     *            LAT(IX1,IY2)+LAT(IX2,IY2))/4.
         MTHET0=(SolZen(IX1,IY1)+SolZen(IX2,IY1)+
     *        SolZen(IX1,IY2)+SolZen(IX2,IY2))/4.

      MTHET =(SatZen(IX1,IY1)+SatZen(IX2,IY1)+
     *        SatZen(IX1,IY2)+SatZen(IX2,IY2))/4.
      
      MHGHT =(Height(IX1,IY1)+Height(IX2,IY1)+
     *        Height(IX1,IY2)+Height(IX2,IY2))/4000.

   
C Otherwise, check for other condisitons that may contain fill value
C
      sumdif = 0.

c.....take longitude difference of pixels 2 through 4 relative to pixel 1
      do i = 2,4
          diff = SolAz_4(i) - SolAz_4(1)

c..........we are working on a sphere, take shortest arc between points
          if ( diff  .gt.  180. ) diff = diff - 360.
          if ( diff  .lt.  -180. ) diff = 360. + diff

           sumdif = sumdif + diff
      enddo

c..... 
c    using shortest arcs

        MPHI0 = SolAz_4(1) + sumdif/4.

      if ( MPHI0  .gt. +180. ) MPHI0 = MPHI0  - 360.
      if ( MPHI0  .lt. -180. ) MPHI0 = MPHI0 + 360. 
        
      
C
      sumdif = 0.

c.....
      do i = 2,4
          diff = SatAz_4(i) - SatAz_4(1)

c..........we are working on a sphere, take shortest arc between points
          if ( diff  .gt.  180. ) diff = diff - 360.
          if ( diff  .lt.  -180. ) diff = 360. + diff

           sumdif = sumdif + diff
      enddo

c..... 
c    using shortest arcs

        MPHI = SatAz_4(1) + sumdif/4.

      if (  MPHI  .gt. +180. )  MPHI =  MPHI  - 360.
      if (  MPHI  .lt. -180. )  MPHI =  MPHI +  360.
      
      
c          MPHI0 =(SolAz(IX1,IY1)+SolAz(IX2,IY1)+
c     *    SolAz(IX1,IY2)+SolAz(IX2,IY2))/4.
c           MPHI  =(satAz(IX1,IY1)+SatAz(IX2,IY1)+
c     *        SatAz(IX1,IY2)+SatAz(IX2,IY2))/4.
     
      
 
       MDPHI = abs(MPHI0 - MPHI -180.0)
       IF(MDPHI.GT.360.0) MDPHI=amod(MDPHI,360.0)
       IF(MDPHI.GT.180.0) MDPHI=360.0-MDPHI

      IF(MTHET0.GT.0.0.AND.MTHET.GT.0.0.AND.MDPHI.GT.0.0) THEN
        MSCATT = -COS(MTHET0*DTR)*COS(MTHET*DTR)
     *          +SIN(MTHET0*DTR)*SIN(MTHET*DTR)
     *          *COS(MDPHI*DTR)
        MSCATT = ACOS(MSCATT)*RTD
      ELSE
        MSCATT=-99.99
      ENDIF 
C endif for lon_max.........     
      ENDIF 
      RETURN
      END
        SUBROUTINE FILLVALUE_missing(SDS_ref,SDS_ref_STD,SDSTAU_best,SDSTAUS_best,SDSTAUB_best,
     *SDSTAU_average,SDSTAUS_average,SDSTAUB_average,SDS_Least_error,SDS_small_weighting,
     *SDS_sol_INDX_small,SDS_sol_INDX_large,SDSASSY_best,SDSASSY_average,SDS_ccn,sds_mass_conc,
     *SDSBACK_best,SDSBACK_average,SDS_effrad,SDS_AOT_model,SDS_RefF_best,SDS_RefF_average,
     *SDS_TranF_best,SDS_TranF_average,SDS_angs_coeff1,SDS_angs_coeff2,SDS_SCAT_ANGLE_OCEAN,
     *SDS_QCONTROL_ocean,SDS_NUMPIXELS,SDS_CLDFRC_ocean,SDS_Tau_Land_Ocean_img,Qcontrol_special,
     *SDS_correc_small_weighting,SDS_Aerosol_Type,SDSTAU_corrected_213,SDS_SCAT_ANGLE_land,
     *SDS_mass_conc_land,SDS_angs_coeff_land,SDS_CLDFRC_land,SDS_dust_weighting,SDS_est_uncer,
     *SDS_RefF_land,SDS_TranF_land,SDS_NUMPIXELS_land,SDSTAU_corrected,SDS_ref_land,
     *SDS_ref_STD_land,SDS_QCONTROL_land,SDSTAU_small_land,SDS_Surface_Reflectance_Land,
     *SDS_Fitting_Error_Land,Qcontrol_special_land,index_wave,SDS_Mean_Reflectance_Land_All,
     *SDS_SDev_Reflectance_Land_All,SDS_Path_Radiance_Land,SDS_Critical_Reflectance_Land,  
     *SDS_Error_Crit_Reflectance_Land,SDS_Error_Path_Radiance_Land,  
     * SDS_QFlag_Critical_Ref_Land,SDS_QFlag_Path_Radiance_Land,SDSLAT,
     *SDSLON,SDS_SST,SDS_MTHET0,SDS_MTHET,SDS_MPHI0,SDS_MPHI,
     *SDS_Scattering_Angle,SDS_CldMskQA,SDS_Tau_Land_Ocean,
     *SDS_ratio_small_Land_Ocean,NumSQ,SDS_Sea_Land_Flag,iscan)

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION: This subroutine Fills the HDF array's with fill values.
C
C !INPUT PARAMETERS:
C IDATA        Index of box number
C !OUTPUT PARAMETERS:
C SDS_ref         Averaged Reflectance array
C SDS_ref_STD     Standard Deviation of Reflectance
C SDSTAU_best     Optical thickness for best solution
C SDSTAUS_best    Optical thickness contribution small particles for best solution
C SDSTAUB_best    Optical thickness contribution large particles for best solution
C SDSTAU_average  Optical thickness for best solution
C SDSTAUS_average Optical thickness contribution small particles for best solution
C SDSTAUB_average Optical thickness contribution large particles for best solution
C SDS_Least_error    Minimm Error function betwwen derived and computed radiances
C SDS_small_weighting Weight factor for large and small mode
C SDS_sol_INDX       Index for solution number
C SDSASSY_best       Assymetry Factor for best solution
C SDSASSY_average    Assymetry Factor for Average solution
C SDSBACK_best       Backscattering Ratio for best solution
C SDSBACK_average    Backscattering Ratio for Average solution
C SDS_effrad         Effective Radiance
C SDS_AOT_model Ration of optical Thickess small
C SDS_RefF_best      Reflected Flux Best solution
C SDS_RefF_average   Reflected Flux Average solution
C SDS_TranF_best     Transmitted Flux Best solution
C SDS_TranF_average  Transmitted Flux Average solution
C QCONTROL           Value for Quality Control
C SDS_SCAT_ANGLE_OCEAN Scattering angle ocean
C SDS_QCONTROL         Quality control SDS array
C SDS_NUMPIXELS        Number of Pixels used
c 
C!DESCRIPTION:     This subroutine stores all the arrays to be written
C                  as output(HDF) file
C
C!INPUT PARAMETERS: all varaiables to be wrriten as output
C   LATM          Latitude of data cell
C   LONM          Longitude of data cell
C   IAER          Aerosol type
C   BARLBLUE      Average reflectance for blue channel
C   IDATA         Data cell number from 1 to NUMDATA
C   NUMDATA       Number of retrieval cells across orbit swath
C   TAUABLUE      Blue channel aerosol optical thickness (Continental)
C   BARLRED       Average reflectance for red channel
C   TAUARED       Red channel aerosol optical thickness (Continental)
C   SDBLUE        STD of blue channel reflectances
C   SDRED         STD of red channel reflectances
C   TAUABLUEC     Blue channel aerosol optical thickness (Corrected)
C   TAUAREDC      Red channel aerosol optical thickness (Corrected)
C   IBLUE         Number of blue channel observations
C   IRED          Number of red channel observations
C   IERROR        Error flag (0-4)
C   WEIGHT        Relative contribution of smoke/sulfate particles to
C                 dust in the computation of the aerosol optical depth
C   ISULP         Sulfate (1)/no sulfate flag (0)
C   ISMK          Smoke (1)/no smoke flag (0)
C   IPROCE        Aerosol retrieval procedure ID (0-4)
C   SAVERATIO     Aerosol path radiance ration (red to blue channel
C                 for continental model)
C   PIXELS        Number of cells across swath (same as NUMDATA)
C   LINES         Number of cells along swath (same as NUMSCAN)
C   NUMXBOX       Not used
C   NUMYBOX       Not used
C
C!OUTPUT PARAMETERS:ARRAY of 13 variables for output to HDF FILE
C  SDS1           HDF array of cell latitudes
C  SDS2           HDF array of cell longitudes
C  SDS3           HDF array of spectral reflectances
C  SDS4           HDF array of aerosol optical thicknesses for
C                     continental model
C  SDS5           HDF array of the STD of spectral reflectances
C  SDS6           HDF array of aerosol optical thicknesses for
C                     corrected model
C  SDS7           HDF array of STD for corrected optical thickness
C                     blue channel for continental model)
C  SDS8           HDF array of aerosol path radiance ratio (red to
C                     blue channel for continental model)
C  SDS9           HDF array of relative aerosol optical depth
C                     (smoke/sulfate particles to dust)
C  SDS10           HDF array of the number of blue channel clear pixels
C  SDS11          HDF array of the number of red channel clear pixels
C  SDS12          HDF array of retrieval procedure ID (0-4)
C  SDS13          HDF array of aerosol type (0-3)
C  SDS14          HDF array of Error flag (0-4)
C
C!REVISION HISTORY: Updated code to comply with most MODIS software
C                   standards.
C
C!TEAM-UNIQUE HEADER: - NOT YET DEFINED
C
C!REFERENCES AND CREDITS
C     WRITTEN BY: Shana Mattoo
C
C!DESIGN NOTES:
C
C!END
C----------------------------------------------------------------------
C
C !REVISION HISTORY:
C 10/15/1999 fhliang
C fixed prolog.
C
C !TEAM-UNIQUE HEADER: - NOT YET DEFINED
C
C !REFERENCES AND CREDITS
C      WRITTEN BY: Shana Mattoo
C
C !DESIGN NOTES:
C
C !END
C----------------------------------------------------------------------

       IMPLICIT NONE
       SAVE

       INCLUDE 'mod04.inc'

       BYTE    SDS_QCONTROL_ocean(QA_LAND,NUMCELLS)
       INTEGER  Numsq,ICASE,IWAV,idata
       INTEGER*2  SDS_ref(NUMCELLS,NWAV),SDS_ref_STD(NUMCELLS,NWAV)
       INTEGER*2  SDSTAU_best(NUMCELLS,NWAV),
     &            SDSTAU_average(NUMCELLS,NWAV),
     &            SDSTAUB_best(NUMCELLS,NWAV),
     &            SDSTAUB_average(NUMCELLS,NWAV),
     &            SDSTAUS_best(NUMCELLS,NWAV),
     &            SDSTAUS_average(NUMCELLS,NWAV),
     &            SDSASSY_best(NUMCELLS,NWAV),
     &            SDSASSY_average(NUMCELLS,NWAV),
     &            SDSBACK_best(NUMCELLS,NWAV),
     &            SDSBACK_average(NUMCELLS,NWAV),
     &            SDS_RefF_best(NUMCELLS,NWAV),
     &            SDS_RefF_average(NUMCELLS,NWAV),
     &            SDS_TranF_best(NUMCELLS,NWAV),
     &            SDS_TranF_average(NUMCELLS,NWAV) 

       INTEGER*2  SDS_small_weighting(NUMCELLS,NUM_solutions),
     &            SDS_correc_small_weighting(NUMCELLS),
     &            SDS_Least_error(NUMCELLS,NUM_solutions),
     &            SDS_effrad(NUMCELLS,NUM_solutions),
     &            SDS_sol_INDX_small(NUMCELLS,NUM_solutions),
     &            SDS_sol_INDX_large(NUMCELLS,NUM_solutions),
     &            SDS_angs_coeff1(NUMCELLS,NUM_solutions),
     &            SDS_angs_coeff2(NUMCELLS,NUM_solutions),
     &            SDS_AOT_model(NUMCELLS,num_model_index),
     &            SDS_CLDFRC_ocean(NUMCELLS),
     &            SDS_Tau_Land_Ocean_img(NUMCELLS),
     &             SDS_Tau_Land_Ocean(NUMCELLS),
     &             SDS_Sea_Land_Flag(ISWATH,ILINE*Tot_scan)
      BYTE         SDS_QCONTROL_land(QA_LAND,NUMCELLS)
      REAL         SDS_mass_conc_land(NUMCELLS)
      INTEGER *2   SDS_Aerosol_Type(NUMCELLS),
     &             SDS_SCAT_ANGLE_land(NUMCELLS),
     &             SDS_angs_coeff_land(NUMCELLS),
     &             SDS_CLDFRC_land(NUMCELLS),
     &             SDS_dust_weighting(NUMCELLS), 
     &             SDS_NUMPIXELS_land(NUMCELLS,Land_Sol1),
     &             SDSTAU_corrected(NUMCELLS,Land_Sol3),
     &             SDS_ref_land(NUMCELLS,Band_land),
     &             SDS_ref_STD_land(NUMCELLS,Band_land), 
     &     SDS_Surface_Reflectance_Land(NUMCELLS,Land_Sol3),
     &     SDS_Fitting_Error_Land(NUMCELLS),
     &     SDSTAU_corrected_213(NUMCELLS),
     &     SDSTAU_small_land(NUMCELLS,Land_Sol4),
     &     SDS_ratio_small_Land_Ocean(NUMCELLS)
       INTEGER index_wave,startwav,endwav
       INTEGER  INT_Fill_value
       INTEGER *2  SDS_Mean_Reflectance_Land_All(NUMCELLS,Land_Sol3), 
     &SDS_SDev_Reflectance_Land_All(NUMCELLS,Land_Sol3),     
     &SDS_Path_Radiance_Land(NUMCELLS,Land_Sol1),   
     & SDS_Critical_Reflectance_Land(NUMCELLS,Land_Sol1),  
     & SDS_Error_Crit_Reflectance_Land(NUMCELLS,Land_Sol1),      
     & SDS_Error_Path_Radiance_Land(NUMCELLS,Land_Sol1),  
     & SDS_QFlag_Critical_Ref_Land(NUMCELLS,Land_Sol1) ,
     & SDS_QFlag_Path_Radiance_Land(NUMCELLS,Land_Sol1)
      REAL SDSLAT(NUMCELLS),SDSLON(NUMCELLS)
      INTEGER*2 SDS_MTHET0(NUMCELLS),SDS_MTHET(NUMCELLS),
     &          SDS_MPHI0(NUMCELLS),SDS_MPHI(NUMCELLS),
     &          SDS_Scattering_Angle(NUMCELLS)
      REAL*8    SDS_SST(NUMCELLS),SST
      BYTE      SDS_CldMskQA(NUMCELLS)

 
      INTEGER *2   
     &            SDS_est_uncer(NUMCELLS,Land_Sol1),
     &            SDS_RefF_land(NUMCELLS,Land_Sol2),
     &            SDS_TranF_land(NUMCELLS,Land_Sol1)


      INTEGER  Qcontrol_special_land
      REAL     FLOAT_Fill_value
       REAL SDS_mass_conc(NUMCELLS,NUM_solutions),
     &      SDS_ccn(NUMCELLS,NUM_solutions)
        INTEGER  Qcontrol_special,iscan
       INTEGER*2 SDS_NUMPIXELS(NUMCELLS),
     &   SDS_SCAT_ANGLE_OCEAN(NUMCELLS)
c       
          INT_Fill_value=-9999
          FLOAT_Fill_value=-999.0
          
           
        Do Idata = 1,NumSq
          SDS_NUMPIXELS(idata)=INT_Fill_value
          SDS_SCAT_ANGLE_OCEAN(IDATA)=INT_Fill_value 
          SDS_CLDFRC_ocean(IDATA)=INT_Fill_value  
           SDS_Tau_Land_Ocean(IDATA)=INT_Fill_value 
           SDS_Tau_Land_Ocean_img(IDATA)=INT_Fill_value 
           SDS_ratio_small_Land_Ocean(IDATA)=INT_Fill_value 
           SDS_Sea_Land_Flag(idata,iscan) =INT_Fill_value
          DO IWAV= 1,NWAV
            SDS_ref(IDATA,IWAV) =INT_Fill_value
            SDS_ref_STD(IDATA,IWAV) =INT_Fill_value
            SDSASSY_best(IDATA,IWAV) =INT_Fill_value
            SDSASSY_average(IDATA,IWAV) =INT_Fill_value
            SDSBACK_BEST(IDATA,IWAV) =INT_Fill_value
            SDSBACK_AVERAGE(IDATA,IWAV) =INT_Fill_value 
            SDSTAU_BEST(IDATA,IWAV) = INT_Fill_value
            SDSTAUS_BEST(IDATA,IWAV) = INT_Fill_value
            SDSTAUB_BEST(IDATA,IWAV) = INT_Fill_value
            SDSTAU_AVERAGE(IDATA,IWAV) = INT_Fill_value
            SDSTAUS_AVERAGE(IDATA,IWAV) = INT_Fill_value
            SDSTAUB_AVERAGE(IDATA,IWAV) = INT_Fill_value
          ENDDO
            
c  Fill with Fill values in all cases
           DO icase=1,NUM_solutions
            SDS_angs_coeff1(IDATA,ICASE) =INT_Fill_value
            SDS_angs_coeff2(IDATA,ICASE) =INT_Fill_value
            SDS_effrad(IDATA,ICASE) =INT_Fill_value
            SDS_Least_error(IDATA,ICASE) = INT_Fill_value
            SDS_small_weighting(IDATA,ICASE) = INT_Fill_value
            SDS_sol_INDX_small(IDATA,ICASE) = INT_Fill_value
            SDS_sol_INDX_large(IDATA,ICASE) = INT_Fill_value
            SDS_ccn(IDATA,ICASE) =FLOAT_Fill_value
            sds_mass_conc(IDATA,ICASE)= FLOAT_Fill_value
          ENDDO 
             DO icase=1,num_model_index
             SDS_AOT_model(IDATA,ICASE) = INT_Fill_value
             Enddo
        ENDDO
   
 
C
C Land
C

      

C Set Fill_Values to SDS arrays
C
      
         Do Idata=1,numsq
         SDS_CLDFRC_land(IDATA)=INT_Fill_value
         SDS_mass_conc_land(IDATA)= FLOAT_Fill_value
        SDS_QCONTROL_land(1,IDATA)=0
        SDS_QCONTROL_land(2,IDATA)=0
        SDS_QCONTROL_land(3,IDATA)=0
        SDS_QCONTROL_land(4,IDATA)=0
        SDS_QCONTROL_land(5,IDATA)=0
         SDS_Aerosol_Type(IDATA)=INT_Fill_value
         SDS_SCAT_ANGLE_land(IDATA)=INT_Fill_value
         SDSTAU_corrected(IDATA,1)=INT_Fill_value
         SDSTAU_corrected(IDATA,2)=INT_Fill_value
         SDSTAU_corrected(IDATA,3)=INT_Fill_value
         SDSTAU_corrected_213(IDATA)=INT_Fill_value
         SDSTAU_small_land(IDATA,1)=INT_Fill_value
         SDSTAU_small_land(IDATA,2)=INT_Fill_value
         SDSTAU_small_land(IDATA,3)=INT_Fill_value
         SDSTAU_small_land(IDATA,4)=INT_Fill_value
         SDS_angs_coeff_land(IDATA)=INT_Fill_value
c         SDS_CLDFRC_land(IDATA)=INT_Fill_value
         SDS_dust_weighting(IDATA)=INT_Fill_value
         SDS_NUMPIXELS_land(IDATA,1)=INT_Fill_value
         SDS_NUMPIXELS_land(IDATA,2)=INT_Fill_value
         SDS_ref_land(IDATA,1) =INT_Fill_value
         SDS_ref_land(IDATA,2) =INT_Fill_value
         SDS_ref_land(IDATA,3) =INT_Fill_value
         SDS_ref_land(IDATA,4) =INT_Fill_value
         SDS_ref_land(IDATA,5) =INT_Fill_value
         SDS_ref_land(IDATA,6) =INT_Fill_value
         SDS_ref_land(IDATA,7) =INT_Fill_value
         SDS_ref_STD_land(IDATA,1)=INT_Fill_value
         SDS_ref_STD_land(IDATA,2)=INT_Fill_value
         SDS_ref_STD_land(IDATA,3)=INT_Fill_value
         SDS_ref_STD_land(IDATA,4)=INT_Fill_value
         SDS_ref_STD_land(IDATA,5)=INT_Fill_value
         SDS_ref_STD_land(IDATA,6)=INT_Fill_value
         SDS_ref_STD_land(IDATA,7)=INT_Fill_value
         SDS_Surface_Reflectance_Land(IDATA,1) =INT_Fill_value
         SDS_Surface_Reflectance_Land(IDATA,2) =INT_Fill_value
         SDS_Surface_Reflectance_Land(IDATA,3) =INT_Fill_value 
         SDS_Fitting_Error_Land(IDATA)=INT_Fill_value 
          Enddo
        
        Do Idata=1,numsq   
       SDSLAT(IDATA) = FLOAT_Fill_value
      SDSLON(IDATA) = FLOAT_Fill_value
      SDS_SST(IDATA)= INT_Fill_value
      SDS_MTHET0(IDATA) =INT_Fill_value
      SDS_MTHET(IDATA) = INT_Fill_value
      SDS_MPHI0(IDATA) =INT_Fill_value
      SDS_MPHI(IDATA) = INT_Fill_value
      SDS_Scattering_Angle(IDATA)=INT_Fill_value
      SDS_CldMskQA(IDATA) = FLOAT_Fill_value
      Enddo 
       
        RETURN
         END
C================================================================================
C
C SUBROUTINE interpol_wspeed_spatial
C
C This routine reads in the wind speeds of each grid cell surrounding the 
C current grid cell, and bilinearly interpolates wind speed between grid cell.
C
C This routine is intended to eliminate the problem of having discontinuities
C in AOD at the edge of a 1x1 lat/lon grid cell due to sharp changes in the 
C wind speed and therefore wind speed correction.
C
C Bilinear interpolation copied from Numerical Recipes in FORTRAN 77 
C Formula 3.6.5 
C
C Introduced by L. Munchak on March 27, 2013. leigh.a.munchak@nasa.gov
C
C================================================================================

      SUBROUTINE interpol_wspeed_spatial(Lat_center,Lon_center,ugrd,vgrd)

      IMPLICIT NONE
      SAVE

      INCLUDE 'mod04.inc'

C     Arguments passed from main
      REAL ugrd,vgrd,lat_center, lon_center       

C     Dummy variables for met I don't want to read
      REAL temp_null, pwat_null, oz_null         

C     FLOOR of lat_center & lon_center
      REAL fl_lat, fl_lon                       

C     CEILING of lat_center & lon_center
      REAL ci_lon, ci_lat                        

C     u-wind component for upper lefthand corner, lower-right, etc.
      REAL ul_u, ur_u, ll_u, lr_u                
 
C     v-wind component for upper lefthand corner, lower-right, etc.
      REAL ul_v, ur_v, ll_v, lr_v

C     Weighting for bilnear interpolation
      REAL bi_t, bi_u         
     
C     Get integer part of lat/lon
      fl_lat = INT(lat_center)
      fl_lon = INT(lon_center) 

      IF(fl_lat .LT. 0 ) THEN
         fl_lat = fl_lat -1
      ENDIF

      IF(fl_lon .LT. 0) THEN
         fl_lon = fl_lon - 1
      ENDIF

C     Get ceiling of lat/lon
      ci_lat = fl_lat + 1
      ci_lon = fl_lon + 1      

C     Read in wind speeds at each corner, put temp/oz/pwat into dummy variables
      CALL GetAncData_PGE04(ci_lat,ci_lon,temp_null,ur_u,ur_v,pwat_null,oz_null)
      CALL GetAncData_PGE04(ci_lat,fl_lon,temp_null,ul_u,ul_v,pwat_null,oz_null)
      CALL GetAncData_PGE04(fl_lat,ci_lon,temp_null,lr_u,lr_v,pwat_null,oz_null)
      CALL GetAncData_PGE04(fl_lat,fl_lon,temp_null,ll_u,ll_v,pwat_null,oz_null)

C     Compute weighting factors for interpolation
      bi_t  = (lon_center-fl_lon)/(ci_lon-fl_lon)    
      bi_u  = (lat_center-fl_lat)/(ci_lat-fl_lat)   

C     Perform interpolation on u & v
      ugrd  = ((1-bi_t)*(1-bi_u)*ll_u)+(bi_t*(1-bi_u)*lr_u)
     * +(bi_t*bi_u*ur_u)+((1-bi_t)*bi_u*ul_u)
      vgrd  = ((1-bi_t)*(1-bi_u)*ll_v)+(bi_t*(1-bi_u)*lr_v)
     * +(bi_t*bi_u*ur_v)+((1-bi_t)*bi_u*ul_v) 

      RETURN
      END
      
C================================================================================
C
C Subroutine indx_wspeed
C
C================================================================================

       Subroutine indx_wspeed(ugrd,vgrd,Indx_wspeed1,Indx_wspeed2,
     *         WSPEED,Wind)
          IMPLICIT NONE
          SAVE

       INCLUDE 'mod04.inc'
        real ugrd,vgrd,WSPEED,Del,Wind(Lut_indx) 
         integer Indx_wspeed1,Indx_wspeed2,ij,inn 
         WSPEED =SQRT(ugrd*ugrd+vgrd*vgrd)  
           do ij=1,lut_indx 
          if( ij .eq.1)WIND(IJ)=2
          if( ij .eq.2)WIND(IJ)=6
          if( ij .eq.3)WIND(IJ)=10
          if( ij .eq.4)WIND(IJ)=14
         enddo 
        Indx_wspeed1=1
        Indx_wspeed2=Indx_wspeed1+1 
         do ij=1,lut_indx 
          if( WIND(ij) .LE. WSPEED)  then
          Indx_wspeed1=ij
          Indx_wspeed2=ij+1
          endif
        enddo 
        if(Indx_wspeed2 .gt. lut_indx)then
        Indx_wspeed1=lut_indx-1
        Indx_wspeed2=lut_indx 
        endif
        RETURN
         END
           
       subroutine MOD04_OUT_Extra(MODFIL,Ncell,Nrow,SDS_Value,name) 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C MOD04_OUT : Interface to two subroutines which populate HDF arrays
C              fort land and ocean aerosol retrival algorithm.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT  NONE
      SAVE

       INCLUDE 'mod04.inc'
      INCLUDE 'mapi.inc'
      INCLUDE 'PGS_MODIS_39500.f'
      INTEGER MODFIL(MODFILLEN),Ncell,Nrow,DIM2_BUF(2),START2(2)
      integer j,i,LL,RTN
      CHARACTER*(*) name
C
C COMMON FOR BOTH OCEAN AND LAND
C

        INTEGER*2   
     &     BUF(2*ISWATH*2*ILINE*Tot_scan),
     &     SDS_Value(ISWATH,ILINE*Tot_scan) 

C Write 500 meter cloud mask to HDF file     
 
C Write Sea_Land_mask from Wisconsin Cloud mask.......  
         DIM2_BUF(1)   = Ncell
         DIM2_BUF(2)   = Nrow
          
          START2(1) = 0 
          START2(2) = 0
            LL=0
       do   j = 1, Nrow 
       do   i = 1, Ncell
            LL=LL+1 
             BUF(LL)= SDS_Value(i,j)
          enddo
          enddo 
        RTN = PMAR(MODFIL,name,' ',START2,DIM2_BUF,BUF)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'Problem picking sds ,MOD04_v2.f')  
           Return
           End   
           
        SUBROUTINE cldMsk_Ocean(Refl_3 ,Refl_4,Refl_1,
     &  Refl_26,Refl_5,CldMsk_500_Ocean,High_Cloud_Flag,
     &  data_size,savecldmask)

C----------------------------------------------------------------------
C!F77
C
C!DESCRIPTION:   This subroutine generates newcloud mask for 250*250
C                pixel resolution.
C
C!INPUT PARAMETERS:
C          CLDMSK_500     cldmask for 500 meter resolution
C          num_resol
C
C!OUTPUT PARAMETERS:
C NCLDMSK_SYN new cloud mask at 250 pixel resolution
C
C
C !REVISION HISTORY:
C $Log: Process_ocean_V2.f,v $
c 10/18/1999 fhliang
c fixed prolog;
C
C !TEAM-UNIQUE HEADER:
C
C !REFERENCES AND CREDITS:
C
C !DESIGN NOTES:
C  Editd 05 March 2012 by Leigh Munchak - lowered the lower bound of
C  1.38 reflectance threshold to apply ratio test from .01 to .005
C  
C
C !END
C-----------------------------------------------------------------------

        IMPLICIT  NONE
        SAVE

        INCLUDE   'mod04.inc'

        REAL refl_4(ISWATH*2,ILINE*2),THRSH
        Real refl_3(ISWATH*2,ILINE*2)
        Real refl_4_1km(ISWATH,ILINE)
        Real refl_26(ISWATH,ILINE)
        Real refl_5(ISWATH*2,ILINE*2)
        Integer High_Cloud_Flag(ISWATH,ILINE)
        integer savecldmask(2*ISWATH,2*ILINE)

c rhucek 08/09/02: resized arrays RMED and RMEDSQ to ISWATH*2 
c       REAL RMED(ISWATH),RMEDSQ(ISWATH)
        REAL RMED(ISWATH*2),RMEDSQ(ISWATH*2)

c rhucek 08/09/02: declared NROWS array of size ISWATH*2
        INTEGER  NROWS(ISWATH*2)
        integer CldMsk_500_Ocean(2*ISWATH,2*ILINE+1)
        REAL     refl_1(4*ISWATH,4*ILINE)
        INTEGER NUMSQ,I,J,N,NUM,ISTART,IEND
        INTEGER SQNUM,M,IX,IXX,IY,IYY
        integer CldMsk_1km(ISWATH,ILINE)
        integer  IBLUE,JBLUE,numaver,data_size(2)
        real AREFW550,AREFW670 ,refl_1_500m(2*ISWATH,2*ILINE)
        integer X2_offset,Y2_offset,l,p,LL,PP,add_clear,JX,JY

C     INTIALIZE THE ARRAY FOR NEW CLOUD MASK

            DO IYY=1,ILINE*2
                DO  IXX=1,data_size(1)*2
change here                
                 CldMsk_500_Ocean(IXX,IYY)=-9999
c                 CldMsk_500_Ocean(IXX,IYY)=-9
                  savecldmask(IXX,IYY)=-9999
                ENDDO
             ENDDO
 
   

           THRSH=.0025  
            CALL Ocean_CldMsk_3by3_pixel(ISWATH,ILINE,refl_4,THRSH,NROWS,
     *                             RMED,RMEDSQ,CldMsk_500_Ocean,data_size)
       
c change 670 into 500 meter resolution

         DO   IYY = 1,ILINE*2  
         IBLUE=2*IYY-1
          DO  IXX=1,data_size(1)*2
          JBLUE=2*IXX-1
           numaver=0
           AREFW670=0
          DO IY = IBLUE,2*IYY 
          DO IX= JBLUE,2*IXX 
         if( refl_1( IX,IY) .gt.0 .and. refl_1( IX,IY) .le.1.0) then
             numaver=numaver+1
             AREFW670=AREFW670+ refl_1( IX,IY)  
        endif
        enddo
        enddo 
       
          IF( numaver .eq.4)THEN
          refl_1_500m(ixx,iyy)=AREFW670/numaver
          else 
          refl_1_500m(ixx,iyy)=-9999.99
           ENDIF
           enddo
           enddo 
           
c recover  overwrite cldmask with clear if ratio condition is satsified to recall thick dust
            DO JY=1,ILINE*2
            DO JX=1,data_size(1)*2
            IF( refl_3(JX,JY) .gt.0 .and.refl_1_500m(JX,JY) .gt.0) then
               IF((refl_3(JX,JY)/refl_1_500m(JX,JY)) .le.0.75) THEN
                   CldMsk_500_Ocean(jx,jy)=1
               ENDIF
             ENDIF
          ENDDO
        ENDDO
       
C overwrite cldmask  with cloudy data if reflectance at 0.47 gt than the threhold
          DO JY=1,ILINE*2
          DO JX=1,data_size(1)*2 
           IF( refl_3(JX,JY).gt.0.40) THEN
          CldMsk_500_Ocean(jx,jy)=0  
             ENDIF
         ENDDO
       ENDDO
       
        
               DO IYY=1,ILINE
               IBLUE=2*IYY-1
        DO  IXX=1,data_size(1)
             JBLUE=2*IXX-1  
         DO IY =IBLUE,2*IYY
         DO IX= JBLUE,2*IXX 
       If (refl_26(IXX,IYY) .gt. 0.03)then
            CldMsk_500_Ocean(IX,IY)=0  
      Else
c if reflactance 1.38 between 0.01 and  and 0.03 apply ratio
C Leigh's modification - changed lower bound of threshold from .01 to .005
      If (refl_26(IXX,IYY).gt. 0.005 .and.refl_26(IXX,IYY).le.0.03) then
c  set ratio 
       IF(refl_5(IX,IY).gt.0)then 
       IF((refl_26(IXX,IYY)/refl_5(IX,IY)).GE.0.30)CldMsk_500_Ocean(IX,IY)=0 
      Endif 
       ENDIF 
      Endif 
         CldMsk_500_Ocean(IX,IY)=
     *  (CldMsk_500_Ocean(IX,IY)*High_Cloud_Flag(IXX,IYY))  
c    ENDDO for 500 * 500 meter 20*20 box
          ENDDO
          ENDDO 
C      ENDDO's for 1km pixels for 10*10 box
         ENDDO
         ENDDO
          DO JY=1,ILINE*2
          DO JX=1,data_size(1)*2 
         savecldmask(jx,jy)= CldMsk_500_Ocean(jx,jy)
         Enddo
         Enddo
         return
         end
      
 
C***************************************************************************
        SUBROUTINE Ocean_CldMsk_3by3_pixel( ISWATH,ILINE,refl_4,THRHLD,
     *                                NROWS,RMED,RMEDSQ,CLDMSK,data_size )
C-----------------------------------------------------------------------
C!F77
C
C!DESCRIPTION:
C              This subroutine derive cloud mask using 3 x 3 spatial  
C              variability based upon MODIS 500 m resolution data
C
C!INPUT PARAMETERS:
C
C             Data_Size   Number of pixels in along and against scan
C                ISWATH   Number of pixels at 1 km resolution along scan
C                 ILINE   Number of pixels at 1 km resolution against scan
C                refl_4   Reflectance at 500 meter resolution
C               THRHLD1   Threshold1 of 3x3 spatial variability
C               THRHLD2   Threshold2 of pixel reflectance value
C
C!OUTPUT PARAMETERS:
C
C                CLDMSK  Cloud mask
C                  RMED  Mean
C                RMEDSQ  Standard deviation
C
C!REVISION HISTORY:
C
C!TEAM-UNIQUE HEADER:
C
C Developed by MODIS Aerosol Team at NASA GSFC, Greenbelt, MD
C
C!DESIGN NOTES:
C
C!END
C-----------------------------------------------------------------------
C
      IMPLICIT NONE
      INTEGER IY,JY,JX,NY,ISWATH,ILINE,Data_Size(2),ISTART
      integer IEND
      INTEGER NROWS(ISWATH*2),CLDMSK(ISWATH*2,ILINE*2),add_numcld
      REAL refl_4(ISWATH*2,ILINE*2),THRHLD
      REAL RMED(ISWATH*2),RMEDSQ(ISWATH*2),STD,VAR
      SAVE

C 
C Initialize work arrays 
      DO JX = 1, data_size(1)*2
        NROWS(JX) = 0
        RMED(JX) = 0
        RMEDSQ(JX) = 0
      ENDDO

C
C Checking 500 meter resolution 

      NY=0
       DO IY=2,ILINE*2-1
        ISTART=IY-1
        IEND=ISTART+2

        DO JY=ISTART,IEND
          NY=NY+1 

          DO JX=2,data_size(1)*2-1

c             IF(refl_4(JX,JY).GT.0.0.AND.refl_4(JX+1,JY).GT.0.0.AND.refl_4(JX-1,JY).GT.0.0) THEN
              NROWS(JX)=NROWS(JX) + 1
              RMED(JX)=RMED(JX)+(refl_4(JX,JY)+refl_4(JX+1,JY)+refl_4(JX-1,JY))
              RMEDSQ(JX)=RMEDSQ(JX)+(refl_4(JX,JY)*refl_4(JX,JY)
     1                             + refl_4(JX+1,JY)*refl_4(JX+1,JY)
     2                             + refl_4(JX-1,JY)*refl_4(JX-1,JY))
c            ENDIF
          
C...........make clear determination where possible and re-initialize work array elements
            IF(NY.EQ.3) THEN
              
C..............make clear/cloud determination only when all 9 pixels in 3x3 array are valid 
               IF(NROWS(JX) .EQ. 3) THEN
                  VAR=9.0/8.0*(RMEDSQ(JX)/9.0-RMED(JX)*RMED(JX)/81.0)

                  IF(VAR.GT.0.0) THEN
                     STD=SQRT(VAR)

                     IF(STD.LT.THRHLD) THEN
                        CLDMSK(JX,IY)=1
                        Else
                        CLDMSK(JX,IY)=0
                     ENDIF
                  ENDIF
               ENDIF

               NROWS(JX)=0
               RMED(JX)=0.0
               RMEDSQ(JX)=0.0
             ENDIF

          ENDDO
        ENDDO

        NY=0

      ENDDO 

      DO JX=1,data_size(1)*2
        CLDMSK(JX,1)=CLDMSK(JX,2)
        CLDMSK(JX,ILINE*2)=CLDMSK(JX,ILINE*2-1)
      ENDDO
   
      DO JY=1,ILINE*2
        CLDMSK(1,JY)=CLDMSK(2,JY)
        CLDMSK(data_size(1)*2,JY)=CLDMSK(data_size(1)*2-1,JY)
      ENDDO
      RETURN
      END  
C*********************************************************************

      SUBROUTINE INTERPOlation_R_eight(M,X1,X,Y,Y1,INDEX)

C----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION: This subroutine is a general purpose routine and
C              interpolates linearly.  Value of y1 is interpolated
C              or extrapolated for x1
C
C !INPUT PARAMETERS:I,M,X1,X,Y
C
C !OUTPUT PARAMETERS:Y1,LOPT
C
C !REVISION HISTORY:
C $Log: Process_ocean_V2.f,v $
c
C
C !TEAM-UNIQUE HEADER:
C
C !END
C----------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER IL,LL,LOPT,M,INDEX
      REAL*8  PINTEN,PPHI,SINTEN,SPHI
      REAL*8  X(M),Y(M),Y1,X1,Diff

      Y1=0.0
      LOPT=0
      LL=M-1
        DO 230 IL=1,LL
C        Extrapolation on lower bound
       IF(X1 .LE.X(1))THEN
          PPHI=X(1)
          SPHI=X(2)
          PINTEN=Y(1)
          SINTEN=Y(2)
           Diff=(SPHI-PPHI)
           if(Diff . Le. 0.0) Diff=1
          Y1=PINTEN+((SINTEN-PINTEN)*((X1-PPHI)/Diff))
               LOPT=1

*/  Modified by JC Guu  01/09/97
*/  "GO TO 290" is changed to RETURN

           RETURN

C        Extrapolation on upper bound
       ELSEIF(X1 .GE.X(LL+1)) THEN
           PPHI=X(LL)
           SPHI=X(LL+1)
         PINTEN=Y(LL)
         SINTEN=Y(LL+1)
         Diff=(SPHI-PPHI)
           if(Diff .Le. 0.0) Diff=1
          Y1=PINTEN+((SINTEN-PINTEN)*((X1-PPHI)/Diff))
              LOPT=1

*/  Modified by JC Guu  01/09/97
*/  "GO TO 290" is changed to RETURN

           RETURN

C       interpolation
       ELSEIF (X1 .GE.X(IL) .AND.X1.LE. X(IL+1)) THEN
         PPHI=X(IL)
         SPHI=X(IL+1)
         PINTEN=Y(IL)
         SINTEN=Y(IL+1)
         Diff=(SPHI-PPHI)
           if(Diff .Le. 0.0) Diff=1
          Y1=PINTEN+((SINTEN-PINTEN)*((X1-PPHI)/Diff))
          LOPT=1

*/  Modified by JC Guu  01/09/97
*/  "GO TO 290" is changed to RETURN
*/  and two lines are commented out.

           RETURN
C        ELSE
C        GO TO 230

         ENDIF
  230    CONTINUE

  290    RETURN
          END
         Subroutine compute_Gascorrection(Total_H2o,Total_O3,
     +    SolZen,SatZen,set_counter_for_Gread,Multi_factor,
     +    RTN_NCEP,handle_gascoeff,QA_Flag_Land,QA_Flag_Ocean)
           IMPLICIT  NONE
           SAVE
            INCLUDE 'mod04.inc'
      real Total_H2o,G_factor,SatZen,SolZen,LOGCON,LOGCON2,Total_O3
      real G_factor_flat, G_factor_H2O, G_factor_O3, G_factor_CO2
      real c_VZ, c_SZ
      real hh, r
      integer iwave,wave_num,ifile,mband,vband,RTN_NCEP,Nfile
      integer Waves_used
      parameter(wave_num=10,Waves_used=10) 
      real EXPONENT,H2o_Coef(wave_num,3), Opt_H2O_Clim(wave_num)
      real O3_Coef(wave_num,2),Opt_O3_Clim(wave_num)
      real Opt_CO2_Clim(wave_num), RTrans_H2O(wave_num)
      real RTrans_O3(wave_num),RTrans_CO2(wave_num)
      real Multi_factor(wave_num),drygas
      real wave,mol(wave_num), DEGRAD
      integer set_counter_for_Gread,Viirs_Table,ik,ij
      INTEGER QA_Flag_Land(19),QA_Flag_Ocean(12)
      character * 1 line(132)
         
           
! read the file  one time only
           if( set_counter_for_Gread .eq.1) then  
           
              Nfile= handle_gascoeff   
        
 
!      Rob's change for C6
!           OPEN( Nfile, FILE=  
!     +  '/atmos/mattoo/PGE_codes/PGE04_C6/Tables/MODIS_GAS_COEFS_v1.dat', 
!     +   status='old', form='formatted') 

!    this one is in C6 paper 
!     +   '/atmos/mattoo/PGE_codes/PGE04_C6/Tables/MODIS_GAS_COEFS_v2.dat', 
!     +   status='old', form='formatted') 
!    Falguni's .............
!              OPEN( Nfile, FILE=  
!     +   '/atmos/mattoo/PGE_codes/PGE04_C6/Tables/MODIS_GAS_COEFS_Linebyline_V1.dat', 
!     +   status='old', form='formatted') 
     
1         format(132 a1)        
               read( Nfile,1)line 
               do Iwave=1,Waves_used
        read(Nfile,*) mband,vband,wave,mol(iwave),
     +  Opt_O3_Clim(iwave),Opt_H2O_Clim(iwave),
     +  Opt_CO2_Clim(iwave),O3_Coef(iwave,1),O3_Coef(iwave,2),
     +  H2o_Coef(iwave,1),H2o_Coef(iwave,2), H2o_Coef(iwave,3)   
               enddo 
          Endif
        
           Do Iwave =1,wave_num 
            RTrans_H2O(Iwave) =1.
            RTrans_O3(Iwave)=1.
            RTrans_CO2(Iwave)=1.
           ENDDO  
           
              Total_H2o=Total_H2o/10.
              
              
C
C Calculate gemoetric factor for 2-way transmission
C
      DEGRAD=ACOS(-1.)/180.
        G_factor=-1.0
        G_factor_flat = -1.0
        G_factor_H2O = -1.0
        G_factor_O3 = -1.0
        G_factor_CO2 = -1.0
        

        IF(SatZen.GT.0.0.AND.SolZen.GT.0.0) THEN
           c_VZ = COS(DEGRAD*SatZen)
           c_SZ = COS(DEGRAD*SolZen)

c  Calculate G_factors

c  plane parallel g_factor
            G_factor_flat = (1./c_VZ) + (1./c_SZ)
c  Spherical geometry g_factor
            hh = 9.  ! 9 km atmos scale height
            r = 6371./9.
            G_factor= (SQRT( (r*c_VZ)**2. + 2.*r + 1) - r*c_VZ)
     +    +  (SQRT( (r*c_SZ)**2. + 2.*r + 1) - r*c_SZ)
c  Kasten and Young g_factors (  1. / cosz + a1*z**a2 * (a3-z)**a4 )
            G_factor_H2O = 
     +         (1./(c_VZ + 0.0311*SatZen**(0.1) * (92.471-SatZen)**(-1.3814)))
     +    +    (1./(c_SZ + 0.0311*SolZen**(0.1) * (92.471-SolZen)**(-1.3814)))
            G_factor_O3 = 
     +         (1./(c_VZ + 268.45*SatZen**(0.5) * (115.42-SatZen)**(-3.2922)))
     +    +    (1./(c_SZ + 268.45*SolZen**(0.5) * (115.42-SolZen)**(-3.2922)))
            G_factor_CO2 = 
     +         (1./(c_VZ + 0.4567*SatZen**(0.07) * (96.4836-SatZen)**(-1.6970)))
     +    +    (1./(c_SZ + 0.4567*SolZen**(0.07) * (96.4836-SolZen)**(-1.6970)))

        ENDIF
! keep like old versio........        
!          G_factor_H2O =G_factor_flat
!          G_factor_O3 =G_factor_flat
!          G_factor_CO2=G_factor_flat
         
         
! If NCEP water is available compute Water transmission else use OptH20 from clim.        
      IF(RTN_NCEP .eq. 0 .and. Total_H2O.GT.0.0.AND.G_factor.GT.0.0) THEN 
            QA_Flag_Land(2) =0
            QA_Flag_Ocean(8)=0
            LOGCON=ALOG(Total_H2O*G_factor_H2O)
            LOGCON2=LOGCON*LOGCON 
             Do Iwave =1,Waves_used
               EXPONENT=H2o_Coef(Iwave,1)+H2o_Coef(Iwave,2)*LOGCON
     *              +H2o_Coef(Iwave,3)*LOGCON2
                RTrans_H2O(Iwave)=EXP(EXP(EXPONENT))  
             Enddo
           Else
            QA_Flag_Land(2) =1
            QA_Flag_Ocean(8)=1
             Do Iwave =1,Waves_used
             RTrans_H2O(Iwave)=EXP(Opt_H2O_Clim(iwave)*G_factor_H2O) 
             Enddo
           Endif
           
            
           
           
! If NCEP Ozone is available compute Ozonetransmission else use OptOzone from clim.             
      IF(RTN_NCEP .eq.0 .and.Total_O3.GT.0.0.AND.G_factor.GT.0.0) THEN
              QA_Flag_Land(2) =0
              QA_Flag_Ocean(8)=0
             Do Iwave =1,Waves_used
             EXPONENT=Total_O3*G_factor_O3
      RTrans_O3(Iwave)=EXP(O3_Coef(Iwave,1)+O3_Coef(Iwave,2)*EXPONENT)
            Enddo 
           Else
            QA_Flag_Land(2) =1
            QA_Flag_Ocean(8)=1
               Do Iwave =1,Waves_used
                 RTrans_O3(Iwave)=EXP(Opt_O3_Clim(iwave)*G_factor_O3)
              Enddo
           Endif
! compute rest of gases from cli.           
            Do Iwave =1,Waves_used
           RTrans_CO2(iwave)=EXP(Opt_CO2_Clim(Iwave)*G_factor_CO2)
           Enddo
! compute total transmission       
         
           
            Do Iwave =1,wave_num 
           Multi_factor(Iwave)=RTrans_H2O(Iwave)*RTrans_O3(Iwave)*
     +     RTrans_CO2(Iwave)  
            
           Enddo  
             
          Return 
           end
           

