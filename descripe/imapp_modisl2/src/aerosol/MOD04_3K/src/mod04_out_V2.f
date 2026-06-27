       SUBROUTINE MOD04_OUT(MODFIL,Ncell,Nrow,SDSLAT,SDSLON,
     *SDS_Tau_Land_Ocean_img,SDS_Scattering_Angle,SDS_SST,SDS_MTHET0,SDS_MTHET,SDS_MPHI0,SDS_MPHI,
     *SDS_CldMskQA,SDS_QCONTROL_ocean,SDS_NUMPIXELS,SDS_SCAT_ANGLE_OCEAN, SDS_CLDFRC_ocean ,
     *SDS_small_weighting,SDS_Least_error,SDS_effrad,SDS_sol_INDX_small,SDS_sol_INDX_large,SDS_ccn,
     *SDS_mass_conc,SDS_angs_coeff1,SDS_angs_coeff2,SDS_AOT_model,SDS_ref,SDS_ref_STD,SDSTAU_best,
     *SDSTAU_average,SDSTAUS_best,SDSTAUS_average,SDSTAUB_best,SDSTAUB_average,SDSASSY_best,
     *SDSASSY_average,SDSBACK_best,SDSBACK_average,SDS_RefF_best,SDS_RefF_average,SDS_TranF_best,
     *SDS_TranF_average,SDS_Aerosol_Type,SDS_SCAT_ANGLE_land,SDS_mass_conc_land,SDS_angs_coeff_land,
     *SDS_CLDFRC_land,SDS_dust_weighting,SDS_est_uncer,SDS_RefF_land,SDS_TranF_land,
     *SDS_NUMPIXELS_land,SDSTAU_corrected,SDS_ref_land,SDS_ref_STD_land,SDS_QCONTROL_land,
     *SDS_Mean_Reflectance_Land_All,SDS_SDev_Reflectance_Land_All,SDS_Path_Radiance_Land, 
     *SDS_Critical_Reflectance_Land, SDS_Error_Crit_Reflectance_Land, SDS_Error_Path_Radiance_Land,
     *SDS_QFlag_Critical_Ref_Land,SDS_QFlag_Path_Radiance_Land,SDS_correc_small_weighting,
     *SDS_QCONTROL_CritRef_land,SDS_ratio_small_Land_Ocean,SDSTAU_corrected_213,SDSTAU_small_land, 
     *SDS_Reflected_flux_Land_Ocean,SDS_Surface_Reflectance_Land,SDS_Fitting_Error_Land,SDS_Tau_Land_Ocean,
     *SDS_Quality_Land_Ocean)
C---------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C
C   Create and write Scientific Data Sets (SDS) for mod04 OCEAN AND LAND
C   This sunroutine will be called every scan line and call MOD04S_OUT
C   for ocean and MOD04L_OUT for land.
C
C !INPUT PARAMETERS:
C
C    INTEGER   MODFIL    File handle structure for HDF files
C    INTEGER   Ncell     Number of cells across scan
C    INTEGER   Nrow      Number of (instrument) scans in granule
C    INTEGER   SDS*     All arrays for land and ocean
c    Common Array's for land and ocean  SDS*
C                       Lat, Lon, Scantime, Solarzenith angle, View angle,
C                       Diff azimuth, Azimuth
C Ocean SDS_ARRAYS..........FOR ocean ONLY
C
C SDS_ref_STD     Standard deviation of reflectances at 7 bands
C SDSTAU_best     Optical thickness for best solution
C SDSTAUS_best    Optical thickness contribution small particles for best solution
C SDSTAUB_best    Optical thickness contribution large particles for best solution
C SDSTAU_average  Optical thickness for best solution
C SDSTAUS_average Optical thickness contribution small particles for best solution
C SDSTAUB_average Optical thickness contribution large particles for best solution
C SDS_Least_error         Least square error estimated
C SDS_small_weighting     Small mode weighting factor
C SDS_sol_INDX_small      Solution Number index small particles
C SDS_sol_INDX_large      Solution Number index large particles
C SDSASSY_best      Asymmetry_Factor at 7 bands for best solution
C SDSASSY_average   Asymmetry_Factor at 7 bands for average solution
C SDSBACK_best      Backscattering ratio at 7 bands of best solution
C SDSBACK_average   Backscattering ratio at 7 bands of average solution
C SDS_effrad         Effective_Radius at 0.55 micron of both solutions
C SDS_AOT_model Ratio of optical depth of small mode vs effective optical depth at 550
C SDS_RefF_best     Normalized reflected_flux at 7 bands of best solution
C SDS_RefF_average  Normalized reflected_flux at 7 bands of average solution
C SDS_TranF_best    Normalized Transmitted_flux at 7 bands of best solution
C SDS_TranF_average Normalized Transmitted_flux at 7 bands of average solution
C SDS_SCAT_ANGLE_OCEAN Scattering angle ocean
C SDS_QCONTROL         Quality control SDS array
C SDS_NUMPIXELS        Number of Pixels used for 0.55 micron
C SDS_ccn              Cloud_Fraction in percentage
C SDS_mass_conc        Mass concentration
C SDS_angs_coeff1      Angstrom Exponent for 0.550 and 0.865 miron
C SDS_angs_coeff2      Angstrom exponent for 0.865 and 2.130 micron

C
C LAND SDS_ARRAYS..........FOR LAND ONLY
C
c SDS_QCONTROL_land       Quality control SDS array
c SDS_Aerosol_Type        Index of Aerosol type
c SDS_SCAT_ANGLE_land     Scattering Angle
c SDS_angs_coeff_land     Angstrom exponent for 0.47 and 0.67 micron
c SDS_CLDFRC_land         Cloud fraction (%)
c SDS_dust_weighting      Dust aerosol weighting factor 
c SDS_est_uncer           Uncertainty of optical thickness at 0.47 and 0.66 micron
c SDS_RefF_land          Normalized reflected flux at 0.47 and 0.66 microns
c SDS_TranF_land         Normalized Transmitted flux at 0.47 and 0.66 microns
c SDS_NUMPIXELS_land     Number of pixels with desired percentile
c SDSTAU_corrected       Corrected optical thickness at 0.47 0.55 and 0.66 micron
c SDS_ref_land           Mean reflectance at five bands
c SDS_ref_STD_land       Standard deviation of reflectance at five bands
c SDS_mass_conc_land     Mass concentration
C
C !OUTPUT PARAMETERS:    NONE
C
C !REVISION HISTORY:
c 01/28/98 fhliang
c fixed prolog.
c
c Revision 1.1  1996/06/04  16:01:42  vlin
c Initial revision
c
c!TEAM-UNIQUE HEADER:
c
c Developed by MODIS Aerosol/Solar Water Vapor Retrieval Team
c GSFC, Greenbelt, MD
c
C !REFERENCES and CREDITS:
C
C    Written by:
C       Shana Mattoo
C       SAC/ARD
C
C       mattoo@climate.gsfc.nasa.gov
C
C !DESIGN NOTES:
C
C    Externals:
C       Functions:
C
C    Internal Subroutine:
C        MOD04S_OUT   writes Ocean HDF arrays
C        MOD04L_OUT   writes Land  HDF arrays
C
C !END
C----------------------------------------------------------------------
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C MOD04_OUT : Interface to two subroutines which populate HDF arrays
C              fort land and ocean aerosol retrival algorithm.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      IMPLICIT  NONE
      SAVE

c rhucek 01/08/98:  replaced include file mapic.inc with mapi.inc
c     Include 'mapic.inc'
      Include 'mapi.inc'
      Include 'mod04.inc'

      INTEGER MODFIL(MODFILLEN),Ncell,Nrow,IFLAG_land,IFLAG_ocean

C
C COMMON FOR BOTH OCEAN AND LAND
C

        INTEGER*2  SDS_MTHET0(NUMCELLS),SDS_MTHET(NUMCELLS),
     &             SDS_MPHI0(NUMCELLS),SDS_MPHI(NUMCELLS),
     &             SDS_Scattering_Angle(NUMCELLS)
        REAL       SDSLAT(NUMCELLS),SDSLON(NUMCELLS)
        REAL*8     SDS_SST(NUMCELLS)
        BYTE       SDS_CldMskQA(NUMCELLS)

C
C         LAND SDS_ARRAYS........
C

        BYTE       SDS_QCONTROL_land(QA_LAND,NUMCELLS),
     &             SDS_QCONTROL_CritRef_land(QA_LAND,NUMCELLS)
        INTEGER*2  SDS_Tau_Land_Ocean_img(NUMCELLS),
     &             SDS_Tau_Land_Ocean(NUMCELLS),   
     &             SDS_Aerosol_Type(NUMCELLS),
     &             SDS_SCAT_ANGLE_land(NUMCELLS),
     &             SDS_angs_coeff_land(NUMCELLS),
     &             SDS_CLDFRC_land(NUMCELLS),
     &             SDS_dust_weighting(NUMCELLS),
     &             SDS_NUMPIXELS_land(NUMCELLS,Land_Sol1),
     &             SDSTAU_corrected(NUMCELLS,Land_Sol3),
     &             SDSTAU_corrected_213(NUMCELLS),
     &             SDS_ref_land(NUMCELLS,Band_land),
     &             SDS_ref_STD_land(NUMCELLS,Band_land),
     &             SDS_ratio_small_Land_Ocean(NUMCELLS),
     &            SDS_Reflected_flux_Land_Ocean(NUMCELLS),
c  9/2005 ( two new SDS's)......     
     &            SDS_Surface_Reflectance_Land(NUMCELLS,Land_Sol3),
     &            SDS_Fitting_Error_Land(NUMCELLS),
     &             SDS_Quality_Land_Ocean(NUMCELLS)
       INTEGER*2  SDSTAU_small_land(NUMCELLS,Land_Sol4)
     
        REAL       SDS_mass_conc_land(NUMCELLS)
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
     &  SDS_QFlag_Path_Radiance_Land(NUMCELLS,Land_Sol1)

C
C Obsolete (02/2006) Land SDS Arrays
C
      INTEGER *2   
     &            SDS_est_uncer(NUMCELLS,Land_Sol1),
     &            SDS_RefF_land(NUMCELLS,Land_Sol2),
     &            SDS_TranF_land(NUMCELLS,Land_Sol1)


C
C         OCEAN SDS_ARRAYS........
C

        BYTE      SDS_QCONTROL_ocean(QA_OCEAN,NUMCELLS)
        REAL      SDS_ccn(NUMCELLS,NUM_solutions),
     &            SDS_mass_conc(NUMCELLS,NUM_solutions)
       INTEGER *2 SDS_ref(NUMCELLS,NWAV_S),
     &            SDS_ref_STD(NUMCELLS,NWAV_S),
     &            SDSTAU_best(NUMCELLS,NWAV_S),
     &            SDSTAU_average(NUMCELLS,NWAV_S),
     &            SDSTAUB_best(NUMCELLS,NWAV_S),
     &            SDSTAUB_average(NUMCELLS,NWAV_S),
     &            SDSTAUS_best(NUMCELLS,NWAV_S),
     &            SDSTAUS_average(NUMCELLS,NWAV_S),
     &            SDSASSY_best(NUMCELLS,NWAV_S),
     &            SDSASSY_average(NUMCELLS,NWAV_S),
     &            SDSBACK_best(NUMCELLS,NWAV_S),
     &            SDSBACK_average(NUMCELLS,NWAV_S),
     &            SDS_RefF_best(NUMCELLS,NWAV_S),
     &            SDS_RefF_average(NUMCELLS,NWAV_S),
     &            SDS_TranF_best(NUMCELLS,NWAV_S),
     &            SDS_TranF_average(NUMCELLS,NWAV_S)
       INTEGER*2  SDS_small_weighting(NUMCELLS,NUM_solutions),
     &       SDS_correc_small_weighting(NUMCELLS),
     &            SDS_Least_error(NUMCELLS,NUM_solutions),
     &            SDS_effrad(NUMCELLS,NUM_solutions),
     &            SDS_sol_INDX_small(NUMCELLS,NUM_solutions),
     &            SDS_sol_INDX_large(NUMCELLS,NUM_solutions),
     &            SDS_angs_coeff1(NUMCELLS,NUM_solutions),
     &            SDS_angs_coeff2(NUMCELLS,NUM_solutions),
     &            SDS_AOT_model(NUMCELLS,num_model_index),
     &            SDS_NUMPIXELS(NUMCELLS), SDS_CLDFRC_ocean(NUMCELLS),
     &            SDS_SCAT_ANGLE_OCEAN(NUMCELLS)
C
C Buffer arrays for latitude, longitude, solar zenith, solar azimuth
C                   satellite zenith, and satellite azimuth
C
      INTEGER  N_Entries,BUF_Line_No
      REAL     BUF_LAT(NUMCELLS*Lines_Per_Write),
     &         BUF_LON(NUMCELLS*Lines_Per_Write)
      INTEGER*2
     &         BUF_MTHET0(NUMCELLS*Lines_Per_Write),
     &         BUF_MTHET(NUMCELLS*Lines_Per_Write),
     &         BUF_MPHI0(NUMCELLS*Lines_Per_Write),
     &         BUF_MPHI(NUMCELLS*Lines_Per_Write)
      DATA      N_Entries /0/

C
C Populate ocean SDS arrays
C

       CALL MOD04S_OUT(MODFIL,Ncell,Nrow,SDSLAT,SDSLON,
     *SDS_SST,SDS_MTHET0,SDS_MTHET,SDS_MPHI0,SDS_MPHI,
     *SDS_Scattering_Angle,SDS_QCONTROL_ocean,SDS_Tau_Land_Ocean_img,
     *SDS_NUMPIXELS,SDS_SCAT_ANGLE_OCEAN, SDS_CLDFRC_ocean ,
     *SDS_small_weighting,SDS_Least_error,SDS_effrad,
     *SDS_sol_INDX_small,SDS_sol_INDX_large,
     *SDS_ccn,SDS_mass_conc,SDS_angs_coeff1,SDS_angs_coeff2,
     *SDS_AOT_model,SDS_ref,SDS_ref_STD,SDSTAU_best,
     *SDSTAU_average,SDSTAUS_best,SDSTAUS_average,SDSTAUB_best,
     *SDSTAUB_average,SDSASSY_best,SDSASSY_average,SDSBACK_best,
     *SDSBACK_average,SDS_RefF_best,SDS_RefF_average,SDS_TranF_best,
     *SDS_TranF_average,SDS_correc_small_weighting,
     *SDS_ratio_small_Land_Ocean,SDS_Reflected_flux_Land_Ocean,
     *SDS_Tau_Land_Ocean,SDS_Quality_Land_Ocean)

C
C Populate land SDS arrays
C

      CALL MOD04L_OUT(MODFIL,Ncell,Nrow,SDSLAT,SDSLON,
     *SDS_SST,SDS_MTHET0,SDS_MTHET,SDS_MPHI0,SDS_MPHI,
     *SDS_Scattering_Angle,
     *SDS_Tau_Land_Ocean_img,SDS_CldMskQA,SDS_Aerosol_Type,
     *SDS_SCAT_ANGLE_land,SDS_mass_conc_land,SDS_angs_coeff_land,
     *SDS_CLDFRC_land,SDS_dust_weighting,
     *SDS_est_uncer,SDS_RefF_land,SDS_TranF_land,SDS_NUMPIXELS_land,
     *SDSTAU_corrected,SDS_ref_land,SDS_ref_STD_land,SDS_QCONTROL_land,
     *SDS_Mean_Reflectance_Land_All,SDS_Tau_Land_Ocean, 
     *SDS_SDev_Reflectance_Land_All,SDSTAU_small_land,     
     *SDS_Path_Radiance_Land,   
     *SDS_Critical_Reflectance_Land,  
     *SDS_Error_Crit_Reflectance_Land,     
     *SDS_Error_Path_Radiance_Land,  
     *SDS_QFlag_Critical_Ref_Land,SDS_QFlag_Path_Radiance_Land,
     *SDS_QCONTROL_CritRef_land,SDS_ratio_small_Land_Ocean,
     *SDS_Reflected_flux_Land_Ocean,SDSTAU_corrected_213,
     *SDS_Surface_Reflectance_Land,SDS_Fitting_Error_Land,
     * SDS_Quality_Land_Ocean)

       RETURN
       END


C*********************************************************************
       SUBROUTINE MOD04S_OUT(MODFIL,Ncell,Nrow,SDSLAT,SDSLON,
     *SDS_SST,SDS_MTHET0,SDS_MTHET,SDS_MPHI0,SDS_MPHI,
     *SDS_Scattering_Angle,SDS_QCONTROL_ocean,SDS_Tau_Land_Ocean_img,
     *SDS_NUMPIXELS,SDS_SCAT_ANGLE_OCEAN, SDS_CLDFRC_ocean ,
     *SDS_small_weighting,SDS_Least_error,SDS_effrad,
     *SDS_sol_INDX_small,SDS_sol_INDX_large,
     *SDS_ccn,SDS_mass_conc,SDS_angs_coeff1,SDS_angs_coeff2,
     *SDS_AOT_model,SDS_ref,SDS_ref_STD,SDSTAU_best,
     *SDSTAU_average,SDSTAUS_best,SDSTAUS_average,SDSTAUB_best,
     *SDSTAUB_average,SDSASSY_best,SDSASSY_average,SDSBACK_best,
     *SDSBACK_average,SDS_RefF_best,SDS_RefF_average,SDS_TranF_best,
     *SDS_TranF_average,SDS_correc_small_weighting,
     *SDS_ratio_small_Land_Ocean,SDS_Reflected_flux_Land_Ocean,
     *SDS_Tau_Land_Ocean,SDS_Quality_Land_Ocean)

C---------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C
C   Create and write Scientific Data Sets (SDS) for mod04 sea.
C   This sunroutine will be called every scan line but write data into
C   HDF file at every Lines_Per_Write.
C
C !INPUT PARAMETERS:
C
C    INTEGER   MODFIL    File handle structure for HDF files
C    INTEGER   Ncell     Number of cells across scan
C    INTEGER   Nrow      Number of (instrument) scans in granule
c
C    All arrays for ocean    SDS*
C                       Lat, Lon, Scantime, Solarzenith angle, View angle,
C                       Diff azimuth, Azimuth
C
C Ocean SDS_ARRAYS..........FOR ocean ONLY
C
C SDS_ref_STD          Standard deviation of reflectances at 7 bands
C SDSTAU_best          Optical thickness for best solution
C SDSTAUS_best         Optical thickness contribution small particles for
C                      best solution
C SDSTAUB_best         Optical thickness contribution large particles for
C                      best solution
C SDSTAU_average       Optical thickness for average solution
C SDSTAUS_average      Optical thickness contribution small particles for
C                      average solution
C SDSTAUB_average      Optical thickness contribution large particles for
C                      average solution
C SDS_Least_error      Least square error estimated
C SDS_small_weighting  Small mode weighting factor
C SDS_sol_INDX_small   Solution Number index small particles
C SDS_sol_INDX_large   Solution Number index large particles
C SDSASSY_best         Asymmetry_Factor at 7 bands for best solution
C SDSASSY_average      Asymmetry_Factor at 7 bands for average solution
C SDSBACK_best         Backscattering ratio at 7 bands of best solution
C SDSBACK_average      Backscattering ratio at 7 bands of average solution
C SDS_effrad           Effective_Radius at 0.55 micron of both solutions
C SDS_AOT_model   Ratio of optical depth of small mode vs effective optical
C                      depth at 550
C SDS_RefF_best        Normalized reflected_flux at 7 bands of best solution
C SDS_RefF_average     Normalized reflected_flux at 7 bands of average solution
C SDS_TranF_best       Normalized Transmitted_flux at 7 bands of best solution
C SDS_TranF_average    Normalized Transmitted_flux at 7 bands of average solution
C SDS_SCAT_ANGLE_OCEAN Scattering angle ocean
C SDS_QCONTROL         Quality control SDS array
C SDS_NUMPIXELS        Number of Pixels used for 0.55 micron
C SDS_ccn              Cloud_Fraction in percentage
C SDS_mass_conc        Mass concentration
C SDS_angs_coeff1      Angstrom Exponent for 0.550 and 0.865 miron
C SDS_angs_coeff2      Angstrom exponent for 0.865 and 2.130 micron
C
C !OUTPUT PARAMETERS:    NONE
C
C !REVISION HISTORY:
c Revision 1.1  1996/06/04  16:01:42  vlin
c Initial revision
c
C !TEAM-UNIQUE HEADER:
C
C    This software is developed by the MODIS Science Data Support
C    Team for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C !REFERENCES and CREDITS:
C
C    Written by:
C       Vicky Lin                    June 1996
C       Research and Data systems Corporation
C       SAIC/GSC MODIS Science Data Support Office
C       7501 Forbes Blvd, Seabrook MD 20706
C       vlin@ltpmail.gsfc.nasa.gov
C
C
C      Modified by shana Mattoo       March 20001
C
C
C !DESIGN NOTES:
C
C    This subroutine checks the return status of all MODIS Application
C    Program Interface (M-API) function calls.  A successful M-API
C    function call is indicated by a return status (RTN) equal to
C    MAPIOK.  If a M-API function call is not successful, a warning
C    message is written to the LogStatus file by passing an mnemonic
C    (.._W_..) to the function MODIS_SMF_SETDYNAMICMSG.
C
C    Externals:
C       Functions:
C          CRMAR   Creates an HDF SDS structure to store a new data
C                  array into a MODIS HDF file. (libmapi.a)
C          PMAR    PMAR places a multi-dimensional array of data
C                  into an HDF SDS array structure previously
C                  created using CRMAR.         (libmapi.a)
C
C       Named Constant:
C          MODIS_W_GENERIC       (PGS_MODIS_39500.f)
C          SDS*_N                (mod04.inc)
C          SCALE*                (mod04.inc)
C          OFFSET*               (mod04.inc)
C
C    Internal Subroutine:
C          MODIS_SetDimName  Write the dimension names of a SDS array.
C
C !END
C----------------------------------------------------------------------

      IMPLICIT NONE
      SAVE

      INCLUDE 'mod04.inc'
      INCLUDE 'mapi.inc'
      INCLUDE 'PGS_MODIS_39500.f'


      INTEGER I,J,Ncell,Nrow,RANK1,RANK3,RANK2,LL,
     &        N_Entries,RTN,BUF_Line_No,Nline,
     &        MODFIL(MODFILLEN),START3(3),START3_2(3),START2(2),
     &        DIM3_BUF(3),DIM3_BUF_2(3),DIM2_BUF(2),
     &        DIM_SDS(3),values(2)

C
C Declare all incoming arrays
C

      REAL*8      SDS_SST(NUMCELLS)
      BYTE        SDS_QCONTROL_ocean(QA_OCEAN,NUMCELLS)
      REAL        SDSLAT(NUMCELLS),SDSLON(NUMCELLS),
     &            SDS_ccn(NUMCELLS,NUM_solutions),
     &            SDS_mass_conc(NUMCELLS,NUM_solutions)
      INTEGER *2  SDS_ref(NUMCELLS,NWAV_S),
     &            SDS_ref_STD(NUMCELLS,NWAV_S),
     &            SDSTAU_best(NUMCELLS,NWAV_S),
     &            SDSTAU_average(NUMCELLS,NWAV_S),
     &            SDSTAUS_best(NUMCELLS,NWAV_S),
     &            SDSTAUS_average(NUMCELLS,NWAV_S),
     &            SDSTAUB_best(NUMCELLS,NWAV_S),
     &            SDSTAUB_average(NUMCELLS,NWAV_S),
     &            SDSASSY_best(NUMCELLS,NWAV_S),
     &            SDSASSY_average(NUMCELLS,NWAV_S),
     &            SDSBACK_best(NUMCELLS,NWAV_S),
     &            SDSBACK_average(NUMCELLS,NWAV_S),
     &            SDS_RefF_best(NUMCELLS,NWAV_S),
     &            SDS_RefF_average(NUMCELLS,NWAV_S),
     &            SDS_TranF_best(NUMCELLS,NWAV_S),
     &            SDS_TranF_average(NUMCELLS,NWAV_S)
c
      INTEGER *2  SDS_small_weighting(NUMCELLS,NUM_solutions),
     &            SDS_correc_small_weighting(NUMCELLS),
     &            SDS_Least_error(NUMCELLS,NUM_solutions),
     &            SDS_effrad(NUMCELLS,NUM_solutions),
     &            SDS_sol_INDX_small(NUMCELLS,NUM_solutions),
     &            SDS_sol_INDX_large(NUMCELLS,NUM_solutions),
     &            SDS_angs_coeff1(NUMCELLS,NUM_solutions),
     &            SDS_angs_coeff2(NUMCELLS,NUM_solutions),
     &            SDS_AOT_model(NUMCELLS,num_model_index),
     &            SDS_NUMPIXELS(NUMCELLS),SDS_MTHET0(NUMCELLS),
     &            SDS_MTHET(NUMCELLS),SDS_MPHI0(NUMCELLS),
     &            SDS_MPHI(NUMCELLS),SDS_Scattering_Angle(NUMCELLS),
     &            SDS_SCAT_ANGLE_OCEAN(NUMCELLS),
     &            SDS_CLDFRC_ocean (NUMCELLS),
     &            SDS_Tau_Land_Ocean_img(NUMCELLS),
     &            SDS_Tau_Land_Ocean(NUMCELLS),
     &            SDS_ratio_small_Land_Ocean(NUMCELLS),
     &            SDS_Reflected_flux_Land_Ocean(NUMCELLS),
     &            SDS_Quality_Land_Ocean(NUMCELLS)
c
c Set up BUFfer for the arrays to write

       BYTE    BUF_QCONTROL(QA_OCEAN*NUMCELLS*Lines_Per_Write)
       REAL    BUF_LAT(NUMCELLS*Lines_Per_Write),
     &         BUF_LON(NUMCELLS*Lines_Per_Write),
     &         BUF_ccn(NUMCELLS*Lines_Per_Write*NUM_solutions),
     &         BUF_mass_conc(NUMCELLS*Lines_Per_Write*NUM_solutions)
       REAL*8  BUF_SST(NUMCELLS*Lines_Per_Write)

      INTEGER*2 BUF_REF(NUMCELLS*Lines_Per_Write*NWAV_S),
     &         BUF_REF_STD(NUMCELLS*Lines_Per_Write*NWAV_S),
     &         BUF_TAU_best(NUMCELLS*Lines_Per_Write*NWAV_S),
     &         BUF_TAU_average(NUMCELLS*Lines_Per_Write*NWAV_S),
     &         BUF_TAUS_best(NUMCELLS*Lines_Per_Write*NWAV_S),
     &         BUF_TAUS_average(NUMCELLS*Lines_Per_Write*NWAV_S),
     &         BUF_TAUB_best(NUMCELLS*Lines_Per_Write*NWAV_S),
     &         BUF_TAUB_average(NUMCELLS*Lines_Per_Write*NWAV_S),
     &         BUF_ASSY_best(NUMCELLS*Lines_Per_Write*NWAV_S),
     &         BUF_ASSY_average(NUMCELLS*Lines_Per_Write*NWAV_S),
     &         BUF_BACK_best(NUMCELLS*Lines_Per_Write*NWAV_S),
     &         BUF_BACK_average(NUMCELLS*Lines_Per_Write*NWAV_S),
     &         BUF_RefF_best(NUMCELLS*Lines_Per_Write*NWAV_S),
     &         BUF_RefF_average(NUMCELLS*Lines_Per_Write*NWAV_S) 
   
       INTEGER*2
     &    BUF_small_weighting(NUMCELLS*Lines_Per_Write*NUM_solutions),
     &    BUF_correc_small_weighting(NUMCELLS*Lines_Per_Write),
     &    BUF_Least_error(NUMCELLS*Lines_Per_Write*NUM_solutions),
     &    BUF_effrad(NUMCELLS*Lines_Per_Write*NUM_solutions),
     &    BUF_sol_INDX_small(NUMCELLS*Lines_Per_Write*NUM_solutions),
     &    BUF_sol_INDX_large(NUMCELLS*Lines_Per_Write*NUM_solutions),
     &    BUF_angs_coeff1(NUMCELLS*Lines_Per_Write*NUM_solutions),
     &    BUF_angs_coeff2(NUMCELLS*Lines_Per_Write*NUM_solutions),
     &    BUF_AOT_model(NUMCELLS*Lines_Per_Write*num_model_index),
     &    BUF_NUMPIXELS(NUMCELLS*Lines_Per_Write)
       INTEGER*2
     &    BUF_MTHET0(NUMCELLS*Lines_Per_Write),
     &    BUF_MTHET(NUMCELLS*Lines_Per_Write),
     &    BUF_MPHI0(NUMCELLS*Lines_Per_Write),
     &    BUF_MPHI(NUMCELLS*Lines_Per_Write),
     &    BUF_SCATT(NUMCELLS*Lines_Per_Write),
     &    BUF_SCATTANG(NUMCELLS*Lines_Per_Write),
     &    BUF_CLDFRC(NUMCELLS*Lines_Per_Write),
     &    BUF_550(NUMCELLS*Lines_Per_Write),
     &    BUF_550_LS(NUMCELLS*Lines_Per_Write),
     &    BUF_small_wgt(NUMCELLS*Lines_Per_Write),
     &    BUF_refflux(NUMCELLS*Lines_Per_Write),
     &    BUF_Quality(NUMCELLS*Lines_Per_Write)

      LOGICAL   Write_Flag
      DATA      N_Entries /0/

C
C Initialization
C
      DIM3_BUF(1)   = Ncell
      DIM3_BUF(2)   = Lines_Per_Write
      DIM3_BUF_2(1) = QA_Ocean
      DIM3_BUF_2(2) = Ncell
      DIM2_BUF(1)   = Ncell
      DIM2_BUF(2)   = Lines_Per_Write
      N_Entries = N_Entries + 1
      Write_Flag = .false.
      START2(1) = 0
    
      DO 10 I = 1, 3
         START3(I)   = 0
         START3_2(I) =0
   10 CONTINUE

C s into dimensions 1 and 3,
C respectively, of 3-D write BUFfers.
C Dimension 2 is BUFfer (scan) line counter.

      BUF_Line_No = mod(N_Entries,Lines_Per_Write)
      if (BUF_Line_No.eq.0) then
          BUF_Line_No = Lines_Per_Write
          Write_Flag = .true.
      end if

       do 15 j=1,Ncell
       do 20 i=1,QA_Ocean
         LL=QA_Ocean*Ncell*(Buf_Line_No-1)+QA_Ocean*(j-1)+i
         Buf_QCONTROL(LL)=SDS_QCONTROL_ocean(i,j)
   20  continue
   15  continue
       
       do 25 i = 1, Ncell
         LL = Ncell*(BUF_Line_No - 1) + i
         
         BUF_LAT(LL)           = SDSLAT(i)
         BUF_LON(LL)           = SDSLON(i)
         BUF_SST(LL)           = SDS_SST(i)
         BUF_MTHET0(LL)        = SDS_MTHET0(i)
         BUF_MTHET(LL)         = SDS_MTHET(i)
         BUF_MPHI0(LL)         = SDS_MPHI0(i)
         BUF_MPHI(LL)          = SDS_MPHI(i)
         BUF_NUMPIXELS(LL)     = SDS_NUMPIXELS(i)
         BUF_SCATT(LL)         = SDS_SCAT_ANGLE_OCEAN(i)
         BUF_SCATTANG(LL)      = SDS_Scattering_Angle(i)
         BUF_CLDFRC(LL)        = SDS_CLDFRC_ocean (i)
         BUF_550(LL)           = SDS_Tau_Land_Ocean_img(i)
         BUF_550_LS(LL)        = SDS_Tau_Land_Ocean(i)
         BUF_small_wgt(LL)=  SDS_ratio_small_Land_Ocean(i) 
         BUF_Quality(LL)=SDS_Quality_Land_Ocean(i)
   25  continue

      do 30 j = 1,NUM_solutions
      do 35 i = 1, Ncell
         LL = Ncell*Lines_Per_Write*(j-1)
     &      + Ncell*(BUF_Line_No - 1) + i
          BUF_small_weighting(LL) = SDS_small_weighting(i,j)
         BUF_Least_error(LL)     = SDS_Least_error(i,j)
         BUF_effrad(LL)          = SDS_effrad(i,j)
         BUF_sol_INDX_small(LL)        = SDS_sol_INDX_small(i,j)
         BUF_sol_INDX_large(LL)        = SDS_sol_INDX_large(i,j)
         BUF_ccn(LL)             = SDS_ccn(i,j)
         BUF_mass_conc(LL)       = SDS_mass_conc(i,j)
         BUF_angs_coeff1(LL)     = SDS_angs_coeff1(i,j)
         BUF_angs_coeff2(LL)     = SDS_angs_coeff2(i,j)
   35  continue
   30  continue
       do 130 j = 1,num_model_index
       do 135 i = 1, Ncell
         LL = Ncell*Lines_Per_Write*(j-1)
     &      + Ncell*(BUF_Line_No - 1) + i
         BUF_AOT_model(LL)     = SDS_AOT_model(i,j)
 135     continue
 130     continue

      do 40 j = 1,NWAV_S
      do 45 i = 1, Ncell
         LL = Ncell*Lines_Per_Write*(j-1)
     &      + Ncell*(BUF_Line_No - 1) + i
c
        BUF_REF(LL)              = SDS_ref(i,j)
        BUF_REF_STD(LL)          = SDS_ref_STD(i,j)
        BUF_TAU_best(LL)         = SDSTAU_best(i,j)
        BUF_TAU_average(LL)      = SDSTAU_average (i,j)
        BUF_TAUS_best(LL)        = SDSTAUS_best(i,j)
        BUF_TAUS_average(LL)     = SDSTAUS_average(i,j)
        BUF_TAUB_best(LL)        = SDSTAUB_best(i,j)
        BUF_TAUB_average(LL)     = SDSTAUB_average(i,j)
        BUF_ASSY_best(LL)        = SDSASSY_best(i,j)
        BUF_ASSY_average(LL)     = SDSASSY_average(i,j)
        BUF_BACK_best(LL)        = SDSBACK_best(i,j)
        BUF_BACK_average(LL)     = SDSBACK_average(i,j)
        BUF_RefF_best(LL)        = SDS_RefF_best(i,j)
        BUF_RefF_average(LL)     = SDS_RefF_average(i,j)
 
   45 continue
   40 continue

CCCCCCCCCCCCCCC   If BUF_Line_No equals to Lines_Per_Write   CCCCCCCCCCCCCC

       IF (Write_Flag) THEN

         Nline = N_Entries/Lines_Per_Write
         START2(2) = (Nline-1) * Lines_Per_Write
         START3(2) = (Nline-1) * Lines_Per_Write
         START3_2(3)= (Nline-1) * Lines_Per_Write
         write(36,*) 'Ocean',Nline,START2(1),START2(2),
     &       START3(1),START3(2), START3_2(3),N_Entries
         
c Dimension (cells along swath,cells across swath)
c

      RTN = PMAR(MODFIL,SDS1_N,' ',START2,DIM2_BUF,BUF_LAT)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for 1st array','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,SDS2_N,' ',START2,DIM2_BUF,BUF_LON)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for 2nd array','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,SDS1_NE1,' ',START2,DIM2_BUF,BUF_SST)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Scan_Start_Time','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,SDS1_NE2,' ',START2,DIM2_BUF,BUF_MTHET0)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Solar_Zenith','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,SDS1_NE4,' ',START2,DIM2_BUF,BUF_MTHET)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Sensor_Zenith','MOD04_SEA_OUT')
C
      RTN = PMAR(MODFIL,SDS1_NE3,' ',START2,DIM2_BUF,BUF_MPHI0)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Solar_Azimuth','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,SDS1_NE5,' ',START2,DIM2_BUF,BUF_MPHI)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Sensor_Azimuth','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_numpix,' ',START2,DIM2_BUF,BUF_NUMPIXELS)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for NUMBER OF PIXELS USED',
     &    'MOD04_SEA_OUT')
c
C      RTN = PMAR(MODFIL,Fid_scatt_angle,' ',START2,DIM2_BUF,BUF_SCATT)
C         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
C     &   (MODIS_W_GENERIC,'PMAR for Scattering Angle',
C     &    'MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_scattering_angle,' ',START2,DIM2_BUF
     &           ,BUF_SCATTANG)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for Scattering Angle',
     &    'MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_cldfrac,' ',START2,DIM2_BUF, BUF_CLDFRC)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for Cld fraction',
     &    'MOD04_SEA_OUT')
C
       RTN = PMAR(MODFIL,Fid_tau_land_and_ocean_img,' ',START2,DIM2_BUF,
     & BUF_550)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for aerosol tau image at 0.55 micron',
     &   'MOD04_SEA_OUT')
     
C 2/2006
       RTN = PMAR(MODFIL,Fid_tau_land_and_ocean,' ',START2,DIM2_BUF,
     & BUF_550_LS)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for aerosol tau at 0.55 micron',
     &   'MOD04_SEA_OUT')
C1/2002
      RTN = PMAR(MODFIL,Fid_Ratio_Small_Land_And_Ocean,' ', 
     & START2,DIM2_BUF,BUF_small_wgt)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for small mode ratio Ocean and land',
     &   'MOD04_SEA_OUT')
 

C 3/2010
       RTN = PMAR(MODFIL,Fid_land_sea_Quality,' ',START2,DIM2_BUF,
     & BUF_Quality)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for land_sea_Quality',
     &   'MOD04_SEA_OUT')

c
c  Dimension (QA_ocean, cells along swath, cells across swath)
c
      DIM3_BUF_2(3) = Lines_Per_Write

      RTN = PMAR(MODFIL,Fid_qa,' ',START3_2,DIM3_BUF_2,BUF_QCONTROL)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for QUALITY CONTROL','MOD04_SEA_OUT')

c
c  Dimension (cells along swath,cells across swath,solution number=2)
c

      DIM3_BUF(3) = NUM_solutions
c
      RTN = PMAR(MODFIL,Fid_smallmode,' ',START3,DIM3_BUF,
     &               BUF_small_weighting)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &  (MODIS_W_GENERIC,'PMAR fo SMALL MODE WEIGHT','MOD04_SEA_OUT')

      RTN = PMAR(MODFIL,Fid_leasterror,' ',START3,DIM3_BUF,
     &               BUF_Least_error)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Least sqare ','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_eff_rad,' ',START3,DIM3_BUF,Buf_effrad)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for EFF RAD','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_sol_small,' ',START3,DIM3_BUF,
     &           BUF_sol_INDX_small)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Indx Solution','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_sol_large,' ',START3,DIM3_BUF,
     &           BUF_sol_INDX_large)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for Indx Solution','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_ccn,' ',START3,DIM3_BUF,BUF_ccn)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Cloud_Condensation_Nuclei',
     &   'MOD04_SEA_OUT')
c
c
         RTN = PMAR(MODFIL,Fid_mass_con,' ',START3,DIM3_BUF,
     &             BUF_mass_conc)
        if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Mass_Concentration_Ocean',
     &   'MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_angs_coff1,' ',START3,DIM3_BUF,
     &             BUF_angs_coeff1)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Angstrom_Exponent_1_Ocean',
     &   'MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_angs_coff2,' ',START3,DIM3_BUF,
     &             BUF_angs_coeff2)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Angstrom_Exponent_2_Ocean',
     &   'MOD04_SEA_OUT')
c
        DIM3_BUF(3) = num_model_index
      RTN = PMAR(MODFIL,Fid_AOT_model,' ',START3,DIM3_BUF,
     &             BUF_AOT_model)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for Optical_Depth_Ratio_Small_Ocean',
     &   'MOD04_SEA_OUT')

c
c Dimension (cells along swath,cells across swath,number of wavelength)
c

         DIM3_BUF(3) = NWAV_S
c
      RTN = PMAR(MODFIL,Fid_ref,' ',START3,DIM3_BUF,BUF_REF)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for MEAN REF','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_refsd,' ',START3,DIM3_BUF,BUF_REF_STD)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for REF_STD','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_AOT_Best,' ',START3,DIM3_BUF,
     &                 BUF_TAU_best)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for AOT_Best','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_AOT_Average,' ',START3,DIM3_BUF,
     &                BUF_TAU_average)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for AOT_Average','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_AOT_best_small,' ',START3,DIM3_BUF,
     &      BUF_TAUS_best)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &  (MODIS_W_GENERIC,'PMAR for AOT_best_small','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_AOT_average_small,' ',START3,DIM3_BUF,
     &     BUF_TAUS_average)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &  (MODIS_W_GENERIC,'PMAR for AOT_Average_small','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_AOT_best_large,' ',START3,DIM3_BUF,
     &         BUF_TAUB_best)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &  (MODIS_W_GENERIC,'PMAR for AOT_best_large','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_AOT_average_large,' ',START3,DIM3_BUF,
     &         BUF_TAUB_average)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &  (MODIS_W_GENERIC,'PMAR for AOT_average_large','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_assy_best,' ',START3,DIM3_BUF,
     &             BUF_ASSY_best)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for ASSY_best','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_assy_average,' ',START3,DIM3_BUF,
     &             BUF_ASSY_average)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for ASSY_average','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_bcscatt_best ,' ',START3,DIM3_BUF,
     &            BUF_BACK_best)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for BCSCATT_best','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_bcscatt_average ,' ',START3,DIM3_BUF,
     &            BUF_BACK_average)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for BCSCATT_average','MOD04_SEA_OUT')
c
c      RTN = PMAR(MODFIL,Fid_refF_best ,' ',START3,DIM3_BUF,
c     &            BUF_RefF_best)
c         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
c     &   (MODIS_W_GENERIC,'PMAR for RefFX_best','MOD04_SEA_OUT')
c
c      RTN = PMAR(MODFIL,Fid_refF_average ,' ',START3,DIM3_BUF,
c     &            BUF_RefF_average)
c         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
c     &   (MODIS_W_GENERIC,'PMAR for RefFX_average','MOD04_SEA_OUT')


CCCCCCCCCCCC   ELSE IF BUF_Line_No equals to the last scanline CCCCCCCCCCCCCCC

      ELSE IF (N_Entries.eq.Nrow) THEN
         Nline = N_Entries/Lines_Per_Write
         START2(2) = Nline * Lines_Per_Write
         START3_2(3)= Nline * Lines_Per_Write
         DIM2_BUF(2) = Buf_Line_No
         write(36,*) 'Ocean-Else',Nline,START2(1),START2(2),
     &       START3(1),START3(2), START3_2(3),N_Entries
c
c Dimensions (cells along swath,cells across swath)
c

      RTN = PMAR(MODFIL,SDS1_N,' ',START2,DIM2_BUF,BUF_LAT)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for 1st array','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,SDS2_N,' ',START2,DIM2_BUF,BUF_LON)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for 2nd array','MOD04_SEA_OUT')
C     
      RTN = PMAR(MODFIL,SDS1_NE1,' ',START2,DIM2_BUF,BUF_SST)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Scan_Start_Time','MOD04_SEA_OUT')
c
c
      RTN = PMAR(MODFIL,SDS1_NE2,' ',START2,DIM2_BUF,BUF_MTHET0)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Solar_Zenith','MOD04_SEA_OUT')
     
C
      RTN = PMAR(MODFIL,SDS1_NE3,' ',START2,DIM2_BUF,BUF_MPHI0)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Solar_Azimuth','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,SDS1_NE4,' ',START2,DIM2_BUF,BUF_MTHET)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Sensor_Zenith','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,SDS1_NE5,' ',START2,DIM2_BUF,BUF_MPHI)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Sensor_Azimuth','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_numpix,' ',START2,DIM2_BUF,BUF_NUMPIXELS)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for NUMBER OF PIXELS USED',
     &    'MOD04_SEA_OUT')
c
C      RTN = PMAR(MODFIL,Fid_scatt_angle,' ',START2,DIM2_BUF,BUF_SCATT)
C         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
C     &   (MODIS_W_GENERIC,'PMAR for Scattering Angle',
C     &    'MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_scattering_angle,' ',START2,DIM2_BUF
     &           ,BUF_SCATTANG)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for Scattering Angle',
     &    'MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_cldfrac,' ',START2,DIM2_BUF, BUF_CLDFRC)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for Cld fraction',
     &    'MOD04_SEA_OUT')
c
       RTN = PMAR(MODFIL,Fid_tau_land_and_ocean_img,' ',START2,DIM2_BUF,
     & BUF_550)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for aerosol tau image at 0.55 micron',
     &   'MOD04_SEA_OUT')
     
C2/2006
       RTN = PMAR(MODFIL,Fid_tau_land_and_ocean,' ',START2,DIM2_BUF,
     & BUF_550_LS)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for aerosol tau at 0.55 micron',
     &   'MOD04_SEA_OUT')
c 1/2002
      RTN = PMAR(MODFIL,Fid_Ratio_Small_Land_And_Ocean,' ',
     & START2,DIM2_BUF,BUF_small_wgt)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for small mode ratio Ocean and land',
     &   'MOD04_SEA_OUT')
     
     
C 3/2010
       RTN = PMAR(MODFIL,Fid_land_sea_Quality,' ',START2,DIM2_BUF,
     & BUF_Quality)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for land_sea_Quality',
     &   'MOD04_SEA_OUT')
c
c Dimensions (QA_ocean, cells across swath, cells along swath)
c

      DIM3_BUF_2(3) = Buf_Line_No

      RTN = PMAR(MODFIL,Fid_qa,' ',START3_2,DIM3_BUF_2,
     &             BUF_QCONTROL)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for QUALITY CONTROL','MOD04_SEA_OUT')

c
c Dimensions (cells along swath,cells across swath,solution number=2)
c

      START3(2) = Nline * Lines_Per_Write
      DIM3_BUF(2) = Buf_Line_No
      DIM3_BUF(3) = 1
c
      do 50 j = 1, NUM_solutions
            LL = Ncell*Lines_Per_Write*(j-1) + 1
            START3(3) = j - 1
      RTN = PMAR(MODFIL,Fid_smallmode,' ',START3,DIM3_BUF,
     &               BUF_small_weighting(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &  (MODIS_W_GENERIC,'PMAR fo SMALL MODE WEIGHT','MOD04_SEA_OUT')
 
      RTN = PMAR(MODFIL,Fid_leasterror,' ',START3,DIM3_BUF,
     &               BUF_Least_error(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Least sqare ','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_eff_rad,' ',START3,DIM3_BUF,
     &          Buf_effrad(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for EFF RAD','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_sol_small,' ',START3,DIM3_BUF,
     &           BUF_sol_INDX_small(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Indx Solution','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_sol_large,' ',START3,DIM3_BUF,
     &           BUF_sol_INDX_large(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Indx Solution','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_ccn,' ',START3,DIM3_BUF,BUF_ccn(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Cloud_Condensation_Nuclei',
     &   'MOD04_SEA_OUT')
c
          RTN = PMAR(MODFIL,Fid_mass_con,' ',START3,DIM3_BUF,
     &              BUF_mass_conc(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Mass_Concentration_Ocean',
     &   'MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_angs_coff1,' ',START3,DIM3_BUF,
     &              BUF_angs_coeff1(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Angstrom_Exponent_1_Ocean',
     &   'MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_angs_coff2,' ',START3,DIM3_BUF,
     &              BUF_angs_coeff2(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Angstrom_Exponent_2_Ocean',
     &   'MOD04_SEA_OUT')

50    continue
        do 150 j = 1, num_model_index
            LL = Ncell*Lines_Per_Write*(j-1) + 1
            START3(3) = j - 1
       RTN = PMAR(MODFIL,Fid_AOT_model,' ',START3,DIM3_BUF,
     &            BUF_AOT_model(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for Optical_Depth_Ratio_Small_Ocean',
     &   'MOD04_SEA_OUT')
150    continue

c
c Dimension (cells along swath,cells across swath,number of wavelength)
C

      do 60 j = 1, NWAV_S
            LL = Ncell*Lines_Per_Write*(j-1) + 1
            START3(3) = j - 1
c
      RTN = PMAR(MODFIL,Fid_ref,' ',START3,DIM3_BUF,BUF_REF(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for MEAN REF','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_refsd,' ',START3,DIM3_BUF,
     &         BUF_REF_STD(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for REF_STD','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_AOT_Best,' ',START3,DIM3_BUF,
     &                 BUF_TAU_best(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for AOT_Best','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_AOT_Average,' ',START3,DIM3_BUF,
     &                BUF_TAU_average(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for AOT_Average','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_AOT_best_small,' ',START3,DIM3_BUF,
     &      BUF_TAUS_best(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &  (MODIS_W_GENERIC,'PMAR for AOT_best_small','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_AOT_average_small,' ',START3,DIM3_BUF,
     &     BUF_TAUS_average(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &  (MODIS_W_GENERIC,'PMAR for AOT_Average_small','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_AOT_best_large,' ',START3,DIM3_BUF,
     &         BUF_TAUB_best(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &  (MODIS_W_GENERIC,'PMAR for AOT_best_large','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_AOT_average_large,' ',START3,DIM3_BUF,
     &         BUF_TAUB_average(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     & (MODIS_W_GENERIC,'PMAR for AOT_average_large','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_assy_best,' ',START3,DIM3_BUF,
     &             BUF_ASSY_best(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for ASSY_best','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_assy_average,' ',START3,DIM3_BUF,
     &             BUF_ASSY_average(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for ASSY_average','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_bcscatt_best ,' ',START3,DIM3_BUF,
     &            BUF_BACK_best(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for BCSCATT_best','MOD04_SEA_OUT')
c
      RTN = PMAR(MODFIL,Fid_bcscatt_average ,' ',START3,DIM3_BUF,
     &            BUF_BACK_average(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for BCSCATT_average','MOD04_SEA_OUT')
c
c      RTN = PMAR(MODFIL,Fid_refF_best ,' ',START3,DIM3_BUF,
c     &            BUF_RefF_best(LL))
c         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
c     &   (MODIS_W_GENERIC,'PMAR for RefFX_best','MOD04_SEA_OUT')
c
c      RTN = PMAR(MODFIL,Fid_refF_average ,' ',START3,DIM3_BUF,
c     &            BUF_RefF_average(LL))
c         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
c     &   (MODIS_W_GENERIC,'PMAR for RefFX_average','MOD04_SEA_OUT')
c
 
60    CONTINUE

C END OF SUBROUTINe MOD04S_OUT

      END IF

      RETURN
      END


C*********************************************************************
      SUBROUTINE MOD04L_OUT(MODFIL,Ncell,Nrow,SDSLAT,SDSLON,SDS_SST,
     *SDS_MTHET0,SDS_MTHET,SDS_MPHI0,SDS_MPHI,SDS_Scattering_Angle,
     *SDS_Tau_Land_Ocean_img,SDS_CldMskQA,SDS_Aerosol_Type,
     *SDS_SCAT_ANGLE_land,SDS_mass_conc_land,SDS_angs_coeff_land,
     *SDS_CLDFRC_land,SDS_dust_weighting,
     *SDS_est_uncer,SDS_RefF_land,SDS_TranF_land,SDS_NUMPIXELS_land,
     *SDSTAU_corrected,SDS_ref_land,SDS_ref_STD_land,SDS_QCONTROL_land,
     *SDS_Mean_Reflectance_Land_All,SDS_Tau_Land_Ocean,
     *SDS_SDev_Reflectance_Land_All,SDSTAU_small_land,     
     *SDS_Path_Radiance_Land,   
     *SDS_Critical_Reflectance_Land,  
     *SDS_Error_Crit_Reflectance_Land,     
     *SDS_Error_Path_Radiance_Land,  
     *SDS_QFlag_Critical_Ref_Land,SDS_QFlag_Path_Radiance_Land,
     *SDS_QCONTROL_CritRef_land,SDS_ratio_small_Land_Ocean,
     *SDS_Reflected_flux_Land_Ocean,SDSTAU_corrected_213,
     *SDS_Surface_Reflectance_Land,SDS_Fitting_Error_Land,
     *SDS_Quality_Land_Ocean)

C---------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C
C   Create and write Scientific Data Sets (SDS) for mod04 sea.
C   This sunroutine will be called every scan line but write data into
C   HDF file at every Lines_Per_Write.
C
C !INPUT PARAMETERS:
C
C    INTEGER   MODFIL    File handle structure for HDF files
C    INTEGER   Ncell     Number of cells across scan
C    INTEGER   Nrow      Number of (instrument) scans in granule
C        All arrays for land   SDS*
C                       Lat, Lon, Scantime, Solarzenith angle, View angle,
C                       Diff azimuth, Azimuth,
C
C   FOR LAND ONLY
C
C   SDS_QCONTROL_land      Quality control SDS array
C   SDS_Aerosol_Type       Index of Aerosol type
C   SDS_SCAT_ANGLE_land    Scattering Angle
C   SDS_angs_coeff_land    Angstrom exponent for 0.47 and 0.67 micron
C   SDS_CLDFRC_land        Cloud fraction (%)
C   SDS_dust_weighting     Dust aerosol weighting factor 
C   SDS_est_uncer          Uncertainty of optical thickness at 0.47 and 0.66 micron
C   SDS_RefF_land          Normalized reflected flux at 0.47 and 0.66 microns
C   SDS_TranF_land         Normalized Transmitted flux at 0.47 and 0.66 microns
C   SDS_NUMPIXELS_land     Number of pixels with desired percentile
C   SDSTAU_corrected       Corrected optical thickness at 0.47 0.55 and 0.66 micron
C   SDS_ref_land           Mean reflectance at five bands
C   SDS_ref_STD_land       Standard deviation of reflectance at five bands
C   SDS_mass_conc_land     Mass concentration
C
C !OUTPUT PARAMETERS:    NONE
C
C !REVISION HISTORY:
C Revision 1.1  1996/06/04  16:01:42  vlin
c 01/29/98 fhliang
c added NCSA acknowledgement.
c
C Initial revision
C
C !TEAM-UNIQUE HEADER:
C
C    This software is developed by the MODIS Science Data Support
C    Team for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
c HDF portions developed at the National Center for Supercomputing
c Applications at the University of Illinois at Urbana-Champaign.
c
C !REFERENCES and CREDITS:
C
C    Written by:
C       Vicky Lin                    June 1996
C       Research and Data systems Corporation
C       SAIC/GSC MODIS Science Data Support Office
C       7501 Forbes Blvd, Seabrook MD 20706
C       vlin@ltpmail.gsfc.nasa.gov
C
C    Modified By
C    Shana Mattoo                    November 1997
C    Allen Chu                       Decemberr 1997 to include
C                                    QA bytes
C
C !DESIGN NOTES:
C
C    This subroutine checks the return status of all MODIS Application
C    Program Interface (M-API) function calls.  A successful M-API
C    function call is indicated by a return status (RTN) equal to
C    MAPIOK.  If a M-API function call is not successful, a warning
C    message is written to the LogStatus file by passing an mnemonic
C    (.._W_..) to the function MODIS_SMF_SETDYNAMICMSG.
C
C    Externals:
C       Functions:
C          CRMAR   Creates an HDF SDS structure to store a new data
C                  array into a MODIS HDF file. (libmapi.a)
C          PMAR    PMAR places a multi-dimensional array of data
C                  into an HDF SDS array structure previously
C                  created using CRMAR.         (libmapi.a)
C
C       Named Constant:
C          MODIS_W_GENERIC       (PGS_MODIS_39500.f)
C          SDS*_N                (mod04.inc)
C          SCALE*                (mod04.inc)
C          OFFSET*               (mod04.inc)
C
C    Internal Subroutine:
C          MODIS_SetDimName  Write the dimension names of a SDS array.
C
C !END
C----------------------------------------------------------------------

      IMPLICIT NONE
      SAVE

      INCLUDE 'mod04.inc'
      INCLUDE 'mapi.inc'
      INCLUDE 'PGS_MODIS_39500.f'


      INTEGER I,J,Ncell,Nrow,RANK1,RANK3,RANK2,LL,LL1,LL2,LL3,
     &        N_Entries,RTN,BUF_Line_No,Nline,
     &        MODFIL(MODFILLEN),START3(3),START3_2(3),START2(2),
     &        DIM3_BUF(3),DIM3_BUF_2(3),DIM2_BUF(2),
     &        DIM_SDS(3),values(2)

C
C  Declare HDF SDS arrays
C
      REAL*8      SDS_SST(NUMCELLS)
      BYTE        SDS_QCONTROL_land(QA_LAND,NUMCELLS),
     &            SDS_CldMskQA(NUMCELLS)
      BYTE        SDS_QCONTROL_CritRef_land(QA_LAND,NUMCELLS)
      REAL        SDSLAT(NUMCELLS),SDSLON(NUMCELLS)
      REAL        SDS_mass_conc_land(NUMCELLS)
      INTEGER *2  SDS_MTHET0(NUMCELLS), SDS_MTHET(NUMCELLS),
     &            SDS_MPHI0(NUMCELLS),SDS_MPHI(NUMCELLS),
     &            SDS_Scattering_Angle(NUMCELLS),
     &            SDS_Tau_Land_Ocean_img(NUMCELLS),
     &            SDS_Tau_Land_Ocean(NUMCELLS),
     &            SDS_Aerosol_Type(NUMCELLS),
     &            SDS_SCAT_ANGLE_land(NUMCELLS),
     &            SDS_angs_coeff_land(NUMCELLS),
     &            SDS_CLDFRC_land(NUMCELLS),
     &            SDS_dust_weighting(NUMCELLS),
C  9/2005 ( two new SDS's)......     
     &            SDS_Surface_Reflectance_Land(NUMCELLS,Land_Sol3),
     &            SDS_Fitting_Error_Land(NUMCELLS)
      INTEGER *2
     &            SDS_NUMPIXELS_land(NUMCELLS,Land_Sol1),
     &            SDSTAU_corrected(NUMCELLS,Land_Sol3),
     &            SDSTAU_small_land(NUMCELLS,Land_Sol4),
     &            SDSTAU_corrected_213(NUMCELLS),
     &            SDS_ref_land(NUMCELLS,Band_land),
     &            SDS_ref_STD_land(NUMCELLS,Band_land),
     &            SDS_ratio_small_Land_Ocean(NUMCELLS),
     &           SDS_Reflected_flux_Land_Ocean(NUMCELLS) 
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
     & SDS_Quality_Land_Ocean(NUMCELLS)

C
C Obsolete (02/2006) Land SDS Arrays
C
      INTEGER *2   
     &            SDS_est_uncer(NUMCELLS,Land_Sol1),
     &            SDS_RefF_land(NUMCELLS,Land_Sol2),
     &            SDS_TranF_land(NUMCELLS,Land_Sol1)

C
C Declare buffer arrays to write to HDF file
C

      REAL*8   BUF_SST(NUMCELLS*Lines_Per_Write)
      BYTE     BUF_QCONTROL(QA_LAND*NUMCELLS*Lines_Per_Write),
     &         BUF_QCONTROL_crit(QA_LAND*NUMCELLS*Lines_Per_Write),
     &         BUF_CldMskQA(NUMCELLS*Lines_Per_Write)
      REAL     BUF_LAT(NUMCELLS*Lines_Per_Write),
     &         BUF_LON(NUMCELLS*Lines_Per_Write),
     &         BUF_mass_conc(NUMCELLS*Lines_Per_Write)

      INTEGER*2
     &         BUF_MTHET0(NUMCELLS*Lines_Per_Write),
     &         BUF_MTHET(NUMCELLS*Lines_Per_Write),
     &         BUF_MPHI0(NUMCELLS*Lines_Per_Write),
     &         BUF_MPHI(NUMCELLS*Lines_Per_Write),
     &         BUF_SCATTANG(NUMCELLS*Lines_Per_Write),
     &         BUF_550(NUMCELLS*Lines_Per_Write),
     &         BUF_550_LS(NUMCELLS*Lines_Per_Write),
     &         BUF_AERSOL(NUMCELLS*Lines_Per_Write),
     &         BUF_SCATT(NUMCELLS*Lines_Per_Write),
     &         BUF_angs_coeff(NUMCELLS*Lines_Per_Write),
     &         BUF_CLDFRC(NUMCELLS*Lines_Per_Write),
     &         BUF_DUSTWGT(NUMCELLS*Lines_Per_Write),
     &         BUF_TAU_CONT(NUMCELLS*Lines_Per_Write*Land_Sol1),
     &         BUF_est_unc(NUMCELLS*Lines_Per_Write*Land_Sol1),
     &         BUF_Fitting_Error(NUMCELLS*Lines_Per_Write),
     &         BUF_Surface_Reflectance(NUMCELLS*Lines_Per_Write*Land_Sol3),
     &         BUF_Quality(NUMCELLS*Lines_Per_Write)
      INTEGER*2
     &         BUF_RefF(NUMCELLS*Lines_Per_Write*Land_Sol2),
     &         BUF_TranF(NUMCELLS*Lines_Per_Write*Land_Sol1),
     &         BUF_NUMPIXELS(NUMCELLS*Lines_Per_Write*Land_Sol1),
     &         BUF_TAU_corrected(NUMCELLS*Lines_Per_Write*Land_Sol3),
     &         BUF_correc_wav213(NUMCELLS*Lines_Per_Write),
     &         BUF_small_TAU(NUMCELLS*Lines_Per_Write*Land_Sol4),
     &         BUF_REF(NUMCELLS*Lines_Per_Write*NWAV_S),
     &         BUF_REF_STD(NUMCELLS*Lines_Per_Write*NWAV_S),
     &         BUF_small_wgt(NUMCELLS*Lines_Per_Write),
     &         BUF_refflux(NUMCELLS*Lines_Per_Write)
      INTEGER *2 
     &BUF_Mean_Reflectance_Land_All
     &(NUMCELLS*Lines_Per_Write*Land_Sol3), 
     &BUF_SDev_Reflectance_Land_All
     &(NUMCELLS*Lines_Per_Write*Land_Sol3),     
     &BUF_Path_Radiance_Land(NUMCELLS*Lines_Per_Write*Land_Sol1),   
     &BUF_Slope_Land(NUMCELLS*Lines_Per_Write*Land_Sol1),  
     &BUF_Error_Slope_Land(NUMCELLS*Lines_Per_Write*Land_Sol1),      
     &BUF_Error_Path_Radiance_Land(NUMCELLS*Lines_Per_Write*Land_Sol1),  
     &BUF_Rsq_Fit_Land(NUMCELLS*Lines_Per_Write*Land_Sol1),
     &BUF_QualityFlag_Fit_Land(NUMCELLS*Lines_Per_Write*Land_Sol1)

      LOGICAL  Write_Flag
      DATA     N_Entries /0/
c          do  j=1,Ncell
c         write(42,*)'inside out',j,
c     *   SDS_Critical_Reflectance_Land(j,1),
c     *   SDS_Error_Crit_Reflectance_Land(j,1),
c     *   SDS_Critical_Reflectance_Land(j,2),
c     *   SDS_Error_Crit_Reflectance_Land(j,2)
c         enddo


C       WRITE(*,*) 'Ncell,Nrow,MODFIL,MODFILLEN'
C       WRITE(*,*) Ncell,Nrow,MODFIL,MODFILLEN
C       WRITE(*,*) 'SDS_Tau_Land_Ocean = in MOD04_OUT_V2.f '
C       WRITE(*,*) (SDS_Tau_Land_Ocean(i),i=1,135)
C
C Initialization
C

      DIM3_BUF(1)   = Ncell
      DIM3_BUF(2)   = Lines_Per_Write
      DIM3_BUF_2(1) = QA_Land
      DIM3_BUF_2(2) = Ncell
      DIM2_BUF(1)   = Ncell
      DIM2_BUF(2)   = Lines_Per_Write
      N_Entries = N_Entries + 1
      Write_Flag = .false.
      START2(1) = 0

      DO 10 I = 1, 3
         START3(I)   = 0
         START3_2(I) = 0
   10 CONTINUE

C
C s into dimensions 1 and 3,
C respectively, of 3-D write BUFfers.
C Dimension 2 is BUFfer (scan) line counter.
C

      BUF_Line_No = mod(N_Entries,Lines_Per_Write)

      if (BUF_Line_No.eq.0) then
          BUF_Line_No = Lines_Per_Write
          Write_Flag = .true.
      end if
      
Change
       do 15 j=1,Ncell
       do 20 i=1,QA_Land
         LL=QA_Land*Ncell*(Buf_Line_No-1)+QA_Land*(j-1)+i
         BUF_QCONTROL(LL)=SDS_QCONTROL_land(i,j)
         BUF_QCONTROL_crit(LL)=SDS_QCONTROL_CritRef_land(i,j)
   20  continue
   15  continue


       do 25 i = 1, Ncell
         LL = Ncell*(BUF_Line_No - 1) + i
         BUF_LAT(LL)           = SDSLAT(i)
         BUF_LON(LL)           = SDSLON(i)
         BUF_SST(LL)           = SDS_SST(i)
         BUF_MTHET0(LL)        = SDS_MTHET0(i)
         BUF_MTHET(LL)         = SDS_MTHET(i)
         BUF_MPHI0(LL)         = SDS_MPHI0(i)
         BUF_MPHI(LL)          = SDS_MPHI(i)
         BUF_550(LL)           = SDS_Tau_Land_Ocean_img(i)
         BUF_550_LS(LL)        = SDS_Tau_Land_Ocean(i)
         BUF_CldMskQA(LL)      = SDS_CldMskQA(i)
         BUF_AERSOL(LL)        = SDS_Aerosol_Type(i)
         BUF_SCATT(LL)         = SDS_SCAT_ANGLE_land(i)
         BUF_SCATTANG(LL)      = SDS_Scattering_Angle(i)
         BUF_mass_conc(LL)     = SDS_mass_conc_land(i)
         BUF_angs_coeff(LL)    = SDS_angs_coeff_land(i)
         BUF_CLDFRC(LL)        = SDS_CLDFRC_land(i)
         BUF_DUSTWGT(LL)       = SDS_dust_weighting(i)
         BUF_small_wgt(LL)=SDS_ratio_small_Land_Ocean(i) 
         BUF_Fitting_Error(LL)=SDS_Fitting_Error_Land(i)
         BUF_correc_wav213(LL)=SDSTAU_corrected_213(i)
         BUF_Quality(LL)=SDS_Quality_Land_Ocean(i)
   25 continue

       LL1=Ncell*(BUF_Line_No - 1)+1
       LL2=Ncell*(BUF_Line_No)

       do 30 j = 1,Land_Sol1
       do 35 i = 1, Ncell
         LL = Ncell*Lines_Per_Write*(j-1)
     &      + Ncell*(BUF_Line_No - 1) + i 
         BUF_NUMPIXELS(LL)     = SDS_NUMPIXELS_land(i,j)
      BUF_Path_Radiance_Land(LL)= SDS_Path_Radiance_Land(i,j)
      BUF_Slope_Land(LL)        = SDS_Critical_Reflectance_Land(i,j)
      BUF_Error_Slope_Land(LL) = SDS_Error_Crit_Reflectance_Land(i,j)
      BUF_Error_Path_Radiance_Land(LL)  = 
     &SDS_Error_Path_Radiance_Land(i,j)
      BUF_Rsq_Fit_Land(LL)=SDS_QFlag_Critical_Ref_Land(i,j)
      BUF_QualityFlag_Fit_Land(LL)=
     &SDS_QFlag_Path_Radiance_Land(i,j)
   35 continue
   30 continue
c
       do 50 j = 1,Land_Sol3
       do 55 i = 1, Ncell
         LL = Ncell*Lines_Per_Write*(j-1)
     &      + Ncell*(BUF_Line_No - 1) + i
c         BUF_RefF(LL)          = SDS_RefF_land(i,j)
         BUF_TAU_corrected(LL)  = SDSTAU_corrected(i,j)
   55 continue
   50 continue
c
      do 40 j = 1,Band_land
      do 45 i = 1, Ncell
         LL = Ncell*Lines_Per_Write*(j-1)
     &      + Ncell*(BUF_Line_No - 1) + i
c
         BUF_REF(LL)         = SDS_ref_land(i,j)
         BUF_REF_STD(LL)     = SDS_ref_STD_land(i,j)
   45 continue
   40 continue
       do 260 j = 1,Land_Sol3
       do 265 i = 1, Ncell
         LL = Ncell*Lines_Per_Write*(j-1)
     &      + Ncell*(BUF_Line_No - 1) + i
        BUF_Mean_Reflectance_Land_All(LL)= 
     & SDS_Mean_Reflectance_Land_All(i,j)
       BUF_SDev_Reflectance_Land_All(LL)= 
     &SDS_SDev_Reflectance_Land_All(i,j)
        
 265     continue
 260     continue
        do 150 j = 1,Land_Sol3
        do 155 i = 1, Ncell
         LL = Ncell*Lines_Per_Write*(j-1)
     &      + Ncell*(BUF_Line_No - 1) + i
      BUF_Surface_Reflectance(LL)=SDS_Surface_Reflectance_Land(i,j)
155       continue
150       continue
       do 350 j = 1,Land_Sol4
       do 355 i = 1, Ncell
         LL = Ncell*Lines_Per_Write*(j-1)
     &      + Ncell*(BUF_Line_No - 1) + i 
         BUF_small_TAU(LL)  = SDSTAU_small_land(i,j)
 355     continue
 350       continue

CCCCCCCCCCCCCCC   If BUF_Line_No equals to Lines_Per_Write   CCCCCCCCCCCCCC

       IF (Write_Flag) THEN

         Nline = N_Entries/Lines_Per_Write

         START2(2) = (Nline-1) * Lines_Per_Write
         START3(2) = (Nline-1) * Lines_Per_Write
         START3_2(3) = (Nline-1) * Lines_Per_Write

c
        write(36,*) 'Land',Nline,START2(1),START2(2),
     &       START3(1),START3(2), START3_2(3)
c Dimensions (cells along swath,cells across swath)
c

        RTN = PMAR(MODFIL,SDS1_N,' ',START2,DIM2_BUF,BUF_LAT)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for 1st array','MOD04_Land_OUT')
c
        RTN = PMAR(MODFIL,SDS2_N,' ',START2,DIM2_BUF,BUF_LON)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for 2nd array','MOD04_Land_OUT')
c
        RTN = PMAR(MODFIL,SDS1_NE1,' ',START2,DIM2_BUF,BUF_SST)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Scan_Start_Time',
     &      'MOD04_Land_OUT')
c
        RTN = PMAR(MODFIL,SDS1_NE2,' ',START2,DIM2_BUF,BUF_MTHET0)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Solar_Zenith','MOD04_Land_OUT')
c
        RTN = PMAR(MODFIL,SDS1_NE4,' ',START2,DIM2_BUF,BUF_MTHET)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Sensor_Zenith','MOD04_Land_OUT')
c
        RTN = PMAR(MODFIL,SDS1_NE3,' ',START2,DIM2_BUF,BUF_MPHI0)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Solar_Azimuth','MOD04_Land_OUT')
c
        RTN = PMAR(MODFIL,SDS1_NE5,' ',START2,DIM2_BUF,BUF_MPHI)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Sensor_Azimuth','MOD04_Land_OUT')
c
        RTN = PMAR(MODFIL,SDS1_NE6,' ',START2,DIM2_BUF,BUF_CldMskQA)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Cloud_Mask_QA','MOD04_Land_OUT')
C
        RTN = PMAR(MODFIL,Fid_tau_land_and_ocean_img,' ',START2,DIM2_BUF,BUF_550)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for aerosol image tau at 0.55 micron',
     &    'MOD04_Land_OUT')
C2/2006
        RTN = PMAR(MODFIL,Fid_tau_land_and_ocean,' ',START2,DIM2_BUF,BUF_550_LS)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for aerosol tau at 0.55 micron',
     &    'MOD04_Land_OUT')
c 1/2002
      RTN = PMAR(MODFIL,Fid_Ratio_Small_Land_And_Ocean,' ',
     & START2,DIM2_BUF,BUF_small_wgt)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for small mode ratio Ocean and land',
     &   'MOD04_SEA_OUT')
     
C   9/2005
      RTN = PMAR(MODFIL,Fid_Fitting_error,' ',
     & START2,DIM2_BUF,BUF_Fitting_Error)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for Fitting Error Land',
     &   'MOD04_SEA_OUT')
     
     
C 7/2009
       RTN = PMAR(MODFIL,Fid_land_sea_Quality,' ',START2,DIM2_BUF,
     & BUF_Quality)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for land_sea_Quality',
     &   'MOD04_SEA_OUT')
     
C   1/2006
        RTN = PMAR(MODFIL,Fid_AOT_correc_land213 ,' ',
     & START2,DIM2_BUF,BUF_correc_wav213)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for corrected tau for 213um',
     &   'MOD04_SEA_OUT')
 
c
c
        RTN = PMAR(MODFIL,Fid_aer_land,' ',START2,DIM2_BUF,BUF_AERSOL)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for AEROsol_land',
     &    'MOD04_Land_OUT')

        RTN = PMAR(MODFIL,Fid_scattering_angle,' ',START2,DIM2_BUF,
     &                BUF_SCATTANG)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for Scattering Angle',
     &    'MOD04_Land_OUT')
c
        RTN = PMAR(MODFIL,Fid_mass_con_land ,' ',START2,DIM2_BUF,
     &          BUF_mass_conc)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for Mass concen Land',
     &    'MOD04_Land_OUT')
c
        RTN = PMAR(MODFIL,Fid_angs_coff_land,' ',START2,DIM2_BUF,
     &          BUF_angs_coeff)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR forAng Coeff Land',
     &    'MOD04_Land_OUT')
c
        RTN = PMAR(MODFIL,Fid_cldfrac_land,' ',START2,DIM2_BUF,
     &          BUF_CLDFRC)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for Cloud Fraction Land',
     &    'MOD04_Land_OUT')
c
        RTN = PMAR(MODFIL,Fid_dust_wght_land,' ',START2,DIM2_BUF,
     &          BUF_DUSTWGT)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for Dust Weighting  Land',
     &    'MOD04_Land_OUT')

c
c Dimension(QA_Land, cells across swath, cells along swath)
c
      DIM3_BUF_2(3) = Lines_Per_Write

        RTN = PMAR(MODFIL,Fid_qa_land,' ',START3_2,DIM3_BUF_2,
     &             BUF_QCONTROL)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for QUALITY CONTROL','MOD04_Land_OUT')
c
c Dimension(QA_crit_ref Land, cells across swath, cells along swath)
c
      DIM3_BUF_2(3) = Lines_Per_Write

        RTN = PMAR(MODFIL,Fid_qa_extra_land,' ',START3_2,DIM3_BUF_2,
     &      BUF_QCONTROL_crit)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for QUALITY CONTROL','MOD04_Land_OUT')

c
c Dimension (cells along swath,cells across swath,SolLand=2)
c

      DIM3_BUF(3) = Land_Sol1

c
   
 
c
 
c
        RTN = PMAR(MODFIL,Fid_numpix_land,' ',START3,DIM3_BUF,
     &             BUF_NUMPIXELS)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR forNum Pixels_land ','MOD04_Land_OUT')
CNew array's Land......
          RTN = PMAR(MODFIL,Fid_path_radiance,' ',START3,DIM3_BUF,
     &               BUF_Path_Radiance_Land)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &(MODIS_W_GENERIC,'PMAR Path Radiance','MOD04_Land_OUT')

      RTN = PMAR(MODFIL,Fid_Critical_Reflectance,' ',START3,DIM3_BUF,
     &               BUF_Slope_Land)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &(MODIS_W_GENERIC,'PMAR Slope radiance','MOD04_Land_OUT')

       RTN = PMAR(MODFIL,Fid_error_Critical_Reflectance,' ',
     &          START3,DIM3_BUF,BUF_Error_Slope_Land)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &(MODIS_W_GENERIC,'PMAR Error Slope radiance','MOD04_Land_OUT')

        RTN = PMAR(MODFIL,Fid_error_path_radiance,' ',START3,DIM3_BUF,
     &    BUF_Error_Path_Radiance_Land)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &(MODIS_W_GENERIC,'PMAR Error Path_Radiance','MOD04_Land_OUT')

      RTN = PMAR(MODFIL,Fid_QWeight_Crit_Reflectance,' ',
     & START3,DIM3_BUF,BUF_Rsq_Fit_Land)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &(MODIS_W_GENERIC,'PMAR Correlation','MOD04_Land_OUT')

        RTN = PMAR(MODFIL,Fid_QWeight_Path_Radiance,' ',
     &  START3,DIM3_BUF, BUF_QualityFlag_Fit_Land)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &(MODIS_W_GENERIC,'PMAR Qyality flag','MOD04_Land_OUT')

       DIM3_BUF(3) = Land_Sol3

       RTN = PMAR(MODFIL,Fid_mean_ref_all,' ',START3,DIM3_BUF,
     &  BUF_Mean_Reflectance_Land_All)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &(MODIS_W_GENERIC,'PMAR Mean All pixels','MOD04_Land_OUT')

       RTN = PMAR(MODFIL,Fid_sd_ref_all,' ',START3,DIM3_BUF,
     &BUF_SDev_Reflectance_Land_All)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &(MODIS_W_GENERIC,'PMAR S Deviation All pixels','MOD04_Land_OUT')
C9/2005     
      RTN = PMAR(MODFIL,Fid_Surface_Reflectance,' ',START3,DIM3_BUF,
     &BUF_Surface_Reflectance)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &(MODIS_W_GENERIC,'PMAR Surface reflactance for Land','MOD04_Land_OUT')

c
c Dimensions (cells along swath,cells across swath,SolLand=3)
c

      DIM3_BUF(3) = Land_Sol3

        RTN = PMAR(MODFIL,Fid_AOT_correc_land,' ',START3,DIM3_BUF,
     &              BUF_TAU_corrected)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &  (MODIS_W_GENERIC,'PMAR fo AER Corrected_land','MOD04_Land_OUT')
     
     c
c Dimensions (cells along swath,cells across swath,SolLand=4)
c

      DIM3_BUF(3) = Land_Sol4

        RTN = PMAR(MODFIL,Fid_AOT_small_land,' ',START3,DIM3_BUF,
     &              BUF_small_TAU)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &  (MODIS_W_GENERIC,'PMAR fo Tau small Land','MOD04_Land_OUT')

 
c Dimension (cells along swath,cells across swath,number of wavelength_land)
c

      DIM3_BUF(3) = Band_land
c
        RTN = PMAR(MODFIL,Fid_ref_land,' ',START3,DIM3_BUF,BUF_REF)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for MEAN REF_land','MOD04_Land_OUT')
c
        RTN = PMAR(MODFIL,Fid_refsd_land,' ',START3,DIM3_BUF,
     &              BUF_REF_STD)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for REF_STD_land','MOD04_Land_OUT')
c


CCCCCCCCCCCC   ELSE IF BUF_Line_No equals to the last scanline CCCCCCCCCCCCCCC

      ELSE IF (N_Entries.eq.Nrow) THEN

         Nline = N_Entries/Lines_Per_Write
         START2(2) = Nline * Lines_Per_Write
         START3_2(3)= Nline * Lines_Per_Write
         DIM2_BUF(2) = Buf_Line_No
         
          write(36,*) 'Land-Else',Nline,START2(1),START2(2),
     &       START3(1),START3(2), START3_2(3),N_Entries
c

c
c Dimensions (cells along swath,cells across swath)
c

        RTN = PMAR(MODFIL,SDS1_N,' ',START2,DIM2_BUF,BUF_LAT)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for 1st array','MOD04_Land_OUT')
c
        RTN = PMAR(MODFIL,SDS2_N,' ',START2,DIM2_BUF,BUF_LON)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for 2nd array','MOD04_Land_OUT')
c
        RTN = PMAR(MODFIL,SDS1_NE1,' ',START2,DIM2_BUF,BUF_SST)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Scan_Start_Time',
     &      'MOD04_Land_OUT')
C
        RTN = PMAR(MODFIL,SDS1_NE2,' ',START2,DIM2_BUF,BUF_MTHET0)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Solar_Zenith','MOD04_Land_OUT')
c
        RTN = PMAR(MODFIL,SDS1_NE4,' ',START2,DIM2_BUF,BUF_MTHET)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Sensor_Zenith','MOD04_Land_OUT')
C
        RTN = PMAR(MODFIL,SDS1_NE3,' ',START2,DIM2_BUF,BUF_MPHI0)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Solar_Azimuth','MOD04_Land_OUT')
c
        RTN = PMAR(MODFIL,SDS1_NE5,' ',START2,DIM2_BUF,BUF_MPHI)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Sensor_Azimuth','MOD04_Land_OUT')
c
        RTN = PMAR(MODFIL,SDS1_NE6,' ',START2,DIM2_BUF,BUF_CldMskQA)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Cloud_Mask_QA','MOD04_Land_OUT')
C
        RTN = PMAR(MODFIL,Fid_tau_land_and_ocean_img,' ',START2,DIM2_BUF,BUF_550)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for aerosol image tau at 0.55 micron',
     &    'MOD04_Land_OUT')
C2/2006
        RTN = PMAR(MODFIL,Fid_tau_land_and_ocean,' ',START2,DIM2_BUF,BUF_550_LS)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for aerosol tau at 0.55 micron',
     &    'MOD04_Land_OUT')
c 1/2002
         RTN = PMAR(MODFIL,Fid_Ratio_Small_Land_And_Ocean,' ',
     & START2,DIM2_BUF,BUF_small_wgt)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for small mode ratio Ocean and land',
     &   'MOD04_SEA_OUT')
     
C   9/2005
         RTN = PMAR(MODFIL,Fid_Fitting_error,' ',
     & START2,DIM2_BUF,BUF_Fitting_Error)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for Fitting Error Land',
     &   'MOD04_SEA_OUT')

C   1/2006
         RTN = PMAR(MODFIL,Fid_AOT_correc_land213,' ',
     & START2,DIM2_BUF, BUF_correc_wav213)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for Corrected Taua 213um Land',
     &   'MOD04_SEA_OUT')
              
c
        RTN = PMAR(MODFIL,Fid_aer_land,' ',START2,DIM2_BUF,BUF_AERSOL)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for AEROsol_land',
     &    'MOD04_Land_OUT')
c
C        RTN = PMAR(MODFIL,Fid_scatt_angle_land,' ',START2,DIM2_BUF,
C     &                BUF_SCATT)
C         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
C     &   (MODIS_W_GENERIC,'PMAR for Scattering Angle_land',
C     &    'MOD04_Land_OUT')
c
        RTN = PMAR(MODFIL,Fid_scattering_angle,' ',START2,DIM2_BUF,
     &                BUF_SCATTANG)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for Scattering Angle_land',
     &    'MOD04_Land_OUT')
c
        RTN = PMAR(MODFIL,Fid_mass_con_land ,' ',START2,DIM2_BUF,
     &          BUF_mass_conc)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for Mass concen Land',
     &    'MOD04_Land_OUT')
c
        RTN = PMAR(MODFIL,Fid_angs_coff_land,' ',START2,DIM2_BUF,
     &          BUF_angs_coeff)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR forAng Coeff Land',
     &    'MOD04_Land_OUT')
c
        RTN = PMAR(MODFIL,Fid_cldfrac_land,' ',START2,DIM2_BUF,
     &          BUF_CLDFRC)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for Cloud Fraction Land',
     &    'MOD04_Land_OUT')
c
        RTN = PMAR(MODFIL,Fid_dust_wght_land,' ',START2,DIM2_BUF,
     &          BUF_DUSTWGT)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for Dust Weighting  Land',
     &    'MOD04_Land_OUT')

         RTN = PMAR(MODFIL,Fid_land_sea_Quality,' ',START2,DIM2_BUF,
     &  BUF_Quality)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for land_sea_Quality',
     &   'MOD04_SEA_OUT')
c Dimensions (QA_Land, cells across swath, cells along swath)
c

        DIM3_BUF_2(3) = Buf_Line_No

        RTN = PMAR(MODFIL,Fid_qa_land,' ',START3_2,DIM3_BUF_2,
     &             BUF_QCONTROL)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for QUALITY CONTROL','MOD04_Land_OUT')

        RTN = PMAR(MODFIL,Fid_qa_extra_land,' ',START3_2,DIM3_BUF_2,
     &     BUF_QCONTROL_crit)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR for QUALITY CONTROL','MOD04_Land_OUT')

c
c Dimensions (cells along swath,cells across swath,landsol1=2)
c
         START3(2) = Nline * Lines_Per_Write
         DIM3_BUF(2) = Buf_Line_No
         DIM3_BUF(3) = 1
c
          do 60 j = 1,Land_Sol1
              LL = Ncell*Lines_Per_Write*(j-1) + 1
             START3(3) = j - 1
c


        RTN = PMAR(MODFIL,Fid_numpix_land,' ',START3,DIM3_BUF,
     &             BUF_NUMPIXELS(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &   (MODIS_W_GENERIC,'PMAR forNum Pixels_land ','MOD04_Land_OUT')

CNew array's Land......
          RTN = PMAR(MODFIL,Fid_path_radiance,' ',START3,DIM3_BUF,
     &               BUF_Path_Radiance_Land(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &(MODIS_W_GENERIC,'PMAR Path Radiance','MOD04_Land_OUT')

       RTN = PMAR(MODFIL,Fid_Critical_Reflectance,' ',START3,DIM3_BUF,
     &               BUF_Slope_Land(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &(MODIS_W_GENERIC,'PMAR Slope radiance','MOD04_Land_OUT')


        RTN = PMAR(MODFIL,Fid_error_path_radiance,' ',START3,DIM3_BUF,
     &    BUF_Error_Path_Radiance_Land(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &(MODIS_W_GENERIC,'PMAR Error Path_Radiance','MOD04_Land_OUT')

      RTN = PMAR(MODFIL,Fid_error_Critical_Reflectance,' ',
     &      START3,DIM3_BUF,BUF_Error_Slope_Land(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &(MODIS_W_GENERIC,'PMAR Error Slope radiance','MOD04_Land_OUT')

      RTN = PMAR(MODFIL,Fid_QWeight_Path_Radiance,' ',
     &  START3,DIM3_BUF,BUF_QualityFlag_Fit_Land(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &(MODIS_W_GENERIC,'PMAR Qyality flag','MOD04_Land_OUT')


       RTN = PMAR(MODFIL,Fid_QWeight_Crit_Reflectance,' ',
     &  START3,DIM3_BUF,BUF_Rsq_Fit_Land(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &(MODIS_W_GENERIC,'PMAR Correlation','MOD04_Land_OUT')

      
60        continue

       do 160 j = 1,Land_Sol3
            LL = Ncell*Lines_Per_Write*(j-1) + 1
            START3(3) = j - 1
       RTN = PMAR(MODFIL,Fid_mean_ref_all,' ',START3,DIM3_BUF,
     &  BUF_Mean_Reflectance_Land_All(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &(MODIS_W_GENERIC,'PMAR Mean All pixels','MOD04_Land_OUT')

       RTN = PMAR(MODFIL,Fid_sd_ref_all,' ',START3,DIM3_BUF,
     &BUF_SDev_Reflectance_Land_All(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &(MODIS_W_GENERIC,'PMAR S Deviation All pixels','MOD04_Land_OUT')
     
       RTN = PMAR(MODFIL,Fid_Surface_Reflectance,' ',START3,DIM3_BUF,
     &BUF_Surface_Reflectance(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &(MODIS_W_GENERIC,'PMAR Surface reflactance for Land','MOD04_Land_OUT')

160    continue
c         write(36,*) 'two',RTN
c
c   Dimensions (cells along swath,cells across swath,landsol2=3)
c

c
         do 70 j = 1,Land_Sol3
            LL = Ncell*Lines_Per_Write*(j-1) + 1
            START3(3) = j - 1
          RTN = PMAR(MODFIL,Fid_AOT_correc_land,' ',START3,DIM3_BUF,
     &          BUF_TAU_corrected(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &  (MODIS_W_GENERIC,'PMAR fo AER Corrected_land','MOD04_Land_OUT')
c
 
c
70        continue
         do 170 j = 1,Land_Sol4
            LL = Ncell*Lines_Per_Write*(j-1) + 1
            START3(3) = j - 1
          RTN = PMAR(MODFIL,Fid_AOT_small_land,' ',START3,DIM3_BUF,
     &          BUF_small_TAU(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &  (MODIS_W_GENERIC,'PMAR fo Tau small land','MOD04_Land_OUT')
c
 
c
170        continue


c
c dimension( cells along swath,cells across swath,number of wavelength_land)
c

         do 80 j = 1, Band_LAND
            LL = Ncell*Lines_Per_Write*(j-1) + 1
            START3(3) = j - 1
          RTN = PMAR(MODFIL,Fid_ref_land,' ',START3,DIM3_BUF,
     &          BUF_REF(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for MEAN REF_land','MOD04_Land_OUT')
c
          RTN = PMAR(MODFIL,Fid_refsd_land,' ',START3,DIM3_BUF,
     &          BUF_REF_STD(LL))
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for REF_STD_land','MOD04_Land_OUT')
   80     continue

C END OF SUBROUTINe MOD04S_OUT

       END IF

      RETURN
      END




