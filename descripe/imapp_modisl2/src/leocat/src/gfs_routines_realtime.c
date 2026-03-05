/*$Id: gfs_routines.c,v 1.1.1.1 2006/10/02 22:40:37 mpav Exp $*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <float.h>
#include <hdf.h>
#include <mfhdf.h>
#include "common_leocat.h"
#include "wgrib.h"
#include "nwp_leocat.h"

int compare(const void *one, const void *two)
{
    return *(int*)one - *(int*)two;
}

void decode_gfs_grib(char *, nwp_params_tmp *, float *);
void get_grid_info(unsigned char *, int, float *, float *,
                   float *, float *, float *, float *, char *);

void main_gfs(char *gfs_dir, int year, int jday, float time, float *nwp_bounds,
              nwp_params *nwp)
{

    char r_name[256], filename1[256], filename2[256], pflg;
    int year1, year2, jday1, jday2, hour1, hour2, month1, month2, day1, day2;
    int ileap1, ileap2;
    int n, ilev, kfirst_rh, klast_rh, kfirst_o3, klast_o3;
    const float forecast_int = 6.0;
    float x, temp;
    double *p_tmp, *t_tmp, *rh_tmp, *w_tmp, *o3_tmp, *z_prof, *p_prof, *t_prof, *w_prof, wo;
    float *o3_prof;
    nwp_params_tmp gfs1, gfs2;

    printf("Inside main_gfs()\n"); // GPC

    strcpy(r_name,"main_gfs");

    kfirst_rh=0;
    klast_rh=0;
    kfirst_o3=0;
    klast_o3=0;

    /*---------------------------------------------------------------------------
    Determine the needed hour and julian day for each gfs file, keeping in mind
    that we are using 12 hour forecasts.
    ---------------------------------------------------------------------------*/

    //printf("year=%d, jday=%d, time=%f\n",year,jday,time); //RMC
    //printf("bounds = %f %f %f %f\n", *(nwp_bounds+0), *(nwp_bounds+1), *(nwp_bounds+2), *(nwp_bounds+3)); //RMC

    get_nwp_time(year, jday, time, 12, &year1, &year2, &ileap1,
                 &ileap2, &jday1, &jday2, &hour1, &hour2);

    /*---------------------------------------------------------------------------
    Determine the correct year, month, and day for each gfs file.
    ---------------------------------------------------------------------------*/

    month1 = jday2month(jday1, ileap1);
    month2 = jday2month(jday2, ileap2);
    day1 = jday2day(jday1, ileap1);
    day2 = jday2day(jday2, ileap2);

    printf("First Forecast Start Time\n");
    printf("year=%d, month=%d, day=%d, hour=%d\n",year1,month1,day1,hour1);
    printf("Second Forecast Start Time\n");
    printf("year=%d, month=%d, day=%d, hour=%d\n",year2,month2,day2,hour2);

    /*---------------------------------------------------------------------------
    Construct the gfs file names.
    ---------------------------------------------------------------------------*/

    sprintf(filename1,"%s/gblav.%2.2d%2.2d%2.2d%2.2d_F012",gfs_dir,year1 - 100*(year1/100),month1,day1,hour1);
    sprintf(filename2,"%s/gblav.%2.2d%2.2d%2.2d%2.2d_F012",gfs_dir,year2 - 100*(year2/100),month2,day2,hour2);

    printf("GFS file #1 = %s\n",filename1);
    printf("GFS file #2 = %s\n",filename2);

    /*---------------------------------------------------------------------------
    Read in gfs from grib files for both times needed to interpolate to satellite
    overpass time.
    ---------------------------------------------------------------------------*/

    decode_gfs_grib(filename1, &gfs1, nwp_bounds);
    decode_gfs_grib(filename2, &gfs2, nwp_bounds);

    /*---------------------------------------------------------------------------
    Check to make sure that NWP at time1 is consistent with NWP at time 2.
    ---------------------------------------------------------------------------*/

    if (gfs1.nlon != gfs2.nlon || gfs1.nlat != gfs2.nlat ||
            gfs1.dlat != gfs2.dlat || gfs1.dlon != gfs2.dlon ||
            gfs1.first_lat != gfs2.first_lat || gfs1.last_lat != gfs2.last_lat ||
            gfs1.first_lon != gfs2.first_lon || gfs1.last_lat != gfs2.last_lat ||
            gfs1.nlevels != gfs2.nlevels || gfs1.nlevels_rh != gfs2.nlevels_rh ||
            gfs1.nlevels_o3 != gfs2.nlevels_o3) {
        fprintf(stderr,"GFS forecast grids are different\n");
        fprintf(stderr,"%f  %f\n",gfs1.first_lat,gfs2.first_lat);
        exit(1);
    }

    printf("%s\n%f  %f  %f\n%f  %f  %f\n",gfs1.grid_name,gfs1.first_lat,gfs1.last_lat,gfs1.dlat,
           gfs1.first_lon,gfs1.last_lon,gfs1.dlon);

    /*---------------------------------------------------------------------------
    Copy geolocation and dimension data into main NWP structure.
    ---------------------------------------------------------------------------*/

    strcpy(nwp->grid_name,gfs1.grid_name);
    nwp->nlon = gfs1.nlon;
    nwp->nlat = gfs1.nlat;
    nwp->first_lat = gfs1.first_lat;
    nwp->last_lat = gfs1.last_lat;
    nwp->first_lon = gfs1.first_lon;
    nwp->last_lon = gfs1.last_lon;
    nwp->dlat = gfs1.dlat;
    nwp->dlon = gfs1.dlon;

    nwp->ncells = (long)gfs1.nlat*(long)gfs1.nlon;
    nwp->nlevels = MAX_LAYERS_CRTM + 1;
    nwp->nlayers = MAX_LAYERS_CRTM;

    /*---------------------------------------------------------------------------
    Allocate memory for temporary low vertical resolution profiles.
    ---------------------------------------------------------------------------*/

    if ((p_tmp = (double *) calloc(gfs1.nlevels+1,sizeof(double))) == NULL)
        error_allo(r_name,"p_tmp");
    if ((t_tmp = (double *) calloc(gfs1.nlevels+1,sizeof(double))) == NULL)
        error_allo(r_name,"t_tmp");
    if ((rh_tmp = (double *) calloc(gfs1.nlevels+1,sizeof(double))) == NULL)
        error_allo(r_name,"rh_tmp");
    if ((w_tmp = (double *) calloc(gfs1.nlevels+1,sizeof(double))) == NULL)
        error_allo(r_name,"w_tmp");
    if ((o3_tmp = (double *) calloc(gfs1.nlevels+1,sizeof(double))) == NULL)
        error_allo(r_name,"o3_tmp");

    /*---------------------------------------------------------------------------
    Allocate memory for temporary high vertical resolution profiles.
    ---------------------------------------------------------------------------*/

    if ((z_prof = (double *) calloc(nwp->nlevels, sizeof(double))) == NULL)
        error_allo(r_name,"z_prof");
    if ((p_prof = (double *) calloc(nwp->nlevels, sizeof(double))) == NULL)
        error_allo(r_name,"p_prof");
    if ((t_prof = (double *) calloc(nwp->nlevels, sizeof(double))) == NULL)
        error_allo(r_name,"t_prof");
    if ((w_prof = (double *) calloc(nwp->nlevels, sizeof(double))) == NULL)
        error_allo(r_name,"w_prof");
    if ((o3_prof = (float *) calloc(nwp->nlevels, sizeof(float))) == NULL)
        error_allo(r_name,"o3_prof");

    /*---------------------------------------------------------------------------
    Allocate memory for nwp grid cell data structure.
    ---------------------------------------------------------------------------*/

    if ((nwp->map = (profile_params *) malloc(sizeof(profile_params) * nwp->ncells)) == NULL)
        error_allo(r_name,"nwp->map*");

    for (n=0; n<nwp->ncells; n++) {
        if ((nwp->map[n].z_lev = (float *) malloc(sizeof(float) * nwp->nlevels)) == NULL)
            error_allo(r_name,"nwp->map.z_lev");
        if ((nwp->map[n].p_lev = (float *) malloc(sizeof(float) * nwp->nlevels)) == NULL)
            error_allo(r_name,"nwp->map.p_lev");
        if ((nwp->map[n].t_lev = (float *) malloc(sizeof(float) * nwp->nlevels)) == NULL)
            error_allo(r_name,"nwp->map.t_lev");
        if ((nwp->map[n].w_lev = (float *) malloc(sizeof(float) * nwp->nlevels)) == NULL)
            error_allo(r_name,"nwp->map.w_lev");
        if ((nwp->map[n].o3_lev = (float *) malloc(sizeof(float) * nwp->nlevels)) == NULL)
            error_allo(r_name,"nwp->map.o3_lev");
    }

    /*---------------------------------------------------------------------------
    Determine the time interpolation weight.
    ---------------------------------------------------------------------------*/

    temp = ((int)((time+forecast_int)/forecast_int)*forecast_int);
    x = (temp < (time+forecast_int)) ? (temp - time)/forecast_int : ((temp-forecast_int)-time)/forecast_int;

    /*---------------------------------------------------------------------------
    Compute weighted average of 2 GFS times for surface and column parameters.
    ---------------------------------------------------------------------------*/

    for (n=0; n<nwp->ncells; n++) {

        nwp->map[n].p_sfc = (gfs1.sfc_pres[n] >= 0.0 && gfs2.sfc_pres[n] >= 0.0) ?
                            (x*(gfs1.sfc_pres[n]/100.0) + (1.0-x)*(gfs2.sfc_pres[n]/100.0)) : MISSING_FLOAT;
        nwp->map[n].p_msl = (gfs1.p_msl[n] >= 0.0 && gfs2.p_msl[n] >= 0.0) ?
                            (x*(gfs1.p_msl[n]/100.0) + (1.0-x)*(gfs2.p_msl[n]/100.0)) : MISSING_FLOAT;
        nwp->map[n].rh_2m = (gfs1.sfc_rh[n] >= 0.0 && gfs2.sfc_rh[n] >= 0.0) ?
                            (x*gfs1.sfc_rh[n] + (1.0-x)*gfs2.sfc_rh[n]) : MISSING_FLOAT;
        nwp->map[n].t_2m = (gfs1.t_2m[n] >= 0.0 && gfs2.t_2m[n] >= 0.0) ?
                           (x*gfs1.t_2m[n] + (1.0-x)*gfs2.t_2m[n]) : MISSING_FLOAT;
        nwp->map[n].z_sfc = (gfs1.sfc_z[n] != MISSING_GFS || gfs2.sfc_z[n] != MISSING_GFS) ?
                            (x*gfs1.sfc_z[n] + (1.0-x)*gfs2.sfc_z[n]) : MISSING_FLOAT;
        nwp->map[n].t_sfc = (gfs1.sfc_t[n] >= 0.0 && gfs2.sfc_t[n] >= 0.0) ?
                            (x*gfs1.sfc_t[n] + (1.0-x)*gfs2.sfc_t[n]) : MISSING_FLOAT;
        nwp->map[n].snow_sfc = (gfs1.sfc_snow[n] >= 0.0 && gfs2.sfc_snow[n] >= 0.0) ?
                               (x*gfs1.sfc_snow[n] + (1.0-x)*gfs2.sfc_snow[n]) : MISSING_FLOAT;
        nwp->map[n].alb_sfc = (gfs1.sfc_alb[n] >= 0.0 && gfs2.sfc_alb[n] >= 0.0) ?
                              (x*gfs1.sfc_alb[n] + (1.0-x)*gfs2.sfc_alb[n]) : MISSING_FLOAT;

        nwp->map[n].t_tropo = (gfs1.tropo_t[n] >= 0.0 && gfs2.tropo_t[n] >= 0.0) ?
                              (x*gfs1.tropo_t[n] + (1.0-x)*gfs2.tropo_t[n]) : MISSING_FLOAT;
        nwp->map[n].p_tropo = (gfs1.tropo_p[n] >= 0.0 && gfs2.tropo_p[n] >= 0.0) ?
                              (x*gfs1.tropo_p[n] + (1.0-x)*gfs2.tropo_p[n]) : MISSING_FLOAT;
        nwp->map[n].z_tropo = (gfs1.tropo_z[n] >= 0.0 && gfs2.tropo_z[n] >= 0.0) ?
                              (x*gfs1.tropo_z[n] + (1.0-x)*gfs2.tropo_z[n]) : MISSING_FLOAT;

        nwp->map[n].o3_col = (gfs1.col_o3[n] >= 0.0 && gfs2.col_o3[n] >= 0.0) ?
                             (x*gfs1.col_o3[n] + (1.0-x)*gfs2.col_o3[n]) : MISSING_FLOAT;
        nwp->map[n].tpw = (gfs1.col_h20[n] >= 0.0 && gfs2.col_h20[n] >= 0.0) ?
                          (x*gfs1.col_h20[n] + (1.0-x)*gfs2.col_h20[n]) : MISSING_FLOAT;

        /*---------------------------------------------------------------------------
          Check to see if the surface pressure is greater than the lowest level
          atmospheric pressure.  If so, additional level is surface; otherwise
          bottom new level is set to MISSING.
          ---------------------------------------------------------------------------*/
        if (nwp->map[n].p_sfc > 0.0 && gfs1.pressure[n][gfs1.nlevels-1] <= nwp->map[n].p_sfc) {
            nwp->map[n].sfc_level = gfs1.nlevels;
        } else {
            nwp->map[n].sfc_level = 0;
        }

        nwp->map[n].lat = gfs1.lat[n];
        nwp->map[n].lon = gfs1.lon[n];
    }

    /*---------------------------------------------------------------------------
    Compute weighted average of 2 GFS times for full profile parameters.
    ---------------------------------------------------------------------------*/

    /*---------------------------------------------------------------------------
    First determine the vertical indices of the RH and O3 parameters with respect
    to the full profile resolution.
    ---------------------------------------------------------------------------*/

    for (ilev=0; ilev<gfs1.nlevels; ilev++) {
        if (gfs1.pressure_rh[0][0] == gfs1.pressure[0][ilev]) kfirst_rh = ilev;
        if (gfs1.pressure_o3[0][0] == gfs1.pressure[0][ilev]) kfirst_o3 = ilev;
        if (gfs1.pressure_rh[0][gfs1.nlevels_rh-1] == gfs1.pressure[0][ilev]) klast_rh = ilev;
        if (gfs1.pressure_o3[0][gfs1.nlevels_o3-1] == gfs1.pressure[0][ilev]) klast_o3 = ilev;
    }

    /*---------------------------------------------------------------------------
    Create the 101 level pressure profile.
    ---------------------------------------------------------------------------*/

    make_profile_101(p_prof);

    /*---------------------------------------------------------------------------
    Start loop over each NWP grid cell.
    ---------------------------------------------------------------------------*/

    pflg = 0;
    printf("GDAS 1: Looping over %ld NWP grid cells...\n",nwp->ncells);
    for (n=0; n<nwp->ncells; n++) {
        /*if (nwp->sfc_level[n] == 0) pflg = 1;*/
        /*if (n == 620) pflg = 1;*/

        /*---------------------------------------------------------------------------
        Loop over original number of gfs levels.
        ---------------------------------------------------------------------------*/

        for (ilev=0; ilev<gfs1.nlevels; ilev++) {

            /*---------------------------------------------------------------------------
            Pressure - atmospheric levels are fixed.
            ---------------------------------------------------------------------------*/

            p_tmp[ilev] = gfs1.pressure[n][ilev];

            /*---------------------------------------------------------------------------
            Temperature - full vertical resolution.
            ---------------------------------------------------------------------------*/
            t_tmp[ilev] = (gfs1.t[n][ilev] >= 0.0 && gfs2.t[n][ilev] >= 0.0) ?
                          (x*gfs1.t[n][ilev] + (1.0-x)*gfs2.t[n][ilev]) : MISSING_FLOAT;
            if (gfs1.pressure[n][ilev] != gfs2.pressure[n][ilev]) {
                fprintf(stderr,"Pressure levels in GFS files differ\n");
                exit(1);
            }

            /*---------------------------------------------------------------------------
            Relative humidity (in g/kg) - partial vertical resolution.
            ---------------------------------------------------------------------------*/
            rh_tmp[ilev] = 0.0;
            w_tmp[ilev] = 0.0;
            if (ilev >= kfirst_rh && ilev <= klast_rh) {
                rh_tmp[ilev] = (gfs1.rh[n][ilev-kfirst_rh] >= 0.0 && gfs2.rh[n][ilev-kfirst_rh] >= 0.0) ?
                               (x*gfs1.rh[n][ilev-kfirst_rh] + (1.0-x)*gfs2.rh[n][ilev-kfirst_rh]) : MISSING_FLOAT;
                if (gfs1.pressure_rh[n][ilev-kfirst_rh] != gfs2.pressure_rh[n][ilev-kfirst_rh]) {
                    fprintf(stderr,"RH pressure levels in GFS files differ\n");
                    exit(1);
                }
                if (rh_tmp[ilev] > 0.0) w_tmp[ilev] = rh_to_mr_wat(rh_tmp[ilev], p_tmp[ilev], t_tmp[ilev]);
            }

            /*---------------------------------------------------------------------------
            Ozone - partial vertical resolution.
            ---------------------------------------------------------------------------*/
            o3_tmp[ilev] = 0.0;
            if (ilev >= kfirst_o3 && ilev <= klast_o3) {
                o3_tmp[ilev] = (gfs1.o3[n][ilev-kfirst_o3] >= 0.0 && gfs2.o3[n][ilev-kfirst_o3] >= 0.0) ?
                               (x*gfs1.o3[n][ilev-kfirst_o3] + (1.0-x)*gfs2.o3[n][ilev-kfirst_o3]) : MISSING_FLOAT;
                if (gfs1.pressure_o3[n][ilev-kfirst_o3] != gfs2.pressure_o3[n][ilev-kfirst_o3]) {
                    fprintf(stderr,"Ozone pressure levels in GFS files differ\n");
                    exit(1);
                }
            }

            /*---------------------------------------------------------------------------
            Check to see if the surface is encountered before lowest non-surface
            pressure level.
            ---------------------------------------------------------------------------*/

            if (nwp->map[n].p_sfc > 0.0 && p_tmp[ilev] > nwp->map[n].p_sfc) {
                if (nwp->map[n].sfc_level == 0) {
                    nwp->map[n].sfc_level = ilev;
                    t_tmp[ilev] = nwp->map[n].t_2m;
                    p_tmp[ilev] = nwp->map[n].p_sfc;
                    o3_tmp[ilev] = 0.0;
                    rh_tmp[ilev] = nwp->map[n].rh_2m;
                    if (rh_tmp[ilev] > 0.0) w_tmp[ilev] = rh_to_mr_wat(rh_tmp[ilev], p_tmp[ilev], t_tmp[ilev]);
                } else {
                    t_tmp[ilev] = nwp->map[n].t_2m;
                    p_tmp[ilev] = gfs1.pressure[n][ilev-1];
                    o3_tmp[ilev] = 0.0;
                    rh_tmp[ilev] = nwp->map[n].rh_2m;
                    w_tmp[ilev] = w_tmp[ilev-1];
                }
            }
            if (pflg) printf("ilev=%d, p=%f, t=%f, rh=%f, w=%f O3=%f\n",ilev,p_tmp[ilev],t_tmp[ilev],rh_tmp[ilev],w_tmp[ilev],o3_tmp[ilev]);
        }

        /*---------------------------------------------------------------------------
        Add on extra level since we are adding the 2-meter parameters to the profile.
        ---------------------------------------------------------------------------*/

        if (nwp->map[n].sfc_level != gfs1.nlevels) {
            p_tmp[gfs1.nlevels] = gfs1.pressure[n][ilev-1];
        } else {
            p_tmp[gfs1.nlevels] = nwp->map[n].p_sfc;
        }

        t_tmp[gfs1.nlevels] = nwp->map[n].t_2m;
        rh_tmp[gfs1.nlevels] = nwp->map[n].rh_2m;
        o3_tmp[gfs1.nlevels] = 0.0;
        if (rh_tmp[ilev] > 0.0) w_tmp[gfs1.nlevels] = rh_to_mr_wat(rh_tmp[gfs1.nlevels], p_tmp[gfs1.nlevels], t_tmp[gfs1.nlevels]);

        if (pflg) {
            printf("ilev=%d, p=%f, t=%f, rh=%f, w=%f\n",ilev,p_tmp[ilev],t_tmp[ilev],rh_tmp[ilev],w_tmp[ilev]);
        }

        /*---------------------------------------------------------------------------
        Extrapolate NWP profile from top level of RH profile to last level of the
        pressure/temperature profile.
        ---------------------------------------------------------------------------*/

        wo = max(w_tmp[kfirst_rh],0.003);
        for (ilev=kfirst_rh; ilev >= 0; ilev--) {
            w_tmp[ilev] = max(wo*pow((p_tmp[ilev]/p_tmp[kfirst_rh]),3),0.003);
        }

        /*---------------------------------------------------------------------------
        Interpolate profile to 101 levels.
        ---------------------------------------------------------------------------*/

        if (profile_to_101(p_tmp,t_tmp,w_tmp,gfs1.nlevels,(double)nwp->map[n].lat,p_prof,t_prof,w_prof)) {
            fprintf(stderr,"%sError interpolating profile in the vertical - aborting\n",EXE_PROMPT);
            exit(EXIT_FAILURE);
        }

        /*---------------------------------------------------------------------------
        Determine the lowest valid level of the new higher resolution profile.
        ---------------------------------------------------------------------------*/

        for (ilev=nwp->nlevels/2; ilev<nwp->nlevels; ilev++) {
            if (p_prof[ilev] >= nwp->map[n].p_sfc) {
                nwp->map[n].sfc_level = ilev;
                break;
            }
        }

        /*---------------------------------------------------------------------------
        Grab an ozone profile from H. Woolf's climatology fortran routine.
        ---------------------------------------------------------------------------*/

        climoz_101_(&nwp->map[n].lat, &month2, o3_prof);


        /*---------------------------------------------------------------------------
        Calculate the height profile of high resolution profile.
        ---------------------------------------------------------------------------*/

        height_profile(p_prof,t_prof,w_prof,z_prof,nwp->nlevels,nwp->map[n].p_msl);

        /*---------------------------------------------------------------------------
        Now store the new profiles in 2D float pointers for each grid cell.
        The sfc_level value now represents the first invalid lower tropospheric
        level which will become a valid level during RTM processing via interpolation.
        ---------------------------------------------------------------------------*/

        for (ilev=0; ilev<nwp->nlevels; ilev++) {
            nwp->map[n].z_lev[ilev] = (float)z_prof[ilev];
            nwp->map[n].p_lev[ilev] = (float)p_prof[ilev];
            nwp->map[n].t_lev[ilev] = (float)t_prof[ilev];
            nwp->map[n].w_lev[ilev] = (float)w_prof[ilev];
            nwp->map[n].o3_lev[ilev] = (float)o3_prof[ilev];
//            if (pflg) printf("ilev=%d, z=%f, p=%f, t=%f, w=%f, o3=%f\n",
//                                 ilev,nwp->map[n].z_lev[ilev],nwp->map[n].p_lev[ilev],nwp->map[n].t_lev[ilev],nwp->map[n].w_lev[ilev],nwp->map[n].o3_lev[ilev]);
            if (nwp->map[n].t_lev[ilev] < 0.){
            	printf("GDAS 1 Temp profile is < 0 at (cell,level) = (%d,%d) = %f\n",n,ilev,
            								nwp->map[n].t_lev[ilev]);
            }
        }

        if (pflg) {
            printf("sfc_lev = %d, PSFC = %f, MLSP = %f\n",nwp->map[n].sfc_level,nwp->map[n].p_sfc,nwp->map[n].p_msl);
            exit(1);
        }

    }

    /*---------------------------------------------------------------------------
    Destroy temporary pointers.
    ---------------------------------------------------------------------------*/

    free(gfs1.lat);
    free(gfs1.lon);
    free(gfs1.sfc_pres);
    free(gfs1.sfc_rh);
    free(gfs1.sfc_z);
    free(gfs1.sfc_t);
    free(gfs1.sfc_snow);
    free(gfs1.sfc_alb);
    free(gfs1.t_2m);
    free(gfs1.tropo_t);
    free(gfs1.tropo_p);
    free(gfs1.tropo_z);
    free(gfs1.col_o3);
    free(gfs1.col_h20);
    free(gfs1.p_msl);
    destroy_2d_float_ptr(nwp->ncells, gfs1.pressure);
    destroy_2d_float_ptr(nwp->ncells, gfs1.pressure_rh);
    destroy_2d_float_ptr(nwp->ncells, gfs1.pressure_o3);
    destroy_2d_float_ptr(nwp->ncells, gfs1.z);
    destroy_2d_float_ptr(nwp->ncells, gfs1.t);
    destroy_2d_float_ptr(nwp->ncells, gfs1.rh);
    destroy_2d_float_ptr(nwp->ncells, gfs1.o3);

    free(gfs2.lat);
    free(gfs2.lon);
    free(gfs2.sfc_pres);
    free(gfs2.sfc_rh);
    free(gfs2.sfc_z);
    free(gfs2.sfc_t);
    free(gfs2.sfc_snow);
    free(gfs2.sfc_alb);
    free(gfs2.t_2m);
    free(gfs2.tropo_t);
    free(gfs2.tropo_p);
    free(gfs2.tropo_z);
    free(gfs2.col_o3);
    free(gfs2.col_h20);
    free(gfs2.p_msl);
    destroy_2d_float_ptr(nwp->ncells, gfs2.pressure);
    destroy_2d_float_ptr(nwp->ncells, gfs2.pressure_rh);
    destroy_2d_float_ptr(nwp->ncells, gfs2.pressure_o3);
    destroy_2d_float_ptr(nwp->ncells, gfs2.z);
    destroy_2d_float_ptr(nwp->ncells, gfs2.t);
    destroy_2d_float_ptr(nwp->ncells, gfs2.rh);
    destroy_2d_float_ptr(nwp->ncells, gfs2.o3);

    free(p_tmp);
    free(t_tmp);
    free(rh_tmp);
    free(w_tmp);
    free(o3_tmp);

    free(z_prof);
    free(p_prof);
    free(t_prof);
    free(w_prof);
    free(o3_prof);

}

extern char gdas_file[];
extern char gdas_file1[];
extern char gdas_file2[];

void adjo3(double *, double *, float *, float *);
void main_gdas(char *gfs_dir, int year, int jday, float time, float *nwp_bounds,
               nwp_params *nwp)
{

    char r_name[256], filename1[256], filename2[256], pflg;
    int year1, year2, jday1, jday2, hour1, hour2, month1, month2, day1, day2;
    int ileap1, ileap2;
    int n=0, ilev=0, kfirst_rh=0, klast_rh=0, kfirst_o3=0, klast_o3=0;
    const float forecast_int = 6.0;
    float x, temp;
    double *p_tmp, *t_tmp, *rh_tmp, *w_tmp, *o3_tmp, *z_prof, *p_prof, *t_prof, *w_prof, wo;
    float *o3_prof, *adjo3_prof;
    nwp_params_tmp gfs1, gfs2;

    strcpy(r_name,"main_gdas");

//  printf("Inside main_gdas()\n"); // GPC

    /*---------------------------------------------------------------------------
    Determine the needed hour and julian day for each gfs file, keeping in mind
    that we are using 0 hour forecasts.
    ---------------------------------------------------------------------------*/
  
    printf("year=%d, jday=%d, time=%f\n",year,jday,time);
    printf("bounds = %f %f %f %f\n", *(nwp_bounds+0), *(nwp_bounds+1), *(nwp_bounds+2), *(nwp_bounds+3));
    get_nwp_time(year, jday, time, 0, &year1, &year2, &ileap1,
                 &ileap2, &jday1, &jday2, &hour1, &hour2);
  

    /*---------------------------------------------------------------------------
    Determine the correct year, month, and day for each gfs file.
    ---------------------------------------------------------------------------*/
  
    month1 = jday2month(jday1, ileap1);
    month2 = jday2month(jday2, ileap2);
    day1 = jday2day(jday1, ileap1);
    day2 = jday2day(jday2, ileap2);
/*
    printf("First Forecast Start Time\n");
    printf("year=%d, month=%d, day=%d, hour=%d\n",year1,month1,day1,hour1);
    printf("Second Forecast Start Time\n");
    printf("year=%d, month=%d, day=%d, hour=%d\n",year2,month2,day2,hour2);
*/
    /*---------------------------------------------------------------------------
    Construct the gfs file names.
    ---------------------------------------------------------------------------*/


    if ( gdas_file[0] == 0 )
    {

//    sprintf(filename1,"%s/gdas1.PGrbF00.%2.2d%2.2d%2.2d.%2.2dz",gfs_dir,year1 - 100*(year1/100),month1,day1,hour1);
//    sprintf(filename2,"%s/gdas1.PGrbF00.%2.2d%2.2d%2.2d.%2.2dz",gfs_dir,year2 - 100*(year2/100),month2,day2,hour2);

      sprintf(filename1,"%s%s",gfs_dir,gdas_file1);
      sprintf(filename2,"%s%s",gfs_dir,gdas_file2); 

    }

    else

    {

       sprintf(filename1,"%s%s", gfs_dir, gdas_file);
       sprintf(filename2,"%s%s", gfs_dir, gdas_file);

    }

    printf("GDAS file #1 = %s\n",filename1);
    printf("GDAS file #2 = %s\n",filename2);

    /*---------------------------------------------------------------------------
    Read in gfs from grib files for both times needed to interpolate to satellite
    overpass time.
    ---------------------------------------------------------------------------*/

    decode_gfs_grib(filename1, &gfs1, nwp_bounds);
    decode_gfs_grib(filename2, &gfs2, nwp_bounds);

    printf("Finished decoding GRIB files...\n\n");

    /*---------------------------------------------------------------------------
    Check to make sure that NWP at time1 is consistent with NWP at time 2.
    ---------------------------------------------------------------------------*/

    //printf("gfs1.nlevels = %i\n",gfs1.nlevels); //RMC 
    //printf("gfs2.nlevels = %i\n",gfs2.nlevels); //RMC

    //printf("gfs1.nlevels_rh = %i\n",gfs1.nlevels_rh); //RMC 
    //printf("gfs2.nlevels_rh = %i\n",gfs2.nlevels_rh); //RMC

    if (gfs1.nlon != gfs2.nlon || gfs1.nlat != gfs2.nlat ||
            gfs1.dlat != gfs2.dlat || gfs1.dlon != gfs2.dlon ||
            gfs1.first_lat != gfs2.first_lat || gfs1.last_lat != gfs2.last_lat ||
            gfs1.first_lon != gfs2.first_lon || gfs1.last_lon != gfs2.last_lon ||
            gfs1.nlevels != gfs2.nlevels || gfs1.nlevels_rh != gfs2.nlevels_rh ||
            gfs1.nlevels_o3 != gfs2.nlevels_o3) {
        fprintf(stderr,"GFS forecast grids are different\n");
        fprintf(stderr,"%d  %d\n",gfs1.nlevels,gfs2.nlevels);
        fprintf(stderr,"%d  %d\n",gfs1.nlevels_rh,gfs2.nlevels_rh);
        fprintf(stderr,"%d  %d\n",gfs1.nlevels_o3,gfs2.nlevels_o3);
        exit(1);
    }

//  printf("%s\n%f  %f  %f\n%f  %f  %f\n",gfs1.grid_name,gfs1.first_lat,gfs1.last_lat,gfs1.dlat,
//         gfs1.first_lon,gfs1.last_lon,gfs1.dlon);

    /*---------------------------------------------------------------------------
    Copy geolocation and dimension data into main NWP structure.
    ---------------------------------------------------------------------------*/

    strcpy(nwp->grid_name,gfs1.grid_name);
    nwp->nlon = gfs1.nlon;
    nwp->nlat = gfs1.nlat;
    nwp->first_lat = gfs1.first_lat;
    nwp->last_lat = gfs1.last_lat;
    nwp->first_lon = gfs1.first_lon;
    nwp->last_lon = gfs1.last_lon;
    nwp->dlat = gfs1.dlat;
    nwp->dlon = gfs1.dlon;

    nwp->ncells = (long)gfs1.nlat*(long)gfs1.nlon;

//  printf("gfs1.nlat: %d\n", gfs1.nlat);
//  printf("gfs1.nlon: %d\n", gfs1.nlon);
//  printf("nwp->ncells: %ld\n", nwp->ncells);

    nwp->nlevels = MAX_LAYERS_CRTM + 1;
    nwp->nlayers = MAX_LAYERS_CRTM;

//  printf("Number of CRTM layers: %d\n",MAX_LAYERS_CRTM); // GPC

    /*---------------------------------------------------------------------------
    Allocate memory for temporary low vertical resolution profiles.
    ---------------------------------------------------------------------------*/

    if ((p_tmp = (double *) calloc(gfs1.nlevels+1,sizeof(double))) == NULL)
        error_allo(r_name,"p_tmp");
    if ((t_tmp = (double *) calloc(gfs1.nlevels+1,sizeof(double))) == NULL)
        error_allo(r_name,"t_tmp");
    if ((rh_tmp = (double *) calloc(gfs1.nlevels+1,sizeof(double))) == NULL)
        error_allo(r_name,"rh_tmp");
    if ((w_tmp = (double *) calloc(gfs1.nlevels+1,sizeof(double))) == NULL)
        error_allo(r_name,"w_tmp");
    if ((o3_tmp = (double *) calloc(gfs1.nlevels+1,sizeof(double))) == NULL)
        error_allo(r_name,"o3_tmp");

    /*---------------------------------------------------------------------------
    Allocate memory for temporary high vertical resolution profiles.
    ---------------------------------------------------------------------------*/
//  printf("Number of levels in z_prof: %d\n",nwp->nlevels);	 // GPC

    if ((z_prof = (double *) calloc(nwp->nlevels, sizeof(double))) == NULL)
        error_allo(r_name,"z_prof");
    if ((p_prof = (double *) calloc(nwp->nlevels, sizeof(double))) == NULL)
        error_allo(r_name,"p_prof");
    if ((t_prof = (double *) calloc(nwp->nlevels, sizeof(double))) == NULL)
        error_allo(r_name,"t_prof");
    if ((w_prof = (double *) calloc(nwp->nlevels, sizeof(double))) == NULL)
        error_allo(r_name,"w_prof");
    if ((o3_prof = (float *) calloc(nwp->nlevels, sizeof(float))) == NULL)
        error_allo(r_name,"o3_prof");
    if ((adjo3_prof = (float *) calloc(nwp->nlevels, sizeof(float))) == NULL)
        error_allo(r_name,"adjo3_prof");

    /*---------------------------------------------------------------------------
    Allocate memory for nwp grid cell data structure.
    ---------------------------------------------------------------------------*/

    if ((nwp->map = (profile_params *) malloc(sizeof(profile_params) * nwp->ncells)) == NULL)
        error_allo(r_name,"nwp->map*");

//  printf("nwp->nlevels inside main_gdas(): %d\n",nwp->nlevels);	//GPC

    for (n=0; n<nwp->ncells; n++) {
        if ((nwp->map[n].z_lev = (float *) malloc(sizeof(float) * nwp->nlevels)) == NULL)
            error_allo(r_name,"nwp->map.z_lev");
        if ((nwp->map[n].p_lev = (float *) malloc(sizeof(float) * nwp->nlevels)) == NULL)
            error_allo(r_name,"nwp->map.p_lev");
        if ((nwp->map[n].t_lev = (float *) malloc(sizeof(float) * nwp->nlevels)) == NULL)
            error_allo(r_name,"nwp->map.t_lev");
        if ((nwp->map[n].w_lev = (float *) malloc(sizeof(float) * nwp->nlevels)) == NULL)
            error_allo(r_name,"nwp->map.w_lev");
        if ((nwp->map[n].o3_lev = (float *) malloc(sizeof(float) * nwp->nlevels)) == NULL)
            error_allo(r_name,"nwp->map.o3_lev");
        if ((nwp->map[n].adjo3_lev = (float *) malloc(sizeof(float) * nwp->nlevels)) == NULL)
            error_allo(r_name,"nwp->map.adjo3_lev");
    }

    /*---------------------------------------------------------------------------
    Determine the time interpolation weight.
    ---------------------------------------------------------------------------*/

    temp = ((int)((time+forecast_int)/forecast_int)*forecast_int);
    x = (temp < (time+forecast_int)) ? (temp - time)/forecast_int : ((temp-forecast_int)-time)/forecast_int;

//  printf("Time interpolation weight x = %f\n",x);

    /*---------------------------------------------------------------------------
    Compute weighted average of 2 GFS times for surface and column parameters.
    ---------------------------------------------------------------------------*/

    for (n=0; n<nwp->ncells; n++) {

        nwp->map[n].p_sfc = (gfs1.sfc_pres[n] >= 0.0 && gfs2.sfc_pres[n] >= 0.0) ?
                            (x*(gfs1.sfc_pres[n]/100.0) + (1.0-x)*(gfs2.sfc_pres[n]/100.0)) : MISSING_FLOAT;
        nwp->map[n].p_msl = (gfs1.p_msl[n] >= 0.0 && gfs2.p_msl[n] >= 0.0) ?
                            (x*(gfs1.p_msl[n]/100.0) + (1.0-x)*(gfs2.p_msl[n]/100.0)) : MISSING_FLOAT;
        nwp->map[n].rh_2m = (gfs1.sfc_rh[n] >= 0.0 && gfs2.sfc_rh[n] >= 0.0) ?
                            (x*gfs1.sfc_rh[n] + (1.0-x)*gfs2.sfc_rh[n]) : MISSING_FLOAT;
        nwp->map[n].t_2m = (gfs1.t_2m[n] >= 0.0 && gfs2.t_2m[n] >= 0.0) ?
                           (x*gfs1.t_2m[n] + (1.0-x)*gfs2.t_2m[n]) : MISSING_FLOAT;
        nwp->map[n].z_sfc = (gfs1.sfc_z[n] != MISSING_GFS || gfs2.sfc_z[n] != MISSING_GFS) ?
                            (x*gfs1.sfc_z[n] + (1.0-x)*gfs2.sfc_z[n]) : MISSING_FLOAT;
        nwp->map[n].t_sfc = (gfs1.sfc_t[n] >= 0.0 && gfs2.sfc_t[n] >= 0.0) ?
                            (x*gfs1.sfc_t[n] + (1.0-x)*gfs2.sfc_t[n]) : MISSING_FLOAT;
        nwp->map[n].snow_sfc = (gfs1.sfc_snow[n] >= 0.0 && gfs2.sfc_snow[n] >= 0.0) ?
                               (x*gfs1.sfc_snow[n] + (1.0-x)*gfs2.sfc_snow[n]) : MISSING_FLOAT;
        nwp->map[n].alb_sfc = (gfs1.sfc_alb[n] >= 0.0 && gfs2.sfc_alb[n] >= 0.0) ?
                              (x*gfs1.sfc_alb[n] + (1.0-x)*gfs2.sfc_alb[n]) : MISSING_FLOAT;

        nwp->map[n].t_tropo = (gfs1.tropo_t[n] >= 0.0 && gfs2.tropo_t[n] >= 0.0) ?
                              (x*gfs1.tropo_t[n] + (1.0-x)*gfs2.tropo_t[n]) : MISSING_FLOAT;
        nwp->map[n].p_tropo = (gfs1.tropo_p[n] >= 0.0 && gfs2.tropo_p[n] >= 0.0) ?
                              (x*gfs1.tropo_p[n] + (1.0-x)*gfs2.tropo_p[n]) : MISSING_FLOAT;
        nwp->map[n].z_tropo = (gfs1.tropo_z[n] >= 0.0 && gfs2.tropo_z[n] >= 0.0) ?
                              (x*gfs1.tropo_z[n] + (1.0-x)*gfs2.tropo_z[n]) : MISSING_FLOAT;

        nwp->map[n].o3_col = (gfs1.col_o3[n] >= 0.0 && gfs2.col_o3[n] >= 0.0) ?
                             (x*gfs1.col_o3[n] + (1.0-x)*gfs2.col_o3[n]) : MISSING_FLOAT;
        nwp->map[n].tpw = (gfs1.col_h20[n] >= 0.0 && gfs2.col_h20[n] >= 0.0) ?
                          (x*gfs1.col_h20[n] + (1.0-x)*gfs2.col_h20[n]) : MISSING_FLOAT;

        nwp->map[n].u_sfc = (abs(gfs1.sfc_u[n]) < 200.0 && abs(gfs2.sfc_u[n]) < 200.0) ?
                            (x*gfs1.sfc_u[n] + (1.0-x)*gfs2.sfc_u[n]) : MISSING_FLOAT;
        nwp->map[n].v_sfc = (abs(gfs1.sfc_v[n]) < 200.0 && abs(gfs2.sfc_v[n]) < 200.0) ?
                            (x*gfs1.sfc_v[n] + (1.0-x)*gfs2.sfc_v[n]) : MISSING_FLOAT;

        /*---------------------------------------------------------------------------
          Check to see if the surface pressure is greater than the lowest level
          atmospheric pressure.  If so, additional level is surface; otherwise
          bottom new level is set to MISSING.
          ---------------------------------------------------------------------------*/
        if (nwp->map[n].p_sfc > 0.0 && gfs1.pressure[n][gfs1.nlevels-1] <= nwp->map[n].p_sfc) {
            nwp->map[n].sfc_level = gfs1.nlevels;
        } else {
            nwp->map[n].sfc_level = 0;
        }

        nwp->map[n].lat = gfs1.lat[n];
        nwp->map[n].lon = gfs1.lon[n];

        nwp->map[n].land = gfs1.land[n];    //RAF added NWP land/sea tag

//      printf("nwp map: %d, %f %f %f\n", n, nwp->map[n].lat, nwp->map[n].lon, nwp->dlat);

    }

    /*---------------------------------------------------------------------------
    Compute weighted average of 2 GFS times for full profile parameters.
    ---------------------------------------------------------------------------*/

    /*---------------------------------------------------------------------------
    First determine the vertical indices of the RH and O3 parameters with respect
    to the full profile resolution.
    ---------------------------------------------------------------------------*/

    for (ilev=0; ilev<gfs1.nlevels; ilev++) {
        if (gfs1.pressure_rh[0][0] == gfs1.pressure[0][ilev]) kfirst_rh = ilev;
        if (gfs1.pressure_o3[0][0] == gfs1.pressure[0][ilev]) kfirst_o3 = ilev;
        if (gfs1.pressure_rh[0][gfs1.nlevels_rh-1] == gfs1.pressure[0][ilev]) klast_rh = ilev;
        if (gfs1.pressure_o3[0][gfs1.nlevels_o3-1] == gfs1.pressure[0][ilev]) klast_o3 = ilev;
    }

    if (gfs1.nlevels_rh == 25) kfirst_rh = 1; // RMC - is there a better way than hardcoding this in for GFS files?
                                              // GDAS files have nlevels_rh=21 while GFS files have nlevels_rh=25

    /*---------------------------------------------------------------------------
    Create the 101 level pressure profile.
    ---------------------------------------------------------------------------*/

    make_profile_101(p_prof);

    /*---------------------------------------------------------------------------
    Start loop over each NWP grid cell.
    ---------------------------------------------------------------------------*/

    pflg = 0;
//  printf("GDAS 1: Looping over %ld NWP grid cells...\n",nwp->ncells);
    for (n=0; n<nwp->ncells; n++) {
        /*if (nwp->sfc_level[n] == 0) pflg = 1;*/
        /*if (n == 620) pflg = 1;*/

        /*---------------------------------------------------------------------------
        Loop over original number of gfs levels.
        ---------------------------------------------------------------------------*/

        for (ilev=0; ilev<gfs1.nlevels; ilev++) {

            /*---------------------------------------------------------------------------
            Pressure - atmospheric levels are fixed.
            ---------------------------------------------------------------------------*/

            p_tmp[ilev] = gfs1.pressure[n][ilev];

            /*---------------------------------------------------------------------------
            Temperature - full vertical resolution.
            ---------------------------------------------------------------------------*/

            t_tmp[ilev] = (gfs1.t[n][ilev] >= 0.0 && gfs2.t[n][ilev] >= 0.0) ?
                          (x*gfs1.t[n][ilev] + (1.0-x)*gfs2.t[n][ilev]) : MISSING_FLOAT;

            if (gfs1.t[n][ilev] < 0. || gfs2.t[n][ilev] < 0.){	// GPC
            	printf("GDAS (gfs1,gfs2,t_tmp) Temp profile is < 0 at (cell,level) = (%d,%d) = (%f , %f, %f)\n",n,ilev,
            			gfs1.t[n][ilev],gfs2.t[n][ilev],t_tmp[ilev]);
            }

            if (gfs1.pressure[n][ilev] != gfs2.pressure[n][ilev]) {
                fprintf(stderr,"Pressure levels in GFS files differ\n");
                exit(1);
            }

            /*---------------------------------------------------------------------------
            Relative humidity - partial vertical resolution.
            ---------------------------------------------------------------------------*/

            rh_tmp[ilev] = 0.0;
            w_tmp[ilev] = 0.0;

            //printf("gfs1.pressure_rh[0][0] = %f gfs1.pressure[0][ilev] = %f \n",gfs1.pressure_rh[0][0],gfs1.pressure[0][ilev]); //RMC
            //printf("gfs1.pressure_rh[n][ilev-kfirst_rh] = %f n = %i ilev = %i kfirst_rh = %i\n ",gfs1.pressure_rh[n][ilev-kfirst_rh],n,ilev,kfirst_rh); //RMC 
            //printf("gfs2.pressure_rh[n][ilev-kfirst_rh] = %f\n",gfs2.pressure_rh[n][ilev-kfirst_rh]); //RMC

            if (ilev >= kfirst_rh && ilev <= klast_rh) {
                rh_tmp[ilev] = (gfs1.rh[n][ilev-kfirst_rh] >= 0.0 && gfs2.rh[n][ilev-kfirst_rh] >= 0.0) ?
                               (x*gfs1.rh[n][ilev-kfirst_rh] + (1.0-x)*gfs2.rh[n][ilev-kfirst_rh]) : MISSING_FLOAT;
                if (gfs1.pressure_rh[n][ilev-kfirst_rh] != gfs2.pressure_rh[n][ilev-kfirst_rh]) {
                    fprintf(stderr,"RH pressure levels in GFS files differ\n");
                    exit(1);
                }
                if (rh_tmp[ilev] > 0.0) w_tmp[ilev] = rh_to_mr_wat(rh_tmp[ilev], p_tmp[ilev], t_tmp[ilev]);
            }

            /*---------------------------------------------------------------------------
            Ozone - partial vertical resolution.
            ---------------------------------------------------------------------------*/

            o3_tmp[ilev] = 0.0;
            if (ilev >= kfirst_o3 && ilev <= klast_o3) {
                o3_tmp[ilev] = (gfs1.o3[n][ilev-kfirst_o3] >= 0.0 && gfs2.o3[n][ilev-kfirst_o3] >= 0.0) ?
                               (x*gfs1.o3[n][ilev-kfirst_o3] + (1.0-x)*gfs2.o3[n][ilev-kfirst_o3]) : MISSING_FLOAT;
                if (gfs1.pressure_o3[n][ilev-kfirst_o3] != gfs2.pressure_o3[n][ilev-kfirst_o3]) {
                    fprintf(stderr,"Ozone pressure levels in GFS files differ\n");
                    exit(1);
                }
            }

            /*---------------------------------------------------------------------------
            Check to see if the surface is encountered before lowest non-surface
            pressure level.
            ---------------------------------------------------------------------------*/

            if (nwp->map[n].p_sfc > 0.0 && p_tmp[ilev] > nwp->map[n].p_sfc) {
                if (nwp->map[n].sfc_level == 0) {
                    nwp->map[n].sfc_level = ilev;
//                       t_tmp[ilev] = nwp->map[n].t_2m;
                    t_tmp[ilev] = nwp->map[n].t_sfc;
                    if( (nwp->map[n].p_sfc - p_tmp[ilev-1]) < 5.0 ||
                        (p_tmp[ilev] - nwp->map[n].p_sfc) < 5.0) {
                      p_tmp[ilev] = (p_tmp[ilev-1] + p_tmp[ilev]) / 2.0;
                    }
                    else {
                      p_tmp[ilev] = nwp->map[n].p_sfc;
                    }
                    o3_tmp[ilev] = 0.0;
                    rh_tmp[ilev] = nwp->map[n].rh_2m;
                    if (rh_tmp[ilev] > 0.0) w_tmp[ilev] = rh_to_mr_wat(rh_tmp[ilev], p_tmp[ilev], t_tmp[ilev]);
                } else {
//                       t_tmp[ilev] = nwp->map[n].t_2m;
                    t_tmp[ilev] = nwp->map[n].t_sfc;
                    p_tmp[ilev] = gfs1.pressure[n][ilev-1];
                    o3_tmp[ilev] = 0.0;
                    rh_tmp[ilev] = nwp->map[n].rh_2m;
                    w_tmp[ilev] = w_tmp[ilev-1];
                }
            }
        }

        /*---------------------------------------------------------------------------
        Add on extra level since we are adding the 2-meter parameters to the profile.
        ---------------------------------------------------------------------------*/

        if (nwp->map[n].sfc_level != gfs1.nlevels) {
            p_tmp[gfs1.nlevels] = gfs1.pressure[n][ilev-1];
        } else {
            p_tmp[gfs1.nlevels] = nwp->map[n].p_sfc;
        }

//      t_tmp[gfs1.nlevels] = nwp->map[n].t_2m;
        t_tmp[gfs1.nlevels] = nwp->map[n].t_sfc;
        rh_tmp[gfs1.nlevels] = nwp->map[n].rh_2m;
        o3_tmp[gfs1.nlevels] = 0.0;
        if (rh_tmp[ilev] > 0.0) w_tmp[gfs1.nlevels] = rh_to_mr_wat(rh_tmp[gfs1.nlevels], p_tmp[gfs1.nlevels], t_tmp[gfs1.nlevels]);

//      printf("t_tmp[gfs1.nlevels] = %f %f %f %f\n",t_tmp[gfs1.nlevels], p_tmp[gfs1.nlevels], nwp->map[n].p_sfc, nwp->map[n].t_sfc);	// GPC

        /*---------------------------------------------------------------------------
        Extrapolate NWP profile from top level of RH profile to last level of the
        pressure/temperature profile.
        ---------------------------------------------------------------------------*/

        wo = max(w_tmp[kfirst_rh],0.003);
        for (ilev=kfirst_rh; ilev >= 0; ilev--) {
            w_tmp[ilev] = max(wo*pow((p_tmp[ilev]/p_tmp[kfirst_rh]),3),0.003);
        }

        /*---------------------------------------------------------------------------
        Interpolate profile to 101 levels.
        ---------------------------------------------------------------------------*/

//      printf("Interpolate profiles from %d levels to 101 levels...\n",gfs1.nlevels);	// GPC

        if (profile_to_101(p_tmp,t_tmp,w_tmp,gfs1.nlevels,(double)nwp->map[n].lat,p_prof,t_prof,w_prof)) {
            fprintf(stderr,"%sError interpolating profile in the vertical - aborting\n",EXE_PROMPT);
            exit(EXIT_FAILURE);
        }

//      printf("gfs1 temperature profile for cell 3694...\n");	// GPC
//      if (n==3694){
        if(pflg) {
          for (ilev=0;ilev<gfs1.nlevels+1;ilev++){
            printf(" %d %f %f %f\n",ilev, p_tmp[ilev], t_tmp[ilev], w_tmp[ilev]);
          }
          for (ilev=0; ilev<101; ilev++) {
            printf("%d %f %f %f %f\n", ilev, (double)nwp->map[n].lat, p_prof[ilev], t_prof[ilev], w_prof[ilev]);
          }
        }

        /*---------------------------------------------------------------------------
        Determine the lowest valid level of the new higher resolution profile.
        ---------------------------------------------------------------------------*/

        for (ilev=nwp->nlevels/2; ilev<nwp->nlevels; ilev++) {
            if (p_prof[ilev] >= nwp->map[n].p_sfc) {
                nwp->map[n].sfc_level = ilev;
                break;
            }
        }

        /*---------------------------------------------------------------------------
        Grab an ozone profile from H. Woolf's climatology fortran routine.
        ---------------------------------------------------------------------------*/

        climoz_101_(&nwp->map[n].lat, &month2, o3_prof);

        /*---------------------------------------------------------------------------
        Adjust climatological ozone profile with stratospheric values from current
        GDAS data.
        ---------------------------------------------------------------------------*/

        adjo3(p_prof, o3_tmp, o3_prof, adjo3_prof);

        /*---------------------------------------------------------------------------
        Calculate the height profile of high resolution profile.
        ---------------------------------------------------------------------------*/

        height_profile(p_prof,t_prof,w_prof,z_prof,nwp->nlevels,nwp->map[n].p_msl);
//      height_profile(p_prof,t_prof,w_prof,z_prof,nwp->nlevels,nwp->map[n].p_sfc);

        /*---------------------------------------------------------------------------
        Now store the new profiles in 2D float pointers for each grid cell.
        The sfc_level value now represents the first invalid lower tropospheric
        level which will become a valid level during RTM processing via interpolation.
        ---------------------------------------------------------------------------*/

        for (ilev=0; ilev<nwp->nlevels; ilev++) {
            nwp->map[n].z_lev[ilev] = (float)z_prof[ilev];
            nwp->map[n].p_lev[ilev] = (float)p_prof[ilev];
            nwp->map[n].t_lev[ilev] = (float)t_prof[ilev];
            nwp->map[n].w_lev[ilev] = (float)w_prof[ilev];
            nwp->map[n].o3_lev[ilev] = (float)o3_prof[ilev];
            nwp->map[n].adjo3_lev[ilev] = (float)adjo3_prof[ilev];
            if (pflg) printf("ilev=%d, z=%f, p=%f, t=%f, w=%f, o3=%f\n",
               ilev,nwp->map[n].z_lev[ilev],nwp->map[n].p_lev[ilev],nwp->map[n].t_lev[ilev],nwp->map[n].w_lev[ilev],nwp->map[n].o3_lev[ilev]);
            if (nwp->map[n].t_lev[ilev] < 0. && ilev <= nwp->map[n].sfc_level) {
            	printf("GDAS Temp profile is < 0 at (cell,level) = (%d,%d) = %f\n",n,ilev,
            								nwp->map[n].t_lev[ilev]);
            	printf("sfc_lev = %d, PSFC = %f, MLSP = %f\n",nwp->map[n].sfc_level,nwp->map[n].p_sfc,nwp->map[n].p_msl);
            	(void)fflush(stdout);
            	exit(1);
            }
        }

      if (pflg) {
        printf("sfc_lev = %d, PSFC = %f, MLSP = %f\n",nwp->map[n].sfc_level,nwp->map[n].p_sfc,nwp->map[n].p_msl);
//        exit(1);
      }

    }

    /*---------------------------------------------------------------------------
    Destroy temporary pointers.
    ---------------------------------------------------------------------------*/

//  printf("\n\nFree the gfs1 pointers... \n"); // GPC
    free(gfs1.lat);
    free(gfs1.lon);
    free(gfs1.land);
    free(gfs1.sfc_pres);
    free(gfs1.sfc_rh);
    free(gfs1.sfc_z);
    free(gfs1.sfc_t);
    free(gfs1.sfc_snow);
    free(gfs1.sfc_alb);
    free(gfs1.t_2m);
    free(gfs1.tropo_t);
    free(gfs1.tropo_p);
    free(gfs1.tropo_z);
    free(gfs1.col_o3);
    free(gfs1.col_h20);
    free(gfs1.p_msl);
    free(gfs1.sfc_u);
    free(gfs1.sfc_v);
    destroy_2d_float_ptr(nwp->ncells, gfs1.pressure);
    destroy_2d_float_ptr(nwp->ncells, gfs1.pressure_rh);
    destroy_2d_float_ptr(nwp->ncells, gfs1.pressure_o3);
    destroy_2d_float_ptr(nwp->ncells, gfs1.z);
    destroy_2d_float_ptr(nwp->ncells, gfs1.t);
    destroy_2d_float_ptr(nwp->ncells, gfs1.rh);
    destroy_2d_float_ptr(nwp->ncells, gfs1.o3);

//  printf("Free the gfs2 pointers... \n"); // GPC
    free(gfs2.lat);
    free(gfs2.lon);
    free(gfs2.land);
    free(gfs2.sfc_pres);
    free(gfs2.sfc_rh);
    free(gfs2.sfc_z);
    free(gfs2.sfc_t);
    free(gfs2.sfc_snow);
    free(gfs2.sfc_alb);
    free(gfs2.t_2m);
    free(gfs2.tropo_t);
    free(gfs2.tropo_p);
    free(gfs2.tropo_z);
    free(gfs2.col_o3);
    free(gfs2.col_h20);
    free(gfs2.p_msl);
    free(gfs2.sfc_u);
    free(gfs2.sfc_v);
    destroy_2d_float_ptr(nwp->ncells, gfs2.pressure);
    destroy_2d_float_ptr(nwp->ncells, gfs2.pressure_rh);
    destroy_2d_float_ptr(nwp->ncells, gfs2.pressure_o3);
    destroy_2d_float_ptr(nwp->ncells, gfs2.z);
    destroy_2d_float_ptr(nwp->ncells, gfs2.t);
    destroy_2d_float_ptr(nwp->ncells, gfs2.rh);
    destroy_2d_float_ptr(nwp->ncells, gfs2.o3);

    free(p_tmp);
    free(t_tmp);
    free(rh_tmp);
    free(w_tmp);
    free(o3_tmp);

    free(z_prof);
    free(p_prof);
    free(t_prof);
    free(w_prof);
    free(o3_prof);

//  printf("Leaving main_gdas()... \n"); // GPC
}

void decode_gfs_grib(char *filename, nwp_params_tmp *gfs, float *gfs_bounds)
{

    char r_name[256];
    unsigned char *buffer, read_flg[400], scale_lon_flg, prime_meridian_flg, dateline_flg;
    float *array, plevels[200], plevels_rh[200], plevels_o3[200];
    float lat_min, lat_max, lon_min, lon_max, f_lat, f_lon, l_lat, l_lon;
    double temp;
    int i, n, nn, nx, ny, ilevel, ncol, nrow, ncol_interp, nrow_interp, ncells_interp;
    int xmin, xmax, ymin, ymax, index;
    long int len_grib, pos = 0, nxny, buffer_size, *isub, n_sub;
    unsigned char *msg, *pds, *gds, *bms, *bds, *pointer;
    int factor = 1;
    FILE *input;

    strcpy(r_name,"decode_gfs_grib");

    gfs->first_lat = -999.0;
    gfs->first_lon = -999.0;
    gfs->last_lat = -999.0;
    gfs->last_lon = -999.0;
    gfs->dlat = -999.0;
    gfs->dlon = -999.0;
    strcpy(gfs->grid_name,"Unknown");

    for (i=0; i<200; i++) {
        plevels[i] = 9999.0;
        plevels_rh[i] = 9999.0;
        plevels_o3[i] = 9999.0;
    }
    for (i=0; i<400; i++) read_flg[i] = 0;

    if ((input = fopen(filename,"rb")) == NULL) {
        fprintf(stderr,"could not open file: %s\n", filename);
        exit(7);
    }

    if ((buffer = (unsigned char *) malloc(BUFF_ALLOC0)) == NULL) {
        fprintf(stderr,"not enough memory\n");
    }
    buffer_size = BUFF_ALLOC0;

    gfs->nlevels = 0;
    gfs->nlevels_rh = 0;
    gfs->nlevels_o3 = 0;
    i = 0;
    for (;;) {

        msg = seek_grib(input, &pos, &len_grib, buffer, MSEEK);
        if (msg == NULL) {
//          printf("Reached end of grib file\n");
            break;
        }

        /* read whole grib record */
        if (len_grib + msg - buffer > buffer_size) {
            buffer_size = len_grib + msg - buffer + 1000;
            buffer = (unsigned char *) realloc((void *) buffer, buffer_size);
            if (buffer == NULL) {
                fprintf(stderr,"ran out of memory\n");
                exit(8);
            }
        }
        read_grib(input, pos, len_grib, buffer);
        pos += len_grib;

        /* parse grib message */
        msg = buffer;
        pds = (msg + 8);
        pointer = pds + PDS_LEN(pds);

        if (PDS_HAS_GDS(pds)) {
            gds = pointer;
            pointer += GDS_LEN(gds);
        } else {
            gds = NULL;
            fprintf(stderr,"There is no GDS in GRIB file.\n");
            exit(8);
        }

        if (PDS_HAS_BMS(pds)) {
            bms = pointer;
            pointer += BMS_LEN(bms);
        } else {
            bms = NULL;
        }

        bds = pointer;
        pointer += BDS_LEN(bds);
        GDS_grid(gds, bds, &nx, &ny, &nxny);

        gfs->nlon = nx;
        gfs->nlat = ny;

        if (nx != gfs->nlon) {
            fprintf(stderr,"x-dimensions of GFS grid are inconsistent: dimx1=%d, dimx2=%d\n",nx,gfs->nlon);
            exit(8);
        }

        if (ny != gfs->nlat) {
            fprintf(stderr,"y-dimensions of GFS grid are inconsistent: dimy1=%d, dimy2=%d\n",ny,gfs->nlat);
            exit(8);
        }

        get_grid_info(gds, nx, &gfs->first_lat, &gfs->last_lat,
                      &gfs->first_lon, &gfs->last_lon, &gfs->dlat, &gfs->dlon, gfs->grid_name);

        if (strcmp(gfs->grid_name,"Equal Angle") != 0) {
            fprintf(stderr,"GFS grid is not Equal Angle, it is %s\n",gfs->grid_name);
            exit(8);
        }

        /*
			For definitions of PDS_PARAM, see http://www.nco.ncep.noaa.gov/pmb/docs/on388/table2.html
			For definitions of PDS_KPDS6, see http://www.nco.ncep.noaa.gov/pmb/docs/on388/table3.html
			For definitions of PDS_Grid,  see http://www.nco.ncep.noaa.gov/pmb/docs/on388/tableb.html
         */

        // Determine the number of height/temp, relative humidity and ozone isobaric levels.

        // If geostrophic hgt and isobaric level, add to layer count.
        if (PDS_PARAM(pds) == 7 && PDS_KPDS6(pds) == 100) {
            plevels[gfs->nlevels] = PDS_KPDS7(pds);
            gfs->nlevels++;
        }

        // If relative humidity and isobaric level, add to RH layer count.
        if (PDS_PARAM(pds) == 52 && PDS_KPDS6(pds) == 100) {
            plevels_rh[gfs->nlevels_rh] = PDS_KPDS7(pds);
            gfs->nlevels_rh++;
        }

        // If ozone mixing ratio and isobaric level, add to ozone layer count.
        if (PDS_PARAM(pds) == 154 && PDS_KPDS6(pds) == 100) {
            plevels_o3[gfs->nlevels_o3] = PDS_KPDS7(pds);
            gfs->nlevels_o3++;
        }

        /*
			These are the grib datasets that are read in. Here we count how
			many there are...
         */

        // RAF:  added 1-10m u,v winds (kpds5=33, 34 and kpds6 = 105)
        // RAF:  added land/sea flag (kpds5=81 and kpds6 = 1)
        if (
                (PDS_PARAM(pds) == 1   && PDS_KPDS6(pds) == 1  ) ||		// Pressure, Ground or water surface
                (PDS_PARAM(pds) == 1   && PDS_KPDS6(pds) == 7  ) ||		// Pressure, Tropopause
                (PDS_PARAM(pds) == 2   && PDS_KPDS6(pds) == 102) ||		// Pressure reduced to MSL, mean sea level
                (PDS_PARAM(pds) == 7   && PDS_KPDS6(pds) == 1  ) ||		// Geopotential height, Ground or water surface
                (PDS_PARAM(pds) == 7   && PDS_KPDS6(pds) == 7  ) ||		// Geopotential height, Tropopause
                (PDS_PARAM(pds) == 7   && PDS_KPDS6(pds) == 100) ||		// Geopotential height, isobaric level
                (PDS_PARAM(pds) == 10  && PDS_KPDS6(pds) == 200) ||		// Total ozone, entire atmosphere (considered as a single layer)
                (PDS_PARAM(pds) == 11  && PDS_KPDS6(pds) == 1  ) ||		// Temperature, Ground or water surface
                (PDS_PARAM(pds) == 11  && PDS_KPDS6(pds) == 7  ) ||		// Temperature, Tropopause
                (PDS_PARAM(pds) == 11  && PDS_KPDS6(pds) == 100) ||		// Temperature, isobaric level
                (PDS_PARAM(pds) == 11  && PDS_KPDS6(pds) == 105) ||		// Temperature, specified height level above ground
                (PDS_PARAM(pds) == 33  && PDS_KPDS6(pds) == 105) ||		// u-component of wind, specified height level above ground
        	(PDS_PARAM(pds) == 34  && PDS_KPDS6(pds) == 105) ||		// v-component of wind, specified height level above ground
                (PDS_PARAM(pds) == 52  && PDS_KPDS6(pds) == 100) ||		// Relative humidity, isobaric level
                (PDS_PARAM(pds) == 52  && PDS_KPDS6(pds) == 105) ||		// Relative humidity, specified height level above ground
                (PDS_PARAM(pds) == 54  && PDS_KPDS6(pds) == 200) ||		// Precipitable water, entire atmosphere (considered as a single layer)
                (PDS_PARAM(pds) == 65  && PDS_KPDS6(pds) == 1  ) ||		// Water equiv. of accum. snow depth, Ground or water surface
                (PDS_PARAM(pds) == 71  && PDS_KPDS6(pds) == 7  ) ||		// Total cloud cover, Tropopause
                (PDS_PARAM(pds) == 81  && PDS_KPDS6(pds) == 1  ) ||             // Land/sea flag (1=land)
                (PDS_PARAM(pds) == 84  && PDS_KPDS6(pds) == 1  ) ||		// Albedo, Ground or water surface
                (PDS_PARAM(pds) == 154 && PDS_KPDS6(pds) == 100) 		// Ozone mixing ratio, isobaric level
            ) {
            if (PDS_Grid(pds) != 3) {
            	// Grid type 3: 65160-point (360x181) global longitude-latitude grid. (1,1) at (0E, 90N),
            	// matrix layout. N.B.: prime meridian not duplicated
                fprintf(stderr,"GFS data are not lat/lon gridded.\n");
                exit(8);
            }
            read_flg[i] = 1;
        }
        i++;
    }

    if ((isub = (long *) malloc(sizeof(long) * nxny)) == NULL)
        error_allo(r_name,"sub_index");

    /*---------------------------------------------------------------------------
    The longitude bounds are expected to be in -180 - 180 format unless near the
    dateline exluding the poles.
    ---------------------------------------------------------------------------*/

    lat_min = gfs_bounds[0];
    lat_max = gfs_bounds[2];
    lon_min = gfs_bounds[1];
    lon_max = gfs_bounds[3];

    /*---------------------------------------------------------------------------
    Determine if the lon bounds are in 0-360 or -180-180 format.
    ---------------------------------------------------------------------------*/

    scale_lon_flg = NO;
    if (lon_max > 180.0 || lon_min > 180.0) scale_lon_flg = YES;

    if (lat_min >= lat_max) {
        fprintf(stderr,"Subset lat min (%f) and max (%f) are invalid.\n",lat_min,lat_max);
        exit(EXIT_FAILURE);
    }

    if (lon_min >= lon_max) {
        fprintf(stderr,"Subset lon min (%f) and max (%f) are invalid.\n",lon_min,lon_max);
        exit(EXIT_FAILURE);
    }

    /*---------------------------------------------------------------------------
    If near the poles, assume that the area crosses over the pole.
    ---------------------------------------------------------------------------*/

    if (lat_min <= -75.0) {
        lon_min = 0.0;
        lon_max = 360.0;
        lat_min = -90.0;
        lat_max = lat_max;
        scale_lon_flg = YES;
    }
    if (lat_max >= 75.0) {
        lon_min = 0.0;
        lon_max = 360.0;
        lat_min = lat_min;
        lat_max = 90.0;
        scale_lon_flg = YES;
    }

    /*---------------------------------------------------------------------------
    Extend the subsampled area by one dlat and dlon to be sure the area is
    covered well and to be sure that any spatial interpolation does not require
    nearest neighbor sampling.
    ---------------------------------------------------------------------------*/

    lat_min = max(lat_min-gfs->dlat,-90.0);
    lat_max = min(lat_max+gfs->dlat,90.0);

    dateline_flg = NO;
    if (lon_min < 0.0 && lon_max > 0.0 && lon_min <= -90.0 && lon_max >= 90.0) dateline_flg = YES;

    if (scale_lon_flg == NO && dateline_flg == NO) {
        lon_min = max(lon_min-gfs->dlon,-179.9);
        lon_max = min(lon_max+gfs->dlon,179.9);
    } else if(dateline_flg == YES && scale_lon_flg == NO) {
    	lon_min = max(lon_min+gfs->dlon,-179.9);
 	    lon_max = min(lon_max-gfs->dlon,179.9);
    } else {
        lon_min = max(lon_min-gfs->dlon,0.0);
        lon_max = min(lon_max+gfs->dlon,359.9);
    }

    /*---------------------------------------------------------------------------
    Check to see if the prime meridian is contained within the bounds.
    ---------------------------------------------------------------------------*/

//  printf("lon min, max: %f %f %f %f \n", lon_min, lon_max, gfs->first_lon, gfs->last_lon);

    prime_meridian_flg = NO;
    if (lon_min < 0.0 && lon_max > 0.0 && lon_min > -90.0 && lon_max < 90.0) prime_meridian_flg = YES;

    /*---------------------------------------------------------------------------
    The original NWP starting and ending lons are always scaled from 0-360.
    ---------------------------------------------------------------------------*/

    if (gfs->first_lon < 0.0) gfs->first_lon += 360.0;
    if (gfs->last_lon < 0.0) gfs->last_lon += 360.0;

    /*---------------------------------------------------------------------------
    Find the index values of the subset boundaries in NWP space.
    ---------------------------------------------------------------------------*/

    index = eqangle_index(lat_min, convert_lon(0,lon_min), gfs->dlat, gfs->dlon,
                          gfs->first_lat, gfs->first_lon, gfs->last_lat, gfs->last_lon, &ymin, &xmin);

    index = eqangle_index(lat_max, convert_lon(0,lon_max), gfs->dlat, gfs->dlon,
                          gfs->first_lat, gfs->first_lon, gfs->last_lat, gfs->last_lon, &ymax, &xmax);


    if (ymin > ymax) swap_int_values(&ymin, &ymax);
    if (xmin == xmax) {
        xmin = 0;
        xmax = nx-1;
    }
    if(xmax == 0) prime_meridian_flg = YES;

//  printf("xmin %d ,xmax %d ,ymin %d ,ymax %d %f %f %f %f\n", xmin, xmax, ymin, ymax, gfs->first_lon, gfs->last_lon,
//  		gfs->first_lat, gfs->last_lat);

    fprintf(stdout,"%sNWP: {minlat=%f,  maxlat=%f,  minlon=%f,  maxlon=%f}\n",
            EXE_PROMPT,lat_min,lat_max,lon_min,lon_max);
    fprintf(stdout,"%sscale_lon_flg = %d, prime_meridian_flg = %d, dateline_flg = %d\n",EXE_PROMPT,
    		       scale_lon_flg,prime_meridian_flg,dateline_flg);

    /*---------------------------------------------------------------------------
    Determine the 1D index array that describes the subsetted area.
    ---------------------------------------------------------------------------*/

    if (prime_meridian_flg == YES) {
        swap_int_values(&xmin, &xmax);
        i = 0;
        for (nn=ymin; nn<=ymax; nn++) {
            for (n=xmax; n<nx; n++) {
                isub[i] = nn*nx + n;
                i++;
            }
            for (n=0; n<=xmin; n++) {
                isub[i] = nn*nx + n;
                i++;
            }
        }
        nrow = (ymax - ymin) + 1;
        ncol = (xmin + 1) + (nx - xmax);
        eqangle_index2latlon(gfs->first_lon, gfs->dlon, gfs->first_lat, gfs->dlat,
                             xmax, ymin, &f_lon, &f_lat);
        eqangle_index2latlon(gfs->first_lon, gfs->dlon, gfs->first_lat, gfs->dlat,
                             xmin, ymax, &l_lon, &l_lat);
    } else {
    	// RAF: Account for new longitude processing.
    	if(dateline_flg == YES) swap_int_values(&xmin, &xmax);
        i = 0;
        for (nn=ymin; nn<=ymax; nn++) {
            for (n=xmin; n<=xmax; n++) {
                isub[i] = nn*nx + n;
                i++;
            }
        }
        nrow = (ymax - ymin) + 1;
        ncol = (xmax - xmin) + 1;
        eqangle_index2latlon(gfs->first_lon, gfs->dlon, gfs->first_lat, gfs->dlat,
                             xmin, ymin, &f_lon, &f_lat);
        eqangle_index2latlon(gfs->first_lon, gfs->dlon, gfs->first_lat, gfs->dlat,
                             xmax, ymax, &l_lon, &l_lat);
    }

    n_sub = i;

    /*---------------------------------------------------------------------------
    Calculate the size of the spatially interpolated and sub-sampled NWP data.
    ---------------------------------------------------------------------------*/

    ncol = n_sub/nrow;

    ncol_interp = factor*ncol;
    nrow_interp = factor*nrow;
    ncells_interp = ncol_interp*nrow_interp;
    printf("ncol=%d, nrow=%d, n_sub=%ld, nx=%d, factor=%d\n",ncol,nrow,n_sub,nx,factor);
    if (ncells_interp != (factor*factor*n_sub)) {
        fprintf(stderr,"%sError in GFS spatial interpolation - inconsistent dimensions\n",EXE_PROMPT);
        exit (EXIT_FAILURE);
    }

    /*---------------------------------------------------------------------------
    Determine the new grid descriptors for the spatially interpolated subsampled
    grid and store in output structure.
    ---------------------------------------------------------------------------*/

    gfs->first_lat = f_lat;
    gfs->last_lat = l_lat;
    gfs->first_lon = f_lon;
    gfs->last_lon = l_lon;
    gfs->nlat = nrow_interp;
    gfs->nlon = ncol_interp;
//  printf("f_lat section: %f %f %f %f\n", f_lat, l_lat, f_lon, l_lon);

    /*---------------------------------------------------------------------------
    Calculate the spatially interpolated lat/lon of the subsetted GFS grid.
    ---------------------------------------------------------------------------*/

    if ((gfs->lat = (float *) malloc(sizeof(float) * ncells_interp)) == NULL)
        error_allo(r_name,"gfs->lat");

    if ((gfs->lon = (float *) malloc(sizeof(float) * ncells_interp)) == NULL)
        error_allo(r_name,"gfs->lon");

//  RAF: Change sign of 'dlat' so that the calculation of lat/lon for each nwp grid cell is correct.
//       Then, change back so that the correct nwp products are stored.  I have no real understanding
//       of this - but it works.
    gfs->dlat = -gfs->dlat;
    i = 0;
    for (nn=0; nn<nrow_interp; nn++) {
        for (n=0; n<ncol_interp; n++) {
            eqangle_index2latlon(gfs->first_lon, gfs->dlon/factor, gfs->first_lat, gfs->dlat/factor,
                                 n, nn, &gfs->lon[i], &gfs->lat[i]);
            i++;
        }
    }
    if (scale_lon_flg == NO && gfs->first_lon > 180.0) gfs->first_lon -= 360.0;
    if (scale_lon_flg == NO && gfs->last_lon > 180.0) gfs->last_lon -= 360.0;
    if (dateline_flg == YES && gfs->last_lon < 0.0) gfs->last_lon += 360.0;          //RAF

//  printf("final gfs: %f %f %f %f\n", gfs->first_lon, gfs->last_lon, gfs->first_lat, gfs->last_lat);
    gfs->dlat = -gfs->dlat;

    qsort(plevels,sizeof(plevels)/sizeof(float),sizeof(float),compare);
    qsort(plevels_rh,sizeof(plevels_rh)/sizeof(float),sizeof(float),compare);
    qsort(plevels_o3,sizeof(plevels_o3)/sizeof(float),sizeof(float),compare);

    if ((array = (float *) malloc(sizeof(float) * nxny)) == NULL)
        error_allo(r_name,"array");

    gfs->z = allocate_2d_float_ptr("gfs->z", ncells_interp, gfs->nlevels);
    gfs->t = allocate_2d_float_ptr("gfs->t", ncells_interp, gfs->nlevels);
    gfs->rh = allocate_2d_float_ptr("gfs->rh", ncells_interp, gfs->nlevels_rh);
    gfs->o3 = allocate_2d_float_ptr("gfs->o3", ncells_interp, gfs->nlevels_o3);
    gfs->pressure = allocate_2d_float_ptr("gfs->pressure", ncells_interp, gfs->nlevels);
    gfs->pressure_rh = allocate_2d_float_ptr("gfs->pressure_rh", ncells_interp, gfs->nlevels_rh);
    gfs->pressure_o3 = allocate_2d_float_ptr("gfs->pressure_o3", ncells_interp, gfs->nlevels_o3);

    if ((gfs->sfc_pres = (float *) malloc(sizeof(float) * ncells_interp)) == NULL)
        error_allo(r_name,"gfs->sfc_pres");
    if ((gfs->p_msl = (float *) malloc(sizeof(float) * ncells_interp)) == NULL)
        error_allo(r_name,"gfs->p_msl");
    if ((gfs->sfc_rh = (float *) malloc(sizeof(float) * ncells_interp)) == NULL)
        error_allo(r_name,"gfs->sfc_rh");
    if ((gfs->sfc_z = (float *) malloc(sizeof(float) * ncells_interp)) == NULL)
        error_allo(r_name,"gfs->sfc_z");
    if ((gfs->sfc_t = (float *) malloc(sizeof(float) * ncells_interp)) == NULL)
        error_allo(r_name,"gfs->sfc_t");
    if ((gfs->sfc_snow = (float *) malloc(sizeof(float) * ncells_interp)) == NULL)
        error_allo(r_name,"gfs->sfc_snow");
    if ((gfs->sfc_alb = (float *) malloc(sizeof(float) * ncells_interp)) == NULL)
        error_allo(r_name,"gfs->sfc_alb");
    if ((gfs->t_2m = (float *) malloc(sizeof(float) * ncells_interp)) == NULL)
        error_allo(r_name,"gfs->t_2m");
    if ((gfs->tropo_t = (float *) malloc(sizeof(float) * ncells_interp)) == NULL)
        error_allo(r_name,"gfs->tropo_t");
    if ((gfs->tropo_p = (float *) malloc(sizeof(float) * ncells_interp)) == NULL)
        error_allo(r_name,"gfs->tropo_p");
    if ((gfs->tropo_z = (float *) malloc(sizeof(float) * ncells_interp)) == NULL)
        error_allo(r_name,"gfs->tropo_z");
    if ((gfs->col_o3 = (float *) malloc(sizeof(float) * ncells_interp)) == NULL)
        error_allo(r_name,"gfs->col_o3");
    if ((gfs->col_h20 = (float *) malloc(sizeof(float) * ncells_interp)) == NULL)
        error_allo(r_name,"gfs->col_h20");
    // Sfc u,v winds
    if ((gfs->sfc_u = (float *) malloc(sizeof(float) * ncells_interp)) == NULL)
            error_allo(r_name,"gfs->sfc_u");
    if ((gfs->sfc_v = (float *) malloc(sizeof(float) * ncells_interp)) == NULL)
            error_allo(r_name,"gfs->sfc_v");
    // Land/sea flag
    if ((gfs->land = (float *) malloc(sizeof(float) * ncells_interp)) == NULL)
            error_allo(r_name,"gfs->land");


    // The pressure levels for the temp (26) and rh (21)
    for (i=0; i<ncells_interp; i++) {
        for (n=0; n<gfs->nlevels; n++) {
            gfs->pressure[i][n] = plevels[n];
            if (n < gfs->nlevels_rh) gfs->pressure_rh[i][n] = plevels_rh[n];
            if (n < gfs->nlevels_o3) gfs->pressure_o3[i][n] = plevels_o3[n];
        }
    }

    pos = 0;
    i = 0;
    for (;;) {

        msg = seek_grib(input, &pos, &len_grib, buffer, MSEEK);
        if (msg == NULL) {
//          printf("Reached end of grib file\n");
            break;
        }

        if (read_flg[i]) {

            /* read whole grib record */
            if (len_grib + msg - buffer > buffer_size) {
                buffer_size = len_grib + msg - buffer + 1000;
                buffer = (unsigned char *) realloc((void *) buffer, buffer_size);
                if (buffer == NULL) {
                    fprintf(stderr,"ran out of memory\n");
                    exit(8);
                }
            }
            read_grib(input, pos, len_grib, buffer);

            /* parse grib message */
            msg = buffer;
            pds = (msg + 8);

            pointer = pds + PDS_LEN(pds);

            if (PDS_HAS_GDS(pds)) {
                gds = pointer;
                pointer += GDS_LEN(gds);
            } else {
                gds = NULL;
            }

            if (PDS_HAS_BMS(pds)) {
                bms = pointer;
                pointer += BMS_LEN(bms);
            } else {
                bms = NULL;
            }

            bds = pointer;
            pointer += BDS_LEN(bds);

            temp = int_power(10.0, - PDS_DecimalScale(pds));

            /* Read in profile variables*/
            if (PDS_KPDS6(pds) == 100) {
                switch (PDS_PARAM(pds)) {
                case 7:
                    ilevel = -999.0;
                    for (n=gfs->nlevels; n>=0; n--) {
                        if (plevels[n] == PDS_KPDS7(pds)) {
                            ilevel = n;
                            break;
                        }
                    }
                    if (ilevel < 0.0 || ilevel >= gfs->nlevels) {
                        fprintf(stderr,"Current level %d does not match any in established profile\n",PDS_KPDS7(pds));
                        exit(1);
                    }
                    BDS_unpack(array, bds, BMS_bitmap(bms), BDS_NumBits(bds), nxny,
                               temp*BDS_RefValue(bds),temp*int_power(2.0, BDS_BinScale(bds)));
                    for (n=0; n<n_sub; n++) gfs->z[n][ilevel] = array[isub[n]];
                    break;
                case 11:
                    ilevel = -999.0;
                    for (n=gfs->nlevels; n>=0; n--) {
                        if (plevels[n] == PDS_KPDS7(pds)) {
                            ilevel = n;
                            break;
                        }
                    }
                    if (ilevel < 0.0 || ilevel >= gfs->nlevels) {
                        fprintf(stderr,"Current level %d does not match any in established profile\n",PDS_KPDS7(pds));
                        exit(1);
                    }
                    BDS_unpack(array, bds, BMS_bitmap(bms), BDS_NumBits(bds), nxny,
                               temp*BDS_RefValue(bds),temp*int_power(2.0, BDS_BinScale(bds)));
                    for (n=0; n<n_sub; n++) gfs->t[n][ilevel] = array[isub[n]];
                    break;
                case 52:
                    ilevel = -999.0;
                    for (n=gfs->nlevels_rh; n>=0; n--) {
                        if (plevels_rh[n] == PDS_KPDS7(pds)) {
                            ilevel = n;
                            break;
                        }
                    }
                    if (ilevel < 0.0 || ilevel >= gfs->nlevels_rh) {
                        fprintf(stderr,"Current level %d does not match any in established profile\n",PDS_KPDS7(pds));
                        exit(1);
                    }
                    BDS_unpack(array, bds, BMS_bitmap(bms), BDS_NumBits(bds), nxny,
                               temp*BDS_RefValue(bds),temp*int_power(2.0, BDS_BinScale(bds)));
                    for (n=0; n<n_sub; n++) gfs->rh[n][ilevel] = array[isub[n]];
                    break;
                case 154:
                    ilevel = -999.0;
                    for (n=gfs->nlevels_o3; n>=0; n--) {
                        if (plevels_o3[n] == PDS_KPDS7(pds)) {
                            ilevel = n;
                            break;
                        }
                    }
                    if (ilevel < 0.0 || ilevel >= gfs->nlevels_o3) {
                        fprintf(stderr,"Current level %d does not match any in established profile\n",PDS_KPDS7(pds));
                        exit(1);
                    }
                    BDS_unpack(array, bds, BMS_bitmap(bms), BDS_NumBits(bds), nxny,
                               temp*BDS_RefValue(bds),temp*int_power(2.0, BDS_BinScale(bds)));
                    for (n=0; n<n_sub; n++) gfs->o3[n][ilevel] = array[isub[n]];
                    break;
                }
            }
            /* Read in surface variables*/
            if (PDS_KPDS6(pds) == 1) {
                switch (PDS_PARAM(pds)) {
                case 1:
                    BDS_unpack(array, bds, BMS_bitmap(bms), BDS_NumBits(bds), nxny,
                               temp*BDS_RefValue(bds),temp*int_power(2.0, BDS_BinScale(bds)));
                    for (n=0; n<n_sub; n++) gfs->sfc_pres[n] = array[isub[n]];
                    break;
                case 7:
                    BDS_unpack(array, bds, BMS_bitmap(bms), BDS_NumBits(bds), nxny,
                               temp*BDS_RefValue(bds),temp*int_power(2.0, BDS_BinScale(bds)));
                    for (n=0; n<n_sub; n++) gfs->sfc_z[n] = array[isub[n]];
                    break;
                case 11:
                    BDS_unpack(array, bds, BMS_bitmap(bms), BDS_NumBits(bds), nxny,
                               temp*BDS_RefValue(bds),temp*int_power(2.0, BDS_BinScale(bds)));
                    for (n=0; n<n_sub; n++) gfs->sfc_t[n] = array[isub[n]];
                    break;
                case 65:
                    BDS_unpack(array, bds, BMS_bitmap(bms), BDS_NumBits(bds), nxny,
                               temp*BDS_RefValue(bds),temp*int_power(2.0, BDS_BinScale(bds)));
                    for (n=0; n<n_sub; n++) gfs->sfc_snow[n] = array[isub[n]];
                    break;
                case 81:
                    BDS_unpack(array, bds, BMS_bitmap(bms), BDS_NumBits(bds), nxny,
                               temp*BDS_RefValue(bds),temp*int_power(2.0, BDS_BinScale(bds)));
                    for (n=0; n<n_sub; n++) gfs->land[n] = array[isub[n]];
                    break;
                case 84:
                    BDS_unpack(array, bds, BMS_bitmap(bms), BDS_NumBits(bds), nxny,
                               temp*BDS_RefValue(bds),temp*int_power(2.0, BDS_BinScale(bds)));
                    for (n=0; n<n_sub; n++) gfs->sfc_alb[n] = array[isub[n]];
                    break;
                }
            }
            /* Read in MSL variables*/
            if (PDS_KPDS6(pds) == 102) {
                switch (PDS_PARAM(pds)) {
                case 2:
                    BDS_unpack(array, bds, BMS_bitmap(bms), BDS_NumBits(bds), nxny,
                               temp*BDS_RefValue(bds),temp*int_power(2.0, BDS_BinScale(bds)));
                    for (n=0; n<n_sub; n++) gfs->p_msl[n] = array[isub[n]];
                    break;
                }
            }
            /* Read in 2 m above ground variables*/
            if (PDS_KPDS6(pds) == 105) {
                switch (PDS_PARAM(pds)) {
                case 11:
                    BDS_unpack(array, bds, BMS_bitmap(bms), BDS_NumBits(bds), nxny,
                               temp*BDS_RefValue(bds),temp*int_power(2.0, BDS_BinScale(bds)));
                    if (PDS_KPDS7(pds) != 2) {
                        printf("WARNING-lowest level air temperature is not 2 m (as expected), it is %d.\n",PDS_KPDS7(pds));
                        printf("Continuing anyway.\n");
                    }
                    for (n=0; n<n_sub; n++) gfs->t_2m[n] = array[isub[n]];
                    break;
                case 52:
                    BDS_unpack(array, bds, BMS_bitmap(bms), BDS_NumBits(bds), nxny,
                               temp*BDS_RefValue(bds),temp*int_power(2.0, BDS_BinScale(bds)));
                    if (PDS_KPDS7(pds) != 2) {
                        printf("WARNING-lowest level RH is not 2 m (as expected), it is %d.\n",PDS_KPDS7(pds));
                        printf("Continuing anyway.\n");
                    }
                    for (n=0; n<n_sub; n++) gfs->sfc_rh[n] = array[isub[n]];
                    break;
                case 33:
                    BDS_unpack(array, bds, BMS_bitmap(bms), BDS_NumBits(bds), nxny,
                               temp*BDS_RefValue(bds),temp*int_power(2.0, BDS_BinScale(bds)));
                    if (PDS_KPDS7(pds) != 10) {
                        printf("WARNING-lowest level u-wind is not 1-10 m (as expected), it is %d.\n",PDS_KPDS7(pds));
                        printf("Continuing anyway.\n");
                    }
                    for (n=0; n<n_sub; n++) gfs->sfc_u[n] = array[isub[n]];
                    break;
                case 34:
                    BDS_unpack(array, bds, BMS_bitmap(bms), BDS_NumBits(bds), nxny,
                               temp*BDS_RefValue(bds),temp*int_power(2.0, BDS_BinScale(bds)));
                    if (PDS_KPDS7(pds) != 10) {
                        printf("WARNING-lowest level v-wind is not 1-10 m (as expected), it is %d.\n",PDS_KPDS7(pds));
                        printf("Continuing anyway.\n");
                    }
                    for (n=0; n<n_sub; n++) gfs->sfc_v[n] = array[isub[n]];
                    break;
                }
            }
            /* Read in tropopause grid variables*/
            if (PDS_KPDS6(pds) == 7) {
                switch (PDS_PARAM(pds)) {
                case 1:
                    BDS_unpack(array, bds, BMS_bitmap(bms), BDS_NumBits(bds), nxny,
                               temp*BDS_RefValue(bds),temp*int_power(2.0, BDS_BinScale(bds)));
                    for (n=0; n<n_sub; n++) gfs->tropo_p[n] = array[isub[n]];
                    break;
                case 7:
                    BDS_unpack(array, bds, BMS_bitmap(bms), BDS_NumBits(bds), nxny,
                               temp*BDS_RefValue(bds),temp*int_power(2.0, BDS_BinScale(bds)));
                    for (n=0; n<n_sub; n++) gfs->tropo_z[n] = array[isub[n]];
                    break;
                case 11:
                    BDS_unpack(array, bds, BMS_bitmap(bms), BDS_NumBits(bds), nxny,
                               temp*BDS_RefValue(bds),temp*int_power(2.0, BDS_BinScale(bds)));
                    for (n=0; n<n_sub; n++) gfs->tropo_t[n] = array[isub[n]];
                    break;
                }
            }
            /* Read in column density variables*/
            if (PDS_KPDS6(pds) == 200) {
                switch (PDS_PARAM(pds)) {
                case 10:
                    BDS_unpack(array, bds, BMS_bitmap(bms), BDS_NumBits(bds), nxny,
                               temp*BDS_RefValue(bds),temp*int_power(2.0, BDS_BinScale(bds)));
                    for (n=0; n<n_sub; n++) gfs->col_o3[n] = array[isub[n]];
                    break;
                case 54:
                    BDS_unpack(array, bds, BMS_bitmap(bms), BDS_NumBits(bds), nxny,
                               temp*BDS_RefValue(bds),temp*int_power(2.0, BDS_BinScale(bds)));
                    for (n=0; n<n_sub; n++) gfs->col_h20[n] = array[isub[n]];
                    break;
                }
            }

        }
        pos += len_grib;
        i++;
    }

    fclose(input);

    free(buffer);
    free(array);
    free(isub);
}

void get_grid_info(unsigned char *gds, int nx,
                   float *first_lat, float *last_lat,
                   float *first_lon, float *last_lon,
                   float *dlat, float *dlon, char *grid_name)

{
    if (gds && GDS_LatLon(gds) && nx != -1) {
        *first_lat = 0.001*GDS_LatLon_La1(gds);
        *first_lon = 0.001*GDS_LatLon_Lo1(gds);
        *last_lat = 0.001*GDS_LatLon_La2(gds);
        *last_lon = 0.001*GDS_LatLon_Lo2(gds);
        *dlat = 0.001*GDS_LatLon_dy(gds);
        *dlon = 0.001*GDS_LatLon_dx(gds);
        strcpy(grid_name,"Equal Angle");
    } else if (gds && GDS_LatLon(gds) && nx == -1) {
        GDS_prt_thin_lon(gds);
        *first_lat = 0.001*GDS_LatLon_La1(gds);
        *first_lon = 0.001*GDS_LatLon_Lo1(gds);
        *last_lat = 0.001*GDS_LatLon_La2(gds);
        *last_lon = 0.001*GDS_LatLon_Lo2(gds);
        *dlat = 0.001*GDS_LatLon_dy(gds);
        *dlon = -999.0;
        strcpy(grid_name,"Thinned Equal Area");
    } else if (gds && GDS_Gaussian(gds) && nx != -1) {
        *first_lat = 0.001*GDS_LatLon_La1(gds);
        *first_lon = 0.001*GDS_LatLon_Lo1(gds);
        *last_lat = 0.001*GDS_LatLon_La2(gds);
        *last_lon = 0.001*GDS_LatLon_Lo2(gds);
        *dlat = -999.0;
        *dlon = 0.001*GDS_LatLon_dx(gds);
        strcpy(grid_name,"Gaussian");
    } else if (gds && GDS_Gaussian(gds) && nx == -1) {
        GDS_prt_thin_lon(gds);
        *first_lat = 0.001*GDS_LatLon_La1(gds);
        *first_lon = 0.001*GDS_LatLon_Lo1(gds);
        *last_lat = 0.001*GDS_LatLon_La2(gds);
        *last_lon = -999.0;
        *dlat = -999.0;
        *dlon = -999.0;
        strcpy(grid_name,"Thinned Gaussian");
    } else if (gds && GDS_Polar(gds)) {
        *first_lat = -999.0;
        *first_lon = -999.0;
        *last_lat = -999.0;
        *last_lon = -999.0;
        *dlat = -999.0;
        *dlon = -999.0;
        strcpy(grid_name,"Polar Stereo Graphic");
    } else if (gds && GDS_Lambert(gds)) {
        *first_lat = -999.0;
        *first_lon = -999.0;
        *last_lat = -999.0;
        *last_lon = -999.0;
        *dlat = -999.0;
        *dlon = -999.0;
        strcpy(grid_name,"Lambert Conformal");
    } else if (gds && GDS_Mercator(gds)) {
        *first_lat = -999.0;
        *first_lon = -999.0;
        *last_lat = -999.0;
        *last_lon = -999.0;
        *dlat = -999.0;
        *dlon = -999.0;
        strcpy(grid_name,"Mercator");
    } else if (gds && GDS_ssEgrid(gds)) {
        *first_lat = -999.0;
        *first_lon = -999.0;
        *last_lat = -999.0;
        *last_lon = -999.0;
        *dlat = -999.0;
        *dlon = -999.0;
        strcpy(grid_name,"Semi-staggered Arakawa E-Grid");
    } else if (gds && GDS_ss2dEgrid(gds)) {
        *first_lat = -999.0;
        *first_lon = -999.0;
        *last_lat = -999.0;
        *last_lon = -999.0;
        *dlat = -999.0;
        *dlon = -999.0;
        strcpy(grid_name,"Semi-staggered Arakawa E-Grid (2D)");
    } else if (gds && GDS_fEgrid(gds)) {
        *first_lat = -999.0;
        *first_lon = -999.0;
        *last_lat = -999.0;
        *last_lon = -999.0;
        *dlat = -999.0;
        *dlon = -999.0;
        strcpy(grid_name,"Filled Arakawa E-Grid");
    } else if (gds && GDS_RotLL(gds)) {
        *first_lat = -999.0;
        *first_lon = -999.0;
        *last_lat = -999.0;
        *last_lon = -999.0;
        *dlat = -999.0;
        *dlon = -999.0;
        strcpy(grid_name,"Rotated LatLon");
    } else if (gds && GDS_Gnomonic(gds)) {
        *first_lat = -999.0;
        *first_lon = -999.0;
        *last_lat = -999.0;
        *last_lon = -999.0;
        *dlat = -999.0;
        *dlon = -999.0;
        strcpy(grid_name,"Gnomonic");
    } else if (gds && GDS_Harmonic(gds)) {
        *first_lat = -999.0;
        *first_lon = -999.0;
        *last_lat = -999.0;
        *last_lon = -999.0;
        *dlat = -999.0;
        *dlon = -999.0;
        strcpy(grid_name,"Harmonic Spectral");
    } else if (gds && GDS_Triangular(gds)) {
        *first_lat = -999.0;
        *first_lon = -999.0;
        *last_lat = -999.0;
        *last_lon = -999.0;
        *dlat = -999.0;
        *dlon = -999.0;
        strcpy(grid_name,"Triangular");
    }
}
