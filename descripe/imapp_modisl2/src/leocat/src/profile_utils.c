/*$Id: profile_utils.c,v 1.1.1.2 2006/12/05 14:27:49 mpav Exp $*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "common_leocat.h"
#include "nwp_leocat.h"

#define GRAV_ACCEL_EARTH	9.80665			/* m / s^2	*/
#define SPEC_GAS_CONST_AIR	287.05			/* J / kg K	*/
#define SPEC_GAS_CONST_H2O	461.51			/* J / kg K	*/

/*******************************************************************************
/
*******************************************************************************/

int profile_to_101(double *p, double *t, double *w, int n, double lat,
                   double *pp, double *tt, double *ww) {
     int i;

     int iy;
     int iz;
     int nlx = 101;
     
     int i_str;
     int i_end;

     double anum;
     double aden;

     double rlogp;

     double delt;
     double delw;

     double wmin = 0.003;
     double wmax;

     double pb[2];
     double tb[2];
     double wb[2];

     if (n < 2) {
          printf("number of levels must be atleast 2\n");
          return -1;
     }

     /*make_profile_101(pp);*/
     int_levels_pp(p, t, w, n, pp, tt, ww, nlx, &i_str, &i_end);

     if (i_str >= 35) {
          printf("ERROR: temperature profile doesn't go high enough\n");
          return -1;
     }

     for (i = 0; i < i_str; ++i) {
          tt[i] = -1.;
          ww[i] = wmin;
     }

     extem101_64_(tt, &lat);
     
     iz = n - 1;

     if (p[iz] < pp[nlx-1]) {
          pb[0] = p[iz];
          tb[0] = t[iz];
          wb[0] = w[iz];

          pb[1] = pp[nlx-1];

          iy = n - 2;
/*
          while (p[iy] >= p[iz])
               --iy;
*/
          anum = log(pb[1] / pb[0]);
//        aden = log(pb[0] / p[iy]);
          aden = log(pb[1] / p[iy]);

          rlogp = anum / aden;

          delt = t[iz] - t[iy];
          tb[1] = t[iz] + delt * rlogp;

          delw = w[iz] - w[iy];
          wb[1] = w[iz] + delw * rlogp;

          if (wb[1] < wmin)
               wb[1] = wmin;
          else {
               wmax = satmix(pb[1], tb[1]);

               if (wb[1] > wmax)
                    wb[1] = wmax;
          }

          i = i_end + 1;

          int_levels_pp(pb, tb, wb, 2, pp+i, tt+i, ww+i, nlx-i, &i_str, &i_end);
     }

     return 0;
}

/*******************************************************************************
/
*******************************************************************************/
int make_profile_101(double *p) {

     int i;

     double l;

     double a = -1.550789414500298e-04;
     double b = -5.593654380586063e-02;
     double c =  7.451622227151780e+00;

     l = 101;
     for (i = 0; i < 101; ++i) {
            p[i] = pow((a*l*l + b*l + c), (7./2.));
            l = l - 1;
     }

     return 0;
}

/*******************************************************************************
/ 
*******************************************************************************/
int int_levels_pp(double *p,  double *t,  double *w, int n,
                  double *pp, double *tt, double *ww, int nn, int *i_str, int *i_end) {

     int i;
     int ii;

     double dl;

     double slope_t;
     double slope_w;

     if (n < 2) {
          *i_str =  0;
          *i_end = -1;

          return 0;
     }

     for (i = 0; i < nn; ++i) {
          if (pp[i] >= p[0])
               break;

          tt[i] = 0.;
          ww[i] = 0.;
     }

     *i_str = i;

     ii = 1;

     dl    = log(p[ii] / p[ii-1]);
     slope_t = (t[ii] - t[ii-1]) / dl;
     slope_w = (w[ii] - w[ii-1]) / dl;

     for (     ; i < nn;    ) {
          if (pp[i] < p[ii]) {
               dl      = log(pp[i] / p[ii-1]);
               tt[i] = t[ii-1] + slope_t * dl;
               ww[i] = w[ii-1] + slope_w * dl;

               ++i;
          }
          else {
               ++ii;
               if (ii >= n) {
                    if (p[ii-1] == pp[i]) {
                         tt[i] = t[ii-1];
                         ww[i] = w[ii-1];

                         ++i;
                    }

                    break;
               }

               dl    = log(p[ii] / p[ii-1]);
               slope_t = (t[ii] - t[ii-1]) / dl;
               slope_w = (w[ii] - w[ii-1]) / dl;
          }
     }

     *i_end = i - 1;

     for (     ; i < nn; ++i) {
          tt[i] = 0.;
          ww[i] = 0.;
     }

     return 0;
}

/*******************************************************************************
/
*******************************************************************************/
int height_profile(double *p, double *t, double *w, double *z, int n, double p0) {

     int i;

     int i_up;
     int i_dn;

     double g;

     double Rd;
     double Rv;

     double epsilon;

     double t0;
     double w0;

     double z_cur;

     double p_last;
     double t_last;
     double w_last;

     double P;
     double T;
     double W;

     double e;

     double rho;

     g = GRAV_ACCEL_EARTH;

     Rd = SPEC_GAS_CONST_AIR;
     Rv = SPEC_GAS_CONST_H2O;

     epsilon = Rd / Rv;


     for (i = 0; i < n; ++i) {
          if (p[i] > p0)
               break;
     }

     if (i == 0) {
          i_up = i;
          i_dn = i+1;
     }
     else {
          i_up = i-1;
          i_dn = i;
     }

     t0 = t[i_up] + lin_int(p[i_up], p[i_dn], p0, t[i_up], t[i_dn]);
     w0 = w[i_up] + lin_int(p[i_up], p[i_dn], p0, w[i_up], w[i_dn]);

     z_cur = 0;

     p_last = p0;
     t_last = t0;
     w_last = w0;
     for (i = i_up; i >= 0; --i) {
          P = (p_last + p[i]) / 2. * 100.;
          T = (t_last + t[i]) / 2.;
          W = (w_last + w[i]) / 2. / 1000.;

          e = W / (W + epsilon) * P;

          rho = (P-e) / (Rd*T) + e / (Rv*T);

          z_cur += (p_last-p[i])*100. / (g*rho);

          z[i] = z_cur;

          p_last = p[i];
          t_last = t[i];
          w_last = w[i];
     }

     z_cur = 0;

     p_last = p0;
     t_last = t0;
     w_last = w0;
     for (i = i_dn; i < n; ++i) {
          P = (p_last + p[i]) / 2. * 100.;
          T = (t_last + t[i]) / 2.;
          W = (w_last + w[i]) / 2. / 1000.;

          e = W / (W + epsilon) * P;

          rho = (P-e) / (Rd*T) + e / (Rv*T);

          z_cur += (p_last-p[i])*100. / (g*rho);

          z[i] = z_cur;

          p_last = p[i];
          t_last = t[i];
          w_last = w[i];
     }

     return 0;
}

void get_nwp_time(int year, int jday, float time,int forecast_len,
                  int *year1, int *year2, int *ileap1, int *ileap2,
		  int *jday1, int *jday2, int *hour1, int *hour2)

{
  char *rout = {"get_nwp_time"};
  int ileap_curr_year=0,ileap_prev_year=0;
  int year_prev;
  
  /*---------------------------------------------------------------------------
  Determine if there is a leap year to deal with.
  ---------------------------------------------------------------------------*/
  
  ileap_curr_year = leap_year_check(year);
      
  /*Check to see if the current (year - 1) was a leap year*/
  year_prev = year - 1;
  ileap_prev_year = leap_year_check(year_prev);
  
  if (forecast_len == 0) {

//  Changed if test formulation in the following three if blocks - see commented lines. (RAF)
    
//  if (time >= 0.00 && time <= 6.00) {
    if (time >= 0.00 && time < 6.00) {
      *jday1 = jday;
      *hour1 = 0;
      *jday2 = jday;
      *hour2 = 6;
    
    }
//  else if (time > 6.00 && time <= 12.00) {
    else if (time >= 6.00 && time < 12.00) {
      *jday1 = jday;
      *hour1 = 6;
      *jday2 = jday;
      *hour2 = 12;
    }
//  else if (time > 12.00 && time <= 18.00) {
    else if (time >= 12.00 && time < 18.00) {
      *jday1 = jday;
      *hour1 = 12;
      *jday2 = jday;
      *hour2 = 18;
    }
//  else if (time > 18.00 && time <= 24.00) {
    else if (time >= 18.00 && time < 24.00) {
      *jday1 = jday;
      *hour1 = 18;
      *jday2 = jday+1;
      *hour2 = 0;
    }
    else {
      fprintf(stderr,"%s%s-Observation time (%f) is out of bounds\n",EXE_PROMPT,rout,time);
      exit(EXIT_FAILURE);
    }
    
  }
  else if (forecast_len == 6) {
    
    if (time >= 0.00 && time < 6.00) {
      *jday1 = jday - 1;
      *hour1 = 18;
      *jday2 = jday;
      *hour2 = 0;
    
    }
    else if (time >= 6.00 && time < 12.00) {
      *jday1 = jday;
      *hour1 = 0;
      *jday2 = jday;
      *hour2 = 6;
    }
    else if (time >= 12.00 && time < 18.00) {
      *jday1 = jday;
      *hour1 = 6;
      *jday2 = jday;
      *hour2 = 12;
    }
    else if (time >= 18.00 && time < 24.00) {
      *jday1 = jday;
      *hour1 = 12;
      *jday2 = jday;
      *hour2 = 18;
    }
    else {
      fprintf(stderr,"%s%s-Observation time (%f) is out of bounds\n",EXE_PROMPT,rout,time);
      exit(EXIT_FAILURE);
    }
  }
  else if (forecast_len == 12) {
    
    if (time >= 0.00 && time < 6.00) {
      *jday1 = jday - 1;
      *hour1 = 12;
      *jday2 = jday - 1;
      *hour2 = 18;
    
    }
    else if (time >= 6.00 && time < 12.00) {
      *jday1 = jday - 1;
      *hour1 = 18;
      *jday2 = jday;
      *hour2 = 0;
    }
    else if (time >= 12.00 && time < 18.00) {
      *jday1 = jday;
      *hour1 = 0;
      *jday2 = jday;
      *hour2 = 6;
    }
    else if (time >= 18.00 && time < 24.00) {
      *jday1 = jday;
      *hour1 = 6;
      *jday2 = jday;
      *hour2 = 12;
    }
    else {
      fprintf(stderr,"%s%s-Observation time (%f) is out of bounds\n",EXE_PROMPT,rout,time);
      exit(EXIT_FAILURE);
    }
    
  }
  
  /*---------------------------------------------------------------------------
  Determine the correct year, month, and day for each gfs file.
  ---------------------------------------------------------------------------*/
  
  *year1 = year;
  *year2 = year;
  *ileap1 = ileap_curr_year;
  *ileap2 = ileap_curr_year;
  if (*jday1 < 1) {
    *ileap1 = ileap_prev_year;
    *year1 = year - 1;
    *jday1 = 365 + (*ileap1);
  }
  if (*jday2 < 1) {
    *ileap2 = ileap_prev_year;
    *year2 = year - 1;
    *jday2 = 365 + (*ileap2);
  }

  if( *jday2 > (365 + (*ileap1))) {
    *jday2 = 1;
    *year2 = *year2 + 1;
  }


}
