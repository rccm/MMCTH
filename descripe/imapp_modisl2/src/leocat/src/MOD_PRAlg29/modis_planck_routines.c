/*$Id: modis_planck_routines.c,v 1.1.1.1 2006/10/02 22:40:37 mpav Exp $*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define h_planck 6.62606876E-34		/* J/s */
#define c_planck 2.99792458E+08		/* m/s */
#define k_planck 1.3806503E-23		/* J/K */
#define c1_planck 2.0 * h_planck * c_planck * c_planck
#define c2_planck (h_planck * c_planck) / k_planck

/******************************************************************************/
/* All MODIS planck constants were obtained from IMAPP                        */
/******************************************************************************/

const float cwn_terra[] = {2.641767E+03, 2.505274E+03, 2.518031E+03, 2.465422E+03,
                           2.235812E+03, 2.200345E+03, 0.0000000000, 1.478026E+03,
		           1.362741E+03, 1.173198E+03, 1.027703E+03, 9.081998E+02,
		           8.315149E+02, 7.483224E+02, 7.309089E+02, 7.188677E+02,
		           7.045309E+02};
  
const float tcs_terra[] = {9.993487E-01, 9.998699E-01, 9.998604E-01, 9.998701E-01,
                           9.998825E-01, 9.998849E-01, 0.0000000000, 9.994942E-01,
		           9.994937E-01, 9.995643E-01, 9.997499E-01, 9.995880E-01,
		           9.997388E-01, 9.999192E-01, 9.999171E-01, 9.999174E-01,
		           9.999264E-01};

const float tci_terra[] = {4.744530E-01, 9.091094E-02, 9.694298E-02, 8.856134E-02,
                           7.287017E-02, 7.037161E-02, 0.0000000000, 2.177889E-01,
		           2.037728E-01, 1.559624E-01, 7.989879E-02, 1.176660E-01,
		           6.856633E-02, 1.903625E-02, 1.902709E-02, 1.859296E-02,
		           1.619453E-02};
			   

const float cwn_aqua[] = {2.647418E+03, 2.511763E+03, 2.517910E+03, 2.462446E+03,
                          2.248296E+03, 2.209550E+03, 0.0000000000, 1.474292E+03,
		          1.361638E+03, 1.169637E+03, 1.028715E+03, 9.076808E+02,
		          8.308397E+02, 7.482977E+02, 7.307761E+02, 7.182089E+02,
		          7.035020E+02};
  
const float tcs_aqua[] = {9.993438E-01, 9.998680E-01, 9.998649E-01, 9.998729E-01,
                          9.998738E-01, 9.998774E-01, 0.0000000000, 9.995732E-01,
		          9.994894E-01, 9.995439E-01, 9.997496E-01, 9.995483E-01,
		          9.997404E-01, 9.999194E-01, 9.999071E-01, 9.999176E-01,
		          9.999211E-01};

const float tci_aqua[] = {4.792821E-01, 9.260598E-02, 9.387793E-02, 8.659482E-02,
                          7.854801E-02, 7.521532E-02, 0.0000000000, 1.833035E-01,
		          2.053504E-01, 1.628724E-01, 8.003410E-02, 1.290129E-01,
		          6.810679E-02, 1.895925E-02, 2.128960E-02, 1.857071E-02,
		          1.733782E-02};

/******************************************************************************/
/* This function is used to compute brightness temperature in Kelvin.  Input  */
/* is as follows: wavelength (wl) in um and radiance in W*m-2*str-1*um-1      */
/******************************************************************************/
 
float planck_btemp_terra (int band, float rad)

{

  float
    planckT,
    wl_meters,
    wl;
  
  wl = 1.0e+4/cwn_terra[(int)band-20];
  wl_meters = 1.0e-6 * wl;
  
  planckT = ((c2_planck / (wl_meters * log(c1_planck / 
	     (1.0e+6 * rad * pow(wl_meters,5)) + 1.0))) - 
	     tci_terra[(int)band-20])/tcs_terra[(int)band-20];
  
  return planckT;
  
}

/******************************************************************************/
/* This function is used to compute brightness temperature in Kelvin.  Input  */
/* is as follows: wavelength (wl) in um and radiance in W*m-2*str-1*um-1      */
/******************************************************************************/
 
float planck_btemp_aqua (int band, float rad)

{

  float
    planckT,
    wl_meters,
    wl;
  
  wl = 1.0e+4/cwn_aqua[(int)band-20];
  wl_meters = 1.0e-6 * wl;
  
  planckT = ((c2_planck / (wl_meters * log(c1_planck / 
	     (1.0e+6 * rad * pow(wl_meters,5)) + 1.0))) - 
	     tci_aqua[(int)band-20])/tcs_aqua[(int)band-20];
  
  return planckT;
  
}

float modis_planck_terra(int band, float temp)

{

  float
    wavelength,
    bt,
    rad;

  wavelength = (1.0e+4 / cwn_terra[band - 20])*1.0e-6;
  bt = temp * tcs_terra[band - 20] + tci_terra[band - 20];
  
  rad = 1.0e-6 * (c1_planck / pow(wavelength,5)) / (exp(c2_planck / (wavelength * bt)) - 1.0);
  
  return rad;
  
}

float modis_planck_aqua(int band, float temp)

{

  float
    wavelength,
    bt,
    rad;

  wavelength = (1.0e+4 / cwn_aqua[band - 20])*1.0e-6;
  bt = temp * tcs_aqua[band - 20] + tci_aqua[band - 20];
  
  rad = 1.0e-6 * (c1_planck / pow(wavelength,5)) / (exp(c2_planck / (wavelength * bt)) - 1.0);
  
  return rad;
  
}

double rad_v2w_aqua(int band, double rad) {

  double lamda, I;
  
  lamda = cwn_aqua[band - 20]/1.0e-2;
  I = rad * 1.0e-5;
  
  return 1.0e-6 * (c1_planck / (exp(log(c1_planck*pow(lamda, 3) / I + 1.)) - 1.) / pow(1. / lamda, 5));

}

double rad_v2w_terra(int band, double rad) {

  double lamda, I;
  
  lamda = (double)cwn_terra[band - 20]/1.0e-2;
  I = rad * 1.0e-5;
  
  return 1.0e-6 * (c1_planck / (exp(log(c1_planck*pow(lamda, 3) / I + 1.)) - 1.) / pow(1. / lamda, 5));

}
  
/******************************************************************************/
/******************************************************************************/
