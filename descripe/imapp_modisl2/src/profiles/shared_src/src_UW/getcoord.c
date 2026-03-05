/*  
!C **********************************************************************
!Description:
    GIHLS211 transforms latitude, longitude coordinates to line, sample 
    Renamed GIHLS211 to GETCOORD !!
    coordinates for an image in the Goode's Interrupted Homolosine 
    projection.  This routine was complied and run using the C compiler on 
    SunOS 4.2 (UNIX).  Results were accurate at the time of testing, but they
    are by no means guaranteed!

!Input Parameters:  
lt      latitude (double precision)
ln      longitude (double precision) 

!Output Parameters:
lne     line number of 1km mask for given latitude/longitude
smp     element number of 1km maks for given latitude/longitude

!Revision History:

!Team-unique Header:

!References and Credits:

1.  Snyder, John P. and Voxland, Philip M., "An Album of Map Projections",
    U.S. Geological Survey Professional Paper 1453 , United State Government
    Printing Office, Washington D.C., 1989.

2.  Snyder, John P., "Map Projections--A Working Manual", U.S. Geological
    Survey Professional Paper 1395 (Supersedes USGS Bulletin 1532), United
    State Government Printing Office, Washington D.C., 1987.

3.  Goode, J.P., 1925,  The Homolosine projection:  a new device for
    portraying the Earth's surface entire:  Assoc. Am. Geographers, Annals,
    v. 15, p. 119-125

4.  Steinwand, Daniel R., "Mapping Raster Imagery to the Interrupted
    Goode Homolosine Projection", 1993.  In press--IJRS.

!Design Notes:
!END************************************************************************/
#include <stdio.h>
#include <math.h>


/*  Modified by JC Guu  08/27/96             */
/*  function prototype is given here for the */
/*  functions in this module.                */

void getcoord_(double *lt, double *ln, double *lne, double *smp);
void bad_input_parms();
void goode_init(double r); 
int  goode_forward(double lon, double lat, double* x, double* y);
/*
void ptitle(char* A);
void radius(double A);
*/
void giherror(char* what, char* where);
void getcoord_sincos(double val, double* sin_val, double* cos_val);
int getcoord_sign(double x);
double getcoord_adjust_lon(double x);
double getcoord_e0fn(double x);
double getcoord_e1fn(double x);
double getcoord_e2fn(double x);
double getcoord_mlfn(double e0,double e1,double e2,double phi);


#define PI      3.141592653589793238
#define HALF_PI PI*0.5
#define TWO_PI  PI*2.0
#define EPSLN   1.0e-10
#define R2D     57.2957795131
#define D2R     0.0174532925199

#define OK      1
#define ERROR  -1
#define IN_BREAK -2

/* Variables common to all subroutines in this code file
  -----------------------------------------------------*/
static double R;                /* Radius of the earth (sphere) */
static double lon_center[12];   /* Central meridians, one for each region */
static double feast[12];        /* False easting, one for each region */

/* Transformation routine
  -------------------------------*/
void getcoord_(double *lt, double *ln, double *lne, double *smp)     
{
float pixsiz;
double x,y;
double lat, lon, line, samp;
/*
int nl,ns;
*/
double ul_x, ul_y;

/*  Modified by JC Guu  08/30/96             */
/*  The variable status is added here to     */
/*  process the return code of               */
/*  goode_forward().                         */

int status;

/*  Modified by JC Guu  08/26/96             */
/*  The type of constant is specified.       */

pixsiz = 1000.0F;

/* Report parameters and image geometric characteristics to the user
  -----------------------------------------------------------------*/
/*printf("Converting latitude, longitude to line, sample coordinates\n\n"); 
printf("Pixel size is %f km\n\n", pixsiz/1000.0);  */

/*  Modified by JC Guu  08/26/96             */
/*  The equality comparision                 */
/*  "pixsiz == 1000.0" is changed to an      */
/*  inequality comparision                   */
/*  "fabs(pixsiz-1000.0F) < 0.000001F".      */

if (fabs(pixsiz-1000.0F) < 0.000001F) {
   ul_x = -20015000.0;
   ul_y =   8673000.0;
/*
   nl = 17347;
   ns = 40031;
*/
   } 
else {
   ul_x = -20011500.0;
   ul_y =   8669500.0;
/*
   nl = 2168;
   ns = 5004;
*/
   } 
/*printf("Image size is %d lines by %d samples with an upper left\n",nl,ns);
printf("corner of UL_X = %lf and UL_Y = %lf meters.\n", ul_x, ul_y);   */

/* Initialize the Interrupted Goode Homolosine projection
  ------------------------------------------------------*/
goode_init(6370997.0);

/* Process point
  ------------------*/
   lat = *lt;
   lon = *ln;
/* printf("%12.6f %12.6f\n",lat, lon);  */
   lat *= D2R;
   lon *= D2R;

/*  Modified by JC Guu  08/30/96           */
/*  The variable status is used to         */
/*  capture the return value form          */
/*  goode_forward().                       */
/*  Valid return code to process the       */
/*  return status could be added in the    */
/*  furture.                               */

   status=goode_forward(lon, lat, &x, &y);
   if(status==ERROR)
       /*printf("Newton-Raphson method failed to converge\n");*/
       ;

   line = (ul_y - y) / pixsiz + 1.0;
   samp = (x - ul_x) / pixsiz + 1.0;
/* printf("%12.6f %12.6f\n",lat, lon); 
   printf("%12.6f %12.6f\n",line, samp);   */

   *lne = line;
   *smp = samp;
   *lt = lat;
   *ln = lon;

   return;
}          



/*  Modified by JC Guu  08/26/96              */
/*  The return type void is given.            */

void bad_input_parms() {

/*  Modified by JC Guu  08/30/96              */
/*  The required prolog components are added  */
/*  here.                                     */

/*
!C
!Description: Function to report bad parameters.
!Input Parameters: none
!Output Parameters: none
!Revision History:
!Team-unique Header:
!END
*/

/*  Modified by JC Guu  09/03/96              */
/*  The prohibited printf() code is           */
/*  commented out.                            */

/* printf("Syntax: gihll2ls pixsize\n");
printf("         pixsize in km = 1 or 8\n\n");*/
}



/*  Modified by JC Guu  08/26/96              */
/*  The old K&R style function definition     */
/*  is changed to the modern form and the     */
/*  return type is given as void.             */

void goode_init(double r) 
{

/*  Modified by JC Guu  08/30/96              */
/*  The required prolog components are added  */
/*  here.                                     */

/*
!C
!Description: Initialize the Goode`s Homolosine projection.
              Place parameters in static storage for common use.
!Input Parameters:   r                Radius of the earth (sphere)
!Output Parameters:  lon_center[12]   Central meridians, one for 
                                      each region
                     feast[12];       False easting, one for each 
                                      region
!Revision History:
!Team-unique Header:
!END
*/

/* Place parameters in static storage for common use
  -------------------------------------------------*/

R = r;

/* Initialize central meridians for each of the 12 regions
  -------------------------------------------------------*/
lon_center[0] = -1.74532925199;         /* -100.0 degrees */
lon_center[1] = -1.74532925199;         /* -100.0 degrees */
lon_center[2] =  0.523598775598;        /*   30.0 degrees */
lon_center[3] =  0.523598775598;        /*   30.0 degrees */
lon_center[4] = -2.79252680319;         /* -160.0 degrees */
lon_center[5] = -1.0471975512;          /*  -60.0 degrees */
lon_center[6] = -2.79252680319;         /* -160.0 degrees */
lon_center[7] = -1.0471975512;          /*  -60.0 degrees */
lon_center[8] =  0.349065850399;        /*   20.0 degrees */
lon_center[9] =  2.44346095279;         /*  140.0 degrees */
lon_center[10] = 0.349065850399;        /*   20.0 degrees */
lon_center[11] = 2.44346095279;         /*  140.0 degrees */

/* Initialize false eastings for each of the 12 regions
  ----------------------------------------------------*/
feast[0] = R * -1.74532925199;
feast[1] = R * -1.74532925199;
feast[2] = R * 0.523598775598;
feast[3] = R * 0.523598775598;
feast[4] = R * -2.79252680319;
feast[5] = R * -1.0471975512;
feast[6] = R * -2.79252680319;
feast[7] = R * -1.0471975512;
feast[8] = R * 0.349065850399;
feast[9] = R * 2.44346095279;
feast[10] = R * 0.349065850399;
feast[11] = R * 2.44346095279;

/* Report parameters to the user
  -----------------------------*/
/*ptitle("Goode`s Homolosine Equal Area"); 
radius(r):   */

/* Modified by JC Guu  08/26/96              */
/* The function should not return anything   */
/* here.                                     */

}



/*  Modified by JC Guu  08/26/96              */
/*  The old K&R style function definition     */
/*  is changed to the modern form and the     */
/*  return type is given as void.             */

int goode_forward(double lon, double lat, double* x, double* y)

/*  Modified by JC Guu  08/30/96              */
/*  The required prolog components are added  */
/*  here.                                     */

/*
!C
!Description:Goode`s Homolosine forward equations--mapping 
             lat,long to x,y
!Input Parameters:  lon      Longitude
                    lat      Latitude
!Output Parameters: x        X projection coordinate
                    y        Y projection coordinate                    
!Revision History:
!Team-unique Header:
!END
*/

{

/*  Modified by JC Guu  08/30/96              */
/*  Function prototype is provided for        */
/*  getcoord_adjust_lon() therefore the following      */
/*  declaration is not needed.                */

/* double getcoord_adjust_lon();*/ /* Function to adjust longitude to -180 - 180 */
double delta_lon;       /* Delta longitude (Given longitude - center */
double theta;
double delta_theta;
double constant;
int i;
int region;

/* Forward equations
  -----------------*/
if (lat >= 0.710987989993)                   /* if on or above 40 44' 11.8" */
   {
   if (lon <= -0.698131700798) region = 0;   /* If to the left of -40 */
   else region = 2;
   }
else if (lat >= 0.0)                         /* Between 0.0 and 40 44' 11.8" */
   {
   if (lon <= -0.698131700798) region = 1;   /* If to the left of -40 */
   else region = 3;
   }
else if (lat >= -0.710987989993)             /* Between 0.0 & -40 44' 11.8" */
   {
   if (lon <= -1.74532925199) region = 4;       /* If between -180 and -100 */
   else if (lon <= -0.349065850399) region = 5; /* If between -100 and -20 */
   else if (lon <= 1.3962634016) region = 8;    /* If between -20 and 80 */
   else region = 9;                             /* If between 80 and 180 */
   }
else                                            /* Below -40 44' */
   {
   if (lon <= -1.74532925199) region = 6;       /* If between -180 and -100 */
   else if (lon <= -0.349065850399) region = 7;     /* If between -100 and -20 */
   else if (lon <= 1.3962634016) region = 10;   /* If between -20 and 80 */
   else region = 11;                            /* If between 80 and 180 */
   }

if (region==1||region==3||region==4||region==5||region==8||region==9)
   {
   delta_lon = getcoord_adjust_lon(lon - lon_center[region]);
   *x = feast[region] + R * delta_lon * cos(lat);
   *y = R * lat;
   }
else
   {
   delta_lon = getcoord_adjust_lon(lon - lon_center[region]);
   theta = lat;
   constant = PI * sin(lat);

/* Iterate using the Newton-Raphson method to find theta
  -----------------------------------------------------*/
   for (i=0;;i++)
      {
      delta_theta = -(theta + sin(theta) - constant) / (1.0 + cos(theta));
      theta += delta_theta;
      if (fabs(delta_theta) < EPSLN) break;
      if (i >= 30) 
         {giherror("Iteration failed to converge","Goode-forward");return(ERROR);}
      }
   theta /= 2.0;
   *x = feast[region] + 0.900316316158 * R * delta_lon * cos(theta);

/*  Modified by JC Guu                        */
/*  The explicit cast is provided.            */

   *y = R * (1.4142135623731 * sin(theta) - 0.0528035274542 * ((double) getcoord_sign(lat)));
   }


return(OK);
}



/*  Modified by JC Guu  08/27/96              */
/*  The function definition style is changed  */
/*  the K&R style to the modern form.         */

/*
void ptitle(char* A) 
*/

/*  Modified by JC Guu  08/30/96              */
/*  The required prolog components are added  */
/*  here.                                     */

/*
!C
!Description: Functions to report projection parameter
!Input Parameters: none
!Output Parameters: none
!Revision History:
!Team-unique Header:
!END
*/

/*  Modified by JC Guu  09/03/96              */
/*  printf() is a prohibited function.        */

/* {printf("\n%s Projection Parameters:\n\n",A);} */



/*
void radius(double A) 
*/

/*  Modified by JC Guu  08/30/96              */
/*  The required prolog components are added  */
/*  here.                                     */

/*
!C
!Description:
!Input Parameters: none
!Output Parameters: none
!Revision History:
!Team-unique Header:
!END
*/

/*  Modified by JC Guu  08/27/96             */
/*  The non-standard format specification    */
/*  %lf is changed to %f.                    */

/*  Modified by JC Guu  09/03/96             */
/*  The prohibited function ptintf() is      */
/*  commented out.                           */

/*{ printf("   Radius of Sphere:     %f meters\n",A); } */



/*  Modified by JC Guu  08/26/96              */
/*  The old K&R style function definition     */
/*  is changed to the modern form and the     */
/*  return type is given as void.             */

void giherror(char* what, char* where) 

/*  Modified by JC Guu  08/30/96              */
/*  The required prolog components are added  */
/*  here.                                     */

/*
!C
!Description: Function to report errors
!Input Parameters: what     cause of the error
                   where    where the error occurs
!Output Parameters: none
!Revision History:
!Team-unique Header:
!END
*/

/*  Modified by JC Guu  09/03/96              */
/*  printf() is a prohibited function.        */
/* commented out by fhliang on 12/08/97 */
{
/* printf("[%s] %s\n",where,what); */
}


#ifndef sun



/*  Modified by JC Guu  08/26/96              */
/*  The old K&R style function definition     */
/*  is changed to the modern form and the     */
/*  return type is given as void.             */

void getcoord_sincos(double val, double* sin_val, double* cos_val)

/*  Modified by JC Guu  08/30/96              */
/*  The required prolog components are added  */
/*  here.                                     */

/*
!C
!Description: Function to calculate the sine and cosine in 
              one call.  Some computer systems have implemented 
              this function, resulting in a faster implementation
              than calling each function separately.  It is 
              provided here for those computer systems which 
              don`t implement this function.

!Input Parameters:  val         the value to be operated on
!Output Parameters: sin_val     the sine of val
                    cos_val     the cosine of val
!Revision History:
!Team-unique Header:
!END
*/

/*  Modified by JC Guu  08/26/96              */
/*  A return statement is not needed here.    */

{ *sin_val = sin(val); *cos_val = cos(val); }
#endif



/*  Modified by JC Guu  08/27/96              */
/*  The old K&R style function definition     */
/*  is changed to the modern form and the     */
/*  return type is given as int.              */

int getcoord_sign(double x) 
{ if (x < 0.0) return(-1); else return(1);}

/*  Modified by JC Guu  08/30/96              */
/*  The required prolog components are added  */
/*  here.                                     */

/*
!C
!Description: Function to return the sign of an argument
!Input Parameters:  x    a floating point number
!Output Parameters: the sign of x
!Revision History:
!Team-unique Header:
!END
*/



/*  Modified by JC Guu  08/27/96              */
/*  The old K&R style function definition     */
/*  is changed to the modern.                 */

double getcoord_adjust_lon(double x)

/*  Modified by JC Guu  08/30/96              */
/*  The required prolog components are added  */
/*  here.                                     */

/*
!C
!Description: Function to adjust longitude to -180 to 180
!Input Parameters: x  
!Output Parameters:
!Revision History:
!Team-unique Header:
!END
*/

/*  Modified by JC Guu  08/27/96              */
/*  The explicit cast is provided.            */

{x=(fabs(x)<PI)?x:(x-(((double) getcoord_sign(x))*TWO_PI));return(x);}



/*  Modified by JC Guu  08/27/96                */
/*  The style of function definition is changed */
/*  from the K&R style to the modern form for   */
/*  the following functions.                    */

double getcoord_e0fn(double x) 

/*  Modified by JC Guu  08/30/96              */
/*  The required prolog components are added  */
/*  here.                                     */

/*
!C
!Description: Function to compute e0
!Input Parameters: x
!Output Parameters: 
!Revision History:
!Team-unique Header:
!END
*/

{return(1.0-0.25*x*(1.0+x/16.0*(3.0+1.25*x)));}



double getcoord_e1fn(double x)

/*  Modified by JC Guu  08/30/96              */
/*  The required prolog components are added  */
/*  here.                                     */

/*
!C
!Description: Function to compute e1
!Input Parameters: x
!Output Parameters: 
!Revision History:
!Team-unique Header:
!END
*/

{return(0.375*x*(1.0+0.25*x*(1.0+0.46875*x)));}

 

double getcoord_e2fn(double x) 

/*  Modified by JC Guu  08/30/96              */
/*  The required prolog components are added  */
/*  here.                                     */

/*
!C
!Description: Function to compute e2
!Input Parameters: x
!Output Parameters:
!Revision History:
!Team-unique Header:
!END
*/

{return(0.05859375*x*x*(1.0+0.75*x));}



double getcoord_mlfn(double e0,double e1,double e2,double phi) 

/*  Modified by JC Guu  08/30/96              */
/*  The required prolog components are added  */
/*  here.                                     */

/*
!C
!Description:
!Input Parameters: e0, e1, e2 and phi
!Output Parameters:
!Revision History:
!Team-unique Header:
!END
*/

{
return(e0*phi-e1*sin(2.0*phi)+e2*sin(4.0*phi));
}
