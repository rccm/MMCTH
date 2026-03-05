/*$Id: bilinear.c,v 1.1.1.1 2006/10/02 22:40:37 mpav Exp $*/

#include <math.h>

void bilinear (int nx, int ny, float *array, 
               int nx_interp, int ny_interp, void **arr_interp)

{
    
  int
    irow,
    icol,
    index,
    index1,
    index2,
    index3,
    index4,
    i,
    j,
    ip,
    jp;  
    
  long
    npts;
    
  float
    ix,
    jy,
    dx,
    dy,
    dx1,
    dy1,
    *arr_interp_f32;
  
  npts = (long) (nx_interp*ny_interp);
  arr_interp_f32 = (float *) malloc(npts*sizeof(float));  
  *arr_interp = arr_interp_f32; 
  
  index = 0;
  for (irow=0; irow < ny_interp; irow++) {
    jy = irow*((float)ny/(ny_interp));
    j = (int)jy;
    jp = j + 1;
    dy = jy - (float)j;
    dy1 = 1.0 - (float)dy;
    for (icol=0; icol < nx_interp; icol++) {
      ix = icol*((float)nx/(nx_interp));
      i = (int)ix;
      ip = i + 1;
      dx = ix - (float)i;
      dx1 = 1.0 - (float)dx;
      index1 = i + (j*(nx));
      index2 = i + (jp*(nx));
      index3 = ip + (j*(nx));
      index4 = ip + (jp*(nx));
      /*If array indices are out of bounds, use nearest neighbor sampling*/
      /*Although, linear extrapolation should be used instead.*/
      if (ip >= nx || jp >= ny) {
	arr_interp_f32[index] = array[index1];
      }
      else {
        arr_interp_f32[index] = array[index1]*dx1*dy1 +
                            array[index2]*dx1*dy  +
		    	    array[index3]*dx*dy1  +
			    array[index4]*dx*dy;
      }
      index++;
    }
  }
  
}
