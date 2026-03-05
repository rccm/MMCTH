/*$Id: destroy_nwp_pointers.c,v 1.1.1.1 2006/10/02 22:40:37 mpav Exp $*/

#include <hdf.h>
#include "common_leocat.h"
#include "nwp_leocat.h"

void destroy_nwp_pointers(profile_params *map)
{  
  free(map->p_lev);
  free(map->z_lev);
  free(map->t_lev);
  free(map->o3_lev);
  free(map->adjo3_lev);
  free(map->w_lev);
  
}

