/*----------------------------------------------------------------------
#    Copyright (C) 2011,  Space Science and Engineering Center,
#    University C  of Wisconsin-Madison, Madison WI.
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#------------------------------------------------------------------------*/
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 Compilation:  See the file README.C for instructions on compiling this file.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/ 
#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <math.h>
#include "globals.h"

/* Converts satellite zenith angle to scan angle*/


float zen2scanner(float sza)

{
float degrad = 0.017453292;          /*; Degrees to radians conversion factor*/
float Rearth = 6357.;                /*; Radius of the earth in km*/
float Hsat = 850.;                   /*; Nominal height of the satellite*/
float ReHsat = Rearth + Hsat;        /*; Rearth + Hsat */
float scan;

scan = asin(Rearth*sin(sza*degrad)/(ReHsat)) / degrad; 

return scan;

}
