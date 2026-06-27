      subroutine init_var(cube,pxldat,pxl250,indat,v250_dat,
     +                    contx_lat,contx_lon,contx_elev,contx_vzen,
     +                    contx_szen,contx_rel,contx_pw,contx_eco,
     +                    contx_topog,contx_ice,contx_snow,contx_sfctmp,
     +                    contx_msl,contx_ugrd,contx_vgrd,contx_sst,
     +                    s_type,s_pixels,n_nadir,nadir,scan_pixels,
     +                    scan_number,scan_type,Sstart_time,max_sol,
     +                    min_sol,ilin,mlin,klin,I_ind,prev_cube)

      implicit none
      save

c----------------------------------------------------------------------
C!F77 
c
c!Description:
C     Routine for initializing variables used in the cloud mask
c     MODIS software.
C
c!Input parameters:
c cube          Scan cube counter
c pxldat        Pixel BT and Reflectance data holder
c pxl250        Pixel 250 m reflectance data holder
c indat         Context of lines of 1km data
c v250_dat      Context of lines of 250m data
c contx_lat     Array containing context of lines of latitude values
c contx_lon     Array containing context of lines of longitude values
c contx_elev    Array containing context of lines of elevation values
c contx_vzen    Array containing context of lines of v. zenith values
c contx_szen    Array containing context of lines of s. zenith values
c contx_rel     Array containing context of lines of rel. angle values
c contx_pw      Array containing context of lines of pw values 
c contx_eco     Array containing context of lines of ecosystem values
c contx_topog   Array containing context of lines of land/sea values
c contx_ice     Array containing context of lines of ice fraction values
c contx_snow    Array containing context of lines of snow pixels
c contx_sfctmp  Array containing context of lines of land surface temp
c contx_msl     Array containing context of lines of mean sea level P
c contx_ugrd    Array containing context of lines of u wind component
c contx_vgrd    Array containing context of lines of v wind component
c contx_sst     Array containing context of lines of sea surface temp
c s_type        Scan type - Day or night for lines in a context
c s_pixels      Number of pixels in this scan for lines in a context
c n_nadir       Nadir pixel number for lines in a context
C nadir         Nadir pixel number
C scan_pixels   Number of pixels in this scan
c scan_number   Scan number in granule
C scan_type     Scan type - Day or night
C Sstart_time   Starting scan time (TAI)
c max_sol       Maximum solar zenith angle of this granule
c min_sol       Minimum solar zenith anlge of this granule
c
c!Output Parameters:
c ilin          Counts total number of lines processed
c mlin          Counts number of lines in current context
c klin          Counts number of lines output to bit file
c I_ind         Current cube line (1-10) that is being processed
c prev_cube     Number of cubes skipped or already processed
c               (Initialized to -999)
c
c!Revision History:
c
c!Team-unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c
c!END
c----------------------------------------------------------------------

      INCLUDE 'global.inc'

c     scalar arguments
      character*1 scan_type
      double precision Sstart_time
      real max_sol,min_sol
      integer ilin,mlin,klin,I_ind,prev_cube,cube,
     +        scan_number,scan_pixels,nadir

c     array arguments
      character*1 s_type(nlcntx)
      integer contx_topog(npixel,nlcntx),s_pixels(nlcntx),
     +        n_nadir(nlcntx),contx_snow(npixel,nlcntx)
     +                   
      real    pxldat(inband),pxl250(vis_band,num250_per_1km),
     +        contx_lat(npixel,nlcntx),contx_lon(npixel,nlcntx),
     +        contx_elev(npixel,nlcntx),contx_sst(npixel,nlcntx),
     +        contx_vzen(npixel,nlcntx),contx_szen(npixel,nlcntx),
     +        contx_rel(npixel,nlcntx),v250_dat(nx*4,nlcntx,vis_band),
     +        indat(npixel,nlcntx,inband),contx_pw(npixel,nlcntx),
     +        contx_ice(npixel,nlcntx),contx_sfctmp(npixel,nlcntx),
     +        contx_msl(npixel,nlcntx),contx_ugrd(npixel,nlcntx),
     +        contx_vgrd(npixel,nlcntx)

      byte    contx_eco(npixel,nlcntx)
     
c     local scalars
      integer i,j,k

c ... Initialize arrays
c ... Array holding current pixel band information
      do 30 i = 1 , inband
        pxldat(i) = 0.0
  30  continue

c     Initialize geolocation and ancillary data variables
      do 40 i = 1 , npixel
        do 45 j = 1 , nlcntx
          contx_lat(i,j) = bad_data
          contx_lon(i,j) = bad_data
          contx_elev(i,j) = bad_data
          contx_vzen(i,j) = bad_data
          contx_szen(i,j) = bad_data
          contx_rel(i,j) = 999.0
          contx_pw(i,j) = bad_data
          contx_eco(i,j) = -1
          contx_topog(i,j) = -1
          contx_ice(i,j) = 0.0
          contx_snow(i,j) = 0
          contx_sfctmp(i,j) = 0.0
          contx_msl(i,j) = 0.0
          contx_ugrd(i,j) = 0.0
          contx_vgrd(i,j) = 0.0
          contx_sst(i,j) = 0.0

c ...     Main array containing nlcntx lines of data for each band
          do 50 k = 1 , inband
            indat(i,j,k) = -99.0
  50      continue
  45    continue
  40  continue

c ... Initalize swath metadata scalar variables
      scan_type = ' '
      scan_number = 0
      scan_pixels = 0
      nadir = 0
      Sstart_time = 0.0d0

c ... Initialize metadata holders needed for each line
      do 55 i = 1 , nlcntx
         s_type(i) = ' '
         s_pixels(i) = 0
         n_nadir(i) = 0
  55  continue

c ... Array holding current 250m res. visible reflectance data
      do 60 i = 1 , vis_band
      do 60 j = 1 , num250_per_1km
           pxl250(i,j) = 0.0
  60  continue

c ... Array holding nlcntx lines of data for each vis band
      do 80 i = 1 , vis_band
      do 80 j = 1 , nlcntx
      do 80 k = 1 , nx * 4
         v250_dat(k,j,i) = -99.0
  80  continue


c     Initialize counters
      cube = 0
      ilin = 0
      mlin = 0
      klin = 0
      I_ind = 0
      prev_cube = -999

c     Initialize max and min 
      max_sol = -999.99
      min_sol = 999.99

      return
      end

