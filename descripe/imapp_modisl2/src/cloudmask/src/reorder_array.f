      subroutine reorder_array(contx_lat,contx_lon,contx_vzen,
     +                         contx_szen,contx_rel,contx_pw,
     +                         contx_eco,contx_topog,contx_ice,
     +                         contx_snow,contx_sfctmp,contx_msl,
     +                         contx_ugrd,contx_vgrd,contx_elev,s_type,
     +                         s_pixels,n_nadir,maxele,indat,
     +                         v250_dat,contx_sst)

      implicit none
      save

c ... Common statement for debug purposes
      common / bug / debug, h_output

c---------------------------------------------------------------------
C!F77 
c
c!Description:
c     Routine for re-ordering array data after completion of
c     processing for a given line.
c
c!Input Parameters:
c contx_lat     Array of scan line latitudes
c contx_lon     Array of scan line longitudes
c contx_vzen    Array of scan line viewing zenith angles
c contx_szen    Array of scan line solar zenith angles
c contx_rel     Array containing context of lines of rel. angle values
c contx_pw      Array containing context of lines of pw values
c contx_eco     Array containing context of lines of ecosystem values
c contx_topog   Array containing context of lines of land/sea values
c contx_ice     Array containing context of lines of ice fraction values
c contx_snow    Array containing context of lines of snow values
c contx_sfctmp  Array containing context of lines of surface temp
c contx_msl     Array containing context of lines of mean sea level P
c contx_ugrd    Array containing context of lines of u wind component
c contx_vgrd    Array containing context of lines of v wind component
c contx_elev    Array containing context of lines of elevation
c s_pixel       Number of pixels in this scan for lines in a context
c s_type        Scan type - Day or night for lines in a context
c n_nadir       Nadir pixel number for lines in a context
c maxele        Maximum element number in context
c indat         Array containing a context of scan line MODIS
c               data
c v250_dat      Array containing a context of scan line MODIS
c               250m visible data
c contx_sst     Array containing context of lines of sst
c
c!Output Parameters:   None
c
c!Revision History:
c
c!Team-unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c
c!END
c---------------------------------------------------------------------

      INCLUDE 'global.inc'

c     scalar arguments
      integer maxele

c     array arguments
      character*1 s_type(nlcntx)
      real contx_lat(npixel,nlcntx),contx_lon(npixel,nlcntx),
     +     contx_vzen(npixel,nlcntx),contx_szen(npixel,nlcntx),
     +     contx_rel(npixel,nlcntx),contx_pw(npixel,nlcntx),
     +     contx_ice(npixel,nlcntx),contx_sfctmp(npixel,nlcntx),
     +     contx_msl(npixel,nlcntx),contx_ugrd(npixel,nlcntx),
     +     contx_vgrd(npixel,nlcntx),indat(npixel,nlcntx,inband),
     +     v250_dat(nx*4,nlcntx,vis_band),contx_elev(npixel,nlcntx),
     +     contx_sst(npixel,nlcntx)
      integer contx_topog(npixel,nlcntx),s_pixels(nlcntx),
     +        contx_snow(npixel,nlcntx),n_nadir(nlcntx)
      byte contx_eco(npixel,nlcntx)

c     local scalars
      integer ie,il,k,ich,debug,h_output

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(10x/,''Subroutine reorder_array  '',/)')
      endif
c ................................................................


c     If not processing in single-line mode, copy context
c     border number of lines into beginning lines to fill 
c     in when reading next line of data
c     First re-roder the 1km variable arrays
c     
      if (nlcntx.gt.1) then
         do 40 il = 2,nlcntx
c              First array are just line holders
               s_type(il-1) = s_type(il)
               s_pixels(il-1) = s_pixels(il)
               n_nadir(il-1) = n_nadir(il)
            do 50 ie = 1,maxele
               contx_lat(ie,il-1) = contx_lat(ie,il)
               contx_lon(ie,il-1) = contx_lon(ie,il)
               contx_vzen(ie,il-1) = contx_vzen(ie,il)
               contx_szen(ie,il-1) = contx_szen(ie,il)
               contx_rel(ie,il-1) = contx_rel(ie,il)
               contx_pw(ie,il-1) = contx_pw(ie,il)
               contx_ice(ie,il-1) = contx_ice(ie,il)
               contx_snow(ie,il-1) = contx_snow(ie,il)
               contx_topog(ie,il-1) = contx_topog(ie,il)
               contx_eco(ie,il-1) = contx_eco(ie,il)
               contx_sfctmp(ie,il-1) = contx_sfctmp(ie,il)
               contx_msl(ie,il-1) = contx_msl(ie,il)
               contx_ugrd(ie,il-1) = contx_ugrd(ie,il)
               contx_vgrd(ie,il-1) = contx_vgrd(ie,il)
               contx_elev(ie,il-1) = contx_elev(ie,il)
               contx_sst(ie,il-1) = contx_sst(ie,il)

               do 60 k = 1,inband
                 indat(ie,il-1,k) = indat(ie,il,k)
  60           continue

  50        continue
  40     continue
      end if

c     Now re-order the 250m visible arrays

      if (nlcntx.gt.1) then
         do 70 il = 2,nlcntx
            do 80 ie = 1,maxele * num250_per_1km
               do 90 ich = 1 , vis_band
                  v250_dat(ie,il-1,ich) = v250_dat(ie,il,ich)
  90           continue
  80        continue
  70     continue
      end if
               
      return
      end
