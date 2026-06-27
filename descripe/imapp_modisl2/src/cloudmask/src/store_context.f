      subroutine store_context(mlin,Iline,maxele,rlat,rlon,elev,
     +                         view_zen,solar_zen,rel_angle,v250_band,
     +                         vis250_band,rad_band,v1km_band,indat,
     +                         v250_dat,cube_pw,cube_eco,
     +                         cube_topog,cube_ice,cube_snow,
     +                         cube_sfctmp,cube_msl,cube_ugrd,cube_vgrd,
     +                         cube_sst,
     +                         contx_lat,contx_lon,contx_elev,contx_vzen,
     +                         contx_szen,contx_rel,contx_pw,contx_eco,
     +                         contx_topog,contx_ice,contx_snow,
     +                         contx_sfctmp,contx_msl,contx_ugrd,
     +                         contx_vgrd,contx_sst)

      implicit none
      save

c ... Common statement for debug purposes
      common / bug / debug, h_output

c---------------------------------------------------------------------
c!F77 
c
c!Description:
c     Routine which takes scan cube's worth of a certain variable
c     and puts it into a context of lines (usually 3 lines worth).
c     Use routine modis_bright to convert radiances to brightness
c     temperatures.
C
c!Input parameters:
c mlin          Context line number
c Iline         Scan cube line number (1-10)
c maxele        Maximum number of elements in region
c rlat          Array containing scan cube of 1km pixel latitude
c               values
c rlon          Array containing scan cube of 1km pixel longitude
c               values
c elev          Array containing scan cube of 1km pixel elevation
c               values
c view_zen      Array containing scan cube of 1km pixel viewing
c               zenith angles
c solar_zen     Array containing scan cube of 1km pixel solar
c               azimulth angles
c rel_angle     Array containing scan cube of 1km pixel relative
c               angles
c v250_band     Byte array containing valid flag for a scan of 250m
c               pixels
c vis250_band   Array containing scan cube of 250m visible band
c               reflectances
c v1km_band     Byte array containing valid flag for a scan of 1km
c               pixels
c rad_band      Array containing the 13 1km resolution channel
c               radiances and reflectances for a scan cube
c cube_pw       Array containing scan cube of 1km pixel precipitable
c               water values (g/cm2)
c cube_eco      Array containing scan cube of 1km pixel ecosystem
c               values (Olson World Ecosystem Map)
c cube_topog    Array containing scan cube of 1km pixel topography
c               values (land/sea values 0 or 1)
c cube_ice      Array containing scan cube of 1km pixel sea ice
c               concentration values (fraction)
c cube_snow     Array containing scan cube of 1km pixel snow cover
c               values 
c cube_sfctmp   Array containing scan cube of 1km pixel surface
c               temperatures - from model (K)
c cube_msl      Array containing scan cube of 1km pixel mean sea level
c               pressure values from model (hPa)
c cube_ugrd     Array containing scan cube of 1km pixel surface wind
c               u component (m/s)
c cube_vgrd     Array containing scan cube of 1km pixel surface wind
c               v component (m/s)
c cube_sst      Array containing scan cube of 1km pixel sea surface
c               temperatures - from Reynolds Blended data (K)
c!Output Parameters:
c
c indat         Array containing context of lines of radiance
c               or reflectance information for each channel
c v250_dat      Array containing context of lines of relfectance
c               data for the 2 250m visible channels 
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
c contx_snow    Array containing context of lines of snow values
c contx_sfctmp  Array containing context of lines of land surface temp
c contx_msl     Array containing context of lines of mean sea level P
c contx_ugrd    Array containing context of lines of u wind component
c contx_vgrd    Array containing context of lines of v wind component
c contx_sst     Array containing context of lines of sea surface temp
c
c!Revision History:
c 06/04 Collection 5  R. Frey
c Added 'contx_sst', 'cube_sst'.
c
c!Team-Unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c
c!END
c---------------------------------------------------------------------

      include 'global.inc'

c     scalar arguments
      integer mlin,Iline,maxele
c     array arguments
      real rlon(npixel,scans_cube),rlat(npixel,scans_cube),
     +       elev(npixel,scans_cube),
     +       view_zen(npixel,scans_cube),solar_zen(npixel,scans_cube),
     +       cube_pw(npixel,scans_cube),cube_ice(npixel,scans_cube),
     +       cube_sfctmp(npixel,scans_cube),cube_msl(npixel,scans_cube),
     +       cube_ugrd(npixel,scans_cube),cube_vgrd(npixel,scans_cube),
     +       cube_sst(npixel,scans_cube),
     +       rel_angle(npixel,scans_cube),vis250_band(nx,ny,vis_band),
     +       rad_band(nx,ny,inband),indat(npixel,nlcntx,inband),
     +       contx_lat(npixel,nlcntx),contx_lon(npixel,nlcntx),
     +       contx_elev(npixel,nlcntx),contx_sst(npixel,nlcntx),
     +       contx_vzen(npixel,nlcntx),contx_szen(npixel,nlcntx),
     +       contx_rel(npixel,nlcntx),v250_dat(nx*4,nlcntx,vis_band),
     +       contx_pw(npixel,nlcntx),contx_ice(npixel,nlcntx),
     +       contx_sfctmp(npixel,nlcntx),contx_msl(npixel,nlcntx),
     +       contx_ugrd(npixel,nlcntx),contx_vgrd(npixel,nlcntx)
      integer cube_topog(npixel,scans_cube),contx_topog(npixel,nlcntx),
     +        cube_snow(npixel,scans_cube),contx_snow(npixel,nlcntx)
      byte cube_eco(npixel,scans_cube),contx_eco(npixel,nlcntx),
     +     v250_band(nx,ny,vis_band),v1km_band(nx,ny,inband)

c     local scalars
      real temperature,rad,rf,cos_ang,cos_vis,pi,dtr
      integer i,j,k,nch,iele,ilin,offset,rfoffset,ie,debug,h_output,lsf

      integer units
      parameter (units=1)

c     external functions
      real modis_bright
      external modis_bright

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(10x/,''Subroutine store_context  '',/)')
      endif
c ................................................................

c     initialize variables
      pi = acos(-1.0)
      dtr = pi / 180.0

c ... FIRST, PUT CHANNELS NEEDED FOR 1KM MASK INTO CONTEXT
c ... Store arrays in nlcntx by necntx array for passing to subroutines
      do 185 nch = 1 , inband
         do 190 iele = 1 , maxele
  
c           Initialize reflectance and BT values holders
            rf = 0.0
            rad = 0.0
            temperature = 0.0

c ...       First read in the geolocation information and
c           ancillary data
            if (nch .eq. 1) then
               contx_lat(iele,mlin) = rlat(iele,Iline)
               contx_lon(iele,mlin) = rlon(iele,Iline)
               contx_elev(iele,mlin) = elev(iele,Iline)
               contx_vzen(iele,mlin) = view_zen(iele,Iline)
               contx_szen(iele,mlin) = solar_zen(iele,Iline)
               contx_rel(iele,mlin) = rel_angle(iele,Iline)
               contx_pw(iele,mlin) = cube_pw(iele,Iline)
               lsf = cube_topog(iele,Iline)
               if (cube_eco(iele,Iline).ge.0 .and. 
     +             cube_eco(iele,Iline).le. 94) then 
c ..               Make sure the ecosystem and land/sea file are
c ...              consistent
                   if (lsf.eq.1 .or. lsf.eq.2 .or. lsf.eq.3 .or.
     +                 lsf.eq.4) then
                      contx_eco(iele,mlin) = cube_eco(iele,Iline)
                   else
                      contx_eco(iele,mlin) = 14
                   endif
               else
                   contx_eco(iele,mlin) = -1
               endif
               if (lsf.ge.0 .and. lsf.le.7) then 
                   contx_topog(iele,mlin) = cube_topog(iele,Iline)
               else
                   contx_topog(iele,mlin) = -1
               endif
               contx_ice(iele,mlin) = cube_ice(iele,Iline)
               contx_snow(iele,mlin) = cube_snow(iele,Iline)
               contx_sfctmp(iele,mlin) = cube_sfctmp(iele,Iline)
               contx_msl(iele,mlin) = cube_msl(iele,Iline)
               contx_ugrd(iele,mlin) = cube_ugrd(iele,Iline)
               contx_vgrd(iele,mlin) = cube_vgrd(iele,Iline)
               contx_sst(iele,mlin) = cube_sst(iele,Iline)
            endif

c           get the cosine of the solar zenith angle for reflectance
c           calculations
c ...       First check for valid data
            if (v1km_band(iele,Iline,nch) .eq. 0) then
               if (nch .le. inband-ir_band .or. nch .eq. 26) then
                 cos_ang = bad_data
                 if (nint(solar_zen(iele,Iline)) .ne. nint(bad_data)) 
     +             cos_ang = cos(dtr * solar_zen(iele,Iline))
c ...            read in the vis channel radiances and convert to refs.
c ...            Then make sure you are not diving by zero
                 if (cos_ang .gt. 0) then 
                   rf=rad_band(iele,Iline,nch) / cos_ang
                   indat(iele,mlin,nch) =  rf
c ...            Else, set reflectance holder to -55.0
                 else
                   indat(iele,mlin,nch) = -55.0
                 endif

               else
c ...            read in the thermal radiances and convert to BTs
                 rad  = rad_band(iele,Iline,nch)
                 temperature = modis_bright(rad,nch,units)
                 indat(iele,mlin,nch) = temperature
               end if

c ...       If data is invalid, set radiance holder to -99.0
            else 
               indat(iele,mlin,nch) = -99.0
            endif

 190     continue
 185  continue

c ... debug statement ............................................
      if (debug .gt. 3) then
        write(h_output,'(15x,'' CONTEXT FOR MLIN '', I4,'' AND INDEX '',
     +                   I4)') mlin, Iline
        write(h_output,'(15x,'' Pixels 832 - 841'')')
        write(h_output,'(2x,'' indat 1 (.66um) '',/,8(5f10.6/))')
     +       (((indat(i,j,k),k=1,1),j=mlin,mlin),i=832,841)
        write(h_output,'(2x,'' indat 2 (.87um) '',/,8(5f10.6/))')
     +       (((indat(i,j,k),k=2,2),j=mlin,mlin),i=832,841)
        write(h_output,'(2x,'' v1km_band 1 (.66um) '',/,8(5i10/))')
     +       (((v1km_band(i,j,k),k=1,1),j=Iline,Iline),i=832,841)
        write(h_output,'(2x,'' v1km_band 2 (.87um) '',/,8(5i10/))')
     +       (((v1km_band(i,j,k),k=2,2),j=Iline,Iline),i=832,841)
        write(h_output,'(2x,'' v1km_band 22 (3.9um) '',/,8(5i10/))')
     +       (((v1km_band(i,j,k),k=22,22),j=Iline,Iline),i=832,841)
        write(h_output,'(2x,'' radband 1 (.66um) '',/,8(5f10.6/))')
     +       (((rad_band(i,j,k),k=1,1),j=Iline,Iline),i=832,841)
        write(h_output,'(2x,'' radband 2 (.87um) '',/,8(5f10.6/))')
     +       (((rad_band(i,j,k),k=2,2),j=Iline,Iline),i=832,841)
        write(h_output,'(2x,'' radband 22 (3.9um) '',/,8(5f10.6/))')
     +       (((rad_band(i,j,k),k=22,22),j=Iline,Iline),i=832,841)
        write(h_output,'(2x,'' indat 22 (3.7um) '',/,8(5f10.4/))')
     +       (((indat(i,j,k),k=22,22),j=mlin,mlin),i=832,841)
        write(h_output,'(2x,'' indat 32 (11um) '',/,8(5f10.4/))')
     +       (((indat(i,j,k),k=31,31),j=mlin,mlin),i=832,841)
      endif
c ................................................................

c     NEXT, PUT BANDS NEEDED FOR 250m CLOUD MASK INTO CONTEXT
c     Decided to put 16 values which correspond to 1km footprint
c     into 16 consecutive element array cells.
      offset = 0
      rfoffset = 0
      do 200 nch = 1 , vis_band
         do 205 iele = 1 , maxele * 4
            ie = (iele - 1) / 4 + 1
            cos_vis = cos(dtr * solar_zen(ie,Iline))
            do 210 ilin = 0 , 3 
              rfoffset = (Iline - 1) * 4 + (ilin + 1)
              if (v250_band(iele,rfoffset,nch) .eq. 0) then
                if (cos_vis .gt. 0) then
                  rf = vis250_band(iele,rfoffset,nch) / cos_vis
                else
                  rf = -55.0
                endif
              else 
                rf = -99.0
              endif
              offset = (iele - 1) * 4 + (ilin + 1)
              v250_dat(offset,mlin,nch) = rf
 210        continue
 205     continue
 200  continue  

c ... debug statement ............................................
      if (debug .gt. 3) then
        write(h_output,'(15x,'' CONTEXT FOR 250M PIXELS '', 
     +                   I4)') Iline
        write(h_output,'(2x,'' v250_dat 1(.66um)'',/,5(8f10.6/))')
     +            (((v250_dat(i,j,k),k=1,1),j=mlin,mlin),i=1,40)
        write(h_output,'(2x,'' v250_dat 2(.87um)'',/,5(8f10.6/))')
     +            (((v250_dat(i,j,k),k=2,2),j=mlin,mlin),i=1,40)
      endif
c ................................................................

      return 
      end
