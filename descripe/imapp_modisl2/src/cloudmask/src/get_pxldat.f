      subroutine get_pxldat(indat,lines_in_edge,contx_lat,
     +                      contx_lon,contx_vzen,contx_szen,contx_elev,
     +                      contx_rel,contx_pw,contx_eco,contx_topog,
     +                      contx_ice,contx_snow,contx_sfctmp,contx_msl,
     +                      contx_ugrd,contx_vgrd,contx_sst,
     +                      s_type,klin,
     +                      ibl,nlin,nc,nbands,nbad_1km,pxldat,
     +                      plat,plon,precip_water,vza,sfctmp,pmsl,
     +                      u_wind,v_wind,eco_type,lsf,geo_flag,polar,
     +                      day,night,land,water,coast,snglnt,visusd,
     +                      ice,vrused,desert,snow,bad_value,bad_geo,
     +                      hi_elev,tbadj,antarctic,sh_ocean,sg_bad_data,
     +                      map_ice,map_snow,refang,sh_lake)

      implicit none
      save

c ... Common statement for debug purposes
      common / bug / debug, h_output

c-----------------------------------------------------------------------
c!F77
c
c!DESCRIPTION:  Routine to extract needed processing values for the
c               given pixels out of the context arrays.
c
c!Input Parameters:
c indat         Array containing nlcltx number of lines of
c               reflectance or BT values for each channel
c lines_in_edge Number of lines outside of processing region
c contx_lat     Array containing context of lines of latitude values
c contx_lon     Array containing context of lines of longitude values
c contx_vzen    Array containing context of lines of v. zenith values
c contx_szen    Array containing context of lines of s. zenith values
c contx_elev    Array containing context of lines of elevation values
c contx_rel     Array containing context of lines of rel. angle values
c contx_pw      Array containing context of lines of pw values
c contx_eco     Array containing context of lines of ecosystem values
c contx_topog   Array containing context of lines of land/sea values
c contx_ice     Array containing context of lines of ice fraction values
c contx_snow    Array containing context of lines of snow values
c contx_sfctmp  Array containing context of lines of land temp values
c contx_msl     Array containing context of lines of surface pressure
c contx_ugrd    Array containing context of lines of u wind comp. values
c contx_vgrd    Array containing context of lines of v wind comp. values
c contx_sst     Array containing context of lines of sst values
c s_type        Scan type - Day or night for lines in a context
c klin          Current processing line
c ibl           Beginning line of of data to process
c nlin          Total number of lines to process
c nc            Current processing element
c
c!Output Parameters:
c nbands        Number of bands with good data for this pixel
c nbad_1km      Number of bands with bad data for this pixel (1km)
c pxldat        Array containing reflectance or brightness temperatures
c               for all bands for a single pixel
c plat          Pixel value latitude
c plon          Pixel value longitude 
c precip_water  value of precipitable water (g/cm2) at pixel location
c vza           Viewing zenith angle at pixel location
c sfctmp        Model or sea surface temperate for current pixel
c pmsl          Model mean sea level pressure for current pixel
c u_wind        Model u wind component for current pixel
c v_wind        Model v wind component for current pixel
c eco_type      Ecosystem type for current pixel
c lsf           Land/Sea flag (integer)
c geo_flag      Integer array containing geolocation good/bad flags
c               1-lat,2-lon,3-szen,4-vzen,5-rel_angle
c polar         Logical variable flagging polar scenes
c day           Logical variable flagging day scenes
c night         Logical variable flagging night scenes
c land          Logical variable flagging land scenes
c water         Logical variable flagging water scenes
c coast         Logical variable flagging coast scenes
c snglnt        Logical variable flagging sunglint contaminated scenes
c visusd        Logical variable flagging scenes where visible
c               data were used
c ice           Logical variable flagging ice background scenes
c vrused        Logical variable flagging scenes where reflectance
c               ration test can be used 
c desert        Logical variable flagging desert defined scenes 
c snow          Logical variable flagging snow background scenes
c bad_value     Logical variable flagging bad pixels
c bad_geo       Logical variable flagging bad lat/long data
c hi_elev       Logical variable indicating high elevation (>2000 meters)
c tbadj         11 um brightness temperature adjustment for high deserts
c antarctic     Logical flag indicating region south of -60 lat
c sh_ocean      Logical flag indicating shallow ocean
c sg_bad_data   Logical flag indicating one or more channels needed for
c               sun-glint processing is bad or missing data
c map_ice       Logical flag indicating ice background (from ancillary
c               data)
c map_snow      Logical flag indicating snow background (from ancillary
c               data)
c sh_lake       Logical flag indicating shallow inland lakes
c
c!Revision History:
c 06/04 Collection 5  R. Frey
c 'sfctmp' may contain either SST or LST, added 'sh_lake', added a
c few eco indecis to 'desert' category, added unique Australian desert
c categories, added logic for snow/ice determination on coastlines
c ("NISE=200"), added logic for snow determination on New Zealand
c 10/04 Collection 5  R. Frey
c Moved definition of 'Greenland' from snow_mask.f to get_pxldat.f.
c Changed definiton of 'hi_elev' in Greenland from 1500m to 200m.
c
c!Team-Unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-35.
c
c!END
c---------------------------------------------------------------------

      INCLUDE 'global.inc'

      character*1 s_type(nlcntx)
      real pxldat(inband),
     +     indat(npixel,nlcntx,inband),
     +     contx_szen(npixel,nlcntx),contx_rel(npixel,nlcntx),
     +     contx_ice(npixel,nlcntx),contx_lat(npixel,nlcntx),
     +     contx_pw(npixel,nlcntx),contx_vzen(npixel,nlcntx),
     +     contx_lon(npixel,nlcntx),contx_sfctmp(npixel,nlcntx),
     +     contx_elev(npixel,nlcntx),contx_sst(npixel,nlcntx),
     +     contx_msl(npixel,nlcntx),contx_ugrd(npixel,nlcntx),
     +     contx_vgrd(npixel,nlcntx)
      integer contx_topog(npixel,nlcntx),contx_snow(npixel,nlcntx),
     +        geo_flag(5)
      byte contx_eco(npixel,nlcntx)

      real vza,precip_water,plat,plon,u_wind,v_wind,sfctmp,pmsl,ndvi,sza,
     +     pelev,tbadj,vrat,refang
      integer i,j,k,klin,kk,lines_in_edge,nlin,nc,lsf,nbands,
     +        nbad_1km,debug,h_output,ibl,band(bands_used),
     +        sg_band(sg_bands_used)
      logical polar,land,day,night,ice,snglnt,visusd,water,bad_value,
     +        coast,desert,vrused,snow,bad_geo,map_ice,map_snow,ndsi_snow,
     +        hi_elev,antarctic,sh_ocean,sg_bad_data,sh_lake,
     +        New_Zealand,Greenland
      byte eco_type

c ... MODIS band numbers used for this algorithm
      data band /4,1,2,19,5,26,6,7,20,22,23,27,29,31,32,33,34,35,36/

c ... MODIS band numbers used for sun-glint regions
      data sg_band /1,2,20,26,27,31,32,35/

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(10x/,''Subroutine get_pxldat '',/)')
      endif
c ...............................................................
      
c     Put pixel data into data array.  If working on edge pixels
c     then take either the first or last of 3 in context...otherwise
c     always take the middle line
      if (klin .le. lines_in_edge) then
         i = klin
      else if (klin .gt. nlin-lines_in_edge) then
         i = nlcntx
      else 
         i = ((nlcntx-1) / 2) + 1
      endif 

c     The element number to take is always just the current
c     pixel being processed
      j = nc

c     First get 1km channels - check for bad data
      do 300 k = 1,inband

        if (k .eq. 26) then

          if (indat(j,i,k) .le. 1.3 .and. indat(j,i,k) .gt. -99.0) then
            pxldat(k) = indat(j,i,k) 
          else
            pxldat(k) = bad_data
          end if

        else if (k .le. inband-ir_band) then

          if (indat(j,i,k).gt. 0.0) then
            pxldat(k) = indat(j,i,k) 
          else
            pxldat(k) = bad_data
          end if

        else

          if (indat(j,i,k).gt. 0.0 .and. indat(j,i,k).lt. 1000.0) then
            pxldat(k) = indat(j,i,k) 
          else
            pxldat(k) = bad_data
          endif

        endif

  300 continue

c ... Loop to count number of good bands out of those used
      do 350 k = 1 , bands_used
        if (nint(pxldat(band(k))) .ne. nint(bad_data)) then
          nbands = nbands + 1
        else
          nbad_1km = nbad_1km + 1
          bad_value = .true.
        endif
  350 continue

c ... Check bands needed for sun-glint processing.  If any are
c     missing or bad, set flag.
      do 360 k = 1 , sg_bands_used
        if (nint(pxldat(sg_band(k))) .eq. nint(bad_data)) then
          sg_bad_data = .true.
        endif
  360 continue

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(15x,'' PXLDAT line  '', I4, '' and element ''
     +       ,I4)') klin+ibl-1, nc
        write(h_output,'(4x,''.55    .66    .87    .936   1.2    1.38'',
     + ''    1.6    2.1'', /, 8f7.2)') (pxldat(band(kk)),kk=1,8)
        write(h_output,'(2x,''  3.7    3.9    4.0    6.7    8.5   11.0'',
     + ''   12.0 '', /, 7f7.2)') (pxldat(band(kk)),kk=9,15)
        write(h_output,'(2x,'' 13.3   13.6   13.9   14.2 '',/,4f7.2,/)')
     +   (pxldat(band(kk)),kk=16,19)
      endif
c ...............................................................

c     get pixel value of precipitable water
      precip_water = contx_pw(j,i)

c     get pixel value of viewing zenith angle 
      vza = contx_vzen(j,i)
      if (vza .lt. 0.0) then
         geo_flag(4) = 1
      endif

c     Get pixel value of reflectance angle.
      refang = contx_rel(j,i)

c     get pixel value of model mean sea level pressure
      pmsl = contx_msl(j,i)

c     get pixel value of model u component of the surface wind
      u_wind = contx_ugrd(j,i)

c     get pixel value of model u component of the surface wind
      v_wind = contx_vgrd(j,i)

c     get pixel latitude and longitude
      plat = contx_lat(j,i)
      plon = contx_lon(j,i)
      if (nint(plat) .eq. nint(bad_data)) then
         geo_flag(1) = 1
         bad_geo = .true.
      endif
      if (nint(plon) .eq. nint(bad_data)) then
         geo_flag(2) = 1
         bad_geo = .true.
      endif

c     get pixel elevation
      pelev = contx_elev(j,i)

c     Get 11 um brightness temperature threshold adjustment for
c     deserts.
      tbadj = (pelev / 1000.0) * 5.0

c     get numeric value of land/sea flag
      lsf = contx_topog(j,i)

c     get ecosystem type
      eco_type = contx_eco(j,i)

c     Define "Antarctic" flags.
      if(plat .lt. -60.0) then
        antarctic = .true.
      end if

c     Determine whether or not current pixel will be processed as
c     desert.

c     Global.
      if(eco_type .eq. 8 .or. eco_type .eq. 46 .or.
     +   eco_type .eq. 50 .or. eco_type .eq. 51 .or.
     +   eco_type .eq. 59 .or. eco_type .eq. 71 .or.
     +   eco_type .eq. 11 .or. eco_type .eq. 9 .or. 
     +   eco_type .eq. 52) then
        desert = .true.

c     High-elevation.
      else if( (pelev .gt. 2000.0) .and. (eco_type .eq. 42) ) then
        desert = .true.
c       Check for locations where there should be no high-elevation desert.
        if( (plat .le. 10.0 .and. plat .ge. -10.0) .and.
     +      (plon .ge. 90.0) ) then
          desert = .false.
        else if( (plat .ge. -30.0 .and. plat .le. -10.0) .and.
     +           (plon .ge. 160.0 .and. plon .le. 180.0) ) then
          desert = .false.
        else if( (plat .ge. 10.0 .and. plat .le. 26.0) .and.
     +           (plon .ge. 120.0 .and. plon .le. 180.0) ) then
          desert = .false.
        end if

c     Africa.
      else if(plat .le. 20.0 .and. plat .ge. -35.0 .and. plon .le. 60.0
     +        .and. plon .ge. -20.0) then
        if(eco_type .eq.  7 .or. eco_type .eq. 41 .or.
     +     eco_type .eq. 43 .or. eco_type .eq. 58 .or.
     +     eco_type .eq. 36 .or. eco_type .eq. 91 .or.
     +     eco_type .eq. 32 .or. eco_type .eq. 29) then
           desert = .true.
        end if

c     Eurasia.
      else if(plat .le. 70.0 .and. plat .ge. -60.0 .and. plon .le. 180.0
     +        .and. plon .ge. -20.0) then
        if(eco_type .eq. 11 .or. eco_type .eq. 2) then
           desert = .true.
        end if
c       Exclude New Zealand.
        if(plat .ge. -50.0 .and. plat .le. -30.0 .and. plon .ge. 160.0
     +     .and. plon .le. 180.0) then
          desert = .false.
        end if

      end if
c     Add in Australia.
      if(plat .le. -11.0 .and. plat .gt. -40.0 .and. plon .le. 155.0
     +        .and. plon .gt. 110.0) then
        if(eco_type .eq. 43 .or. eco_type .eq. 41 .or. 
     +     eco_type .eq. 91) then
          desert = .true.
        end if
      end if

c     Determine whether or not visible ratio test may be used over
c     land surfaces.
      if(eco_type .eq. 2 .or. eco_type .eq. 8 .or.
     *   eco_type .eq. 11 .or. eco_type .eq. 40 .or.
     *   eco_type .eq. 41 .or. eco_type .eq. 46 .or.
     *   eco_type .eq. 51 .or. eco_type .eq. 52 .or.
     *   eco_type .eq. 59 .or. eco_type .eq. 71 .or.
     *   eco_type .eq. 50) then
        vrused = .false.
      endif

c     Now we turn our attention to the ancillary data sets
c     get logical flag variables for given pixel
c     Day/night flag (current definition is > 85 degrees)
      sza = contx_szen(j,i) 
      if (sza .lt. 0.0) then
        bad_geo = .true.
        geo_flag(3) = 1
      else if (sza .gt. 85.0 .or. s_type(i) .ne. 'D') then
         night = .true.
      else
         day = .true.
      endif

c     set polar flag (if lat is higher then +/-60)
      if (contx_lat(j,i) .gt. 60.0 .or. contx_lat(j,i) .lt. -60.0) then
         polar = .true.
      endif

c     set visusd flag (currently set to true if szen < 85 degrees)
      if (contx_szen(j,i) .gt. 85.0 .and. s_type(i) .ne. 'D') then
         visusd = .false.
      endif

c     set the sunglnt flag
      if (nint(contx_rel(j,i)) .eq. 999.0) then
         geo_flag(5) = 1
      elseif (contx_rel(j,i) .le. 36.0) then
         snglnt = .true.
      endif

c     set land/sea flag
c ... First make sure that it is not a missing value
      if (lsf .ne. -1) then
        if (lsf .eq. 1 .or. lsf .eq. 4) then
          land = .true.
          if(eco_type .eq. 14) then
c           Fix-up for missing ecosystem data in eastern Greenland and
c           north-eastern Siberia.  Without this, these regions become
c           completely "coast".
            if( (plat .lt. 64.0) .or. (plat .ge. 67.5 .and. ((plon
     *                .lt. -40.0 .and. plon .gt. -168.5) .or.
     *                                 plon .gt. -12.5)) .or.
     *          ((plat .ge. 64.0 .and. plat .lt. 67.5) .and.
     *          ((plon .lt. -40.0 .and. plon .gt. -168.5) .or.
     *                                  plon .gt. -30.0)) ) then
              coast = .true.
            end if
          end if
        elseif (lsf .eq. 2) then
          coast = .true.
          land = .true.
        elseif (lsf .eq. 3) then
          land = .true.
          sh_lake = .true.
c         Need shallow lakes to be processed as "coast" for day, but
c         not night.
          if(day) then
            coast = .true.
          end if
        else 
          water = .true.
          if(lsf .eq. 0) then
            sh_ocean = .true.
          end if
        endif

c     If land/sea flag is missing, then calculate visible ratio to 
c     determine if land or water.

      elseif (nint(pxldat(1)) .ne. nint(bad_data) .and.    
     +  nint(pxldat(2)) .ne. nint(bad_data)) then 
        vrat = pxldat(2) / pxldat(1)
        if (vrat .gt. 0.9) then 
          land = .true.
        else
          water = .true.
        endif

c ... If all else fails, call it land.

      else
        land = .true.
        water = .false.
      endif

c     Get pixel value of SST or model surface temp
      if(land) then
        sfctmp = contx_sfctmp(j,i)
      else
        sfctmp = contx_sst(j,i)
      end if

c     Set high elevation flag.
c     First, define "Greenland".
      Greenland = .false.
      if(land) then
        if(plat .ge. 60.0 .and. plat .lt. 67.0) then
          if(plon .ge. -60.0 .and. plon .lt. -30.0) then
            Greenland = .true.
          end if
        else if(plat .ge. 67.0 .and. plat .lt. 75.0) then
          if(plon .ge. -60.0 .and. plon .lt. -10.0) then
            Greenland = .true.
          end if
        else if(plat .ge. 75.0) then
          if(plon .ge. -70.0 .and. plon .lt. -10.0) then
            Greenland = .true.
          end if
        end if
      end if
      if( (pelev .gt. 2000.0) .or. (pelev .gt. 200.0 .and. Greenland
     +    .and. land) .or. (plat .ge. 75.7 .and. plat .le. 79.0
     +    .and. plon .ge. -73.0 .and. plon .le. -50.0 .and. land) ) then
        hi_elev = .true.
      end if

c ... Calculate raw NDVI value.
      if(nint(pxldat(1)) .ne. nint(bad_data) .and. nint(pxldat(2))
     +                   .ne. nint(bad_data)) then
        ndvi = (pxldat(2)-pxldat(1)) / (pxldat(2)+pxldat(1))
      end if

c ... Set the ancillary ice flag.
      if ( (contx_ice(j,i) .gt. 0.5 .and. contx_ice(j,i) .le. 1.0) 
     +    .or. (contx_snow(j,i) .gt. 25 .and. contx_snow(j,i) .lt. 102
     +    .and. water) ) then
        map_ice = .true.
      endif
      if(contx_snow(j,i) .eq. 200 .and. water .and. plat .lt. -60.0)
     +          map_ice = .true.

c ... Set the ancillary snow flag.
      if (contx_snow(j,i) .eq. 103 .or. contx_snow(j,i) .eq. 104 .or.
     +    contx_snow(j,i) .eq. 101) then
        map_snow = .true.
      endif
      if(contx_snow(j,i) .gt. 25 .and. contx_snow(j,i) .lt. 101 .and.
     +   land) map_snow = .true.
      if( (contx_snow(j,i) .eq. 200 .and. land) .and. (plat .ge. 44.0
     +  .or. plat .le. -40.0) ) map_snow = .true.

c ****************************************************************

      if(day) then

c       Run quick version of D.Hall's snow detection algorithm.
        call snow_mask(pxldat,plat,land,snglnt,water,hi_elev,
     +                 Greenland,ndsi_snow)

        if(water) then

          if(plat .ge. -60.0 .and. plat .le. 25.0) then
            if(map_ice .and. ndsi_snow) then
              ice = .true.
            end if
          else if( (contx_snow(j,i).eq.252 .or. contx_snow(j,i).eq.200)
     +        .and. (plat .ge. 44.0 .or. plat .le. -40.0) ) then
            if(ndsi_snow) then
              ice = .true.
            end if
          else if((lsf .eq. 3 .or. lsf .eq. 5) .and. ndsi_snow) then
            ice = .true.
          else if(map_ice .and. ndsi_snow) then
            ice = .true.
          endif

        else if(land) then

          if(plat .ge. -60.0 .and. plat .le. 25.0) then
c           Define New Zealand region which receives snow but snow
c           map does not show it.
            if( (plat .ge. -48.0 .and. plat .le. -34.0) .and.
     +          (plon .ge. 165.0 .and. plon .le. 180.0) ) then
              New_Zealand = .true.
            else
              New_Zealand = .false.
            end if
            if(map_snow .and. ndsi_snow) then
              snow = .true.
            else if(New_Zealand .and. ndsi_snow) then
              snow = .true.
            end if
          else if(plat .lt. -60.0) then
            snow = .true.
          else
            if(ndsi_snow) then
              snow = .true.
            end if
          endif

        endif

      else

        ice = map_ice
        snow = map_snow

      endif

c ****************************************************************

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(15x,'' PXLDAT variables '',/)')
        write(h_output,'(5x,'' # bad bands, bad_value, bad_geo '',
     +        i10,2l4/)') nbad_1km,bad_value,bad_geo
        write(h_output,'(15x,'' ncep ice, nsidc ice '',f6.3,i5,3l5/)')
     +        contx_ice(j,i),contx_snow(j,i),map_ice,map_snow,ndsi_snow
        write(h_output,'(2x,'' Day Polar Desert Visusd Snglnt Land Water
     + Coast Ice'',/,3L5,2x,12L6)') day,polar,desert,visusd,snglnt,land,
     + water,coast,map_ice,map_snow,ice,snow,hi_elev
        write(h_output,'(2x,''  pw     vza    sfctmp  pmsl    u_wind  v_
     +wind eco_type lsf'',/,7f8.2,2i5)') precip_water,vza,sfctmp,pmsl,
     + u_wind,v_wind,ndvi,eco_type,lsf
      endif
c ...............................................................

      
      return
      end
