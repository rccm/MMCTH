      subroutine snow_mask(pxldat,plat,land,snglnt,water,hi_elev,
     +                     Greenland,ndsi_snow)

      implicit none
      save

c ... Common statement for debug purposes
      common / bug / debug, h_output

c---------------------------------------------------------------------
c!F77
c
c!Description:
c       Subroutine which implements a quick and dirty version of
c       the snow algorithm by Dorothy Hall, George Riggs and 
c       Vince Salmonson.
c
c!Input Parameters:
c pxldat        Array containing reflectance or brightness temperatures
c               for all bands for a single pixel
c plat          Latitude of current pixel
c land          Flag indicating a pixel to be processed as land
c snglnt        Logical variable flagging background as sunglint 
c               contaminated
c water         Logical variable indicating water surface for current 
c               pixel
c hi_elev       Logical flag indicating elevation > 2000 meters
c Greenland     Logical flag indicating location is "Greenland"
c
c!Output Parameters:
c ndsi-snow     logical variable indicating the prescence of snow in a
c               FOV according to the NDSI test
c
c!Revision History:
c 10/04  Collection 5b  R. Frey
c Removed definition of 'Greenland'.
c Removed tests on 1.38 um, 'sm_mnir', 3.7-11 um BTD in Greenland.
c Added 1.5K to 8.5-11 um BTD test threshold in Greenland.
c
c!Team-unique Header:
c
c!References and Credits:
c See snow product ATBD.
c
c!END
c---------------------------------------------------------------------

      include 'global.inc'
      include 'snow_mask.inc'
      include 'platform_name.inc'
c
c ... scalar arguments ..
      logical snglnt,water,ndsi_snow,hi_elev,land,Greenland
      real plat
c ...
c ... array arguments ..
      real pxldat(inband)
c ...
c ... local scalars ..
      real ndsi,masv88,masv55,masv188,masir13,masir11,
     +     masir85,masir37,diff,sth37_11,masnir,diff2,sth85_11
      integer h_output,debug,I_bad
c
c ... identification 
      masv55 = pxldat(4)
      masv88 = pxldat(2) 
      masv188 = pxldat(26)
      masir37 = pxldat(20)
      masir85 = pxldat(29)
      masir11 = pxldat(31)
      masir13 = pxldat(35)
c     Select NIR channel for NDSI based on platform name.
c     Band 6 for Terra, band 7 for Aqua.
      if(platform_name .eq. 'Terra') then
        masnir = pxldat(6)
      else if(platform_name .eq. 'Aqua ') then
        masnir = pxldat(7)
      else
        call message( 'snow_mask', 'Platform name not recognized' //
     &    ' [OPERATOR ACTION: Contact SDST]', 0, 2)
      endif

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(10x/,''Subroutine Snow_Mask '',/)')
      endif
c ................................................................

c ... Initialize
      ndsi = 0.0
      ndsi_snow = .false.
      I_bad = nint(bad_data)

c ................................................................

c ... First, check the 11 micron brightness temperature.
      if (nint(masir11) .ne. I_bad .and. masir11 .le. sm_bt11(1)) then

c ...    perform the NDSI on the current pixel
         if (nint(masv55) .ne. I_bad .and. nint(masnir).ne.I_bad) then
           ndsi  = (masv55 - masnir) / (masv55 + masnir)

           if ( (ndsi .gt. sm_ndsi(1)) .and. (masv88 .gt. sm_ref2(1)) ) 
     +     then

             ndsi_snow = .true.

c            Check for false snow detection.

c            Now, make sure NDSI is not flagging a thin cirrus
             if( .not. (Greenland .and. land) ) then
               if(nint(masv188).ne.I_bad.and.nint(masir13).ne.I_bad) then
                 if (masv188 .gt. sm_ref3(1) .and. masir13 .lt. sm_co2(1)) then
                   ndsi_snow = .false.
                 endif
               endif
             endif
c            If in sunglint region, disregard if between -60 and 50 lat.
             if(water .and. snglnt .and. plat .le. 50.0 .and.
     +          plat .ge. -60.0) then
               ndsi_snow = .false.
             end if
c            Check for ice clouds mis-identified as snow. Note modified
c            BTD thresholds for Greenland.
             if(masir85 .ne. I_bad .and. masir11 .ne. I_bad) then
               diff = masir85 - masir11
               sth85_11 = sm85_11(1)
               if(masir37 .ne. I_bad .and. (.not. (Greenland .and. land)) ) then
                 diff2 = masir37 - masir11
                 sth37_11 = sm37_11(1)
                 if(diff .ge. sth85_11 .and. diff2 .ge. sth37_11) then
                   ndsi_snow = .false.
                 end if
               else
                 if(Greenland .and. land) then
                   sth85_11 = sm85_11(1) + 1.5
                 else
                   sth85_11 = sm85_11(1)
                 end if
                 if(diff .ge. sth85_11) then
                   ndsi_snow = .false.
                 end if
               end if
             end if
c            Check for water clouds mis-identified as snow.
             if(masir37 .ne. I_bad .and. masir11 .ne. I_bad) then
               if( .not. (land .and. Greenland) ) then
                 diff = masir37 - masir11
                 if(hi_elev) then
                   sth37_11 = sm37_11hel(1)
                 else
                   sth37_11 = sm37_11(1)
                 end if
                 if(diff .ge. sth37_11) then
                   ndsi_snow = .false.
                 end if
               end if
             end if

             if( (.not. hi_elev) .and. (.not. (Greenland .and. land)) ) then
               if(masnir .gt. sm_mnir(1)) then
                 ndsi_snow = .false.
               end if
             end if

           endif 
         endif
      endif

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(15x,'' snow mask variables '',L4)') ndsi_snow
        write(h_output,'(2x,''masv55 masnir ndsi masv188 masv88 masir13
     +      snglnt hi_elev ndsi_snow '',/,6f8.4,f7.2,3L8)') masv55,masnir,ndsi,
     +                   masv188,masv88,sth85_11,masir13,snglnt,hi_elev,ndsi_snow
        write(h_output,'(''IR vals '',3f10.2)') masir85,masir11,masir37
      endif
c ................................................................

      return
      end
