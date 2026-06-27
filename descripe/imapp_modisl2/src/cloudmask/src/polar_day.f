      subroutine polar_day(pxldat,vza,snglnt,visusd,refang,vrused,
     *                     cirrus_vis,land,ice,snow,desert,coast,eco_type,
     *                     uniform,hi_elev,kele,indat,nmtests,testbits,tbadj,
     *                     antarctic,sh_ocean,sfctmp,qa_bits,confdnc)

      implicit none
      save

C----------------------------------------------------------------------
C!F77 
C
C!Description:
C     Routine for providing conditional input parameter pertaining to 
C     polar daytime processing.
C
c!Input Parameters:
c pxldat        Array containing reflectance or brightness temperatures
c               for all bands for a single pixel
c vza           Current pixel viewing angle
c snglnt        Logical variable flagging sunglint pixels
c visusd        Logical variable indicating whether vis data used or not
c visusd        Logical variable indicating when reflect ratio test
c               can be implemented
c refang        reflectance angle for current pixel
c cirrus_vis    Logical variable flagging thin cirrus contaminated
c               scenes in the visible
c land          Logical variable flagging land backgrounds
c ice           Logical variable flagging ice backgrounds
c snow          Logical variable flagging snow backgrounds
c desert        Logical variable flagging desert backgrounds
c coast         Logical variable flagging coast backgrounds
c uniform       Logical variable flagging contexts with similar surface properties
c hi_elev       Logical flag indicating elevations > 2000 meters.
c kele          Current pixel being processed
c indat         array containing 'nlcntx' lines of data
c tbadj         11 um brightness temperature threshold adjustment for 
c               deserts (based on elevation)
c eco_type      Byte variable containing ecosystem index for current pixel
c antarctic     Logical flag indicating Antarctic region (< -60 latitude)
c sh_ocean      Logical variable indicating shallow ocean
c sfctmp        SST for current pixel
c
c!Output Parameters:
c nmtests       Number of tests applied to this pixel
c testbits      Byte array containing cloud mask results
c qa_bits       Byte array containing qa bit results
c confdnc       Current pixel unobstructed confidence
c
c!Revision History:
c 06/04 Collection 5  R. Frey
c Updated calling arguments to PolarDay_snow and Antarctic_day.
c 10/04 Collection 5b R. Frey
c Updated calling arguments to PolarDay_ocean and Antarctic_day.
c
c!Team-Unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c
c!END
C----------------------------------------------------------------------

      include 'global.inc'

c     scalar arguments
      real vza,confdnc,tbadj,refang,sfctmp
      logical visusd,land,snow,snglnt,ice,coast,desert,vrused,cirrus_vis,
     *        uniform,hi_elev,antarctic,sh_ocean
      integer nmtests,kele
      byte eco_type

c     array arguments
      real pxldat(inband),indat(npixel,nlcntx,inband)
      byte testbits(6),qa_bits(10)

c     local scalars
      integer debug,h_output

c     external subroutines
      external PolarDay_snow,PolarDay_land,PolarDay_ocean,PolarDay_coast,
     *         PolarDay_desert,PolarDay_desert_c,chk_land,chk_coast,
     *         chk_spatial_var,chk_shallow_water

c     Common statement for debug purposes
      common / bug / debug, h_output

c----------------------------------------------------------------------

c     Debug statement 
      if (debug .gt. 0) then
        write(h_output,'(10x/,''Using polar_day processing path'',/)')
      endif

c----------------------------------------------------------------------

c     Polar processing further defined based on ecosystem map, land/sea
c     tag file, and snow/ice data.

      if(antarctic .and. land) then
c        Antarctica
         call Antarctic_day(pxldat,visusd,testbits,qa_bits,
     *                      nmtests,confdnc)

      else if(snow .or. ice) then
c        snow or ice covered surfaces

         call PolarDay_snow(pxldat,vza,visusd,cirrus_vis,hi_elev,
     *                      testbits,qa_bits,nmtests,confdnc)

      else if (land) then
c        Land surfaces

         if (desert .and. coast) then

            call PolarDay_desert_c(pxldat,vza,visusd,cirrus_vis,hi_elev,
     *                             testbits,qa_bits,nmtests,confdnc)

         else if (coast) then

            call PolarDay_coast(pxldat,vza,visusd,cirrus_vis,hi_elev,
     *                          testbits,qa_bits,nmtests,confdnc)

         else if (desert) then

            call PolarDay_desert(pxldat,vza,visusd,cirrus_vis,hi_elev,
     *                           testbits,qa_bits,nmtests,confdnc)

         else 

            call PolarDay_land(pxldat,vza,visusd,vrused,cirrus_vis,
     *                         hi_elev,testbits,qa_bits,nmtests,
     *                         confdnc)

         endif

      else

c        water surface

         call PolarDay_ocean(pxldat,vza,snglnt,visusd,cirrus_vis,
     *                       refang,sfctmp,sh_ocean,testbits,qa_bits,
     *                       nmtests,confdnc)

      end if

c----------------------------------------------------------------------

c     Debug statement 
      if (debug .gt. 0) then
        write(h_output,'(10x/,'' Polar day premlinary confidence '',
     +        f10.2,3l5/)') confdnc,land,coast,sh_ocean
      endif

c----------------------------------------------------------------------

c     Perform clear sky confidence confirmation tests.

c     Under certain conditions, apply spatial variability test.
      if(confdnc .le. 0.95 .and. confdnc .gt. 0.05 .and. uniform .and.
     *   (.not. land)) then
        call chk_spatial_var(indat,kele,confdnc,qa_bits,testbits)
      end if

c     Check land.
      if(land .and. (.not. (snow .or. ice))) then
        if(confdnc .le. 0.95) then
          call chk_land(pxldat,eco_type,desert,tbadj,confdnc,qa_bits,
     *                  testbits)
        end if
      end if

c     Check coastal regions.
      if(coast .and. (.not. (snow .or. ice))) then
        call chk_coast(pxldat,confdnc,qa_bits,testbits)
      end if

c     Perform "final test" in shallow ocean conditions.
      if(sh_ocean .and. (.not. ice)) then
         call chk_shallow_water(pxldat,confdnc,qa_bits,testbits)
      end if

c-----------------------------------------------------------------------

      return
      end
