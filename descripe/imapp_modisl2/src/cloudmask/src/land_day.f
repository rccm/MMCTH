      subroutine land_day(pxldat,vza,visusd,vrused,cirrus_vis,
     *                    desert,coast,snow,ice,hi_elev,tbadj,eco_type,
     *                    testbits,qa_bits,nmtests,confdnc)

      implicit none
      save

c-----------------------------------------------------------------------
C!F77 
c
c!Description:
c     Routine for setting appropriate flags and processing path
c     for daytime observations over land surfaces.
c
c!Input Parameters:
c pxldat        Array containing reflectance or brightness temperatures
c               for all bands for a single pixel
c vza           Current pixel viewing angle
c visusd        Logical variable indicating whether vis data used or not
c vrused        Logical variable indicating when vrat test can be used
c cirrus_vis    Logical variable flagging thin cirrus contaminated
c               scenes in the visible
c hi_elev       Logical variable indicating high elevation (> 2000 meters)
c tbadj         11 um brightness temperature threshold adjustment for 
c               deserts
c desert        Logical variable indicating desert ecosystems
c coast         Logical variable indicating coast ecosystems
c snow          Logical variable flagging snow backgrounds
c ice           Logical variable flagging ice backgrounds
c eco_type      Byte variable containing ecosystem index for current pixel
c
c!Output Parameters:
c testbits      Byte array containing cloud mask results
c qa_bits       Byte array containing qa bit results
c nmtests       Number of tests actually applied to the given pixel
c confdnc       Current pixel unobstructed confidence
c
c!Revision History:
c 06/04 Collection 5  R. Frey
c Updated calling arguments to Day_snow.
c
c!Team-Unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-35.
c
c!END
c-----------------------------------------------------------------------

      include 'global.inc'

c     scalar arguments
      real vza,confdnc,tbadj
      integer nmtests
      logical visusd,vrused,snow,ice,desert,coast,cirrus_vis,hi_elev
      byte eco_type

c     array arguments
      real pxldat(inband)
      byte testbits(6),qa_bits(10)

c     local scalars
      integer h_output,debug

c     external subroutines
      external LandDay_desert,LandDay,Day_snow,LandDay_desert_c,
     +         LandDay_coast,chk_land,chk_coast

c     Common statement for debug purposes
      common / bug / debug, h_output

c-----------------------------------------------------------------------

c     debug statement 
      if (debug .gt. 0) then
        write(h_output,'(10x/,''Using land_day processing path '',/)')
      endif

c-----------------------------------------------------------------------

c     Land processing further defined based on ecosystem map, land/sea
c     tag and snow/ice data.

      if (snow .or. ice) then

         call Day_snow(pxldat,vza,visusd,cirrus_vis,hi_elev,
     +                 testbits,qa_bits,nmtests,confdnc)

      else if (desert .and. coast) then

         call LandDay_desert_c(pxldat,vza,visusd,cirrus_vis,
     +                         hi_elev,testbits,qa_bits,
     +                         nmtests,confdnc)

      else if (coast) then

         call LandDay_coast(pxldat,vza,visusd,cirrus_vis,
     +                      hi_elev,testbits,qa_bits,nmtests,
     +                      confdnc)
 
      else if (desert) then


         call LandDay_desert(pxldat,vza,visusd,cirrus_vis,
     +                       hi_elev,testbits,qa_bits,nmtests,
     +                       confdnc)
 
      else
 
         call LandDay(pxldat,vza,visusd,vrused,cirrus_vis,
     +                hi_elev,testbits,qa_bits,nmtests,confdnc)

      endif

c-----------------------------------------------------------------------

c     Perform clear sky confidence confirmation tests.

c     Check non-snow covered land regions.
      if(.not. (snow .or. ice)) then
       if(confdnc .le. 0.95) then
         call chk_land(pxldat,eco_type,desert,tbadj,confdnc,qa_bits,
     +                 testbits)
       end if
      end if

c     Check coastal regions.
      if(coast .and. (.not. (snow .or. ice))) then
        call chk_coast(pxldat,confdnc,qa_bits,testbits)
      end if
       
c-----------------------------------------------------------------------

      return
      end
