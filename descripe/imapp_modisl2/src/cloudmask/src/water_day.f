      subroutine water_day(pxldat,vza,snglnt,visusd,refang,cirrus_vis,
     *                     sfctmp,hi_elev,uniform,ice,maxele,klin,
     *                     line_edge,sh_ocean,indat,kele,nmtests,
     *                     testbits,qa_bits,confdnc)

      implicit none
      save

c---------------------------------------------------------------------
C!F77 
c
c!Description:
c     Routine for setting appropriate flags and processing path
c     for daytime observations over water surfaces.
c
c     If the confidence determined from spectral tests is
c     uncertain, then other tests may be applied. 
c
c!Input parameters:
c pxldat        Array containing reflectance or brightness temperatures
c               for all bands for a single pixel
c vza           Current pixel viewing angle
c snglnt        Logical variable flagging sunglint pixels
c visusd        Logical variable indicating whether vis data used or not
c refang        Reflectance angle for current pixel
c cirrus_vis    Logical variable flagging thin cirrus contaminated
c               scenes in the visible
c sfctmp        SST from ancillary data
c hi_elev       Logical variable indicating elevations > 2000 meters
c uniform       Logical variable indicating uniformity of context
c ice           Logical variable flagging ice backgrounds
c sh_ocean      Logical variable indicating shallow ocean
c indat         Array containing nlcntx lines of data
c kele          Current granule element number being processed
c
c!Output Parameters:
c nmtests       Number of tests applied to this pixel
c testbits      Byte array containing cloud mask results
c qa_bits       Byte array containing qa bit results
c confdnc       Current pixel unobstructed confidence
c
c!Revision History:
c 06/04 Collection 5  R. Frey
c Updated argument list to Day_snow
c 10/04 Collection 5b R. Frey
c Updates to calling arguments lists.
c
c!Team-unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-35.
c
c!END
c---------------------------------------------------------------------

      include 'global.inc'

c     scalar arguments
      real vza,confdnc,refang,sfctmp
      integer kele,nmtests,maxele,klin
      logical visusd,snglnt,uniform,ice,cirrus_vis,hi_elev,sh_ocean,
     *        line_edge

c     array arguments
      real pxldat(inband),indat(npixel,nlcntx,inband)
      byte testbits(6),qa_bits(10)

c     local scalars
      integer debug,h_output

c     external subroutines
      external ocean_day,Day_snow,chk_spatial_var,chk_sunglint,
     *         chk_shallow_water

c     Common statement for debug purposes
      common / bug / debug, h_output

c---------------------------------------------------------------------

c     debug statement 
      if (debug .gt. 0) then
        write(h_output,'(10x/,''Using water_day processing path '',/)')
      endif

c---------------------------------------------------------------------

c     Further define processing path.

      if (ice) then
c        Processing for ice-covered scenes.

         call Day_snow(pxldat,vza,visusd,cirrus_vis,hi_elev,testbits,
     +                 qa_bits,nmtests,confdnc)

      else 
  
c        Normal ocean processing
         call ocean_day(pxldat,vza,snglnt,visusd,cirrus_vis,sfctmp,
     +                  refang,sh_ocean,testbits,qa_bits,nmtests,
     +                  confdnc)

      endif

c---------------------------------------------------------------------

c     debug statement 
      if (debug .gt. 0) then
        write(h_output,'(10x/,'' Water day premlinary confidence '',
     +        f10.2)') confdnc
        write(h_output,'(10x/,'' uniform,snglnt,sh_ocean,ice '',
     +        4l5/)') uniform,snglnt,sh_ocean,ice
      endif

c---------------------------------------------------------------------

c     Perform clear sky confidence confirmation tests.

c     Under certain conditions, apply the spatial variability test.
      if(confdnc .le. 0.95 .and. confdnc .gt. 0.05 .and. uniform) then
         call chk_spatial_var(indat,kele,confdnc,qa_bits,testbits)
      end if

c     Perform clear sky restoral tests in sun-glint regions.
      if(snglnt .and. uniform .and. confdnc.lt.0.95) then
         call chk_sunglint(indat,pxldat,kele,confdnc,maxele,klin,
     +                     line_edge,qa_bits,testbits)
      end if

c     Perform clear sky restoral tests in shallow ocean conditions.
      if(sh_ocean .and. (.not. ice)) then
         call chk_shallow_water(pxldat,confdnc,qa_bits,testbits)
      end if

c---------------------------------------------------------------------

      return
      end
