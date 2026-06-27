      subroutine polar_nite(pxldat,vza,land,ice,snow,desert,hi_elev,
     +                      sfctmp,eco_type,nmtests,testbits,uniform,
     +                      indat,kele,antarctic,sh_ocean,qa_bits,confdnc)

      implicit none
      save

c----------------------------------------------------------------------
C!F77 
c
c!Description:
C     Routine for providing conditional input parameters pertaining
C     to polar nighttime processing.
C
c!Input parameters:
c pxldat        Array containing reflectance or brightness temperatures
c               for all bands for a single pixel
c vza           Viewing zenith angle
c land          Logical variable flagging land backgrounds
c ice           Logical variable flagging ice backgrounds
c snow          Logical variable flagging snow backgrounds
c desert        Logical variable flagging desert backgrounds
c hi_elev       Logical variable flagging high elevations (> 2000 meters)
c sfctmp        SST from Reynolds Blended data set
c eco_type      Ecosystem index from Olson ecosystem data set
c uniform       Logical variable flagging uniform background
c kele          Current element number being processed
c indat         Array containing 'nlcntx' lines of radiance data
c antarctic     Logical variable indicating pixel is south of -60 lat
c sh_ocean      Logical variable indicating ocean depths < 50 m or within
c               5 km of shoreline
c
c!Output Parameters:
c nmtests       Number of tests applied to this pixel
c testbits      Byte array containing cloud mask results
c qa_bits       Byte array containing qa bit results
c confdnc       Current pixel unobstructed confidence
c
c!Revision History:
c 06/04 Collection 5  R. Frey
c Removed 11 micron standard deviation calculation 
c Added 'desert','hi_elev','antarctic','eco_type','vza','sfctmp','indat'
c   needed for Collection 5 processing
c 10/04 Collection 5b R. Frey
c Updated calling arguments for PolarNite_ocean.
c
c!Team-Unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c
c!END
c----------------------------------------------------------------------

      include 'global.inc'

c     scalar arguments
      real confdnc,sfctmp
      logical land,snow,ice,uniform,desert,hi_elev,antarctic,sh_ocean
      integer nmtests,kele
      byte eco_type

c     array arguments
      real pxldat(inband),vza,indat(npixel,nlcntx,inband)
      byte testbits(6),qa_bits(10)

c     local scalars
      integer debug,h_output

c     local arrays

c     external subroutines
      external PolarNite_snow,PolarNite_land,PolarNite_ocean,
     *         chk_spatial_var

c     Common statement for debug purposes
      common / bug / debug, h_output

c----------------------------------------------------------------------

c     Debug statement.
      if (debug .gt. 0) then
        write(h_output,'(10x/,''Using polar_nite processing path '',/)')
      endif

c----------------------------------------------------------------------

      if (snow .or. ice) then

c        snow or ice covered surfaces
         call PolarNite_snow(pxldat,vza,hi_elev,land,antarctic,testbits,
     +                       qa_bits,nmtests,confdnc)

      else if (land) then

c        Land surfaces
         call PolarNite_land(pxldat,vza,desert,hi_elev,sfctmp,
     +                       eco_type,testbits,qa_bits,nmtests,confdnc)

      else

c        Water surfaces.
         call PolarNite_ocean(indat,kele,pxldat,vza,sfctmp,sh_ocean,
     +                        uniform,testbits,qa_bits,nmtests,confdnc)

      end if

c----------------------------------------------------------------------

c     Debug statement.
      if (debug .gt. 0) then
        write(h_output,'(10x/,'' Polar Nite preliminary confidence '',
     +        f10.2,/)') confdnc
      endif

c----------------------------------------------------------------------

c     Perform clear sky confidence confirmation tests. 

      if(confdnc .le. 0.95 .and. confdnc .gt. 0.05 .and. uniform .and.
     *  (.not. land)) then
        call chk_spatial_var(indat,kele,confdnc,qa_bits,testbits)
      end if

c----------------------------------------------------------------------

      return
      end
