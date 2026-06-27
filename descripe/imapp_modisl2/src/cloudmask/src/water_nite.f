      subroutine water_nite(pxldat,vza,uniform,ice,indat,sfctmp,sh_ocean,
     *                      kele,nmtests,testbits,qa_bits,confdnc)

      implicit none
      save

c---------------------------------------------------------------------
C!F77 
c
c!Description:
c     Routine for setting appropriate flags and processing path
c     for nighttime observations over water surfaces.
c
c     If the confidence determined from the water background is
c     uncertain, than a spatial variability test is applied if
c     the pixels within the processing context are all from the
c     same ecosystem (not including a coastline).
c
c!Input parameters:
c pxldat        Array containing reflectance or brightness temperatures
c               for all bands for a single pixel
c uniform       Logical variable flagging context of similar ecosystem
c ice           Locical variable flagging ice backgrounds
c indat         Array containing nlcntx lines of data
c sfctmp        Surface temperature for current pixel
c sh_ocean      Logical flag indicating ocean depths < 50 m or within
c               5 km of shorlines
c kele          Current granule element number being processed
c vza           Viewing zenith angle
c
c!Output Parameters:
c nmtests       Number of tests applied for this pixel
c testbits      Byte array containing cloud mask results
c qa_bits       10-byte array containing qa bit results
c confdnc       Current pixel unobstructed confidence
c
c!Revision History:
c Removed 11 micron standard deviation calculation
c Added 'lnd' and 'sfctmp' for Collection 5 processing.
c 10/04  Collection 5b  R. Frey
c Updated calling argument list for ocean_nite.
c
c!Team-Unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-35.
c
c!END
c---------------------------------------------------------------------
c
      include 'global.inc'
    
c     scalar arguments
      real vza,confdnc,sfctmp
      integer kele,nmtests
      logical uniform,ice,sh_ocean
c
c     array arguments
      real pxldat(inband),indat(npixel,nlcntx,inband)
      byte testbits(6),qa_bits(10)

c     local scalars
      integer debug,h_output
      logical lnd

c     local arrays

c     external subroutines
      external ocean_nite,Nite_snow,chk_spatial_var

c     Common statement for debug purposes
      common / bug / debug, h_output

c---------------------------------------------------------------------
c
c     Debug statement.
      if (debug .gt. 0) then
        write(h_output,'(10x/,''Using water_nite processing path '',/)')
      endif

c---------------------------------------------------------------------

      if (ice) then

c        Ice processing

         lnd = .false.
         call Nite_snow(pxldat,vza,lnd,testbits,qa_bits,nmtests,
     +                  confdnc)

      else

c        Normal nighttime ocean processing

         call ocean_nite(indat,kele,pxldat,vza,sfctmp,sh_ocean,
     +                   uniform,testbits,qa_bits,nmtests,confdnc)

      endif

c---------------------------------------------------------------------

c     Debug statement.
      if (debug .gt. 0) then
        write(h_output,'(10x/,'' Water nite preliminary confidence '',
     +        f10.2,/)') confdnc
      endif

c---------------------------------------------------------------------

c     Perform clear sky confidence confirmation tests.

      if(confdnc .le. 0.95 .and. confdnc .gt. 0.05 .and. uniform) then
        call chk_spatial_var(indat,kele,confdnc,qa_bits,testbits)
      end if

c---------------------------------------------------------------------

      return
      end
