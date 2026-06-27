      subroutine land_nite(pxldat,plat,vza,ice,snow,coast,tbadj,desert,
     +                     hi_elev,sh_lake,sfctmp,eco_type,nmtests,
     +                     testbits,qa_bits,confdnc)

      implicit none
      save

c----------------------------------------------------------------------
c!F77 
c
c!Description:
c     Routine for setting appropriate flags and processing path
c     for nighttime observations over land surfaces.
c
c!Input Parameters:
c pxldat        Array containing reflectance or brightness temperatures
c               for all bands for a single pixel
c vza           viewing zenith angle
c ice           Logical variable flagging ice backgrounds
c snow          Logical variable flagging snow backgrounds
c coast         Logical variable indicating coast processing
c tbadj         11 um brightness temperature elevation adjustment
c desert        Logical variable flagging arid ecosystem
c hi_elev       Logical variable flagging high elevation (> 2000 m)
c sh_lake       Logical flag indicating shallow inland lakes
c sfctmp        Surface air temperature from model output 
c eco_type      Ecosystem index
c
c!Output Parameters:
c nmtests       Number of tests applied to this pix
c testbits      Byte array containing cloud mask results
c qa_bits       Byte array containing qa bit results
c confdnc       Current pixel unobstructed confidence
c
c!Revision History:
c 06/04 Collection 5  R. Frey
c Updated calling argument lists for Nite_snow and LandNite.
c Added 'desert','hi_elev','sfctmp','eco_type','lnd', removed 'plat'.
c
c!Team-unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-35.
c
c!END
c----------------------------------------------------------------------

      include 'global.inc'

c     scalar arguments
      real confdnc,vza,tbadj,plat,sfctmp
      logical snow,ice,coast,desert,hi_elev,sh_lake
      integer nmtests
      byte eco_type

c     array arguments
      real pxldat(inband)
      byte testbits(6),qa_bits(10)

c     local scalars
      integer debug,h_output
      logical lnd

c     external subroutines
      external LandNite,Nite_snow,chk_land_nite

c     Common statement for debug purposes
      common / bug / debug, h_output

c----------------------------------------------------------------------

c     Debug statement.
      if (debug .gt. 0) then
        write(h_output,'(10x/,''Using land_nite processing path '',/)')
      endif

c----------------------------------------------------------------------

      if (snow .or. ice) then
c        snow or ice covered surfaces

         lnd = .true.
         call Nite_snow(pxldat,vza,lnd,testbits,qa_bits,nmtests,confdnc)
      
      else
c        Standard processing
        
         call LandNite(pxldat,plat,vza,coast,desert,hi_elev,sh_lake,
     +                 sfctmp,eco_type,testbits,qa_bits,nmtests,
     +                 confdnc)
 
      endif

c----------------------------------------------------------------------

c     Debug statement.
      if (debug .gt. 0) then
        write(h_output,'(10x/,'' Land nite confidence '',
     +        f10.2,/)') confdnc
      endif

c----------------------------------------------------------------------

c     Perform clear-sky restoral tests.

      if( .not. (snow .or. ice)) then
        if(confdnc .le. 0.95) then
          call chk_land_nite(pxldat,tbadj,confdnc,qa_bits,
     *                       testbits)
        end if
      end if

      return
      end
