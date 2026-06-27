      subroutine store_cm_scan(nc,confdnc,day,visusd,snglnt,snow,ice,coast,
     *                         desert,land,process,scan_confdnc,
     *                         scan_visusd,scan_snglnt,scan_snow,scan_ice,
     *                         scan_coast,scan_desert,scan_land,
     *                         scan_process,scan_day)

      implicit none
      save

c-----------------------------------------------------------------------
c!F77
c
c!DESCRIPTION:  Routine to save 1-km values for each pixel on current scan.
c
c!Input Parameters:
c nc            Current pixel number
c confdnc       Clear-sky confidence of current pixel
c day           Logical variable flagging day scenes
c visusd        Logical variable indicating vis. data used
c snglnt        Logical variable flagging sunglint regions
c snow          Logical variable flagging snow processing 
c ice           Logical variable flagging ice processing 
c coast         Logical variable flagging coast processing 
c desert        Logical variable flagging desert processing 
c land          Logical variable flagging land processing 
c process       Logical variable indicating current pixel processed
c
c!Output Parameters:
c scan_confdnc  Clear-sky confidence array for pixels on current scan
c scan_visusd   Indicates vis. data used for pixels on current scan
c scan_snglnt   Indicates sunglint for pixels on current scan
c scan_snow     Indicates snow processing for pixels on current scan
c scan_ice      Indicates ice processing for pixels on current scan
c scan_coast    Indicates coast processing for pixels on current scan
c scan_desert   Indicates desert processing for pixels on current scan
c scan_land     Indicates land processing for pixels on current scan
c scan_process  Indicates pixel processed for pixels on current scan
c scan_day      Indicates day processing for pixels on current scan
c
c!Revision History:
c
c!Team-Unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-35.
c
c!END
c---------------------------------------------------------------------

      INCLUDE 'global.inc'

c     Array arguments.
      real confdnc,scan_confdnc(npixel)
      logical scan_day(npixel),scan_visusd(npixel),scan_snglnt(npixel),
     *        scan_snow(npixel),scan_ice(npixel),scan_coast(npixel),
     *        scan_desert(npixel),scan_land(npixel),
     *        scan_process(npixel)

c     Array scalars.
      integer nc,debug,h_output
      logical day,visusd,snglnt,snow,ice,coast,desert,land,process

c ... Common statement for debug purposes
      common / bug / debug, h_output

c---------------------------------------------------------------------

c     Fill arrays.
 
      scan_confdnc(nc) = confdnc
      scan_visusd(nc) = visusd
      scan_snglnt(nc) = snglnt
      scan_snow(nc) = snow
      scan_ice(nc) = ice
      scan_coast(nc) = coast
      scan_desert(nc) = desert
      scan_land(nc) = land
      scan_process(nc) = process
      scan_day(nc) = day

c---------------------------------------------------------------------

      return
      end
