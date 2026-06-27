      subroutine mod06uw_stats_out(pcnnl, pcnnw, pcnni, ng_co2, ng_irp,
     +                     nmpix, nl, nm, nh, nn, nk,nco2,nbias, np, nc,
     +                     ni, nw, nx, nu, success_co2, success_irp,
     +                     plow, pmid, phigh, pthin, pthick, popq, 
     +                     pcirrus, pice, pwater, pmixed, punc, pocean,
     +                     pland, psnow )


      implicit none
      save

c----------------------------------------------------------------------
C!F77 
c
c!Description:
c     Routine for calculating granule level output statistics
c
c!Input parameters:
c pcnnl         Percentage of land processed pixels per granule
c pcnnw         Percentage of water processed pixels per granule
c pcnni         Percentage of snow/ice processed pixels per granule
c ng_co2        Number of good co2 retrievals
c ng_irp        Number of good IRPHASE retrievals
c nmpix         Total number of pixels used to produce stats
c nl            Number of low cloud pixels found
c nm            Number of mid cloud pixels found
c nh            Number of high cloud pixels found
c nn            Number of thin cloud pixels found
c nk            Number of thick cloud pixels found
c np            Number of opaque cloud pixels found
c nc            Number of cirrus cloud pixels found
c ni            Number of pixels with ice cloud found
c nw            Number of pixels with water cloud found
c nx            Number of pixels with mixed cloud found
c nu            Number of pixels with uncertain clouds found
c nco2          Number of co2 retrival 
c nbias         Number of bias applied
c
c!Output Parameters: None
c
c success_co2   Percentage of pixels with successful CO2 retrievals
c success_irp   Percentage of pixels with successful irp retrievals
c plow          Percentage of pixels found to have low cloud
c pmid          Percentage of pixels found to have mid cloud
c phigh         Percentage of pixels found to have high cloud
c pthin         Percentage of pixels found to have thin cloud
c pthick        Percentage of pixels found to have thick cloud
c popq          Percentage of pixels found to have opaque cloud
c pcirrus       Percentage of pixels found to have cirrus cloud
c pice          Percentage of pixels found to have ice cloud
c pwater        Percentage of pixels found to have water cloud
c pmixed        Percentage of pixels found to have mixed phase cloud
c punc          Percentage of pixels found to have uncertain phase cloud
c pocean        Percentage of pixels with ocean background
c pland         Percentage of pixels with land background
c psnow         Percentage of pixels with snow background
c pco2          Percentage of co2 retrival 
c pbias         Percentage of bias applied
c
c!Revision History:
c $Id: mod06uw_stats_out.f,v 1.1.1.1 2005/02/22 17:15:54 gumley Exp $
c
c!Team-Unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c
c!END
c----------------------------------------------------------------------
      include 'mod06uw_debug.inc'

c     scalar arguments
      integer ng_co2, ng_irp, nmpix, nl, nm, nh, nn, nk, np, nc,
     +        ni, nw, nx, nu,nco2,nbias
      real pcnnl, pcnnw, pcnni, success_co2, success_irp,
     +     plow, pmid, phigh, pthin, pthick, popq,
     +     pcirrus, pice, pwater, pmixed, punc, pocean,
     +     pland, psnow, pco2, pbias 

c     parameters - fill value for metadata
      real Meta_misg
      parameter (Meta_misg = -99999.0)


c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(10x/,''Within stats_out routine '',/)')
      endif
c ...............................................................

      if(nmpix .gt. 0) then
         success_co2 = (float(ng_co2) / nmpix) * 100.0
         success_irp = (float(ng_irp) / nmpix) * 100.0
         plow   = (float(nl) / nmpix) * 100.0
         pmid   = (float(nm) / nmpix) * 100.0
         phigh  = (float(nh) / nmpix) * 100.0
         pthin  = (float(nn) / nmpix) * 100.0
         pthick = (float(nk) / nmpix) * 100.0
         popq   = (float(np) / nmpix) * 100.0
         pcirrus= (float(nc) / nmpix) * 100.0
         pice   = (float(ni) / nmpix) * 100.0
         pwater = (float(nw) / nmpix) * 100.0
         pmixed = (float(nx) / nmpix) * 100.0
         punc   = (float(nu) / nmpix) * 100.0
         pco2   = (float(nco2) / nmpix) * 100.0
         pbias   = (float(nbias) / nmpix) * 100.0
         pocean =  pcnnw
         pland  =  pcnnl
         psnow  =  pcnni
      else
c        write error message to log file
         call message( 'stats_out',
     &   'Warning - No good pixels were found to process ',
     &   0, 3 )

         success_co2 = 0.0
         success_irp = 0.0
         plow   =  Meta_misg
         pmid   =  Meta_misg
         phigh  =  Meta_misg
         pthin  =  Meta_misg
         pthick =  Meta_misg
         popq   =  Meta_misg
         pcirrus=  Meta_misg
         pice   =  Meta_misg
         pwater =  Meta_misg
         pmixed =  Meta_misg
         punc   =  Meta_misg
         pocean =  Meta_misg
         pland  =  Meta_misg
         psnow  =  Meta_misg
      end if

      write (h_output,fmt='(/,30x,''Granule Statstics: '',/)')
      write (h_output,fmt='(5x,''Parameter         Total Number '',
     +       ''    Percent '',/)')
      write (h_output,fmt='(5x,''Good SS CTP rets. '',I10,f15.2)')
     +       ng_co2,success_co2
      write (h_output,fmt='(5x,''Good CO2 CTP rets. '',I10,f15.2)')
     +       nco2,pco2
      write (h_output,fmt='(5x,''Good bias applied. '',I10,f15.2)')
     +       nbias,pbias
      write (h_output,fmt='(5x,''Good IRP rets. '',I10,f15.2)')
     +       ng_irp,success_irp
      write (h_output,fmt='(5x,''Low Cloud      '',I10,f15.2)')
     +       nl,plow
      write (h_output,fmt='(5x,''Mid Cloud      '',I10,f15.2)')
     +       nm,pmid
      write (h_output,fmt='(5x,''High Cloud     '',I10,f15.2)')
     +       nh,phigh
      write (h_output,fmt='(5x,''Thin Cloud     '',I10,f15.2)')
     +       nn,pthin
      write (h_output,fmt='(5x,''Thick Cloud    '',I10,f15.2)')
     +       nk,pthick
      write (h_output,fmt='(5x,''Cirrus Cloud   '',I10,f15.2)')
     +       nc,pcirrus
      write (h_output,fmt='(5x,''Opaque Cloud   '',I10,f15.2)')
     +       np,popq
      write (h_output,fmt='(5x,''Ice Cloud      '',I10,f15.2)')
     +       ni,pice
      write (h_output,fmt='(5x,''Water Cloud    '',I10,f15.2)')
     +       nw,pwater
      write (h_output,fmt='(5x,''Mixed Phase    '',I10,f15.2)')
     +       nx,pmixed
      write (h_output,fmt='(5x,''Uncertain Phase'',I10,f15.2,/)')
     +       nu,punc
      write (h_output,fmt='(1x,''pct ocean background '',f9.2)')pocean
      write (h_output,fmt='(1x,''pct land background  '',f9.2)')pland
      write (h_output,fmt='(1x,''pct snow background  '',f9.2)')psnow
      write (h_output,fmt='(4x)')
      write (h_output,fmt='(1x,''Total pixels: '',i10)')nmpix

      return
      end
