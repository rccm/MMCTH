      subroutine stats_out(nmpix,num_bad,n1s,n2s,n3s,n4s,no,ns,
     +                     nd,nt,ng,ni,nl,nw,nr,nv,n250,ncr,
     +                     ncd,nu,nn,nz,na,ne,pcbad,pcn1s,pcn2s,
     +                     pcn3s,pcn4s,pcnnc,pcnns,pcnnd,pcnnt,pcnng,
     +                     pcnni,pcnnl,pcnnw,pcnnr,pcnnv,pcnncr,
     +                     pcnncd,pcnu,pcnn,pcnz,pcna,pcne)

      implicit none
      save

c ... Common statement for debug purposes
      common / bug / debug, h_output


c----------------------------------------------------------------------
C!F77 
c
c!Description:
c     Routine for calculating granule level output statistics
c
c!Input parameters:
c npix          Total number of pixels used to produce stats
c num_bad       Sum of bad pixels
c n1s           Confidence category 1 sum
c n2s           Confidence category 2 sum
c n3s           Confidence category 3 sum
c n4s           Confidence category 4 sum
c no            Sum of non-cloud obstruction pixels
c ns            Sum of shadow pixels found
c nd            Sum of day processed pixels
c nt            Sum of night processed pixels
c ng            Sum of sunglint found pixles
c ni            Sum of snow/ice processed pixels
c nl            Sum of land processed pixels
c nw            Sum of water processed pixels
c nr            Sum of thin cirrus (ir) found pixels
c nv            Sum of thin cirrus (vis) found pixels
c n250          Counter for number of pixels included in stats
c               (includes all possible pixels)
c ncr           Number of 250m pixels found to be clear
c               (out of all possible pixels)
c ncd           Number of 250m pixels found to be cloudy
c               (out of all possible pixels)
c nu            Sum of bad lat pixels
c nn            Sum of bad lon pixels
c nz            Sum of bad szen pixels
c na            Sum of bad vzen pixels
c ne            Sum of bad rel pixels
c
c!Output Parameters: None
c
c pcbad         Percentage of bad pixels
c pcn1s         Percentage of pixels in category 1 per granule
c pcn2s         Percentage of pixels in category 2 per granule
c pcn3s         Percentage of pixels in category 3 per granule
c pcn4s         Percentage of pixels in category 4 per granule
c pcnnc         Percentage of non-cloud obstruction pixels per granule
c pcnns         Percentage of shadow pixels found per granule
c pcnnd         Percentage of day processed pixels per granule
c pcnnt         Percentage of night processed pixels per granule
c pcnng         Percentage of sunglint found pixles per granule
c pcnni         Percentage of snow/ice processed pixels per granule
c pcnnl         Percentage of land processed pixels per granule
c pcnnw         Percentage of water processed pixels per granule
c pcnnr         Percentage of thin cirrus (ir) found pixels per granule
c pcnnv         Percentage of thin cirrus (vis) found pixels per granule
c pcnncr        Number of 250m pixels found to be clear per granule
c               (out of all possible pixels)
c pcnncd        Number of 250m pixels found to be cloudy per granule
c               (out of all possible pixels)
c pcnu          Percentage of bad lat pixels in granule
c pcnn          Percentage of bad lon pixels in granule
c pcnz          Percentage of bad sza pixels in granule
c pcna          Percentage of bad vza pixels in granule
c pcne          Percentage of bad rel angle  pixels in granule
c
c!Revision History:
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
      integer nmpix,n1s,n2s,n3s,n4s,num_bad,no,ns,nd,nt,ng,ni,nl,
     +        nw,nr,nv,n250,ncr,ncd,nu,nn,nz,na,ne
      real pcn1s,pcn2s,pcn3s,pcn4s,pcbad,pcnnc,pcnns,pcnnd,pcnnt,
     +     pcnng,pcnni,pcnnl,pcnnw,pcnnr,pcnnv,pcnncr,pcnncd,
     +     pcnu,pcnn,pcnz,pcna,pcne

c     local arguments
      integer h_output,debug

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(10x/,''Within stats_out routine '',/)')
      endif
c ...............................................................

      if(nmpix .gt. 0) then
         pcn1s = (float(n1s) / nmpix) * 100.0
         pcn2s = (float(n2s) / nmpix) * 100.0
         pcn3s = (float(n3s) / nmpix) * 100.0
         pcn4s = (float(n4s) / nmpix) * 100.0
         pcbad = (float(num_bad) / nmpix) * 100.0
         pcnnc = (float(no) / nmpix) * 100.0
         pcnns = (float(ns) / nmpix) * 100.0
         pcnnd = (float(nd) / nmpix) * 100.0
         pcnnt = (float(nt) / nmpix) * 100.0
         pcnng = (float(ng) / nmpix) * 100.0
         pcnni = (float(ni) / nmpix) * 100.0
         pcnnl = (float(nl) / nmpix) * 100.0
         pcnnw = (float(nw) / nmpix) * 100.0
         pcnnr = (float(nr) / nmpix) * 100.0
         pcnnv = (float(nv) / nmpix) * 100.0
         pcnncr = (float(ncr) / n250) * 100.0
         pcnncd = (float(ncd) / n250) * 100.0
         pcnu = (float(nu) / nmpix) * 100.0
         pcnn = (float(nn) / nmpix) * 100.0
         pcnz = (float(nz) / nmpix) * 100.0
         pcna = (float(na) / nmpix) * 100.0
         pcne = (float(ne) / nmpix) * 100.0
      else
c        write error message to log file
         call message( 'stats_out',
     &   'Warning - No good pixels were found to process' //
     &   ' [OPERATOR ACTION:  Notify SDST]',
     &   0, 1 )

         pcn1s = Meta_misg
         pcn2s = Meta_misg
         pcn3s = Meta_misg
         pcn4s = Meta_misg
         pcbad = 100.0
         pcnnc = Meta_misg
         pcnns = Meta_misg
         pcnnd = Meta_misg
         pcnnt = Meta_misg
         pcnng = Meta_misg
         pcnni = Meta_misg
         pcnnl = Meta_misg
         pcnnw = Meta_misg
         pcnnr = Meta_misg
         pcnnv = Meta_misg
         pcnncr = Meta_misg
         pcnncd = Meta_misg
         pcnu = Meta_misg
         pcnn = Meta_misg
         pcnz = Meta_misg
         pcna = Meta_misg
         pcne = Meta_misg
      end if

      write (h_output,fmt='(/,30x,''Granule Statstics: '',/)')
      write (h_output,fmt='(5x,''Parameter         Total Number '', 
     +       ''    Percent '',/)')
      write (h_output,fmt='(5x,''Bad data       '',I10,f15.2)')
     +       num_bad,pcbad
      write (h_output,fmt='(5x,''NCO            '',I10,f15.2)')
     +       no,pcnnc
      write (h_output,fmt='(5x,''Shadow         '',I10,f15.2)')
     +       ns,pcnns
      write (h_output,fmt='(5x,''Day            '',I10,f15.2)')
     +       nd,pcnnd
      write (h_output,fmt='(5x,''Night          '',I10,f15.2)')
     +       nt,pcnnt
      write (h_output,fmt='(5x,''Sunglint       '',I10,f15.2)')
     +       ng,pcnng
      write (h_output,fmt='(5x,''Snow/ice       '',I10,f15.2)')
     +                           ni,pcnni
      write (h_output,fmt='(5x,''Land           '',I10,f15.2)')
     +       nl,pcnnl
      write (h_output,fmt='(5x,''Water          '',I10,f15.2)')
     +       nw,pcnnw
      write (h_output,fmt='(5x,''Thin cirrus sol'',I10,f15.2)')
     +       nv,pcnnv
      write (h_output,fmt='(5x,''Thin cirrus ir '',I10,f15.2)')
     +                           nr,pcnnr
      write (h_output,fmt='(5x,''250m clear     '',I10,f15.2)')
     +                           ncr,pcnncr
      write (h_output,fmt='(5x,''250m cloudy    '',I10,f15.2)')
     +                           ncd,pcnncd
      write (h_output,fmt='(5x,''Bad Latitude   '',I10,f15.2)')
     +       nu,pcnu
      write (h_output,fmt='(5x,''Bad Longitude  '',I10,f15.2)')
     +       nn,pcnn
      write (h_output,fmt='(5x,''Bad Solar Zen. '',I10,f15.2)')
     +       nz,pcnz
      write (h_output,fmt='(5x,''Bad View. Zen. '',I10,f15.2)')
     +       na,pcna
      write (h_output,fmt='(5x,''Bad Rel. Angle '',I10,f15.2)')
     +       ne,pcne
      write (h_output,fmt='(4x)')
      write (h_output,fmt='(1x,''# with confidence > 99% '',I9)')n1s
      write (h_output,fmt='(1x,''# with confidence > 95% '',I9)')n2s
      write (h_output,fmt='(1x,''# with confidence > 66% '',I9)')n3s
      write (h_output,fmt='(1x,''# with confidence <  1% '',I9,/)')n4s
      write (h_output,fmt='(1x,''% with confidence > 99% '',f9.2)')pcn1s
      write (h_output,fmt='(1x,''% with confidence > 95% '',f9.2)')pcn2s
      write (h_output,fmt='(1x,''% with confidence > 66% '',f9.2)')pcn3s
      write (h_output,fmt='(1x,''% with confidence <  1% '',f9.2)')pcn4s
      write (h_output,fmt='(4x)')
      write (h_output,fmt='(1x,''% of bad pixels '',f9.2)')pcbad
      write (h_output,fmt='(/,1x,''Total number of pixels '',I10)')nmpix

      return
      end
