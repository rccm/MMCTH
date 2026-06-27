       subroutine ocean_day(pxldat,vza,snglnt,visusd,cirrus_vis,sfctmp,
     +                      refang,sh_ocean,testbits,qa_bits,nmtests,
     +                      confdnc)

       implicit none
       save 

c---------------------------------------------------------------------
C!F77 
c
c!Description:
c Performs clear-sky spectral tests for water surfaces during daylight
c conditions.
c
c Each spectral test is placed in one of five test groups. The groups
c represented in this routine are:
c
c      
c          Group 1: High thick cloud
c                   11 micron bt test 
c                   13.9 micron bt test  
c                   6.75 micron bt test 
c
c          Group 2: Low cloud - thick
c                   8-11 micron and 11-12 micron bt tests
c                   11-4 micron bt tests
c                  
c          Group 3: Thick cloud
c                   .87 micron reflectance test
c                   .87/.66 micron reflectance ratio test
c 
c          Group 4: Thin cirrus test
c                   1.38 micron reflectance test 
c
c
c A "confidence of clear sky" is computed for each spectral test.
c Confidences from single-threshold tests are calculated in subroutine
c 'conf_test'.  Those from double-threshold tests (where clear-sky
c radiance data falls in a range between two thresholds or lies on
c either side of them) are generated in 'conf_test_2val'.
c
c The minimum confidence value in each group is defined as the group
c confidence.  Final confidence is defined as the nth root of the
c product of the group confidences, where n is the number of groups.
c
c A "qa bit" (in array 'qa_bits') is set for each test performed and
c a "test bit" is set (in array 'testbits') for each clear-sky test
c passed.
c
c!Input Parameters:
c pxldat        Array containing reflectance or brightness temperatures
c               for all bands for a single pixel
c vza           Current pixel viewing angle
c snglnt        Logical variable indicating sun glint contamination
c visusd        Logical variable indicating whether vis data used or not
c cirrus_vis 	Logical variable flagging thin cirrus contaminated
c 		scenes in the visible
c sfctmp        SST from ancillary data
c refang        Reflectance angle
c sh_ocean      Logical flag indicating ocean depths < 50 m or within 5 km
c               of shoreline
c
c!Output Parameters:
c testbits      6 byte array containing bit results
c qa_bits       10 byte array contining qa bit results
c nmtests 	 Counts number of inidividual tests applied
c confdnc       product of all applied individual confidences
c
c!Revision History:
c 06/04 Collection 5  R. Frey
c Implemented new 11-12 um thin cirrus test (Key version)
c 10/04 Collection 5b R. Frey
c Added SST test.
c
c!Team-Unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c
c!Design Notes:
c    Externals:
c        Subroutines conf_test_2val,conf_test,tview,set_bit,
c                    clear_bit,check_bits,set_qa_bit
c        Functions rega,regb
c
c!END
c-------------------------------------------------------------------

c     Declarations.

      include 'global.inc'
      include 'ocean_day_thr.inc'
      include 'snglntr_thr.inc'

c     Scalar arguments. 
      real confdnc,vza,refang,sfctmp
      logical visusd,snglnt,cirrus_vis,sh_ocean
      integer nmtests

c     Array arguments. 
      real pxldat(inband)
      byte testbits(6),qa_bits(10)
        
c     Local scalars.
      real c1,c2,c3,c4,c6,c7,c8,c9,c11,cosvza,dfthrsh,diftemp,
     +     dtr,m31_22,m31_32,m29_31,m31,m32,m35,m22,m29,m26,
     +     pi,schi,vrat,m27,Rel_equality_EPS,m01,m02,
     +     cmin1,cmin2,cmin3,cmin4,tri_thres,fac,groups,
     +     pre_confdnc,locut,hicut,midpt,power,max_vza,a,corr,
     +     sfcdif,c10,sst_thrsh
      integer nptests,rtn,debug,h_output,kk

c     Local arrays.
      real hicuta(2),locuta(2),midpta(2)
      integer ngtests(5)
        
c     Parameter statements.
      parameter (Rel_equality_EPS = 0.000001)
      parameter(max_vza = 65.49)

c     External functions.
      real trispc 
      external trispc
         
c     External subroutines.
      external conf_test,tview,set_bit,clear_bit,conf_test_2val,
     +         set_qa_bit,get_sg_thresholds
c    
c     Intrinsic functions.
      intrinsic acos,cos

c     Common statement for debug purposes
      common / bug / debug, h_output

c-------------------------------------------------------------------

c     Initialize variables.

      pi = acos(-1.0)
      dtr = pi/180.0

c     'nmtests' counts the number of tests applied to this pixel.
      nmtests = 0

c     'nptests' counts the number of tests which found no evidience 
c     of cloud.
      nptests = 0

c     Place reflectance and brightness temperature values into easy-to-
c     identify variables.
      m01 = pxldat(1)
      m02 = pxldat(2)
      m26 = pxldat(26)
      m22 = pxldat(22)
      m27 = pxldat(27)
      m29 = pxldat(29)
      m31 = pxldat(31)
      m32 = pxldat(32)
      m35 = pxldat(35)

c     Initialize test group confidences.
      cmin1 = 1.0
      cmin2 = 1.0
      cmin3 = 1.0
      cmin4 = 1.0
 
c     Initialize array containing number of tests in each test group.
      do 10 kk = 1 , 5
         ngtests(kk) = 0 
  10  continue

c-------------------------------------------------------------------

c     debug statement 
      if (debug .gt. 0) then
        write(h_output,'(10x/,''Processing subroutine ocean_day '',/)')
      endif

c-------------------------------------------------------------------

c     Begin clear sky tests.

c-------------------------------------------------------------------

c     GROUP 1 TESTS

c     11 micron brightness temperature threshold test.

      if (nint(m31) .ne. nint(bad_data)) then
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,13)
        if (m31 .ge. dobt11(2)) then
          call set_bit(testbits,13)
          nptests = nptests + 1
        end if
        call conf_test(m31,dobt11(1),dobt11(3),dobt11(4),
     *                 dobt11(2),1,c1)
        cmin1 = min(cmin1,c1)
        ngtests(1) = ngtests(1) + 1
      endif

c-------------------------------------------------------------------

c     debug statement
      if (debug .gt. 0) then
        write(h_output,'(1x,''m31: '',5f10.2)') m31,dobt11(1),
     +          dobt11(2),dobt11(3),dobt11(4)
      endif

c-------------------------------------------------------------------
 
c     co2 high cloud test

      if (nint(m35) .ne. nint(bad_data)) then
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,14)
        if (m35 .gt. doco2(2)) then
          call set_bit(testbits,14)
          nptests = nptests + 1
        end if
        call conf_test(m35,doco2(1),doco2(3),doco2(4),
     *                doco2(2),1,c2)
        cmin1 = min(cmin1,c2)
        ngtests(1) = ngtests(1) + 1
      endif

c-------------------------------------------------------------------
 
c     debug statement 
      if (debug .gt. 0) then
        write(h_output,'(1x,''m35: '',5f10.2)') m35,doco2(1),
     +          doco2(2),doco2(3),doco2(4)
      endif

c-------------------------------------------------------------------

c     H20 vapor channel (6.7 micron) high cloud test 

      if (nint(m27) .ne. nint(bad_data)) then
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,15)
        if (m27 .gt. doh20(2)) then
          call set_bit(testbits,15)
          nptests = nptests + 1
        end if
        call conf_test(m27,doh20(1),doh20(3),doh20(4),
     *                doh20(2),1,c3)
        cmin1 = min(cmin1,c3)
        ngtests(1) = ngtests(1) + 1
      endif

c-------------------------------------------------------------------

c     debug statement 
      if (debug .gt. 0) then
        write(h_output,'(1x,''m27: '',5f10.2)') m27,doh20(1),
     +          doh20(2),doh20(3),doh20(4)
      endif

c-------------------------------------------------------------------

c ... SST test

      if ( (nint(m31) .ne. nint(bad_data)) .and.
     +     (nint(m32) .ne. nint(bad_data)) ) then

       if (sfctmp .gt. 0.0 .and. sfctmp .lt. 350.0) then

        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,27)

        if(sh_ocean) then
          sst_thrsh = 10.0
        else
          sst_thrsh = 6.0
        end if

        m31_32 = m31 - m32
        if(m31_32 .ge. 1.0) then
          midpt = sst_thrsh + (2.0 * nint(m31_32))
        else
          midpt = sst_thrsh
        end if

        a = vza / max_vza
        corr = (a**4) * 3.0
        midpt = midpt + corr
        locut = midpt + 1.0
        hicut = midpt - 2.0

        sfcdif = sfctmp - m31

        if( sfcdif .lt. midpt ) then
          call set_bit(testbits,27)
          nptests = nptests + 1
        end if

        call conf_test(sfcdif,locut,hicut,1.0,midpt,1,c10)
        cmin1 = min(cmin1,c10)
        ngtests(1) = ngtests(1) + 1

       endif

      endif

c-----------------------------------------------------------------

c     debug statement
      if (debug .gt. 0) then
        write(h_output,'(1x,''sfctmp: '',9f10.3)') sfctmp,m31,
     +           sfcdif,a,corr,locut,midpt,hicut,c10
      endif

c-----------------------------------------------------------------

c     GROUP 2 TESTS

c     tri-spectral tests - 8, 11 and 12 micron BTDIF's

      if (nint(m31) .ne. nint(bad_data)
     +                 .and.
     +    nint(m32) .ne. nint(bad_data)
     +                 .and.
     +    nint(m29)  .ne. nint(bad_data)) then

        m29_31 = m29 - m31
        m31_32 = m31 - m32

c       Get clear sky 8-11 micron clear sky thresholds based
c       upon 11-12 difference and compare to global regressions
c       determined from global HIRS data
        tri_thres = trispc(m31_32) 
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,18) 
        if (m29_31 .lt. tri_thres) then
          nptests = nptests + 1
          call set_bit(testbits,18)
        end if
        locut = tri_thres + .5
        hicut = tri_thres - .5
        call conf_test(m29_31,locut,hicut,1.0,tri_thres,1,c4)
        cmin2 = min(cmin2,c4)
        ngtests(2) = ngtests(2) + 1
      endif
 
c-------------------------------------------------------------------

c     debug statement 
      if (debug .gt. 0) then
        write(h_output,'(1x,''trispec: '',5f10.2)') m31_32,m29_31,
     +          locut,tri_thres,hicut
      endif

c-------------------------------------------------------------------

c       11-12um brightness temperature difference test
c       for thin cirrus.

      if (nint(m31) .ne. nint(bad_data)
     +                 .and.
     +    nint(m32) .ne. nint(bad_data)
     +                  .and.
     +              vza .gt. 0.0) then

        m31_32 = m31 - m32

c       calculate secant of viewing zenith angle.
        cosvza = cos(vza*dtr)
        if (abs(cosvza).gt.Rel_equality_EPS) then
          schi = 1.0/cosvza
        else
          schi = 99.0
        end if
 
c       Interpolate look-up table values of 11 - 12 micron bt
c       difference thresholds (function of viewing zenith
c       and 11 micron brightness temperature).
        call tview(1,schi,m31,diftemp)
 
        if (diftemp.lt.0.1 .or. abs(schi-99.0).lt.0.0001) then
          dfthrsh = do11_12hi(1)
        else
          dfthrsh = diftemp
        end if
 
c       Since the IR BTDIF bit has possibly been set already,
c       change the bit only if the current test failed.
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,18) 
        if (m31_32 .le. dfthrsh) then
          nptests = nptests + 1
        else
          call check_bits(testbits,18,rtn)
          if (rtn .eq. 1) then
            call clear_bit(testbits,18)
          end if
        endif
        locut = dfthrsh + (0.3 * dfthrsh)
        hicut = dfthrsh - 1.25
        call conf_test(m31_32,locut,hicut,1.0,dfthrsh,1,c6)
        cmin2 = min(cmin2,c6)
        ngtests(2) = ngtests(2) + 1
      endif

c-------------------------------------------------------------------

c     debug statement 
      if (debug .gt. 0) then
        write(h_output,'(1x,''APOLLO m31_32: '',4f10.2)') m31_32,
     +          dfthrsh,locut,hicut
      endif

c-------------------------------------------------------------------

c     11 minus 4 micron BTDIF fog and low cloud test.

      if (visusd .and. .not. snglnt) then
        if (nint(m31) .ne. nint(bad_data)
     +                   .and.
     +      nint(m22) .ne. nint(bad_data)) then

          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,19) 
          m31_22 = m31 - m22

          if (m31_22 .ge. do11_4lo(2)) then
            call set_bit(testbits,19)
            nptests = nptests + 1
          end if

          call conf_test(m31_22,do11_4lo(1),do11_4lo(3),do11_4lo(4),
     +                  do11_4lo(2),1,c7)

          cmin2 = min(cmin2,c7)
          ngtests(2) = ngtests(2) + 1

        endif

c-------------------------------------------------------------------

c       debug statement 
        if (debug .gt. 0) then
         write(h_output,'(1x,''m31_22: '',5f10.2)')m31_22,do11_4lo(1),
     +            do11_4lo(2),do11_4lo(3),do11_4lo(4)
        endif

c-------------------------------------------------------------------

      end if

c     GROUP 3 TESTS 

c     NIR reflectance threshold test.

      if (visusd) then
        if (nint(m02) .ne. nint(bad_data)) then

c         Take into account sunglint problems
          if(snglnt) then
            call get_sg_thresholds(refang,locut,hicut,midpt,power)
          else
            locut = doref2(1)
            hicut = doref2(3)
            midpt = doref2(2)
            power = doref2(4)
          end if

          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,20) 
          if (m02.le.midpt) then
            call set_bit(testbits,20)
            nptests = nptests + 1
          end if
          call conf_test(m02,locut,hicut,power,
     +                   midpt,1,c8)
          cmin3 = min(cmin3,c8)
          ngtests(3) = ngtests(3) + 1
        endif

c-------------------------------------------------------------------

c       debug statement 
        if (debug .gt. 0) then
          write(h_output,'(1x,''m02: '',6f10.4)') m02,locut,
     +            hicut,midpt,power,refang
        endif

c-------------------------------------------------------------------

      end if

c     Visible channel ratio test 

      if (visusd) then
        if (nint(m01) .ne. nint(bad_data)
     +                   .and.
     +      nint(m02) .ne. nint(bad_data)) then

c         Account for sun glint contamination
          if (snglnt) then
            locuta(1) = snglntvcl(1)
            locuta(2) = snglntvcl(2)
            hicuta(1) = snglntvch(1)
            hicuta(2) = snglntvch(2)
            midpta(1) = snglntv(1)
            midpta(2) = snglntv(2)
          else
            locuta(1) = dovratlo(1)
            locuta(2) = dovrathi(1)
            hicuta(1) = dovratlo(3)
            hicuta(2) = dovrathi(3)
            midpta(1) = dovratlo(2)
            midpta(2) = dovrathi(2)
          end if

          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,21) 
          vrat = m02 / m01
          if (vrat .lt. midpta(1) .or. vrat .gt. midpta(2)) then
            call set_bit(testbits,21)
            nptests = nptests + 1
          end if
          call conf_test_2val(vrat,locuta,hicuta,1.0,midpta,2,c9)
          cmin3 = min(cmin3,c9)
          ngtests(3) = ngtests(3) + 1
        endif

c-------------------------------------------------------------------

c       debug statement 
        if (debug .gt. 0) then
          write(h_output,'(1x,''vrat: '',7f10.2)') vrat,locuta(1),
     +            locuta(2),hicuta(1),hicuta(2),midpta(1),midpta(2)
        endif

c-------------------------------------------------------------------

      end if

c     GROUP 4 TESTS 

c ... Near-infrared high cloud test.

      if (visusd) then
        if (nint(m26) .ne. nint(bad_data)) then
          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,16) 
          if (m26 .le. doref3(2)) then
            call set_bit(testbits,16)
            nptests = nptests + 1
          end if
          call conf_test(m26,doref3(1),doref3(3),doref3(4),
     +                  doref3(2),1,c11)
          cmin4 = min(cmin4,c11)
          ngtests(4) = ngtests(4) + 1
        endif

c-------------------------------------------------------------------

c       debug statement 
        if (debug .gt. 0) then
           write(h_output,'(1x,''m26: '',6f10.4)')m26,doref3(1),
     +                doref3(2),doref3(3),doref3(4)
        endif

c-------------------------------------------------------------------

      end if

c     Thin cirrus test.

      if (visusd) then
        if (nint(m26) .ne. nint(bad_data)) then
          call set_qa_bit(qa_bits,9)
          if (m26 .lt. dotci(1) .and. m26 .ge. dotci(2)) then
            call clear_bit(testbits,9)
            cirrus_vis = .true.
          endif

c-------------------------------------------------------------------

c         debug statement 
          if (debug .gt. 0) then
             write(h_output,'(1x,''NIR Thin cirrus: '',3f10.4)') m26,
     +             dotci(1),dotci(2)
          endif

c-------------------------------------------------------------------

        endif
      endif
 
c-------------------------------------------------------------------

c     Determine initial confidence based on group values.
      pre_confdnc = cmin1 * cmin2 * cmin3 * cmin4
      
c     Find the number of test groups used for current pixel.
      groups = 0
      do kk = 1,5
        if(ngtests(kk) .gt. 0) then
          groups = groups + 1.0
        end if
      enddo
      fac = 1.0
      if (groups .gt. 0) fac = 1.0 / groups

c     Final pixel confidence is nth root of group confidence product.
      confdnc = pre_confdnc**fac

c-------------------------------------------------------------------

c     debug statement 
      if (debug .gt. 0) then
        write(h_output,'(1x,''tests '',6i10)') nmtests,nptests,ngtests
        write(h_output,'(1x,''confdnc '',7f8.5/,8f8.5)') c1,c2,c3,c4,
     +         c6,c7,c8,c9,c11,cmin1,cmin2,cmin3,cmin4,fac,confdnc
      endif

c-------------------------------------------------------------------

      return
      end
