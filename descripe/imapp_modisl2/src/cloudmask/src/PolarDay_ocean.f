       subroutine PolarDay_ocean(pxldat,vza,snglnt,visusd,cirrus_vis,
     +                           refang,sfctmp,sh_ocean,testbits,
     +                           qa_bits,nmtests,confdnc)
       implicit none
       save 

c---------------------------------------------------------------------
C!F77 
c
c!Description:
c      Routine for performing clear sky tests over polar water
c      surfaces during daylight hours.
c
c      For daytime polar ocean the groups are:
c          Group 1: High thick cloud
c                   11 micron bt test (masir11)
c                   6.75 micron bt test (not is use with mas)
c
c          Group 2: Low cloud - thick
c                   8-11 micron and 11-12 micron bt tests
c                   11-4 micron bt tests
c                  
c          Group 3: Thick cloud
c                   .87 micron reflectance test (masv88)
c                   .87/.66 micron reflectance ratio test
c 
c          Group 4: Thin cirrus test
c                   1.38 micron reflectance test 
c
c
c!Input Parameters:
c pxldat        Array containing reflectance or brightness temperatures
c               for all bands for a single pixel
c vza           Current pixel viewing angle
c snglnt        Logical variable indicating sun glint contamination
c visusd        Logical variable indicating whether vis data used or not
c cirrus_vis 	Logical variable flagging thin cirrus contaminated
c 		scenes in the visible
c refang        reflectance angle
c sfctmp        SST for current pixel
c sh_ocean      Logical flag indicating ocean depths < 50 m or within 5 km
c
c!Output Parameters:
c testbits      6 byte array containing bit results
c qa_bits       10 byte array contining qa bit results
c nmtests 	 Counts number of inidividual tests applied
c confdnc       product of all applied individual confidences
c
c!Revision History:
c 06/04 Collection 5  R. Frey:
c Implemented new version of 11-12 um thin cirrus test (J. Key version)
c 10/04 Collection 5b  R. Frey:
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

      include 'global.inc'
      include 'PolarDay_ocean_thr.inc'
      include 'snglntr_thr.inc'

c ...
c ... scalar arguments ..
      real confdnc,vza,refang,sfctmp
      logical visusd,snglnt,cirrus_vis,sh_ocean
      integer nmtests
c ...
c ... array arguments ..
      real pxldat(inband)
      byte testbits(6),qa_bits(10)
c ...
c ... local scalars ..
      real c1,c3,c4,c6,c7,cosvza,dfthrsh,diftemp,
     +     dtr,mas11_4,masdf1,masdf2,masir11,masir12,
     +     masir4,masir8,masv188,masv66,masv88,
     +     pi,schi,vrat,masir65,c8,c9,Rel_equality_EPS,
     +     cmin1,cmin2,cmin3,cmin4,tri_thres,fac,groups,
     +     pre_confdnc,c11,locut,hicut,midpt,power,c10,sst_thrsh,
     +     a,max_vza,corr,sfcdif
      integer nptests,rtn,debug,h_output,kk
      parameter (Rel_equality_EPS = 0.000001)
      parameter(max_vza = 65.49)
         
c ... local arrays ..
      real hicuta(2),locuta(2),midpta(2)
      integer ngtests(4)
c ...
c ... external functions ..
      real trispc 
      external trispc
c ...
c ... external subroutines ..
      external conf_test,tview,set_bit,clear_bit,conf_test_2val,
     +         set_qa_bit
c ...
c ... intrinsic functions ..
      intrinsic acos,cos

c ... Common statement for debug purposes
      common / bug / debug, h_output

c ...
c ... initialize variables
      pi = acos(-1.0)
      dtr = pi/180.0

c ... nmtests counts the number of tests applied to this pixel
      nmtests = 0
c ... nptests counts the number of tests passed
      nptests = 0
c ... set confidence to 1.0 to begin with
      confdnc = 1.0
c ... place band values into individual variables for easy
c ... identification 
      masv66 = pxldat(1)
      masv88 = pxldat(2)
      masv188 = pxldat(26)
      masir4 = pxldat(22)
      masir65 = pxldat(27)
      masir8 = pxldat(29)
      masir11 = pxldat(31)
      masir12 = pxldat(32)

c ...
      masdf2 = 0.0
      masdf1 = 0.0 
      cosvza = 0.0
      schi = 0.0
      diftemp = 0.0
      dfthrsh = 0.0
      vrat = 0.0
      mas11_4 = 0.0
      rtn = 0

c ... the c suffix variables represent individual test confidences
      c1 = 0.0
      c3 = 0.0
      c4 = 0.0
      c6 = 0.0
      c7 = 0.0
      c8 = 0.0
      c9 = 0.0
      c10 = 0.0
      c11 = 0.0
 
c     cmin variables represent group confidences
      cmin1 = 1.0
      cmin2 = 1.0
      cmin3 = 1.0
      cmin4 = 1.0
 
c ... initialize group number holder
      do 10 kk = 1 , 4
         ngtests(kk) = 0 
  10  continue

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(10x/,''Subroutine PolarDay_ocean '',/)')
      endif
c ................................................................


C     **** GROUP 1 TESTS *************************************
c     11 micron brightness temperature threshold test
      if (nint(masir11) .ne. nint(bad_data)) then
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,13)
c ...   compare to daytime ocean threshold, set clear bit if passed
        if (masir11 .ge. pdobt11(2)) then
          call set_bit(testbits,13)
          nptests = nptests + 1
        end if
c ...   calculate confidence compared to low and high confidence cutoffs
        call conf_test(masir11,pdobt11(1),pdobt11(3),pdobt11(4),
     *                 pdobt11(2),1,c1)
        cmin1 = min(cmin1,c1)
        ngtests(1) = ngtests(1) + 1
      endif

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(1x,''masir11: '',5f10.2)') masir11,pdobt11(1),
     +          pdobt11(2),pdobt11(3),pdobt11(4)
      endif
c ................................................................

 
c ... H20 vapor channel (6.7 micron) high cloud test 
      if (nint(masir65) .ne. nint(bad_data)) then
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,15)
        if (masir65 .gt. pdoh20(2)) then
          call set_bit(testbits,15)
          nptests = nptests + 1
        end if
        call conf_test(masir65,pdoh20(1),pdoh20(3),pdoh20(4),
     *                pdoh20(2),1,c3)
        cmin1 = min(cmin1,c3)
        ngtests(1) = ngtests(1) + 1
      endif

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(1x,''masir65: '',5f10.2)') masir65,pdoh20(1),
     +          pdoh20(2),pdoh20(3),pdoh20(4)
      endif
c ................................................................

c ... SST test

      if ( (nint(masir11) .ne. nint(bad_data)) .and.
     +     (nint(masir12) .ne. nint(bad_data)) ) then

       if (sfctmp .gt. 0.0 .and. sfctmp .lt. 350.0) then

        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,27)

        if(sh_ocean) then
          sst_thrsh = 10.0
        else
          sst_thrsh = 6.0
        end if

        masdf1 = masir11 - masir12
        if(masdf1 .ge. 1.0) then
          midpt = sst_thrsh + (2.0 * nint(masdf1))
        else
          midpt = sst_thrsh
        end if

        a = vza / max_vza
        corr = (a**4) * 3.0
        midpt = midpt + corr
        locut = midpt + 1.0
        hicut = midpt - 2.0

        sfcdif = sfctmp - masir11

        if( sfcdif .lt. midpt ) then
          call set_bit(testbits,27)
          nptests = nptests + 1
        end if

        call conf_test(sfcdif,locut,hicut,1.0,midpt,1,c10)
        cmin1 = min(cmin1,c10)
        ngtests(1) = ngtests(1) + 1

       endif

      endif

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(1x,''sfctmp: '',9f10.3)') sfctmp,masir11,
     +           sfcdif,a,corr,locut,midpt,hicut,c10
      endif
c ................................................................

c     *****  END OF GROUP 1 TESTS  *************************** 
 
 
 
c     ****  GROUP 2 TESTS  ***********************************
c ... tri-spectral tests - 8, 11 and 12 micron BTDIF's
c ... calculate 8 minus 11 and 11 minus 12 micron BTDIFs
      if (nint(masir11) .ne. nint(bad_data)
     +                 .and.
     +    nint(masir12) .ne. nint(bad_data)
     +                 .and.
     +    nint(masir8)  .ne. nint(bad_data)) then

        masdf2 = masir8 - masir11
        masdf1 = masir11 - masir12

c       Get clear sky 8-11 micron clear sky thresholds based
c       upon 11-12 difference and compare to global regressions
c       determined from global HIRS data
        tri_thres = trispc(masdf1) 
c       tri_thres = 100.0
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,18) 
        if (masdf2.lt.tri_thres) then
          nptests = nptests + 1
          call set_bit(testbits,18)
        end if
        locut = tri_thres + .5
        hicut = tri_thres - .5
        call conf_test(masdf2,locut,hicut,1.0,tri_thres,1,c4)
        cmin2 = min(cmin2,c4)
        ngtests(2) = ngtests(2) + 1
      endif
 
c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(1x,''trispec: '',5f10.2)') masdf1,masdf2,
     +          locut,tri_thres,hicut
      endif
c ................................................................

      if (nint(masir11) .ne. nint(bad_data)
     +                 .and.
     +    nint(masir12) .ne. nint(bad_data)
     +                  .and.
     +              vza .gt. 0.0) then

        masdf1 = masir11 - masir12

c ...   11-12um brightness temperature difference test
c ...   for thin cirrus).
c ...   added apollo viewing angle/av4t regressed threshold.
c ...   calculate secant of viewing zenith angle.
        cosvza = cos(vza*dtr)
        if (abs(cosvza).gt.Rel_equality_EPS) then
          schi = 1.0/cosvza
        else
          schi = 99.0
        end if
 
c ...   interpolate look-up table values of 11 - 12 micron bt
c ...   difference thresholds (function of viewing zenith
c ...   and 11 micron brightness temperature).
        call tview(1,schi,masir11,diftemp)
 
c ...   if a threshold was determined by apollo, then use this
c ...   as the thin cirrus test, otherwise use a standard threshold
        if (diftemp.lt.0.1 .or. abs(schi-99.0).lt.0.0001) then
          dfthrsh = pdo11_12hi(1)
        else
          dfthrsh = diftemp
        end if
 
c ...   Since the IR BTDIF testbit was already potentially set,
c ...   then only change the bit if the current test failed
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,18) 
        if (masdf1.le.dfthrsh) then
          nptests = nptests + 1
        else
          call check_bits(testbits,18,rtn)
          if (rtn .eq. 1) then
            call clear_bit(testbits,18)
          end if
        endif
        locut = dfthrsh + (0.3 * dfthrsh)
        hicut = dfthrsh - 1.25
        call conf_test(masdf1,locut,hicut,1.0,dfthrsh,1,c6)
        cmin2 = min(cmin2,c6)
        ngtests(2) = ngtests(2) + 1
      endif

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(1x,''APOLLO masdf1: '',4f10.2)') masdf1,
     +          dfthrsh,locut,hicut
      endif
c ................................................................


c ... 11 minus 4 micron BTDIF fog and low cloud test.
c ... for now placing in the SWIR bit place holder

      if (visusd .and. .not. snglnt) then
        if (nint(masir11) .ne. nint(bad_data)
     +                   .and.
     +      nint(masir4) .ne. nint(bad_data)) then
          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,19) 
          mas11_4 = masir11 - masir4
          if (mas11_4 .ge. pdo11_4lo(2)) then
            call set_bit(testbits,19)
            nptests = nptests + 1
          end if
          call conf_test(mas11_4,pdo11_4lo(1),pdo11_4lo(3),pdo11_4lo(4),
     +                  pdo11_4lo(2),1,c7)
          cmin2 = min(cmin2,c7)
          ngtests(2) = ngtests(2) + 1
        endif

c ...   debug statement ............................................
        if (debug .gt. 0) then
         write(h_output,'(1x,''mas11_4: '',5f10.2)')mas11_4,pdo11_4lo(1),
     +            pdo11_4lo(2),pdo11_4lo(3),pdo11_4lo(4)
        endif
c ................................................................

      end if
c *******     END OF GROUP 2 TESTS ****************************
c
c
c
c ********  START OF GROUP 3 TESTS ****************************
c ... visible (channel 2) reflectance threshold test.
      if (visusd) then
        if (nint(masv88) .ne. nint(bad_data)) then

c         Take into account sunglint problems
          if(snglnt) then
            call get_sg_thresholds(refang,locut,hicut,midpt,power)
          else
            locut = pdoref2(1)
            hicut = pdoref2(3)
            midpt = pdoref2(2)
            power = pdoref2(4)
          end if

          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,20) 
          if (masv88.le.midpt) then
            call set_bit(testbits,20)
            nptests = nptests + 1
          end if
          call conf_test(masv88,locut,hicut,power,
     +                   midpt,1,c8)
          cmin3 = min(cmin3,c8)
          ngtests(3) = ngtests(3) + 1
        endif

c ...   debug statement ............................................
        if (debug .gt. 0) then
          write(h_output,'(1x,''masv88: '',6f10.4)') masv88,locut,
     +            hicut,midpt,power,refang
        endif
c ................................................................
      end if
 

c ... visible channel ratio test (channel 2 / channel 1)
      if (visusd) then
        if (nint(masv66) .ne. nint(bad_data)
     +                   .and.
     +      nint(masv88) .ne. nint(bad_data)) then

c         Account for sun glint contamination
          if (snglnt) then
            locuta(1) = snglntvcl(1)
            locuta(2) = snglntvcl(2)
            hicuta(1) = snglntvch(1)
            hicuta(2) = snglntvch(2)
            midpta(1) = snglntv(1)
            midpta(2) = snglntv(2)
          else
            locuta(1) = pdovratlo(1)
            locuta(2) = pdovrathi(1)
            hicuta(1) = pdovratlo(3)
            hicuta(2) = pdovrathi(3)
            midpta(1) = pdovratlo(2)
            midpta(2) = pdovrathi(2)
          end if

          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,21) 
          vrat = masv88/masv66
          if (vrat .lt. midpta(1) .or. vrat .gt. midpta(2)) then
            call set_bit(testbits,21)
            nptests = nptests + 1
          end if
          call conf_test_2val(vrat,locuta,hicuta,1.0,midpta,2,c9)
          cmin3 = min(cmin3,c9)
          ngtests(3) = ngtests(3) + 1
        endif

c ...   debug statement ............................................
        if (debug .gt. 0) then
          write(h_output,'(1x,''vrat: '',7f10.2)') vrat,locuta(1),
     +            locuta(2),hicuta(1),hicuta(2),midpta(1),midpta(2)
        endif
c ................................................................
      end if
c ******       END OF GROUP 3 TESTS   ****************************
c
c
c
c
c ***********   START OF GROUP 4 TESTS  ************************* 
c ... near infrared high cloud test
      if (visusd) then
        if (nint(masv188) .ne. nint(bad_data)) then
          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,16) 
          if (masv188 .le. pdoref3(2)) then
            call set_bit(testbits,16)
            nptests = nptests + 1
          end if
          call conf_test(masv188,pdoref3(1),pdoref3(3),pdoref3(4),
     +                  pdoref3(2),1,c11)
          cmin4 = min(cmin4,c11)
          ngtests(4) = ngtests(4) + 1
        endif

c ...   debug statement ............................................
        if (debug .gt. 0) then
           write(h_output,'(1x,''masv188: '',6f10.4)')masv188,pdoref3(1),
     +                pdoref3(2),pdoref3(3),pdoref3(4)
        endif
c ................................................................
      end if
c ************   END OF GROUP 4 TESTS   ****************************
c
c
c     Check to see if thin cirrus bit should be set
      if (visusd) then
        if (nint(masv188) .ne. nint(bad_data)) then
          call set_qa_bit(qa_bits,9)
          if (masv188 .lt. pdotci(1) .and. masv188 .ge. pdotci(2)) then
            call clear_bit(testbits,9)
            cirrus_vis = .true.
          endif
c ...     debug statement ............................................
          if (debug .gt. 0) then
             write(h_output,'(1x,''NIR Thin cirrus: '',3f10.4)')masv188,
     +                  pdotci(1),pdotci(2)
          endif
c ..................................................................
        endif
      endif
c
c     Determine final confidence based on group values
      pre_confdnc = cmin1 * cmin2 * cmin3 * cmin4
      
c     Next, make sure you have all groups covered
      groups = 0
      do kk = 1,4
        if(ngtests(kk) .gt. 0) then
          groups = groups + 1.0
        end if
      enddo
      if (groups .gt. 0) fac = 1.0 / groups
c     Find final pixel confidence as nth root of group tests
      confdnc = pre_confdnc**fac

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(1x,''tests '',6i10)') nmtests,nptests,ngtests
        write(h_output,'(1x,''confdnc '',7f8.5/,8f8.5)') c1,c3,c4,
     +         c6,c7,c8,c9,c10,c11,cmin1,cmin2,cmin3,cmin4,fac,confdnc
      endif
c ................................................................

      return
      end
