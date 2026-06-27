      subroutine LandNite(pxldat,plat,vza,coast,desert,hi_elev,sh_lake,
     +                    sfctmp,eco_type,testbits,qa_bits,nmtests,
     +                    confdnc)

      implicit none
      save

c---------------------------------------------------------------------
c!F77 
c
c!Description:
c      Routine for performing clear sky tests over land
c      surfaces during nightime hours.
c
c      For nighttime land the groups are:
c          Group 1: High thick cloud
c                   13.9 micron bt test (masir13) 
c                   6.75 micron bt test 
c                   Surface temperature test 
c
c          Group 2: Low cloud - thick
c                   11-12 micron bt tests
c                   11-4 micron bt tests
c                   7.3-11 micron bt test (thick mid-level)
c
c          Group 5: High cloud - thin
c                   3.7-12 micron bt test
c
c!Input Parameters:
c pxldat        Array containing reflectance or brightness temperatures
c               for all bands for a single pixel
c vza           viewing zenith angle
c coast         flag indicating coast processing
c desert        flag indicating desert processing
c hi_elev       flag indicating high elevation processing (> 2000 m)
c sh_lake       Logical flag indicating shallow inland lakes
c sfctmp        Surface air temperature from model data
c eco_type      Ecosystem index
c
c!Output Parameters:
c testbits      four-byte integer containing bit results
c qa_bits       10 byte array containing QA bit results
c nmtests       Number of tests actually applied in this routine
c confdnc       product of all applied individual confidences
c
c!Revision History:
c
c 07-15-02 raf:
c Modified 3.9-11 um test so that the test threshold is a function of
c the 11-12 um BTD.  Thresholds are calculated in get_nl_thresholds.f
c Added the 7.3-11 um mid-level cloud test.
c Added 2K to 3.9-11 um test threshold for coastal areas.
c 06/04 Collection 5   raf:
c Added surface temperature test    
c Implemented new 11-12 um thin cirrus test (J. Key version)
c 10/04 Collection 5   raf:
c Added 11-12 and 3.9-11 um BTD conditions to choice of LST thresholds.
c Changed basic LST threshold from 10K to 12K.
c
c!Team-Unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c
c!Design Notes:
c
c Externals:
c       Subroutines: conf_test,set_bit,set_qa_bit,tview,get_nl_thresholds
c
c!END
c-----------------------------------------------------------------------

      include 'global.inc'
      include 'LandNite_thr.inc'
c ...
c ... scalar arguments ..
      real confdnc,vza,plat,sfctmp
      integer nmtests
      logical coast,desert,hi_elev,sh_lake
      byte eco_type
c ...
c ... array arguments ..
      real pxldat(inband)
      byte testbits(6),qa_bits(10)
c ...
c ... local scalars ..
      real c1,c2,mas4_12,masir11,masir12,masir13,masir4,
     +     c3,masir65,mas11_4,c4,cmin1,cmin2,cmin5,groups,
     +     fac,pre_confdnc,c5,masdf1,schi,cosvza,dtr,pi,diftemp,
     +     dfthrsh,locut,hicut,Rel_equality_EPS,masir73,mas7_11,c6,
     +     power,midpt,a,c9,sfcdif,max_vza,corr,lst_thrsh,masdf2
      integer nptests,debug,h_output,kk

c ... local arrays
      integer ngtests(3)
c ...
      parameter(Rel_equality_EPS = 0.000001)
      parameter(max_vza = 65.49)

c ... external subroutines ..
      external conf_test,set_bit,set_qa_bit,tview,get_nl_thresholds

c ... Common statement for debug purposes
      common / bug / debug, h_output

c ... initialize variables
      pi = acos(-1.0)
      dtr = pi/180.0
c ... nmtests counts the number of tests applied to this pixel
      nmtests = 0
c ... nptests counts the number of tests passed
      nptests = 0
c ... confidence to 1.0 to begin with
      confdnc = 1.0

c ... place band values into individual variables for easy
c ... identification
      masir4 = pxldat(22)
      masir65 = pxldat(27)
      masir73 = pxldat(28)
      masir11 = pxldat(31)
      masir12 = pxldat(32)
      masir13 = pxldat(35)

c ...
      mas4_12 = 0.0
      mas11_4 = 0.0
      masdf1 = 0.0
      schi = 0.0

c ... the c suffix variables represent individual test confidences
      c1 = 0.0
      c2 = 0.0
      c3 = 0.0
      c4 = 0.0
      c5 = 0.0
      c6 = 0.0
      c9 = 0.0
      cmin1 = 1.0
      cmin2 = 1.0
      cmin5 = 1.0

c ... initialize group number holder
      do 10 kk = 1 , 3
         ngtests(kk) = 0
  10  continue

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(10x/,''Processing subroutine LandNite '',/)')
      endif
c ................................................................


C     **** GROUP 1 TESTS *************************************
c ... co2 high cloud test
      if (nint(masir13) .ne. nint(bad_data)) then
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,14)
        if (masir13.gt.nlco2(2)) then
          call set_bit(testbits,14)
          nptests = nptests + 1
        end if
        call conf_test(masir13,nlco2(1),nlco2(3),nlco2(4),
     *                nlco2(2),1,c1)
        cmin1 = min(cmin1,c1)
        ngtests(1) = ngtests(1) + 1
      endif

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(1x,''masir13: '',5f10.2)') masir13,nlco2(1),
     +          nlco2(2),nlco2(3),nlco2(4)
      endif
c ................................................................

c     H20 vapor channel (6.7 micron) high cloud test
      if (nint(masir65) .ne. nint(bad_data)) then
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,15)
        if (masir65.gt.nlh20(2)) then
          call set_bit(testbits,15)
          nptests = nptests + 1
        end if
        call conf_test(masir65,nlh20(1),nlh20(3),nlh20(4),
     *                nlh20(2),1,c2)
        cmin1 = min(cmin1,c2)
        ngtests(1) = ngtests(1) + 1
      endif

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(1x,''masir65: '',5f10.2)') masir65,nlh20(1),
     +          nlh20(2),nlh20(3),nlh20(4)
      endif

c ................................................................

c ... Surface Temperature Test

      if ( nint(masir11) .ne. nint(bad_data) .and. (.not. hi_elev) .and.
     *     nint(masir12) .ne. nint(bad_data) .and. eco_type .ne. 8) then

       if (sfctmp .gt. 0.0 .and. sfctmp .lt. 350.0) then

        masdf1 = masir11 - masir12
        masdf2 = masir11 - masir4

        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,27)

        if(desert) then
          lst_thrsh = 20.0
        else if(masdf1 .ge. -0.2 .or.
     *  (masdf1 .lt. -0.2 .and. (masdf2 .le. -0.5 .or. masdf2 .ge. 1.0))) then
          lst_thrsh = 12.0
        else
          lst_thrsh = 20.0
        end if

        if(masdf1 .ge. 1.0) then
          midpt = lst_thrsh + (2.0 * nint(masdf1))
        else
          midpt = lst_thrsh 
        end if
        a = vza / max_vza
        corr = (a**4) * 3.0
        midpt = midpt + corr
        locut = midpt + 2.0
        hicut = midpt - 2.0

        sfcdif = sfctmp - masir11

        if( sfcdif .lt. midpt ) then
          call set_bit(testbits,27)
          nptests = nptests + 1
        end if

        call conf_test(sfcdif,locut,hicut,1.0,midpt,1,c9)
        cmin1 = min(cmin1,c9)
        ngtests(1) = ngtests(1) + 1

       endif

      endif

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(1x,''sfctmp: '',9f9.3)') masdf1,masdf2,
     +   sfctmp,masir11,sfcdif,locut,midpt,hicut,c9
      end if

c ................................................................
c     *****  END OF GROUP 1 TESTS  ***************************
c
c 
c     ****  GROUP 2 TESTS  ***********************************
c ... 11-12um brightness temperature difference test
c ... for thin cirrus).
      if (nint(masir11) .ne. nint(bad_data) 
     +                 .and.
     +    nint(masir12) .ne. nint(bad_data)
     +                  .and.
     +              vza .gt. 0.0) then
      
        masdf1 = masir11 - masir12
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
          dfthrsh = nl11_12hi(1)
        else
          dfthrsh = diftemp
        end if

c...    Set flags if test passed
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,18)
        if (masdf1.le.dfthrsh) then
          call set_bit(testbits,18)
          nptests = nptests + 1
        end if

        locut = dfthrsh
        midpt = dfthrsh - (0.3 * dfthrsh)
        if(masir11 .lt. 270.0) then
          if(abs(plat) .le. 30.0) then
            hicut = midpt - 1.25
          else
            a = (90.0 - abs(plat)) / 60.0
            hicut = -0.1 - ((a**4) * 1.15)
          end if
        else
          hicut = midpt - 1.25
        end if

        call conf_test(masdf1,locut,hicut,1.0,midpt,1,c5)
        cmin2 = min(cmin2,c5)
        ngtests(2) = ngtests(2) + 1

      endif

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(1x,''APOLLO masdf1: '',5f10.2)') masdf1,
     +          nl11_12hi(1),dfthrsh,locut,hicut
      endif
c ................................................................

c ... 11 minus 4 micron BTDIF fog and low cloud test.
      if (nint(masir11) .ne. nint(bad_data) .and.
     +     nint(masir4) .ne.  nint(bad_data) .and.
     +     nint(masir12) .ne. nint(bad_data)) then

        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,19)
        mas11_4 = masir11 - masir4

        masdf1 = masir11 - masir12
        call get_nl_thresholds(masdf1,locut,hicut,midpt,power)

        if(sh_lake) then
          locut = locut + 2.0
          midpt = midpt + 2.0
          hicut = hicut + 2.0
        end if

        if (mas11_4 .le. midpt) then
          call set_bit(testbits,19)
          nptests = nptests + 1
        end if

        call conf_test(mas11_4,locut,hicut,power,midpt,1,c3)

        cmin2 = min(cmin2,c3)
        ngtests(2) = ngtests(2) + 1

      endif

c ... debug statement ............................................
      if (debug .gt. 0) then
         write(h_output,'(1x,''mas11_4: '',l4,6f10.2)')coast,mas11_4,
     +       masdf1,locut,midpt,hicut,c3
      endif
c ................................................................


c ... 7.3-11um brightness temperature difference test
c ... for thick, mid-level clouds
      if (nint(masir11) .ne. nint(bad_data)
     +                  .and.
     +     nint(masir73) .ne.  nint(bad_data) .and.
     +     nint(masir4) .ne. nint(bad_data)) then

        mas11_4 = masir11 - masir4

        if(mas11_4 .le. -2.0) then

          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,23)
          mas7_11 = masir73 - masir11

          if ( mas7_11 .le. nl7_11s(2) ) then
            nptests = nptests + 1
            call set_bit(testbits,23)
          end if

          call conf_test(mas7_11,nl7_11s(1),nl7_11s(3),nl7_11s(4),
     +                   nl7_11s(2),1,c6)

          cmin2 = min(cmin2,c6)
          ngtests(2) = ngtests(2) + 1

        end if

      endif

c ... debug statement ............................................
      if (debug .gt. 0) then
         write(h_output,'(1x,''mas7_11: '',6f10.2)')mas7_11,masir11,
     +            nl7_11s(1),nl7_11s(2),nl7_11s(3),c6
      endif
c ................................................................


c *******     END OF GROUP 2 TESTS ****************************



c *******    START OF GROUP 5 TESTS  **************************
c ... 4-12um brightness temperature difference test
c ... for thin cirrus)
      if (nint(masir12) .ne. nint(bad_data)
     +                  .and.
     +     nint(masir4) .ne.  nint(bad_data)) then
        mas4_12 = masir4 - masir12
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,17)
        if (mas4_12.le.nl4_12hi(2)) then
          nptests = nptests + 1
          call set_bit(testbits,17)
        end if
        call conf_test(mas4_12,nl4_12hi(1),nl4_12hi(3),nl4_12hi(4),
     +                nl4_12hi(2),1,c4)
        cmin5 = min(cmin5,c4)
        ngtests(3) = ngtests(3) + 1
      endif

c ... debug statement ............................................
      if (debug .gt. 0) then
         write(h_output,'(1x,''mas4_12: '',5f10.2)')mas4_12,nl4_12hi(1),
     +            nl4_12hi(2),nl4_12hi(3),nl4_12hi(4)
      endif
c ................................................................
c ********    END OF GROUP 5 TESTS  *****************************

c     Determine final confidence based on group values
      pre_confdnc = cmin1 * cmin2 * cmin5

c     Next, make sure you have all groups covered
      groups = 0
      do kk = 1,3
        if(ngtests(kk) .gt. 0) then
          groups = groups + 1.0
        end if
      enddo
      if (groups .gt. 0) fac = 1.0 / groups
c     Find final pixel confidence as nth root of group tests
      confdnc = pre_confdnc**fac


c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(1x,''tests '',5i10)') nmtests,nptests,ngtests
        write(h_output,'(1x,''confdnc '',10f8.5/,f8.5)') c1,c2,c3,c4,c5,
     +      c6,cmin1,cmin2,cmin5,fac,confdnc
      endif
c ................................................................

      return
      end
