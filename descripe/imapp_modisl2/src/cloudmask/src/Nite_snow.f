      subroutine Nite_snow(pxldat,vza,lnd,testbits,qa_bits,nmtests,
     +                     confdnc)

      implicit none
      save
 
c--------------------------------------------------------------------
C!F77 
c
c!Description:
c      Routine for performing clear sky tests over snow 
c      surfaces during nighttime hours.
c
c      The cloud test groups are:
c          Group 1: High thick cloud
c                   13.9 micron bt test (masir13) 
c                   6.75 micron bt test (not is use with mas)
c
c          Group 2: Low cloud - thick
c                   11-4 micron bt tests
c                   11-12 micron bt test
c                   7.3-11 micron bt test
c
c          Group 5: High cloud - thin
c                   3.7-12 micron bt test
c
c!Input Parameters:
c pxldat        Array containing reflectance or brightness temperatures
c               for all bands for a single pixel
c vza           Viewing zenith angle for current pixel
c lnd           Indicates land surface beneath snow
c
c!Output Parameters:
c testbits      six byte integer containing bit results
c qa_bits       ten byte array containing QA bit results
c nmtests       Acutal number of tests applied in this subroutine
c confdnc       product of all applied individual confidences
c
c!Revision History:
c 06/04 Collection 5  R. Frey:
c Added 11-12 um thin cirrus test (J. Key version)
c Added 7.3-11 thin cloud test (Y. Liu)
c
c!Team-unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c
c!Design Notes:
c    Externals:
c       Subroutines conf_test,set_bit,set_qa_bit,tview
c
c!END
c--------------------------------------------------------------------

      include 'global.inc'
      include 'Nite_snow_thr.inc'
      include 'PolarNite_snow_thr.inc'
c ...
c ... scalar arguments ..
      real confdnc,vza
      integer nmtests
      logical lnd
c ...
c ... array arguments ..
      real pxldat(inband)
      byte testbits(6),qa_bits(10)
c ...
c ... local scalars ..
      real c1,c2,mas4_12,masir11,masir12,masir13,masir4,
     +     c3,masir65,mas11_4,c4,cmin1,cmin2,cmin5,fac,
     +     pre_confdnc,groups,pi,dtr,masdf1,cosvza,schi,diftemp,
     +     dfthrsh,locut,hicut,c5,midpt,masir7,mas7_11,c6,power
      integer nptests,debug,h_output,kk

c     local arrays
      integer ngtests(3)

      real Rel_equality_EPS
      parameter(Rel_equality_EPS = 0.000001)

c ... external subroutines ..
      external conf_test,set_bit,set_qa_bit,tview

c ... Common statement for debug purposes
      common / bug / debug, h_output

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
      masir4 = pxldat(22)
      masir65 = pxldat(27)
      masir7 = pxldat(28)
      masir11 = pxldat(31)
      masir12 = pxldat(32)
      masir13 = pxldat(35)

      mas4_12 = 0.0
      mas11_4 = 0.0
      groups = 0.0
      pre_confdnc = 0.0
      fac = 0

c ... the c suffix variables represent individual test confidences
      c1 = 0.0
      c2 = 0.0
      c3 = 0.0
      c4 = 0.0
      c5 = 0.0
      c6 = 0.0
      cmin1 = 1.0
      cmin2 = 1.0
      cmin5 = 1.0

c ... initialize group number holder
      do 10 kk = 1 , 3
         ngtests(kk) = 0
  10  continue

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(10x/,''Processing subroutine Nite_snow '',
     +                   /)')
      endif
c ................................................................

C     **** GROUP 1 TESTS *************************************
c ... co2 high cloud test
      if (nint(masir13) .ne. nint(bad_data)) then
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,14)
        if (masir13 .gt. nsco2(2)) then
          call set_bit(testbits,14)
          nptests = nptests + 1
        end if
        call conf_test(masir13,nsco2(1),nsco2(3),nsco2(4),
     *                nsco2(2),1,c1)
        cmin1 = min(cmin1,c1)
        ngtests(1) = ngtests(1) + 1
      endif

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(1x,''masir13: '',5f10.2)') masir13,nsco2(1),
     +          nsco2(2),nsco2(3),nsco2(4)
      endif
c ................................................................


c     H20 vapor channel (6.7 micron) high cloud test
      if (nint(masir65) .ne. nint(bad_data)) then
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,15)
        if (masir65 .gt. nsh20(2)) then
          call set_bit(testbits,15)
          nptests = nptests + 1
        end if
        call conf_test(masir65,nsh20(1),nsh20(3),nsh20(4),
     *                nsh20(2),1,c2)
        cmin1 = min(cmin1,c2)
        ngtests(1) = ngtests(1) + 1
      endif

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(1x,''masir65: '',5f10.2)') masir65,nsh20(1),
     +          nsh20(2),nsh20(3),nsh20(4)
      endif
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
          dfthrsh = ns11_12hi(1)
        else
c         Add adjustment for snow cover.
          dfthrsh = diftemp + ns11_12adj(1)
        end if

c...    Set flags if test passed
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,18)
        if (masdf1 .le. dfthrsh) then
          call set_bit(testbits,18)
          nptests = nptests + 1
        end if

        locut = dfthrsh
        midpt = dfthrsh - (0.3 * dfthrsh)
        if(masir11 .lt. 270.0) then
          hicut = midpt - (0.2 * dfthrsh)
        else
          hicut = midpt - 1.25
        end if

        call conf_test(masdf1,locut,hicut,1.0,dfthrsh,1,c5)
        cmin2 = min(cmin2,c5)
        ngtests(2) = ngtests(2) + 1

      endif

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(1x,''APOLLO masdf1: '',8f10.2)') masdf1,
     +          ns11_12hi(1),ns11_12adj(1),
     +          masir11,schi,dfthrsh,locut,hicut
      endif
c ................................................................

c ... 11 minus 4 micron BTDIF fog and low cloud test.
      if (nint(masir11) .ne. nint(bad_data)
     +                 .and.
     +    nint(masir4) .ne.  nint(bad_data)) then
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,19)
        mas11_4 = masir11 - masir4
        if (mas11_4 .le. ns11_4lo(2)) then
          call set_bit(testbits,19)
          nptests = nptests + 1
        end if
        call conf_test(mas11_4,ns11_4lo(1),ns11_4lo(3),ns11_4lo(4),
     +                ns11_4lo(2),1,c3)
        cmin2 = min(cmin2,c3)
        ngtests(2) = ngtests(2) + 1
      endif

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(1x,''mas11_4: '',5f10.2)')mas11_4,ns11_4lo(1),
     +            ns11_4lo(2),ns11_4lo(3),ns11_4lo(4)
      endif
c ................................................................


c ... 7.3 minus 11 micron cloud test
      if (nint(masir11) .ne. nint(bad_data)
     +                  .and.
     +    nint(masir7) .ne.  nint(bad_data)) then

        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,23)
        mas7_11 = masir7 - masir11

        if(lnd) then
c         Get thresholds for land (snow) surface.
          call get_pn_thresholds(masir11,bt_11_bnds2,pn_7_11l,pn_7_11m1,
     *                          pn_7_11m2,pn_7_11m3,pn_7_11h,locut,hicut,
     *                          midpt,power)
        else
c         Get thresholds for water (ice) surface.
          call get_pn_thresholds(masir11,bt_11_bnds2,pn_7_11lw,pn_7_11m1w,
     *                          pn_7_11m2w,pn_7_11m3w,pn_7_11hw,locut,hicut,
     *                          midpt,power)
        end if

        if (mas7_11 .gt. midpt) then
          call set_bit(testbits,23)
          nptests = nptests + 1
        end if
        call conf_test(mas7_11,locut,hicut,power,midpt,1,c6)
        cmin2 = min(cmin2,c6)
        ngtests(2) = ngtests(2) + 1
      endif

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(1x,''mas7_11: '',7f10.2)')mas7_11,locut,
     +            hicut,midpt,power,masir11,c6
      end if
c ................................................................
c *******     END OF GROUP 2 TESTS ****************************


c *******    START OF GROUP 5 TESTS  **************************
c ... 4-12um brightness temperature difference test
c ... for thin cirrus).
      if (nint(masir12) .ne. nint(bad_data)
     +                 .and.
     +    nint(masir4) .ne.  nint(bad_data)) then
        mas4_12 = masir4 - masir12
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,17)
        if (mas4_12 .le. ns4_12hi(2)) then
          nptests = nptests + 1
          call set_bit(testbits,17)
        end if
        call conf_test(mas4_12,ns4_12hi(1),ns4_12hi(3),ns4_12hi(4),
     +                ns4_12hi(2),1,c4)
        cmin5 = min(cmin5,c4)
        ngtests(3) = ngtests(3) + 1
      endif

c ... debug statement ............................................
      if (debug .gt. 0) then
      write(h_output,'(1x,''mas4_12: '',5f10.2)')mas4_12,ns4_12hi(1),
     +            ns4_12hi(2),ns4_12hi(3),ns4_12hi(4)
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
      if(groups .gt. 0) fac = 1.0 / groups
c     Find final pixel confidence as nth root of group tests
      confdnc = pre_confdnc**fac


c     One last test.  If the 6.5 micron brightness temperature is
c     greater than the 11 micron value, then there is an inversion
c     and is very, very likely to be clear
      if (nint(masir11) .ne. nint(bad_data)
     +                 .and.
     +    nint(masir65) .ne.  nint(bad_data)) then
        call set_qa_bit(qa_bits,26)
        if ((masir65 - masir11) .gt. n65_11(1)) then 
              confdnc = 1.0
              call set_bit(testbits,26)
        endif
c ...   debug statement ............................................
        if (debug .gt. 0) then
           write(h_output,'(10x,''Final inversion test:'',/,4f10.2,/)')
     +                masir65,masir11,masir11-masir65,confdnc
        endif
c ...................................................................
      endif


c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(1x,''tests '',6i10)') nmtests,nptests,ngtests
        write(h_output,'(1x,''confdnc '',9f8.5/,2f8.5)') c1,c2,c3,c4,
     +         c5,c6,cmin1,cmin2,cmin5,fac,confdnc
      endif
c ................................................................


      return
      end
