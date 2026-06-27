      subroutine PolarNite_land(pxldat,vza,desert,hi_elev,sfctmp,
     +                          eco_type,testbits,qa_bits,nmtests,
     +                          confdnc)

      implicit none
      save

c---------------------------------------------------------------------
c!F77 
c
c!Description:
c      Routine for performing clear sky tests over polar land
c      surfaces during nightime hours.
c
c      For nighttime polar land the groups are:
c          Group 1: High thick cloud
c                   6.75 micron bt test 
c                   surface temperature test
c
c          Group 2: Low cloud - thick
c                   11-12 micron bt test
c                   11-4 micron bt test
c                   7.3-11 micron bt test
c
c          Group 5: High cloud - thin
c                   3.7-12 micron bt test
c
c!Input Parameters:
c pxldat        Array containing reflectance or brightness temperatures
c               for all bands for a single pixel
c vza           viewing zenith angle
c desert        flag indicating desert processing
c hi_elev       flag indicating high elevation processing (> 2000 m)
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
c 06/04 Collection 5  R. Frey
c Added surface temperature cloud test (GDAS sfc air vs 11 micron bt)
c Added 7.3-11 micron cloud test
c Implemented new 11-12 micron thin cirrus test (Key version)
c Modified 3.9-12 micron test (dynamic thresholds based on 11 micron bt)
c 10/04 Collection 5  R. Frey
c Added 11-12 and 3.9-11 um BTD conditions on choice of LST threshold
c Changed basic LST threshold from 10K to 12K
c
c!Team-Unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-35.
c
c!Design Notes:
c
c Externals:
c       Subroutines: conf_test,set_bit,set_qa_bit,tview,get_pn_thresholds
c
c!END
c-----------------------------------------------------------------------

      include 'global.inc'
      include 'PolarNite_land_thr.inc'
      include 'PolarNite_snow_thr.inc'
c ...
c ... scalar arguments ..
      real confdnc,vza,sfctmp
      integer nmtests
      logical desert,hi_elev
      byte eco_type
c ...
c ... array arguments ..
      real pxldat(inband)
      byte testbits(6),qa_bits(10)
c ...
c ... local scalars ..
      real c2,mas4_12,masir11,masir12,masir4,masir7,mas7_11,
     +     c3,masir65,mas11_4,c4,cmin1,cmin2,cmin5,groups,
     +     fac,pre_confdnc,c5,c6,masdf1,schi,cosvza,dtr,pi,diftemp,
     +     dfthrsh,locut,hicut,Rel_equality_EPS,midpt,power,a,
     +     lst_thrsh,corr,max_vza,sfcdif,c7,masdf2
      integer nptests,debug,h_output,kk

c ... local arrays
      integer ngtests(3)
c ...
      parameter(Rel_equality_EPS = 0.000001)
      parameter(max_vza = 65.49)

c ... external subroutines ..
      external conf_test,set_bit,set_qa_bit,tview,get_pn_thresholds

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
c ... confidence to 1.0 to begin with
      confdnc = 1.0

c ... place band values into individual variables for easy
c ... identification
      masir4 = pxldat(22)
      masir65 = pxldat(27)
      masir7 = pxldat(28)
      masir11 = pxldat(31)
      masir12 = pxldat(32)

c ...
      mas4_12 = 0.0
      mas11_4 = 0.0
      masdf1 = 0.0
      schi = 0.0

c ... the c suffix variables represent individual test confidences
      c2 = 0.0
      c3 = 0.0
      c4 = 0.0
      c5 = 0.0
      c6 = 0.0
      c7 = 0.0
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

c     H20 vapor channel (6.7 micron) high cloud test
      if (nint(masir65) .ne. nint(bad_data)) then
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,15)
        if (masir65 .gt. pnlh20(2)) then
          call set_bit(testbits,15)
          nptests = nptests + 1
        end if
        call conf_test(masir65,pnlh20(1),pnlh20(3),pnlh20(4),
     *                pnlh20(2),1,c2)
        cmin1 = min(cmin1,c2)
        ngtests(1) = ngtests(1) + 1
      endif

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(1x,''masir65: '',5f10.2)') masir65,pnlh20(1),
     +          pnlh20(2),pnlh20(3),pnlh20(4)
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
        else if(masdf1 .ge. 0.0 .or.
     *  (masdf1 .lt. 0.0 .and. (masdf2 .le. -0.5 .or. masdf2 .ge. 1.0))) then
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

        call conf_test(sfcdif,locut,hicut,1.0,midpt,1,c7)
        cmin1 = min(cmin1,c7)
        ngtests(1) = ngtests(1) + 1

       endif

      endif

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(1x,''sfctmp: '',9f9.3)') masdf1,masdf2,
     +   sfctmp,masir11,sfcdif,locut,midpt,hicut,c7
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
          dfthrsh = pnl11_12hi(1)
        else
c         Add 0.2 for likely snow cover.
          dfthrsh = diftemp + 0.2
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
        write(h_output,'(1x,''APOLLO masdf1: '',5f10.2)') masdf1,
     +          pnl11_12hi(1),dfthrsh,locut,hicut
      endif
c ................................................................

c ... 11 minus 4 micron BTDIF fog and low cloud test.
      if (nint(masir11) .ne. nint(bad_data)
     +                  .and.
     +     nint(masir4) .ne.  nint(bad_data)) then

        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,19)
        mas11_4 = masir11 - masir4

c       Use polar night snow thresholds here.  Logic is that land
c       surfaces poleward of 60 latitude are snow covered most of
c       the year and many times ancillary snow map is incomplete.

        call get_pn_thresholds(masir11,bt_11_bounds,pn_11_4l,pn_11_4m1,
     *                         pn_11_4m2,pn_11_4m3,pn_11_4h,locut,hicut,
     *                         midpt,power)

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
        write(h_output,'(1x,''mas11_4: '',6f10.2)')mas11_4,locut,
     +            hicut,midpt,power,masir11
      end if
c ................................................................

c ... 7.3 minus 11 micron cloud test 
      if (nint(masir11) .ne. nint(bad_data)
     +                  .and.
     +    nint(masir7) .ne.  nint(bad_data)) then

c       Check 11 um brightness temperature.  This to guard against
c       false cloud retrievals during polar summer.
        if(masir11 .lt. 270.0) then

          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,23)
          mas7_11 = masir7 - masir11

c         Use polar night snow thresholds here.  Logic is that land
c         surfaces poleward of 60 latitude are snow covered most of
c         the year and many times ancillary snow map is incomplete.

          call get_pn_thresholds(masir11,bt_11_bnds2,pn_7_11l,pn_7_11m1,
     *                           pn_7_11m2,pn_7_11m3,pn_7_11h,locut,hicut,
     *                           midpt,power)

          if (mas7_11 .gt. midpt) then
            call set_bit(testbits,23)
            nptests = nptests + 1
          end if
          call conf_test(mas7_11,locut,hicut,power,midpt,1,c6)
          cmin2 = min(cmin2,c6)
          ngtests(2) = ngtests(2) + 1

        end if
      endif

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(1x,''mas7_11: '',6f10.2)')mas7_11,locut,
     +            hicut,midpt,power,masir11
      end if
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

c       Use polar night snow thresholds here.  Logic is that land
c       surfaces poleward of 60 latitude are snow covered most of
c       the year and many times ancillary snow map is incomplete.

        call get_pn_thresholds(masir11,bt_11_bounds,pn_4_12l,pn_4_12m1,
     *                         pn_4_12m2,pn_4_12m3,pn_4_12h,locut,hicut,
     *                         midpt,power)

        if (mas4_12 .le. midpt) then
          nptests = nptests + 1
          call set_bit(testbits,17)
        end if
        call conf_test(mas4_12,locut,hicut,power,midpt,1,c4)
        cmin5 = min(cmin5,c4)
        ngtests(3) = ngtests(3) + 1
      endif

c ... debug statement ............................................
      if (debug .gt. 0) then
      write(h_output,'(1x,''mas4_12: '',6f10.2)')mas4_12,locut,
     +            hicut,midpt,power,masir11
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
        write(h_output,'(1x,''confdnc '',8f8.5/,2f8.5)') c2,c3,c4,c5,
     +         c7,cmin1,cmin2,cmin5,fac,confdnc
      endif
c ................................................................

      return
      end
