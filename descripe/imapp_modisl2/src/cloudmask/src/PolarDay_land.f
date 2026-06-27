      subroutine PolarDay_land(pxldat,vza,visusd,vrused,cirrus_vis,
     +                         hi_elev,testbits,qa_bits,nmtests,
     +                         confdnc)
 
c---------------------------------------------------------------------
C!F77 
c
c!Description:
c      Routine for performing clear sky tests over polar land 
c      surfaces during daylight hours.
c
c      For daytime polar land the groups are:
c          Group 1: High thick cloud
c                   6.75 micron bt test 
c
c          Group 2: Low cloud - thick
c                   8-11 micron and 11-12 micron bt tests
c                   11-4 micron bt tests
c        
c          Group 3: Thick cloud
c                   .66 micron reflectance test (masv66)
c                   .87/.66 micron reflectance ratio test
c 
c          Group 4: Thin cirrus test
c                   1.38 micron reflectance test 
c
c!Input Parameters:
c pxldat        Array containing reflectance or brightness temperatures
c               for all bands for a single pixel
c vza           Current pixel viewing angle
c visusd        Logical variable indicating whether vis data used or not
c vrused        Logical variable defining when reflectance ratio
c               test can be used.
c cirrus_vis 	Logical variable flagging thin cirrus contaminated
c               scenes in the visible
c hi_elev       Logical flag indicating elevation > 2000 meters
c
c
c!Output Parameters:
c testbits      6 byte array containing cloud mask bit results
c qa_bits       10 byte array containing QA bit results
c confdnc       product of all applied individual confidences
c
c!Revision History:
c 06/04 Collection 5  R. Frey:
c Implemented new version of 11-12 um thin cirrus test (J. Key version)
c
c!Team-unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-35.
c
c!Design Notes:
c
c   Externals:
c      Subroutines: 
c      conftest,tview,set_qa_bit,set_bit,clear_bit
c
c!End
c---------------------------------------------------------------------

      implicit none
      save

      include 'global.inc'
      include 'PolarDay_land_thr.inc'


c ... scalar arguments ..
      real confdnc,vza
      integer nmtests
      logical visusd,vrused,cirrus_vis,hi_elev
c ...
c ... array arguments ..
      real pxldat(inband)
      byte testbits(6),qa_bits(10)
c ...
c ... local scalars ..
      real c2,c3,c4,c5,cosvza,dfthrsh,diftemp,dtr,mas11_4,masdf1,
     +     masir11,masir12,masir4,masv188,
     +     masv66,masv88,pi,schi,vrat,c6,
     +     masir65,c7,cmin1,cmin2,cmin3,cmin4,
     +     Rel_equality_EPS,groups,fac,pre_confdnc,
     +     eta,etad,etan,locut,hicut,s1,s2
      parameter (Rel_equality_EPS = 0.000001)
 
      integer nptests,debug,h_output,kk
c ...
c ... local arrays ..
      integer ngtests(4)
c ...
c ... external subroutines ..
      external conf_test,tview,set_bit,clear_bit,set_qa_bit

c     Common statement for debug purposes
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
      masv66 = pxldat(1)
      masv88 = pxldat(2)
      masv188 = pxldat(26)
      masir4 = pxldat(22)
      masir65 = pxldat(27)
      masir11 = pxldat(31)
      masir12 = pxldat(32)

c ... Initialization
      masdf1 = 0.0
      cosvza = 0.0
      schi = 0.0
      diftemp = 0.0
      dfthrsh = 0.0
      vrat = 0.0
      mas11_4 = 0.0

c ... the c suffix variables represent individual test confidences
      c2 = 0.0
      c3 = 0.0
      c4 = 0.0
      c5 = 0.0
      c6 = 0.0
      c7 = 0.0
c ... cmin are the group confidences
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
        write(h_output,'(10x/,''Processing subroutine LandDay '',/)') 
      endif
c ................................................................

 
C     **** GROUP 1 TESTS *************************************

c     H20 vapor channel (6.7 micron) high cloud test
      if (nint(masir65) .ne. nint(bad_data)) then
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,15)
        if (masir65 .gt. pdlh20(2)) then
          call set_bit(testbits,15)
          nptests = nptests + 1
        end if
        call conf_test(masir65,pdlh20(1),pdlh20(3),pdlh20(4),
     *                pdlh20(2),1,c2)
        cmin1 = min(cmin1,c2)
        ngtests(1) = ngtests(1) + 1
      endif

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(1x,''masir65: '',5f10.2)') masir65,pdlh20(1),
     +          pdlh20(2),pdlh20(3),pdlh20(4)
      endif
c ................................................................
c     *****  END OF GROUP 1 TESTS  ***************************
c
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
          dfthrsh = pdl11_12hi(1)
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
        locut = dfthrsh + (0.3 * dfthrsh)
        hicut = dfthrsh - (0.3 * dfthrsh)
        call conf_test(masdf1,locut,hicut,1.0,dfthrsh,1,c3)
        cmin2 = min(cmin2,c3)
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
      if (visusd) then
        if (nint(masir11) .ne. nint(bad_data) 
     +                   .and.
     +      nint(masir4) .ne.  nint(bad_data)) then
          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,19)
          mas11_4 = masir11 - masir4
          if (mas11_4 .ge. pdl11_4lo(2)) then
            call set_bit(testbits,19)
            nptests = nptests + 1
          end if
          call conf_test(mas11_4,pdl11_4lo(1),pdl11_4lo(3),pdl11_4lo(4),
     +                  pdl11_4lo(2),1,c4)
          cmin2 = min(cmin2,c4)
          ngtests(2) = ngtests(2) + 1
        endif

c ... debug statement ............................................
      if (debug .gt. 0) then
         write(h_output,'(1x,''mas11_4: '',5f10.2)')mas11_4,pdl11_4lo(1),
     +            pdl11_4lo(2),pdl11_4lo(3),pdl11_4lo(4)
      endif
c ................................................................
      end if
c *******     END OF GROUP 2 TESTS ****************************
c
c
c
c ********  START OF GROUP 3 TESTS ****************************
c ... visible (channel 1) reflectance threshold test.
      if (visusd) then
        if (nint(masv66) .ne. nint(bad_data)) then 
          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,20)
          if (masv66 .le. pdlref1(2)) then
            call set_bit(testbits,20)
            nptests = nptests + 1
          end if
          call conf_test(masv66,pdlref1(1),pdlref1(3),pdlref1(4),
     +                   pdlref1(2),1,c5)
          cmin3 = min(cmin3,c5)
          ngtests(3) = ngtests(3) + 1
        end if

c ...   debug statement ............................................
        if (debug .gt. 0) then
          write(h_output,'(1x,''masv66: '',5f10.2)') masv66,pdlref1(1),
     +            pdlref1(2),pdlref1(3),pdlref1(4)
        endif
c ................................................................
      end if

c ... visible channel ratio test (channel 2 / channel 1)
c ... Changed to implement GEMI test instead of straight ratio
c ... Apply only to scenes without sunglint and certain
c     ecosystem types
      if (visusd .and. vrused) then
        if (nint(masv66) .ne. nint(bad_data) 
     +                   .and.
     +      nint(masv88) .ne. nint(bad_data)) then
          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,21)
c ...     Scale values by 100 to make consistent with MAS version
          s1 = masv66 * 100.
          s2 = masv88 * 100.
          etan = 2.0 * (s2-s1) + 1.5*s2 + 0.5*s1
          etad = s2 + s1 + 0.5
          eta = etan / etad
          vrat=eta * (1.0-0.25*eta) - ((s1-0.125) / (1.0-s1))
          if(vrat .gt. pdlvrat(2)) then
            nptests = nptests + 1
            call set_bit(testbits,21)
          end if
          call conf_test(vrat,pdlvrat(1),pdlvrat(3),pdlvrat(4),
     +           pdlvrat(2),1,c6)
          cmin3 = min(cmin3,c6)
          ngtests(3) = ngtests(3) + 1
        end if

c ...   debug statement ............................................
        if (debug .gt. 0) then
           write(h_output,'(1x,''GEMI: '',7f10.2)')vrat,masv88,
     +                masv66,pdlvrat(1),pdlvrat(2),pdlvrat(3),pdlvrat(4)
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
      if (visusd .and. (.not. hi_elev) ) then
        if (nint(masv188) .ne. nint(bad_data)) then 
          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,16)
          if (masv188 .le. pdlref3(2)) then
            call set_bit(testbits,16)
            nptests = nptests + 1
          end if
          call conf_test(masv188,pdlref3(1),pdlref3(3),pdlref3(4),
     +                  pdlref3(2),1,c7)
          cmin4 = min(cmin4,c7)
          ngtests(4) = ngtests(4) + 1
        endif

c ...   debug statement ............................................
        if (debug .gt. 0) then
           write(h_output,'(1x,''masv188: '',6f10.4)')masv188,pdlref3(1),
     +                pdlref3(2),pdlref3(3),pdlref3(4)
        endif
c ................................................................
      end if
c ************   END OF GROUP 4 TESTS   ****************************
 
 
 
c     Check to see if thin cirrus bit should be set
      if (visusd .and. (.not. hi_elev) ) then
        if (nint(masv188) .ne. nint(bad_data)) then 
          call set_qa_bit(qa_bits,9)
          if (masv188 .lt. pdltci(1) .and. masv188 .ge. pdltci(2)) then
            call clear_bit(testbits,9)
            cirrus_vis = .true.
          endif
c ...     debug statement ............................................
          if (debug .gt. 0) then
             write(h_output,'(1x,''NIR Thin cirrus: '',3f10.4)')masv188,
     +                         pdltci(1),pdltci(2)
          endif
c ................................................................
        endif
      endif

c     Determine intermediate confidence based on group values
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
        write(h_output,'(1x,''confdnc '',7f8.5/,5f8.5)') c2,c3,c4,c5,
     +         c6,c7,cmin1,cmin2,cmin3,cmin4,fac,confdnc
      endif
c ................................................................

      return
      end
