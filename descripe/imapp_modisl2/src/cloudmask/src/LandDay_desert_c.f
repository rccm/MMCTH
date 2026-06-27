      subroutine LandDay_desert_c(pxldat,vza,visusd,cirrus_vis,
     +                            hi_elev,testbits,qa_bits,
     +                            nmtests,confdnc)
 
      implicit none
      save

c--------------------------------------------------------------------
C!F77 
c
c!Description:
c
c      Routine for performing clear sky tests over coastal desert 
c      surfaces during daylight hours.
c
c      For daytime coastal desert the groups are:
c          Group 1: High thick cloud
c                   13.9 micron bt test (masir13) 
c                   6.75 micron bt test 
c
c          Group 2: Mainly Low cloud - thick
c                   11-4 micron bt tests
c                   11-12 Thin cirrus tests
c        
c          Group 3: Thick cloud
c                   .87 micron reflectance test (masv87)
c 
c          Group 4: Thin cirrus test
c                   1.38 micron reflectance test 
c
c!Input Parameters:
c real    pxldat Array containing reflectance or brightness temperatures
c                for all bands for a single pixel
c real    vza    Current pixel viewing angle
c logical visusd Logical variable indicating whether vis data used 
c logical cirrus_vis Logical variable flagging thin cirrus contaminated
c               scenes in the visible
c logical hi_elev Logical variable indicating elevation > 2000 meters
c
c!Output Parameters:
c byte    testbits six byte integer containing bit results
c byte    qa_bits  ten byte array containing QA bit results
c integer nmtests  Counts number of inidividual tests applied
c real    confdnc  product of all applied individual confidences
c
c!Revision History:
c Implemented new version of 11-12 um thin cirrus test (J. Key version)
c 06/04 Collection 5   R. Frey
c
c!Team-Unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c 1996/4/3 13:55:25
c K. Strabala (kathys@ssec.wisc.edu)
c Original version.
c
c!Design Notes:
c
c   Externals:
c      Subroutines:  
c      conf_test_2val,conftest,tview,set_bit,clear_bit
c      set_qa_bit
c
c!End
c----------------------------------------------------------------------

      include 'global.inc'
      include 'LandDay_desert_c_thr.inc'

c ... scalar arguments ..
      real confdnc,vza
      integer nmtests
      logical visusd,cirrus_vis,hi_elev
c ...
c ... array arguments ..
      real pxldat(inband)
      byte testbits(6),qa_bits(10)
c ...
c ... local scalars ..
      real c1,c2,c3,c4,c5,cosvza,dfthrsh,diftemp,dtr,mas11_4,masdf1,
     +     masir11,masir12,masir13,masir4,masv188,masv88,pi,schi,
     +     masir65,c6,cmin1,cmin2,cmin3,cmin4,pre_confdnc,
     +     Rel_equality_EPS,groups,fac,locut,hicut
      parameter (Rel_equality_EPS = 0.000001)

      integer nptests,kk,debug,h_output
c ...
c ... local arrays ..
      real hicuta(2),locuta(2),midpta(2)
      integer ngtests(4)
c ...
c ... external subroutines ..
      external conf_test,conf_test_2val,tview,set_bit,
     *         clear_bit,set_qa_bit

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
      masv88 = pxldat(2)
      masv188 = pxldat(26)
      masir4 = pxldat(22)
      masir65 = pxldat(27)
      masir11 = pxldat(31)
      masir12 = pxldat(32)
      masir13 = pxldat(35)

c ...
      masdf1 = 0.0
      cosvza = 0.0
      schi = 0.0
      diftemp = 0.0
      dfthrsh = 0.0
      mas11_4 = 0.0

c ... the c suffix variables represent individual test confidences
      c1 = 0.0
      c2 = 0.0
      c3 = 0.0
      c4 = 0.0
      c5 = 0.0
      c6 = 0.0
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
        write(h_output,'(10x/,''Subroutine LandDay_desert_c '',
     +                      /)')
      endif
c ................................................................

 
c ... perform tests.  note that some tests are not used
c ... in sun glint conditions.
c
C     **** GROUP 1 TESTS *************************************
c ... co2 high cloud test
      if (nint(masir13) .ne. nint(bad_data)) then
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,14)
        if (masir13 .gt. ldsco2_c(2)) then
          call set_bit(testbits,14)
          nptests = nptests + 1
        end if
        call conf_test(masir13,ldsco2_c(1),ldsco2_c(3),ldsco2_c(4),
     *                ldsco2_c(2),1,c1)
        cmin1 = min(cmin1,c1)
        ngtests(1) = ngtests(1) + 1
      endif

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(1x,''masir13: '',5f10.2)') masir13,ldsco2_c(1),
     +          ldsco2_c(2),ldsco2_c(3),ldsco2_c(4)
      endif
c ................................................................


c     H20 vapor channel (6.7 micron) high cloud test
      if (nint(masir65) .ne. nint(bad_data)) then
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,15)
        if (masir65 .gt. ldsh20_c(2)) then
          call set_bit(testbits,15)
          nptests = nptests + 1
        end if
        call conf_test(masir65,ldsh20_c(1),ldsh20_c(3),ldsh20_c(4),
     *                 ldsh20_c(2),1,c2)
        cmin1 = min(cmin1,c2)
        ngtests(1) = ngtests(1) + 1
      endif

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(1x,''masir65: '',5f10.2)') masir65,ldsh20_c(1),
     +          ldsh20_c(2),ldsh20_c(3),ldsh20_c(4)
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
          dfthrsh = lds11_12hi_c(2)
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
        hicut = dfthrsh - 1.25
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
      if (visusd) then
        if (nint(masir11) .ne. nint(bad_data)
     +                   .and.
     +      nint(masir4) .ne.  nint(bad_data)) then
          if(masir11 .le. 320.0) then
            nmtests = nmtests + 1
            mas11_4 = masir11 - masir4
            call set_qa_bit(qa_bits,19)
            if (mas11_4 .ge. lds11_4lo_c(2) .and. 
     +          mas11_4 .le. lds11_4hi_c(2)) then
              call set_bit(testbits,19)
              nptests = nptests + 1
            end if
            locuta(1) = lds11_4lo_c(1)
            locuta(2) = lds11_4hi_c(1)
            hicuta(1) = lds11_4lo_c(3)
            hicuta(2) = lds11_4hi_c(3)
            midpta(1) = lds11_4lo_c(2)
            midpta(2) = lds11_4hi_c(2)
  
            call conf_test_2val(mas11_4,locuta,hicuta,1.0,midpta,2,c4)
            cmin2 = min(cmin2,c4)
            ngtests(2) = ngtests(2) + 1
          endif
        end if

c ..... debug statement ............................................
        if (debug .gt. 0) then
        write(h_output,'(1x,''mas11_4a: '',5f10.2)')mas11_4,
     +        lds11_4lo_c(1),lds11_4lo_c(2),lds11_4lo_c(3),lds11_4lo_c(4)
        write(h_output,'(1x,''mas11_4b: '',5f10.2)')mas11_4,
     +        lds11_4hi_c(1),lds11_4hi_c(2),lds11_4hi_c(3),lds11_4hi_c(4)
        endif
c ..................................................................
      endif


c *******     END OF GROUP 2 TESTS ****************************
c
c
c
c ********  START OF GROUP 3 TESTS ****************************
c ... visible (.88 micron) reflectance threshold test.
      if (visusd) then
        if (nint(masv88) .ne. nint(bad_data)) then
          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,20)
          if (masv88 .le. ldsref2_c(2)) then
            call set_bit(testbits,20)
            nptests = nptests + 1
          end if
          call conf_test(masv88,ldsref2_c(1),ldsref2_c(3),ldsref2_c(4),
     +                  ldsref2_c(2),1,c5)
          cmin3 = min(cmin3,c5)
          ngtests(3) = ngtests(3) + 1
        end if

c  ...  debug statement ............................................
        if (debug .gt. 0) then
          write(h_output,'(1x,''masv88: '',5f10.2)') masv88,ldsref2_c(1),
     +            ldsref2_c(2),ldsref2_c(3),ldsref2_c(4)
        endif
c ..................................................................
      endif
c ******       END OF GROUP 3 TESTS   ****************************
 
 
 
 
c ***********   START OF GROUP 4 TESTS  *************************
c ... near infrared high cloud test
      if ((.not. hi_elev) .and. visusd) then
        if (nint(masv188) .ne. nint(bad_data)) then
          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,16)
          if (masv188 .le. ldsref3_c(2)) then
            call set_bit(testbits,16)
            nptests = nptests + 1
          end if
          call conf_test(masv188,ldsref3_c(1),ldsref3_c(3),ldsref3_c(4),
     +                  ldsref3_c(2),1,c6)
          cmin4 = min(cmin4,c6)
          ngtests(4) = ngtests(4) + 1
        end if
 
c ...   debug statement ............................................
        if (debug .gt. 0) then
           write(h_output,'(1x,''masv188: '',6f10.4)')masv188,
     +           ldsref3_c(1),ldsref3_c(2),ldsref3_c(3),ldsref3_c(4)
        endif
c ................................................................
      endif
c ************   END OF GROUP 4 TESTS   ****************************
 
 

c     Check to see if thin cirrus bit should be set
      if ((.not. hi_elev) .and. visusd) then
        if (nint(masv188) .ne. nint(bad_data)) then
          call set_qa_bit(qa_bits,9)
          if (masv188.lt.ldstci_c(1) .and. masv188.ge.ldstci_c(2)) then
            call clear_bit(testbits,9)
            cirrus_vis = .true.
          endif
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
        write(h_output,'(1x,''confdnc '',8f8.5/,4f8.5)') c1,c2,c3,c4,c5,
     +         c6,cmin1,cmin2,cmin3,cmin4,fac,confdnc
      endif
c ................................................................

      return
      end
