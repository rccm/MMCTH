      subroutine Day_snow(pxldat,vza,visusd,cirrus_vis,hi_elev,
     +                    testbits,qa_bits,nmtests,confdnc)

      implicit none
      save
 
c---------------------------------------------------------------------
c!F77 
c
c!Description:
c      Routine for performing clear sky tests over snow 
c      surfaces during daylight hours.
c
c      For daytime snow the groups are:
c          Group 1: High thick cloud
c                   13.9 micron bt test
c                   6.75 micron bt test 
c
c          Group 2: Low cloud - thick
c                   11-4 micron bt test
c                   11-12 micron bt test
c        
c          Group 4: Thin cirrus test
c                   1.38 micron reflectance test 
c
c!Input Parameters:
c pxldat        Array containing reflectance or brightness temperatures
c               for all bands for a single pixel
c vza           Viewing zenith angle for current pixel
c visusd        Logical variable indicating whether vis data used or not
c cirrus_vis 	Logical variable flagging thin cirrus contaminated
c               scenes in the visible
c hi_elev       Logical variable indicating high elevation (> 2000 meters)

c
c!Output Parameters:
c testbits      six byte integer containing bit results
c qa_bits       ten byte integer containing qa bit results
c confdnc       product of all applied individual confidences
c
c!Revision History:
c
c Added 11-12 thin cirrus test 
c 06/04 Collection 5   R. Frey
c
c!Team-unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-35.
c
c!Design Notes:
c    Externals:
c       Subroutines conf_test,set_bit,clear_bit,set_qa_bit
c
c!END
c-----------------------------------------------------------------------

      include 'global.inc'
      include 'Day_snow_thr.inc'

c ...
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
      real c1,c2,c3,c4,c6,mas11_4,cmin1,cmin2,cmin4,locut,hicut,
     +     masir11,masir13,masir4,masv188,masir65,masir12,
     +     groups,fac,pre_confdnc,schi,cosvza,masdf1,pi,dtr,diftemp,
     +     dfthrsh
      integer nptests,kk,debug,h_output

c ... local arrays
      integer ngtests(3)
      real sn4_11(4)
      
      real Rel_equality_EPS
      parameter(Rel_equality_EPS = 0.000001)

c ... external subroutines ..
      external conf_test,set_bit,clear_bit,set_qa_bit

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
      masv188 = pxldat(26)
      masir4 = pxldat(22)
      masir65 = pxldat(27)
      masir11 = pxldat(31)
      masir12 = pxldat(32)
      masir13 = pxldat(35)

c ...
      mas11_4 = 0.0

c ... the c suffix variables represent individual test confidences
      c1 = 0.0
      c2 = 0.0
      c3 = 0.0
      c4 = 0.0
      cmin1 = 1.0
      cmin2 = 1.0
      cmin4 = 1.0

c ... initialize group number holder
      do 10 kk = 1 , 3
         ngtests(kk) = 0
  10  continue

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(10x/,''Processing subroutine Day_snow '',
     +                   /)')
      endif
c ................................................................
 
c ... perform tests.  note that some tests are not used
c ... in sunglint conditions.
c
c ... co2 high cloud test
      if (nint(masir13) .ne. nint(bad_data)) then
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,14)
        if (masir13 .gt. dsco2(2)) then
          call set_bit(testbits,14)
          nptests = nptests + 1
        end if
        call conf_test(masir13,dsco2(1),dsco2(3),dsco2(4),
     *                dsco2(2),1,c1)
        cmin1 = min(cmin1,c1)
        ngtests(1) = ngtests(1) + 1
      endif

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(1x,''masir13: '',5f10.2)') masir13,dsco2(1),
     +          dsco2(2),dsco2(3),dsco2(4)
      endif
c ................................................................


c     H20 vapor channel (6.7 micron) high cloud test
      if (nint(masir65) .ne. nint(bad_data)) then
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,15)
        if (masir65 .gt. dsh20(2)) then
          call set_bit(testbits,15)
          nptests = nptests + 1
        end if
        call conf_test(masir65,dsh20(1),dsh20(3),dsh20(4),
     *               dsh20(2),1,c2)
        cmin1 = min(cmin1,c2)
        ngtests(1) = ngtests(1) + 1
      endif

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(1x,''masir65: '',5f10.2)') masir65,dsh20(1),
     +          dsh20(2),dsh20(3),dsh20(4)
      endif
c ................................................................
c     *****  END OF GROUP 1 TESTS  ***************************
 
 
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
          dfthrsh = ds11_12hi(1)
        else
c         Add adjustment for snow cover.
          dfthrsh = diftemp + ds11_12adj(1)
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
c       hicut = dfthrsh - 1.25
        call conf_test(masdf1,locut,hicut,1.0,dfthrsh,1,c6)
        cmin2 = min(cmin2,c6)
        ngtests(2) = ngtests(2) + 1
      endif
 
c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(1x,''APOLLO masdf1: '',8f10.2)') masdf1,
     +          ds11_12hi(1),ds11_12adj(1),
     +          masir11,schi,dfthrsh,locut,hicut
      endif
c ................................................................

c ... 11 minus 4 micron BTDIF fog and low cloud test.
      if (visusd) then
        if (nint(masir11) .ne. nint(bad_data)
     +                   .and.
     +      nint(masir4) .ne.  nint(bad_data)) then
          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,19)
          mas11_4 = masir4 - masir11

          if(hi_elev) then
            sn4_11(1) = ds4_11hel(1)
            sn4_11(2) = ds4_11hel(2)
            sn4_11(3) = ds4_11hel(3)
            sn4_11(4) = ds4_11hel(4)
          else
            sn4_11(1) = ds4_11(1)
            sn4_11(2) = ds4_11(2)
            sn4_11(3) = ds4_11(3)
            sn4_11(4) = ds4_11(4)
          end if

          if (mas11_4 .le. sn4_11(2)) then
            call set_bit(testbits,19)
            nptests = nptests + 1
          end if
          call conf_test(mas11_4,sn4_11(1),sn4_11(3),sn4_11(4),
     +                  sn4_11(2),1,c3)
          cmin2 = min(cmin2,c3)
          ngtests(2) = ngtests(2) + 1
        end if

c ..... debug statement ............................................
        if (debug .gt. 0) then
         write(h_output,'(1x,''mas11_4: '',5f10.2)')mas11_4,sn4_11(1),
     +            sn4_11(2),sn4_11(3),sn4_11(4)
        endif
c ..................................................................
      end if

c *******     END OF GROUP 2 TESTS ****************************
c
c
c ***********   START OF GROUP 4 TESTS  *************************
c ... near infrared high cloud test
      if ((.not. hi_elev) .and. visusd) then
        if (nint(masv188) .ne. nint(bad_data)) then
          nmtests = nmtests + 1
          call set_qa_bit(qa_bits,16)
          if (masv188 .le. dsref3(2)) then
            call set_bit(testbits,16)
            nptests = nptests + 1
          end if
          call conf_test(masv188,dsref3(1),dsref3(3),dsref3(4),
     +                  dsref3(2),1,c4)
          cmin4 = min(cmin4,c4)
          ngtests(3) = ngtests(3) + 1
        end if

c ...   debug statement ............................................
        if (debug .gt. 0) then
          write(h_output,'(1x,''masv188: '',6f10.4)')masv188,dsref3(1),
     +                dsref3(2),dsref3(3),dsref3(4)
        endif
c ................................................................
      endif

c ************   END OF GROUP 4 TESTS   ****************************
c
c     Check to see if thin cirrus bit should be set
      if ((.not. hi_elev) .and. visusd) then
        if (nint(masv188) .ne. nint(bad_data)) then
          call set_qa_bit(qa_bits,9)
          if(masv188 .lt. dstci(1) .and. masv188 .ge. dstci(2)) then
            call clear_bit(testbits,9)
            cirrus_vis = .true.
          endif
c ...     debug statement ............................................
          if (debug .gt. 0) then
            write(h_output,'(1x,''NIR Thin cirrus: '',3f10.4)')masv188,
     +                    dstci(1),dstci(2)
          endif
c ................................................................
        endif
      endif
c
c     Determine final confidence based on group values
      pre_confdnc = cmin1 * cmin2 * cmin4

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

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(1x,''tests '',6i10)') nmtests,nptests,ngtests
        write(h_output,'(1x,''confdnc '',9f8.5/,2f8.5)') c1,c2,c3,c4,c6,
     +         cmin1,cmin2,cmin4,fac,confdnc
      endif
c ................................................................

      return
      end
