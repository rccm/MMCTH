      subroutine Antarctic_day(pxldat,visusd,testbits,qa_bits,
     +                         nmtests,confdnc)

      implicit none
      save
 
c---------------------------------------------------------------------
c!F77 
c
c!Description:
c      Routine for performing clear sky tests over Antarctic snow 
c      surfaces during daylight hours.
c
c      For daytime land type 1 the groups are:
c          Group 1: High thick cloud
c                   6.7 micron bt test
c
c          Group 2: Low cloud - thick
c                   11-4 micron bt tests
c        
c          Group 4: Thin cirrus test
c                   1.38 micron reflectance test 
c
c!Input Parameters:
c pxldat        Array containing reflectance or brightness temperatures
c               for all bands for a single pixel
c visusd        Logical variable indicating whether vis data used or not
c
c
c!Output Parameters:
c testbits      six byte integer containing bit results
c qa_bits       ten byte integer containing qa bit results
c nmtests       number spectral tests performed
c confdnc       product of all applied individual confidences
c
c!Revision History:
c
c Added 11-12 um thin cirrus test
c 06/04 Collection 5    R. Frey
c Removed 11-12 um thin cirrus test.
c Added 11 um BT-dependent 3.0-11 um BTD test.  Replaces test with static
c thresholds.
c 10/04 Collection 5b   R. Frey
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
      include 'Antarctic_day_thr.inc'

c ...
c ... scalar arguments ..
      real confdnc
      integer nmtests
      logical visusd
c ...
c ... array arguments ..
      real pxldat(inband)
      byte testbits(6),qa_bits(10)
c ...
c ... local scalars ..
      real c2,c3,mas11_4,cmin1,cmin2,
     +     masir11,masir4,locut,hicut,
     +     masir65,groups,fac,pre_confdnc,midpt,power
      integer nptests,kk,debug,h_output

c ... local arrays
      integer ngtests(3)
     
      real Rel_equality_EPS
      parameter(Rel_equality_EPS = 0.000001)

c ... external subroutines ..
      external conf_test,set_bit,clear_bit,set_qa_bit

c     Common statement for debug purposes
      common / bug / debug, h_output

c ...
c ... nmtests counts the number of tests applied to this pixel
      nmtests = 0

c ... initialize variables

c ... nptests counts the number of tests passed
      nptests = 0
c ... set confidence to 1.0 to begin with
      confdnc = 1.0
c ... place band values into individual variables for easy
c ... identification
      masir4 = pxldat(22)
      masir65 = pxldat(27)
      masir11 = pxldat(31)
c ...
      mas11_4 = 0.0

c ... the c suffix variables represent individual test confidences
      c2 = 0.0
      c3 = 0.0
      cmin1 = 1.0
      cmin2 = 1.0

c ... initialize group number holder
      do 10 kk = 1 , 2
         ngtests(kk) = 0
  10  continue

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(10x/,''Processing subroutine Antarctic_day '',
     +                   /)')
      endif
c ................................................................
 
c ... perform tests.  
 
c     H20 vapor channel (6.7 micron) high cloud test
      if (nint(masir65) .ne. nint(bad_data)) then
        nmtests = nmtests + 1
        call set_qa_bit(qa_bits,15)
        if (masir65 .gt. anth20(2)) then
          call set_bit(testbits,15)
          nptests = nptests + 1
        end if
        call conf_test(masir65,anth20(1),anth20(3),anth20(4),
     *               anth20(2),1,c2)
        cmin1 = min(cmin1,c2)
        ngtests(1) = ngtests(1) + 1
      endif

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(1x,''masir65: '',5f10.2)') masir65,anth20(1),
     +          anth20(2),anth20(3),anth20(4)
      endif
c ................................................................
c     *****  END OF GROUP 1 TESTS  ***************************
 
 
 
c     ****  GROUP 2 TESTS  ***********************************

c ... 11 minus 3.9 micron BTDIF fog and low cloud test.
      if (visusd) then
        if (nint(masir11) .ne. nint(bad_data)
     +                   .and.
     +      nint(masir4) .ne.  nint(bad_data)) then

          if(masir11 .gt. 230.0) then

            nmtests = nmtests + 1
            call set_qa_bit(qa_bits,19)

            mas11_4 = masir4 - masir11

            call get_pn_thresholds(masir11,bt_11_bnds4,ant4_11l,ant4_11m1,
     *                             ant4_11m2,ant4_11m3,ant4_11h,locut,hicut,
     *                             midpt,power)

            if (mas11_4 .le. midpt) then
              call set_bit(testbits,19)
              nptests = nptests + 1
            end if
            call conf_test(mas11_4,locut,hicut,power,midpt,1,c3)
            cmin2 = min(cmin2,c3)
            ngtests(2) = ngtests(2) + 1

          end if
        end if

c ..... debug statement ............................................
        if (debug .gt. 0) then
         write(h_output,'(1x,''mas11_4: '',6f9.3)')masir11,mas11_4,
     +         locut,midpt,hicut,c3
        endif
c ..................................................................
      end if

c *******     END OF GROUP 2 TESTS ****************************
 
c     Determine final confidence based on group values
      pre_confdnc = cmin1 * cmin2

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
        write(h_output,'(1x,''tests '',3i10)') nmtests,nptests,ngtests
        write(h_output,'(1x,''confdnc '',6f8.5)') c2,c3,
     +         cmin1,cmin2,fac,confdnc
      endif
c ................................................................

      return
      end
