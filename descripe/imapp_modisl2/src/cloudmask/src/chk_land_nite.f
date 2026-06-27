      subroutine chk_land_nite(pxldat,tbadj,confdnc,qa_bits,
     *                         testbits)

      implicit none
      save

c-----------------------------------------------------------------------
C!F77
c
c!Description:
c     Perform final clear-sky confidence check on land pixels.
c     If confidence of clear sky is low but temp is warm and if the
c     IR clear sky tests all passed, then use the 11 um brightness
c     temperature to assign a final confidence.
c
c!Input Parameters:
c pxldat        Array containing reflectance or brightness temperatures
c               for all bands for a single pixel
c tbadj         11 um brightness temperature threshold adjustment for
c               deserts
c
c!Output Parameters:
c testbits      Byte array containing cloud mask results
c qa_bits       Byte array containing qa bit results
c confdnc       Current pixel unobstructed confidence
c
c!Revision History:
c
c!Team-Unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-35.
c
c!END
c-----------------------------------------------------------------------

      include 'global.inc'
      include 'land_restoral.inc'

c     scalar arguments
      real confdnc,tbadj

c     array arguments
      real pxldat(inband)
      byte testbits(6),qa_bits(10)

c     local scalars
      integer j,h_output,debug
      real m31
      logical bit_test

c     local arrays
      integer bitno(4),rtn(4),rtnqa(4)
      real hds11(3)

c     external subroutines
      external check_bits,check_qa_bits,set_qa_bit

c     Common statement for debug purposes
      common / bug / debug, h_output

c     Set bit numbers to test in "final test".
      data bitno /14,15,17,23/

c-----------------------------------------------------------------------

      m31 = pxldat(31)

      if (nint(m31) .ne. nint(bad_data)) then

c       Check IR clear sky tests.
        bit_test = .true.
        do j = 1,4
          call check_qa_bits(qa_bits,bitno(j),rtnqa(j))
          call check_bits(testbits,bitno(j),rtn(j))
          if(rtnqa(j) .eq. 1) then
            if(rtn(j) .eq. 0) then
              bit_test = .false.
            end if
          end if
        enddo

        if(bit_test) then

c         Get elevation-adjusted 11 micron brightness temperature threshold.

          call set_qa_bit(qa_bits,26)

          hds11(1) = lnbt11(1) - tbadj
          hds11(2) = lnbt11(2) - tbadj
          hds11(3) = lnbt11(3) - tbadj

c         Check for hot scene.
          if(m31 .gt. hds11(1)) then

c           Assign confidence level based on 11 micron Tbb.
            if(m31 .gt. hds11(3)) then
c             Assign pixel to confident clear, set bit #26.
              confdnc = 1.0
              call set_bit(testbits,26)
            else if(m31 .gt. hds11(2)) then
c             Assign pixel to probably clear, set bit #26.
              confdnc = 0.96
              call set_bit(testbits,26)
            else
c             Assign pixel to uncertain, do not set bit #26.
              confdnc = 0.95
            end if

          end if

c-----------------------------------------------------------------------

c         debug statement
          if (debug .gt. 0) then
            write(h_output,'(''Night land restoral test: '',/,6f10.2,l5,f10.5,
     *            /)') lnbt11(1),tbadj,hds11(1),hds11(2),hds11(3),
     *            m31,bit_test,confdnc
          endif

c-----------------------------------------------------------------------

        end if

      end if

c-----------------------------------------------------------------------

      return
      end
