      subroutine chk_coast(pxldat,confdnc,qa_bits,testbits)


c---------------------------------------------------------------------
C!F77 
c
c!Description:
c
c     Routine which checks for extreme values of NDVI in coastal 
c     regions.  Set confidence to "confident clear" if very low or
c     very high values are found.
c
c!Input parameters:
c pxldat        Array containing reflectance or brightness temperatures
c               for all bands for current pixel
c confdnc       Current pixel unobstructed confidence
c
c!Output Parameters:
c qa_bits       Byte array containing qa bits
c testbits      Byte array containing test results
c
c!Revision History:
c
c!Team-unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-35.
c
c!END
c---------------------------------------------------------------------

      save 

      include 'global.inc'
      include 'swc_ndvi.inc'

c     Scalar arguments
      real confdnc

c     Array arguments
      real pxldat(inband)
      byte testbits(6),qa_bits(10)

c     Local scalars
      integer rtn,debug,h_output
      real ndvi
      logical irclr

c     External subroutines
      external check_bits,set_bit,set_qa_bit

c     Common statement for debug purposes
      common / bug / debug, h_output

c---------------------------------------------------------------------

c     Determine the logical flag 'irclr' - true if ir cloud tests below
c     have all been passed. APOLLO test makes final decision for bit 18.

      irclr = .true.
      call check_bits(testbits,14,rtn)
      if(rtn .eq. 0) irclr = .false.
      call check_bits(testbits,15,rtn)
      if(rtn .eq. 0) irclr = .false.
      call check_bits(testbits,18,rtn)
      if(rtn .eq. 0) irclr = .false.

      if(irclr) then
        if (nint(pxldat(2)) .ne. nint(bad_data) .and.
     +           nint(pxldat(1)) .ne. nint(bad_data)) then

c         Check for very low or very high ndvi values.
          call set_qa_bit(qa_bits,22)
          ndvi = (pxldat(2) - pxldat(1)) / (pxldat(2) + pxldat(1))
          if(ndvi .le. swc_ndvi(1) .or. ndvi .ge. swc_ndvi(2)) then
            confdnc = 1.0
            call set_bit(testbits,22)
          end if

        end if
      end if

c---------------------------------------------------------------------

c     Debug statement.
      if(debug .gt. 0) then
        write(h_output,'(10x,'' Coastal NDVI test results: '')')
        write(h_output,'(10x,'' NDVI thresholds: '',2f10.5)') swc_ndvi
        write(h_output,'(10x,''irclr,ndvi,confdnc: '',l5,2f10.5/)')
     *         irclr,ndvi,confdnc
      end if

c---------------------------------------------------------------------

      return
      end
