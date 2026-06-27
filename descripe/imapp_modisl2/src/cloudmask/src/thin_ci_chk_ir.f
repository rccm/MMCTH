      subroutine thin_ci_chk_ir(pxldat,vza,cirrus_ir,qa_bits,testbits)


      implicit none
      save

c---------------------------------------------------------------------
C!F77
c
c!Description:
c ... Routine to test for thin cirrus using IR channels.  This 
c ... will indicate whether the we believe the cirrus is thin
c ... enough for most tests to be applied without affecting
c ... results.  It will allow PI's with algorithms which are
c ... very sensitive to thin cirrus contamination to see if
c ... it might be there without affecting the final cloud mask
c ... confidence.
c
c!Input parameters:
c pxldat        Array containing reflectance or brightness temperatures
c               for all bands for a single pixel
c vza           Current pixel viewing angle
c
c!Output Parameters:
c cirrus_ir     Logical variable flagging thin cirrus contaminated
c               scenes in the infrared
c qa_bits       10 byte array contining qa bit results
c testbits      6 byte array containing bit results
c
c!Revision History:
c 10/04  Collection 5b   R. Frey
c Corrected thresholds.
c
c!Team-Unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c
c!END
c---------------------------------------------------------------------

      include 'global.inc'

c ... scalar arguments ..
      real vza
      logical cirrus_ir

c ... array arguments ..
      real pxldat(inband)
      byte testbits(6),qa_bits(10)

c ... local scalars
      real cosvza,pi,dtr,masdf1,masir11,masir12,schi,dfthrsh,ci1,ci2,
     +     diftemp,Rel_equality_EPS
      logical code
      integer debug,h_output
      parameter (Rel_equality_EPS = 0.000001)

c ... external subroutines ..
      external tview,clear_bit,set_qa_bit

c ... intrinsic functions
      intrinsic cos

c ... Common statement for debug purposes
      common / bug / debug, h_output

c     Routine which checks for the presence of thin cirrus. This check 
c     is made independently of other spectral tests which may check
c     for similar conditions.
c
c     Check to see if IR thin cirrus bit should be set.
c     Right now this is based upon the APOLLO thin cirrus 11-12 BTDIF.
c     This test has been fairly robust over all but snow covered 
c     regions.

c ... assignment statements
      masir11 = pxldat(31)
      masir12 = pxldat(32)

c ... initialize variables
      pi = acos(-1.0)
      dtr = pi/180.0
      masdf1 = 0.0
      cosvza = 0.0
      diftemp = 0.0
      schi = 0.0 
      dfthrsh = 0.0
      code = .true.
      ci1 = 0.0
      ci2 = 0.0

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(10x/,''Processing thin_ci_chk_ir routine'',/)')
      endif
c ................................................................


c ... 11-12um brightness temperature difference test
c ... for low clouds).
      if (nint(masir11) .ne. nint(bad_data)
     +                 .and.
     +    nint(masir12) .ne. nint(bad_data)
     +                  .and.
     +              vza .gt. 0.0) then

        masdf1 = masir11 - masir12

c ...   11-12um brightness temperature difference test
c ...   for thin cirrus).
c ...   added apollo viewing angle/masir11 regressed threshold.
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
c ...   else don't use this threshold
        if (diftemp.lt.0.1 .or. abs(schi-99.0).lt.0.0001) then
          code = .false.
        else
          dfthrsh = diftemp
        endif
      
c ...   Want to use a threshold range of very thin cirrus.
        if (code) then
          call set_qa_bit(qa_bits,11)
          ci1 = dfthrsh
          ci2 = dfthrsh + (0.3 * dfthrsh) 
          if (masdf1 .gt. ci1 .and. masdf1 .le. ci2) then
             call clear_bit(testbits,11)
             cirrus_ir = .true.
          endif
        endif

      endif

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(10x,'' thin_cirrus_ir vars:'',/,6f10.2,/)')
     *        masir11,masir12,masdf1,vza,schi,diftemp
        write(h_output,'(10x,'' more variables:'',/,l4,3f10.2,l4,/)')
     *        code,dfthrsh,ci1,ci2,cirrus_ir
      endif
c .................................................................

      return
      end

