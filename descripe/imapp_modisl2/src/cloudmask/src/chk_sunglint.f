      subroutine chk_sunglint(indat,pxldat,kele,confdnc,maxele,klin,
     *                        line_edge,qa_bits,testbits)


c---------------------------------------------------------------------
C!F77 
c
c!Description:
c
c     Performs clear sky restoral tests in sun-glint conditions.
c
c!Input parameters:
c indat         Array containing nlcntx lines of data
c pxldat        Array containing reflectance or brightness temperatures
c               for all bands for current pixel
c kele          Current granule element number being processed
c confdnc       Current pixel unobstructed confidence
c
c!Output Parameters:
c qa_bits       Byte array containing qa bits
c testbits      Byte array containing test results
c
c!Revision History:
c $Id: chk_sunglint.f,v 1.1.2.6 2004/10/25 17:13:35 raf Exp $
c 06/04 Collection 5  R. Frey
c Updated argument list in call to spatial_var
c 10/04 Collection 5  R. Frey
c Added band 2 "sigma*mean" clear-sky restoral test
c
c!Team-unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD35.
c
c!END
c---------------------------------------------------------------------

      save
   
      include 'global.inc'
      include 'snglntr_thr.inc'

c     Scalar arguments
      integer kele,klin,maxele
      real confdnc
      logical line_edge

c     Array arguments
      real indat(npixel,nlcntx,inband),pxldat(inband)
      byte testbits(6),qa_bits(10)

c     Local scalars
      integer h_output,debug,varslt,rtn,ipt
      real modir37,modir11,d37_11,modv895,modv935,modv443,rat,sigma,
     +     mean
      logical irclr

      integer band

c     Local Arrays
      real diff(var_band,8),rgdata(nlcntx,necntx,var_band)

c     External subroutines
      external spatial_var,get_regdif,check_bits,set_qa_bit,set_bit

c ... Common statement for debug purposes
      common / bug / debug, h_output

c---------------------------------------------------------------------

c     Get brightness temperature differences between pixel of interest
c     and the 8 surrounding it.
      call get_regdif(indat,kele,rgdata,diff)

c     Check variation in the region.
      call spatial_var(diff,ipt,varslt)
        
      if(varslt .eq. 1) then

c       Region is uniform in the 11 micron IR window.
c       Determine the logical flag 'irclr' - true if ir cloud tests below
c       have all been passed. APOLLO test makes final decision for bit 18.
        irclr = .true.
        call check_bits(testbits,13,rtn)
        if(rtn .eq. 0) irclr = .false.
        call check_bits(testbits,14,rtn)
        if(rtn .eq. 0) irclr = .false.
        call check_bits(testbits,15,rtn)
        if(rtn .eq. 0) irclr = .false.
        call check_bits(testbits,18,rtn)
        if(rtn .eq. 0) irclr = .false.
        call check_bits(testbits,27,rtn)
        if(rtn .eq. 0) irclr = .false.

        if(irclr) then

          modir37 = pxldat(20)
          modir11 = pxldat(31)
          modv895 = pxldat(17)
          modv935 = pxldat(18)
          modv443 = pxldat(9)

          if(modir37.ne.nint(bad_data) .and. modir11.ne.nint(bad_data) .and.
     +       modv895.ne.nint(bad_data) .and. modv935.ne.nint(bad_data)) then

c           Set bit which indicates this test was applied.
            call set_qa_bit(qa_bits,26)

            d37_11 = modir37 - modir11

c           Set bit if tests passed.
            if(d37_11 .ge. sg_tbdfl(1) ) then

              confdnc = 0.67

              rat = modv895 / modv935

              if(rat.gt.snglrat(1) .and. modv443.ne.nint(bad_data)) then

                call set_bit(testbits,26)
                confdnc = 0.96

              else

                band = 2
                call get_regstd(indat,kele,maxele,line_edge,klin,band,sigma,mean)

                if(mean .ne. bad_data) then
                  if( (sigma * mean) .lt. 0.001) then
                      call set_bit(testbits,26)
                      confdnc = 0.96
                  end if
                end if

              end if

            end if

          end if

c ......  debug statement ............................................
          if(debug .gt. 0) then
            write(h_output,'(10x,'' Sun-glint clear sky restoral: '')')
            write(h_output,'(10x,'' 3.7-11 um difference '',i5,f8.2)') varslt,d37_11
            write(h_output,'(10x,'' Channel 17-18 ratio  '',f8.2)') rat
            write(h_output,'(10x,'' Sigma test  '',3f10.5)') mean,sigma,sigma*mean
            write(h_output,'(10x,'' Confidence after sun-glint CSR tests
     +         is  '',f10.2,/)') confdnc
          end if
c ....................................................................

        end if

      end if

      return
      end
