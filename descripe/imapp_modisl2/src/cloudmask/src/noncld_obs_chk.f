      subroutine noncld_obs_chk(indat,pxldat,confdnc,kele,maxele,
     +                          line_edge,klin,qa_bits,testbits,
     +                          smoke)

      implicit none
      save

C!F77 ************************************************************
c!Description:
c ... Routine which checks for the presence of a non-cloud obstruction.
c ... Currently checks for smoke.
c
c!Input Parameters:
c indat         Array containing nlcltx number of lines of
c               reflectance or BT values for each channel
c pxldat        Array containing reflectance or brightness temperatures
c               for all bands for a single pixel
c confdnc       clear-sky confidence for current pixel
c kele          First element pixel in the region (context)
c maxele        Maximum number of pixels in a scan
c line_edge     Logical variable indicating first or last line of
c               data in a granule
c klin          Counter indicating number of lines processed
c
c!Output Parameters:
c qa_bits       Byte array containing qa results
c testbits      Byte array containing test results
c smoke         logical variable indicating whether smoke is present
c
c!Revision History:
c 10/04  Collection 5b    R. Frey
c Changed call to get_regstd.
c
c!Team-unique Header:
c
c    This software is developed by the MODIS Science Data Support Team
c    for the National Aeronautics and Space Administration,
c    Goddard Space Flight Center, under contract NAS5-32373.
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c
c!Design Notes:
c
c!END****************************************************************
c
      include 'global.inc'
      include 'noncld_obs_chk.inc'

c ... scalar arguments ..
      logical smoke,line_edge
      real confdnc
      integer kele,maxele,klin

c ... array arguments
      real pxldat(inband),indat(npixel,nlcntx,inband)
      byte qa_bits(10),testbits(6)

c ... local scalars ..
      logical thk_smoke,fire,bit_test
      real masv21,masv66,masir3,masir11,tdif,coef,masir12,tdiff,
     +     masv55,masv47,masv86,smkrat,vrat,sigma,masir3f,mean
      integer band,debug,h_output,j

c ... local arrays
      integer bitno(4),rtn(4),rtnqa(4)

C ... Common statement for debug purposes
      common / bug / debug, h_output

c     Routine which checks for the possible presence of smoke.

c     Set bit numbers to test.
      data bitno /15,16,18,19/

c     Initializations.
      fire = .false.
      thk_smoke = .false.

      masv47 = pxldat(3)
      masv86 = pxldat(2)
      masv66 = pxldat(1)
      masv55 = pxldat(4)
      masv21 = pxldat(7)
      masir3 = pxldat(20)
      masir3f = pxldat(21)
      masir11 = pxldat(31)
      masir12 = pxldat(32)

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(10x/,''Within nco bit testing routine '',/)')
      endif
c ................................................................

c     Check clear sky tests.
      bit_test = .true.
      do j = 1,4
        call check_qa_bits(qa_bits,bitno(j),rtnqa(j))
        call check_bits(testbits,bitno(j),rtn(j))
        if(rtnqa(j) .eq. 1) then
          if(rtn(j) .eq. 0) then
            bit_test = .false.
          end if
        else
          bit_test = .false.
        end if
      enddo

c     Do not perform fire or smoke tests if tests reported in bits
c     15, 16, 18, or 19 reported cloud.

      if(bit_test) then

c       Set bit to say that we did attempt nco smoke test
        call set_qa_bit(qa_bits,8)

c       Check for fires (hot spots).

        if(masir3f .ne. nint(bad_data) .and. masir11 .ne. nint(bad_data)) then
          tdif = masir3f - masir11
          if (masir3f .gt. nc_bt37(1) .and. tdif .gt. nc37_11(1)) then
            fire = .true.
          end if 
        end if 

c       Test for thick smoke. 

        if (nint(masv66) .ne. nint(bad_data) .and.
     +    nint(masv21) .ne. nint(bad_data) .and.
     +    nint(masv47) .ne. nint(bad_data) .and.
     +    nint(masv47) .ne. nint(bad_data) .and.
     +    nint(masir3) .ne. nint(bad_data) .and.
     +    nint(masir11) .ne. nint(bad_data) ) then

          if((masv21*100.0) .lt. nc21(1)) then

            coef = 6.0 + masv21*100.0
            if((masv66*100.0) .gt. coef) then

              smkrat = masv47 / masv66
              if(smkrat .ge. ncrat(1)) then


                vrat = masv86 / masv66
                if(vrat .ge. ncvrat(1)) then

                  band = 1
                  call get_regstd(indat,kele,maxele,line_edge,klin,
     *                            band,sigma,mean)
                  if(sigma .le. ncsig(1)) then

                    thk_smoke = .true.

                  end if
                end if
              end if
            end if
          end if
        end if

      end if

      if(thk_smoke .or. fire) smoke = .true.

c ...   debug statement ............................................
        if (debug .gt. 0) then
          write(h_output,'(10x,'' NCO test '',/,11f8.2,3l4/)')
     *    masir3,masir11,tdif,masv21,coef,masv66,masv47,smkrat,masv86,
     *    vrat,sigma,bit_test,fire,thk_smoke
        endif
c ..................................................................

c     Perform dust test.
      if (nint(masir11) .ne. nint(bad_data) .and.
     +    nint(masir12) .ne. nint(bad_data)) then

        if(confdnc .gt. 0.67) then

c         Set bit to say that we attempted the nco dust test
          call set_qa_bit(qa_bits,28)

          tdiff = masir11 - masir12

          if(tdiff .lt. nc11_12(1)) then
            call clear_bit(testbits,28)
          end if

        end if

      end if

c ..................................................................

c ....debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(10x,'' Dust test '',/,4f8.2/)')
     *    masir11,masir12,tdiff,nc11_12(1)
      endif
c ................................................................

      return
      end
