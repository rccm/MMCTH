      subroutine shadows(pxldat,shadow,visusd,qa_bits)

      implicit none
      save

C!F77 ************************************************************
c!Description:
c determines the presence of shadows and sets the appropriate bit flag.
c
c!Input Parameters:
c pxldat        Array containing reflectance or brightness temperatures
c               for all bands for a single pixel
c visusd        Logical variable indicating whether visible data
c               was used or not
c 
c!Output Parameters:
c shadow        logical variable indicating shadow is present if .true.
c qa_bits       Byte array holding current pixel qa results
c
c!Revision History:
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
      include 'shadows_thr.inc'
      include 'global.inc'

c ... scalar arguments ..
      logical shadow,visusd
c ...
c ... array arguments ..
      real pxldat(inband)
      byte qa_bits(10)
c ...
c ... local scalars ..
      real masv88,masv66,masv945,vrat,masv124
      integer debug,h_output

c ... Common statement for debug purposes
      common / bug / debug, h_output

c ...
      masv66 = pxldat(1)
      masv88 = pxldat(2)
      masv945 = pxldat(19)
      masv124 = pxldat(5)

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(10x/,''Within shadows testing routine '',/)')
      endif
c ................................................................

 
c     Reflectance at 0.945 um must be ge 12% and "visible ratio" gt the
c     ocean threshold or else we're seeing a shadow.  The test on the 
c     ratio is to assure that the scene is not clear-sky conditions over
c     a sub-grid scale water body. 
      if (visusd) then
        if (nint(masv66) .ne. nint(bad_data)
     +                   .and.
     +      nint(masv88) .ne. nint(bad_data)
     +                   .and.
     +      nint(masv945) .ne. nint(bad_data)
     +                   .and.
     +      nint(masv124) .ne. nint(bad_data)) then

          vrat = (masv88 - masv66) / (masv66 + masv88)

c         Set qa bits which show that we actually did test for shadows
          call set_qa_bit(qa_bits,10)

c ...     debug statement ............................................
          if (debug .gt. 0) then
             write(h_output,'(10x,'' Shadows test '',5f10.2,/)') 
     *           masv66,masv88,vrat,masv945,masv124
          endif
c ....................................................................

          if(masv945 .lt. shadnir(1) .and. vrat .gt. shavrat(1) .and.
     *       masv945 .gt. shadnir(2) .and. masv124 .lt. shad124(1)) then
            shadow = .true.
          else
            shadow = .false.
          end if

c ....    debug statement ............................................
          if (debug .gt. 0) then
            if (shadow) then
              write(h_output,'(10x,'' Shadow found ''/)') 
            else 
              write(h_output,'(10x,'' Shadow not found ''/)') 
            endif
          endif
c ....................................................................
        endif
      endif

      return
      end
