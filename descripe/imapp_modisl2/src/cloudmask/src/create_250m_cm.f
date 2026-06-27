      subroutine create_250m_cm(nc,pxl250,cm1km,ice_250,snow_250,
     *                          desert_250,coast_250,land_250,
     *                          bitarray,qa_bitarray,n250,ncr,ncd,
     *                          day_250,proc_250,visusd_250,qa_250,
     *                          snglnt_250)

      implicit none
      save

c-----------------------------------------------------------------------
c!F77
c
c!DESCRIPTION:  Create the 250-m cloud mask for current 4x4 region.
c
c!Input Parameters:
c nc            Current 1-km pixel number
c pxl250        Array containing 250m relfectances in 4x4 region
c cm1km         1-km cloud mask result associated with each 250-m pixel
c ice_250       Indicates ice processing for 250-m pixels
c snow_250      Indicates snow processing for 250-m pixels
c desert_250    Indicates desert processing for 250-m pixels
c coast_250     Indicates coast processing for 250-m pixels
c land_250      Indicates land processing for 250-m pixels
c scan_day      Indicates day processing for pixels in current scan
c scan_process  Indicates 1-km processing began for pixels in current scan
c scan_visusd   Indicates visible data used for pixels in current scan
c qa_bitarray   QA 1-km info for pixels on current scan
c snglnt_250    Indicates sunglint processing for corresponding 1-km pixels
c
c!Output Parameters:
c bitarray      Array containing line of 48 bit test results
c qa_bitarray   Array containing line of 10 byte qa results
c n250          Counter for number of pixels included in stats
c               (includes all possible pixels)
c ncr           Number of 250m pixels found to be clear
c ncd           Number of 250m pixels found to be cloudy
c
c!Revision History:
c
c!Team-Unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-35.
c
c!END
c---------------------------------------------------------------------

      INCLUDE 'global.inc'

c     Array arguments.
      real pxl250(vis_band,num250_per_1km)
      logical ice_250(num250_per_1km),snow_250(num250_per_1km),
     *        desert_250(num250_per_1km),coast_250(num250_per_1km),
     *        land_250(num250_per_1km),
     *        cm1km(num250_per_1km),day_250(num250_per_1km),
     *        proc_250(num250_per_1km),visusd_250(num250_per_1km),
     *        qa_250(num250_per_1km),snglnt_250(num250_per_1km)
      byte bitarray(npixel,6),qa_bitarray(npixel,10)

c     Scalar arguments.
      integer nc,n250,ncr,ncd

c     Local arrays.
      integer bits(48),qa(80)
      byte temp_bits(6),temp_qa(10)

c     Local scalars.
      integer i,npx,n250pix,nclr,ncld,ival,ipos,ibyte,j,debug,
     *        h_output,rtn,kk,jj
      logical lval

c     External routines.
      external land_day_250m,ocean_day_250m,check_bits,check_qa_bits

c ... Common statement for debug purposes
      common / bug / debug, h_output

c     Initialize counters.
      data n250pix /0/, nclr /0/ , ncld /0/

c---------------------------------------------------------------------

c     Get bit arrays for current 1-km pixel.

      do i = 1, 10
        if(i .le. 6) temp_bits(i) = bitarray(nc,i)
        temp_qa(i) = qa_bitarray(nc,i)
      enddo

c---------------------------------------------------------------------

c     Process 250-m pixels.

      do npx = 1, num250_per_1km

        if(land_250(npx)) then

          call land_day_250m(npx,pxl250,desert_250,coast_250,ice_250,
     *                       snow_250,day_250,proc_250,visusd_250,
     *                       qa_250,cm1km,temp_qa,temp_bits)

        else

          call ocean_day_250m(npx,pxl250,snglnt_250,ice_250,cm1km,
     *                        day_250,proc_250,visusd_250,qa_250,
     *                        temp_qa,temp_bits)

        end if

      enddo

c---------------------------------------------------------------------

c     Fill in 250-m information.

      do i = 1, 10
        if(i .le. 6) bitarray(nc,i) = temp_bits(i)
        qa_bitarray(nc,i) = temp_qa(i)
      enddo

c---------------------------------------------------------------------

c     Collect info for output statistics.

      do i = 1 , num250_per_1km
        n250pix = n250pix + 1
        call check_qa_bits(temp_qa,31+i,rtn)
        if (rtn .eq. 1) then
          call check_bits(temp_bits,31+i,rtn)
          if (rtn .eq. 1) then
            nclr = nclr + 1
          else
            ncld = ncld + 1
          endif
        endif
      enddo
      n250 = n250pix
      ncr = nclr
      ncd = ncld

c---------------------------------------------------------------------

c     Debug statement.
      if (debug .gt. 0) then
        write(h_output,'(10x,'' create_250m_cm:'',/,3i10,/)')
     *        n250,ncr,ncd
      end if

c---------------------------------------------------------------------

c     Debug section.
      if (debug .gt. 0) then
c        strip out and print bit values for current pixel
         ival = 0
         ipos = 0
         do i = 1 , 48
            bits(i) = 0
         enddo
         do ibyte = 1,6
            ival = bitarray(nc,ibyte)
            do j = 1,8
               ipos = ipos + 1
               lval = btest(ival,j-1)
               if(lval) then
                  bits(ipos) = 1
               else
                  bits(ipos) = 0
              end if
            enddo
         enddo

         write(h_output,'(1x,3x,'' CM Bits Results'')')
         write(h_output,'(1x,3x,''FOV No'',3x,16(i2,2x))') (kk,kk=0,15)
         write(h_output,'(1x)')
         write(h_output,'(1x,i8,2x,16i4)') nc,(bits(jj),jj=1,16)
         write(h_output,'(1x,i8,2x,16i4)') nc,(bits(jj),jj=17,32)
         write(h_output,'(1x,i8,2x,16i4,/)') nc,(bits(jj),jj=33,48)

         ipos = 0
         do i = 1 , 80
            qa(i) = 0
         enddo
         do ibyte = 1,10
            ival = qa_bitarray(nc,ibyte)
            do j = 1,8
               ipos = ipos + 1
               lval = btest(ival,j-1)
               if(lval) then
                  qa(ipos) = 1
               else
                  qa(ipos) = 0
              end if
            enddo
         enddo

         write(h_output,'(1x,3x,'' QA Bits Results'')')
         write(h_output,'(1x,3x,''FOV No'',3x,16(i2,2x))') (kk,kk=0,15)
         write(h_output,'(1x)')
         write(h_output,'(1x,i8,2x,16i4)') nc,(qa(jj),jj=1,16)
         write(h_output,'(1x,i8,2x,16i4)') nc,(qa(jj),jj=17,32)
         write(h_output,'(1x,i8,2x,16i4)') nc,(qa(jj),jj=33,48)
         write(h_output,'(1x,i8,2x,16i4)') nc,(qa(jj),jj=49,64)
         write(h_output,'(1x,i8,2x,16i4,/)') nc,(qa(jj),jj=65,80)

        write(h_output,'(10x,'' testbits pixel results '',/,6i5,/,6i5)')
     *        temp_bits(1),temp_bits(2),temp_bits(3),temp_bits(4),
     *        temp_bits(5),temp_bits(6),bitarray(nc,1),bitarray(nc,2),
     *        bitarray(nc,3),bitarray(nc,4),bitarray(nc,5),
     *        bitarray(nc,6)
        write(h_output,'(10x,'' qa_bits pixel results ''/10i5,/10i5)')
     *        temp_qa(1),temp_qa(2),temp_qa(3),temp_qa(4),
     *        temp_qa(5),temp_qa(6),temp_qa(7),temp_qa(8),
     *        temp_qa(9),temp_qa(10),qa_bitarray(nc,1),
     *        qa_bitarray(nc,2),qa_bitarray(nc,3),qa_bitarray(nc,4),
     *        qa_bitarray(nc,5),qa_bitarray(nc,6),qa_bitarray(nc,7),
     *        qa_bitarray(nc,8),qa_bitarray(nc,9),qa_bitarray(nc,10)
      endif

c---------------------------------------------------------------------

      return
      end
