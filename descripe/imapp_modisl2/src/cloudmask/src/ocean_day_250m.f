      subroutine ocean_day_250m(npx,pxl250,snglnt_250,ice_250,cm1km,
     *                          day_250,proc_250,visusd_250,qa_250,
     *                          temp_qa,temp_bits)

      implicit none
      save

c---------------------------------------------------------------------
C!F77 
c
c!Description:
c     Routine for generating a cloud mask result for each of 16
c     pixels corresponding to the current 1-km FOV in daytime 
c     conditions over water surfaces.
c
c!Input parameters:
c npx           250m pixel number (1-16)
c pxl250        Array containing reflectances for 250m visible
c               fov's within the 1km footprint
c snglnt_250    Logical variable flagging sunglint pixels
c ice_250       Logical variable indicating ice on water surface
c cm1km         Confidence of clear sky based on closest 1-km FOV
c scan_day      Indicates day processing for pixels in current scan
c scan_process  Indicates 1-km processing began for pixels in current scan
c scan_visusd   Indicates visible data used for pixels in current scan
c qa_bitarray   QA 1-km info for pixels on current scan
c
c!Output Parameters:
c temp_bits     Byte array containing cloud mask results
c temp_qa       Byte array containing qa bit results
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

      include 'global.inc'
      include 'ocean_day_thr.inc'

c     scalar arguments
      integer i,npx

c     array arguments
      real pxl250(vis_band,num250_per_1km)
      logical snglnt_250(num250_per_1km),ice_250(num250_per_1km),
     *        cm1km(num250_per_1km),day_250(num250_per_1km),
     *        proc_250(num250_per_1km),visusd_250(num250_per_1km),
     *        qa_250(num250_per_1km)
      byte temp_bits(6),temp_qa(10)

c     local scalars
      integer debug,h_output
      real mqkm_88

c     External subroutines
      external set_qa_bit,set_bit,check_bits

c     Common statement for debug purposes
      common / bug / debug, h_output

c---------------------------------------------------------------------
c     Debug statement
      if (debug .gt. 0) then
        write(h_output,'(1x,''Subroutine ocean_day_250m.f: '',i4/)')
     *        npx
      endif
c---------------------------------------------------------------------

      i = npx

c     Create 250m mask.
        
      if(day_250(i) .and. proc_250(i) .and. visusd_250(i) .and.
     *   qa_250(i)) then

        mqkm_88 = pxl250(2,i)
        if (nint(mqkm_88) .ne. nint(bad_data)) then

          call set_qa_bit(temp_qa,31+i)

          if((mqkm_88 .ge. doref2(2) .and. mqkm_88 .le. doref2(1)) .or.
     *       (snglnt_250(i) .and. mqkm_88 .ge. doref2(2)) .or.
     *        ice_250(i)) then

c           Don't perform 250 m tests - copy 1 km results.
            if (cm1km(i)) then
              call set_bit(temp_bits,31+i)
            end if

          else if(mqkm_88 .lt. doref2(2)) then

c           Current 250 m pixel is clear - set bit.
            call set_bit(temp_bits,31+i)

          end if

        end if

      end if

c---------------------------------------------------------------------

c     Debug statement
      if (debug .gt. 0) then
        write(h_output,'(1x,''Data from ocean_day_250m.f: '',i4/)')
     *        npx
        write(h_output,'(6l5)') day_250(i),proc_250(i),visusd_250(i),
     *        qa_250(i),snglnt_250(i),ice_250(i)
        write(h_output,'(3f6.3)') doref2(1),doref2(3),mqkm_88
      endif

c---------------------------------------------------------------------


      return
      end
