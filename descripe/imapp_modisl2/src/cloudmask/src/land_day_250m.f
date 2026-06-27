      subroutine land_day_250m(npx,pxl250,desert_250,coast_250,ice_250,
     *                         snow_250,day_250,proc_250,visusd_250,
     *                         qa_250,cm1km,temp_qa,temp_bits)

      implicit none
      save

c---------------------------------------------------------------------
C!F77 
c
c!Description:
c     Routine for generating a cloud mask result for each of 16
c     pixels corresponding to the current 1-km FOV in daytime 
c     conditions over land surfaces.
c
c!Input parameters:
c npx           250-m pixel number (1-16)
c pxl250        Array containing reflectances for 250m fovs
c desert_250    Logical variable indicating desert processing
c snow_250      Logical variable indicating snow on surface
c ice_250       Logical variable indicating ice on surface
c coast_250     Logical variable indicating coastal region
c cm1km         Confidence of clear sky for closest 1-km FOV
c scan_day      Indicates day processing for pixels in current scan
c scan_process  Indicates 1-km processing began for pixels in current scan
c scan_visusd   Indicates visible data used for pixels in current scan
c qa_bitarray   QA 1-km info for pixels on current scan
c
c!Output Parameters:
c temp_bits      Byte array containing cloud mask results
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
      include 'LandDay_thr.inc'
      include 'ocean_day_thr.inc'

c     scalar arguments
      integer npx

c     array arguments
      real pxl250(vis_band,num250_per_1km)
      logical desert_250(num250_per_1km),snow_250(num250_per_1km),
     *        ice_250(num250_per_1km),coast_250(num250_per_1km),
     *        cm1km(num250_per_1km),day_250(num250_per_1km),
     *        proc_250(num250_per_1km),visusd_250(num250_per_1km),
     *        qa_250(num250_per_1km)
      byte temp_bits(6),temp_qa(10)

c     local scalars
      integer i,debug,h_output,vis_indx
      real mqkm_vis,vis1,vis2,s1,s2,vrat,etan,etad,eta

c     External subroutines
      external set_qa_bit,set_bit,conf_test

c     Common statement for debug purposes
      common / bug / debug, h_output

c---------------------------------------------------------------------
c     Debug statement 
      if (debug .gt. 0) then
        write(h_output,'(1x,''Subroutine land_day_250m.f: '',i4/)') 
     *        npx
      endif
c---------------------------------------------------------------------

      i = npx

c     Create 250m mask.

      if(day_250(i) .and. proc_250(i) .and. visusd_250(i) .and.
     *   qa_250(i)) then

        if(desert_250(i)) then
          vis_indx = 2
        else
          vis_indx = 1
        end if

        if (cm1km(i)) then

c         Assume 250 m pixel is clear - fill bit with 1 km result.
          mqkm_vis = pxl250(vis_indx,i)
          if (nint(mqkm_vis) .ne. nint(bad_data)) then
            call set_qa_bit(temp_qa,31+i)
            call set_bit(temp_bits,31+i)
          end if

        else if (ice_250(i) .or. snow_250(i)) then

c         Assume 250 m pixel is cloudy - but do not perform individual
c         pixel tests, so qa bit is set and test bit stays unset (0).
c         Bit agrees with 1-km result.
          mqkm_vis = pxl250(vis_indx,i)
          if (nint(mqkm_vis) .ne. nint(bad_data)) then
            call set_qa_bit(temp_qa,31+i)
          end if

        else 

c         Perform tests on current 250-m pixel. Look for clear sky.

          vis1 = pxl250(1,i) 
          vis2 = pxl250(2,i)
          s1 = pxl250(1,i) * 100.0
          s2 = pxl250(2,i) * 100.0

          if ((nint(vis1) .ne. nint(bad_data)) .and. (nint(vis2) .ne.
     *         nint(bad_data))) then

c           GEMI visible ratio test.

            call set_qa_bit(temp_qa,31+i)

            etan = 2.0 * (s2-s1) + 1.5*s2 + 0.5*s1
            etad = s2 + s1 + 0.5
            eta = etan / etad
            vrat=eta * (1.0-0.25*eta) - ((s1-0.125) / (1.0-s1))

c---------------------------------------------------------------------

c           Debug statement 
            if (debug .gt. 3) then
              write(h_output,'(1x,''250 m GEMI: '',i10,5f10.3)') i,
     *            vis1,vis2,vrat,dlvrat(3),doref2(3)
            endif

c---------------------------------------------------------------------

            if(vrat .gt. dlvrat(3)) then
c             Clear land.
              call set_bit(temp_bits,31+i)
            else if((.not. coast_250(i))) then
c             Can't be sure of clear, so copy 1-km mask to current bit.
              if(cm1km(i)) then
                call set_bit(temp_bits,31+i)
              end if
            else
c             If coast, perform visible threshold test using water threshold.
              if(vis2 .lt. doref2(3)) then
                call set_bit(temp_bits,31+i)
              else
c               Can't be sure of clear, so copy 1-km mask to current bit.
                if(cm1km(i)) then
                  call set_bit(temp_bits,31+i)
                end if
              end if
            end if
 
          end if

        end if

      end if

c--------------------------------------------------------------------

c     Debug statement 
      if (debug .gt. 3) then
          write(h_output,'(''Data from land_day_250m'')')
          write(h_output,'(9l5)') day_250(i),proc_250(i),visusd_250(i),
     *      qa_250(i),cm1km(i),desert_250(i),ice_250(i),snow_250(i),
     *      coast_250(i)
      endif

c---------------------------------------------------------------------

      return
      end
