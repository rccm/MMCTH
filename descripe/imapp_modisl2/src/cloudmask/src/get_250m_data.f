      subroutine get_250m_data(nlin,lines_in_edge,klin,nc,maxele,
     *                         v250_dat,scan_confdnc,scan_ice,
     *                         scan_snow,scan_desert,scan_coast,
     *                         scan_land,pxl250,cm1km,
     *                         ice_250,snow_250,desert_250,coast_250,
     *                         land_250,scan_snglnt,snglnt_250,
     *                         scan_day,scan_process,scan_visusd,
     *                         qa_bitarray,day_250,proc_250,visusd_250,
     *                         qa_250)

      implicit none
      save

c-----------------------------------------------------------------------
c!F77
c
c!DESCRIPTION:  Routine to acquire data for 250-m pixel processing in
c               current 4x4 region.
c
c!Input Parameters:
c nlin          Total number of lines to process in granule
c lines_in_edge Number of lines outside of processing region
c klin          Current line being processed (1-km)
c nc            Current 1-km pixel number
c maxele        Number 1-km pixels in scan line
c v250_dat      Array containing 250m visible reflectance values
c scan_confdnc  Clear-sky confidence of all 1-km pixels in current scan
c scan_ice      Indicates ice processing for pixels in current scan
c scan_snow     Indicates snow processing for pixels in current scan
c scan_desert   Indicates desert processing for pixels in current scan
c scan_coast    Indicates coast processing for pixels in current scan
c scan_land     Indicates land processing for pixels in current scan
c scan_snglnt   Indicates sunglint processing for pixels in current scan
c scan_day      Indicates day processing for pixels in current scan
c scan_process  Indicates 1-km processing began for pixels in current scan
c scan_visusd   Indicates visible data used for pixels in current scan
c qa_bitarray   QA 1-km info for pixels on current scan
c
c!Output Parameters:
c pxl250        Array containing 250m relfectances in 4x4 region 
c cm1km         1-km cloud mask result associated with each 250-m pixel
c ice_250       Indicates ice processing for corresponding 1-km pixels
c snow_250      Indicates snow processing for corresponding 1-km pixels
c desert_250    Indicates desert processing for corresponding 1-km pixels
c coast_250     Indicates coast processing for corresponding 1-km pixels
c land_250      Indicates land processing for corresponding 1-km pixels
c snglnt_250    Indicates sunglint processing for corresponding 1-km pixels
c day_250       Indicates day processing for corresponding 1-km pixels
c proc_250      Indicates 1-km cm processed for corresponding 1-km pixels
c visusd_250    Indicates visible data used for corresponding 1-km pixels
c qa_250        QA info for corresponding 1-km pixels
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
      real v250_dat(nx*4,nlcntx,vis_band),scan_confdnc(npixel),
     *     pxl250(vis_band,num250_per_1km)
      logical scan_ice(npixel),scan_snow(npixel),scan_desert(npixel),
     *        scan_coast(npixel),scan_land(npixel),
     *        ice_250(num250_per_1km),snow_250(num250_per_1km),
     *        desert_250(num250_per_1km),coast_250(num250_per_1km),
     *        land_250(num250_per_1km),
     *        cm1km(num250_per_1km),snglnt_250(num250_per_1km),
     *        scan_snglnt(npixel),scan_day(npixel),scan_process(npixel),
     *        scan_visusd(npixel),day_250(num250_per_1km),
     *        proc_250(num250_per_1km),visusd_250(num250_per_1km),
     *        qa_250(num250_per_1km)
      byte qa_bitarray(npixel,10)

c     Scalar arguments.
      integer nc,nlin,klin,lines_in_edge,maxele

c     Local scalars.
      integer i,jj,k,kk,j250,debug,h_output

c ... Common statement for debug purposes
      common / bug / debug, h_output

c---------------------------------------------------------------------

c     Put 250-m pixel data into data array.  If working on edge pixels
c     then take either the first or last of 3 in context...otherwise
c     always take the middle line

      if (klin .le. lines_in_edge) then
         i = klin
      else if (klin .gt. nlin-lines_in_edge) then
         i = nlcntx
      else
         i = ((nlcntx-1) / 2) + 1
      endif

      do 500 k = 1 , vis_band
         do 400 jj = 1 , num250_per_1km

           j250 = (nc - 1) * num250_per_1km  + jj
           if (v250_dat(j250,i,k) .gt. 0.0) then
              pxl250(k,jj) = v250_dat(j250,i,k)
           else
              pxl250(k,jj) = bad_data
           end if

  400    continue
  500 continue

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(1x,'' pxldat 250m'',i5/,4(4f8.3/))') nc,
     +       ((pxl250(jj,kk),kk=1,num250_per_1km),jj=1,vis_band)
      endif
c ...............................................................

c     Get information about nearest 1-km pixel.

      do i = 1, num250_per_1km

        if(i .le. 12 .or. nc .eq. maxele) then
        
          if(scan_confdnc(nc) .gt. 0.95) then
            cm1km(i) = .true.
          else
            cm1km(i) = .false.
          end if
          if(qa_bitarray(nc,1) .eq. 0) then
            qa_250(i) = .false.
          else
            qa_250(i) = .true.
          end if
          ice_250(i) = scan_ice(nc)
          snow_250(i) = scan_snow(nc)
          desert_250(i) = scan_desert(nc)
          coast_250(i) = scan_coast(nc)
          land_250(i) = scan_land(nc)
          snglnt_250(i) = scan_snglnt(nc)
          day_250(i) = scan_day(nc)
          proc_250(i) = scan_process(nc)
          visusd_250(i) = scan_visusd(nc)

        else

          if(scan_confdnc(nc+1) .gt. 0.95) then
            cm1km(i) = .true.
          else
            cm1km(i) = .false.
          end if
          if(qa_bitarray(nc+1,1) .eq. 0) then
            qa_250(i) = .false.
          else
            qa_250(i) = .true.
          end if
          ice_250(i) = scan_ice(nc+1)
          snow_250(i) = scan_snow(nc+1)
          desert_250(i) = scan_desert(nc+1)
          coast_250(i) = scan_coast(nc+1)
          land_250(i) = scan_land(nc+1)
          snglnt_250(i) = scan_snglnt(nc+1)
          day_250(i) = scan_day(nc+1)
          proc_250(i) = scan_process(nc+1)
          visusd_250(i) = scan_visusd(nc+1)

        end if

      enddo

c---------------------------------------------------------------------

      return
      end
