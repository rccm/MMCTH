      subroutine put_proc_var_irphase(npts,line,pixel,ncloud,nmissing,
     +                                ngood_msk,indat_bt,var)

C!F77-------------------------------------------------------------
c!Description:
c	Routine which go through a scan cube of data and
c       gets the brightness temperature differences and scene
c       variance for a given box size.
c
c!Input parameters:
c npts		Counts number of boxes in scan cube
c line          Current processing line
c pixel         Current processing pixel
c ncloud        Number of cloudy pixels within box
c nmissing      Number of missing pixels within box (bad radiances)
c ngood_msk     Number of valid cloud mask pixels within box
c
c!Output Parameters:
c indat_bt	Array containing boxes of btdif values
c var		Array containing boxes of scene variance values
c
c!Revision History:
c
c!Team-unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c
c!Design Notes:
c!END-------------------------------------------------------------
      implicit none
      save

c ..  main include file for mod06uw
      include 'mod06uw_data.inc'
      include 'mod06uw_debug.inc'

c ... scalar arguments
      integer npts, line, pixel, ncloud, nmissing, ngood_msk

c ... array arguments      
c rhucek 5/8/01:  FORTRAN PARAMETER dayphase_band replaced with irphase_band 
c so that dimensions of dummy argument "var" agree with the array size passed by 
c caller (MOD_PR06UW).
c     real var(max_samp_pixel*max_samp_line,dayphase_band)
      real var(max_samp_pixel*max_samp_line,irphase_band)
      integer indat_bt(max_samp_pixel*max_samp_line,irphase_band+5)

c     local scalars
      real rad_avg, storad, sumsq, s2, rad_band
      integer ii, jj, nch, units, i,
     +        j, num_tot

c     external functions
      real modis_bright
      external modis_bright

c     intrinsic functions
      intrinsic real, nint

c ...   now process the data in 5x5 segments.  if more than
c ...   5 of the 25 pixels are found to have cloud in them, then
c ...   average  the radiances, convert to brightness temperatures and
c ...   save in array.
        units = 1
 
c ... debug statement ............................................
      if (debug .gt. 2) then
        write(h_output,'(2x/,''Subroutine put_proc_var  '',/)')
        write(h_output,'(2x,''line,pixel,npts '', 3i10)')
     +        line, pixel, npts 
        write(h_output,'(2x,''ncloud,nmissing '', 2i10,/)')
     +        ncloud, nmissing
      endif
c ................................................................

 
      if ( (view(pixel+isamp/2,line+isamp/2) .lt. Zen_Ang_Lim .and.
     +      view(pixel+isamp/2,line+isamp/2) .gt. 0.0) .and. 
     +     (ngood_msk .eq. (isamp*isamp)) ) then

        if(ncloud .ge. ncprct) then

c ...     Do for each channel
          do 155 nch = 1 , irphase_band 

            storad = 0.0
            num_tot = 0
            sumsq = 0.0

c ...       Do for each pixel inside box
            do 160 jj = 0,4
              do 165 ii = 0,4
                rad_band = radiance1(pixel+ii,line+jj,ir_aband(nch))

	        if (rad_band .gt. 0.) then
                  num_tot = num_tot + 1
                  storad = storad + rad_band
c ...             store radiance square sum for standard dev
                  sumsq = sumsq + rad_band**2
	        endif
	
  165         continue
  160       continue

            if(num_tot .gt. 0) then
              rad_avg = storad / real(num_tot)
              indat_bt(npts,nch) = 
     +          nint(100. * modis_bright(rad_avg,ir_mband(nch),units))
              s2 = ((real(num_tot) * sumsq) - (storad**2))
     +             / (real(num_tot) ** 2)
              if (s2 .lt. 0.0) then
                var(npts,nch) = 0.0
              else
                var(npts,nch) = sqrt(s2)
              endif
            endif

  155     continue

          indat_bt(npts,4) = indat_bt(npts,1) - indat_bt(npts,2)
          indat_bt(npts,5) = indat_bt(npts,2) - indat_bt(npts,3)
          indat_bt(npts,6) = nint(100.*var(npts,1))
          indat_bt(npts,7) = nint(100.*var(npts,2))
          indat_bt(npts,8) = nint(100.*var(npts,3))

        else

	  do i=1,irphase_band+5
            indat_bt(npts,i) = clr_data
	  enddo

        end if

c ....end of zenith angle limitation
      else

        do 190 i = 1 , irphase_band+5
          indat_bt(npts,i) = out_misg
 190    continue

      endif


c ... debug statement ............................................
      if (debug .gt. 2) then
        write(h_output,'(15x,'' INDAT_BT Line'', I4,'' AND PIXEL '',
     +                I4, '' POINTS SO FAR '',I10)') line, pixel, npts
        write(h_output,'(2x,'' indat_bt '',/,(8i10/))')
     +       ((indat_bt(i,j),j=1,8),i=npts,npts)
        write(h_output,'(2x,'' var  '',/,(3f10.2/))')
     +       ((var(i,j),j=1,3),i=npts,npts)
      endif
c ................................................................



      return
      end
