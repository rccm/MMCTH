      subroutine mod06CT_compute_products( line, pixel, ngood_msk,
     &      ncloudy, ncloud_co2, ngood_co2, met_date, qual_flag, conf_flag, 
     &      hc_flag, ci_flag, os_top_flag, cldhgt_cat, nearnad_flag,
     &      cldhgtmet_flag, ssctpr, slcd, smcd, shcd, sthncd, sthkcd,
     &      sopcd, scicd, sco2, sbias)

c-----------------------------------------------------------------------
c
c!F77
c
c!Description:
c    Compute MOD06CT retrieval products.
c
c!Input Parameters:
c    LINE             Line number within a swath (1-10)
c    PIXEL            Pixel number within a 1km scan (1-1500)
c    NGOOD_MSK        Number pixels where cloud mask was determined
c    NCLOUDY          Number cloudy pixels in box
c                     (from valid cloud mask pixels)
c    NCLOUD_CO2       Number cloudy pixels in box 
c                     (must be valid pixel for CO2-slicing retrieval)
c    MET_DATE         Year, month, day, hour of gridded met data
c
c    The following in COMMON /MOD06_DATA/ are used:
c    RADIANCE1        Radiances for IR bands
c    BOX_FLAG_CO2     Valid data flag array for box (0=bad,1=good)
c    BOX_CLOUD        Clear/cloud flag array for box (0=clear,1=cloud)
c    NGOOD_CO2        Number valid numbers
c
c!Output Parameters:
c    QUAL_FLAG        Cloud top pressure quality flag (0=unusable,1=usable)
c    CONF_FLAG        Cloud top pressure confidence flag (0-3, low to high)
c    HC_FLAG          Indicates high cloud found (0=missing, 1=no, 2=yes)
c    CI_FLAG          Indicates cirrus cloud found (0=missing, 1=no, 2=yes)
c    OS_TOP_FLAG      Indicates over-shooting thunderstorm top (0=missing, 1=no, 2=yes)
c    CLDHGT_CAT       Indicates low, mid, hi clouds for QA output
c    NEARNAD_FALG     Indicates near-nadir view (1=yes, 2=no)
c    CLDHGTMET_FLAG   Indicates CTP retrieval method
c    SSCTPR           Sum of good cloud pressure retrievals in granule
c    SLCD             Sum of low cloud retrievals in granule
c    SMCD             Sum of mid-level cloud retrievals in granule
c    SHCD             Sum of high cloud retrievals in granule
c    STHNCN           Sum of "thin" cloud retrievals in granule
c    STHKCD           Sum of "thick" cloud "retrievals in granule
c    SOPCD            Sum of opaque cloud retrievals in granule
c    SCICD            Sum of cirrus cloud retrievals in granule
c                              detected (0=missing, 1=no, 2=yes)
c
c    The following are output to COMMON /MOD06_DATA/:
c    CLOUD_FRACTION   Cloud fraction determined from cloud mask output
c    PROCESSING_FLAG  Flag indicating which technique was used to 
c                      retrieve cloud parameters
c
c!Revision History:
c
c calls:
c function modis_bright       Compute brightness temperature from radiances
c subroutine ialign           Order radiances from coldest to warmest
c subroutine get_co2cld       compute CO2-slicing cloud top pressure
c subroutine get_ct_stats     sum number occurrances of various cloud types
c
c!Team-unique Header:
c
c!End
c
c-----------------------------------------------------------------------
      
      implicit none
      save
      
      include 'mod06uw_data.inc'
      include 'mod06uw_debug.inc'
      
c ... arguments

      integer line, pixel, ngood_msk, ncloudy, qual_flag, conf_flag,
     &        hc_flag, ci_flag, ssctpr, slcd, smcd, shcd, sthncd,
     &        sthkcd, sopcd, scicd, met_date(4), ngood_co2,bias_flag,
     &        co2_flag, sco2, sbias, os_top_flag, cldhgt_cat, nearnad_flag,
     &        cldhgtmet_flag, ncloud_co2

c ... local variables

      integer mbnds(ntbct), abnds(ntbct), n, ipx, iout, proc_flag,
     &        jout, k, j, i, ii, jj, iord11(isamp*isamp),
     &        obs_clr, tmprad(isamp*isamp,nbct), ibuf11(isamp*isamp),
     &        mid_px, mid_ln
      real sum, rwarm(ntbct), rcold(ntbct), twarm(ntbct), tcold(ntbct),
     &     cld_frac, sumc, sumw, pct, eca, rad
      real sum2

c ... external functions

      real modis_bright
      external modis_bright

c-----------------------------------------------------------------------

c ... Data statements

c ... Define MODIS IR bands used
      data mbnds /36,35,34,33,31,29,32/

c ... Define indices of input radiance array corresponding to the needed
c ... MODIS IR bands
      data abnds /15,14,13,12,10,8,11/

c-----------------------------------------------------------------------

c ... Define output array indices for the current line and pixel.
      iout = (pixel / isamp) + 1
      jout = (line / isamp) + 1

c ... Initialize mean "cold" and "warm" brightness temperatures and
c ... radiances
      do k = 1,ntbct
        twarm(k) = bad_value
        rwarm(k) = bad_value
        tcold(k) = bad_value
        rcold(k) = bad_value
      end do
      proc_flag = 0

c-----------------------------------------------------------------------

c ... Check number of cloudy pixels in the box
      if(ncloudy .ge. ncprct) then

c ...   Compute box-average cloudy and/or clear radiances based on 
c ...   the processing option flag
        proc_flag = mod06ct_proc_opt

c-----------------------------------------------------------------------

c ...   Write debug information about input data.
        if(debug .gt. 0) then
c        if(line .eq. 6 .and. pixel .eq. 676) then
          write(h_output,'(1x)')
          write(h_output,'(''Routine mod06CT_compute_products'')')
          write(h_output,'(''Line, pixel are '',2i5)') line,pixel
          write(h_output,'(''Process option is '',i3)')
     &        mod06ct_proc_opt
          write(h_output,'(''Pixel quality (box_flag_co2)'')')
          write(h_output,'(5i5)') box_flag_co2
          write(h_output,'(1x)')
          write(h_output,'(''Input radiances (radiance1)'')')
          do k = 1,ntbct
            write(h_output,'(''Ch '',i5,'' Pos '',i5)')
     &        mbnds(k),abnds(k)
            write(h_output,'(1x,5f10.3)') ((radiance1(ii,jj,abnds(k)),
     &        ii=pixel,pixel+isamp-1),jj=line,line+isamp-1)
          end do
c        end if
        end if

c-----------------------------------------------------------------------

c ...   Get average radiances over entire box for the necessary bands.
c ...  6/12/04 HZ average cloud radiance only 

        do k = 1, ntbct
          n = 0
          sum = 0.0
          sum2 = 0.0
          do j = 1, isamp
            do i = 1, isamp
              if(box_flag_co2(i,j) .eq. 1) then
                if(box_cloud(i,j) .eq. 1) then
c               if(box_cloud(i,j) .eq. 0) then
                  n = n + 1
                  sum = sum + radiance1(pixel+i-1,line+j-1,abnds(k))
c                 write(*,'(''summing: '',4i5, 2f10.3)') k,i,j,n,radiance1(pixel+i-1,line+j-1,abnds(k)),sum
                end if
              end if
              sum2 = sum2 + radiance1(pixel+i-1,line+j-1,abnds(k))
            end do
          end do
          if(n .gt. 0) then
            rcold(k) = sum / n
            tcold(k) = modis_bright(rcold(k),mbnds(k),1)
          end if
c         write(*,'(''sum: '',2i5,2f10.3)') k, n, sum, rcold(k)
        end do

c-----------------------------------------------------------------------

        if(mod06ct_proc_opt .eq. 2) then

c ...     Get average clear-sky radiance over box for 5 CO2 bands
          do k = 1, nbct
            n = 0
            sum = 0.0
            do j = 1, isamp
              do i = 1, isamp
                if(box_flag_co2(i,j).eq.1 .and. box_cloud(i,j).eq.0)then
                  n = n + 1
                  sum = sum + radiance1(pixel+i-1,line+j-1,abnds(k))
                end if
              end do
            end do
            if(n .gt. 0) then
              obs_clr = 1
              rwarm(k) = sum / n
              twarm(k) = modis_bright(rwarm(k),mbnds(k),1)
            else
              obs_clr = 0
            end if
          end do

c-----------------------------------------------------------------------

        else if(mod06ct_proc_opt .eq. 3) then

c-----------------------------------------------------------------------

c ...     Get "warm" and "cold" values from the observed radiances

c ...     First, get a 1-D array of all usable radiances
          do k = 1,nbct
            ipx = 0
            do j = 1, isamp
              do i = 1, isamp
                if(box_flag_co2(i,j) .eq. 1) then
                  ipx = ipx + 1
                  rad = radiance1(pixel+i-1,line+j-1,abnds(k))
                  tmprad(ipx,k) = nint(rad * 100.0)
                end if
              end do
            end do
          end do

c ...     Order the 11 micron radiances from coldest to warmest
          call ialign(tmprad(1,kwc),ibuf11,iord11,ipx)

c ...     Get average of "warm" and "cold" radiances.
          do k = 1,nbct
            sumw = 0.0
            sumc = 0.0
            do i = 1,n_avgct
              sumc = sumc + tmprad(iord11(i+1), k) * 0.01
              sumw = sumw + tmprad(iord11(ipx-i), k) * 0.01
            enddo
            rcold(k) = sumc / n_avgct
            rwarm(k) = sumw / n_avgct
            tcold(k) = modis_bright(rcold(k),mbnds(k),1)
            twarm(k) = modis_bright(rwarm(k),mbnds(k),1)
          end do
             
c-----------------------------------------------------------------------

        end if

c-----------------------------------------------------------------------

c ...   Write debug information about CO2-slicing input data
        if(debug .gt. 0) then
c        if(line .eq. 6 .and. pixel .eq. 676) then
          write(h_output,'(1x)')
          write(h_output,'(''Input data to get_co2cld'')')
          write(h_output,'(''Cold and warm t(bb), radiances:'')')
          write(h_output,'(7f10.3)') (tcold(ii),ii=1,ntbct)
          write(h_output,'(7f10.3)') (twarm(ii),ii=1,ntbct)
          write(h_output,'(7f10.3)') (rcold(ii),ii=1,ntbct)
          write(h_output,'(7f10.3)') (rwarm(ii),ii=1,ntbct)
c        end if
        end if

c-----------------------------------------------------------------------

c ...   Compute CO2-slicing cloud top height for the box
        call get_co2cld ( line, pixel, obs_clr, tcold, twarm, 
     &    met_date, pct, eca, qual_flag, conf_flag, co2_flag, bias_flag,
     &    hc_flag, ci_flag, os_top_flag, cldhgt_cat, nearnad_flag,
     &    cldhgtmet_flag, ncloud_co2, ngood_co2 )

c-----------------------------------------------------------------------

c ...   Determine cloud category and compute sums for output metadata
        call get_ct_stats ( line, pixel, qual_flag, co2_flag, bias_flag,
     &    pct, eca, ssctpr, 
     &    slcd, smcd, shcd, sthncd, sthkcd, sopcd, scicd,sco2,sbias )

c-----------------------------------------------------------------------


      else  
       
c ...   Box is too clear - can't do cloud top height retrieval

c ...   Set cloud flags to indicate clear skies.
        hc_flag = 3
        ci_flag = 3
        cldhgt_cat = 1
        cldhgtmet_flag = 7
        if(lat1(pixel,line) .lt. 50.0) then
          os_top_flag = 1
        end if

      end if

c-----------------------------------------------------------------------

c ... Compute cloud fraction (areal coverage) based on the cloud mask.
      if(ngood_msk .eq. (isamp*isamp)) then
        cld_frac = float(ncloudy) / ngood_msk
        cloud_fraction(iout,jout) = cld_frac
        mid_px = pixel + isamp / 2
        mid_ln = line + isamp / 2
        if(view(mid_px,mid_ln) .le. near_nadir_vza_limit) then
          cloud_fraction_nearnad(iout,jout) = cloud_fraction(iout,jout)
        end if
      end if

c-----------------------------------------------------------------------

c ... Write debug information about output
      if(debug .gt. 0) then
c      if(line .eq. 6 .and. pixel .eq. 676) then
        write(h_output,'(1x)')
        write(h_output,'(''Routine mod06CT_compute_products'')')
        write(h_output,'(''Cloud fraction: '',3f10.3)') view(mid_px,mid_ln),
     *         cld_frac, cloud_fraction_nearnad(iout,jout)
c      end if
      endif

c-----------------------------------------------------------------------

      return
      end
