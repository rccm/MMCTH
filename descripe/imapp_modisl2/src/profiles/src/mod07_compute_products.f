      SUBROUTINE MOD07_COMPUTE_PRODUCTS(LINE, PIXEL, MONTH, LANDFRAC, 
     &  FLAG)

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C     Compute MOD07 retrieval products.
C
C !INPUT PARAMETERS:
C     LINE                  Line number within a 1km scan (1-10)
C     PIXEL                 Pixel number within a 1km scan (1-1500)
C     MONTH                 Month of observation (1-12)
C     LANDFRAC              Fraction of land in observation (0.0-1.0)
C
C     The following arrays in COMMON /MOD07_DATA/ are used:
C     RADIANCE1             Radiances for IR bands
C     PRE1                  Surface pressure
C     ZEN1                  Senzor zenith angle
C     LAND1                 Land/water flag
C     LAT1                  Latitude
C     ELV1                  Surface elevation
C     MASK1                 Cloudmask
C
C !OUTPUT PARAMETERS:
C     FLAG                  Success flag (0=Success, -1=Failure)
C
C     The following arrays in COMMON /MOD07_DATA/ are filled:
C     BRIGHTNESS_TEMP       Brightness temperature
C     RETR_TEMP_PROFILE     Retrieved temperature profile
C     RETR_DEWP_PROFILE     Retrieved dewpoint profile
C     RETR_WVMR_PROFILE     Retrieved water vapour mixing ratio profile
C     RETR_HITE_PROFILE     Retrieved geopotential height profile
C     RETR_OZONE_PROFILE     Retrieved ozone profile
C     SFC_TEMP              Surface temperature
C     SFC_PRES              Surface pressure
C     SFC_ELEV              Surface elevation
C     WATER_VAPOR           Water vapor for entire column,
C                           by integrating regression profiles
C     WATER_VAPOR_DIRECT    Water vapor for entire column,
C                           by direct regression
C     WATER_VAPOR_LOW       Water vapor for surface to 900 hPa
C     WATER_VAPOR_HIGH      Water vapor for 700 hPa to 300 hPa
C     TOTAL_OZONE           Total column ozone
C     TOTAL_TOTALS          Total totals index
C     LIFTED_INDEX          Lifted index
C     K_INDEX               K index
C
C !REVISION HISTORY:
C     26 February 2002: SWS made updates for delivery:
C         - Included global biases, read in from namelist file bias.nl
C           Biases are separated into three latitude zones (tropical, 
C           midlatitude, polar) each for land and ocean.  bias_version
C           character string is read in with bias file but currently is
C           not used anywhere. 
C         - Included new variable WATER_VAPOR_DIRECT
C         - Updated call to modis_ges101 to receive output tpw_ret, swemis
C         - Limit maximum BT for any band to 330K
C         - Skip retrieval for first ten pixel in each scan line
C
C     February & March 2003: SWS updated
C         - Used retrieved skin temperature instead of GDAS sfc T in sfc_temp
C         - Only skip first 10 pixels if for Terra (include platform_name.inc)
C         - Read bias from file
C         - Read in detector flags, skip detectors if flagged 0 
C         - Use average zenith angle, latitude, sfcP, elev, pwat of "good" pixels 
C           in 5x5 instead of center pixel.  This is important when some pixels 
C           in a 5x5 retrieval area are not retrieved (e.g., cloudy).
C         - Brightness temperature or TPW dependent biases are used instead of the
C           old latitude dependent ones
C         - init_bias parameter set so only read bias, detector flags once
c     May 2005: SWS fixed geopotential height profile calculation,
c               and fixed output of t, q profiles below ground
c     May 2006: EvaBorbas added ozone profiles to the outputs and changed the pressure level to 101
c 
c     March 2007 EvaBorbas changing back to 20 level output
c     Jan 2009: dewpoint temperature and mixing ratio profiles are in the output now 
c     Jan 2010: EvaBorbas removed ozone profiles from the outputs 
c     May 2010: EvaBorbas changed modis_bright to modis_bright_shift
c     Aug  2010: EvaB, Ln 524 & 530 added: mixr and dewT is written out if mixr > 0.001 
c                               this results the top 4 levels filled with missing value
c
C !TEAM-UNIQUE HEADER:
C     Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
C !END
C-----------------------------------------------------------------------
      
      IMPLICIT NONE
      
      include 'mod07_data.inc'
      include 'mod07_pcfnum.inc'
      include 'debug.inc'
      include 'platform_name.inc'

c ... Arguments
      integer line, pixel, month, flag
      real landfrac
      
c ... Local variables
      real rad, rsum, ravg, psfc, elev, tsfc, wsfc, zsfc
      integer i, j, k, n, units, iok
      integer iout, jout, outindex(outlevels), ntemp, idir
      logical goodbox

      real tbrt(36)
      integer nlevels, bias_flag, dummy_flag
      parameter (nlevels = 101)
      real pres(nlevels), temp(nlevels), wvmr(nlevels)
      real z101(nlevels), ozon(nlevels)
      real alat, azen, amonth, tpw, toz, swemis
      integer lsfc, lsfc40, errflg
      character*100 errtxt
      character*20 bias_version

c ... New for bias, detector flag read 
      integer ncols_bias, ncols_det 
      integer num_bands, num_cols
      integer pcf_num_bias, pcf_num_det
      parameter (num_cols = 10)
      character*100 header_bias, header_det  
      character*4 platform
      character*5 platform_coeff
      real bias_data(num_cols, nbands), bias 
      real det_flag(num_cols, nbands)
      integer band_data(nbands)
      logical init_bias
      real gdas_pwat

c ... new for zenith angle, lat, sfcP, elev average
      real zensum, alatsum, sfcpsum, elevsum, pwatsum

      real tpw_ret, skint
      
      integer nstd
      parameter (nstd = 40)
      integer iuc
      parameter (iuc = 100)
      real pref(nstd), pstd(nstd), tstd(nstd), wstd(nstd), 
     &  zstd(nstd), pwvp(nstd), ostd(nstd)
      
      integer level
      real index_lift, index_tt, index_ki, tpw_lo, tpw_hi
      
c ... External functions
c      real modis_bright, dewtem
c      external modis_bright, dewtem
      real modis_bright_shift, dewtem
      external modis_bright_shift, dewtem
      integer get_coeff
      external get_coeff

c ... Save statements

      SAVE init_bias, bias_flag, bias_data, ncols_bias, det_flag

c ... Data statements

      DATA init_bias / .true. /

c ... Standard pressure profile (hPa) at 40 levels
      data pref /.1,.2,.5,1.,1.5,2.,3.,4.,5.,7.,10.,15.,20.,25.,30.,50.,
     &  60.,70.,85.,100.,115.,135.,150.,200.,250.,300.,350.,400.,430.,
     &  475.,500.,570.,620.,670.,700.,780.,850.,920.,950.,1000./

c ... Indices of output levels in standard pressure profile
      data outindex / 9, 11, 13, 15, 16, 18, 20, 23, 24, 25, 26,
     &  28, 31, 33, 35, 36, 37, 38, 39, 40 /

      if ( init_bias ) then

c-----------------------------------------------------------------------
c     Read in bias values and detector flags - Only Once
c-----------------------------------------------------------------------

c ... Assign pcf number and platform string 

c rhucek 06/11/03: replaced if-elseif block with similar logic that 
c        assumes common pcfnum's for Terra and Aqua detector flags and 
c        radiance bias tables.
c
c     if (platform_name(1:5) .eq. 'Terra') then
c       pcf_num_bias = biasT_pcfnum
c       pcf_num_det = detT_pcfnum
c       platform = 'Terr'
c     elseif (platform_name(1:4) .eq. 'Aqua') then
c       pcf_num_bias = biasA_pcfnum
c       pcf_num_det = detA_pcfnum
c       platform = 'Aqua'
c     endif

      pcf_num_bias = radbias_pcfnum
      pcf_num_det  = detflgs_pcfnum
      platform     = 'Terr'
      if (platform_name(1:4) .eq. 'Aqua') platform = 'Aqua'

c ... Read bias file
      errflg = get_coeff(pcf_num_bias, nbands, num_cols, 
     &  header_bias, platform_coeff, bias_flag, num_bands, 
     &  ncols_bias, band_data, bias_data, errtxt) 

c ...   Check error status
        if (errflg .ne. 0) call message('mod07_compute_products',
     &    'Error in get_coeff.f: ' // errtxt // char(10) //
     &    '[OPERATOR ACTION: Verify name and format of bias ' //
     &    ' input files. If error persists, contact SDST]',
     &    errflg, 2)

      errflg = 0
      if (nbands. ne. num_bands) then
         errflg = -1 
      elseif (platform .ne. platform_coeff(1:4)) then
         errflg = -2 
      else 
         do j = 1,nbands
            if (band_data(j) .ne. bands(j)) then
               errflg = -3
            endif
         enddo
      endif
 
        if (errflg .ne. 0) call message('mod07_compute_products',
     &    'Error in number of bands, band list, or platform ' //
     &    ' name in bias input file [OPERATOR ACTION: Verify ' //
     &    ' format of bias files. If error persists, contact SDST]',
     &    errflg, 2)

      errflg = get_coeff(pcf_num_det, nbands, num_cols, header_det, 
     &  platform_coeff, dummy_flag, num_bands, ncols_det, 
     &  band_data, det_flag, errtxt) 

c ...   Check error status
        if (errflg .ne. 0) call message('mod07_compute_products',
     &    'Error in get_coeff.f: ' // errtxt // char(10) //
     &    '[OPERATOR ACTION: Verify name and format of detector ' //
     &    ' flag files. If error persists, contact SDST]',
     &    errflg, 2)

      errflg = 0
      if (nbands. ne. num_bands) then
         errflg = -1 
      elseif (platform .ne. platform_coeff(1:4)) then
         errflg = -2 
      else
         do j = 1,nbands
            if (band_data(j) .ne. bands(j)) then
               errflg = -3
            endif
         enddo
      endif
 
        if (errflg .ne. 0) call message('mod07_compute_products',
     &    'Error in number of bands, band list, or platform ' //
     &    ' name in detector  flag file [OPERATOR ACTION: Verify ' //
     &    ' format of detector flag files. If error persists, contact SDST]',
     &    errflg, 2)

c ...   Unset initialization flag

        init_bias = .false.

      endif

c-----------------------------------------------------------------------
c     Compute box average radiance and brightness temp in 1km ir bands
c-----------------------------------------------------------------------

c ... Set default return flag
      flag = -1
      
c ... Initialize brightness temperature array for this box
      do i = 1, 36
        tbrt(i) = -1.0
      end do

c ... Loop over ir bands
      do k = 1, nbands

c ...   Sum the good radiances within this box
        rsum = 0.0

c .. march 03: new parameters for zen, lat, sfcp, elev average
        zensum = 0.0
        alatsum = 0.0
        sfcpsum = 0.0
        elevsum = 0.0
        pwatsum = 0.0

        n = 0
        do j = 1, isamp
          do i = 1, isamp

c ... 19feb03: added screening for bad detectors
            if (det_flag(line + j - 1,k) .eq. 1) then

            rad = radiance1(pixel + i - 1, line + j - 1, k)
            if (box_flag(i, j) .eq. 1 .and. rad .gt. 0.0) then
              rsum = rsum + rad
              n = n + 1

c ... march 03: added zenith angle sum to use average instead of middle
c ...          also use average latitude, sfcp, elev, pwat
                if (k.eq.nbands) then
                   zensum = zensum + zen1(pixel+i-1,line+j-1)
                   alatsum = alatsum + lat1(pixel+i-1,line+j-1)
                   sfcpsum = sfcpsum + pre1(pixel+i-1,line+j-1)
                   elevsum = elevsum + elv1(pixel+i-1,line+j-1)
                   pwatsum = pwatsum + pwat1(pixel+i-1,line+j-1)
                endif
              endif   
            endif
          end do
        end do

c ...   If sufficient good pixels were found,
c ...   compute average radiance and convert to brightness temp
        if (n .ge. min_ngood .and. n .gt. 0) then
          ravg = rsum / real(n)
          units = 1
          tbrt(bands(k)) = modis_bright_shift(ravg, bands(k), units)
        endif

      end do

c ... march 03: Compute average zenith angle, psfc, pwat, etc.
c ...   This uses n from the last band (nbands); it will
c ...   only be different from other bands if different
c ...   pixels had rad < 0
      if (n .ge. min_ngood .and. n .gt. 0) then
         azen = zensum / real(n)
         psfc = sfcpsum / real(n)
         gdas_pwat = pwatsum / real(n)
         elev = elevsum / real(n)
         alat = alatsum / real(n)
      endif

c ... Check the resulting brightness temperatures for this box
      goodbox = .true.

c ... 4feb02: added check for BT > 330

      do k = 1, nbands
        if (tbrt(bands(k)) .le.   0.0) goodbox = .false.
        if (tbrt(bands(k)) .gt. 330.0) goodbox = .false.
      end do

c ... 26 feb02: do not retrieve first ten pixels in each scan line
c ... 14 jan03: do this only for terra
      if (platform_name(1:5) .eq. 'Terra') then
         if (pixel.lt.11) goodbox = .false.
      endif

c-----------------------------------------------------------------------
c     If data in this box are good, compute retrieval
c-----------------------------------------------------------------------

      if (goodbox) then

c ...   Get surface ancillary data values

c ...   psfc, elev, azen, alat now computed above (average instead of middle)
c        psfc = pre1(pixel + isamp / 2, line + isamp / 2)
c        elev = elv1(pixel + isamp / 2, line + isamp / 2)
c        azen = zen1(pixel + isamp / 2, line + isamp / 2)
c        alat = lat1(pixel + isamp / 2, line + isamp / 2)

        amonth = real(month)
                
c-----------------------------------------------------------------------
c  march 03: TPW or BT dependent biases
c-----------------------------------------------------------------------
      do k = 1,nbands
        if (bias_flag .eq. 1) then
          bias = bias_data(1,k) * tbrt(bands(k)) + bias_data(2,k) 
        elseif (bias_flag .eq. 2) then
          bias = bias_data(1,k) * gdas_pwat + bias_data(2,k)
        endif
 
c ...  If there are 4 columns in input bias file, check that bias 
c ...  is within min (3rd column) and max (4th column)
        if (ncols_bias .eq. 4) then
          if (bias .lt. bias_data(3,k)) then
            bias = bias_data(3,k)
          endif
          if (bias .gt. bias_data(4,k)) then
            bias = bias_data(4,k)
          endif
        endif

        tbrt(bands(k)) = tbrt(bands(k)) - bias
      enddo

c ...   Temperature and moisture regression retrieval
c ...   (includes new output parameters tpw_ret, swemis, skint)
        call modis_ges101(tbrt, psfc, landfrac, azen, alat, amonth, 
     &    iuctw, iuang, pres, temp, wvmr, ozon,tpw, toz, lsfc, errflg, 
     &    errtxt, swemis, tpw_ret, skint)

c ...   Check error status
        if (errflg .ne. 0) call message('mod07_compute_products',
     &    'Error in modis_ges101.f: ' // errtxt // char(10) //
     &    '[OPERATOR ACTION: Verify name of view angle and ' //
     &    ' coefficient files. If error persists, contact SDST]',
     &    errflg, 2)

c ...   Set return flag to indicate successful retrieval
        flag = 0

      endif

c-----------------------------------------------------------------------
c     Compute derived products (surface level handling could be better)
c-----------------------------------------------------------------------

      if (goodbox) then
      
c ...   Interpolate profiles to the 40 standard levels
        do i = 1, nstd
          pstd(i) = pref(i)
        end do
        call interp(nlevels, pres, temp, nstd, pstd, tstd)
        call interp(nlevels, pres, wvmr, nstd, pstd, wstd)
        call interp(nlevels, pres, ozon, nstd, pstd, ostd)

c ...   Compute stability indices
        level = 40
        call indices(pstd, tstd, wstd, psfc, temp(lsfc), wvmr(lsfc),
     &    level, index_lift, index_tt, index_ki)
 
c ...   Compute geopotential heights
        idir = 1

c ... SWS May 05: Find surface level in 40 level pressure profile
c ...  require that surface pressure be >= surface level (i.e., not below) 
      lsfc40 = 1
      do i = 1, nstd
         if (psfc .ge. pstd(i) ) lsfc40 = i
      end do

c ... SWS May 05: changed nstd to lsfc40 to fix bug in geopotential height
c ...  note that now zstd will only be filled with values to lsfc40
c ...  level.
        call gphite(pstd, tstd, wstd, elev, lsfc40, idir, zstd)
        call gphite(pres, temp, wvmr, elev, lsfc, idir, z101)

c ...   Precipitable water; low level (920 hPa to surface)
c        do i = 1, nstd
c          pwvp(i) = 0.0
c        end do
c        n = 3
c        call precwv(pstd(38), wstd(38), n, pwvp)
c        tpw_lo = pwvp(n)

c ...   Precipitable water; high level (300 hPa to 700 hPa)
c        do i = 1, nstd
c          pwvp(i) = 0.0
c        end do
c        n = 10
c        call precwv(pstd(26), wstd(26), n, pwvp)
c        tpw_hi = pwvp(n)

c Nov 2009 EvaB: computation of tpw_lo and tpw_hi is computed now from the 101 level
c                profile data and the new tprecw is used as it is used for tpw now

c ...  EvaB Precipitable water; low level (680 hPa to surface)
c                     top level of integration is 85 (683.667 hPa)
       call tprecw(pres, wvmr, tpw_lo, nlevels, 85, lsfc, psfc)

c ...   Precipitable water; high level (10 hPa to 440 hPa)
c                     top level of integration is 21 (9.5 hPa)
c                     bottom level of integration is 73 (441.882 hPa)  (64 (300.00 hPa))
      call tprecw(pres, wvmr, tpw_hi, nlevels, 21, 73, pres(73))

      endif
      
c-----------------------------------------------------------------------
c ... Store results in product arrays
c-----------------------------------------------------------------------

c ... Pixel and line coordinates in input 1km array (center of box)
      i = pixel + isamp / 2
      j = line + isamp / 2

c ... Pixel and line coordinates in output 5x5 sampled array
      iout = (pixel / isamp) + 1
      jout = (line / isamp) + 1

c ... Copy the brightness temperatures
      do k = 1, nbands
        brightness_temp(iout, jout, k) = tbrt(bands(k))
      end do

c ... Copy the products
      if (goodbox) then
      
c ...   Parameters copied from input ancillary data
        sfc_pres(iout, jout) = psfc
        sfc_elev(iout, jout) = elev
c 1 oct 02: sws changed sfc_temp from GDAS input temp(lsfc) to skint retrieved
c        sfc_temp(iout, jout) = temp(lsfc)
        sfc_temp(iout, jout) = skint

c ...   26 feb 02: added water_vapor_direct 
c ...   Retrieval products
        water_vapor(iout, jout)        = tpw
        water_vapor_direct(iout, jout) = tpw_ret
        water_vapor_low(iout, jout)    = tpw_lo
        water_vapor_high(iout, jout)   = tpw_hi
        total_totals(iout, jout)       = index_tt
c rhucek 05/20/02: commented next line at request of Liam Gumley
c       total_totals(iout, jout)       = tpw_ret
        lifted_index(iout, jout)       = index_lift
        k_index(iout, jout)            = index_ki
        total_ozone(iout, jout)        = toz

c ...   Retrieval profile at levels above surface
        do k = 1, outlevels
c ... SWS changed 12 May 2005 to set values to only those levels
c ...  above the surface.
          if (outindex(k) .le. lsfc40) then
            retr_temp_profile(iout, jout, k) = 
     &        tstd(outindex(k))

c   evab added on 08.26.2010 to fix dew point temperature and mix ratio profile above 50 hPa
	    if (wstd(outindex(k)).gt.0.001) then
	    	 retr_dewp_profile(iout, jout, k) =
     &             dewtem(pstd(outindex(k)), tstd(outindex(k)),
     &             wstd(outindex(k)))
                 retr_wvmr_profile(iout, jout, k) =
     &               wstd(outindex(k))     
           endif
	    
            retr_hite_profile(iout, jout, k) = 
     &        zstd(outindex(k))
c            retr_ozone_profile(iout, jout, k) = 
c     &        ostd(outindex(k))

          endif
c
c may 2006 Eva Borbas modification for 101 level and ozone profile is added
c
c          if (pres(k) .le. pres(lsfc)) then
c            retr_temp_profile(iout, jout, k) = 
c     &        temp(k)
c            retr_wvmr_profile(iout, jout, k) =
c     &        dewtem(pres(k), temp(k), wvmr(k))
c           retr_wvmr_profile(iout, jout, k) =
c     &        wvmr(k)
c            retr_hite_profile(iout, jout, k) = 
c     &        z101(k)
c            retr_ozone_profile(iout, jout, k) = 
c     &        ozon(k)
c          endif
        end do

c ...   Variables not computed in this version of the algorithm
c ...   (set to missing value on output)
c ...   HEIGHT_TROPOPAUSE, GUESS_TEMP_PROFILE, GUESS_DEWP_PROFILE

      endif
      
c-----------------------------------------------------------------------
c ... Write debug information for nadir retrieval box
c-----------------------------------------------------------------------

      if (debug .eq. 1) then

        if (line .eq. 6 .and. pixel .eq. 676) then

          write(h_output, '(/,a,/)') 'MOD07_COMPUTE_PRODUCTS DEBUG INFO'

          write(h_output, '(''Input Data for Nadir Box'')')
          write(h_output, '(''IR Bands: '',20i7)')
     &      (bands(i), i = 1, nbands)
          write(h_output, '(''Avg Temp: '',20f7.2)')
     &      (tbrt(bands(i)), i = 1, nbands)

          if (goodbox) then
          
            write(h_output, '(/,''Retrieval for Nadir Box'')')
            write(h_output, '(''Pressure Levels: '',20f8.1)')
     &        (pstd(i), i = 40, 26, -1)
            write(h_output, '(''Retrieval Temp.: '',20f8.1)')
     &        (tstd(i), i = 40, 26, -1)
            write(h_output, '(''Retrieval Mixr.: '',20f8.3)')
     &        (wstd(i), i = 40, 26, -1)
            write(h_output, '(''Retrieval Hite.: '',20f8.1)')
     &        (zstd(i), i = 40, 26, -1)
c           write(h_output, '(''Retrieval Ozon.: '',20f7.1)')
c    &        (ostd(i), i = 40, 26, -1)
c
c may 2006 Eva Borbas modification for 101 level and ozone profile is added
c
c            write(h_output, '(''Pressure Levels: '',101f7.1)')
c     &        (pres(i), i = 101, 1, -1)
c            write(h_output, '(''Retrieval Temp.: '',101f7.1)')
c     &        (temp(i), i = 101, 1, -1)
c            write(h_output, '(''Retrieval Mixr.: '',101f7.1)')
c     &        (wvmr(i), i = 101, 1, -1)
c            write(h_output, '(''Retrieval Hite.: '',101f7.1)')
c     &        (z101(i), i = 101, 1, -1)
c            write(h_output, '(''Retrieval Ozon.: '',101f7.1)')
c     &        (ozon(i), i = 101, 1, -1)
          
            write(h_output, '(''Lifted Index          ='',f7.2)') index_lift
            write(h_output, '(''Total-totals Index    ='',f7.2)') index_tt
            write(h_output, '(''K Index               ='',f7.2)') index_ki
            write(h_output, '(''Tot. Precip. Water    ='',f7.2)') tpw
            write(h_output, '(''Tot. Precip. Water Int='',f7.2)') tpw_ret
            write(h_output, '(''Tot. Precip. Water Lo ='',f7.2)') tpw_lo
            write(h_output, '(''Short Wave Emissivity ='',f7.2)') swemis
            write(h_output, '(''Tot. Precip. Water Hi ='',f7.2)') tpw_hi
            write(h_output, '(''Total Ozone           ='',f7.2)') toz
            write(h_output, '(''Surface Skin Temp. ='',f7.2)') skint

c may 2006 Eva Borbas modification for 101 level and ozone profile is added

             write(h_output, '(''Surface Level: '',2i6)') lsfc,lsfc40


            write(h_output, '(/,''Output profile for Nadir Box'')')
            write(h_output, '(''Pressure Levels: '',20f8.1)')
     &        (pstd(outindex(k)), k = outlevels, 1, -1)
c     &        (pres(k), k = outlevels, 1, -1)
            write(h_output, '(''Retrieval Temp.: '',20f8.1)')
     &        (retr_temp_profile(iout, jout, k), k = outlevels, 1, -1)
            write(h_output, '(''Retrieval Dewp.: '',20f8.1)')
     &        (retr_dewp_profile(iout, jout, k), k = outlevels, 1, -1)
            write(h_output, '(''Retrieval Mixr.: '',20f8.3)')
     &        (retr_wvmr_profile(iout, jout, k), k = outlevels, 1, -1)
            write(h_output, '(''Retrieval Hite.: '',20i8)')
     &        (nint(retr_hite_profile(iout, jout, k)), k = outlevels, 1, -1)
c            write(h_output, '(''Retrieval Ozon.: '',20f7.1)')
c     &        (retr_ozone_profile(iout, jout, k), k = outlevels, 1, -1)

          else
          
            write(h_output, '(''Input Data Bad: No Retrieval'')')
            
          endif       
               
        endif

      endif

c-----------------------------------------------------------------------
c ... Set success flag (0=Success, -1=Failure)
c-----------------------------------------------------------------------

      if (goodbox) flag = 0
      
      END
