      SUBROUTINE MODIS_GES101(BRT, PSURF, PLAND, THETA, ALAT, AMONTH,
     &  IUCTW, IUANG, PRES, TEMP, WVMR, OZON,TPW, TOZ, LSFC, ERRFLG, 
     &  ERRTXT, SWEMIS, TPW_RET, SKINT)

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C     MODIS regression retrieval of atmospheric profiles of
C     temperature, moisture, and ozone at 101 levels. In addition,
C     surface skin temperature, microwave emissivity, and IR surface
C     emissivities and reflectivities are retrieved with 1km resolution.
C     Algorithm developed by Jun Li, CIMSS/SSEC.
C
C !INPUT PARAMETERS:
C     BRT      Brightness temperatures (K) for all 36 MODIS bands
C              (caller should set to 0.0 for bands 1-19 and 26).
C              Local array KUSE determines which bands are used.
C     PSURF    Surface pressure (hPa)
C     PLAND    Fraction of land (0.0-1.0)
C     THETA    Local view angle (aka senzor zenith angle, degrees)
C     ALAT     Latitude (degrees, -90S to +90N)
C     AMONTH   Month of year (1-12)
C     IUCTW    FORTRAN unit number for regression coefficient file
C              MODIS_REGCOEF.bin (caller must open this file
C              in direct access mode with reclen=1340 bytes - was 1276)
C     IUANG    FORTRAN unit number for view angle file
C              MODIS_senzen.bin (caller must open this file
C              in direct access mode with reclen=2720 bytes)
C
C !OUTPUT PARAMETERS:
C     PRES     Pressure levels for retrieved profiles (hPa)
C     TEMP     Retrieved profile of atmospheric temperature (K)
C     WVMR     Retrieved profile of water vapor mixing ratio (g/kg)
C     O    Retrieved profile of ozone mixing ratio (g/kg)
C     TPW      Total precipitable water vapor (cm)
C     TOZ      Total ozone (Dobsons)
C     LSFC     Level index of surface in pressure profile
C     ERRFLG   Error flag (0=success, 1=failure)
C     ERRTXT   Error text
C     SWEMIS   Short wave emissivity used in retrieval
C     TPW_RET  Total precipitable water vapor (cm) retrieved directly
C              instead of integrated from a retrieved profile
C
C !REVISION HISTORY:
C     February 2002: SWS made updates for delivery:
C        - Now uses 11 predictors instead of 12.  Instead of individual
C          bands 24,25,27-36, predictors are brightness temperature 
C          difference between 25-24, and individual bands 27 through 36
C        - Included 7 BT classes: perform retrieval using only that portion
C          of the training data set with band 31 BT in the same class.
C          Classes are < 245, 245-269, 269-285, 285-294, 294-300, 300-310, > 310 
C        - Included output parameters swemis and tpw_ret. 
C
C    May 2003: SWS updated
C        - Call to tprecw_new2.f instead of tprecw.f
C        - Added retrieved surface skin temperature variable skint  
C        - Removed band 24 from the regression (kuse), and removed band 25-24 difference
C        - Updated version number to agree with new coefficients, 2003073
C
C    April 2004: SWS updated to accommodate coefficients with 15 emissivities (reclen 1340)
C        - Updated version number to agree with new coefficients, 2004122 (May 1, 2004) 
C
C    April 20, 2004: SWS modified for land/ocean partitioning
C    April 21, 2003: SWS included new zones for land (4) and ocean (3)
C           because after partitioning there was not enough training
C           in all zones.
C       
C    May 2006 delivery:
C     SWS updated to accommodate coefficients with 8 emissivities
C     inflection point wavelengths, [3.7 4.3 4.5 7.6 8.3 9.3 10.8 14.3]
C        - reclen 1284
C        - Updated version number to agree with new coefficients, 2006135
C        - swemis indexing
C
C    May 2006 Eva Borbas: adding ozone profile to the output, changing the output level for 101
c
c    Nov 2006 Eva Borbas: no Ocean Zones, changing emis from 8 to 10 values 
c    July 2007 SWS: putting back in ocean zones
C    Sept 2009 EvaB  Updated version number to agree with new coefficients,2009260
C    Jan 2010 EvaB  updated the reclen =1308 and nemis=11 (for the MODIS used bands)
c    Aug 2010 EvaB  min value (0.001) is added to mixr profile (ln 450). 

C !TEAM-UNIQUE HEADER:
C     Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.C
C !END
C-----------------------------------------------------------------------

      implicit none
      
c ... Include files
      include 'platform_name.inc'
            
c ... Parameters
      integer nl, nemis, nall
c sws May 06 delivery : changed nemis from 15 to 8
c EvaB Nov 2006 : changed nemis from 8 to 10
      parameter (nl = 101, nemis = 11, nall = 36)
      integer leng
      parameter (leng = 3 * nl + 2 + 2 * nemis)
      integer nxx,nx

C --------------------
C	Number of predictors, nxx = 11 for bands 25, 27-36
C --------------------
      parameter (nxx=11, nx = 2*nxx + 4)
      integer ntwx, ntwy
      parameter (ntwx = nx, ntwy = leng)
c sws 22 april 04: added separate land/ocean numcls
      integer numang, numLandOcean
      integer numcls_land, numcls_ocean

      parameter (numang = 680)
c EvaB Nov 2006 : changed numcls_ocean from 3 to 1
c SWS 7/06: changing back
      parameter (numcls_land = 4, numcls_ocean = 3)
      integer nx1, nrx
      parameter (nx1 = nx + 1, nrx = nx1 * numang * 
     &    (numcls_land+numcls_ocean))

c ... Array arguments
      real brt(nall), pres(nl), temp(nl), wvmr(nl)

c ... Scalar arguments
      real psurf, pland, theta, alat, amonth, tpw, toz
      integer iuctw, iuang, lsfc, errflg
      character*(*) errtxt

c ... new variables for directly retrieved TPW, skin temperature
      real tpw_ret, skint

      integer land_flag

c ... Local scalars 
      integer i, ii, ixx, j, l, lang, init, ios, reclen, recnum
      real sum, mindiff, thisdiff, wmax
      logical opened

c ... Local arrays 
      real angp(numang), coef(ntwx, ntwy), coeftw(nrx, ntwy),
     &  xtw(ntwx), ytw(ntwy), ytwb(ntwy), header(ntwy)
      real rtv(leng), ozon(nl), profile(nl)
      real xx
      integer kuse(nall)
      integer ncls
      real swemis

c ... External functions
      real satmix
      external satmix

c ... Save statement
      save
      
c ... Data statements

C------------------------------------------------------
c ... Band use flags (1=use band, 0=don't use band)
C------------------------------------------------------
      data kuse/23*0, 0, 1, 0, 1, 9*1/

c ... Initialization flag (1=initialize, 0=don't initialize)      
      data init/1/

c ... Pressure profile at 101 levels (hPa)
      data profile/0.005,.016,.038,.077,.137,.224,.345,.506,.714,
     &  .975,1.297,1.687,2.153,2.701,3.340,4.077,4.920,
     &  5.878,6.957,8.165,9.512,11.004,12.649,14.456,16.432,
     &  18.585,20.922,23.453,26.183,29.121,32.274,35.651,39.257,
     &  43.100,47.188,51.528,56.126,60.989,66.125,71.540,77.240,
     &  83.231,89.520,96.114,103.017,110.237,117.777,125.646,133.846,
     &  142.385,151.266,160.496,170.078,180.018,190.320,200.989,212.028,
     &  223.441,235.234,247.408,259.969,272.919,286.262,300.000,314.137,
     &  328.675,343.618,358.966,374.724,390.893,407.474,424.470,441.882,
     &  459.712,477.961,496.630,515.720,535.232,555.167,575.525,596.306,
     &  617.511,639.140,661.192,683.667,706.565,729.886,753.628,777.790,
     &  802.371,827.371,852.788,878.620,904.866,931.524,958.591,986.067,
     &  1013.948,1042.232,1070.917,1100.000/
     
c-----------------------------------------------------------------------
c     INITIALIZATION
c-----------------------------------------------------------------------

c ... Intialize the output array
      do i = 1, leng
        rtv(i) = 0.0
      end do

c ... Set error flag and text
      errflg = 0
      errtxt = ' '
      
c ... Read view angle and regression coefficient files if required
      if (init .eq. 1) then

c ...   Check that view angle file was opened with
c ...   correct record length (2720 bytes)
        inquire(iuang, opened=opened, recl=reclen)
        if (.not. opened) then
          errflg = 1
          errtxt = 'View angle file is not open'
        endif
        if (reclen .ne. 2720) then
          errflg = 1
          errtxt = 'View angle file record length is incorrect'
        endif
        
c ...   Read view angle file
        recnum = 1
        read (iuang, rec=recnum, iostat=ios) angp
        if (ios .ne. 0) then
          errflg = 1
          errtxt = 'Error reading view angle file'
          return
        endif

c ...   Check that regression coefficient file was opened with
c ...  NEW SWS, May 2006 delivery: record length is 1284 for 8 emis
c ...  EvaB, Nov 2006: record length is 1300 for 10 emis
c ...  EvaB, Jan 2010: record length is 1308 for 11 emis
        inquire(iuctw, opened=opened, recl=reclen)
        if (.not. opened) then
          errflg = 1
          errtxt = 'Regression coefficient file is not open'
        endif
        if (reclen .ne. 1308) then
          errflg = 1
          errtxt = 'Regression coefficient file record length is incorrect'
        endif

c ...   Read the header record
        recnum = 1
        read(iuctw, rec=recnum, iostat=ios) header
        if (ios .ne. 0) then
          errflg = 1
          errtxt = 'Error reading regression coefficient file'
          return
        endif

c ...   Check the following items in the header record
c ...   (1) Satellite number (1=Terra, 2=Aqua)
c ...   (2) Creation date (YYYYDDD, Current version=2004122)
c ...   (NOTE: 'platform_name' is defined in platform_name.inc')
c ...   May 2006 delivery: sws updated creation date

        if (nint(header(1)) .ne. 1 .and.
     &      nint(header(1)) .ne. 2) then
          errflg = 1
          errtxt = 'Invalid satellite number in regr coefficient file'
          return
        endif
        if ((platform_name(1:5) .eq. 'Terra' .and. nint(header(1)) .eq. 2) .or.
     &      (platform_name(1:4) .eq. 'Aqua'  .and. nint(header(1)) .eq. 1)) then
          errflg = 1
          errtxt = 'Mismatch between platform name and satellite number in regression coefficient file'
          return
        endif          
c        if (nint(header(2)) .ne. 2006135) then
        if (nint(header(2)) .ne. 2009260) then
          errflg = 1
          errtxt = 'Invalid version number in regression coefficient file'
          return
        endif          
        
c ...   Read the regression coefficients
        do i = 1, nrx

c ...     Read this record
          recnum = i + 1
          read (iuctw, rec=recnum, iostat=ios) ytw
          if (ios .ne. 0) then
            errflg = 1
            errtxt = 'Error reading regression coefficient file'
            return
          endif

c ...     Save the coefficients          
          do j = 1, ntwy
            coeftw(i, j) = ytw(j)
          end do
          
        end do

c ...   Unset the initialization flag
        init = 0

      end if

c-----------------------------------------------------------------------
c     RETRIEVAL
c-----------------------------------------------------------------------

c ... Get coefficient index for this view angle
      lang = 1
      mindiff = abs(angp(1) - theta)
      do i = 2, numang
        thisdiff = abs(angp(i) - theta)
        if (thisdiff .lt. mindiff) then
          mindiff = thisdiff
          lang = i
        endif
      end do

C Add BT Classes

c ... Get the class number: ncls
      xx = brt(31) 

C Set land_flag to 0 (ocean) or 1 (land)
C Different zones for land and ocean

      if (pland.gt.0.5) then
        land_flag = 1

        if(xx .lt. 272.0) ncls = 1
        if(xx .ge. 272.0 .and. xx .lt. 287.0) ncls = 2
        if(xx .ge. 287.0 .and. xx .lt. 296.0) ncls = 3
        if(xx .ge. 296.0) ncls = 4

      end if

      if (pland.le.0.5) then
        land_flag = 0

        if(xx .lt. 283.5) ncls = 1
        if(xx .ge. 283.5 .and. xx .lt. 293.0) ncls = 2
        if(xx .ge. 293.0) ncls = 3

      end if
    
C ----------------------------------------------------------
c     get coefficients for this view angle
C ----------------------------------------------------------
      do i = 1,ntwx

C original:
c        ii = (lang - 1) * nx1 + i
C Add BT Classes (12 nov 01)
c        ii = (ncls-1)*numang*nx1 + (lang - 1) * nx1 + i
c Added land/ocean partitioning
         ii = land_flag*numcls_ocean*numang*nx1 + 
     &         (ncls-1)*numang*nx1 + (lang - 1) * nx1 + i
        do j = 1, ntwy
          coef(i, j) = coeftw(ii, j)
        end do
      end do

C Add BT Classes (12 nov 01)
c      ii = lang * nx1
      ii = land_flag*numcls_ocean*numang*nx1 +
     &      (ncls-1)*numang*nx1 + lang * nx1

      do j = 1, ntwy
        ytwb(j) = coeftw(ii, j)
      end do

c ... linear term of regression equation
      ixx=0

      do l = 1, nall
        if (kuse(l) . ne . 0) then
          ixx = ixx + 1
          xtw(ixx) = brt(l)
        end if
      end do

c ... Non-linear term of regression equation
      do l = 1, nall
        if (kuse(l) . ne . 0) then
          ixx = ixx + 1
          xtw(ixx) = brt(l)**2 / 250.0
        end if
      end do

c ... Additional predictors

c ... Surface pressure
      ixx = ixx + 1
      xtw(ixx) = psurf

c ... Percentage of land
      ixx = ixx + 1
      xtw(ixx) = pland

c ... Latitude
      ixx = ixx + 1
      xtw(ixx) = alat

c ... Month
      ixx = ixx + 1
      xtw(ixx) = amonth

c ... Evaluate regression equation
      do i = 1, ntwy

        sum = 0.0
        do j = 1, ntwx
          sum = sum + xtw(j) * coef(j, i)
        end do
        rtv(i) = sum + ytwb(i)

c ...   Water vapor mixing ratios and ozone partial pressures
c ...   are retrieved as logarithms, so convert to actual values
        if (i .le. 3 * nl .and. i .gt. nl) rtv(i) = exp(rtv(i))

c ...   22 jan 02: also take the exponent of retrieved tpw
c ...   tpw because we saved this as a logarithm in the coefficients
        if (i .eq. nl*3+2) rtv(i) = exp(rtv(i))

c ...   Limit minimum retrieved value for any parameter to 0.001
        if (rtv(i) .le. 0.001) rtv(i) = 0.001

      end do

c ... Contents of array RTV (319 elements) are as follows:
c ...   Elements   1-101 are atmospheric temperatures (K)
c ...   Elements 102-202 are water vapor mixing ratios (g/kg)
c ...   Elements 203-303 are ozone partial pressures (ppmv)
c ...   Element  304     is surface skin temperature (K)
c ...   Element  305    IS NOW integrated TPW (cm)

C for 10 emissivities
C wavelengths [3.7 4.3 5.0 5.8 7.6 8.3 9.3 10.8 12.1 14.3]
c ...   Elements 306-315 are IR surface emissivities
c ...   Elements 316-325 are IR surface reflectivities

c ... Extract retrieved profiles
      do i = 1, nl
        pres(i) = profile(i)
        temp(i) = rtv(i)
        wvmr(i) = rtv(i + nl)
        ozon(i) = rtv(i + (nl * 2))
      end do

C swemis (element 307, 4.3 microns, nl*3+4)
      swemis = rtv(nl*3 + 4)

C 17jan02: save directly retrieved tpw (element 305, nl*3 + 2)
      tpw_ret = rtv(nl*3 + 2)

C 1oct02: save skin temperature (element 304, nl*3 + 1)
      skint = rtv(nl*3 + 1)

c ... Do not allow retrieved water vapor mixing ratios greater than
c ... saturation mixing ratio
      do i = 1, nl
        wmax = satmix(pres(i), temp(i))
        wvmr(i) = min(wvmr(i), wmax)
C evab added on 08.26.2010 to fix dew point temperature and mixing ratio profile
c                          above 50 hPa
        wvmr(i) = max(wvmr(i), 0.001)
      end do

c-----------------------------------------------------------------------
c     DERIVED PRODUCTS
c-----------------------------------------------------------------------

c ... Find surface level in pressure profile
      lsfc = 1
      do i = 1, nl - 1
        if (psurf .ge. (pres(i) + 5.0)) lsfc = i + 1
      end do
      
c ... Compute total precipitable water (centimeters)
c      call tprecw_new2(pres, wvmr, tpw, nl, lsfc, psurf)

c ...  Nov 2009 EvaB: new algorithm is used for tpw_lo and tpw_hi too
c                     top level of integration is 21 (9.5 hPa)
      call tprecw(pres, wvmr, tpw, nl,21, lsfc, psurf)
      
c ... Compute total ozone (Dobson units)
      call total_ozone(pres, ozon, nl, lsfc, toz)

      END
