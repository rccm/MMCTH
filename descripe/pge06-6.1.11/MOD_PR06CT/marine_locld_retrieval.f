      subroutine marine_locld_retrieval(rlat, met_month, isp,
     *           ltrp, ttpp, z, tp, tw, lwin, lmlc, lapse_rate)
 
c!F77-------------------------------------------------------------------
c!Description: Calculate marine low cloud height using a lapse rate
c              method.
c
c!Input parameters:
c debug               Debug write flag
c rlat                Latitude
c met_month           Month                     
c isp                 Level of surface pressure (1-101)
c ltrp                Level of tropopause (1-101)
c ttpp                Brightness temperatures for each profile level
c z                   Geopotential height profile
c tp                  Atm. temperature profile
c tw                  Observed 11 micron brightness temperature         
c lwin                IRW cloud top pressure level
c
c!Output parameters:
c lmlc                Level of cloud top pressure retrieval (1-101)
c
c!Revision history:
c                     12/10       R. Frey 
c                     Modified for MODIS from HIRS CO2-slicing code
c
c!Calls:
c None
c
c!End-------------------------------------------------------------------

      implicit none
      save

      include 'mod06uw_debug.inc'

      integer plev
      parameter (plev = 101)

c     Scalar arguments.
      integer met_month, isp, ltrp, lwin, lmlc
      real    rlat, tw

c     Array arguments.
      real    ttpp(plev), z(plev), tp(plev)

c     Local scalars.
      double precision dplat, c0, c1, c2, c3, c4
      real             lapse_rate, dt, zcld
      integer          jm, llwin, k, zpl, zph

c     Local arrays.
      double precision month_coeffs(5,12,3), breakpts(2,12)
      double precision a0(12), a1(12), a2(12)
     
c-----------------------------------------------------------------------

c     Lapse rate coefficients, by month of year.
      data a0 /3.9345, 3.7507, 3.7440, 3.7810, 3.7956, 3.9949,
     *         4.0161, 4.1329, 3.9528, 4.0972, 3.9939, 4.1510/
      data a1 / 0.0104,  0.0099,  0.0066, -0.0004, -0.0093, -0.0144,
     *         -0.0185, -0.0225, -0.0189, -0.0081,  0.0003,  0.0080/
      data a2 /0.0008, 0.0008, 0.0009, 0.0009, 0.0008, 0.0007,
     *         0.0006, 0.0005, 0.0007, 0.0008, 0.0008, 0.0008/
      
c     Lapse rate coefficients, by month of year.

      data month_coeffs /
     *  2.9769800876, -0.0515871084,  0.0027409105,  0.0001135740,  0.00000113040,
     *  3.3483238557,  0.1372575458,  0.0133258850,  0.0003043608,  0.00000218650,
     *  2.4060295675,  0.0372001609,  0.0096472724,  0.0002334206,  0.00000165450,
     *  2.6522386726,  0.0325728824,  0.0100892620,  0.0002601226,  0.00000198560,
     *  1.9578262599, -0.2112028966, -0.0057943564, -0.0001050464, -0.00000074313,
     *  2.7659753980, -0.1186500984,  0.0011626989,  0.0000936998,  0.00000101060,
     *  2.1106811602, -0.3073665907, -0.0090862456, -0.0000889596,  0.00000003552,
     *  3.0982173723, -0.1629588435, -0.0020384299,  0.0000286274,  0.00000060283,
     *  3.0760551826, -0.2043463270, -0.0053969994, -0.0000541329, -0.00000001768,
     *  3.6377215316, -0.0857783614,  0.0024313179,  0.0001495010,  0.00000170850,
     *  3.3206165420, -0.1411094234, -0.0026068389,  0.0000057937,  0.00000042113,
     *  3.0526632533, -0.1121521836, -0.0009912556,  0.0000179690,  0.00000027070,

     *  2.9426577089, -0.0510674066,  0.0052420293,  0.0001097927, -0.00000372380,
     *  2.6499605646, -0.0105152229,  0.0042895903,  0.0000719741, -0.00000066735,
     *  2.3652046763,  0.0141129341,  0.0059242144, -0.0000158816, -0.00000265790,
     *  2.5433158163, -0.0046876415,  0.0059325140,  0.0000143938, -0.00000346360,
     *  2.4994027830, -0.0364706332,  0.0082001522,  0.0000843577, -0.00000768780,
     *  2.7641495752, -0.0728625243,  0.0088877822,  0.0001767765, -0.00001168390,
     *  3.1202042743, -0.1002374752,  0.0064054059,  0.0002620230, -0.00001078950,
     *  3.4331195144, -0.1021765880,  0.0010498850,  0.0001614861,  0.00000510150,
     *  3.4539389485, -0.1158261776,  0.0015449592,  0.0001711651,  0.00000248080,
     *  3.6013336912, -0.0775800028,  0.0041940388,  0.0000941307, -0.00000408720,
     *  3.1947419143, -0.1045316345,  0.0049986486,  0.0001910731, -0.00000505500,
     *  3.1276377012, -0.0707628268,  0.0055532926,  0.0001549833, -0.00000570980,

     *  1.9009562748,  0.0236905223,  0.0086504022, -0.0002167013,  0.00000151230,
     *  2.4878735828, -0.0076514110,  0.0079443995, -0.0001773726,  0.00000114730,
     *  3.1251275103, -0.1214572133,  0.0146488407, -0.0003187508,  0.00000210290,
     * 13.3931706579, -1.2206947755,  0.0560380539, -0.0009873591,  0.00000598210,
     *  1.6432070460,  0.1151206937,  0.0033130967, -0.0001458434,  0.00000128610,
     * -5.2366360253,  1.0105574562, -0.0355440449,  0.0005187964, -0.00000262410,
     * -4.7396480830,  0.9625734101, -0.0355846807,  0.0005522497, -0.00000299860,
     * -1.4424842734,  0.4769307320, -0.0139027010,  0.0001758823, -0.00000079846,
     * -3.7140186247,  0.6720953861, -0.0210548327,  0.0002974491, -0.00000149380,
     *  8.2237401369, -0.5127532741,  0.0205285436, -0.0003015662,  0.00000157680,
     * -0.4502046794,  0.2629679617, -0.0018419395, -0.0000368887,  0.00000048223,
     *  9.3930897423, -0.8836682302,  0.0460453172, -0.0008450362,  0.00000517810/

      data breakpts /
     * - 3.8,  22.1,
     * -21.5,  12.8,
     * - 2.8,  10.7,
     * -23.4,  29.4,
     * -12.3,  14.9,
     * - 7.0,  16.8,
     * -10.5,  15.0,
     * - 7.8,  19.5,
     * - 8.6,  17.4,
     * - 7.0,  27.0,
     * - 9.2,  22.0,
     * - 3.7,  19.0/

c-----------------------------------------------------------------------

c     Define month of year as an index,
      jm = met_month

c     Get monthly mean lapse rates.  There are three separate sets of
c     coefficients, delineated by 'breakpts' (latitudes).

c     Determine which coefficients to use based on latitude.
      if(rlat .lt. breakpts(1,jm)) then
        c0 = month_coeffs(1,jm,1)
        c1 = month_coeffs(2,jm,1)
        c2 = month_coeffs(3,jm,1)
        c3 = month_coeffs(4,jm,1)
        c4 = month_coeffs(5,jm,1)
      else if(rlat .ge. breakpts(1,jm) .and. rlat .le. breakpts(2,jm)) then
        c0 = month_coeffs(1,jm,2)
        c1 = month_coeffs(2,jm,2)
        c2 = month_coeffs(3,jm,2)
        c3 = month_coeffs(4,jm,2)
        c4 = month_coeffs(5,jm,2)
      else
        c0 = month_coeffs(1,jm,3)
        c1 = month_coeffs(2,jm,3)
        c2 = month_coeffs(3,jm,3)
        c3 = month_coeffs(4,jm,3)
        c4 = month_coeffs(5,jm,3)
      end if

      dplat = dble(rlat)

c     lapse_rate = sngl ( a0(met_month) + (a1(met_month) * dplat) +
c    *             (a2(met_month) * dplat**2.d0) )

      lapse_rate = c0 + c1*rlat + c2*rlat**2 + c3*rlat**3 + c4*rlat**4
      if(lapse_rate .lt. 2.0) lapse_rate = 2.0
      if(lapse_rate .gt. 10.0) lapse_rate = 10.0

c     Get temperature difference between surface (calculated) and cloud
c     (observed).
      dt = ttpp(isp) - tw

c     Calculate cloud height.
      if(dt .gt. 0.0) then
        zcld = dt / lapse_rate
      else
        zcld = 0.0
      end if

c     Find closest z-level.
      if(zcld .eq. 0.0) then

        llwin = isp

      else

c       Initialize. Sometimes 'zcld' < 'z(isp)'.
        zpl = isp
        zph = isp - 1

        do k = isp, ltrp, -1
          if (zcld .ge. z(k)) then
            zpl = k
            zph = k - 1
          end if
        enddo

        if( (abs(zcld-z(zph))) .le. (abs(zcld-z(zpl))) ) then
          llwin = zph
        else
          llwin = zpl
        end if

      end if

      lmlc = llwin

c-----------------------------------------------------------------------

      if(debug .gt. 0) then
        write(h_output,'(/,''Ocean low cloud algorithm data:'')')
        write(h_output,'(''Month, lat., coeffs.: '', i5, f10.2, 5f15.10)')
     *       jm, rlat, c0, c1, c2, c3, c4
        write(h_output,'(''Calc. sfc. BT, Obs. BT, dt, lapse rate, '',
     *            ''cld hgt, z(sfc)'')')
        write(h_output,'(7f10.3)') ttpp(isp),tw,dt,lapse_rate,zcld,
     *                             z(isp)  
        write(h_output,'(1x)')
        write(h_output,'(''Interpolation:'')')
        write(h_output,'(''trop and sfc levels, zph, zpl, cld lev'')')
        write(h_output,'(6i10)') ltrp,isp,zph,zpl,lmlc
        write(h_output,'(''zlo, zhi, z: '')')
        write(h_output,'(3f10.3)') z(zpl),z(zph),z(lmlc)
        write(h_output,'(1x)')
      end if

c-----------------------------------------------------------------------

      return
      end
