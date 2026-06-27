      REAL FUNCTION MODIS_PLANCK(platform_name, TEMP, BAND, UNITS)

C-----------------------------------------------------------------------
C!F77
C
C!DESCRIPTION:
C    Compute Planck radiance for a MODIS infrared band on Terra or Aqua.
C
C    Spectral responses for each IR detector were obtained from MCST:
C    ftp://ftp.mcst.ssai.biz/pub/permanent/MCST in directories
C    PFM_L1B_LUT_4-30-99 (Terra) and FM1_RSR_LUT_07-10-01 (Aqua).
C
C    An average spectral response for all detectors in each band was
C    computed. The detector-averaged spectral response data were used
C    to compute the effective central wavenumbers and temperature
C    correction coefficients included in this module.
C
C    NOTE: The plaform name ('Terra' or 'Aqua') is passed to this
C    function via the common block defined in 'platform_name.inc'.
C
C!INPUT PARAMETERS:
C    TEMP (REAL)     Brightness temperature (Kelvin)
C    BAND (LONG)     MODIS IR band number (20-25, 27-36)
C    UNITS (LONG)    Flag defining radiance units
C                    0 => milliWatts per square meter per
C                         steradian per inverse centimeter
C                    1 => Watts per square meter per
C                         steradian per micron
C
C!OUTPUT PARAMETERS:
C    MODIS_PLANCK  Planck radiance (units are determined by UNITS)
C                  Note that a value of -1.0 is returned if
C                  TEMP .LE. 0.0, or BAND is not in range 20-25, 27-36.
C
C!REVISION HISTORY:
C    Liam.Gumley@ssec.wisc.edu
C
C!TEAM-UNIQUE HEADER:
C    Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
C
C!END
C-----------------------------------------------------------------------
      
      IMPLICIT NONE

c ... Arguments
      real temp
      integer band, units
      character*(*) platform_name

c ... Local variables
      real cwn_terra(16), tcs_terra(16), tci_terra(16)
      real cwn_aqua(16), tcs_aqua(16), tci_aqua(16)
      real cwn_mas(6), tcs_mas(6), tci_mas(6)
      real cwn_master(6), tcs_master(6), tci_master(6)
      real cwn_seviri(6), tcs_seviri(6), tci_seviri(6)
      real cwn_viirs(5), tcs_viirs(5), tci_viirs(5)
      real cwn_ams(2), tcs_ams(2), tci_ams(2)

      real cwn, tcs, tci
      integer index
      
c ... External functions
      real planck_m, planc_m
      external planck_m, planc_m
            
c ... Data statements

c-----------------------------------------------------------------------

c     TERRA MODIS DETECTOR-AVERAGED SPECTRAL RESPONSE
c     (LIAM GUMLEY 2003/06/05)

c     BAND 20 TEMPERATURE RANGE WAS  180.00 K TO  350.00 K
c     BAND 21 TEMPERATURE RANGE WAS  180.00 K TO  400.00 K
c     BAND 22 TEMPERATURE RANGE WAS  180.00 K TO  350.00 K
c     BAND 23 TEMPERATURE RANGE WAS  180.00 K TO  350.00 K
c     BAND 24 TEMPERATURE RANGE WAS  180.00 K TO  320.00 K
c     BAND 25 TEMPERATURE RANGE WAS  180.00 K TO  320.00 K
c     BAND 27 TEMPERATURE RANGE WAS  180.00 K TO  320.00 K
c     BAND 28 TEMPERATURE RANGE WAS  180.00 K TO  320.00 K
c     BAND 29 TEMPERATURE RANGE WAS  180.00 K TO  340.00 K
c     BAND 30 TEMPERATURE RANGE WAS  180.00 K TO  340.00 K
c     BAND 31 TEMPERATURE RANGE WAS  180.00 K TO  340.00 K
c     BAND 32 TEMPERATURE RANGE WAS  180.00 K TO  340.00 K
c     BAND 33 TEMPERATURE RANGE WAS  180.00 K TO  330.00 K
c     BAND 34 TEMPERATURE RANGE WAS  180.00 K TO  320.00 K
c     BAND 35 TEMPERATURE RANGE WAS  180.00 K TO  310.00 K
c     BAND 36 TEMPERATURE RANGE WAS  180.00 K TO  310.00 K

c     BANDS
c      20,  21,  22,  23,
c      24,  25,  27,  28,
c      29,  30,  31,  32,
c      33,  34,  35,  36,

c ... Effective central wavenumbers (inverse centimeters)
      data cwn_terra/
     &  2.641767E+03, 2.505274E+03, 2.518031E+03, 2.465422E+03,
     &  2.235812E+03, 2.200345E+03, 1.478026E+03, 1.362741E+03,
     &  1.173198E+03, 1.027703E+03, 9.081998E+02, 8.315149E+02,
     &  7.483224E+02, 7.309089E+02, 7.188677E+02, 7.045309E+02/

c ... Temperature correction slopes (no units)
      data tcs_terra/ 
     &  9.993487E-01, 9.998699E-01, 9.998604E-01, 9.998701E-01,
     &  9.998825E-01, 9.998849E-01, 9.994942E-01, 9.994937E-01,
     &  9.995643E-01, 9.997499E-01, 9.995880E-01, 9.997388E-01,
     &  9.999192E-01, 9.999171E-01, 9.999174E-01, 9.999264E-01/

c ... Temperature correction intercepts (Kelvin)
      data tci_terra/
     &  4.744530E-01, 9.091094E-02, 9.694298E-02, 8.856134E-02,
     &  7.287017E-02, 7.037161E-02, 2.177889E-01, 2.037728E-01,
     &  1.559624E-01, 7.989879E-02, 1.176660E-01, 6.856633E-02,
     &  1.903625E-02, 1.902709E-02, 1.859296E-02, 1.619453E-02/

c-----------------------------------------------------------------------

c     AQUA MODIS DETECTOR-AVERAGED SPECTRAL RESPONSE
c     (LIAM GUMLEY 2003/06/05)

c     BAND 20 TEMPERATURE RANGE WAS  180.00 K TO  350.00 K
c     BAND 21 TEMPERATURE RANGE WAS  180.00 K TO  400.00 K
c     BAND 22 TEMPERATURE RANGE WAS  180.00 K TO  350.00 K
c     BAND 23 TEMPERATURE RANGE WAS  180.00 K TO  350.00 K
c     BAND 24 TEMPERATURE RANGE WAS  180.00 K TO  320.00 K
c     BAND 25 TEMPERATURE RANGE WAS  180.00 K TO  320.00 K
c     BAND 27 TEMPERATURE RANGE WAS  180.00 K TO  320.00 K
c     BAND 28 TEMPERATURE RANGE WAS  180.00 K TO  320.00 K
c     BAND 29 TEMPERATURE RANGE WAS  180.00 K TO  340.00 K
c     BAND 30 TEMPERATURE RANGE WAS  180.00 K TO  340.00 K
c     BAND 31 TEMPERATURE RANGE WAS  180.00 K TO  340.00 K
c     BAND 32 TEMPERATURE RANGE WAS  180.00 K TO  340.00 K
c     BAND 33 TEMPERATURE RANGE WAS  180.00 K TO  330.00 K
c     BAND 34 TEMPERATURE RANGE WAS  180.00 K TO  320.00 K
c     BAND 35 TEMPERATURE RANGE WAS  180.00 K TO  310.00 K
c     BAND 36 TEMPERATURE RANGE WAS  180.00 K TO  310.00 K

c     BANDS
c       20,  21,  22,  23,
c       24,  25,  27,  28,
c       29,  30,  31,  32,
c       33,  34,  35,  36,

c ... Effective central wavenumbers (inverse centimeters)
      data cwn_aqua/
     &  2.647418E+03, 2.511763E+03, 2.517910E+03, 2.462446E+03,
     &  2.248296E+03, 2.209550E+03, 1.474292E+03, 1.361638E+03,
     &  1.169637E+03, 1.028715E+03, 9.076808E+02, 8.308397E+02,
     &  7.482977E+02, 7.307761E+02, 7.182089E+02, 7.035020E+02/

c ... Temperature correction slopes (no units)
      data tcs_aqua/ 
     &  9.993438E-01, 9.998680E-01, 9.998649E-01, 9.998729E-01,
     &  9.998738E-01, 9.998774E-01, 9.995732E-01, 9.994894E-01,
     &  9.995439E-01, 9.997496E-01, 9.995483E-01, 9.997404E-01,
     &  9.999194E-01, 9.999071E-01, 9.999176E-01, 9.999211E-01/

c ... Temperature correction intercepts (Kelvin)
      data tci_aqua/
     &  4.792821E-01, 9.260598E-02, 9.387793E-02, 8.659482E-02,
     &  7.854801E-02, 7.521532E-02, 1.833035E-01, 2.053504E-01,
     &  1.628724E-01, 8.003410E-02, 1.290129E-01, 6.810679E-02,
     &  1.895925E-02, 2.128960E-02, 1.857071E-02, 1.733782E-02/

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     MAS BANDS
c       30, 42, 45, 46, 48, 49

c ... Effective central wavenumbers (inverse centimeters)
      data cwn_mas/
     &  2652.519, 1170.411, 912.99, 830.77, 751.16, 723.96/

c ... Temperature correction slopes (no units)
      data tcs_mas/ 
     &  0.99937, 0.99937, 0.99955, 0.99967, 0.99978, 0.99977/

c ... Temperature correction intercepts (Kelvin)
      data tci_mas/
     &  0.44728, 0.22513, 0.13007, 0.087251, 0.05295, 0.053137/
c-----------------------------------------------------------------------
c     MASTER BANDS
c       30, 43, 47, 49 plus SHIS MAS-equivalent bands 48 and 49

c ... Effective central wavenumbers (inverse centimeters)
      data cwn_master/
     &  2670.227, 1155.4015, 941.17, 820.51, 751.16, 723.96 /

c ... Temperature correction slopes (no units)
      data tcs_master/ 
     &  0.99944, 0.99934, 0.99937, 0.99963, 0.99978, 0.99977/

c ... Temperature correction intercepts (Kelvin)
      data tci_master/
     &  0.41362, 0.23203, 0.18665, 0.094882, 0.05295, 0.053137/
c-----------------------------------------------------------------------
c     SEVIRI BANDS
c       3.9, 7.3, 8.5, 11., 12. and 13.4 um

c ... Effective central wavenumbers (inverse centimeters)
      data cwn_seviri/
     &  2569.094, 1362.142, 1149.083, 930.659, 839.661, 746.27 /

c ... Temperature correction slopes (no units)
      data tcs_seviri/ 
     &  0.9959, 0.9991, 0.9996, 0.9983, 0.9988, 0.9981/

c ... Temperature correction intercepts (Kelvin)
      data tci_seviri/
     &  3.471, 0.485, 0.181, 0.627, 0.397, 0.576/
	 
	 
c-----------------------------------------------------------------------
c     VIIRS BANDS
c       3., 7.3, 8.5, 11., 12. and 13.4 um

c ... Effective central wavenumbers (inverse centimeters)
      data cwn_viirs/
     &  2708.3865, 2460.8193, 1166.1845, 935.10476, 845.79752  /

c ... Temperature correction slopes (no units)
      data tcs_viirs/ 
     &  0.999354, 0.999623,  0.999698, 0.998273, 0.998778/

c ... Temperature correction intercepts (Kelvin)
      data tci_viirs/
     &  0.593537, 0.325879, 0.146942, 0.650338,  0.421701 /
	 
c-----------------------------------------------------------------------
c     AMS BANDS
c       3.7 and 10.2

c ... Effective central wavenumbers (inverse centimeters)
      data cwn_ams/
     &  2656.04241, 979.9118 /

c ... Temperature correction slopes (no units)
      data tcs_ams/ 
     &  0.99942, 0.99778/

c ... Temperature correction intercepts (Kelvin)
      data tci_ams/
     &  0.42667, 0.65987 /



c ... Set default return value
      modis_planck = -1.0


		if (platform_name(1:6) .eq. 'Master' .or.
     &         platform_name(1:6) .eq. 'master' .or.
     &         platform_name(1:6) .eq. 'MASTER') then

		if (band .eq. 30) index = 1
		if (band .eq. 43) index = 2
		if (band .eq. 47) index = 3
		if (band .eq. 49) index = 4
		if (band .eq. 60) index = 5
		if (band .eq. 61) index = 6

		else if (platform_name(1:3) .eq. 'Mas' .or.
     &         platform_name(1:3) .eq. 'mas' .or.
     &         platform_name(1:3) .eq. 'MAS') then

		if (band .eq. 30) index = 1
		if (band .eq. 42) index = 2
		if (band .eq. 45) index = 3
		if (band .eq. 46) index = 4
		if (band .eq. 48 .or. band .eq. 60) index = 5
		if (band .eq. 49 .or. band .eq. 61) index = 6

		else if (platform_name(1:6) .eq. 'SEVIRI' .or.
     &         platform_name(1:6) .eq. 'seviri' .or.
     &         platform_name(1:6) .eq. 'Seviri') then

		if (band .eq. 4) index = 1
		if (band .eq. 6) index = 2
		if (band .eq. 7) index = 3
		if (band .eq. 9) index = 4
		if (band .eq. 10) index = 5
		if (band .eq. 11) index = 6

		else if (platform_name(1:11) .eq. 'NPP_:_VIIRS') then

			index = band - 11
		else if (platform_name(1:3) .eq. 'AMS') then 
			if (band .eq. 11) index = 1
			if (band .eq. 12) index = 2

		else

c ... Check input parameters and return if they are bad
      if (temp .le. 0.0 .or.
     &    band .lt. 20 .or.
     &    band .gt. 36 .or.
     &    band .eq. 26) return

c ... Get index into coefficient arrays
      if (band .le. 25) then
        index = band - 19
      else
        index = band - 20
      endif

		endif

      
c ... Get the coefficients for Terra or Aqua
      if (platform_name(1:5) .eq. 'Terra' .or.
     &    platform_name(1:5) .eq. 'terra' .or.
     &    platform_name(1:5) .eq. 'TERRA') then
        cwn = cwn_terra(index)
        tcs = tcs_terra(index)
        tci = tci_terra(index)
      else if (platform_name(1:4) .eq. 'Aqua' .or.
     &         platform_name(1:4) .eq. 'aqua' .or.
     &         platform_name(1:4) .eq. 'AQUA') then
        cwn = cwn_aqua(index)
        tcs = tcs_aqua(index)
        tci = tci_aqua(index)
      else if (platform_name(1:3) .eq. 'Mas' .or.
     &         platform_name(1:3) .eq. 'mas' .or.
     &         platform_name(1:3) .eq. 'MAS') then
        cwn = cwn_mas(index)
        tcs = tcs_mas(index)
        tci = tci_mas(index)
      else if (platform_name(1:6) .eq. 'Master' .or.
     &         platform_name(1:6) .eq. 'master' .or.
     &         platform_name(1:6) .eq. 'MASTER') then
        cwn = cwn_master(index)
        tcs = tcs_master(index)
        tci = tci_master(index)
      else if (platform_name(1:6) .eq. 'Seviri' .or.
     &         platform_name(1:6) .eq. 'seviri' .or.
     &         platform_name(1:6) .eq. 'SEVIRI') then
        cwn = cwn_seviri(index)
        tcs = tcs_seviri(index)
        tci = tci_seviri(index)
      else if (platform_name(1:11) .eq. 'NPP_:_VIIRS') then
        cwn = cwn_viirs(index)
        tcs = tcs_viirs(index)
        tci = tci_viirs(index)
      else if (platform_name(1:3) .eq. 'AMS') then 
      	cwn = cwn_ams(index)
      	tcs = tcs_ams(index)
      	tci = tci_ams(index)
      else
      endif
                
c ... Compute Planck radiance
      if (units .eq. 1) then

c ...   Radiance units are
c ...   Watts per square meter per steradian per micron
        modis_planck = planck_m(1.0e+4 / cwn, temp * tcs + tci)

      else

c ...   Radiance units are
c ...   milliWatts per square meter per steradian per wavenumber
        modis_planck = planc_m(cwn, temp * tcs + tci)

      endif

      END
