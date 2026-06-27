      REAL FUNCTION MODIS_BRIGHT(platform_name,RAD, BAND, UNITS)

C-----------------------------------------------------------------------
C!F77
C
C!DESCRIPTION:
C    Compute brightness temperature for a MODIS infrared band
C    on Terra or Aqua.
C
C    Spectral responses for each IR detector were obtained from MCST:
C    ftp://mcstftp.gsfc.nasa.gov/pub/permanent/MCST in
C    directories PFM_L1B_LUT_4-30-99 and FM1_RSR_LUT_07-10-01
C
C    An average spectral response for each infrared band was
C    computed. The band-averaged spectral response data were used
C    to compute the effective central wavenumbers and temperature
C    correction coefficients included in this module.
C
C    NOTE: The plaform name ('Terra' or 'Aqua') is passed to this
C    function via the common block defined in 'platform_name.inc'.
C
C!INPUT PARAMETERS:
C    RAD (REAL)      Planck radiance (units are determined by UNITS)
C    BAND (LONG)     MODIS IR band number (20-25, 27-36)
C    UNITS (LONG)    Flag defining radiance units
C                    0 => milliWatts per square meter per
C                         steradian per inverse centimeter
C                    1 => Watts per square meter per
C                         steradian per micron
C
C!OUTPUT PARAMETERS:
C    MODIS_BRIGHT  Brightness temperature (Kelvin)
C                  Note that a value of -1.0 is returned if
C                  RAD .LE. 0.0, or BAND is not in range 20-25, 27-36.
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
      real rad
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
	  real cwn_emas(8), tcs_emas(8), tci_emas(8)
	  		
      real cwn, tcs, tci
      integer index

c ... External functions
      real bright_m, brite_m
      external bright_m, brite_m
            
c ... Data statements

c-----------------------------------------------------------------------
c     BAND-AVERAGED MODIS SPECTRAL RESPONSE FUNCTIONS FOR TERRA
c     (LG 08-JAN-2002)
c     TEMPERATURE RANGE FOR FIT WAS  180.00 K TO  320.00 K
c
c     BANDS
c      20,  21,  22,  23,
c      24,  25,  27,  28,
c      29,  30,  31,  32,
c      33,  34,  35,  36,

c ... Effective central wavenumbers (inverse centimeters)
      data cwn_terra/
     & 2.641775E+03, 2.505277E+03, 2.518028E+03, 2.465428E+03,
     & 2.235815E+03, 2.200346E+03, 1.477967E+03, 1.362737E+03,
     & 1.173190E+03, 1.027715E+03, 9.080884E+02, 8.315399E+02,
     & 7.483394E+02, 7.308963E+02, 7.188681E+02, 7.045367E+02/

c ... Temperature correction slopes (no units)
      data tcs_terra/ 
     &  9.993411E-01, 9.998646E-01, 9.998585E-01, 9.998682E-01,
     &  9.998820E-01, 9.998845E-01, 9.994878E-01, 9.994918E-01,
     &  9.995496E-01, 9.997399E-01, 9.995607E-01, 9.997256E-01,
     &  9.999161E-01, 9.999167E-01, 9.999192E-01, 9.999282E-01/

c ... Temperature correction intercepts (Kelvin)
      data tci_terra/
     &  4.770522E-01, 9.264053E-02, 9.756834E-02, 8.928794E-02,
     &  7.309468E-02, 7.061112E-02, 2.204659E-01, 2.045902E-01,
     &  1.599076E-01, 8.249532E-02, 1.302885E-01, 7.181662E-02,
     &  1.970605E-02, 1.912743E-02, 1.816222E-02, 1.579983E-02/

c-----------------------------------------------------------------------
c     BAND-AVERAGED MODIS SPECTRAL RESPONSE FUNCTIONS FOR AQUA
c     (LG 08-JAN-2002)
c     TEMPERATURE RANGE FOR FIT WAS  180.00 K TO  320.00 K
c
c     BANDS
c      20,  21,  22,  23,
c      24,  25,  27,  28,
c      29,  30,  31,  32,
c      33,  34,  35,  36,

c ... Effective central wavenumbers (inverse centimeters)
      data cwn_aqua/
     &  2.647409E+03, 2.511760E+03, 2.517908E+03, 2.462442E+03,
     &  2.248296E+03, 2.209547E+03, 1.474262E+03, 1.361626E+03,
     &  1.169626E+03, 1.028740E+03, 9.076813E+02, 8.308411E+02,
     &  7.482978E+02, 7.307766E+02, 7.182094E+02, 7.035007E+02/

c ... Temperature correction slopes (no units)
      data tcs_aqua/ 
     &  9.993363E-01, 9.998626E-01, 9.998627E-01, 9.998707E-01,
     &  9.998737E-01, 9.998770E-01, 9.995694E-01, 9.994867E-01,
     &  9.995270E-01, 9.997382E-01, 9.995270E-01, 9.997271E-01,
     &  9.999173E-01, 9.999070E-01, 9.999198E-01, 9.999233E-01/

c ... Temperature correction intercepts (Kelvin)
      data tci_aqua/
     &  4.818401E-01, 9.426663E-02, 9.458604E-02, 8.736613E-02,
     &  7.873285E-02, 7.550804E-02, 1.848769E-01, 2.064384E-01,
     &  1.674982E-01, 8.304364E-02, 1.343433E-01, 7.135051E-02,
     &  1.948513E-02, 2.131043E-02, 1.804156E-02, 1.683156E-02/
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

c-----------------------------------------------------------------------
c     EMAS BANDS
c       3.7, 7.3, 8.5, 11, 12, 13.3, 13.6 and 13.9
c        26,  28,  30, 33, 34,   36,   37,      38
c ... Effective central wavenumbers (inverse centimeters)
      data cwn_emas/
     &  2683.8433, 1372.4951, 1167.6786, 903.9956, 832.9168, 
     &  749.9062, 733.5681, 716.1784 /

c ... Temperature correction slopes (no units)
      data tcs_emas/ 
     &  0.99932, 0.99925, 0.99977, 0.99985, 0.99974, 0.99985,
     &  0.99993, 0.99990/

c ... Temperature correction intercepts (Kelvin)
      data tci_emas/
     &  0.50705, 0.32168, 0.082686, 0.041587, 0.069156, 
     &  0.036514, 0.017001, 0.022627 /



c ... Set default return value
      modis_bright = -1.0


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

		else if (platform_name(1:11) .eq. 'NPP_:_VIIRS' .or. 
     &			platform_name(1:3) .eq. 'npp' ) then
			index = band - 11
	
		else if (platform_name(1:3) .eq. 'AMS') then 
			if (band .eq. 11) index = 1
			if (band .eq. 12) index = 2

		else if (platform_name(1:4) .eq. 'EMAS') then 
		
			if (band .eq. 26) index = 1
			if (band .eq. 28) index = 2
			if (band .eq. 30) index = 3
			if (band .eq. 33) index = 4
			if (band .eq. 34) index = 5
			if (band >= 36) index = band-30
			
c MODIS
		else 

c ... Check input parameters and return if they are bad
      if (rad .le. 0.0 .or.
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
      else if (platform_name(1:11) .eq. 'NPP_:_VIIRS' .or. 
     & 		   platform_name(1:3) .eq. 'npp' ) then
        cwn = cwn_viirs(index)
        tcs = tcs_viirs(index)
        tci = tci_viirs(index)
      else if (platform_name(1:3) .eq. 'AMS') then 
      	cwn = cwn_ams(index)
      	tcs = tcs_ams(index)
      	tci = tci_ams(index)
      else if (platform_name(1:4) .eq. 'EMAS') then 
      	cwn = cwn_emas(index)
      	tcs = tcs_emas(index)
      	tci = tci_emas(index)
      else
      endif
     
c ... Compute brightness temperature
      if (units .eq. 1) then

c ...   Radiance units are
c ...   Watts per square meter per steradian per micron
        modis_bright = (bright_m(1.0e+4 / cwn, rad) - tci) / tcs

      else

c ...   Radiance units are
c ...   milliWatts per square meter per steradian per wavenumber
        modis_bright = (brite_m(cwn, rad) - tci) / tcs

      endif

      END
