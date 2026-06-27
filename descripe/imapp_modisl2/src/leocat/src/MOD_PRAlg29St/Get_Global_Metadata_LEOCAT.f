      SUBROUTINE GET_GLOBAL_METADATA_LEOCAT( MODFIL, PCGOOD, 
     $ PCN1S, PCN2S, PCN3S, PCN4S, PCNNCD, PCNNCR, PCNND,
     $ PCNNT, PCNNG, PCNNI, PCNNL, PCNNW, PCNNV, PCNNR,
     $ PCNNC, MAX_SOL, MIN_SOL, NDVIFIL, THRESFIL)
C-----------------------------------------------------------------------
C !F77
C
C !PURPOSE:
C    Extract metadata from Collection 6 processed cloud mask granule.
C
C !DESCRIPTION:
C    Extracts Global metadata from the Collection 6 Processed granule which
C    is needed for UW product generation.
C    Please note that calling this subroutine will cause
C    a warning message to be inserted into the LogStatus
C    file.  Please see description of message MTYPEf2c().
C
C !INPUT PARAMETERS:
C    MODFIL       File information array containing
C                 HDF SD_ID (index 1)
C                 VS file identifier (index 2)
C                 File access mode (index 3)
C
C !OUTPUT PARAMETERS:
C    PCBAD        Percentage of good data
C    PCN1S        Very high confident clear percentage
C    PCN2S        High confident clear percentage
C    PCN3S        Uncertain confident clear percentage
C    PCN4S        Low confident clear percentage
C
C !REVISION HISTORY:
C
C !TEAM-UNIQUE HEADER:
C    Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
C
C !END
C-----------------------------------------------------------------------

      IMPLICIT NONE

      INCLUDE 'mapi.inc'
      INCLUDE 'PGS_SMF.f'

c ... Arguments

      integer modfil( MODFILLEN )
      real pcgood, pcn1s, pcn2s, pcn3s, pcn4s, pcnnd, pcnnt,
     &  pcnng, pcnni, pcnnl, pcnnw, pcnnv, pcnnr, pcnnc, pcnncd, pcnncr,
     &  max_sol, min_sol
      character*34 ndvifil
      character*23 thresfil

c ... Local variables

      character*72 dtype, attrib
      integer nelmnt, rtn
      real value, pcbad
      character*35 charval
      character*256 cline

c ... Get number of successful retrievals for this granule

      dtype = '  '
      attrib = 'SuccessfulRetrievalPct'
      nelmnt = 0
      value = 0.0
      rtn = 0

c ... The error message
c ... 'MTYPEf2c():MAPI_E_ERR:324431360 (PID=26923)
c ... ERROR: MTYPEf2c unable to use empty input parameter.'
c ... appears in LogStatus when empty inputs for variables 'dtype'
c ... and 'nelmnt' are passed to GMFIN().
c ... This is an unfortunate message since the M-API library
c ... was designed to accept empty inputs, and to return their actual
c ... values.

c ... Make initial call to get correct input parameters
c ... without hardwiring

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )
      
c ... Second call will replace variables with correct inputs

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )
      
      if (rtn .ne. 0) then
        call message( 'Get_Global_Metadata_Leocat',
     &    'Error getting Successful Retrieval Pct from C6 file' //
     &    ' [OPERATOR ACTION: Stage correct C6, rerun PGE]', 0, 2 )
      endif

c ... Set return value

      pcgood = value 

c ... Get number of very high percentage clear for this granule

      dtype = '  '
      attrib = 'VeryHighConfidentClearPct'
      nelmnt = 0
      value = 0.0
      rtn = 0

c ... Make initial call to get correct input parameters 
c ... without hardwiring

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )

c ... Second call will replace variables with correct inputs

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )
      
      if (rtn .ne. 0) then
        call message( 'Get_Global_Metadata_Leocat',
     &    'Error getting Very High Confident Clear from C6 file' //
     &    ' [OPERATOR ACTION: Stage correct C6, rerun PGE]', 0, 2 )
      end if

c ... Set return value

      pcn1s = value 

c ... Get number of high percentage clear for this granule

      dtype = '  '
      attrib = 'HighConfidentClearPct'
      nelmnt = 0
      value = 0.0
      rtn = 0

c ... Make initial call to get correct input parameters 
c ... without hardwiring

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )

c ... Second call will replace variables with correct inputs

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )
      
      if (rtn .ne. 0) then
        call message( 'Get_Global_Metadata_Leocat',
     &    'Error getting High Confident Clear from C6 file' //
     &    ' [OPERATOR ACTION: Stage correct C6, rerun PGE]', 0, 2 )
      end if

c ... Set return value

      pcn2s = value 

c ... Get number of uncertain percentage clear for this granule

      dtype = '  '
      attrib = 'UncertainConfidentClearPct'
      nelmnt = 0
      value = 0.0
      rtn = 0

c ... Make initial call to get correct input parameters 
c ... without hardwiring

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )

c ... Second call will replace variables with correct inputs

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )
      
      if (rtn .ne. 0) then
        call message( 'Get_Global_Metadata_Leocat',
     &    'Error getting Uncertain Confident Clear from C6 file' //
     &    ' [OPERATOR ACTION: Stage correct C6, rerun PGE]', 0, 2 )
      end if

c ... Set return value

      pcn3s = value 

c ... Get number of low percentage clear for this granule

      dtype = '  '
      attrib = 'LowConfidentClearPct'
      nelmnt = 0
      value = 0.0
      rtn = 0

c ... Make initial call to get correct input parameters 
c ... without hardwiring

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )

c ... Second call will replace variables with correct inputs

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )
      
      if (rtn .ne. 0) then
        call message( 'Get_Global_Metadata_Leocat',
     &    'Error getting Low Confident Clear from C6 file' //
     &    ' [OPERATOR ACTION: Stage correct C6, rerun PGE]', 0, 2 )
      end if

c ... Set return value

      pcn4s = value 

c ... Get number of cloud cover (250 m) for this granule

      dtype = '  '
      attrib = 'CloudCoverPct250m'
      nelmnt = 0
      value = 0.0
      rtn = 0

c ... Make initial call to get correct input parameters 
c ... without hardwiring

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )

c ... Second call will replace variables with correct inputs

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )
      
      if (rtn .ne. 0) then
        call message( 'Get_Global_Metadata_Leocat',
     &    'Error getting Cloud Cover 250 m from C6 file' //
     &    ' [OPERATOR ACTION: Stage correct C6, rerun PGE]', 0, 2 )
      end if

c ... Set return value

      pcnncd = value 

c ... Get number of clear (250 m) for this granule

      dtype = '  '
      attrib = 'ClearPct250m'
      nelmnt = 0
      value = 0.0
      rtn = 0

c ... Make initial call to get correct input parameters 
c ... without hardwiring

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )

c ... Second call will replace variables with correct inputs

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )
      
      if (rtn .ne. 0) then
        call message( 'Get_Global_Metadata_Leocat',
     &    'Error getting Clear 250 m from C6 file' //
     &    ' [OPERATOR ACTION: Stage correct C6, rerun PGE]', 0, 2 )
      end if

c ... Set return value

      pcnncr = value 

c ... Get number of day processed for this granule

      dtype = '  '
      attrib = 'DayProcessedPct'
      nelmnt = 0
      value = 0.0
      rtn = 0

c ... Make initial call to get correct input parameters 
c ... without hardwiring

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )

c ... Second call will replace variables with correct inputs

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )
      
      if (rtn .ne. 0) then
        call message( 'Get_Global_Metadata_Leocat',
     &    'Error getting Day Processed from C6 file' //
     &    ' [OPERATOR ACTION: Stage correct C6, rerun PGE]', 0, 2 )
      end if

c ... Set return value

      pcnnd = value 

c ... Get number of night processed for this granule

      dtype = '  '
      attrib = 'NightProcessedPct'
      nelmnt = 0
      value = 0.0
      rtn = 0

c ... Make initial call to get correct input parameters 
c ... without hardwiring

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )

c ... Second call will replace variables with correct inputs

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )
      
      if (rtn .ne. 0) then
        call message( 'Get_Global_Metadata_Leocat',
     &    'Error getting Night Processed from C6 file' //
     &    ' [OPERATOR ACTION: Stage correct C6, rerun PGE]', 0, 2 )
      end if

c ... Set return value

      pcnnt = value 

c ... Get number of sunglint processed for this granule

      dtype = '  '
      attrib = 'SunglintProcessedPct'
      nelmnt = 0
      value = 0.0
      rtn = 0

c ... Make initial call to get correct input parameters 
c ... without hardwiring

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )

c ... Second call will replace variables with correct inputs

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )
      
      if (rtn .ne. 0) then
        call message( 'Get_Global_Metadata_Leocat',
     &    'Error getting Sunglint Processed from C6 file' //
     &    ' [OPERATOR ACTION: Stage correct C6, rerun PGE]', 0, 2 )
      end if

c ... Set return value

      pcnng = value 

c ... Get number of snow/ice processed for this granule

      dtype = '  '
      attrib = 'Snow_IceSurfaceProcessedPct'
      nelmnt = 0
      value = 0.0
      rtn = 0

c ... Make initial call to get correct input parameters 
c ... without hardwiring

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )

c ... Second call will replace variables with correct inputs

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )
      
      if (rtn .ne. 0) then
        call message( 'Get_Global_Metadata_Leocat',
     &    'Error getting Snow/Ice Processed from C6 file' //
     &    ' [OPERATOR ACTION: Stage correct C6, rerun PGE]', 0, 2 )
      end if

c ... Set return value

      pcnni = value 

c ... Get number of land processed for this granule

      dtype = '  '
      attrib = 'LandProcessedPct'
      nelmnt = 0
      value = 0.0
      rtn = 0

c ... Make initial call to get correct input parameters 
c ... without hardwiring

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )

c ... Second call will replace variables with correct inputs

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )
      
      if (rtn .ne. 0) then
        call message( 'Get_Global_Metadata_Leocat',
     &    'Error getting Land processed from C6 file' //
     &    ' [OPERATOR ACTION: Stage correct C6, rerun PGE]', 0, 2 )
      end if

c ... Set return value

      pcnnl = value 

c ... Get number of water processed for this granule

      dtype = '  '
      attrib = 'WaterProcessedPct'
      nelmnt = 0
      value = 0.0
      rtn = 0

c ... Make initial call to get correct input parameters 
c ... without hardwiring

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )

c ... Second call will replace variables with correct inputs

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )
      
      if (rtn .ne. 0) then
        call message( 'Get_Global_Metadata_Leocat',
     &    'Error getting Water processed from C6 file' //
     &    ' [OPERATOR ACTION: Stage correct C6, rerun PGE]', 0, 2 )
      end if

c ... Set return value

      pcnnw = value 

c ... Get number of thin cirrus solar processed for this granule

      dtype = '  '
      attrib = 'ThinCirrusSolarFoundPct'
      nelmnt = 0
      value = 0.0
      rtn = 0

c ... Make initial call to get correct input parameters 
c ... without hardwiring

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )

c ... Second call will replace variables with correct inputs

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )
      
      if (rtn .ne. 0) then
        call message( 'Get_Global_Metadata_Leocat',
     &    'Error getting Thin Cirrus Solar processed from C6 file' //
     &    ' [OPERATOR ACTION: Stage correct C6, rerun PGE]', 0, 2 )
      end if

c ... Set return value

      pcnnv = value 

c ... Get number of thin cirrus infrared processed for this granule

      dtype = '  '
      attrib = 'ThinCirrusIR_FoundPct'
      nelmnt = 0
      value = 0.0
      rtn = 0

c ... Make initial call to get correct input parameters 
c ... without hardwiring

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )

c ... Second call will replace variables with correct inputs

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )
      
      if (rtn .ne. 0) then
        call message( 'Get_Global_Metadata_Leocat',
     &    'Error getting Thin Cirrus IR processed from C6 file' //
     &    ' [OPERATOR ACTION: Stage correct C6, rerun PGE]', 0, 2 )
      end if

c ... Set return value

      pcnnr = value 

c ... Get number of non cloud obstruction processed for this granule

      dtype = '  '
      attrib = 'NonCloudObstructionFoundPct'
      nelmnt = 0
      value = 0.0
      rtn = 0

c ... Make initial call to get correct input parameters 
c ... without hardwiring

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )

c ... Second call will replace variables with correct inputs

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )
      
      if (rtn .ne. 0) then
        call message( 'Get_Global_Metadata_Leocat',
     &    'Error getting Noncloud obstruction processed from C6 file' //
     &    ' [OPERATOR ACTION: Stage correct C6, rerun PGE]', 0, 2 )
      end if

c ... Set return value

      pcnnc = value 

c ... Get value of maximum solar zenith angle

      dtype = '  '
      attrib = 'MaxSolarZenithAngle'
      nelmnt = 0
      value = 0.0
      rtn = 0

c ... Make initial call to get correct input parameters 
c ... without hardwiring

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )

c ... Second call will replace variables with correct inputs

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )
      
      if (rtn .ne. 0) then
        call message( 'Get_Global_Metadata_Leocat',
     &    'Error getting max solar zenith processed from C6 file' //
     &    ' [OPERATOR ACTION: Stage correct C6, rerun PGE]', 0, 2 )
      end if

c ... Set return value

      max_sol = value 

c ... Get value of minimum solar zenith angle

      dtype = '  '
      attrib = 'MinSolarZenithAngle'
      nelmnt = 0
      value = 0.0
      rtn = 0

c ... Make initial call to get correct input parameters 
c ... without hardwiring

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )

c ... Second call will replace variables with correct inputs

      rtn = gmfin( modfil, attrib, dtype, nelmnt, value )
      
      if (rtn .ne. 0) then
        call message( 'Get_Global_Metadata_Leocat',
     &    'Error getting min solar zenith processed from C6 file' //
     &    ' [OPERATOR ACTION: Stage correct C6, rerun PGE]', 0, 2 )
      end if

c ... Set return value

      min_sol = value 

c ... Get value of NDVI file

      dtype = '  '
      attrib = 'NDVIFile'
      nelmnt = 0
      charval = '  '
      rtn = 0

c ... Make initial call to get correct input parameters 
c ... without hardwiring

      rtn = gmfin( modfil, attrib, dtype, nelmnt, charval )

c ... Second call will replace variables with correct inputs
c ... Add one to  nelmnt to account for difference in C/FORTRAN strings
      rtn = gmfin( modfil, attrib, dtype, nelmnt+1, charval )
      
      if (rtn .ne. 0) then
        call message( 'Get_Global_Metadata_Leocat',
     &    'Error getting NDVI filename processed from C6 file' //
     &    ' [OPERATOR ACTION: Stage correct C6, rerun PGE]', 0, 2 )
      end if

c ... Set return value

      ndvifil = charval

c ... Get value of Thresholds file

      dtype = '  '
      attrib = 'ThresholdFile'
      nelmnt = 0
      charval = '  '
      rtn = 0

c ... Make initial call to get correct input parameters 
c ... without hardwiring

      rtn = gmfin( modfil, attrib, dtype, nelmnt, charval )

c ... Second call will replace variables with correct inputs
c ... Add one to  nelmnt to account for difference in C/FORTRAN strings
      rtn = gmfin( modfil, attrib, dtype, nelmnt+1, charval )
      
      if (rtn .ne. 0) then
        call message( 'Get_Global_Metadata_Leocat',
     &    'Error getting Threshold filename processed from C6 file' //
     &    ' [OPERATOR ACTION: Stage correct C6, rerun PGE]', 0, 2 )
      end if

c ... Set return value

      thresfil = charval

      END
