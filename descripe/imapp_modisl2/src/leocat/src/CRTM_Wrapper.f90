!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

MODULE CRTM_Wrapper

CONTAINS

  SUBROUTINE Run_Forward_Model (S_ID, &
                                nlevels_actual, &
				SURFACE_TEMPERATURE, &
				SECANT_VIEW_ANGLE, &
				SECANT_SOLAR_ANGLE, &
				Surface_Emissivity, &
				Surface_Reflectivity, &
				LEVEL_HEIGHT_IN, &
				LEVEL_PRESSURE_IN, &
				LEVEL_TEMPERATURE_IN, &
				LEVEL_WATER_VAPOR_IN, &
				LEVEL_OZONE_IN, &
				Z_TOP, &
				P_TOP, &
				T_TOP, &
				W_TOP, &
				o3_TOP, &
				LAYER_TEMPERATURE_OUT, &
				Brightness_Temperature, &
				Upwelling_Radiance, &
				tau)


    ! ------------
    ! Module usage
    ! ------------

    ! -- Utility modules
    USE Type_Kinds
    USE Error_Handler


    ! -- pCRTM modules
    USE Initialize
    USE Parameters
    USE Forward_Model
    
    USE Units_Conversion
    USE Level_Layer_Conversion


    ! ---------------------------
    ! Disable all implicit typing
    ! ---------------------------

    IMPLICIT NONE

    ! ----------
    ! Parameters
    ! ----------

    CHARACTER( * ), PARAMETER :: PROGRAM_NAME   = 'CRTM_WRAPPER'
    INTEGER, PARAMETER :: NCHANNELS_CRTM = 16
    INTEGER, PARAMETER :: MAX_LAYERS_CRTM = 100
    INTEGER, PARAMETER :: NPROF = 1
    
    ! --------------
    ! Input
    ! --------------
    
    INTEGER, INTENT( IN ) :: S_ID
    INTEGER, INTENT( IN ) :: nlevels_actual
    REAL( fp_kind ), INTENT ( IN ), DIMENSION( NPROF ) :: SURFACE_TEMPERATURE
    REAL( fp_kind ), INTENT ( IN ), DIMENSION( NPROF ) :: SECANT_VIEW_ANGLE
    REAL( fp_kind ), INTENT ( IN ), DIMENSION( NPROF ) :: SECANT_SOLAR_ANGLE
    REAL( fp_kind ), INTENT ( IN ), DIMENSION( NCHANNELS_CRTM*NPROF ) :: Surface_Emissivity
    REAL( fp_kind ), INTENT ( IN ), DIMENSION( NCHANNELS_CRTM*NPROF ) :: Surface_Reflectivity
    REAL( fp_kind ), INTENT ( IN ), DIMENSION( MAX_LAYERS_CRTM ) :: LEVEL_HEIGHT_IN
    REAL( fp_kind ), INTENT ( IN ), DIMENSION( MAX_LAYERS_CRTM ) :: LEVEL_PRESSURE_IN
    REAL( fp_kind ), INTENT ( IN ), DIMENSION( MAX_LAYERS_CRTM ) :: LEVEL_TEMPERATURE_IN
    REAL( fp_kind ), INTENT ( IN ), DIMENSION( MAX_LAYERS_CRTM ) :: LEVEL_WATER_VAPOR_IN
    REAL( fp_kind ), INTENT ( IN ), DIMENSION( MAX_LAYERS_CRTM ) :: LEVEL_OZONE_IN
    
    REAL( fp_kind ), DIMENSION( MAX_LAYERS_CRTM+1 ) :: LEVEL_HEIGHT_TMP
    REAL( fp_kind ), DIMENSION( MAX_LAYERS_CRTM+1 ) :: LEVEL_PRESSURE_TMP
    REAL( fp_kind ), DIMENSION( MAX_LAYERS_CRTM+1 ) :: LEVEL_TEMPERATURE_TMP
    REAL( fp_kind ), DIMENSION( MAX_LAYERS_CRTM+1 ) :: LEVEL_WATER_VAPOR_TMP
    REAL( fp_kind ), DIMENSION( MAX_LAYERS_CRTM+1 ) :: LEVEL_OZONE_TMP
    REAL( fp_kind ), DIMENSION( MAX_LAYERS_CRTM+1 ) :: LEVEL_e
    
    REAL( fp_kind ), INTENT ( IN ), DIMENSION( NPROF ) :: Z_TOP
    REAL( fp_kind ), INTENT ( IN ), DIMENSION( NPROF ) :: P_TOP
    REAL( fp_kind ), INTENT ( IN ), DIMENSION( NPROF ) :: T_TOP
    REAL( fp_kind ), INTENT ( IN ), DIMENSION( NPROF ) :: W_TOP
    REAL( fp_kind ), INTENT ( IN ), DIMENSION( NPROF ) :: o3_TOP
    
    REAL( fp_kind ), DIMENSION( MAX_LAYERS_CRTM ) :: LAYER_PRESSURE_OUT
    REAL( fp_kind ), DIMENSION( MAX_LAYERS_CRTM ) :: LAYER_TEMPERATURE_OUT
    REAL( fp_kind ), DIMENSION( MAX_LAYERS_CRTM ) :: LAYER_WATER_VAPOR_OUT
    REAL( fp_kind ), DIMENSION( MAX_LAYERS_CRTM ) :: LAYER_OZONE_OUT
    
    REAL( fp_kind ), DIMENSION( :, : ), ALLOCATABLE :: LEVEL_PRESSURE
    REAL( fp_kind ), DIMENSION( :, : ), ALLOCATABLE :: LAYER_PRESSURE
    REAL( fp_kind ), DIMENSION( :, : ), ALLOCATABLE :: LAYER_TEMPERATURE
    REAL( fp_kind ), DIMENSION( :, : ), ALLOCATABLE :: LAYER_WATER_VAPOR
    REAL( fp_kind ), DIMENSION( :, : ), ALLOCATABLE :: LAYER_OZONE


    ! ---------
    ! Variables
    ! ---------

    ! -- Variable dimension determined during initialisation
    INTEGER :: n_Channels  ! L dimension

    ! -- Other forward model inputs
    !REAL( fp_kind ), DIMENSION( NCHANNELS_CRTM*NPROF ) :: Surface_Emissivity         ! L*M  
    !REAL( fp_kind ), DIMENSION( NCHANNELS_CRTM*NPROF ) :: Surface_Reflectivity       ! L*M  
    INTEGER,         DIMENSION( NPROF )     :: n_Channels_Per_Profile     ! M  
    INTEGER,         DIMENSION( NCHANNELS_CRTM*NPROF ) :: Channel_Index              ! L*M 


    ! -- Forward outputs                                                                                
    REAL( fp_kind ), DIMENSION( MAX_LAYERS_CRTM, NCHANNELS_CRTM*NPROF ) :: Tau                     ! K x L*M
    REAL( fp_kind ), DIMENSION( MAX_LAYERS_CRTM, NCHANNELS_CRTM*NPROF ) :: Flux_Tau                ! K x L*M
    REAL( fp_kind ), DIMENSION( MAX_LAYERS_CRTM, NCHANNELS_CRTM*NPROF ) :: Solar_Tau               ! K x L*M
    REAL( fp_kind ), DIMENSION( NCHANNELS_CRTM*NPROF ) :: Upwelling_Radiance      ! L*M  
    REAL( fp_kind ), DIMENSION( NCHANNELS_CRTM*NPROF ) :: Brightness_Temperature  ! L*M  


    ! -- Optional inputs
    REAL( fp_kind ), DIMENSION( NCHANNELS_CRTM*NPROF ) :: Solar_Reflectivity         ! L*M  
    REAL( fp_kind ), DIMENSION( NPROF )     :: Secant_Flux_Angle          ! M


    ! -- Error status variables
    INTEGER :: Error_Status
    INTEGER :: Allocate_Status


    ! -- Some integers to play around with
    ! -- the number of channels to process
    INTEGER :: n_Channel_Skip
    INTEGER :: n_Channels_to_Process


    ! -- Loop counters
    INTEGER :: l, m, lm, k


    ! -- Coefficient filenames
    CHARACTER( 256 ) :: SpcCoeff_File
    CHARACTER( 256 ) :: TauCoeff_File

    !#----------------------------------------------------------------------------#
    !#              -- GET THE REQUIRED COEFFICIENT DATAFILE NAMES --             #
    !#----------------------------------------------------------------------------#

    ! ----------------------
    ! The SpcCoeff data file
    ! ----------------------

    if (S_ID == 0) then
      !print*,"CRTM  - using Terra MODIS"
      SpcCoeff_File = TRIM( "../data/modis_terra.SpcCoeff.bin" )
      TauCoeff_File = TRIM( "../data/modis_terra.TauCoeff.bin" )
    elseif (S_ID == 1) then
      !print*,"CRTM  - using Aqua MODIS"
      SpcCoeff_File = TRIM( "../data/modis_aqua.SpcCoeff.bin" )
      TauCoeff_File = TRIM( "../data/modis_aqua.TauCoeff.bin" )
    endif


    !#----------------------------------------------------------------------------#
    !#                      -- INITIALISE THE pCRTM MODEL --                      # 
    !#                                                                            #
    !#                 This function is in the INITIALIZE module                  #
    !#----------------------------------------------------------------------------#

    Error_Status = Initialize_RTM( Spectral_File = ADJUSTL( SpcCoeff_File ), &
                                   Tau_File      = ADJUSTL( TauCoeff_File ), &
				   quiet         = 1  )

    IF ( Error_Status /= SUCCESS ) THEN 
       CALL Display_Message( PROGRAM_NAME, &
                             'Error initializing the pCRTM', & 
                              Error_Status)  
      STOP
    END IF



    !#----------------------------------------------------------------------------#
    !#                -- ALLOCATE THE CHANNEL DEPENDENT ARRAYS --                 #
    !#                                                                            #
    !# Rather than hard-wire the code for a particular number of channels, the    #
    !# following allocations are done so that:                                    #
    !#   1) the channel dimension can be determined dynamically based on the      #
    !#      channel dimension of the input coefficient data files, and            #
    !#   2) there is an illustration of the use of the input arguments            #
    !#        n_Channels_per_Profile                                              #
    !#      and                                                                   #
    !#        Channel_Index                                                       #
    !#                                                                            #
    !# So, if the next block of code seems overly obtuse, you can always set the  #
    !# number of channels to some fixed value, declare the channel-dependent      #
    !# arrays accordingly, and avoid the allocations.                             #
    !#----------------------------------------------------------------------------#

    ! ----------------------------------------------
    ! Retrieve the number of channels defined during
    ! the initialization.
    !
    ! This subroutine is in the PARAMETERS module.
    ! ----------------------------------------------

    CALL Get_Max_n_Channels( n_Channels )
    if (n_Channels /= NCHANNELS_CRTM) then
      CALL Display_Message( PROGRAM_NAME, &
                            'The number of channels from spectral files does not match MODIS', & 
                             Error_Status)  
      STOP
    END IF


    ! ----------------------------------------------------
    ! Get the number of channels to skip in the processing
    !
    ! This is just to illustrate how what is
    ! contained in the input arrays
    !   n_Channels_per_Profile
    ! and
    !   Channel_Index
    ! controls the channel processing
    ! ----------------------------------------------------

    n_Channel_Skip = 1

    ! -- So... let's skip some channels instead of doing them all
    n_Channels_to_Process = n_Channels / n_Channel_Skip
    IF ( MOD( n_Channels, n_Channel_Skip ) /= 0 ) n_Channels_to_Process = n_Channels_to_Process + 1



    !#----------------------------------------------------------------------------#
    !#                   -- FILL THE REMAINING INPUT ARRAYS --                    #
    !#----------------------------------------------------------------------------#
    
    ALLOCATE( LEVEL_PRESSURE( nlevels_actual, NPROF ), &
              LAYER_PRESSURE( nlevels_actual, NPROF ), &
	      LAYER_TEMPERATURE( nlevels_actual, NPROF), &
	      LAYER_WATER_VAPOR( nlevels_actual, NPROF), &
	      LAYER_OZONE( nlevels_actual, NPROF), &
	      STAT = Allocate_Status )
	     
    IF ( Allocate_Status /= 0 ) THEN
      CALL Display_Message( PROGRAM_NAME, &
                            'Error allocating profile arrays', & 
                             FAILURE )  
      STOP
    END IF
    
    LEVEL_HEIGHT_TMP(2:MAX_LAYERS_CRTM+1) = LEVEL_HEIGHT_IN
    LEVEL_HEIGHT_TMP(1) = Z_TOP(1)
    
    LEVEL_PRESSURE_TMP(2:MAX_LAYERS_CRTM+1) = LEVEL_PRESSURE_IN
    LEVEL_PRESSURE_TMP(1) = P_TOP(1)
    
    LEVEL_TEMPERATURE_TMP(2:MAX_LAYERS_CRTM+1) = LEVEL_TEMPERATURE_IN
    LEVEL_TEMPERATURE_TMP(1) = T_TOP(1)
    
    LEVEL_WATER_VAPOR_TMP(2:MAX_LAYERS_CRTM+1) = LEVEL_WATER_VAPOR_IN
    LEVEL_WATER_VAPOR_TMP(1) = W_TOP(1)
    
    LEVEL_OZONE_TMP(2:MAX_LAYERS_CRTM+1) = LEVEL_OZONE_IN
    LEVEL_OZONE_TMP(1) = o3_TOP(1)
    
    LEVEL_e = MR_to_PP(LEVEL_PRESSURE_TMP, LEVEL_WATER_VAPOR_TMP)
    Error_Status = Effective_Layer_TP(LEVEL_HEIGHT_TMP, &
                                      LEVEL_PRESSURE_TMP, &
                                      LEVEL_TEMPERATURE_TMP, &
                                      LEVEL_e,        & 
                                      LAYER_PRESSURE_OUT, &
                                      LAYER_TEMPERATURE_OUT)	
    
    DO k=1, nlevels_actual
      LEVEL_PRESSURE(k,1) = LEVEL_PRESSURE_IN(k)
      LAYER_WATER_VAPOR_OUT(k) = (LEVEL_WATER_VAPOR_TMP(k) + LEVEL_WATER_VAPOR_TMP(k+1))/2.0
      LAYER_OZONE_OUT(k) = (LEVEL_OZONE_TMP(k) + LEVEL_OZONE_TMP(k+1))/2.0
      
      LAYER_PRESSURE(k,1) = LAYER_PRESSURE_OUT(k)
      LAYER_TEMPERATURE(k,1) = LAYER_TEMPERATURE_OUT(k)
      LAYER_WATER_VAPOR(k,1) = LAYER_WATER_VAPOR_OUT(k)
      LAYER_OZONE(k,1) = LAYER_OZONE_OUT(k)
      
     !write (unit=6,fmt="(5(f10.4,2x))") LEVEL_PRESSURE(k,1),LAYER_PRESSURE(k,1),LAYER_TEMPERATURE(k,1), &
     !LAYER_WATER_VAPOR(k,1),LAYER_OZONE(k,1)
    END DO    
    
    ! -----------------------------
    ! Surface and flux angle inputs
    ! -----------------------------

    !Surface_Emissivity   = 1.0_fp_kind
    !Surface_Reflectivity = ONE - Surface_Emissivity

    Solar_Reflectivity = ZERO
    Secant_Flux_Angle  = SECANT_DIFFUSIVITY_ANGLE


    ! --------------------------------------
    ! Process this many channels per profile
    ! --------------------------------------

    n_Channels_per_Profile = n_Channels_to_Process


    ! ---------------------------------
    ! The index of required channels in
    ! the coefficient data structures.
    ! ---------------------------------

    Channel_Index( 1:n_Channels_to_Process*NPROF ) = &
      (/ (( l, l = 1, n_Channels, n_Channel_Skip ), m = 1, NPROF ) /)



    !#----------------------------------------------------------------------------#
    !#                         -- CALL THE FORWARD MODEL --                       #
    !#                                                                            #
    !#                This function is in the FORWARD_MODEL module                #
    !#----------------------------------------------------------------------------#
    
    Error_Status = Compute_RTM( LEVEL_PRESSURE,         &  ! Input, K x M
                                LAYER_PRESSURE,         &  ! Input, K x M 
                                LAYER_TEMPERATURE,      &  ! Input, K x M 
                                LAYER_WATER_VAPOR,      &  ! Input, K x M 
                                LAYER_OZONE,            &  ! Input, K x M    
                                SURFACE_TEMPERATURE,    &  ! Input, M     
                                Surface_Emissivity,     &  ! Input, L*M     
                                Surface_Reflectivity,   &  ! Input, L*M       
                                SECANT_VIEW_ANGLE,      &  ! Input, M
                                SECANT_SOLAR_ANGLE,     &  ! Input, M   
                                n_Channels_Per_Profile, &  ! Input, M
                                Channel_Index,          &  ! Input, L*M   
                                Tau,                    &  ! Output, K x L*M  
                                Flux_Tau,               &  ! Output, K x L*M   
                                Solar_Tau,              &  ! Output, K x L*M  
                                Upwelling_Radiance,     &  ! Output, L*M  
                                Brightness_Temperature, &  ! Output, L*M
                                Solar_Reflectivity = Solar_Reflectivity, &  ! Optional input, L*M 
                                Secant_Flux_Angle  = Secant_Flux_Angle   )  ! Optional input, M

    IF ( Error_Status /= SUCCESS ) THEN 
       CALL Display_Message( PROGRAM_NAME, &
                             'Error in Compute_RTM call', & 
                              Error_Status )                           
     STOP
    END IF



    !#----------------------------------------------------------------------------#
    !#                        -- OUTPUT RESULTS AND CLEAN UP --                   #
    !#----------------------------------------------------------------------------#

    DEALLOCATE( LEVEL_PRESSURE, &
                LAYER_PRESSURE, &
		LAYER_TEMPERATURE, &
		LAYER_WATER_VAPOR, &
		LAYER_OZONE, &
		STAT = Allocate_Status )
		
    IF ( Allocate_Status /= 0 ) THEN
      CALL Display_Message( PROGRAM_NAME, &
                            'Error deallocating profile arrays', & 
                             WARNING )  
    END IF

    !#----------------------------------------------------------------------------#
    !#                        -- DESTROY THE pCRTM SPACE --                       # 
    !#                                                                            #
    !#                 This function is in the INITIALIZE module                  #
    !#----------------------------------------------------------------------------#

    !WRITE( *, '( /5x, "Destroying the pCRTM space..." )' )

    Error_Status = Destroy_RTM()

    IF ( Error_Status /= SUCCESS ) THEN 
       CALL Display_Message( PROGRAM_NAME, &
                             'Error destroying the pCRTM space', & 
                              Error_Status )
     STOP
    END IF
  
  END SUBROUTINE Run_Forward_Model

END MODULE CRTM_Wrapper
