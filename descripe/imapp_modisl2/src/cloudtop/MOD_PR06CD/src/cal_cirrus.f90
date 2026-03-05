        SUBROUTINE CAL_CIRRUS(N_PARTITION_X, N_PARTITION_Y, &
                              REFL_BAND_1,REFL_BAND_26, &
                              CENT_SOLAR_ANG, CENT_VIEW_ANG, &
                              CIRRUS_REFL, N_X, N_Y, &
                        Max_Frames_Per_Line,Max_Lines_Per_Img)
!
!...................................................................
!
!        This is the core code for deriving cirrus reflectance using MODIS
!        visible and 1.38 micron (band 26) data.
! 
!
!       History:
!
!        first version (three slopes for each sub-image)  4/23/2001
!
!         Written by 
!           Dr. Ping Yang 
!           GEST, UMBC / Code 913, NASA-GSFC
!           email: pyang@climate.gsfc.nasa.gov
!           phone: 301-614-6127
!           
!
!
!       second version (one slope for each subimage. If a sub-image is
!       cirrus free, the derived slope may not be correct. In this case, the slope
!       is calculated using the transmittance of h2O based on an analytical expression
!       derived from fitting of  Gao's look-up table of the transmittance)
!
!       written by 
!           Dr. Ping Yang
!           Department of Atmospheric Sciences
!           Texas A&M University
!           Tel:979-845-4923
!           email: pyang@ariel.met.tamu.edu
!                            Nov. 2, 2001 
!
!       This subroutine first calls subroutine "REMOVE_SURFACE_EFFECT" that
!       returns the slope information for subimages. Then the cirrus 
!       reflectance is calculated for each subimages. Note that we have
!       applied a four-points interpolation scheme in calculating the reflectance 
!       to avoid the "chess board" effect.
!......................................................................


        REAL, DIMENSION(Max_Frames_Per_Line,Max_Lines_Per_Img), &
                       INTENT(IN) :: REFL_BAND_1,REFL_BAND_26
        INTEGER, INTENT(IN) :: N_X, N_Y, Max_Frames_Per_Line, &
                               Max_Lines_Per_Img
        INTEGER, INTENT(IN) :: N_PARTITION_X, N_PARTITION_Y
        REAL, DIMENSION(N_PARTITION_X,N_PARTITION_Y) :: CENT_SOLAR_ANG, CENT_VIEW_ANG
        REAL, DIMENSION(Max_Frames_Per_Line,Max_Lines_Per_Img), &
                              INTENT(OUT) :: CIRRUS_REFL
!

        INTERFACE

       SUBROUTINE COMPUTE_CIRRUS_REFLECTANCE(N_PARTITION_X, &
                                             N_PARTITION_Y, &
                   DATA_SOLAR_ZENITH, DATA_VIEW_ZENITH,     &
                   REFL_BAND_1,REFL_BAND_26, &
                DERIVED_CIRRUS_REFLECTION, N_X, N_Y, &
                   Max_Frames_Per_Line,Max_Lines_Per_Img)


        INTRINSIC SELECTED_INT_KIND, KIND, SIZE, SUM, MINVAL, MAXVAL

        INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
        INTEGER, PARAMETER :: SP = KIND(1.0)

        INTEGER, INTENT(IN) :: N_PARTITION_X, N_PARTITION_Y,  &
                               N_X, N_Y, Max_Frames_Per_Line, &
                               Max_Lines_Per_Img

        REAL, DIMENSION(Max_Frames_Per_Line,Max_Lines_Per_Img), &
                            INTENT(IN) :: REFL_BAND_1,REFL_BAND_26

        REAL, DIMENSION(N_PARTITION_X,N_PARTITION_Y), &
                INTENT(IN) :: DATA_SOLAR_ZENITH, DATA_VIEW_ZENITH

        REAL, DIMENSION(Max_Frames_Per_Line,Max_Lines_Per_Img), &
                 INTENT(OUT) :: DERIVED_CIRRUS_REFLECTION


!        INTEGER(I4B), PARAMETER ::  N_X=1354, N_Y=2030
!                          N_X and N_Y are the x and y dimensions of input image
!
!        INTEGER, PARAMETER :: N_PARTITION_X = 4, N_PARTITION_Y = 4 
!
!             N_PARTITION_X= number of sub-images along x-direction
!             N_PARTITION_Y= number of sub-images along Y-direction


        INTEGER(I4B), PARAMETER :: N_FOR_BAR_AVERAGE = 10, N_SLOPE=20


!          N_SLOPE is the number of binned slices for 1.38 um reflectance
!          The number of points  in each slice to get average values

        REAL(SP), PARAMETER :: X_LOWER_BOUND = 0.0,  &
                               X_UPPER_BOUND = 2000.0, &
                               Y_LOWER_BOUND = 0.0,  &
                               Y_UPPER_BOUND = 1000.0

        REAL(SP), PARAMETER :: REFL_138_LOW_PERCENTAGE = 0.7, &
                               REFL_138_HIGH_PERCENTAGE = 0.999, &
                               REFL_138_L_REFLECTION = 30.0, &
                               REFL_138_H_REFLECTION = 1000.0

        REAL(SP) :: OUT_SLOPE         ! the output of subroutine REMOVE_SURFACE_EFFECT

        INTEGER, DIMENSION(N_PARTITION_X + 1) ::  I_X  ! nodes for partition along x
        INTEGER, DIMENSION(N_PARTITION_Y + 1) ::  I_Y  ! nodes for partition along y

        REAL(SP), DIMENSION(N_PARTITION_X,N_PARTITION_Y) :: &
                                              SLOPE_AT_CENTER
!                                These are slopes at the centers of sub-images

        REAL(SP), DIMENSION(N_PARTITION_X+1,N_PARTITION_Y+1) :: &
                                             SLOPE_AT_NODE
!                                These are slopes at the corners of subimages


        REAL(SP), DIMENSION(N_X,N_Y) :: DATA_FLOAT_138, &
                                        DATA_FLOAT_066_OR_124


        REAL(SP), DIMENSION(:), ALLOCATABLE :: SUB_DATA_138, &
                                               SUB_DATA_066_OR_124

        INTEGER(I4B) :: N_SIZE_SUB_DATA          !size of dynamic arrays SUB_DATA_138, 
!                        and SUB_DATA_066_OR_124, which is determined in running time.

        INTEGER(I4B) :: INDEX_QA_SLOPE   


       END SUBROUTINE COMPUTE_CIRRUS_REFLECTANCE

!..........................................................................

        END INTERFACE


        CALL COMPUTE_CIRRUS_REFLECTANCE(N_PARTITION_X, &
                                        N_PARTITION_Y, &
                   CENT_SOLAR_ANG, CENT_VIEW_ANG,      &
                   REFL_BAND_1,REFL_BAND_26, &
                   CIRRUS_REFL, N_X, N_Y, &
                   Max_Frames_Per_Line,Max_Lines_Per_Img)

        END SUBROUTINE CAL_CIRRUS
!
!
       SUBROUTINE COMPUTE_CIRRUS_REFLECTANCE(N_PARTITION_X, &
                                             N_PARTITION_Y, &
                   DATA_SOLAR_ZENITH, DATA_VIEW_ZENITH,     &
                   REFL_BAND_1,REFL_BAND_26, &
                DERIVED_CIRRUS_REFLECTION, N_X, N_Y, &
                   Max_Frames_Per_Line,Max_Lines_Per_Img)


        INTRINSIC SELECTED_INT_KIND, KIND, SIZE, SUM, MINVAL, MAXVAL

        INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
        INTEGER, PARAMETER :: SP = KIND(1.0)

        INTEGER, INTENT(IN) :: N_PARTITION_X, N_PARTITION_Y,  &
                               N_X, N_Y, Max_Frames_Per_Line, &
                               Max_Lines_Per_Img

        REAL, DIMENSION(Max_Frames_Per_Line,Max_Lines_Per_Img), &
                            INTENT(IN) :: REFL_BAND_1,REFL_BAND_26

        REAL, DIMENSION(N_PARTITION_X,N_PARTITION_Y), &
                INTENT(IN) :: DATA_SOLAR_ZENITH, DATA_VIEW_ZENITH

        REAL, DIMENSION(Max_Frames_Per_Line,Max_Lines_Per_Img), &
                 INTENT(OUT) :: DERIVED_CIRRUS_REFLECTION


!        INTEGER(I4B), PARAMETER ::  N_X=1354, N_Y=2030
!                          N_X and N_Y are the x and y dimensions of input image
!
!        INTEGER, PARAMETER :: N_PARTITION_X = 4, N_PARTITION_Y = 4 
!
!             N_PARTITION_X= number of sub-images along x-direction
!             N_PARTITION_Y= number of sub-images along Y-direction


        INTEGER(I4B), PARAMETER :: N_FOR_BAR_AVERAGE = 10, N_SLOPE=20


!          N_SLOPE is the number of binned slices for 1.38 um reflectance
!          The number of points  in each slice to get average values

        REAL(SP), PARAMETER :: X_LOWER_BOUND = 0.0,  &
                               X_UPPER_BOUND = 2000.0, &
                               Y_LOWER_BOUND = 0.0,  &
                               Y_UPPER_BOUND = 1000.0

        REAL(SP), PARAMETER :: REFL_138_LOW_PERCENTAGE = 0.7, &
                               REFL_138_HIGH_PERCENTAGE = 0.999, &
                               REFL_138_L_REFLECTION = 30.0, &
                               REFL_138_H_REFLECTION = 1000.0

        REAL(SP) :: OUT_SLOPE         ! the output of subroutine REMOVE_SURFACE_EFFECT

        INTEGER, DIMENSION(N_PARTITION_X + 1) ::  I_X  ! nodes for partition along x
        INTEGER, DIMENSION(N_PARTITION_Y + 1) ::  I_Y  ! nodes for partition along y

        REAL(SP), DIMENSION(N_PARTITION_X,N_PARTITION_Y) :: &
                                              SLOPE_AT_CENTER
!                                These are slopes at the centers of sub-images

        REAL(SP), DIMENSION(N_PARTITION_X+1,N_PARTITION_Y+1) :: &
                                             SLOPE_AT_NODE
!                                These are slopes at the corners of subimages


        REAL(SP), DIMENSION(N_X,N_Y) :: DATA_FLOAT_138, &
                                        DATA_FLOAT_066_OR_124


        REAL(SP), DIMENSION(:), ALLOCATABLE :: SUB_DATA_138, &
                                               SUB_DATA_066_OR_124

        INTEGER(I4B) :: N_SIZE_SUB_DATA          !size of dynamic arrays SUB_DATA_138, 
!                        and SUB_DATA_066_OR_124, which is determined in running time.

        INTEGER(I4B) :: INDEX_QA_SLOPE   



!..........................................................................

        INTERFACE

               
        SUBROUTINE REMOVE_SURFACE_EFFECT(INDEX_QA_SLOPE, &
                                         DATA_FLOAT_138, &
           DATA_FLOAT_066_OR_124, N_SLOPE,                  &
           N_FOR_BAR_AVERAGE, REFL_138_LOW_PERCENTAGE,      &
                              REFL_138_HIGH_PERCENTAGE,     &
           REFL_138_L_REFLECTION, REFL_138_H_REFLECTION,    & 
           OUT_SLOPE )

!
! DESCRIPTION: This subroutine is first to select the useful data of 1.38 and 
!   0.66 (or 1.24) um reflectance. Further, it removes the surface effect in 
!   the co-relation plot of the two channels.  
!
        INTRINSIC SELECTED_INT_KIND, KIND, SIZE, SUM 

        INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
        INTEGER, PARAMETER :: SP = KIND(1.0)
        REAL(SP), PARAMETER :: X_REFL_H_REFLECTION = 1000.0 
        REAL(SP), PARAMETER :: X_REFL_L_REFLECTION = 0.0

        INTEGER(I4B), INTENT(IN) :: N_SLOPE
        INTEGER(I4B), INTENT(OUT) :: INDEX_QA_SLOPE
!          N_SLOPE is the number of binned slices for 1.38 um reflectance

        INTEGER(I4B), INTENT(IN) :: N_FOR_BAR_AVERAGE ! the number of points 
                                            ! in each slice to get average values

        REAL(SP), INTENT(IN) :: REFL_138_LOW_PERCENTAGE,  &
                                REFL_138_HIGH_PERCENTAGE, &
                                REFL_138_L_REFLECTION,    &
                                REFL_138_H_REFLECTION

        REAL(SP), INTENT(OUT) :: OUT_SLOPE

!           These are the output of this subroutine. OUT_X_RAYLEIGH is the 
!       x-axis interception of the lower part of the derived 0.66-1.38 curve.
!       OUT_Y_DIVISION determines the division of the two parts. 
   

        INTEGER(I4B), DIMENSION(N_SLOPE+1) :: N_SUM
        REAL(SP), DIMENSION(N_SLOPE+1) :: Y_138_NODE

        REAL(SP), DIMENSION(N_SLOPE) :: &
                          X_066_OR_124_BAR, Y_138_BAR



        REAL(SP), DIMENSION(:), INTENT(INOUT) :: DATA_FLOAT_138, &
                     DATA_FLOAT_066_OR_124


        REAL(SP), DIMENSION(:), ALLOCATABLE ::  SAMPLE_138, &
                                                SAMPLE_066_OR_124

        CHARACTER(2) :: FILE_NAME_NUMBER_X,  FILE_NAME_NUMBER_Y

        END SUBROUTINE REMOVE_SURFACE_EFFECT
        END INTERFACE

!..................................................................

! rhucek 8/26/05: Changed initialization from 0.0 to 1999.80 on
!    request from Bill Ridgway. 
!      DERIVED_CIRRUS_REFLECTION = 0.0
       DERIVED_CIRRUS_REFLECTION = -1999.80 

       DATA_FLOAT_066_OR_124(1:N_X,1:N_Y) = &
                    1000.0 * REFL_BAND_1(1:N_X,1:N_Y)
       DATA_FLOAT_138(1:N_X,1:N_Y) = &
                    1000.0 * REFL_BAND_26(1:N_X,1:N_Y)


!
!              determine partition nodes

        I_X(1) = 0
        DO I = 2, N_PARTITION_X 
          I_X(I) = IFIX( FLOAT(I-1) &
             * (FLOAT(N_X) / FLOAT(N_PARTITION_X) ) )
        END DO
        I_X(N_PARTITION_X+1) = N_X
        
        I_Y(1) = 0
        DO I = 2, N_PARTITION_Y 
          I_Y(I) = IFIX(FLOAT(I-1) &
             * (FLOAT(N_Y) / FLOAT(N_PARTITION_Y) ) )
        END DO
        I_Y(N_PARTITION_Y+1) = N_Y


!......... determine the envelope at the centers of subimages  ........
!                ******* important *********
!          The reflectance has been multplied by a factor of 1000 for
!          slicing the data in deriving the slope in 1.38 and 0.66 (1.24) scatter-plot
!
!       ............. slopes  at centers of subimages ..............
!
      DO J_PARTITION_INDEX = 1, N_PARTITION_Y
         J0 = I_Y(J_PARTITION_INDEX)+1
         J1 = I_Y(J_PARTITION_INDEX+1)


      DO I_PARTITION_INDEX = 1, N_PARTITION_X
         I0 = I_X(I_PARTITION_INDEX)+1
         I1 = I_X(I_PARTITION_INDEX+1)


      N_SIZE_SUB_DATA = (I1 - I0 + 1) * ( J1 - J0 + 1)

      ALLOCATE(SUB_DATA_138(N_SIZE_SUB_DATA), &
               SUB_DATA_066_OR_124(N_SIZE_SUB_DATA))

        L1 = 1
      DO I = I0, I1
        L2 = L1 + J1 - J0
       SUB_DATA_138(L1:L2) = DATA_FLOAT_138(I,J0:J1)
       SUB_DATA_066_OR_124(L1:L2) = DATA_FLOAT_066_OR_124(I,J0:J1)
       L1 = L2 + 1
      END DO


        CALL REMOVE_SURFACE_EFFECT(INDEX_QA_SLOPE, &
                                   SUB_DATA_138, &
           SUB_DATA_066_OR_124, N_SLOPE,                  &
           N_FOR_BAR_AVERAGE, REFL_138_LOW_PERCENTAGE,      &
                              REFL_138_HIGH_PERCENTAGE,     &
           REFL_138_L_REFLECTION, REFL_138_H_REFLECTION,    & 
           OUT_SLOPE)

      DEALLOCATE(SUB_DATA_138,SUB_DATA_066_OR_124) 

      IF(INDEX_QA_SLOPE == 1) THEN    

!            ......   good slope 

      SLOPE_AT_CENTER(I_PARTITION_INDEX,J_PARTITION_INDEX) = &
            OUT_SLOPE

      ELSE

!**************** temporary output ********************************
!
!        print *, 'sub_images (NX,NX) that needs  look-up table: ( ', &
!               I_PARTITION_INDEX, J_PARTITION_INDEX, ')', &
!                ',  slope(original)=', OUT_SLOPE
!
!
!          ..... use default slope obtained from look-up table

      ANG_SOLAR = &
         DATA_SOLAR_ZENITH(I_PARTITION_INDEX,J_PARTITION_INDEX)

           IF(ANG_SOLAR > 88.0) ANG_SOLAR = 88.0
      ANG_VIEW = &
         DATA_VIEW_ZENITH(I_PARTITION_INDEX,J_PARTITION_INDEX)
           IF(ANG_VIEW > 88.0) ANG_VIEW = 88.0
     
      GEOM_FACTOR = 1.0 / COS(0.0174533 * ANG_SOLAR) &
                  + 1.0 / COS(0.0174533 * ANG_VIEW)

      IF(GEOM_FACTOR < 6.0) THEN
        SLOPE_DEFAULT = 0.8325104 &
           * EXP( - 5.719236E-2 * GEOM_FACTOR) 
      ELSE
        SLOPE_DEFAULT = 0.7378554 &
           * EXP( - 3.705750E-2 * GEOM_FACTOR) 
      END IF

      IF(SLOPE_DEFAULT < 0.25 ) SLOPE_DEFAULT = 0.25  ! extreme case (e.g.,
!                     nearly horizontal view when parallel plane model for 
!                     transmittance breaks down)
  
      SLOPE_AT_CENTER(I_PARTITION_INDEX,J_PARTITION_INDEX) = &
             SLOPE_DEFAULT 
      END IF

      END DO
      END DO

!                     *** get the slope at the corners of subimages ****

      SLOPE_AT_NODE(1,1) = SLOPE_AT_CENTER(1,1)
      SLOPE_AT_NODE(N_PARTITION_X+1,1) = SLOPE_AT_CENTER(N_PARTITION_X,1)
      SLOPE_AT_NODE(N_PARTITION_X+1,N_PARTITION_Y+1) = &
                               SLOPE_AT_CENTER(N_PARTITION_X,N_PARTITION_Y)
      SLOPE_AT_NODE(1,N_PARTITION_Y+1) = SLOPE_AT_CENTER(1,N_PARTITION_X)

      
      DO I = 1, N_PARTITION_X - 1
      SLOPE_AT_NODE(I+1,1) = 0.5*(SLOPE_AT_CENTER(I,1) &
                                  + SLOPE_AT_CENTER(I+1,1))
      SLOPE_AT_NODE(I+1,N_PARTITION_Y+1) = 0.5*(SLOPE_AT_CENTER(I,N_PARTITION_Y) &
                                  + SLOPE_AT_CENTER(I+1,N_PARTITION_Y))
      END DO
   
      DO J = 1, N_PARTITION_Y - 1
      SLOPE_AT_NODE(1,J+1) = 0.5*(SLOPE_AT_CENTER(1,J) &
                                  + SLOPE_AT_CENTER(1,J+1))
      SLOPE_AT_NODE(N_PARTITION_X+1,J+1) = 0.5*(SLOPE_AT_CENTER(N_PARTITION_X,J) &
                                  + SLOPE_AT_CENTER(N_PARTITION_X,J+1))
      END DO

      DO I = 2, N_PARTITION_X
      DO J = 2, N_PARTITION_Y
      SLOPE_AT_NODE(I,J) = 0.25*(SLOPE_AT_CENTER(I-1,J-1) &
          + SLOPE_AT_CENTER(I,J-1) + SLOPE_AT_CENTER(I-1,J) &
          + SLOPE_AT_CENTER(I,J) )
           
      END DO
      END DO     




!**************************temperary output *********************************8
!
!         WRITE(*,"(1x, (10X, '*** slopes at centers (4x4) of subimage ***  '))")
!
!
!      DO J = N_PARTITION_Y,1,-1
!         WRITE(*,"(1x, 4(E12.6,2x))") SLOPE_AT_CENTER(1:N_PARTITION_X,J)
!      END DO
!         WRITE(*,"(1x, 4(E12.6,2x))")
!         WRITE(*,"(1x, 4(E12.6,2x))")
!
!         WRITE(*,"(1x, (10X, '*** slopes at nodes (5x5) of subimage *** '))")
!
!      DO J = N_PARTITION_Y+1,1,-1
!         WRITE(*,"(1x, 5(E12.6,2x))") SLOPE_AT_NODE(1:N_PARTITION_X+1,J)
!      END DO
!*****************************************************************************




! ................ handle bad data (just in case) ....................

       DO I = 1, N_X
       DO J = 1, N_Y

        X = DATA_FLOAT_066_OR_124(I,J)
        Y = DATA_FLOAT_138(I,J)

        IF( X < X_LOWER_BOUND .OR. Y < Y_LOWER_BOUND &
         .OR. X > X_UPPER_BOUND .OR. Y > Y_UPPER_BOUND) THEN

           DATA_FLOAT_066_OR_124(I,J) = - 1999.80

        END IF

       END DO
       END DO

! ................ .......................................................


   
      DO I = 1, N_X
          I1 = 1
          I2 = N_PARTITION_X + 1
 
          loop_x_check: DO
             K = ( I1 + I2 ) / 2
             IF( I < I_X(K)) THEN
                I2 = K
             ELSE
                I1 = K
             END IF
             IF(I2- I1 <= 1) EXIT loop_x_check
          END DO loop_x_check
          S_X = FLOAT(I- I_X(I1))/FLOAT(I_X(I2) - I_X(I1))

      loop_j_index_check: DO J = 1, N_Y

          IF( DATA_FLOAT_066_OR_124(I,J) < - 1999.80) &
                                 CYCLE loop_j_index_check
          J1 = 1
          J2 = N_PARTITION_Y + 1
 
          loop_y_check: DO
             K = ( J1 + J2 ) / 2
             IF( J < I_Y(K)) THEN
                J2 = K
             ELSE
                J1 = K
             END IF
             IF(J2- J1 <= 1) EXIT loop_y_check
          END DO loop_y_check

          S_Y = FLOAT(J- I_Y(J1))/FLOAT(I_Y(J2) - I_Y(J1))

          SLOPE_AT_PIXEL =  (1.0-S_X) * (1.0 - S_Y) * SLOPE_AT_NODE(I1,J1) &
             + (1.0-S_X) * S_Y * SLOPE_AT_NODE(I1,J2) &
             + S_X * ( 1.0 - S_Y) * SLOPE_AT_NODE(I2,J1) &
             + S_X * S_Y * SLOPE_AT_NODE(I2,J2)

        DERIVED_CIRRUS_REFLECTION(I,J) =  DATA_FLOAT_138(I,J) /  SLOPE_AT_PIXEL

        END DO loop_j_index_check
        END DO


      DO I = 1, N_X
      DO J = 1, N_Y
      IF( DATA_FLOAT_066_OR_124(I,J) > - 1999.80 .AND. &
         DERIVED_CIRRUS_REFLECTION(I,J) > &
          DATA_FLOAT_066_OR_124(I,J) ) &
      DERIVED_CIRRUS_REFLECTION(I,J) = DATA_FLOAT_066_OR_124(I,J)

!!! The following statement is commented out in view of the new initialization
!!!     statement earlier (DERIVED_CIRRUS_REFLECTION = -1999.80)
!!!       IF(DERIVED_CIRRUS_REFLECTION(I,J) < 0.0) &
!!!           DERIVED_CIRRUS_REFLECTION(I,J) = 0.0

      END DO
      END DO

        DERIVED_CIRRUS_REFLECTION = 1.0E-3 * DERIVED_CIRRUS_REFLECTION

      END  SUBROUTINE COMPUTE_CIRRUS_REFLECTANCE



!---------------------------------------------------------------------------
        SUBROUTINE REMOVE_SURFACE_EFFECT(INDEX_QA_SLOPE, &
                                         DATA_FLOAT_138, &
           DATA_FLOAT_066_OR_124, N_SLOPE,                  &
           N_FOR_BAR_AVERAGE, REFL_138_LOW_PERCENTAGE,      &
                              REFL_138_HIGH_PERCENTAGE,     &
           REFL_138_L_REFLECTION, REFL_138_H_REFLECTION,    & 
           OUT_SLOPE )

!
! DESCRIPTION: This subroutine is first to select the useful data of 1.38 and 
!   0.66 (or 1.24) um reflectance. Further, it removes the surface effect in 
!   the co-relation plot of the two channels.  
!
        INTRINSIC SELECTED_INT_KIND, KIND, SIZE, SUM 

        INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
        INTEGER, PARAMETER :: SP = KIND(1.0)
        REAL(SP), PARAMETER :: X_REFL_H_REFLECTION = 1000.0 
        REAL(SP), PARAMETER :: X_REFL_L_REFLECTION = 0.0
                               
!           X_REFL_H_REFLECTION is the maximum visible reflectance
!           X_REFL_L_REFLECTION is the minimum visible reflectance

        INTEGER(I4B), INTENT(IN) :: N_SLOPE
        INTEGER(I4B), INTENT(OUT) :: INDEX_QA_SLOPE
!          N_SLOPE is the number of binned slices for 1.38 um reflectance

        INTEGER(I4B), INTENT(IN) :: N_FOR_BAR_AVERAGE ! the number of points 
                                            ! in each slice to get average values

        REAL(SP), INTENT(IN) :: REFL_138_LOW_PERCENTAGE,  &
                                REFL_138_HIGH_PERCENTAGE, &
                                REFL_138_L_REFLECTION,    &
                                REFL_138_H_REFLECTION

        REAL(SP), INTENT(OUT) :: OUT_SLOPE

!           These are the output of this subroutine. OUT_X_RAYLEIGH is the 
!       x-axis interception of the lower part of the derived 0.66-1.38 curve.
!       OUT_Y_DIVISION determines the division of the two parts. 
     

        INTEGER(I4B), DIMENSION(N_SLOPE+1) :: N_SUM
        REAL(SP), DIMENSION(N_SLOPE+1) :: Y_138_NODE

        REAL(SP), DIMENSION(N_SLOPE) :: &
                          X_066_OR_124_BAR, Y_138_BAR


        REAL(SP), DIMENSION(:), INTENT(INOUT) :: DATA_FLOAT_138, &
                     DATA_FLOAT_066_OR_124


        REAL(SP), DIMENSION(:), ALLOCATABLE ::  SAMPLE_138, &
                                                SAMPLE_066_OR_124

        CHARACTER(2) :: FILE_NAME_NUMBER_X,  FILE_NAME_NUMBER_Y
!......................................................................
        INTERFACE

      FUNCTION locate(xx,x)
        INTRINSIC SELECTED_INT_KIND, KIND, SIZE

        INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
        INTEGER, PARAMETER :: SP = KIND(1.0)

      REAL(SP), INTENT(IN) :: x
      REAL(SP), DIMENSION(:), INTENT(IN) :: xx

      INTEGER(I4B) :: n, jl,jm,ju
      INTEGER(I4B) :: locate
     
      LOGICAL :: ascnd
      END FUNCTION locate


        SUBROUTINE sort2(arr,brr)
        INTRINSIC SELECTED_INT_KIND, KIND, SIZE
        INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
        INTEGER, PARAMETER :: SP = KIND(1.0)

        REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr, brr
        INTEGER(I4B), PARAMETER :: M = 15, NSTACK = 200
        INTEGER(I4B) :: n, i, ir, j, jstack, k, l
        INTEGER(I4B), DIMENSION(NSTACK) :: istack
        REAL(SP) ::  a, b, temp
        END SUBROUTINE sort2

      END INTERFACE 
!......................................................................

!- Step 1: Sort the image array DATA_FLOAT_138, then get DATA_FLOAT_138
!         in ascending order and re-arranged array DATA_FLOAT_066_OR_124,
!         accordingly.

!         note: DATA_FLOAT_066_OR_124 indicates 0.66 um reflectance in
!         the case of overland, and indicates 1.24 um reflectance in
!         the case of overocean.
 


           N = SIZE(DATA_FLOAT_138)
 

	   CALL SORT2(DATA_FLOAT_066_OR_124, &
                      DATA_FLOAT_138 )

           IF(DATA_FLOAT_066_OR_124(N) < &
               X_REFL_H_REFLECTION ) THEN
             N_MAX_X = N
           ELSE
             N_MAX_X = LOCATE(DATA_FLOAT_066_OR_124, &
               X_REFL_H_REFLECTION )
           END IF

           IF(DATA_FLOAT_066_OR_124(1) >= &
               X_REFL_L_REFLECTION ) THEN
             N_MIN_X = 1
           ELSE
             N_MIN_X = LOCATE(DATA_FLOAT_066_OR_124, &
               X_REFL_L_REFLECTION ) + 1
           END IF


!                N_X_GOOD = N_MAX_X - N_MIN_X + 1
!             The number of good points checked with 0.66 or 1.24-um band


	   CALL SORT2(DATA_FLOAT_138(N_MIN_X:N_MAX_X), &
                 DATA_FLOAT_066_OR_124(N_MIN_X:N_MAX_X) )


        IF(DATA_FLOAT_138(N_MAX_X) < REFL_138_H_REFLECTION) THEN
            N_MAX_GOOD = N_MAX_X
        ELSE
           N_MAX_GOOD = LOCATE(DATA_FLOAT_138(N_MIN_X:N_MAX_X), &
                    REFL_138_H_REFLECTION) + N_MIN_X - 1
        END IF

!               Note that the return of locate( ) is 
!                                between 1 and (N_MAX_X - N_MIN_X)


        IF( DATA_FLOAT_138(N_MIN_X) >= REFL_138_L_REFLECTION) THEN
            N_MIN_GOOD = N_MIN_X
        ELSE
            N_MIN_GOOD = LOCATE(DATA_FLOAT_138(N_MIN_X:N_MAX_X), &
            REFL_138_L_REFLECTION) + N_MIN_X   
        END IF


!..... Cut off some percentages at both ends ............

	N_MIN_138 =   N_MIN_GOOD &
         + IFIX(  REFL_138_LOW_PERCENTAGE * FLOAT(N_MAX_GOOD - N_MIN_GOOD))

	N_MAX_138 =   N_MIN_GOOD &
         + IFIX(  REFL_138_HIGH_PERCENTAGE * FLOAT(N_MAX_GOOD - N_MIN_GOOD))


!****************************************************************************

        N_SAMPLE = N_MAX_138 - N_MIN_138 + 1
        NON_ZERO_BIN = 0
       
        IF( N_SAMPLE >= 600) THEN 
!           Note that we think the results are statistically meaningless if
!           the good points are less than 600
        
        ALLOCATE(SAMPLE_138(N_SAMPLE), &
                 SAMPLE_066_OR_124(N_SAMPLE) )

        SAMPLE_138(1:N_SAMPLE) = DATA_FLOAT_138(N_MIN_138:N_MAX_138)
        SAMPLE_066_OR_124(1:N_SAMPLE) = &
            DATA_FLOAT_066_OR_124(N_MIN_138:N_MAX_138)

!- Step 2: bin 1.38-micron image array DATA_FLOAT_138 btween 
!         REFL_TH_MIN_138 and REFL_TH_MAX_138. The node for the bins
!         are Y_138_NODE. The number of point in the bins are given
!         by N_SUM

! The bin-width is not uniform. Currently, we use two values.

        N_SLOPE_REGION_1 = IFIX((3.0/7.0) * FLOAT(N_SLOPE)) 

        RATIO_STEP1_TO_STEP2 = 2.0 / 3.0     !  ratio for the
                       ! steps for the two regions (STEP_1/STEP_2)

        STEP_1 = RATIO_STEP1_TO_STEP2 &
           * (SAMPLE_138(N_SAMPLE) - SAMPLE_138(1)) & 
            / ( FLOAT(N_SLOPE) &
            - (1.0 - RATIO_STEP1_TO_STEP2) * N_SLOPE_REGION_1)

        STEP_2 =( (SAMPLE_138(N_SAMPLE) - SAMPLE_138(1)) &
             - FLOAT(N_SLOPE_REGION_1) * STEP_1 ) &
               / FLOAT(N_SLOPE - N_SLOPE_REGION_1)
      

        DO I = 1, N_SLOPE_REGION_1 + 1
          Y_138_NODE(I) = SAMPLE_138(1) + FLOAT(I-1) * STEP_1
	END DO

        DO I = N_SLOPE_REGION_1 + 2, N_SLOPE + 1
          Y_138_NODE(I) = Y_138_NODE(N_SLOPE_REGION_1 + 1) &
             + FLOAT(I- N_SLOPE_REGION_1 - 1) * STEP_2
	END DO


        N_SUM = 0  
	DO I = 1, N_SAMPLE
           J = LOCATE(Y_138_NODE, SAMPLE_138(I))
           N_SUM(J+1)  = N_SUM(J+1) + 1 
	END DO

        DO I = 3, N_SLOPE 
         N_SUM(I) = N_SUM(I) + N_SUM(I-1) 
        END DO

         N_SUM(N_SLOPE + 1) = N_SAMPLE 


!- Step 3: Get Y_138_BAR and X_066_OR_124_BAR for each bin

        slice_loop: DO  I = 2, N_SLOPE + 1

         IF(N_SUM(I)-N_SUM(I-1) <= N_FOR_BAR_AVERAGE+3) &
                                              CYCLE slice_loop

           NON_ZERO_BIN = NON_ZERO_BIN + 1

           N_LOWER = N_SUM(I-1) + 1
           N_UPPER = N_SUM(I)
 

          CALL SORT2(SAMPLE_066_OR_124(N_LOWER:N_UPPER), &
                    SAMPLE_138(N_LOWER:N_UPPER) )

          N_LOWER_2 = N_LOWER 
          N_UPPER_2 = N_LOWER_2 + N_FOR_BAR_AVERAGE - 1


          N_POINT_TEST = 15000          ! this number will be in header
          ERR_FUNC_0 = 1.0E20

          SLICE_X_MIN = MINVAL(SAMPLE_066_OR_124(N_LOWER_2:N_UPPER_2))
          SLICE_X_MAX = MAXVAL(SAMPLE_066_OR_124(N_LOWER_2:N_UPPER_2))

          SLICE_Y_MIN = MINVAL(SAMPLE_138(N_LOWER_2:N_UPPER_2))
          SLICE_Y_MAX = MAXVAL(SAMPLE_138(N_LOWER_2:N_UPPER_2))

           IY = 67107

          loop_points_test : DO I_POINT_TEST = 1, N_POINT_TEST


            IY = 125 * IY
            IY = IY - IY / 2796203 * 2796203
            HARVEST = float(IY)
            HARVEST1 = HARVEST / 2796203.0

            IY = 125 * IY
            IY = IY - IY / 2796203 * 2796203
            HARVEST = float(IY)
            HARVEST2 = HARVEST / 2796203.0


           SLICE_X = SLICE_X_MIN &
             + HARVEST1 * ( SLICE_X_MAX - SLICE_X_MIN )          

           SLICE_Y = SLICE_Y_MIN &
             + HARVEST2 * ( SLICE_Y_MAX - SLICE_Y_MIN )          

           ERR_FUNC_1 = SUM(                             &
            + ABS(SLICE_X                                    &
               - SAMPLE_066_OR_124(N_LOWER_2:N_UPPER_2))  &
            + ABS(  SLICE_Y                                  &
               - SAMPLE_138(N_LOWER_2:N_UPPER_2))  ) 

          IF(ERR_FUNC_1 > ERR_FUNC_0 ) CYCLE loop_points_test
            
             ERR_FUNC_0 = ERR_FUNC_1
          
          X_066_OR_124_BAR(NON_ZERO_BIN) =  SLICE_X

          Y_138_BAR(NON_ZERO_BIN) =  SLICE_Y

          END DO  loop_points_test
 
         END DO slice_loop       ! end the loop for a slice

        DEALLOCATE(SAMPLE_138,SAMPLE_066_OR_124)

        END IF    ! end slicing calculation



!-step 4: 
!        find the left envelope that is free of cirrus 
!     contamination

      IF(NON_ZERO_BIN < 2) THEN

        OUT_SLOPE = 0.0
   
      ELSE

      x0 = MINVAL(X_066_OR_124_BAR(1:NON_ZERO_BIN))
      y0 = MINVAL(Y_138_BAR(1:NON_ZERO_BIN))

      x_max = MAXVAL(X_066_OR_124_BAR(1:NON_ZERO_BIN))
      y_max = MAXVAL(Y_138_BAR(1:NON_ZERO_BIN))


      N_TEST = 250000


       err_func_line = 1.0e20

       IY = 67107

      monte_carlo_1: DO I_TEST = 1, N_TEST


       IY = 125 * IY
       IY = IY - IY / 2796203 * 2796203
       HARVEST = float(IY)
       HARVEST1 = HARVEST / 2796203.0

       IY = 125 * IY
       IY = IY - IY / 2796203 * 2796203
       HARVEST = float(IY)
       HARVEST2 = HARVEST / 2796203.0


       x1 = x0 + HARVEST1 *   (x_max - x0 ) 
       y1 = y0 + HARVEST2 * (y_max - y0 ) 

       IF(IY == 0) IY = 67107
       IY = 125 * IY
       IY = IY - IY / 2796203 * 2796203
       HARVEST = float(IY)
       a = HARVEST / 2796203.0

        b = y1 - a * x1

            err_test = SUM(      &
                abs( a * X_066_OR_124_BAR(1:NON_ZERO_BIN) + b &
                - Y_138_BAR(1:NON_ZERO_BIN)  ) )

           if( err_test > err_func_line) CYCLE  monte_carlo_1

                 err_func_line = err_test

                 good_a = a


           end do  monte_carlo_1

        OUT_SLOPE = good_a

         END IF

            IF( OUT_SLOPE < 0.3 ) THEN
                 INDEX_QA_SLOPE = 0
            ELSE
                 INDEX_QA_SLOPE = 1
            END IF


        END SUBROUTINE REMOVE_SURFACE_EFFECT


!
!---------------------------------------------------------------------------
!
      FUNCTION locate(xx,x)
!
!---------------------------------------------------------------------------
        INTRINSIC SELECTED_INT_KIND, KIND, SIZE

        INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
        INTEGER, PARAMETER :: SP = KIND(1.0)

      REAL(SP), INTENT(IN) :: x
      REAL(SP), DIMENSION(:), INTENT(IN) :: xx

      INTEGER(I4B) :: n, jl,jm,ju
      INTEGER(I4B) :: locate
     
      LOGICAL :: ascnd

!
!       Given an array xx(1:N), and given a value x, returns
!      a value j such that x is between xx(j) and xx(j+1). 
!      xx must be monotonic, either increasing or decreasing.
!      j=0 or j=N is retrurned to indicate that x is out of range.

       n = SIZE(xx)
       ascnd = ( xx(n) >= xx(1) )
      jl=1
      ju=n
      do 
      if(ju-jl <= 1) exit
        jm=(ju+jl)/2
        if(ascnd .eqv. (x >= xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      end do

        locate=jl
      END FUNCTION locate


!------------------------------------------------------------------------------
!
	SUBROUTINE sort2(arr,brr)
!
!-------------------------------------------------------------------------------
        INTRINSIC SELECTED_INT_KIND, KIND, SIZE
        INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
        INTEGER, PARAMETER :: SP = KIND(1.0)

	REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr, brr
	INTEGER(I4B), PARAMETER :: M = 15, NSTACK = 200
	INTEGER(I4B) :: n, i, ir, j, jstack, k, l
	INTEGER(I4B), DIMENSION(NSTACK) :: istack
	REAL(SP) ::  a, b, temp

!	Sorts an array arr(1:n) into ascending order using Quicksort, while 
!	making the corresponding rearrangement of the array brr(1:n).
!       Both arr and brr are replaced on output by their
!       sorted rearrangements.

!       Parameters: M is the size of subarrys sorted by straight
!       insertion and NSTACK is the required auxiliary storge.


        n=SIZE(arr)
	jstack = 0
	l = 1
	ir = n

        do 

       	if (ir-l < M) then
		do  j = l+1,ir
			a = arr(j)
			b = brr(j)
			do  i = j-1, 1, -1
				if (arr(i) <= a) exit
				arr(i+1) = arr(i)
				brr(i+1) = brr(i)
			end do 
			arr(i+1) = a
			brr(i+1) = b
		end do 
		if (jstack == 0) return
		ir     = istack(jstack)
		l      = istack(jstack-1)
		jstack = jstack-2
	else
		k = (l+ir)/2

! swap the data

		temp     = arr(k)
		arr(k)   = arr(l+1)
		arr(l+1) = temp

		temp     = brr(k)
		brr(k)   = brr(l+1)
		brr(l+1) = temp

		if(arr(l) > arr(ir)) then
			temp    = arr(l)
			arr(l)  = arr(ir)
			arr(ir) = temp

			temp    = brr(l)
			brr(l)  = brr(ir)
			brr(ir) = temp
		endif

		if (arr(l+1) > arr(ir)) then
			temp     = arr(l+1)
			arr(l+1) = arr(ir)
			arr(ir)  = temp

			temp     = brr(l+1)
			brr(l+1) = brr(ir)
			brr(ir)  = temp
		endif

		if (arr(l) > arr(l+1)) then
			temp     = arr(l+1)    
			arr(l+1) = arr(l)
			arr(l)   = temp

			temp     = brr(l+1)
			brr(l+1) = brr(l)
			brr(l)   = temp
		endif

		i = l+1
		j = ir
		a = arr(l+1)
		b = brr(l+1)

           do             
          
             do
			i = i+1
		if (arr(i) >= a) exit
             end do

             do
			j = j-1
		if (arr(j) <= a) exit
              end do 

		if (j < i) exit
		temp   = arr(i)
		arr(i) = arr(j)
		arr(j) = temp

		temp   = brr(i)
		brr(i) = brr(j)
		brr(j) = temp
                end do

		arr(l+1) = arr(j)
		arr(j) = a
		brr(l+1) = brr(j)
		brr(j) = b
		jstack = jstack+2
!	write(*,*) 'jstack=',jstack,'NSTACK=',NSTACK
		if(jstack > NSTACK) pause 'NSTACK too small in sort2'
		if(ir-i+1 >= j-1) then
			istack(jstack)   = ir
			istack(jstack-1) = i
			ir               = j-1
		else
			istack(jstack)   = j-1
			istack(jstack-1) = l
			l                = i
		endif	
	endif
	end do

	END SUBROUTINE SORT2


