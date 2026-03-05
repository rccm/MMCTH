!-----------------------------------------------------------------------
!
      SUBROUTINE Find_Contrail_Pixels(Xdim,Ydim,Image_Brightness,Data_Mask)
!
!-----------------------------------------------------------------------
! !F90
!
! !DESCRIPTION: This subroutine derives an index of aircraft
!              contrails from the image of the 1.375-micron MODIS channel.
!
! !Input Parameters:
!   INTEGER   Image_Brightness  MODIS reflectance at channel 26
!
! !OUTPUT PARAMETERS:
!   BYTE      Data_Mask         Two dimensional (2-D) array for passing
!                               contrail index data. Index 0 refers to 
!                               non-contrail pixel, and index 1 refers to 
!                               contral pixel.
!
! !REVISION HISTORY:
! Modified by Liqun Ma   Feb.  1998
! Temporary outputs were moved out.
!
! !TEAM-UNIQUE HEADER: - NOT YET DEFINED
!
! !REFERENCES AND CREDITS
!    First version created by Dr. Bill Ridgway
!    (ridgway@climate.gsfc.nasa.gov)
!                             1997 11 03
!
! !DESIGN NOTES: Must be compiled with f90, free format.
!
!   Functions and Subroutines:
!      DrawLine
!      Hough
!
! !END--------------------------------------------------------------------

      IMPLICIT NONE

      INTEGER :: Xdim, Ydim                          ! Image dimensions (input)
      INTEGER, PARAMETER :: Rdim = 300               ! Hough image R dimension
      INTEGER, PARAMETER :: Tdim = 360               ! Hough image T dimension
      INTEGER, PARAMETER :: Hthresh = 80             ! hough threshold for edge detection
      INTEGER, PARAMETER :: Lthresh = 40             ! minimum line length in pixels
      INTEGER, PARAMETER :: Gthresh = 40             ! gradient threshold ridge detection

      REAL :: Image_Brightness(0:Xdim-1,0:Ydim-1)    ! image brightness (input)
      INTEGER :: Data_Mask(0:Xdim-1,0:Ydim-1)        ! contrail mask (output)

      INTEGER, DIMENSION (:,:), ALLOCATABLE :: Data_Sobel
      INTEGER, DIMENSION (:,:), ALLOCATABLE :: Temp_Lines
      INTEGER, DIMENSION (:,:), ALLOCATABLE :: Data_Hough
      INTEGER, DIMENSION (:,:), ALLOCATABLE :: Data_Hplus
      INTEGER, DIMENSION (:,:), ALLOCATABLE :: Temp_Hplus
      INTEGER, DIMENSION (:,:), ALLOCATABLE :: Image_Integer_Copy

      CHARACTER(1), DIMENSION (:,:), ALLOCATABLE :: c1

      INTEGER :: x,y,r,t,irepeat
      INTEGER :: Max_Hough,Max_Sobel

      ALLOCATE (Data_Sobel (0:Xdim-1,0:Ydim-1) )     ! integer Sobel gradient field
      ALLOCATE (Temp_Lines (0:Xdim-1,0:Ydim-1) )     ! integer temporary
      ALLOCATE (Data_Hough (0:Tdim-1,0:Rdim-1) )     ! integer raw Hough field
      ALLOCATE (Data_Hplus (0:Tdim-1,0:Rdim-1) )     ! integer raw Hough field
      ALLOCATE (Temp_Hplus (0:Tdim-1,0:Rdim-1) )     ! integer raw Hough field
      ALLOCATE (Image_Integer_Copy(0:Xdim-1,0:Ydim-1))  ! Work array, array destroyed by Hough
      ALLOCATE (c1 (0:Xdim-1,0:Ydim-1))              ! byte image array

! read original image

      Data_Sobel         = 0
      Temp_Lines         = 0
      Data_Hough         = 0
      Data_Hplus         = 0
      Temp_Hplus         = 0
      Image_Integer_Copy = Image_Brightness

! write binary_original image
! commented out by Liqun Ma 02/17/98

!     OPEN(1,file="binary_original.out5",form='unformatted',access='direct',recl=Xdim*Ydim)
!     WRITE(1,rec=1) char(Image_Integer_Copy) ; CLOSE(1)

      CALL Hough(Lthresh, Gthresh, Hthresh, &
                 Ydim, Xdim, Rdim, Tdim, Image_Integer_Copy, &
                 Data_Sobel, Data_Mask, Data_Hplus, Data_Hough)

! make Data_Mask fatter

      DO irepeat=1,1
      Temp_Lines = Data_Mask
      DO x=1,Xdim-2
      DO y=1,Ydim-2
         Data_Mask(x,y) = MAX(Data_Mask(x,y),Temp_Lines(x-1,y))
         Data_Mask(x,y) = MAX(Data_Mask(x,y),Temp_Lines(x+1,y))
         Data_Mask(x,y) = MAX(Data_Mask(x,y),Temp_Lines(x,y-1))
         Data_Mask(x,y) = MAX(Data_Mask(x,y),Temp_Lines(x,y+1))
      END DO
      END DO
      END DO

! make Data_Hplus fatter (points that exceed Hough threshold)

      DO irepeat=1,3
      Temp_Hplus = Data_Hplus
      DO t=1,Tdim-2
      DO r=1,Rdim-2
         Data_Hplus(t,r) = MAX(Data_Hplus(t,r),Temp_Hplus(t-1,r))
         Data_Hplus(t,r) = MAX(Data_Hplus(t,r),Temp_Hplus(t+1,r))
         Data_Hplus(t,r) = MAX(Data_Hplus(t,r),Temp_Hplus(t,r-1))
         Data_Hplus(t,r) = MAX(Data_Hplus(t,r),Temp_Hplus(t,r+1))
      END DO
      END DO
      END DO

      IF(.TRUE.) THEN       ! THIS IF/THEN BLOCK IS FOR DIAGNOSTICS ONLY

! write Sobel pixel image

! cmnmented out by Liqun Ma 02/17/98
!
!     OPEN(1,file="binary_sobel.out6",form='unformatted',access='direct',recl=Xdim*Ydim)
!     Max_Sobel=MAXVAL(Data_Sobel)
!     WRITE(1,rec=1) char(Data_Sobel*255/Max_Sobel) ; CLOSE(1)

! write Hplus pixel image


!     OPEN(1,file="binary_hplus.out7",form='unformatted',access='direct',recl=Tdim*Rdim)
!     WRITE(1,rec=1) char(255*Data_Hplus) ; CLOSE(1)

! write hough image

!     OPEN(1,file="binary_hough.out8",form='unformatted',access='direct',recl=Tdim*Rdim)
!     Max_Hough=MAXVAL(Data_Hough)
!     WRITE(1,rec=1) char(Data_Hough*255/Max_Hough) ; CLOSE(1)

! write overlay image

!     OPEN(1,file="binary_overlay.out9",form='unformatted',access='direct',recl=Xdim*Ydim)
!     WRITE(1,rec=1) char(max(int(Image_Brightness),255*Data_Mask)) ; CLOSE(1)

      ENDIF

      DEALLOCATE (Data_Sobel) 
      DEALLOCATE (Temp_Lines)
      DEALLOCATE (Data_Hough)   
      DEALLOCATE (Data_Hplus) 
      DEALLOCATE (Temp_Hplus) 
      DEALLOCATE (Image_Integer_Copy)

      END SUBROUTINE Find_Contrail_Pixels

      SUBROUTINE Hough(Lthresh, Gthresh, Hthresh, &
                       Ydim, Xdim, Rdim, Tdim, Data_Image, &
                       Data_Sobel, Data_Mask, Data_Hplus, Data_Hough)

!     Version 971110 (November 10, 1997)
!     Questions: Bill.Ridgway@gsfc.nasa.gov
!     Must be compiled with f90, free format

      IMPLICIT NONE

      INTEGER :: Rdim                ! Size of the output image R dimension
      INTEGER :: Tdim                ! Size of the output image T dimension

      INTEGER :: Ydim, Xdim          ! Image dimensions

      INTEGER :: x,y,r,t             ! Image indices
      REAL :: Radius, Theta          ! Hough transform values
      REAL :: Max_Radius             ! Maximum Radius range
      REAL :: Max_Theta              ! Maximum Theta range
      INTEGER :: dr,dt               ! Hough offset indices
      INTEGER :: Max_Flag            ! Boolean flag

      INTEGER :: Data_Image(0:Xdim-1,0:Ydim-1)
      INTEGER :: Data_Sobel(0:Xdim-1,0:Ydim-1)       ! integer sobel gradient field
      INTEGER :: Data_Mask(0:Xdim-1,0:Ydim-1)
      INTEGER :: Data_Hough(0:Tdim-1,0:Rdim-1)       ! integer raw Hough field
      INTEGER :: Data_Hplus(0:Tdim-1,0:Rdim-1)       ! integer raw Hough field

      INTEGER :: Max_Hough

      REAL :: table_sin(0:Tdim-1)    ! Lookup tables for sin and cos
      REAL :: table_cos(0:Tdim-1)    ! Lookup tables for sin and cos

      INTEGER :: last                ! Last value of r

      INTEGER :: Lthresh             ! Minimum line length to draw
      INTEGER :: Gthresh             ! Gradient threshold, edge detection
      INTEGER :: Hthresh             ! Hough threshold for edge detection

      INTEGER :: j,k

      REAL :: pi
      pi = acos(-1.)

! Calculate range of Radius and Theta
      Max_Radius = sqrt(REAL(Ydim)*REAL(Ydim)+REAL(Xdim)*REAL(Xdim))
      Max_Theta = 2*pi

! Initialize SIN and COS lookup tables from 0 to 2 PI
      DO t=0,Tdim-1
        Theta = (REAL(t)-REAL(Tdim)/2.) * Max_Theta / REAL(Tdim)
        table_sin(t) = sin(Theta)
        table_cos(t) = cos(Theta)
      END DO

! Process input image to calculate Sobel gradient

      Data_Sobel = 0

      DO y=24,Ydim-25
      DO x=24,Xdim-25
         Data_Sobel(x,y) = 0
         DO j = 0,24,6
         DO k = 0,24,6
              IF(2*Data_Image(x,y).gt.3*Data_Image(x+j,y+k).and. &
                 2*Data_Image(x,y).gt.3*Data_Image(x-j,y-k).and. &
                 Data_Image(x,y).gt.20+Data_Image(x+j,y+k).and. &
                 Data_Image(x,y).gt.20+Data_Image(x-j,y-k)) then
                 Data_Sobel(x,y) = Data_Sobel(x,y) + 1
              ENDIF
         end do
         end do
      END DO
      END DO

! rhucek 06/29/98: replaced 1 line on request from SCF
!     Data_Sobel = 100*Data_Sobel/MAXVAL(Data_Sobel)
      if (MAXVAL(Data_Sobel) .ne. 0) Data_Sobel = 100*Data_Sobel/MAXVAL(Data_Sobel)

      if(.true.) then

! Data_Mask used temporarily to mark pixels exceeding threshold
      where (Data_Sobel > Gthresh)
        Data_Mask = 1
      elsewhere
        Data_Mask = 0
      end where

! Record line of (r,t) values for strong edge pixels in Hough array
      DO y=1,Ydim-2
      DO x=1,Xdim-2
      IF (Data_Sobel(x,y).gt.Gthresh) THEN

! Draw first point on Hough image
         t = 1
         Radius = (REAL(y) - REAL(Ydim)/2.)*table_cos(t) &
                + (REAL(x) - REAL(Xdim)/2.)*table_sin(t)
         r = 0.5+(Radius/Max_Radius + 0.5)*Rdim
         IF ((r.ge.0).and.(r.lt.Rdim)) THEN
            Data_Hough(t,r)=Data_Hough(t,r)+1
            Data_Hough(t+Tdim/2,Rdim-r)=Data_Hough(t+Tdim/2,Rdim-r)+1
         ENDIF
         last = r

! Draw remaining points on Hough image
      DO t=2,Tdim/2
            Radius = (REAL(y) - REAL(Ydim)/2)*table_cos(t) &
                   + (REAL(x) - REAL(Xdim)/2)*table_sin(t)
            r = 0.5+(Radius/Max_Radius + 0.5)*Rdim
            IF ((r.ge.0).and.(r.lt.Rdim)) THEN
! Draw next point on Hough image
               Data_Hough(t,r)=Data_Hough(t,r)+1
               Data_Hough(t+Tdim/2,Rdim-r)=Data_Hough(t+Tdim/2,Rdim-r)+1

   1  continue
      IF((last-1.gt.r).and.(last.lt.Rdim)) THEN
        last = last - 1
        Data_Hough(t,last) = Data_Hough(t,last) + 1
        Data_Hough(t+Tdim/2,Rdim-last) = Data_Hough(t+Tdim/2,Rdim-last) + 1
        GOTO 1
      ENDIF

! Fill in any missing points on Hough image

   2  continue
      IF((last+1.lt.r).and.(last.gt.0)) THEN
        last = last + 1
        Data_Hough(t,last) = Data_Hough(t,last) + 1
        Data_Hough(t+Tdim/2,Rdim-last) = Data_Hough(t+Tdim/2,Rdim-last) + 1
        GOTO 2
      ENDIF

      END IF

      last = r
      END DO

      END IF
      END DO
      END DO

! Convert original Data_Image to pixels exceeding gradient threshold
! Define boundary frame

      Data_Image(:,0)      = 0
      Data_Image(:,Ydim-1) = 0
      Data_Image(0,:)      = 0
      Data_Image(Xdim-1,:) = 0

! Make adjacent pixels non-zero

      DO y=1,Ydim-2
      DO x=1,Xdim-2
         Data_Image(x,y) =  Data_Mask(x+1,y+1) + Data_Mask(x,y+1) + Data_Mask(x-1,y+1) &
                         +  Data_Mask(x+1,y)   + Data_Mask(x,y)   + Data_Mask(x-1,y)   &
                         +  Data_Mask(x+1,y-1) + Data_Mask(x,y-1) + Data_Mask(x-1,y-1)
      END DO
      END DO

! Clear the output image

      Data_Mask(:,:) = 0

! Process Hough array to find strong maxima

      Max_Hough=MAXVAL(Data_Hough)

      DO t=0,Tdim/2-1
      DO r=0,Rdim-1
      IF (Data_Hough(t,r).gt.Max_Hough*Hthresh/100) THEN

! Check 5x5 neighborhood for higher values, set flag
      Max_Flag = 1
      DO dt=MAX(0,t-2),min(t+2,Tdim-1)
      DO dr=MAX(0,r-2),min(r+2,Rdim-1)
         IF(Data_Hough(dt,dr).gt.Data_Hough(t,r)) THEN
            Max_Flag = 0
            GOTO 3
         ENDIF
      END DO
      END DO
   3  CONTINUE

! Process line array to plot line segments
      IF (Max_Flag.eq.1) THEN
         Data_Hplus(t,r)=1
         Radius = (REAL(r)/REAL(Rdim) - 0.5)*Max_Radius
         Theta  = (REAL(t)/REAL(Tdim) - 0.5)*Max_Theta
         CALL DrawLine(Lthresh, Radius, Theta, Ydim, Xdim, &
                          Data_Image, Data_Mask)
      ENDIF
      ENDIF
      END DO
      END DO

      ENDIF

      END SUBROUTINE Hough

      SUBROUTINE DrawLine(Lthresh, radius, theta, Ydim, Xdim, &
                          Data_Image, Data_Mask)

      IMPLICIT NONE

      INTEGER :: Lthresh             ! Minimum line length to draw
      REAL :: radius,theta           ! Polar radius and angle
      INTEGER :: Ydim, Xdim          ! Image dimensions

      INTEGER :: Data_Image(0:Xdim-1,0:Ydim-1)
      INTEGER :: Data_Mask(0:Xdim-1,0:Ydim-1)

      REAL :: M,B                    ! Euclidean coordinate line parameters
      INTEGER :: line_start          ! Line start point
      INTEGER :: line_end            ! Line end point
      INTEGER :: x,y,x2,y2           ! Image indices

      REAL :: pi
      pi = acos(-1.)

! Handle lines with y changing faster than x

      IF (((theta.ge.PI/4).and.(theta.le.3*PI/4)) .or. &
          ((-theta.ge.PI/4).and.(-theta.le.3*PI/4))) THEN

        M = cos(theta)/sin(theta)
        B = radius / sin(theta) + (Ydim/2)*M + Xdim/2
        line_start = -1
        line_end = -1

        DO y = 0, Ydim-1

! Get line coordinates

          x = 0.5 + B - y*M

          IF ((x.ge.0).and.(x.lt.Xdim)) THEN

            IF ((Data_Image(x,y).gt.0).and.(line_start.lt.0)) THEN
               line_start = y
            ELSE IF ((Data_Image(x,y).eq.0).and.(line_start.gt.0)) THEN
               line_end = y  ! Draw only long lines
               IF (line_end-line_start.gt.Lthresh) THEN
                  DO y2 = line_start, line_end-1
                     x2 = (0.5 + B - y2*M)
                     Data_Mask(x2,y2) = 1
                  END DO
               ENDIF
               line_start = -1
               line_end = -1
            ENDIF

          ENDIF

        END DO

! Handle lines with x changing faster than y

      ELSE

        M = sin(theta)/cos(theta)
        B = radius / cos(theta) + (Xdim/2)*M + Ydim/2
        line_start = -1
        line_end = -1

        DO x = 0, Xdim-1

! Get line coordinates

          y = 0.5 + B - x*M

          IF ((y.ge.0).and.(y.lt.Ydim)) THEN

            IF ((Data_Image(x,y).gt.0).and.(line_start.lt.0)) THEN
               line_start = x
            ELSE IF ((Data_Image(x,y).eq.0).and.(line_start.gt.0)) THEN
               line_end = x  ! Draw only long lines
               IF (line_end-line_start.gt.Lthresh) THEN
                  DO x2 = line_start, line_end-1
                     y2 = (0.5 + B - x2*M)
                     Data_Mask(x2,y2) = 1
                  END DO
               ENDIF
               line_start = -1
               line_end = -1
            ENDIF

          ENDIF

        END DO

      ENDIF
      END SUBROUTINE DrawLine
