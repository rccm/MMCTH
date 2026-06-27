      subroutine read_data(error_file,xnvalm5)
C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C
C !INPUT PARAMETERS:
C
C    Type       Name             Description
C    ====       ====             ===========
C    REAL*4     xnvalm5          MODIS reflectivities to be converted
C                                 to n-values
C
C !OUTPUT PARAMETERS:
C
C    Type       Name             Description
C    ====       ====             ===========
C    LOGICAL*4  error_file       returned status flag
C
C !REVISION HISTORY:
C
C    Initial Version by Jeremy Warner   12/01/2006
C     Modified from legacy code by Christina Hsu.
C
C !TEAM-UNIQUE HEADER:
C
C    This software is developed by the Deep Blue Science Team
C    for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-02041.
C
C !REFERENCES AND CREDITS
C
C !DESIGN NOTES:
C
C   Externals:
C
C     xnvalm                     (sample.inc)
C     xlat                       (sample.inc)
C     xlong                      (sample.inc)
C     sfref412                   (sample.inc)
C     sfref470                   (sample.inc)
C     sfref650                   (sample.inc)
C
C   Functions:  none
C
C !END
C-----------------------------------------------------------------------
c
      use seawifs_surface_pressure, only: get_surface_pressure, 
     1                                    get_elevation
      
      implicit none
      logical error_file
      integer*4 i,j,k, ilat5, ilon5, itmp, ilat4, ilon4, status
      real*4    xnvalm5(3), psi, cc, s1, s2, ss
      real*4    rr, xmu, xtmp
      include 'sample.inc'
      include 'aottbl.inc'
 
c     -- first undo sza normalization in the L1B code
 
      xmu = cos(sza * 3.14159/180.)
 
      do i = 1,3
      rr   = xnvalm5(i) * xmu/3.14159 ! use this for data from Wei
      if (rr .gt. 0.0) then
      xnvalm5(i) =  -100.*alog10(rr)
      	else
      	  print *, 'err - xvalm5(i) le 0, i= ',i
      	  go to 300
      endif
      enddo
 
      do i = 1,2
      xnvalm(i) = xnvalm5(1)
      enddo
      do i = 3,4
      xnvalm(i) = xnvalm5(2)
      enddo
 
      xnvalm(5) = xnvalm5(3)
      xnvalm(6) = xnvalm5(1)

c     -- find terrain pressure
      ilat5 = (-1.* xlat + 90.) *2.
      ilat5 = ilat5 + 1
      if (ilat5.gt.360) ilat5 = 360
      if (ilat5.lt.1)   ilat5 = 1

      if (xlong.ge.0.0) then
      ilon5 = xlong *2.
      ilon5 = ilon5 + 1
        else
         ilon5 = (xlong + 360.) *2.
         ilon5 = ilon5 + 1
      endif
      if (ilon5.gt.720) ilon5 = 720
      if (ilon5.lt.1)   ilon5 = 1

c      pteran = sfcprs(ilat5,ilon5)/1013.
      
      ilat4 = (-1.* xlat + 90.) *12.
      ilat4 = ilat4 + 1
      if (ilat4.gt.2160) ilat4 = 2160
      if (ilat4.lt.1)   ilat4 = 1

      ilon4 = (xlong + 180.) *12.
      ilon4 = ilon4 + 1
      if (ilon4.gt.4320) ilon4 = 4320
      if (ilon4.lt.1)   ilon4 = 1

      pteran = get_surface_pressure(xlat, xlong, status)
      pelev = get_elevation(xlat, xlong, status)

      isnow = 0
      do i = 1,6
         xnvalm(i) = xnvalm(i)/100.
      enddo
      error_file = .false.
      return

300   continue
      error_file = .true.
      return
      end
      
