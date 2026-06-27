       subroutine check_reg_uniformity(line_edge,ele_edge,
     +                                 contx_eco,contx_topog,
     +                                 contx_snow,n_nadir,kele,
     +                                 eco_type,day,land,water,
     +                                 coast,snow,ice,uniform)

      implicit none
      save


c ... Common statement for debug purposes
      common / bug / debug, h_output

c--------------------------------------------------------------------
c!F77
c
c!DESCRIPTION:
c This entire routine checks is used to decide of the regional 
c processing box is consistent enough for running the uniformity check.
c The diffent backgroun variables are checked against the center
c pixel, and if they are not all the same, then the variable "uniform"
c is set to false.  All variables must be consistent for uniform to
c be set to true.
c
c!Input parameters:
c line_edge     Logical variable - true if processing border line
c ele_edge      Logical variable - true if processing border pixel
c contx_eco     Array containing context of lines of ecosystem values
c contx_topog   Array containing context of lines of land/sea values
c contx_snow    Array containing context of lines of snow values
c n_nadir       Nadir pixel number for lines in a context
c kele          Current processing element
c eco_type      Ecosystem type for current pixel
c day           Logical variable flagging day scenes
c land          Logical variable flagging land scenes
c water         Logical variable flagging water scenes
c coast         Logical variable flagging coast scenes
c snow          Logical variable flagging snow background scenes
c ice           Logical variable flagging ice background scenes
c
c!Output Parameters:
c uniform       Logical variable flagging uniform scenes
c               (Places where all pixels in context are similar)
c
c!Revision History:
c
c!Team-unique Header:
c
c!References and Credits:
c See Cloud Mask ATBD-MOD-06.
c
c!END
c----------------------------------------------------------------------

      include 'global.inc'
c
c     scalar arguments
      integer kele
      logical day,uniform,line_edge,ele_edge,
     +        coast,land,water,snow,ice
      byte eco_type

c     array arguments
      integer contx_topog(npixel,nlcntx),n_nadir(nlcntx),
     +        contx_snow(npixel,nlcntx)
      byte contx_eco(npixel,nlcntx)

c     local scalars
      integer i,j,imv,ide,itotal,nland,nwater,ncoast,h_output,debug,
     *        jj,kk

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(10x/,''Subroutine check_reg_uniformity '')')
        write(h_output,'(10x/,''Line Edge or Element Edge? '',2L5/)')
     +        line_edge, ele_edge
      endif
c ................................................................

c ... initialize variables
      imv = ((necntx - 1) / 2) + 1
      itotal = nlcntx * necntx
      nland = 0
      nwater = 0
      ncoast = 0


c ... Check all surface variables for consistency
c     If any of the checks fail, then set uniformity to zero

c ... First check if line or pixel is in a border region
      if (ele_edge .or. line_edge) then
         uniform = .false.

c ... Check if middle pixel has been flagged as snow or ice.
      else if(snow .or. ice) then
         uniform = .false.

      else

c ...    Check pixels in nlcntx by necntx region for consistency
         do 100 i = 1 , nlcntx
            do 200 j = 1 , necntx
               ide = kele + (j - imv)
c ...          First check consistency of ecosystem
               if ((eco_type - contx_eco(ide,i)) .ne. 0) 
     +              uniform = .false.
c ...          Next, check land/water consistency
               if (contx_topog(ide,i) .eq. 1   .or.  
     +             contx_topog(ide,i) .eq. 4) then
                  nland = nland + 1
               else if(contx_topog(ide,i) .eq. 2 .or. 
     +                 contx_topog(ide,i) .eq. 3) then
                  ncoast = ncoast + 1
               else 
                  nwater = nwater + 1
               end if
c ...          Check snow consistency
               if (snow) then
                 if (contx_snow(ide,i).ne.103
     +             .and.contx_snow(ide,i).ne.104)
     +             uniform = .false.
               elseif (.not. snow) then
                 if (contx_snow(ide,i).eq.103 
     +             .or. contx_snow(ide,i).eq.104)
     +             uniform = .false.
               endif
200         continue
c ...       One more check for a common nadir pixel.  Don't want to
c ...       invoke uniformity test if nadir pixels area not matched up.
            if (n_nadir(i) .ne. n_nadir(nlcntx/2 + 1)) then
               uniform = .false.
            endif
100      continue

c        At the end now, we want to decide if we have
c        a coastline in our region.
         if(nwater + ncoast .eq. itotal ) then
           if (nwater .ne. 9) then
             uniform = .false.
c            Provide "double coastlines".
             coast = .true.
             land = .true.
             water = .false.
           endif
         else if(nland + ncoast .eq. itotal) then
           if (nland .ne. 9) then 
             uniform = .false.
           endif
         else 
           uniform = .false.
           coast = .true.
           land = .true.
           water = .false.
         end if

      endif

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(15x,'' Final Uniformity '',L4)') uniform
        write(h_output,'(2x,''nwater nland ncoast itotal Day Land Water
     +  Coast'',/,4I6,2x,4L5)') nwater,nland,ncoast,itotal,day,land,
     + water,coast
        write(h_output,'(15x,9i5)') ((contx_eco(jj,kk),jj=kele-1,
     *         kele+1),kk=1,3)
      endif
c ................................................................

      return
      end
