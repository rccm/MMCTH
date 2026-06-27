      subroutine mod06ct_set_qa( line, pixel, nmissing, ncloudy,
     &  qual_flag, conf_flag, hc_flag, ci_flag, os_top_flag,
     &  cldhgt_cat, nearnad_flag, cldhgtmet_flag )

c-----------------------------------------------------------------------
c
c!F77
c
c!Description:
c    Set the appropriate MOD06 QA bits for a box of pixels for the
c    cloud top properties.
c
c!Input Parameters:
c    line            Line number at corner of this box
c    pixel           Pixel number at corner of this box
c    nmissing        Number of missing (bad radiance) pixels in this box
c    ncloudy         Number of cloudy pixels in this box
c    qual_flag       quality of product flag
c    conf_flag       confidence in product flag
c    ci_flag         cirrus cloud found flag
c    hc_flag         high cloud found flag
c    os_top_flag     overshooting thunderstorm top flag
c    cldhgt_cat      Indicates lo, mid, hi clouds for QA output
c    nearnad_flag    Indicates near-nadir or oblique view
c    cldhgtmet_flag  Indicates CTP retrieval method
c
c
c!Output Parameters: None
c
c!Revision History:
c $Id: mod06ct_set_qa.f,v 1.5 1999/06/11 22:39:30 kis Exp $
c
c!Team-unique Header:
c
c!End
c
c-----------------------------------------------------------------------
      
      implicit none

      
      include 'mod06uw_data.inc'
      include 'mod06uw_debug.inc'

c ... arguments

      integer line, pixel, nmissing, ncloudy, 
     +        qual_flag, conf_flag, hc_flag, ci_flag, os_top_flag,
     +        cldhgt_cat, nearnad_flag, cldhgtmet_flag

c ... local variables

      integer iout, jout, nclear
      
c whuang 05/07/01: Added save statement for common / mod06_data / and
c                  / mod06_debug /
      save / mod06_data /, / mod06_debug /

c ... Pixel and line coordinates in output 5x5 sampled array
      
      iout = ( pixel / isamp ) + 1
      jout = ( line / isamp ) + 1

c ... If retrieval was good, set product run time QA flags      

      if ( qual_flag .eq. 1 ) then
      
c ...   Set the quality bits first
c ...   Cloud Top Pressure 
        call set_qa_bit( product_qa( 1, iout, jout ), 0 )

c ...   Cloud Top Temperature
        call set_qa_bit( product_qa( 1, iout, jout ), 4 )

c ...   Cloud Fraction 
        call set_qa_bit( product_qa( 1, iout, jout ), 8 )

c ...   Cloud Effective Emissivity
        call set_qa_bit( product_qa( 1, iout, jout ), 12 )

c ...   Cloud Height
        call set_qa_bit( product_qa( 1, iout, jout ), 48 )

c ...   Now the confidence bits.  Note that cloud fraction always
c ...   has a confidence of "3" (C6 change 5-14-10, R. Frey).

c ...   Cloud Fraction 
        call set_qa_bit( product_qa( 1, iout, jout ), 9 )
        call set_qa_bit( product_qa( 1, iout, jout ), 10 )

        if (conf_flag .eq. 1) then

c ...     Cloud Top Pressure 
          call set_qa_bit( product_qa( 1, iout, jout ), 1 )

c ...     Cloud Top Temperature
          call set_qa_bit( product_qa( 1, iout, jout ), 5 )

c ...     Cloud Fraction 
c         call set_qa_bit( product_qa( 1, iout, jout ), 9 )

c ...     Cloud Effective Emissivity
          call set_qa_bit( product_qa( 1, iout, jout ), 13 )
 
c ...     Cloud Height
          call set_qa_bit( product_qa( 1, iout, jout ), 49 )
 
        else if (conf_flag .eq. 2) then

c ...     Cloud Top Pressure 
          call set_qa_bit( product_qa( 1, iout, jout ), 2 )

c ...     Cloud Top Temperature
          call set_qa_bit( product_qa( 1, iout, jout ), 6 )

c ...     Cloud Fraction 
c         call set_qa_bit( product_qa( 1, iout, jout ), 10 )

c ...     Cloud Effective Emissivity
          call set_qa_bit( product_qa( 1, iout, jout ), 14 )

c ...     Cloud Height
          call set_qa_bit( product_qa( 1, iout, jout ), 50 )
 
        else if (conf_flag .eq. 3) then

c ...     Cloud Top Pressure 
          call set_qa_bit( product_qa( 1, iout, jout ), 1 )
          call set_qa_bit( product_qa( 1, iout, jout ), 2 )

c ...     Cloud Top Temperature
          call set_qa_bit( product_qa( 1, iout, jout ), 5 )
          call set_qa_bit( product_qa( 1, iout, jout ), 6 )

c ...     Cloud Fraction 
c         call set_qa_bit( product_qa( 1, iout, jout ), 9 )
c         call set_qa_bit( product_qa( 1, iout, jout ), 10 )

c ...     Cloud Effective Emissivity
          call set_qa_bit( product_qa( 1, iout, jout ), 13 )
          call set_qa_bit( product_qa( 1, iout, jout ), 14 )

c ...     Cloud Height
          call set_qa_bit( product_qa( 1, iout, jout ), 49 )
          call set_qa_bit( product_qa( 1, iout, jout ), 50 )
 
        endif
      endif

c ... Don't forget to set the cirrus and high cloud flags
      if (ci_flag .eq. 2 ) then
        call set_qa_bit( product_qa( 1, iout, jout ), 21 )
      else if  (ci_flag .eq. 1 ) then
        call set_qa_bit( product_qa( 1, iout, jout ), 20 )
      else if  (ci_flag .eq. 3 ) then
        call set_qa_bit( product_qa( 1, iout, jout ), 20 )
        call set_qa_bit( product_qa( 1, iout, jout ), 21 )
      endif
      if (hc_flag .eq. 2 ) then
        call set_qa_bit( product_qa( 1, iout, jout ), 23 )
      else if (hc_flag .eq. 1 ) then
        call set_qa_bit( product_qa( 1, iout, jout ), 22 )
      else if (hc_flag .eq. 3 ) then
        call set_qa_bit( product_qa( 1, iout, jout ), 22 )
        call set_qa_bit( product_qa( 1, iout, jout ), 23 )
      endif
c ... And the over-shooting thunderstorm top flag
      if (os_top_flag .eq. 2) then
        call set_qa_bit( product_qa( 1, iout, jout ), 53 )
      else if( os_top_flag .eq. 1) then
        call set_qa_bit( product_qa( 1, iout, jout ), 52 )
      end if


c ... Set optional run time QA flags

      product_qa( 4, iout, jout ) = ncloudy
      nclear = isamp*isamp - nmissing - ncloudy
      product_qa( 5, iout, jout ) = nclear
      product_qa( 6, iout, jout ) = nmissing

c ... Finally, set the input data resource flags - Setting
c ...  here because I really don't know where else to
c ...  set them.
c ...  Clear radiance origin
       if (mod06ct_proc_opt .eq. 1 ) then
         call set_qa_bit( product_qa( 1, iout, jout ), 54 )
       else if (mod06ct_proc_opt .eq. 3 ) then
         call set_qa_bit( product_qa( 1, iout, jout ), 55 )
       endif

c ...  Moisture profile, Temp profile, surface temp over land
       if (tem1(pixel,line) .lt. 1.0) then
          call set_qa_bit( product_qa( 1, iout, jout ), 56 )
          call set_qa_bit( product_qa( 1, iout, jout ), 57 )
          call set_qa_bit( product_qa( 1, iout, jout ), 58 )
          call set_qa_bit( product_qa( 1, iout, jout ), 59 )
          call set_qa_bit( product_qa( 1, iout, jout ), 60 )
          call set_qa_bit( product_qa( 1, iout, jout ), 61 )
          call set_qa_bit( product_qa( 1, iout, jout ), 62 )
          call set_qa_bit( product_qa( 1, iout, jout ), 63 )
       endif

c ...  surface pressure
       if (pre1(pixel,line) .lt. 1.0) then
         call set_qa_bit( product_qa( 1, iout, jout ), 64 )
         call set_qa_bit( product_qa( 1, iout, jout ), 65 )
       endif
        
c ...  Surface emissivity currently not used
       call set_qa_bit( product_qa( 1, iout, jout ), 68 )
       call set_qa_bit( product_qa( 1, iout, jout ), 69 )

c ...  Ecosystem type
       if (eco(pixel,line) .lt. 0) then
         call set_qa_bit( product_qa( 1, iout, jout ), 70 )
         call set_qa_bit( product_qa( 1, iout, jout ), 71 )
       endif

c ...  Cloud Height Category
       if(cldhgt_cat .eq. 1) then
         call set_qa_bit( product_qa( 1, iout, jout ), 72 )
       else if(cldhgt_cat .eq. 2) then
         call set_qa_bit( product_qa( 1, iout, jout ), 73 )
       else if(cldhgt_cat .eq. 3) then
         call set_qa_bit( product_qa( 1, iout, jout ), 72 )
         call set_qa_bit( product_qa( 1, iout, jout ), 73 )
       else if(cldhgt_cat .eq. 4) then
         call set_qa_bit( product_qa( 1, iout, jout ), 74 )
       else if(cldhgt_cat .eq. 5) then
         call set_qa_bit( product_qa( 1, iout, jout ), 72 )
         call set_qa_bit( product_qa( 1, iout, jout ), 74 )
       end if

c ...  Nadir View Angle Flag
       if(nearnad_flag .eq. 1) then
         call set_qa_bit( product_qa( 1, iout, jout ), 75 )
       else if(nearnad_flag .eq. 2) then
         call set_qa_bit( product_qa( 1, iout, jout ), 76 )
       end if
       
c ...  Cloud Height Method 
       if(cldhgtmet_flag .eq. 1) then
         call set_qa_bit( product_qa( 1, iout, jout ), 77 )
       else if (cldhgtmet_flag .eq. 2) then
         call set_qa_bit( product_qa( 1, iout, jout ), 78 )
       else if (cldhgtmet_flag .eq. 3) then
         call set_qa_bit( product_qa( 1, iout, jout ), 77 )
         call set_qa_bit( product_qa( 1, iout, jout ), 78 )
       else if (cldhgtmet_flag .eq. 4) then
         call set_qa_bit( product_qa( 1, iout, jout ), 79 )
       else if (cldhgtmet_flag .eq. 5) then
         call set_qa_bit( product_qa( 1, iout, jout ), 77 )
         call set_qa_bit( product_qa( 1, iout, jout ), 79 )
       else if (cldhgtmet_flag .eq. 6) then
         call set_qa_bit( product_qa( 1, iout, jout ), 78 )
         call set_qa_bit( product_qa( 1, iout, jout ), 79 )
       else if (cldhgtmet_flag .eq. 7) then
         call set_qa_bit( product_qa( 1, iout, jout ), 77 )
         call set_qa_bit( product_qa( 1, iout, jout ), 78 )
         call set_qa_bit( product_qa( 1, iout, jout ), 79 )
       end if
       
c ... debug statement ............................................
      if (debug .gt. 2) then
         write(h_output,'(1x,3x,'' CT QA Bits set '')')
         write(h_output,'(''qual_flag,conf_flag,hc_flag,ci_flag,os_top_flag,'',
     &    ''cldhgt_cat,nearnad_flag,cldhgtmet_flag '',8i3)') qual_flag,conf_flag,
     &    hc_flag,ci_flag,os_top_flag,cldhgt_cat,nearnad_flag,cldhgtmet_flag
      endif
c ................................................................



      
      return
      end
