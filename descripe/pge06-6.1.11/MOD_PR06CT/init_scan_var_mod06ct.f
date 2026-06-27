      subroutine init_scan_var_mod06ct( line, pixel, qual_flag, conf_flag,
     *  ci_flag, hc_flag, os_top_flag, cldhgt_cat, nearnad_flag,
     *  cldhgtmet_flag )
      
c-----------------------------------------------------------------------
c
c!F77
c
c!Description:
c    Initialize output product arrays for MOD06CT 5x5 sampling as well
c    as QA flags.
c
c!Input Parameters:
c
c    line                        first line of current 5x5
c    pixel                       first pixel of current 5x5
c    QA flags:
c    qual_flag                   quality of product
c    conf_flag                   confidence in product
c    ci_flag                     cirrus cloud found
c    hc_flag                     high cloud
c    os_top_flag                 over-shooting top flag
c    cldhgt_cat                  cloud height category flag
c    nearnad_flag                near-nadir view flag
c    cldhgtmet_flag              index indicating CTP retrieval method
c
c!Output Parameters:
c
c    QA flags:
c    qual_flag                   quality of product
c    conf_flag                   confidence in product
c    ci_flag                     cirrus cloud found
c    hc_flag                     high cloud
c    os_top_flag                 over-shooting top flag
c    cldhgt_cat                  cloud height category flag
c    nearnad_flag                near-nadir view flag
c    cldhgtmet_flag              index indicating CTP retrieval method
c
c!Revision History:
c $Id: init_scan_var_mod06ct.f,v 1.4 1999/04/16 22:55:21 kis Exp $
c
c    Added 'os_top_flag'         R. Frey  (5/10)
c    Added 'cldhgt_cat'          R. Frey  (5/10)
c    Added 'nearnad_flag'        R. Frey  (5/10)
c    Added 'cldhgtmet_flag'      R. Frey  (5/10)
c
c!Team-unique Header:
c
c!End
c
c-----------------------------------------------------------------------

      implicit none
      save

      include 'mod06uw_data.inc'
          
c ... arguments

      integer qual_flag,conf_flag,hc_flag,ci_flag,os_top_flag,
     *        cldhgt_cat,nearnad_flag,cldhgtmet_flag,line,pixel

c ... local

      integer mid_px, mid_ln

c     Initialize flags
      qual_flag = 0
      conf_flag = 0
      hc_flag = 0
      ci_flag = 0
      os_top_flag = 0
      cldhgt_cat = 0
      cldhgtmet_flag = 0

      mid_px = pixel + isamp / 2
      mid_ln = line + isamp / 2
      if(view(mid_px,mid_ln) .le. near_nadir_vza_limit) then
        nearnad_flag = 1
      else
        nearnad_flag = 2
      end if

      return
      end
