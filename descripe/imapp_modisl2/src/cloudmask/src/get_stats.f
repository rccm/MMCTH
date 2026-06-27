      subroutine get_stats(day,land,water,snglnt,snow,coast,desert,
     +                     ice,shadow,smoke,cirrus_ir,cirrus_vis,
     +                     confdnc,nmtests,nbands,bad_value,bad_geo,
     +                     geo_flag,no,ns,nd,nt,ng,ni,
     +                     nl,nw,nr,nv,nu,nn,nz,na,ne,
     +                     npix,n1sm,n2sm,n3sm,n4sm,num_bad)

      implicit none
      save

c----------------------------------------------------------------------
C!F77 
c
c!Description:
c Routine for summing up the various qa granule statistics to be
c filed as PSA's.
c
c!Input Parameters:
c day           Logical variable flagging day scenes
c land          Logical variable flagging land scenes
c water         Logical variable flagging water scenes
c snglnt        Logical variable flagging sunglint contaminated scenes
c snow          Logical variable flagging snow background scenes
c coast         Logical variable flagging coastal scenes
c desert        Logical variable flagging desert background scenes
c ice           Logical variable flagging ice background scenes
c shadow        Logical variable flagging shadow contaminated scenes
c smoke		Logical variable flagging smoke contaminated scenes
c cirrus_ir     Logical variable flagging thin cirrus contaminated
c               scenes in the infrared
c cirrus_vis    Logical variable flagging thin cirrus contaminated
c               scenes in the visible
c confdnc       Unobstructed fov confidence
c nmtests       Number of tests applied to this pixel
c nbands        Number of bands successfully read for this pixel
c bad_value     Logical variable indicating bad pixel radiance or 
c               reflectance bits      
c bad_geo       Logical variable flagging bad lat/long data
c geo_flag      Integer array containing geolocation good/bad flags
c               1-lat,2-lon,3-szen,4-vzen,5-rel_angle
c
c!Output Parameters:
c no		Sum of non-cloud obstruction pixels
c ns		Sum of shadow pixels found
c nd		Sum of day processed pixels
c nt		Sum of night processed pixels
c ng		Sum of sunglint found pixles
c ni		Sum of snow/ice processed pixels 
c nl		Sum of land processed pixels
c nw		Sum of water processed pixels
c nr		Sum of thin cirrus (ir) found pixels
c nv		Sum of thin cirrus (vis) found pixels
c nu            Sum of bad lat pixels
c nn            Sum of bad lon pixels
c nz            Sum of bad szen pixels
c na            Sum of bad vzen pixels
c ne            Sum of bad rel pixels
c npix          Counter for number of pixels included in stats
c n1sm          Confidence category 1 sum
c n2sm          Confidence category 2 sum
c n3sm          Confidence category 3 sum
c n4sm          Confidence category 4 sum
c num_bad       Bad data sum
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
    
c     scalar arguments
      integer npix,n1sm,n2sm,n3sm,n4sm,num_bad,no,ns,nd,nt,
     +        ng,ni,nl,nw,nr,nv,nmtests,nbands,
     +        nu,nn,nz,na,ne
      real confdnc
      logical bad_value,smoke,shadow,day,snglnt,snow,ice,land,
     +        coast,desert,water,cirrus_ir,cirrus_vis,bad_geo
c ... array arguments
      integer geo_flag(5)
      
c      
c     local scalars
      integer nmpix,n1s,n2s,n3s,n4s,nbad,nco,nshadow,nday,nite,
     +        nsnglnt,nice,nland,nwater,ncir,ncvis,
     +        debug,h_output,nlat,nlon,nsza,nvza,nrel

c ... Common statement for debug purposes
      common / bug / debug, h_output

c ... Data statement for holding variable initialization
      data nmpix /0/ , n1s /0/ , n2s /0/ , n3s /0/ , n4s /0/
      data nbad /0/ , nday /0/ , nshadow /0/ , nite /0/ , nice /0/
      data nland /0/ , nwater /0/ , nco /0/ , nsnglnt /0/
      data ncir /0/ , ncvis /0/
      data nlat /0/ , nlon /0/ , nsza /0/ , nvza /0/ , nrel /0/

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(10x/,''Within get_stats  routine '',/)')
      endif
c ...............................................................

      nmpix = nmpix + 1

c     if no lat/lon, or no tests or bands were used than don't count
      if (bad_geo .or. 
     +   (bad_value .and. (nmtests .eq. 0 .or. nbands .eq. 0))) then
         nbad = nbad + 1
c     else add to confidence sum
      else
         if (smoke) then
            nco = nco + 1
         endif
         if (shadow) then
            nshadow = nshadow + 1
         endif
         if (day) then
            nday = nday + 1
         else
            nite = nite + 1
         endif
         if (snglnt) then
            nsnglnt = nsnglnt + 1
         endif
         if (snow .or. ice) then
            nice = nice + 1
         endif

         if (land .or. coast .or. desert) then
            nland = nland + 1
         endif
           
         if (water) then
             nwater = nwater + 1
         endif

         if (cirrus_ir) then
             ncir = ncir + 1
         endif
        
         if (cirrus_vis) then
             ncvis = ncvis + 1
         endif
 
         if (geo_flag(4) .eq. 1) nvza = nvza + 1
         if (geo_flag(5) .eq. 1) nrel = nrel + 1

         if(confdnc .gt. 0.99) then
           n1s = n1s + 1
         else if(confdnc .gt. 0.95) then
           n2s = n2s + 1
         else if(confdnc .gt. 0.66) then
           n3s = n3s + 1
         else if(confdnc .le. 0.66) then
           n4s = n4s + 1
         end if
      endif
      
c ... Sum up bad geo flags what would cause bad_geo to be true
      if (bad_geo) then
         if (geo_flag(1) .eq. 1) nlat = nlat + 1
         if (geo_flag(2) .eq. 1) nlon = nlon + 1
         if (geo_flag(3) .eq. 1) nsza = nsza + 1
         if (geo_flag(4) .eq. 1) nvza = nvza + 1
         if (geo_flag(5) .eq. 1) nrel = nrel + 1
      endif

      npix = nmpix
      num_bad = nbad
      n1sm = n1s
      n2sm = n2s
      n3sm = n3s
      n4sm = n4s

      no = nco
      ns = nshadow 
      nd = nday
      nt = nite
      ng = nsnglnt
      ni = nice
      nl = nland
      nw = nwater
      nr = ncir
      nv = ncvis
      nu = nlat
      nn = nlon
      nz = nsza
      na = nvza
      ne = nrel

c ... debug statement ............................................
      if (debug .gt. 0) then
        write(h_output,'(10x,'' get_stats confidence vars:'',/,6i10,/)')
     *        npix,num_bad,n1sm,n2sm,n3sm,n4sm
        write(h_output,'(10x,'' get_stats psa vars:'',/,6i10,/)')
     *        no,ns,nd,nt,ng,ni
        write(h_output,'(10x,'' get_stats psa vars2:'',/,4i10,/)')
     *        nl,nw,nr,nv
        write(h_output,'(10x,'' get_stats psa vars3:'',/,5i10,/)')
     *        nu,nn,nz,na,ne
      endif
c ................................................................

      return
      end
