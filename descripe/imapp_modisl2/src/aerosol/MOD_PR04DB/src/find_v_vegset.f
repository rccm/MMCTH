      subroutine find_v_veg_test(itmp,jtmp,realbuf,tmpvg,
     1           outbufvg)
!c.  1           outbufvg,qa_flag)

c      use lut_arrays   

      include 'aottbl.inc'
      include 'newaottbl.inc'

      real realbuf(26), outbufvg(20)
c      real realbuf(13), outbuf(20)
      integer qa_flag(4)
      integer month
c...
c... Driver for MODIS aerosol retrieval over vegetated surfaces
c... This code takes input of MODIS reflectance (binary form -
c... cloud-screened data), and calls main-subroutine for aerosol
c... retrieval over vegetated surfaces (subroutine myd_aero_veg).
c... It utilizes inter-bands sfc refl relationships between 2.1um and 0.47 & 0.65um.
c...
c... A version to process MODIS data on granule-basis
c... Inherited from "find_aero_470_v4_appl_all_10x10avgrelatn_whlgrnl.f"
c... --> copied from "find_aero_modis_veg_v042_stdalone_AM01-04.f"
c...
c...    Written by Myeong-Jae Jeong (MJ)
c...    Last modified Aug 08, 2011
c...
c...

      parameter(nxgrd=1354,nygrd=2030)
c.    parameter(nx2=900, ny2=500)
      parameter(nx2=3600, ny2=1800) ! landcover data

c      real   realbuf(13)
      real*4 xlat, xlon, tmpvg(6)
      integer*4 begday,endday
      integer*4 ihh, imm, itmp, jtmp !, lsf
      character*48 fninp
      character*55 fninp2
      character*29 fnbinin  ! MODIS binary refl file
c      character*29 pathbin  ! binary refl file location
c      character*56 pnbinin  ! full path of binary refl file
      character*31 fnout,fnout2  ! daily bin/asc refl files
c     character*33 pnout,pnout2  ! full path bin/asc files
c      character*42 pnout,pnout2  ! full path bin/asc files
c     character*31 sitename(nsites)
      character*4  swstr
      character*27 fnpwbin,fnlatbin,fnlonbin
      character*51 pathpw
      character*78 pnpwbin,pnlatbin,pnlonbin
c----------------<-- for variables associated with binary files
c
c     ---input tables
c.    real nvalx412(10,46,30,10,8,20)
c      real nvalx470(10,46,30,10,4,24),nvalx650(10,46,30,10,24)
c.    real sfc_ref412(20),sfc_ref470(24), sfc_ref650(24)
c      real*4 xlcvr_2(nx2,ny2)

c     ---input parameters
      real tmp(14), xnvalm6(6)

c     ---intermediate parameters
c.    character*4 amidstr
c      integer*8 naccrow
c      character*32  pathout

c    added by ces beginning 29 jun 2010
      integer*4 doy
      common    /xday/ doy

c.    common /angle_node/ theta0, theta, phi
c.    common /sfcref_node/ sfc_ref412, sfc_ref470, sfc_ref650
c.    common /fname_node/ aer_tab, w0_name470
c.    data pi      /3.14159/

c... ****************************************************************
c... ****************************************************************
c... ****************************************************************
      
      xday = real(doy)

c      print *, 'HERE', xday
c      print *, itmp, jtmp
c... Input/output data paths
c... -------------------------
c/dust/csalustro/vegDB/resulting/
c      pathout='/dust/csalustro/vegDB/resulting/'     ! output path, including temporary output
c      pathbin='/dust/csalustro/vegDB/binary/'
      
c... ****************************************************************
c... ****************************************************************
c... ****************************************************************


c-----------------------------------------------------------------------
c      Load Tables
c-----------------------------------------------------------------------
c... read LUTs of TOA reflectance and Global Land Cover
c      call get_lut4vegaer(nvalx470,nvalx650,xlcvr_2)
c      call get_lut4vegaer(xlcvr_2)
c-----------------------------------------------------------------------


c... ==================================================
c...  output binary file info
c      open(2,file=pathout//'numrow_info_s1.txt',status='unknown')
c... ==================================================

c... different aerosol models (SSA) for 470nm band
c... ##########
c... ##########
c... ##########
c.    imod=3  ! test
c... ##########
c... ##########
c... ##########
c.    amidstr='AM01'
c.    do 2500 imod=1,4
c.    write(amidstr(3:4),'(i2.2)') imod

c-----------------------------------------------------------------------
c      Start processing the data
c-----------------------------------------------------------------------
c...  get the list of input MODIS binary file(s)
C         open(11,file='/dust/csalustro/vegDB/Clare.dir/grnl_listtmp.dat',status='old')
C         read(11,*) nf_bin   ! number of granules per day
C
Cc        noob=0
C         do 1010 ij=1,nf_bin ! different granules in a day
C            read(11,211) fnbinin
Cc!... need to get obs time info from "fnbinin"
C            read(fnbinin,214) iyr_in, doy
C            read(fnbinin,212) ihh, imm
C            print *, 'processing ... ', iyr_in, doy, ihh, imm
C212   format(20x,2i2,4x)     ! @@@ tmp --> should match with "fnbinin"
C214   format(12x,i4,i3,10x)  ! @@@ tmp --> should match with "fnbinin"
C
C!/dust/csalustro/vegDB/binary/MYDL1B_TmpA2004153.2130.bin
Cc            print *, 'pathbin: ', pathbin
C            print *, 'fnbinin: ',fnbinin
Cc            pnbinin=pathbin//fnbinin(2:28)
Cc            call system('ls -l '//pnbinin)
C
C            fnout= 'MYD_outx.A2004001.0005.'//amidstr//'.bin'  ! @@@
Cc            fnout2='MYD_outx.A2004001.0005.'//amidstr//'.asc'
C            write(fnout(11:14),'(i4.4)') iyr_in
C            write(fnout(15:17),'(i3.3)') doy
C            write(fnout(19:22),'(2i2.2)') ihh, imm
CC            write(fnout2(11:14),'(i4.4)') iyr_in
CC            write(fnout2(15:17),'(i3.3)') xday
CC            write(fnout2(19:22),'(2i2.2)') ihh, imm
Cc            pnout=pathout//fnout
Cc            pnout2=pathout//fnout2
C
Cc... open and read data for a granule
C            
CC            open(3,err=51,file=pathbin//fnbinin(2:28),access='sequential',
CC     &           form='unformatted',status='old')
C            open(15,file=pathout//fnout,access='sequential',
C     &           form='unformatted',status='unknown')
Cc...        open(16,file=pnout2,status='unknown')

c-------------------------------------------------------------
c            naccrow=0
c   10       read(3,end=200) itmp, jtmp, realbuf, lsf
c            if(itmp.le.525.and.jtmp.le.525.and.imod.eq.1) print *, itmp, jtmp, realbuf
c   10       continue

c-------------------------------------------------------------

c-----------------------------------------------------------------------
c...  assign season ID ...
c-----------------------------------------------------------------------
      if(xday.ge.60.and.xday.lt.152) then      ! MAM
         iopss = 1
      elseif(xday.ge.152.and.xday.lt.244) then
         iopss = 2                                 ! JJA
      elseif(xday.ge.244.and.xday.lt.335) then
         iopss = 3   ! SON, iopss def. change (4-->3; as of Jun 15, 2010)
      elseif(xday.ge.335.and.xday.lt.367.or.xday.ge.0) then
         iopss = 1   ! --> slot for 4, DJF @@@@@ temporary set-up @@@@@
c        iopss = 4
      else
         print *, 'Invalid iopss value ...'
         go to 10
         stop
      endif
      
      month = 0
      if (xday .ge. 1 .AND. xday .le. 31)      month = 1
      if (xday .ge. 32 .AND. xday .le. 59)     month = 2
      if (xday .ge. 60 .AND. xday .le. 90)     month = 3
      if (xday .ge. 91 .AND. xday .le. 120)    month = 4
      if (xday .ge. 121 .AND. xday .le. 151)   month = 5
      if (xday .ge. 152 .AND. xday .le. 181)   month = 6
      if (xday .ge. 182 .AND. xday .le. 212)   month = 7
      if (xday .ge. 213 .AND. xday .le. 243)   month = 8
      if (xday .ge. 244 .AND. xday .le. 273)   month = 9
      if (xday .ge. 274 .AND. xday .le. 304)   month = 10
      if (xday .ge. 305 .AND. xday .le. 334)   month = 11
      if (xday .ge. 335 .AND. xday .le. 366)   month = 12
      
c-----------------------------------------------------------------------
c           if(itmp.ge.100.and.itmp.lt.110.and.
c    1         jtmp.ge.200.and.jtmp.lt.210)
c    2      print *, itmp, jtmp, realbuf, lsf, iopss, xday


c-----------------------------------------------------------------------
c... aerosol retrieval over vegetated surface
c-----------------------------------------------------------------------
      call myd_aero_veg(month,iopss,itmp,jtmp,realbuf,tmpvg,outbufvg)
c.    call myd_aero_veg(imod,iopss,itmp,jtmp,realbuf,tmpvg,outbufvg)
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c... write retrieval results on a binary file
c.     write(15) itmp,jtmp,outbufvg
c-----------------------------------------------------------------------
c.           if(itmp.ge.100.and.itmp.lt.110.and.
c.   1         jtmp.ge.200.and.jtmp.lt.210)
c.   2 write(*,811) itmp,jtmp,(outbufvg(kk),kk=1,20)
c.    if(itmp.ge.100.and.itmp.lt.103.and.
c.   1 jtmp.ge.200.and.jtmp.lt.202)
c.   2 write(*,*) 'writing fort.16 and fort.23 for test...'
c.    write(16,811) itmp,jtmp,(outbufvg(kk),kk=1,20)
c      naccrow=naccrow+1


c  51  print *, 'Error in reading binary files' 
c     stop                                     

c200   close (3)  ! close input binary files (indiv. granule)
c      close(15)  ! close output binary file, "pnout" (indiv. granule)

c      write(2,322) naccrow, fnout, amidstr
c1010  continue   ! different granules
c      close(11)  ! file with a list of granules to process
c     print *, noob, ' lcvr data were out of bound'

c.2500  continue ! aerosol model
c      close(2)

      go to 10

c-------------------------------------------------------------------
c... format statements
c-------------------------------------------------------------------
  322 format(i10,2x,a31,2x,a4)
  211 format(a29)
  811 format(2(i4,1x),8(1x,f9.4),6(2x,f8.3),3x,f7.1,2(2x,f10.5),
     1       2(2x,f9.4),2x,f7.1)
c-------------------------------------------------------------------

c... ****************************************
c... The end of main
c... ****************************************

10    continue
      return
      end
c... ****************************************
c
c-------------------------------------------------------------------
c


c... ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine myd_aero_veg(month, iopss,itmp,jtmp,realbuf,tmpvg,
     1           outbufvg)
c.    subroutine myd_aero_veg(imod,iopss,itmp,jtmp,realbuf,tmpvg,
c.   1           outbufvg)

c...  Main subroutines for aerosol retrieval over vegetated surfaces.
c...     Written by Myeong-Jae Jeong (MJ)
c...     Last modified Aug 08, 2011
c...
c...  Inputs:
c...     iopss:    season ID (1:MAM; 2:JJA; 3:SON; 4:DJF)
c...     imod:     aerosol model index for 470nm --> no longer an input
c...     nvalx470: TOA refl LUT for 470nm
c...     nvalx650: TOA refl LUT for 650nm
c...     xlcvr_2:  MODIS land cover (IGBP)
c...     regid_2:  Region Index (to choose a set of aerosol models
c...               depending on regions and land cover)
c...     realbuf:  see below for definitions
c...     tmpvg:
c...
c...  Output:
c...     outbufvg: see the end of this subroutine for definitions
c...
      include 'aottbl.inc'
      include 'newaottbl.inc'

      parameter(nx2=3600, ny2=1800)  ! landcover data dimension
      
      integer month
c      real nvalx470(10,46,30,10,4,24),nvalx650(10,46,30,10,24)
      real nval(10,46,30), yy(10), yyw(8) !tau(10), 
c.    real xnvalm6(6), realbuf(13), outbuf(20)
      real xnvalm6(6), realbuf(26), outbufvg(20), tmpvg(6)
c      real*4 xlcvr_2(nx2,ny2)
      integer tau_x470_flag, tau_x650_flag

c      real theta0(10), theta(46), phi(30)
c      real sfc_ref412(20), sfc_ref470(24), sfc_ref650(24)
      real*4 r412db,r470db,r650db,cl_flag
      real xtau(3),ssa(3),qa_flag(4),aot_mod(6) !,w0_470(4)
      character*4 w0_name470(4)
      character*12 aer_tab(10)
      real*4 ctharr(3)     ! 08/05/2011
      integer imodarr(3)   ! 08/05/2011

c      common /angle_node/ theta0, theta, phi
c      common /sfcref_node/ sfc_ref412, sfc_ref470, sfc_ref650
      common /fname_node/ aer_tab, w0_name470
      data pi      /3.14159/
C      data theta0  /0.0,8.0,16.0,24.0,32.0,40.0,48.0,56.0,64.0,72.0/
C      data tau     /0.0, 0.1, 0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5/
C      data w0_470  /0.91, 0.94, 0.96, 0.99/
C      data sfc_ref412 /1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,
C     1                 11.,12.,13.,14.,15.,16.,17.,18.,19.,20./
C      data sfc_ref470 /1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,
C     1                 11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,
C     2                 21.,22.,23.,24./
C      data sfc_ref650 /1.,3.,5.,7.,9.,11.,13.,15.,17.,19.,21.,
C     1                 23.,25.,27.,29.,31.,33.,35.,37.,39.,41.,
C     2                 43.,45.,47./

c.    data aer_tab /'table_0.0aot', 'table_0.1aot', 'table_0.3aot',
c.   1              'table_0.5aot', 'table_1.0aot', 'table_1.5aot',
c.   2              'table_2.0aot', 'table_2.5aot', 'table_3.0aot',
c.   3              'table_3.5aot'/
      data w0_name470 /'0.91','0.94','0.96','0.99'/
  
c-------------------------------------------
c... define nodes for angles
      do i = 1, 46
         theta(i) = 2.*float(i-1)
      enddo

      do i = 1, 30
         phi(i) = 6. + 6.*float(i-1)
      enddo

      mm = 10     ! solar zenith
      nn = 46     ! satellite zenith
      ll = 30     ! rel. azimuth
      ma = 10     ! tau
      mw =  8     ! ssa
c-------------------------------------------

c-------------------------------------------
c... IGBP Landcover data
      xlatbeg=-90.0
      xlonbeg=-180.0
      xyintv2=0.10
c-------------------------------------------

c... ************************************************
c... ************************************************
c... come back later to set these threshold values, if necessary...(MJ) Jun 15, 2010
      xmnsfc1=1.0  ! allowed min. sfc refl. for B1
      xmxsfc1=47.0 ! allowed max. sfc refl. for B1
      xmnsfc3=1.0  ! allowed min. sfc refl. for B3
      xmxsfc3=24.0 ! allowed max. sfc refl. for B3
c... ************************************************
c... ************************************************

c.?   xlat5  = realbuf(1)
c.?   xlong5 = realbuf(2)
      xlat   = realbuf(1)
      xlong  = realbuf(2)
      sza    = realbuf(3)
      xthet  = realbuf(4)
      xphi   = realbuf(5)
      xnvalm6(1)=tmpvg(1)     ! (B1)
      xnvalm6(2)=tmpvg(2)     ! (B3)
      xnvalm6(3)=tmpvg(3)     ! (B8)
      xnvalm6(4)=tmpvg(4)     ! 2.1um (B7)
      xnvalm6(5)=tmpvg(5)     ! 850nm (B2)
      xnvalm6(6)=tmpvg(6)     ! 1.2um (B5)
c.    band26    =realbuf(10)  ! 1.38um(B26)
      band02    =xnvalm6(5)
      band05    =xnvalm6(6)
      band07    =xnvalm6(4)
      band08    =xnvalm6(3)   ! (B8)
      band03    =xnvalm6(2)   ! (B3)
      band01    =xnvalm6(1)   ! (B1)
c.    cl_flag   =realbuf(13)
      r412db    =realbuf(24)*100.
      r470db    =realbuf(25)*100.
      r650db    =realbuf(26)*100.

c     print *, itmp, jtmp, realbuf, lsf   ! test

c... ######################################
      lcvr=1  ! initialization
      ioprg=5 ! ititialization
c...  get lcvr (IGBP landcover type; integer) data here
      idx=int((xlong-xlonbeg)/xyintv2) + 1
      idy=int((xlat-xlatbeg)/xyintv2) + 1
            
      if(idx.ge.1.and.idx.le.nx2.and.idy.ge.1.and.idy.le.ny2) then
         sfc_typ = xlcvr_2(idx,idy)
         xreg_id = regid_2(idx,idy)  ! 08/05/2011
      else
c        sfc_typ = 1.0*lcvr  ! @@@@@ just for test
         print *, 'lcvr data out of bound: idx,nx2,idy,ny2: ',idx,nx2,idy,ny2,xlat,xlong
c        noob=noob+1
         stop ! @@@@@
      endif
      lcvr=int(sfc_typ)
      ioprg=int(xreg_id)  ! 08/05/2011
c.    if(itmp.ge.100.and.itmp.lt.110.and.jtmp.ge.200.and.
c.   1   jtmp.lt.210) print *, sfc_typ, lcvr, xreg_id, ioprg
c... ######################################


c-----------------------------------------------------------------------

      xmu = cos(sza * 3.14159/180.)

c... refl unit conversion (pi*L/F/xmu --> L/F)
      do i = 1, 6
         xnvalm6(i) = xnvalm6(i) * xmu/3.14159
      enddo

      x1 = sza
      x2 = xthet
      x3 = xphi

c... ============================================================
c... ============================================================
      refl1= xnvalm6(3)    ! 412 nm
      refl6= xnvalm6(1)    ! 650 nm
      refl3= xnvalm6(2)    ! 470 nm
      refl21=band07        ! 2.1um 
c... ============================================================
c... ============================================================


c-----------------------------------------------------------------------
c...  estimate surface reflectance at 470nm & 650nm
c-----------------------------------------------------------------------
      xeAs_B1 = -999.0  ! initialization
      xeAs_B3 = -999.0
      call get_sfcrfl_veg(iopss,refl21,xnvalm6,lcvr,xeAs_B1,xeAs_B3
     1                   ,sirndvi,s_ndvi)

c-----------------------------------------------------------------------

c           if(itmp.ge.100.and.itmp.lt.200.and.
c    1         jtmp.ge.200.and.jtmp.lt.300)
c      if(xeAs_B1.gt.0.and.xeAs_B3.gt.0)
c    2 print *, 'sfc refl= ',xeAs_B1,xeAs_B3,r650db,r470db
c... write outputs on fort.23 for tests !MJ@@@@@
c.    if(xeAs_B1.gt.0.and.xeAs_B1.lt.47.and.xeAs_B3.gt.0.and.
c.   1   xeAs_B3.lt.24)
c.   2write(23,2323) itmp,jtmp,xlat,xlong,sza,xthet,xphi,lcvr,
c.   3  (tmpvg(kk),kk=1,6),xeAs_B1,xeAs_B3,r650db,r470db,r412db,
c.   4  sirndvi,s_ndvi,ioprg,lcvr
2323  format(2(i4,2x),2(f9.4,2x),3(f8.3,2x),i3,6(2x,e14.6),
     1       5(2x,f8.3),2(2x,f9.4),2(2x,i3))

c... Begin aerosol retrieval .....
c... getting thresholds for selecting a set of aerosol models
c... for retrieval depending on regions (08/05/2011)
!... below 37 lines --> @@new@@
      ctharr(:) = -999.0
      imodarr(:) = -999
      if(ioprg.eq.1) then      ! N. America
         ctharr(1)=0.2  ! 
         ctharr(2)=1.0
         ctharr(3)=2.0
         imodarr(1)=3   ! w0=0.96
         imodarr(2)=2   ! w0=0.94
         imodarr(3)=1   ! w0=0.91
      elseif(ioprg.eq.2) then  ! China
         ctharr(1)=0.5  ! cth0
         ctharr(2)=1.5  ! cth1
         ctharr(3)=3.5  ! cth2
         imodarr(1)=2   ! w0=0.94
         imodarr(2)=1   ! w0=0.91
         imodarr(3)=1   ! w0=0.89*-->0.91
      elseif(ioprg.eq.3) then  ! S. America
         ctharr(1)=0.1
         ctharr(2)=0.2
         ctharr(3)=0.5
         imodarr(1)=2   ! w0=0.94
         imodarr(2)=2   ! w0=0.94
         imodarr(3)=2   ! w0=0.94
      elseif(ioprg.eq.6) then  ! S. Africa
         if (month .ge. 6 .AND. month .le. 11) then     ! June through November, use more absorbing model.
           ctharr(1)=0.2
           ctharr(2)=0.5
           ctharr(3)=1.0
           imodarr(1)=1   ! w0=0.91
           imodarr(2)=1   ! w0=0.91
           imodarr(3)=1   ! w0=0.91
         else
           ctharr(1)=0.2
           ctharr(2)=0.5
           ctharr(3)=1.0
           imodarr(1)=2   ! w0=0.94
           imodarr(2)=2   ! w0=0.94
           imodarr(3)=2   ! w0=0.94
         end if
      elseif(ioprg.eq.7) then  ! SE Asia
         if (month .eq. 12 .OR. month .le. 2) then   ! Winter
           ctharr(1)=0.2
           ctharr(2)=0.5
           ctharr(3)=1.0
           imodarr(1)=1   ! w0=0.91
           imodarr(2)=1   ! w0=0.91
           imodarr(3)=1   ! w0=0.91
         else
           ctharr(1)=0.2
           ctharr(2)=0.5
           ctharr(3)=1.0
           imodarr(1)=2   ! w0=0.94
           imodarr(2)=1   ! w0=0.91
           imodarr(3)=1   ! w0=0.91
         end if
      elseif(ioprg.gt.3) then ! default
         ctharr(1)=0.2
         ctharr(2)=0.5
         ctharr(3)=1.0
         imodarr(1)=2   ! w0=0.94
         imodarr(2)=1   ! w0=0.91
         imodarr(3)=1   ! w0=0.91
      elseif(ioprg.le.0) then ! no veg-retrieval 
         ctharr(1)=0.2
         ctharr(2)=0.5
         ctharr(3)=1.0
         imodarr(1)=3   ! w0=0.96
         imodarr(2)=3   ! w0=0.96
         imodarr(3)=3   ! w0=0.96
      endif
c---------------------------------------------------------------
c     Initialization
c---------------------------------------------------------------
c     --intermediate parameters
      w0_x       = -999.
      w0_int     = -999.
      tau_x412   = -999.
      tau_x412_91 = -999.
      tau_x470   = -999.
      tau_x650   = -999.
      xxrat      = -999.
      xxrat2     = -999.
      aot        = -999.
c     -- output parameters
      tau550     = -999.
      alpha      = -999.

      do i = 1, 3
      xtau(i) = -999.
      ssa(i) = -999.
      enddo

      do i = 1, 4
      qa_flag(i) = 0
      enddo

      do i = 1, 6
      aot_mod(i) = -999.
      enddo

      do i = 1,20
        outbufvg(i) = -999.
      enddo

c     sfc_typ = -999.

c      if(xeAs_B1.lt.xmnsfc1.or.xeAs_B1.gt.xmxsfc1.or.
c     1   xeAs_B3.lt.xmnsfc3.or.xeAs_B3.gt.xmxsfc3) return  ! return to main w/o retrieval
c     print *, 'xeAs_B1,xeAs_B3',xeAs_B1,xeAs_B3,refl21,lcvr
      if (xeAs_B1.lt.xmnsfc1.and.xeAs_B1.gt.-900.0) xeAs_B1=xmnsfc1
      if (xeAs_B3.lt.xmnsfc3.and.xeAs_B3.gt.-900.0) xeAs_B3=xmnsfc3
      
c---------------------------------------------------------------
c     Screen for pixels outside reasonable ranges of reflectance
c---------------------------------------------------------------
c... will need to re-open the following lines .... (MJ)
      if (refl1.gt.0.0.and.refl1.lt.0.09.and.
     1    refl6.gt.0.0.and.refl6.lt.0.14) go to 11
c      if (refl1.gt.0.09.and.refl1.lt.0.50.and.
c     1    res.gt.6.0) go to 11
c...   go to 10 ! replaced by "return" below
      if(s_ndvi.gt.0.1.and.refl21.gt.0.01.and.refl21.le.0.25) go to 11 ! remove water/bright sfc contamination
      return   ! return to main w/o retrieval
            
11    continue

c     if (xphi.gt.179.99) go to 10
c     if (xphi.gt.179.99) return   ! choose this or below (MJ)
      if (xphi.gt.179.99) xphi=179.99  ! choose this or above (MJ)
      if (xphi.lt.6.0) xphi = 6.

c     -- sun glint mask
c
      cc     = 3.14159/180.
      psi    = acos(cos(sza*cc)*cos(xthet*cc) +
     1         sin(sza*cc)*sin(xthet*cc)*cos(xphi*cc))
      glint_ang = psi/cc

c      if (abs(psi/cc).lt.35.0) go to 10
c
c     -- scattering angle (scat_ang)
c
      cc     = 3.14159/180.
      psi    = acos(cos(sza*cc)*cos(xthet*cc) -
     1         sin(sza*cc)*sin(xthet*cc)*cos(xphi*cc))
      scat_ang = 180. - psi/cc

c      if (scat_ang .gt. 158.) go to 10
c      if (scat_ang .gt. 175.) go to 10

c... added here on Feb 26, 2010
      x1 = sza
      x2 = xthet
      x3 = xphi

c--------------------------------------------------------
c   Input surface reflectance at 470 nm
c--------------------------------------------------------
c   Input surface reflectance at 470 nm
      r470 = xeAs_B3
c
c   Input toa reflectance (L/F) at 470 nm
      refl = refl3
c     x3 = xphi                 ! ori. position (moved up)

c     Retrieving 470 nm AOT
c     imod = 1                  ! w0 = 0.91
c     imod = 2                  ! w0 = 0.94
c     imod = 3                  ! w0 = 0.96
c     imod = 4                  ! w0 = 0.995
!c... below 40 lines --> AOT-dependent aerosol model selection (new@@@)
      As_21=-999.0
      tau_x470    = -999.
      tau_x4701   = -999.
      tau_x4702   = -999.
      tau_x4703   = -999.
      tau_x470_flag=0
      tau_x650_flag=0
      tau_x470_flag1=0
      tau_x470_flag2=0
      tau_x470_flag3=0
      cth0=ctharr(1)
      cth1=ctharr(2)
      cth2=ctharr(3)
c     print *, ctharr
c     print *, cth0, cth1, cth2
      imod=imodarr(1)
      call aero_470veg(refl,x1,x2,x3,mm,nn,ll,ma,
     1    imod,r470,tau_x4701,tau_x470_flag1)

c... **********************************************
c... **********************************************
c.... new addition: calculate 2.1um sfc refl. using
c...  first guess of 470nm AOT
c...  Aug. 12, 2011
      tau_x4701_00=tau_x4701
      xeAs_B3_00=xeAs_B3
      tref21in=xnvalm6(4)
      alpest=1.5
      aot21est=tau_x4701*(466.0/2110.)**alpest
      if(aot21est.lt.0.0.or.aot21est.gt.0.5) return  ! no retr. for too-heavy aerosol
      if(aot21est.ge.0.0.and.aot21est.lt.0.06) go to 333
c.
      call calc_sfc21(x1,x2,x3,tref21in,aot21est,As_21)
c.    
c...  refl21r=3.14159*As_21/100.0/xmu
      refl21r=As_21/100.0
      if(refl21r.lt.0) return
      call get_sfcrfl_veg(iopss,refl21r,xnvalm6,lcvr,xeAs_B1,xeAs_B3
     1                   ,sirndvi,s_ndvi)
      if(xeAs_B1.lt.xmnsfc1.or.xeAs_B1.gt.xmxsfc1.or.
     1   xeAs_B3.lt.xmnsfc3.or.xeAs_B3.gt.xmxsfc3) return  ! return to main w/o retrieval
      r470=xeAs_B3
c... **********************************************
c... **********************************************
      call aero_470veg(refl,x1,x2,x3,mm,nn,ll,ma,
     1    imod,r470,tau_x4701,tau_x470_flag1)

333   continue
c... ********
c      if(itmp.ge.800.and.itmp.lt.820.and.
c     1         jtmp.ge.900.and.jtmp.lt.920)
c     2  print *,100.*refl21,As_21,aot21est,tau_x4701_00,tau_x4701,
c     3          xeAs_B3_00, xeAs_B3 
c... ********

      if(tau_x4701.ge.0.0.and.tau_x4701.lt.cth0) then
         tau_x470=tau_x4701
         tau_x470_flag=tau_x470_flag1
      elseif(tau_x4701.ge.cth0.and.tau_x4701.lt.cth1) then
         imod=imodarr(2)
         call aero_470veg(refl,x1,x2,x3,mm,nn,ll,ma,
     1      imod,r470,tau_x4702,tau_x470_flag2)
         f_1=(tau_x4701-cth0)/(cth1-cth0)
         tau_x470=f_1*tau_x4702+(1.0-f_1)*tau_x4701
         tau_x470_flag=tau_x470_flag2
      elseif(tau_x4701.ge.cth1.and.tau_x4701.lt.cth2) then
         imod=imodarr(2)
         call aero_470veg(refl,x1,x2,x3,mm,nn,ll,ma,
     1      imod,r470,tau_x4702,tau_x470_flag2)
         imod=imodarr(3)
         call aero_470veg(refl,x1,x2,x3,mm,nn,ll,ma,
     1      imod,r470,tau_x4703,tau_x470_flag3)
         f_2=(tau_x4701-cth1)/(cth2-cth1)
         tau_x470=f_2*tau_x4703+(1.0-f_2)*tau_x4702
         tau_x470_flag=tau_x470_flag3
      elseif(tau_x4701.ge.cth2) then
         imod=imodarr(3)
         call aero_470veg(refl,x1,x2,x3,mm,nn,ll,ma,
     1      imod,r470,tau_x4703,tau_x470_flag3)
         tau_x470=tau_x4703
         tau_x470_flag=tau_x470_flag3
      else
         tau_x470=tau_x4701
         tau_x470_flag=tau_x470_flag1
c        tau_x470=-999.0
      endif
c;--------------------------------------------<
c.   original call for 470nm AOT retrieval (before 08/05/2011)
c.       imod=3
c.       call aero_470veg(refl,x1,x2,x3,mm,nn,ll,ma, ! commented out
c.   1       imod,r470,tau_x470,tau_x470_flag)       ! 08/05/2011

!c. @@@test@@@--> 5lines
      if(tau_x470.gt.5.or.tau_x470.lt.0) then
          tau_x470=-999.
          tau_x650=-999.
          return
      endif
c      print *,'tau_x470 =', tau_x470


c--------------------------------------------------------
c          Retrieving 650 nm AOT --> added Feb 26, 2010
c--------------------------------------------------------
c   Input surface reflectance at 650 nm
      r650 = xeAs_B1
c
c   Input toa reflectance (L/F) at 650 nm
      refl = refl6

c     Retrieving 650 nm AOT
      call aero_650veg(refl,x1,x2,x3,mm,nn,ll,ma,
     1    r650,tau_x650,tau_x650_flag,
     2    tau_x470_flag,tau_x470)

c      print *,'tau_x650 =', tau_x650, r650

!c. @@@test@@@--> 5lines
      if(tau_x650.gt.5.or.tau_x650.lt.0) then
          tau_x470=-999.
          tau_x650=-999.
          return
      endif
c--------------------------------------------------------

c   -- tweak AOT over S.Africa, Mongu site.
      if (ioprg .eq. 6) then
        if (month .ge. 6 .AND. month .le. 8) then
          tau_x470 = tau_x470 * 1.5
          tau_x650 = tau_x650 * 1.5
        end if
        if (month .ge. 9 .AND. month .le. 11) then
          tau_x470 = tau_x470 * 1.25
          tau_x650 = tau_x650 * 1.25
        end if
      end if

c   -- tweak over Mukdahan, SE Asia site.
      if (ioprg .eq. 7) then
        if (month .ge. 3 .AND. month .le. 5) then
          if (tau_x470 .gt. 0.7) then
            tau_x470 = tau_x470 * 1.1
            tau_x650 = tau_x650 * 1.1
          end if
        end if
        
        if (month .eq. 12 .OR. month .le. 2) then
          tau_x470 = tau_x470 * 1.25
          tau_x650 = tau_x650 * 1.25
        end if
        if (month .ge. 9 .AND. month .le. 11) then
          tau_x470 = tau_x470 * 1.1
          tau_x650 = tau_x650 * 1.1
        end if
      end if
      
c... original set up (temporary)
c...==============================================
c... MJ re-defines output
      alpha=alog10(tau_x470/tau_x650)/alog10(650./466.)
c     if(tau_x470.lt.0.0) tau_x470=-999.0
c     if(tau_x650.lt.0.0) tau_x650=-999.0
c     if(tau_x470.lt.0.0.and.tau_x650.lt.0.0) alpha=-999.0
c     if(alpha.ge.-0.5.and.alpha.le.3.5.and.tau_x470.gt.0.and.
c    &   tau_x650.gt.0.) then
      if(alpha.ge.-0.5.and.alpha.le.3.5) then
         tau550=tau_x470*(466.0/550.0)**alpha
c        xtau(1)=tau_x470*(466.0/412.0)**alpha  ! AOT 412nm (extrapol.)
      else
         if((tau_x470_flag.eq.1.or.tau_x470_flag.eq.2).and.
     1      (tau_x650_flag.eq.1.or.tau_x650_flag.eq.2)) then    
            alpha=1.0
            tau550=(tau_x470+tau_x650)/2.0
            tau_x470=tau550*(550./466.)**alpha
            tau_x650=tau550*(550./650.)**alpha
         elseif((tau_x470_flag.eq.1.or.tau_x470_flag.eq.2).and.
     1           tau_x650_flag.eq.0) then
            alpha=1.0
            tau550=tau_x650*(650./550.)**alpha
            tau_x470=tau_x650*(650./466.)**alpha
         elseif(tau_x470_flag.eq.0.and.
     1          (tau_x650_flag.eq.1.or.tau_x650_flag.eq.2)) then
            alpha=1.0
            tau550=tau_x470*(466./550.)**alpha
            tau_x650=tau_x470*(466./650.)**alpha
         else
            alpha=1.0
            tau550=(tau_x470+tau_x650)/2.0
            tau_x470=tau550*(550./466.)**alpha
            tau_x650=tau550*(550./650.)**alpha
         endif
      endif
c...==============================================

c...==========================================================
c... MJ re-defines output
c     if(tau_x470.lt.0.0) tau_x470=-999.0
c     if(tau_x650.lt.0.0) tau_x650=-999.0
c     if(tau_x470.lt.0.0.or.tau_x650.lt.0.0) return
c-------------------------------------------------------------
c Set output buffer --> based on find_v.f
c-------------------------------------------------------------
      read(w0_name470(imod),'(f4.2)') ssatmp
      ssa(2)=ssatmp
      ssa(3)=0.976
      xtau(1)=-999.0  ! fillvalue at 412nm for over-veg retr.
      xtau(2)=tau_x470
      xtau(3)=tau_x650
      do i=1,3
         outbufvg(i) = xtau(i)
         outbufvg(i+3) = ssa(i)
      enddo
      outbufvg(7) = tau550
      outbufvg(8) = alpha
      outbufvg(9) =  1.0*tau_x470_flag  ! used to be "r412"
      outbufvg(10) = 1.0*tau_x650_flag
      outbufvg(11) = r470   ! over-veg. sfc. refl.
      outbufvg(12) = r650   ! over-veg. sfc. refl.
      outbufvg(13) = xthet
      outbufvg(14) = scat_ang
      outbufvg(15) = sfc_typ

      outbufvg(16) = xlat
      outbufvg(17) = xlong
      outbufvg(18) = sirndvi
      outbufvg(19) = s_ndvi
      outbufvg(20) = xreg_id

c     print *,'outbufvg(1,2,3) ',outbufvg(1),outbufvg(2),outbufvg(3)
      return
      end
c... --------------------------------------
c... The end of subroutine myd_aero_veg
c... (main subroutine for aerosol retr. over vegetation)
c... ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



c... ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine get_sfcrfl_veg(iopss,refl21,xnvalm6,lcvr,
     1                          xeAs_B1,xeAs_B3,sirndvi,s_ndvi)
c...
c...  subroutine to estimate surface reflectance at 470nm & 650nm over
c...  vegetated surfaces using 2.1um TOA reflectance
c...     program written by Myeong-Jae Jeong (MJ)
c...     last modified Jun 15, 2010
c...
c...     Inputs
c...        iopss:   season ID (1: MAM; 2: JJA; 3: SON; 4: DJF)
c...        refl21:  2.1um TOA reflectance
c...        xnvalm6: TOA reflectances at other bands
c...        lcvr:    MODIS land cover (IGBP)
c...
c...     Output
c...        xeAs_B1: estimated surface refl. for band 1 (LER in %)
c...        xeAs_B3: estimated surface refl. for band 3 (LER in %)
c...

      real xnvalm6(6)

      sirndvi=(xnvalm6(6)-xnvalm6(4))/(xnvalm6(6)+xnvalm6(4))
      s_ndvi =(xnvalm6(5)-xnvalm6(1))/(xnvalm6(5)+xnvalm6(1))

c...  initialization
      xeAs_B1=-999.0
      xeAs_B3=-999.0

c      print *, 'lcvr: ',lcvr
c... ----------------------------------------------------------
c... Get surface reflectance over non-arid surfaces (MAM2004)
c... ----------------------------------------------------------
      sb7=100.0*refl21
      if(refl21.lt.0.01) return    ! skip very small and (-)ve B7refl, no retr.
      if(s_ndvi.lt.0.10) return    ! NDVI' <0.1, no retr.
c...
c...  ---------------------- MAM 2004 ------------------------
      if(iopss.eq.1) then  ! MAM 2004
      if((lcvr.gt.0.and.lcvr.lt.12).or.lcvr.eq.14) then
         xeAs_B1=0.552625 + 0.480117*sb7 + 0.003773*sb7**2
         sb1=xeAs_B1
         xeAs_B3=-0.330536 + 0.482954*sb1 + 0.000000*sb1**2
      elseif(lcvr.eq.12) then
         if(sirndvi.lt.0.35) then
            xeAs_B1=6.282806 + 0.165798*sb7 + 0.000000*sb7**2
            sb1=xeAs_B1
            xeAs_B3=2.688409 + 0.275116*sb1 + 0.000000*sb1**2
         else
            xeAs_B1=-0.976568 + 0.621287*sb7 + 0.000000*sb7**2
            sb1=xeAs_B1
            xeAs_B3=0.912642 + 0.398208*sb1 + 0.000000*sb1**2
         endif
      elseif(lcvr.eq.13) then
         xeAs_B1=-3.685365 + 0.941037*sb7 + 0.000000*sb7**2
         sb1=xeAs_B1
         xeAs_B3=-0.349138 + 0.634389*sb1 + 0.000000*sb1**2
      else
         return
      endif
c...
c...  ---------------------- JJA 2004 ------------------------
      elseif(iopss.eq.2) then  ! JJA 2004
      if((lcvr.gt.0.and.lcvr.lt.12).or.lcvr.eq.14) then
         xeAs_B1=0.441294 + 0.460608*sb7 + 0.004524*sb7**2
         sb1=xeAs_B1
         xeAs_B3=-0.584121 + 0.496080*sb1 + 0.000000*sb1**2
      elseif(lcvr.eq.12) then
         if(sirndvi.lt.0.35) then
            xeAs_B1= 5.239523+ 0.207669*sb7 + 0.000000*sb7**2  !sg1
            sb1=xeAs_B1
            xeAs_B3= 0.245139 + 0.544203*sb1 + 0.000000*sb1**2
         else
            xeAs_B1=-0.118660 + 0.503617*sb7 + 0.000000*sb7**2 !sg2
            sb1=xeAs_B1
            xeAs_B3=-0.073580 + 0.534508*sb1 + 0.000000*sb1**2
         endif
      elseif(lcvr.eq.13) then
         xeAs_B1=-1.702984 + 0.830312*sb7 + 0.000000*sb7**2
         sb1=xeAs_B1
         xeAs_B3=-0.050284 + 0.623932*sb1 + 0.000000*sb1**2
      else
         return
      endif

c...
c...  ---------------------- SON 2004 ------------------------
      elseif(iopss.eq.3) then  ! SON 2004

      if((lcvr.gt.0.and.lcvr.lt.12).or.lcvr.eq.14) then
         xeAs_B1=1.174910 + 0.356037*sb7 + 0.006684*sb7**2
         sb1=xeAs_B1
         xeAs_B3= 0.004815 + 0.442943*sb1 + 0.000000*sb1**2
      elseif(lcvr.eq.12) then
         if(sirndvi.lt.0.35) then
            xeAs_B1=-2.264160 + 0.678071*sb7 + 0.000000*sb7**2  !sg1
            sb1=xeAs_B1
            xeAs_B3= 1.249338 + 0.357649*sb1 + 0.000000*sb1**2
         else
            xeAs_B1=-1.279875 + 0.616125*sb7 + 0.000000*sb7**2  !sg2
            sb1=xeAs_B1
            xeAs_B3= 1.272351 + 0.203916*sb1 + 0.000000*sb1**2
         endif
      elseif(lcvr.eq.13) then
         xeAs_B1=-1.516383 + 0.801742*sb7 + 0.000000*sb7**2
         sb1=xeAs_B1
         xeAs_B3=-0.058548 + 0.604228*sb1 + 0.000000*sb1**2
      else
         return
      endif

c...  ---------------------- DJF 2004 ------------------------
      elseif(iopss.eq.4) then  ! DJF 2004

c        ! data not ready... use MAM relations for the time being

      else
c
         return
c
      endif  !...

c     if(xeAs_B1.lt.0.0.or.xeAs_B3.lt.0.0) go to 10  ! no retr.

      return
      end
c... --------------------------------------
c... The end of subroutine get_sfcrfl_veg
c... ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


c... ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c--------------------------------------------------------
      subroutine aero_470veg(refl,x1,x2,x3,mm,nn,ll,ma,
     1    imod,r470,tau_x470,tau_x470_flag)

      include 'aottbl.inc'
      include 'newaottbl.inc'

c      common /angle_node/ theta0, theta, phi
c      common /sfcref_node/ sfc_ref412, sfc_ref470, sfc_ref650
c      real theta0(10), theta(46), phi(30), tau(10)
c      real sfc_ref412(20), sfc_ref470(24), sfc_ref650(24)
      real nnvalx(4,4,2,10), yy(10), yy2(8)
c      real nvalx470(10,46,30,10,4,24)
      integer tau_x470_flag, imod
      data pi   /3.14159/
c      data tau  /0.0, 0.1, 0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5/
      
c      print *, 'aero_470veg in: ', refl,x1,x2,x3,mm,nn,ll,ma,imod,r470,tau_x470,tau_x470_flag
      index_ii = r470
      if (index_ii.lt.0) return  ! orig.

      frac = (r470-sfc_ref470(index_ii))/
     1       (sfc_ref470(index_ii+1)-sfc_ref470(index_ii))

      if (index_ii.lt.1.or.index_ii.gt.24)
     1    print *,'index_iir470 = ', index_ii,xlat,xlong
      if (frac.lt.0.0.or.frac.gt.1.0)
     1    print *,'frac on sfc470=', frac
      if (index_ii.lt.1.or.index_ii.gt.24) then  ! MJ added
        tau_x470_flag = 9
        tau_x470=-999.0
        return             
      endif


      call search(dflag2,x3,phi,ll,ii)
      xfrac = (x3-phi(ii))/(phi(ii+1)-phi(ii))

      nsm=1
      dif=x1-theta0(1)
      do i=1,mm
        dift=x1-theta0(i)
        if (dift.gt.0. .and. dift.lt.dif) then
          nsm=i
          dif=dift
        endif
      enddo
      mbeg = nsm - 2
      if (mbeg.le.0) then
        mbeg = 0
      else if (mbeg.gt.mm-4) then
        mbeg = mm-4
      endif

      nsn=1
      dif=x2-theta(1)
      do i=1,nn
        dift=x2-theta(i)
        if (dift.gt.0. .and. dift.lt.dif) then
          nsn=i
          dif=dift
        endif
      enddo
      nbeg = nsn - 2
      if (nbeg.le.0) then
        nbeg = 0
      else if (nbeg.gt.nn-4) then
        nbeg = nn-4
      endif

c      write(6,*) "frac =",frac
      do ia = 1, 10
       do i = 1, 4
        do j = 1, 4
       nnvalx(i,j,1,ia) = nvalx470(mbeg+i,nbeg+j,ii,ia,imod,index_ii)*
     1  (1.-frac) + nvalx470(mbeg+i,nbeg+j,ii,ia,imod,index_ii+1)*frac
       nnvalx(i,j,2,ia) = nvalx470(mbeg+i,nbeg+j,ii+1,ia,imod,index_ii)*
     1  (1.-frac) + nvalx470(mbeg+i,nbeg+j,ii+1,ia,imod,index_ii+1)*frac
        enddo
       enddo
      enddo

c---     interpolating AOT tables

      do 105 ia = 1, 10

      call new_intep(theta0, theta, phi, nnvalx, mm, nn, ll, ia,
     1          x1,x2,x3,y,dy,mbeg,nbeg,xfrac)

      yy(ia) = y/pi
c      print *,'tau, i/f=', tau(ia), y/pi, dy
 105  continue
      if (refl.le.yy(1)) then
       tau_x470 = 0.02
       tau_x470_flag = 2 ! MJ added
       return
      endif

      if (refl.ge.yy(10)) then
       xxrat = 0.8
       tau_x470_flag = 1
       return
      endif
c
c
c     Check if the reflectance increase with AOT
c

      if (yy(1).lt.yy(2)) go to 650

      if (refl.lt.yy(4)) return

      do i = 1, 7
      yy2(i) = yy(i+3)
      enddo

      if (yy2(2).lt.yy2(1)) return
      call search2(dflag2,refl,yy2,7,index_ii,frac)
c      print *,'after 1st search'
c      if (yy2(1).gt.yy2(2)) print *,'yy=',refl,(yy(i),I=1,10)
      tau_x = frac*tau(index_ii+1+3) + (1.-frac)*tau(index_ii+3)
      tau_x470 = tau_x
      return

650   continue
c
c     Pass the monotonic order check
c
      call search2(dflag2,refl,yy,10,index_ii,frac)
c      print *,'after 2nd search'
      tau_x = frac*tau(index_ii+1) + (1.-frac)*tau(index_ii)

      tau_x470 = tau_x

c      print *,'tau_x470 =', tau_x470

      return
      end
c... --------------------------------------
c... The end of subroutine aero_470veg
c... ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



c... ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c--------------------------------------------------------
      subroutine aero_650veg(refl,x1,x2,x3,mm,nn,ll,ma,
     1    r650,tau_x650,tau_x650_flag,
     2    tau_x470_flag,tau_x470)

      include 'aottbl.inc'
      include 'newaottbl.inc'

c      common /angle_node/ theta0, theta, phi
c      common /sfcref_node/ sfc_ref412, sfc_ref470, sfc_ref650

c      real theta0(10), theta(46), phi(30), tau(10)
c      real sfc_ref412(20), sfc_ref470(24), sfc_ref650(24)
      real nnvalx(4,4,2,10), yy(10), yy2(8), yy3(3), yy5(6)
      real tau_x650, tau_x412, tau_x470
      real refl,x1,x2,x3,r650
c      real nvalx650(10,46,30,10,24)
      integer tau_x650_flag, tau_x412_flag_91
      logical dflag2
      data pi   /3.14159/
c      data tau  /0.0, 0.1, 0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5/

      index_ii = (r650+1.)/2.

      frac = (r650-sfc_ref650(index_ii))/
     1       (sfc_ref650(index_ii+1)-sfc_ref650(index_ii))

c      if (index_ii.lt.1.or.index_ii.gt.24)
c     1    print *,'index_iir650 = ', index_ii,xlat,xlong
c      if (frac.lt.0.0.or.frac.gt.1.0)
c     1    print *,'frac on sfc=', frac

c... --------------------------------------------------
c...  MJ added 12Feb2010 @@@@@
c... --------------------------------------------------
c.    if(index_ii.lt.2.or.index_ii.gt.12) then ! ok ver.
      if(index_ii.lt.1.or.index_ii.gt.23) then
         tau_x650_flag = 9 
         tau_x650 = -999.0
         return
      endif ! if(index_ii.lt.1.or.index_ii.gt.23) then
c... --------------------------------------------------

      if (index_ii.lt.1) then
        index_ii = 1
        frac = 0.
      endif

      call search(dflag2,x3,phi,ll,ii)
      xfrac = (x3-phi(ii))/(phi(ii+1)-phi(ii))

      nsm=1
      dif=x1-theta0(1)
      do i=1,mm
        dift=x1-theta0(i)
        if (dift.gt.0. .and. dift.lt.dif) then
          nsm=i
          dif=dift
        endif
      enddo
      mbeg = nsm - 2
      if (mbeg.le.0) then
        mbeg = 0
      else if (mbeg.gt.mm-4) then
        mbeg = mm-4
      endif

      nsn=1
      dif=x2-theta(1)
      do i=1,nn
        dift=x2-theta(i)
        if (dift.gt.0. .and. dift.lt.dif) then
          nsn=i
          dif=dift
        endif
      enddo
      nbeg = nsn - 2
      if (nbeg.le.0) then
        nbeg = 0
      else if (nbeg.gt.nn-4) then
        nbeg = nn-4
      endif

c      write(6,*) "frac =",frac
      do ia = 1, 10
       do i = 1, 4
        do j = 1, 4
          nnvalx(i,j,1,ia) = nvalx650(mbeg+i,nbeg+j,ii,ia,index_ii)*
     1       (1.-frac) + nvalx650(mbeg+i,nbeg+j,ii,ia,index_ii+1)*frac
          nnvalx(i,j,2,ia) = nvalx650(mbeg+i,nbeg+j,ii+1,ia,index_ii)*
     1       (1.-frac) + nvalx650(mbeg+i,nbeg+j,ii+1,ia,index_ii+1)*frac
        enddo
       enddo
      enddo

c---     interpolating AOT tables

      do 600 ia = 1, 10

      call new_intep(theta0, theta, phi, nnvalx, mm, nn, ll, ia,
     1          x1,x2,x3,y,dy,mbeg,nbeg,xfrac)

      yy(ia) = y/pi
c      print *,'tau, i/f=', tau(ia), y/pi, dy
 600  continue


      if (refl.le.yy(1).and.yy(1).lt.yy(2)) then
       tau_x650_flag = 2  ! MJ added
       tau_x650 = 0.02
       return
      endif

      if (refl.ge.yy(10)) then
      tau_x650 = 4.0
      w0_x = -999.
      tau_x650_flag = 1
      return
      endif

c!!!  if (tau_x470_flag.gt.0) go to 670 ! ori open 20100216 @@@@@
c      if (tau_x412_flag_91.gt.0) go to 680
c
c
c     Check if the reflectance increase with AOT
c
      if (yy(1).lt.yy(2)) go to 650

c... *****************************
c... *****************************
c... *****************************
        tau_x650_flag = 4  ! mj added 08/05/2011
        return             ! mj added 20100216 @@@@@
                           ! no retrieval for non-monotonic
                           ! changes in LUT reflectance
c... needs investigation to comment this out!
c... *****************************
c... *****************************
c... *****************************

      if (refl.lt.yy(4)) return

      do i = 1, 7
      yy2(i) = yy(i+3)
      enddo

      if (yy2(2).lt.yy2(1)) return
      call search2(dflag2,refl,yy2,7,index_ii,frac)
c      print *,'after 1st search'
c      if (yy2(1).gt.yy2(2)) print *,'yy=',refl,(yy(i),I=1,10)
      tau_x = frac*tau(index_ii+1+3) + (1.-frac)*tau(index_ii+3)
      tau_x650 = tau_x
      return

670   continue
      if (refl.lt.yy(8)) return

      do i = 1, 3
      yy3(i) = yy(i+7)
      enddo

      if (yy3(2).lt.yy3(1)) return
      call search2(dflag2,refl,yy3,3,index_ii,frac)
c      print *,'after 1st search'
c      if (yy3(1).gt.yy3(2)) print *,'yy=',refl,(yy(i),I=1,10)
      tau_x = frac*tau(index_ii+1+7) + (1.-frac)*tau(index_ii+7)
      tau_x650 = tau_x
      return

680   continue
      if (refl.lt.yy(5)) return

      do i = 1, 6
      yy5(i) = yy(i+4)
      enddo

      if (yy5(2).lt.yy5(1)) return
      call search2(dflag2,refl,yy5,6,index_ii,frac)
c      print *,'after 1st search'
c      if (yy3(1).gt.yy3(2)) print *,'yy=',refl,(yy(i),I=1,10)
      tau_x = frac*tau(index_ii+1+4) + (1.-frac)*tau(index_ii+4)
      tau_x650 = tau_x
      return

650   continue
c
c     Pass the monotonic order check
c
      call search2(dflag2,refl,yy,10,index_ii,frac)
c      print *,'after 2nd search'
      tau_x = frac*tau(index_ii+1) + (1.-frac)*tau(index_ii)

      tau_x650 = tau_x

      return
      end
c... --------------------------------------
c... The end of subroutine aero_650veg
c... ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

c **********************************************************************************************************


c **********************************************************************************************************
      subroutine cmb_stddb_veg(outbuf, outbufvg, tmpvg, idalg)
c... 
c...  A subroutine to combine standard Deep Blue and over-vegetation
c...  retrievals, depending on ndvi, ndvi_swir, refl at 2.1um, etc.
c...     Written by Myeong-Jae Jeong (MJ)
c...     Ver. 0.5 Aug 8, 2011
c...
c.... Inputs:
c...     outbuf(20)   - standard DeepBlue output
c...     outbufvg(20) - over-vegetation retr. output
c...     tmpvg(6)     - TOA reflectances
c...
c...  Output:
c...     outfub(20)   - combined output
c...     idalg        - index for algorithm in effect
c...
c...     idalg=0 : standard Deep Blue 
c...     idalg=1 : over-vegetation retrieval
c...     idalg=2 : mixture of standard DB and over-veg. retrievals

      real outbuf(20), outbufvg(20), tmpvg(6)
      real outbuftmp(20), xfrdb
      real c_ndsir1,c_ndvi1,c_rf21_1,c_rf21_2,c_rf21_3
      integer idalg

      c_ndsir1=0.1  ! NDVI_swir minimum threshold  
      c_ndvi1=0.1   ! NDVI minimum threshold 
      c_rf21_1=0.01 ! 2.1um refl. threshold 1
      c_rf21_2=0.25 ! 2.1um refl. threshold 2
      c_rf21_3=0.35 ! 2.1um refl. threshold 3

c...  Note. Also try to employ NDVI_swir thresholds later to combine
c...        standard DeepBlue and over-vegeation retrievals... (MJ)

c...  Initialization ... set to standard DeepBlue retrieval results
      do i=1,20
         outbuftmp(i)=outbuf(i)
      enddo
      idalg=0 ! standard DeepBlue 


c      if(outbufvg(18).ge.c_ndsir1.and.outbufvg(19).ge.c_ndvi1) then
c         if(tmpvg(4).ge.c_rf21_1.and.tmpvg(4).lt.c_rf21_2) then  ! veg.
c            do i=1,8
c               outbuftmp(i)=outbufvg(i)
c            enddo
c            do i=11,14
c               outbuftmp(i)=outbufvg(i)
c            enddo
c            outbuftmp(15) = 100.0
c            idalg=1
c         elseif(tmpvg(4).ge.c_rf21_2.and.tmpvg(4).lt.c_rf21_3) then
cc            xfrdb=(tmpvg(4)-c_rf21_2)/(c_rf21_3-c_rf21_2)     
cc            if(outbufvg(7).ge.0.and.outbuf(7).ge.0) then      ! veg+DB
cc               do i=2,3
cc                  outbuftmp(i)=(1.0-xfrdb)*outbufvg(i)+xfrdb*outbuf(i)
cc                  outbuftmp(i+3)=((1.0-xfrdb)*outbufvg(i+3)*outbufvg(i) 
cc     1                           + xfrdb*outbuf(i+3)*outbuf(i)) 
cc     2                           /(outbuftmp(i))
cc               enddo
cc               outbuftmp(7)=(1.0-xfrdb)*outbufvg(7)+xfrdb*outbuf(7)
cc               outbuftmp(1)=-999.0  ! set aot412nm to -999.0
cc               outbuftmp(4)=-999.0  ! set ssa412nm to -999.0
cc               outbuftmp(8)=alog(outbuftmp(2)/outbuftmp(3))/
cc     1                           alog(650./466.)
cc               if(outbufvg(11).gt.0.and.outbuf(11).gt.0) then
cc                  outbuftmp(11)=((1.0-xfrdb)*outbufvg(11)*outbufvg(2)
cc     1                          +xfrdb*outbuf(11)*outbuf(2))
cc     2                          /(outbuftmp(2))       
cc                  outbuftmp(12)=((1.0-xfrdb)*outbufvg(12)*outbufvg(3)
cc     1                          +xfrdb*outbuf(12)*outbuf(3))
cc     2                          /(outbuftmp(3))       
cc               else 
cc                  outbuftmp(11)=-999.0
cc                  outbuftmp(12)=-999.0
cc               endif
cc               idalg=2
c            if(outbufvg(7).ge.0.and.outbuf(7).lt.0) then
c               do i=1,8
c                  outbuftmp(i)=outbufvg(i)
c               enddo
c               do i=11,14
c                  outbuftmp(i)=outbufvg(i)
c               enddo
c               outbuftmp(15) = 100.0
c               idalg=1
c            else
c               do i=1,8
c                  outbuftmp(i)=outbuf(i)
c               enddo
c               do i=11,14
c                  outbuftmp(i)=outbuf(i)
c               enddo
c               idalg=0
c            endif
c         endif
c      endif 

c...  Finalization ... push the results to outbuf
      do i=1,20
         outbuf(i)=outbuftmp(i)
      enddo

      return
      end
c... --------------------------------------
c... The end of subroutine cmb_stddb_veg
c... ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c **********************************************************************************************************


c **********************************************************************************************************
      subroutine calc_sfc21(sza,xthet,xphi,refl21,aot21,As_21)
c...
c...  ----------------------------------------------------------
c...  written by Myeong-Jae Jeong (MJ) at GWNU, Korea
c...  ver. 1.0 Aug 12, 2011
c...
c...  subroutine to calculate surface reflectance at 2.1um.
c...  code inherited from "sfc470_5wav.f" by N.C. Hsu 
c...
c...  Inputs
c...     sza   : solar zenith angle
c...     xthet : viewing zenith angle
c...     xphi  : relative azimuth angle
c...     refl21: 2.1um TOA reflectance
c...     aot21 : 2.1um AOT
c... 
c...  Output
c...     As_21 : 2.1 surface reflectance (%)
c...
c...  ----------------------------------------------------------
c...     refl21  = refl21  * xmu/3.14159

      include 'aottbl.inc'
      include 'sfc21tbl.inc'

c     real tau2(4), theta0(10), theta(46), phi(30)
      real tau2(4)
      real nnvalx(4,4,2), rr0x(4,4,2)
      real ttx(4,4,2), ssx(4,4,2)
      logical dflag3
      integer index_ii

      data tau2    /0.0, 0.1, 0.3, 0.5/
      data pi      /3.14159/
c     data theta0  /0.0,8.0,16.0,24.0,32.0,40.0,48.0,56.0,64.0,72.0/

      x1 = sza
      x2 = xthet
      x3 = xphi
      refl7 = refl21        ! 2.1 micron

c     do i = 1, 46
c        theta(i) = 2.*float(i-1)
c     enddo

c     do i = 1, 30
c        phi(i) = 6. + 6.*float(i-1)
c     enddo

      mm = 10     ! solar zenith
      nn = 46     ! satellite zenith
      ll = 30     ! rel. azimuth
      ma2=  4     ! tau2 (aot at 2.1um)

c     -- Derive SFC

c--   interpolating for 2.1 micron channel

      nsm=1
      dif=x1-theta0(1)
      do i=1,mm
        dift=x1-theta0(i)
        if (dift.gt.0. .and. dift.lt.dif) then
          nsm=i
          dif=dift
        endif
      enddo
      mbeg = nsm - 2
      if (mbeg.le.0) then
        mbeg = 0
      else if (mbeg.gt.mm-4) then
        mbeg = mm-4
      endif
 
      nsn=1
      dif=x2-theta(1)
      do i=1,nn
        dift=x2-theta(i)
        if (dift.gt.0. .and. dift.lt.dif) then
          nsn=i
          dif=dift
        endif
      enddo
      nbeg = nsn - 2
      if (nbeg.le.0) then
        nbeg = 0
      else if (nbeg.gt.nn-4) then
        nbeg = nn-4
      endif
c... 
      aot2=aot21
      call search(dflag3,aot2,tau2,ma2,index_ii)
      frac  = (aot2-tau2(index_ii))
     1        /(tau2(index_ii+1)-tau2(index_ii))
      if(frac.lt.0.or.frac.gt.1)
     1   print *,'aot2, frac, index_ii=',aot2, frac, index_ii
 
      call search(dflag3,x3,phi,ll,ii)
      xfrac = (x3-phi(ii))/(phi(ii+1)-phi(ii))

c     print *, 'nsm, nsn, ii, aot2, frac, index_ii'
c     print *, nsm, nsn, ii, aot2, frac, index_ii

       do 123 i = 1, 4
        do 124 j = 1, 4
          nnvalx(i,j,1) = nvalx21(mbeg+i,nbeg+j,ii,index_ii)*
     1       (1.-frac) + nvalx21(mbeg+i,nbeg+j,ii,index_ii+1)*frac   
          nnvalx(i,j,2) = nvalx21(mbeg+i,nbeg+j,ii+1,index_ii)*
     1       (1.-frac) + nvalx21(mbeg+i,nbeg+j,ii+1,index_ii+1)*frac   
 
          rr0x(i,j,1) = r0x_21(mbeg+i,nbeg+j,ii,index_ii)*
     1       (1.-frac) + r0x_21(mbeg+i,nbeg+j,ii,index_ii+1)*frac
          rr0x(i,j,2) = r0x_21(mbeg+i,nbeg+j,ii+1,index_ii)*
     1       (1.-frac) + r0x_21(mbeg+i,nbeg+j,ii+1,index_ii+1)*frac
 
          ttx(i,j,1) = tx_21(mbeg+i,nbeg+j,ii,index_ii)*
     1       (1.-frac) + tx_21(mbeg+i,nbeg+j,ii,index_ii+1)*frac
          ttx(i,j,2) = tx_21(mbeg+i,nbeg+j,ii+1,index_ii)*
     1       (1.-frac) + tx_21(mbeg+i,nbeg+j,ii+1,index_ii+1)*frac
 
          ssx(i,j,1) = sx_21(mbeg+i,nbeg+j,ii,index_ii)*
     1       (1.-frac) + sx_21(mbeg+i,nbeg+j,ii,index_ii+1)*frac
          ssx(i,j,2) = sx_21(mbeg+i,nbeg+j,ii+1,index_ii)*
     1       (1.-frac) + sx_21(mbeg+i,nbeg+j,ii+1,index_ii+1)*frac

124    continue
123   continue
c       enddo
c      enddo
 
c---     
 
        call new_intepsf(theta0, theta, phi, nnvalx, mm, nn, ll,
     1                 x1,x2,x3,y,dy,mbeg,nbeg,xfrac)
 
        RR = y/pi
 
        call new_intepsf(theta0, theta, phi, rr0x, mm, nn, ll,
     1                 x1,x2,x3,y,dy,mbeg,nbeg,xfrac)
 
        RR0 = y/pi
 
        call new_intepsf(theta0, theta, phi, ttx, mm, nn, ll,
     1                 x1,x2,x3,y,dy,mbeg,nbeg,xfrac)
 
        TT = y/pi
 
        call new_intepsf(theta0, theta, phi, ssx, mm, nn, ll,
     1                 x1,x2,x3,y,dy,mbeg,nbeg,xfrac)
 
c!      SS = y/pi
        SS = y    ! chg20090818
 
        dff = refl7 - RR0
        As_21  = 100.* dff / (TT + SS*dff)

c        print *,nday,sza,xthet,xphi,refl3,refl7,As,As_21
c        print *,'refl7,RR0,dff,TT,SS,As_21= ',
c     1     refl7,RR0,dff,TT,SS,As_21

      return
      end
c... --------------------------------------
c... The end of subroutine calc_sfc21
c... ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c **********************************************************************************************************

c
      subroutine new_intepsf(x1a,x2a,x3a,ya,m,n,l,x1,x2,x3,y,dy,
     1                     mbeg,nbeg,frac)
      parameter (nmax=46,mmax=10,lmax=30)
      dimension x1a(m),x2a(n),x3a(l),ya(4,4,2)
      dimension xx2a(4), xx1a(4)
      dimension yntmp(4),ymtmp(4),yltmp(2)
 
      do 12 j=1,4
        do 11 k=1,4
          yltmp(1)=ya(j,k,1)
          yltmp(2)=ya(j,k,2)
          yntmp(k) = yltmp(1)*(1.-frac) + yltmp(2)*frac
          xx2a(k) = x2a(k+nbeg)
11      continue
        call polint(xx2a,yntmp,4,x2,ymtmp(j),dy)
        xx1a(j) = x1a(j+mbeg)
12    continue
      call polint(xx1a,ymtmp,4,x1,y,dy)
      return
      end

