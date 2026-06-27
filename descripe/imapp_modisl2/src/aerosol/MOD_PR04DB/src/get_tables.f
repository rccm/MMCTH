      subroutine get_tables(fhandle)
c
c     Load AOT Tables
c
      implicit none

      include 'aottbl.inc'
      include 'newaottbl.inc'
      include 'table.inc'
      include 'PGS_MODIS_39500.f'
      INCLUDE 'PGS_SMF.f'
      INCLUDE 'PGS_PC.f'
      INCLUDE 'PGS_PC_9.f'
      INCLUDE 'PGS_IO.f'
      INCLUDE 'PGS_IO_1.f'
      
      integer LUN_TERPRS
      integer LUN_MODIS_TABLE
      integer LUN_MODIS_FILE
C      integer LUN_CORR_MODIS_FILE
      integer LUN_412_AER_TBL
      integer LUN_466_AER_TBL
      integer LUN_645_AER_TBL
      integer LUN_DOY

      parameter (LUN_MODIS_TABLE=412000)
      parameter (LUN_TERPRS=412001)
      parameter (LUN_412_AER_TBL=412002)
      parameter (LUN_466_AER_TBL=412003)
      parameter (LUN_645_AER_TBL=412004)
      parameter (LUN_MODIS_FILE=412200)
C      parameter (LUN_CORR_MODIS_FILE=412201)
      parameter (LUN_DOY=405037)
			real*4 ps(360,720)
c      real*4 gref412(3600,1800), gref470(3600,1800), 
c      real*4 gref650(3600,1800)
      integer*4 i, j, k,fv, reclen, fh, fhandle, fhandle_corr
      integer*4 prtn, pgs_io_gen_openf, pgs_io_gen_closef
      integer*4 PGS_PC_getconfigdata
      integer*4 gap_beg(3), gap_end(3), doy
      real*4 drat
      
      character pathname*256, msg*256, CheckStr*512, input_file
            
      common   /surpre/ ps
      common   /xday/ doy
      data theta0  /0.0,8.0,16.0,24.0,32.0,40.0,48.0,56.0,64.0,72.0/
      data w0      /0.82, 0.85, 0.87, 0.91, 0.94, 0.96, 0.98, 0.995/
      data w0_470  /0.91, 0.94, 0.96, 0.99/
      data tau     /0.0, 0.1, 0.3, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5/
      data sfc_ref412 /1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,
     1                 11.,12.,13.,14.,15.,16.,17.,18.,19.,20./
      data sfc_ref470 /1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,
     1                 11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,
     2                 21.,22.,23.,24./
      data sfc_ref650 /1.,3.,5.,7.,9.,11.,13.,15.,17.,19.,21.,
     1                 23.,25.,27.,29.,31.,33.,35.,37.,39.,41.,
     2                 43.,45.,47./

      data gap_beg /60, 182, 305/
      data gap_end /90, 212, 334/
      
c-----------------------------------------------------------------------
c     Get Day of Year from PCF and calculate day ratio for interpolation
c-----------------------------------------------------------------------

      prtn = PGS_PC_getconfigdata(LUN_DOY, CheckStr)
      if (prtn .ne. 0) then
         WRITE (msg,'(A)') "Failed to read Day of Year from PCF"
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'get_tables')
      ENDIF
                                                                                
      READ (CheckStr,'(I3)') doy

      drat = 0.0
      do i = 1,3
         if (doy .ge. gap_beg(i) .and. doy .le. gap_end(i)) then
            drat = (1.0*(doy - gap_beg(i)))/(gap_end(i) - gap_beg(i))
         endif
      enddo
      
c      fv=1
c      prtn=pgs_io_gen_openf(LUN_TERPRS,PGSd_IO_Gen_RseqUnf,reclen,fh,fv)
c      if (prtn .ne. 0) then
c         msg = "Error opening LUN_TERPRS lun."
c         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'get_tables')
c      endif
c      read (fh, end=900, err=900) sfcprs
c      prtn = pgs_io_gen_closef(fh)
      
c      do i = 1,360
c         do j = 1,720
c            ps(i,j) = sfcprs(i,j)
c         enddo
c      enddo

c-----------------------------------------------------------------------
c      Load Tables
c-----------------------------------------------------------------------

      fv=1
      reclen = 0
      prtn=pgs_io_gen_openf(LUN_MODIS_TABLE,PGSd_IO_Gen_RseqUnf,reclen,fh,fv)
      if (prtn .ne. 0) then
         msg = "Error opening LUN_MODIS_TABLE lun."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'get_tables')
      endif
      read (fh) logi0, z1i0, z2i0, ti0, sb
      read (fh) li0r, z1i0r, z2i0r, ti0r, sbr
      prtn = pgs_io_gen_closef(fh)

      fv = 1
      prtn=pgs_io_gen_openf(LUN_412_AER_TBL,PGSd_IO_Gen_RseqUnf,reclen,fh,fv)
      if (prtn .ne. 0) then
         msg = "Error opening LUN_412_AER_TBL lun."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'get_tables')
      endif
      msg = "opened LUN_412_AER_TBL lun."
      call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC, msg, 'get_tables')
      read (fh, end=900, err=900) nvalx412
      prtn = pgs_io_gen_closef(fh)

      fv = 1
      prtn=pgs_io_gen_openf(LUN_466_AER_TBL,PGSd_IO_Gen_RseqUnf,reclen,fh,fv)
      if (prtn .ne. 0) then
         msg = "Error opening LUN_466_AER_TBL lun."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'get_tables')
      endif
      msg = "opened LUN_466_AER_TBL lun."
      call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC, msg, 'get_tables')
      read (fh, end=900, err=900) nvalx470
      prtn = pgs_io_gen_closef(fh)

      fv = 1
      prtn=pgs_io_gen_openf(LUN_645_AER_TBL,PGSd_IO_Gen_RseqUnf,reclen,fh,fv)
      if (prtn .ne. 0) then
         msg = "Error opening LUN_645_AER_TBL lun."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'get_tables')
      endif
      msg = "opened LUN_645_AER_TBL lun."
      call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC, msg, 'get_tables')
      read (fh, end=900, err=900) nvalx650
      prtn = pgs_io_gen_closef(fh)

c-----------------------------------------------------------------------
c     End Loading Tables
c-----------------------------------------------------------------------

      do i = 1, 46
         theta(i) = 2.*float(i-1)
      enddo

      do i = 1, 30
         phi(i) = 6. + 6.*float(i-1)
      enddo

c-----------------------------------------------------------------------
c     Open intermediate data file
c-----------------------------------------------------------------------

      fv=1
      prtn=pgs_io_gen_openf(LUN_MODIS_FILE,PGSd_IO_Gen_WseqUnf,reclen,fhandle,fv)
      if (prtn .ne. 0) then
         msg = "Error opening LUN_MODIS_FILE lun."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'get_tables')
      endif
			
      return

900   continue
      msg = "Error reading file."
      call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'get_tables')
      end 

c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------

      subroutine interp_file(LUN, drat, first)
      implicit none

      include 'PGS_MODIS_39500.f'
      INCLUDE 'PGS_SMF.f'
      INCLUDE 'PGS_PC.f'
      INCLUDE 'PGS_PC_9.f'
      INCLUDE 'PGS_IO.f'
      INCLUDE 'PGS_IO_1.f'

      real    first(3600,1800)
      integer LUN
      integer reclen, fh
      integer fv, i, j, prtn, pgs_io_gen_openf, pgs_io_gen_closef, pgs_pc_getreference
      real    second(3600,1800), drat
      character msg*256
      
      reclen = 0
      fv = 1
      prtn = pgs_io_gen_openf(LUN, PGSd_IO_Gen_RseqUnf, reclen, fh, fv)
      if (prtn .ne. 0) then
         msg = "Error opening first file."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'interp_file')
      endif
      read (fh) first
      prtn = pgs_io_gen_closef(fh)

      if (drat .gt. 0) then
         fv = 2
         prtn = pgs_io_gen_openf(LUN, PGSd_IO_Gen_RseqUnf, reclen, fh, fv)
         if (prtn .ne. 0) then
            msg = "Error opening second file."
            call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'interp_file')
         endif
         read (fh) second
         prtn = pgs_io_gen_closef(fh)

         do j = 1,1800
            do i = 1,3600
               first(i,j) = first(i,j) + (second(i,j)-first(i,j))*drat
            enddo
         enddo
      endif

      end
