      subroutine get_newtables
c
c     Load AOT Tables added/changed for Col 5.1
c
      implicit none

      include 'newaottbl.inc'
      include 'table.inc'
      include 'PGS_MODIS_39500.f'
      INCLUDE 'PGS_SMF.f'
      INCLUDE 'PGS_PC.f'
      INCLUDE 'PGS_PC_9.f'
      INCLUDE 'PGS_IO.f'
      INCLUDE 'PGS_IO_1.f'

C      integer LUN_R412_REF
C      integer LUN_R470_REF
C      integer LUN_R470b_REF
C      integer LUN_R470_REF2
C      integer LUN_R470b_REF2
C      integer LUN_R650_REF
C      integer LUN_R650_REF2
c      integer LUN_TERRAIN
c      integer LUN_TERRAIN_NEW      
c      integer LUN_R865
      integer LUN_RAYL412
      integer LUN_RAYL466
      integer LUN_RAYL647
            
C      parameter (LUN_R412_REF=412008)
C      parameter (LUN_R470_REF=412009)
C      parameter (LUN_R470b_REF=412010)
C      parameter (LUN_R470_REF2=412013)
C      parameter (LUN_R470b_REF2=412014)
C      parameter (LUN_R650_REF=412011)
C      parameter (LUN_R650_REF2=412012)
c      parameter (LUN_TERRAIN=412006)
c      parameter (LUN_TERRAIN_NEW=412007)
c      parameter (LUN_R865=412005)
      parameter (LUN_RAYL412=412300)
      parameter (LUN_RAYL466=412302)
      parameter (LUN_RAYL647=412303)

      integer*4 fv, reclen, fh, file_version
      integer*4 prtn, pgs_io_gen_openf, pgs_io_gen_closef, pgs_pc_getreference
      character pathname*256, msg*256, CheckStr*512, raylfile*256

c-----------------------------------------------------------------------
c      Load surface reflectivity data base
c-----------------------------------------------------------------------

      reclen = 0
      
c      fv=1
C      prtn=pgs_io_gen_openf(LUN_TERRAIN,PGSd_IO_Gen_RseqUnf,reclen,fh,fv)
C      if (prtn .ne. 0) then
C         msg = "Error opening LUN_TERRAIN lun."
C         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'get_newtables')
C      endif
C      read (fh) terrain_flag
C      prtn = pgs_io_gen_closef(fh)
C
C      fv=1
C      prtn=pgs_io_gen_openf(LUN_TERRAIN_NEW,PGSd_IO_Gen_RseqUnf,reclen,fh,fv)
C      if (prtn .ne. 0) then
C         msg = "Error opening LUN_TERRAIN_NEW lun."
C         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'get_newtables')
C      endif
C      read (fh) terrain_flag_new
C      prtn = pgs_io_gen_closef(fh)
      
C      fv=1
C      prtn=pgs_io_gen_openf(LUN_R412_REF,PGSd_IO_Gen_RseqUnf,reclen,fh,fv)
C      if (prtn .ne. 0) then
C         msg = "Error opening LUN_R412_REF lun."
C         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'get_newtables')
C      endif
C      read (fh) xsfc412_bk
C      prtn = pgs_io_gen_closef(fh)
      
C      fv=1
C      prtn=pgs_io_gen_openf(LUN_R470_REF,PGSd_IO_Gen_RseqUnf,reclen,fh,fv)
C      if (prtn .ne. 0) then
C         msg = "Error opening LUN_R470_REF lun."
C         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'get_newtables')
C      endif
C      read (fh) xsfc470_bk
C      prtn = pgs_io_gen_closef(fh)
      
C      fv=1
C      prtn=pgs_io_gen_openf(LUN_R470b_REF,PGSd_IO_Gen_RseqUnf,reclen,fh,fv)
C      if (prtn .ne. 0) then
C         msg = "Error opening LUN_R470b_REF lun."
C         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'get_newtables')
C      endif
C      read (fh) xsfc470b_bk
C      prtn = pgs_io_gen_closef(fh)

C      fv=1
C      prtn=pgs_io_gen_openf(LUN_R470_REF2,PGSd_IO_Gen_RseqUnf,reclen,fh,fv)
C      if (prtn .ne. 0) then
C         msg = "Error opening LUN_R470_REF2 lun."
C         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'get_newtables')
C      endif
C      read (fh) xsfc470_bk2
C      prtn = pgs_io_gen_closef(fh)
      
C      fv=1
C      prtn=pgs_io_gen_openf(LUN_R470b_REF2,PGSd_IO_Gen_RseqUnf,reclen,fh,fv)
C      if (prtn .ne. 0) then
C         msg = "Error opening LUN_R470b_REF2 lun."
C         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'get_newtables')
C      endif
C      read (fh) xsfc470b_bk2
C      prtn = pgs_io_gen_closef(fh)
      
C      fv=1
C      prtn=pgs_io_gen_openf(LUN_R650_REF,PGSd_IO_Gen_RseqUnf,reclen,fh,fv)
C      if (prtn .ne. 0) then
C         msg = "Error opening LUN_R650_REF lun."
C         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'get_newtables')
C      endif
C      read (fh) xsfc650_bk
C      prtn = pgs_io_gen_closef(fh)

C      fv=1
C      prtn=pgs_io_gen_openf(LUN_R650_REF2,PGSd_IO_Gen_RseqUnf,reclen,fh,fv)
C      if (prtn .ne. 0) then
C         msg = "Error opening LUN_R650_REF2 lun."
C         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'get_newtables')
C      endif
C      read (fh) xsfc650_bk2
C      prtn = pgs_io_gen_closef(fh)

c-----------------------------------------------------------------------
c      Also Rayleigh file for polarization correction
c-----------------------------------------------------------------------

C      fv = 1
C      prtn=pgs_io_gen_openf(LUN_RAYL,PGSd_IO_Gen_RseqUnf,reclen,fh,fv)
C      if (prtn .ne. 0) then
C         msg = "Error opening LUN_RAYL lun."
C         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'get_newtables')
C      endif
C      msg = "opened LUN_RAYL lun."
C      call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC, msg, 'get_newtables')
C      read (fh, end=900, err=900) nvalc
C      prtn = pgs_io_gen_closef(fh)

      file_version = 1
      prtn=pgs_pc_getreference(LUN_RAYL412,file_version,raylfile)
      if (prtn .lt. 0) then
         msg = "Error opening LUN_RAYL412 lun."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      endif
      msg = "Will use RAYL file "//raylfile
      call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC, msg, 'DeepBlue')

      open (1,file=raylfile, status='old') ! ori.
      read (1,*) nvalc
      close (1)

      file_version = 1
      prtn=pgs_pc_getreference(LUN_RAYL466,file_version,raylfile)
      if (prtn .lt. 0) then
         msg = "Error opening LUN_RAYL466 lun."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      endif
      msg = "Will use RAYL file "//raylfile
      call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC, msg, 'DeepBlue')

      open (1,file=raylfile, status='old') ! ori.
      read (1,*) nvalc2
      close (1)

      file_version = 1
      prtn=pgs_pc_getreference(LUN_RAYL647,file_version,raylfile)
      if (prtn .lt. 0) then
         msg = "Error opening LUN_RAYL647 lun."
         call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'DeepBlue')
      endif
      msg = "Will use RAYL file "//raylfile
      call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC, msg, 'DeepBlue')

      open (1,file=raylfile, status='old') ! ori.
      read (1,*) nvalc3
      close (1)
      
      return

900   continue
      msg = "Error reading file."
      call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC, msg, 'get_newtables')
      end
      
