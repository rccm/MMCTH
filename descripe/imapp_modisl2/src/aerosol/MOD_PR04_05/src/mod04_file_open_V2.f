         subroutine file_open(choice_flag,error_flag,
     &           modfil_L1B_1km,modfil_L1B_500,modfil_L1B_250,
     &           FN_L1B_1km,FN_L1B_500,FN_L1B_250,
     &           modfil_Geo,modfil_CldMsk,modfil_mod04,modfil_mod05,
     &           modfil_mod07,modfil_Anc,
     &           handle_LUT466,handle_LUT553,
     &           handle_LUT644,handle_LUT213,handle_LUTMAP,
     &           handle_S,handle_L,
     &           handle_INSCI,handle_OUTSCI,handle_QC)
      implicit none
      include 'mapi.inc'
      include 'PGS_SMF.f'
      include 'PGS_PC.f'
      include 'PGS_IO.f'
      include 'PGS_MODIS_39500.f'
      include 'mod04.inc'
      SAVE
C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION: Open files according to value of choice_flag.
C
C !INPUT PARAMETERS:   
C               
C    character  choice_flag  The flag used to indicate which algorithm 
C                            to run.
C                            'MA': Run MODIS MOD04 land & Sea
C                            'ML': Run MODIS MOD04 land
C                            'MS': Run MODIS MOD04 Sea
C                            'NL': Run non-MODIS Land
C                            'NS': Run non-MODIS Sea
C
C !OUTPUT PARAMETERS: 
C
C    logical    error_flag The flag used to indicate code status.
C    character  FN_L1B     L1B file name
C    integer    modfil_*   File handle structure for HDF files
C    integer    handle_*   File handle used to manipulate files opened
C                          for read or write.
C
C !REVISION HISTORY:  
C $Log$
C !TEAM-UNIQUE HEADER: 
C
C    This software was developed by the MODIS Science Data Support Team
C    (SDST) for the National Aeronautics and Space Administration, 
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C !REFERENCES and CREDITS:
C
C    Written by:
C    Vicky Lin                       03/03/97
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd, Seabrook MD 20706
C
C    vlin@ltpmail.gsfc.nasa.gov
C
C !END
C----------------------------------------------------------------------

      character*(*) FUNCNAME
      parameter    (FUNCNAME = 'file_open')
 
      character*(PGSd_PC_FILE_PATH_MAX) HDF_FILENAME
*/ start : modified by Allen Chu
      character*(*) FN_L1B_1km,FN_L1B_500,FN_L1B_250
*/ end
      character*(*)  choice_flag
      character*1024 msg
      integer i,rtn,file_access,rec_len,file_version,
     &        modfil_L1B_1km(*),modfil_L1B_500(*),modfil_L1B_250(*),
     &        modfil_Geo(*),modfil_CldMsk(*),modfil_mod04(*),
     &        modfil_mod05(*),modfil_mod07(*),modfil_Anc(*)
      integer pgs_io_gen_openf,pgs_pc_getreference
      logical error_flag

      integer num_args
      integer FlagRA
      character FlagBuff*10
      integer iargc

      num_args = iargc ( )
      
      if(num_args .eq. 1) then
         call getarg ( 1, FlagBuff )
         read (FlagBuff,*) FlagRA
      else
C      This is the default value
         FlagRA = 0
      endif

c ... initialization

      FN_L1B_1km = ' '
      FN_L1B_500 = ' '
      FN_L1B_250 = ' '
      error_flag = .false.
      handle_LUT466 = -1
      handle_LUT553 = -1
      handle_LUT644 = -1
      handle_LUT213 = -1
      handle_LUTMAP = -1
      handle_INSCI = -1
      handle_OUTSCI = -1
      handle_QC = -1
      file_version = 1

      do i = 1, lut_indx
         handle_S(i) = -1
         handle_L(i) = -1
      end do

      do i = 1, MODFILLEN
         modfil_L1B_1km(i) = 0
         modfil_L1B_500(i) = 0
         modfil_L1B_250(i) = 0
         modfil_Geo(i) = 0
         modfil_CldMsk(i) = 0
         modfil_mod04(i) = 0
         modfil_mod05(i) = 0
         modfil_mod07(i) = 0
         modfil_Anc(i) = 0
      end do
      
c ... input check

      if (choice_flag.ne.'MA' .and. choice_flag.ne.'ML' .and.
     &    choice_flag.ne.'MS' .and. choice_flag.ne.'NS' .and.
     &    choice_flag.ne.'NL') then
          error_flag = .true.
          call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     &         'Invalid input for choice_flag','file_open')
          return
      end if

c------------------------------------------------------------
c ... Open MODIS MOD_PR04 Land Lookup table files:
c------------------------------------------------------------

      if (choice_flag.eq.'MA' .or. choice_flag.eq.'ML') then
          file_access = PGSd_IO_Gen_RSeqFrm

          call MODIS_IO_GEN_OPENF(LRN_LUT466,file_access,rec_len,
     &         handle_LUT466,file_version)
          call MODIS_IO_GEN_OPENF(LRN_LUT553,file_access,rec_len,
     &         handle_LUT553,file_version)
          call MODIS_IO_GEN_OPENF(LRN_LUT644,file_access,rec_len,
     &         handle_LUT644,file_version)
          call MODIS_IO_GEN_OPENF(LRN_LUT213,file_access,rec_len,
     &         handle_LUT213,file_version)
          call MODIS_IO_GEN_OPENF(LRN_LUTMAP,file_access,rec_len,
     &         handle_LUTMAP,file_version)

      end if

c rhucek 02/09/98:  Supplemented status messaging when opening QC file.

c.....try to open QC file for update
      file_version = 1
      file_access = PGSd_IO_Gen_USeqFrm ! Assume file exist

      rtn = pgs_io_gen_openf(LRN_QC,file_access,rec_len,handle_QC,file_version)

      If (rtn .NE. PGS_S_SUCCESS) Then
         msg =  'PGS_IO_GEN_OpenF unable to open MOD04_QC file in update mode, '
     1   // char(10) // 'i.e. file does not exit.  Next try opening file in "create" mode.' 

         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,MSG,FUNCNAME)

c........Try to open file for create (PGSd_IO_Gen_WSeqFrm)
         FILE_ACCESS = PGSd_IO_Gen_WSeqFrm 

         rtn = pgs_io_gen_openf(LRN_QC,FILE_ACCESS,rec_len,handle_QC,FILE_VERSION)

         If (rtn .NE. PGS_S_SUCCESS) Then
            msg = 'PGS_IO_GEN_OpenF able to open MOD05_QC in "create" mode. '
     1      // char(10) // 'Continue processing without access to MOD04_QC.'

            CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,MSG,FUNCNAME)
         Else
            msg = 'PGS_IO_GEN_OpenF successfully opened MOD04_QC in "create" mode. '
            CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,MSG,FUNCNAME)
         End If

      End If


c------------------------------------------------------------
c ... Open MODIS MOD_PR04 Sea files:
c ... 'modaeroceanv0.ext','smallcase1-6', and 'largecase1-6'.
c------------------------------------------------------------

      if (choice_flag.eq.'MA' .or. choice_flag.eq.'MS') then
          
          file_access = PGSd_IO_Gen_RSeqFrm
          call MODIS_IO_GEN_OPENF(S1,file_access,rec_len,
     &         handle_S(1),file_version)
          call MODIS_IO_GEN_OPENF(S2,file_access,rec_len,
     &         handle_S(2),file_version)
          call MODIS_IO_GEN_OPENF(S3,file_access,rec_len,
     &         handle_S(3),file_version)
           call MODIS_IO_GEN_OPENF(S4,file_access,rec_len,
     &         handle_S(4),file_version)
          call MODIS_IO_GEN_OPENF(L1,file_access,rec_len,
     &         handle_L(1),file_version)
          call MODIS_IO_GEN_OPENF(L2,file_access,rec_len,
     &         handle_L(2),file_version)
          call MODIS_IO_GEN_OPENF(L3,file_access,rec_len,
     &         handle_L(3),file_version)
             call MODIS_IO_GEN_OPENF(L4,file_access,rec_len,
     &         handle_L(4),file_version)  
          call MODIS_IO_GEN_OPENF(INSCI,file_access,rec_len,
     &                            handle_INSCI,file_version)
          

              end if

c------------------------------------------------------------
c ... Open common MODIS files: 
c ... L1B, Geolocation, Cloud Mask, MOD_PRANC,
c ... MOD_PR05, MOD_PR07, and MOD_PR04 files.
c------------------------------------------------------------

      if (choice_flag(1:1).eq.'M') then

c ... Open L1B - 1km resolution file

      file_version = 1
      rtn = -1
      if( FlagRA .eq. 1) then
         rtn = pgs_pc_getreference(LRN_L1B_1km_RA,file_version,FN_L1B_1km)
      else
         rtn = pgs_pc_getreference(LRN_L1B_1km,file_version,FN_L1B_1km)
      endif
      if (rtn.ne.PGS_S_SUCCESS) then
          error_flag = .true.
         msg=
     1   'pgs_pc_getreference for L1B 1km file'
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msg,
     &        'mod04_file_open_V2.f')
      else
          rtn = -1 
          rtn = opmfil(FN_L1B_1km,'r',modfil_L1B_1km)
          if (rtn.ne.mapiok) then
            error_flag = .true.
            msg=
     1      'Failed to open L1B 1km file'
     2      // char(10) // 'Operator Action:  Check system resources/environment, '
     3      // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4      // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
            CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msg,
     &           'mod04_file_open_V2.f')
          end if
      end if
      if (error_flag) return

c ... Open L1B - 500m resolution file

      file_version = 1
      rtn = -1
      rtn = pgs_pc_getreference(LRN_L1B_500,file_version,FN_L1B_500)
      if (rtn.ne.PGS_S_SUCCESS) then
          error_flag = .true.
          msg=
     1    'pgs_pc_getreference for L1B 500m file'
     2    // char(10) // 'Operator Action:  Check system resources/environment, '
     3    // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4    // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msg,
     &         'mod04_file_open_V2.f')
      else
          rtn = -1 
          rtn = opmfil(FN_L1B_500,'r',modfil_L1B_500)
          if (rtn.ne.mapiok) then
            error_flag = .true.
            msg=
     1      'Failed to open L1B 500m file'
     2      // char(10) // 'Operator Action:  Check system resources/environment, '
     3      // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4      // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
            CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msg,
     &           'mod04_file_open_V2.f')
          end if
      end if
      if (error_flag) return

c ... Open L1B - 250m resolution file

      file_version = 1
      rtn = -1
      rtn = pgs_pc_getreference(LRN_L1B_250,file_version,FN_L1B_250)
      if (rtn.ne.PGS_S_SUCCESS) then
          error_flag = .true.
          msg=
     1    'pgs_pc_getreference for L1B 250m file'
     2    // char(10) // 'Operator Action:  Check system resources/environment, '
     3    // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4    // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
          CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msg,
     &         'mod04_file_open_V2.f')
      else
          rtn = -1 
          rtn = opmfil(FN_L1B_250,'r',modfil_L1B_250)
          if (rtn.ne.mapiok) then
            error_flag = .true.
            msg=
     1      'Failed to open L1B 250m file'
     2      // char(10) // 'Operator Action:  Check system resources/environment, '
     3      // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4      // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
            CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msg,
     &           'mod04_file_open_V2.f')
          end if
      end if
      if (error_flag) return

c ... Open Geolocation file

         file_version = 1
         HDF_FILENAME = ' '
         rtn = -1
         rtn = pgs_pc_getreference(LRN_Geo,file_version,HDF_FILENAME)
         if (rtn.ne.PGS_S_SUCCESS) then
           error_flag = .true.
           msg=
     1     'pgs_pc_getreference for Geo file'
     2     // char(10) // 'Operator Action:  Check system resources/environment, '
     3     // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4     // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
             CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msg,
     &       'mod04_file_open_V2.f')
         else
           rtn = -1
           rtn = opmfil(HDF_FILENAME,'r',modfil_Geo)
           if (rtn.ne.mapiok) then
            error_flag = .true.
            msg=
     1      'Failed to open Geolocation file'
     2      // char(10) // 'Operator Action:  Check system resources/environment, '
     3      // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4      // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
            CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msg,
     &      'mod04_file_open_V2.f')
             end if
         end if
   
c ... Open Cloud Mask file
   
         file_version = 1
         HDF_FILENAME = ' '
         rtn = -1
         rtn = pgs_pc_getreference(LRN_CldMsk,file_version,HDF_FILENAME)
         if (rtn.ne.PGS_S_SUCCESS) then
           error_flag = .true.
           msg=
     1     'pgs_pc_getreference for cloud mask file'
     2     // char(10) // 'Operator Action:  Check system resources/environment, '
     3     // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4     // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
           CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msg,
     &          'mod04_file_open_V2.f')
         else
           rtn = -1
           rtn = opmfil(HDF_FILENAME,'r',modfil_CldMsk)
           if (rtn.ne.mapiok) then
             error_flag = .true.
             msg=
     1       'Failed to open cloud mask file'
     2       // char(10) // 'Operator Action:  Check system resources/environment, '
     3       // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4       // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
             CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msg,
     &       'mod04_file_open_V2.f')
           end if
         end if

c ... Open MOD_PR04 file

         file_version = 1
         HDF_FILENAME = ' '
         rtn = -1
         rtn = pgs_pc_getreference(LRN_MOD04,FILE_VERSION,HDF_FILENAME)
         if (rtn.ne.PGS_S_SUCCESS) then
           error_flag = .true.
           msg=
     1     'pgs_pc_getreference for MOD_PR04 file'
     2     // char(10) // 'Operator Action:  Check system resources/environment, '
     3     // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4     // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
           CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msg,
     &     'mod04_file_open_V2.f')
         else
            rtn = -1 
            rtn = OPMFIL(HDF_FILENAME,"a",MODFIL_MOD04)
            IF (RTN.EQ.mapiok) THEN
                CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,
     &            'Successfully opened "mod04.hdf" file',
     &            'file_open')
            ELSE
             error_flag = .true.
             msg=
     1       'Failed to open MOD_PR04 file'
     2       // char(10) // 'Operator Action:  Check system resources/environment, '
     3       // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4       // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
             CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msg,
     &       'mod04_file_open_V2.f')
            END IF
         end if
       
C   
C ... Open MOD_PR05 file
C
         file_version = 1
         HDF_FILENAME = ' '
         rtn = -1
         rtn = pgs_pc_getreference(LRN_MOD05,file_version,HDF_FILENAME)
         if (rtn.ne.PGS_S_SUCCESS) then
           error_flag = .true.
           msg=
     1     'pgs_pc_getreference for MOD_PR05 file'
     2     // char(10) // 'Operator Action:  Check system resources/environment, '
     3     // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4     // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
           CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msg,
     &     'mod04_file_open_V2.f')
           
         else
            rtn = -1 
            rtn = OPMFIL(HDF_FILENAME,"r",MODFIL_MOD05)
            IF (RTN.EQ.mapiok) THEN
                CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,
     &            'Successfully opened "mod05.hdf" file',
     &            'file_open')
            ELSE
             error_flag = .true.
             msg=
     1       'Failed to open MOD_PR05 file'
     2       // char(10) // 'Operator Action:  Check system resources/environment, '
     3       // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4       // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
             CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msg,
     &       'mod04_file_open_V2.f')
            END IF
         end if
C
C ... Open MOD_PR07 file
C   
         file_version = 1
         HDF_FILENAME = ' '
         rtn = -1
         rtn = pgs_pc_getreference(LRN_MOD07,file_version,HDF_FILENAME)
         if (rtn.ne.PGS_S_SUCCESS) then
           error_flag = .true.
           msg=
     1     'pgs_pc_getreference for MOD_PR05 file'
     2     // char(10) // 'Operator Action:  Check system resources/environment, '
     3     // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4     // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
           CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msg,
     &     'mod04_file_open_V2.f')
        else
            rtn = -1 
            rtn = OPMFIL(HDF_FILENAME,"r",MODFIL_MOD07)
            IF (RTN.EQ.mapiok) THEN
                CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,
     &            'Successfully opened "mod07.hdf" file',
     &            'file_open')
            ELSE
             error_flag = .true.
             msg=
     1       'Failed to open MOD_PR07 file'
     2       // char(10) // 'Operator Action:  Check system resources/environment, '
     3       // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4       // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
             CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msg,
     &       'mod04_file_open_V2.f')
            END IF
         end if
C   
C ... Open WISCONSIN Ancillary file
C   
C         file_version = 1
C         HDF_FILENAME = ' '
C         rtn = -1
C         rtn = pgs_pc_getreference(LRN_WISC_ANC,file_version,
C     &                             HDF_FILENAME)
C         if (rtn.ne.PGS_S_SUCCESS) then
C             error_flag = .true.
C             call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,
C     &         'pgs_pc_getreference for "wis_anci"','file_open')
C         else
C             rtn = -1
C             rtn = opmfil(HDF_FILENAME,'r',modfil_ANC)
C              IF (RTN.EQ.mapiok) THEN
C                CALL MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,
C     &            'Successfully opened "wis_anci" file',
C     &            'file_open')
C            ELSE
C                error_flag = .true.
C                CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,
C     &       'OPMFIL for "wis_anci"','file_open')
C            END IF
C         end if
      end if

c------------------------------------------------------------
c ... Open non-MODIS files:
c ... MOD_PR04 land :'modaerland.input' and 'modaerland.output'
c ... MOD_PR04 sea :'modaeroceanv2.input' and 'modaeroceanv2.output'
c------------------------------------------------------------

 
C          file_access = PGSd_IO_Gen_WSeqFrm 
C                                 ! Assume file does not exist
C          rtn = pgs_io_gen_openf(OUTSCI,file_access,rec_len,
C     &                           handle_OUTSCI,file_version)
C          if (rtn .ne. PGS_S_SUCCESS) then
C              file_access = PGSd_IO_Gen_USeqFrm 
C                                 ! open file for update
C              call MODIS_IO_GEN_OPENF(OUTSCI,file_access,rec_len,
C     &                                handle_OUTSCI,file_version)
C          end if

      

      return
      end
