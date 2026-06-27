             subroutine file_close(choice_flag,error_flag,
     &     modfil_L1B_1km,modfil_L1B_500,modfil_L1B_250,
     &     MODFIL_Geo,MODFIL_Cldmsk,modfil_MOD04,MODFIL_MOD05,
     &     MODFIL_MOD07,MODFIL_ANC,
     &     handle_LUT466,handle_LUT553,handle_LUT644,handle_LUT213,
     &     handle_LUTMAP,
     &     HANDLE_S,HANDLE_L,HANDLE_INSCI,HANDLE_OUTSCI,handle_QC,
     &     groups,HDFattNms,NumHandles)


      implicit none
       SAVE
      include 'mapi.inc'
      include 'PGS_MODIS_39500.f'
      include 'mod04.inc'

C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION: Close files according to value of 'choice_flag'.
C
C !INPUT PARAMETERS:
C
C    character  choice_flag  
C                        The flag used to indicate which algorithm to run.
C                        'MA': Run MODIS MOD04 Land & Sea
C                        'ML': Run MODIS MOD04 Land
C                        'MS': Run MODIS MOD04 Sea
C                        'NL': Run non-MODIS Land
C                        'NS': Run non-MODIS Sea
C
C    integer   modfil_*  File handle structure for HDF files
C    integer   handle_*  File handle used to manipulate files opened
C                        for read or write.
C    character  groups   A character array with size of
C                        [pGSd_MET_NUM_OF_GROUPS][PGS_MET_GROUP_NAME_L],
C                        where PGSd_MET_NUM_OF_GROUPS is 20 and
C                        PGSd_MET_GROUP_NAME_L is 50.  Each row in the 
C                        array stores a handles to an internal ODL tree 
C                        structure which will be written out a an 
C                        ECS PVL attribute.
C
C !OUTPUT PARAMETERS:
C
C    logical    error_flag The flag used to indicate code status.
C
C !REVISION HISTORY:
C
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

      character*(*) choice_flag
      CHARACTER*(PGSd_MET_GROUP_NAME_L) groups(PGSd_MET_NUM_OF_GROUPS)
      CHARACTER*(MAX_ECS_NAME_L-1) HDFAttNms(PGSd_MET_NUM_OF_GROUPS)
      character*1024 msg

      integer rtn,INVENTORYMETADATA,
     &        modfil_L1B_1km(*),modfil_L1B_500(*),modfil_L1B_250(*),
     &        modfil_Geo(*),modfil_CldMsk(*),modfil_mod04(*),
     &        modfil_mod05(*),modfil_mod07(*),modfil_Anc(*),
     &        NumHandles
      logical error_flag

C      parameter(INVENTORYMETADATA = 2)

c ... initialization

C      HDFAttNms(INVENTORYMETADATA) = MECS_CORE

c rhucek 02/11/98:  for testing only
       HANDLE_OUTSCI = 5
        
c ... input check

      if (choice_flag.ne.'MA' .and. choice_flag.ne.'ML' .and.
     &    choice_flag.ne.'MS' .and. choice_flag.ne.'NS' .and.
     &    choice_flag.ne.'NL') then
          error_flag = .true.
          call MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     &    'Invalid input for choice_flag','mod04_file_close_V2.f')
          return
      end if

c--------------------------------------------------------------
c ... Close MODIS MOD_PR04 Land Lookup table files:
c--------------------------------------------------------------

      call Modis_IO_Gen_Closef(LRN_QC,handle_QC)

      if (choice_flag.eq.'MA' .or. choice_flag.eq.'ML') then
          call Modis_IO_Gen_Closef(LRN_LUT466,handle_LUT466)
          call Modis_IO_Gen_Closef(LRN_LUT553,handle_LUT553)
          call Modis_IO_Gen_Closef(LRN_LUT644,handle_LUT644)
          call Modis_IO_Gen_Closef(LRN_LUT213,handle_LUT213)
          call Modis_IO_Gen_Closef(LRN_LUTMAP,handle_LUTMAP)
      end if

c--------------------------------------------------------------
c ... Close MODIS MOD_PR04 Sea files:
c ... ','smallcase1-3', and 'largecase1-3'.
c--------------------------------------------------------------

      if (choice_flag.eq.'MA' .or. choice_flag.eq.'MS') then
          call Modis_IO_Gen_Closef(S1,handle_S(1))
          call Modis_IO_Gen_Closef(S2,handle_S(2))
          call Modis_IO_Gen_Closef(S3,handle_S(3))
          call Modis_IO_Gen_Closef(S4,handle_S(4))
          call Modis_IO_Gen_Closef(L1,handle_L(1))
          call Modis_IO_Gen_Closef(L2,handle_L(2))
          call Modis_IO_Gen_Closef(L3,handle_L(3))
          call Modis_IO_Gen_Closef(L4,handle_L(4))
      end if

c--------------------------------------------------------------
c ... Close common MODIS files:
c ... L1B, Geolocation, Cloud Mask, MOD_PRANC,
c ... MOD_PR05, MOD_PR07, and MOD_PR04 files.
c--------------------------------------------------------------

      if (choice_flag(1:1).eq.'M') then

c ... Close L1B - 1km resolution file

       rtn = clmfil(modfil_L1B_1km)
       if (rtn.eq.mapiok) then
           call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,
     &     'close L1B 1km file','mod04_file_close_V2.f')
       else
         error_flag = .true.
         msg=
     1   'Failed to close L1B 1km file'
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msg,
     &   'mod04_file_close_V2.f')
       end if
    
c ... Close L1B - 500m resolution file

       rtn = clmfil(modfil_L1B_500)
       if (rtn.eq.mapiok) then
           call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,
     &     'close L1B 500m file','mod04_file_close_V2.f')
       else
         error_flag = .true.
         msg=
     1   'Failed to close L1B 500m file'
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'

c rhucek 10/07/99:  added "msg" buffer to MODIS_SMF_SETDYNAMICMSG to argument list
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msg,'mod04_file_close_V2.f')
       end if
    
c ... Close L1B - 250m resolution file

       rtn = clmfil(modfil_L1B_250)
       if (rtn.eq.mapiok) then
           call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,
     &     'close L1B 250m file','mod04_file_close_V2.f')
       else
         error_flag = .true.
         msg=
     1   'Failed to close L1B 250m file'
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msg,
     &   'mod04_file_close_V2.f')
       end if
    
c ... Close Geolocation file

       rtn = clmfil(modfil_Geo)
       if (rtn.eq.mapiok) then
           call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,
     &     'close Geolocation file','mod04_file_close_V2.f')
       else
         error_flag = .true.
         msg=
     1   'Failed to close geolocation file'
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msg,
     &   'mod04_file_close_V2.f')
       end if
    
c ... Close Cloud Mask file

       rtn = clmfil(modfil_CldMsk)
       if (rtn.eq.mapiok) then
           call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,
     &     'close cloud mask file','mod04_file_close_V2.f')
       else
         error_flag = .true.
         msg=
     1   'Failed to close geolocation file'
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msg,
     &   'mod04_file_close_V2.f')
       end if
    
c ... Close ancillary file

C          rtn = clmfil(modfil_ANC)
C          if (rtn.eq.mapiok) then
C              call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,
C     &             'close ancillary file','file_close')
C          else
C              error_flag = .true.
C              call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,
C     &             'Failed to close ancillary file','file_close')
C          end if
    
c ... Close MOD_PR05 file

       rtn = clmfil(modfil_mod05)
       if (rtn.eq.mapiok) then
         call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,
     &   'close mod05 product file','mod04_file_close_V2.f')
       else
         error_flag = .true.
         msg=
     1   'Failed to close MOD_PR05 file'
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msg,
     &   'mod04_file_close_V2.f')
       end if
    
c ... Close MOD_PR07 file

       rtn = clmfil(modfil_mod07)
       if (rtn.eq.mapiok) then
         call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,
     &   'close mod07 product file','mod04_file_close_V2.f')
       else
         error_flag = .true.
         msg=
     1   'Failed to close MOD_PR05 file'
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msg,
     &   'mod04_file_close_V2.f')
       end if
    
    
C ... Write HDF attributes & close MOD_PR04 product file

*/  cpmfil -> clmfil

C*     rtn = cpmfil(modfil_MOD04,groups,HDFAttNms,1)

       rtn = cpmfil(modfil_MOD04,groups,HDFAttNms,NumHandles)
       if (rtn.eq.mapiok) then
         call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,
     &   'close mod04 product file','mod04_file_close_V2.f')
       else
         error_flag = .true.
         msg=
     1   'Failed to close MOD_PR04 file while writing metadata'
     2   // char(10) // 'Operator Action:  Check system resources/environment, '
     3   // char(10) // 'PCF, and SDPTK configuration.  If a fault is identified, '
     4   // char(10) // 'correct and rerun PGE.  Otherwise, notify SDST.'
         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msg,
     &   'mod04_file_close_V2.f')
       end if

      end if

c--------------------------------------------------------------
c ... Close non-MODIS files:
c ... MOD_PR04 land :'modaerland.input' and 'modaerland.output'
c ... MOD_PR04 sea :'modaeroceanv2.input' and 'modaeroceanv2.output'
c--------------------------------------------------------------

      if (choice_flag(1:1).eq.'N') then
          call Modis_IO_Gen_Closef(INSCI,handle_insci)
C          call Modis_IO_Gen_Closef(OUTSCI,handle_outsci)
      end if

      return
      end
