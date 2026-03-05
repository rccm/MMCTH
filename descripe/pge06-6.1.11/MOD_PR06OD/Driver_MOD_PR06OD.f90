      program Driver_MOD_PR06OD

!-------------------------------------------------------------------
!
!!F77
!
!!Description:
! Driver program for MODIS cloud property algorithm.
!
!!Input Parameters:
! None
!
!!Output Parameters:
! None
!
!!Revision History:
! Revision 3.0  2004/08/24 mag
! the new driver
!
!!Team-unique Header:
! Cloud Retrieval Group, NASA/GSFC, Code 913, Greenbelt, Maryland, USA
!
!!Credit and Reference:
! Programmed by:
! Mark Gray (mag)
! L-3 GSI
! NASA/GSFC, code 913
! Greenbelt, Maryland, USA
!
!!End
!---------------------------------------------------------------------
	use names
	use mod06_run_settings

      implicit none
 
#ifdef USE_PGS
      include 'PGS_SMF.f'
      include 'PGS_PC.f'
      include 'PGS_MODIS_39500.f'
      include 'mapi.f90.inc'
#endif

#ifdef USE_PGS
!.....ECS inventory metadata variables
	real :: statistics(7)
	integer :: ErrorFlag

!--------------------------------------------------------------------------------
! Begin executable code by writing "Begin" message to LogStatus file 
!--------------------------------------------------------------------------------

         integer NUM_Input_LRN, NUM_PSA, NUM_MeasParam, NUM_ArchivePSA_SC, rtn, Vrsn_No, LUN

         character*(*) FUNCNAME
         PARAMETER    (FUNCNAME = 'Driver_MOD_PR06OD')

         integer       MAX_NUM_INPUT_PNTRS
         parameter    (MAX_NUM_INPUT_PNTRS = 50)

         integer       MAX_NUM_PSA
         parameter    (MAX_NUM_PSA = 20)

         integer       MAX_NUM_ArchivePSA_SC
         parameter    (MAX_NUM_ArchivePSA_SC = 8)

         integer       FAIL,    SUCCEED
         parameter    (FAIL=-1, SUCCEED=0)

         integer       MAX_NUM_MEAS_PARM
         parameter    (MAX_NUM_MEAS_PARM = 10)

         integer       INVENTORYMETADATA,     ARCHIVEMETADATA
         parameter    (INVENTORYMETADATA = 2, ARCHIVEMETADATA=3)

         character*(MAX_ECS_NAME_L-1)      HDFAttNms(PGSd_MET_NUM_OF_GROUPS)


         character*(PGSd_MET_GROUP_NAME_L) Met_Handles(PGSd_MET_NUM_OF_GROUPS)
         
         real    PSAValue(MAX_NUM_PSA) 
         

         integer Input_LRN     (MAX_NUM_INPUT_PNTRS)  /MAX_NUM_INPUT_PNTRS*-99999/, &
              Input_Vrsn    (MAX_NUM_INPUT_PNTRS)  /MAX_NUM_INPUT_PNTRS*1/, &
              OUTPUT_LRN, DATETIME_LRN
         integer QA_Miss(MAX_NUM_MEAS_PARM)     / MAX_NUM_MEAS_PARM*0 /
         
         character*100 Name_ArchivePSA_SC(MAX_NUM_ArchivePSA_SC)
         real Value_ArchivePSA_SC(MAX_NUM_ArchivePSA_SC)

         character*6 LUNSTR


         !.....Function declarations
         integer  Set_ArchMet_MOD06, pgs_pc_getreference, pgs_pc_getconfigdata, set_invmet_mod06, pgs_td_asciitime_atob

         logical error_flag
         logical ExtGeoPntr_Flag 
         character*2048                    msgbuf
         character*60   Name_MeasParm(MAX_NUM_MEAS_PARM) / MAX_NUM_MEAS_PARM*' ' /
         character*30   Auto_QA_Flag(MAX_NUM_MEAS_PARM) / MAX_NUM_MEAS_PARM*' ' /
         character*30   Auto_QA_Expl(MAX_NUM_MEAS_PARM) / MAX_NUM_MEAS_PARM*' ' /
         character*30   PSAName(MAX_NUM_PSA) / MAX_NUM_PSA*' ' /
         integer modfil_out(MODFILLEN)
         integer NUM_of_HDFAttrNms, i
         character(len=500) datetime, jdatetime, temp, tempDateTime1, tempDateTime2, tempString
         integer iDateTime1, iDateTime2
         integer datePos
         
        integer  num_args
        integer  FlagRA
        character FlagBuff*10
        integer  iargc

        logical lErrorFlag
        integer  Set_PSA
        character*4 pcf_satid
        character*32  doi_char
        character*128 proctype
        integer LUN_Rep_Actual
        parameter (LUN_Rep_Actual = 800504)
        integer LUN_Sat_Instrument
        parameter (LUN_Sat_Instrument = 800510)

        integer DFACC_WRITE, DFNT_CHAR8
        parameter   (DFACC_WRITE = 2, &
             DFNT_CHAR8  = 4)

        integer sfstart, sfscatt
        integer sfend
        integer doiLen

#else
!.....Parameter Declarations


!*****************************************************
!     Variable Declarations
!*****************************************************


!.....ECS inventory metadata variables

	  character(len = 200) :: par_file, dummy
	character(len=300) :: dirname, dayname, yearname, timename, outdirname
	character(len=500) :: temp
	real :: statistics(7)
	integer :: ErrorFlag
	integer :: mylun, i, idx
#endif




#ifdef USE_PGS

!........add MODIS products to list of inputs
#ifdef USE_TOAST
         NUM_Input_LRN = 31
#else
         NUM_Input_LRN = 30 
#endif


         num_args = iargc ( )
      
         if(num_args .eq. 1) then
            call getarg ( 1, FlagBuff )
            read (FlagBuff,*) FlagRA
         else
            !     This is the default value
            FlagRA = 0
         endif
      
         if( FlagRA .eq. 1) then
            Input_LRN(1)  = 430001
         else  
            Input_LRN(1)  = 700002
         endif
         Input_LRN(2)  = 422500
         Input_LRN(3)  = 600000

! rhucek 08/26/02: revise list of library files  
         !........add ancillary products to list of inputs
         Input_LRN(4) = 900000
#ifdef USE_TOAST
         Input_LRN(5) = 900020
#endif
         Input_LRN(6) = 900040
         Input_LRN(7) = 900100

         Input_LRN(8)  = 417502
         Input_LRN(9)  = 417501
         Input_LRN(10) = 417504
         Input_LRN(11) = 417505
         Input_LRN(12) = 417506
         Input_LRN(13) = 417507
         Input_LRN(14) = 417508
         Input_LRN(15) = 417509
         Input_LRN(16) = 417510
         Input_LRN(17) = 417511
         Input_LRN(18) = 417512
         Input_LRN(19) = 417513
         Input_LRN(20) = 417514
         Input_LRN(21) = 417515
         Input_LRN(22) = 417503
         Input_LRN(23) = 417521
         Input_LRN(24) = 417540
         Input_LRN(25) = 417527
         Input_LRN(26) = 417522
         Input_LRN(27) = 417523
         Input_LRN(28) = 417524
         Input_LRN(29) = 417525
         Input_LRN(30) = 417526
         Input_LRN(31) = 417516

         OUTPUT_LRN = 412500
         
         DATETIME_LRN = 10258

         NUM_PSA     = 0 
         NUM_MeasParam = 0
         
!........No Archive metadata set by MOD_PR06OD 
         NUM_ArchivePSA_SC       = 0
         temp = ''
         LUN = DATETIME_LRN
         write(LUNSTR,'(I6)') LUN
         rtn = pgs_pc_getconfigdata(LUN, datetime)
         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.
            
            msgbuf = &
                 'pgs_pc_getconfigdata '// trim(LUNSTR) // ' from PCF. ' &
                 // char(10) // 'Cloud Optical Properties Retrieval code to exit.' &
                 // char(10) // 'Operator Action:  Notify SDST. '
            
            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         end if
         
         rtn = PGS_TD_ASCIItime_AtoB(datetime,jdatetime)
         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.
            
            msgbuf = &
                 'problem in PGS_TD_UTCtoUTCjd function '
            
            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         end if
         
         temp = ''
         temp = jdatetime(1:4)
         read(temp,*) MYYEAR
         
         temp = ''
         temp = jdatetime(6:8)
         read(temp, *) MYDAY
         
         temp = ''
         temp = jdatetime(10:11) // jdatetime(13:14)
         read(temp, *) MYTIME
         
         
         
         
         temp = ''
         Vrsn_No = 1
         LUN = Input_LRN(1)
         write(LUNSTR,'(I6)') LUN
         !         ........retrieve name of hdf product from pcf
         rtn = pgs_pc_getreference(LUN,Vrsn_No,temp)
         
         
         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.
            
            msgbuf = &
                 'pgs_pc_getreference detected error retrieving file with LUN# '// trim(LUNSTR) // ' from PCF. ' &
                 // char(10) // 'Cloud Optical Properties Retrieval code to exit.' &
                 // char(10) // 'Operator Action:  Notify SDST. '
            
            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         end if
         Alevel1b_name(1) = temp
         
         temp = ''
         Vrsn_No = 1
         LUN = Input_LRN(2)
         write(LUNSTR,'(I6)') LUN
         !         ........retrieve name of hdf product from pcf
         rtn = pgs_pc_getreference(LUN,Vrsn_No,temp)
         
         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.
            
            msgbuf = &
                 'pgs_pc_getreference detected error retrieving file with LUN# '// trim(LUNSTR) // ' from PCF. ' &
                 // char(10) // 'Cloud Optical Properties Retrieval code to exit.' &
                 // char(10) // 'Operator Action:  Notify SDST. '
            
            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         end if
         Acloudmask_name = temp
         
         temp = ''
         Vrsn_No = 1
         LUN = Input_LRN(3)
         write(LUNSTR,'(I6)') LUN
         !         ........retrieve name of hdf product from pcf
         rtn = pgs_pc_getreference(LUN,Vrsn_No,temp)
         
         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.
            
            msgbuf = &
                 'pgs_pc_getreference detected error retrieving file with LUN# '// trim(LUNSTR) // ' from PCF. ' &
                 // char(10) // 'Cloud Optical Properties Retrieval code to exit.' &
                 // char(10) // 'Operator Action:  Notify SDST. '
            
            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         end if
         Ageolocation_name = temp

         temp = ''
         Vrsn_No = 1
         LUN = OUTPUT_LRN
         write(LUNSTR,'(I6)') LUN
         !         ........retrieve name of hdf product from pcf
         rtn = pgs_pc_getreference(LUN,Vrsn_No,temp)
         
         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.
            
            msgbuf = &
                 'pgs_pc_getreference detected error retrieving file with LUN# '// trim(LUNSTR) // ' from PCF. ' &
                 // char(10) // 'Cloud Optical Properties Retrieval code to exit.' &
                 // char(10) // 'Operator Action:  Notify SDST. '
            
            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         end if
         Amod06_name = temp
         
         temp = ''
         Vrsn_No = 1
         LUN = Input_LRN(4)
         write(LUNSTR,'(I6)') LUN
         !         ........retrieve name of hdf product from pcf
         rtn = pgs_pc_getreference(LUN,Vrsn_No,temp)
         
         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.
            
            msgbuf = &
                 'pgs_pc_getreference detected error retrieving 1st file with LUN# '// trim(LUNSTR) // ' from PCF. ' &
                 // char(10) // 'Cloud Optical Properties Retrieval code to exit.' &
                 // char(10) // 'Operator Action:  Notify SDST. '
            
            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         end if
         Agdas_name = temp         
         tempString = Agdas_name
         
         !Get the time of the gdas file
         datePos = SCAN(tempString, '.', .true.)
         tempDateTime1 = tempString(datePos+1:datePos+2)
         
         tempString = tempString(1:datePos-1)
         datePos = SCAN(tempString, '.', .true.)
         tempDateTime1 = tempString(datePos+1:datePos+6) // tempDateTime1
         read (tempDateTime1,*) iDateTime1
                  
         temp = ''
         Vrsn_No = 2
         LUN = Input_LRN(4)
         write(LUNSTR,'(I6)') LUN
         !         ........retrieve name of hdf product from pcf
         rtn = pgs_pc_getreference(LUN,Vrsn_No,temp)
         
         If (rtn .ne. PGS_S_SUCCESS) Then
            ! By gbritzol on 01/06/2014: Changed the logic so it runs with only one gdas file
            ! without going to run error for cases that only 1 gdas is provided because of 
            ! gdas structure file inconsistencies (Same implementation as the MOD_PRLCAT Alg17/Alg29).
            ! In case 1 gdas is provided the 2nd gdas file will be the same as the 1st one giving
            ! a warning in the LogStatus file that the 2nd file was not present and so the 2nd one
            ! will point to the 1st one.
            error_flag = .FALSE.
            
            msgbuf = &
                 'pgs_pc_getreference detected error retrieving 2nd file with LUN# '// trim(LUNSTR) // ' from PCF. ' &
                 // char(10) // 'Continuing process using the previously fetched gdas file' &
                 // char(10) // 'as the 2nd gdas file. '
            
            Call modis_smf_setdynamicmsg(MODIS_W_GENERIC,msgbuf,FUNCNAME)
            
            Agdas_name2 = Agdas_name
         else
            Agdas_name2 = temp
            tempString = Agdas_name2
            
            datePos = SCAN(tempString, '.', .true.)
            tempDateTime2 = tempString(datePos+1:datePos+2)
            
            tempString = tempString(1:datePos-1)
            datePos = SCAN(tempString, '.', .true.)
            tempDateTime2 = tempString(datePos+1:datePos+6) // tempDateTime2
            read (tempDateTime2,*) iDateTime2
            
            if(iDateTime2<iDatetime1) then
               temp = Agdas_name
               Agdas_name = Agdas_name2
               Agdas_name2 = temp
            end if
         end if
                  
         
#ifdef USE_TOAST
         temp = ''
         Vrsn_No = 1
         LUN = Input_LRN(5)
         write(LUNSTR,'(I6)') LUN
!         ........retrieve name of hdf product from pcf
         rtn = pgs_pc_getreference(LUN,Vrsn_No,temp)

         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.

            msgbuf = &
            'pgs_pc_getreference detected error retrieving file with LUN# '// trim(LUNSTR) // ' from PCF. ' &
            // char(10) // 'Cloud Optical Properties Retrieval code to exit.' &
            // char(10) // 'Operator Action:  Notify SDST. '
            
            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         end if
         Aozone_name = temp
#else
!	Aozone_name = "none"
#endif	

         temp = ''
         Vrsn_No = 1
         LUN = Input_LRN(6)
         write(LUNSTR,'(I6)') LUN
         !         ........retrieve name of hdf product from pcf
         rtn = pgs_pc_getreference(LUN,Vrsn_No,temp)
         
         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.
            
            msgbuf = &
                 'pgs_pc_getreference detected error retrieving file with LUN# '// trim(LUNSTR) // ' from PCF. ' &
                 // char(10) // 'Cloud Optical Properties Retrieval code to exit.' &
                 // char(10) // 'Operator Action:  Notify SDST. '
            
            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         end if
         Ancepice_name = temp
         
         temp = ''
         Vrsn_No = 1
         LUN = Input_LRN(7)
         write(LUNSTR,'(I6)') LUN
         !         ........retrieve name of hdf product from pcf
         rtn = pgs_pc_getreference(LUN,Vrsn_No,temp)
         
         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.
            
            msgbuf = &
                 'pgs_pc_getreference detected error retrieving file with LUN# '// trim(LUNSTR) // ' from PCF. ' &
                 // char(10) // 'Cloud Optical Properties Retrieval code to exit.' &
                 // char(10) // 'Operator Action:  Notify SDST. '
            
            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         end if
         Anise_name = temp
         
         temp = ''
         Vrsn_No = 1
         LUN = Input_LRN(8)
         write(LUNSTR,'(I6)') LUN
!         ........retrieve name of hdf product from pcf
         rtn = pgs_pc_getreference(LUN,Vrsn_No,temp)

         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.

            msgbuf = &
            'pgs_pc_getreference detected error retrieving file with LUN# '// trim(LUNSTR) // ' from PCF. ' &
            // char(10) // 'Cloud Optical Properties Retrieval code to exit.' &
            // char(10) // 'Operator Action:  Notify SDST. '
            
            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         end if
         Aice_library = temp

         temp = ''
         Vrsn_No = 1
         LUN = Input_LRN(9)
         write(LUNSTR,'(I6)') LUN
         !         ........retrieve name of hdf product from pcf
         rtn = pgs_pc_getreference(LUN,Vrsn_No,temp)
         
         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.
            
            msgbuf = &
                 'pgs_pc_getreference detected error retrieving file with LUN# '// trim(LUNSTR) // ' from PCF. ' &
                 // char(10) // 'Cloud Optical Properties Retrieval code to exit.' &
                 // char(10) // 'Operator Action:  Notify SDST. '
            
            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         end if
         Awater_library = temp
         
         
         temp = ''
         Vrsn_No = 1
         LUN = Input_LRN(10)
         write(LUNSTR,'(I6)') LUN
         !         ........retrieve name of hdf product from pcf
         rtn = pgs_pc_getreference(LUN,Vrsn_No,temp)
         
         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.
            
            msgbuf = &
                 'pgs_pc_getreference detected error retrieving file with LUN# '// trim(LUNSTR) // ' from PCF. ' &
                 // char(10) // 'Cloud Optical Properties Retrieval code to exit.' &
                 // char(10) // 'Operator Action:  Notify SDST. '
            
            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         end if
         Alibnames_ice(1) = temp
         
         temp = ''
         Vrsn_No = 1
         LUN = Input_LRN(11)
         write(LUNSTR,'(I6)') LUN
         !         ........retrieve name of hdf product from pcf
         rtn = pgs_pc_getreference(LUN,Vrsn_No,temp)
         
         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.
            
            msgbuf = &
                 'pgs_pc_getreference detected error retrieving file with LUN# '// trim(LUNSTR) // ' from PCF. ' &
                 // char(10) // 'Cloud Optical Properties Retrieval code to exit.' &
                 // char(10) // 'Operator Action:  Notify SDST. '
            
            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         end if
         Alibnames_ice_sdev(1) = temp
         
         temp = ''
         Vrsn_No = 1
         LUN = Input_LRN(12)
         write(LUNSTR,'(I6)') LUN
         !         ........retrieve name of hdf product from pcf
         rtn = pgs_pc_getreference(LUN,Vrsn_No,temp)
         
         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.
            
            msgbuf = &
                 'pgs_pc_getreference detected error retrieving file with LUN# '// trim(LUNSTR) // ' from PCF. ' &
                 // char(10) // 'Cloud Optical Properties Retrieval code to exit.' &
                 // char(10) // 'Operator Action:  Notify SDST. '
            
            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         end if
         Alibnames_ice(2) = temp
         
         
         temp = ''
         Vrsn_No = 1
         LUN = Input_LRN(13)
         write(LUNSTR,'(I6)') LUN
         !         ........retrieve name of hdf product from pcf
         rtn = pgs_pc_getreference(LUN,Vrsn_No,temp)
         
         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.
            
            msgbuf = &
                 'pgs_pc_getreference detected error retrieving file with LUN# '// trim(LUNSTR) // ' from PCF. ' &
                 // char(10) // 'Cloud Optical Properties Retrieval code to exit.' &
                 // char(10) // 'Operator Action:  Notify SDST. '
            
            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         end if
         Alibnames_ice_sdev(2) = temp
         
         temp = ''
         Vrsn_No = 1
         LUN = Input_LRN(14)
         write(LUNSTR,'(I6)') LUN
         !         ........retrieve name of hdf product from pcf
         rtn = pgs_pc_getreference(LUN,Vrsn_No,temp)
         
         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.
            
            msgbuf = &
                 'pgs_pc_getreference detected error retrieving file with LUN# '// trim(LUNSTR) // ' from PCF. ' &
                 // char(10) // 'Cloud Optical Properties Retrieval code to exit.' &
                 // char(10) // 'Operator Action:  Notify SDST. '
            
            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         end if
         Alibnames_ice(3) =  temp
         
         
         temp = ''
         Vrsn_No = 1
         LUN = Input_LRN(15)
         write(LUNSTR,'(I6)') LUN
         !         ........retrieve name of hdf product from pcf
         rtn = pgs_pc_getreference(LUN,Vrsn_No,temp)
         
         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.
            
            msgbuf = &
                 'pgs_pc_getreference detected error retrieving file with LUN# '// trim(LUNSTR) // ' from PCF. ' &
                 // char(10) // 'Cloud Optical Properties Retrieval code to exit.' &
                 // char(10) // 'Operator Action:  Notify SDST. '
            
            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         end if
         Alibnames_ice_sdev(3) = temp
         
         
         temp = ''
         Vrsn_No = 1
         LUN = Input_LRN(16)
         write(LUNSTR,'(I6)') LUN
         !         ........retrieve name of hdf product from pcf
         rtn = pgs_pc_getreference(LUN,Vrsn_No,temp)
         
         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.
            
            msgbuf = &
                 'pgs_pc_getreference detected error retrieving file with LUN# '// trim(LUNSTR) // ' from PCF. ' &
                 // char(10) // 'Cloud Optical Properties Retrieval code to exit.' &
                 // char(10) // 'Operator Action:  Notify SDST. '
            
            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         end if
         Alibnames_water(1) =  temp
         
         temp = ''
         Vrsn_No = 1
         LUN = Input_LRN(17)
         write(LUNSTR,'(I6)') LUN
         !         ........retrieve name of hdf product from pcf
         rtn = pgs_pc_getreference(LUN,Vrsn_No,temp)
         
         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.
            
            msgbuf = &
                 'pgs_pc_getreference detected error retrieving file with LUN# '// trim(LUNSTR) // ' from PCF. ' &
                 // char(10) // 'Cloud Optical Properties Retrieval code to exit.' &
                 // char(10) // 'Operator Action:  Notify SDST. '
            
            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         end if
         Alibnames_water_sdev(1) =  temp
         
         temp = ''
         Vrsn_No = 1
         LUN = Input_LRN(18)
         write(LUNSTR,'(I6)') LUN
         !         ........retrieve name of hdf product from pcf
         rtn = pgs_pc_getreference(LUN,Vrsn_No,temp)
         
         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.
            
            msgbuf = &
                 'pgs_pc_getreference detected error retrieving file with LUN# '// trim(LUNSTR) // ' from PCF. ' &
                 // char(10) // 'Cloud Optical Properties Retrieval code to exit.' &
                 // char(10) // 'Operator Action:  Notify SDST. '
            
            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         end if
         Alibnames_water(2) =  temp
         
         
         temp = ''
         Vrsn_No = 1
         LUN = Input_LRN(19)
         write(LUNSTR,'(I6)') LUN
         !         ........retrieve name of hdf product from pcf
         rtn = pgs_pc_getreference(LUN,Vrsn_No,temp)
         
         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.
            
            msgbuf = &
                 'pgs_pc_getreference detected error retrieving file with LUN# '// trim(LUNSTR) // ' from PCF. ' &
                 // char(10) // 'Cloud Optical Properties Retrieval code to exit.' &
                 // char(10) // 'Operator Action:  Notify SDST. '
            
            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         end if
         Alibnames_water_sdev(2) = temp
         
         temp = ''
         Vrsn_No = 1
         LUN = Input_LRN(20)
         write(LUNSTR,'(I6)') LUN
         !         ........retrieve name of hdf product from pcf
         rtn = pgs_pc_getreference(LUN,Vrsn_No,temp)
         
         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.
            
            msgbuf = &
                 'pgs_pc_getreference detected error retrieving file with LUN# '// trim(LUNSTR) // ' from PCF. ' &
                 // char(10) // 'Cloud Optical Properties Retrieval code to exit.' &
                 // char(10) // 'Operator Action:  Notify SDST. '
            
            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         end if
         Alibnames_water(3) = temp
         
         temp = ''
         Vrsn_No = 1
         LUN = Input_LRN(21)
         write(LUNSTR,'(I6)') LUN
         !         ........retrieve name of hdf product from pcf
         rtn = pgs_pc_getreference(LUN,Vrsn_No,temp)
         
         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.
            
            msgbuf = &
                 'pgs_pc_getreference detected error retrieving file with LUN# '// trim(LUNSTR) // ' from PCF. ' &
                 // char(10) // 'Cloud Optical Properties Retrieval code to exit.' &
                 // char(10) // 'Operator Action:  Notify SDST. '
            
            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         end if
         Alibnames_water_sdev(3) =  temp
         
         temp = ''
         Vrsn_No = 1
         LUN = Input_LRN(22)
         write(LUNSTR,'(I6)') LUN
         !         ........retrieve name of hdf product from pcf
         rtn = pgs_pc_getreference(LUN,Vrsn_No,temp)
         
         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.
            
            msgbuf = &
                 'pgs_pc_getreference detected error retrieving file with LUN# '// trim(LUNSTR) // ' from PCF. ' &
                 // char(10) // 'Cloud Optical Properties Retrieval code to exit.' &
                 // char(10) // 'Operator Action:  Notify SDST. '
            
            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         end if
         Aphase_library =  temp
         
         temp = ''
         Vrsn_No = 1
         LUN = Input_LRN(23)
         write(LUNSTR,'(I6)') LUN
         !         ........retrieve name of hdf product from pcf
         rtn = pgs_pc_getreference(LUN,Vrsn_No,temp)
         
         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.
            
            msgbuf = &
                 'pgs_pc_getreference detected error retrieving file with LUN# '// trim(LUNSTR) // ' from PCF. ' &
                 // char(10) // 'Cloud Optical Properties Retrieval code to exit.' &
                 // char(10) // 'Operator Action:  Notify SDST. '
            
            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         end if
         Atransmittance_library = temp

         temp = ''
         Vrsn_No = 1
         LUN = Input_LRN(24)
         write(LUNSTR,'(I6)') LUN
         !         ........retrieve name of hdf product from pcf
         rtn = pgs_pc_getreference(LUN,Vrsn_No,temp)
         
         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.
            
            msgbuf = &
                 'pgs_pc_getreference detected error retrieving file with LUN# '// trim(LUNSTR) // ' from PCF. ' &
                 // char(10) // 'Cloud Optical Properties Retrieval code to exit.' &
                 // char(10) // 'Operator Action:  Notify SDST. '
            
            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         end if
         Aecosystem_data_name = temp
         
         temp = ''
         Vrsn_No = 1
         LUN = Input_LRN(25)
         write(LUNSTR,'(I6)') LUN
         !         ........retrieve name of hdf product from pcf
         rtn = pgs_pc_getreference(LUN,Vrsn_No,temp)
         
         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.
            
            msgbuf = &
                 'pgs_pc_getreference detected error retrieving file with LUN# '// trim(LUNSTR) // ' from PCF. ' &
                 // char(10) // 'Cloud Optical Properties Retrieval code to exit.' &
                 // char(10) // 'Operator Action:  Notify SDST. '
            
            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         end if
         Asnowicealbedo_data_name =  temp
         
         temp = ''
         Vrsn_No = 1
         LUN = Input_LRN(26)
         write(LUNSTR,'(I6)') LUN
         !         ........retrieve name of hdf product from pcf
         rtn = pgs_pc_getreference(LUN,Vrsn_No,temp)
         
         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.
            
            msgbuf = &
                 'pgs_pc_getreference detected error retrieving file with LUN# '// trim(LUNSTR) // ' from PCF. ' &
                 // char(10) // 'Cloud Optical Properties Retrieval code to exit.' &
                 // char(10) // 'Operator Action:  Notify SDST. '
            
            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         end if
         Asurfacealbedo_lib_659 = temp
         
         temp = ''
         Vrsn_No = 1
         LUN = Input_LRN(27)
         write(LUNSTR,'(I6)') LUN
         !         ........retrieve name of hdf product from pcf
         rtn = pgs_pc_getreference(LUN,Vrsn_No,temp)
         
         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.
            
            msgbuf = &
                 'pgs_pc_getreference detected error retrieving file with LUN# '// trim(LUNSTR) // ' from PCF. ' &
                 // char(10) // 'Cloud Optical Properties Retrieval code to exit.' &
                 // char(10) // 'Operator Action:  Notify SDST. '
            
            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         end if
         Asurfacealbedo_lib_858 = temp
         
         temp = ''
         Vrsn_No = 1
         LUN = Input_LRN(28)
         write(LUNSTR,'(I6)') LUN
         !         ........retrieve name of hdf product from pcf
         rtn = pgs_pc_getreference(LUN,Vrsn_No,temp)
         
         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.
            
            msgbuf = &
                 'pgs_pc_getreference detected error retrieving file with LUN# '// trim(LUNSTR) // ' from PCF. ' &
                 // char(10) // 'Cloud Optical Properties Retrieval code to exit.' &
                 // char(10) // 'Operator Action:  Notify SDST. '
            
            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         end if
         Asurfacealbedo_lib_124 =  temp
         
         temp = ''
         Vrsn_No = 1
         LUN = Input_LRN(29)
         write(LUNSTR,'(I6)') LUN
         !         ........retrieve name of hdf product from pcf
         rtn = pgs_pc_getreference(LUN,Vrsn_No,temp)
         
         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.
            
            msgbuf = &
                 'pgs_pc_getreference detected error retrieving file with LUN# '// trim(LUNSTR) // ' from PCF. ' &
                 // char(10) // 'Cloud Optical Properties Retrieval code to exit.' &
                 // char(10) // 'Operator Action:  Notify SDST. '
            
            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         end if
         Asurfacealbedo_lib_164 = temp
         
         temp = ''
         Vrsn_No = 1
         LUN = Input_LRN(30)
         write(LUNSTR,'(I6)') LUN
         !         ........retrieve name of hdf product from pcf
         rtn = pgs_pc_getreference(LUN,Vrsn_No,temp)
         
         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.
            
            msgbuf = &
                 'pgs_pc_getreference detected error retrieving file with LUN# '// trim(LUNSTR) // ' from PCF. ' &
                 // char(10) // 'Cloud Optical Properties Retrieval code to exit.' &
                 // char(10) // 'Operator Action:  Notify SDST. '
            
            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         end if
         Asurfacealbedo_lib_21 = temp
         
         
         LUN = 413155
         temp = ''
         write(LUNSTR,'(I6)') LUN
         !         ........retrieve name of hdf product from pcf
         rtn = pgs_pc_getconfigdata(LUN,temp)
         
         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.
            
            msgbuf = &
                 'pgs_pc_getconfigdata detected error with LUN# '// trim(LUNSTR) // ' from PCF. ' &
                 // char(10) // 'Cloud Optical Properties Retrieval code to exit.' &
                 // char(10) // 'Operator Action:  Notify SDST. '
            
            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         end if
         ACT_lib_path = temp
         
         temp = ''
         Vrsn_No = 1
         LUN = Input_LRN(31)
         write(LUNSTR,'(I6)') LUN
         !         ........retrieve name of hdf product from pcf
         rtn = pgs_pc_getreference(LUN,Vrsn_No,temp)
         
         If (rtn .ne. PGS_S_SUCCESS) Then
            error_flag = .TRUE.
            
            msgbuf = &
                 'pgs_pc_getreference detected error retrieving file with LUN# '// trim(LUNSTR) // ' from PCF. ' &
                 // char(10) // 'Cloud Optical Properties Retrieval code to exit.' &
                 // char(10) // 'Operator Action:  Notify SDST. '
            
            Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
         end if
         Aemissivity_name = temp


#else
!--------------------------------------------------------------------------------
! Begin executable code by writing "Begin" message to LogStatus file 
!--------------------------------------------------------------------------------
	
	mylun = 8


	call getarg(1, par_file)
	print*, par_file

	open (unit=mylun, file=trim(par_file), status='old', form='formatted')
	read(mylun, *) dummy
	read(mylun, *) yearname, dayname, timename
	
	read(yearname,*) MYYEAR
	read(dayname, *) MYDAY
	read(timename, *) MYTIME
	dirname = './inputs-' // trim(yearname) // trim(dayname) // '.' // trim(timename) // '/'
	outdirname = './outputs-' // trim(yearname) // trim(dayname) // '.' // trim(timename) // '/OD/'
	print*, trim(dirname) 
	
	read(mylun, *) dummy
	read(mylun, *) temp
	Alevel1b_name(1) = trim(dirname) // trim(temp)
	print*, trim(Alevel1b_name(1))

	read(mylun, *) dummy
	read(mylun, *) temp
	Acloudmask_name = trim(dirname) // trim(temp)
	print*, trim(Acloudmask_name)
	
	read(mylun, *) dummy
	read(mylun, *) temp
	Ageolocation_name = trim(dirname) // trim(temp)
	print*, trim(Ageolocation_name)
	
	
	read(mylun, *) dummy
	read(mylun, *) temp
	Amod06_name = trim(outdirname) // trim(temp)
	print*, trim(Amod06_name)
!	MY_TEXT_FILE = trim(Amod06_name) // ".txt"
	
	read(mylun, *) dummy
	read(mylun, *) temp
	Agdas_name = trim(dirname) // trim(temp)
	read(mylun, *) temp
	Agdas_name2 = trim(dirname) // trim(temp)
	print*, trim(Agdas_name)
	print*, trim(Agdas_name2)

#ifdef USE_TOAST
	read(mylun, *) dummy
	read(mylun, *) temp
	Aozone_name = trim(dirname) // trim(temp)
	print*, trim(Aozone_name)
#else
!	Aozone_name = "none"
#endif	
	
	read(mylun, *) dummy
	read(mylun, *) temp
	Ancepice_name = trim(dirname) // trim(temp)
	print*, trim(Ancepice_name)
	
	read(mylun, *) dummy
	read(mylun, *) temp
	Anise_name = trim(dirname) // trim(temp)
	print*, trim(Anise_name)

! skip the CSR bias file name if present
	
	read(mylun, *) dummy
	print*, trim(dummy)
	idx = index (dummy, "CSR")
	if (idx /= 0) then 
		read(mylun,*) dummy
		read(mylun,*) dummy
	endif


	read(mylun, '(a)') temp
	Aice_library =  trim(temp)
	read(mylun, '(a)') temp
	Awater_library = trim(temp)
	
	do i=1, 3

		read(mylun, '(a)') temp
		Alibnames_ice(i) =  trim(temp)
		read(mylun, '(a)') temp
		Alibnames_ice_sdev(i) =  trim(temp)
	
	end do
	do i=1, 3

		read(mylun, '(a)') temp
		Alibnames_water(i) =  trim(temp)
		read(mylun, '(a)') temp
		Alibnames_water_sdev(i) =  trim(temp)
	
	end do


	read(mylun, '(a)') temp
	Aphase_library = trim(temp)
    print*, "ICE: ", trim(Aice_library)
	print*, trim(Awater_library)
	do i=1, 3
		print*, trim(Alibnames_ice(i))
		print*, trim(Alibnames_ice_sdev(i))
		print*, trim(Alibnames_water(i))
		print*, trim(Alibnames_water_sdev(i))
	end do

	print*, trim(Aphase_library)
	
	read(mylun, *) dummy
	read(mylun, '(a)') temp
	Atransmittance_library =  trim(temp)
	print*, trim(Atransmittance_library)
	
	read(mylun, *) dummy
	read(mylun, '(a)') temp
	Aecosystem_data_name =  trim(temp)
	print*, trim(Aecosystem_data_name)

	read(mylun, *) dummy
	read(mylun, '(a)') temp
	Asnowicealbedo_data_name =  trim(temp)
	print*, trim(Asnowicealbedo_data_name)

	read(mylun, *) dummy
	read(mylun, '(a)') temp
	Asurfacealbedo_lib_659 =  trim(temp)
	read(mylun, '(a)') temp
	Asurfacealbedo_lib_858 =  trim(temp)
	read(mylun, '(a)') temp
	Asurfacealbedo_lib_124 =  trim(temp)
	read(mylun, '(a)') temp
	Asurfacealbedo_lib_164 =  trim(temp)
	read(mylun, '(a)') temp
	Asurfacealbedo_lib_21 =  trim(temp)

	read(mylun, *) dummy
	read(mylun, '(a)') temp
	ACT_lib_path = trim(temp)
	read(mylun, *) dummy
	read(mylun, '(a)') temp
	Aemissivity_name = trim(temp)

	
	print*, trim(Asurfacealbedo_lib_659)
	print*, trim(Asurfacealbedo_lib_858)
	print*, trim(Asurfacealbedo_lib_124)
	print*, trim(Asurfacealbedo_lib_164)
	print*, trim(Asurfacealbedo_lib_21)
	
	
	close(mylun)

#endif


!...........Perform cloud optical properties retrieval


               call MOD_PR06OD(statistics,ErrorFlag)

#ifdef USE_PGS
        ExtGeoPntr_Flag = .TRUE.
        Input_Vrsn(4) = 2
        
        rtn = opmfil(Amod06_name, 'a', Modfil_Out)
        
        If (rtn .ne. MAPIOK) Then
           error_flag = .TRUE.
           
           msgbuf = &
                'OPMFIL detected error reopening Cloud Product.  Halt Processing.' &
                // char(10) // 'Operator Action:  Notify SDST. '
           
           Call modis_smf_setdynamicmsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)
        Else
           
           
           !..............write ECS inventory metadata
           rtn = Set_InvMet_MOD06(ExtGeoPntr_Flag, &
                NUM_Input_LRN,Input_LRN,Input_Vrsn, &
                NUM_MeasParam,Name_MeasParm,QA_Miss, &
                Auto_QA_Flag,Auto_QA_Expl, &
                NUM_PSA,PSAName,PSAValue, &
                Met_Handles) 
           
           If (rtn .eq. FAIL) Then
              error_flag = .TRUE.
              
              msgbuf = &
                   'Set_InvMet_MOD06 detected error attempting to ' &
                   // char(10) // 'update ECS inventory metadata. ' &
                   // char(10) // 'Operator Action:  Notify SDST. ' 
              
              Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
           end if


           !  -----------------------------------------------------------------------
           !  Set the doi Attributes (PSAs).
           !  -----------------------------------------------------------------------

           NUM_PSA = 17
           !     Get satellite instrument name.
           rtn = pgs_pc_getconfigdata(LUN_Sat_Instrument,pcf_satid)
           if (rtn .ne. 0) then
              msgbuf = &
                   'Error reading instrument name from pcf LUN 800510.' // &
                   ' [OPERATOR ACTION: Check PCF file, fix if corrupt.]'
              Call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,msgbuf,FUNCNAME)
           endif
!          Get ReprocessingActual value. 
           rtn = pgs_pc_getconfigdata(LUN_Rep_Actual,proctype)   
           if (rtn .ne. 0) then
              msgbuf = &
                  'Error reading ReprocessingActual from pcf LUN 800504.' // &
                  '[OPERATOR ACTION: Check PCF file, fix if corrupt.]'
              Call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,msgbuf,FUNCNAME)
           endif

           if (pcf_satid .eq. 'AM1M') then
             if (proctype(1:1) .eq. 'N') then
                doi_char = '10.5067/MODIS/MOD06_L2.NRT.061'
             else
                doi_char = '10.5067/MODIS/MOD06_L2.061'
             endif
           else
             if (proctype(1:1) .eq. 'N') then
                doi_char = '10.5067/MODIS/MYD06_L2.NRT.061'
             else
                doi_char = '10.5067/MODIS/MYD06_L2.061'
             endif
           endif
           
           rtn = Set_PSA(Met_Handles, 'identifier_product_doi', NUM_PSA+1, doi_char)
           If (rtn .eq. FAIL) Then
              lErrorFlag = .true.
              
              msgbuf = 'Set_PSA detected error setting DOI PSAs.'  // &
                   'Operator Action:  Refer to prior low level LogStatus error ' // &
                   'messages originating in function Set_InvPSA_Atmos.'
              
              Call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,msgbuf,FUNCNAME)
           EndIf
           
           rtn = Set_PSA(Met_Handles, 'identifier_product_doi_authority', NUM_PSA+2, 'https://doi.org')
           If (rtn .eq. FAIL) Then
              lErrorFlag = .true.
              
              msgbuf = 'Set_PSA detected error setting DOI PSAs.' // &
                   'Operator Action:  Refer to prior low level LogStatus error ' // &
                   'messages originating in function Set_InvPSA_Atmos.'
              
              Call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,msgbuf,FUNCNAME)
           EndIf
           
           Num_of_HDfAttrNms = 1
           HDFAttNms(INVENTORYMETADATA) = MECS_CORE
           !               HDFAttNms(ARCHIVEMETADATA)   = MECS_ARCHIVE
           
           rtn = CPMFIL(Modfil_Out,Met_Handles,HDFAttNms,Num_Of_HDFAttrNms)
           
           If (rtn .ne. MAPIOK) Then 
              error_flag = .TRUE.
              
              msgbuf = &
                   'CPMFIL detected error closing Cloud Product HDF file. ' &
                   // char(10) // 'Operator Action:  Notify SDST. '
              
              Call MODIS_SMF_SETDYNAMICMSG(MODIS_E_GENERIC,msgbuf,FUNCNAME)
           Else
              msgbuf = 'CPMFIL successfully closed Cloud Product HDF file. '
              
              Call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,msgbuf,FUNCNAME)
           End If



           rtn = -1
           rtn = opmfil( Amod06_name, 'a', Modfil_Out )
           if ( rtn .ne. MAPIOK ) then
              msgbuf = 'Error opening output MOD06_L2 file.' // &
                   ' [OPERATOR ACTION: Check system resources]'
              Call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,msgbuf,FUNCNAME)
           endif
           
           doiLen = 26
           if (proctype(1:1) .eq. 'N') then
             doiLen = 30
           endif
           rtn = sfscatt(Modfil_Out, 'identifier_product_doi', DFNT_CHAR8, doiLen, doi_char) 
           if (rtn .eq. -1) then
              msgbuf = 'Problem writting the global attribute identifier_product_doi'
              Call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,msgbuf,FUNCNAME)
           endif
           
           rtn = sfscatt(Modfil_Out, 'identifier_product_doi_authority', DFNT_CHAR8, 17,'https://doi.org') 
           if (rtn .eq. -1) then
              msgbuf = 'Problem writting the global attribute identifier_product_doi_authority'
              Call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,msgbuf,FUNCNAME)
           endif
           
           ! ... Close the output file
           
           rtn = clmfil(Modfil_Out)
           if( rtn .eq. 0 ) then
              msgbuf = 'Success closing output MOD_PR06 file '
              Call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,msgbuf,FUNCNAME)
           else
              msgbuf = 'Error closing output MOD_PR06 file.' // &
                   '[OPERATOR ACTION: Check system resources]'
              Call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC,msgbuf,FUNCNAME)
           endif
           
        end if

#endif

      End program Driver_MOD_PR06OD
