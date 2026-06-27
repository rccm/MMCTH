      Function Initial_Input_Checks(L1B_MAPI_fhandle,L1B_fname,
     &                              Band_No,Gain,Cal_Type,Resol)
      implicit none
      include 'mapi.inc'
      include 'hdf.inc'
      include 'PGS_MODIS_39500.f'
      include 'L1B_Reader.inc'
  
C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION: 
C
C  Perform quality assurance testing on the input arguments passed
C  to function Initial_Input_Checks.  These tests come mainly in the 
C  form of out of bounds and valid range checking.
C
C  Function Initial_Input_Checks separately compares the input 
C  arguments L1B_MAPI_fhandle, Gain, Cal_type and Resol against a 
C  discrete set of permissible values.  If for each case no match is 
C  found, an error message is reported and an error flag is set.  
C  Out of bounds checking is performed on arguments L1B_fname and 
C  Band_No.  If either satisfies the out of bounds criteria, an error 
C  message is reported to the SDP toolkit LogStatus file. 
C
C !INPUT PARAMETERS:
C
C  integer L1B_MAPI_fhandle  Array containing HDF SD (element 1) and 
C                            VS (element 2) interface IDs, and 
C                            the L1B file access mode (element 3).  
C                            Additional elements carrying pointers to 
C                            C language structures may also be passed, 
C                            and are used internally by the MAPI 
C                            library to maximize CPU efficiency.
C  character*(*) Cal_type    Calibration type applied to stored integers.
C                            Retrieve data as reflectances ("refl"),
C                            radiances ("rad") or 
C                            corrected counts ("counts").
C  character*(*) Gain        Band radiometric gain factor, 
C                            high ("H") or low ("L")
C  character*(*) L1B_fname   Variable containing full path name of 
C                            L1B data file
C  integer   Band_No         MODIS spectral band number from 1 to 36
C  integer   Resol           Data spatial resolution factor relative to 
C                            1 km.  Possible values are 16 (250 m), 
C                            4 (500 m) and 1 (1 km).
C
C !OUTPUT PARAMETERS:      NONE
C
C !REVISION HISTORY:
C Revision 1.5  1998/06/05  15:28:15  rhucek
C Modified error messages and description of input parameters
C
C Revision 1.4  1998/05/20  14:09:07  rhucek
C Update error messages with "Operator Action" strings.
C
C Revision 1.3  1997/12/15  17:01:23  fhliang
C changed 'mapic.inc' to 'mapi.inc', and added 'hdf.inc'.
C
C Revision 1.2  1997/08/12  13:13:50  rhucek
C Added explicitly type declarations for dummy arguments
C  Gain and Cal_type [character*(*)], and Band_No and Resol
C (integer).
C
C Changed declaration of dummy argument L1B_MAPI_fhandle
C from fixed to assumed-size array.
C
C Revision 1.1  1997/04/15  14:13:56  vlin
C Initial revision
C
C !TEAM-UNIQUE HEADER:
C
C    This software is developed by the MODIS Science Data Support
C    Team for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C !REFERENCES AND CREDITS:
C
C    Written by    
C    S. Vicky Lin                  March 1997 
C    *** Under guidance of Richard Hucek  ***
C
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd, Seabrook MD 20706
C    rhucek@ltpmail.gsfc.nasa.gov
C    vlin@modis1.gsfc.nasa.gov
C
C !DESIGN NOTES:
C  
C  Return Value:     0 for success; -1 for fail
C
C  Externals:
C    Named Constants:
C        Band_Is_*         (L1B_Reader.inc)
C        Cal_Is_*          (L1B_Reader.inc)
C        DFACC_READ        (hdf.inc)
C        FAIL              (hdf.inc)
C        First_*_Band      (L1B_Reader.inc)
C        High_Gain         (L1B_Reader.inc)
C        Last_*_Band       (L1B_Reader.inc)
C        Low_Gain          (L1B_Reader.inc)
C        P_SDID            (mapi.inc)
C        P_ACCESS          (mapi.inc)
C        Resol_Is_*        (L1B_Reader.inc)
C        SUCCEED           (hdf.inc)
C
C  Internals:
C    Subroutines:
C        MODIS_SMF_SetDynamicMsg
C
C !END
C---------------------------------------------------------------------

c PARAMETER declarations
      character*(*) FUNCNAME
      PARAMETER    (FUNCNAME = 'Initial_Input_Checks')

      character*(*) Cal_Type,Gain,L1B_fname
      character*6   msg6
      character*135 msgbuf

      integer Band_No,L1B_MAPI_fhandle(*),Resol
      integer string_length,fbyte,lbyte
      integer Initial_Input_Checks

c ... Initialization

      Initial_Input_Checks = succeed
      string_length = -1

c ... Check file access to L1B dataset

      if (L1B_MAPI_fhandle(P_SDID).le.0) then
          Initial_Input_Checks = fail
          msgbuf = 'L1B file HDF SD interface index is not set. ' 
     1             // char(10) //  'READ_L1B cannot access file.'
     2             // char(10) // 'Operator Action:  Notify SDST.'


          call Modis_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

      else if (L1B_MAPI_fhandle(P_ACCESS).ne.DFACC_READ) then
          Initial_Input_Checks = fail
          msgbuf = 'L1B file not opened for "read-only" access.' 
     1             // char(10) // 'Operator Action:  Notify SDST.'

          call Modis_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

      end if

c ... Check for file name out of bounds

      call CUT_NAME_L1B(L1B_fname,fbyte,lbyte)
      string_length = lbyte - fbyte
      if (string_length.le.0) then
          Initial_Input_Checks = fail
          msgbuf = 'L1B file name is a blank string, READ_L1B unable to continue.'
     1             // char(10) // 'Operator Action:  Notify SDST.'

          call Modis_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

      end if

c ... Check for acceptable value of L1B spatial resolution factor.

      if (Resol.eq.Resol_Is_250m) then

         if (Band_No.lt.First_250m_Band .or.
     1       Band_No.gt.Last_250m_Band) then
             Initial_Input_Checks = fail
             write(msg6,'(i6)') Band_No
             msgbuf = 'MODIS L1B 250M band set does not include band' // msg6 
     1                // char(10) // 'Operator Action:  Notify SDST.'

             call Modis_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

         end if

      else if (Resol.eq.Resol_Is_500m) then

         if (Band_No.lt.First_250m_Band .or.
     1       Band_No.gt.Last_500m_Band) then
             Initial_Input_Checks = fail
             write(msg6,'(i6)') Band_No
             msgbuf = 'MODIS L1B 500M band set does not include band' // msg6 
     1                // char(10) // 'Operator Action:  Notify SDST.'

             call Modis_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

         end if

      else if (Resol.eq.Resol_Is_1km) then

         if (Band_No.lt.First_Modis_Band .or.
     1       Band_No.gt.Last_Modis_Band) then
             Initial_Input_Checks = fail
             write(msg6,'(i6)') Band_No
             msgbuf = 'MODIS L1B 1KM band set does not include band' // msg6 
     1                // char(10) // 'Operator Action:  Notify SDST.'

             call Modis_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

         else if (Band_No.eq.Band_Is_13.or.Band_No.eq.Band_Is_14) then 
              if (Gain.ne.High_Gain.and.Gain.ne.Low_Gain) then
                 Initial_Input_Checks = fail
                 write(msg6,'(i6)') Band_No
                 msgbuf = 'Gain is not one of "H" or "L" for MODIS L1B band' // msg6
     1                    // char(10) // 'Operator Action:  Notify SDST.'

                 call Modis_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

              end if
         end if
 
      else

         Initial_Input_Checks = fail
         write(msg6,'(i6)') Resol
         msgbuf = 'The resolution factor is not set to 16, 4, or 1: = ' // msg6
     1            // char(10) // 'Operator Action:  Notify SDST.'

         call Modis_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

      end if

c ... Check for acceptable calibration type

      if (Cal_Type.eq.Cal_Is_Refl) then

          if (Band_No.ge.First_Thermal_Band.and.
     1        Band_No.le.Last_Thermal_Band .and.
     2        Band_No.ne.Last_Solar_Band) then
              Initial_Input_Checks = fail
              write(msg6,'(i6)') Band_No

              msgbuf = 'MODIS L1B reflectance band set does not include band' // msg6 
     1                 // char(10) // 'Operator Action:  Notify SDST.'

              call Modis_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

          end if

      else if (Cal_Type.ne.Cal_Is_Rad   .and.
     1         Cal_Type.ne.Cal_Is_Count) then

              Initial_Input_Checks = fail
              msgbuf = 'Calibration type is not one of "Refl", "Rad" or "Count"'
     1                 // char(10) // 'Operator Action:  Notify SDST.'

              call Modis_SMF_SetDynamicMsg(MODIS_E_GENERIC,msgbuf,FUNCNAME)

      end if

      return
      end

c***********************************************************************

      SUBROUTINE CUT_NAME_L1B(string, fbyte, lbyte)
      IMPLICIT NONE

C----------------------------------------------------------------------
C !F77
C
C !Description:  
C
C    CUT_NAME_L1B finds the position (first and last bytes) of the file 
C    name within a string buffer which contains both path (file location) 
C    and file name. 
C
C !INPUT PARAMETERS:
C
C    character string  A string variable of arbitrary length which 
C                      shall consist of a unix path and file name.
C
C !OUTPUT PARAMETERS:
C
C    integer fbyte     The byte location of the first nonblank character 
C                      of the input string.
C    integer lbyte     The byte location of the last nonblank character 
C                      of the input string.
C
C !REVISION HISTORY:
c Revision 1.2  1997/02/25  18:55:14  vlin
c Changed from function to a subroutine
C
C !TEAM-UNIQUE HEADER:
C
C    This software was developed by the MODIS Science Data Support Team
C    (SDST) for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-32373.
C
C !REFERENCES AND CREDITS
C
C    Written by Richard Hucek
C
C    Research and Data systems Corporation
C    SAIC/GSC MODIS Science Data Support Office
C    7501 Forbes Blvd, Seabrook MD 20706
C
C    rhucek@ltpmail.gsfc.nasa.gov
C
C !DESIGN NOTES:
C
C    There is no check to identify the unpredictable value of an 
C    undefined string.  Consequently, users must take care to initialize 
C    all string variables before passing them to CUT_NAME.
C
C !END
C-----------------------------------------------------------------------

      CHARACTER*(*) string
      INTEGER fbyte, lbyte, string_len

C Initialize variables
      string_len=len(string)
      fbyte=1
      lbyte = string_len

C Determine byte position of last non-blank, non-slash character.
C This is last character in file name
      DO WHILE ( (string(lbyte:lbyte).eq.' ').and.(lbyte.ge.1) )
         lbyte=lbyte-1
      END DO

C Determine byte position of last slash (/) in character string
      fbyte = lbyte - 1

      DO WHILE ( (fbyte.ge.1) .and. (string(fbyte:fbyte).ne.'/') )
         fbyte=fbyte-1
      END DO

      fbyte = fbyte + 1 

      RETURN
      END
