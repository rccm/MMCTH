      SUBROUTINE MESSAGE( ROUTINE, SCRIPT, CODE, LEVEL )

C-----------------------------------------------------------------------
C !F77
C 
C !DESCRIPTION:
C     Subroutine for writing error messages to the PGS Toolkit log file.
C
C !INPUT PARAMETERS:
C     ROUTINE     Name of the routine where message originated
C     SCRIPT      Message string
C     CODE        Error code number
C     LEVEL       Error level number
C                 0 = Warning
C                 1 = Recoverable error
C                 2 = Fatal error (execution will stop)
c                 3 = Notice - General information you wish
c                     to have recorded in the log files
C
C !OUTPUT PARAMETERS:
C     None
C
C !REVISION HISTORY:
C
C !TEAM-UNIQUE HEADER:
C     Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
C
C !DESIGN NOTES:
C MESSAGE should be used directly in routines where errors occur.
C Routines should also include the error code as the last argument in
C the call arguments, so that the calling program can attempt to recover
C from the error and continue. The error level 'Fatal error' should be
C used sparingly, as it will cause execution to stop in ERRMSG. The
C calling program should always be given the chance to recover from
C any errors that occur, unless the error is truly fatal (e.g. an
C infinite loop will result).
C
C !END
C-----------------------------------------------------------------------

      implicit none

!      include 'PGS_MODIS_39500.f'

c ... Arguments
     
      character*(*) routine, script
      integer code, level

c ... Local variables

      character errorlevel*30, errorcode*30, outstring*255
      logical remove_all
      integer length
      
c ... Set error level string

      errorlevel = 'Warning'
      if( level .eq. 1 ) errorlevel = 'Recoverable error'
      if( level .eq. 2 ) errorlevel = 'Fatal Error'
      if( level .eq. 3 ) errorlevel = 'Notice'

c ... Set error code string

      write( errorcode, '( ''Error code = '', i6 )' ) code

c ... Construct output message string

      if( level .eq. 3 ) then
        outstring = errorlevel // char( 10 ) //
     &              script
      else
        outstring = errorlevel // char( 10 ) //
     &              script // char( 10 ) //
     &              errorcode
      endif
            
c ... Compress output message string

      remove_all = .false.
      outstring = trim(outstring)
	  length = len(outstring)
      
c ... Write message to log file
 
      if ( level .eq. 0 .or. level .eq. 1 ) then

c ...   Error is non-fatal

	    print*, outstring(1:length), routine

!        call modis_smf_setdynamicmsg( MODIS_E_GENERIC, 
!     &    outstring( 1 : length ), routine )

      else if ( level .eq. 2 ) then

c ...   Error was fatal

		print*, "STOPPING BECAUSE OF FATAL ERROR: ", outstring(1:length), routine
!        call modis_smf_setdynamicmsg( MODIS_E_GENERIC,
!     &    'STOPPING BECAUSE OF FATAL ERROR: ' //
!     &    outstring( 1 : length ), routine )
        call exit(1)

      else 

c ...   Message is notice only

!        call modis_smf_setdynamicmsg( MODIS_A_GENERIC,
!     &    outstring( 1 : length ), routine )
 
	   print*, outstring(1:length), routine
		     
      end if

      END
