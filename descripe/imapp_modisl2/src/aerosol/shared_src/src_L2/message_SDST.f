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
C     20-SEP-1995 Liam Gumley, CIMSS/SSEC
C                 Created
C     07-MAR-1996 Liam Gumley, CIMSS/SSEC
C                 Changed messages to lowercase
C     06-JUN-1996 Liam Gumley, CIMSS/SSEC
C                 Modified to work with PGSTK
C     23-APR-1997 Kathy Strabala, CIMSS/SSEC
C                 Added status message option
C     27-OCT-1998 Liam Gumley, CIMSS/SSEC
C                 Simplified this routine
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

      include 'PGS_MODIS_39500.f'

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
      call strcompress( outstring, remove_all, length )
      
c ... Write message to log file
 
c ... Error is non-fatal
      if ( level.eq.0 .or. level.eq.1) then
         call modis_smf_setdynamicmsg( MODIS_E_GENERIC, 
     &        outstring( 1 : length ), routine )

c ... Report error "script" and stop because error is fatal
      else if ( level .eq. 2 ) then
         call modis_smf_setdynamicmsg( MODIS_E_GENERIC,
     &        outstring( 1 : length ), routine )

         call modis_smf_setdynamicmsg( MODIS_E_GENERIC,
     &        'Stopping in message.f because of fatal error ' //
     &        '[OPERATOR ACTION: Contact SDST]', 'message' )

         call exit(1)

c ... Message is notice only
      else if (level .eq. 3) then 
         call modis_smf_setdynamicmsg( MODIS_A_GENERIC,
     &        outstring( 1 : length ), routine )

c ...Unrecognized message level
      else 
         call modis_smf_setdynamicmsg( MODIS_E_GENERIC, 
     &        outstring( 1 : length ), routine )

         call modis_smf_setdynamicmsg( MODIS_E_GENERIC,
     &        'Stopping in message.f because message level is not one of 0,1,2, or 3. ' //
     &        '[OPERATOR ACTION: Contact SDST]', 'message' )

         call exit(1)
      end if


      END
