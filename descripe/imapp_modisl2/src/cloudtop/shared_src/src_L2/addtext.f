      integer function addtext(newtext,oldtext)

      implicit none
      include 'PGS_MODIS_39500.f'

C-----------------------------------------------------------------------
C !F77
C
C !Description:  This routine simply appends extra text to a text buffer.
C
C !Input Parameters:
C
C     character*(*)      newtext        Text to be appended
C     character*(*)      oldtext        Buffer to be appended to
C
C !Output Parameters: None
C
C !Revision History:
C
c Revision 1.1  1997/11/05  23:35:43  rhucek
c Initial revision
c
C
C !Team-Unique Header:
C
C   This software was developed by the MODIS Science Data Support Team
C   (SDST) for the National Aeronautics and Space Administration,
C   Goddard Space Flight Center, under contract NAS5-32373.
C
C !References and Credits:
C
C   Developer: JC Guu 06/23/97 jguu@ltpmail.gsfc.nasa.gov
C
C !Externals:
C
C   Function:
C      strlen             strlen.f
C
C !Internals:
C
C      integer       buflen    Nominal length of the buffer
C      integer       sln       Actual length of the buffer
C      integer       sl        Actual length of the text to be appended
C
C !END
C-----------------------------------------------------------------------

C-----Declaration of PARAMETERs 
      character*1 BLANK
      PARAMETER (BLANK = ' ')

      character*(*) FUNCNAME
      parameter (FUNCNAME = 'addtext')

      integer FAIL, SUCCEED
      PARAMETER (FAIL=-1, SUCCEED=0)

C-----Function argument declarations
      character*(*) newtext,oldtext

C-----Local variable and array declarations
      character*512 msgbuf

      integer strlen
      integer sl_old,sl_new,buflen

C-----------------------------------------------------------------------
C If "newtext" is blank, return without update.  Otherwise attempt
C to update "oldtext".  If "oldtext" buffer length is too short, 
C report this with LogStatus message and return without update.
C-----------------------------------------------------------------------
 
      If (newtext .EQ. BLANK) Then
         addtext = SUCCEED

      Else
         buflen = len(oldtext)
         sl_old = strlen(oldtext)
         sl_new = strlen(newtext)

         If (buflen .LT. (sl_new+sl_old)) Then
            addtext = FAIL

            msgbuf = 
     1         'Old text buffer too small to contain '
     2         // 'sum of old and new '
     3         // char(10) 
     4         // 'text strings.  '
     5         // 'Old text buffer not updated by AddText'

            Call modis_smf_setdynamicmsg(MODIS_W_GENERIC,msgbuf,
     2           FUNCNAME)
         Else 
            oldtext(sl_old+1:sl_new+sl_old) = newtext(1:sl_new)
            addtext = SUCCEED
         EndIf
      EndIf

      Return
      End
