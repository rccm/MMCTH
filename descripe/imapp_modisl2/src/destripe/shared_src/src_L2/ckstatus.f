
      subroutine ckstatus_f(rtn,usrlog,funcname,runmode)

      implicit none

      include 'ckstatus.inc'
      include 'PGS_MODIS_39500.f'
      include 'hdf.inc'

C----------------------------------------------------------------------------
C !F77
C
C !Description:
C
C   This routine is used to check the return code of functions. If
C   a function only returns SUCCEED(0) or FAIL(-1), then this routine
C   should be used instead of ckstatus_s(). 
C
C !Input Parameters: 
C
C   integer         rtn        The return code form a function
C   character*(*)   usrlog     The message to be written to the LogStatus
C   character*(*)   funcname   The name of the function which returns.
C   integer         runmode    The flag for turning on or off the writing
C                              of SUCCESS messages to LogStatus 
C
C !Output Parameters: None
C
C !Revision History:
C
C
c Revision 1.4  1997/11/11  19:55:48  rhucek
c Removed reference to hdf.inc file in ckstatus_s.
c
c Revision 1.3  1997/08/26  15:57:45  rhucek
c removed functions int2str and colseq
c
c Revision 1.2  1997/06/18  13:45:07  vlin
c used to run V2 Metadata writer
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
C   Developer: JC Guu 03/10/97 jguu@ltpmail.gsfc.nasa.gov
C
C !Externals:
C
C   FAIL                        hdf.inc
C   MODIS_E_GENERIC             PGS_MODIS_39500.f
C   MODIS_S_GENERIC             PGS_MODIS_39500.f
C   DEBUG                       ckstatus.inc
C
C !Internals:
C
C !END
C----------------------------------------------------------------------------

      integer rtn,msgcode,runmode
      character*(*) funcname,usrlog

*---------------------------------------------------------------------------*


      if(rtn.eq.FAIL) then
          msgcode=MODIS_E_GENERIC
          call modis_smf_setdynamicmsg(msgcode,usrlog,funcname)
      elseif(runmode.eq.DEBUG) then
          msgcode=MODIS_S_GENERIC
          call modis_smf_setdynamicmsg(msgcode,usrlog,funcname)
      end if
      

      end



******/---------------------------------------------------------------/******

      subroutine ckstatus_s(rtn,usrlog,funcname,runmode)

      implicit none

      include 'PGS_SMF.f'
      include 'PGS_MODIS_39500.f'
      include 'ckstatus.inc'

C----------------------------------------------------------------------------
C !F77
C
C !Description:
C
C   This routine is used to check the return code of functions. If
C   a function returns PGS_S_SUCCESS for success, this routine should
C   be used instead of ckstatus_f().
C
C !Input Parameters:
C
C   integer         rtn        The return code form a function
C   character*(*)   usrlog     The message to be written to the LogStatus
C   character*(*)   funcname   The name of the function which returns.
C   integer         runmode    The flag for turning on or off the writing
C                              of SUCCESS messages to LogStatus
C
C !Output Parameters: None
C
C !Revision History:
C
C
c Revision 1.4  1997/11/11  19:55:48  rhucek
c Removed reference to hdf.inc file in ckstatus_s.
c
c Revision 1.3  1997/08/26  15:57:45  rhucek
c removed functions int2str and colseq
c
c Revision 1.2  1997/06/18  13:45:07  vlin
c used to run V2 Metadata writer
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
C   Developer: JC Guu 03/10/97 jguu@ltpmail.gsfc.nasa.gov
C
C !Externals:
C
C   MODIS_E_GENERIC             PGS_MODIS_39500.f
C   MODIS_S_GENERIC             PGS_MODIS_39500.f
C   DEBUG                       ckstatus.inc
C
C !Internals:
C
C !END
C----------------------------------------------------------------------------
      integer rtn,msgcode,runmode
      character*(*) funcname,usrlog

*---------------------------------------------------------------------------*


      if(rtn.eq.PGS_S_SUCCESS) then
          if(runmode.eq.DEBUG) then
              msgcode=MODIS_S_GENERIC
              call modis_smf_setdynamicmsg(msgcode,usrlog,funcname)
          end if
      else
          msgcode=MODIS_E_GENERIC
          call modis_smf_setdynamicmsg(msgcode,usrlog,funcname)
      end if


      end
