      INTEGER FUNCTION GET_PLATFORM_NAME(FILEID, VERSION, STATUS, NAME)
      
c-----------------------------------------------------------------------
c !F77
c
c !DESCRIPTION:
c     Given the process control file (PCF) logical unit number and
c     version number of a MODIS MOD021KM or MYD021KM product file,
c     return the name of the satellite platform (Terra or Aqua).
c
c !INPUT PARAMETERS:
c     FILEID     PCF logical unit number of a MOD021KM or MYD021KM file
c     VERSION    PCF version number of a MOD021KM or MYD021KM file
c           
c !OUTPUT PARAMETERS:
c     GET_PLATFORM_NAME  Return status
c                        (0=success, -1=failure)
c     STATUS             PGS toolkit error status
c                        (see Note 1)
c     NAME               Platform name
c                        ('Terra' or 'Aqua' if successful)
c                        (' ' otherwise)
c
c NOTES:
c     (1) For a description of STATUS values, see the SDP toolkit
c         documentation for 'PGS_MET_GetPCAttr'. The status values
c         are defined in the include file PGS_MET_13.f.
c
c !REVISION HISTORY:
c
c !TEAM-UNIQUE HEADER:
c     Developed by the MODIS Group, CIMSS/SSEC, UW-Madison.
c
c !END
c-----------------------------------------------------------------------

      IMPLICIT NONE

c ... Include files
      include 'PGS_SMF.f'
      
c ... Arguments
      integer fileid, version, status
      character*(*) name
      
c ... Local variables
      character*80  hdfattrname, parmname, parmvalue
      integer rtn

c ... External functions      
      integer pgs_met_getpcattr_s
      external pgs_met_getpcattr_s

c ... Get the metadata attribute which contains the platform name
c ... ('.1' must be appended to the parameter name because it
c ... is the first instance of the name in a group)
      hdfattrname = 'CoreMetadata.0'
      parmname = 'ASSOCIATEDPLATFORMSHORTNAME.1'
      parmvalue = ' '
      rtn = pgs_met_getpcattr_s(fileid, version, hdfattrname, parmname,
     &  parmvalue)

c ... Set values to return to caller
c ... (PGS_S_SUCCESS is defined in PGS_SMF.f)
      if (rtn .eq. PGS_S_SUCCESS) then

c ...   Success
        get_platform_name = 0      
        status = 0
        name = parmvalue
        
      else

c ...   Failure      
        get_platform_name = -1
        status = rtn
        name = ' '
        
      endif
                  
      END
