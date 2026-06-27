c     Subroutine Copy_cldmask_M(Modfil_cldmsk,cldmsk_lrn,MCF_cmlrn,
      Subroutine Copy_cldmask_M(Modfil_cldmsk,cldmsk_lrn,
     +                          pcn4s,pcnnd,
     +                          pcnnt,pcnng,pcnni,pcnnl,
     +                          pcnnw,pcnns,pcnnv,pcnnr,
     +                          pcnnc,max_sol,min_sol)

      implicit none
      save

      include 'mod06uw_debug.inc'
c 02/17/98 fhliang: changed 'mapic.inc' to 'mapi.inc' and 'hdf.inc'
c     include 'mapic.inc'
      include 'mapi.inc'
      include 'hdf.inc'
      include 'PGS_PC.f'
      include 'PGS_MODIS_39500.f'
c 02/18/98 fhliang: added following include file:
      include 'PGS_SMF.f'


c-----------------------------------------------------------------------
c!F77
c
c!Purpose:
c   To extract the Cloud Mask Product Specific Attributes from the
c   the cloud mask product HDF file that are needed by
c   the atmospheric profiles code.
c
c!Description:
c   The program uses information found in the cloud mask MCF
c   file and extracts using PGS calls.  
c
c!Input Parameters:
c  integer  Modfil_cldmsk   Array containing HDF SD (index 1) and VS 
c                           (index 2) file identifiers and 
c                           access mode (index 3) for Cloud Mask file.
c  integer  cldmsk_lrn      Variable containing the PCF logical 
c                           reference number (LRN) to the cloud mask
c                           file.
cc  integer  MCF_cmlrn       Variable containing the PCF logical 
cc                           reference number (LRN) to the cloud 
cc                           mask metadata configuration file (MCF)
c!Output Parameters:
c
c pcn4s         Percentage of pixels in category 4 per granule
c pcnnd         Percentage of day processed pixels per granule
c pcnnt         Percentage of night processed pixels per granule
c pcnng         Percentage of sunglint found pixles per granule
c pcnni         Percentage of snow/ice processed pixels per granule
c pcnnl         Percentage of land processed pixels per granule
c pcnnw         Percentage of water processed pixels per granule
c pcnns         Percentage of shadow pixels found per granule
c pcnnv         Percentage of thin cirrus (vis) found pixels per granule
c pcnnr         Percentage of thin cirrus (ir) found pixels per granule
c pcnnc         Percentage of non-cloud obstruction pixels per granule
c max_sol       Maximum solar zenith angle for this granule
c min_sol       Minimum solar zenith angle for this granule
c
c !REVISION HISTORY:
c $Id: Copy_cldmask_M.f,v 1.6 1999/04/16 22:43:44 kis Exp $
c
c!Team Unique Header:
c HDF portions developed at the National Center for Supercomputing
c Applications at the University of Illinois at Urbana-Champaign.
c
c !REFERENCES AND CREDITS:
c
c!Design Notes:
c
c  Externals:
c
c    Functions and Subroutines:
c        pgs_met_init          (libPGSTK.a)
c        pgs_met_setattr_*     (libPGSTK.a)
c        message  writes error message to LogStatus file
c
c    Named Constant:
c        DFACC_READ            (hdf.inc)
c        MAPIOK                (mapi.inc)
c        MCORE_*               (mapi.inc)
c        MECS_CORE             (mapi.inc)
c        MFAIL                 (mapi.inc)
c        MODFILLEN             (mapi.inc)
c        MODIS_E_GENERIC       (PGS_MODIS_39500.f)
c        MODIS_W_GENERIC       (PGS_MODIS_39500.f)
c        PGSd_MET_NUM_OF_GROUPS(PGS_MET.f:included in "mapi.inc")
c        PGSd_PC_FILE_PATH_MAX (PGS_PC.f)
c        PGS_S_SUCCESS         (PGS_SMF.f)
c        P_SDID                (mapi.inc)
c        P_ACCESS              (mapi.inc)
c
c!END
c----------------------------------------------------------------------

c     Scalar Arguments
c     integer cldmsk_lrn,MCF_cmlrn
      integer cldmsk_lrn
      real pcn4s,pcnnd,pcnnt,pcnng,pcnni,pcnnl,
     +     pcnnw,pcnns,pcnnv,pcnnr,pcnnc,max_sol,min_sol


      integer modfil_cldmsk(MODFILLEN)

c     parameters
      real Meta_misg 
      parameter (Meta_misg = -99999.0)

c     Local scalars
c     character*49  Master_Grp_Hndls(PGSd_MET_NUM_OF_GROUPS)
      integer    i,rtn,version

c     Local arrays
      character*36  attr,Field_Name(13)
      character*70 buf_char(13)

c ... set program name for error messaging
      character*32 routine_name
      parameter ( routine_name = 'mod06_close_files' )
      
c     external routines
      integer    pgs_met_init,pgs_met_getpcattr_s
      external   pgs_met_init,pgs_met_getpcattr_s,message

c ----------------------------------------------------------------------

C     Initialize holding arrays
      do 10 i = 1, 13
         buf_char(i) = '  '
   10 continue

c     do 11 i = 1 , PGSd_MET_NUM_OF_GROUPS
c        Master_Grp_Hndls(i) = '  '
c  11 continue

c ... Initialize output variables
      pcn4s = Meta_misg
      pcnnd = Meta_misg
      pcnnt = Meta_misg
      pcnng = Meta_misg
      pcnni = Meta_misg
      pcnnl = Meta_misg
      pcnnw = Meta_misg
      pcnns = Meta_misg
      pcnnv = Meta_misg
      pcnnr = Meta_misg
      pcnnc = Meta_misg
      max_sol = Meta_misg
      min_sol = Meta_misg

c     Check inputs
      if (modfil_cldmsk(P_SDID).le.0) then
        call message( routine_name,
     &  'Failed to get cldmask metadata.' //
     &  ' [OPERATOR ACTION: Contact SDST]', 0, 2 )

      else if (modfil_cldmsk(P_ACCESS).ne.DFACC_READ) then
         call message( routine_name,
     &  'Failed to get cldmask metadata.' //
     &  ' [OPERATOR ACTION: Contact SDST]', 0, 2 )
      endif

c 03/13/98 fhliang commented out following block.
cC     Initialize metadata tool defining Master_Grp_Hndls, and set 
cC     CoreMetadata.0 attribute fields that have values provided in MCF.
c      rtn = pgs_met_init(MCF_cmlrn, Master_Grp_Hndls)
c      if (rtn.ne.PGS_S_SUCCESS) then
c        call message( routine_name,
c     &  'Failed to initialize cldmask metadata file ',
c     &  0, 2 )
c      endif

c ... Initialize output variables
      

c----------------------------------------------------------------
c     Set CoreMetadata.0 attribute fields that have values retrieved 
c     from mask product.  Note that the use of M-API parameters are 
c     defined in "mapi.inc".
c----------------------------------------------------------------

      Field_Name(1) = 'PARAMETERVALUE.5'
      Field_Name(2) = 'PARAMETERVALUE.8'
      Field_Name(3) = 'PARAMETERVALUE.9'
      Field_Name(4) = 'PARAMETERVALUE.10'
      Field_Name(5) = 'PARAMETERVALUE.11'
      Field_Name(6) = 'PARAMETERVALUE.12'
      Field_Name(7) = 'PARAMETERVALUE.13'
      Field_Name(8) = 'PARAMETERVALUE.14'
      Field_Name(9) = 'PARAMETERVALUE.15'
      Field_Name(10) = 'PARAMETERVALUE.16'
      Field_Name(11) = 'PARAMETERVALUE.17'
      Field_Name(12) = 'PARAMETERVALUE.18'
      Field_Name(13) = 'PARAMETERVALUE.19'

      do 20 i = 1, 13
         attr = Field_Name(i)
         version = 1
         rtn = pgs_met_getpcattr_s(cldmsk_lrn,version,
     +                             'CoreMetadata.0',attr,buf_char(i))
         if (rtn.ne.PGS_S_SUCCESS) then
           call message( routine_name,
     +      'Failed to extract PSA from Cloud Mask File.' //
     +      ' [OPERATOR ACTION: Contact SDST]', 0, 1 )
            buf_char(i) = '-9999900'
         endif
   20 continue

c ... Place extracted values into correct (real) variable names
      read(buf_char(1),'(f8.2)') pcn4s
      read(buf_char(2),'(f8.2)') pcnnd
      read(buf_char(3),'(f8.2)') pcnnt
      read(buf_char(4),'(f8.2)') pcnng
      read(buf_char(5),'(f8.2)') pcnni
      read(buf_char(6),'(f8.2)') pcnnl
      read(buf_char(7),'(f8.2)') pcnnw
      read(buf_char(8),'(f8.2)') pcnns
      read(buf_char(9),'(f8.2)') pcnnv
      read(buf_char(10),'(f8.2)') pcnnr
      read(buf_char(11),'(f8.2)') pcnnc
      read(buf_char(12),'(f8.2)') max_sol
      read(buf_char(13),'(f8.2)') min_sol

c----------------------------------------------------------------
c ....................................................................
      if (debug .gt. 0) then

      write (h_output, '(/,72(''-''),/)' )
      write (h_output, '(''Copy_cldmask_meta debug info'',/)' )
      write (h_output,fmt='(1x,''Metadata Cloud Mask stats: '')')
      write (h_output,fmt='(1x,''% with confidence <  1% '',f9.2)')pcn4s
      write (h_output,fmt='(1x,''percent NCO found '',f9.2)')pcnnc
      write (h_output,fmt='(1x,''percent shadow found '',f9.2)')pcnns
      write (h_output,fmt='(1x,''percent day processed '',f9.2)')pcnnd
      write (h_output,fmt='(1x,''percent night processed '',f9.2)')pcnnt
      write (h_output,fmt='(1x,''percent sunglint found '',f9.2)')pcnng
      write (h_output,fmt='(1x,''percent snow/ice processed '',f9.2)')
     +                           pcnni
      write (h_output,fmt='(1x,''percent land processed '',f9.2)')pcnnl
      write (h_output,fmt='(1x,''percent water processed '',f9.2)')pcnnw
      write (h_output,fmt='(1x,''percent thin cirrus solar found '',
     +                    f9.2)')pcnnr
      write (h_output,fmt='(1x,''percent thin cirrus ir found '',f9.2)')
     +                           pcnnv
      write (h_output,fmt='(1x,''max solar zenith angle '',f9.2)')
     +                           max_sol
      write (h_output,fmt='(1x,''min solar zenith angle '',f9.2)')
     +                           min_sol

      endif
c ....................................................................


      return
      end
