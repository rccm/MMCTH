       SUBROUTINE write_db(dims,QAflags,DB550_aot,DB550_aotbe,DB_aex,DB_aot,DB_ssa,DB_sr,DB550_sdaot,
     1                    DB_ref,DB_naot,DB_cldfrac,DB_alg,DB_confflg,DTDB_aot,DTDB_QA, DTDB_flag,DB_uncert)
C-----------------------------------------------------------------------
C !F77
C
C !DESCRIPTION:
C
C    write_db writes Deep Blue data to new SDSs in MOD04 product.
C
C !INPUT PARAMETERS:
C
C    Type       Name             Description
C    ====       ====             ===========
C    INTEGER*4  dims             array of array dimensions
C    BYTE       QAflags          3-D array of bit flags
C    INTEGER*2  DB550_aot        2-D array of scaled integer values for
C                                 550nm tau.
C    INTEGER*2  DB_aex           2-D array of scaled integer values for
C                                 alpha.
C    INTEGER*2  DB_aot           3-D array of scaled integer values for
C                                 tau at 412, 470, and 650nm.
C    INTEGER*2  DB_ssa           3-D array of scaled integer values for
C                                 omega at 412, 470, and 650nm.
C    INTEGER*2  DB_sr            3-D array of scaled integer values for
C                                 surface reflectivity at at 412, 470, 
C                                 and 650nm.
C    INTEGER*2  DB550_sdaot      2-D array of scaled integer values for
C                                 standard deviation of 550nm tau.
C    INTEGER*2  DB_ref           3-D array of scaled integer values for
C                                 reflectivity at at 412, 470, and 650nm.
C    INTEGER*2  DB_naot          2-D array of integer values for number of						                                
C                                 pixels used in retrieval of AOT at 412,
C                                 470, and 650nm.
C    INTEGER*2  DB_cldfrac       2-D array of scaled integer values for 
C                                 cloud fraction
C    INTEGER*2  DB_alg           flag indicated algorithm used to retrieve
C                                 the majority of pixels per cell. 0=DB, 1=VEG, 2=MIXED
C    INTEGER*2  DB_confflg       Confidence flag for DB AOT retrievals
C
C !OUTPUT PARAMETERS:  none
C
C !REVISION HISTORY:
C
C    Initial Version by Jeremy Warner   12/01/2006
C    Clare Salustro added STD,          05/02/2008
C        Number_Pixels_Used, Reflectance_Land
C
C !TEAM-UNIQUE HEADER:
C
C    This software is developed by the Deep Blue Science Team
C    for the National Aeronautics and Space Administration,
C    Goddard Space Flight Center, under contract NAS5-02041.
C
C !REFERENCES AND CREDITS
C
C !DESIGN NOTES:
C
C   Externals:
C
C     MODIS_W_GENERIC            (MODIS_39500.f)
C
C   Functions:
C
C !END
C-----------------------------------------------------------------------
      IMPLICIT  NONE
      SAVE

c rhucek 01/08/98:  replaced include file mapic.inc with mapi.inc
c     Include 'mapic.inc'
      Include 'mapi.inc'
      Include 'mapi_hdfeos.inc'
      include 'PGS_MODIS_39500.f'
      include 'hdf.inc'
c     Include 'mod04.inc'

      INTEGER Ncell,Nrow

      CHARACTER*512 usrlog
      CHARACTER*256 fname, funcname, msg
      parameter (funcname="write_db")
      CHARACTER*26 istr
      INTEGER      prtn, file_version, pgs_pc_getreference
      INTEGER      LOGFLAG
      INTEGER      rtn, string_loc, fbyte1, lbyte1, dims(3)
      INTEGER      No_byte,Fmax,Lmax,Num_WL
      INTEGER      i,j,k,l,m,n,p, START2(3), DIM2_BUF(2), DIM3_BUF(3)
      PARAMETER    (Num_WL=3,No_byte=6,Fmax=150,Lmax=550)
      BYTE         QAflags(No_byte,Fmax,Lmax)
      INTEGER*2    DB_aot(Fmax,Lmax,Num_WL), DB_aex(Fmax,Lmax), DB550_aot(Fmax,Lmax)
      INTEGER*2    DB550_aotbe(Fmax,Lmax), DB_cldfrac(Fmax,Lmax)
      INTEGER*2    DB_ssa(Fmax,Lmax,Num_WL), DB_sr(Fmax,Lmax,Num_WL)
      INTEGER*2    DB550_sdaot(Fmax,Lmax)
      INTEGER*2    DB_ref(Fmax,Lmax,Num_WL), DB_naot(Fmax,Lmax)
      INTEGER*2    DB_alg(Fmax,Lmax), DB_uncert(Fmax,Lmax)
      INTEGER*4    DB_confflg(Fmax,Lmax)
      INTEGER*2    DTDB_aot(Fmax,Lmax)
      INTEGER*2    DTDB_flag(Fmax,Lmax)
      INTEGER*4    DTDB_QA(Fmax,Lmax) 

      INTEGER      modfil(MODFILLEN),ierr
C
C Buffer arrays
C
      INTEGER*2 AOT_BUF(Num_WL*Fmax*Lmax), AE_BUF(Fmax*Lmax), AOT550_BUF(Fmax*Lmax)
      INTEGER*2 AOT550BE_BUF(Fmax*Lmax), CLDFRAC_BUF(Fmax*Lmax)
      INTEGER*2 SSA_BUF(Num_WL*Fmax*Lmax), SR_BUF(Num_WL*Fmax*Lmax)
      INTEGER*2 SDAOT_BUF(Num_WL*Fmax*Lmax), SD550_BUF(Fmax*Lmax)
      INTEGER*2 REF_BUF(Num_WL*Fmax*Lmax), NAOT_BUF(Fmax*Lmax)
      INTEGER*2 ALG_BUF(Fmax*Lmax), CONF_BUF(Fmax*Lmax)
      INTEGER*2 DTDBFLG_BUF(Fmax*Lmax), DTDB_BUF(Fmax*Lmax)
      INTEGER*2 DTDBQA_BUF(Fmax*Lmax)
      INTEGER*2 UNCERT_BUF(Fmax*Lmax)
      BYTE      BBUF(No_byte*Fmax*Lmax)
      INTEGER*2 TMP_BUF(2*Fmax*Lmax)
      Ncell = dims(1)
      Nrow = dims(2)

      START2(1) = 0
      START2(2) = 0
      START2(3) = 0

C Open file

*/  Retrieve the filename of the output swath file.
      file_version=1
      prtn=pgs_pc_getreference(LUN_sw,file_version,fname)

      WRITE (istr, '(I25)') LUN_sw
      rtn = string_loc(istr, fbyte1, lbyte1)
      usrlog = "Retrieving filename for mod04, LUN "
     +          //istr(fbyte1:lbyte1)// " - pgs_pc_getreference"
      CALL ckstatus_s(prtn,usrlog,funcname,LOGFLAG)

      IF (OPMFIL(fname, 'a', MODFIL) .ne. MAPIOK) THEN
        msg = "Error opening MOD04 file "//fname
        call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,msg,funcname)
      ENDIF

      msg = "Opened MOD04 file "//fname
      call MODIS_SMF_SETDYNAMICMSG(MODIS_S_GENERIC, msg, funcname)

C Check for valid file and band numbers
      IF (modfil(1).le.0.OR.modfil(3).ne.DFACC_RDWR) THEN
        call MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,
     & 'Invalid SD_ID or invalid file access type',funcname)
      END IF

c Store data in buffers
      DIM2_BUF(1) = Ncell
      DIM2_BUF(2) = Nrow

      k = 0
      l = 0
      DO 10 i=1,Nrow
        DO 15 j=1,Ncell
          k = k + 1
          AE_BUF(k) = DB_aex(j,i)
          AOT550_BUF(k) = DB550_aot(j,i)
          AOT550BE_BUF(k) = DB550_aotbe(j,i)
          NAOT_BUF(k) = DB_naot(j,i)
          SD550_BUF(k) = DB550_sdaot(j,i)
          CLDFRAC_BUF(k) = DB_cldfrac(j,i)
          ALG_BUF(k)  = DB_alg(j,i)
          CONF_BUF(k) = DB_confflg(j,i)
          DTDB_BUF(k) = DTDB_aot(j,i)
          DTDBQA_BUF(k) = DTDB_QA(j,i)
          DTDBFLG_BUF(k) = DTDB_flag(j,i)
          UNCERT_BUF(k) = DB_uncert(j,i)
          DO 20 m=1,No_byte
             l = l + 1
             BBUF(l) = QAflags(m,j,i)
  20      CONTINUE  
  15    CONTINUE
  10  CONTINUE

      p = 0
      DO 25 n=1,Num_WL
        DO 25 i=1,Nrow
          DO 25 j=1,Ncell
             p = p + 1
             AOT_BUF(p) = DB_aot(j,i,n)
             SSA_BUF(p) = DB_ssa(j,i,n)
             SR_BUF(p) = DB_sr(j,i,n)
             REF_BUF(p) = DB_ref(j,i,n)

  25  CONTINUE

c Write data to file

        RTN = PMAR(MODFIL,'Deep_Blue_Angstrom_Exponent_Land',
     &             ' ',START2,DIM2_BUF,AE_BUF) 
        if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG 
     &      (MODIS_W_GENERIC,'PMAR for Deep_Blue_Angstrom_Exponent_Land','write_db')

        RTN = PMAR(MODFIL,'Deep_Blue_Aerosol_Optical_Depth_550_Land',
     &             ' ',START2,DIM2_BUF,AOT550_BUF) 
        if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG 
     &      (MODIS_W_GENERIC,'PMAR for Deep_Blue_Aerosol_Optical_Depth_550_Land','write_db')
     
        RTN = PMAR(MODFIL,'Deep_Blue_Aerosol_Optical_Depth_550_Land_Best_Estimate',
     &             ' ',START2,DIM2_BUF,AOT550BE_BUF) 
        if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG 
     &      (MODIS_W_GENERIC,'PMAR for Deep_Blue_Aerosol_Optical_Depth_550_Land_Best_Estimate','write_db')

        RTN = PMAR(MODFIL,'Deep_Blue_Aerosol_Optical_Depth_550_Land_STD',
     &             ' ',START2,DIM2_BUF,SD550_BUF)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Deep_Blue_Aerosol_Optical_Depth_550_Land_STD','write_db')
      
        RTN = PMAR(MODFIL,'Deep_Blue_Algorithm_Flag_Land',
     &             ' ',START2,DIM2_BUF,ALG_BUF) 
        if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG 
     &      (MODIS_W_GENERIC,'PMAR for Deep_Blue_Algorithm_Flag','write_db')

        RTN = PMAR(MODFIL,'Deep_Blue_Aerosol_Optical_Depth_550_Land_QA_Flag',
     &             ' ',START2,DIM2_BUF,CONF_BUF) 
        if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG 
     &      (MODIS_W_GENERIC,'PMAR for Deep_Blue_Aerosol_Optical_Depth_550_Land','write_db')
        
        RTN = PMAR(MODFIL,'Deep_Blue_Aerosol_Optical_Depth_550_Land_Estimated_Uncertainty',
     &             ' ',START2,DIM2_BUF,UNCERT_BUF) 
        if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG 
     &      (MODIS_W_GENERIC,'PMAR for Deep_Blue_Aerosol_Optical_Depth_550_Land_Estimated_Uncertainty','write_db')
     
        RTN = PMAR(MODFIL,'Deep_Blue_Cloud_Fraction_Land',
     &             ' ',START2,DIM2_BUF,CLDFRAC_BUF) 
        if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG 
     &      (MODIS_W_GENERIC,'PMAR for Deep_Blue_Cloud_Fraction_Land','write_db')
        
c        RTN = PMAR(MODFIL,'Dark_Target_Deep_Blue_Optical_Depth_550_Land_And_Ocean',
        RTN = PMAR(MODFIL,'AOD_550_Dark_Target_Deep_Blue_Combined',
     &             ' ',START2,DIM2_BUF,DTDB_BUF)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'AOD_550_Dark_Target_Deep_Blue_Combined','write_db')

c        RTN = PMAR(MODFIL,'Dark_Target_Deep_Blue_Optical_Depth_550_Land_And_Ocean_QA',
        RTN = PMAR(MODFIL,'AOD_550_Dark_Target_Deep_Blue_Combined_QA_Flag',
     &             ' ',START2,DIM2_BUF,DTDBQA_BUF)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Dark_Target_Deep_Blue_Optical_Depth_550_Land_And_Ocean_QA','write_db')
     
c        RTN = PMAR(MODFIL,'Dark_Target_Deep_Blue_Optical_Depth_550_Land_And_Ocean_AlgFlag',
        RTN = PMAR(MODFIL,'AOD_550_Dark_Target_Deep_Blue_Combined_Algorithm_Flag',
     &             ' ',START2,DIM2_BUF,DTDBFLG_BUF)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Dark_Target_Deep_Blue_Optical_Depth_550_Land_And_Ocean_AlgFlag','write_db')
        
        RTN = PMAR(MODFIL,'Deep_Blue_Number_Pixels_Used_550_Land',
     &             ' ',START2,DIM2_BUF,NAOT_BUF)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Deep_Number_Pixels_Used_550_Land','write_db')


      DIM3_BUF(1) = Ncell
      DIM3_BUF(2) = Nrow
      DIM3_BUF(3) = Num_WL

        RTN = PMAR(MODFIL,'Deep_Blue_Spectral_Aerosol_Optical_Depth_Land',
     &             ' ',START2,DIM3_BUF,AOT_BUF)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Deep_Blue_Aerosol_Optical_Depth_Land','write_db')

c        RTN = PMAR(MODFIL,'Deep_Blue_Aerosol_Optical_Depth_Land_STD',
c     &             ' ',START2,DIM3_BUF,SDAOT_BUF)
c         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
c     &      (MODIS_W_GENERIC,'PMAR for Deep_Blue_Aerosol_Optical_Depth_Land_STD','write_db')
     
        RTN = PMAR(MODFIL,'Deep_Blue_Spectral_Single_Scattering_Albedo_Land',
     &             ' ',START2,DIM3_BUF,SSA_BUF)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Deep_Blue_Single_Scattering_Albedo_Land','write_db')

        RTN = PMAR(MODFIL,'Deep_Blue_Spectral_Surface_Reflectance_Land',
     &             ' ',START2,DIM3_BUF,SR_BUF)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Deep_Blue_Surface_Reflectance_Land','write_db')
     
        RTN = PMAR(MODFIL,'Deep_Blue_Spectral_TOA_Reflectance_Land',
     &             ' ',START2,DIM3_BUF,REF_BUF)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Deep_Blue_Reflectance_Land','write_db')

      DIM3_BUF(1) = No_byte
      DIM3_BUF(2) = Ncell
      DIM3_BUF(3) = Nrow

        RTN = PMAR(MODFIL,'Quality_Assurance_Land',
     &             ' ',START2,DIM3_BUF,BBUF)
         if (RTN.NE.MAPIOK) CALL MODIS_SMF_SETDYNAMICMSG
     &      (MODIS_W_GENERIC,'PMAR for Quality_Assurance_Land','write_db')

c Close file and go home

      ierr = CLMFIL(MODFIL)

      RETURN
      END
