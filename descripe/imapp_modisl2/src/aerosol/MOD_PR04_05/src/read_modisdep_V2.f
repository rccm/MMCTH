      SUBROUTINE READ_AVERAGE_MODIS_DEPEN(MODFIL_MOD05,MODFIL_MOD07,
     *           ISCAN,MODFIL_L1B_1KM)

C-----------------------------------------------------------------------
C!F77
C
C!DESCRIPTION:
C                        **This subroutine reads the ancillary data
C                         MOD05 and MOD07 In HDF format
C
C
C!INPUT PARAMETERS:
C         MODFIL_*       value for reading the filemod05,mod07
C         ISCAN         scan number in a scan cube 1- 100 typically
C!OUTPUT PARAMETERS:
C
c
C
C!REVISION HISTORY:
C $Log: read_modisdep_V2.f,v $
c 10/19/1999 fhliang
c fixed prolog;
C
C!TEAM-UNIQUE HEADER: - NOT YET DEFINED
C
C!REFERENCES AND CREDITS
C     WRITTEN BY: Shana Mattoo
C
C!DESIGN NOTES:
C
C!End
C-----------------------------------------------------------------------

      IMPLICIT NONE
      SAVE

      Include 'PGS_MODIS_39500.f'
      Include 'mapi.inc'
      Include 'mod04.inc'

      CHARACTER*80 DATA_TYPE,ARRNM,GRPNM
      INTEGER DIM_SIZE(3),RANK,NUMSCAN
      INTEGER RESULT, MODFIL_MOD05(MODFILLEN), START( 3 ),RTN
      INTEGER  EDGE_water( 3 ),EDGE_ozone( 3 )
      INTEGER MODFIL_L1B_1KM(MODFILLEN)
      INTEGER MODFIL_MOD07(MODFILLEN),ISCAN
      INTEGER * 2 WATER(ISWATH*ILINE)
      INTEGER * 2 TOTOZONE((ISWATH/5)*(ILINE/5))
      REAL  NEWWATER_VAPOR(ISWATH,ILINE),ADD_OFFSET,SCALE_FACTOR
      REAL  NEWTOT_OZONE((ISWATH/5),(ILINE/5))
      REAL  WATER_VAPOR(NUMCELLS),fillvalue_water,fillvalue_ozone
      REAL  OZONE_COL(NUMCELLS)
      INTEGER IK,IS,L
      EXTERNAL GMAR,GMARDM

C Call GMFIN to get 'Number of Scans' from L1B dataset

       arrnm= 'Number of Scans'
       data_type = 'INTEGER*4'
       Rank = 1
       RTN = GMFIN(Modfil_L1B_1km,arrnm,data_type,Rank,NUMSCAN)
       IF (rtn.ne.0) CALL MODIS_SMF_SETDYNAMICMSG(MODIS_W_GENERIC,
     &   'MAPI function GMFIN failed',
     &'READ_AVERAGE_MODIS_DEPEN Get_Metatable for number of scans')


            grpnm=' '
            arrnm='Water_Vapor_Near_Infrared'

c retrieve the rank, dimensions and data type of sds data.

        rank=2
        add_offset=0.0
        scale_factor=1000.0
        fillvalue_water=-99

C RETRIEVE THE RANK, DIMENSIONS AND DATA TYPE OF SDS DATA.


       If (GMARDM(MODFIL_MOD05,arrnm,grpnm,data_type,Rank,
     &  Dim_Size).ne.MAPIOK)
     &   CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,'GMARDM failed',
     &   'mod05 water ')

         START( 1 ) = 0
         START( 2 ) = (ISCAN-1)*10
         EDGE_water(1)=DIM_SIZE(1)
c             ******* 10 lines and 200  in one scan cube
         EDGE_water(2)=DIM_SIZE(2)/NUMSCAN

       RESULT = GMAR(MODFIL_MOD05,arrnm,grpnm,start,edge_water,WATER)


         IF( RESULT .NE. MAPIOK)
     &         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,
     &       'GMAR  failed for water','READ_MODISDEP')

c                    ******** save data in new array newwater
           L=0
          DO IK=1,EDGE_water(2)
          DO IS=1,EDGE_water(1)
           L=L+1
          NEWWATER_VAPOR(IS,IK)=(REAL(WATER(L))-ADD_OFFSET)
     *    /SCALE_FACTOR
          ENDDO
          ENDDO

c
C          READ Total Column Ozone DOBSON UNITS

                grpnm=' '
            arrnm='Total_Ozone'

c retrieve the rank, dimensions and data type of sds data.
c
        rank=2
        add_offset=0.0
        scale_factor=10.0
       fillvalue_ozone=-10000
C
C RETRIEVE THE RANK, DIMENSIONS AND DATA TYPE OF SDS DATA.
         If (GMARDM(MODFIL_MOD07,arrnm,grpnm,data_type,Rank,Dim_Size)
     &                    .ne.MAPIOK)
     &   CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,
     &'GMARDM failed dim_info failed FOR TOTAL OZONE','READ_MODISDEP')

C
         START( 1 ) = 0
         START( 2 ) = (ISCAN-1)*2
         EDGE_Ozone(1)=DIM_SIZE(1)
c             ******* 2 lines and 200  in one scan cube
         EDGE_Ozone(2)=DIM_SIZE(2)/NUMSCAN
c
      RESULT = GMAR(MODFIL_MOD07,arrnm,grpnm,start,edge_ozone,TOTOZONE)
        IF( RESULT .NE. MAPIOK)
     &         CALL MODIS_SMF_SETDYNAMICMSG(MODIS_F_GENERIC,
     &       'GMAR  failed for TOTOZONE','READ_MODISDEP')
c                    ******** save data in new array newwater
           L=0
          DO IK=1,EDGE_Ozone(2)
          DO IS=1,EDGE_Ozone(1)
           L=L+1
           NEWTOT_OZONE(IS,IK)=(REAL(TOTOZONE(L))-ADD_OFFSET)
     *    /SCALE_FACTOR
          ENDDO
          ENDDO
1115       continue

c            DO Is= 1,EDGE_Ozone(1)
c         WRITE(36,*)ISCAN,is,(NEWTOT_OZONE(is,ik ),ik=1,EDGE_Ozone(2))
c           ENDDO
         CALL  AVERAGE_MODIS(EDGE_water,EDGE_ozone,
     *         NEWWATER_VAPOR,NEWTOT_OZONE,WATER_VAPOR,OZONE_COL,
     *         fillvalue_water,fillvalue_ozone)

       RETURN
       END


C***********************************************************************

       SUBROUTINE AVERAGE_MODIS(EDGE_water,EDGE_ozone,
     * NEWWATER_VAPOR,NEWTOT_OZONE,WATER_VAPOR,OZONE_COL,
     * fillvalue_water,fillvalue_ozone)

C-----------------------------------------------------------------------
C!F77
C
C!DESCRIPTION:
C                        **This subroutine averages on 10*10 pixel
C                          mod05 and mod07  data
C
C
C!INPUT PARAMETERS:
C         ISWATH        number of pixels across the swath
C         ILINE         number of lines
C         NUMCELLS      maximum dimension for 10*10 sqaures in 1354*10
c        NEWWATER_VAPOR Total_Ozone in "Dobson"
c        NEWTOT_OZONE   Column_Water_Vapor in "cm"
C!OUTPUT PARAMETERS:
C         WATER_VAPOR    Total_Ozone in "Dobson" averaged for 10*10 pixels
C                         from 1*1 km resolution
c         OZONE_COL      Column_Water_Vapor in "cm"averageg on 10 *10 pixels
C                        from 5*5  pixel resolution
C
C!REVISION HISTORY:
C $Log: read_modisdep_V2.f,v $
c 10/19/1999 fhliang
c fixed prolog;
C
C!TEAM-UNIQUE HEADER: - NOT YET DEFINED
C
C!REFERENCES AND CREDITS
C     WRITTEN BY: Shana Mattoo
C
C!DESIGN NOTES:
C
C!End
C-----------------------------------------------------------------------

      IMPLICIT NONE
      SAVE

      Include 'mod04.inc'

      INTEGER IA,ICOL,IDIMX,IE,IFACTOR,IGRIDXTOT
      INTEGER EDGE_water(3),EDGE_ozone(3)
      INTEGER NUMSQ
      INTEGER IGRIDYTOT,IGX,IP,IROW,IS,ISQ,JDIMY,JE,JGY,JP,JS
      REAL WATER_VAPOR1,WATER_VAPOR(NUMCELLS)
      REAL OZONE_COL1,OZONE_COL(NUMCELLS)
      REAL NEWWATER_VAPOR(ISWATH,ILINE),fillvalue_ozone,fillvalue_water
      REAL NEWTOT_OZONE((ISWATH/5),(ILINE/5))

*
C ********AVERAGEWATER_VAPOR  FOR EACH 10*10 SQUARE

      IGRIDXTOT=EDGE_water(1)
      IGRIDYTOT=ILINE
      NUMSQ=(IGRIDXTOT/10)*(IGRIDYTOT/10)
      IFACTOR=1
      IDIMX=10*IFACTOR
      JDIMY=10*IFACTOR
      ISQ=0
      JS=1
      JE=JDIMY
      DO 200 IROW=1,IGRIDYTOT*IFACTOR,JDIMY
        IS=1
        IE=IDIMX
        DO 210 ICOL=1,IGRIDXTOT*IFACTOR,IDIMX
          IGX=0
          JGY=0
          ISQ=ISQ+1
            WATER_VAPOR1=0.0
C           Check for _FillValues
          DO 220 IP=IS,IE
            IGX=IGX+1
            DO 230 JP=JS,JE
            IF(NEWWATER_VAPOR(IP,JP).GT.0.0) THEN
              JGY=JGY+1
              WATER_VAPOR1=WATER_VAPOR1+NEWWATER_VAPOR(IP,JP)
            ENDIF
  230       CONTINUE
  220     CONTINUE
          IF(JGY .GE.1) THEN
            WATER_VAPOR(ISQ)= WATER_VAPOR1/JGY
          ELSE
            WATER_VAPOR(ISQ)=fillvalue_water
          ENDIF
          IS=IS+IDIMX
          IE=IE+IDIMX
  210   CONTINUE
        JS=JS+JDIMY
        JE=JE+JDIMY
  200 CONTINUE
  250 CONTINUE
C            ********AVERAGE OZONE_COLUM FOR EACH 10*10 SQUARE
       IGRIDXTOT=EDGE_ozone(1)
       IGRIDYTOT=ILINE/5
       NUMSQ=IGRIDXTOT/IGRIDYTOT
      IFACTOR=1
      IDIMX=2
      JDIMY=2
      ISQ=0
      JS=1
      JE=JDIMY
      DO 300 IROW=1,IGRIDYTOT*IFACTOR,JDIMY
        IS=1
        IE=IDIMX
        DO 310 ICOL=1,IGRIDXTOT*IFACTOR,IDIMX
          IGX=0
          JGY=0
          ISQ=ISQ+1
          OZONE_COL1=0.0
          DO 320 IP=IS,IE
            IGX=IGX+1
            DO 330 JP=JS,JE
          IF(NEWTOT_OZONE(IP,JP).GT.0.0) THEN
             JGY=JGY+1
             OZONE_COL1= OZONE_COL1+NEWTOT_OZONE(IP,JP)
          ENDIF
  330       CONTINUE
  320     CONTINUE
          IF(JGY .GE.1) THEN
            OZONE_COL(ISQ)= OZONE_COL1/JGY
          ELSE
           OZONE_COL(ISQ)=fillvalue_ozone
          ENDIF
          IS=IS+IDIMX
          IE=IE+IDIMX
  310   CONTINUE
        JS=JS+JDIMY
        JE=JE+JDIMY
  300 CONTINUE
  350 CONTINUE
c            DO Ia = 1,NUMSQ
c          WRITE(55,*)NUMSQ,Ia,WATER_VAPOR(Ia), OZONE_COL(Ia)
c         ENDDO

       RETURN
      END
