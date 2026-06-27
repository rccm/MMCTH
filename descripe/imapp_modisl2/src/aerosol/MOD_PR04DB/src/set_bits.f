      SUBROUTINE set_bits(QA_V,Bit_SP,QA_Byte)
      IMPLICIT NONE
      SAVE
C-------------------------------------------------------------------
C
C !F77
C
C !DESCRIPTION:
C                 This subroutine is to set QA bit given nth byte
C                 Note: the QA byte setting starts from the rightmost
C                       (not the leftmost) bit
C
C        1st word       7th  6th  5th  4th  3rd  2nd  1st  0th  (bit)
C        2nd word       7th  6th  5th  4th  3rd  2nd  1st  0th  (bit)
C        3rd word       7th  6th  5th  4th  3rd  2nd  1st  0th  (bit)
C        4th word       7th  6th  5th  4th  3rd  2nd  1st  0th  (bit)
C        5th word       7th  6th  5th  4th  3rd  2nd  1st  0th  (bit)
C
C !INPUT PARAMETERS:
C
C        QA_V       QA parameter value
C        Bit_SP     Bit starting position of Ith QA parameter
C                    (see MODIS atmosphere QA plan)
C
C !OUTPUT PARAMETERS:
C
C        QA_Byte    Byte set for quality control
C
C !REVISION HISTORY:
C
C        WRITTEN BY
C        Dr. Allen Chu          11/25/97
C        Code 913
C        NASA Goddard Space Flight Center
C        Greenbelt, MD 20771
C
C !TEAM-UNIQUE HEADER:
C
C   This software was developed by the MODIS Atmosphere Science Team
C   for the National Aeronautics and Space Administration at
C   Goddard Space Flight Center.
C
C !END
C--------------------------------------------------------------------
C
      Intrinsic ibset,ibclr
      INTEGER QA_V,Bit_SP,Byte_Temp
      BYTE QA_Byte
      Byte_Temp=QA_Byte

      IF(QA_V.EQ.0) THEN

        Byte_Temp = ibclr(Byte_Temp,Bit_SP)

      ELSE IF(QA_V.EQ.1) THEN

        Byte_Temp = ibset(Byte_Temp,Bit_SP)

      ELSE IF(QA_V.EQ.2) THEN

        Byte_Temp = ibclr(Byte_Temp,Bit_SP)
        Byte_Temp = ibset(Byte_Temp,Bit_SP+1)

      ELSE IF(QA_V.EQ.3) THEN

        Byte_temp = ibset(Byte_Temp,Bit_SP)
        Byte_temp = ibset(Byte_Temp,Bit_SP+1)

      ENDIF
      QA_Byte=Byte_Temp

      RETURN
      END


