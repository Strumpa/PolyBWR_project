*DECK OUTKIN
      SUBROUTINE OUTKIN (IPMAC1,IPMAC2,NBMIX,NGRP,NEL,NUN,MAT,VOL,IDL,
     1 EVECT,ADECT,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform direct-adjoint homogenization and condensation for producing
* point kinetics parameters.
*
*Copyright:
* Copyright (C) 2026 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPMAC1  L_MACROLIB pointer to the input macrolib.
* IPMAC2  L_MACROLIB pointer to the output extended macrolib.
* NBMIX   number of material mixtures.
* NGRP    total number of energy groups.
* NEL     number of volumes.
* NUN     number of unknowns per energy group.
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* IDL     position of the average flux component associated with
*         each volume.
* EVECT   unknowns.
* ADECT   adjoint flux unknowns.
* IMPX    print parameter (equal to zero for no print).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAC1,IPMAC2
      PARAMETER(NREAC=11)
      INTEGER NBMIX,NBFIS,NGRP,NEL,NUN,MAT(NEL),IDL(NEL),IMPX
      REAL VOL(NEL),EVECT(NUN,NGRP),ADECT(NUN,NGRP)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPMAC1,KPMAC1
      PARAMETER(NSTATE=40)
      INTEGER ISTATE(NSTATE)
      DOUBLE PRECISION TGENRS,DNUM,DDEN
      CHARACTER TEXT12*12
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, DIMENSION(:), ALLOCATABLE :: RBETA,RLAMB,ZOVER
      REAL, DIMENSION(:,:), ALLOCATABLE :: ZUFIS,DLAMB,DBETA
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ACCUM
*----
*  RECOVER STATE-VECTOR INFORMATION
*----
      CALL LCMGET(IPMAC1,'STATE-VECTOR',ISTATE)
      NBFIS=ISTATE(4)
      NPRC=ISTATE(7)
      IF(NBFIS.EQ.0) THEN
        CALL XABORT('OUTKIN: NO FISSILE ISOTOPES IN INPUT MACROLIB.')
      ELSE IF(NPRC.EQ.0) THEN
        CALL XABORT('OUTKIN: NO DELAYED NEUTRONS IN INPUT MACROLIB.')
      ENDIF
      ALLOCATE(DLAMB(NPRC,NBFIS),DBETA(NPRC,NBFIS))
      CALL LCMGET(IPMAC1,'LAMBDA-D',DLAMB)
      CALL LCMLEN(IPMAC1,'BETA-D',ILBETA,ITYLCM)
      IF(ILBETA.GT.0) CALL LCMGET(IPMAC1,'BETA-D',DBETA)
*----
*  COMPUTE LAMBDA AND BETA VALUES
*----
      ALLOCATE(ACCUM(0:NPRC,NBFIS),ZOVER(NBMIX),ZUFIS(NBMIX,NBFIS))
      JPMAC1=LCMGID(IPMAC1,'GROUP')
      TGENRS=0.0
      ACCUM(0:NPRC,:NBFIS)=0.0D0
      DO IGR=1,NGRP
        KPMAC1=LCMGIL(JPMAC1,IGR)
        CALL LCMGET(KPMAC1,'OVERV',ZOVER)
        TGENRS=0.0D0
        DO K=1,NEL
          IBM=MAT(K)
          IPFL=IDL(K)
          IF((IBM.NE.0).AND.(IPFL.NE.0)) THEN
             TGENRS=TGENRS+VOL(K)*ADECT(IPFL,IGR)*EVECT(IPFL,IGR)*
     1       ZOVER(IBM)
          ENDIF
        ENDDO
        CALL LCMGET(KPMAC1,'NUSIGF',ZUFIS)
        DO IFIS=1,NBFIS
          DO K=1,NEL
            IBM=MAT(K)
            IPFL=IDL(K)
            IF((IBM.NE.0).AND.(IPFL.NE.0)) THEN
              ACCUM(0,IFIS)=ACCUM(0,IFIS)+VOL(K)*ADECT(IPFL,IGR)*
     1        EVECT(IPFL,IGR)*ZUFIS(IBM,IFIS)
              IF(ILBETA.NE.0) THEN
                DO IPRC=1,NPRC
                  ACCUM(IPRC,IFIS)=ACCUM(IPRC,IFIS)+DBETA(IPRC,IFIS)*
     1            VOL(K)*ADECT(IPFL,IGR)*EVECT(IPFL,IGR)*ZUFIS(IBM,IFIS)
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDDO ! IFIS
        IF(ILBETA.EQ.0) THEN
          DO IPRC=1,NPRC
            WRITE(TEXT12,'(A6,I2.2)') 'NUSIGF',IPRC
            CALL LCMLEN(KPMAC1,TEXT12,ILONG,ITYLCM)
            IF(ILONG.EQ.0) CALL XABORT('OUTKIN: MISSING NUSIGF01.')
            CALL LCMGET(KPMAC1,TEXT12,ZUFIS)
            DO IFIS=1,NBFIS
              DO K=1,NEL
                IBM=MAT(K)
                IPFL=IDL(K)
                IF((IBM.NE.0).AND.(IPFL.NE.0)) THEN
                  ACCUM(IPRC,IFIS)=ACCUM(IPRC,IFIS)+VOL(K)*
     1            ADECT(IPFL,IGR)*EVECT(IPFL,IGR)*ZUFIS(IBM,IFIS)
                ENDIF
              ENDDO
            ENDDO ! IFIS
          ENDDO ! IPRC
        ENDIF
      ENDDO ! IGR
      DEALLOCATE(ZUFIS,ZOVER)
*----
*  AVERAGE PRECURSOR DECAY CONSTANTS. IT IS RECOMMENDED TO USE 8-GROUP
*  PRECURSOR DECAY CONSTANTS, IDENTICAL FOR ALL ISOTOPES.
*----
      ALLOCATE(RLAMB(NPRC),RBETA(NPRC))
      DO IPRC=1,NPRC
        DNUM=0.0D0
        DDEN=0.0D0
        DO IFIS=1,NBFIS
          DNUM=DNUM+DLAMB(IPRC,IFIS)*ACCUM(IPRC,IFIS)
          DDEN=DDEN+ACCUM(IPRC,IFIS)
        ENDDO
        RLAMB(IPRC)=REAL(DNUM/DDEN)
      ENDDO
*----
*  STORE POINT KINETICS INFORMATION.
*----
      CALL LCMSIX(IPMAC2,'P-KINETICS',1)
      DO IPRC=1,NPRC
        RBETA(IPRC)=REAL(SUM(ACCUM(IPRC,:NBFIS))/SUM(ACCUM(0,:NBFIS)))
      ENDDO
      CALL LCMPUT(IPMAC2,'BETAI',NPRC,2,RBETA)
      CALL LCMPUT(IPMAC2,'LAMBDAI',NPRC,2,DLAMB(1,1))
      RGENRS=REAL(TGENRS/SUM(ACCUM(0,:NBFIS)))
      CALL LCMPUT(IPMAC2,'LAMBDA',1,2,RGENRS)
      CALL LCMSIX(IPMAC2,' ',2)
      DEALLOCATE(ACCUM,DBETA,DLAMB)
*----
*  PRINT POINT KINETICS INFORMATION.
*----
      IF(IMPX.GT.0) THEN
        WRITE(6,'(/36H OUTKIN: POINT KINETICS INFORMATION:)')
        WRITE(6,'(17X,15HEFFECTIVE BETA:,1P,E12.4)') SUM(RBETA(:NPRC))
        WRITE(6,'(8X,24HPROMPT-NEUTRON LIFETIME:,1P,E12.4)') RGENRS
        WRITE(6,'(6X,26HPRECURSOR DECAY CONSTANTS:,1P,10E12.4)') RLAMB
        WRITE(6,'(6X,26HDELAYED NEUTRON FRACTIONS:,1P,10E12.4)') RBETA
      ENDIF
      DEALLOCATE(RBETA,RLAMB)
      RETURN
      END
