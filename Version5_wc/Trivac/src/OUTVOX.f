*DECK OUTVOX
      SUBROUTINE OUTVOX (IPMAC1,IPVAL,IPMAC2,NBMIX,NL,NBFIS,NGRP,NALBP,
     1 TITR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Driver for the post-treatment of reactor calculation results using
* voxelized flux information.
*
*Copyright:
* Copyright (C) 2026 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A Hebert
*
*Parameters: input
* IPMAC1  L_MACROLIB pointer to the nuclear properties.
* IPVAL   L_FVIEW pointer to the interpflux data structure.
* IPMAC2  L_MACROLIB pointer to the edition information.
* NBMIX   number of material mixtures.
* NL      scattering anisotropy.
* NBFIS   number of fissionable isotopes.
* NGRP    total number of energy groups.
* NALBP   number of physical albedos.
* TITR    title.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAC1,IPMAC2,IPVAL
      CHARACTER TITR*72
      INTEGER NBMIX,NL,NBFIS,NGRP,NALBP
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=20)
      INTEGER ISTATE(NSTATE)
      TYPE(C_PTR) JPMAC1,KPMAC1,JPVAL
      CHARACTER TEXT4*4
      REAL CXYZ(3)
      LOGICAL LVAL
      DOUBLE PRECISION DFLOTT,ZNORM
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, DIMENSION(:), ALLOCATABLE :: IGCOND
      INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: MATXYZ
      REAL, DIMENSION(:), ALLOCATABLE :: SGD
      REAL, DIMENSION(:,:), ALLOCATABLE :: ZUFIS
      REAL, DIMENSION(:,:,:), ALLOCATABLE :: FLU
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IGCOND(NGRP))
*
      TKR=0.0
      IMPX=1
      NGCOND=NGRP
      DO IGR=1,NGRP
        IGCOND(IGR)=IGR
      ENDDO
      LMOD=0
      LVAL=.FALSE.
      ZNORM=1.0D0
      CALL KDRCPU(TK1)
*
   30 CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
      IF(INDIC.NE.3) CALL XABORT('OUTVOX: CHARACTER DATA EXPECTED(2).')
*
   40 IF(TEXT4.EQ.'EDIT') THEN
         CALL REDGET(INDIC,IMPX,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('OUTVOX: INTEGER DATA EXPECTED(1).')
      ELSE IF(TEXT4.EQ.'POWR') THEN
*        NORMALIZATION TO A GIVEN FISSION POWER.
         CALL REDGET (INDIC,NITMA,POWER,TEXT4,DFLOTT)
         IF(INDIC.NE.2) CALL XABORT('OUTVOX: REAL DATA EXPECTED.')
         ZNORM=0.0D0
         JPMAC1=LCMGID(IPMAC1,'GROUP')
         CALL LCMGET(IPVAL,'STATE-VECTOR',ISTATE)
         IF(ISTATE(1).NE.NGRP) CALL XABORT('OUTVOX: INVALID NGRP.')
         CALL LCMGET(IPVAL,'CXYZ',CXYZ)
         VOL=CXYZ(1)*CXYZ(2)*CXYZ(3)
         IXLG=ISTATE(2)
         IYLG=ISTATE(3)
         IZLG=ISTATE(4)
         ALLOCATE(MATXYZ(IXLG,IYLG,IZLG),FLU(IXLG,IYLG,IZLG),SGD(NBMIX))
         CALL LCMGET(IPVAL,'MATXYZ',MATXYZ)
         JPVAL=LCMGID(IPVAL,'FLUX')
         DO IGR=1,NGRP
           KPMAC1=LCMGIL(JPMAC1,IGR)
           CALL LCMLEN(KPMAC1,'H-FACTOR',LENGT,ITYLCM)
           IF(LENGT.GT.0) THEN
             CALL LCMGET(KPMAC1,'H-FACTOR',SGD)
           ELSE
             WRITE(6,'(/43H OUTVOX: *** WARNING *** NO H-FACTOR FOUND ,
     1       28HON LCM. USE NU*SIGF INSTEAD.)')
             ALLOCATE(ZUFIS(NBMIX,NBFIS))
             SGD(:NBMIX)=0.0
             CALL LCMGET(KPMAC1,'NUSIGF',ZUFIS)
             DO IBM=1,NBMIX
               DO IFISS=1,NBFIS
                 SGD(IBM)=SGD(IBM)+ZUFIS(IBM,IFISS)
               ENDDO
             ENDDO
             DEALLOCATE(ZUFIS)
           ENDIF
           CALL LCMGDL(JPVAL,IG,FLU)
           DO K=1,IZLG
             DO J=1,IYLG
               DO I=1,IXLG
                 L=MATXYZ(I,J,K)
                 IF(L.GT.0) ZNORM=ZNORM+FLU(I,J,K)*VOL*SGD(L)
               ENDDO
             ENDDO
           ENDDO
         ENDDO ! IGR
         ZNORM=POWER/ZNORM
         WRITE(6,300) ' DIRECT',ZNORM
         DEALLOCATE(SGD,FLU,MATXYZ)
      ELSE IF(TEXT4.EQ.'FISS') THEN
*        NORMALIZATION TO A GIVEN FISSION SECONDARY NEUTRON PRODUCTION.
         CALL REDGET (INDIC,NITMA,POWER,TEXT4,DFLOTT)
         IF(INDIC.NE.2) CALL XABORT('OUTVOX: REAL DATA EXPECTED.')
         ZNORM=0.0D0
         JPMAC1=LCMGID(IPMAC1,'GROUP')
         CALL LCMGET(IPVAL,'STATE-VECTOR',ISTATE)
         IF(ISTATE(1).NE.NGRP) CALL XABORT('OUTVOX: INVALID NGRP.')
         CALL LCMGET(IPVAL,'CXYZ',CXYZ)
         VOL=CXYZ(1)*CXYZ(2)*CXYZ(3)
         IXLG=ISTATE(2)
         IYLG=ISTATE(3)
         IZLG=ISTATE(4)
         ALLOCATE(MATXYZ(IXLG,IYLG,IZLG),FLU(IXLG,IYLG,IZLG),SGD(NBMIX))
         CALL LCMGET(IPVAL,'MATXYZ',MATXYZ)
         JPVAL=LCMGID(IPVAL,'FLUX')
         DO IGR=1,NGRP
           KPMAC1=LCMGIL(JPMAC1,IGR)
           CALL LCMLEN(KPMAC1,'NUSIGF',LENGT,ITYLCM)
           IF(LENGT.EQ.0) THEN
             CALL LCMLIB(KPMAC1)
             CALL XABORT('OUTVOX: NUSIGF RECORD MISSING IN MACROLIB.')
           ENDIF
           ALLOCATE(ZUFIS(NBMIX,NBFIS))
           SGD(:NBMIX)=0.0
           CALL LCMGET(KPMAC1,'NUSIGF',ZUFIS)
           DO IBM=1,NBMIX
             DO IFISS=1,NBFIS
               SGD(IBM)=SGD(IBM)+ZUFIS(IBM,IFISS)
             ENDDO
           ENDDO
           DEALLOCATE(ZUFIS)
           CALL LCMGDL(JPVAL,IG,FLU)
           DO K=1,IZLG
             DO J=1,IYLG
               DO I=1,IXLG
                 L=MATXYZ(I,J,K)
                 IF(L.GT.0) ZNORM=ZNORM+FLU(I,J,K)*VOL*SGD(L)
               ENDDO
             ENDDO
           ENDDO
         ENDDO ! IGR
         ZNORM=POWER/ZNORM
         WRITE(6,300) ' DIRECT',ZNORM
         DEALLOCATE(SGD,FLU,MATXYZ)
      ELSE IF(TEXT4.EQ.'SOUR') THEN
*        NORMALIZATION TO A GIVEN SOURCE INTENSITY.
         CALL REDGET (INDIC,NITMA,SNUMB,TEXT4,DFLOTT)
         IF(INDIC.NE.2) CALL XABORT('OUTVOX: REAL DATA EXPECTED.')
         ZNORM=0.0D0
         JPMAC1=LCMGID(IPMAC1,'GROUP')
         CALL LCMGET(IPVAL,'STATE-VECTOR',ISTATE)
         IF(ISTATE(1).NE.NGRP) CALL XABORT('OUTVOX: INVALID NGRP.')
         IXLG=ISTATE(2)
         IYLG=ISTATE(3)
         IZLG=ISTATE(4)
         CALL LCMGET(IPVAL,'CXYZ',CXYZ)
         VOL=CXYZ(1)*CXYZ(2)*CXYZ(3)
         ALLOCATE(MATXYZ(IXLG,IYLG,IZLG),SGD(NBMIX))
         CALL LCMGET(IPVAL,'MATXYZ',MATXYZ)
         DO IGR=1,NGRP
           KPMAC1=LCMGIL(JPMAC1,IGR)
           CALL LCMLEN(KPMAC1,'FIXE',LENGT,ITYLCM)
           IF(LENGT.EQ.0) THEN
             CALL LCMLIB(KPMAC1)
             CALL XABORT('OUTVOX: SOURCE RECORD MISSING IN MACROLIB.')
           ENDIF
           CALL LCMGET(KPMAC1,'FIXE',SGD)
           DO K=1,IZLG
             DO J=1,IYLG
               DO I=1,IXLG
                 L=MATXYZ(I,J,K)
                 IF(L.GT.0) ZNORM=ZNORM+VOL*SGD(L)
               ENDDO
             ENDDO
           ENDDO
         ENDDO ! IGR
         ZNORM=SNUMB/ZNORM
         WRITE(6,305) ' DIRECT',ZNORM
         DEALLOCATE(SGD,MATXYZ)
      ELSE IF(TEXT4.EQ.'COND') THEN
         NGCOND=0
         CALL REDGET (INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.EQ.3) THEN
           IF(TEXT4.EQ.'NONE') THEN
             NGCOND=NGRP
             DO IGR=1,NGRP
               IGCOND(IGR)=IGR
             ENDDO
             GO TO 30
           ENDIF
           NGCOND=1
           IGCOND(NGCOND)=NGRP
           GO TO 40
         ELSE IF(INDIC.EQ.1) THEN
  170      IF(NITMA.GT.NGRP) NITMA=NGRP
           NGCOND=NGCOND+1
           IGCOND(NGCOND)=NITMA
           CALL REDGET (INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
           IF(INDIC.EQ.1) THEN
             GO TO 170
           ELSE IF(INDIC.EQ.3) THEN
             GO TO 40
           ELSE
             CALL XABORT('OUTVOX: INTEGER OR CHARACTER DATA EXPECTED(3'
     1       //').')
           ENDIF
         ELSE
           CALL XABORT('OUTVOX: INTEGER OR CHARACTER DATA EXPECTED(4).')
         ENDIF
      ELSE IF(TEXT4.EQ.'INTG') THEN
*        READ THE MERGE INDICES.
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
         IF((INDIC.EQ.3).AND.(TEXT4.EQ.'VAL')) THEN
           IF(IMPX.GT.0) WRITE(6,320) TITR
           LVAL=.TRUE.
           IF(IMPX.GT.0) WRITE(6,330) (IGCOND(IG),IG=1,NGCOND)
           CALL OUTVAL(IPMAC1,IPMAC2,IPVAL,NBMIX,NL,NBFIS,NGRP,NALBP,
     1     NGCOND,IGCOND,ZNORM,IMPX)
           GO TO 180
         ELSE
           CALL XABORT('OUTVOX: INVALID KEY WORD.')
         ENDIF
      ELSE
         CALL XABORT('OUTVOX: '//TEXT4//' IS AN INVALID KEY WORD.')
      ENDIF
      GO TO 30
*----
*  SCRATCH STORAGE DEALLOCATION
*----
  180 DEALLOCATE(IGCOND)
      CALL KDRCPU(TK2)
      TKR=TK2-TK1
      WRITE(6,310) TKR
      RETURN
*
  300 FORMAT(/9H OUTVOX: ,A7,28H FLUX NORMALIZATION FACTOR =,1P,E13.5)
  305 FORMAT(/9H OUTVOX: ,A7,30H SOURCE NORMALIZATION FACTOR =,1P,E13.5)
  310 FORMAT(/49H OUTVOX: CPU TIME FOR REACTION RATE CALCULATION =,F7.3)
  320 FORMAT(/12H OUTVOX: ***,A72,3H***)
  330 FORMAT(/20H CONDENSATION INDEX:/(1X,14I5))
      END
