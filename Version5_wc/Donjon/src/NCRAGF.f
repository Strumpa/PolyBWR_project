*DECK NCRAGF
      SUBROUTINE NCRAGF(IPMAC,IPCPO,IACCS,NMIL,NMIX,NGRP,NGFF,NALBP,
     1 IMPX,NCAL,TERP,MIXC,IDF,NTYPE,NFINF)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Build the macrolib by scanning the NCAL elementary calculations and
* weighting them with TERP factors. ADF, GFF and physical albedos part.
*
*Copyright:
* Copyright (C) 2015 Ecole Polytechnique de Montreal
*
*Author(s): 
* R. Chambon, A. Hebert
*
*Parameters: input
* IPMAC   address of the output macrolib LCM object.
* IPCPO   address of the multicompo object.
* IACCS   =0 macrolib is created; =1 ... is updated.
* NMIL    number of material mixtures in the multicompo.
* NMIX    maximum number of material mixtures in the macrolib.
* NGRP    number of energy groups.
* NGFF    number of group form factors per energy group.
* NALBP   number of physical albedos per energy group.
* IMPX    print parameter (equal to zero for no print).
* NCAL    number of elementary calculations in the multicompo.
* TERP    interpolation factors.
* MIXC    mixture index in the multicompo corresponding to each macrolib
*         mixture. Equal to zero if a macrolib mixture is not updated.
* IDF     ADF type, 0 = none, 1 = Albedo, 2 = FD_B/FD_C/..., 3 = ADF.
* NTYPE   number of ADF.
* NFINF   number of 'enriched' flux (for pin power reconstruction in
*         NAP:).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAC,IPCPO
      INTEGER IACCS,NMIL,NMIX,NGRP,NGFF,NALBP,IMPX,NCAL,MIXC(NMIX),IDF,
     1 NTYPE,NFINF
      REAL TERP(NCAL,NMIX)
*----
*  LOCAL VARIABLES
*----
      INTEGER, PARAMETER::IOUT=6
      INTEGER, PARAMETER::MAXIFX=5
      INTEGER, PARAMETER::NSTATE=40
      INTEGER FINF(MAXIFX),NITMA
      REAL WEIGHT,FACTOR,ZZZ
      CHARACTER FINFN*8,HSMG*131
      TYPE(C_PTR) JPCPO,KPCPO,LPCPO,MPCPO
      INTEGER IKEFF,IKINF,I,IBM,IBMOLD,ICAL,IGR,JGR,IGFF,ILONG,ITYLCM,
     1 ITYPE,ITYP2,JTYPE,IAL,NTYPE2
      INTEGER ISTATE(NSTATE)
      DOUBLE PRECISION GAR1,GAR2
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: GAR4,VOL,ZKINF,ZKEFF
      REAL, ALLOCATABLE, DIMENSION(:,:) :: GAR6,ALBP
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: GAR5,ADF2,ALBP2
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: GFF,ADF2M
      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:) :: HADF,HADF2
*----
*  SCRATCH STORAGE ALLOCATION
*----
        PRINT *,'NCRAGF: NMIX=',NMIX
      ALLOCATE(GAR4(NGRP*NGRP),GAR6(NGRP,2),GFF(NMIX,NGFF,NGRP,2+NFINF),
     1 GAR5(NGFF,NGRP,2+MAXIFX),ALBP(NALBP,NGRP),ALBP2(NMIX,NALBP,NGRP),
     2 ZKINF(NMIX),ZKEFF(NMIX),HADF(NTYPE),ADF2(NMIX,NGRP,NTYPE),
     3 ADF2M(NMIX,NGRP,NGRP,NTYPE))
*----
*  OVERALL MULTICOMPO MIXTURE LOOP
*----
      IKINF=0
      IKEFF=0
      JPCPO=LCMGID(IPCPO,'MIXTURES')
      IF(NALBP.NE.0) ALBP2(:NMIX,:NALBP,:NGRP)=0.0
      ZKINF(:NMIX)=0.0
      ZKEFF(:NMIX)=0.0
      DO 500 IBMOLD=1,NMIL
      IF(IMPX.GT.0) WRITE(IOUT,'(/33H NCRAGF: PROCESS MULTICOMPO MIXTU,
     1 2HRE,I5)') IBMOLD
      KPCPO=LCMGIL(JPCPO,IBMOLD)
      LPCPO=LCMGID(KPCPO,'CALCULATIONS')
*----
*  READ EXISTING MACROLIB INFORMATION
*----
      MPCPO=LCMGIL(LPCPO,1)
      CALL LCMGET(MPCPO,'STATE-VECTOR',ISTATE)
      IF(ISTATE(1).NE.1) CALL XABORT('NCRAGF: THE NUMBER OF MIXTURE SH'
     1 //'OULD ALWAYS BE EQUAL TO 1 IN A MULTICOMPO MICROLIB BRANCH.')
      IF(IACCS.EQ.0) THEN  !IACCS
        IF((IDF.NE.0).OR.(NGFF.NE.0)) CALL LCMSIX(MPCPO,'MACROLIB',1)
        IF(IDF.NE.0) THEN
          !copy ADF names from multicompo
          CALL LCMSIX(IPMAC,'ADF',1)
          CALL LCMLEN(MPCPO,'ADF',ILONG,ITYLCM)
          IF(ILONG.EQ.0) CALL XABORT('NCRAGF: MISSING ADF DIRECTORY I'
     1    //'N MULTICOMPO OBJECT.')
          CALL LCMSIX(MPCPO,'ADF',1)
          CALL LCMEQU(MPCPO,IPMAC)
          IF(IDF.EQ.1) THEN
            CALL LCMLEN(IPMAC,'ALBS00',ILONG,ITYLCM)
            IF(ILONG.GT.0) CALL LCMDEL(IPMAC,'ALBS00')
            ADF2(:NMIX,:NGRP,:NTYPE)=0.0
          ELSE IF((IDF.EQ.2).OR.(IDF.EQ.3)) THEN
            CALL LCMGET(MPCPO,'NTYPE',NITMA)
            IF(NITMA.NE.NTYPE) CALL XABORT('NCRAGF: INVALID NTYPE(1).')
            IF(NTYPE.GT.0) THEN
              CALL LCMGTC(MPCPO,'HADF',8,NTYPE,HADF)
              DO ITYPE=1,NTYPE
                CALL LCMLEN(IPMAC,HADF(ITYPE),ILONG,ITYLCM)
                IF(ILONG.GT.0) CALL LCMDEL(IPMAC,HADF(ITYPE))
              ENDDO
            ENDIF
            ADF2(:NMIX,:NGRP,:NTYPE)=0.0
          ELSE IF(IDF.EQ.4) THEN
            CALL LCMGET(MPCPO,'NTYPE',NITMA)
            IF(NITMA.NE.NTYPE) CALL XABORT('NCRAGF: INVALID NTYPE(2).')
            IF(NTYPE.GT.0) THEN
              CALL LCMGTC(MPCPO,'HADF',8,NTYPE,HADF)
              DO ITYPE=1,NTYPE
                CALL LCMLEN(IPMAC,HADF(ITYPE),ILONG,ITYLCM)
                IF(ILONG.GT.0) CALL LCMDEL(IPMAC,HADF(ITYPE))
              ENDDO
            ENDIF
            ADF2M(:NMIX,:NGRP,:NGRP,:NTYPE)=0.0
          ENDIF
          CALL LCMSIX(MPCPO,' ',2)
          CALL LCMSIX(IPMAC,' ',2)
        ENDIF
        IF(NGFF.NE.0) THEN
          !copy GFF geom and FINF names from multicompo
          CALL LCMSIX(IPMAC,'GFF',1)
          CALL LCMLEN(MPCPO,'GFF',ILONG,ITYLCM)
          IF(ILONG.EQ.0) CALL XABORT('NCRAGF: MISSING GFF DIRECTORY I'
     1    //'N MULTICOMPO OBJECT.')
          CALL LCMSIX(MPCPO,'GFF',1)
          CALL LCMEQU(MPCPO,IPMAC)
          IF(NFINF.GT.0) THEN
            CALL LCMGET(IPMAC,'FINF_NUMBER ',FINF)
            DO I=1,NFINF
              WRITE(FINFN,'(5HFINF_,I3.3)') FINF(I)
              CALL LCMLEN(IPMAC,FINFN,ILONG,ITYLCM)
              IF(ILONG.GT.0) CALL LCMDEL(IPMAC,FINFN)
            ENDDO
          ENDIF
          CALL LCMSIX(MPCPO,' ',2)
          CALL LCMSIX(IPMAC,' ',2)
          GFF(:NMIX,:NGFF,:NGRP,:2+NFINF)=0.0
        ENDIF
        IF((IDF.NE.0).OR.(NGFF.NE.0)) CALL LCMSIX(MPCPO,' ',2)
        CALL LCMGET(IPMAC,'STATE-VECTOR',ISTATE)
        ISTATE(8)=NALBP
        ISTATE(12)=IDF
        ISTATE(16)=NGFF
        CALL LCMPUT(IPMAC,'STATE-VECTOR',NSTATE,1,ISTATE)
        IACCS=1
      ELSE !IACCS
* Recover ADF, GFF and physical albedos previously computed
        IF(NGFF.NE.0) THEN
          CALL LCMSIX(IPMAC,'GFF',1)
          CALL LCMGET(IPMAC,'NWT0',GFF(1,1,1,1))
          CALL LCMGET(IPMAC,'H-FACTOR',GFF(1,1,1,2))
          IF(NFINF.GT.0) THEN
            CALL LCMGET(IPMAC,'FINF_NUMBER ',FINF)
            DO I=1,NFINF
              WRITE(FINFN,'(5HFINF_,I3.3)') FINF(I)
              CALL LCMGET(IPMAC,FINFN,GFF(1,1,1,2+I))
            ENDDO
          ENDIF
          DO IBM=1,NMIX
            IF(MIXC(IBM).EQ.IBMOLD) GFF(IBM,:NGFF,:NGRP,:NFINF+2)=0.0
          ENDDO
          CALL LCMSIX(IPMAC,' ',2)
        ENDIF
        IF(IDF.NE.0) THEN
          CALL LCMSIX(IPMAC,'ADF',1)
          IF(IDF.EQ.1) THEN
            DO IBM=1,NMIX
              IF(MIXC(IBM).EQ.IBMOLD) ADF2(IBM,:NGRP,1)=0.0
            ENDDO
          ELSE IF((IDF.EQ.2).OR.(IDF.EQ.3)) THEN
            CALL LCMGTC(IPMAC,'HADF',8,NTYPE,HADF)
            DO ITYPE=1,NTYPE
              CALL LCMGET(IPMAC,HADF(ITYPE),ADF2(1,1,ITYPE))
              DO IBM=1,NMIX
                IF(MIXC(IBM).EQ.IBMOLD) ADF2(IBM,:NGRP,ITYPE)=0.0
              ENDDO
            ENDDO
          ELSE IF(IDF.EQ.4) THEN
            CALL LCMGTC(IPMAC,'HADF',8,NTYPE,HADF)
            DO ITYPE=1,NTYPE
              CALL LCMGET(IPMAC,HADF(ITYPE),ADF2M(1,1,1,ITYPE))
              DO IBM=1,NMIX
                IF(MIXC(IBM).EQ.IBMOLD) ADF2M(IBM,:NGRP,:NGRP,ITYPE)=0.0
              ENDDO
            ENDDO
          ENDIF
          CALL LCMSIX(IPMAC,' ',2)
        ENDIF
        DO IBM=1,NMIX
          IF(MIXC(IBM).EQ.IBMOLD) THEN
            IF(NALBP.NE.0) ALBP2(IBM,:NALBP,:NGRP)=0.0
            ZKINF(IBM)=0.0
            ZKEFF(IBM)=0.0
          ENDIF
        ENDDO
        CALL LCMGET(IPMAC,'STATE-VECTOR',ISTATE)
        ISTATE(8)=NALBP
        ISTATE(12)=IDF
        ISTATE(16)=NGFF
        CALL LCMPUT(IPMAC,'STATE-VECTOR',NSTATE,1,ISTATE)
*
      ENDIF  !IACCS
*----
*  OVERALL ELEMENTARY CALCULATION LOOP
*----
      DO 210 ICAL=1,NCAL
      MPCPO=LCMGIL(LPCPO,ICAL)
      DO 200 IBM=1,NMIX
      WEIGHT=TERP(ICAL,IBM)
      IF((MIXC(IBM).NE.IBMOLD).OR.(WEIGHT.EQ.0.0)) GO TO 200
*----
*  PERFORM INTERPOLATION
*----
*----
*  PROCESS GROUP FORM FACTOR (GFF) INFORMATION
*----
      IF(NGFF.NE.0) THEN
        CALL LCMSIX(MPCPO,'MACROLIB',1)
        CALL LCMLEN(MPCPO,'GFF',ILONG,ITYLCM)
        IF(ILONG.NE.0) THEN
          CALL LCMSIX(MPCPO,'GFF',1)
          CALL LCMLEN(MPCPO,'NWT0',ILONG,ITYLCM)
          IF(ILONG.GT.NGFF*NGRP*(2+MAXIFX)) THEN
            CALL LCMLIB(MPCPO)
            WRITE(6,'(6H NGFF=,I6,6H NGRP=,I6,11H LEN(NWT0)=,I6)')
     >      NGFF,NGRP,ILONG
            CALL XABORT('NCRAGF: MAXIFX OVERFLOW.')
          ENDIF
          CALL LCMGET(MPCPO,'NWT0',GAR5(1,1,1))
          CALL LCMGET(MPCPO,'H-FACTOR',GAR5(1,1,2))
          CALL LCMLEN(MPCPO,'FINF_NUMBER ',NFINF,ITYLCM)
          IF(NFINF.GT.0) THEN
            CALL LCMGET(MPCPO,'FINF_NUMBER ',FINF)
            DO I=1,NFINF
              WRITE(FINFN,'(5HFINF_,I3.3)') FINF(I)
              CALL LCMGET(MPCPO,FINFN,GAR5(1,1,2+I))
            ENDDO
          ENDIF
          DO IGFF=1,NGFF
            DO IGR=1,NGRP
              GFF(IBM,IGFF,IGR,1)=GFF(IBM,IGFF,IGR,1)
     1                            +WEIGHT*GAR5(IGFF,IGR,1)
              GFF(IBM,IGFF,IGR,2)=GFF(IBM,IGFF,IGR,2)
     1                            +WEIGHT*GAR5(IGFF,IGR,2)
              DO I=1,NFINF
                GFF(IBM,IGFF,IGR,2+I)=GFF(IBM,IGFF,IGR,2+I)
     1                           +WEIGHT*GAR5(IGFF,IGR,2+I)
              ENDDO
            ENDDO
          ENDDO
          CALL LCMSIX(MPCPO,' ',2)
        ENDIF
        CALL LCMSIX(MPCPO,' ',2)
      ENDIF
*----
*  PROCESS ADF INFORMATION
*----
      IF(IDF.NE.0) THEN
        CALL LCMSIX(MPCPO,'MACROLIB',1)
        CALL LCMLEN(MPCPO,'ADF',ILONG,ITYLCM)
        IF(ILONG.NE.0) THEN
          CALL LCMSIX(MPCPO,'ADF',1)
          IF(IDF.EQ.1) THEN
            GAR6(:NGRP,:2)=0.0
            CALL LCMGET(MPCPO,'ALBS00',GAR6)
            DO IGR=1,NGRP
              ADF2(IBM,IGR,:2)=ADF2(IBM,IGR,:2)+WEIGHT*GAR6(IGR,:2)
            ENDDO
          ELSE IF((IDF.EQ.2).OR.(IDF.EQ.3)) THEN
            CALL LCMGET(MPCPO,'NTYPE',NTYPE2)
            ALLOCATE(HADF2(NTYPE2))
            CALL LCMGTC(MPCPO,'HADF',8,NTYPE2,HADF2)
            IF(NTYPE2.EQ.1) THEN
*             assign the same ADF to all sides.
              CALL LCMLEN(MPCPO,HADF2(1),ILONG,ITYLCM)
              IF(ILONG.NE.NGRP) CALL XABORT('NCRAGF: INVALID ADF LENGT'
     1        //'H(1).')
              CALL LCMGET(MPCPO,HADF2(1),GAR4)
              DO ITYPE=1,NTYPE
                DO IGR=1,NGRP
                  ADF2(IBM,IGR,ITYPE)=ADF2(IBM,IGR,ITYPE)+WEIGHT*
     1            GAR4(IGR)
                ENDDO
              ENDDO
            ELSE
              IF(NTYPE2.GT.NTYPE) CALL XABORT('NCRAGF: NTYPE OVERFLOW.')
              DO ITYP2=1,NTYPE2
                ITYPE=0
                DO JTYPE=1,NTYPE
                  IF(HADF2(ITYP2).EQ.HADF(JTYPE)) THEN
                    ITYPE=JTYPE
                    GO TO 180
                  ENDIF
                ENDDO
                WRITE(HSMG,'(18HNCRAGF: ADF NAMED ,A,11H NOT FOUND.)')
     1          TRIM(HADF2(ITYP2))
                CALL XABORT(HSMG)
  180           CALL LCMLEN(MPCPO,HADF2(ITYP2),ILONG,ITYLCM)
                IF(ILONG.NE.NGRP) CALL XABORT('NCRAGF: INVALID ADF LEN'
     1          //'GTH(2).')
                CALL LCMGET(MPCPO,HADF2(ITYP2),GAR4)
                DO IGR=1,NGRP
                  ADF2(IBM,IGR,ITYPE)=ADF2(IBM,IGR,ITYPE)+WEIGHT*
     1            GAR4(IGR)
                ENDDO
              ENDDO
            ENDIF
            DEALLOCATE(HADF2)
          ELSE IF(IDF.EQ.4) THEN
            CALL LCMGET(MPCPO,'NTYPE',NTYPE2)
            ALLOCATE(HADF2(NTYPE2))
            CALL LCMGTC(MPCPO,'HADF',8,NTYPE2,HADF2)
            IF(NTYPE2.EQ.1) THEN
*             assign the same MADF to all sides.
              CALL LCMLEN(MPCPO,HADF2(1),ILONG,ITYLCM)
              IF(ILONG.NE.NGRP*NGRP) CALL XABORT('NCRAGF: INVALID ADFM'
     1        //'LENGTH(1).')
              CALL LCMGET(MPCPO,HADF2(1),GAR4)
              DO ITYPE=1,NTYPE
                DO JGR=1,NGRP
                  DO IGR=1,NGRP
                    ADF2M(IBM,IGR,JGR,ITYPE)=ADF2M(IBM,IGR,JGR,ITYPE)+
     1              WEIGHT*GAR4((JGR-1)*NGRP+IGR)
                  ENDDO
                ENDDO
              ENDDO
            ELSE
              IF(NTYPE2.GT.NTYPE) CALL XABORT('NCRAGF: NTYPE OVERFLOW.')
              DO ITYP2=1,NTYPE2
                ITYPE=0
                DO JTYPE=1,NTYPE
                  IF(HADF2(ITYP2).EQ.HADF(JTYPE)) THEN
                    ITYPE=JTYPE
                    GO TO 190
                  ENDIF
                ENDDO
                WRITE(HSMG,'(19HNCRAGF: ADFM NAMED ,A,11H NOT FOUND.)')
     1          TRIM(HADF2(ITYP2))
                CALL XABORT(HSMG)
                CALL LCMLEN(MPCPO,HADF2(ITYP2),ILONG,ITYLCM)
  190           IF(ILONG.NE.NGRP*NGRP) CALL XABORT('NCRAGF: INVALID AD'
     1          //'FM LENGTH(2).')
                CALL LCMGET(MPCPO,HADF2(ITYP2),GAR4)
                DO JGR=1,NGRP
                  DO IGR=1,NGRP
                    ADF2M(IBM,IGR,JGR,ITYPE)=ADF2M(IBM,IGR,JGR,ITYPE)+
     1              WEIGHT*GAR4((JGR-1)*NGRP+IGR)
                  ENDDO
                ENDDO
              ENDDO
            ENDIF
            DEALLOCATE(HADF2)
          ENDIF
          CALL LCMSIX(MPCPO,' ',2)
        ENDIF
        CALL LCMSIX(MPCPO,' ',2)
      ENDIF
*----
*  PROCESS PHYSICAL ALBEDO INFORMATION
*----
      IF(NALBP.NE.0) THEN
        CALL LCMSIX(MPCPO,'MACROLIB',1)
        CALL LCMGET(MPCPO,'ALBEDO',ALBP)
        DO IGR=1,NGRP
          DO IAL=1,NALBP
            FACTOR=(1.0-ALBP(IAL,IGR))/(1.0+ALBP(IAL,IGR))
            ALBP2(IBM,IAL,IGR)=ALBP2(IBM,IAL,IGR)+WEIGHT*FACTOR
          ENDDO
        ENDDO
        CALL LCMSIX(MPCPO,' ',2)
      ENDIF
*----
*  PROCESS KINF
*----
      CALL LCMLEN(MPCPO,'K-INFINITY',IKINF,ITYLCM)
      IF(IKINF.EQ.1) THEN
        CALL LCMGET(MPCPO,'K-INFINITY',ZZZ)
        ZKINF(IBM)=ZKINF(IBM)+WEIGHT*ZZZ
      ENDIF
*----
*  PROCESS KEFF
*----
      CALL LCMLEN(MPCPO,'K-EFFECTIVE',IKEFF,ITYLCM)
      IF(IKEFF.EQ.1) THEN
        CALL LCMGET(MPCPO,'K-EFFECTIVE',ZZZ)
        ZKEFF(IBM)=ZKEFF(IBM)+WEIGHT*ZZZ
      ENDIF
  200 CONTINUE
  210 CONTINUE
*----
*  WRITE INTERPOLATED MACROLIB INFORMATION
*----
      IF(IDF.EQ.1) THEN
        CALL LCMSIX(IPMAC,'ADF',1)
        CALL LCMPUT(IPMAC,'ALBS00',NMIX*NGRP*2,2,ADF2(1,1,1))
        CALL LCMSIX(IPMAC,' ',2)
      ELSE IF((IDF.EQ.2).OR.(IDF.EQ.3)) THEN
        CALL LCMSIX(IPMAC,'ADF',1)
        CALL LCMGTC(IPMAC,'HADF',8,NTYPE,HADF)
        DO ITYPE=1,NTYPE
          CALL LCMPUT(IPMAC,HADF(ITYPE),NMIX*NGRP,2,
     1    ADF2(1,1,ITYPE))
        ENDDO
        CALL LCMSIX(IPMAC,' ',2)
        IF(IMPX.GT.1) THEN
          DO IBM=1,NMIX
            IF(MIXC(IBM).EQ.0) CYCLE
            WRITE(6,'(/40H NCRAGF: DISCONTINUITY FACTORS - MIXTURE,I5)')
     1      IBM
            DO ITYPE=1,NTYPE
              WRITE(6,'(1X,A,1H:,1P,(5X,10E12.4))') TRIM(HADF(ITYPE)),
     1        (ADF2(IBM,IGR,ITYPE),IGR=1,NGRP)
            ENDDO
          ENDDO
        ENDIF
      ELSE IF(IDF.EQ.4) THEN
        CALL LCMSIX(IPMAC,'ADF',1)
        CALL LCMGTC(IPMAC,'HADF',8,NTYPE,HADF)
        DO ITYPE=1,NTYPE
          CALL LCMPUT(IPMAC,HADF(ITYPE),NMIX*NGRP*NGRP,2,
     1    ADF2M(1,1,1,ITYPE))
        ENDDO
        CALL LCMSIX(IPMAC,' ',2)
        IF(IMPX.GT.1) THEN
          DO IBM=1,NMIX
            IF(MIXC(IBM).EQ.0) CYCLE
            WRITE(6,'(/40H NCRAGF: DISCONTINUITY FACTORS - MIXTURE,I5)')
     1      IBM
            DO ITYPE=1,NTYPE
              WRITE(6,'(1X,A,1H:,1P,(5X,10E12.4))') TRIM(HADF(ITYPE)),
     1        ((ADF2M(IBM,IGR,JGR,ITYPE),IGR=1,NGRP),JGR=1,NGRP)
            ENDDO
          ENDDO
        ENDIF
      ENDIF
      IF(NGFF.NE.0) THEN
        CALL LCMSIX(IPMAC,'GFF',1)
        CALL LCMPUT(IPMAC,'NWT0',NMIX*NGFF*NGRP,2,GFF(1,1,1,1))
        CALL LCMPUT(IPMAC,'H-FACTOR',NMIX*NGFF*NGRP,2,GFF(1,1,1,2))
        IF(NFINF.GT.0) THEN
          CALL LCMGET(IPMAC,'FINF_NUMBER ',FINF)
          DO I=1,NFINF
            WRITE(FINFN,'(5HFINF_,I3.3)') FINF(I)
            CALL LCMPUT(IPMAC,FINFN,NMIX*NGFF*NGRP,2,GFF(1,1,1,2+I))
          ENDDO
        ENDIF
        CALL LCMSIX(IPMAC,' ',2)
      ENDIF
      IACCS=1
*----
*  END OF OVERALL MULTICOMPO MIXTURE LOOP
*----
      IF(IMPX.GT.0) WRITE(IOUT,'(/33H NCRAGF: PROCESS MULTICOMPO MIXTU,
     1 6HRE-OUT,I5)') IBMOLD
  500 CONTINUE
*----
*  AVERAGE PHYSICAL ALBEDO INFORMATION
*----
      IF(NALBP.NE.0) THEN
        ALLOCATE(VOL(NMIX))
        CALL LCMGET(IPMAC,'VOLUME',VOL)
        DO IGR=1,NGRP
          DO IAL=1,NALBP
            GAR1=0.0D0
            GAR2=0.0D0
            DO IBM=1,NMIX
              GAR1=GAR1+ALBP2(IBM,IAL,IGR)*VOL(IBM)
              GAR2=GAR2+VOL(IBM)
            ENDDO
            ALBP(IAL,IGR)=REAL((1.0D0-GAR1/GAR2)/(1.0D0+GAR1/GAR2))
          ENDDO
        ENDDO
        DEALLOCATE(VOL)
        CALL LCMPUT(IPMAC,'ALBEDO',NALBP*NGRP,2,ALBP(1,1))
      ENDIF
*----
*  AVERAGE KINF
*----
      IF(IKINF.EQ.1) THEN
        ALLOCATE(VOL(NMIX))
        CALL LCMGET(IPMAC,'VOLUME',VOL)
        GAR1=0.0D0
        GAR2=0.0D0
        DO IBM=1,NMIX
          GAR1=GAR1+ZKINF(IBM)*VOL(IBM)
          GAR2=GAR2+VOL(IBM)
        ENDDO
        ZZZ=REAL(GAR1/GAR2)
        DEALLOCATE(VOL)
        CALL LCMPUT(IPMAC,'K-INFINITY',1,2,ZZZ)
      ENDIF
*----
*  AVERAGE KEFF
*----
      IF(IKEFF.EQ.1) THEN
        ALLOCATE(VOL(NMIX))
        CALL LCMGET(IPMAC,'VOLUME',VOL)
        GAR1=0.0D0
        GAR2=0.0D0
        DO IBM=1,NMIX
          GAR1=GAR1+ZKEFF(IBM)*VOL(IBM)
          GAR2=GAR2+VOL(IBM)
        ENDDO
        ZZZ=REAL(GAR1/GAR2)
        DEALLOCATE(VOL)
        CALL LCMPUT(IPMAC,'K-EFFECTIVE',1,2,ZZZ)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(ADF2M,ADF2,HADF,ZKEFF,ZKINF,ALBP2,ALBP,GAR5,GFF,GAR6,
     1 GAR4)
      RETURN
      END
