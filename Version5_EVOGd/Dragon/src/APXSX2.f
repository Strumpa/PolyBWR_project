*DECK APXSX2
      SUBROUTINE APXSX2(IPAPX,IPTEMP,NGRP,NL,NMAC,NISO,NMIL,IMIL,ITRANC,
     1 RECNAM,NOMMAC,TYPISO,NOMREA,IPERM,CONCES,B2)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover the cross sections of an elementary calculation and single
* mixture in the edit structure and copy them in the Apex file.
*
*Copyright:
* Copyright (C) 2025 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPAPX   pointer to the Apex file.
* IPTEMP  pointer to the edit structure.
* NGRP    number of energy groups in the Apex file.
* NL      number of Legendre orders.
* NMAC    number of macroscopic sets in the Apex file.
* NISO    number of particularized isotopes in the Apex file.
* NMIL    number of mixtures in the Apex file.
* ITRANC
* IMIL    mixture index.
* RECNAM  character identification of calculation.
* NOMMAC  names of the macroscopic sets.
* TYPISO  types of the particularized isotopes.
* NOMREA  name of the Apex reaction.
* IPERM   pointer to the particularized isotopes in the edit structure.
* CONCES  number densities of particularized isotopes.
* B2      buckling.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE hdf5_wrap
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPAPX,IPTEMP,IPERM(NISO)
      CHARACTER*80 RECNAM
      CHARACTER*4 TYPISO(NISO)
      CHARACTER*8 NOMMAC(NMAC)
      CHARACTER*12 NOMREA
      INTEGER NGRP,NL,NMAC,NISO,IMIL,ITRANC
      REAL B2,CONCES(NISO)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPEDIT,KPEDIT
      INTEGER RANK,TYPE,NBYTE,DIMSR(5)
      CHARACTER RECNAM2*80,CM*2,TEXT12*12
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: JSO,ITYPR,IPOS,NJJ,IJJ
      REAL, ALLOCATABLE, DIMENSION(:) :: WORK1D,WO1D
      REAL, ALLOCATABLE, DIMENSION(:,:) :: WORK2D,WP2D,WF2D,WO2D
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: WORK3D,WP3D,WF3D,WO3D
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: WP4D,WF4D,WO4D
*----
*  FIND NISOP, NISOF AND NISOO
*----
      NISOF=0
      NISOP=0
      NISOO=0
      IF(NISO.EQ.0) GO TO 10
      ALLOCATE(JSO(NISO))
      DO ISO=1,NISO
        IF(TYPISO(ISO).EQ.'FISS') THEN
          NISOF=NISOF+1
          JSO(ISO)=NISOF
        ELSE IF(TYPISO(ISO).EQ.'F.P.') THEN
          NISOP=NISOP+1
          JSO(ISO)=NISOP
        ELSE IF(TYPISO(ISO).EQ.'OTHE') THEN
          NISOO=NISOO+1
          JSO(ISO)=NISOO
        ENDIF
      ENDDO
      IF(NISOF.GT.0) THEN
        CALL hdf5_create_group(IPAPX,TRIM(RECNAM)//"mic/fiss/")
      ENDIF
      IF(NISOP.GT.0) THEN
        CALL hdf5_create_group(IPAPX,TRIM(RECNAM)//"mic/f.p. /")
      ENDIF
      IF(NISOO.GT.0) THEN
        CALL hdf5_create_group(IPAPX,TRIM(RECNAM)//"mic/othe/")
      ENDIF
*----
*  RECOVER DIFF AND SCAT OF PARTICULARIZED ISOTOPES
*----
      IF((NOMREA.EQ.'DIFF').OR.(NOMREA.EQ.'SCAT')) THEN
        ALLOCATE(WF3D(NGRP,NL,NISOF),WP3D(NGRP,NL,NISOP),
     1  WO3D(NGRP,NL,NISOO),WF4D(NGRP,NGRP,NL,NISOF),
     2  WP4D(NGRP,NGRP,NL,NISOP),WO4D(NGRP,NGRP,NL,NISOO))
        WF3D(:NGRP,:NL,:NISOF)=0.0
        WP3D(:NGRP,:NL,:NISOP)=0.0
        WO3D(:NGRP,:NL,:NISOO)=0.0
        WF4D(:NGRP,:NGRP,:NL,:NISOF)=0.0
        WP4D(:NGRP,:NGRP,:NL,:NISOP)=0.0
        WO4D(:NGRP,:NGRP,:NL,:NISOO)=0.0
        DO ISO=1,NISO
          IF(.NOT.C_ASSOCIATED(IPERM(ISO))) CYCLE
          DO IL=1,NL
            WRITE(CM,'(I2.2)') IL-1
            CALL LCMLEN(IPERM(ISO),'SCAT'//CM,ILONG,ITYLCM)
            IF(ILONG.EQ.0) CYCLE
            FACT=2.0*REAL(IL)-1.0
            ALLOCATE(WORK2D(NGRP,NL),WORK3D(NGRP,NGRP,NL),ITYPR(NL))
            CALL XDRLGS(IPERM(ISO),-1,0,0,NL-1,1,NGRP,WORK2D,WORK3D,
     1      ITYPR)
            IF(TYPISO(ISO).EQ.'FISS') THEN
              WF3D(:,:,JSO(ISO))=WORK2D(:,:)
              WF4D(:,:,:,JSO(ISO))=WORK3D(:,:,:)*FACT
            ELSE IF(TYPISO(ISO).EQ.'F.P.') THEN
              WP3D(:,:,JSO(ISO))=WORK2D(:,:)
              WP4D(:,:,:,JSO(ISO))=WORK3D(:,:,:)*FACT
            ELSE IF(TYPISO(ISO).EQ.'OTHE') THEN
              WO3D(:,:,JSO(ISO))=WORK2D(:,:)
              WO4D(:,:,:,JSO(ISO))=WORK3D(:,:,:)*FACT
            ENDIF
            DEALLOCATE(ITYPR,WORK3D,WORK2D)
          ENDDO
        ENDDO
        ! remove (n,2n) from 'DIFF'
        DO ISO=1,NISO
          IF(.NOT.C_ASSOCIATED(IPERM(ISO))) CYCLE
          CALL LCMLEN(IPERM(ISO),'N2N',ILONG,ITYLCM)
          IF(ILONG.GT.0) THEN
            ALLOCATE(WORK1D(NGRP))
            CALL LCMGET(IPERM(ISO),'N2N',WORK1D)
            IF(TYPISO(ISO).EQ.'FISS') THEN
              WF3D(:,1,JSO(ISO))=WF3D(:,1,JSO(ISO))-WORK1D(:)
            ELSE IF(TYPISO(ISO).EQ.'F.P.') THEN
              WP3D(:,1,JSO(ISO))=WF3D(:,1,JSO(ISO))-WORK1D(:)
            ELSE IF(TYPISO(ISO).EQ.'OTHE') THEN
              WO3D(:,1,JSO(ISO))=WF3D(:,1,JSO(ISO))-WORK1D(:)
            ENDIF
            DEALLOCATE(WORK1D)
          ENDIF
        ENDDO
        ! remove (n,3n) from 'DIFF'
        DO ISO=1,NISO
          IF(.NOT.C_ASSOCIATED(IPERM(ISO))) CYCLE
          CALL LCMLEN(IPERM(ISO),'N3N',ILONG,ITYLCM)
          IF(ILONG.GT.0) THEN
            ALLOCATE(WORK1D(NGRP))
            CALL LCMGET(IPERM(ISO),'N3N',WORK1D)
            IF(TYPISO(ISO).EQ.'FISS') THEN
              WF3D(:,1,JSO(ISO))=WF3D(:,1,JSO(ISO))-2.0*WORK1D(:)
            ELSE IF(TYPISO(ISO).EQ.'F.P.') THEN
              WP3D(:,1,JSO(ISO))=WF3D(:,1,JSO(ISO))-2.0*WORK1D(:)
            ELSE IF(TYPISO(ISO).EQ.'OTHE') THEN
              WO3D(:,1,JSO(ISO))=WF3D(:,1,JSO(ISO))-2.0*WORK1D(:)
            ENDIF
            DEALLOCATE(WORK1D)
          ENDIF
        ENDDO
        IF(NOMREA.EQ.'DIFF') THEN
          IF(NISOF.GT.0) THEN
            WRITE(RECNAM2,'(A,9Hmic/fiss/,A)') TRIM(RECNAM),'DIFF'
            CALL hdf5_write_data(IPAPX,RECNAM2,WF3D)
          ENDIF
          IF(NISOP.GT.0) THEN
            WRITE(RECNAM2,'(A,10Hmic/f.p. /,A)') TRIM(RECNAM),'DIFF'
            CALL hdf5_write_data(IPAPX,RECNAM2,WP3D)
          ENDIF
          IF(NISOO.GT.0) THEN
            WRITE(RECNAM2,'(A,9Hmic/othe/,A)') TRIM(RECNAM),'DIFF'
            CALL hdf5_write_data(IPAPX,RECNAM2,WO3D)
          ENDIF
        ELSE IF(NOMREA.EQ.'SCAT') THEN
          IF(NISOF.GT.0) THEN
            WRITE(RECNAM2,'(A,9Hmic/fiss/,A)') TRIM(RECNAM),'SCAT'
            CALL hdf5_write_data(IPAPX,RECNAM2,WF4D)
          ENDIF
          IF(NISOP.GT.0) THEN
            WRITE(RECNAM2,'(A,10Hmic/f.p. /,A)') TRIM(RECNAM),'SCAT'
            CALL hdf5_write_data(IPAPX,RECNAM2,WP4D)
          ENDIF
          IF(NISOO.GT.0) THEN
            WRITE(RECNAM2,'(A,9Hmic/othe/,A)') TRIM(RECNAM),'SCAT'
            CALL hdf5_write_data(IPAPX,RECNAM2,WO4D)
          ENDIF
        ENDIF
        DEALLOCATE(WO4D,WP4D,WF4D,WO3D,WP3D,WF3D)
        GO TO 10
      ENDIF
*----
*  RECOVER OTHER REACTIONS OF PARTICULARIZED ISOTOPES
*----
      ALLOCATE(WF2D(NGRP,NISOF),WP2D(NGRP,NISOP),WO2D(NGRP,NISOO))
      WF2D(:NGRP,:NISOF)=0.0
      WP2D(:NGRP,:NISOP)=0.0
      WO2D(:NGRP,:NISOO)=0.0
      IF(NOMREA.EQ.'ABSO') THEN
        DO ISO=1,NISO
          IF(.NOT.C_ASSOCIATED(IPERM(ISO))) CYCLE
          CALL LCMLEN(IPERM(ISO),'NTOT0',ILONG,ITYLCM)
          IF(ILONG.EQ.0) CYCLE
          ALLOCATE(WORK1D(NGRP))
          CALL LCMGET(IPERM(ISO),'NTOT0',WORK1D)
          IF(TYPISO(ISO).EQ.'FISS') THEN
            WF2D(:,JSO(ISO))=WORK1D(:)
          ELSE IF(TYPISO(ISO).EQ.'F.P.') THEN
            WP2D(:,JSO(ISO))=WORK1D(:)
          ELSE IF(TYPISO(ISO).EQ.'OTHE') THEN
            WO2D(:,JSO(ISO))=WORK1D(:)
          ENDIF
          DEALLOCATE(WORK1D)
        ENDDO
        ! remove 'DIFF' from 'TOTA'
        DO ISO=1,NISO
          IF(.NOT.C_ASSOCIATED(IPERM(ISO))) CYCLE
          IF(TYPISO(ISO).EQ.'FISS') THEN
            WRITE(RECNAM2,'(A,9Hmic/fiss/,A)') TRIM(RECNAM),'DIFF'
            CALL hdf5_info(IPAPX,RECNAM2,RANK,TYPE,NBYTE,DIMSR)
            IF(TYPE.EQ.99) CALL XABORT('APXSX2: MISSING DIFF INFO(1).')
            CALL hdf5_read_data(IPAPX,RECNAM2,WF3D)
            WF2D(:,JSO(ISO))=WF2D(:,JSO(ISO))-WF3D(:,1,JSO(ISO))
            DEALLOCATE(WF3D)
          ELSE IF(TYPISO(ISO).EQ.'F.P.') THEN
            WRITE(RECNAM2,'(A,10Hmic/f.p. /,A)') TRIM(RECNAM),'DIFF'
            CALL hdf5_info(IPAPX,RECNAM2,RANK,TYPE,NBYTE,DIMSR)
            IF(TYPE.EQ.99) CALL XABORT('APXSX2: MISSING DIFF INFO(2).')
            CALL hdf5_read_data(IPAPX,RECNAM2,WP3D)
            WP2D(:,JSO(ISO))=WP2D(:,JSO(ISO))-WP3D(:,1,JSO(ISO))
            DEALLOCATE(WP3D)
          ELSE IF(TYPISO(ISO).EQ.'OTHE') THEN
            WRITE(RECNAM2,'(A,9Hmic/othe/,A)') TRIM(RECNAM),'DIFF'
            CALL hdf5_info(IPAPX,RECNAM2,RANK,TYPE,NBYTE,DIMSR)
            IF(TYPE.EQ.99) CALL XABORT('APXSX2: MISSING DIFF INFO(3).')
            CALL hdf5_read_data(IPAPX,RECNAM2,WO3D)
            WO2D(:,JSO(ISO))=WO2D(:,JSO(ISO))-WO3D(:,1,JSO(ISO))
            DEALLOCATE(WO3D)
          ENDIF
        ENDDO
      ELSE
        IF(NOMREA.EQ.'TOTA') THEN
          TEXT12='NTOT0'
        ELSE IF(NOMREA.EQ.'TOT1') THEN
          TEXT12='NTOT1'
        ELSE IF(NOMREA.EQ.'NUFI') THEN
          TEXT12='NUSIGF'
        ELSE IF(NOMREA.EQ.'FISS') THEN
          TEXT12='NFTOT'
        ELSE IF(NOMREA.EQ.'ENER') THEN
          TEXT12='H-FACTOR'
        ELSE IF((NOMREA.EQ.'CORR').AND.(ITRANC.EQ.1).AND.(NL.GE.2)) THEN
          TEXT12='SIGS01'
        ELSE IF((NOMREA.EQ.'CORR').AND.(ITRANC.EQ.2)) THEN
          TEXT12='TRANC'
        ELSE
          TEXT12=NOMREA
        ENDIF
        DO ISO=1,NISO
          IF(.NOT.C_ASSOCIATED(IPERM(ISO))) CYCLE
          CALL LCMLEN(IPERM(ISO),TEXT12,ILONG,ITYLCM)
          IF(ILONG.EQ.0) CYCLE
          ALLOCATE(WORK1D(NGRP))
          CALL LCMGET(IPERM(ISO),TEXT12,WORK1D)
          IF(TYPISO(ISO).EQ.'FISS') THEN
            WF2D(:,JSO(ISO))=WORK1D(:)
          ELSE IF(TYPISO(ISO).EQ.'F.P.') THEN
            WP2D(:,JSO(ISO))=WORK1D(:)
          ELSE IF(TYPISO(ISO).EQ.'OTHE') THEN
            WO2D(:,JSO(ISO))=WORK1D(:)
          ENDIF
          DEALLOCATE(WORK1D)
        ENDDO
        IF(NOMREA.EQ.'ENER') THEN
          WF2D(:,:)=WF2D(:,:)*1.0E-6
          WP2D(:,:)=WP2D(:,:)*1.0E-6
          WO2D(:,:)=WO2D(:,:)*1.0E-6
        ELSE IF(NOMREA.EQ.'LEAK') THEN
          WF2D(:,:)=WF2D(:,:)*B2
          WP2D(:,:)=WP2D(:,:)*B2
          WO2D(:,:)=WO2D(:,:)*B2
        ENDIF
      ENDIF
      IF(NISOF.GT.0) THEN
        WRITE(RECNAM2,'(A,9Hmic/fiss/,A)') TRIM(RECNAM),TRIM(NOMREA)
        CALL hdf5_write_data(IPAPX,RECNAM2,WF2D)
      ENDIF
      IF(NISOP.GT.0) THEN
        WRITE(RECNAM2,'(A,10Hmic/f.p. /,A)') TRIM(RECNAM),TRIM(NOMREA)
        CALL hdf5_write_data(IPAPX,RECNAM2,WP2D)
      ENDIF
      IF(NISOO.GT.0) THEN
        WRITE(RECNAM2,'(A,9Hmic/othe/,A)') TRIM(RECNAM),TRIM(NOMREA)
        CALL hdf5_write_data(IPAPX,RECNAM2,WO2D)
      ENDIF
      DEALLOCATE(WO2D,WP2D,WF2D)
*----
*  RECOVER DIFF AND SCAT OF MACROSCOPIC SETS
*----
   10 CALL LCMSIX(IPTEMP,'MACROLIB',1)
      JPEDIT=LCMGID(IPTEMP,'GROUP')
      IF(NMAC.GT.0) THEN
        CALL hdf5_create_group(IPAPX,TRIM(RECNAM)//"mac/TOTAL/")
      ENDIF
      DO IMAC=1,NMAC
        IF(NOMMAC(IMAC).EQ.'TOTAL') THEN
          IF((NOMREA.EQ.'DIFF').OR.(NOMREA.EQ.'SCAT')) THEN
            ALLOCATE(WO2D(NGRP,NL),WO3D(NGRP,NGRP,NL))
            WO2D(:NGRP,:NL)=0.0
            WO3D(:NGRP,:NGRP,:NL)=0.0
            DO IGR=1,NGRP
              KPEDIT=LCMGIL(JPEDIT,IGR)
              ALLOCATE(IJJ(NMIL),NJJ(NMIL),IPOS(NMIL),WORK1D(NGRP*NMIL))
              DO IL=1,NL
                WRITE(CM,'(I2.2)') IL-1
                CALL LCMLEN(KPEDIT,'SCAT'//CM,ILONG,ITYLCM)
                IF(ILONG.EQ.0) CYCLE
                CALL LCMGET(KPEDIT,'IJJS'//CM,IJJ)
                CALL LCMGET(KPEDIT,'NJJS'//CM,NJJ)
                CALL LCMGET(KPEDIT,'IPOS'//CM,IPOS)
                CALL LCMGET(KPEDIT,'SCAT'//CM,WORK1D)
                IPO=IPOS(IMIL)
                J2=IJJ(IMIL)
                J1=IJJ(IMIL)-NJJ(IMIL)+1
                DO JGR=J2,J1,-1
                  WO2D(JGR,IL)=WO2D(JGR,IL)+WORK1D(IPO)
                  WO3D(IGR,JGR,IL)=WORK1D(IPO)*REAL(2*IL-1)
                  IPO=IPO+1
                ENDDO
              ENDDO ! IL
              DEALLOCATE(WORK1D,IPOS,NJJ,IJJ)
              ! remove (n,2n) from 'DIFF'
              CALL LCMLEN(KPEDIT,'N2N',ILONG,ITYLCM)
              IF(ILONG.GT.0) THEN
                ALLOCATE(WORK1D(NMIL))
                CALL LCMGET(KPEDIT,'N2N',WORK1D)
                WO2D(IGR,1)=WO2D(IGR,1)-WORK1D(IMIL)
                DEALLOCATE(WORK1D)
              ENDIF
              ! remove (n,2n) from 'DIFF'
              CALL LCMLEN(KPEDIT,'N2N',ILONG,ITYLCM)
              IF(ILONG.GT.0) THEN
                ALLOCATE(WORK1D(NMIL))
                CALL LCMGET(KPEDIT,'N2N',WORK1D)
                WO2D(IGR,1)=WO2D(IGR,1)-2.0*WORK1D(IMIL)
                DEALLOCATE(WORK1D)
              ENDIF
            ENDDO ! IGR
            IF(NOMREA.EQ.'DIFF') THEN
              WRITE(RECNAM2,'(A,10Hmac/TOTAL/,A)') TRIM(RECNAM),'DIFF'
              CALL hdf5_write_data(IPAPX,RECNAM2,WO2D)
            ELSE IF(NOMREA.EQ.'SCAT') THEN
              WRITE(RECNAM2,'(A,10Hmac/TOTAL/,A)') TRIM(RECNAM),'SCAT'
              CALL hdf5_write_data(IPAPX,RECNAM2,WO3D)
            ENDIF
            DEALLOCATE(WO3D,WO2D)
*----
*  RECOVER OTHER REACTIONS OF MACROSCOPIC SETS
*----
          ELSE IF(NOMREA.EQ.'ABSO') THEN
            ALLOCATE(WO1D(NGRP))
            WO1D(:NGRP)=0.0
            DO IGR=1,NGRP
              KPEDIT=LCMGIL(JPEDIT,IGR)
              CALL LCMLEN(KPEDIT,'NTOT0',ILONG,ITYLCM)
              IF(ILONG.EQ.0) CYCLE
              ALLOCATE(WORK1D(NMIL))
              CALL LCMGET(KPEDIT,'NTOT0',WORK1D)
              WO1D(IGR)=WORK1D(IMIL)
              DEALLOCATE(WORK1D)
            ENDDO
            ! remove 'DIFF' from 'TOTA'
            WRITE(RECNAM2,'(A,10Hmac/TOTAL/,A)') TRIM(RECNAM),'DIFF'
            CALL hdf5_info(IPAPX,RECNAM2,RANK,TYPE,NBYTE,DIMSR)
            IF(TYPE.EQ.99) CALL XABORT('APXSX2: MISSING DIFF INFO(4).')
            CALL hdf5_read_data(IPAPX,RECNAM2,WO2D)
            WO1D(:)=WO1D(:)-WO2D(:,1)
            DEALLOCATE(WO2D)
            WRITE(RECNAM2,'(A,10Hmac/TOTAL/,A)') TRIM(RECNAM),'ABSO'
            CALL hdf5_write_data(IPAPX,RECNAM2,WO1D)
            DEALLOCATE(WO1D)
          ELSE
            IF(NOMREA.EQ.'TOTA') THEN
              TEXT12='NTOT0'
            ELSE IF(NOMREA.EQ.'TOT1') THEN
              TEXT12='NTOT1'
            ELSE IF(NOMREA.EQ.'NUFI') THEN
              TEXT12='NUSIGF'
            ELSE IF(NOMREA.EQ.'FISS') THEN
              TEXT12='NFTOT'
            ELSE IF(NOMREA.EQ.'ENER') THEN
              TEXT12='H-FACTOR'
            ELSE IF(NOMREA.EQ.'LEAK') THEN
              TEXT12='DIFF'
            ELSE
              TEXT12=NOMREA
            ENDIF
            ALLOCATE(WO1D(NGRP))
            WO1D(:NGRP)=0.0
            DO IGR=1,NGRP
              KPEDIT=LCMGIL(JPEDIT,IGR)
              CALL LCMLEN(KPEDIT,TEXT12,ILONG,ITYLCM)
              IF(ILONG.EQ.0) CYCLE
              ALLOCATE(WORK1D(NMIL))
              CALL LCMGET(KPEDIT,TEXT12,WORK1D)
              WO1D(IGR)=WORK1D(IMIL)
              DEALLOCATE(WORK1D)
            ENDDO
            IF(NOMREA.EQ.'ENER') THEN
              WO1D(:)=WO1D(:)*1.0E-6
            ELSE IF(NOMREA.EQ.'LEAK') THEN
              WO1D(:)=WO1D(:)*B2
            ENDIF
            WRITE(RECNAM2,'(A,10Hmac/TOTAL/,A)') TRIM(RECNAM),
     1      TRIM(NOMREA)
            CALL hdf5_write_data(IPAPX,RECNAM2,WO1D)
            DEALLOCATE(WO1D)
          ENDIF
        ELSE IF(NOMMAC(IMAC).EQ.'RESIDUAL') THEN
          ! substract particularized contributions
          CALL hdf5_create_group(IPAPX,TRIM(RECNAM)//"mac/RESIDUAL/")
          IF(NOMREA.EQ.'DIFF') THEN
            WRITE(RECNAM2,'(A,10Hmac/TOTAL/,A)') TRIM(RECNAM),'DIFF'
            CALL hdf5_read_data(IPAPX,RECNAM2,WORK2D)
            IF(NISOF.GT.0) THEN
              WRITE(RECNAM2,'(A,9Hmic/fiss/,A)') TRIM(RECNAM),'DIFF'
              CALL hdf5_read_data(IPAPX,RECNAM2,WF3D)
            ENDIF
            IF(NISOP.GT.0) THEN
              WRITE(RECNAM2,'(A,10Hmic/f.p. /,A)') TRIM(RECNAM),'DIFF'
              CALL hdf5_read_data(IPAPX,RECNAM2,WP3D)
            ENDIF
            IF(NISOO.GT.0) THEN
              WRITE(RECNAM2,'(A,9Hmic/othe/,A)') TRIM(RECNAM),'DIFF'
              CALL hdf5_read_data(IPAPX,RECNAM2,WO3D)
            ENDIF
            DO ISO=1,NISO
              CONC=CONCES(ISO)
              IF(TYPISO(ISO).EQ.'FISS') THEN
                WORK2D(:,:)=WORK2D(:,:)-CONC*WF3D(:,:,JSO(ISO))
               ELSE IF(TYPISO(ISO).EQ.'F.P.') THEN
                WORK2D(:,:)=WORK2D(:,:)-CONC*WP3D(:,:,JSO(ISO))
              ELSE IF(TYPISO(ISO).EQ.'OTHE') THEN
                WORK2D(:,:)=WORK2D(:,:)-CONC*WO3D(:,:,JSO(ISO))
              ENDIF
            ENDDO
            IF(NISOF.GT.0) DEALLOCATE(WF3D)
            IF(NISOP.GT.0) DEALLOCATE(WP3D)
            IF(NISOO.GT.0) DEALLOCATE(WO3D)
            WRITE(RECNAM2,'(A,13Hmac/RESIDUAL/,A)') TRIM(RECNAM),'DIFF'
            CALL hdf5_write_data(IPAPX,RECNAM2,WORK2D)
            DEALLOCATE(WORK2D)
          ELSE IF(NOMREA.EQ.'SCAT') THEN
            WRITE(RECNAM2,'(A,10Hmac/TOTAL/,A)') TRIM(RECNAM),'SCAT'
            CALL hdf5_read_data(IPAPX,RECNAM2,WORK3D)
            IF(NISOF.GT.0) THEN
              WRITE(RECNAM2,'(A,9Hmic/fiss/,A)') TRIM(RECNAM),'SCAT'
              CALL hdf5_read_data(IPAPX,RECNAM2,WF4D)
            ENDIF
            IF(NISOP.GT.0) THEN
              WRITE(RECNAM2,'(A,10Hmic/f.p. /,A)') TRIM(RECNAM),'SCAT'
              CALL hdf5_read_data(IPAPX,RECNAM2,WP4D)
            ENDIF
            IF(NISOO.GT.0) THEN
              WRITE(RECNAM2,'(A,9Hmic/othe/,A)') TRIM(RECNAM),'SCAT'
              CALL hdf5_read_data(IPAPX,RECNAM2,WO4D)
            ENDIF
            DO ISO=1,NISO
              CONC=CONCES(ISO)
              IF(TYPISO(ISO).EQ.'FISS') THEN
                WORK3D(:,:,:)=WORK3D(:,:,:)-CONC*WF4D(:,:,:,JSO(ISO))
              ELSE IF(TYPISO(ISO).EQ.'F.P.') THEN
                WORK3D(:,:,:)=WORK3D(:,:,:)-CONC*WP4D(:,:,:,JSO(ISO))
              ELSE IF(TYPISO(ISO).EQ.'OTHE') THEN
                WORK3D(:,:,:)=WORK3D(:,:,:)-CONC*WO4D(:,:,:,JSO(ISO))
              ENDIF
            ENDDO
            IF(NISOF.GT.0) DEALLOCATE(WF4D)
            IF(NISOP.GT.0) DEALLOCATE(WP4D)
            IF(NISOO.GT.0) DEALLOCATE(WO4D)
            WRITE(RECNAM2,'(A,13Hmac/RESIDUAL/,A)') TRIM(RECNAM),'SCAT'
            CALL hdf5_write_data(IPAPX,RECNAM2,WORK3D)
            DEALLOCATE(WORK3D)
          ELSE
            WRITE(RECNAM2,'(A,10Hmac/TOTAL/,A)') TRIM(RECNAM),
     1      TRIM(NOMREA)
            CALL hdf5_read_data(IPAPX,RECNAM2,WORK1D)
            IF(NISOF.GT.0) THEN
              WRITE(RECNAM2,'(A,9Hmic/fiss/,A)') TRIM(RECNAM),
     1        TRIM(NOMREA)
              CALL hdf5_read_data(IPAPX,RECNAM2,WF2D)
            ENDIF
            IF(NISOP.GT.0) THEN
              WRITE(RECNAM2,'(A,10Hmic/f.p. /,A)') TRIM(RECNAM),
     1        TRIM(NOMREA)
              CALL hdf5_read_data(IPAPX,RECNAM2,WP2D)
            ENDIF
            IF(NISOO.GT.0) THEN
              WRITE(RECNAM2,'(A,9Hmic/othe/,A)') TRIM(RECNAM),
     1        TRIM(NOMREA)
              CALL hdf5_read_data(IPAPX,RECNAM2,WO2D)
            ENDIF
            DO ISO=1,NISO
              CONC=CONCES(ISO)
              IF(TYPISO(ISO).EQ.'FISS') THEN
                WORK1D(:)=WORK1D(:)-CONC*WF2D(:,JSO(ISO))
              ELSE IF(TYPISO(ISO).EQ.'F.P.') THEN
                WORK1D(:)=WORK1D(:)-CONC*WP2D(:,JSO(ISO))
              ELSE IF(TYPISO(ISO).EQ.'OTHE') THEN
                WORK1D(:)=WORK1D(:)-CONC*WO2D(:,JSO(ISO))
              ENDIF
            ENDDO
            IF(NISOF.GT.0) DEALLOCATE(WF2D)
            IF(NISOP.GT.0) DEALLOCATE(WP2D)
            IF(NISOO.GT.0) DEALLOCATE(WO2D)
            WRITE(RECNAM2,'(A,13Hmac/RESIDUAL/,A)') TRIM(RECNAM),
     1      TRIM(NOMREA)
            CALL hdf5_write_data(IPAPX,RECNAM2,WORK1D)
            DEALLOCATE(WORK1D)
          ENDIF
        ENDIF
      ENDDO
      CALL LCMSIX(IPTEMP,' ',2)
      IF(NISO.GT.0) DEALLOCATE(JSO)
      RETURN
      END
