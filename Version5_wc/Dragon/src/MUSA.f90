!
!-----------------------------------------------------------------------
!
!Purpose:
! Calculation of cellwise scattering-reduced collision, escape and
! transmission probabilities for the current iteration method in the
! multicell surfacic approximation.
!
!Copyright:
! Copyright (C) 2025 Ecole Polytechnique de Montreal
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version
!
!Author(s): A. Hebert
!
!Parameters: input
! IPSYS   pointer to the system matrices.
! IPTRK   pointer to the tracking (L_TRACK signature).
! IFTRAK  tracking file unit.
! IMPX    print flag (equal to zero for no print).
! NREG    total number of merged regions for which specific values
!         of the neutron flux and reactions rates are required.
! NBMIX   number of mixtures.
! SIGT0   total macroscopic cross sections ordered by mixture.
! SIGW0   within-group scattering macroscopic cross section ordered
!         by mixture.
! NBATCH  number of tracks dispatched in each OpenMP core.
! TITREC  title.
! NALBP   number of multigroup physical albedos.
! ALBP    multigroup physical albedos.
!
!-----------------------------------------------------------------------
!
SUBROUTINE MUSA(IPSYS,IPTRK,IFTRAK,IMPX,NREG,NBMIX,SIGT0,SIGW0,NBATCH, &
  & TITREC,NALBP,ALBP)
  USE GANLIB
  !----
  !  SUBROUTINE ARGUMENTS
  !----
  TYPE(C_PTR) IPSYS,IPTRK
  INTEGER IFTRAK,IMPX,NREG,NBMIX,NBATCH,NALBP
  REAL SIGT0(0:NBMIX),SIGW0(0:NBMIX),ALBP(NALBP)
  CHARACTER TITREC*72
  !----
  !  LOCAL VARIABLES
  !----
  PARAMETER (EPS1=1.0E-4,NSTATE=40)
  TYPE(C_PTR) JPTRK,KPTRK
  INTEGER ISTATT(NSTATE),NNPSYS(1)
  CHARACTER TITRE2*72
  logical LSKIP
  !----
  !  ALLOCATABLE ARRAYS
  !----
  INTEGER, ALLOCATABLE, DIMENSION(:) :: MATALB,NMC_NODE,NMC_SURF,MAT2,IGEN,INUM
  REAL, ALLOCATABLE, DIMENSION(:) :: SIGT2,SIGW2,WORK
  REAL, ALLOCATABLE, DIMENSION(:,:) :: VOLSUR
  REAL, POINTER, DIMENSION(:) :: PSSW,PSJW,PISW,PIJW
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: DPROB,DPROBX
  TYPE(C_PTR) :: PSSW_PTR,PSJW_PTR,PISW_PTR,PIJW_PTR
  !
  IND(I,J) = MAX(I+J3+1,J+J3+1)*(MAX(I+J3+1,J+J3+1)-1)/2 &
  & + MIN(I+J3+1,J+J3+1)
  !
  WPR(I,J)= REAL(DPROB( IND(I,J),1 ) / DPROB( IND(I,0),1 ))
  !----
  !  BICKLEY FLAG
  !----
  SAVE IBICKL
  DATA IBICKL/0/
  !----
  !  RECOVER BICKLEY TABLES
  !----
  IF(IBICKL.EQ.0) THEN
    CALL XDRTA2
    IBICKL=1
  ENDIF
  !----
  !  RECOVER SALT SPECIFIC PARAMETERS
  !----
  REWIND IFTRAK
  CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATT)
  IF(NREG.NE.ISTATT(1)) THEN
     CALL XABORT('MUSA: STATE VECTOR HAS INVALID # OF ZONES.')
  ENDIF
  NMACRO=ISTATT(24) ! NGEN
  NMCEL=NMACRO
  NMERGE=NMACRO
  NGEN=NMACRO
  ALLOCATE(IGEN(NMERGE),INUM(NMCEL))
  DO IK=1,NMERGE
    IGEN(IK)=IK
  ENDDO
  DO IK=1,NMCEL
    INUM(IK)=IK
  ENDDO
  IF(NMACRO.EQ.0) CALL XABORT('MUSA: MUST MODULE TRACKING IS MANDATORY.')
  ALLOCATE(NMC_NODE(NMACRO+1),NMC_SURF(NMACRO+1))
  JPTRK=LCMGID(IPTRK,'MACRO-TRACK')
  CALL LCMGET(IPTRK,'NMC_NODE',NMC_NODE)
  CALL LCMGET(IPTRK,'NMC_SURF',NMC_SURF)
  NMIX=NMC_SURF(NMACRO+1)
  NIFR=NMC_SURF(NMACRO+1)
  !----
  ! LOOP OVER MACRO GEOMETRIES AND COMPUTE PIJ MATRICES USING EXCELP
  !----
  J1=0
  NMIX=0
  NPIJ=0
  NPIS=0
  NPSS=0
  DO IMACRO=1,NMACRO
    J2=NMC_NODE(IMACRO+1)-NMC_NODE(IMACRO)
    J3=NMC_SURF(IMACRO+1)-NMC_SURF(IMACRO)
    J1=J1+J2
    NMIX=NMIX+J3
    NPIJ=NPIJ+J2*J2
    NPIS=NPIS+J2*J3
    NPSS=NPSS+J3*J3
  ENDDO
  IF(J1.NE.NREG) CALL XABORT('MUSA: INVALID NREG.')
  IF(NMIX.NE.NMC_SURF(NMACRO+1)) CALL XABORT('MUSA: INVALID NMIX.')
  PIJW_PTR=LCMARA(NPIJ)
  PISW_PTR=LCMARA(NPIS)
  PSJW_PTR=LCMARA(NPIS)
  PSSW_PTR=LCMARA(NPSS)
  CALL C_F_POINTER(PIJW_PTR,PIJW,(/ NPIJ /))
  CALL C_F_POINTER(PISW_PTR,PISW,(/ NPIS /))
  CALL C_F_POINTER(PSJW_PTR,PSJW,(/ NPIS /))
  CALL C_F_POINTER(PSSW_PTR,PSSW,(/ NPSS /))
  J1=0
  IPIJ=0
  IPIS=0
  IPSS=0
  DO IMACRO=1,NMACRO
    J2=NMC_NODE(IMACRO+1)-NMC_NODE(IMACRO)
    J3=NMC_SURF(IMACRO+1)-NMC_SURF(IMACRO)
    N2PRO=(J2+J3+1)**2
    WRITE(TITRE2,'(A,9H -- MACRO,I5.5)') TRIM(TITREC),IMACRO
    KPTRK=LCMGIL(JPTRK,IMACRO)
    KNORM=4 ! use HELIOS-type normalization
    NNPSYS(1)=1
    ALLOCATE(MAT2(J2),SIGT2(J2),SIGW2(J2))
    ALLOCATE(MATALB(-J3:J2),VOLSUR(-J3:J2,1),DPROB(N2PRO,1),DPROBX(N2PRO,1))
    CALL LCMGET(KPTRK,'MATCOD',MAT2)
    CALL EXCELP(KPTRK,IFTRAK,IMPX,J3,J2,NBMIX,MAT2,KNORM,SIGT0,1,N2PRO, &
    & 1,NNPSYS(1),NBATCH,TITRE2,NALBP,ALBP,MATALB,VOLSUR,DPROB,DPROBX)
    !----
    !  CHECK IF SCATTERING REDUCTION IS REQUIRED 
    !----
    DO I=1,J2
      SIGT2(I)=SIGT0(MAT2(I)) ! sigt by node
      SIGW2(I)=SIGW0(MAT2(I)) ! sigw by node
    ENDDO
    LSKIP=.TRUE.
    DO I=1,J2
      LSKIP=LSKIP.AND.(SIGW2(I).EQ.0.0)
    ENDDO
    !----
    !  SCATTERING REDUCTION IF LSKIP=.FALSE.
    !----
    IF(LSKIP) THEN
      ! DO NOT PERFORM SCATTERING REDUCTION.
      DO I=1,J2
        DO J=1,J2
          IF(SIGT2(J).EQ.0.0) THEN
            PIJW(IPIJ+(J-1)*J2+I)=WPR(I,J)
          ELSE
            PIJW(IPIJ+(J-1)*J2+I)=WPR(I,J)/SIGT2(J)
          ENDIF
        ENDDO
      ENDDO
      DO I=1,J2
        DO JC=1,J3
          PISW(IPIS+(JC-1)*J2+I)=WPR(I,-JC)
          IF(SIGT2(I).EQ.0.0) THEN
            PSJW(IPIS+(I-1)*J3+JC)=WPR(-JC,I)
          ELSE
            PSJW(IPIS+(I-1)*J3+JC)=WPR(-JC,I)/SIGT2(I)
          ENDIF
        ENDDO
      ENDDO
      DO IC=1,J3
        DO JC=1,J3
          PSSW(IPSS+(JC-1)*J3+IC)=WPR(-IC,-JC)
        ENDDO
      ENDDO
    ELSE
      ! COMPUTE THE SCATTERING-REDUCED COLLISION AND ESCAPE MATRICES.
      DO I=1,J2
        DO J=1,J2
          IF(SIGT2(J).EQ.0.0) THEN
            PIJW(IPIJ+(J-1)*J2+I)=0.0
          ELSE
            PIJW(IPIJ+(J-1)*J2+I)=-WPR(I,J)*SIGW2(J)/SIGT2(J)
          ENDIF
        ENDDO
        PIJW(IPIJ+(I-1)*J2+I)=1.0+PIJW(IPIJ+(I-1)*J2+I)
      ENDDO
      CALL ALINV(J2,PIJW(IPIJ+1),J2,IER)
      IF(IER.NE.0) CALL XABORT('MUSA: SINGULAR MATRIX.')
      ALLOCATE(WORK(J2))
      DO I=1,J2
        DO K=1,J2
          WORK(K)=PIJW(IPIJ+(K-1)*J2+I)
        ENDDO
        DO J=1,J2
          WGAR=0.0
          DO K=1,J2
            IF(SIGT2(J).EQ.0.0) THEN
              WGAR=WGAR+WORK(K)*WPR(K,J)
            ELSE
              WGAR=WGAR+WORK(K)*WPR(K,J)/SIGT2(J)
            ENDIF
          ENDDO
          PIJW(IPIJ+(J-1)*J2+I)=WGAR
        ENDDO
        DO JC=1,J3
          WGAR=0.0
          DO K=1,J2
            WGAR=WGAR+WORK(K)*WPR(K,-JC)
          ENDDO
          PISW(IPIS+(JC-1)*J2+I)=WGAR
        ENDDO
      ENDDO
      DEALLOCATE(WORK)
      !
      ! COMPUTE THE SCATTERING-REDUCED COLLISION PROBABILITY MATRIX
      ! FOR INCOMING NEUTRONS.
      DO IC=1,J3
        DO J=1,J2
          IF(SIGT2(J).EQ.0.0) THEN
            WGAR=WPR(-IC,J)
          ELSE
            WGAR=WPR(-IC,J)/SIGT2(J)
          ENDIF
          DO K=1,J2
            IF(SIGT2(K).NE.0.0) THEN
              WGAR=WGAR+WPR(-IC,K)*PIJW(IPIJ+(J-1)*J2+K)*SIGW2(K)/SIGT2(K)
            ENDIF
          ENDDO
          PSJW(IPIS+(J-1)*J3+IC)=WGAR
        ENDDO
      ENDDO
      !
      ! COMPUTE THE SCATTERING-REDUCED TRANSMISSION PROBABILITY MATRIX.
      DO IC=1,J3
        DO JC=1,J3
          WGAR=WPR(-IC,-JC)
          DO K=1,J2
            IF(SIGT2(K).NE.0.0) THEN
              WGAR=WGAR+WPR(-IC,K)*PISW(IPIS+(JC-1)*J2+K)*SIGW2(K)/SIGT2(K)
            ENDIF
          ENDDO
          PSSW(IPSS+(JC-1)*J3+IC)=WGAR
        ENDDO
      ENDDO
    ENDIF
    DEALLOCATE(DPROBX,DPROB,VOLSUR,MATALB)
    IF(IMPX.GE.8) THEN
      IF(LSKIP) THEN
        IN=1
      ELSE
        IN=2
      ENDIF
      CALL SYBPRX(IN,J3,J2,IMACRO,SIGT2,SIGW2,PIJW(IPIJ+1),PISW(IPIS+1), &
      & PSJW(IPIS+1),PSSW(IPSS+1))
    ENDIF
    DEALLOCATE(SIGW2,SIGT2,MAT2)
    J1=J1+J2
    IPIJ=IPIJ+J2*J2
    IPIS=IPIS+J2*J3
    IPSS=IPSS+J3*J3
  ENDDO
  ! end of SYB004 equivalent
  CALL LCMPPD(IPSYS,'PSSW$SYBIL',NPSS,2,PSSW_PTR)
  CALL LCMPPD(IPSYS,'PSJW$SYBIL',NPIS,2,PSJW_PTR)
  CALL LCMPPD(IPSYS,'PISW$SYBIL',NPIS,2,PISW_PTR)
  CALL LCMPPD(IPSYS,'PIJW$SYBIL',NPIJ,2,PIJW_PTR)
  IF(IMPX.GT.1) THEN
    WRITE(6,'(/31H MUSA: PIJ INFORMATION IN GROUP)')
    CALL LCMLIB(IPSYS)
  ENDIF
  RETURN
END SUBROUTINE MUSA
