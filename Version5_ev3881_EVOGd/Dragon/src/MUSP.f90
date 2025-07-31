!
!-----------------------------------------------------------------------
!
!Purpose:
! Calculation of the collision probabilities for the multicell
! surfacic approximation.
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
! IPTRK   pointer to the tracking (L_TRACK signature).
! IFTRAK  tracking file unit.
! IMPX    print flag (equal to zero for no print).
! NREGIO  total number of merged blocks for which specific values
!         of the neutron flux and reactions rates are required.
! NBMIX   number of mixtures (NBMIX=max(MAT(i))).
! MAT     index-number of the mixture type assigned to each volume.
! VOL     volumes.
! SIGT0   total macroscopic cross sections ordered by mixture.
! SIGW0   P0 within-group scattering macroscopic cross sections
!         ordered by mixture.
! NELPIJ  number of elements in pij matrix.
! ILK     leakage flag (=.true. if neutron leakage through external
!         boundary is present).
! NBATCH  number of tracks dispached in eack OpenMP core.
! TITREC  title.
!
!Parameters: output
! PIJ     reduced and symmetrized collision probabilities.
!
!-----------------------------------------------------------------------
!
SUBROUTINE MUSP(IPTRK,IFTRAK,IMPX,NREGIO,NBMIX,MAT,VOL,SIGT0,SIGW0,NELPIJ, &
  & ILK,NBATCH,TITREC,PIJ)
  USE GANLIB
  !----
  !  SUBROUTINE ARGUMENTS
  !----
  LOGICAL ILK
  TYPE(C_PTR) IPTRK
  INTEGER IFTRAK,IMPX,NREGIO,NBMIX,MAT(NREGIO),NELPIJ,NBATCH
  REAL VOL(NREGIO),SIGT0(0:NBMIX),SIGW0(0:NBMIX),PIJ(NELPIJ)
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
  INTEGER, ALLOCATABLE, DIMENSION(:) :: MATALB,NMC_NODE,NMC_SURF,MAT2,IGEN,INUM,IFR,MIX
  REAL, ALLOCATABLE, DIMENSION(:) :: SIGT2,SIGW2,PIJW,PISW,PSJW,PSSW,WORK,ALB,DVX
  REAL, ALLOCATABLE, DIMENSION(:,:) :: VOLSUR,PP,PSSB
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: DPROB,DPROBX
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
  !  RECOVER MUST SPECIFIC PARAMETERS
  !----
  CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATT)
  IF(NREGIO.NE.ISTATT(1)) THEN
     CALL XABORT('MUSP: STATE VECTOR HAS INVALID # OF ZONES.')
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
  IF(NMACRO.EQ.0) CALL XABORT('MUSP: MUST MODULE TRACKING IS MANDATORY.')
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
  IJAT=0
  NPIJ=0
  NPIS=0
  NPSS=0
  DO IMACRO=1,NMACRO
    J2=NMC_NODE(IMACRO+1)-NMC_NODE(IMACRO)
    J3=NMC_SURF(IMACRO+1)-NMC_SURF(IMACRO)
    J1=J1+J2
    IJAT=IJAT+J3
    NPIJ=NPIJ+J2*J2
    NPIS=NPIS+J2*J3
    NPSS=NPSS+J3*J3
  ENDDO
  IF(J1.NE.NREGIO) CALL XABORT('MUSP: INVALID NREGIO.')
  IF(IJAT.NE.NMC_SURF(NMACRO+1)) CALL XABORT('MUSP: INVALID IJAT.')
  ALLOCATE(SIGT2(NREGIO),SIGW2(NREGIO))
  DO IREGIO=1,NREGIO
    SIGT2(IREGIO)=SIGT0(MAT(IREGIO)) ! sigt by node
    SIGW2(IREGIO)=SIGW0(MAT(IREGIO)) ! sigw by node
  ENDDO
  ALLOCATE(PIJW(NPIJ),PISW(NPIS),PSJW(NPIS),PSSW(NPSS))
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
    ALLOCATE(MAT2(J2),MATALB(-J3:J2),VOLSUR(-J3:J2,1),DPROB(N2PRO,1), &
    & DPROBX(N2PRO,1))
    CALL LCMGET(KPTRK,'MATCOD',MAT2)
    CALL EXCELP(KPTRK,IFTRAK,IMPX,J3,J2,NBMIX,MAT2,KNORM,SIGT0,1,N2PRO, &
    & 1,NNPSYS(1),NBATCH,TITRE2,NALBP,ALBP,MATALB,VOLSUR,DPROB,DPROBX)
    !----
    !  CHECK IF SCATTERING REDUCTION IS REQUIRED 
    !----
    LSKIP=.TRUE.
    DO I=1,J2
      LSKIP=LSKIP.AND.(SIGW2(J1+I).EQ.0.0)
    ENDDO
    !----
    !  SCATTERING REDUCTION IF LSKIP=.FALSE.
    !----
    IF(LSKIP) THEN
      ! DO NOT PERFORM SCATTERING REDUCTION.
      DO I=1,J2
        DO J=1,J2
          IF(SIGT2(J1+J).EQ.0.0) THEN
            PIJW(IPIJ+(J-1)*J2+I)=WPR(I,J)
          ELSE
            PIJW(IPIJ+(J-1)*J2+I)=WPR(I,J)/SIGT2(J1+J)
          ENDIF
        ENDDO
      ENDDO
      DO I=1,J2
        DO JC=1,J3
          PISW(IPIS+(JC-1)*J2+I)=WPR(I,-JC)
          IF(SIGT2(J1+I).EQ.0.0) THEN
            PSJW(IPIS+(I-1)*J3+JC)=WPR(-JC,I)
          ELSE
            PSJW(IPIS+(I-1)*J3+JC)=WPR(-JC,I)/SIGT2(J1+I)
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
          IF(SIGT2(J1+J).EQ.0.0) THEN
            PIJW(IPIJ+(J-1)*J2+I)=0.0
          ELSE
            PIJW(IPIJ+(J-1)*J2+I)=-WPR(I,J)*SIGW2(J1+J)/SIGT2(J1+J)
          ENDIF
        ENDDO
        PIJW(IPIJ+(I-1)*J2+I)=1.0+PIJW(IPIJ+(I-1)*J2+I)
      ENDDO
      CALL ALINV(J2,PIJW(IPIJ+1),J2,IER)
      IF(IER.NE.0) CALL XABORT('MUSP: SINGULAR MATRIX.')
      ALLOCATE(WORK(J2))
      DO I=1,J2
        DO K=1,J2
          WORK(K)=PIJW(IPIJ+(K-1)*J2+I)
        ENDDO
        DO J=1,J2
          WGAR=0.0
          DO K=1,J2
            IF(SIGT2(J1+J).EQ.0.0) THEN
              WGAR=WGAR+WORK(K)*WPR(K,J)
            ELSE
              WGAR=WGAR+WORK(K)*WPR(K,J)/SIGT2(J1+J)
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
          IF(SIGT2(J1+J).EQ.0.0) THEN
            WGAR=WPR(-IC,J)
          ELSE
            WGAR=WPR(-IC,J)/SIGT2(J1+J)
          ENDIF
          DO K=1,J2
            IF(SIGT2(J1+K).NE.0.0) THEN
              WGAR=WGAR+WPR(-IC,K)*PIJW(IPIJ+(J-1)*J2+K)*SIGW2(J1+K)/SIGT2(J1+K)
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
            IF(SIGT2(J1+K).NE.0.0) THEN
              WGAR=WGAR+WPR(-IC,K)*PISW(IPIS+(JC-1)*J2+K)*SIGW2(J1+K)/SIGT2(J1+K)
            ENDIF
          ENDDO
          PSSW(IPSS+(JC-1)*J3+IC)=WGAR
        ENDDO
      ENDDO
    ENDIF
    DEALLOCATE(DPROBX,DPROB,VOLSUR,MATALB,MAT2)
    IF(IMPX.GE.8) THEN
      IF(LSKIP) THEN
        IN=1
      ELSE
        IN=2
      ENDIF
      CALL SYBPRX(IN,J3,J2,IMACRO,SIGT2(J1+1),SIGW2(J1+1),PIJW(IPIJ+1), &
      & PISW(IPIS+1),PSJW(IPIS+1),PSSW(IPSS+1))
    ENDIF
    J1=J1+J2
    IPIJ=IPIJ+J2*J2
    IPIS=IPIS+J2*J3
    IPSS=IPSS+J3*J3
  ENDDO
  ! end of SYB004 equivalent
  !----
  ! COMPUTE THE GLOBAL SCATTERING-REDUCED COLLISION PROBABILITY MATRIX
  !----
  ALLOCATE(PP(NREGIO,NREGIO))
  PP(:NREGIO,:NREGIO)=0.0
  IPIJ=0
  DO JKG=1,NGEN
    J2=NMC_NODE(JKG+1)-NMC_NODE(JKG)
    I1=0
    DO IKK=1,NMERGE
      IKG=IGEN(IKK)
      I2=NMC_NODE(IKG+1)-NMC_NODE(IKG)
      IF(IKG.EQ.JKG) THEN
        DO J=1,J2
          DO I=1,J2
            PP(I1+I,I1+J)=PIJW(IPIJ+(J-1)*J2+I)
          ENDDO
        ENDDO
      ENDIF
      I1=I1+I2
    ENDDO
    IPIJ=IPIJ+J2*J2
  ENDDO
  !----
  ! COMPUTE PSSB=A*(I-PSS*A)**-1
  !----
  ALLOCATE(PSSB(IJAT,2*IJAT),IFR(NIFR),ALB(NIFR),MIX(NMIX),DVX(NMIX))
  CALL LCMGET(IPTRK,'IFR',IFR)
  CALL LCMGET(IPTRK,'ALB',ALB)
  CALL LCMGET(IPTRK,'MIX',MIX)
  CALL LCMGET(IPTRK,'DVX',DVX)
  PSSB(:IJAT,:2*IJAT)=0.0
  DO I=1,IJAT
    PSSB(I,I)=1.0
  ENDDO
  DO ICEL=1,NMCEL
    IKK=INUM(ICEL)
    IKG=IGEN(IKK)
    J3=NMC_SURF(IKG+1)-NMC_SURF(IKG)
    IT=0
    DO IK=1,IKK-1
      IT=IT+(NMC_SURF(IGEN(IK)+1)-NMC_SURF(IGEN(IK)))
    ENDDO
    IS=0
    DO IK=1,ICEL-1
      IS=IS+(NMC_SURF(IGEN(INUM(IK))+1)-NMC_SURF(IGEN(INUM(IK))))
    ENDDO
    IPSS=0
    DO IK=1,IKG-1
      IPSS=IPSS+(NMC_SURF(IK+1)-NMC_SURF(IK))**2
    ENDDO
    DO JC=1,J3
      J1=IFR(IS+JC)
      J2=MIX(IT+JC)
      ALBEDO=ALB(IS+JC)
      PSSB(J1,IJAT+J2)=PSSB(J1,IJAT+J2)+ALBEDO*DVX(IT+JC)
      DO IC=1,J3
        J2=MIX(IT+IC)
        PSSB(J1,J2)=PSSB(J1,J2)-PSSW(IPSS+(JC-1)*J3+IC)*ALBEDO*DVX(IT+IC)
      ENDDO
    ENDDO
  ENDDO
  CALL ALSB(IJAT,IJAT,PSSB,IER,IJAT)
  IF(IER.NE.0) CALL XABORT('MUSP: SINGULAR MATRIX.')
  !----
  !  COMPUTATION OF PISW*PSSB*PSJW
  !----
  I1=0
  DO IKK=1,NMERGE
    IKG=IGEN(IKK)
    I2=NMC_NODE(IKG+1)-NMC_NODE(IKG)
    I3=NMC_SURF(IKG+1)-NMC_SURF(IKG)
    IT=0
    DO IK=1,IKK-1
      IT=IT+(NMC_SURF(IGEN(IK)+1)-NMC_SURF(IGEN(IK)))
    ENDDO
    IPIS=0
    DO IK=1,IKG-1
      IPIS=IPIS+(NMC_NODE(IK+1)-NMC_NODE(IK))*(NMC_SURF(IK+1)-NMC_SURF(IK))
    ENDDO
    DO I=1,I2
      DO IC=1,I3
        ICC=MIX(IT+IC)
        ZZZ=PISW(IPIS+(IC-1)*I2+I)*SIGN(1.0,DVX(IT+IC))
        J1=0
        DO JKK=1,NMERGE
          JKG=IGEN(JKK)
          J2=NMC_NODE(JKG+1)-NMC_NODE(JKG)
          J3=NMC_SURF(JKG+1)-NMC_SURF(JKG)
          JT=0
          DO IK=1,JKK-1
            JT=JT+(NMC_SURF(IGEN(IK)+1)-NMC_SURF(IGEN(IK)))
          ENDDO
          IPSJ=0
          DO IK=1,JKG-1
            IPSJ=IPSJ+(NMC_NODE(IK+1)-NMC_NODE(IK))*(NMC_SURF(IK+1)-NMC_SURF(IK))
          ENDDO
          DO J=1,J2
            DO JC=1,J3
              JCC=MIX(JT+JC)
              PBJ=PSJW(IPSJ+(J-1)*J3+JC)
              PP(I1+I,J1+J)=PP(I1+I,J1+J)+ZZZ*DVX(JT+JC)*PSSB(JCC,IJAT+ICC)*PBJ
            ENDDO
          ENDDO
          J1=J1+J2
        ENDDO
      ENDDO
    ENDDO
    I1=I1+I2
  ENDDO
  ! end of SYBRX3 equivalent
  DEALLOCATE(INUM,IGEN,DVX,MIX,ALB,IFR)
  DEALLOCATE(PSSW,PSJW,PISW,PIJW,PSSB)
  DEALLOCATE(NMC_SURF,NMC_NODE)
  !
  IF(IMPX.GE.7) THEN
    WRITE (6,170) (J,J=1,NREGIO)
    DO I=1,NREGIO
      WRITE (6,180) I,(PP(I,J),J=1,NREGIO)
    ENDDO
    WRITE (6,'(//)')
  ENDIF
  IF((IMPX.GE.-10).OR.(IMPX.LT.0)) THEN
    ! CHECK THE RECIPROCITY CONDITIONS.
    VOLTOT=0.0
    DO I=1,NREGIO
      VOLTOT=VOLTOT+VOL(I)
    ENDDO
    VOLTOT=VOLTOT/REAL(NREGIO)
    WRK=0.0
    DO I=1,NREGIO
      DO J=1,NREGIO
        AAA=PP(I,J)*VOL(I)
        BBB=PP(J,I)*VOL(J)
        WRK=MAX(WRK,ABS(AAA-BBB)/VOLTOT)
      ENDDO
    ENDDO
    IF(WRK.GE.EPS1) WRITE (6,150) WRK
    ! CHECK THE CONSERVATION CONDITIONS.
    IF(.NOT.ILK) THEN
      WRK=0.0
      DO I=1,NREGIO
        F1=1.0
        DO J=1,NREGIO
          AAA=PP(I,J)
          F1=F1-AAA*(SIGT2(J)-SIGW2(J))
        ENDDO
        WRK=AMAX1(WRK,ABS(F1))
      ENDDO
      IF(WRK.GE.EPS1) WRITE (6,160) WRK
    ENDIF
  ENDIF
  !
  IC=0
  DO IKK=1,NREGIO
    IOF=(IKK-1)*NREGIO
    DO JKK=1,IKK
      IC=IC+1
      PIJ(IC)=PP(JKK,IKK)*VOL(JKK)
    ENDDO
  ENDDO
  DEALLOCATE(PP,SIGT2,SIGW2)
  RETURN
  !
  150 FORMAT (/56H MUSP: THE SCATTERING-REDUCED PIJ DO NOT MEET THE RECIPR, &
  & 25HOCITY CONDITIONS. RECIP =,1P,E10.3/)
  160 FORMAT (/56H MUSP: THE SCATTERING-REDUCED PIJ DO NOT MEET THE CONSER, &
  & 25HVATION CONDITIONS. LEAK =,1P,E10.3/)
  170 FORMAT (//47H MUSP: SCATTERING-REDUCED COLLISION PROBABILITY, &
  & 9H MATRIX ://(11X,2HJ=,I4,:,5X,2HJ=,I4,:,5X,2HJ=,I4,:,5X,2HJ=, &
  & I4,:,5X,2HJ=,I4,:,5X,2HJ=,I4,:,5X,2HJ=,I4,:,5X,2HJ=,I4,:,5X,2HJ=, &
  & I4,:,5X,2HJ=,I4,:,5X,2HJ=,I4))
  180 FORMAT (3H I=,I4,2H: ,1P,11E11.3/(9X,11E11.3))
END SUBROUTINE MUSP
