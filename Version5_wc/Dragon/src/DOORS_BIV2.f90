SUBROUTINE DOORS_BIV2(IPTRK,NANIS,NREG,NMAT,NUN,MATCOD,VOL,SIGG,SUNKNO,FUNKNO)
  !
  !-----------------------------------------------------------------------
  !
  !Purpose:
  ! Compute the source for the solution of diffusion or PN equations.
  ! Use a BIVAC tracking without a flat flux approximation
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
  ! IPTRK   pointer to the tracking LCM object.
  ! NANIS   maximum cross section Legendre order (=0: isotropic).
  ! NREG    number of regions.
  ! NMAT    number of mixtures.
  ! NUN     number of unknowns per energy group including net current.
  ! MATCOD  mixture indices.
  ! VOL     volumes. Volumes are included in SUNKNO.
  ! SIGG    cross section.
  ! FUNKNO  unknown vector.
  !
  !Parameters: input/output
  ! SUNKNO  integrated sources.
  !
  !-----------------------------------------------------------------------
  !
  USE GANLIB
  !----
  !  SUBROUTINE ARGUMENTS
  !----
  TYPE(C_PTR) IPTRK
  INTEGER NANIS,NREG,NMAT,NUN,MATCOD(NREG)
  REAL VOL(NREG),SIGG(0:NMAT,NANIS+1),SUNKNO(NUN),FUNKNO(NUN)
  !----
  !  LOCAL VARIABLES
  !----
  PARAMETER(NSTATE=40)
  INTEGER JPAR(NSTATE)
  !----
  !  RECOVER BIVAC SPECIFIC PARAMETERS.
  !----
  CALL LCMGET(IPTRK,'STATE-VECTOR',JPAR)
  IF(JPAR(1).NE.NREG) CALL XABORT('BIV2: INCONSISTENT NREG.')
  IF(JPAR(2).NE.NUN) CALL XABORT('BIV2: INCONSISTENT NUN.')
  ITYPE=JPAR(6)
  IELEM=JPAR(8)
  ICOL=JPAR(9)
  NLF=JPAR(14)
  IF((ITYPE.EQ.2).OR.(ITYPE.EQ.5)) THEN
    ! Cartesian 1D or 2D geometry
    IF((IELEM.GT.0).AND.(ICOL.LE.3)) THEN
      ! Raviart-Thomas / diffusion or SPN
      CALL BIV2GSO(IPTRK,NANIS,NREG,NMAT,NUN,MATCOD,VOL,SIGG,SUNKNO,FUNKNO)
    ELSE IF((IELEM.LT.0).AND.(NLF.EQ.0)) THEN
      ! Lagrange / diffusion
      CALL BIV2FSO(IPTRK,NREG,NMAT,NUN,MATCOD,VOL,SIGG,SUNKNO,FUNKNO)
    ELSE
      CALL XABORT('BIV2: DISCRETIZATION TYPE NOT AVAILABLE(1).')
    ENDIF
  ELSE IF(ITYPE.EQ.8) THEN
    ! Hexagonal 2D geometry
    IF((IELEM.GT.0).AND.(ICOL.LE.3)) THEN
      ! Raviart-Thomas / diffusion or SPN
      CALL BIV2GSO(IPTRK,NANIS,NREG,NMAT,NUN,MATCOD,VOL,SIGG,SUNKNO,FUNKNO)
    ELSE IF((IELEM.LT.0).AND.(NLF.EQ.0)) THEN
      ! Lagrange / diffusion
      CALL BIV2FSH(IPTRK,NREG,NMAT,NUN,MATCOD,VOL,SIGG,SUNKNO,FUNKNO)
    ELSE
      CALL XABORT('BIV2: DISCRETIZATION TYPE NOT AVAILABLE(2).')
    ENDIF
  ELSE
    CALL XABORT('BIV2: GEOMETRY TYPE NOT AVAILABLE.')
  ENDIF
  RETURN
CONTAINS
  SUBROUTINE BIV2FSO(IPTRK,NREG,NMAT,NUN,MATCOD,VOL,SIGG,SUNKNO,FUNKNO)
    !
    !-----------------------------------------------------------------------
    !
    !Purpose:
    ! Source term calculation for finite element or mesh corner finite
    ! differences in Cartesian geometry.
    !
    !-----------------------------------------------------------------------
    !
    USE GANLIB
    !----
    !  SUBROUTINE ARGUMENTS
    !----
    TYPE(C_PTR) IPTRK
    INTEGER NREG,NMAT,NUN,MATCOD(NREG)
    REAL VOL(NREG),SIGG(0:NMAT),SUNKNO(NUN),FUNKNO(NUN)
    !----
    !  LOCAL VARIABLES
    !----
    PARAMETER(NSTATE=40)
    INTEGER JPAR(NSTATE),IJ1(25),IJ2(25)
    LOGICAL CYLIND
    !----
    !  ALLOCATABLE ARRAYS
    !----
    INTEGER, ALLOCATABLE, DIMENSION(:) :: KN,IDL
    REAL, ALLOCATABLE, DIMENSION(:) :: XX,DD,T,TS
    REAL, ALLOCATABLE, DIMENSION(:,:) :: R,RS
    !----
    !  RECOVER BIVAC SPECIFIC PARAMETERS.
    !----
    CALL LCMGET(IPTRK,'STATE-VECTOR',JPAR)
    ITYPE=JPAR(6)
    IELEM=JPAR(8)
    ICOL=JPAR(9)
    L4=JPAR(11)
    LX=JPAR(12)
    NLF=JPAR(14)
    ISPN=JPAR(15)
    ISCAT=JPAR(16)
    CYLIND=(ITYPE.EQ.3).OR.(ITYPE.EQ.6)
    IF(IELEM.GT.0) CALL XABORT('BIV2FSO: LAGRANGE METHOD EXPECTED.')
    CALL LCMSIX(IPTRK,'BIVCOL',1)
    CALL LCMLEN(IPTRK,'T',LC,ITYLCM)
    ALLOCATE(R(LC,LC),RS(LC,LC),T(LC),TS(LC))
    CALL LCMGET(IPTRK,'R',R)
    CALL LCMGET(IPTRK,'RS',RS)
    CALL LCMGET(IPTRK,'T',T)
    CALL LCMGET(IPTRK,'TS',TS)
    CALL LCMSIX(IPTRK,' ',2)
    !
    CALL LCMLEN(IPTRK,'KN',MAXKN,ITYLCM)
    ALLOCATE(XX(NREG),DD(NREG),KN(MAXKN),IDL(NREG))
    CALL LCMGET(IPTRK,'XX',XX)
    CALL LCMGET(IPTRK,'DD',DD)
    CALL LCMGET(IPTRK,'KN',KN)
    CALL LCMGET(IPTRK,'KEYFLX',IDL)
    !----
    !  COMPUTE VECTORS IJ1 AND IJ2
    !----
    LL=LC*LC
    DO I=1,LL
      IJ1(I)=1+MOD(I-1,LC)
      IJ2(I)=1+(I-IJ1(I))/LC
    ENDDO
    !----
    !  COMPUTE THE SOURCE
    !----
    NUM1=0
    DO K=1,NREG
      IBM=MATCOD(K)
      IF(IBM.LE.0) CYCLE
      IF(VOL(K).EQ.0.0) GO TO 10
      DO I=1,LL
        IND1=KN(NUM1+I)
        IF(IND1.EQ.0) CYCLE
        I1=IJ1(I)
        I2=IJ2(I)
        DO J=1,LL
          IND2=KN(NUM1+J)
          IF(IND2.EQ.0) CYCLE
          IF(CYLIND) THEN
            RR=(R(I1,IJ1(J))+RS(I1,IJ1(J))*XX(K)/DD(K))*R(I2,IJ2(J))*VOL(K)
          ELSE
            RR=R(I1,IJ1(J))*R(I2,IJ2(J))*VOL(K)
          ENDIF
          SUNKNO(IND1)=SUNKNO(IND1)+RR*FUNKNO(IND2)*SIGG(IBM)
        ENDDO ! J
      ENDDO ! I
      10 NUM1=NUM1+LL
    ENDDO ! K
    !----
    !  APPEND THE INTEGRATED VOLUMIC SOURCES
    !----
    NUM1=0
    DO K=1,NREG
      IBM=MATCOD(K)
      IF(IBM.LE.0) CYCLE
      SUNKNO(IDL(K))=SUNKNO(IDL(K))+FUNKNO(IDL(K))*VOL(K)*SIGG(IBM)
    ENDDO
    DEALLOCATE(IDL,KN,DD,XX)
    DEALLOCATE(TS,T,RS,R)
    RETURN
  END SUBROUTINE BIV2FSO
  !
  SUBROUTINE BIV2GSO(IPTRK,NANIS,NREG,NMAT,NUN,MATCOD,VOL,SIGG,SUNKNO,FUNKNO)
    !
    !-----------------------------------------------------------------------
    !
    !Purpose:
    ! Source term calculation for a mixed-dual formulation of the finite
    ! element technique in a 2-D Cartesian or hexagonal geometry.
    !
    !-----------------------------------------------------------------------
    !
    USE GANLIB
    !----
    !  SUBROUTINE ARGUMENTS
    !----
    TYPE(C_PTR) IPTRK
    INTEGER NANIS,NREG,NMAT,NUN,MATCOD(NREG)
    REAL VOL(NREG),SIGG(0:NMAT,NANIS+1),SUNKNO(NUN),FUNKNO(NUN)
    !----
    !  LOCAL VARIABLES
    !----
    PARAMETER(NSTATE=40)
    INTEGER JPAR(NSTATE)
    LOGICAL LHEX
    !----
    !  ALLOCATABLE ARRAYS
    !----
    INTEGER, ALLOCATABLE, DIMENSION(:) :: KN,IPERT
    REAL, ALLOCATABLE, DIMENSION(:,:) :: RR
    !----
    !  RECOVER BIVAC SPECIFIC PARAMETERS.
    !----
    CALL LCMGET(IPTRK,'STATE-VECTOR',JPAR)
    ITYPE=JPAR(6)
    IELEM=JPAR(8)
    ICOL=JPAR(9)
    L4=JPAR(11)
    LX=JPAR(12)
    NLF=JPAR(14)
    ISPN=JPAR(15)
    ISCAT=JPAR(16)
    LHEX=(ITYPE.EQ.8)
    IF((IELEM.LT.0).OR.(ICOL.GT.3)) CALL XABORT('BIV2GSO: RAVIA' &
    & //'RT-THOMAS METHOD EXPECTED.')
    CALL LCMLEN(IPTRK,'KN',MAXKN,ITYLCM)
    ALLOCATE(KN(MAXKN))
    CALL LCMGET(IPTRK,'KN',KN)
    NBLOS=0
    SIDE=0.0
    IF(LHEX) THEN
      ! Raviart-Thomas-Schneider method
      NBLOS=LX/3
      ALLOCATE(IPERT(NBLOS))
      CALL LCMGET(IPTRK,'IPERT',IPERT)
      CALL LCMGET(IPTRK,'SIDE',SIDE)
    ENDIF
    !----
    !  RECOVER THE FINITE ELEMENT UNIT STIFFNESS MATRIX.
    !----
    IF(NLF.GT.0) THEN
      CALL LCMSIX(IPTRK,'BIVCOL',1)
      CALL LCMLEN(IPTRK,'T',LC,ITYLCM)
      ALLOCATE(RR(LC,LC))
      CALL LCMGET(IPTRK,'R',RR)
      CALL LCMSIX(IPTRK,' ',2)
    ENDIF
    !----
    !  COMPUTE THE SOURCE
    !----
    IF(NLF.EQ.0) THEN
      !----
      !  CARTESIAN 2D DUAL (RAVIART-THOMAS) CASE.
      !----
      IF(.NOT.LHEX) THEN
        NUM1=0
        DO IR=1,NREG
          IBM=MATCOD(IR)
          IF(IBM.LE.0) CYCLE
          IF(VOL(IR).EQ.0.0) GO TO 10
          DO I0=1,IELEM*IELEM
            IND=KN(NUM1+1)+I0-1
            SUNKNO(IND)=SUNKNO(IND)+FUNKNO(IND)*VOL(IR)*SIGG(IBM,1)
          ENDDO ! I0
          10 NUM1=NUM1+5
        ENDDO ! IR
      ELSE IF(LHEX) THEN
        TTTT=0.5*SQRT(3.0)*SIDE*SIDE
        NUM1=0
        DO KEL=1,NBLOS
          IF(IPERT(KEL).EQ.0) CYCLE
          NUM1=NUM1+1
          IBM=MATCOD((IPERT(KEL)*3-1)+1)
          IF(IBM.LE.0) CYCLE
          GARS=SIGG(IBM,1)
          DO K2=0,IELEM-1
            DO K1=0,IELEM-1
              JND1=KN(NUM1)+K2*IELEM+K1
              JND2=KN(NBLOS+NUM1)+K2*IELEM+K1
              JND3=KN(2*NBLOS+NUM1)+K2*IELEM+K1
              SUNKNO(JND1)=SUNKNO(JND1)+FUNKNO(JND1)*TTTT*GARS
              SUNKNO(JND2)=SUNKNO(JND2)+FUNKNO(JND2)*TTTT*GARS
              SUNKNO(JND3)=SUNKNO(JND3)+FUNKNO(JND3)*TTTT*GARS
            ENDDO ! K1
          ENDDO ! K2
        ENDDO ! KEL
      ENDIF
    ELSE
      !----
      !  CARTESIAN 2D SPN CASE.
      !----
      IF(LHEX) CALL XABORT('BIV2GSO: HEXAGONAL SPN IS NOT IMPLEMENTED.')
      DO IL=0,MIN(ABS(ISCAT)-1,NANIS)
        FACT=REAL(2*IL+1)
        NUM1=0
        DO IR=1,NREG
          IBM=MATCOD(IR)
          IF(IBM.LE.0) CYCLE
          IF(VOL(IR).EQ.0.0) GO TO 20
          IF(MOD(IL,2).EQ.0) THEN
            DO I0=1,IELEM*IELEM
              IND=(IL/2)*L4+KN(NUM1+1)+I0-1
              SUNKNO(IND)=SUNKNO(IND)+FACT*FUNKNO(IND)*VOL(IR)*SIGG(IBM,IL+1)
            ENDDO ! I0
          ELSE
            DO I0=1,IELEM
              DO IC=1,2
                IIC=1+(IC-1)*IELEM
                IND1=(IL/2)*L4+ABS(KN(NUM1+1+IC))+I0-1
                S1=REAL(SIGN(1,KN(NUM1+1+IC)))
                DO JC=1,2
                  JJC=1+(JC-1)*IELEM
                  IND2=(IL/2)*L4+ABS(KN(NUM1+1+JC))+I0-1
                  IF((KN(NUM1+1+IC).NE.0).AND.(KN(NUM1+1+JC).NE.0)) THEN
                    S2=REAL(SIGN(1,KN(NUM1+1+JC)))
                    AUXX=S1*S2*FACT*RR(IIC,JJC)*VOL(IR)
                    SUNKNO(IND1)=SUNKNO(IND1)-AUXX*FUNKNO(IND2)*SIGG(IBM,IL+1)
                  ENDIF
                ENDDO ! JC
              ENDDO ! IC
              DO IC=3,4
                IIC=1+(IC-3)*IELEM
                IND1=(IL/2)*L4+ABS(KN(NUM1+1+IC))+I0-1
                S1=REAL(SIGN(1,KN(NUM1+1+IC)))
                DO JC=3,4
                  JJC=1+(JC-3)*IELEM
                  IND2=(IL/2)*L4+ABS(KN(NUM1+1+JC))+I0-1
                  IF((KN(NUM1+1+IC).NE.0).AND.(KN(NUM1+1+JC).NE.0)) THEN
                    S2=REAL(SIGN(1,KN(NUM1+1+JC)))
                    AUXX=S1*S2*FACT*RR(IIC,JJC)*VOL(IR)
                    SUNKNO(IND1)=SUNKNO(IND1)-AUXX*FUNKNO(IND2)*SIGG(IBM,IL+1)
                  ENDIF
                ENDDO ! JC
              ENDDO ! IC
            ENDDO ! I0
          ENDIF
          20 NUM1=NUM1+5
        ENDDO ! IR
      ENDDO ! IL
    ENDIF
    IF(LHEX) DEALLOCATE(IPERT)
    IF(NLF.GT.0) DEALLOCATE(RR)
    DEALLOCATE(KN)
    RETURN
  END SUBROUTINE BIV2GSO
  !
  SUBROUTINE BIV2FSH(IPTRK,NREG,NMAT,NUN,MATCOD,VOL,SIGG,SUNKNO,FUNKNO)
    !
    !-----------------------------------------------------------------------
    !
    !Purpose:
    ! Source term calculation for finite element or mesh corner finite
    ! differences in hexagonal geometry.
    !
    !-----------------------------------------------------------------------
    !
    USE GANLIB
    !----
    !  SUBROUTINE ARGUMENTS
    !----
    TYPE(C_PTR) IPTRK
    INTEGER NREG,NMAT,NUN,MATCOD(NREG)
    REAL VOL(NREG),SIGG(0:NMAT),SUNKNO(NUN),FUNKNO(NUN)
    !----
    !  LOCAL VARIABLES
    !----
    PARAMETER(NSTATE=40)
    INTEGER JPAR(NSTATE)
    INTEGER ISR(6,2),ISRH(6,2),ISRT(3,2)
    REAL TH(6),RH2(6,6),RH(6,6),RT(3,3)
    !----
    !  ALLOCATABLE ARRAYS
    !----
    INTEGER, ALLOCATABLE, DIMENSION(:) :: KN,IDL
    REAL, ALLOCATABLE, DIMENSION(:) :: QFR
    !----
    !  DATA STATEMENTS
    !----
    DATA ISRH/2,1,4,5,6,3,1,4,5,6,3,2/
    DATA ISRT/1,2,3,2,3,1/
    !----
    !  RECOVER BIVAC SPECIFIC PARAMETERS.
    !----
    CALL LCMGET(IPTRK,'STATE-VECTOR',JPAR)
    ITYPE=JPAR(6)
    IELEM=JPAR(8)
    ISPLH=JPAR(10)
    L4=JPAR(11)
    LX=JPAR(12)
    NLF=JPAR(14)
    ISPN=JPAR(15)
    ISCAT=JPAR(16)
    IF(ISPLH.EQ.1) THEN
      NELEM=MAXKN/7
    ELSE
      NELEM=MAXKN/4
    ENDIF
    IF(IELEM.GT.0) CALL XABORT('BIV2FSH: LAGRANGE METHOD EXPECTED.')
    IF(NLF.GT.0) CALL XABORT('BIV2FSH: SPN NOT IMPLEMENTED.')
    CALL LCMSIX(IPTRK,'BIVCOL',1)
    CALL LCMGET(IPTRK,'RH',RH)
    CALL LCMGET(IPTRK,'RT',RT)
    CALL LCMSIX(IPTRK,' ',2)
    !
    CALL LCMLEN(IPTRK,'KN',MAXKN,ITYLCM)
    CALL LCMLEN(IPTRK,'QFR',MAXQF,ITYLCM)
    ALLOCATE(KN(MAXKN),QFR(MAXQF),IDL(NREG))
    CALL LCMGET(IPTRK,'KN',KN)
    CALL LCMGET(IPTRK,'QFR',QFR)
    CALL LCMGET(IPTRK,'KEYFLX',IDL)
    CALL LCMGET(IPTRK,'SIDE',SIDE)
    IF(ISPLH.EQ.1) THEN
      NELEM=MAXKN/7
    ELSE
      NELEM=MAXKN/4
    ENDIF
    !----
    !  RECOVER THE HEXAGONAL MASS (RH2) MATRICES.
    !----
    IF(ISPLH.EQ.1) THEN
      ! hexagonal basis
      LH=6
      DO I=1,6
      DO J=1,2
          ISR(I,J)=ISRH(I,J)
        ENDDO ! J
      ENDDO ! I
      DO I=1,6
        TH(I)=0.0
        DO J=1,6
          RH2(I,J)=RH(I,J)
          TH(I)=TH(I)+RH(I,J)
        ENDDO ! J
      ENDDO ! I
      CONST=1.5*SQRT(3.0)
      AA=SIDE
    ELSE
      ! triangular basis
      LH=3
      DO I=1,3
        DO J=1,2
          ISR(I,J)=ISRT(I,J)
        ENDDO ! J
      ENDDO ! I
      DO I=1,3
        TH(I)=0.0
        DO J=1,3
          RH2(I,J)=RT(I,J)
          TH(I)=TH(I)+RT(I,J)
        ENDDO ! J
      ENDDO ! I
      CONST=0.25*SQRT(3.0)
      AA=SIDE/REAL(ISPLH-1)
    ENDIF
    !----
    !  COMPUTE THE SOURCE
    !----
    NUM1=0
    DO K=1,NELEM
      KHEX=KN(NUM1+LH+1)
      IF(VOL(KHEX).EQ.0.0) GO TO 10
      IBM=MATCOD(KHEX)
      VOL0=QFR(NUM1+LH+1)
      GARS=SIGG(IBM)
      DO I=1,LH
        IND1=KN(NUM1+I)
        IF(IND1.EQ.0) CYCLE
        DO J=1,LH
          IND2=KN(NUM1+J)
          IF(IND2.EQ.0) CYCLE
          SUNKNO(IND1)=SUNKNO(IND1)+RH2(I,J)*FUNKNO(IND2)*VOL0*GARS
        ENDDO ! J
      ENDDO ! I
      10 NUM1=NUM1+LH+1
    ENDDO ! K
    !----
    !  APPEND THE INTEGRATED VOLUMIC SOURCES
    !----
    NUM1=0
    DO K=1,NREG
      IBM=MATCOD(K)
      IF(IBM.LE.0) CYCLE
      SUNKNO(IDL(K))=SUNKNO(IDL(K))+FUNKNO(IDL(K))*VOL(K)*SIGG(IBM)
    ENDDO
    DEALLOCATE(IDL,QFR,KN)
    RETURN
  END SUBROUTINE BIV2FSH
END SUBROUTINE DOORS_BIV2
