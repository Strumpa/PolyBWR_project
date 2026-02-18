SUBROUTINE DOORS_TRI2(IPTRK,NANIS,NREG,NMAT,NUN,MATCOD,VOL,SIGG,SUNKNO,FUNKNO)
  !
  !-----------------------------------------------------------------------
  !
  !Purpose:
  ! Compute the source for the solution of diffusion or PN equations.
  ! Use a TRIVAC tracking without a flat flux approximation
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
  IF(JPAR(1).NE.NREG) CALL XABORT('DOORS_TRI2: INCONSISTENT NREG.')
  IF(JPAR(2).NE.NUN) CALL XABORT('DOORS_TRI2: INCONSISTENT NUN.')
  IF(NANIS.NE.0) CALL XABORT('DOORS_TRI2: SPN NOT IMPLEMENTED.')
  ITYPE=JPAR(6)
  IF((ITYPE.EQ.2).OR.(ITYPE.EQ.5).OR.(ITYPE.EQ.7)) THEN
    ! Cartesian 1D, 2D or 3D geometry
    CALL TRI2GSO(IPTRK,NREG,NMAT,NUN,MATCOD,VOL,SIGG,SUNKNO,FUNKNO)
  ELSE IF((ITYPE.EQ.8).OR.(ITYPE.EQ.9))  THEN
    ! Hexagonal 2D or 3D geometry
    CALL TRI2GSR(IPTRK,NREG,NMAT,NUN,MATCOD,SIGG,SUNKNO,FUNKNO)
  ELSE
    CALL XABORT('DOORS_TRI2: GEOMETRY TYPE NOT AVAILABLE.')
  ENDIF
  RETURN
CONTAINS
  SUBROUTINE TRI2GSO(IPTRK,NREG,NMAT,NUN,MATCOD,VOL,SIGG,SUNKNO,FUNKNO)
    !
    !-----------------------------------------------------------------------
    !
    !Purpose:
    ! Source term calculation for a mixed-dual formulation of the finite
    ! element technique in a 3-D Cartesian geometry.
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
    INTEGER IPAR(NSTATE)
    CHARACTER HSMG*131
    !----
    !  ALLOCATABLE ARRAYS
    !----
    INTEGER, ALLOCATABLE, DIMENSION(:) :: KN
    !----
    !  RECOVER TRIVAC SPECIFIC PARAMETERS.
    !----
    CALL LCMGET(IPTRK,'STATE-VECTOR',IPAR)
    ITYPE=IPAR(6)
    IELEM=IPAR(9)
    ICOL=IPAR(10)
    LX=IPAR(14)
    LY=IPAR(15)
    LZ=IPAR(16)
    ISCAT=IPAR(32)
    IF((ITYPE.NE.2).AND.(ITYPE.NE.5).AND.(ITYPE.NE.7)) THEN
      CALL XABORT('TRI2GSO: INVALID CARTESIAN GEOMETRY.')
    ELSE IF((IELEM.LT.0).OR.(ICOL.GT.3)) THEN
      CALL XABORT('TRI2GSO: RAVIART-THOMAS METHOD EXPECTED.')
    ELSE IF(ISCAT.GT.1) THEN
      WRITE(HSMG,'(54HTRI2GSO: MACRO-CALCULATION WITH ANISOTROPIC SCATTERING, &
      & 63H CURRENTLY NOT IMPLEMENTED; USE SCAT 1 KEYWORD IN TRIVAT: DATA.)')
      CALL XABORT(HSMG)
    ELSE IF(LX*LY*LZ.NE.NREG) THEN
      CALL XABORT('TRI2GSO: INVALID NREG.')
    ENDIF
    CALL LCMLEN(IPTRK,'KN',MAXKN,ITYLCM)
    ALLOCATE(KN(MAXKN))
    CALL LCMGET(IPTRK,'KN',KN)
    NUM1=0
    DO K=1,NREG
      L=MATCOD(K)
      IF(L.LE.0) CYCLE
      IF(VOL(K).EQ.0.0) GO TO 10
      GARS=SIGG(L)
      DO I0=1,IELEM**3
        IND1=KN(NUM1+1)+I0-1
        SUNKNO(IND1)=SUNKNO(IND1)+FUNKNO(IND1)*VOL(K)*GARS
      ENDDO ! I0
      10 NUM1=NUM1+1+6*IELEM**2
    ENDDO ! K
    DEALLOCATE(KN)
    RETURN
  END SUBROUTINE TRI2GSO
  !
  SUBROUTINE TRI2GSR(IPTRK,NREG,NMAT,NUN,MATCOD,SIGG,SUNKNO,FUNKNO)
    !
    !-----------------------------------------------------------------------
    !
    !Purpose:
    ! Source term calculation for a Thomas-Raviart-Schneider formulation of
    ! the finite element technique in a 3-D hexagonal geometry.
    !
    !-----------------------------------------------------------------------
    !
    USE GANLIB
    !----
    !  SUBROUTINE ARGUMENTS
    !----
    TYPE(C_PTR) IPTRK
    INTEGER NREG,NMAT,NUN,MATCOD(3,NREG/3)
    REAL SIGG(0:NMAT),SUNKNO(NUN),FUNKNO(NUN)
    !----
    !  LOCAL VARIABLES
    !----
    PARAMETER(NSTATE=40)
    INTEGER IPAR(NSTATE)
    CHARACTER HSMG*131
    !----
    !  ALLOCATABLE ARRAYS
    !----
    INTEGER, ALLOCATABLE, DIMENSION(:) :: IPERT
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: KN
    REAL, ALLOCATABLE, DIMENSION(:) :: FRZ
    REAL, ALLOCATABLE, DIMENSION(:,:) :: ZZ
    !----
    !  RECOVER TRIVAC SPECIFIC PARAMETERS.
    !----
    CALL LCMGET(IPTRK,'STATE-VECTOR',IPAR)
    ITYPE=IPAR(6)
    IELEM=IPAR(9)
    ICOL=IPAR(10)
    LX=IPAR(14)
    LZ=IPAR(16)
    ISCAT=IPAR(32)
    NBLOS=LX*LZ/3
    IF((ITYPE.NE.8).AND.(ITYPE.NE.9)) THEN
      CALL XABORT('TRI2GSR: INVALID HEXAGONAL GEOMETRY.')
    ELSE IF((IELEM.LT.0).OR.(ICOL.GT.3)) THEN
      CALL XABORT('TRI2GSR: RAVIART-THOMAS METHOD EXPECTED.')
    ELSE IF(ISCAT.GT.1) THEN
      WRITE(HSMG,'(54HTRI2GSR: MACRO-CALCULATION WITH ANISOTROPIC SCATTERING, &
      & 63H CURRENTLY NOT IMPLEMENTED; USE SCAT 1 KEYWORD IN TRIVAT: DATA.)')
      CALL XABORT(HSMG)
    ELSE IF(3*NBLOS.NE.NREG) THEN
      CALL XABORT('TRI2GSR: INVALID NREG.')
    ENDIF
    CALL LCMLEN(IPTRK,'KN',MAXKN,ITYLCM)
    ALLOCATE(ZZ(3,NBLOS),KN(NBLOS,MAXKN/NBLOS),IPERT(NBLOS),FRZ(NBLOS))
    CALL LCMGET(IPTRK,'SIDE',SIDE)
    CALL LCMGET(IPTRK,'ZZ',ZZ)
    CALL LCMGET(IPTRK,'KN',KN)
    CALL LCMGET(IPTRK,'IPERT',IPERT)
    CALL LCMGET(IPTRK,'FRZ',FRZ)
    NELEM=IELEM*(IELEM+1)
    TTTT=0.5*SQRT(3.0)*SIDE*SIDE
    NUM=0
    DO KEL=1,NBLOS
      IF(IPERT(KEL).EQ.0) CYCLE
      NUM=NUM+1
      L=MATCOD(1,IPERT(KEL))
      IF(L.EQ.0) CYCLE
      VOL0=TTTT*ZZ(1,IPERT(KEL))*FRZ(KEL)
      GARS=SIGG(L)
      DO K3=0,IELEM-1
        DO K2=0,IELEM-1
          DO K1=0,IELEM-1
            JND1=(NUM-1)*IELEM**3+K3*IELEM**2+K2*IELEM+K1+1
            JND2=(KN(NUM,1)-1)*IELEM**3+K3*IELEM**2+K2*IELEM+K1+1
            JND3=(KN(NUM,2)-1)*IELEM**3+K3*IELEM**2+K2*IELEM+K1+1
            SUNKNO(JND1)=SUNKNO(JND1)+FUNKNO(JND1)*VOL0*GARS
            SUNKNO(JND2)=SUNKNO(JND2)+FUNKNO(JND2)*VOL0*GARS
            SUNKNO(JND3)=SUNKNO(JND3)+FUNKNO(JND3)*VOL0*GARS
          ENDDO ! K1
        ENDDO ! K2
      ENDDO ! K3
    ENDDO ! KEL
    DEALLOCATE(FRZ,IPERT,KN,ZZ)
  END SUBROUTINE TRI2GSR
END SUBROUTINE DOORS_TRI2
