SUBROUTINE OUTVAL(IPMAC1,IPMAC2,IPVAL,NBMIX,NL,NBFIS,NGRP,NALBP,NGCOND, &
& IGCOND,ZNORM,IMPX)
  !
  !-----------------------------------------------------------------------
  !
  !Purpose:
  ! Perform an homogenization based on a voxelized flux.
  !
  !Copyright:
  ! Copyright (C) 2026 Ecole Polytechnique de Montreal
  ! This library is free software; you can redistribute it and/or
  ! modify it under the terms of the GNU Lesser General Public
  ! License as published by the Free Software Foundation; either
  ! version 2.1 of the License, or (at your option) any later version
  !
  !Author(s): A. Hebert
  !
  !Parameters: input
  ! IPMAC1  L_MACROLIB pointer to the input macrolib.
  ! IPMAC2  L_MACROLIB pointer to the output extended macrolib.
  ! IPVAL   L_FVIEW pointer to the interpflux data structure.
  ! NBMIX   number of material mixtures.
  ! NL      scattering anisotropy.
  ! NBFIS   number of fissionable isotopes.
  ! NGRP    total number of energy groups.
  ! NALBP   number of physical albedos.
  ! NGCOND  number of macrogroups after energy condensation.
  ! IGCOND  limit of condensed groups.
  ! ZNORM   flux normalization factor.
  ! IMPX    print parameter (equal to zero for no print).
  !
  !-----------------------------------------------------------------------
  !
  USE GANLIB
  IMPLICIT NONE
  !----
  !  Subroutine arguments
  !----
  TYPE(C_PTR) IPMAC1,IPMAC2,IPVAL
  INTEGER NBMIX,NL,NBFIS,NGRP,NALBP,NGCOND,IGCOND(NGCOND),IMPX
  DOUBLE PRECISION ZNORM
  !----
  !  Local variables
  !----
  TYPE(C_PTR) JPVAL
  INTEGER, PARAMETER :: NSTATE=40
  INTEGER I,J,K,IS,JS,KS,IXLG,IYLG,IZLG,NZS,NZSOLD,IM,II,IUNK,IBM,IG,IM_KEEP, &
  & ITYPLU,INTLIR,ISTATE(NSTATE),ILONG,ITYLCM
  REAL REALIR,SXYZ(3)
  DOUBLE PRECISION DBLLIR
  CHARACTER CARLIR*12
  LOGICAL L3D,FV1,FV2
  !----
  !  Allocatable arrays
  !----
  INTEGER, DIMENSION(:), ALLOCATABLE :: ITNODE,IREMIX,MAT_VOX,IDL_VOX,IHOM_VOX
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: MATXYZ
  REAL, DIMENSION(:), ALLOCATABLE :: MXI,MYI,MZI,VOL_VOX
  REAL, DIMENSION(:,:), ALLOCATABLE :: FLU_VOX
  REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: FLUX
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: NODE
  !----
  !  Read homogeneous node definitions
  !----
  L3D=.FALSE.
  CALL REDGET(ITYPLU,NZS,REALIR,CARLIR,DBLLIR)
  IF(ITYPLU.NE.1) CALL XABORT('OUTVAL: INTEGER VARIABLE EXPECTED.')
  IF(NZS.LE.0) CALL XABORT('OUTVAL: INVALID VALUE OF NZS.')
  ALLOCATE(NODE(8,NZS),ITNODE(NZS))
  CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
  IF(ITYPLU.NE.3) CALL XABORT('OUTVAL: CHARACTER VARIABLE EXPECTED(1).')
  DO IM=1,NZS
    IF(CARLIR.EQ.'RECT') THEN
      ITNODE(IM)=1
      DO I=1,4
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.2) CALL XABORT('OUTVAL: REAL VARIABLE EXPECTED(1).')
        NODE(I,IM)=REALIR
      ENDDO
      NODE(5:6,IM)=0.0D0
    ELSE IF(CARLIR.EQ.'TRIA') THEN
      ITNODE(IM)=2
      DO I=1,6
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.2) CALL XABORT('OUTVAL: REAL VARIABLE EXPECTED(2).')
        NODE(I,IM)=REALIR
      ENDDO
    ELSE IF(CARLIR.EQ.'QUAD') THEN
      ITNODE(IM)=3
      DO I=1,8
        CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
        IF(ITYPLU.NE.2) CALL XABORT('OUTVAL: REAL VARIABLE EXPECTED(3).')
        NODE(I,IM)=REALIR
      ENDDO
    ELSE
      CALL XABORT('OUTVAL: *RECT*, *TRIA* OR *QUAD* KEYWORD EXPECTED.')
    ENDIF
    CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
    IF(ITYPLU.NE.3) CALL XABORT('OUTVAL: CHARACTER DATA EXPECTED(2).')
    IF(CARLIR.EQ.'AXIAL') THEN
      L3D=.TRUE.
      CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
      IF(ITYPLU.NE.2) CALL XABORT('OUTVAL: REAL VARIABLE EXPECTED(4).')
      NODE(1,IM)=REALIR
      CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
      IF(ITYPLU.NE.2) CALL XABORT('OUTVAL: REAL VARIABLE EXPECTED(5).')
      NODE(2,IM)=REALIR
      CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
      IF(ITYPLU.NE.3) CALL XABORT('OUTVAL: CHARACTER DATA EXPECTED(3).')
    ENDIF
  ENDDO
  !----
  !  Recover voxelized information
  !----
  CALL LCMGET(IPVAL,'STATE-VECTOR',ISTATE)
  IF(ISTATE(1).NE.NGRP) CALL XABORT('OUTVAL: invalid number of groups.')
  IXLG=ISTATE(2)
  IYLG=ISTATE(3)
  IZLG=ISTATE(4)
  ALLOCATE(MAT_VOX(IXLG*IYLG*IZLG),IHOM_VOX(IXLG*IYLG*IZLG),VOL_VOX(IXLG*IYLG*IZLG), &
  & FLU_VOX(IXLG*IYLG*IZLG,NGRP))
  ALLOCATE(MXI(IXLG),MYI(IYLG),MZI(IZLG),MATXYZ(IXLG,IYLG,IZLG),FLUX(IXLG,IYLG,IZLG,NGRP))
  CALL LCMGET(IPVAL,'SXYZ',SXYZ)
  CALL LCMGET(IPVAL,'MXI',MXI)
  CALL LCMLEN(IPVAL,'MYI',ILONG,ITYLCM)
  IF(ILONG.GT.0) THEN
    CALL LCMGET(IPVAL,'MYI',MYI)
  ELSE
    MYI(1)=0.0
  ENDIF
  CALL LCMLEN(IPVAL,'MZI',ILONG,ITYLCM)
  IF(ILONG.GT.0) THEN
    CALL LCMGET(IPVAL,'MZI',MZI)
  ELSE
    MZI(1)=0.0
  ENDIF
  CALL LCMGET(IPVAL,'MATXYZ',MATXYZ)
  JPVAL=LCMGID(IPVAL,'FLUX')
  DO IG=1,NGRP
    CALL LCMGDL(JPVAL,IG,FLUX(1,1,1,IG))
  ENDDO
  IUNK=0
  DO K=1,IZLG
    DO J=1,IYLG
      DO I=1,IXLG
        IBM=MATXYZ(I,J,K)
        IF(IBM.EQ.0) CYCLE
        IM_KEEP=0
        IS=0 ; JS=0 ; KS=0
        DO IM=1,NZS
          IF(ITNODE(IM).EQ.1) THEN
            IS=I
            IF((MXI(I).GE.NODE(1,IM)).AND.(MXI(I).LE.NODE(2,IM))) GO TO 10
            CYCLE
            10 JS=J
            IF((MYI(J).GE.NODE(3,IM)).AND.(MYI(J).LE.NODE(4,IM))) GO TO 20
            CYCLE
          ELSE IF(ITNODE(IM).EQ.2) THEN
            FV1=OUT_CONTAINS(NODE(1,IM), NODE(3,IM), NODE(5,IM), MXI(I), MYI(J))
            IF(.NOT.FV1) CYCLE
            IS=I ; JS=J
          ELSE IF(ITNODE(IM).EQ.3) THEN
            FV1=OUT_CONTAINS(NODE(1,IM), NODE(3,IM), NODE(5,IM), MXI(I), MYI(J))
            FV2=OUT_CONTAINS(NODE(1,IM), NODE(5,IM), NODE(7,IM), MXI(I), MYI(J))
            IF((.NOT.FV1).AND.(.NOT.FV2)) CYCLE
            IS=I ; JS=J
          ENDIF
          20 IF(L3D) THEN
            KS=K
            IF((MZI(K).GE.NODE(1,IM)).AND.(MZI(K).LE.NODE(2,IM))) GO TO 30
            CYCLE
          ELSE
            KS=K
          ENDIF
          30 IM_KEEP=IM
          EXIT
        ENDDO ! IM
        IF((IS.EQ.0).OR.(JS.EQ.0.).OR.(KS.EQ.0)) CYCLE
        IUNK=IUNK+1
        MAT_VOX(IUNK)=IBM
        IF(IM_KEEP.GT.NZS) call XABORT('OUTVAL: IM overflow.')
        IHOM_VOX(IUNK)=IM_KEEP
        VOL_VOX(IUNK)=SXYZ(1)*SXYZ(2)*SXYZ(3)
        IF((I.EQ.1).OR.(I.EQ.IXLG)) VOL_VOX(IUNK)=VOL_VOX(IUNK)/2.0
        IF((J.EQ.1).OR.(J.EQ.IYLG)) VOL_VOX(IUNK)=VOL_VOX(IUNK)/2.0
        IF((IZLG.GT.1).AND.((K.EQ.1).OR.(K.EQ.IZLG))) VOL_VOX(IUNK)=VOL_VOX(IUNK)/2.0
        FLU_VOX(IUNK,:)=FLUX(IS,JS,KS,:)*REAL(ZNORM)
      ENDDO ! I
    ENDDO ! J
  ENDDO ! K
  DEALLOCATE(FLUX,MATXYZ,MZI,MYI,MXI)
  !----
  !  Remix homogenized indices
  !----
  IF(CARLIR.EQ.'REMIX') THEN
    NZSOLD=NZS
    NZS=0
    ALLOCATE(IREMIX(NZSOLD))
    DO II=1,NZSOLD
      CALL REDGET(ITYPLU,IREMIX(II),REALIR,CARLIR,DBLLIR)
      IF(ITYPLU.NE.1) CALL XABORT('OUTVAL: INTEGER DATA EXPECTED(4).')
    ENDDO
    DO K=1,IUNK
      IM=IHOM_VOX(K)
      IF(IM.GT.0) THEN
        IF(IM.GT.NZSOLD) CALL XABORT('OUTVAL: IHOM_VOX OVERFLOW.')
        IHOM_VOX(K)=IREMIX(IM)
        NZS=MAX(NZS,IHOM_VOX(K))
      ENDIF
    ENDDO
    CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
    IF(ITYPLU.NE.3) CALL XABORT('OUTVAL: CHARACTER DATA EXPECTED(4).')
    DEALLOCATE(IREMIX)
  ENDIF
  IF(CARLIR.NE.';') CALL XABORT('OUTVAL: ; expected.')
  !----
  !  Perform homogenization and condensation
  !----
  IF(IUNK.GT.0) THEN
    ALLOCATE(IDL_VOX(IUNK))
    DO I=1,IUNK
      IDL_VOX(I)=I
    ENDDO
    CALL OUTAUX(IPMAC1,IPMAC2,NBMIX,NL,NBFIS,NGRP,IUNK,IXLG*IYLG*IZLG,NALBP, &
    & NZS,NGCOND,MAT_VOX,VOL_VOX,IDL_VOX,FLU_VOX,IHOM_VOX,IGCOND,IMPX)
    DEALLOCATE(IDL_VOX)
  ENDIF
  DEALLOCATE(FLU_VOX,VOL_VOX,IHOM_VOX,MAT_VOX)
  RETURN
  !
  CONTAINS
    FUNCTION OUT_CONTAINS(AXY, BXY, CXY, X, Y) RESULT(FVALUE)
      DOUBLE PRECISION, INTENT(IN) :: AXY(2),BXY(2),CXY(2)
      REAL, INTENT(IN) :: X,Y
      LOGICAL FVALUE
      DOUBLE PRECISION DET
      !
      DET = (BXY(1) - AXY(1)) * (CXY(2) - AXY(2)) - (BXY(2) - AXY(2)) * (CXY(1) - AXY(1))
      FVALUE = DET * ((BXY(1) - AXY(1)) * (Y - AXY(2)) - (BXY(2) - AXY(2)) * (X - AXY(1))) >= 0.0D0 .AND. &
             & DET * ((CXY(1) - BXY(1)) * (Y - BXY(2)) - (CXY(2) - BXY(2)) * (X - BXY(1))) >= 0.0D0 .AND. &
             & DET * ((AXY(1) - CXY(1)) * (Y - CXY(2)) - (AXY(2) - CXY(2)) * (X - CXY(1))) >= 0.0D0
    END FUNCTION OUT_CONTAINS
END SUBROUTINE OUTVAL
