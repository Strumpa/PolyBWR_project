SUBROUTINE OUTVAL(IPGEO2,IPMAC1,IPMAC2,IPVAL,NBMIX,NL,NBFIS,NGRP,NALBP, &
& NGCOND,IGCOND,ZNORM,IMPX)
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
  ! IPGEO2  L_GEOM pointer to the macrogeometry containing INTG data.
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
  TYPE(C_PTR) IPGEO2,IPMAC1,IPMAC2,IPVAL
  INTEGER NBMIX,NL,NBFIS,NGRP,NALBP,NGCOND,IGCOND(NGCOND),IMPX
  DOUBLE PRECISION ZNORM
  !----
  !  Local variables
  !----
  TYPE(C_PTR) JPVAL
  INTEGER, PARAMETER :: NSTATE=40
  INTEGER ISTATE(NSTATE),I,J,K,IT,IS,JS,KS,IXLG,IYLG,IZLG,NZS,IM,IUNK,IBM, &
  & IG,IM_KEEP,NMERGE,NZS_NEW,IDIM,ILONG,ITYLCM
  REAL SXYZ(3)
  LOGICAL L3D,LTEST
  !----
  !  Allocatable arrays
  !----
  INTEGER, DIMENSION(:), ALLOCATABLE :: ITNODE,IREMIX,MAT_VOX,IDL_VOX,IHOM_VOX
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: MATXYZ
  REAL, DIMENSION(:), ALLOCATABLE :: MXI,MYI,MZI,VOL_VOX
  REAL, DIMENSION(:,:), ALLOCATABLE :: FLU_VOX
  REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: FLUX
  REAL, DIMENSION(:,:), ALLOCATABLE :: NODE
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: TLAMB
  !----
  !  Recover macrogeometry
  !----
  IF(.NOT.C_ASSOCIATED(IPGEO2)) THEN
    CALL XABORT('OUTVAL: MACROGEOMETRY IS NOT DEFINED.')
  ENDIF
  CALL LCMGET(IPGEO2,'STATE-VECTOR',ISTATE)
  IF(ISTATE(1).NE.31) CALL XABORT('OUTVAL: INVALID TYPE OF MACROGEOMETRY.')
  CALL LCMSIX(IPGEO2,'MACRO-GEOM',1)
    CALL LCMGET(IPGEO2,'STATE-VECTOR',ISTATE)
    IDIM=ISTATE(1)
    NZS=ISTATE(2)
    NMERGE=ISTATE(4)
    L3D=(ISTATE(5).EQ.3)
    ALLOCATE(NODE(IDIM,NZS),ITNODE(NZS))
    CALL LCMGET(IPGEO2,'NODE',NODE)
    CALL LCMGET(IPGEO2,'ITNODE',ITNODE)
    IF(NMERGE.GT.0) THEN
      ALLOCATE(IREMIX(NZS))
      CALL LCMGET(IPGEO2,'IREMIX',IREMIX)
    ENDIF
    IF(IMPX.GT.0) WRITE(6,200) ISTATE(:5)
  CALL LCMSIX(IPGEO2,' ',2)
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
            IF(IDIM.LT.4) CALL XABORT('OUTVAL: NODE OVERFLOW(1).')
            IS=I
            IF((MXI(I).GE.NODE(1,IM)).AND.(MXI(I).LE.NODE(2,IM))) GO TO 10
            CYCLE
            10 JS=J
            IF((MYI(J).GE.NODE(3,IM)).AND.(MYI(J).LE.NODE(4,IM))) GO TO 20
            CYCLE
          ELSE IF(ITNODE(IM).EQ.2) THEN
            CALL XABORT('OUTVAL: INVALID ITNODE VALUE.')
          ELSE IF(ITNODE(IM).GE.3) THEN
            IF(IDIM.LT.2*ITNODE(IM)) CALL XABORT('OUTVAL: NODE OVERFLOW(2).')
            ALLOCATE(TLAMB(ITNODE(IM)))
            TLAMB=OUTBAR(ITNODE(IM),NODE(1,IM),MXI(I), MYI(J))
            LTEST=.TRUE.
            DO IT=1,ITNODE(IM)
              LTEST=LTEST.AND.(TLAMB(IT).GE.-1.0E-6)
            ENDDO
            DEALLOCATE(TLAMB)
            IF(.NOT.LTEST) CYCLE
            IS=I ; JS=J
          ENDIF
          20 IF(L3D) THEN
            IF(IDIM.LT.14) CALL XABORT('OUTVAL: NODE OVERFLOW(3).')
            KS=K
            IF((MZI(K).GE.NODE(13,IM)).AND.(MZI(K).LE.NODE(14,IM))) GO TO 30
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
        IF(IM_KEEP.GT.NZS) call XABORT('OUTVAL: IM OVERFLOW.')
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
  IF(NMERGE.GT.0) THEN
    NZS_NEW=0
    DO K=1,IUNK
      IM=IHOM_VOX(K)
      IF(IM.GT.0) THEN
        IF(IM.GT.NZS) CALL XABORT('OUTVAL: IHOM_VOX OVERFLOW.')
        IHOM_VOX(K)=IREMIX(IM)
        NZS_NEW=MAX(NZS_NEW,IHOM_VOX(K))
      ENDIF
    ENDDO
    DEALLOCATE(IREMIX)
    IF(NZS_NEW.NE.NMERGE) CALL XABORT('OUTVAL: INVALID NMERGE.')
  ELSE
    NMERGE=NZS
  ENDIF
  !----
  !  Perform homogenization and condensation
  !----
  IF(IUNK.GT.0) THEN
    ALLOCATE(IDL_VOX(IUNK))
    DO I=1,IUNK
      IDL_VOX(I)=I
    ENDDO
    CALL OUTAUX(IPMAC1,IPMAC2,NBMIX,NL,NBFIS,NGRP,IUNK,IXLG*IYLG*IZLG,NALBP, &
    & NMERGE,NGCOND,MAT_VOX,VOL_VOX,IDL_VOX,FLU_VOX,IHOM_VOX,IGCOND,IMPX)
    DEALLOCATE(IDL_VOX)
  ENDIF
  DEALLOCATE(FLU_VOX,VOL_VOX,IHOM_VOX,MAT_VOX)
  200 FORMAT(/22H MACROGEOMETRY OPTIONS/1X,21(1H-)/ &
  & 7H IDIM  ,I8,34H   (FIRST DIMENSION OF ARRAY NODE)/ &
  & 7H NZS   ,I8,23H   (NUMBER OF MIXTURES)/ &
  & 7H ISYM  ,I8,39H   (=0/1/2/3: UNDEFINED/NONE/REFL/ROTA)/ &
  & 7H NMERGE,I8,33H   (=0/NUMBER OF MERGED MIXTURES)/ &
  & 7H IDIM  ,I8,16H   (=2/3: 2D/3D))
  RETURN
  !
  CONTAINS
    FUNCTION OUTBAR(N, NODE, CX, CY) RESULT(TLAMB)
      ! Compute the barycentric coordinates of (CX,CY)
      INTEGER, INTENT(IN)   :: N             ! number of vertices
      REAL, INTENT(IN)  :: NODE(2*N)         ! polygon vertex coordinates
      REAL, INTENT(IN)  :: CX, CY            ! target point coordinates
      DOUBLE PRECISION :: TLAMB(N)           ! resulting normalized weights

      DOUBLE PRECISION, ALLOCATABLE :: R(:), ALPHA(:), W(:)
      DOUBLE PRECISION :: DX, DY, DX_NEXT, DY_NEXT, DOT, CROSS, TOTAL_WEIGHT
      INTEGER  :: I, IM1, IP1

      ALLOCATE(R(N), ALPHA(N), W(N))
      TLAMB = 0.0D0
      W = 0.0D0

      ! step 1: precompute distances (R) from (CX,CY) to all vertices
      DO I = 1, N
        R(I) = SQRT((NODE(2*I-1) - CX)**2 + (NODE(2*I) - CY)**2)
            
        ! ROBUSTNESS CHECK: IF POINT IS EXACTLY ON CORNER I
        IF(R(I) < 1.0D-12) THEN
          TLAMB = 0.0D0
          TLAMB(I) = 1.0D0
          RETURN
        ENDIF
      ENDDO

      ! step 2: compute signed angles (ALPHA) between adjacent vectors using atan2
      DO I = 1, N
        IP1 = I + 1; IF(IP1 > N) IP1 = 1 ! wrap around to first vertex
            
        ! vectors from (CX,CY) to CORNER(I) and CORNER(I+1)
        DX = NODE(2*I-1) - CX
        DY = NODE(2*I) - CY
        DX_NEXT = NODE(2*IP1-1) - CX
        DY_NEXT = NODE(2*IP1) - CY

        ! dot product and 2d cross product determinant
        DOT = DX * DX_NEXT + DY * DY_NEXT
        CROSS = DX * DY_NEXT - DY * DX_NEXT

        ! ATAN2 yields the exact signed angle between vectors in [-pi, pi]
        ALPHA(I) = ATAN2(CROSS, DOT)

        ! robustness check: if point lies precisely on the segment between I and I+1
        IF(ABS(CROSS) < 1.0D-12 .AND. DOT < 0.0D0) THEN
          TLAMB = 0.0D0
          TLAMB(I) = R(IP1) / (R(I) + R(IP1))
          TLAMB(IP1) = 1.0D0 - TLAMB(I)
          RETURN
        ENDIF
      ENDDO

      ! step 3: compute unnormalized TLAMB using the tangent of half-angles
      DO I = 1, N
        IM1 = I - 1; IF(IM1 < 1) IM1 = N ! wrap around to last vertex
        W(I) = (TAN(ALPHA(IM1) / 2.0D0) + TAN(ALPHA(I) / 2.0D0)) / R(I)
      ENDDO

      ! step 4: normalize TLAMB to sum up to 1.0
      TOTAL_WEIGHT = SUM(W)
      IF(ABS(TOTAL_WEIGHT) > 1.0D-12) THEN
        TLAMB = W / TOTAL_WEIGHT
      ELSE
        TLAMB = 0.0D0
      ENDIF
    END FUNCTION OUTBAR
END SUBROUTINE OUTVAL
