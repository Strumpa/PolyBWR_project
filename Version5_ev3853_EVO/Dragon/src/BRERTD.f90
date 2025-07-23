SUBROUTINE BRERTD(IELEM,ICOL,NGRP,DELX,DIFF,SIGT,SUNXS,JXM,JXP,IMPX, &
& FHOMM,FHOMP)
!
!-----------------------------------------------------------------------
!
!Purpose:
! Compute the Raviart-Thomas boundary fluxes for a single node.
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
! IELEM   Raviart-Thomas polynomial order.
! ICOL    Raviart-Thomas polynomial integration type.
! NGRP    number of energy groups.
! DELX    node width along X-axis.
! DIFF    diffusion coefficient array (cm).
! SIGT    total XS (cm^-1).
! SUNXS   production (fission + scattering) cross sections. The second
!         dimension is for primary neutrons.
! JXM     left boundary currents.
! JXP     right boundary currents.
! IMPX    print flag.
!
!Parameters: output
! FHOMM   left boundary fluxes.
! FHOMP   right boundary fluxes.
!
!-----------------------------------------------------------------------
  USE GANLIB
  !
  !----
  !  subroutine arguments
  !----
  INTEGER, INTENT(IN) :: IELEM,ICOL,NGRP,IMPX
  REAL, INTENT(IN) :: DELX
  REAL, DIMENSION(NGRP), INTENT(IN) :: DIFF,JXM,JXP
  REAL, DIMENSION(NGRP), INTENT(IN) :: SIGT
  REAL, DIMENSION(NGRP,NGRP), INTENT(IN) :: SUNXS
  REAL, DIMENSION(NGRP), INTENT(OUT) :: FHOMM,FHOMP
  !----
  !  LOCAL VARIABLES
  !----
  TYPE(C_PTR) IPTRK
  REAL QQ(5,5)
  !----
  !  ALLOCATABLE ARRAYS
  !----
  REAL, DIMENSION(:,:), ALLOCATABLE :: R,V,SYS,FUNKNO
  !----
  !  RECOVER FINITE ELEMENT UNIT MATRICES.
  !----
  CALL LCMOP(IPTRK,'***DUMMY***',0,1,0)
  CALL BIVCOL(IPTRK,IMPX,IELEM,ICOL)
  ALLOCATE(V(IELEM+1,IELEM),R(IELEM+1,IELEM+1))
  CALL LCMSIX(IPTRK,'BIVCOL',1)
  CALL LCMGET(IPTRK,'V',V)
  CALL LCMGET(IPTRK,'R',R)
  CALL LCMSIX(IPTRK,' ',2)
  CALL LCMCL(IPTRK,2)
  !----
  !  LEAKAGE-REMOVAL SYSTEM MATRIX ASSEMBLY.
  !----
  ALLOCATE(SYS(IELEM*NGRP,IELEM*NGRP+1))
  SYS(:IELEM*NGRP,:IELEM*NGRP+1)=0.0
  DO I0=1,IELEM
    DO J0=1,IELEM
      QQ(I0,J0)=0.0
      DO K0=2,IELEM
        QQ(I0,J0)=QQ(I0,J0)+V(K0,I0)*V(K0,J0)/R(K0,K0)
      ENDDO
    ENDDO
  ENDDO
  INX=IELEM*NGRP+1
  DO IG=1,NGRP
    DO J0=1,IELEM
      JND1=(IG-1)*IELEM+J0
      SYS(JND1,JND1)=SYS(JND1,JND1)+DELX*SIGT(IG)
      DO K0=1,IELEM
        IF(QQ(J0,K0).EQ.0.0) CYCLE
        KND1=(IG-1)*IELEM+K0
        SYS(JND1,KND1)=SYS(JND1,KND1)+QQ(J0,K0)*DIFF(IG)/DELX
      ENDDO
      SYS(JND1,INX)=SYS(JND1,INX)-V(1,J0)*JXM(IG)-V(IELEM+1,J0)*JXP(IG)
    ENDDO
  ENDDO
  !----
  !  SOURCE SYSTEM MATRIX ASSEMBLY.
  !----
  DO IG=1,NGRP
    DO JG=1,NGRP
      DO J0=1,IELEM
        JND1=(IG-1)*IELEM+J0
        JND2=(JG-1)*IELEM+J0
        SYS(JND1,JND2)=SYS(JND1,JND2)-DELX*SUNXS(IG,JG) ! IG <-- JG
      ENDDO
    ENDDO
  ENDDO
  !----
  !  SOLVE A ONE-NODE PROBLEM.
  !----
  CALL ALSB(IELEM*NGRP,1,SYS,IER,IELEM*NGRP)
  IF(IER.NE.0) CALL XABORT('BRERTD: SINGULAR MATRIX.')
  ALLOCATE(FUNKNO(IELEM,NGRP))
  DO IG=1,NGRP
    DO J0=1,IELEM
      FUNKNO(J0,IG)=SYS((IG-1)*IELEM+J0,IELEM*NGRP+1)
    ENDDO
  ENDDO
  DEALLOCATE(SYS,R,V)
  !----
  !  COMPUTE NODAL SURFACE FLUXES
  !----
  DENOM=REAL(IELEM*(IELEM+1))
  FHOMM(:NGRP)=0.0
  FHOMP(:NGRP)=0.0
  DO IG=1,NGRP
    IF(ICOL.EQ.1) THEN
      ! NEM relations
      IF(IELEM.EQ.1) THEN
        FHOMM(IG)=FUNKNO(1,IG)+((1./3.)*JXM(IG)+(1./6.)*JXP(IG))*DELX/DIFF(IG)
        FHOMP(IG)=FUNKNO(1,IG)-((1./6.)*JXM(IG)+(1./3.)*JXP(IG))*DELX/DIFF(IG)
      ELSE IF(IELEM.EQ.2) THEN
        FHOMM(IG)=FUNKNO(1,IG)-(5.0*SQRT(3.0)/6.0)*FUNKNO(2,IG)+((1./8.)*JXM(IG)- &
        & (1./24.)*JXP(IG))*DELX/DIFF(IG)
        FHOMP(IG)=FUNKNO(1,IG)+(5.0*SQRT(3.0)/6.0)*FUNKNO(2,IG)-(-(1./24.)*JXM(IG)+ &
        & (1./8.)*JXP(IG))*DELX/DIFF(IG)
      ELSE IF(IELEM.EQ.3) THEN
        FHOMM(IG)=FUNKNO(1,IG)-(5.0*SQRT(3.0)/6.0)*FUNKNO(2,IG)+(7.0*SQRT(5.0)/10.0)*FUNKNO(3,IG)+ &
        & ((1./15.)*JXM(IG)+(1./60.)*JXP(IG))*DELX/DIFF(IG)
        FHOMP(IG)=FUNKNO(1,IG)+(5.0*SQRT(3.0)/6.0)*FUNKNO(2,IG)+(7.0*SQRT(5.0)/10.0)*FUNKNO(3,IG)- &
        & ((1./60.)*JXM(IG)+(1./15.)*JXP(IG))*DELX/DIFF(IG)
      ELSE
        CALL XABORT('BRERTD: IELEM OVERFLOW.')
      ENDIF
    ELSE IF(ICOL.EQ.2) THEN
      ! MCFD relations
      FHOMM(IG)=JXM(IG)*DELX/(DIFF(IG)*DENOM)
      FHOMP(IG)=-JXP(IG)*DELX/(DIFF(IG)*DENOM)
      DO J0=1,IELEM
        FP=SQRT(REAL(2*J0-1))*(1.0-REAL(J0*(J0-1))/DENOM)
        FM=FP*(-1.0)**(J0-1)
        FHOMM(IG)=FHOMM(IG)+FP*(-1.0)**(J0-1)*FUNKNO(J0,IG)
        FHOMP(IG)=FHOMP(IG)+FP*FUNKNO(J0,IG)
      ENDDO
    ELSE IF(ICOL.EQ.3) THEN
      IF(IELEM.EQ.1) THEN
        FHOMM(IG)=FUNKNO(1,IG)+((1./4.)*JXM(IG)+(1./4.)*JXP(IG))*DELX/DIFF(IG)
        FHOMP(IG)=FUNKNO(1,IG)-((1./4.)*JXM(IG)+(1./4.)*JXP(IG))*DELX/DIFF(IG)
      ELSE IF(IELEM.EQ.2) THEN
        FHOMM(IG)=FUNKNO(1,IG)-SQRT(3.0)*FUNKNO(2,IG)+((1./12.)*JXM(IG)- &
        & (1./12.)*JXP(IG))*DELX/DIFF(IG)
        FHOMP(IG)=FUNKNO(1,IG)+SQRT(3.0)*FUNKNO(2,IG)-(-(1./12.)*JXM(IG)+ &
        & (1./12.)*JXP(IG))*DELX/DIFF(IG)
      ELSE IF(IELEM.EQ.3) THEN
        FHOMM(IG)=FUNKNO(1,IG)-(5.0*SQRT(3.0)/6.0)*FUNKNO(2,IG)+SQRT(5.0)*FUNKNO(3,IG)+ &
        & ((1./24.)*JXM(IG)+(1./24.)*JXP(IG))*DELX/DIFF(IG)
        FHOMP(IG)=FUNKNO(1,IG)+(5.0*SQRT(3.0)/6.0)*FUNKNO(2,IG)+SQRT(5.0)*FUNKNO(3,IG)- &
        & ((1./24.)*JXM(IG)+(1./24.)*JXP(IG))*DELX/DIFF(IG)
      ELSE
        CALL XABORT('BRERTD: IELEM OVERFLOW.')
      ENDIF
    ELSE
      CALL XABORT('BRERTD: ICOL OVERFLOW.')
    ENDIF
    IF(IMPX.GT.0) THEN
      WRITE(6,'(46H BRERTD: RAVIART-THOMAS FLUX UNKNOWNS IN GROUP,I4,1H:/(1P,12E12.4))') &
      & IG,FUNKNO(:IELEM,IG)
      WRITE(6,'(48H BRERTD: RAVIART-THOMAS BOUNDARY FLUXES IN GROUP,I4,1H:/(1P,12E12.4))') &
      & IG,FHOMM(IG),FHOMP(IG)
    ENDIF
  ENDDO
  DEALLOCATE(FUNKNO)
END SUBROUTINE BRERTD
