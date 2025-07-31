SUBROUTINE BRERTS(IELEM,ICOL,NGRP,NLF,DELX,RCAT,JXM,JXP,IMPX,FHOMM,FHOMP)
!
!-----------------------------------------------------------------------
!
!Purpose:
! Compute the Raviart-Thomas boundary fluxes for a single node in SPN
! theory.
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
! NLF     (NLF-1) is the SPN order (NLF is an even integer).
! DELX    node width along X-axis.
! RCAT    removal matrix (total minus scattering cross sections). The
!         second dimension is for primary neutrons. The (2*IL+1) factor
!         is included.
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
  !  SUBROUTINE ARGUMENTS
  !----
  INTEGER, INTENT(IN) :: IELEM,ICOL,NGRP,NLF,IMPX
  REAL, INTENT(IN) :: DELX
  REAL, DIMENSION(NGRP,NLF/2), INTENT(IN) :: JXM,JXP
  REAL(KIND=8), DIMENSION(NGRP,NGRP,NLF), INTENT(IN) :: RCAT
  REAL, DIMENSION(NGRP,NLF/2), INTENT(OUT) :: FHOMM,FHOMP
  !----
  !  LOCAL VARIABLES
  !----
  TYPE(C_PTR) IPTRK
  REAL QQ(5,5)
  REAL(KIND=8) FHOMM_IG,FHOMP_IG
  !----
  !  ALLOCATABLE ARRAYS
  !----
  INTEGER, DIMENSION(:), ALLOCATABLE :: IND
  REAL, DIMENSION(:,:), ALLOCATABLE :: R,V
  REAL, DIMENSION(:,:,:), ALLOCATABLE :: FUNKNO
  REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: DXM,DXP,SYS
  REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: RCATI
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
  DO I0=1,IELEM
    DO J0=1,IELEM
      QQ(I0,J0)=0.0
      DO K0=2,IELEM
        QQ(I0,J0)=QQ(I0,J0)+V(K0,I0)*V(K0,J0)/R(K0,K0)
      ENDDO
    ENDDO
  ENDDO
  !----
  !  INVERT THE RESIDUAL MATRIX.
  !----
  IF(MOD(NLF,2).NE.0) CALL XABORT('BRERTS: EVEN NLF EXPECTED.')
  ALLOCATE(RCATI(NGRP,NGRP,NLF),IND(NGRP))
  DO IL=2,NLF,2
    RCATI(:NGRP,:NGRP,IL)=RCAT(:NGRP,:NGRP,IL)
    CALL ALINVD(NGRP,RCATI(1,1,IL),NGRP,IER,IND)
    IF(IER.NE.0) CALL XABORT('BRERTS: SINGULAR MATRIX(1).')
  ENDDO
  DEALLOCATE(IND)
  !----
  !  LEAKAGE-REMOVAL SYSTEM MATRIX ASSEMBLY FOR THE RAVIART-THOMAS METHOD.
  !----
  FHOMM(:NGRP,:NLF/2)=0.0
  FHOMP(:NGRP,:NLF/2)=0.0
  L4=IELEM*NGRP
  INX=L4*NLF/2+1
  ALLOCATE(SYS(L4*NLF/2,INX))
  SYS(:L4*NLF/2,:INX)=0.0D0
  DO IL=0,NLF-1
    IF(MOD(IL,2).EQ.0) THEN
      !----
      !  EVEN PARITY EQUATION.
      !----
      DO IG=1,NGRP
        DO JG=1,NGRP
          DO J0=1,IELEM
            JND1=(IL/2)*L4+(IG-1)*IELEM+J0
            JND2=(IL/2)*L4+(JG-1)*IELEM+J0
            SYS(JND1,JND2)=SYS(JND1,JND2)+DELX*RCAT(IG,JG,IL+1) ! IG <-- JG
          ENDDO
        ENDDO
      ENDDO
    ELSE
      !----
      !  ODD PARITY EQUATION.
      !----
      DO IG=1,NGRP
        IF(IELEM.GT.1) THEN
          ! get rid of net current collocation points inside finite elements.
          DO JG=1,NGRP
            DO J0=1,IELEM
              JND1=((IL-1)/2)*L4+(IG-1)*IELEM+J0
              DO K0=1,IELEM
                IF(QQ(J0,K0).EQ.0.0) CYCLE
                KND2=((IL-1)/2)*L4+(JG-1)*IELEM+K0
                SYS(JND1,KND2)=SYS(JND1,KND2)+(REAL(IL)**2)*QQ(J0,K0)*RCATI(IG,JG,IL+1)/DELX
                IF(IL.LE.NLF-3) THEN
                  KND2=((IL+1)/2)*L4+(JG-1)*IELEM+K0
                  SYS(JND1,KND2)=SYS(JND1,KND2)+REAL(IL*(IL+1))*QQ(J0,K0)*RCATI(IG,JG,IL+1)/DELX
                ENDIF
                IF(IL.GE.3) THEN
                  KND2=((IL-3)/2)*L4+(JG-1)*IELEM+K0
                  SYS(JND1,KND2)=SYS(JND1,KND2)+REAL((IL-1)*(IL-2))*QQ(J0,K0)*RCATI(IG,JG,IL-1)/DELX
                  KND2=((IL-1)/2)*L4+(JG-1)*IELEM+K0
                  SYS(JND1,KND2)=SYS(JND1,KND2)+(REAL(IL-1)**2)*QQ(J0,K0)*RCATI(IG,JG,IL-1)/DELX
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDIF
        DO J0=1,IELEM
          JND1=((IL-1)/2)*L4+(IG-1)*IELEM+J0
          GJXM=JXM(IG,1)
          GJXP=JXP(IG,1)
          IF(IL.EQ.3) THEN
            GJXM=2.0*JXM(IG,1)+3.0*JXM(IG,2)
            GJXP=2.0*JXP(IG,1)+3.0*JXP(IG,2)
          ELSE IF(IL.EQ.5) THEN
            GJXM=4.0*JXM(IG,2)+5.0*JXM(IG,3)
            GJXP=4.0*JXP(IG,2)+5.0*JXP(IG,3)
          ELSE IF(IL.EQ.7) THEN
            GJXM=6.0*JXM(IG,3)+7.0*JXM(IG,4)
            GJXP=6.0*JXP(IG,3)+7.0*JXP(IG,4)
          ELSE IF(IL.GE.9) THEN
            CALL XABORT('BRERTS: SPN ORDER NOT IMPLEMENTED(1).')
          ENDIF
          SYS(JND1,INX)=SYS(JND1,INX)-(V(1,J0)*GJXM+V(IELEM+1,J0)*GJXP)
        ENDDO
      ENDDO
    ENDIF
  ENDDO
  !----
  !  SOLVE A ONE-NODE PROBLEM.
  !----
  CALL ALSBD(L4*NLF/2,1,SYS,IER,L4*NLF/2)
  IF(IER.NE.0) CALL XABORT('BRERTS: SINGULAR MATRIX(2).')
  ALLOCATE(FUNKNO(IELEM,NGRP,NLF/2))
  DO IL=1,NLF/2
    DO IG=1,NGRP
      DO J0=1,IELEM
        FUNKNO(J0,IG,IL)=REAL(SYS((IL-1)*L4+(IG-1)*IELEM+J0,INX))
      ENDDO
    ENDDO
  ENDDO
  DEALLOCATE(SYS,RCATI,R,V)
  !----
  !  COMPUTE SURFACE FLUX GRADIENTS USING ODD PARITY EQUATIONS
  !----
  ALLOCATE(DXM(NGRP,NLF/2),DXP(NGRP,NLF/2))
  IF(NLF.EQ.2) THEN
    DXM(:NGRP,1)=MATMUL(RCAT(:,:,2),JXM(:,1))*DELX
    DXP(:NGRP,1)=MATMUL(RCAT(:,:,2),JXP(:,1))*DELX
  ELSE IF(NLF.EQ.4) THEN
    DXM(:NGRP,2)=MATMUL(RCAT(:,:,4),JXM(:,2))*DELX/3.0
    DXP(:NGRP,2)=MATMUL(RCAT(:,:,4),JXP(:,2))*DELX/3.0
    DXM(:NGRP,1)=MATMUL(RCAT(:,:,2),JXM(:,1))*DELX-2.0*DXM(:NGRP,2)
    DXP(:NGRP,1)=MATMUL(RCAT(:,:,2),JXP(:,1))*DELX-2.0*DXP(:NGRP,2)
  ELSE IF(NLF.EQ.6) THEN
    DXM(:NGRP,3)=MATMUL(RCAT(:,:,6),JXM(:,3))*DELX/5.0
    DXP(:NGRP,3)=MATMUL(RCAT(:,:,6),JXP(:,3))*DELX/5.0
    DXM(:NGRP,2)=MATMUL(RCAT(:,:,4),JXM(:,2))*DELX/3.0-4.0*DXM(:NGRP,3)/3.0
    DXP(:NGRP,2)=MATMUL(RCAT(:,:,4),JXP(:,2))*DELX/3.0-4.0*DXP(:NGRP,3)/3.0
    DXM(:NGRP,1)=MATMUL(RCAT(:,:,2),JXM(:,1))*DELX-2.0*DXM(:NGRP,2)
    DXP(:NGRP,1)=MATMUL(RCAT(:,:,2),JXP(:,1))*DELX-2.0*DXP(:NGRP,2)
  ELSE
    CALL XABORT('BRERTS: SPN ORDER NOT IMPLEMENTED(2).')
  ENDIF
  !----
  !  COMPUTE NODAL SURFACE FLUXES
  !----
  DENOM=REAL(IELEM*(IELEM+1))
  DO IG=1,NGRP
    DO IL=1,NLF/2
      FHOMM_IG=0.0D0
      FHOMP_IG=0.0D0
      IF(ICOL.EQ.1) THEN
        ! NEM relations
        IF(IELEM.EQ.1) THEN
          FHOMM_IG=FUNKNO(1,IG,IL)+((1./3.)*DXM(IG,IL)+(1./6.)*DXP(IG,IL))
          FHOMP_IG=FUNKNO(1,IG,IL)-((1./6.)*DXM(IG,IL)+(1./3.)*DXP(IG,IL))
        ELSE IF(IELEM.EQ.2) THEN
          FHOMM_IG=FUNKNO(1,IG,IL)-(5.0*SQRT(3.)/6.)*FUNKNO(2,IG,IL)+((1./8.)*DXM(IG,IL)- &
          & (1./24.)*DXP(IG,IL))
          FHOMP_IG=FUNKNO(1,IG,IL)+(5.0*SQRT(3.)/6.)*FUNKNO(2,IG,IL)-(-(1./24.)*DXM(IG,IL)+ &
          & (1./8.)*DXP(IG,IL))
        ELSE IF(IELEM.EQ.3) THEN
          FHOMM_IG=FUNKNO(1,IG,IL)-(5.0*SQRT(3.)/6.)*FUNKNO(2,IG,IL)+(7.0*SQRT(5.)/10.)* &
          & FUNKNO(3,IG,IL)+((1./15.)*DXM(IG,IL)+(1./60.)*DXP(IG,IL))
          FHOMP_IG=FUNKNO(1,IG,IL)+(5.0*SQRT(3.)/6.)*FUNKNO(2,IG,IL)+(7.0*SQRT(5.)/10.)* &
          & FUNKNO(3,IG,IL)-((1./60.)*DXM(IG,IL)+(1./15.)*DXP(IG,IL))
        ELSE
          CALL XABORT('BRERTS: IELEM OVERFLOW.')
        ENDIF
      ELSE IF(ICOL.EQ.2) THEN
        ! MCFD relations
        FHOMM_IG=DXM(IG,IL)/DENOM
        FHOMP_IG=-DXP(IG,IL)/DENOM
        DO J0=1,IELEM
          FP=SQRT(REAL(2*J0-1))*(1.0-REAL(J0*(J0-1))/DENOM)
          FM=FP*(-1.0)**(J0-1)
          FHOMM_IG=FHOMM_IG+FP*(-1.0)**(J0-1)*FUNKNO(J0,IG,IL)
          FHOMP_IG=FHOMP_IG+FP*FUNKNO(J0,IG,IL)
        ENDDO
      ELSE IF(ICOL.EQ.3) THEN
        IF(IELEM.EQ.1) THEN
          FHOMM_IG=FUNKNO(1,IG,IL)+((1./4.)*DXM(IG,IL)+(1./4.)*DXP(IG,IL))
          FHOMP_IG=FUNKNO(1,IG,IL)-((1./4.)*DXM(IG,IL)+(1./4.)*DXP(IG,IL))
        ELSE IF(IELEM.EQ.2) THEN
          FHOMM_IG=FUNKNO(1,IG,IL)-SQRT(3.0)*FUNKNO(2,IG,IL)+((1./12.)*DXM(IG,IL)- &
          & (1./12.)*DXP(IG,IL))
          FHOMP_IG=FUNKNO(1,IG,IL)+SQRT(3.0)*FUNKNO(2,IG,IL)-(-(1./12.)*DXM(IG,IL)+ &
          & (1./12.)*DXP(IG,IL))
        ELSE IF(IELEM.EQ.3) THEN
          FHOMM_IG=FUNKNO(1,IG,IL)-(5.0*SQRT(3.)/6.0)*FUNKNO(2,IG,IL)+SQRT(5.)*FUNKNO(3,IG,IL)+ &
          & ((1./24.)*DXM(IG,IL)+(1./24.)*DXP(IG,IL))
          FHOMP_IG=FUNKNO(1,IG,IL)+(5.0*SQRT(3.)/6.0)*FUNKNO(2,IG,IL)+SQRT(5.)*FUNKNO(3,IG,IL)- &
          & ((1./24.)*DXM(IG,IL)+(1./24.)*DXP(IG,IL))
        ELSE
          CALL XABORT('BRERTS: IELEM OVERFLOW.')
        ENDIF
      ELSE
        CALL XABORT('BRERTS: ICOL OVERFLOW.')
      ENDIF
      FHOMM(IG,IL)=REAL(FHOMM_IG)
      FHOMP(IG,IL)=REAL(FHOMP_IG)
    ENDDO
    IF(IMPX.GT.0) THEN
      DO IL=1,NLF/2
        WRITE(6,'(14H BRERTS: ORDER,I2,38H RAVIART-THOMAS FLUX UNKNOWNS IN GROUP,I4,1H:/ &
        & (1P,12E12.4))') 2*IL-2,IG,FUNKNO(:IELEM,IG,IL)
        WRITE(6,'(14H BRERTS: ORDER,I2,40H RAVIART-THOMAS BOUNDARY FLUXES IN GROUP,I4,1H:/ &
        & (1P,12E12.4))') 2*IL-2,IG,FHOMM(IG,IL),FHOMP(IG,IL)
      ENDDO
    ENDIF
  ENDDO
  DEALLOCATE(DXP,DXM,FUNKNO)
END SUBROUTINE BRERTS
