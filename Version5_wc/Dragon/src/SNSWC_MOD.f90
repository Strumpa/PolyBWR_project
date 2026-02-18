!-----------------------------------------------------------------------
!
! Purpose:
!   Compute the flux solution along a direction, over the spatial domain
!   in 1D/2D/3D Cartesian geometry, within the discrete ordinates (SN)
!   framework. Boltzmann and Boltzmann Fokker-Planck solvers. Albedo 
!   boundary conditions.
!
! Copyright:
!  Copyright (C) 2025 Ecole Polytechnique de Montreal
!  This library is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License as published by the Free Software Foundation; either
!  version 2.1 of the License, or (at your option) any later version
!
! Author(s): A. Hebert, A. A. Calloo and C. Bienvenue
!
!---------------------------------------------------------------------
!
MODULE SNSWC_MOD

  USE SNBTE_MOD ! Boltzmann
  USE SNBFP_MOD ! Boltzmann Fokker-Planck

CONTAINS

  SUBROUTINE SNSWC1(NMAT,LX,NM,NMX,NSCT,ISCHM,IELEM,MAT,DX,SIGMA,V,CST,WX,BS,QEXT,FLUX,XNI,MN,DN,IS,ISBS, &
  & ISBSM,MAXL,NBS,LFIXUP,LSHOOT,IBFP,EELEM,NME,ESCHM,ESTOPW,DELTAE,WE,FEP,FEP0)
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    !   Compute the flux solution along a direction, over the spatial domain
    !   in 1D Cartesian geometry, within the discrete ordinates (SN)
    !   framework. Boltzmann and Boltzmann Fokker-Planck solvers. Albedo 
    !   boundary conditions.
    !
    ! Parameters: input
    !   NMAT    total number of materials.
    !   LX      total number of spatial cells.
    !   NM      total number of moments of the flux in space and in energy.
    !   NMX     number of incoming/outgoing boundary flux moments in space.
    !   NSCT    total number of scattering cross-section moments.
    !   ISCHM   spatial discretization scheme index.
    !   IELEM   (order+1) of the spatial approximation polynomial.
    !   MAT     material index in each spatial cell.
    !   DX      factor containing first direction cosines (mu).
    !   SIGMA   macroscopic total cross section in the cell.
    !   V       cell volume.
    !   CST     Legendre coefficients for the polynomial approximations.
    !   WX      spatial closure relation weighting factors.
    !   BS      intensities of boundary fixed sources.
    !   QEXT    angular external source term.
    !   MN      Moment-to-discrete matrix.
    !   DN      Discrete-to-moment matrix.
    !   IS      sweeping iteration.
    !   ISBS    flag for boundary sources sources.
    !   ISBSM   flag array for boundary fixed sources in each unit surface.
    !   MAXL    maximum number of unit surface in the boundary source.
    !   NBS     total number of boundary sources.
    !   LFIXUP  flag to enable negative flux fixup.
    !   LSHOOT  flag to enable shooting method.
    !   IBFP    type of energy proparation relation.
    !   EELEM   (order+1) of the energy approximation polynomial.
    !   NME     number of incoming/outgoing boundary flux moments in energy.
    !   ESCHM   energy discretization scheme index.
    !   ESTOPW  stopping powers at the upper and lower group boundaries.
    !   DELTAE  energy group width.
    !   WE      energy closure relation weighting factors.
    !   FEP0    initial incoming/outgoing boundary flux in energy.
    !
    ! Parameters: input/output
    !   FLUX    Legendre moments of the flux in space and in energy.
    !   XNI     incoming boundary flux in space.
    !   FEP     outgoing boundary flux in energy.
    !
    !-----------------------------------------------------------------------

    ! VARIABLE DECLARATION
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: LFIXUP,LSHOOT
    INTEGER, INTENT(IN) :: LX,NM,NMX,NSCT,ISCHM,IELEM,MAT(LX),NMAT,MAXL,NBS,ISBS,ISBSM(2),IBFP,EELEM,NME,ESCHM,IS
    REAL, INTENT(IN) :: DX,SIGMA(0:NMAT),V(LX),CST(IELEM),WX(IELEM+1),BS(MAXL,NBS),MN(NSCT),DN(NSCT),QEXT(NM,NSCT,LX)
    DOUBLE PRECISION, INTENT(INOUT) :: XNI(NMX)
    REAL, INTENT(INOUT) :: FLUX(NM,NSCT,LX)
    REAL, INTENT(IN), OPTIONAL :: ESTOPW(0:NMAT,2),DELTAE,WE(EELEM+1),FEP0(NME,LX)
    DOUBLE PRECISION, INTENT(INOUT), OPTIONAL :: FEP(NME,LX)

    ! LOCAL VARIABLES
    REAL :: SX
    INTEGER :: I0,I,IEL,IBM,P
    DOUBLE PRECISION :: QN(NM),Q2(NM,NM+1)

    ! SWEEPING DIRECTION
    SX = SIGN(1.0,DX)

    ! BOUNDARY FIXED SOURCES
    IF(DX.GT.0) THEN
        IF(ISBS.EQ.1.AND.ISBSM(1).NE.0) THEN
            XNI(1)=XNI(1)+BS(1,ISBSM(1))
        ENDIF
    ELSE
        IF(ISBS.EQ.1.AND.ISBSM(2).NE.0) THEN
            XNI(1)=XNI(1)+BS(1,ISBSM(2))
        ENDIF
    ENDIF

    ! SWEEPING OVER ALL VOXELS
    DO I0=1,LX
        I=I0
        IF(DX.LT.0) I=LX+1-I

        ! DATA
        IBM=MAT(I)
        IF(IBM.EQ.0) CYCLE

        ! SOURCE DENSITY TERM
        DO IEL=1,NM
            QN(IEL)=0.0
            DO P=1,NSCT
                QN(IEL)=QN(IEL)+QEXT(IEL,P,I)*MN(P)
            ENDDO
        ENDDO

        ! UPPER ENERGY GROUP BOUNDARY INCIDENT FLUX
        IF(IBFP.NE.0) FEP(:,I)=FEP0(:,I)

        ! FLUX CALCULATION IN THE VOXEL
        Q2=0.0D0
        IF(IBFP.EQ.0) THEN
            ! BOLTZMANN SOLVER
            CALL SNBTE1(DX,SIGMA(IBM),V(I),QN,Q2,XNI(1),NM,ISCHM,IELEM,CST,WX,SX,LFIXUP)
        ELSE
            ! BOLTZMANN FOKKER-PLANCK SOLVER
            CALL SNBFP1(DX,SIGMA(IBM),ESTOPW(IBM,1)/DELTAE,ESTOPW(IBM,2)/DELTAE,V(I),QN,Q2,XNI,FEP(:,I), &
            & NM,NMX,NME,ISCHM,ESCHM,IELEM,EELEM,CST,WX,WE,SX,LFIXUP,IBFP,DELTAE)
        ENDIF
        IF(LSHOOT.AND.IS.LT.5) CYCLE ! SHOOTING METHOD

        ! SAVE LEGENDRE MOMENT OF THE FLUX
        DO P=1,NSCT
            DO IEL=1,NM
                FLUX(IEL,P,I)=FLUX(IEL,P,I)+REAL(Q2(IEL,NM+1))*DN(P)
            ENDDO
        ENDDO
    ENDDO ! END OF X-LOOP

    RETURN
  END SUBROUTINE SNSWC1

  SUBROUTINE SNSWC2(NMAT,IXM,IYM,IXP,IYP,LX,LY,NM,NMX,NMY,NSCT,ISCHM,IELEM,MAT,DX,DY,SIGMA,V,CST,WX,WY,BS,QEXT, &
  & FLUX,XNI0,XNJ0,MN,DN,ZCODE,ISBS,ISBSM,MAXL,NBS,LFIXUP,IBFP,EELEM,NME,ESCHM,ESTOPW,DELTAE,WE,FEP,FEP0)
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    !   Compute the flux solution along a direction, over the spatial domain
    !   in 2D Cartesian geometry, within the discrete ordinates (SN)
    !   framework. Boltzmann and Boltzmann Fokker-Planck solvers. Albedo 
    !   boundary conditions.
    !
    !
    ! Parameters: input
    !   NMAT    total number of materials.
    !   IXM     starting index for sweeping over x-axis.
    !   IYM     starting index for sweeping over y-axis.
    !   IXP     ending index for sweeping over x-axis.
    !   IYP     ending index for sweeping over y-axis.
    !   LX      total number of spatial cells.
    !   LY      total number of spatial cells.
    !   NM      total number of moments of the flux in space and in energy.
    !   NMX     number of incoming/outgoing boundary flux moments in space.
    !   NMY     number of incoming/outgoing boundary flux moments in space.
    !   NSCT    total number of scattering cross-section moments.
    !   ISCHM   spatial discretization scheme index.
    !   IELEM   (order+1) of the spatial approximation polynomial.
    !   MAT     material index in each spatial cell.
    !   DX      factor containing first direction cosines (mu).
    !   DY      factor containing second direction cosines (eta).
    !   SIGMA   macroscopic total cross section in the cell.
    !   V       cell volume.
    !   CST     Legendre coefficients for the polynomial approximations.
    !   WX      spatial closure relation weighting factors.
    !   WY      spatial closure relation weighting factors.
    !   BS      intensities of boundary fixed sources.
    !   QEXT    angular external source term.
    !   MN      Moment-to-discrete matrix.
    !   DN      Discrete-to-moment matrix.
    !   ZCODE   albedo boundary condition coefficients.
    !   ISBS    flag for boundary sources sources.
    !   ISBSM   flag array for boundary fixed sources in each unit surface.
    !   MAXL    maximum number of unit surface in the boundary source.
    !   NBS     total number of boundary sources.
    !   LFIXUP  flag to enable negative flux fixup.
    !   IBFP    type of energy proparation relation.
    !   EELEM   (order+1) of the energy approximation polynomial.
    !   NME     number of incoming/outgoing boundary flux moments in energy.
    !   ESCHM   energy discretization scheme index.
    !   ESTOPW  stopping powers at the upper and lower group boundaries.
    !   DELTAE  energy group width.
    !   WE      energy closure relation weighting factors.
    !   FEP0    initial incoming boundary flux in energy.
    !
    ! Parameters: input/output
    !   FLUX    Legendre moments of the flux in space and in energy.
    !   XNI0    incoming/outgoing boundary flux in space.
    !   XNJ0    incoming/outgoing boundary flux in space.
    !   FEP     outgoing boundary flux in energy.
    !
    !-----------------------------------------------------------------------

    ! VARIABLE DECLARATION
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: LFIXUP
    INTEGER, INTENT(IN) :: IXM,IYM,IXP,IYP,LX,LY,NM,NMX,NMY,NSCT,ISCHM,IELEM,MAT(LX,LY),NMAT,MAXL,NBS,ISBS, &
                           & ISBSM(4),IBFP,EELEM,NME,ESCHM
    REAL, INTENT(IN) :: DX(LX,LY),DY(LX),SIGMA(0:NMAT),V(LX,LY),CST(IELEM),WX(IELEM+1),WY(IELEM+1),BS(MAXL,NBS), &
                        & MN(NSCT),DN(NSCT),ZCODE(6),QEXT(NM,NSCT,LX,LY)
    REAL, INTENT(INOUT) :: XNI0(NMX,LY),XNJ0(NMY,LX),FLUX(NM,NSCT,LX,LY)
    REAL, INTENT(IN), OPTIONAL :: ESTOPW(0:NMAT,2),DELTAE,WE(EELEM+1),FEP0(NME,LX,LY)
    REAL, INTENT(INOUT), OPTIONAL :: FEP(NME,LX,LY)

    ! LOCAL VARIABLES
    REAL :: SX,SY
    INTEGER :: I0,J0,I,J,P,IEL,IBM
    DOUBLE PRECISION :: XNI(NMX,LY),XNJ(NMY),QN(NM),Q2(NM,NM+1)

    ! SWEEPING DIRECTIONS
    SX = SIGN(1.0,DX(1,1))
    SY = SIGN(1.0,DY(1))

    ! SWEEPING OVER X-AXIS
    DO I0=IXM,IXP
        I=I0
        IF(SX.LT.0) I=LX+1-I

        ! SWEEPING OVER Y-AXIS
        DO J0=IYM,IYP
            J=J0
            IF(SY.LT.0) J=LY+1-J

            ! X-BOUNDARIES CONDITIONS
            IF(I0.EQ.1) THEN
                XNI(1:NMX,J)=0.0
                DO IEL=1,NMX
                    IF(SX.GT.0) THEN
                        XNI(IEL,J)=XNI0(IEL,J)*ZCODE(1)
                    ELSE
                        XNI(IEL,J)=XNI0(IEL,J)*ZCODE(2)
                    ENDIF
                ENDDO

                ! X-BOUNDARIES FIXED SOURCES
                IF(ISBS.EQ.1) THEN
                    IF((SX.LT.0).AND.ISBSM(2).NE.0) THEN
                        XNI(1,J)=XNI(1,J)+BS(J,ISBSM(2)) 
                    ELSEIF((SX.GT.0).AND.ISBSM(1).NE.0) THEN
                        XNI(1,J)=XNI(1,J)+BS(J,ISBSM(1)) 
                    ENDIF
                ENDIF
            ELSEIF(I0.EQ.IXM) THEN
                DO IEL=1,NMX
                    XNI(IEL,J)=XNI0(IEL,J)
                ENDDO
            ENDIF

            ! Y-BOUNDARIES CONDITIONS
            IF(J0.EQ.1) THEN
                XNJ(1:NMY)=0.0
                DO IEL=1,NMY
                    IF(SY.GT.0) THEN
                        XNJ(IEL)=XNJ0(IEL,I)*ZCODE(3)
                    ELSE
                        XNJ(IEL)=XNJ0(IEL,I)*ZCODE(4)
                    ENDIF
                ENDDO

                !Y-BOUNDARIES FIXED SOURCES
                IF(ISBS.EQ.1) THEN
                    IF((SY.LT.0).AND.ISBSM(4).NE.0) THEN
                        XNJ(1)=XNJ(1)+BS(I,ISBSM(4)) 
                    ELSEIF((SY.GT.0).AND.ISBSM(3).NE.0) THEN
                        XNJ(1)=XNJ(1)+BS(I,ISBSM(3))
                    ENDIF
                ENDIF
            ELSEIF(J0.EQ.IYM) THEN
                DO IEL=1,NMY
                    XNJ(IEL)=XNJ0(IEL,I)
                ENDDO
            ENDIF

            ! CHECK IF THE CELL IS VOID
            IBM=MAT(I,J)
            IF(IBM.EQ.0) CYCLE

            ! SOURCE DENSITY TERM
            QN(:NM)=0.0D0
            DO P=1,NSCT
                DO IEL=1,NM
                    QN(IEL)=QN(IEL)+QEXT(IEL,P,I,J)*MN(P)
                ENDDO
            ENDDO

            ! UPPER ENERGY GROUP BOUNDARY INCIDENT FLUX
            IF(IBFP.NE.0) FEP(:,I,J)=FEP0(:,I,J)

            ! FLUX CALCULATION
            Q2(:NM,:NM+1)=0.0D0
            IF(IBFP.EQ.0) THEN
                ! BOLTZMANN SOLVER
                CALL SNBTE2(DX(I,J),DY(I),SIGMA(IBM),V(I,J),QN,Q2,XNI(1:NMX,J),XNJ(1:NMY), &
                & NM,NMX,NMY,ISCHM,IELEM,CST,WX,WY,SX,SY,LFIXUP)
            ELSE
                ! BOLTZMANN FOKKER-PLANCK SOLVER
                CALL SNBFP2(DX(I,J),DY(I),SIGMA(IBM),ESTOPW(IBM,1)/DELTAE,ESTOPW(IBM,2)/DELTAE,V(I,J),QN,Q2,XNI(1:NMX,J), &
                & XNJ(1:NMY),FEP(:,I,J),NM,NMX,NMY,NME,ISCHM,ESCHM,IELEM,EELEM,CST,WX,WY,WE,SX,SY,LFIXUP,IBFP,DELTAE)
            ENDIF

            ! SAVE LEGENDRE MOMENT OF THE FLUX
            DO P=1,NSCT
            DO IEL=1,NM
                FLUX(IEL,P,I,J)=FLUX(IEL,P,I,J)+REAL(Q2(IEL,NM+1))*DN(P)
            ENDDO
            ENDDO

            ! SAVE Y-BOUNDARY CONDITIONS
            DO IEL=1,NMY
                XNJ0(IEL,I)=REAL(XNJ(IEL))
            ENDDO
        ENDDO

        ! SAVE X-BOUNDARY CONDITIONS
        DO J0=IYM,IYP
            J=J0
            IF(SY.LT.0) J=LY+1-J
            DO IEL=1,NMX
                XNI0(IEL,J)=REAL(XNI(IEL,J))
            ENDDO
        ENDDO

    ENDDO
    RETURN
  END SUBROUTINE SNSWC2

  SUBROUTINE SNSWC3(NMAT,IXM,IYM,IZM,IXP,IYP,IZP,LX,LY,LZ,NM,NMX,NMY,NMZ,NSCT,ISCHM,IELEM,MAT,DX,DY,DZ,SIGMA,V,CST,WX,WY,WZ, &
  & BS,QEXT,FLUX,XNI0,XNJ0,XNK0,MN,DN,ZCODE,ISBS,ISBSM,MAXL,NBS,LFIXUP,IBFP,EELEM,NME,ESCHM,ESTOPW,DELTAE,WE,FEP,FEP0)
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    !   Compute the flux solution along a direction, over the spatial domain
    !   in 3D Cartesian geometry, within the discrete ordinates (SN)
    !   framework. Boltzmann and Boltzmann Fokker-Planck solvers. Albedo 
    !   boundary conditions.
    !
    ! Parameters: input
    !   NMAT    total number of materials.
    !   IXM     starting index for sweeping over x-axis.
    !   IYM     starting index for sweeping over y-axis.
    !   IZM     starting index for sweeping over z-axis.
    !   IXP     ending index for sweeping over x-axis.
    !   IYP     ending index for sweeping over y-axis.
    !   IZP     ending index for sweeping over z-axis.
    !   LX      total number of spatial cells.
    !   LY      total number of spatial cells.
    !   LZ      total number of spatial cells.
    !   NM      total number of moments of the flux in space and in energy.
    !   NMX     number of incoming/outgoing boundary flux moments in space.
    !   NMY     number of incoming/outgoing boundary flux moments in space.
    !   NMZ     number of incoming/outgoing boundary flux moments in space.
    !   NSCT    total number of scattering cross-section moments.
    !   ISCHM   spatial discretization scheme index.
    !   IELEM   (order+1) of the spatial approximation polynomial.
    !   MAT     material index in each spatial cell.
    !   DX      factor containing first direction cosines (mu).
    !   DY      factor containing second direction cosines (eta).
    !   DZ      factor containing third direction cosines (xi).
    !   SIGMA   macroscopic total cross section in the cell.
    !   V       cell volume.
    !   CST     Legendre coefficients for the polynomial approximations.
    !   WX      spatial closure relation weighting factors.
    !   WY      spatial closure relation weighting factors.
    !   WZ      spatial closure relation weighting factors.
    !   BS      intensities of boundary fixed sources.
    !   QEXT    angular external source term.
    !   MN      Moment-to-discrete matrix.
    !   DN      Discrete-to-moment matrix.
    !   ZCODE   albedo boundary condition coefficients.
    !   ISBS    flag for boundary sources sources.
    !   ISBSM   flag array for boundary fixed sources in each unit surface.
    !   MAXL    maximum number of unit surface in the boundary source.
    !   NBS     total number of boundary sources.
    !   LFIXUP  flag to enable negative flux fixup.
    !   IBFP    type of energy proparation relation.
    !   EELEM   (order+1) of the energy approximation polynomial.
    !   NME     number of incoming/outgoing boundary flux moments in energy.
    !   ESCHM   energy discretization scheme index.
    !   ESTOPW  stopping powers at the upper and lower group boundaries.
    !   DELTAE  energy group width.
    !   WE      energy closure relation weighting factors.
    !   FEP0    initial incoming boundary flux in energy.
    !
    ! Parameters: input/output
    !   FLUX    Legendre moments of the flux in space and in energy.
    !   XNI0    incoming/outgoing boundary flux in space.
    !   XNJ0    incoming/outgoing boundary flux in space.
    !   XNK0    incoming/outgoing boundary flux in space.
    !   FEP     outgoing boundary flux in energy.
    !
    !-----------------------------------------------------------------------

    ! VARIABLE DECLARATION
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: LFIXUP
    INTEGER, INTENT(IN) :: IXM,IYM,IZM,IXP,IYP,IZP,LX,LY,LZ,NM,NMX,NMY,NMZ,NSCT,ISCHM,IELEM,MAT(LX,LY,LZ),NMAT,MAXL, &
    & NBS,ISBS,ISBSM(6),IBFP,EELEM,NME,ESCHM
    REAL, INTENT(IN) :: DX(LY,LZ),DY(LX,LZ),DZ(LX,LY),SIGMA(0:NMAT),V(LX,LY,LZ),CST(IELEM),WX(IELEM+1),WY(IELEM+1), &
    & WZ(IELEM+1),BS(MAXL,NBS),MN(NSCT),DN(NSCT),ZCODE(6), QEXT(NM,NSCT,LX,LY,LZ)
    REAL, INTENT(INOUT) :: XNI0(NMX,LY,LZ),XNJ0(NMY,LX,LZ),XNK0(NMZ,LX,LY),FLUX(NM,NSCT,LX,LY,LZ)
    REAL, INTENT(IN), OPTIONAL :: ESTOPW(0:NMAT,2),DELTAE,WE(EELEM+1),FEP0(NME,LX,LY,LZ)
    REAL, INTENT(INOUT), OPTIONAL :: FEP(NME,LX,LY,LZ)

    ! LOCAL VARIABLES
    REAL :: SX,SY,SZ
    INTEGER :: I0,J0,K0,I,J,K,P,IEL,IBM
    DOUBLE PRECISION :: XNI(NMX,LY,LZ),XNJ(NMY,LZ),XNK(NMZ),QN(NM),Q2(NM,NM+1)

    ! SWEEPING DIRECTIONS
    SX = SIGN(1.0,DX(1,1))
    SY = SIGN(1.0,DY(1,1))
    SZ = SIGN(1.0,DZ(1,1))

    ! SWEEPING OVER X-AXIS
    DO I0=IXM,IXP
        I=I0
        IF(SX.LT.0) I=LX+1-I

        ! SWEEPING OVER Y-AXIS
        DO J0=IYM,IYP
            J=J0
            IF(SY.LT.0) J=LY+1-J

            ! SWEEPING OVER Z-AXIS
            DO K0=IZM,IZP
                K=K0
                IF(SZ.LT.0) K=LZ+1-K

                ! X-BOUNDARIES CONDITIONS
                IF(I0.EQ.1) THEN
                    XNI(:NMX,J,K)=0.0
                    DO IEL=1,NMX
                        IF(SX.GT.0) THEN
                            XNI(IEL,J,K)=XNI0(IEL,J,K)*ZCODE(1)
                        ELSE
                            XNI(IEL,J,K)=XNI0(IEL,J,K)*ZCODE(2)
                        ENDIF
                    ENDDO

                    ! X-BOUNDARIES FIXED SOURCES
                    IF(ISBS.EQ.1) THEN
                        IF((SX.LT.0).AND.ISBSM(2).NE.0) THEN
                            XNI(1,J,K)=XNI(1,J,K)+BS((J-1)*LZ+K,ISBSM(2)) 
                        ELSEIF((SX.GT.0).AND.ISBSM(1).NE.0) THEN
                            XNI(1,J,K)=XNI(1,J,K)+BS((J-1)*LZ+K,ISBSM(1)) 
                        ENDIF
                    ENDIF
                ELSEIF(I0.EQ.IXM) THEN
                    DO IEL=1,NMX
                        XNI(IEL,J,K)=XNI0(IEL,J,K)
                    ENDDO
                ENDIF

                ! Y-BOUNDARIES CONDITIONS
                IF(J0.EQ.1) THEN
                    XNJ(:NMY,K)=0.0
                    DO IEL=1,NMY
                        IF(SY.GT.0) THEN
                            XNJ(IEL,K)=XNJ0(IEL,I,K)*ZCODE(3)
                        ELSE
                            XNJ(IEL,K)=XNJ0(IEL,I,K)*ZCODE(4)
                        ENDIF
                    ENDDO

                    !Y-BOUNDARIES FIXED SOURCES
                    IF(ISBS.EQ.1) THEN
                        IF((SY.LT.0).AND.ISBSM(4).NE.0) THEN
                            XNJ(1,K)=XNJ(1,K)+BS((I-1)*LZ+K,ISBSM(4)) 
                        ELSEIF((SY.GT.0).AND.ISBSM(3).NE.0) THEN
                            XNJ(1,K)=XNJ(1,K)+BS((I-1)*LZ+K,ISBSM(3)) 
                        ENDIF
                    ENDIF
                ELSEIF(J0.EQ.IYM) THEN
                    DO IEL=1,NMY
                        XNJ(IEL,K)=XNJ0(IEL,I,K)
                    ENDDO
                ENDIF

                ! Z-BOUNDARIES CONDITIONS
                IF(K0.EQ.1) THEN
                    XNK(:NMZ)=0.0
                    DO IEL=1,NMZ
                        IF(SZ.GT.0) THEN
                            XNK(IEL)=XNK0(IEL,I,J)*ZCODE(5)
                        ELSE
                            XNK(IEL)=XNK0(IEL,I,J)*ZCODE(6)
                        ENDIF
                    ENDDO

                    ! Z-BOUNDARIES FIXED SOURCES
                    IF(ISBS.EQ.1) THEN
                        IF((SZ.LT.0).AND.ISBSM(6).NE.0) THEN
                            XNK(1)=XNK(1)+BS((I-1)*LY+J,ISBSM(6)) 
                        ELSEIF((SZ.GT.0).AND.ISBSM(5).NE.0) THEN
                            XNK(1)=XNK(1)+BS((I-1)*LY+J,ISBSM(5)) 
                        ENDIF
                    ENDIF
                ELSEIF(K0.EQ.IZM) THEN
                    DO IEL=1,NMZ
                        XNK(IEL)=XNK0(IEL,I,J)
                    ENDDO
                ENDIF

                ! CHECK IF THE CELL IS VOID
                IBM=MAT(I,J,K)
                IF(IBM.EQ.0) CYCLE

                ! SOURCE DENSITY TERM
                QN(:NM)=0.0D0
                DO P=1,NSCT
                DO IEL=1,NM
                    QN(IEL)=QN(IEL)+QEXT(IEL,P,I,J,K)*MN(P)
                ENDDO
                ENDDO

                ! UPPER ENERGY GROUP BOUNDARY INCIDENT FLUX
                IF(IBFP.NE.0) FEP(:NME,I,J,K)=FEP0(:NME,I,J,K)

                ! FLUX CALCULATION
                Q2(:NM,:NM+1)=0.0D0
                IF(IBFP.EQ.0) THEN
                    ! BOLTZMANN SOLVER
                    CALL SNBTE3(DX(J,K),DY(I,K),DZ(I,J),SIGMA(IBM),V(I,J,K),QN,Q2,XNI(:NMX,J,K),XNJ(:NMY,K),XNK(:NMZ), &
                    & NM,NMX,NMY,NMZ,ISCHM,IELEM,CST,WX,WY,WZ,SX,SY,SZ,LFIXUP)
                ELSE
                    ! BOLTZMANN FOKKER-PLANCK SOLVER
                    CALL SNBFP3(DX(J,K),DY(I,K),DZ(I,J),SIGMA(IBM),ESTOPW(IBM,1)/DELTAE,ESTOPW(IBM,2)/DELTAE,V(I,J,K),QN,Q2, &
                    & XNI(:NMX,J,K),XNJ(:NMY,K),XNK(:NMZ),FEP(:NME,I,J,K),NM,NMX,NMY,NMZ,NME,ISCHM,ESCHM,IELEM,EELEM,CST,WX,WY, &
                    & WZ,WE,SX,SY,SZ,LFIXUP,IBFP,DELTAE)
                ENDIF

                ! SAVE LEGENDRE MOMENT OF THE FLUX
                DO P=1,NSCT
                DO IEL=1,NM
                    FLUX(IEL,P,I,J,K)=FLUX(IEL,P,I,J,K)+REAL(Q2(IEL,NM+1))*DN(P)
                ENDDO
                ENDDO
            ENDDO

            ! SAVE Z-BOUNDARY CONDITIONS
            DO IEL=1,NMZ
                XNK0(IEL,I,J)=REAL(XNK(IEL))
            ENDDO
        ENDDO

        ! SAVE Y-BOUNDARY CONDITIONS
        DO K0=IZM,IZP
            K=K0
            IF(SZ.LT.0) K=LZ+1-K
            DO IEL=1,NMY
                XNJ0(IEL,I,K)=REAL(XNJ(IEL,K))
            ENDDO
        ENDDO
    ENDDO

    ! SAVE X-BOUNDARY CONDITIONS
    DO K0=IZM,IZP
        K=K0
        IF(SZ.LT.0) K=LZ+1-K
        DO J0=IYM,IYP
            J=J0
            IF(SY.LT.0) J=LY+1-J
            DO IEL=1,NMX
                XNI0(IEL,J,K)=REAL(XNI(IEL,J,K))
            ENDDO
        ENDDO
    ENDDO
    RETURN
  END SUBROUTINE SNSWC3
END MODULE SNSWC_MOD
