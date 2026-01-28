!-----------------------------------------------------------------------
!
! Purpose:
!   Compute the flux solution along a direction, over the spatial domain
!   in 2D/3D hexagonal geometry, within the discrete ordinates (SN)
!   framework. Boltzmann solvers. Albedo boundary conditions.
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
!-----------------------------------------------------------------------
!
MODULE SNSWH_MOD

  USE SNBTE_MOD

CONTAINS

  SUBROUTINE SNSWH2(IND,JAC,LOZSWP,ISPLH,TMPXNI,TMPXNJ,TMPXND,FLUX,MAT,TOTAL,VOL,DA,DB,QEXT,NUN,NM,NMX,NMY,NSCT, &
  & ISCHM,IELEM,CST,WX,WY,NCOL,NMAT,IHEX,I,J,DCOORD,LFIXUP,MN,DN)
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    !   Compute the flux solution along a direction, over the spatial domain
    !   in 2D hexagonal geometry, within the discrete ordinates (SN)
    !   framework. Boltzmann solvers. Albedo boundary conditions.
    !
    ! Parameters: input
    !   IND     angle index for sweeping direction.
    !   JAC     Jacobian matrix elements for the affine transformation.
    !   LOZSWP  lozenge sweeping order for the given angle.
    !   ISPLH   number of subcells in each dimension within a lozenge.
    !   MAT     material index in each spatial cell.
    !   TOTAL   macroscopic total cross section in each material.
    !   VOL     cell volume.
    !   DX      factor containing first direction cosines (mu).
    !   DY      factor containing second direction cosines (eta).
    !   QEXT    angular external source term.
    !   NUN     total number of unknowns.
    !   NM      total number of moments of the flux in space.
    !   NMX     number of incoming/outgoing boundary flux moments in space.
    !   NMY     number of incoming/outgoing boundary flux moments in space.
    !   NSCT    total number of scattering cross-section moments.
    !   ISCHM   spatial discretization scheme index.
    !   IELEM   (order+1) of the spatial approximation polynomial.
    !   CST     Legendre coefficients for the polynomial approximations.
    !   WX      spatial closure relation weighting factors.
    !   WY      spatial closure relation weighting factors.
    !   NCOL    number of columns in the temporary boundary flux arrays.
    !   NMAT    total number of materials.
    !   IHEX    hexagonal mesh index.
    !   I       I-index related to the current direction.
    !   J       J-index related to the current direction.
    !   DCOORD  D-index related to the current direction.
    !   LFIXUP  flag to enable negative flux fixup.
    !   MN      Moment-to-discrete matrix.
    !   DN      Discrete-to-moment matrix.
    !
    ! Parameters: input/output
    !   TMPXNI  incoming/outgoing boundary flux in space (I-axis).
    !   TMPXNJ  incoming/outgoing boundary flux in space (J-axis).
    !   TMPXND  incoming/outgoing boundary flux in space (D-axis).
    !   FLUX    Legendre moments of the flux in space and in energy.
    !
    ! Comments:
    !   1. The direction of the axes I, J and D for the surface boundary 
    !      fluxes are shown in the diagram below. This means that 
    !      a) lozenge A has I- and D-boundaries (instead of I and J)
    !      b) lozenge B has I- and J-boundaries
    !      c) lozenge C has D- and J-boundaries (instead of I and J)
    !
    !                                  ^
    !                         j-axis   |
    !      ^  y-axis                   |          ^
    !      |                       _________     /    d-axis
    !      |                      /       / \   /
    !      |                     /   B   /   \
    !      | - - - ->           /       /     \
    !           x-axis         (-------(   A   )
    !                           \       \     /
    !                            \  C    \   / 
    !                             \_______\_/   \
    !                                            \   i-axis
    !                                             ^
    !
    !-----------------------------------------------------------------------

    ! SUBROUTINE ARGUMENTS
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: LFIXUP
    INTEGER, INTENT(IN) :: IND,ISPLH,NUN,NM,NMX,NMY,NSCT,ISCHM,IELEM,NCOL,LOZSWP(3,6),NMAT,IHEX,I,J,DCOORD,MAT(ISPLH,ISPLH,3)
    REAL, INTENT(IN) :: JAC(2,2,3),WX(IELEM+1),WY(IELEM+1),CST(IELEM),MN(NSCT),DN(NSCT),TOTAL(0:NMAT),VOL(ISPLH,ISPLH,3), &
    & DA(ISPLH,ISPLH,3),DB(ISPLH,ISPLH,3),QEXT(NUN)
    DOUBLE PRECISION, INTENT(INOUT) :: FLUX(NM,NSCT,3*ISPLH**2),TMPXNI(IELEM,ISPLH,NCOL),TMPXNJ(IELEM,ISPLH,NCOL), &
    & TMPXND(IELEM,ISPLH,NCOL)

    ! LOCAL VARIABLES
    INTEGER :: ILOZLOOP,ILOZ,IL,JL,I2,J2,IY,IX,IBM,I_FETCH,IOF,IEL,P
    REAL :: SX,SY,SIGMA,V,MUH,ETAH,AAA,BBB,CCC,DDD
    DOUBLE PRECISION :: Q(NM),Q2(NM,NM+1),XNI(NMX),XNJ(NMY)

    ! LOOP OVER LOZENGES
    DO ILOZLOOP=1,3
        ILOZ=LOZSWP(ILOZLOOP,IND)

        ! Get Jacobian elements values
        AAA = JAC(1,1,ILOZ)
        BBB = JAC(1,2,ILOZ)
        CCC = JAC(2,1,ILOZ)
        DDD = JAC(2,2,ILOZ)

        ! LOOP OVER SUBMESH WITHIN EACH LOZENGE
        DO IL=1,ISPLH
            I2=IL
            ! Account for different sweep direction depending on angle 
            IF((ILOZ.EQ.1).OR.(ILOZ.EQ.2))THEN
                IF((IND.EQ.2).OR.(IND.EQ.3).OR.(IND.EQ.4)) I2=ISPLH+1-I2
            ELSEIF(ILOZ.EQ.3)THEN
                IF((IND.EQ.3).OR.(IND.EQ.4).OR.(IND.EQ.5)) I2=ISPLH+1-I2
            ENDIF
            DO JL=1,ISPLH
                J2=JL
                ! Account for different sweep direction depending on angle 
                IF((ILOZ.EQ.2).OR.(ILOZ.EQ.3))THEN
                    IF((IND.EQ.4).OR.(IND.EQ.5).OR.(IND.EQ.6)) J2=ISPLH+1-J2
                ELSEIF(ILOZ.EQ.1)THEN
                    IF((IND.EQ.3).OR.(IND.EQ.4).OR.(IND.EQ.5)) J2=ISPLH+1-J2
                ENDIF

                ! READ IN XNI AND XNJ DEPENDING ON LOZENGE
                I_FETCH=0

                IF((ILOZ.EQ.1))THEN
                    ! Read boundary fluxes in reverse for lozenge A since affine transformation of lozenges causes the D and
                    ! I directions of lozenges C and A respectively to be reversed
                    I_FETCH=ISPLH+1-I2
                    XNI(:) = TMPXNI(:,J2,J)
                    XNJ(:) = TMPXND(:,I_FETCH,DCOORD)
                ELSEIF((ILOZ.EQ.2))THEN
                    XNI(:) = TMPXNI(:,J2,J)
                    XNJ(:) = TMPXNJ(:,I2,I)
                ELSEIF((ILOZ.EQ.3))THEN
                    XNI(:) = TMPXND(:,J2,DCOORD)
                    XNJ(:) = TMPXNJ(:,I2,I)
                ENDIF

                ! DATA
                IBM=MAT(I2,J2,ILOZ)
                IF(IBM.EQ.0) CYCLE ! Skip loop if virtual element 
                SIGMA=TOTAL(IBM)
                V=VOL(I2,J2,ILOZ)

                ! COMPUTE EFFECTIVE DIRECTION COSINES
                MUH=(DA(I2,J2,ILOZ)*DDD)-(DB(I2,J2,ILOZ)*BBB)
                ETAH=(-DA(I2,J2,ILOZ)*CCC)+(DB(I2,J2,ILOZ)*AAA)

                ! SOURCE DENSITY TERM
                DO IEL=1,NM
                    Q(IEL)=0.0D0
                    DO P=1,NSCT
                        Q(IEL)=Q(IEL)+QEXT(((((((IHEX-1)*3+(ILOZ-1))*ISPLH+(J2-1))*ISPLH+(I2-1))*NSCT+(P-1))*NM)+IEL)*MN(P)
                    ENDDO
                ENDDO

                ! FLUX CALCULATION - CARTESIAN BOLTZMANN SOLVER
                Q2(:NM,:NM+1)=0.0D0
                SX=SIGN(1.0,MUH)
                SY=SIGN(1.0,ETAH)
                CALL SNBTE2(MUH,ETAH,SIGMA,REAL(V),Q,Q2,XNI,XNJ,NM,NMX,NMY,ISCHM,IELEM,CST,WX,WY,SX,SY,LFIXUP)

                ! Assign I-boundary fluxes if lozenges A or B
                IF((ILOZ.EQ.1).OR.(ILOZ.EQ.2)) THEN
                    TMPXNI(:,J2,J)=REAL(XNI(:))
                ENDIF

                ! Assign J-boundary fluxes if lozenges B or C (write back to J-boundary)
                IF((ILOZ.EQ.2).OR.(ILOZ.EQ.3)) THEN
                    TMPXNJ(:,I2,I)=REAL(XNJ(:))
                ENDIF

                ! Assign D-boundary fluxes if lozenge A or C (map from the other face)
                IF((ILOZ.EQ.1)) THEN
                    TMPXND(:,I_FETCH,DCOORD)=REAL(XNJ(:))
                ELSEIF((ILOZ.EQ.3)) THEN
                    TMPXND(:,J2,DCOORD)=REAL(XNI(:))
                ENDIF

                ! Flip gradient components on D-boundary when orientation reverses
                DO IY=1,IELEM
                    IF((MOD(IY,2).EQ.0).AND.(ILOZ.EQ.3).AND.(IL.EQ.ISPLH)) TMPXND(IY,J2,DCOORD)=TMPXND(IY,J2,DCOORD)*(-1)
                ENDDO
                DO IX=1,IELEM
                    IF((MOD(IX,2).EQ.0).AND.(ILOZ.EQ.1).AND.(JL.EQ.ISPLH)) TMPXND(IX,I_FETCH,DCOORD)=TMPXND(IX,I_FETCH,DCOORD)*(-1)
                ENDDO

                ! SAVE LEGENDRE MOMENT OF THE FLUX
                IOF=((ILOZ-1)*ISPLH+(J2-1))*ISPLH+I2 ! lozenge-local subcell index
                DO P=1,NSCT
                    DO IEL=1,NM
                        FLUX(IEL,P,IOF)=FLUX(IEL,P,IOF)+Q2(IEL,NM+1)*DN(P)
                    ENDDO
                ENDDO

            ENDDO ! END OF LOZANGE SUBMESH J-LOOP
        ENDDO ! END OF LOZANGE SUBMESH I-LOOP
    ENDDO ! END OF LOZENGES LOOP

    RETURN
  END SUBROUTINE SNSWH2

  SUBROUTINE SNSWH3(IND,IND_XY,JAC,LOZSWP,ISPLH,LZ,TMPXNI,TMPXNJ,TMPXND,TMPXNK,FLUX,MAT,TOTAL,VOL,DA,DB,DC,QEXT,NUN,NM, &
  & NMX,NMY,NMZ,NSCT,ISCHM,IELEM,CST,WX,WY,WZ,NCOL,NHEX,NMAT,IHEX_XY,I,J,DCOORD,LFIXUP,MN,DN,K_MC,NCELLZ)
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    !   Compute the flux solution along a direction, over the spatial domain
    !   in 3D hexagonal geometry (i.e. hexagonal in the xy-plane, cartesian 
    !   in the z-axis), within the discrete ordinates (SN) framework. 
    !   Boltzmann solvers. Albedo boundary conditions.
    !
    ! Parameters: input
    !   IND     angle index for sweeping direction.
    !   IND_XY  angle index in the XY-plane for sweeping direction.
    !   JAC     Jacobian matrix elements for the affine transformation.
    !   LOZSWP  lozenge sweeping order for the given angle.
    !   ISPLH   number of subcells in each dimension within a lozenge.
    !   LZ      total number of planes in the z-axis within a macrocell.
    !   MAT     material index in each spatial cell.
    !   TOTAL   macroscopic total cross section in each material.
    !   VOL     cell volume.
    !   DX      factor containing first direction cosines (mu).
    !   DY      factor containing second direction cosines (eta).
    !   DZ      factor containing third direction cosines (xi).
    !   QEXT    angular external source term.
    !   NUN     total number of unknowns.
    !   NM      total number of moments of the flux in space.
    !   NMX     number of incoming/outgoing boundary flux moments in space.
    !   NMY     number of incoming/outgoing boundary flux moments in space.
    !   NMZ     number of incoming/outgoing boundary flux moments in space.
    !   NSCT    total number of scattering cross-section moments.
    !   ISCHM   spatial discretization scheme index.
    !   IELEM   (order+1) of the spatial approximation polynomial.
    !   CST     Legendre coefficients for the polynomial approximations.
    !   WX      spatial closure relation weighting factors.
    !   WY      spatial closure relation weighting factors.
    !   WZ      spatial closure relation weighting factors.
    !   NCOL    number of columns in the temporary boundary flux arrays.
    !   NHEX    total number of hexagonal macrocells in the xy-plane.
    !   NMAT    total number of materials.
    !   IHEX_XY hexagonal mesh index in the xy-plane.
    !   I       I-index related to the current direction.
    !   J       J-index related to the current direction.
    !   DCOORD  D-index related to the current direction.
    !   LFIXUP  flag to enable negative flux fixup.
    !   MN      Moment-to-discrete matrix.
    !   DN      Discrete-to-moment matrix.
    !   K_MC    macrocell z-index.
    !   NCELLZ  number of z-axis planes in each macrocell.
    !
    ! Parameters: input/output
    !   TMPXNI  incoming/outgoing boundary flux in space (I-axis).
    !   TMPXNJ  incoming/outgoing boundary flux in space (J-axis).
    !   TMPXND  incoming/outgoing boundary flux in space (D-axis).
    !   TMPXNK  incoming/outgoing boundary flux in space (K-axis).
    !   FLUX    Legendre moments of the flux in space and in energy.
    !
    ! Comments:
    !   1. The direction of the axes I, J and D for the surface boundary 
    !      fluxes are shown in the diagram below. This means that 
    !      a) lozenge A has I- and D-boundaries (instead of I and J)
    !      b) lozenge B has I- and J-boundaries
    !      c) lozenge C has D- and J-boundaries (instead of I and J)
    !
    !                                  ^
    !                         j-axis   |
    !      ^  y-axis                   |          ^
    !      |                       _________     /    d-axis
    !      |                      /       / \   /
    !      |                     /   B   /   \
    !      | - - - ->           /       /     \
    !           x-axis         (-------(   A   )
    !                           \       \     /
    !                            \  C    \   / 
    !                             \_______\_/   \
    !                                            \   i-axis
    !                                             ^
    !
    !-----------------------------------------------------------------------

    ! SUBROUTINE ARGUMENTS
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: LFIXUP
    INTEGER, INTENT(IN) :: IND,IND_XY,ISPLH,NUN,NM,NMX,NMY,NMZ,NSCT,ISCHM,IELEM,NCOL,NHEX,LOZSWP(3,6),NMAT,IHEX_XY,I,J,DCOORD, &
    & MAT(ISPLH,ISPLH,3,NHEX,LZ),LZ,K_MC,NCELLZ
    REAL, INTENT(IN) :: JAC(2,2,3),WX(IELEM+1),WY(IELEM+1),WZ(IELEM+1),CST(IELEM),MN(NSCT),DN(NSCT),TOTAL(0:NMAT), &
    & VOL(ISPLH,ISPLH,3,NHEX,LZ),DC(ISPLH*ISPLH*3*NHEX,1),DB(ISPLH*ISPLH*3*NHEX,LZ),DA(1,LZ),QEXT(NUN)
    DOUBLE PRECISION, INTENT(INOUT) :: FLUX(NM,NSCT,3*ISPLH**2,NHEX,LZ),TMPXNI(NMX,ISPLH,NCOL,LZ),TMPXNJ(NMY,ISPLH,NCOL,LZ), &
    & TMPXND(NMX,ISPLH,NCOL,LZ),TMPXNK(NMZ,ISPLH,ISPLH,3)

    ! LOCAL VARIABLES
    INTEGER :: ILOZLOOP,ILOZ,IL,JL,I2,J2,IY,IX,IBM,I_FETCH,IOF,IEL,P,IMZ,K,IZ
    REAL :: SX,SY,SZ,SIGMA,V,MUH,ETAH,XIH,AAA,BBB,CCC,DDD
    DOUBLE PRECISION :: Q(NM),Q2(NM,NM+1),XNI(NMX),XNJ(NMY),XNK(NMZ)

    ! LOOP OVER Z-AXIS PLANES IN MACROCELL
    DO IMZ=1,MIN(NCELLZ,LZ-(K_MC-1)*NCELLZ)
        K=(K_MC-1)*NCELLZ+IMZ
        IF(IND.GE.7) K=LZ+1-K

        ! LOOP OVER LOZENGES
        DO ILOZLOOP=1,3
            ILOZ=LOZSWP(ILOZLOOP,IND_XY)

            ! Get Jacobian elements values
            AAA = JAC(1,1,ILOZ)
            BBB = JAC(1,2,ILOZ)
            CCC = JAC(2,1,ILOZ)
            DDD = JAC(2,2,ILOZ)

            ! LOOP OVER SUBMESH WITHIN EACH LOZENGE
            DO IL=1,ISPLH
                I2=IL
                ! Account for different sweep direction depending on angle 
                IF((ILOZ.EQ.1).OR.(ILOZ.EQ.2))THEN
                    IF((IND_XY.EQ.2).OR.(IND_XY.EQ.3).OR.(IND_XY.EQ.4)) I2=ISPLH+1-I2
                ELSEIF(ILOZ.EQ.3)THEN
                    IF((IND_XY.EQ.3).OR.(IND_XY.EQ.4).OR.(IND_XY.EQ.5)) I2=ISPLH+1-I2
                ENDIF
                DO JL=1,ISPLH
                    J2=JL
                    ! Account for different sweep direction depending on angle 
                    IF((ILOZ.EQ.2).OR.(ILOZ.EQ.3))THEN
                        IF((IND_XY.EQ.4).OR.(IND_XY.EQ.5).OR.(IND_XY.EQ.6)) J2=ISPLH+1-J2
                    ELSEIF(ILOZ.EQ.1)THEN
                        IF((IND_XY.EQ.3).OR.(IND_XY.EQ.4).OR.(IND_XY.EQ.5)) J2=ISPLH+1-J2
                    ENDIF

                    ! READ IN XNI AND XNJ DEPENDING ON LOZENGE
                    I_FETCH=0

                    IF((ILOZ.EQ.1))THEN
                        ! Read boundary fluxes in reverse for lozenge A since affine transformation of lozenges
                        ! causes the D and I directions of lozenges C and A respectively to be reversed
                        I_FETCH=ISPLH+1-I2
                        XNI(:) = TMPXNI(:,J2,J,K)
                        XNJ(:) = TMPXND(:,I_FETCH,DCOORD,K)
                    ELSEIF((ILOZ.EQ.2))THEN
                        XNI(:) = TMPXNI(:,J2,J,K)
                        XNJ(:) = TMPXNJ(:,I2,I,K)
                    ELSEIF((ILOZ.EQ.3))THEN
                        XNI(:) = TMPXND(:,J2,DCOORD,K)
                        XNJ(:) = TMPXNJ(:,I2,I,K)
                    ENDIF
                    XNK(:) = TMPXNK(:,I2,J2,ILOZ)

                    ! DATA
                    IBM=MAT(I2,J2,ILOZ,IHEX_XY,K)
                    IF(IBM.EQ.0) CYCLE ! Skip loop if virtual element 
                    SIGMA=TOTAL(IBM)
                    V=VOL(I2,J2,ILOZ,IHEX_XY,K)

                    ! COMPUTE EFFECTIVE DIRECTION COSINES
                    MUH=(DA(1,K)*DDD)-(DB(1,K)*BBB)
                    ETAH=(-DA(1,K)*CCC)+(DB(1,K)*AAA)
                    XIH=DC(1,1)

                    ! SOURCE DENSITY TERM
                    DO IEL=1,NM
                        Q(IEL)=0.0D0
                        DO P=1,NSCT
                            Q(IEL)=Q(IEL)+QEXT((((((((K-1)*NHEX+(IHEX_XY-1))*3+(ILOZ-1))*ISPLH+(J2-1))*ISPLH+(I2-1))* &
                            & NSCT+(P-1))*NM)+IEL)*MN(P)
                        ENDDO
                    ENDDO

                    ! FLUX CALCULATION - CARTESIAN BOLTZMANN SOLVER
                    Q2(:NM,:NM+1)=0.0D0
                    SX=SIGN(1.0,MUH)
                    SY=SIGN(1.0,ETAH)
                    SZ=SIGN(1.0,XIH)
                    CALL SNBTE3(MUH,ETAH,XIH,SIGMA,REAL(V),Q,Q2,XNI,XNJ,XNK,NM,NMX,NMY,NMZ, &
                    & ISCHM,IELEM,CST,WX,WY,WZ,SX,SY,SZ,LFIXUP)

                    ! Assign I-boundary fluxes if lozenges A or B
                    IF((ILOZ.EQ.1).OR.(ILOZ.EQ.2)) THEN
                        TMPXNI(:,J2,J,K)=REAL(XNI(:))
                    ENDIF

                    ! Assign J-boundary fluxes if lozenges B or C (write back to J-boundary)
                    IF((ILOZ.EQ.2).OR.(ILOZ.EQ.3)) THEN
                        TMPXNJ(:,I2,I,K)=REAL(XNJ(:))
                    ENDIF

                    ! Assign D-boundary fluxes if lozenge A or C (map from the other face)
                    IF((ILOZ.EQ.1)) THEN
                        TMPXND(:,I_FETCH,DCOORD,K)=REAL(XNJ(:))
                    ELSEIF((ILOZ.EQ.3)) THEN
                        TMPXND(:,J2,DCOORD,K)=REAL(XNI(:))
                    ENDIF

                    ! Z-SPACE
                    TMPXNK(:,I2,J2,ILOZ)=XNK(:)

                    ! Flip gradient components on D-boundary when orientation reverses
                    ! Match the exact logic used in SNFKH3
                    DO IZ=1,IELEM
                        DO IY=1,IELEM
                            IF((MOD(IY,2).EQ.0).AND.(ILOZ.EQ.3).AND.(IL.EQ.ISPLH)) THEN
                                TMPXND(IELEM*(IZ-1)+IY,J2,DCOORD,K)=TMPXND(IELEM*(IZ-1)+IY,J2,DCOORD,K)*(-1)
                            ENDIF
                        ENDDO
                    ENDDO
                    DO IZ=1,IELEM
                        DO IX=1,IELEM
                            IF((MOD(IX,2).EQ.0).AND.(ILOZ.EQ.1).AND.(JL.EQ.ISPLH)) THEN
                                TMPXND(IELEM*(IZ-1)+IX,I_FETCH,DCOORD,K)=TMPXND(IELEM*(IZ-1)+IX,I_FETCH,DCOORD,K)*(-1)
                            ENDIF
                        ENDDO
                    ENDDO

                    ! SAVE LEGENDRE MOMENT OF THE FLUX
                    IOF=((ILOZ-1)*ISPLH+(J2-1))*ISPLH+I2 ! lozenge-local subcell index
                    DO P=1,NSCT
                        DO IEL=1,NM
                            FLUX(IEL,P,IOF,IHEX_XY,K)=FLUX(IEL,P,IOF,IHEX_XY,K)+Q2(IEL,NM+1)*DN(P)
                        ENDDO
                    ENDDO

                ENDDO ! END OF LOZANGE SUBMESH J-LOOP
            ENDDO ! END OF LOZANGE SUBMESH I-LOOP
        ENDDO ! END OF LOZENGES LOOP
    ENDDO ! END OF Z-AXIS PLANES LOOP
    RETURN
  END SUBROUTINE SNSWH3
END MODULE SNSWH_MOD
