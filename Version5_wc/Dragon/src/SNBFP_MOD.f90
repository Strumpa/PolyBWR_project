!
!---------------------------------------------------------------------
!
!Purpose:
! Compute flux solution in a cell in 1D/2D/3D Cartesian geometry for the 
! discrete ordinates (SN) Boltzmann Fokker-Planck equation.
!
! Copyright:
!  Copyright (C) 2025 Ecole Polytechnique de Montreal
!  This library is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License as published by the Free Software Foundation; either
!  version 2.1 of the License, or (at your option) any later version
!
! Author(s): A. Hebert and C. Bienvenue
!
!---------------------------------------------------------------------
!
MODULE SNBFP_MOD

  USE SNADPT_MOD

CONTAINS

  SUBROUTINE SNBFP1(DX,SIGMA,BM,BP,V,QN,Q2,XNI,FEP,NM,NMX,NME,ISCHM,ESCHM,IELEM,EELEM,CST,WX,WE,SX,LFIXUP,IBFP,DELTAE)
    !
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    !   Compute flux solution in a cell in 1D Cartesian geometry for the 
    !   discrete ordinates (SN) Boltzmann Fokker-Planck equation.
    !
    ! Parameters: input
    !   DX      factor containing first direction cosines (mu).
    !   SIGMA   macroscopic total cross section in the cell.
    !   BM      stopping power at the upper group boundary.
    !   BP      stopping power at the lower group boundary.
    !   V       cell volume.
    !   QN      Legendre components of the angular source.
    !   NM      total number of moments of the flux in space and in energy.
    !   NMX     number of incoming/outgoing boundary flux moments in space.
    !   NME     number of incoming/outgoing boundary flux moments in energy.
    !   ISCHM   spatial discretization scheme index.
    !   ESCHM   energy discretization scheme index.
    !   IELEM   (order+1) of the spatial approximation polynomial.
    !   EELEM   (order+1) of the energy approximation polynomial.
    !   CST     Legendre coefficients for the polynomial approximations.
    !   WX      spatial closure relation weighting factors.
    !   WE      energy closure relation weighting factors.
    !   SX      sweeping direction.
    !   LFIXUP  flag to enable negative flux fixup.
    !   IBFP    type of energy proparation relation.
    !   DELTAE  energy group width.
    !
    ! Parameters: input/output
    !   Q2      matrix of the linear system to solve for the flux moments.
    !   XNI     incoming/outgoing boundary flux in space.
    !   FEP     incoming/outgoing boundary flux in energy.
    !
    !-----------------------------------------------------------------------

    ! VARIABLE DECLARATION
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: LFIXUP
    INTEGER, INTENT(IN) :: NM,NMX,NME,ISCHM,ESCHM,IELEM,EELEM,IBFP
    REAL, INTENT(IN) :: SIGMA,BM,BP,V,CST(IELEM),WX(IELEM+1),WE(EELEM+1),DX,SX,DELTAE
    DOUBLE PRECISION, INTENT(IN) :: QN(NM)
    DOUBLE PRECISION, INTENT(INOUT) :: XNI(NMX),Q2(NM,NM+1),FEP(NME)

    ! LOCAL VARIABLES
    INTEGER :: IER,II,JJ,IIX,IIE,IE,IX,JX,JE,KE
    REAL, PARAMETER :: RLOG=1.0E-8
    REAL :: TB
    DOUBLE PRECISION :: FACT
    REAL :: WX0(2),WE0(2)

    !----
    ! TYPE OF ENERGY PROPAGATION FACTOR
    !----
    IF(IBFP.EQ.1) THEN
        ! GALERKIN TYPE
        TB=BM/BP
    ELSE
        ! PRZYBYLSKI AND LIGOU TYPE
        TB=1.0
    ENDIF

    !----
    ! MATRIX OF LEGENDRE MOMENT COEFFICIENTS OF THE FLUX
    !----
    
    ! CONSTANT EXPANSION OF THE FLUX
    IF(IELEM.EQ.1.AND.EELEM.EQ.1) THEN

        IF(ISCHM.EQ.3.OR.ESCHM.EQ.3) THEN
            IF(IELEM.NE.1.OR.EELEM.NE.1) CALL XABORT('SNBFP1: AWD SCHEME NOT IMPLEMENTED FOR IELEM>1.')
            CALL SNADPT2D(WX0,WE0,1.0,TB,XNI(1),DBLE(FEP(1)),QN(1),DX,V*BP,SIGMA,V)
            IF(ISCHM.NE.3) WX0=WX
            IF(ESCHM.NE.3) WE0=WE
            Q2(1,1)=SIGMA*V+((WE0(2)-1.0)*TB+1.0)*BP*V+WX0(2)*DX*SX
            Q2(1,2)=QN(1)*V-(BP*WE0(1)*TB-BM)*V*FEP(1)-DX*SX*(WX0(1)-1.0)*XNI(1)
        ELSE
            Q2(1,1)=SIGMA*V+((WE(2)-1.0)*TB+1.0)*BP*V+WX(2)*DX*SX
            Q2(1,2)=QN(1)*V-(BP*WE(1)*TB-BM)*V*FEP(1)-DX*SX*(WX(1)-1.0)*XNI(1)
        ENDIF

    ! GENERAL POLYNOMIAL EXPANSION OF THE FLUX
    ELSE

        IF(ISCHM.EQ.3.AND.ESCHM.EQ.3) CALL XABORT('SNBFP1: AWD SCHEME NOT IMPLEMENTED FOR IELEM>1 OR EELEM>1.')

        ! COLLISION TERM
        DO IX=1,IELEM
        DO IE=1,EELEM
            II=EELEM*(IX-1)+IE
            Q2(II,II)=Q2(II,II)+SIGMA*V
        ENDDO
        ENDDO

        ! CSD TERM
        DO IE=1,EELEM
        DO JE=1,EELEM
            FACT=CST(IE)*CST(JE)*((WE(JE+1)-1.0)*TB+1.0)*BP*(-1.0)**(IE+JE)*V
            DO KE=1,IE-1
                IF(MOD(IE-KE,2).EQ.1) THEN ! LINEAR STOPPING POWER
                    IF(JE.EQ.KE) FACT=FACT+CST(IE)*CST(JE)*CST(KE)*(BP+BM)*V*1/(2*JE-1)
                    IF(JE+1.EQ.KE) FACT=FACT+CST(IE)*CST(JE)*CST(KE)*(BM-BP)*JE/((2*JE-1)*(2*JE+1))*V
                    IF(JE-1.EQ.KE) FACT=FACT+CST(IE)*CST(JE)*CST(KE)*(BM-BP)*(JE-1)/((2*JE-3)*(2*JE-1))*V
                ENDIF
            ENDDO
            DO IX=1,IELEM
                II=EELEM*(IX-1)+IE
                JJ=EELEM*(IX-1)+JE
                Q2(II,JJ)=Q2(II,JJ)+FACT
            ENDDO
        ENDDO
        ENDDO

        ! X-STREAMING TERM
        DO IX=1,IELEM
        DO JX=1,IELEM
            FACT=CST(IX)*CST(JX)*WX(JX+1)*DX*SX**(IX+JX-1)
            IF(IX.GE.JX+1.AND.MOD(IX-JX,2).EQ.1) THEN
                FACT=FACT-2.0*CST(IX)*CST(JX)*DX
            ENDIF
            DO IE=1,EELEM
                II=EELEM*(IX-1)+IE
                JJ=EELEM*(JX-1)+IE
                Q2(II,JJ)=Q2(II,JJ)+FACT
            ENDDO
        ENDDO
        ENDDO

        !----
        ! SOURCE TERM
        !---- 
        DO IX=1,IELEM
        DO IE=1,EELEM
            II=EELEM*(IX-1)+IE
            Q2(II,NM+1)=Q2(II,NM+1)+QN(II)*V
        ENDDO
        ENDDO
        DO IE=1,EELEM
            FACT=CST(IE)*((-1.0)**(IE-1)*BP*WE(1)*TB-BM)*V
            DO IX=1,IELEM
                II=EELEM*(IX-1)+IE
                IIE=IX
                Q2(II,NM+1)=Q2(II,NM+1)-FACT*FEP(IIE)
            ENDDO
        ENDDO
        DO IX=1,IELEM
            FACT=CST(IX)*DX*SX**IX*(WX(1)-(-1.0)**(IX-1))
            DO IE=1,EELEM
                II=EELEM*(IX-1)+IE
                IIX=IE
                Q2(II,NM+1)=Q2(II,NM+1)-FACT*XNI(IIX)
            ENDDO
        ENDDO
    ENDIF

    !----
    ! SOLVE THE LINEAR SYSTEM
    !----
    CALL ALSBD(NM,1,Q2,IER,NM)
    IF(IER.NE.0) CALL XABORT('SNBFP1: SINGULAR MATRIX.')

    !----
    ! CLOSURE RELATIONS
    !----

    ! CONSTANT EXPANSION OF THE FLUX
    IF(IELEM.EQ.1.AND.EELEM.EQ.1) THEN

        IF(LFIXUP.AND.(Q2(1,NM+1).LE.RLOG)) Q2(1,NM+1)=0.0 ! FLUX FIX-UP
        IF(ISCHM.EQ.3.OR.ESCHM.EQ.3) THEN
            IF(ISCHM.NE.3) WX0=WX
            IF(ESCHM.NE.3) WE0=WE
            FEP(1)=WE0(1)*TB*FEP(1)+((WE0(2)-1.0)*TB+1.0)*REAL(Q2(1,2))
            XNI(1)=WX0(1)*XNI(1)+WX0(2)*Q2(1,2)
        ELSE
            FEP(1)=WE(1)*TB*FEP(1)+((WE(2)-1.0)*TB+1.0)*REAL(Q2(1,2))
            XNI(1)=WX(1)*XNI(1)+WX(2)*Q2(1,2)
        ENDIF
        IF(LFIXUP) THEN ! FLUX FIX-UP
            IF(FEP(1).LE.RLOG) FEP(1)=0.0
            IF(XNI(1).LE.RLOG) XNI(1)=0.0
        ENDIF

    ! GENERAL POLYNOMIAL EXPANSION OF THE FLUX
    ELSE
        FEP(:NME)=WE(1)*TB*FEP(:NME)
        XNI(:NMX)=WX(1)*XNI(:NMX)
        DO IE=1,EELEM
            FACT=CST(IE)*(-1.0)**(IE-1)*((WE(IE+1)-1.0)*TB+1.0)
            DO IX=1,IELEM
                II=EELEM*(IX-1)+IE
                IIE=IX
                FEP(IIE)=FEP(IIE)+REAL(FACT*Q2(II,NM+1))
            ENDDO
        ENDDO
        DO IX=1,IELEM
            FACT=CST(IX)*SX**(IX-1)*WX(IX+1)
            DO IE=1,EELEM
                II=EELEM*(IX-1)+IE
                IIX=IE
                XNI(IIX)=XNI(IIX)+FACT*Q2(II,NM+1)
            ENDDO
        ENDDO
    ENDIF
    FEP(:NME)=FEP(:NME)/DELTAE

    ! RETURN FLUX MOMENTS
    RETURN
  END SUBROUTINE SNBFP1
  !
  SUBROUTINE SNBFP2(DX,DY,SIGMA,BM,BP,V,QN,Q2,XNI,XNJ,FEP,NM,NMX,NMY,NME, &
  & ISCHM,ESCHM,IELEM,EELEM,CST,WX,WY,WE,SX,SY,LFIXUP,IBFP,DELTAE)
    !
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    !   Compute flux solution in a cell in 2D Cartesian geometry for the 
    !   discrete ordinates (SN) Boltzmann Fokker-Planck equation.
    !
    ! Parameters: input
    !   DX      factor containing first direction cosines (mu).
    !   DY      factor containing second direction cosines (eta).
    !   SIGMA   macroscopic total cross section in the cell.
    !   BM      stopping power at the upper group boundary.
    !   BP      stopping power at the lower group boundary.
    !   V       cell volume.
    !   QN      Legendre components of the angular source.
    !   NM      total number of moments of the flux in space and in energy.
    !   NMX     number of incoming/outgoing boundary flux moments in x.
    !   NMY     number of incoming/outgoing boundary flux moments in y.
    !   NME     number of incoming/outgoing boundary flux moments in energy.
    !   ISCHM   spatial discretization scheme index.
    !   ESCHM   energy discretization scheme index.
    !   IELEM   (order+1) of the spatial approximation polynomial.
    !   EELEM   (order+1) of the energy approximation polynomial.
    !   CST     Legendre coefficients for the polynomial approximations.
    !   WX      spatial closure relation weighting factors.
    !   WY      spatial closure relation weighting factors.
    !   WE      energy closure relation weighting factors.
    !   SX      sweeping direction.
    !   SY      sweeping direction.
    !   LFIXUP  flag to enable negative flux fixup.
    !   IBFP    type of energy proparation relation.
    !   DELTAE  energy group width.
    !
    ! Parameters: input/output
    !   Q2      matrix of the linear system to solve for the flux moments.
    !   XNI     incoming/outgoing boundary flux in space.
    !   XNJ     incoming/outgoing boundary flux in space.
    !   FEP     incoming/outgoing boundary flux in energy.
    !
    !-----------------------------------------------------------------------

    ! VARIABLE DECLARATION
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: LFIXUP
    INTEGER, INTENT(IN) :: NM,NMX,NMY,NME,ISCHM,ESCHM,IELEM,EELEM,IBFP
    REAL, INTENT(IN) :: SIGMA,BM,BP,V,CST(IELEM),WX(IELEM+1),WY(IELEM+1), &
    & WE(EELEM+1),DX,DY,SX,SY,DELTAE
    REAL, INTENT(INOUT) :: FEP(NME)
    DOUBLE PRECISION, INTENT(IN) :: QN(NM)
    DOUBLE PRECISION, INTENT(INOUT) :: XNI(NMX),XNJ(NMY),Q2(NM,NM+1)

    ! LOCAL VARIABLES
    INTEGER :: IER,II,JJ,IIX,IIY,IIE,IY,IE,IX,JX,JY,JE,KE
    REAL, PARAMETER :: RLOG=1.0E-8
    REAL :: TB
    DOUBLE PRECISION :: FACT
    REAL :: WX0(2),WY0(2),WE0(2)

    !----
    ! TYPE OF ENERGY PROPAGATION FACTOR
    !----
    IF(IBFP.EQ.1) THEN
        ! GALERKIN TYPE
        TB=BM/BP
    ELSE
        ! PRZYBYLSKI AND LIGOU TYPE
        TB=1.0
    ENDIF

    !----
    ! MATRIX OF LEGENDRE MOMENT COEFFICIENTS OF THE FLUX
    !----
    
    ! CONSTANT EXPANSION OF THE FLUX
    IF(IELEM.EQ.1.AND.EELEM.EQ.1) THEN

        IF(ISCHM.EQ.3.OR.ESCHM.EQ.3) THEN
            IF(IELEM.NE.1.OR.EELEM.NE.1) CALL XABORT('SNBFP2: AWD SCHEME NOT IMPLEMENTED FOR IELEM>1.')
            CALL SNADPT3D(WX0,WY0,WE0,1.0,1.0,TB,XNI(1),XNJ(1),DBLE(FEP(1)),QN(1),DX,DY,V*BP,SIGMA,V)
            IF(ISCHM.NE.3) THEN
                WX0=WX
                WY0=WY
            ENDIF
            IF(ESCHM.NE.3) WE0=WE
            Q2(1,1)=SIGMA*V+((WE0(2)-1.0)*TB+1.0)*BP*V+WX0(2)*DX*SX+WY0(2)*DY*SY
            Q2(1,2)=Q2(1,2)+QN(1)*V-(BP*WE0(1)*TB-BM)*V*FEP(1)-DX*SX*(WX0(1)-1.0)*XNI(1)-DY*SY*(WY0(1)-1.0)*XNJ(1)
        ELSE
            Q2(1,1)=SIGMA*V+((WE(2)-1.0)*TB+1.0)*BP*V+WX(2)*DX*SX+WY(2)*DY*SY
            Q2(1,2)=Q2(1,2)+QN(1)*V-(BP*WE(1)*TB-BM)*V*FEP(1)-DX*SX*(WX(1)-1.0)*XNI(1)-DY*SY*(WY(1)-1.0)*XNJ(1)
        ENDIF

    ! GENERAL POLYNOMIAL EXPANSION OF THE FLUX
    ELSE

        IF(ISCHM.EQ.3.AND.ESCHM.EQ.3) CALL XABORT('SNBFP2: AWD SCHEME NOT IMPLEMENTED FOR IELEM>1 OR EELEM>1.')

        ! COLLISION TERM
        DO IY=1,IELEM
        DO IX=1,IELEM
        DO IE=1,EELEM
            II=EELEM*IELEM*(IY-1)+EELEM*(IX-1)+IE
            Q2(II,II)=Q2(II,II)+SIGMA*V
        ENDDO
        ENDDO
        ENDDO

        ! CSD TERM
        DO IE=1,EELEM
        DO JE=1,EELEM
            FACT=CST(IE)*CST(JE)*((WE(JE+1)-1.0)*TB+1.0)*BP*(-1.0)**(IE+JE)*V
            DO KE=1,IE-1
                IF(MOD(IE-KE,2).EQ.1) THEN ! LINEAR STOPPING POWER
                    IF(JE.EQ.KE) FACT=FACT+CST(IE)*CST(JE)*CST(KE)*(BP+BM)*V*1/(2*JE-1)
                    IF(JE+1.EQ.KE) FACT=FACT+CST(IE)*CST(JE)*CST(KE)*(BM-BP)*JE/((2*JE-1)*(2*JE+1))*V
                    IF(JE-1.EQ.KE) FACT=FACT+CST(IE)*CST(JE)*CST(KE)*(BM-BP)*(JE-1)/((2*JE-3)*(2*JE-1))*V
                ENDIF
            ENDDO
            DO IY=1,IELEM
            DO IX=1,IELEM
                II=EELEM*IELEM*(IY-1)+EELEM*(IX-1)+IE
                JJ=EELEM*IELEM*(IY-1)+EELEM*(IX-1)+JE
                Q2(II,JJ)=Q2(II,JJ)+FACT
            ENDDO
            ENDDO
        ENDDO
        ENDDO

        ! X-STREAMING TERM
        DO IX=1,IELEM
        DO JX=1,IELEM
            FACT=CST(IX)*CST(JX)*WX(JX+1)*DX*SX**(IX+JX-1)
            IF(IX.GE.JX+1.AND.MOD(IX-JX,2).EQ.1) THEN
                FACT=FACT-2.0*CST(IX)*CST(JX)*DX
            ENDIF
            DO IY=1,IELEM
            DO IE=1,EELEM
                II=EELEM*IELEM*(IY-1)+EELEM*(IX-1)+IE
                JJ=EELEM*IELEM*(IY-1)+EELEM*(JX-1)+IE
                Q2(II,JJ)=Q2(II,JJ)+FACT
            ENDDO
            ENDDO
        ENDDO
        ENDDO

        ! Y-STREAMING TERM
        DO IY=1,IELEM
        DO JY=1,IELEM
            FACT=CST(IY)*CST(JY)*WY(JY+1)*DY*SY**(IY+JY-1)
            IF(IY.GE.JY+1.AND.MOD(IY-JY,2).EQ.1) THEN
                FACT=FACT-2.0*CST(IY)*CST(JY)*DY
            ENDIF
            DO IX=1,IELEM
            DO IE=1,EELEM
                II=EELEM*IELEM*(IY-1)+EELEM*(IX-1)+IE
                JJ=EELEM*IELEM*(JY-1)+EELEM*(IX-1)+IE
                Q2(II,JJ)=Q2(II,JJ)+FACT
            ENDDO
            ENDDO
        ENDDO
        ENDDO

        !----
        ! SOURCE TERM
        !---- 
        DO IY=1,IELEM
        DO IX=1,IELEM
        DO IE=1,EELEM
            II=EELEM*IELEM*(IY-1)+EELEM*(IX-1)+IE
            Q2(II,NM+1)=Q2(II,NM+1)+QN(II)*V
        ENDDO
        ENDDO
        ENDDO
        DO IE=1,EELEM
            FACT=CST(IE)*((-1.0)**(IE-1)*BP*WE(1)*TB-BM)*V
            DO IX=1,IELEM
            DO IY=1,IELEM
                II=EELEM*IELEM*(IY-1)+EELEM*(IX-1)+IE
                IIE=IELEM*(IY-1)+IX
                Q2(II,NM+1)=Q2(II,NM+1)-FACT*FEP(IIE)
            ENDDO
            ENDDO
        ENDDO
        DO IX=1,IELEM
            FACT=CST(IX)*DX*SX**IX*(WX(1)-(-1.0)**(IX-1))
            DO IE=1,EELEM
            DO IY=1,IELEM
                II=EELEM*IELEM*(IY-1)+EELEM*(IX-1)+IE
                IIX=EELEM*(IY-1)+IE
                Q2(II,NM+1)=Q2(II,NM+1)-FACT*XNI(IIX)
            ENDDO
            ENDDO
        ENDDO
        DO IY=1,IELEM
            FACT=CST(IY)*DY*SY**IY*(WY(1)-(-1.0)**(IY-1))
            DO IE=1,EELEM
            DO IX=1,IELEM
                II=EELEM*IELEM*(IY-1)+EELEM*(IX-1)+IE
                IIY=EELEM*(IX-1)+IE
                Q2(II,NM+1)=Q2(II,NM+1)-FACT*XNJ(IIY)
            ENDDO
            ENDDO
        ENDDO
    ENDIF

    !----
    ! SOLVE THE LINEAR SYSTEM
    !----
    CALL ALSBD(NM,1,Q2,IER,NM)
    IF(IER.NE.0) CALL XABORT('SNBFP2: SINGULAR MATRIX.')

    !----
    ! CLOSURE RELATIONS
    !----

    ! CONSTANT EXPANSION OF THE FLUX
    IF(IELEM.EQ.1.AND.EELEM.EQ.1) THEN

        IF(LFIXUP.AND.(Q2(1,NM+1).LE.RLOG)) Q2(1,NM+1)=0.0 ! FLUX FIX-UP
        IF(ISCHM.EQ.3.OR.ESCHM.EQ.3) THEN
            IF(ISCHM.NE.3) THEN
                WX0=WX
                WY0=WY
            ENDIF
            IF(ESCHM.NE.3) WE0=WE
            FEP(1)=WE0(1)*TB*FEP(1)+((WE0(2)-1.0)*TB+1.0)*REAL(Q2(1,2))
            XNI(1)=WX0(1)*XNI(1)+WX0(2)*Q2(1,2)
            XNJ(1)=WY0(1)*XNJ(1)+WY0(2)*Q2(1,2)
        ELSE
            FEP(1)=WE(1)*TB*FEP(1)+((WE(2)-1.0)*TB+1.0)*REAL(Q2(1,2))
            XNI(1)=WX(1)*XNI(1)+WX(2)*Q2(1,2)
            XNJ(1)=WY(1)*XNJ(1)+WY(2)*Q2(1,2)
        ENDIF
        IF(LFIXUP) THEN ! FLUX FIX-UP
            IF(FEP(1).LE.RLOG) FEP(1)=0.0
            IF(XNI(1).LE.RLOG) XNI(1)=0.0
            IF(XNJ(1).LE.RLOG) XNJ(1)=0.0
        ENDIF

    ! GENERAL POLYNOMIAL EXPANSION OF THE FLUX
    ELSE
        FEP(:NME)=WE(1)*TB*FEP(:NME)
        XNI(:NMX)=WX(1)*XNI(:NMX)
        XNJ(:NMY)=WY(1)*XNJ(:NMY)
        DO IE=1,EELEM
            FACT=CST(IE)*(-1.0)**(IE-1)*((WE(IE+1)-1.0)*TB+1.0)
            DO IX=1,IELEM
            DO IY=1,IELEM
                II=EELEM*IELEM*(IY-1)+EELEM*(IX-1)+IE
                IIE=IELEM*(IY-1)+IX
                FEP(IIE)=FEP(IIE)+REAL(FACT*Q2(II,NM+1))
            ENDDO
            ENDDO
        ENDDO
        DO IX=1,IELEM
            FACT=CST(IX)*SX**(IX-1)*WX(IX+1)
            DO IE=1,EELEM
            DO IY=1,IELEM
                II=EELEM*IELEM*(IY-1)+EELEM*(IX-1)+IE
                IIX=EELEM*(IY-1)+IE
                XNI(IIX)=XNI(IIX)+FACT*Q2(II,NM+1)
            ENDDO
            ENDDO
        ENDDO
        DO IY=1,IELEM
            FACT=CST(IY)*SY**(IY-1)*WY(IY+1)
            DO IE=1,EELEM
            DO IX=1,IELEM
                II=EELEM*IELEM*(IY-1)+EELEM*(IX-1)+IE
                IIY=EELEM*(IX-1)+IE
                XNJ(IIY)=XNJ(IIY)+FACT*Q2(II,NM+1)
            ENDDO
            ENDDO
        ENDDO
    ENDIF

    FEP(:NME)=FEP(:NME)/DELTAE

    ! RETURN FLUX MOMENTS
    RETURN
  END SUBROUTINE SNBFP2
  !
  SUBROUTINE SNBFP3(DX,DY,DZ,SIGMA,BM,BP,V,QN,Q2,XNI,XNJ,XNK,FEP,NM,NMX,NMY,NMZ,NME,ISCHM, &
  & ESCHM,IELEM,EELEM,CST,WX,WY,WZ,WE,SX,SY,SZ,LFIXUP,IBFP,DELTAE)
    !
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    !   Compute flux solution in a cell in 3D Cartesian geometry for the 
    !   discrete ordinates (SN) Boltzmann Fokker-Planck equation.
    !
    ! Parameters: input
    !   DX      factor containing first direction cosines (mu).
    !   DY      factor containing second direction cosines (eta).
    !   DZ      factor containing third direction cosines (xi).
    !   SIGMA   macroscopic total cross section in the cell.
    !   BM      stopping power at the upper group boundary.
    !   BP      stopping power at the lower group boundary.
    !   V       cell volume.
    !   QN      Legendre components of the angular source.
    !   NM      total number of moments of the flux in space and in energy.
    !   NMX     number of incoming/outgoing boundary flux moments in x.
    !   NMY     number of incoming/outgoing boundary flux moments in y.
    !   NMZ     number of incoming/outgoing boundary flux moments in z.
    !   NME     number of incoming/outgoing boundary flux moments in energy.
    !   ISCHM   spatial discretization scheme index.
    !   ESCHM   energy discretization scheme index.
    !   IELEM   (order+1) of the spatial approximation polynomial.
    !   EELEM   (order+1) of the energy approximation polynomial.
    !   CST     Legendre coefficients for the polynomial approximations.
    !   WX      spatial closure relation weighting factors.
    !   WY      spatial closure relation weighting factors.
    !   WZ      spatial closure relation weighting factors.
    !   WE      energy closure relation weighting factors.
    !   SX      sweeping direction.
    !   SY      sweeping direction.
    !   SZ      sweeping direction.
    !   LFIXUP  flag to enable negative flux fixup.
    !   IBFP    type of energy proparation relation.
    !   DELTAE  energy group width.
    !
    ! Parameters: input/output
    !   Q2      matrix of the linear system to solve for the flux moments.
    !   XNI     incoming/outgoing boundary flux in space.
    !   XNJ     incoming/outgoing boundary flux in space.
    !   XNK     incoming/outgoing boundary flux in space.
    !   FEP     incoming/outgoing boundary flux in energy.
    !
    !-----------------------------------------------------------------------

    ! VARIABLE DECLARATION
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: LFIXUP
    INTEGER, INTENT(IN) :: NM,NMX,NMY,NMZ,NME,ISCHM,ESCHM,IELEM,EELEM,IBFP
    REAL, INTENT(IN) :: SIGMA,BM,BP,V,CST(IELEM),WX(IELEM+1),WY(IELEM+1),WZ(IELEM+1), &
    & WE(EELEM+1),DX,DY,DZ,SX,SY,SZ,DELTAE
    REAL, INTENT(INOUT) :: FEP(NME)
    DOUBLE PRECISION, INTENT(IN) :: QN(NM)
    DOUBLE PRECISION, INTENT(INOUT) :: XNI(NMX),XNJ(NMY),XNK(NMZ),Q2(NM,NM+1)

    ! LOCAL VARIABLES
    INTEGER :: IER,II,JJ,IIX,IIY,IIZ,IIE,IY,IZ,IE,IX,JX,JY,JZ,JE,KE
    REAL, PARAMETER :: RLOG=1.0E-8
    REAL :: TB
    DOUBLE PRECISION :: FACT
    REAL :: WX0(2),WY0(2),WZ0(2),WE0(2)

    !----
    ! TYPE OF ENERGY PROPAGATION FACTOR
    !----
    IF(IBFP.EQ.1) THEN
        ! GALERKIN TYPE
        TB=BM/BP
    ELSE
        ! PRZYBYLSKI AND LIGOU TYPE
        TB=1.0
    ENDIF

    !----
    ! MATRIX OF LEGENDRE MOMENT COEFFICIENTS OF THE FLUX
    !----
    
    ! CONSTANT EXPANSION OF THE FLUX
    IF(IELEM.EQ.1.AND.EELEM.EQ.1) THEN

        IF(ISCHM.EQ.3.OR.ESCHM.EQ.3) THEN
            IF(IELEM.NE.1.OR.EELEM.NE.1) CALL XABORT('SNBFP2: AWD SCHEME NOT IMPLEMENTED FOR IELEM>1.')
            CALL SNADPT4D(WX0,WY0,WZ0,WE0,1.0,1.0,1.0,TB,XNI(1),XNJ(1),XNJ(1),DBLE(FEP(1)),QN(1),DX,DY,DZ,V*BP,SIGMA,V)
            IF(ISCHM.NE.3) THEN
                WX0=WX
                WY0=WY
                WZ0=WZ
            ENDIF
            IF(ESCHM.NE.3) WE0=WE
            Q2(1,1)=SIGMA*V+((WE0(2)-1.0)*TB+1.0)*BP*V+WX0(2)*DX*SX+WY0(2)*DY*SY+WZ0(2)*DZ*SZ
            Q2(1,2)=Q2(1,2)+QN(1)*V-(BP*WE0(1)*TB-BM)*V*FEP(1)-DX*SX*(WX0(1)-1.0)*XNI(1)-DY*SY*(WY0(1)-1.0)*XNJ(1)- &
            & DZ*SZ*(WZ0(1)-1.0)*XNK(1)
        ELSE
            Q2(1,1)=SIGMA*V+((WE(2)-1.0)*TB+1.0)*BP*V+WX(2)*DX*SX+WY(2)*DY*SY+WZ(2)*DZ*SZ
            Q2(1,2)=Q2(1,2)+QN(1)*V-(BP*WE(1)*TB-BM)*V*FEP(1)-DX*SX*(WX(1)-1.0)*XNI(1)-DY*SY*(WY(1)-1.0)*XNJ(1)- &
            & DZ*SZ*(WZ(1)-1.0)*XNK(1)
        ENDIF

    ! GENERAL POLYNOMIAL EXPANSION OF THE FLUX
    ELSE

        ! COLLISION TERM
        DO IZ=1,IELEM
        DO IY=1,IELEM
        DO IX=1,IELEM
        DO IE=1,EELEM
            II=EELEM*IELEM**2*(IZ-1)+EELEM*IELEM*(IY-1)+EELEM*(IX-1)+IE
            Q2(II,II)=Q2(II,II)+SIGMA*V
        ENDDO
        ENDDO
        ENDDO
        ENDDO

        ! CSD TERM
        DO IE=1,EELEM
        DO JE=1,EELEM
            FACT=CST(IE)*CST(JE)*((WE(JE+1)-1.0)*TB+1.0)*BP*(-1.0)**(IE+JE)*V
            DO KE=1,IE-1
                IF(MOD(IE-KE,2).EQ.1) THEN ! LINEAR STOPPING POWER
                    IF(JE.EQ.KE) FACT=FACT+CST(IE)*CST(JE)*CST(KE)*(BP+BM)*V*1/(2*JE-1)
                    IF(JE+1.EQ.KE) FACT=FACT+CST(IE)*CST(JE)*CST(KE)*(BM-BP)*JE/((2*JE-1)*(2*JE+1))*V
                    IF(JE-1.EQ.KE) FACT=FACT+CST(IE)*CST(JE)*CST(KE)*(BM-BP)*(JE-1)/((2*JE-3)*(2*JE-1))*V
                ENDIF
            ENDDO
            DO IZ=1,IELEM
            DO IY=1,IELEM
            DO IX=1,IELEM
                II=EELEM*IELEM**2*(IZ-1)+EELEM*IELEM*(IY-1)+EELEM*(IX-1)+IE
                JJ=EELEM*IELEM**2*(IZ-1)+EELEM*IELEM*(IY-1)+EELEM*(IX-1)+JE
                Q2(II,JJ)=Q2(II,JJ)+FACT
            ENDDO
            ENDDO
            ENDDO
        ENDDO
        ENDDO

        ! X-STREAMING TERM
        DO IX=1,IELEM
        DO JX=1,IELEM
            FACT=CST(IX)*CST(JX)*WX(JX+1)*DX*SX**(IX+JX-1)
            IF(IX.GE.JX+1.AND.MOD(IX-JX,2).EQ.1) THEN
                FACT=FACT-2.0*CST(IX)*CST(JX)*DX
            ENDIF
            DO IZ=1,IELEM
            DO IY=1,IELEM
            DO IE=1,EELEM
                II=EELEM*IELEM**2*(IZ-1)+EELEM*IELEM*(IY-1)+EELEM*(IX-1)+IE
                JJ=EELEM*IELEM**2*(IZ-1)+EELEM*IELEM*(IY-1)+EELEM*(JX-1)+IE
                Q2(II,JJ)=Q2(II,JJ)+FACT
            ENDDO
            ENDDO
            ENDDO
        ENDDO
        ENDDO

        ! Y-STREAMING TERM
        DO IY=1,IELEM
        DO JY=1,IELEM
            FACT=CST(IY)*CST(JY)*WY(JY+1)*DY*SY**(IY+JY-1)
            IF(IY.GE.JY+1.AND.MOD(IY-JY,2).EQ.1) THEN
                FACT=FACT-2.0*CST(IY)*CST(JY)*DY
            ENDIF
            DO IZ=1,IELEM
            DO IX=1,IELEM
            DO IE=1,EELEM
                II=EELEM*IELEM**2*(IZ-1)+EELEM*IELEM*(IY-1)+EELEM*(IX-1)+IE
                JJ=EELEM*IELEM**2*(IZ-1)+EELEM*IELEM*(JY-1)+EELEM*(IX-1)+IE
                Q2(II,JJ)=Q2(II,JJ)+FACT
            ENDDO
            ENDDO
            ENDDO
        ENDDO
        ENDDO

        ! Z-STREAMING TERM
        DO IZ=1,IELEM
        DO JZ=1,IELEM
            FACT=CST(IZ)*CST(JZ)*WZ(JZ+1)*DZ*SZ**(IZ+JZ-1)
            IF(IZ.GE.JZ+1.AND.MOD(IZ-JZ,2).EQ.1) THEN
                FACT=FACT-2.0*CST(IZ)*CST(JZ)*DZ
            ENDIF
            DO IY=1,IELEM
            DO IX=1,IELEM
            DO IE=1,EELEM
                II=EELEM*IELEM**2*(IZ-1)+EELEM*IELEM*(IY-1)+EELEM*(IX-1)+IE
                JJ=EELEM*IELEM**2*(JZ-1)+EELEM*IELEM*(IY-1)+EELEM*(IX-1)+IE
                Q2(II,JJ)=Q2(II,JJ)+FACT
            ENDDO
            ENDDO
            ENDDO
        ENDDO
        ENDDO

        !----
        ! SOURCE TERM
        !---- 
        DO IZ=1,IELEM
        DO IY=1,IELEM
        DO IX=1,IELEM
        DO IE=1,EELEM
            II=EELEM*IELEM**2*(IZ-1)+EELEM*IELEM*(IY-1)+EELEM*(IX-1)+IE
            Q2(II,NM+1)=Q2(II,NM+1)+QN(II)*V
        ENDDO
        ENDDO
        ENDDO
        ENDDO
        DO IE=1,EELEM
            FACT=CST(IE)*((-1.0)**(IE-1)*BP*WE(1)*TB-BM)*V
            DO IX=1,IELEM
            DO IY=1,IELEM
            DO IZ=1,IELEM
                II=EELEM*IELEM**2*(IZ-1)+EELEM*IELEM*(IY-1)+EELEM*(IX-1)+IE
                IIE=IELEM**2*(IZ-1)+IELEM*(IY-1)+IX
                Q2(II,NM+1)=Q2(II,NM+1)-FACT*FEP(IIE)
            ENDDO
            ENDDO
            ENDDO
        ENDDO
        DO IX=1,IELEM
            FACT=CST(IX)*DX*SX**IX*(WX(1)-(-1.0)**(IX-1))
            DO IE=1,EELEM
            DO IY=1,IELEM
            DO IZ=1,IELEM
                II=EELEM*IELEM**2*(IZ-1)+EELEM*IELEM*(IY-1)+EELEM*(IX-1)+IE
                IIX=EELEM*IELEM*(IZ-1)+EELEM*(IY-1)+IE
                Q2(II,NM+1)=Q2(II,NM+1)-FACT*XNI(IIX)
            ENDDO
            ENDDO
            ENDDO
        ENDDO
        DO IY=1,IELEM
            FACT=CST(IY)*DY*SY**IY*(WY(1)-(-1.0)**(IY-1))
            DO IE=1,EELEM
            DO IX=1,IELEM
            DO IZ=1,IELEM
                II=EELEM*IELEM**2*(IZ-1)+EELEM*IELEM*(IY-1)+EELEM*(IX-1)+IE
                IIY=EELEM*IELEM*(IZ-1)+EELEM*(IX-1)+IE
                Q2(II,NM+1)=Q2(II,NM+1)-FACT*XNJ(IIY)
            ENDDO
            ENDDO
            ENDDO
        ENDDO
        DO IZ=1,IELEM
            FACT=CST(IZ)*DZ*SZ**IZ*(WZ(1)-(-1.0)**(IZ-1))
            DO IE=1,EELEM
            DO IX=1,IELEM
            DO IY=1,IELEM
                II=EELEM*IELEM**2*(IZ-1)+EELEM*IELEM*(IY-1)+EELEM*(IX-1)+IE
                IIZ=EELEM*IELEM*(IY-1)+EELEM*(IX-1)+IE
                Q2(II,NM+1)=Q2(II,NM+1)-FACT*XNK(IIZ)
            ENDDO
            ENDDO
            ENDDO
        ENDDO
    ENDIF

    !----
    ! SOLVE THE LINEAR SYSTEM
    !----
    CALL ALSBD(NM,1,Q2,IER,NM)
    IF(IER.NE.0) CALL XABORT('SNBFP3: SINGULAR MATRIX.')

    !----
    ! CLOSURE RELATIONS
    !----

    ! CONSTANT EXPANSION OF THE FLUX
    IF(IELEM.EQ.1.AND.EELEM.EQ.1) THEN

        IF(LFIXUP.AND.(Q2(1,NM+1).LE.RLOG)) Q2(1,NM+1)=0.0 ! FLUX FIX-UP
        IF(ISCHM.EQ.3.OR.ESCHM.EQ.3) THEN
            IF(ISCHM.NE.3) THEN
                WX0=WX
                WY0=WY
                WZ0=WZ
            ENDIF
            IF(ESCHM.NE.3) WE0=WE
            FEP(1)=WE0(1)*TB*FEP(1)+((WE0(2)-1.0)*TB+1.0)*REAL(Q2(1,2))
            XNI(1)=WX0(1)*XNI(1)+WX0(2)*Q2(1,2)
            XNJ(1)=WY0(1)*XNJ(1)+WY0(2)*Q2(1,2)
            XNK(1)=WZ0(1)*XNK(1)+WZ0(2)*Q2(1,2)
        ELSE
            FEP(1)=WE(1)*TB*FEP(1)+((WE(2)-1.0)*TB+1.0)*REAL(Q2(1,2))
            XNI(1)=WX(1)*XNI(1)+WX(2)*Q2(1,2)
            XNJ(1)=WY(1)*XNJ(1)+WY(2)*Q2(1,2)
            XNK(1)=WZ(1)*XNK(1)+WZ(2)*Q2(1,2)
        ENDIF
        IF(LFIXUP) THEN ! FLUX FIX-UP
            IF(FEP(1).LE.RLOG) FEP(1)=0.0
            IF(XNI(1).LE.RLOG) XNI(1)=0.0
            IF(XNJ(1).LE.RLOG) XNJ(1)=0.0
            IF(XNK(1).LE.RLOG) XNK(1)=0.0
        ENDIF

    ! GENERAL POLYNOMIAL EXPANSION OF THE FLUX
    ELSE
        FEP(:NME)=WE(1)*TB*FEP(:NME)
        XNI(:NMX)=WX(1)*XNI(:NMX)
        XNJ(:NMY)=WY(1)*XNJ(:NMY)
        XNK(:NMZ)=WZ(1)*XNK(:NMZ)
        DO IE=1,EELEM
            FACT=CST(IE)*(-1.0)**(IE-1)*((WE(IE+1)-1.0)*TB+1.0)
            DO IX=1,IELEM
            DO IY=1,IELEM
            DO IZ=1,IELEM
                II=EELEM*IELEM**2*(IZ-1)+EELEM*IELEM*(IY-1)+EELEM*(IX-1)+IE
                IIE=IELEM**2*(IZ-1)+IELEM*(IY-1)+IX
                FEP(IIE)=FEP(IIE)+REAL(FACT*Q2(II,NM+1))
            ENDDO
            ENDDO
            ENDDO
        ENDDO
        DO IX=1,IELEM
            FACT=CST(IX)*SX**(IX-1)*WX(IX+1)
            DO IE=1,EELEM
            DO IY=1,IELEM
            DO IZ=1,IELEM
                II=EELEM*IELEM**2*(IZ-1)+EELEM*IELEM*(IY-1)+EELEM*(IX-1)+IE
                IIX=EELEM*IELEM*(IZ-1)+EELEM*(IY-1)+IE
                XNI(IIX)=XNI(IIX)+FACT*Q2(II,NM+1)
            ENDDO
            ENDDO
            ENDDO
        ENDDO
        DO IY=1,IELEM
            FACT=CST(IY)*SY**(IY-1)*WY(IY+1)
            DO IE=1,EELEM
            DO IX=1,IELEM
            DO IZ=1,IELEM
                II=EELEM*IELEM**2*(IZ-1)+EELEM*IELEM*(IY-1)+EELEM*(IX-1)+IE
                IIY=EELEM*IELEM*(IZ-1)+EELEM*(IX-1)+IE
                XNJ(IIY)=XNJ(IIY)+FACT*Q2(II,NM+1)
            ENDDO
            ENDDO
            ENDDO
        ENDDO
        DO IZ=1,IELEM
            FACT=CST(IZ)*SZ**(IZ-1)*WZ(IZ+1)
            DO IE=1,EELEM
            DO IX=1,IELEM
            DO IY=1,IELEM
                II=EELEM*IELEM**2*(IZ-1)+EELEM*IELEM*(IY-1)+EELEM*(IX-1)+IE
                IIZ=EELEM*IELEM*(IY-1)+EELEM*(IX-1)+IE
                XNK(IIZ)=XNK(IIZ)+FACT*Q2(II,NM+1)
            ENDDO
            ENDDO
            ENDDO
        ENDDO
    ENDIF

    FEP(:NME)=FEP(:NME)/DELTAE

    ! RETURN FLUX MOMENTS
    RETURN
  END SUBROUTINE SNBFP3
END MODULE SNBFP_MOD
