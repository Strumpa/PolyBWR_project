!-----------------------------------------------------------------------
!
! Purpose:
!   Compute flux solution in a cell in 1D/2D/3D Cartesian geometry for the 
!   discrete ordinates (SN) Boltzmann transport equation.
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
MODULE SNBTE_MOD

  USE SNADPT_MOD

CONTAINS

  SUBROUTINE SNBTE1(DX,SIGMA,V,QN,Q2,XNI,NM,ISCHM,IELEM,CST,WX,SX,LFIXUP)
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    !   Compute flux solution in a cell in 1D Cartesian geometry for the 
    !   discrete ordinates (SN) Boltzmann transport equation.
    !
    ! Parameters: input
    !   DX      factor containing first direction cosines (mu).
    !   SIGMA   macroscopic total cross section in the cell.
    !   V       cell volume.
    !   QN      Legendre components of the angular source.
    !   NM      total number of moments of the flux in space and in energy.
    !   ISCHM   spatial discretization scheme index.
    !   IELEM   (order+1) of the spatial approximation polynomial.
    !   CST     Legendre coefficients for the polynomial approximations.
    !   WX      spatial closure relation weighting factors.
    !   SX      sweeping direction.
    !   LFIXUP  flag to enable negative flux fixup.
    !
    ! Parameters: input/output
    !   Q2      matrix of the linear system to solve for the flux moments.
    !   XNI     incoming/outgoing boundary flux in space.
    !
    !-----------------------------------------------------------------------

    ! VARIABLE DECLARATION
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: LFIXUP
    INTEGER, INTENT(IN) :: NM,ISCHM,IELEM
    REAL, INTENT(IN) :: SIGMA,V,CST(IELEM),WX(IELEM+1),DX,SX
    DOUBLE PRECISION, INTENT(IN) :: QN(NM)
    DOUBLE PRECISION, INTENT(INOUT) :: XNI,Q2(NM,NM+1)

    ! LOCAL VARIABLES
    INTEGER :: IER,IX,JX
    REAL, PARAMETER :: RLOG=1.0E-8
    DOUBLE PRECISION :: FACT
    REAL :: WX0(2)

    !----
    ! MATRIX OF LEGENDRE MOMENT COEFFICIENTS OF THE FLUX AND SOURCE TERMS
    !----
    IF(ISCHM.EQ.3.AND.IELEM.NE.1) CALL XABORT('SNBTE2: AWD SCHEME NOT IMPLEMENTED FOR IELEM>1.')

    ! CONSTANT EXPANSION OF THE FLUX
    IF(IELEM.EQ.1.AND.ISCHM.LE.3) THEN

        ! DD0 SCHEME
        IF(ISCHM.EQ.1) THEN

            Q2(1,1)=V*SIGMA+2.0D0*ABS(DX)
            Q2(1,2)=V*QN(1)+2.0D0*ABS(DX)*XNI

        ! DG0 SCHEME
        ELSEIF(ISCHM.EQ.2) THEN

            Q2(1,1)=V*SIGMA+ABS(DX)
            Q2(1,2)=V*QN(1)+ABS(DX)*XNI

        ! AWD0 SCHEME
        ELSEIF(ISCHM.EQ.3) THEN

            CALL SNADPT1D(WX0,1.0,XNI,QN(1),DX,SIGMA,V)
            Q2(1,1)=V*SIGMA+WX0(2)*ABS(DX)
            Q2(1,2)=V*QN(1)+(1.0-WX0(1))*ABS(DX)*XNI

        ENDIF

    ! LINEAR EXPANSION OF THE FLUX
    ELSEIF(IELEM.EQ.2.AND.ISCHM.LE.2) THEN

        ! DD1 SCHEME
        IF(ISCHM.EQ.1) THEN

            Q2(1,1)=V*SIGMA
            Q2(1,2)=2.0D0*CST(2)*DX
            Q2(2,1)=-2.0D0*CST(2)*DX
            Q2(2,2)=V*SIGMA+6.0D0*ABS(DX)

            Q2(1,3)=V*QN(1)
            Q2(2,3)=V*QN(2)-2.0D0*CST(2)*DX*XNI

        ! DG1 SCHEME
        ELSEIF(ISCHM.EQ.2) THEN

            Q2(1,1)=V*SIGMA+ABS(DX)
            Q2(1,2)=CST(2)*DX
            Q2(2,1)=-CST(2)*DX
            Q2(2,2)=V*SIGMA+3*ABS(DX)

            Q2(1,3)=V*QN(1)+DX*SX*XNI
            Q2(2,3)=V*QN(2)-CST(2)*DX*XNI

        ENDIF

    ! GENERAL POLYNOMIAL EXPANSION OF THE FLUX
    ELSE

        ! COLLISION TERM
        DO IX=1,IELEM
            Q2(IX,IX)=Q2(IX,IX)+SIGMA*V
        ENDDO

        ! X-STREAMING TERM
        DO IX=1,IELEM
        DO JX=1,IELEM
            FACT=CST(IX)*CST(JX)*WX(JX+1)*DX*SX**(IX+JX-1)
            IF(IX.GE.JX+1.AND.MOD(IX-JX,2).EQ.1) THEN
                FACT=FACT-2.0*CST(IX)*CST(JX)*DX
            ENDIF
            Q2(IX,JX)=Q2(IX,JX)+FACT
        ENDDO
        ENDDO

        ! SOURCE TERM
        DO IX=1,IELEM
            Q2(IX,NM+1)=Q2(IX,NM+1)+QN(IX)*V-CST(IX)*DX*SX**IX*(WX(1)-(-1.0)**(IX-1))*XNI
        ENDDO

    ENDIF

    !----
    ! SOLVE THE LINEAR SYSTEM
    !----
    CALL ALSBD(NM,1,Q2,IER,NM)
    IF(IER.NE.0) CALL XABORT('SNBTE: SINGULAR MATRIX.')

    !----
    ! CLOSURE RELATIONS
    !----

    ! CONSTANT EXPANSION OF THE FLUX
    IF(IELEM.EQ.1.AND.ISCHM.LE.3) THEN

        IF(LFIXUP.AND.(Q2(1,NM+1).LE.RLOG)) Q2(1,2)=0.0 ! FLUX FIX-UP

        ! DD SCHEME
        IF(ISCHM.EQ.1) THEN

            XNI=2.0D0*Q2(1,2)-XNI

        ! DG SCHEME
        ELSEIF(ISCHM.EQ.2) THEN

            XNI=Q2(1,2) 

        ! AWD SCHEME
        ELSEIF(ISCHM.EQ.3) THEN

            XNI=WX0(1)*XNI+WX0(2)*Q2(1,2)

        ENDIF

        IF(LFIXUP.AND.XNI.LE.RLOG) XNI=0.0 ! FLUX FIX-UP

    ! LINEAR EXPANSION OF THE FLUX
    ELSEIF(IELEM.EQ.2.AND.ISCHM.LE.2) THEN

        ! DD SCHEME
        IF(ISCHM.EQ.1) THEN

            XNI=XNI+2.0D0*SX*CST(2)*Q2(2,3)

        ! DG SCHEME
        ELSEIF(ISCHM.EQ.2) THEN

            XNI=Q2(1,3)+SX*CST(2)*Q2(2,3)

        ENDIF

    ! GENERAL POLYNOMIAL EXPANSION OF THE FLUX
    ELSE
        XNI=WX(1)*XNI
        DO IX=1,IELEM
            XNI=XNI+CST(IX)*SX**(IX-1)*WX(IX+1)*Q2(IX,NM+1)
        ENDDO
    ENDIF

    RETURN
  END SUBROUTINE SNBTE1

  SUBROUTINE SNBTE2(DX,DY,SIGMA,V,QN,Q2,XNI,XNJ,NM,NMX,NMY,ISCHM,IELEM,CST,WX,WY,SX,SY,LFIXUP)
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    !   Compute flux solution in a cell in 2D Cartesian geometry for the 
    !   discrete ordinates (SN) Boltzmann transport equation.
    !
    ! Parameters: input
    !   DX      factor containing first direction cosines (mu).
    !   DY      factor containing second direction cosines (eta).
    !   SIGMA   macroscopic total cross section in the cell.
    !   V       cell volume.
    !   QN      Legendre components of the angular source.
    !   NM      total number of moments of the flux in space and in energy.
    !   ISCHM   spatial discretization scheme index.
    !   IELEM   (order+1) of the spatial approximation polynomial.
    !   CST     Legendre coefficients for the polynomial approximations.
    !   WX      spatial closure relation weighting factors.
    !   WY      spatial closure relation weighting factors.
    !   SX      sweeping direction.
    !   SY      sweeping direction.
    !   LFIXUP  flag to enable negative flux fixup.
    !
    ! Parameters: input/output
    !   Q2      matrix of the linear system to solve for the flux moments.
    !   XNI     incoming/outgoing boundary flux in space.
    !   XNJ     incoming/outgoing boundary flux in space.
    !
    !-----------------------------------------------------------------------

    ! VARIABLE DECLARATION
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: LFIXUP
    INTEGER, INTENT(IN) :: NM,NMX,NMY,ISCHM,IELEM
    REAL, INTENT(IN) :: SIGMA,V,CST(IELEM),WX(IELEM+1),WY(IELEM+1),DX,DY,SX,SY
    DOUBLE PRECISION, INTENT(IN) :: QN(NM)
    DOUBLE PRECISION, INTENT(INOUT) :: XNI(NMX),XNJ(NMY),Q2(NM,NM+1)

    ! LOCAL VARIABLES
    INTEGER :: IER,II,JJ,IIX,IIY,IX,IY,JX,JY
    REAL, PARAMETER :: RLOG=1.0E-8
    DOUBLE PRECISION :: FACT
    REAL :: WX0(2),WY0(2)

    !----
    ! MATRIX OF LEGENDRE MOMENT COEFFICIENTS OF THE FLUX AND SOURCE TERMS
    !----
    IF(ISCHM.EQ.3.AND.IELEM.NE.1) CALL XABORT('SNBTE2: AWD SCHEME NOT IMPLEMENTED FOR IELEM>1.')

    ! CONSTANT EXPANSION OF THE FLUX
    IF(IELEM.EQ.1.AND.ISCHM.LE.3) THEN

        ! DD0 SCHEME
        IF(ISCHM.EQ.1) THEN

            Q2(1,1)=V*SIGMA+2.0D0*ABS(DX)+2.0D0*ABS(DY)
            Q2(1,2)=V*QN(1)+2.0D0*ABS(DX)*XNI(1)+2.0D0*ABS(DY)*XNJ(1)

        ! DG0 SCHEME
        ELSEIF(ISCHM.EQ.2) THEN

            Q2(1,1)=V*SIGMA+ABS(DX)+ABS(DY) 
            Q2(1,2)=V*QN(1)+ABS(DX)*XNI(1)+ABS(DY)*XNJ(1)

        ! AWD0 SCHEME
        ELSEIF(ISCHM.EQ.3) THEN

            CALL SNADPT2D(WX0,WY0,1.0,1.0,XNI(1),XNJ(1),QN(1),DX,DY,SIGMA,V)
            Q2(1,1)=V*SIGMA+WX0(2)*ABS(DX)+WY0(2)*ABS(DY)
            Q2(1,2)=V*QN(1)+(1.0-WX0(1))*ABS(DX)*XNI(1)+(1.0-WY0(1))*ABS(DY)*XNJ(1)

        ENDIF

    ! LINEAR EXPANSION OF THE FLUX
    ELSEIF(IELEM.EQ.2.AND.ISCHM.LE.2) THEN

        ! DD1 SCHEME
        IF(ISCHM.EQ.1) THEN

            Q2(1,1)=V*SIGMA
            Q2(1,2)=2.0D0*CST(2)*DX
            Q2(1,3)=2.0D0*CST(2)*DY
            Q2(2,1)=-2.0D0*CST(2)*DX
            Q2(2,2)=V*SIGMA+6.0D0*ABS(DX)
            Q2(2,4)=2.0D0*CST(2)*DY
            Q2(3,1)=-2.0D0*CST(2)*DY
            Q2(3,3)=V*SIGMA+6.0D0*ABS(DY)
            Q2(3,4)=2.0D0*CST(2)*DX
            Q2(4,2)=-2.0D0*CST(2)*DY
            Q2(4,3)=-2.0D0*CST(2)*DX
            Q2(4,4)=V*SIGMA+6.0D0*ABS(DX)+6.0D0*ABS(DY)

            Q2(1,5)=V*QN(1)
            Q2(2,5)=V*QN(2)-2.0D0*CST(2)*DX*XNI(1)
            Q2(3,5)=V*QN(3)-2.0D0*CST(2)*DY*XNJ(1)
            Q2(4,5)=V*QN(4)-2.0D0*CST(2)*DX*XNI(2)-2.0D0*CST(2)*DY*XNJ(2)
            

        ! DG1 SCHEME
        ELSEIF(ISCHM.EQ.2) THEN

            Q2(1,1)=V*SIGMA+ABS(DX)+ABS(DY)
            Q2(1,2)=CST(2)*DX 
            Q2(1,3)=CST(2)*DY 
            Q2(2,1)=-CST(2)*DX 
            Q2(2,2)=V*SIGMA+3*ABS(DX)+ABS(DY) 
            Q2(2,4)=CST(2)*DY 
            Q2(3,1)=-CST(2)*DY 
            Q2(3,3)=V*SIGMA+ABS(DX)+3*ABS(DY)
            Q2(3,4)=CST(2)*DX 
            Q2(4,2)=-CST(2)*DY 
            Q2(4,3)=-CST(2)*DX 
            Q2(4,4)=V*SIGMA+3*ABS(DX)+3*ABS(DY) 

            Q2(1,5)=V*QN(1)+DX*SX*XNI(1)+DY*SY*XNJ(1)
            Q2(2,5)=V*QN(2)-CST(2)*DX*XNI(1)+DY*SY*XNJ(2) 
            Q2(3,5)=V*QN(3)+DX*SX*XNI(2)-CST(2)*DY*XNJ(1) 
            Q2(4,5)=V*QN(4)-CST(2)*DY*XNJ(2)-CST(2)*DX*XNI(2) 

        ENDIF

    ! GENERAL POLYNOMIAL EXPANSION OF THE FLUX
    ELSE

        ! COLLISION TERM
        DO IY=1,IELEM
        DO IX=1,IELEM
            II=IELEM*(IY-1)+IX
            Q2(II,II)=Q2(II,II)+SIGMA*V
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
                II=IELEM*(IY-1)+IX
                JJ=IELEM*(IY-1)+JX
                Q2(II,JJ)=Q2(II,JJ)+FACT
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
                II=IELEM*(IY-1)+IX
                JJ=IELEM*(JY-1)+IX
                Q2(II,JJ)=Q2(II,JJ)+FACT
            ENDDO
        ENDDO
        ENDDO

        ! SOURCE TERMS
        DO IY=1,IELEM
        DO IX=1,IELEM
            II=IELEM*(IY-1)+IX
            Q2(II,NM+1)=Q2(II,NM+1)+QN(II)*V
        ENDDO
        ENDDO
        DO IX=1,IELEM
            FACT=CST(IX)*DX*SX**IX*(WX(1)-(-1.0)**(IX-1))
            DO IY=1,IELEM
                II=IELEM*(IY-1)+IX
                IIX=IY
                Q2(II,NM+1)=Q2(II,NM+1)-FACT*XNI(IIX)
            ENDDO
        ENDDO
        DO IY=1,IELEM
            FACT=CST(IY)*DY*SY**IY*(WY(1)-(-1.0)**(IY-1))
            DO IX=1,IELEM
                II=IELEM*(IY-1)+IX
                IIY=IX
                Q2(II,NM+1)=Q2(II,NM+1)-FACT*XNJ(IIY)
            ENDDO
        ENDDO

    ENDIF

    !----
    ! SOLVE THE LINEAR SYSTEM
    !----
    CALL ALSBD(NM,1,Q2,IER,NM)
    IF(IER.NE.0) CALL XABORT('SNBTE: SINGULAR MATRIX.')

    !----
    ! CLOSURE RELATIONS
    !----

    ! CONSTANT EXPANSION OF THE FLUX
    IF(IELEM.EQ.1.AND.ISCHM.LE.3) THEN

        IF(LFIXUP.AND.(Q2(1,NM+1).LE.RLOG)) Q2(1,2)=0.0 ! FLUX FIX-UP

        ! DD SCHEME
        IF(ISCHM.EQ.1) THEN

            XNI(1)=2.0D0*Q2(1,2)-XNI(1)
            XNJ(1)=2.0D0*Q2(1,2)-XNJ(1)

        ! DG SCHEME
        ELSEIF(ISCHM.EQ.2) THEN

            XNI(1)=Q2(1,2) 
            XNJ(1)=Q2(1,2)

        ! AWD SCHEME
        ELSEIF(ISCHM.EQ.3) THEN

            XNI(1)=WX0(1)*XNI(1)+WX0(2)*Q2(1,2)
            XNJ(1)=WY0(1)*XNJ(1)+WY0(2)*Q2(1,2)

        ENDIF

        IF(LFIXUP) THEN ! FLUX FIX-UP
            IF(XNI(1).LE.RLOG) XNI(1)=0.0
            IF(XNJ(1).LE.RLOG) XNJ(1)=0.0
        ENDIF

    ! LINEAR EXPANSION OF THE FLUX
    ELSEIF(IELEM.EQ.2.AND.ISCHM.LE.2) THEN

        ! DD SCHEME
        IF(ISCHM.EQ.1) THEN

            XNI(1)=XNI(1)+2.0D0*SX*CST(2)*Q2(2,5)
            XNI(2)=XNI(2)+2.0D0*SX*CST(2)*Q2(4,5)
            XNJ(1)=XNJ(1)+2.0D0*SY*CST(2)*Q2(3,5)
            XNJ(2)=XNJ(2)+2.0D0*SY*CST(2)*Q2(4,5)

        ! DG SCHEME
        ELSEIF(ISCHM.EQ.2) THEN

            XNI(1)=Q2(1,5)+SX*CST(2)*Q2(2,5)
            XNI(2)=Q2(3,5)+SX*CST(2)*Q2(4,5)
            XNJ(1)=Q2(1,5)+SY*CST(2)*Q2(3,5)
            XNJ(2)=Q2(2,5)+SY*CST(2)*Q2(4,5)

        ENDIF

    ! GENERAL POLYNOMIAL EXPANSION OF THE FLUX
    ELSE
        XNI(:NMX)=WX(1)*XNI(:NMX)
        XNJ(:NMY)=WY(1)*XNJ(:NMY)
        DO IX=1,IELEM
            FACT=CST(IX)*SX**(IX-1)*WX(IX+1)
            DO IY=1,IELEM
                II=IELEM*(IY-1)+IX
                IIX=IY
                XNI(IIX)=XNI(IIX)+FACT*Q2(II,NM+1)
            ENDDO
        ENDDO
        DO IY=1,IELEM
            FACT=CST(IY)*SY**(IY-1)*WY(IY+1)
            DO IX=1,IELEM
                II=IELEM*(IY-1)+IX
                IIY=IX
                XNJ(IIY)=XNJ(IIY)+FACT*Q2(II,NM+1)
            ENDDO
        ENDDO
    ENDIF
    RETURN
  END SUBROUTINE SNBTE2

  SUBROUTINE SNBTE3(DX,DY,DZ,SIGMA,V,QN,Q2,XNI,XNJ,XNK,NM,NMX,NMY,NMZ,ISCHM,IELEM,CST,WX,WY,WZ,SX,SY,SZ,LFIXUP)
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    !   Compute flux solution in a cell in 3D Cartesian geometry for the 
    !   discrete ordinates (SN) Boltzmann transport equation.
    !
    ! Parameters: input
    !   DX      factor containing first direction cosines (mu).
    !   DY      factor containing second direction cosines (eta).
    !   DZ      factor containing third direction cosines (xi).
    !   SIGMA   macroscopic total cross section in the cell.
    !   V       cell volume.
    !   QN      Legendre components of the angular source.
    !   NM      total number of moments of the flux in space and in energy.
    !   ISCHM   spatial discretization scheme index.
    !   IELEM   (order+1) of the spatial approximation polynomial.
    !   CST     Legendre coefficients for the polynomial approximations.
    !   WX      spatial closure relation weighting factors.
    !   WY      spatial closure relation weighting factors.
    !   WZ      spatial closure relation weighting factors.
    !   SX      sweeping direction.
    !   SY      sweeping direction.
    !   SZ      sweeping direction.
    !   LFIXUP  flag to enable negative flux fixup.
    !
    ! Parameters: input/output
    !   Q2      matrix of the linear system to solve for the flux moments.
    !   XNI     incoming/outgoing boundary flux in space.
    !   XNJ     incoming/outgoing boundary flux in space.
    !   XNK     incoming/outgoing boundary flux in space.
    !
    !-----------------------------------------------------------------------

    ! VARIABLE DECLARATION
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: LFIXUP
    INTEGER, INTENT(IN) :: NM,NMX,NMY,NMZ,ISCHM,IELEM
    REAL, INTENT(IN) :: SIGMA,V,CST(IELEM),WX(IELEM+1),WY(IELEM+1),WZ(IELEM+1),DX,DY,DZ,SX,SY,SZ
    DOUBLE PRECISION, INTENT(IN) :: QN(NM)
    DOUBLE PRECISION, INTENT(INOUT) :: XNI(NMX),XNJ(NMY),XNK(NMZ),Q2(NM,NM+1)

    ! LOCAL VARIABLES
    INTEGER :: IER,II,JJ,IIX,IIY,IIZ,IX,IY,IZ,JX,JY,JZ
    REAL, PARAMETER :: RLOG=1.0E-8
    DOUBLE PRECISION :: FACT
    REAL :: WX0(2),WY0(2),WZ0(2)

    !----
    ! MATRIX OF LEGENDRE MOMENT COEFFICIENTS OF THE FLUX AND SOURCE TERMS
    !----
    IF(ISCHM.EQ.3.AND.IELEM.NE.1) CALL XABORT('SNBTE3: AWD SCHEME NOT IMPLEMENTED FOR IELEM>1.')

    ! CONSTANT EXPANSION OF THE FLUX
    IF(IELEM.EQ.1.AND.ISCHM.LE.3) THEN

        ! DD0 SCHEME
        IF(ISCHM.EQ.1) THEN

            Q2(1,1)=V*SIGMA+2.0D0*ABS(DX)+2.0D0*ABS(DY)+2.0D0*ABS(DZ)
            Q2(1,2)=V*QN(1)+2.0D0*ABS(DX)*XNI(1)+2.0D0*ABS(DY)*XNJ(1)+2.0D0*ABS(DZ)*XNK(1)

        ! DG0 SCHEME
        ELSEIF(ISCHM.EQ.2) THEN

            Q2(1,1)=V*SIGMA+ABS(DX)+ABS(DY)+ABS(DZ) 
            Q2(1,2)=V*QN(1)+ABS(DX)*XNI(1)+ABS(DY)*XNJ(1)+ABS(DZ)*XNK(1)

        ! AWD0 SCHEME
        ELSEIF(ISCHM.EQ.3) THEN

            CALL SNADPT3D(WX0,WY0,WZ0,1.0,1.0,1.0,XNI(1),XNJ(1),XNK(1),QN(1),DX,DY,DZ,SIGMA,V)
            Q2(1,1)=V*SIGMA+WX0(2)*ABS(DX)+WY0(2)*ABS(DY)+WZ0(2)*ABS(DZ)
            Q2(1,2)=V*QN(1)+(1.0-WX0(1))*ABS(DX)*XNI(1)+(1.0-WY0(1))*ABS(DY)*XNJ(1)+(1.0-WZ0(1))*ABS(DZ)*XNK(1)

        ENDIF

    ! LINEAR EXPANSION OF THE FLUX
    ELSEIF(IELEM.EQ.2.AND.ISCHM.LE.2) THEN

        ! DD1 SCHEME
        IF(ISCHM.EQ.1) THEN

            Q2(1,1)=V*SIGMA
            Q2(1,2)=2.0D0*CST(2)*DX
            Q2(1,3)=2.0D0*CST(2)*DY
            Q2(1,5)=2.0D0*CST(2)*DZ
            Q2(2,1)=-2.0D0*CST(2)*DX
            Q2(2,2)=V*SIGMA+6.0D0*ABS(DX)
            Q2(2,4)=2.0D0*CST(2)*DY
            Q2(2,6)=2.0D0*CST(2)*DZ
            Q2(3,1)=-2.0D0*CST(2)*DY
            Q2(3,3)=V*SIGMA+6.0D0*ABS(DY)
            Q2(3,4)=2.0D0*CST(2)*DX
            Q2(3,7)=2.0D0*CST(2)*DZ
            Q2(4,2)=-2.0D0*CST(2)*DY
            Q2(4,3)=-2.0D0*CST(2)*DX
            Q2(4,4)=V*SIGMA+6.0D0*ABS(DX)+6.0D0*ABS(DY)
            Q2(4,8)=2.0D0*CST(2)*DZ
            Q2(5,1)=-2.0D0*CST(2)*DZ
            Q2(5,5)=V*SIGMA+6.0D0*ABS(DZ)
            Q2(5,6)=2.0D0*CST(2)*DX
            Q2(5,7)=2.0D0*CST(2)*DY
            Q2(6,2)=-2.0D0*CST(2)*DZ
            Q2(6,5)=-2.0D0*CST(2)*DX
            Q2(6,6)=V*SIGMA+6.0D0*ABS(DX)+6.0D0*ABS(DZ)
            Q2(6,8)=2.0D0*CST(2)*DY
            Q2(7,3)=-2.0D0*CST(2)*DZ
            Q2(7,5)=-2.0D0*CST(2)*DY
            Q2(7,7)=V*SIGMA+6.0D0*ABS(DY)+6.0D0*ABS(DZ)
            Q2(7,8)=2.0D0*CST(2)*DX
            Q2(8,4)=-2.0D0*CST(2)*DZ
            Q2(8,6)=-2.0D0*CST(2)*DY
            Q2(8,7)=-2.0D0*CST(2)*DX
            Q2(8,8)=V*SIGMA+6.0D0*ABS(DX)+6.0D0*ABS(DY)+6.0D0*ABS(DZ)

            Q2(1,9)=V*QN(1)
            Q2(2,9)=V*QN(2)-2.0D0*CST(2)*DX*XNI(1)
            Q2(3,9)=V*QN(3)-2.0D0*CST(2)*DY*XNJ(1)
            Q2(4,9)=V*QN(4)-2.0D0*CST(2)*DX*XNI(2)-2.0D0*CST(2)*DY*XNJ(2)
            Q2(5,9)=V*QN(5)-2.0D0*CST(2)*DZ*XNK(1)
            Q2(6,9)=V*QN(6)-2.0D0*CST(2)*DX*XNI(3)-2.0D0*CST(2)*DZ*XNK(2)
            Q2(7,9)=V*QN(7)-2.0D0*CST(2)*DY*XNJ(3)-2.0D0*CST(2)*DZ*XNK(3)
            Q2(8,9)=V*QN(8)-2.0D0*CST(2)*DX*XNI(4)-2.0D0*CST(2)*DY*XNJ(4)-2.0D0*CST(2)*DZ*XNK(4)

        ! DG1 SCHEME
        ELSEIF(ISCHM.EQ.2) THEN

            Q2(1,1)=V*SIGMA+ABS(DX)+ABS(DY)+ABS(DZ)
            Q2(1,2)=CST(2)*DX
            Q2(1,3)=CST(2)*DY
            Q2(1,5)=CST(2)*DZ
            Q2(2,1)=-CST(2)*DX
            Q2(2,2)=V*SIGMA+3*ABS(DX)+ABS(DY)+ABS(DZ)
            Q2(2,4)=CST(2)*DY
            Q2(2,6)=CST(2)*DZ
            Q2(3,1)=-CST(2)*DY
            Q2(3,3)=V*SIGMA+ABS(DX)+3*ABS(DY)+ABS(DZ)
            Q2(3,4)=CST(2)*DX
            Q2(3,7)=CST(2)*DZ
            Q2(4,2)=-CST(2)*DY
            Q2(4,3)=-CST(2)*DX
            Q2(4,4)=V*SIGMA+3*ABS(DX)+3*ABS(DY)+ABS(DZ)
            Q2(4,8)=CST(2)*DZ
            Q2(5,1)=-CST(2)*DZ
            Q2(5,5)=V*SIGMA+ABS(DX)+ABS(DY)+3*ABS(DZ)
            Q2(5,6)=CST(2)*DX
            Q2(5,7)=CST(2)*DY
            Q2(6,2)=-CST(2)*DZ
            Q2(6,5)=-CST(2)*DX
            Q2(6,6)=V*SIGMA+3*ABS(DX)+ABS(DY)+3*ABS(DZ)
            Q2(6,8)=CST(2)*DY
            Q2(7,3)=-CST(2)*DZ
            Q2(7,5)=-CST(2)*DY
            Q2(7,7)=V*SIGMA+ABS(DX)+3*ABS(DY)+3*ABS(DZ)
            Q2(7,8)=CST(2)*DX
            Q2(8,4)=-CST(2)*DZ
            Q2(8,6)=-CST(2)*DY
            Q2(8,7)=-CST(2)*DX
            Q2(8,8)=V*SIGMA+3*ABS(DX)+3*ABS(DY)+3*ABS(DZ)
            
            Q2(1,9)=QN(1)*V+XNK(1)*DZ*SZ+DY*XNJ(1)*SY+DX*XNI(1)*SX
            Q2(2,9)=QN(2)*V+XNK(2)*DZ*SZ+DY*XNJ(2)*SY-CST(2)*DX*XNI(1)
            Q2(3,9)=QN(3)*V+XNK(3)*DZ*SZ+DX*XNI(2)*SX-CST(2)*DY*XNJ(1)
            Q2(4,9)=QN(4)*V+XNK(4)*DZ*SZ-CST(2)*DY*XNJ(2)-CST(2)*DX*XNI(2)
            Q2(5,9)=QN(5)*V+DY*XNJ(3)*SY+DX*XNI(3)*SX-CST(2)*XNK(1)*DZ
            Q2(6,9)=QN(6)*V+DY*XNJ(4)*SY-CST(2)*DX*XNI(3)-CST(2)*XNK(2)*DZ
            Q2(7,9)=QN(7)*V+DX*XNI(4)*SX-CST(2)*DY*XNJ(3)-CST(2)*XNK(3)*DZ
            Q2(8,9)=QN(8)*V-CST(2)*DY*XNJ(4)-CST(2)*DX*XNI(4)-CST(2)*XNK(4)*DZ

        ENDIF

    ! GENERAL POLYNOMIAL EXPANSION OF THE FLUX
    ELSE

        ! COLLISION TERM
        DO IZ=1,IELEM
        DO IY=1,IELEM
        DO IX=1,IELEM
            II=IELEM**2*(IZ-1)+IELEM*(IY-1)+IX
            Q2(II,II)=Q2(II,II)+SIGMA*V
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
                II=IELEM**2*(IZ-1)+IELEM*(IY-1)+IX
                JJ=IELEM**2*(IZ-1)+IELEM*(IY-1)+JX
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
            DO IZ=1,IELEM
            DO IX=1,IELEM
                II=IELEM**2*(IZ-1)+IELEM*(IY-1)+IX
                JJ=IELEM**2*(IZ-1)+IELEM*(JY-1)+IX
                Q2(II,JJ)=Q2(II,JJ)+FACT
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
                II=IELEM**2*(IZ-1)+IELEM*(IY-1)+IX
                JJ=IELEM**2*(JZ-1)+IELEM*(IY-1)+IX
                Q2(II,JJ)=Q2(II,JJ)+FACT
            ENDDO
            ENDDO
        ENDDO
        ENDDO

        ! SOURCE TERMS
        DO IZ=1,IELEM
        DO IY=1,IELEM
        DO IX=1,IELEM
            II=IELEM**2*(IZ-1)+IELEM*(IY-1)+IX
            Q2(II,NM+1)=Q2(II,NM+1)+QN(II)*V
        ENDDO
        ENDDO
        ENDDO
        DO IX=1,IELEM
            FACT=CST(IX)*DX*SX**IX*(WX(1)-(-1.0)**(IX-1))
            DO IY=1,IELEM
            DO IZ=1,IELEM
                II=IELEM**2*(IZ-1)+IELEM*(IY-1)+IX
                IIX=IELEM*(IZ-1)+IY
                Q2(II,NM+1)=Q2(II,NM+1)-FACT*XNI(IIX)
            ENDDO
            ENDDO
        ENDDO
        DO IY=1,IELEM
            FACT=CST(IY)*DY*SY**IY*(WY(1)-(-1.0)**(IY-1))
            DO IX=1,IELEM
            DO IZ=1,IELEM
                II=IELEM**2*(IZ-1)+IELEM*(IY-1)+IX
                IIY=IELEM*(IZ-1)+IX
                Q2(II,NM+1)=Q2(II,NM+1)-FACT*XNJ(IIY)
            ENDDO
            ENDDO
        ENDDO
        DO IZ=1,IELEM
            FACT=CST(IZ)*DZ*SZ**IZ*(WZ(1)-(-1.0)**(IZ-1))
            DO IX=1,IELEM
            DO IY=1,IELEM
                II=IELEM**2*(IZ-1)+IELEM*(IY-1)+IX
                IIZ=IELEM*(IY-1)+IX
                Q2(II,NM+1)=Q2(II,NM+1)-FACT*XNK(IIZ)
            ENDDO
            ENDDO
        ENDDO

    ENDIF

    !----
    ! SOLVE THE LINEAR SYSTEM
    !----
    CALL ALSBD(NM,1,Q2,IER,NM)
    IF(IER.NE.0) CALL XABORT('SNBTE: SINGULAR MATRIX.')

    !----
    ! CLOSURE RELATIONS
    !----

    ! CONSTANT EXPANSION OF THE FLUX
    IF(IELEM.EQ.1.AND.ISCHM.LE.3) THEN

        IF(LFIXUP.AND.(Q2(1,NM+1).LE.RLOG)) Q2(1,2)=0.0 ! FLUX FIX-UP

        ! DD SCHEME
        IF(ISCHM.EQ.1) THEN

            XNI(1)=2.0D0*Q2(1,2)-XNI(1)
            XNJ(1)=2.0D0*Q2(1,2)-XNJ(1)
            XNK(1)=2.0D0*Q2(1,2)-XNK(1)

        ! DG SCHEME
        ELSEIF(ISCHM.EQ.2) THEN

            XNI(1) = Q2(1,2) 
            XNJ(1) = Q2(1,2) 
            XNK(1) = Q2(1,2)

        ! AWD SCHEME
        ELSEIF(ISCHM.EQ.3) THEN

            XNI(1)=WX0(1)*XNI(1)+WX0(2)*Q2(1,2)
            XNJ(1)=WY0(1)*XNJ(1)+WY0(2)*Q2(1,2)
            XNK(1)=WZ0(1)*XNK(1)+WZ0(2)*Q2(1,2)

        ENDIF

        IF(LFIXUP) THEN ! FLUX FIX-UP
            IF(XNI(1).LE.RLOG) XNI(1)=0.0
            IF(XNJ(1).LE.RLOG) XNJ(1)=0.0
            IF(XNK(1).LE.RLOG) XNK(1)=0.0
        ENDIF

    ! LINEAR EXPANSION OF THE FLUX
    ELSEIF(IELEM.EQ.2.AND.ISCHM.LE.2) THEN

        ! DD SCHEME
        IF(ISCHM.EQ.1) THEN

            XNI(1)=XNI(1)+SX*2*CST(2)*Q2(2,9)
            XNI(2)=XNI(2)+SX*2*CST(2)*Q2(4,9)
            XNI(3)=XNI(3)+SX*2*CST(2)*Q2(6,9)
            XNI(4)=XNI(4)+SX*2*CST(2)*Q2(8,9)
            XNJ(1)=XNJ(1)+SY*2*CST(2)*Q2(3,9)
            XNJ(2)=XNJ(2)+SY*2*CST(2)*Q2(4,9)
            XNJ(3)=XNJ(3)+SY*2*CST(2)*Q2(7,9)
            XNJ(4)=XNJ(4)+SY*2*CST(2)*Q2(8,9)
            XNK(1)=XNK(1)+SZ*2*CST(2)*Q2(5,9)
            XNK(2)=XNK(2)+SZ*2*CST(2)*Q2(6,9)
            XNK(3)=XNK(3)+SZ*2*CST(2)*Q2(7,9)
            XNK(4)=XNK(4)+SZ*2*CST(2)*Q2(8,9)

        ! DG SCHEME
        ELSEIF(ISCHM.EQ.2) THEN

            XNI(1)=Q2(1,9)+SX*CST(2)*Q2(2,9)
            XNI(2)=Q2(3,9)+SX*CST(2)*Q2(4,9)
            XNI(3)=Q2(5,9)+SX*CST(2)*Q2(6,9)
            XNI(4)=Q2(7,9)+SX*CST(2)*Q2(8,9)
            XNJ(1)=Q2(1,9)+SY*CST(2)*Q2(3,9)
            XNJ(2)=Q2(2,9)+SY*CST(2)*Q2(4,9)
            XNJ(3)=Q2(5,9)+SY*CST(2)*Q2(7,9)
            XNJ(4)=Q2(6,9)+SY*CST(2)*Q2(8,9)
            XNK(1)=Q2(1,9)+SZ*CST(2)*Q2(5,9)
            XNK(2)=Q2(2,9)+SZ*CST(2)*Q2(6,9)
            XNK(3)=Q2(3,9)+SZ*CST(2)*Q2(7,9)
            XNK(4)=Q2(4,9)+SZ*CST(2)*Q2(8,9)

        ENDIF

    ! GENERAL POLYNOMIAL EXPANSION OF THE FLUX
    ELSE
        XNI(:NMX)=WX(1)*XNI(:NMX)
        XNJ(:NMY)=WY(1)*XNJ(:NMY)
        XNK(:NMZ)=WZ(1)*XNK(:NMZ)
        DO IX=1,IELEM
            FACT=CST(IX)*SX**(IX-1)*WX(IX+1)
            DO IY=1,IELEM
            DO IZ=1,IELEM
                II=IELEM**2*(IZ-1)+IELEM*(IY-1)+IX
                IIX=IELEM*(IZ-1)+IY
                XNI(IIX)=XNI(IIX)+FACT*Q2(II,NM+1)
            ENDDO
            ENDDO
        ENDDO
        DO IY=1,IELEM
            FACT=CST(IY)*SY**(IY-1)*WY(IY+1)
            DO IX=1,IELEM
            DO IZ=1,IELEM
                II=IELEM**2*(IZ-1)+IELEM*(IY-1)+IX
                IIY=IELEM*(IZ-1)+IX
                XNJ(IIY)=XNJ(IIY)+FACT*Q2(II,NM+1)
            ENDDO
            ENDDO
        ENDDO
        DO IZ=1,IELEM
            FACT=CST(IZ)*SZ**(IZ-1)*WZ(IZ+1)
            DO IX=1,IELEM
            DO IY=1,IELEM
                II=IELEM**2*(IZ-1)+IELEM*(IY-1)+IX
                IIZ=IELEM*(IY-1)+IX
                XNK(IIZ)=XNK(IIZ)+FACT*Q2(II,NM+1)
            ENDDO
            ENDDO
        ENDDO
    ENDIF
    RETURN
  END SUBROUTINE SNBTE3
END MODULE SNBTE_MOD
