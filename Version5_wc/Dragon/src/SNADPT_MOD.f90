!
!---------------------------------------------------------------------
!
!Purpose:
!  Compute the weighting factors for the adaptive weighted diamond
!  difference (AWDD) scheme in a 1D/2D/3D Cartesian cell within the
!  discrete ordinates (SN) framework.
!
! Copyright:
!  Copyright (C) 2025 Ecole Polytechnique de Montreal
!  This library is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License as published by the Free Software Foundation; either
!  version 2.1 of the License, or (at your option) any later version
!
! Author(s): C. Bienvenue
!
!---------------------------------------------------------------------
!
MODULE SNADPT_MOD

CONTAINS

  SUBROUTINE SNADPT1D(W1,TB1,XN1,Q,D1,SIGMA,V)
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    !   Compute the weighting factors for the adaptive weighted diamond
    !   difference (AWDD) scheme in a 1D Cartesian cell (Boltzmann transport
    !   equation).
    !
    ! Parameters: input
    !   TB1     ratio of stopping powers or 1.0 for spatial dimension #1.
    !   XN1     incoming/outgoing boundary flux in space along axis #1.
    !   Q       angular source term.
    !   D1      factor with first direction cosines (mu) along axis #1.
    !   SIGMA   macroscopic total cross section in the cell.
    !   V       cell volume.
    !
    ! Parameters: output
    !   W1      weighting factors for the AWDD scheme along axis #1.
    !
    !-----------------------------------------------------------------------

    ! SUBROUTINE ARGUMENTS
    IMPLICIT NONE
    REAL, INTENT(IN) :: TB1,D1,SIGMA,V
    DOUBLE PRECISION, INTENT(IN) :: Q,XN1
    REAL, INTENT(OUT) :: W1(2)

    ! LOCAL VARIABLES
    REAL :: P1,U1,F1
    DOUBLE PRECISION :: PHI_DD
    REAL, PARAMETER :: B=2.0

    ! DIAMOND DIFFERENCE FLUX
    PHI_DD = (Q*V+2*TB1*ABS(D1)*XN1)/((TB1+1)*ABS(D1)+SIGMA*V)

    ! CALCULATION OF P1
    IF(PHI_DD.NE.0) THEN
        U1 = REAL((XN1-PHI_DD)/PHI_DD)
        F1=2.0*B*ABS(U1)
        IF(F1.LE.1) THEN
            P1=1.0
        ELSE
            P1=1.0/F1
        ENDIF
    ELSE
        P1=0.0
    ENDIF

    ! COMPUTE WEIGHTING FACTORS
    W1(1)=-P1
    W1(2)=1+P1

    RETURN
  END SUBROUTINE SNADPT1D

  SUBROUTINE SNADPT2D(W1,W2,TB1,TB2,XN1,XN2,Q,D1,D2,SIGMA,V)
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    !   Compute the weighting factors for the adaptive weighted diamond
    !   difference (AWDD) scheme in a 2D Cartesian cell (Boltzmann transport
    !   equation) or 1D Cartesian cell (Boltzmann Fokker-Planck transport
    !   equation)
    !
    ! Parameters: input
    !   TB1     ratio of stopping powers or 1.0 for spatial dimension #1.
    !   TB2     ratio of stopping powers or 1.0 for spatial dimension #2.
    !   XN1     incoming/outgoing boundary flux in space along axis #1.
    !   XN2     incoming/outgoing boundary flux in space along axis #2.
    !   Q       angular source term.
    !   D1      factor with first direction cosines (mu) along axis #1.
    !   D2      factor with second direction cosines (eta) along axis #2.
    !   SIGMA   macroscopic total cross section in the cell.
    !   V       cell volume.
    !
    ! Parameters: output
    !   W1      weighting factors for the AWDD scheme along axis #1.
    !   W2      weighting factors for the AWDD scheme along axis #2.
    !
    !-----------------------------------------------------------------------

    ! SUBROUTINE ARGUMENTS
    IMPLICIT NONE
    REAL, INTENT(IN) :: TB1,TB2,D1,D2,SIGMA,V
    DOUBLE PRECISION, INTENT(IN) :: Q,XN1,XN2
    REAL, INTENT(OUT) :: W1(2),W2(2)

    ! LOCAL VARIABLES
    REAL :: P1,P2,U1,U2,F1,F2
    DOUBLE PRECISION :: PHI_DD
    REAL, PARAMETER :: B=2.0

    ! DIAMOND DIFFERENCE FLUX
    PHI_DD = (Q*V+2*TB1*ABS(D1)*XN1+2*TB2*ABS(D2)*XN2)/((TB1+1)*ABS(D1)+(TB2+1)*ABS(D2)+SIGMA*V)

    ! CALCULATION OF P1, P2
    IF(PHI_DD.NE.0) THEN
        U1 = REAL((XN1-PHI_DD)/PHI_DD)
        U2 = REAL((XN2-PHI_DD)/PHI_DD)
        F1=2.0*B*ABS(U1)
        F2=2.0*B*ABS(U2)
        IF(F1.LE.1) THEN
            P1=1.0
        ELSE
            P1=1.0/F1
        ENDIF
        IF(F2.LE.1) THEN
            P2=1.0
        ELSE
            P2=1.0/F2
        ENDIF
    ELSE
        P1=0.0
        P2=0.0
    ENDIF

    ! COMPUTE WEIGHTING FACTORS
    W1(1)=-P1
    W1(2)=1+P1
    W2(1)=-P2
    W2(2)=1+P2

    RETURN
  END SUBROUTINE SNADPT2D

  SUBROUTINE SNADPT3D(W1,W2,W3,TB1,TB2,TB3,XN1,XN2,XN3,Q,D1,D2,D3,SIGMA,V)
    !
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    !   Compute the weighting factors for the adaptive weighted diamond
    !   difference (AWDD) scheme in a 3D Cartesian cell (Boltzmann transport
    !   equation) or 2D Cartesian cell (Boltzmann Fokker-Planck transport
    !   equation)
    !
    ! Parameters: input
    !   TB1     ratio of stopping powers or 1.0 for spatial dimension #1.
    !   TB2     ratio of stopping powers or 1.0 for spatial dimension #2.
    !   TB3     ratio of stopping powers or 1.0 for spatial dimension #3.
    !   XN1     incoming/outgoing boundary flux in space along axis #1.
    !   XN2     incoming/outgoing boundary flux in space along axis #2.
    !   XN3     incoming/outgoing boundary flux in space along axis #3.
    !   Q       angular source term.
    !   D1      factor with first direction cosines (mu) along axis #1.
    !   D2      factor with second direction cosines (eta) along axis #2.
    !   D3      factor with third direction cosines (xi) along axis #3.
    !   SIGMA   macroscopic total cross section in the cell.
    !   V       cell volume.
    !
    ! Parameters: output
    !   W1      weighting factors for the AWDD scheme along axis #1.
    !   W2      weighting factors for the AWDD scheme along axis #2.
    !   W3      weighting factors for the AWDD scheme along axis #3.
    !
    !-----------------------------------------------------------------------

    ! SUBROUTINE ARGUMENTS
    IMPLICIT NONE
    REAL, INTENT(IN) :: TB1,TB2,TB3,D1,D2,D3,SIGMA,V
    DOUBLE PRECISION, INTENT(IN) :: Q,XN1,XN2,XN3
    REAL, INTENT(OUT) :: W1(2),W2(2),W3(2)

    ! LOCAL VARIABLES
    REAL :: P1,P2,P3,U1,U2,U3,F1,F2,F3
    DOUBLE PRECISION :: PHI_DD
    REAL, PARAMETER :: B=2.0

    ! DIAMOND DIFFERENCE FLUX
    PHI_DD = (Q*V+2*TB1*ABS(D1)*XN1+2*TB2*ABS(D2)*XN2+2*TB3*ABS(D3)*XN3)/((TB1+1)*ABS(D1)+(TB2+1)*ABS(D2)+(TB3+1)*ABS(D3)+SIGMA*V)

    ! CALCULATION OF P1, P2, P3
    IF(PHI_DD.NE.0) THEN
        U1 = REAL((XN1-PHI_DD)/PHI_DD)
        U2 = REAL((XN2-PHI_DD)/PHI_DD)
        U3 = REAL((XN3-PHI_DD)/PHI_DD)
        F1=2.0*B*ABS(U1)
        F2=2.0*B*ABS(U2)
        F3=2.0*B*ABS(U3)
        IF(F1.LE.1) THEN
            P1=1.0
        ELSE
            P1=1.0/F1
        ENDIF
        IF(F2.LE.1) THEN
            P2=1.0
        ELSE
            P2=1.0/F2
        ENDIF
        IF(F3.LE.1) THEN
            P3=1.0
        ELSE
            P3=1.0/F3
        ENDIF
    ELSE
        P1=0.0
        P2=0.0
        P3=0.0
    ENDIF

    ! COMPUTE WEIGHTING FACTORS
    W1(1)=-P1
    W1(2)=1+P1
    W2(1)=-P2
    W2(2)=1+P2
    W3(1)=-P3
    W3(2)=1+P3

    RETURN
  END SUBROUTINE SNADPT3D

  SUBROUTINE SNADPT4D(W1,W2,W3,W4,TB1,TB2,TB3,TB4,XN1,XN2,XN3,XN4,Q,D1,D2,D3,D4,SIGMA,V)
    !
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    !   Compute the weighting factors for the adaptive weighted diamond
    !   difference (AWDD) scheme in a 3D Cartesian cell (Boltzmann Fokker-
    !   Planck transport equation)
    !
    ! Parameters: input
    !   TB1     ratio of stopping powers or 1.0 for spatial dimension #1.
    !   TB2     ratio of stopping powers or 1.0 for spatial dimension #2.
    !   TB3     ratio of stopping powers or 1.0 for spatial dimension #3.
    !   TB4     ratio of stopping powers or 1.0 for spatial dimension #4.
    !   XN1     incoming/outgoing boundary flux in space along axis #1.
    !   XN2     incoming/outgoing boundary flux in space along axis #2.
    !   XN3     incoming/outgoing boundary flux in space along axis #3.
    !   XN4     incoming/outgoing boundary flux in space along axis #4.
    !   Q       angular source term.
    !   D1      factor with first direction cosines (mu) along axis #1.
    !   D2      factor with second direction cosines (eta) along axis #2.
    !   D3      factor with third direction cosines (xi) along axis #3.
    !   D4      factor with fourth direction cosines along axis #4.
    !   SIGMA   macroscopic total cross section in the cell.
    !   V       cell volume.
    !
    ! Parameters: output
    !   W1      weighting factors for the AWDD scheme along axis #1.
    !   W2      weighting factors for the AWDD scheme along axis #2.
    !   W3      weighting factors for the AWDD scheme along axis #3.
    !   W4      weighting factors for the AWDD scheme along axis #4.
    !
    !-----------------------------------------------------------------------

    ! SUBROUTINE ARGUMENTS
    IMPLICIT NONE
    REAL, INTENT(IN) :: TB1,TB2,TB3,TB4,D1,D2,D3,D4,SIGMA,V
    DOUBLE PRECISION, INTENT(IN) :: Q,XN1,XN2,XN3,XN4
    REAL, INTENT(OUT) :: W1(2),W2(2),W3(2),W4(2)

    ! LOCAL VARIABLES
    REAL :: P1,P2,P3,P4,U1,U2,U3,U4,F1,F2,F3,F4
    DOUBLE PRECISION :: PHI_DD
    REAL, PARAMETER :: B=2.0

    ! DIAMOND DIFFERENCE FLUX
    PHI_DD = (Q*V+2*TB1*ABS(D1)*XN1+2*TB2*ABS(D2)*XN2+2*TB3*ABS(D3)*XN3+2*TB4*ABS(D4)*XN4)/ &
    & ((TB1+1)*ABS(D1)+(TB2+1)*ABS(D2)+(TB3+1)*ABS(D3)+(TB4+1)*ABS(D4)+SIGMA*V)

    ! CALCULATION OF P1, P2, P3, P4
    IF(PHI_DD.NE.0) THEN
        U1 = REAL((XN1-PHI_DD)/PHI_DD)
        U2 = REAL((XN2-PHI_DD)/PHI_DD)
        U3 = REAL((XN3-PHI_DD)/PHI_DD)
        U4 = REAL((XN4-PHI_DD)/PHI_DD)
        F1=2.0*B*ABS(U1)
        F2=2.0*B*ABS(U2)
        F3=2.0*B*ABS(U3)
        F4=2.0*B*ABS(U4)
        IF(F1.LE.1) THEN
            P1=1.0
        ELSE
            P1=1.0/F1
        ENDIF
        IF(F2.LE.1) THEN
            P2=1.0
        ELSE
            P2=1.0/F2
        ENDIF
        IF(F3.LE.1) THEN
            P3=1.0
        ELSE
            P3=1.0/F3
        ENDIF
        IF(F4.LE.1) THEN
            P4=1.0
        ELSE
            P4=1.0/F4
        ENDIF
    ELSE
        P1=0.0
        P2=0.0
        P3=0.0
        P4=0.0
    ENDIF

    ! COMPUTE WEIGHTING FACTORS
    W1(1)=-P1
    W1(2)=1+P1
    W2(1)=-P2
    W2(2)=1+P2
    W3(1)=-P3
    W3(2)=1+P3
    W4(1)=-P4
    W4(2)=1+P4

    RETURN
  END SUBROUTINE SNADPT4D
END MODULE SNADPT_MOD
