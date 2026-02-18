SUBROUTINE THMPV(SPEED, POULET, VCOOL, DCOOL, PCOOL, TCOOL, MUT, XFL, HD, NZ, HZ, EPS, RHOL, RHOG, VGJ, IDFM, ACOOL)
!
!-----------------------------------------------------------------------
!
!Purpose:
! Update the pressure and velocity vectors in the THM model to model the
! pressure drop and the velocity of the fluid in the channel
!
!Copyright:
! Copyright (C) 2025 Ecole Polytechnique de Montreal
!
!Author(s): C. Huet
! 02/2025: C. Huet - Creation
!
!Parameters: input
! SPEED   inlet velocity of the fluid in the channel
! POULET  Pressure at the outlet
! VCOOL   velocity of the fluid in the channel
! DCOOL   density of the fluid in the channel
! PCOOL   pressure of the fluid in the channel
! TCOOL   temperature of the fluid in the channel
! MUT     dynamic viscosity of the fluid in the channel
! XFL     quality of the fluid in the channel
! HD      hydraulic diameter of the channel
! NZ      number of nodes in the channel
! HZ      height of the channel
! EPS     coolant void fraction in the channel
! RHOL    density of the liquid fraction
! RHOG    density of the vapour fraction
! VGJ     drift velocity in the channel
! IDFM    flag for the use of the drift flux model
! ACOOL   cross-sectional area of the channel
!
!Parameters: output
! VCOOL   velocity of the fluid in the channel
! PCOOL   pressure of the fluid in the channel
!
!-----------------------------------------------------------------------
!
    USE GANLIB
    IMPLICIT NONE
!----
!   SUBROUTINE ARGUMENTS
!----
    INTEGER NZ, IDFM
    REAL SPEED, POULET, VCOOL(NZ), DCOOL(NZ), PCOOL(NZ), TCOOL(NZ), MUT(NZ), XFL(NZ)
    REAL HZ(NZ),VGJ(NZ),RHOL(NZ), RHOG(NZ), EPS(NZ), HD(NZ), ACOOL(NZ) 
!----
!   LOCAL VARIABLES
!----
    REAL g
    REAL(kind=8), ALLOCATABLE, DIMENSION(:,:) :: A

    INTEGER K, I, J, IER
    REAL PHIL0, TPMULT, TPMULT0
    REAL REY, REY0, FRIC,FRIC0,DELTA, UL
    REAL CP11, H11, K11, RHO11, MUL

    g = 9.81 !gravity
    ALLOCATE(A(2*NZ,2*NZ+1))
    FORALL (I=1:2*NZ, J=1:2*NZ+1) A(I, J) = 0.0

!----
!   MATRIX FILLING FOR THE PRESSURE AND VELOCITY CALCULATION
!----
!   BOTTOM OF THE CHANNEL
!----
    PRINT *, 'THMPV: Filling the matrix for pressure and velocity calculation'
    PRINT *, 'THMPV: NZ = ', NZ
    PRINT *, 'POULET = ', POULET
    DO K = 1, NZ 
        IF (K .EQ. 1) THEN
            IF(IDFM.GT.0) THEN
    !       COMPUTE MUL, UL and Reynolds AT K
                CALL THMTX(TCOOL(K), 0.0, RHO11, H11, K11, MUL, CP11)
                UL = VCOOL(K) - (EPS(K) / (1.0 - EPS(K)))*RHOG(K)/DCOOL(K) * VGJ(K)
                REY0 = ABS(UL*RHOL(K)) * HD(K) / MUL
    !       COMPUTE MUL, UL and Reynolds AT K+1
                CALL THMTX(TCOOL(K+1), 0.0, RHO11, H11, K11, MUL, CP11)
                UL = VCOOL(K+1) - (EPS(K+1) / (1.0 - EPS(K+1)))*RHOG(K+1)/DCOOL(K+1) * VGJ(K+1)
                REY = ABS(UL*RHOL(K+1)) * HD(K+1) / MUL
            ELSE 
                REY = ABS(VCOOL(K+1)*DCOOL(K+1)) * (1.0 - XFL(K+1)) * HD(K+1) / MUT(K+1)
                REY0 = ABS(VCOOL(K)*DCOOL(K)) * (1.0 - XFL(K)) * HD(K) / MUT(K)
            ENDIF
            
            
            CALL THMFRI(REY,EPS(K+1),HD(K+1),FRIC)!MUT à isoler vapeur/liquide : passer par THMTX(TCOOL, X=0)
            CALL THMFRI(REY0,EPS(K),HD(K),FRIC0)

            IF (XFL(K) .GT. 0.0) THEN
                CALL THMPLO(PCOOL(K), XFL(K), PHIL0)
                TPMULT0 = PHIL0
                CALL THMPLO(PCOOL(K+1), XFL(K+1), PHIL0)
                TPMULT = PHIL0
            ELSE
                TPMULT = 1.0
                TPMULT0 = 1.0
            ENDIF
            A(1,1) = 1.0
!   MOMENTUM CONSERVATION EQUATION
            IF (IDFM .GT. 0) THEN
                DELTA = ((EPS(K)/1-EPS(K))*RHOL(K)*RHOG(K)/DCOOL(K)*VGJ(K)**2) - &
            ((EPS(K+1)/1-EPS(K+1))*RHOL(K+1)*RHOG(K+1)/DCOOL(K+1)*VGJ(K+1)* &
            ACOOL(K+1)/ACOOL(K)**2)
            ELSE
                DELTA = 0.0
            ENDIF
            A(K+NZ,K) = - (VCOOL(K)*DCOOL(K))*(1.0 - (TPMULT0*FRIC0*HZ(K))/(2.0*HD(K)))
            A(K+NZ,K+1) = (VCOOL(K+1)*DCOOL(K+1))*(1.0 + (TPMULT*FRIC*HZ(K))/ &
             (2.0*HD(K+1)))*ACOOL(K+1)/ACOOL(K)
            A(K+NZ, 2*NZ+1) =  - ((DCOOL(K+1)* HZ(K+1)*ACOOL(K+1)/ACOOL(K) + DCOOL(K)* HZ(K)) &
               * g ) /2 + DELTA
            A(K+NZ,K-1+NZ) = 0.0
            A(K+NZ,K+NZ) = -1.0
            A(K+NZ,K+1+NZ) = ACOOL(K+1)/ACOOL(K)

!    MASS CONSERVATION EQUATION
            A(1, 2*NZ+1) = SPEED

!----
!   TOP OF THE CHANNEL
!----
        ELSE IF (K .EQ. NZ) THEN
!   MASS CONSERVATION EQUATION
            A(K,K-1) = - DCOOL(K-1)*ACOOL(K-1)/ACOOL(K)
            A(K,K) = DCOOL(K)
!   MOMENTUM CONSERVATION EQUATION
            A(K, 2*NZ+1) = 0.0
            A(2*NZ, 2*NZ+1) = POULET
            A(2*NZ, 2*NZ) = 1.0
!----
!   MIDDLE OF THE CHANNEL
!----
        ELSE
            IF (IDFM.GT.0) THEN
!       COMPUTE MUL, UL and Reynolds AT K
                CALL THMTX(TCOOL(K), 0.0, RHO11, H11, K11, MUL, CP11)
                UL = VCOOL(K) - (EPS(K) / (1.0 - EPS(K)))*RHOG(K)/DCOOL(K) * VGJ(K)
                REY0 = ABS(UL*RHOL(K)) * HD(K) / MUL
!       COMPUTE MUL, UL and Reynolds AT K+1
                CALL THMTX(TCOOL(K+1), 0.0, RHO11, H11, K11, MUL, CP11)
                UL = VCOOL(K+1) - (EPS(K+1) / (1.0 - EPS(K+1)))*RHOG(K+1)/DCOOL(K+1) * VGJ(K+1)
                REY = ABS(UL*RHOL(K+1)) * HD(K+1) / MUL
            ELSE
                REY = ABS(VCOOL(K+1)*DCOOL(K+1)) * (1.0 - XFL(K+1)) * HD(K+1) / MUT(K+1)
                REY0 = ABS(VCOOL(K)*DCOOL(K)) * (1.0 - XFL(K)) * HD(K) / MUT(K)
            ENDIF        
            CALL THMFRI(REY,EPS(K+1),HD(K+1),FRIC)
            CALL THMFRI(REY0,EPS(K),HD(K),FRIC0)

            IF (XFL(K) .GT. 0.0) THEN
                CALL THMPLO(PCOOL(K+1), XFL(K+1), PHIL0)
                TPMULT = PHIL0
                CALL THMPLO(PCOOL(K), XFL(K), PHIL0)
                TPMULT0 = PHIL0
            ELSE
                TPMULT = 1.0
                TPMULT0 = 1.0
            ENDIF
!   MASS CONSERVATION EQUATION
            A(K,K-1) = - DCOOL(K-1)*ACOOL(K-1)/ACOOL(K)
            A(K,K) = DCOOL(K)
            A(K,K+1) = 0.0
            A(K, 2*NZ+1) = 0.0 
!----
!   MOMENTUM CONSERVATION EQUATION  
!----
            IF (IDFM .GT. 0) THEN
                DELTA = ((EPS(K)/1-EPS(K))*RHOL(K)*RHOG(K)/DCOOL(K)*VGJ(K)**2) - &
            ((EPS(K+1)/1-EPS(K+1))*RHOL(K+1)*RHOG(K+1)/DCOOL(K+1)*VGJ(K+1)**2*ACOOL(K+1) &
            /ACOOL(K))
            ELSE
                DELTA = 0.0
            ENDIF
            A(K+NZ,K) = - (VCOOL(K)*DCOOL(K))*(1.0 - (TPMULT0*FRIC0*HZ(K))/(2.0*HD(K)))
            A(K+NZ,K+1) = (VCOOL(K+1)*DCOOL(K+1))*(1.0 + (TPMULT*FRIC*HZ(K))/ &
             (2.0*HD(K+1)))*ACOOL(K+1)/ACOOL(K)
            A(K+NZ, 2*NZ+1) =  - ((DCOOL(K+1)* HZ(K+1)*ACOOL(K+1)/ACOOL(K) + DCOOL(K)* &
             HZ(K)) * g ) /2 + DELTA
            A(K+NZ,K-1+NZ) = 0.0
            A(K+NZ,K+NZ) = -1.0
            A(K+NZ,K+1+NZ) = ACOOL(K+1)/ACOOL(K)
        ENDIF
    END DO
!----
!   SOLVING THE LINEAR SYSTEM
!----
    call ALSBD(2*NZ, 1, A, IER, 2*NZ)

    if (IER /= 0) CALL XABORT('THMPV: SINGULAR MATRIX.')
!----
!   RECOVER THE PRESSURE AND VELOCITY VECTORS
!----
    DO K = 1, NZ
        VCOOL(K) = REAL(A(K, 2*NZ+1))
        PCOOL(K) = REAL(A(K+NZ, 2*NZ+1))
    END DO

    DEALLOCATE(A)
    
    RETURN
    END