!DECK THMDFM
      SUBROUTINE THMDFM(PCOOL,VCOOL,HMAVG,HD,TL,TSAT,IDFM,EPS,XFL,RHO,RHOL,RHOG, VGJ, VGJprime, C0, HLV)
!
!-----------------------------------------------------------------------
!
! Purpose:
! Drift-flux Model for the computation of thermohydraulics parameters in two-phase flow
!
!Copyright:
! Copyright (C) 2025 Ecole Polytechnique de Montreal.
!
!Author(s): 
! M. Bellier
!
!Parameters: input
! PCOOL   pressure in Pascal
! VCOOL   coolant velocity in m/s
! HMAVG   averaged enthalpy
! HD      hydraulic diameter in m
! TL      liquid temperature in K
! TSAT    saturation temperature in K
! IDFM    flag indicating if the drift flux model is to be used 
!         (0=HEM1(no drift velocity)/1=EPRI/2=MODEBSTION/3=GERAMP/4=CHEXAL) 
! EPS     input coolant void fraction
!
!
!Parameters: output
! XFL     coolant flow quality
! RHO     coolant density in Kg/m^3
! RHOL    liquid density in kg/m^3
! RHOG    vapour density in kg/m^3
! VGJ     drift velocity
! C0      concentration parameter
! VGJprime 
! HLV     delta between liquid and vapour enthaply
!
!-----------------------------------------------------------------------
!
!----
!  SUBROUTINE ARGUMENTS
!----
      REAL PCOOL,VCOOL,HMAVG,HD,TL,TSAT,EPS,XFL,RHO,RHOL,RHOG, VGJ, VGJprime, C0, HLV
      INTEGER IDFM
!----
!  LOCAL VARIABLES
!----
      REAL EPSold, ERREPS, VLIQ, VVAP, TCALO, HLSAT, HGSAT, ZMUL, ZMUG, CPL, CPG, ZKL, ZKG, ZMU, REY
      INTEGER NITER
!----
! INITIALIZE VARIABLES
!----
      VGJ = 0
      C0 = 1
      VGJprime = 0

!----
!  MAIN LOOP
!----
     NITER=0
     ERREPS=1

   10 CONTINUE
!----
!  SAVE THE OLD EPSILON VALUE
!----
      EPSold = EPS 
      NITER = NITER+1

!----
! TEST ON ERR EPS
!----
      IF (NITER .GT. 150) GOTO 20
      IF (ERREPS .LT. 1E-8) GOTO 20 

!----
!  COMPUTE DENSITIES
!----
      TCALO=EPS*TSAT+(1.0-EPS)*TL
      CALL THMTX(TCALO,0.0,RHOL,HLSAT,ZKL,ZMUL,CPL)
      CALL THMTX(TCALO,1.0,RHOG,HGSAT,ZKG,ZMUG,CPG)

      RHO = RHOL*(1 - EPS)+ EPS*RHOG
    
!----
!  COMPUTE PHASES VELOCITIES AND REYNOLDS
!----
      VLIQ = VCOOL - (1.0/(1.0- EPS) - RHOL/RHO) *VGJprime
      VVAP = VCOOL + RHOL/RHO *VGJprime
      ZMU = (ZMUL*ZMUG/ (ZMUL*(1.0-EPS) + ZMUG*EPS))
      REY = RHO * ABS(VCOOL) * HD / ZMU

!----
!  COMPUTE FLOW QUALITY
!----

      IF (HLSAT .GT. HMAVG) THEN 
        XFL = 0
      ELSE IF (HMAVG .GT. HGSAT) THEN
        XFL = 1
      ELSE     
        XFL = (HMAVG - HLSAT)/(HGSAT - HLSAT)
      ENDIF

!----
!  COMPUTE VGJ, VGJprime AND C0 AFTER CHOSEN CORRELATION
!----
      CALL THMVGJ(VCOOL, RHO, PCOOL, ZMU, XFL, HD, RHOG, RHOL, EPS, IDFM, VGJ, C0)
      VGJprime = VGJ + (C0-1)*VCOOL

!----
!  COMPUTE HLV
!----
      HLV=HGSAT-HLSAT
!----
!  COMPUTE NEW EPS VALUE
!----
      IF (XFL.EQ.0) THEN 
        EPS = 0
      ELSE IF (XFL.EQ.1) THEN
        EPS = 1
      ELSE 
        EPS = XFL / (C0 * (XFL + (RHOG/RHOL) * (1 - XFL)) + (RHOG * VGJ) / (RHOL * VCOOL))
      ENDIF
!----
!  COMPUTE DELTA BETWEEN EPSold AND EPS
!----
    ERREPS = ABS(EPSold - EPS)
   GOTO 10


!----
! EXIT LOOP
!----
   20 CONTINUE

      IF (NITER.GT.150) THEN
        PRINT *, 'THMDFM: Maximum number of iterations reached (150)'
      ELSE
        PRINT *, 'THMDFM: Convergence reached in I = ', NITER, 'iterations'
      ENDIF
END
