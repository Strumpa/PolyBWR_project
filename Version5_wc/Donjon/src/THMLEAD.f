*DECK THMLEAD
      SUBROUTINE THMLEAD(ITIME,I,J,K,K0,MFLOW,HMAVG,ENT,HD,IHCONV,
     > KHCONV,ISUBM,RADCL,ZF,PHI,TCALO,RHO,RHOLAV,TSCLAD)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Adaptation of THMH2O.f for convection of liquid lead using Gnielinski
* correlation.
*
*Copyright:
* Copyright (C) 2026 Ecole Polytechnique de Montreal.
*
*Author(s): Francesco Pepe
* 
*Parameters: input
* ITIME   type of calculation  (0=steady-state; 1=transient.
* I       position of channel alon X-axis
* J       position of channel alon Y-axis
* K       position of channel alon Z-axis
* K0      onset of nuclear boiling point
* MFLOW   massic coolant flow rate in Kg/m^2/s
* HMAVG   averaged enthalpy
* ENT     four values of enthalpy in J/Kg to be used in Gaussian
*         integration
* HD      hydraulic diameter in m
* IHCONV  flag indicating HCONV chosen (0=default/1=user-provided).
* KHCONV  fixed user-provided HCONV value in W/m^2/K.
* ISUBM   subcooling model (0: one-phase; 1: Jens-Lottes model;
*         2: Saha-Zuber model). [ERROR IF NOT 0]
* RADCL   outer clad radius in m
* ZF      parameters used to compute heat flux on clad surface in
*         transient cases.
* PHI     heat flow exchanged between clad and fluid in W/m^2.
*         Given in steady-state cases.
*
*Parameters: output
* PHI     heat flow exchanged between clad and fluid in W/m^2.
*         Computed in transient cases.
* TCALO   coolant temperature in K
* RHO     coolant density in Kg/m^3
* RHOLAV  liquid density in kg/m^3
* TSCLAD  clad temperature in K
*
*-----------------------------------------------------------------------
*
      USE FREELFR, only: THMLT, THMLH, CHECKLIM, THMLMBT 
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER I,J,K,K0,IHCONV,ISUBM
      REAL MFLOW,HMAVG,ENT(4),HD,KHCONV,RADCL,ZF(2),PHI,TCALO,RHO,
     > RHOLAV,TSCLAD
*----
*  LOCAL VARIABLES
*----
      REAL W(4),HL(4)
      CHARACTER HSMG*131
*----
*  SAVE VARIABLES
*----
      SAVE W
      DATA W /0.347855,0.652145,0.652145,0.347855/
*----
*  COMPUTE THE DENSITY AND TEMPERATURE OF THE LIQUID
*----
      IF(HMAVG.LT.0.0) CALL XABORT('THMLEAD: NEGATIVE INPUT ENTHALPY.')
      IF(ISUBM.NE.0) CALL XABORT('THMLEAD: NOT A ONE PHASE FLOW.')
      CALL THMLMBT(TM,TSAT,TMAX)
      HL(1)=ENT(1)
      HL(2)=ENT(2)
      HL(3)=ENT(3)
      HL(4)=ENT(4)
      CALL THMLH(HL(1),R11,TL1)
      CALL THMLH(HL(2),R11,TL2)
      CALL THMLH(HL(3),R11,TL3)
      CALL THMLH(HL(4),R11,TL4)
      CALL THMLT(TL1,RHO1,R2,R3,R4,CP1)
      CALL THMLT(TL2,RHO2,R2,R3,R4,CP2)
      CALL THMLT(TL3,RHO3,R2,R3,R4,CP3)
      CALL THMLT(TL4,RHO4,R2,R3,R4,CP4)
      TL=0.5*(W(1)*TL1+W(2)*TL2+W(3)*TL3+W(4)*TL4)
      CALL CHECKLIM(TL)
      RHOLAV=0.5*(W(1)*RHO1+W(2)*RHO2+W(3)*RHO3+W(4)*RHO4)
      CPLAV=0.5*(W(1)*CP1+W(2)*CP2+W(3)*CP3+W(4)*CP4)
*----
*  COMPUTE THE FLUID PROPERTIES
*  RHO: fluid density
*  REL: Reynolds number of liquid phase
*  PRL: Prandtl number of liquid phase
*----
*     One phase liquid
      TB=TSAT
      IF(TL.LT.TB) THEN
        TCALO=TL
      ELSE
        TCALO=TB
      ENDIF
      CALL THMLT(TCALO,R1,R2,ZKONE,ZMUONE,CPONE)
      RHO=RHOLAV
      REL=MFLOW*HD/ZMUONE
      PRL=ZMUONE*CPONE/ZKONE
      ZKL=ZKONE
*----
*  THERMAL EXCHANGE BETWEEN CLAD AND FLUID USING THE DITTUS AND BOELTER
*  CORRELATION (SINGLE PHASE) 
*----
      IF(IHCONV.EQ.0) THEN
        ITER=0
        TSCLAD=TCALO
        DO
          ITER=ITER+1
          IF(ITER.GT.50) CALL XABORT('THMLEAD: HCONV FAILURE.')
          IF(TSCLAD.GT.TSAT) THEN
            CALL XABORT('THMLEAD: INVALID HEAT TRANSFER REGIME')
          ENDIF
          HA=0.023*(ZKONE/HD)*REL**0.8*PRL**0.4
          F=1.0
          S=1.0
          HB=0.0
          K0=0
          HCONV=F*HA+S*HB
          IF(HCONV.LE.0.0) THEN
            WRITE(HSMG,'(35HTHMLEAD: DRY OUT REACHED IN CHANNEL,3I5)')
     >      I,J,K
            CALL XABORT(HSMG)
          ENDIF
          IF(ITIME.EQ.0) THEN
            TWAL=(PHI+S*HB*TSAT+F*HA*TCALO)/(S*HB+F*HA)
          ELSE
            TWAL=0.0
            CALL XABORT('THMLEAD: TRANSIENT NOT IMPLEMENTED YET')
          ENDIF
          IF(ABS(TSCLAD-TWAL).GT.1.0E-5*TSCLAD) THEN
            TSCLAD=TWAL
          ELSE
            EXIT
          ENDIF
        ENDDO
      ELSE IF(IHCONV.EQ.1) THEN
        IF(ITIME.EQ.0) THEN
          TSCLAD=TCALO+PHI/KHCONV
        ELSE
          RCHC=RADCL*KHCONV
          TSCLAD=MAX(273.15,(ZF(1)+RCHC*TCALO)/(ZF(2)+RCHC))
          PHI=(ZF(1)-TSCLAD*ZF(2))/RADCL
        ENDIF
      ENDIF
      RETURN
      END
