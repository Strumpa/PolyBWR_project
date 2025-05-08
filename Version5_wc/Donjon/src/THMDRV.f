*DECK THMDRV
      SUBROUTINE THMDRV(MPTHM,IMPX,IX,IY,NZ,XBURN,VOLXY,HZ,CFLUX,POROS,
     > FNFU,NFD,NDTOT,IFLUID,SNAME,SCOMP,IGAP,IFUEL,FNAME,FCOMP,FCOOL,
     > FFUEL,ACOOL,HD,PCH,RAD,MAXIT1,MAXITL,ERMAXT,SPEED,TINLET,PINLET,
     > FRACPU,ICONDF,NCONDF,KCONDF,UCONDF,ICONDC,NCONDC,KCONDC,UCONDC,
     > IHGAP,KHGAP,IHCONV,KHCONV,WTEFF,IFRCDI,ISUBM,FRO,POW,IPRES,TCOMB,
     > DCOOL,TCOOL,TSURF,HCOOL,PCOOL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Driver of the steady-state thermal-hydraulics calculation.
*
*Copyright:
* Copyright (C) 2012 Ecole Polytechnique de Montreal.
*
*Author(s): 
* A. Hebert
* C. Garrido 
* 08/2023: Modifications to include Molten Salt heat transfer in coolant
* C. Garrido 
* 07/2024: Modifications to include Molten Salt heat transfer in static
*          fuel
* C. Huet
* 02/2025: Modifications to include pressure drop calculation
*
*Parameters: input
* MPTHM   directory of the THM object containing steady-state
*         thermohydraulics data.
* IMPX    printing index (=0 for no print).
* IX      position of mesh along X direction.
* IY      position of mesh along Y direction.
* NZ      number of meshes along Z direction (channel direction).
* XBURN   burnup distribution in MWday/tonne.
* VOLXY   mesh area in the radial plane.
* HZ      Z-directed mesh widths.
* CFLUX   critical heat flux in W/m^2.
* POROS   oxyde porosity.
* FNFU    number of active fuel rods in the fuel bundle.
* NFD     number of discretization points in fuel region.
* NDTOT   number of total discretization points in the the fuel
*         pellet and the cladding.
* IFLUID  type of fluid (0=H2O; 1=D2O; 2=SALT).
* SNAME   Name of the molten salt (e.g. "LiF-BeF2")
* SCOMP   Composition of the molten salt (e.g. "0.66-0.34")
* FCOOL   power density fraction in coolant.
* FFUEL   power density fraction in fuel.
* ACOOL   coolant cross section area in m^2.
* HD      hydraulic diameter of one assembly in m.
* PCH     heating perimeter in m.
* RAD     fuel and clad radii in m.
* MAXIT1  maximum number of conduction iterations.
* MAXITL  maximum number of center-pellet iterations.
* ERMAXT  convergence criterion.
* SPEED   inlet flow velocity in m/s.
* TINLET  inlet temperature in K.
* PINLET  inlet pressure in Pa.
* FRACPU  plutonium fraction in fuel.
* ICONDF  fuel conductivity flag (0=Stora-Chenebault or COMETHE/
*         1=user-provided polynomial + inverse term).
* NCONDF  degree of user-provided fuel conductivity polynomial.
* KCONDF  polynomial coefficients for fuel conductivity in W/m/K^(k+1)
*         (except for the two last coefficients which belongs to the
*         inverse term).
* UCONDF  required unit of temperature in polynomial for fuel
*         conductivity (KELVIN or CELSIUS).
* ICONDC  clad conductivity flag (0=default/1=user-provided
*         polynomial).
* NCONDC  degree of user-provided clad conductivity polynomial.
* KCONDC  polynomial coefficients for clad conductivity in W/m/K^(k+1).
* UCONDC  required unit of temperature in polynomial for clad
*         conductivity (KELVIN or CELSIUS).
* IHGAP   flag indicating HGAP chosen (0=default/1=user-provided).
* KHGAP   fixed user-provided HGAP value in W/m^2/K.
* IHCONV  flag indicating HCONV chosen (0=default/1=user-provided).
* KHCONV  fixed user-provided HCONV value in W/m^2/K.
* WTEFF   surface temperature's weighting factor in effective fuel
*         temperature.
* IFRCDI  flag indicating if average approximation is forced during
*         fuel conductivity evaluation (0=default/1=average
*         approximation forced).
* ISUBM   subcooling model (0: one-phase; 1: Jens-Lottes model;
*         2: Saha- Zuber model).
* FRO     radial power form factors.
* POW     power distribution in W.
* IGAP    Flag indicating if the gap is considered (0=gap/1=no gap)
* IFUEL   type of fuel (0=UO2/MOX; 1=SALT).
* FNAME   Name of the molten salt (e.g. "LiF-BeF2")
* FCOMP   Composition of the molten salt (e.g. "0.66-0.34")
* IPRES   flag indicating if pressure is to be computed (0=nonstant/
*         1=variable).
*
*Parameters: output
* TCOMB   averaged fuel temperature distribution in K.
* DCOOL   coolant density distribution in g/cc.
* TCOOL   coolant temperature distribution in K.
* TSURF   surface fuel temperature distribution in K.
* HCOOL   coolant enthalpty distribution in J/kg.
* PCOOL   coolant pressure distribution in Pa.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE t_saltdata
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) MPTHM
      INTEGER IMPX,IX,IY,NZ,NFD,NDTOT,IFLUID,MAXIT1,MAXITL,IHGAP,IGAP,
     > IFUEL,IPRES
      REAL XBURN(NZ),VOLXY,HZ(NZ),CFLUX,POROS,FNFU,FCOOL,FFUEL,ACOOL,
     > HD,PCH,RAD(NDTOT-1,NZ),ERMAXT,SPEED,TINLET,PINLET,FRACPU,
     > KCONDF(NCONDF+3),KCONDC(NCONDC+1),KHGAP,KHCONV,WTEFF,FRO(NFD-1),
     > POW(NZ),TCOMB(NZ),DCOOL(NZ),TCOOL(NZ),TSURF(NZ),HCOOL(NZ),
     > PCOOL(NZ),MUT(NZ)
      CHARACTER UCONDF*12,UCONDC*12
*----
*  LOCAL VARIABLES
*----
      TYPE(tpdata) STP,FTP
      PARAMETER (KMAXO=100,MAXNPO=40)
      REAL TRE11(MAXNPO),RADD(MAXNPO),ENT(4),MFLOW,TLC(NZ)
      CHARACTER HSMG*131,SNAME*32,SCOMP*32,FNAME*32,FCOMP*32
      REAL XS(4),TC1,PC(NZ),TP(NZ),RHOL,XFL(NZ),EPS(NZ),HINLET,
     > TCLAD(NZ),ENTH(NZ),SLIP(NZ),AGM(NZ),QFUEL(NZ),QCOOL(NZ),K11
      INTEGER KWA(NZ)
      REAL XX2(MAXNPO),XX3(MAXNPO),ZF(2)
      DATA XS/-0.861136,-0.339981,0.339981,0.861136/

      REAL TBUL(NZ)

      INTEGER I
      REAL ERRV, ERRP
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: VCOOL,TCENT,DLCOOL
      REAL, ALLOCATABLE, DIMENSION(:,:) :: TEMPT

      REAL, ALLOCATABLE, DIMENSION(:) :: PTEMP, VTEMP
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(VCOOL(NZ),TEMPT(NDTOT,NZ),TCENT(NZ),DLCOOL(NZ))
      ALLOCATE(PTEMP(NZ), VTEMP(NZ))
*----
*  COMPUTE THE INLET FLOW ENTHALPY AND VELOCITY
*----
      IF(NDTOT.GT.MAXNPO) CALL XABORT('THMDRV: MAXNPO OVERFLOW.')
      IF(IFLUID.EQ.0) THEN
         CALL THMSAT(PINLET,TSAT)
      ELSE IF(IFLUID.EQ.1) THEN
         CALL THMHST(PINLET,TSAT)
*CGT TODO: GET ALSO FREEZING??
      ELSE IF(IFLUID.EQ.2) THEN
         CALL THMSGT(SNAME,SCOMP,STP,IMPX)
         CALL THMSST(STP,TSAT,IMPX)
*CGT
      ENDIF
      IF (IFUEL.EQ.1) THEN
         CALL THMSGT(FNAME,FCOMP,FTP,IMPX)
      ENDIF

      IF(TINLET.GT.TSAT) THEN
         WRITE(HSMG,'(27HTHMDRV: INLET TEMPERATURE (,1P,E12.4,
     >   40H K) GREATER THAN SATURATION TEMPERATURE.)') TINLET
         CALL XABORT(HSMG)
      ENDIF
      IF(IFLUID.EQ.0) THEN
        CALL THMPT(PINLET,TINLET,RHOIN,HINLET,R3,R4,R5)
      ELSE IF(IFLUID.EQ.1) THEN
        CALL THMHPT(PINLET,TINLET,RHOIN,HINLET,R3,R4,R5)
      ELSE IF(IFLUID.EQ.2) THEN
        CALL THMSPT(STP,TINLET,RHOIN,HINLET,R3,R4,R5,IMPX)
      ENDIF
      MFLOW=SPEED*RHOIN
      HMSUP=HINLET
*----
*  INITIALIZE VALUES OF STEAM QUALITIES, VOID FRACTION AND DENSITY
*  PRESSURE, VELOCITY AND TEMPERATURE OF THE COOLANT ALONG THE CHANNEL.
*---
      DO K=1,NZ
         EPS(K)=0.0
         XFL(K)=0.0
         SLIP(K)=1.0
         KWA(K)=0
         MUT(K)=0.0
         QFUEL(K)=0.0

         PCOOL(K)=PINLET
         VCOOL(K)=MFLOW/RHOIN
         DCOOL(K)=RHOIN
         DLCOOL(K)=RHOIN
         TCOOL(K)=TINLET
*----
*  COMPUTE THE SATURATION TEMPERATURE AND THE THERMODYNAMIC PROPERTIES
*  IF THE PRESSURE DROP IS COMPUTED
*---

        IF (IPRES.EQ.1) THEN
          IF(POW(K).EQ.0.0) CYCLE
          IF(IFLUID.EQ.0) THEN
            CALL THMSAT(PCOOL(K),TSAT)
          ELSE IF(IFLUID.EQ.1) THEN
            CALL THMHST(PCOOL(K),TSAT)
          ENDIF

          TB=TSAT-0.1
          IF(TCOOL(K).LT.TB) THEN
            IF(IFLUID.EQ.0) THEN
              CALL THMPT(PCOOL(K),TCOOL(K),RHOIN,H11,K11,MUT(K),C11)
            ELSE IF(IFLUID.EQ.1) THEN
              CALL THMHPT(PCOOL(K),TCOOL(K),RHOIN,H11,K11,MUT(K),C11)
            ELSE IF(IFLUID.EQ.2) THEN
              CALL THMSPT(STP,TCOOL(K),R11,H11,K11,MUT(K),C11,IMPX)
            ENDIF
          ELSE
            IF(IFLUID.EQ.0) THEN
              CALL THMPT(PCOOL(K),TB,R11,H11,K11,MUT(K),C11)
            ELSE IF(IFLUID.EQ.1) THEN
              CALL THMHPT(PCOOL(K),TB,R11,H11,K11,MUT(K),C11)
            ELSE IF(IFLUID.EQ.2) THEN
              CALL THMSPT(STP,TB,R11,H11,K11,MUT(K),C11,IMPX)
            ENDIF
          ENDIF
        ENDIF
      ENDDO
*----
*  MAIN LOOP ALONG THE 1D CHANNEL.
*----
      ERRV = 1.0
      ERRP = 1.0
      I=0
      IF (IPRES .EQ. 0) GOTO 30
   10 CONTINUE
*----
*  UPDATE HINLET FUNCTION OF INLET PRESSURE AND TEMPERATURE
*----
       HMSUP=HINLET
       SPEED=MFLOW/DCOOL(1)
*----
*  WHILE LOOP FOR PRESSURE AND VELOCITY CONVERGENCE
*  CHECK FOR CONVERGENCE
*----
        IF (I .GT. 1000) GOTO 20
        IF ((ERRP < 1E-3) .AND. (ERRV < 1E-3)) GOTO 20

          I = I + 1

          PTEMP = PCOOL
          VTEMP = VCOOL
          CALL THMPV(SPEED, PINLET, VCOOL, DCOOL, 
     >              PCOOL, MUT, XFL, HD, NZ, HZ)
   30 CONTINUE
*----
*  MAIN LOOP ALONG THE 1D CHANNEL.
*----
      K0=0 ! onset of nuclear boiling point
      DO K=1,NZ
        IF(POW(K).EQ.0.0) CYCLE
        IF(IFLUID.EQ.0) THEN
          CALL THMSAT(PCOOL(K),TSAT)
        ELSE IF(IFLUID.EQ.1) THEN
          CALL THMHST(PCOOL(K),TSAT)
        ENDIF
        TBUL(K)=TSAT
*----
*  COMPUTE THE LINEAR POWER, THE VOLUMIC POWER AND THE THERMAL EXCHANGE
*  COEFFICIENT OF THE GAP
*----
        DV=VOLXY*HZ(K)
*       linear power in W/m
        POWLIN=(POW(K)/DV)*VOLXY/FNFU
*       volumic power in W/m^3
        QFUEL(K)=POW(K)*FFUEL/DV
        QCOOL(K)=POW(K)*FCOOL/DV
*----
*  INITIALIZATION OF PINCELL TEMPERATURES
*----
        IF(POW(K).EQ.0.0) CYCLE
        IF(IMPX.GT.4) WRITE(6,190) K
        DO L=1,NDTOT
           TRE11(L)=TCOMB(K)
        ENDDO
        DO L=1,NDTOT-1
           RADD(L)=RAD(L,K)
        ENDDO
*----
*  COMPUTE THE POWER DENSITY AND HEAT FLOW ALONG THE CHANNEL
*----
*       out-of-clad heat flow in W/m2
        IF(IMPX.GT.5) WRITE(6,'(15H THMDRV: QFUEL(,I5,2H)=,1P,E12.4,
     >  6H W/m2.)') K,QFUEL(K)
        PHI2=0.5*QFUEL(K)*RAD(NFD,K)**2/RAD(NDTOT-1,K)
        IF(PHI2.GT.CFLUX) THEN
          WRITE(HSMG,'(23HTHMDRV: THE HEAT FLUX (,1P,E12.4,5H) IS ,
     >    37HGREATER THAN THE CRITICAL HEAT FLUX (,E12.4,2H).)')
     >    PHI2,CFLUX
          CALL XABORT(HSMG)
        ENDIF
*----
*  COMPUTE FOUR VALUES OF ENTHALPY IN J/KG TO BE USED IN GAUSSIAN
*  INTEGRATION. DELTH1 IS THE ENTHALPY INCREASE IN EACH AXIAL MESH.
*----
        DELTH1=(PCH/ACOOL*PHI2+QCOOL(K))*HZ(K)/MFLOW
        DO I1=1,4
          POINT=(1.0+XS(I1))/2.0
          ENT(I1)=HMSUP+POINT*DELTH1
        ENDDO
        HMSUP=HMSUP+DELTH1
*----
*  COMPUTE THE VALUE OF THE DENSITY AND THE CLAD-COOLANT HEAT TRANSFER
*  COEFFICIENT
*----
        IF(K.GT.1) THEN
          XFL(K)=XFL(K-1)
          EPS(K)=EPS(K-1)
          SLIP(K)=SLIP(K-1)
        ENDIF
*CGT
        IF ((IFLUID.EQ.0).OR.(IFLUID.EQ.1)) THEN
          CALL THMH2O(0,IX,IY,K,K0,PCOOL(K),MFLOW,HMSUP,ENT,HD,IFLUID,
     >    IHCONV,KHCONV,ISUBM,RAD(NDTOT-1,K),ZF,PHI2,XFL(K),EPS(K),
     >    SLIP(K),ACOOL,PCH,HZ(K),TCALO,RHO,RHOL,TRE11(NDTOT),KWA(K))
        ELSEIF (IFLUID.EQ.2) THEN
          CALL THMSAL(IMPX,0,IX,IY,K,K0,MFLOW,HMSUP,ENT,HD,STP,
     >    IHCONV,KHCONV,ISUBM,RAD(NDTOT-1,K),ZF,PHI2,XFL(K),
     >    EPS(K),SLIP(K),HZ(K),TCALO,RHO,RHOL,TRE11(NDTOT),
     >    KWA(K))
        ENDIF
*CGT
*----
*  STEADY-STATE SOLUTION OF THE CONDUCTION EQUATIONS IN A FUEL PIN.
*----
        DTINV=0.0
        IF(IGAP.EQ.0) THEN
          CALL THMROD(IMPX,NFD,NDTOT-1,MAXIT1,MAXITL,ERMAXT,DTINV,RADD,
     >    TRE11,TRE11,QFUEL(K),FRO,TRE11(NDTOT),POWLIN,XBURN(K),
     >    POROS,FRACPU,ICONDF,NCONDF,KCONDF,UCONDF,ICONDC,NCONDC,
     >    KCONDC,UCONDC,IHGAP,KHGAP,IFRCDI,TC1,XX2,XX3,ZF)
        ELSE
          CALL THMRNG(IMPX,NFD,NDTOT-1,MAXIT1,MAXITL,ERMAXT,DTINV,RADD,
     >    TRE11,TRE11,QFUEL(K),FRO,TRE11(NDTOT),XBURN(K),
     >    POROS,FRACPU,ICONDF,NCONDF,KCONDF,UCONDF,ICONDC,NCONDC,
     >    KCONDC,UCONDC,IFRCDI,IFUEL,FTP,TC1,XX2,XX3,ZF)
        ENDIF
*
        DO K1=1,NDTOT-1
          TRE11(K1)=XX2(K1)+TRE11(NDTOT)*XX3(K1)
        ENDDO
*----
*  RECOVER MESHWISE TEMPERATURES AND FLUID DENSITY. BY DEFAULT, USE THE
*  ROWLANDS FORMULA TO COMPUTE THE EFFECTIVE FUEL TEMPERATURE, OTHERWISE
*  USE USER-SPECIFIED WEIGHTING FACTOR.
*----  
        TCOMB(K)=(1.0-WTEFF)*TC1+WTEFF*TRE11(NFD)
        TCENT(K)=TC1
        TSURF(K)=TRE11(NFD)
        TCLAD(K)=TRE11(NDTOT)
        TCOOL(K)=TCALO
        DCOOL(K)=RHO
        DLCOOL(K)=RHOL
        HCOOL(K)=HMSUP
        !PCOOL(K)=PINLET
        PC(K)=PINLET
        TP(K)=TCLAD(K)
        TLC(K)=TCOOL(K)
        ENTH(K)=HCOOL(K)
        AGM(K)=MFLOW ! constant flow rate
        DO K2=1,NDTOT
          TEMPT(K2,K)=TRE11(K2)
        ENDDO
        IF (IPRES .EQ. 0) THEN
          PCOOL(K)=PINLET
          VCOOL(K)=MFLOW/DCOOL(K)
        ENDIF
*----
*  COMPUTE THE SATURATION TEMPERATURE AND THE THERMODYNAMIC PROPERTIES
*  IF THE PRESSURE DROP IS COMPUTED
*---
        IF (IPRES.EQ.1) THEN
          IF(POW(K).EQ.0.0) CYCLE
          IF(IFLUID.EQ.0) THEN
            CALL THMSAT(PCOOL(K),TSAT)
          ELSE IF(IFLUID.EQ.1) THEN
            CALL THMHST(PCOOL(K),TSAT)
          ENDIF
 
          TB=TSAT-0.1
          IF(TCOOL(K).LT.TB) THEN
            IF(IFLUID.EQ.0) THEN
             CALL THMPT(PCOOL(K),TCOOL(K),RHOIN,H11,K11,MUT(K),C11)
            ELSE IF(IFLUID.EQ.1) THEN
             CALL THMHPT(PCOOL(K),TCOOL(K),RHOIN,H11,K11,MUT(K),C11)
            ELSE IF(IFLUID.EQ.2) THEN
             CALL THMSPT(STP,TCOOL(K),R11,H11,K11,MUT(K),C11,IMPX)
            ENDIF
          ELSE
            IF(IFLUID.EQ.0) THEN
             CALL THMPT(PCOOL(K),TB,R11,H11,K11,MUT(K),C11)
            ELSE IF(IFLUID.EQ.1) THEN
             CALL THMHPT(PCOOL(K),TB,R11,H11,K11,MUT(K),C11)
            ELSE IF(IFLUID.EQ.2) THEN
             CALL THMSPT(STP,TB,R11,H11,K11,MUT(K),C11,IMPX)
            ENDIF
          ENDIF
         ENDIF
      ENDDO
*----
* IF THE PRESSURE DROP IS COMPUTED, COMPUTE THE 
* THE PRESSURE AND VELOCITY RESIDUALS
*----
      IF (IPRES .EQ. 0) GOTO 20
      ERRV = 0
      ERRP = 0
      DO K=1,NZ
        ERRV = ERRV + ABS(VCOOL(K) - VTEMP(K))/VTEMP(K)
        ERRP = ERRP + ABS(PCOOL(K) - PTEMP(K))/PTEMP(K)
      ENDDO
      ERRV = ERRV/NZ
      ERRP = ERRP/NZ
      GO TO 10

   20 CONTINUE

      IF (I.GE.1000) THEN
        CALL XABORT('THMDRV: MAXIMUM NB OF ITERATIONS REACHED.')
      ELSE IF(IMPX.GT.0) THEN
        WRITE(6,'(37H THMDRV: CONVERGENCE REACHED AT ITER=,I5,1H.)') I
      ENDIF
*----
* PRINT THE THERMOHYDRAULICAL PARAMETERS
*----
      IF(IMPX.GT.4) THEN
        WRITE(6,250) 'POW', POW(:NZ)
        WRITE(6,250) 'PCOOL', PCOOL(:NZ)
        WRITE(6,250) 'VCOOL', VCOOL(:NZ)
        WRITE(6,250) 'DCOOL', DCOOL(:NZ)
        WRITE(6,250) 'TCOOL', TCOOL(:NZ)
        WRITE(6,250) 'EPS', EPS(:NZ)
        WRITE(6,250) 'XFL', XFL(:NZ)
        WRITE(6,250) 'TSAT', TBUL(:NZ)
        WRITE(6,250) 'MUT', MUT(:NZ)
      ENDIF
*----
*  PRINT THE OUTLET THERMOHYDRAULICAL PARAMETERS
*----
      IF(IMPX.GT.3) THEN
        WRITE(6,'(/16H THMDRV: CHANNEL,2I6/1X,27(1H-))') IX,IY
        WRITE(6,210) ' ____________________________________________',
     >          '_____________________________________________________',
     >          '_____________________________________________________',
     >          '______________'
        WRITE(6,210) '|     |   TCOMB    |   TSURF    |    DCOOL ',
     >          '   |    TCOOL    |    PCOOL    |    HCOOL    |    ',
     >          'QFUEL    |    QCOOL    |    VOID   |     QUAL    |',
     >          '     SLIP    |  FLOW  |',
     >          '|     |     K      |     K      |    Kg/m3    |   ',
     >          '   K      |     Pa      |    J/Kg     |    W/m3   ',
     >          '  |    W/m3     |           |             |       ',
     >          '      | REGIME |'
        WRITE(6,210) '|_____|____________|____________|____________',
     >          '_|_____________|_____________|_____________|_________',
     >          '____|_____________|___________|_____________|________',
     >          '_____|________|'
        DO L=NZ,1,-1
          IF(L.EQ.1) THEN
            WRITE(6,220) '| BOT |',TCOMB(L),' |',TSURF(L),
     >            ' |',DCOOL(L),' |',TCOOL(L),' |',PCOOL(L),
     >            ' |',HCOOL(L),' |',QFUEL(L),' |',QCOOL(L),' |',
     >            EPS(L),' |',XFL(L),' |',SLIP(L),' |',KWA(L),' |'
          ELSEIF(L.EQ.NZ) THEN
            WRITE(6,220) '| TOP |',TCOMB(L),' |',TSURF(L),
     >            ' |',DCOOL(L),' |',TCOOL(L),' |',PCOOL(L),
     >            ' |',HCOOL(L),' |',QFUEL(L),' |',QCOOL(L),' |',
     >            EPS(L),' |',XFL(L),' |',SLIP(L),' |',KWA(L),' |'
          ELSE
            WRITE(6,230) '| ',L,' |',TCOMB(L),' |',TSURF(L),
     >            ' |',DCOOL(L),' |',TCOOL(L),' |',PCOOL(L),
     >            ' |',HCOOL(L),' |',QFUEL(L),' |',QCOOL(L),' |',
     >            EPS(L),' |',XFL(L),' |',SLIP(L),' |',KWA(L),' |'
          ENDIF
        ENDDO
        WRITE(6,210) '|_____|____________|____________|____________',
     >          '_|_____________|_____________|_____________|_________',
     >          '____|_____________|___________|_____________|________',
     >          '_____|________|'
        WRITE(6,240) MFLOW
      ENDIF
*----
*  MODIFICATION OF THE VECTORS TO FIT THE GEOMETRY OF THE CHANNELS AND
*  THE BUNDLES AND WRITE THE DATA IN LCM OBJECT THM 
*----
      CALL LCMPUT(MPTHM,'PRESSURE',NZ,2,PCOOL)
      CALL LCMPUT(MPTHM,'DENSITY',NZ,2,DCOOL)
      CALL LCMPUT(MPTHM,'LIQUID-DENS',NZ,2,DLCOOL)
      CALL LCMPUT(MPTHM,'ENTHALPY',NZ,2,HCOOL)
      CALL LCMPUT(MPTHM,'VELOCITIES',NZ,2,VCOOL)
      CALL LCMPUT(MPTHM,'CENTER-TEMPS',NZ,2,TCENT)
      CALL LCMPUT(MPTHM,'COOLANT-TEMP',NZ,2,TCOOL)
      CALL LCMPUT(MPTHM,'POWER',NZ,2,POW)
      CALL LCMPUT(MPTHM,'TEMPERATURES',NDTOT*NZ,2,TEMPT)
      CALL LCMPUT(MPTHM,'PINLET',1,2,PINLET)
      CALL LCMPUT(MPTHM,'TINLET',1,2,TINLET)
      CALL LCMPUT(MPTHM,'VINLET',1,2,SPEED)
      CALL LCMPUT(MPTHM,'POULET',1,2,PINLET)
      CALL LCMPUT(MPTHM,'RADII',(NDTOT-1)*NZ,2,RAD)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(DLCOOL,TCENT,TEMPT,VCOOL)
      DEALLOCATE(PTEMP, VTEMP)
      RETURN
*
  190 FORMAT(/21H THMDRV: AXIAL SLICE=,I5)
  210 FORMAT(1X,A,A,A,A)
  220 FORMAT(1X,A,F11.2,A,F11.2,A,F12.4,A,F12.2,A,3P,E12.4,
     >       A,1P,E12.4,A,1P,E12.4,A,1P,E12.4,A,0P,F10.4,A,E12.4,A,
     >       E12.4,A,I5,2X,A)
  230 FORMAT(1X,A,I3,A,F11.2,A,F11.2,A,F12.4,A,F12.2,A,3P,E12.4,
     >       A,1P,E12.4,A,1P,E12.4,A,1P,E12.4,A,0P,F10.4,A,E12.4,A,
     >       E12.4,A,I5,2X,A)
  240 FORMAT(7H MFLXT=,1P,E12.4,8H Kg/m2/s)
  250 FORMAT(9H THMDRV: ,A6,1H:,1P,11E12.4/(4X,12E12.4))
      END
