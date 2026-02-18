*DECK LIBPRQ
      SUBROUTINE LIBPRQ(MAXTRA,DELI,AWR,E0,Q,IALTER,IL,N,PRI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the PRI array for various Legendre orders using Gaussian
* integration. Inelastic scattering case.
*
*Copyright:
* Copyright (C) 2025 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* MAXTRA  allocated dimension of array PRI.
* DELI    elementary lethargy width of the equi-width lethargy mesh.
* AWR     mass ratio for current isotope.
* E0      energy corresponding to the upper limit of primary group.
* Q       Q-value (negative value) for an inelastic diffusion.
* IALTER  type of approximation (=0: use exponentials; =1: use Taylor
*         expansions).
* IL      Legendre order (=0: isotropic kernel in LAB).
*
*Parameters: output
* N       exact dimension of array PRI.
* PRI     array containing the slowing-down probabilities defined on
*         an equi-width lethargy mesh.
*
*Reference:
* M. Grandotto-Biettoli, "AUTOSECOL, un calcul automatique de
* l'auto-protection des resonances des isotopes lourds,"
* Note CEA-N-1961, Commissariat a l'Energie Atomique, 1977. 
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER MAXTRA,IALTER,IL,N
      REAL DELI,AWR,E0,Q,PRI(MAXTRA)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NGPT=6,MAXNL=50)
      DOUBLE PRECISION AWRB,ALP,T0,FACT,COEF0,GAM,UM,UMIN,UMAX
      REAL UI(NGPT),WI(NGPT),UJ(NGPT),WJ(NGPT)
      REAL POLY(0:MAXNL),CALC(0:MAXNL,0:2)
      CHARACTER HSMG*131
      ZMU(AWR,RGAM,U)=0.5*(AWR+1.0)*EXP(-0.5*U)-0.5*(RGAM/(AWR+1.0)+
     1 AWR-1.0)*EXP(0.5*U)
*----
*  COMPUTE THE LEGENDRE POLYNOMIAL OF ORDER IL.
*----
      IF(IL.GT.MAXNL) CALL XABORT('LIBPRQ: IL OVERFLOW.')
      IF(IL.EQ.0) THEN
        POLY(0)=1.0
      ELSE IF(IL.EQ.1) THEN
        POLY(0)=0.0
        POLY(1)=1.0
      ELSE
        CALC(0:MAXNL,0:1)=0.0
        CALC(0,0)=1.0
        CALC(1,1)=1.0
        DO J=2,IL
          DO I=0,IL
            T0=-REAL(J-1)*CALC(I,MOD(J-2,3))
            IF(I.GT.0) T0=T0+(2.0*REAL(J-1)+1.0)*CALC(I-1,MOD(J-1,3))
            CALC(I,MOD(J,3))=REAL(T0)/REAL(J)
          ENDDO
        ENDDO
        POLY(0:IL)=CALC(0:IL,MOD(IL,3))
      ENDIF
*
      AWRB=AWR
      IF(AWR.LT.1.0001) AWRB=1.0001
      GAM=AWRB*(AWRB+1.D0)*Q/E0
      UMIN=LOG((AWRB+1.D0)**2/(AWRB**2+1.D0+GAM+2.D0*SQRT(AWRB**2+GAM)))
      IF(-GAM.GT.AWRB**2) CALL XABORT('LIBPRQ: NEGATIVE SQRT ARGUMENT.')
      IF(UMIN.LT.DELI) THEN
        ! pathological case where the treshold is negligible
        CALL LIBPRI(MAXTRA,DELI,AWR,IALTER,IL,N,PRI)
        RETURN
      ENDIF
      PRI(:MAXTRA)=0.0
      ALP=((AWRB-1.D0)/(AWRB+1.D0))**2
      CALL ALGPT(NGPT,0.0,DELI,UI,WI)
      N=0
      DO I=1,NGPT
        ! primary group base point
        WII=WI(I)
        UII=UI(I)
        EN=E0*EXP(-UI(I)*DELI)
        GAM=AWRB*(AWRB+1.D0)*Q/EN
        UM=AWRB**2+1.D0+GAM
        UMIN=LOG((AWRB+1.D0)**2/(UM+2.D0*SQRT(AWRB**2+GAM)))
        UMAX=LOG((AWRB+1.D0)**2/(UM-2.D0*SQRT(AWRB**2+GAM)))
        NMIN=1+INT((UMIN-1.D-6)/DELI)
        NMAX=1+INT((UMAX-1.D-6)/DELI)
        IF(NMAX.GT.MAXTRA) THEN
          WRITE(HSMG,'(25HLIBPRQ: MAXTRA OVERFLOW (,I8,2H >,I8,2H).)')
     1    NMAX,MAXTRA
          CALL XABORT(HSMG)
        ENDIF
        COEF0=AWRB/SQRT(AWRB**2+GAM)
        WII=WII*REAL(COEF0)
        IF(NMAX.EQ.NMIN) THEN
          CALL ALGPT(NGPT,REAL(UMIN),REAL(UMAX),UJ,WJ)
          DO J=1,NGPT
            FACT=POLY(0)
            T0=1.0D0
            DO K=1,IL
              T0=T0*ZMU(AWR,REAL(GAM),UJ(J)-UII)
              FACT=FACT+POLY(K)*T0
            ENDDO
            IF(IALTER.EQ.0) THEN
              PRI(NMIN)=PRI(NMIN)+WII*WJ(J)*EXP(UII-UJ(J))*REAL(FACT)
            ELSE
              PRI(NMIN)=PRI(NMIN)+WII*WJ(J)*REAL(FACT/(UMAX-UMIN))
            ENDIF
          ENDDO ! J
        ELSE
          CALL ALGPT(NGPT,REAL(UMIN),NMIN*DELI,UJ,WJ)
          DO J=1,NGPT
            FACT=POLY(0)
            T0=1.0D0
            DO K=1,IL
              T0=T0*ZMU(AWR,REAL(GAM),UJ(J)-UII)
              FACT=FACT+POLY(K)*T0
            ENDDO
            IF(IALTER.EQ.0) THEN
              PRI(NMIN)=PRI(NMIN)+WII*WJ(J)*EXP(UII-UJ(J))*REAL(FACT)
            ELSE
              PRI(NMIN)=PRI(NMIN)+WII*WJ(J)*REAL(FACT/(UMAX-UMIN))
            ENDIF
          ENDDO ! J
          DO N=NMIN+1,NMAX-1
            CALL ALGPT(NGPT,(N-1)*DELI,N*DELI,UJ,WJ)
            DO J=1,NGPT
              FACT=POLY(0)
              T0=1.0D0
              DO K=1,IL
                T0=T0*ZMU(AWR,REAL(GAM),UJ(J)-UII)
                FACT=FACT+POLY(K)*T0
              ENDDO
              IF(IALTER.EQ.0) THEN
                PRI(N)=PRI(N)+WII*WJ(J)*EXP(UII-UJ(J))*REAL(FACT)
              ELSE
                PRI(N)=PRI(N)+WII*WJ(J)*REAL(FACT/(UMAX-UMIN))
              ENDIF
            ENDDO ! J
          ENDDO ! N
          CALL ALGPT(NGPT,(NMAX-1)*DELI,REAL(UMAX),UJ,WJ)
          DO J=1,NGPT
            FACT=POLY(0)
            T0=1.0D0
            DO K=1,IL
              T0=T0*ZMU(AWR,REAL(GAM),UJ(J)-UII)
              FACT=FACT+POLY(K)*T0
            ENDDO
            IF(IALTER.EQ.0) THEN
              PRI(NMAX)=PRI(NMAX)+WII*WJ(J)*EXP(UII-UJ(J))*REAL(FACT)
            ELSE
              PRI(NMAX)=PRI(NMAX)+WII*WJ(J)*REAL(FACT/(UMAX-UMIN))
            ENDIF
          ENDDO ! J
        ENDIF
        N=MAX(N,NMAX)
      ENDDO ! I
      IF(IALTER.EQ.0) THEN
        PRI(:N)=PRI(:N)/DELI/REAL(1.0D0-ALP)
      ELSE
        PRI(:N)=PRI(:N)/DELI
      ENDIF
*----
*  NORMALIZATION
*----
      IF(IL.EQ.0) THEN
        FACT=0.0D0
        DO I=1,N
          FACT=FACT+PRI(I)
        ENDDO
        PRI(:N)=PRI(:N)/REAL(FACT)
      ENDIF
      RETURN
      END
