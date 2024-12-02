*DECK SNQU07
      SUBROUTINE SNQU07(NLF,X,W)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To define Gauss-Lobatto points and weights (1D quadrature).
*
*Copyright:
* Copyright (C) 2008 Ecole Polytechnique de Montreal.
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): 
* C. Bienvenue
*
*Parameters: input
* NLF  order of the SN approximation (even number).
*
*Parameters: output
* X     base points in $\\mu$ of the quadrature.
* W     weights of the quadrature.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NLF
      REAL X(NLF),W(NLF)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(PI=3.141592654,NMAX=100,EPS=1.0E-6)
      INTEGER N
      REAL XOLD(NLF),P(NLF,NLF)
*----
*  COMPUTE QUADRATURE PARAMETERS
*----

      ! INITIAL GUESS BASE POINTS
      DO 10 I=1,NLF
         X(I)=COS(PI*(I-1)/(NLF-1))
         XOLD(I)=2.0
   10 ENDDO

      ! NEWTON-RAPHSON METHOD TO COMPUTE BASE POINTS
      P(:NLF,:NLF)=0.0
      DO 20 I=1,NLF
        N=1
        DO WHILE (ABS(X(I)-XOLD(I)).GT.EPS)
           XOLD(I) = X(I)
           P(I,1)=1
           P(I,2)=X(I)
           DO 30 J=2,NLF-1
              P(I,J+1)=((2*J-1)*X(I)*P(I,J)-(J-1)*P(I,J-1))/J 
   30      ENDDO
           X(I)=XOLD(I)-(X(I)*P(I,NLF)-P(I,NLF-1))/(NLF*P(I,NLF))
           IF(N.GT.NMAX) CALL XABORT('SNQU07: CONVERGENCE ISSUE.')
           N=N+1
        ENDDO
        W(I)=2/((NLF-1)*NLF*P(I,NLF)**2)
   20 ENDDO
      
      ! COMPUTE WEIGHTS
      DO 40 I=1,NLF
         W(I)=2/((NLF-1)*NLF*P(I,NLF)**2)
   40 ENDDO

      RETURN
      END