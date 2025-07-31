*DECK XDREXP
      SUBROUTINE XDREXP(DX,NBX,PARAM,E00,E01,E10,E11)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Construct exponential tables for linear interpolation.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy, R. Le Tellier 
*
*Parameters: input
* DX      step for tables (here, DX=0.02d0).
* NBX     order of tables (here, NBX=7936).
*
*Parameters: output
* PARAM   exponential table characteristics.
* E00     exponential table.
* E01     exponential table.
* E10     exponential table.
* E11     exponential table.
*
*Comments:
*  Modified in order to tabulate (1-exp(-x))/x instead of (1-exp(-x))).
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER          NBX
      DOUBLE PRECISION DX
      REAL             PARAM(3), E00(0:NBX), E01(0:NBX), E10(0:NBX),
     >                 E11(0:NBX)
*----
*  LOCAL VARIABLES
*----
      INTEGER          I
      DOUBLE PRECISION X0, X1, X2, X3, X4, PAS, XLIM, EX00, EX01,
     >                 EX10, EX11
      DOUBLE PRECISION DEPS, DREF, DZERO, DONE, TAUDMIN
      PARAMETER      ( DEPS= 1.D-10, DREF= 1.D0/512.D0,
     >                 DZERO= 0.D0, DONE=1.D0, TAUDMIN=2.D-2 )
*----
*  SPACE AND VALUES FOR EXPONENTIAL TABLES
*----
      INTEGER          MEX1
      PARAMETER      ( MEX1=7936 )
*
      IF( NBX.NE.  MEX1 ) GO TO 97
      IF( DX .GT. DREF+DEPS .OR. DX .LT. DREF-DEPS  ) GO TO 98
      PAS=   DONE/DX
      XLIM=  DBLE(NBX)*DX
*----
*  WE CONSTRUCT THE LINEAR TABLES USING ACCURATE *EXP* VALUES
*----
      X0=  DZERO
      EX00= DZERO
      EX10= DONE
      DO 20 I= 0, NBX-1
*
*        STORE CONSTANT VALUE:
         X1= X0 + DX
         IF (X1.LE.TAUDMIN) THEN
            X2=0.5D0*X1
            X3=X1/3.D0
            X4=0.5D0*X2
            EX11=DONE-X2*(DONE-X3*(DONE-X4))
         ELSE
            EX11=(DONE - EXP(-X1))/X1
         ENDIF
         EX01=(DONE - EXP(-X1))
*
*        STORE STEP AND CONSTANT VALUES:
         X2=(EX01-EX00)/(X1-X0)
         X3=(EX11-EX10)/(X1-X0)
         E01(I)= REAL(X2)
         E00(I)= REAL(EX00-X0*X2)
         E11(I)= REAL(X3)
         E10(I)= REAL(EX10-X0*X3)
         X0= X1
         EX00= EX01
         EX10= EX11
   20 CONTINUE
      PARAM(1)=REAL(PAS)
      PARAM(2)=REAL(DX)
      PARAM(3)=REAL(XLIM)
      E10(NBX)= REAL(1.D0/XLIM)
      E11(NBX)= 0.0
      E00(NBX)= 1.0
      E01(NBX)= 0.0
      RETURN
*----
*  ERROR SECTION
*----
   97 CALL XABORT('XDREXP: EXP LINEAR TABLES HAVE MORE THAN 7936 ELEM'
     > //'ENTS')
   98 CALL XABORT('XDREXP: EXP LINEAR TABLES HAVE A STEP OF 1/512')
      END
