*DECK B1DIF
      SUBROUTINE B1DIF(OPTION,TYPE,NGRO,ST,SFNU,XHI,IJJ0,IJJ1,NJJ0,NJJ1,
     1 SCAT0,SCAT1,REFKEF,LFISSI,IMPX,DHOM,GAMMA,B2,ALAM1,CAET,A2,PHI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Solve the B-n equations and perform a buckling search if required.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* OPTION  type of leakage coefficient. Can be 'LKRD', 'RHS', 'B0', 'P0',
*         'B1', 'P1', 'B0TR' or 'P0TR'. 'LKRD' and 'RHS' are used to
*         impose a leakage coefficient.
* TYPE    type of buckling search. Can be 'DIFF', 'K', 'B' or 'L'.
* NGRO    number of energy groups.
* ST      macroscopic total cross sections.
* SFNU    nu*macroscopic fission cross sections.
* XHI     fission spectrum normalized to one.
* IJJ0    most thermal group in band for P0 scattering.
* NJJ0    number of groups in band for P0 scattering.
* IJJ1    most thermal group in band for P1 scattering.
* NJJ1    number of groups in band for P1 scattering.
* SCAT0   packed diffusion P0 macroscopic cross sections.
* SCAT1   packed diffusion P1 macroscopic cross sections.
* REFKEF  target K-effective for type B or type L calculations.
* LFISSI  fissile isotope flag (=.TRUE. if present).
* IMPX    print flag.
*
*Parameters: input/output
* PHI     homogeneous flux from heterogeneous calculation on input and
*         fundamental flux at output.
*
*Parameters: output
* DHOM    homogeneous leakage coefficients.
* GAMMA   gamma factors.
* B2      buckling.
* ALAM1   effective multiplication factor.
* CAET    infinite multiplication factor.
* A2      migration area.
*
*-----------------------------------------------------------------------
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER OPTION*4,TYPE*4
      INTEGER NGRO,IJJ0(NGRO),IJJ1(NGRO),NJJ0(NGRO),NJJ1(NGRO),IMPX
      REAL ST(NGRO),SFNU(NGRO),XHI(NGRO),SCAT0(*),SCAT1(*),DHOM(NGRO),
     > GAMMA(NGRO)
      DOUBLE PRECISION B2,ALAM1,CAET,A2,PHI(NGRO),REFKEF
      LOGICAL LFISSI
*----
*  LOCAL VARIABLES
*----
      PARAMETER (EPS=1.0D-6,MAXIT=50)
      DOUBLE PRECISION FFITX(MAXIT),B2ITX(MAXIT)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: CSTOC,SA
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: ASTOC,BSTOC
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(ASTOC(NGRO,NGRO+1),BSTOC(NGRO,NGRO),CSTOC(NGRO),SA(NGRO))
*
      IF((IMPX.GT.0).AND.(TYPE.EQ.'DIFF')) THEN
         WRITE (6,400)
      ELSE IF(IMPX.GT.0) THEN
         WRITE (6,410) OPTION,TYPE
      ENDIF
      IAPROX=2
      IF((OPTION.EQ.'P0').OR.(OPTION.EQ.'P1').OR.(OPTION.EQ.'P0TR'))
     > IAPROX=1
      IF((OPTION.EQ.'LKRD').OR.(OPTION.EQ.'RHS')) IAPROX=0
      BIL1=0.0D0
      DO 5 I=1,NGRO
      BIL1=BIL1+XHI(I)
      SA(I)=ST(I)
    5 CONTINUE
      IF((BIL1.GT.0.9D0).AND.(ABS(BIL1-1.0D0).GT.EPS)) THEN
        IF(IMPX.GT.0)
     >  WRITE(6,'(46H B1DIF: WARNING INCONSISTENT FISSION SPECTRUM.)')
      ENDIF
      IGAR=0
      DO 11 I=1,NGRO
      DO 10 J=IJJ0(I),IJJ0(I)-NJJ0(I)+1,-1
      IGAR=IGAR+1
      SA(J)=SA(J)-SCAT0(IGAR)
   10 CONTINUE
   11 CONTINUE
      IF(TYPE.EQ.'DIFF') THEN
         DO 15 I=1,NGRO
         ST2=DBLE(ST(I))
         GAMMA(I)=REAL(B1GAMA(IAPROX,B2,ST2))
         IF(IAPROX.NE.0) DHOM(I)=REAL(1.0D0/(3.0D0*GAMMA(I)*ST2))
   15    CONTINUE
         RETURN
      ELSE IF((TYPE.EQ.'B').OR.(TYPE.EQ.'L')) THEN
*        COMPUTE THE INITIAL BUCKLING.
         BIL1=0.0D0
         BIL2=0.0D0
         IF(LFISSI) THEN
            DO 20 I=1,NGRO
            BIL1=BIL1+(SFNU(I)/REFKEF)*PHI(I)
   20       CONTINUE
         ENDIF
         DO 21 I=1,NGRO
         BIL1=BIL1-SA(I)*PHI(I)
         BIL2=BIL2+DHOM(I)*PHI(I)
   21    CONTINUE
         B2=BIL1/BIL2
         DO 25 I=1,NGRO
         ST2=-0.7D0*(ST(I)**2)
         IF(B2.LT.ST2) THEN
            IF(IMPX.GT.0) WRITE (6,415) B2,ST2
            B2=ST2
         ENDIF
         PHI(I)=1.0D0
   25    CONTINUE
      ENDIF
      IF(IMPX.GT.1) WRITE(6,420) B2
      IF(TYPE.EQ.'L') GO TO 160
*----
*  COMPUTE THE FUNDAMENTAL FLUX WITH TYPE K OR TYPE B
*----
      ITEX=0
   30 ITEX=ITEX+1
      IF(ITEX.GT.MAXIT) CALL XABORT('B1DIF: UNABLE TO CONVERGE(1).')
      IGAR=0
      DO 55 I=1,NGRO
      ST2=DBLE(ST(I))
      DD=DBLE(DHOM(I))
      BETA=B1BETA(IAPROX,B2,ST2,DD)
      DO 40 J=1,NGRO
      ASTOC(I,J)=0.0D0
      BSTOC(I,J)=0.0D0
   40 CONTINUE
      ASTOC(I,I)=ST2
      ASTOC(I,NGRO+1)=(1.0D0-B2*BETA)*XHI(I)
      DO 50 J=IJJ0(I),IJJ0(I)-NJJ0(I)+1,-1
      IGAR=IGAR+1
      ASTOC(I,J)=ASTOC(I,J)-(1.0D0-B2*BETA)*SCAT0(IGAR)
      BSTOC(I,J)=SCAT0(IGAR)
   50 CONTINUE
   55 CONTINUE
      IF((OPTION.EQ.'P1').OR.(OPTION.EQ.'B1')) THEN
         DO 72 J=1,NGRO
         DO 60 K=1,NGRO
         CSTOC(K)=BSTOC(K,J)
         BSTOC(K,J)=0.0D0
   60    CONTINUE
         IGAR=0
         DO 71 I=1,NGRO
         DO 70 K=IJJ1(I),IJJ1(I)-NJJ1(I)+1,-1
         IGAR=IGAR+1
         BSTOC(I,J)=BSTOC(I,J)+3.0D0*SCAT1(IGAR)*CSTOC(K)
   70    CONTINUE
   71    CONTINUE
   72    CONTINUE
         IGAR=0
         DO 95 I=1,NGRO
         ST2=DBLE(ST(I))
         DD=DBLE(DHOM(I))
         BETA=B1BETA(IAPROX,B2,ST2,DD)*ST2
         DO 80 J=1,NGRO
         ASTOC(I,J)=ASTOC(I,J)+BETA*BSTOC(I,J)
   80    CONTINUE
         DO 90 J=IJJ1(I),IJJ1(I)-NJJ1(I)+1,-1
         IGAR=IGAR+1
         ASTOC(I,NGRO+1)=ASTOC(I,NGRO+1)-3.0D0*BETA*SCAT1(IGAR)*XHI(J)
         ASTOC(I,J)=ASTOC(I,J)-3.0D0*BETA*SCAT1(IGAR)*ST(J)
   90    CONTINUE
   95    CONTINUE
      ENDIF
      CALL B1SOL(NGRO,ASTOC,IER)
      IF(IER.NE.0) CALL XABORT('B1DIF: SINGULAR MATRIX(1).')
      ALAM1=0.0D0
      CAET=0.0D0
      DO 130 I=1,NGRO
      ALAM1=ALAM1+SFNU(I)*ASTOC(I,NGRO+1)
      CAET=CAET+SA(I)*ASTOC(I,NGRO+1)
  130 CONTINUE
      IF(IMPX.GT.1) WRITE (6,430) ITEX,ALAM1,B2
      B2ITX(ITEX)=B2
      FFITX(ITEX)=REFKEF-ALAM1
      IF(TYPE.EQ.'K') THEN
         DO 140 I=1,NGRO
         PHI(I)=ASTOC(I,NGRO+1)/ALAM1
  140    CONTINUE
      ELSE IF(TYPE.EQ.'B') THEN
*        COMPUTE THE EXTRAPOLATED BUCKLING.
         IF(ITEX.LE.5) THEN
*          USE A BALANCE RELATION.
           B2=B2*(ALAM1/REFKEF-CAET)/(1.0D0-CAET)
         ELSE
           IF(ITEX.EQ.6) THEN
*            SORT THE ROOT CONVERGENCE HISTORY.
             DO I=2,ITEX-1
               WORKF=FFITX(ITEX-I)
               WORKB=B2ITX(ITEX-I)
               J=I
               DO WHILE((J.GT.0).AND.
     >                  (ABS(FFITX(ITEX-J+1)).GT.ABS(WORKF)))
                 FFITX(ITEX-J)=FFITX(ITEX-J+1)
                 B2ITX(ITEX-J)=B2ITX(ITEX-J+1)
                 J=J-1
               ENDDO
               FFITX(ITEX-J)=WORKF
               B2ITX(ITEX-J)=WORKB
             ENDDO
           ENDIF
           J=0
           DO I=ITEX-1,2,-1
             IF(FFITX(I)*FFITX(ITEX).LT.0.0) THEN
               J=I
               EXIT
             ENDIF
           ENDDO
           IF(J.NE.0) THEN
*            USE A BISSECTION METHOD.
             B2=0.5D0*(B2ITX(J)+B2ITX(ITEX))
           ELSE
*            USE THE SECANT METHOD.
             AA=FFITX(ITEX)-FFITX(ITEX-1)
             B2=(B2ITX(ITEX-1)*FFITX(ITEX)-B2ITX(ITEX)*FFITX(ITEX-1))/AA
           ENDIF
         ENDIF
*        CHECK THE CONVERGENCE.
         BIL1=0.0D0
         BIL2=0.0D0
         DO 150 I=1,NGRO
         ST2=ST(I)**2
         BIL1=MAX(BIL1,ABS(ASTOC(I,NGRO+1)/ALAM1))
         BIL2=MAX(BIL2,ABS(PHI(I)-ASTOC(I,NGRO+1)/ALAM1))
         PHI(I)=ASTOC(I,NGRO+1)/ALAM1
  150    CONTINUE
         ERR3=ABS(REFKEF-ALAM1)
         IF((BIL2.GE.10*EPS*BIL1).OR.(ERR3.GE.EPS)) GO TO 30
      ENDIF
      GO TO 300
*----
*  COMPUTE THE FUNDAMENTAL FLUX WITH TYPE L
*----
  160 ITEX=0
  170 ITEX=ITEX+1
      IF(ITEX.GT.MAXIT) CALL XABORT('B1DIF: UNABLE TO CONVERGE(2).')
      IF((OPTION.EQ.'P1').OR.(OPTION.EQ.'B1')) THEN
         IGAR=0
         DO 200 I=1,NGRO
         ST2=DBLE(ST(I))
         DD=DBLE(DHOM(I))
         BETA=B1BETA(IAPROX,B2,ST2,DD)
         BETA=BETA*ST2/(1.0D0-B2*BETA)
         DO 180 J=1,NGRO
         ASTOC(I,J)=0.0D0
  180    CONTINUE
         ASTOC(I,I)=1.0D0
         ASTOC(I,NGRO+1)=BETA
         DO 190 J=IJJ1(I),IJJ1(I)-NJJ1(I)+1,-1
         IGAR=IGAR+1
         ASTOC(I,J)=ASTOC(I,J)-3.0D0*BETA*SCAT1(IGAR)*PHI(J)/PHI(I)
  190    CONTINUE
  200    CONTINUE
         CALL B1SOL(NGRO,ASTOC,IER)
         IF(IER.NE.0) CALL XABORT('B1DIF: SINGULAR MATRIX(2).')
         DO 210 I=1,NGRO
         DHOM(I)=REAL(ASTOC(I,NGRO+1))
  210    CONTINUE
      ELSE IF((OPTION.NE.'LKRD').AND.(OPTION.NE.'RHS')) THEN
         DO 220 I=1,NGRO
         ST2=ST(I)
         GAMMA(I)=REAL(B1GAMA(IAPROX,B2,ST2))
         IF(IAPROX.NE.0) DHOM(I)=REAL(1.0D0/(3.0D0*GAMMA(I)*ST2))
  220    CONTINUE
      ENDIF
      ASTOC(:NGRO,:NGRO)=0.0D0
      BSTOC(:NGRO,:NGRO)=0.0D0
      ASTOC(:NGRO,NGRO+1)=PHI(:NGRO)
      IGAR=0
      DO 250 I=1,NGRO
      ASTOC(I,I)=ASTOC(I,I)+ST(I)
      BSTOC(I,I)=BSTOC(I,I)-DHOM(I)
      DO 230 J=IJJ0(I),IJJ0(I)-NJJ0(I)+1,-1
      IGAR=IGAR+1
      ASTOC(I,J)=ASTOC(I,J)-SCAT0(IGAR)
  230 CONTINUE
      IF(LFISSI) THEN
         DO 240 J=1,NGRO
         ASTOC(I,J)=ASTOC(I,J)-XHI(I)*SFNU(J)/REFKEF
  240    CONTINUE
      ENDIF
  250 CONTINUE
      B2OLD=B2
      CALL ALEIGD(ASTOC,BSTOC,NGRO,B2,ASTOC(1,NGRO+1),EPS,IT)
      IF(IMPX.GT.1) WRITE (6,440) ITEX,IT,B2
      BIL1=SQRT(DOT_PRODUCT(ASTOC(:NGRO,NGRO+1),ASTOC(:NGRO,NGRO+1)))
      BIL2=0.0D0
      DO 260 I=1,NGRO
      BIL2=MAX(BIL2,ABS(PHI(I)-ASTOC(I,NGRO+1)/BIL1))
      PHI(I)=0.5*(PHI(I)+ASTOC(I,NGRO+1)/BIL1)
  260 CONTINUE
      ERR3=ABS(B2-B2OLD)
      IF((BIL2.GE.10.0*EPS*BIL1).OR.(ERR3.GE.EPS)) GO TO 170
*----
*  COMPUTE THE LEAKAGE COEFFICIENTS
*----
  300 IF((OPTION.EQ.'P1').OR.(OPTION.EQ.'B1')) THEN
         IGAR=0
         DO 325 I=1,NGRO
         ST2=DBLE(ST(I))
         DD=DBLE(DHOM(I))
         BETA=B1BETA(IAPROX,B2,ST2,DD)
         BETA=BETA*ST2/(1.0D0-B2*BETA)
         DO 310 J=1,NGRO
         ASTOC(I,J)=0.0D0
  310    CONTINUE
         ASTOC(I,I)=1.0D0
         ASTOC(I,NGRO+1)=BETA
         DO 320 J=IJJ1(I),IJJ1(I)-NJJ1(I)+1,-1
         IGAR=IGAR+1
         ASTOC(I,J)=ASTOC(I,J)-3.0D0*BETA*SCAT1(IGAR)*PHI(J)/PHI(I)
  320    CONTINUE
  325    CONTINUE
         CALL B1SOL(NGRO,ASTOC,IER)
         IF(IER.NE.0) CALL XABORT('B1DIF: SINGULAR MATRIX(2).')
         DO 330 I=1,NGRO
         DHOM(I)=REAL(ASTOC(I,NGRO+1))
  330    CONTINUE
      ELSE IF((OPTION.NE.'LKRD').AND.(OPTION.NE.'RHS')) THEN
         DO 340 I=1,NGRO
         ST2=ST(I)
         GAMMA(I)=REAL(B1GAMA(IAPROX,B2,ST2))
         IF(IAPROX.NE.0) DHOM(I)=REAL(1.0D0/(3.0D0*GAMMA(I)*ST2))
  340    CONTINUE
      ENDIF
      A2=0.0D0
      CAET=0.0D0
      ZXC=0.0D0
      DO 350 I=1,NGRO
      A2=A2+DHOM(I)*PHI(I)
      CAET=CAET+SA(I)*PHI(I)
      IF(LFISSI) ZXC=ZXC+SFNU(I)*PHI(I)/REFKEF
      ST2=DBLE(ST(I))
      GAMMA(I)=REAL(B1GAMA(IAPROX,B2,ST2))
  350 CONTINUE
      A2=A2/CAET
      CAET=REFKEF*ZXC/CAET
      IF(CAET.EQ.0.0) THEN
        ALAM2=0.0
      ELSE
        ALAM2=CAET/(1.0D0+A2*B2)
      ENDIF
      IF(TYPE.EQ.'L') ALAM1=ALAM2
      IF(IMPX.GT.0) WRITE (6,450) ITEX,B2,ALAM1,ALAM2,CAET,A2
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(SA,CSTOC,BSTOC,ASTOC)
      RETURN
*
  400 FORMAT(/42H B1DIF: DIFFUSION COEFFICIENT CALCULATION.)
  410 FORMAT(/20H B1DIF: SOLUTION OF ,A4,21H EQUATIONS WITH TYPE ,A4)
  415 FORMAT(47H B1DIF: THE INITIAL BUCKLING WAS INCREASED FROM,1P,
     1 E13.5,3H TO,E13.5)
  420 FORMAT(26H B1DIF: INITIAL BUCKLING =,1P,E13.5)
  430 FORMAT(33H B1DIF: K-EFFECTIVE ITERATION NO.,I3,13H. K-EFFECTIVE,
     1 2H =,F10.6,11H BUCKLING =,1P,E13.5)
  440 FORMAT(30H B1DIF: BUCKLING ITERATION NO.,I3,13H CONVERGED IN,I5,
     1 29H INNER ITERATIONS. BUCKLING =,1P,E13.5)
  450 FORMAT(8X,22HNUMBER OF ITERATIONS =,I3/8X,10HBUCKLING =,1P,E13.5,
     1 0P/8X,13HK-EFFECTIVE =,F10.6,3H  (,F10.6,2H )/8X,12HK-INFINITE =,
     2 F10.6/8X,16HMIGRATION AREA =,1P,E13.5/)
      END
