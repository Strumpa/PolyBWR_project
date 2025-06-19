*DECK XDRKIN
      SUBROUTINE XDRKIN(DX,NBX,MLOG,BIV,PASV,XLIMV)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Construct Bickley tables for KI1(X), KI2(X), KI3(X), KI4(X), KI5(X),
* taking into account logarithmic singularities for KI1(X) and KI2(X).
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* DX      step for tables (here, DX=0.02d0).
* NBX     nb element in tables (here, NBX=600).
* MLOG    interval for logarithmic singularities (suggested values:
*         MLOG(1)=30, MLOG(2)=15, MLOG(3)= 0, MLOG(4)= 0, MLOG(5)= 0).
*
*Parameters: output
* BIV     elements of quadratic BICKLEY table.
* PASV    step quadratic of BICKLEY table.
* XLMV    upper limit of quadratic BICKLEY table.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
***** CALLS:    *AKIN10* ROUTINE FOR ACCURATE KIN(X) BICKLEY VALUES
*               *AK0BES* ROUTINE FOR ACCURATE K0(X)   BESSEL VALUES
*               *AK1BES* ROUTINE FOR ACCURATE K1(X)   BESSEL VALUES
*----
*  SUBROUTINE ARGUMENTS
*----
      DOUBLE PRECISION DX
      INTEGER NBX,MLOG(5)
      REAL BIV(0:NBX,3,5),PASV(5),XLIMV(5)
*----
*  LOCAL VARIABLES
*----
      INTEGER  I,KI,IORD
      DOUBLE PRECISION X, AKIN(-1:10), AK0BES, AK1BES
      DOUBLE PRECISION GAMMA, PIO2, PAS, XLIM
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: C0, C1, C2
      PARAMETER ( GAMMA=0.57721566490153D0,
     >            PIO2= 1.57079632679490D0 )
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(C0(5,0:NBX), C1(5,0:NBX), C2(5,0:NBX))
*
      IF( MLOG(3).GT.0.OR.MLOG(4).GT.0.OR.MLOG(5).GT.0 )GOTO 99
*
      IF( NBX .NE.   600 )GOTO 97
      IF( DX  .NE.0.02D0 )GOTO 98
      PAS=   1.0D0/DX
      XLIM=  DBLE(NBX)*DX
      DO 10 I=1,5
      XLIMV(I)=REAL(XLIM)
      PASV(I)=REAL(PAS)
   10 CONTINUE
*----
*  FIRST, WE CONSTRUCT THE TABLES USING ACCURATE *AKIN10* VALUES
*----
      X= 0.D0
      CALL AKIN10(X,AKIN(1))
      AKIN( 0)= 0.D0
      AKIN(-1)= 0.D0
      DO 30 I= 0, NBX-1
         DO 20 KI= 1, 5
            C2(KI,I)=   0.5D0 * AKIN(KI-2)
            C1(KI,I)= -(AKIN(KI-1)+X*AKIN(KI-2))
            C0(KI,I)=   AKIN(KI)+X*(AKIN(KI-1)+X*C2(KI,I))
   20    CONTINUE
         X= X + DX
         CALL AKIN10(X,AKIN(1))
         AKIN( 0)= AK0BES(X)
         AKIN(-1)= AK1BES(X)
   30 CONTINUE
      DO 40 KI= 1, 5
         C0(KI,NBX)= 0.D0
         C1(KI,NBX)= 0.D0
         C2(KI,NBX)= 0.D0
   40 CONTINUE
*----
*  KI1(X) ADJUSTMENTS
*----
      X= 0.D0
      DO 50 I= 1, MLOG(1)-1
         X= X + DX
         C0(1,I)= C0(1,I) + 0.5D0*X
         C1(1,I)= C1(1,I) - LOG(X)
         C2(1,I)= C2(1,I) - 0.5D0/X
   50 CONTINUE
*
      C1(1,0)= GAMMA-(LOG(2.D0)+1.D0)
*----
*  KI2(X) ADJUSTMENTS
*----
      X= 0.D0
      DO 60 I= 1, MLOG(2)-1
         X= X + DX
         C0(2,I)= C0(2,I) + 0.25D0*X*X
         C1(2,I)= C1(2,I) - X
         C2(2,I)= C2(2,I) + 0.5D0*LOG(X) +0.75D0
   60 CONTINUE
*
      C1(2,0)=  -PIO2
      C2(2,0)=   0.5D0*(LOG(2.D0)+1.5D0-GAMMA)
*----
*  OUTPUT VALUES
*----
      DO 80 I= 0, NBX
        DO 70 IORD=1,5
          BIV(I,1,IORD)= REAL(C0(IORD,I))
          BIV(I,2,IORD)= REAL(C1(IORD,I))
          BIV(I,3,IORD)= REAL(C2(IORD,I))
   70   CONTINUE
   80 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(C2, C1, C0)
      RETURN
*----
*  ERROR SECTION
*----
   97 CALL XABORT('XDRKIN: KIN TABLES HAVE MORE THAN 600 ELEMENTS')
   98 CALL XABORT('XDRKIN: KIN TABLES HAVE A STEP OF 0.02')
   99 CALL XABORT('XDRKIN: NO LOG SINGULARITY TAKEN FOR KI345')
      END
