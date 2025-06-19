*DECK SNFT13
      SUBROUTINE SNFT13(NUN,NGEFF,IMPX,INCONV,NGIND,LX,LY,LZ,IELEM,
     1 NMAT,NPQ,NSCT,MAT,VOL,TOTAL,NCODE,ZCODE,QEXT,LFIXUP,DU,DE,DZ,
     2 W,MRMX,MRMY,MRMZ,DC,DB,DA,FUNKNO,ISBS,NBS,ISBSM,BS,MAXL,MN,DN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform one inner iteration for solving SN equations in 3D Cartesian
* geometry for the HODD method. Energy-angle multithreading. Albedo
* boundary conditions.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): N. Martin
*
*Parameters: input
* NUN     total number of unknowns in vector FUNKNO.
* NGEFF   number of energy groups processed in parallel.
* IMPX    print flag (equal to zero for no print).
* INCONV  energy group convergence flag (set to .FALSE. if converged).
* NGIND   energy group indices assign to the NGEFF set.
* LX      number of meshes along X axis.
* LY      number of meshes along Y axis.
* LZ      number of meshes along Z axis.
* IELEM   degree of the spatial approximation (=1 diamond scheme).
* NMAT    number of material mixtures.
* NPQ     number of SN directions in height octants.
* NSCT    maximum number of spherical harmonics moments of the flux.
* MAT     material mixture index in each region.
* VOL     volumes of each region.
* TOTAL   macroscopic total cross sections.
* NCODE   boundary condition indices.
* ZCODE   albedos.
* QEXT    Legendre components of the fixed source.
* LFIXUP  flag to enable negative flux fixup.
* DU      first direction cosines ($\\mu$).
* DE      second direction cosines ($\\eta$).
* DZ      third direction cosines ($\\xi$).
* W       weights.
* MRMX    quadrature index.
* MRMY    quadrature index.
* MRMZ    quadrature index.
* DC      diamond-scheme parameter.
* DB      diamond-scheme parameter.
* DA      diamond-scheme parameter.
* MN      moment-to-discrete matrix.
* DN      discrete-to-moment matrix.
*
*Parameters: input/output
* FUNKNO  Legendre components of the flux and boundary fluxes.
*
*-----------------------------------------------------------------------
*
#if defined(_OPENMP)
      USE omp_lib
#endif
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NUN,NGEFF,IMPX,NGIND(NGEFF),LX,LY,LZ,IELEM,NMAT,NPQ,
     1 NSCT,MAT(LX,LY,LZ),NCODE(6),MRMX(NPQ),MRMY(NPQ),MRMZ(NPQ),ISBS,
     2 NBS,ISBSM(6*ISBS,NPQ*ISBS,NGEFF*ISBS),MAXL
      LOGICAL INCONV(NGEFF)
      REAL VOL(LX,LY,LZ),TOTAL(0:NMAT,NGEFF),ZCODE(6),QEXT(NUN,NGEFF),
     1 DU(NPQ),DE(NPQ),DZ(NPQ),W(NPQ),DC(LX,LY,NPQ),DB(LX,LZ,NPQ),
     2 DA(LY,LZ,NPQ),FUNKNO(NUN,NGEFF),BS(MAXL*ISBS,NBS*ISBS),
     3 MN(NPQ,NSCT),DN(NSCT,NPQ)
      LOGICAL LFIXUP
*----
*  LOCAL VARIABLES
*----
      INTEGER NPQD(8),IIND(8),P
      PARAMETER(IUNOUT=6,RLOG=1.0E-8,PI=3.141592654)
      DOUBLE PRECISION VT,CONST0,CONST1,CONST2,Q(IELEM**3),
     1 Q2(IELEM**3,IELEM**3+1),XNK(IELEM,IELEM)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: INDANG
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: FLUX
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: FLUX_G
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: XNI
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: XNJ
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(INDANG(NPQ,8))
      ALLOCATE(XNI(IELEM,IELEM,LY,LZ),XNJ(IELEM,IELEM,LZ))
      ALLOCATE(FLUX(IELEM**3,NSCT,LX,LY,LZ),
     1 FLUX_G(IELEM**3,NSCT,LX,LY,LZ,NGEFF))
*----
*  DEFINITION OF CONSTANTS.
*----
      L4=IELEM*IELEM*IELEM*LX*LY*LZ*NSCT
      IOFY=IELEM*IELEM*LY*LZ*NPQ
      IOFZ=IELEM*IELEM*(LX+LY)*LZ*NPQ
      CONST0=2.0D0*DSQRT(3.0D0)
      CONST1=2.0D0*DSQRT(5.0D0)
      CONST2=2.0D0*DSQRT(15.0D0)
*----
*  PARAMETER VALIDATION.
*----
      IF(IELEM.GT.3) CALL XABORT('SNFT13: INVALID IELEM (DIAM) VALUE. '
     1 //'CHECK INPUT DATA FILE.')
      FLUX_G(:IELEM**3,:NSCT,:LX,:LY,:LZ,:NGEFF)=0.0D0
*----
*  SET OCTANT SWAPPING ORDER.
*----
      NPQD(:8)=0
      INDANG(:NPQ,:8)=0
      IIND(:)=0
      DO 10 M=1,NPQ
        VU=DU(M)
        VE=DE(M)
        VZ=DZ(M)
        IF(W(M).EQ.0) CYCLE
        IF((VU.GE.0.0).AND.(VE.GE.0.0).AND.(VZ.GE.0.0)) THEN
          IND=1
          JND=8
        ELSE IF((VU.LE.0.0).AND.(VE.GE.0.0).AND.(VZ.GE.0.0)) THEN
          IND=2
          JND=7
        ELSE IF((VU.LE.0.0).AND.(VE.LE.0.0).AND.(VZ.GE.0.0)) THEN
          IND=3
          JND=5
        ELSE IF((VU.GE.0.0).AND.(VE.LE.0.0).AND.(VZ.GE.0.0)) THEN
          IND=4
          JND=6
        ELSE IF((VU.GE.0.0).AND.(VE.GE.0.0).AND.(VZ.LE.0.0)) THEN
          IND=5
          JND=4
        ELSE IF((VU.LE.0.0).AND.(VE.GE.0.0).AND.(VZ.LE.0.0)) THEN
          IND=6
          JND=3
        ELSE IF((VU.LE.0.0).AND.(VE.LE.0.0).AND.(VZ.LE.0.0)) THEN
          IND=7
          JND=1
        ELSE
          IND=8
          JND=2
        ENDIF
        IIND(JND)=IND
        NPQD(IND)=NPQD(IND)+1
        INDANG(NPQD(IND),IND)=M
   10 CONTINUE
*----
*  MAIN LOOP OVER OCTANTS.
*----
      DO 420 JND=1,8
      IND=IIND(JND)
*----
*  PRELIMINARY LOOPS FOR SETTING BOUNDARY CONDITIONS.
*----
*$OMP  PARALLEL DO
*$OMP1 PRIVATE(M,IG,VU,VE,VZ,M1,E1,IOF,JOF,IEL,JEL,I,J,K)
*$OMP2 SHARED(FUNKNO) COLLAPSE(2)
      DO 150 IG=1,NGEFF
      DO 140 IPQD=1,NPQD(IND)
      IF(.NOT.INCONV(IG)) GO TO 140
      M=INDANG(IPQD,IND)
      VU=DU(M)
      VE=DE(M)
      VZ=DZ(M)
*
      IF(VU.GT.0.0)THEN
         M1=MRMX(M)
         IF(NCODE(1).NE.4)THEN
         DO JEL=1,IELEM
         DO IEL=1,IELEM
            DO J=1,LY
            DO K=1,LZ
             IOF=((((M-1)*LZ+(K-1))*LY+(J-1))*IELEM+(JEL-1))*IELEM+IEL
             JOF=((((M1-1)*LZ+(K-1))*LY+(J-1))*IELEM+(JEL-1))*IELEM+IEL
             FUNKNO(L4+IOF,IG)=FUNKNO(L4+JOF,IG)
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ENDIF
      ELSEIF(VU.LT.0.0)THEN
         M1=MRMX(M)
         IF(NCODE(2).NE.4)THEN
         DO JEL=1,IELEM
         DO IEL=1,IELEM
            DO J=1,LY
            DO K=1,LZ
             IOF=((((M-1)*LZ+(K-1))*LY+(J-1))*IELEM+(JEL-1))*IELEM+IEL
             JOF=((((M1-1)*LZ+(K-1))*LY+(J-1))*IELEM+(JEL-1))*IELEM+IEL
             FUNKNO(L4+IOF,IG)=FUNKNO(L4+JOF,IG)
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ENDIF
      ENDIF
      IF(VE.GT.0.0)THEN
         M1=MRMY(M)
         IF(NCODE(3).NE.4)THEN
         DO IEL=1,IELEM
         DO JEL=1,IELEM
            DO I=1,LX
            DO K=1,LZ
             IOF=((((M-1)*LZ+(K-1))*LX+(I-1))*IELEM+(JEL-1))*IELEM+IEL
             JOF=((((M1-1)*LZ+(K-1))*LX+(I-1))*IELEM+(JEL-1))*IELEM+IEL
             FUNKNO(L4+IOFY+IOF,IG)=FUNKNO(L4+IOFY+JOF,IG)
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ENDIF
      ELSEIF(VE.LT.0.0)THEN
         M1=MRMY(M)
         IF(NCODE(4).NE.4)THEN
         DO IEL=1,IELEM
         DO JEL=1,IELEM
            DO I=1,LX
            DO K=1,LZ
             IOF=((((M-1)*LZ+(K-1))*LX+(I-1))*IELEM+(JEL-1))*IELEM+IEL
             JOF=((((M1-1)*LZ+(K-1))*LX+(I-1))*IELEM+(JEL-1))*IELEM+IEL
             FUNKNO(L4+IOFY+IOF,IG)=FUNKNO(L4+IOFY+JOF,IG)
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ENDIF
      ENDIF
      IF(VZ.GT.0.0)THEN
         M1=MRMZ(M)
         IF(NCODE(5).NE.4)THEN
         DO IEL=1,IELEM
         DO JEL=1,IELEM
            DO I=1,LX
            DO J=1,LY
             IOF=((((M-1)*LY+(J-1))*LX+(I-1))*IELEM+(JEL-1))*IELEM+IEL
             JOF=((((M1-1)*LY+(J-1))*LX+(I-1))*IELEM+(JEL-1))*IELEM+IEL
             E1=FUNKNO(L4+IOFZ+IOF,IG)
             FUNKNO(L4+IOFZ+IOF,IG)=FUNKNO(L4+IOFZ+JOF,IG)
             FUNKNO(L4+IOFZ+JOF,IG)=E1
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ENDIF
      ELSEIF(VZ.LT.0.0)THEN
         M1=MRMZ(M)
         IF(NCODE(6).NE.4)THEN
         DO IEL=1,IELEM
         DO JEL=1,IELEM
            DO I=1,LX
            DO J=1,LY
             IOF=((((M-1)*LY+(J-1))*LX+(I-1))*IELEM+(JEL-1))*IELEM+IEL
             JOF=((((M1-1)*LY+(J-1))*LX+(I-1))*IELEM+(JEL-1))*IELEM+IEL
             E1=FUNKNO(L4+IOFZ+IOF,IG)
             FUNKNO(L4+IOFZ+IOF,IG)=FUNKNO(L4+IOFZ+JOF,IG)
             FUNKNO(L4+IOFZ+JOF,IG)=E1
            ENDDO
            ENDDO
         ENDDO
         ENDDO
         ENDIF
      ENDIF
*
  140 CONTINUE
  150 CONTINUE
*$OMP END PARALLEL DO
*----
*  MAIN SWAPPING LOOPS FOR SN FLUX CALCULATION
*----
*$OMP  PARALLEL DO
*$OMP1 PRIVATE(ITID,FLUX,M,IG,XNI,XNJ,XNK,Q,Q2,IOF,IER,IX,IY,IZ,IEL,JEL)
*$OMP2 PRIVATE(I,J,K,L,VT) SHARED(FUNKNO) REDUCTION(+:FLUX_G)
*$OMP3 COLLAPSE(2)
      DO 410 IG=1,NGEFF
      DO 400 IPQD=1,NPQD(IND)
#if defined(_OPENMP)
        ITID=omp_get_thread_num()
#else
        ITID=0
#endif
      IF(IMPX.GT.5) WRITE(IUNOUT,500) ITID,NGIND(IG),IPQD
      IF(.NOT.INCONV(IG)) GO TO 400
      M=INDANG(IPQD,IND)
      FLUX(:IELEM**3,:NSCT,:LX,:LY,:LZ)=0.0D0
      IF(W(M).EQ.0.0) GO TO 400
*----
*  LOOP OVER X-, Y- AND Z-DIRECTED AXES.
*----
      DO 350 IX=1,LX
      I=IX
      IF((IND.EQ.2).OR.(IND.EQ.3).OR.(IND.EQ.6).OR.(IND.EQ.7)) I=LX+1-I
      DO 310 IY=1,LY
      J=IY
      IF((IND.EQ.3).OR.(IND.EQ.4).OR.(IND.EQ.7).OR.(IND.EQ.8)) J=LY+1-J
      DO 170 JEL=1,IELEM
      DO 160 IEL=1,IELEM
      IOF=((((M-1)*LY+(J-1))*LX+(I-1))*IELEM+(JEL-1))*IELEM+IEL
      IF((IND.EQ.1).OR.(IND.EQ.2).OR.(IND.EQ.3).OR.(IND.EQ.4)) THEN
        XNK(IEL,JEL)=FUNKNO(L4+IOFZ+IOF,IG)*ZCODE(5)
      ELSE
        XNK(IEL,JEL)=FUNKNO(L4+IOFZ+IOF,IG)*ZCODE(6)
      ENDIF
  160 CONTINUE
  170 CONTINUE
      IF(ISBS.EQ.1) THEN
        IF(((IND.EQ.5).OR.(IND.EQ.6).OR.(IND.EQ.7).OR.(IND.EQ.8))
     1  .AND.ISBSM(6,M,IG).NE.0) THEN
          XNK(1,1)=XNK(1,1)+BS((I-1)*LY+J,ISBSM(6,M,IG)) 
        ELSEIF(((IND.EQ.1).OR.(IND.EQ.2).OR.(IND.EQ.3).OR.(IND.EQ.4))
     1  .AND.ISBSM(5,M,IG).NE.0) THEN
          XNK(1,1)=XNK(1,1)+BS((I-1)*LY+J,ISBSM(5,M,IG)) 
        ENDIF
      ENDIF
      DO 280 IZ=1,LZ
      K=IZ
      IF((IND.EQ.5).OR.(IND.EQ.6).OR.(IND.EQ.7).OR.(IND.EQ.8)) K=LZ+1-K
      IF(IY.EQ.1) THEN
        DO 190 JEL=1,IELEM
        DO 180 IEL=1,IELEM
        IOF=((((M-1)*LZ+(K-1))*LX+(I-1))*IELEM+(JEL-1))*IELEM+IEL
        IF((IND.EQ.1).OR.(IND.EQ.2).OR.(IND.EQ.5).OR.(IND.EQ.6)) THEN
          XNJ(IEL,JEL,K)=FUNKNO(L4+IOFY+IOF,IG)*ZCODE(3)
        ELSE
          XNJ(IEL,JEL,K)=FUNKNO(L4+IOFY+IOF,IG)*ZCODE(4)
        ENDIF
  180   CONTINUE
  190   CONTINUE
        IF(ISBS.EQ.1) THEN
          IF(((IND.EQ.3).OR.(IND.EQ.4).OR.(IND.EQ.7).OR.(IND.EQ.8))
     1    .AND.ISBSM(4,M,IG).NE.0) THEN
            XNJ(1,1,K)=XNJ(1,1,K)+BS((I-1)*LZ+K,ISBSM(4,M,IG)) 
          ELSEIF(((IND.EQ.1).OR.(IND.EQ.2).OR.(IND.EQ.5).OR.(IND.EQ.6))
     1    .AND.ISBSM(3,M,IG).NE.0) THEN
            XNJ(1,1,K)=XNJ(1,1,K)+BS((I-1)*LZ+K,ISBSM(3,M,IG)) 
          ENDIF
        ENDIF
      ENDIF
      IF(IX.EQ.1) THEN
        DO 210 JEL=1,IELEM
        DO 200 IEL=1,IELEM
        IOF=((((M-1)*LZ+(K-1))*LY+(J-1))*IELEM+(JEL-1))*IELEM+IEL
        IF((IND.EQ.1).OR.(IND.EQ.4).OR.(IND.EQ.5).OR.(IND.EQ.8)) THEN
          XNI(IEL,JEL,J,K)=FUNKNO(L4+IOF,IG)*ZCODE(1)
        ELSE
          XNI(IEL,JEL,J,K)=FUNKNO(L4+IOF,IG)*ZCODE(2)
        ENDIF
  200   CONTINUE
  210   CONTINUE
        IF(ISBS.EQ.1) THEN
          IF(((IND.EQ.2).OR.(IND.EQ.3).OR.(IND.EQ.6).OR.(IND.EQ.7))
     1    .AND.ISBSM(2,M,IG).NE.0) THEN
            XNI(1,1,J,K)=XNI(1,1,J,K)+BS((J-1)*LZ+K,ISBSM(2,M,IG)) 
          ELSEIF(((IND.EQ.1).OR.(IND.EQ.4).OR.(IND.EQ.5).OR.(IND.EQ.8))
     1    .AND.ISBSM(1,M,IG).NE.0) THEN
            XNI(1,1,J,K)=XNI(1,1,J,K)+BS((J-1)*LZ+K,ISBSM(1,M,IG)) 
          ENDIF
        ENDIF
      ENDIF
      IF(MAT(I,J,K).EQ.0) GO TO 280
      DO 230 IEL=1,IELEM**3
      Q(IEL)=0.0D0
      DO 220 P=1,NSCT
      IOF=((((K-1)*LY+(J-1))*LX+(I-1))*NSCT+(P-1))*IELEM**3+IEL
      Q(IEL)=Q(IEL)+QEXT(IOF,IG)*MN(M,P)
  220 CONTINUE
  230 CONTINUE
      VT=VOL(I,J,K)*TOTAL(MAT(I,J,K),IG)
      Q2(:IELEM**3,:(IELEM**3)+1)=0.0D0
*---------------------------------------
      IF(IELEM.EQ.1) THEN
        Q2(1,1)=2.0D0*ABS(DA(J,K,M))+2.0D0*ABS(DB(I,K,M))+
     1        2.0D0*ABS(DC(I,J,M))+VT
        Q2(1,2)=2.0D0*ABS(DA(J,K,M))*XNI(1,1,J,K)+
     1        2.0D0*ABS(DB(I,K,M))*XNJ(1,1,K)+
     2        2.0D0*ABS(DC(I,J,M))*XNK(1,1)+VOL(I,J,K)*Q(1)
      ELSE IF(IELEM.EQ.2) THEN
      Q2(1,1) = VT
      Q2(1,2) = 2*3**(0.5D0)*DA(J,K,M)
      Q2(1,3) = 2*3**(0.5D0)*DB(I,K,M)
      Q2(1,5) = 2*3**(0.5D0)*DC(I,J,M)
      Q2(2,1) = 2*3**(0.5D0)*DA(J,K,M)
      Q2(2,2) = - VT - 6*ABS(DA(J,K,M))
      Q2(2,4) = -2*3**(0.5D0)*DB(I,K,M)
      Q2(2,6) = -2*3**(0.5D0)*DC(I,J,M)
      Q2(3,1) = 2*3**(0.5D0)*DB(I,K,M)
      Q2(3,3) = - VT - 6*ABS(DB(I,K,M))
      Q2(3,4) = -2*3**(0.5D0)*DA(J,K,M)
      Q2(3,7) = -2*3**(0.5D0)*DC(I,J,M)
      Q2(4,2) = -2*3**(0.5D0)*DB(I,K,M)
      Q2(4,3) = -2*3**(0.5D0)*DA(J,K,M)
      Q2(4,4) = VT + 6*ABS(DB(I,K,M)) + 6*ABS(DA(J,K,M))
      Q2(4,8) = 2*3**(0.5D0)*DC(I,J,M)
      Q2(5,1) = 2*3**(0.5D0)*DC(I,J,M)
      Q2(5,5) = - VT - 6*ABS(DC(I,J,M))
      Q2(5,6) = -2*3**(0.5D0)*DA(J,K,M)
      Q2(5,7) = -2*3**(0.5D0)*DB(I,K,M)
      Q2(6,2) = -2*3**(0.5D0)*DC(I,J,M)
      Q2(6,5) = -2*3**(0.5D0)*DA(J,K,M)
      Q2(6,6) = VT + 6*ABS(DA(J,K,M)) + 6*ABS(DC(I,J,M))
      Q2(6,8) = 2*3**(0.5D0)*DB(I,K,M)
      Q2(7,3) = -2*3**(0.5D0)*DC(I,J,M)
      Q2(7,5) = -2*3**(0.5D0)*DB(I,K,M)
      Q2(7,7) = VT + 6*ABS(DB(I,K,M)) + 6*ABS(DC(I,J,M))
      Q2(7,8) = 2*3**(0.5D0)*DA(J,K,M)
      Q2(8,4) = 2*3**(0.5D0)*DC(I,J,M)
      Q2(8,6) = 2*3**(0.5D0)*DB(I,K,M)
      Q2(8,7) = 2*3**(0.5D0)*DA(J,K,M)
      Q2(8,8) = - VT - 6*ABS(DB(I,K,M)) - 6*ABS(DA(J,K,M)) -
     > 6*ABS(DC(I,J,M))
      Q2(1,9) = Q(1)*VOL(I,J,K)
      Q2(2,9) = -Q(2)*VOL(I,J,K)
      Q2(3,9) = -Q(3)*VOL(I,J,K)
      Q2(4,9) = Q(4)*VOL(I,J,K)
      Q2(5,9) = -Q(5)*VOL(I,J,K)
      Q2(6,9) = Q(6)*VOL(I,J,K)
      Q2(7,9) = Q(7)*VOL(I,J,K)
      Q2(8,9) = -Q(8)*VOL(I,J,K)
      Q2(1,9) = Q2(1,9) + 0
      Q2(2,9) = Q2(2,9) + 2*3**(0.5D0)*DA(J,K,M)*XNI(1,1,J,K)
      Q2(3,9) = Q2(3,9) + 2*3**(0.5D0)*DB(I,K,M)*XNJ(1,1,K)
      Q2(4,9) = Q2(4,9) - 2*3**(0.5D0)*DB(I,K,M)*XNJ(2,1,K) -
     > 2*3**(0.5D0)*DA(J,K,M)*XNI(2,1,J,K)
      Q2(5,9) = Q2(5,9) + 2*3**(0.5D0)*XNK(1,1)*DC(I,J,M)
      Q2(6,9) = Q2(6,9) - 2*3**(0.5D0)*DA(J,K,M)*XNI(1,2,J,K) -
     > 2*3**(0.5D0)*XNK(2,1)*DC(I,J,M)
      Q2(7,9) = Q2(7,9) - 2*3**(0.5D0)*DB(I,K,M)*XNJ(1,2,K) -
     >  2*3**(0.5D0)*XNK(1,2)*DC(I,J,M)
      Q2(8,9) = Q2(8,9) + 2*3**(0.5D0)*DB(I,K,M)*XNJ(2,2,K) +
     >  2*3**(0.5D0)*DA(J,K,M)*XNI(2,2,J,K) + 
     >  2*3**(0.5D0)*XNK(2,2)*DC(I,J,M)
      ELSE IF(IELEM.EQ.3) THEN
      Q2(1,1) = - VT - 2*ABS(DB(I,K,M)) - 2*ABS(DA(J,K,M)) -
     > 2*ABS(DC(I,J,M))
      Q2(1,3) = -2*5**(0.5D0)*ABS(DA(J,K,M))
      Q2(1,7) = -2*5**(0.5D0)*ABS(DB(I,K,M))
      Q2(1,19) = -2*5**(0.5D0)*ABS(DC(I,J,M))
      Q2(2,2) = - VT - 2*ABS(DB(I,K,M)) - 2*ABS(DC(I,J,M))
      Q2(2,3) = -2*15**(0.5D0)*DA(J,K,M)
      Q2(2,8) = -2*5**(0.5D0)*ABS(DB(I,K,M))
      Q2(2,20) = -2*5**(0.5D0)*ABS(DC(I,J,M))
      Q2(3,1) = 2*5**(0.5D0)*ABS(DA(J,K,M))
      Q2(3,2) = -2*15**(0.5D0)*DA(J,K,M)
      Q2(3,3) = VT + 2*ABS(DB(I,K,M)) + 10*ABS(DA(J,K,M)) +
     > 2*ABS(DC(I,J,M))
      Q2(3,9) = 2*5**(0.5D0)*ABS(DB(I,K,M))
      Q2(3,21) = 2*5**(0.5D0)*ABS(DC(I,J,M))
      Q2(4,4) = - VT - 2*ABS(DA(J,K,M)) - 2*ABS(DC(I,J,M))
      Q2(4,6) = -2*5**(0.5D0)*ABS(DA(J,K,M))
      Q2(4,7) = -2*15**(0.5D0)*DB(I,K,M)
      Q2(4,22) = -2*5**(0.5D0)*ABS(DC(I,J,M))
      Q2(5,5) = - VT - 2*ABS(DC(I,J,M))
      Q2(5,6) = -2*15**(0.5D0)*DA(J,K,M)
      Q2(5,8) = -2*15**(0.5D0)*DB(I,K,M)
      Q2(5,23) = -2*5**(0.5D0)*ABS(DC(I,J,M))
      Q2(6,4) = 2*5**(0.5D0)*ABS(DA(J,K,M))
      Q2(6,5) = -2*15**(0.5D0)*DA(J,K,M)
      Q2(6,6) = VT + 10*ABS(DA(J,K,M)) + 2*ABS(DC(I,J,M))
      Q2(6,9) = 2*15**(0.5D0)*DB(I,K,M)
      Q2(6,24) = 2*5**(0.5D0)*ABS(DC(I,J,M))
      Q2(7,1) = 2*5**(0.5D0)*ABS(DB(I,K,M))
      Q2(7,4) = -2*15**(0.5D0)*DB(I,K,M)
      Q2(7,7) = VT + 10*ABS(DB(I,K,M)) + 2*ABS(DA(J,K,M)) +
     > 2*ABS(DC(I,J,M))
      Q2(7,9) = 2*5**(0.5D0)*ABS(DA(J,K,M))
      Q2(7,25) = 2*5**(0.5D0)*ABS(DC(I,J,M))
      Q2(8,2) = 2*5**(0.5D0)*ABS(DB(I,K,M))
      Q2(8,5) = -2*15**(0.5D0)*DB(I,K,M)
      Q2(8,8) = VT + 10*ABS(DB(I,K,M)) + 2*ABS(DC(I,J,M))
      Q2(8,9) = 2*15**(0.5D0)*DA(J,K,M)
      Q2(8,26) = 2*5**(0.5D0)*ABS(DC(I,J,M))
      Q2(9,3) = -2*5**(0.5D0)*ABS(DB(I,K,M))
      Q2(9,6) = 2*15**(0.5D0)*DB(I,K,M)
      Q2(9,7) = -2*5**(0.5D0)*ABS(DA(J,K,M))
      Q2(9,8) = 2*15**(0.5D0)*DA(J,K,M)
      Q2(9,9) = - VT - 10*ABS(DB(I,K,M)) - 10*ABS(DA(J,K,M)) -
     > 2*ABS(DC(I,J,M))
      Q2(9,27) = -2*5**(0.5D0)*ABS(DC(I,J,M))
      Q2(10,10) = VT + 2*ABS(DB(I,K,M)) + 2*ABS(DA(J,K,M))
      Q2(10,12) = 2*5**(0.5D0)*ABS(DA(J,K,M))
      Q2(10,16) = 2*5**(0.5D0)*ABS(DB(I,K,M))
      Q2(10,19) = 2*15**(0.5D0)*DC(I,J,M)
      Q2(11,11) = VT + 2*ABS(DB(I,K,M))
      Q2(11,12) = 2*15**(0.5D0)*DA(J,K,M)
      Q2(11,17) = 2*5**(0.5D0)*ABS(DB(I,K,M))
      Q2(11,20) = 2*15**(0.5D0)*DC(I,J,M)
      Q2(12,10) = -2*5**(0.5D0)*ABS(DA(J,K,M))
      Q2(12,11) = 2*15**(0.5D0)*DA(J,K,M)
      Q2(12,12) = - VT - 2*ABS(DB(I,K,M)) - 10*ABS(DA(J,K,M))
      Q2(12,18) = -2*5**(0.5D0)*ABS(DB(I,K,M))
      Q2(12,21) = -2*15**(0.5D0)*DC(I,J,M)
      Q2(13,13) = VT + 2*ABS(DA(J,K,M))
      Q2(13,15) = 2*5**(0.5D0)*ABS(DA(J,K,M))
      Q2(13,16) = 2*15**(0.5D0)*DB(I,K,M)
      Q2(13,22) = 2*15**(0.5D0)*DC(I,J,M)
      Q2(14,14) = VT
      Q2(14,15) = 2*15**(0.5D0)*DA(J,K,M)
      Q2(14,17) = 2*15**(0.5D0)*DB(I,K,M)
      Q2(14,23) = 2*15**(0.5D0)*DC(I,J,M)
      Q2(15,13) = -2*5**(0.5D0)*ABS(DA(J,K,M))
      Q2(15,14) = 2*15**(0.5D0)*DA(J,K,M)
      Q2(15,15) = - VT - 10*ABS(DA(J,K,M))
      Q2(15,18) = -2*15**(0.5D0)*DB(I,K,M)
      Q2(15,24) = -2*15**(0.5D0)*DC(I,J,M)
      Q2(16,10) = -2*5**(0.5D0)*ABS(DB(I,K,M))
      Q2(16,13) = 2*15**(0.5D0)*DB(I,K,M)
      Q2(16,16) = - VT - 10*ABS(DB(I,K,M)) - 2*ABS(DA(J,K,M))
      Q2(16,18) = -2*5**(0.5D0)*ABS(DA(J,K,M))
      Q2(16,25) = -2*15**(0.5D0)*DC(I,J,M)
      Q2(17,11) = -2*5**(0.5D0)*ABS(DB(I,K,M))
      Q2(17,14) = 2*15**(0.5D0)*DB(I,K,M)
      Q2(17,17) = - VT - 10*ABS(DB(I,K,M))
      Q2(17,18) = -2*15**(0.5D0)*DA(J,K,M)
      Q2(17,26) = -2*15**(0.5D0)*DC(I,J,M)
      Q2(18,12) = 2*5**(0.5D0)*ABS(DB(I,K,M))
      Q2(18,15) = -2*15**(0.5D0)*DB(I,K,M)
      Q2(18,16) = 2*5**(0.5D0)*ABS(DA(J,K,M))
      Q2(18,17) = -2*15**(0.5D0)*DA(J,K,M)
      Q2(18,18) = VT + 10*ABS(DB(I,K,M)) + 10*ABS(DA(J,K,M))
      Q2(18,27) = 2*15**(0.5D0)*DC(I,J,M)
      Q2(19,1) = 2*5**(0.5D0)*ABS(DC(I,J,M))
      Q2(19,10) = -2*15**(0.5D0)*DC(I,J,M)
      Q2(19,19) = VT + 2*ABS(DB(I,K,M)) + 2*ABS(DA(J,K,M)) +
     > 10*ABS(DC(I,J,M))
      Q2(19,21) = 2*5**(0.5D0)*ABS(DA(J,K,M))
      Q2(19,25) = 2*5**(0.5D0)*ABS(DB(I,K,M))
      Q2(20,2) = 2*5**(0.5D0)*ABS(DC(I,J,M))
      Q2(20,11) = -2*15**(0.5D0)*DC(I,J,M)
      Q2(20,20) = VT + 2*ABS(DB(I,K,M)) + 10*ABS(DC(I,J,M))
      Q2(20,21) = 2*15**(0.5D0)*DA(J,K,M)
      Q2(20,26) = 2*5**(0.5D0)*ABS(DB(I,K,M))
      Q2(21,3) = -2*5**(0.5D0)*ABS(DC(I,J,M))
      Q2(21,12) = 2*15**(0.5D0)*DC(I,J,M)
      Q2(21,19) = -2*5**(0.5D0)*ABS(DA(J,K,M))
      Q2(21,20) = 2*15**(0.5D0)*DA(J,K,M)
      Q2(21,21) = - VT - 2*ABS(DB(I,K,M)) - 10*ABS(DA(J,K,M)) -
     > 10*ABS(DC(I,J,M))
      Q2(21,27) = -2*5**(0.5D0)*ABS(DB(I,K,M))
      Q2(22,4) = -2*5**(0.5D0)*ABS(DC(I,J,M))
      Q2(22,13) = 2*15**(0.5D0)*DC(I,J,M)
      Q2(22,22) = - VT - 2*ABS(DA(J,K,M)) - 10*ABS(DC(I,J,M))
      Q2(22,24) = -2*5**(0.5D0)*ABS(DA(J,K,M))
      Q2(22,25) = -2*15**(0.5D0)*DB(I,K,M)
      Q2(23,5) = -2*5**(0.5D0)*ABS(DC(I,J,M))
      Q2(23,14) = 2*15**(0.5D0)*DC(I,J,M)
      Q2(23,23) = - VT - 10*ABS(DC(I,J,M))
      Q2(23,24) = -2*15**(0.5D0)*DA(J,K,M)
      Q2(23,26) = -2*15**(0.5D0)*DB(I,K,M)
      Q2(24,6) = 2*5**(0.5D0)*ABS(DC(I,J,M))
      Q2(24,15) = -2*15**(0.5D0)*DC(I,J,M)
      Q2(24,22) = 2*5**(0.5D0)*ABS(DA(J,K,M))
      Q2(24,23) = -2*15**(0.5D0)*DA(J,K,M)
      Q2(24,24) = VT + 10*ABS(DA(J,K,M)) + 10*ABS(DC(I,J,M))
      Q2(24,27) = 2*15**(0.5D0)*DB(I,K,M)
      Q2(25,7) = -2*5**(0.5D0)*ABS(DC(I,J,M))
      Q2(25,16) = 2*15**(0.5D0)*DC(I,J,M)
      Q2(25,19) = -2*5**(0.5D0)*ABS(DB(I,K,M))
      Q2(25,22) = 2*15**(0.5D0)*DB(I,K,M)
      Q2(25,25) = - VT - 10*ABS(DB(I,K,M)) - 2*ABS(DA(J,K,M)) -
     > 10*ABS(DC(I,J,M))
      Q2(25,27) = -2*5**(0.5D0)*ABS(DA(J,K,M))
      Q2(26,8) = 2*5**(0.5D0)*ABS(DC(I,J,M))
      Q2(26,17) = -2*15**(0.5D0)*DC(I,J,M)
      Q2(26,20) = 2*5**(0.5D0)*ABS(DB(I,K,M))
      Q2(26,23) = -2*15**(0.5D0)*DB(I,K,M)
      Q2(26,26) = VT + 10*ABS(DB(I,K,M)) + 10*ABS(DC(I,J,M))
      Q2(26,27) = 2*15**(0.5D0)*DA(J,K,M)
      Q2(27,9) = -2*5**(0.5D0)*ABS(DC(I,J,M))
      Q2(27,18) = 2*15**(0.5D0)*DC(I,J,M)
      Q2(27,21) = -2*5**(0.5D0)*ABS(DB(I,K,M))
      Q2(27,24) = 2*15**(0.5D0)*DB(I,K,M)
      Q2(27,25) = -2*5**(0.5D0)*ABS(DA(J,K,M))
      Q2(27,26) = 2*15**(0.5D0)*DA(J,K,M)
      Q2(27,27) = - VT - 10*ABS(DB(I,K,M)) - 10*ABS(DA(J,K,M)) -
     > 10*ABS(DC(I,J,M))
      Q2(1,28) = -Q(01)*VOL(I,J,K)
      Q2(2,28) = -Q(02)*VOL(I,J,K)
      Q2(3,28) = Q(03)*VOL(I,J,K)
      Q2(4,28) = -Q(04)*VOL(I,J,K)
      Q2(5,28) = -Q(05)*VOL(I,J,K)
      Q2(6,28) = Q(06)*VOL(I,J,K)
      Q2(7,28) = Q(07)*VOL(I,J,K)
      Q2(8,28) = Q(08)*VOL(I,J,K)
      Q2(9,28) = -Q(09)*VOL(I,J,K)
      Q2(10,28) = Q(10)*VOL(I,J,K)
      Q2(11,28) = Q(11)*VOL(I,J,K)
      Q2(12,28) = -Q(12)*VOL(I,J,K)
      Q2(13,28) = Q(13)*VOL(I,J,K)
      Q2(14,28) = Q(14)*VOL(I,J,K)
      Q2(15,28) = -Q(15)*VOL(I,J,K)
      Q2(16,28) = -Q(16)*VOL(I,J,K)
      Q2(17,28) = -Q(17)*VOL(I,J,K)
      Q2(18,28) = Q(18)*VOL(I,J,K)
      Q2(19,28) = Q(19)*VOL(I,J,K)
      Q2(20,28) = Q(20)*VOL(I,J,K)
      Q2(21,28) = -Q(21)*VOL(I,J,K)
      Q2(22,28) = -Q(22)*VOL(I,J,K)
      Q2(23,28) = -Q(23)*VOL(I,J,K)
      Q2(24,28) = Q(24)*VOL(I,J,K)
      Q2(25,28) = -Q(25)*VOL(I,J,K)
      Q2(26,28) = Q(26)*VOL(I,J,K)
      Q2(27,28) = -Q(27)*VOL(I,J,K)
      Q2(1,28) = Q2(1,28) - 2*XNK(1,1)*ABS(DC(I,J,M)) -
     > 2*ABS(DB(I,K,M))*XNJ(1,1,K) - 2*ABS(DA(J,K,M))*XNI(1,1,J,K)
      Q2(2,28) = Q2(2,28) - 2*XNK(2,1)*ABS(DC(I,J,M)) -
     > 2*ABS(DB(I,K,M))*XNJ(2,1,K)
      Q2(3,28) = Q2(3,28) + 2*XNK(3,1)*ABS(DC(I,J,M)) +
     > 2*ABS(DB(I,K,M))*XNJ(3,1,K) + 
     > 2*5**(0.5D0)*ABS(DA(J,K,M))*XNI(1,1,J,K)
      Q2(4,28) = Q2(4,28) - 2*XNK(1,2)*ABS(DC(I,J,M)) -
     > 2*ABS(DA(J,K,M))*XNI(2,1,J,K)
      Q2(5,28) = Q2(5,28) -2*XNK(2,2)*ABS(DC(I,J,M))
      Q2(6,28) = Q2(6,28) + 2*XNK(3,2)*ABS(DC(I,J,M)) +
     > 2*5**(0.5D0)*ABS(DA(J,K,M))*XNI(2,1,J,K)
      Q2(7,28) = Q2(7,28) + 2*XNK(1,3)*ABS(DC(I,J,M)) +
     > 2*ABS(DA(J,K,M))*XNI(3,1,J,K) +
     > 2*5**(0.5D0)*ABS(DB(I,K,M))*XNJ(1,1,K)
      Q2(8,28) = Q2(8,28) + 2*XNK(2,3)*ABS(DC(I,J,M)) +
     > 2*5**(0.5D0)*ABS(DB(I,K,M))*XNJ(2,1,K)
      Q2(9,28) = Q2(9,28) - 2*XNK(3,3)*ABS(DC(I,J,M)) -
     > 2*5**(0.5D0)*ABS(DB(I,K,M))*XNJ(3,1,K) -
     > 2*5**(0.5D0)*ABS(DA(J,K,M))*XNI(3,1,J,K)
      Q2(10,28) = Q2(10,28) + 2*ABS(DB(I,K,M))*XNJ(1,2,K) +
     > 2*ABS(DA(J,K,M))*XNI(1,2,J,K)
      Q2(11,28) = Q2(11,28) + 2*ABS(DB(I,K,M))*XNJ(2,2,K)
      Q2(12,28) = Q2(12,28) - 2*ABS(DB(I,K,M))*XNJ(3,2,K) -
     > 2*5**(0.5D0)*ABS(DA(J,K,M))*XNI(1,2,J,K)
      Q2(13,28) = Q2(13,28) + 2*ABS(DA(J,K,M))*XNI(2,2,J,K)
      Q2(14,28) = Q2(14,28) + 0
      Q2(15,28) = Q2(15,28) -2*5**(0.5D0)*ABS(DA(J,K,M))*XNI(2,2,J,K)
      Q2(16,28) = Q2(16,28) - 2*ABS(DA(J,K,M))*XNI(3,2,J,K) -
     > 2*5**(0.5D0)*ABS(DB(I,K,M))*XNJ(1,2,K)
      Q2(17,28) = Q2(17,28) -2*5**(0.5D0)*ABS(DB(I,K,M))*XNJ(2,2,K)
      Q2(18,28) = Q2(18,28) + 2*5**(0.5D0)*ABS(DB(I,K,M))*XNJ(3,2,K) +
     > 2*5**(0.5D0)*ABS(DA(J,K,M))*XNI(3,2,J,K)
      Q2(19,28) = Q2(19,28) + 2*ABS(DB(I,K,M))*XNJ(1,3,K) +
     > 2*ABS(DA(J,K,M))*XNI(1,3,J,K) +
     > 2*5**(0.5D0)*XNK(1,1)*ABS(DC(I,J,M))
      Q2(20,28) = Q2(20,28) + 2*ABS(DB(I,K,M))*XNJ(2,3,K) +
     > 2*5**(0.5D0)*XNK(2,1)*ABS(DC(I,J,M))
      Q2(21,28) = Q2(21,28) - 2*ABS(DB(I,K,M))*XNJ(3,3,K) -
     > 2*5**(0.5D0)*ABS(DA(J,K,M))*XNI(1,3,J,K) -
     > 2*5**(0.5D0)*XNK(3,1)*ABS(DC(I,J,M))
      Q2(22,28) = Q2(22,28) - 2*ABS(DA(J,K,M))*XNI(2,3,J,K) -
     > 2*5**(0.5D0)*XNK(1,2)*ABS(DC(I,J,M))
      Q2(23,28) = Q2(23,28) - 2*5**(0.5D0)*XNK(2,2)*ABS(DC(I,J,M))
      Q2(24,28) = Q2(24,28) + 2*5**(0.5D0)*ABS(DA(J,K,M))*XNI(2,3,J,K) +
     > 2*5**(0.5D0)*XNK(3,2)*ABS(DC(I,J,M))
      Q2(25,28) = Q2(25,28) - 2*ABS(DA(J,K,M))*XNI(3,3,J,K) -
     > 2*5**(0.5D0)*ABS(DB(I,K,M))*XNJ(1,3,K) -
     > 2*5**(0.5D0)*XNK(1,3)*ABS(DC(I,J,M))
      Q2(26,28) = Q2(26,28) + 2*5**(0.5D0)*ABS(DB(I,K,M))*XNJ(2,3,K) +
     > 2*5**(0.5D0)*XNK(2,3)*ABS(DC(I,J,M))
      Q2(27,28) = Q2(27,28) - 2*5**(0.5D0)*ABS(DB(I,K,M))*XNJ(3,3,K) -
     > 2*5**(0.5D0)*ABS(DA(J,K,M))*XNI(3,3,J,K) -
     > 2*5**(0.5D0)*XNK(3,3)*ABS(DC(I,J,M))
      ENDIF
*
      CALL ALSBD(IELEM**3,1,Q2,IER,IELEM**3)
      IF(IER.NE.0) CALL XABORT('SNFT13: SINGULAR MATRIX.')
*
      IF(IELEM.EQ.1) THEN
      IF(LFIXUP.AND.(Q2(1,2).LE.RLOG)) Q2(1,2)=0.0
      XNI(1,1,J,K)=2.0D0*Q2(1,2)-XNI(1,1,J,K)
      XNJ(1,1,K)=2.0D0*Q2(1,2)-XNJ(1,1,K)
      XNK(1,1)=2.0D0*Q2(1,2)-XNK(1,1)
      IF(LFIXUP.AND.(XNI(1,1,J,K).LE.RLOG)) XNI(1,1,J,K)=0.0
      IF(LFIXUP.AND.(XNJ(1,1,K).LE.RLOG)) XNJ(1,1,K)=0.0
      IF(LFIXUP.AND.(XNK(1,1).LE.RLOG)) XNK(1,1)=0.0
      ELSE IF(IELEM.EQ.2) THEN
      XNI(1,1,J,K) = XNI(1,1,J,K) + SIGN(1.0,DU(M))*2*3**(0.5D0)*Q2(2,9)
      XNI(2,1,J,K) = XNI(2,1,J,K) + SIGN(1.0,DU(M))*2*3**(0.5D0)*Q2(4,9)
      XNI(1,2,J,K) = XNI(1,2,J,K) + SIGN(1.0,DU(M))*2*3**(0.5D0)*Q2(6,9)
      XNI(2,2,J,K) = XNI(2,2,J,K) + SIGN(1.0,DU(M))*2*3**(0.5D0)*Q2(8,9)
      XNJ(1,1,K) = XNJ(1,1,K) + SIGN(1.0,DE(M))*2*3**(0.5D0)*Q2(3,9)
      XNJ(2,1,K) = XNJ(2,1,K) + SIGN(1.0,DE(M))*2*3**(0.5D0)*Q2(4,9)
      XNJ(1,2,K) = XNJ(1,2,K) + SIGN(1.0,DE(M))*2*3**(0.5D0)*Q2(7,9)
      XNJ(2,2,K) = XNJ(2,2,K) + SIGN(1.0,DE(M))*2*3**(0.5D0)*Q2(8,9)
      XNK(1,1) = XNK(1,1) + SIGN(1.0,DZ(M))*2*3**(0.5D0)*Q2(5,9)
      XNK(2,1) = XNK(2,1) + SIGN(1.0,DZ(M))*2*3**(0.5D0)*Q2(6,9)
      XNK(1,2) = XNK(1,2) + SIGN(1.0,DZ(M))*2*3**(0.5D0)*Q2(7,9)
      XNK(2,2) = XNK(2,2) + SIGN(1.0,DZ(M))*2*3**(0.5D0)*Q2(8,9)
      ELSE IF(IELEM.EQ.3) THEN
      XNI(1,1,J,K) = 2*Q2(01,28) - XNI(1,1,J,K) + 2*5**(0.5D0)*Q2(03,28)
      XNI(2,1,J,K) = 2*Q2(04,28) - XNI(2,1,J,K) + 2*5**(0.5D0)*Q2(06,28)
      XNI(3,1,J,K) = 2*Q2(07,28) - XNI(3,1,J,K) + 2*5**(0.5D0)*Q2(09,28)
      XNI(1,2,J,K) = 2*Q2(10,28) - XNI(1,2,J,K) + 2*5**(0.5D0)*Q2(12,28)
      XNI(2,2,J,K) = 2*Q2(13,28) - XNI(2,2,J,K) + 2*5**(0.5D0)*Q2(15,28)
      XNI(3,2,J,K) = 2*Q2(16,28) - XNI(3,2,J,K) + 2*5**(0.5D0)*Q2(18,28)
      XNI(1,3,J,K) = 2*Q2(19,28) - XNI(1,3,J,K) + 2*5**(0.5D0)*Q2(21,28)
      XNI(2,3,J,K) = 2*Q2(22,28) - XNI(2,3,J,K) + 2*5**(0.5D0)*Q2(24,28)
      XNI(3,3,J,K) = 2*Q2(25,28) - XNI(3,3,J,K) + 2*5**(0.5D0)*Q2(27,28)
      XNJ(1,1,K) = 2*Q2(01,28) - XNJ(1,1,K) + 2*5**(0.5D0)*Q2(07,28)
      XNJ(2,1,K) = 2*Q2(02,28) - XNJ(2,1,K) + 2*5**(0.5D0)*Q2(08,28)
      XNJ(3,1,K) = 2*Q2(03,28) - XNJ(3,1,K) + 2*5**(0.5D0)*Q2(09,28)
      XNJ(1,2,K) = 2*Q2(10,28) - XNJ(1,2,K) + 2*5**(0.5D0)*Q2(16,28)
      XNJ(2,2,K) = 2*Q2(11,28) - XNJ(2,2,K) + 2*5**(0.5D0)*Q2(17,28)
      XNJ(3,2,K) = 2*Q2(12,28) - XNJ(3,2,K) + 2*5**(0.5D0)*Q2(18,28)
      XNJ(1,3,K) = 2*Q2(19,28) - XNJ(1,3,K) + 2*5**(0.5D0)*Q2(25,28)
      XNJ(2,3,K) = 2*Q2(20,28) - XNJ(2,3,K) + 2*5**(0.5D0)*Q2(26,28)
      XNJ(3,3,K) = 2*Q2(21,28) - XNJ(3,3,K) + 2*5**(0.5D0)*Q2(27,28)
      XNK(1,1) = 2*Q2(01,28) - XNK(1,1) + 2*5**(0.5D0)*Q2(19,28)
      XNK(2,1) = 2*Q2(02,28) - XNK(2,1) + 2*5**(0.5D0)*Q2(20,28)
      XNK(3,1) = 2*Q2(03,28) - XNK(3,1) + 2*5**(0.5D0)*Q2(21,28)
      XNK(1,2) = 2*Q2(04,28) - XNK(1,2) + 2*5**(0.5D0)*Q2(22,28)
      XNK(2,2) = 2*Q2(05,28) - XNK(2,2) + 2*5**(0.5D0)*Q2(23,28)
      XNK(3,2) = 2*Q2(06,28) - XNK(3,2) + 2*5**(0.5D0)*Q2(24,28)
      XNK(1,3) = 2*Q2(07,28) - XNK(1,3) + 2*5**(0.5D0)*Q2(25,28)
      XNK(2,3) = 2*Q2(08,28) - XNK(2,3) + 2*5**(0.5D0)*Q2(26,28)
      XNK(3,3) = 2*Q2(09,28) - XNK(3,3) + 2*5**(0.5D0)*Q2(27,28)
      ENDIF
      DO 270 P=1,NSCT
      DO 260 IEL=1,IELEM**3
      FLUX(IEL,P,I,J,K)=FLUX(IEL,P,I,J,K)+Q2(IEL,IELEM**3+1)*DN(P,M)
  260 CONTINUE
  270 CONTINUE
*--
  280 CONTINUE
      DO 300 JEL=1,IELEM
      DO 290 IEL=1,IELEM
      IOF=((((M-1)*LY+(J-1))*LX+(I-1))*IELEM+(JEL-1))*IELEM+IEL
      FUNKNO(L4+IOFZ+IOF,IG)=REAL(XNK(IEL,JEL))
  290 CONTINUE
  300 CONTINUE
*--
  310 CONTINUE
      DO 340 K=1,LZ
      DO 330 JEL=1,IELEM
      DO 320 IEL=1,IELEM
      IOF=((((M-1)*LZ+(K-1))*LX+(I-1))*IELEM+(JEL-1))*IELEM+IEL
      FUNKNO(L4+IOFY+IOF,IG)=REAL(XNJ(IEL,JEL,K))
  320 CONTINUE
  330 CONTINUE
  340 CONTINUE
*--
  350 CONTINUE
      DO 390 K=1,LZ
      DO 380 J=1,LY
      DO 370 JEL=1,IELEM
      DO 360 IEL=1,IELEM
      IOF=((((M-1)*LZ+(K-1))*LY+(J-1))*IELEM+(JEL-1))*IELEM+IEL
      FUNKNO(L4+IOF,IG)=REAL(XNI(IEL,JEL,J,K))
  360 CONTINUE
  370 CONTINUE
  380 CONTINUE
  390 CONTINUE
      FLUX_G(:,:,:,:,:,IG)=FLUX_G(:,:,:,:,:,IG)+FLUX(:,:,:,:,:)
  400 CONTINUE
  410 CONTINUE
*$OMP END PARALLEL DO
  420 CONTINUE
      DO 430 IG=1,NGEFF
        IF(.NOT.INCONV(IG)) GO TO 430
        FUNKNO(:L4,IG)=
     1  RESHAPE(REAL(FLUX_G(:IELEM**3,:NSCT,:LX,:LY,:LZ,IG)), (/L4/))
  430 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(FLUX_G,FLUX,XNJ,XNI,INDANG)
      RETURN
  500 FORMAT(16H SNFT13: thread=,I8,12H --->(group=,I4,7H angle=,I4,1H))
      END
