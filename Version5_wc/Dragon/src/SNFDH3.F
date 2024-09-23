*DECK SNFDH3
      SUBROUTINE SNFDH3(NMZ,NUN,NGEFF,IMPX,INCONV,NGIND,NHEX,LZ,NBC,
     1 IELEM,ISPLH,SIDE,NMAT,NPQ,NSCT,MAT,VOL,NCODE,ZCODE,TOTAL,QEXT,
     2 DU,DE,DZ,W,MRMZ,DC,DB,DA,PL,CONNEC,CONFROM,TASKLST,TASKSINWAVE,
     3 CON_WHOAMI,CON_TOWHOM,FUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform one inner iteration for solving SN equations in 3D hexagonal
* geometry for the HODD method. VOID boundary conditions on sides, and,
* VOID or REFL boundary conditions on top and bottom. KBA algorithm and 
* angle parallelisation. 
*
*Copyright:
* Copyright (C) 2022 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. A. Calloo
*
*Parameters: input
* NMZ     number of macrocells along z-axis.
* NUN     total number of unknowns in vector FUNKNO.
* NGRP    number of energy groups.
* NGEFF   dimension of arrays INCONV and INGIND.
* IMPX    print flag (equal to zero for no print).
* INCONV  energy group convergence flag (set to .FALSE. if converged).
* INGIND  energy group index assign to 1:NGEFF arrays.
* NHEX    number of hexagons in X-Y plane.
* LZ      number of meshes along Z axis.
* IELEM   measure of order of the spatial approximation polynomial:
*         =1 constant (classical diamond scheme - default for HODD);
*         =2 linear (default for DG);
*         =3 parabolic;
*         =4 cubic (only for DG).
* SIDE    side of an hexagon.
* NMAT    number of material mixtures.
* ISPLH   splitting option for hexagons.
* NPQ     number of SN directions in 12 dodecants (including zero-weight
*         directions).
* NSCT    maximum number of spherical harmonics moments of the flux.
* MAT     material mixture index in each region.
* VOL     volumes of each region.
* TOTAL   macroscopic total cross sections.
* NCODE   boundary condition indices.
* ZCODE   albedos.
* QEXT    Legendre components of the fixed source.
* DU      first direction cosines ($\\mu$).
* DE      second direction cosines ($\\eta$).
* DZ      third direction cosines ($\\xi$).
* W       weights.
* MRMZ    quadrature index.
* DC      diamond-scheme parameter.
* DB      diamond-scheme parameter.
* DA      diamond-scheme parameter.
* PL      discrete values of the spherical harmonics corresponding
*         to the 2D SN quadrature.
* CONNEC  connectivity matrix for flux swapping -- which lozenges is the
*         lozenge under consideration connected to; in order to pass the
*         flux along. This is dependent on direction.
* CONFROM matrix for incoming flux -- which lozenges are feeding into
*         the lozenge under consideration. This is dependent on
*         direction.
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
      INTEGER NMZ,NUN,NGEFF,IMPX,NGIND(NGEFF),NHEX,LZ,IELEM,ISPLH,
     1 NMAT,NPQ,NSCT,MAT(ISPLH,ISPLH,3,NHEX,LZ),NBC,NCODE(6),MRMZ(NPQ),
     2 CONNEC(3,(NHEX*3)*2,6),CONFROM(2,3,6),
     3 TASKLST(2,NBC*NMZ,((NBC-1)*4 +NMZ),6),
     4 TASKSINWAVE(((NBC-1)*4 +NMZ))
      LOGICAL INCONV(NGEFF)
      REAL SIDE,VOL(ISPLH,ISPLH,3,NHEX,LZ),ZCODE(6),TOTAL(0:NMAT,NGEFF),
     1 QEXT(NUN,NGEFF),DU(NPQ),DE(NPQ),DZ(NPQ),W(NPQ),
     2 DC(ISPLH*ISPLH*3*NHEX,1,NPQ),DB(ISPLH*ISPLH*3*NHEX,LZ,NPQ),
     3 DA(1,LZ,NPQ),PL(NSCT,NPQ),FUNKNO(NUN,NGEFF)
*----
*  LOCAL VARIABLES
*----
      INTEGER NPQD(12),IIND(12)
      DOUBLE PRECISION Q(IELEM**3), Q2(IELEM**3,(IELEM**3)+1),
     1 VT, CONST0, CONST1, CONST2
      PARAMETER(IUNOUT=6,PI=3.141592654)
      INTEGER :: ILOZSWP(3,6), IFROMI, IFROMJ
      REAL :: JAC(2,2,3), MUH, ETAH, XIH, AAA, BBB, CCC, DDD, MUHTEMP,
     1 ETAHTEMP
      DOUBLE PRECISION :: THETA, XNI(IELEM,IELEM,ISPLH,LZ),
     1 XNJ(IELEM,IELEM,LZ),XNK(IELEM,IELEM),C1
      INTEGER MAXWAVECNT,NCOL
      INTEGER :: CON_WHOAMI(2,(NHEX*3)*2,6), CON_TOWHOM(2,(NHEX*3)*2,6)
      INTEGER CHEXx,CSIDEx,CHEXy,CSIDEy
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: INDANG
      DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:,:) :: FLUX
      DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: FLUX_G
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: BFLUX

      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: BFLUZ
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: TMPBZ

      REAL, ALLOCATABLE, DIMENSION(:,:) :: TMP,TMP2
*----
*  DEFINITION OF CONSTANTS.
*----
      L4=3*(IELEM**3)*(ISPLH**2)*NHEX*LZ*NSCT
      L6=NHEX*3*ISPLH*ISPLH*NPQ*IELEM*IELEM
      CONST0=2.0D0*DSQRT(3.0D0)
      CONST1=2.0D0*DSQRT(5.0D0)
      CONST2=2.0D0*DSQRT(15.0D0)
      MAXWAVECNT = ((NBC-1)*4 +NMZ)
      NCOL=2*NBC -1
      NCELZ=CEILING(REAL(LZ)/NMZ)
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(INDANG(NPQ,12))
      ALLOCATE(FLUX(IELEM**3,NSCT))
      ALLOCATE(FLUX_G(IELEM**3,NSCT,3*ISPLH**2,NHEX,LZ,NGEFF))
      ALLOCATE(BFLUX(3,NCOL,ISPLH,IELEM,IELEM,LZ))

      ALLOCATE(BFLUZ(3,NHEX,ISPLH,ISPLH,IELEM,IELEM))
      ALLOCATE(TMPBZ(3,NHEX,ISPLH,ISPLH,IELEM,IELEM))

      ALLOCATE(TMP(NHEX*3*ISPLH*ISPLH*NPQ*IELEM*IELEM,NGEFF))
      ALLOCATE(TMP2(NHEX*3*ISPLH*ISPLH*NPQ*IELEM*IELEM,NGEFF))

      FLUX_G(:IELEM**3,:NSCT,:3*ISPLH**2,:NHEX,:LZ,:NGEFF)=0.0D0
*----
*  CONSTRUCT JACOBIAN MATRIX FOR EACH LOZENGE
*----
      ILOZSWP = RESHAPE((/ 3, 2, 1, 3, 1, 2, 1, 3, 2, 1, 2, 3, 2, 1,
     1   3, 2, 3, 1 /), SHAPE(ILOZSWP))
      JAC = RESHAPE((/ 1., -SQRT(3.), 1., SQRT(3.), 2., 0., 1.,
     1    SQRT(3.), 2., 0., -1., SQRT(3.) /), SHAPE(JAC))
      JAC = (SIDE/2.)*JAC

*----
*  INITIALISE ARRAYS AND STORE PREVIOUS ITERATION OUTGOING FLUXES
*----
      TMP(:NHEX*3*ISPLH*ISPLH*NPQ*IELEM*IELEM,:NGEFF)=0.0
      TMP2(:NHEX*3*ISPLH*ISPLH*NPQ*IELEM*IELEM,:NGEFF)=0.0
      IF((NCODE(5).NE.1).or.(NCODE(6).NE.1))THEN
         DO 776 IG=1,NGEFF
            IF(.NOT.INCONV(IG)) GO TO 776
            TMP2(1:L6,IG)=FUNKNO(L4+1:L4+L6,IG)
  776    CONTINUE
      ENDIF

************************************************************************
************************************************************************
************************************************************************

*----
*  PARAMETER VALIDATION.
*----
      IF((IELEM.LT.1).OR.(IELEM.GT.3))
     1   CALL XABORT('SNFDH3: INVALID IELEM (DIAM) VALUE. '
     2   //'CHECK INPUT DATA FILE.')
*----
*  SET DODECANT SWAPPING ORDER.
*----
      NPQD(:12)=0
      INDANG(:NPQ,:12)=0
      IIND(:12)=0
      DO 10 M=1,NPQ
        VU=DU(M)
        VE=DE(M)
        VZ=DZ(M)
        IF(W(M).EQ.0) CYCLE
        THETA=0.0D0
        IF(VE.GT.0.0)THEN
          IF(VU.EQ.0.0)THEN
            THETA = PI/2
          ELSEIF(VU.GT.0.0)THEN
            THETA = ATAN(ABS(VE/VU))
          ELSEIF(VU.LT.0.0)THEN
            THETA = PI - ATAN(ABS(VE/VU))
          ENDIF
        ELSEIF(VE.LT.0.0)THEN
          IF(VU.EQ.0.0)THEN
            THETA = 3*PI/2
          ELSEIF(VU.LT.0.0)THEN
            THETA = PI + ATAN(ABS(VE/VU))
          ELSEIF(VU.GT.0.0)THEN
            THETA = 2.*PI - ATAN(ABS(VE/VU))
          ENDIF
        ENDIF
*---- UNFOLD DODECANTS
        IND=0
        IF(VZ.GE.0.0)THEN
          IF((THETA.GT.0.0).AND.(THETA.LT.(PI/3.)))THEN
            IND=1
          ELSEIF((THETA.GT.(PI/3.)).AND.(THETA.LT.(2.*PI/3.)))THEN
            IND=2
          ELSEIF((THETA.GT.(2.*PI/3.)).AND.(THETA.LT.(PI)))THEN
            IND=3
          ELSEIF((THETA.GT.(PI)).AND.(THETA.LT.(4.*PI/3.)))THEN
            IND=4
          ELSEIF((THETA.GT.(4.*PI/3.)).AND.(THETA.LT.(5.*PI/3.)))THEN
            IND=5
          ELSEIF((THETA.GT.(5.*PI/3.)).AND.(THETA.LT.(2.*PI)))THEN
            IND=6
          ENDIF
        ELSEIF(VZ.LT.0.0)THEN
          IF((THETA.GT.0.0).AND.(THETA.LT.(PI/3.)))THEN
            IND=7
          ELSEIF((THETA.GT.(PI/3.)).AND.(THETA.LT.(2.*PI/3.)))THEN
            IND=8
          ELSEIF((THETA.GT.(2.*PI/3.)).AND.(THETA.LT.(PI)))THEN
            IND=9
          ELSEIF((THETA.GT.(PI)).AND.(THETA.LT.(4.*PI/3.)))THEN
            IND=10
          ELSEIF((THETA.GT.(4.*PI/3.)).AND.(THETA.LT.(5.*PI/3.)))THEN
            IND=11
          ELSEIF((THETA.GT.(5.*PI/3.)).AND.(THETA.LT.(2.*PI)))THEN
            IND=12
          ENDIF
        ENDIF
        IIND(IND)=IND ! assume IIND(I)=I in hexagonal geometry
        NPQD(IND)=NPQD(IND)+1
        INDANG(NPQD(IND),IND)=M
   10 CONTINUE
*----
*  MAIN LOOP OVER DODECANTS.
*----
      DO 310 JND=1,12
      IND=IIND(JND)
      IND2=MOD(IND-1,6)+1
      IF(IND.EQ.0) GO TO 310 ! needed because of S2 LS (8 dir. for 12 dodecants)
*----
*  PRELIMINARY LOOPS FOR SETTING BOUNDARY CONDITIONS.
*----
      IF((NCODE(5).NE.1).or.(NCODE(6).NE.1))THEN
      DO 150 IG=1,NGEFF
      DO 140 IPQD=1,NPQD(IND)
      IF(.NOT.INCONV(IG)) GO TO 140
      M=INDANG(IPQD,IND)
      VZ=DZ(M)
*
      IF(VZ.GT.0.0)THEN
         M1=MRMZ(M)
         IF(NCODE(5).NE.4)THEN
         DO IEL=1,IELEM
         DO JEL=1,IELEM
           DO I=1,NHEX
           DO J=1,3
           DO IL=1,ISPLH
           DO JL=1,ISPLH
           IOF=((((((M-1)*NHEX+(I-1))*3+(J-1))*ISPLH+(JL-1))*ISPLH+
     >      (IL-1))*IELEM+(JEL-1))*IELEM + IEL
           JOF=((((((M1-1)*NHEX+(I-1))*3+(J-1))*ISPLH+(JL-1))*ISPLH+
     >      (IL-1))*IELEM+(JEL-1))*IELEM+IEL
           FUNKNO(L4+IOF,IG)=FUNKNO(L4+JOF,IG)
           TMP(IOF,IG)=TMP(JOF,IG)
           ENDDO
           ENDDO
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
           DO I=1,NHEX
           DO J=1,3
           DO IL=1,ISPLH
           DO JL=1,ISPLH
           IOF=((((((M-1)*NHEX+(I-1))*3+(J-1))*ISPLH+(JL-1))*ISPLH+
     >      (IL-1))*IELEM+(JEL-1))*IELEM+IEL
           JOF=((((((M1-1)*NHEX+(I-1))*3+(J-1))*ISPLH+(JL-1))*ISPLH+
     >      (IL-1))*IELEM+(JEL-1))*IELEM+IEL
           FUNKNO(L4+IOF,IG)=FUNKNO(L4+JOF,IG)
           TMP(IOF,IG)=TMP(JOF,IG)
           ENDDO
           ENDDO
           ENDDO
           ENDDO
         ENDDO
         ENDDO
         ENDIF
      ENDIF
*
  140 CONTINUE
  150 CONTINUE
      ENDIF

*----
*  KBA DIAGONAL LOOP
*----
      BFLUX(:3,:NCOL,:ISPLH,:IELEM,:IELEM,:LZ)=0.0D0
      BFLUZ(:3,:NHEX,:ISPLH,:ISPLH,:IELEM,:IELEM)=0.0D0

      NHEXPROCD = 0

      DO 305 IDI=1,MAXWAVECNT

      TMPBZ(:3,:NHEX,:ISPLH,:ISPLH,:IELEM,:IELEM)=0.0D0

      IF(IDI.NE.1) NHEXPROCD = NHEXPROCD + TASKSINWAVE(IDI-1)

* Loop over macrocells // hexes in same sweep front
      DO 304 ICEL=1,TASKSINWAVE(IDI)
      I=TASKLST(1,ICEL,IDI,IND2)
      KK=TASKLST(2,ICEL,IDI,IND2)
      KK2=TASKLST(2,ICEL,IDI,IND2)
*----
*  MAIN SWAPPING LOOPS FOR SN FLUX CALCULATION
*----
*$OMP  PARALLEL DO
      DO 300 IG=1,NGEFF
      DO 290 IPQD=1,NPQD(IND)
#if defined(_OPENMP)
        ITID=omp_get_thread_num()
#else
        ITID=0
#endif
      IF(IMPX.GT.5) WRITE(IUNOUT,500) ITID,NGIND(IG),IPQD
      IF(.NOT.INCONV(IG)) GO TO 290
      M=INDANG(IPQD,IND)
      IF(W(M).EQ.0.0) GO TO 290

* Loop over z-levels in z-macrocell
      DO 280 ICELZ=1,NCELZ

      K=((KK-1)*NCELZ)+ICELZ
      K2=K
      IF(IND.GE.7) K=LZ+1-K
      IF(K.GT.LZ) GO TO 280

* Loop over lozenges
      DO 270 JJ=1,3
      J=ILOZSWP(JJ,IND2)
*
      AAA = JAC(1,1,J)
      BBB = JAC(1,2,J)
      CCC = JAC(2,1,J)
      DDD = JAC(2,2,J)

* Identify DOMAIN-WISE indices of connected sides 
      IHEXI  = CONNEC(1,((I-1)*3*2) + ((J-1)*2) +1,IND2)
      ILOZI  = CONNEC(2,((I-1)*3*2) + ((J-1)*2) +1,IND2)
      IHEXJ  = CONNEC(1,((I-1)*3*2) + ((J-1)*2) +2,IND2)
      ILOZJ  = CONNEC(2,((I-1)*3*2) + ((J-1)*2) +2,IND2)
      IFROMI = CONFROM(1,J,IND2)
      IFROMJ = CONFROM(2,J,IND2)

* Identify BF-ARRAY-WISE indices of sides of current lozenge
      MHEXx  = CON_WHOAMI(1,((I-1)*3*2) + ((J-1)*2) +1,IND2)
      MSIDEx = CON_WHOAMI(2,((I-1)*3*2) + ((J-1)*2) +1,IND2)
      MHEXy  = CON_WHOAMI(1,((I-1)*3*2) + ((J-1)*2) +2,IND2)
      MSIDEy = CON_WHOAMI(2,((I-1)*3*2) + ((J-1)*2) +2,IND2)

* Identify BF-ARRAY-WISE indices of connected sides
      CHEXx  = CON_TOWHOM(1,((I-1)*3*2) + ((J-1)*2) +1,IND2)
      CSIDEx = CON_TOWHOM(2,((I-1)*3*2) + ((J-1)*2) +1,IND2)
      CHEXy  = CON_TOWHOM(1,((I-1)*3*2) + ((J-1)*2) +2,IND2)
      CSIDEy = CON_TOWHOM(2,((I-1)*3*2) + ((J-1)*2) +2,IND2)

* Loop over internal split of lozenges
      DO 260 IL=1,ISPLH
      I2=IL
      IF(CONFROM(1,J,IND2).EQ.3) I2=ISPLH+1-I2
      I4=I2

!     Needed here because otherwise I will overwrite the boundary flux 
!     I will need. Assuming split of 3 for each lozenge, BF 1 will go in 
!     the spot of BF 3, which I have not taken yet... 
!     TO BE CHECKED
      IF(J.EQ.1) I4=ISPLH+1-I2
      DO IEL=1,IELEM
      DO JEL=1,IELEM
         XNJ(IEL,JEL,K) = BFLUX(MSIDEy,MHEXy,I4,IEL,JEL,K)
         BFLUX(MSIDEy,MHEXy,I4,IEL,JEL,K) = 0.0D0
      ENDDO
      ENDDO

      DO 250 JL=1,ISPLH
      J2=JL
      IF(CONFROM(2,J,IND2).EQ.4) J2=ISPLH+1-J2
      J4=J2

      IF(IL.EQ.1)THEN
         DO IEL=1,IELEM
         DO JEL=1,IELEM
            XNI(IEL,JEL,J2,K) = BFLUX(MSIDEx,MHEXx,J4,IEL,JEL,K)
            BFLUX(MSIDEx,MHEXx,J4,IEL,JEL,K) = 0.0D0
         ENDDO
         ENDDO
      ENDIF
*-----
* Pick up z-axis boundary conditions.
      IF(K2.EQ.1)THEN
         DO IEL=1,IELEM
         DO JEL=1,IELEM
            IF((NCODE(5).NE.1).or.(NCODE(6).NE.1))THEN
            IOF=((((((M-1)*NHEX+(I-1))*3+(J-1))*ISPLH+(J2-1))*ISPLH+
     1         (I2-1))*IELEM+(JEL-1))*IELEM+IEL
            IF(IND.LT.7) THEN
               BFLUZ(J,I,I2,J2,IEL,JEL)=FUNKNO(L4+IOF,IG)*ZCODE(5)
            ELSE
               BFLUZ(J,I,I2,J2,IEL,JEL)=FUNKNO(L4+IOF,IG)*ZCODE(6)
            ENDIF
            ELSE
               BFLUZ(J,I,I2,J2,IEL,JEL)=0.0D0
            ENDIF
         ENDDO
         ENDDO
      ENDIF

      DO IEL=1,IELEM
      DO JEL=1,IELEM
         XNK(IEL,JEL) = BFLUZ(J,I,I2,J2,IEL,JEL)
         BFLUZ(J,I,I2,J2,IEL,JEL) = 0
      ENDDO
      ENDDO
*
      FLUX(:IELEM**3,:NSCT)=0.0D0
*
      MUHTEMP  =  DA(1,K,M)
      ETAHTEMP =  DB(1,K,M)
      MUH  = (MUHTEMP*DDD) - (ETAHTEMP*BBB)
      ETAH = (-MUHTEMP*CCC) + (ETAHTEMP*AAA)
      XIH  = DC(1,1,M)
*
      IF(MAT(I2,J2,J,I,K).EQ.0) GO TO 250

*     -----------------------------------------------------
      DO 170 IEL=1,IELEM**3
      Q(IEL)=0.0D0
      DO 160 ISCT=1,NSCT
      IOF=(((((((K-1)*NHEX+(I-1))*3+(J-1))*ISPLH+(J2-1))*ISPLH+
     1 (I2-1))*NSCT+(ISCT-1))*IELEM**3)+IEL
      Q(IEL)=Q(IEL)+QEXT(IOF,IG)*PL(ISCT,M)/(4.0D0*PI)
  160 CONTINUE
  170 CONTINUE
*     -----------------------------------------------------
*
      VT=VOL(I2,J2,J,I,K)*TOTAL(MAT(I2,J2,J,I,K),IG)
      CALL XDDSET(Q2,(IELEM**3)*((IELEM**3)+1),0.0D0)
*
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
      IF(IELEM.EQ.1) THEN
      Q2(1,1)=2.0D0*ABS(MUH)+2.0D0*ABS(ETAH)+
     1        2.0D0*ABS(XIH)+VT 
      Q2(1,2)=2.0D0*ABS(MUH)*XNI(1,1,J2,K)+
     1        2.0D0*ABS(ETAH)*XNJ(1,1,K)+
     2        2.0D0*ABS(XIH)*XNK(1,1)+VOL(I2,J2,J,I,K)*Q(1)
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
      ELSE IF(IELEM.EQ.2) THEN
      Q2(1,1) = VT 
      Q2(1,2) = 2*3**(0.5D0)*MUH 
      Q2(1,3) = 2*3**(0.5D0)*ETAH 
      Q2(1,5) = 2*3**(0.5D0)*XIH 
      Q2(2,1) = 2*3**(0.5D0)*MUH 
      Q2(2,2) = - VT - 6*ABS(MUH) 
      Q2(2,4) = -2*3**(0.5D0)*ETAH 
      Q2(2,6) = -2*3**(0.5D0)*XIH 
      Q2(3,1) = 2*3**(0.5D0)*ETAH 
      Q2(3,3) = - VT - 6*ABS(ETAH) 
      Q2(3,4) = -2*3**(0.5D0)*MUH 
      Q2(3,7) = -2*3**(0.5D0)*XIH 
      Q2(4,2) = -2*3**(0.5D0)*ETAH 
      Q2(4,3) = -2*3**(0.5D0)*MUH 
      Q2(4,4) = VT + 6*ABS(ETAH) + 6*ABS(MUH) 
      Q2(4,8) = 2*3**(0.5D0)*XIH 
      Q2(5,1) = 2*3**(0.5D0)*XIH 
      Q2(5,5) = - VT - 6*ABS(XIH) 
      Q2(5,6) = -2*3**(0.5D0)*MUH 
      Q2(5,7) = -2*3**(0.5D0)*ETAH 
      Q2(6,2) = -2*3**(0.5D0)*XIH 
      Q2(6,5) = -2*3**(0.5D0)*MUH 
      Q2(6,6) = VT + 6*ABS(MUH) + 6*ABS(XIH) 
      Q2(6,8) = 2*3**(0.5D0)*ETAH 
      Q2(7,3) = -2*3**(0.5D0)*XIH 
      Q2(7,5) = -2*3**(0.5D0)*ETAH 
      Q2(7,7) = VT + 6*ABS(ETAH) + 6*ABS(XIH) 
      Q2(7,8) = 2*3**(0.5D0)*MUH 
      Q2(8,4) = 2*3**(0.5D0)*XIH 
      Q2(8,6) = 2*3**(0.5D0)*ETAH 
      Q2(8,7) = 2*3**(0.5D0)*MUH 
      Q2(8,8) = - VT - 6*ABS(ETAH) - 6*ABS
     >   (MUH) - 6*ABS(XIH) 
      Q2(1,9) = Q(1)*VOL(I2,J2,J,I,K) 
      Q2(2,9) = -Q(2)*VOL(I2,J2,J,I,K) 
      Q2(3,9) = -Q(3)*VOL(I2,J2,J,I,K) 
      Q2(4,9) = Q(4)*VOL(I2,J2,J,I,K) 
      Q2(5,9) = -Q(5)*VOL(I2,J2,J,I,K) 
      Q2(6,9) = Q(6)*VOL(I2,J2,J,I,K) 
      Q2(7,9) = Q(7)*VOL(I2,J2,J,I,K) 
      Q2(8,9) = -Q(8)*VOL(I2,J2,J,I,K) 
      Q2(1,9) = Q2(1,9) + 0 
      Q2(2,9) = Q2(2,9) + 2*3**(0.5D0)*MUH
     >   *XNI(1,1,J2,K) 
      Q2(3,9) = Q2(3,9) + 2*3**(0.5D0)*ETAH
     >   *XNJ(1,1,K) 
      Q2(4,9) = Q2(4,9) - 2*3**(0.5D0)*ETAH
     >   *XNJ(2,1,K) - 2*3**(0.5D0)*MUH*X
     >   NI(2,1,J2,K) 
      Q2(5,9) = Q2(5,9) + 2*3**(0.5D0)*XNK(1,1)
     >   *XIH 
      Q2(6,9) = Q2(6,9) - 2*3**(0.5D0)*MUH
     >   *XNI(1,2,J2,K) - 2*3**(0.5D0)*XNK(2,1)*
     >   XIH 
      Q2(7,9) = Q2(7,9) - 2*3**(0.5D0)*ETAH
     >   *XNJ(1,2,K) - 2*3**(0.5D0)*XNK(1,2)*XIH
      Q2(8,9) = Q2(8,9) + 2*3**(0.5D0)*ETAH
     >   *XNJ(2,2,K) + 2*3**(0.5D0)*MUH*XNI
     >   (2,2,J2,K) + 2*3**(0.5D0)*XNK(2,2)*XIH
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
      ELSE IF(IELEM.EQ.3) THEN
      Q2(1,1) = - VT - 2*ABS(ETAH) - 2*ABS(MUH)
     >    - 2*ABS(XIH) 
      Q2(1,2) = 0 
      Q2(1,3) = -2*5**(0.5D0)*ABS(MUH) 
      Q2(1,4) = 0 
      Q2(1,5) = 0 
      Q2(1,6) = 0 
      Q2(1,7) = -2*5**(0.5D0)*ABS(ETAH) 
      Q2(1,8) = 0 
      Q2(1,9) = 0 
      Q2(1,10) = 0 
      Q2(1,11) = 0 
      Q2(1,12) = 0 
      Q2(1,13) = 0 
      Q2(1,14) = 0 
      Q2(1,15) = 0 
      Q2(1,16) = 0 
      Q2(1,17) = 0 
      Q2(1,18) = 0 
      Q2(1,19) = -2*5**(0.5D0)*ABS(XIH) 
      Q2(1,20) = 0 
      Q2(1,21) = 0 
      Q2(1,22) = 0 
      Q2(1,23) = 0 
      Q2(1,24) = 0 
      Q2(1,25) = 0 
      Q2(1,26) = 0 
      Q2(1,27) = 0 
      Q2(2,1) = 0 
      Q2(2,2) = - VT - 2*ABS(ETAH) - 2*ABS(XIH)
      Q2(2,3) = -2*15**(0.5D0)*MUH 
      Q2(2,4) = 0 
      Q2(2,5) = 0 
      Q2(2,6) = 0 
      Q2(2,7) = 0 
      Q2(2,8) = -2*5**(0.5D0)*ABS(ETAH) 
      Q2(2,9) = 0 
      Q2(2,10) = 0 
      Q2(2,11) = 0 
      Q2(2,12) = 0 
      Q2(2,13) = 0 
      Q2(2,14) = 0 
      Q2(2,15) = 0 
      Q2(2,16) = 0 
      Q2(2,17) = 0 
      Q2(2,18) = 0 
      Q2(2,19) = 0 
      Q2(2,20) = -2*5**(0.5D0)*ABS(XIH) 
      Q2(2,21) = 0 
      Q2(2,22) = 0 
      Q2(2,23) = 0 
      Q2(2,24) = 0 
      Q2(2,25) = 0 
      Q2(2,26) = 0 
      Q2(2,27) = 0 
      Q2(3,1) = 2*5**(0.5D0)*ABS(MUH) 
      Q2(3,2) = -2*15**(0.5D0)*MUH 
      Q2(3,3) = VT + 2*ABS(ETAH) + 10*ABS(MUH) 
     >   + 2*ABS(XIH) 
      Q2(3,4) = 0 
      Q2(3,5) = 0 
      Q2(3,6) = 0 
      Q2(3,7) = 0 
      Q2(3,8) = 0 
      Q2(3,9) = 2*5**(0.5D0)*ABS(ETAH) 
      Q2(3,10) = 0 
      Q2(3,11) = 0 
      Q2(3,12) = 0 
      Q2(3,13) = 0 
      Q2(3,14) = 0 
      Q2(3,15) = 0 
      Q2(3,16) = 0 
      Q2(3,17) = 0 
      Q2(3,18) = 0 
      Q2(3,19) = 0 
      Q2(3,20) = 0 
      Q2(3,21) = 2*5**(0.5D0)*ABS(XIH) 
      Q2(3,22) = 0 
      Q2(3,23) = 0 
      Q2(3,24) = 0 
      Q2(3,25) = 0 
      Q2(3,26) = 0 
      Q2(3,27) = 0 
      Q2(4,1) = 0 
      Q2(4,2) = 0 
      Q2(4,3) = 0 
      Q2(4,4) = - VT - 2*ABS(MUH) - 2*ABS(XIH) 
      Q2(4,5) = 0 
      Q2(4,6) = -2*5**(0.5D0)*ABS(MUH) 
      Q2(4,7) = -2*15**(0.5D0)*ETAH 
      Q2(4,8) = 0 
      Q2(4,9) = 0 
      Q2(4,10) = 0 
      Q2(4,11) = 0 
      Q2(4,12) = 0 
      Q2(4,13) = 0 
      Q2(4,14) = 0 
      Q2(4,15) = 0 
      Q2(4,16) = 0 
      Q2(4,17) = 0 
      Q2(4,18) = 0 
      Q2(4,19) = 0 
      Q2(4,20) = 0 
      Q2(4,21) = 0 
      Q2(4,22) = -2*5**(0.5D0)*ABS(XIH) 
      Q2(4,23) = 0 
      Q2(4,24) = 0 
      Q2(4,25) = 0 
      Q2(4,26) = 0 
      Q2(4,27) = 0 
      Q2(5,1) = 0 
      Q2(5,2) = 0 
      Q2(5,3) = 0 
      Q2(5,4) = 0 
      Q2(5,5) = - VT - 2*ABS(XIH) 
      Q2(5,6) = -2*15**(0.5D0)*MUH 
      Q2(5,7) = 0 
      Q2(5,8) = -2*15**(0.5D0)*ETAH 
      Q2(5,9) = 0 
      Q2(5,10) = 0 
      Q2(5,11) = 0 
      Q2(5,12) = 0 
      Q2(5,13) = 0 
      Q2(5,14) = 0 
      Q2(5,15) = 0 
      Q2(5,16) = 0 
      Q2(5,17) = 0 
      Q2(5,18) = 0 
      Q2(5,19) = 0 
      Q2(5,20) = 0 
      Q2(5,21) = 0 
      Q2(5,22) = 0 
      Q2(5,23) = -2*5**(0.5D0)*ABS(XIH) 
      Q2(5,24) = 0 
      Q2(5,25) = 0 
      Q2(5,26) = 0 
      Q2(5,27) = 0 
      Q2(6,1) = 0 
      Q2(6,2) = 0 
      Q2(6,3) = 0 
      Q2(6,4) = 2*5**(0.5D0)*ABS(MUH) 
      Q2(6,5) = -2*15**(0.5D0)*MUH 
      Q2(6,6) = VT + 10*ABS(MUH) + 2*ABS(XIH) 
      Q2(6,7) = 0 
      Q2(6,8) = 0 
      Q2(6,9) = 2*15**(0.5D0)*ETAH 
      Q2(6,10) = 0 
      Q2(6,11) = 0 
      Q2(6,12) = 0 
      Q2(6,13) = 0 
      Q2(6,14) = 0 
      Q2(6,15) = 0 
      Q2(6,16) = 0 
      Q2(6,17) = 0 
      Q2(6,18) = 0 
      Q2(6,19) = 0 
      Q2(6,20) = 0 
      Q2(6,21) = 0 
      Q2(6,22) = 0 
      Q2(6,23) = 0 
      Q2(6,24) = 2*5**(0.5D0)*ABS(XIH) 
      Q2(6,25) = 0 
      Q2(6,26) = 0 
      Q2(6,27) = 0 
      Q2(7,1) = 2*5**(0.5D0)*ABS(ETAH) 
      Q2(7,2) = 0 
      Q2(7,3) = 0 
      Q2(7,4) = -2*15**(0.5D0)*ETAH 
      Q2(7,5) = 0 
      Q2(7,6) = 0 
      Q2(7,7) = VT + 10*ABS(ETAH) + 2*ABS(MUH) 
     >   + 2*ABS(XIH) 
      Q2(7,8) = 0 
      Q2(7,9) = 2*5**(0.5D0)*ABS(MUH) 
      Q2(7,10) = 0 
      Q2(7,11) = 0 
      Q2(7,12) = 0 
      Q2(7,13) = 0 
      Q2(7,14) = 0 
      Q2(7,15) = 0 
      Q2(7,16) = 0 
      Q2(7,17) = 0 
      Q2(7,18) = 0 
      Q2(7,19) = 0 
      Q2(7,20) = 0 
      Q2(7,21) = 0 
      Q2(7,22) = 0 
      Q2(7,23) = 0 
      Q2(7,24) = 0 
      Q2(7,25) = 2*5**(0.5D0)*ABS(XIH) 
      Q2(7,26) = 0 
      Q2(7,27) = 0 
      Q2(8,1) = 0 
      Q2(8,2) = 2*5**(0.5D0)*ABS(ETAH) 
      Q2(8,3) = 0 
      Q2(8,4) = 0 
      Q2(8,5) = -2*15**(0.5D0)*ETAH 
      Q2(8,6) = 0 
      Q2(8,7) = 0 
      Q2(8,8) = VT + 10*ABS(ETAH) + 2*ABS(XIH) 
      Q2(8,9) = 2*15**(0.5D0)*MUH 
      Q2(8,10) = 0 
      Q2(8,11) = 0 
      Q2(8,12) = 0 
      Q2(8,13) = 0 
      Q2(8,14) = 0 
      Q2(8,15) = 0 
      Q2(8,16) = 0 
      Q2(8,17) = 0 
      Q2(8,18) = 0 
      Q2(8,19) = 0 
      Q2(8,20) = 0 
      Q2(8,21) = 0 
      Q2(8,22) = 0 
      Q2(8,23) = 0 
      Q2(8,24) = 0 
      Q2(8,25) = 0 
      Q2(8,26) = 2*5**(0.5D0)*ABS(XIH) 
      Q2(8,27) = 0 
      Q2(9,1) = 0 
      Q2(9,2) = 0 
      Q2(9,3) = -2*5**(0.5D0)*ABS(ETAH) 
      Q2(9,4) = 0 
      Q2(9,5) = 0 
      Q2(9,6) = 2*15**(0.5D0)*ETAH 
      Q2(9,7) = -2*5**(0.5D0)*ABS(MUH) 
      Q2(9,8) = 2*15**(0.5D0)*MUH 
      Q2(9,9) = - VT - 10*ABS(ETAH) - 10*ABS(MU
     >   H) - 2*ABS(XIH) 
      Q2(9,10) = 0 
      Q2(9,11) = 0 
      Q2(9,12) = 0 
      Q2(9,13) = 0 
      Q2(9,14) = 0 
      Q2(9,15) = 0 
      Q2(9,16) = 0 
      Q2(9,17) = 0 
      Q2(9,18) = 0 
      Q2(9,19) = 0 
      Q2(9,20) = 0 
      Q2(9,21) = 0 
      Q2(9,22) = 0 
      Q2(9,23) = 0 
      Q2(9,24) = 0 
      Q2(9,25) = 0 
      Q2(9,26) = 0 
      Q2(9,27) = -2*5**(0.5D0)*ABS(XIH) 
      Q2(10,1) = 0 
      Q2(10,2) = 0 
      Q2(10,3) = 0 
      Q2(10,4) = 0 
      Q2(10,5) = 0 
      Q2(10,6) = 0 
      Q2(10,7) = 0 
      Q2(10,8) = 0 
      Q2(10,9) = 0 
      Q2(10,10) = VT + 2*ABS(ETAH) + 2*ABS(MUH)
      Q2(10,11) = 0 
      Q2(10,12) = 2*5**(0.5D0)*ABS(MUH) 
      Q2(10,13) = 0 
      Q2(10,14) = 0 
      Q2(10,15) = 0 
      Q2(10,16) = 2*5**(0.5D0)*ABS(ETAH) 
      Q2(10,17) = 0 
      Q2(10,18) = 0 
      Q2(10,19) = 2*15**(0.5D0)*XIH 
      Q2(10,20) = 0 
      Q2(10,21) = 0 
      Q2(10,22) = 0 
      Q2(10,23) = 0 
      Q2(10,24) = 0 
      Q2(10,25) = 0 
      Q2(10,26) = 0 
      Q2(10,27) = 0 
      Q2(11,1) = 0 
      Q2(11,2) = 0 
      Q2(11,3) = 0 
      Q2(11,4) = 0 
      Q2(11,5) = 0 
      Q2(11,6) = 0 
      Q2(11,7) = 0 
      Q2(11,8) = 0 
      Q2(11,9) = 0 
      Q2(11,10) = 0 
      Q2(11,11) = VT + 2*ABS(ETAH) 
      Q2(11,12) = 2*15**(0.5D0)*MUH 
      Q2(11,13) = 0 
      Q2(11,14) = 0 
      Q2(11,15) = 0 
      Q2(11,16) = 0 
      Q2(11,17) = 2*5**(0.5D0)*ABS(ETAH) 
      Q2(11,18) = 0 
      Q2(11,19) = 0 
      Q2(11,20) = 2*15**(0.5D0)*XIH 
      Q2(11,21) = 0 
      Q2(11,22) = 0 
      Q2(11,23) = 0 
      Q2(11,24) = 0 
      Q2(11,25) = 0 
      Q2(11,26) = 0 
      Q2(11,27) = 0 
      Q2(12,1) = 0 
      Q2(12,2) = 0 
      Q2(12,3) = 0 
      Q2(12,4) = 0 
      Q2(12,5) = 0 
      Q2(12,6) = 0 
      Q2(12,7) = 0 
      Q2(12,8) = 0 
      Q2(12,9) = 0 
      Q2(12,10) = -2*5**(0.5D0)*ABS(MUH) 
      Q2(12,11) = 2*15**(0.5D0)*MUH 
      Q2(12,12) = - VT - 2*ABS(ETAH) - 10*ABS(M
     >   UH) 
      Q2(12,13) = 0 
      Q2(12,14) = 0 
      Q2(12,15) = 0 
      Q2(12,16) = 0 
      Q2(12,17) = 0 
      Q2(12,18) = -2*5**(0.5D0)*ABS(ETAH) 
      Q2(12,19) = 0 
      Q2(12,20) = 0 
      Q2(12,21) = -2*15**(0.5D0)*XIH 
      Q2(12,22) = 0 
      Q2(12,23) = 0 
      Q2(12,24) = 0 
      Q2(12,25) = 0 
      Q2(12,26) = 0 
      Q2(12,27) = 0 
      Q2(13,1) = 0 
      Q2(13,2) = 0 
      Q2(13,3) = 0 
      Q2(13,4) = 0 
      Q2(13,5) = 0 
      Q2(13,6) = 0 
      Q2(13,7) = 0 
      Q2(13,8) = 0 
      Q2(13,9) = 0 
      Q2(13,10) = 0 
      Q2(13,11) = 0 
      Q2(13,12) = 0 
      Q2(13,13) = VT + 2*ABS(MUH) 
      Q2(13,14) = 0 
      Q2(13,15) = 2*5**(0.5D0)*ABS(MUH) 
      Q2(13,16) = 2*15**(0.5D0)*ETAH 
      Q2(13,17) = 0 
      Q2(13,18) = 0 
      Q2(13,19) = 0 
      Q2(13,20) = 0 
      Q2(13,21) = 0 
      Q2(13,22) = 2*15**(0.5D0)*XIH 
      Q2(13,23) = 0 
      Q2(13,24) = 0 
      Q2(13,25) = 0 
      Q2(13,26) = 0 
      Q2(13,27) = 0 
      Q2(14,1) = 0 
      Q2(14,2) = 0 
      Q2(14,3) = 0 
      Q2(14,4) = 0 
      Q2(14,5) = 0 
      Q2(14,6) = 0 
      Q2(14,7) = 0 
      Q2(14,8) = 0 
      Q2(14,9) = 0 
      Q2(14,10) = 0 
      Q2(14,11) = 0 
      Q2(14,12) = 0 
      Q2(14,13) = 0 
      Q2(14,14) = VT 
      Q2(14,15) = 2*15**(0.5D0)*MUH 
      Q2(14,16) = 0 
      Q2(14,17) = 2*15**(0.5D0)*ETAH 
      Q2(14,18) = 0 
      Q2(14,19) = 0 
      Q2(14,20) = 0 
      Q2(14,21) = 0 
      Q2(14,22) = 0 
      Q2(14,23) = 2*15**(0.5D0)*XIH 
      Q2(14,24) = 0 
      Q2(14,25) = 0 
      Q2(14,26) = 0 
      Q2(14,27) = 0 
      Q2(15,1) = 0 
      Q2(15,2) = 0 
      Q2(15,3) = 0 
      Q2(15,4) = 0 
      Q2(15,5) = 0 
      Q2(15,6) = 0 
      Q2(15,7) = 0 
      Q2(15,8) = 0 
      Q2(15,9) = 0 
      Q2(15,10) = 0 
      Q2(15,11) = 0 
      Q2(15,12) = 0 
      Q2(15,13) = -2*5**(0.5D0)*ABS(MUH) 
      Q2(15,14) = 2*15**(0.5D0)*MUH 
      Q2(15,15) = - VT - 10*ABS(MUH) 
      Q2(15,16) = 0 
      Q2(15,17) = 0 
      Q2(15,18) = -2*15**(0.5D0)*ETAH 
      Q2(15,19) = 0 
      Q2(15,20) = 0 
      Q2(15,21) = 0 
      Q2(15,22) = 0 
      Q2(15,23) = 0 
      Q2(15,24) = -2*15**(0.5D0)*XIH 
      Q2(15,25) = 0 
      Q2(15,26) = 0 
      Q2(15,27) = 0 
      Q2(16,1) = 0 
      Q2(16,2) = 0 
      Q2(16,3) = 0 
      Q2(16,4) = 0 
      Q2(16,5) = 0 
      Q2(16,6) = 0 
      Q2(16,7) = 0 
      Q2(16,8) = 0 
      Q2(16,9) = 0 
      Q2(16,10) = -2*5**(0.5D0)*ABS(ETAH) 
      Q2(16,11) = 0 
      Q2(16,12) = 0 
      Q2(16,13) = 2*15**(0.5D0)*ETAH 
      Q2(16,14) = 0 
      Q2(16,15) = 0 
      Q2(16,16) = - VT - 10*ABS(ETAH) - 2*ABS(M
     >   UH) 
      Q2(16,17) = 0 
      Q2(16,18) = -2*5**(0.5D0)*ABS(MUH) 
      Q2(16,19) = 0 
      Q2(16,20) = 0 
      Q2(16,21) = 0 
      Q2(16,22) = 0 
      Q2(16,23) = 0 
      Q2(16,24) = 0 
      Q2(16,25) = -2*15**(0.5D0)*XIH 
      Q2(16,26) = 0 
      Q2(16,27) = 0 
      Q2(17,1) = 0 
      Q2(17,2) = 0 
      Q2(17,3) = 0 
      Q2(17,4) = 0 
      Q2(17,5) = 0 
      Q2(17,6) = 0 
      Q2(17,7) = 0 
      Q2(17,8) = 0 
      Q2(17,9) = 0 
      Q2(17,10) = 0 
      Q2(17,11) = -2*5**(0.5D0)*ABS(ETAH) 
      Q2(17,12) = 0 
      Q2(17,13) = 0 
      Q2(17,14) = 2*15**(0.5D0)*ETAH 
      Q2(17,15) = 0 
      Q2(17,16) = 0 
      Q2(17,17) = - VT - 10*ABS(ETAH) 
      Q2(17,18) = -2*15**(0.5D0)*MUH 
      Q2(17,19) = 0 
      Q2(17,20) = 0 
      Q2(17,21) = 0 
      Q2(17,22) = 0 
      Q2(17,23) = 0 
      Q2(17,24) = 0 
      Q2(17,25) = 0 
      Q2(17,26) = -2*15**(0.5D0)*XIH 
      Q2(17,27) = 0 
      Q2(18,1) = 0 
      Q2(18,2) = 0 
      Q2(18,3) = 0 
      Q2(18,4) = 0 
      Q2(18,5) = 0 
      Q2(18,6) = 0 
      Q2(18,7) = 0 
      Q2(18,8) = 0 
      Q2(18,9) = 0 
      Q2(18,10) = 0 
      Q2(18,11) = 0 
      Q2(18,12) = 2*5**(0.5D0)*ABS(ETAH) 
      Q2(18,13) = 0 
      Q2(18,14) = 0 
      Q2(18,15) = -2*15**(0.5D0)*ETAH 
      Q2(18,16) = 2*5**(0.5D0)*ABS(MUH) 
      Q2(18,17) = -2*15**(0.5D0)*MUH 
      Q2(18,18) = VT + 10*ABS(ETAH) + 10*ABS(MU
     >   H) 
      Q2(18,19) = 0 
      Q2(18,20) = 0 
      Q2(18,21) = 0 
      Q2(18,22) = 0 
      Q2(18,23) = 0 
      Q2(18,24) = 0 
      Q2(18,25) = 0 
      Q2(18,26) = 0 
      Q2(18,27) = 2*15**(0.5D0)*XIH 
      Q2(19,1) = 2*5**(0.5D0)*ABS(XIH) 
      Q2(19,2) = 0 
      Q2(19,3) = 0 
      Q2(19,4) = 0 
      Q2(19,5) = 0 
      Q2(19,6) = 0 
      Q2(19,7) = 0 
      Q2(19,8) = 0 
      Q2(19,9) = 0 
      Q2(19,10) = -2*15**(0.5D0)*XIH 
      Q2(19,11) = 0 
      Q2(19,12) = 0 
      Q2(19,13) = 0 
      Q2(19,14) = 0 
      Q2(19,15) = 0 
      Q2(19,16) = 0 
      Q2(19,17) = 0 
      Q2(19,18) = 0 
      Q2(19,19) = VT + 2*ABS(ETAH) + 2*ABS(MUH)
     >    + 10*ABS(XIH) 
      Q2(19,20) = 0 
      Q2(19,21) = 2*5**(0.5D0)*ABS(MUH) 
      Q2(19,22) = 0 
      Q2(19,23) = 0 
      Q2(19,24) = 0 
      Q2(19,25) = 2*5**(0.5D0)*ABS(ETAH) 
      Q2(19,26) = 0 
      Q2(19,27) = 0 
      Q2(20,1) = 0 
      Q2(20,2) = 2*5**(0.5D0)*ABS(XIH) 
      Q2(20,3) = 0 
      Q2(20,4) = 0 
      Q2(20,5) = 0 
      Q2(20,6) = 0 
      Q2(20,7) = 0 
      Q2(20,8) = 0 
      Q2(20,9) = 0 
      Q2(20,10) = 0 
      Q2(20,11) = -2*15**(0.5D0)*XIH 
      Q2(20,12) = 0 
      Q2(20,13) = 0 
      Q2(20,14) = 0 
      Q2(20,15) = 0 
      Q2(20,16) = 0 
      Q2(20,17) = 0 
      Q2(20,18) = 0 
      Q2(20,19) = 0 
      Q2(20,20) = VT + 2*ABS(ETAH) + 10*ABS(XIH
     >   ) 
      Q2(20,21) = 2*15**(0.5D0)*MUH 
      Q2(20,22) = 0 
      Q2(20,23) = 0 
      Q2(20,24) = 0 
      Q2(20,25) = 0 
      Q2(20,26) = 2*5**(0.5D0)*ABS(ETAH) 
      Q2(20,27) = 0 
      Q2(21,1) = 0 
      Q2(21,2) = 0 
      Q2(21,3) = -2*5**(0.5D0)*ABS(XIH) 
      Q2(21,4) = 0 
      Q2(21,5) = 0 
      Q2(21,6) = 0 
      Q2(21,7) = 0 
      Q2(21,8) = 0 
      Q2(21,9) = 0 
      Q2(21,10) = 0 
      Q2(21,11) = 0 
      Q2(21,12) = 2*15**(0.5D0)*XIH 
      Q2(21,13) = 0 
      Q2(21,14) = 0 
      Q2(21,15) = 0 
      Q2(21,16) = 0 
      Q2(21,17) = 0 
      Q2(21,18) = 0 
      Q2(21,19) = -2*5**(0.5D0)*ABS(MUH) 
      Q2(21,20) = 2*15**(0.5D0)*MUH 
      Q2(21,21) = - VT - 2*ABS(ETAH) - 10*ABS(M
     >   UH) - 10*ABS(XIH) 
      Q2(21,22) = 0 
      Q2(21,23) = 0 
      Q2(21,24) = 0 
      Q2(21,25) = 0 
      Q2(21,26) = 0 
      Q2(21,27) = -2*5**(0.5D0)*ABS(ETAH) 
      Q2(22,1) = 0 
      Q2(22,2) = 0 
      Q2(22,3) = 0 
      Q2(22,4) = -2*5**(0.5D0)*ABS(XIH) 
      Q2(22,5) = 0 
      Q2(22,6) = 0 
      Q2(22,7) = 0 
      Q2(22,8) = 0 
      Q2(22,9) = 0 
      Q2(22,10) = 0 
      Q2(22,11) = 0 
      Q2(22,12) = 0 
      Q2(22,13) = 2*15**(0.5D0)*XIH 
      Q2(22,14) = 0 
      Q2(22,15) = 0 
      Q2(22,16) = 0 
      Q2(22,17) = 0 
      Q2(22,18) = 0 
      Q2(22,19) = 0 
      Q2(22,20) = 0 
      Q2(22,21) = 0 
      Q2(22,22) = - VT - 2*ABS(MUH) - 10*ABS(XI
     >   H) 
      Q2(22,23) = 0 
      Q2(22,24) = -2*5**(0.5D0)*ABS(MUH) 
      Q2(22,25) = -2*15**(0.5D0)*ETAH 
      Q2(22,26) = 0 
      Q2(22,27) = 0 
      Q2(23,1) = 0 
      Q2(23,2) = 0 
      Q2(23,3) = 0 
      Q2(23,4) = 0 
      Q2(23,5) = -2*5**(0.5D0)*ABS(XIH) 
      Q2(23,6) = 0 
      Q2(23,7) = 0 
      Q2(23,8) = 0 
      Q2(23,9) = 0 
      Q2(23,10) = 0 
      Q2(23,11) = 0 
      Q2(23,12) = 0 
      Q2(23,13) = 0 
      Q2(23,14) = 2*15**(0.5D0)*XIH 
      Q2(23,15) = 0 
      Q2(23,16) = 0 
      Q2(23,17) = 0 
      Q2(23,18) = 0 
      Q2(23,19) = 0 
      Q2(23,20) = 0 
      Q2(23,21) = 0 
      Q2(23,22) = 0 
      Q2(23,23) = - VT - 10*ABS(XIH) 
      Q2(23,24) = -2*15**(0.5D0)*MUH 
      Q2(23,25) = 0 
      Q2(23,26) = -2*15**(0.5D0)*ETAH 
      Q2(23,27) = 0 
      Q2(24,1) = 0 
      Q2(24,2) = 0 
      Q2(24,3) = 0 
      Q2(24,4) = 0 
      Q2(24,5) = 0 
      Q2(24,6) = 2*5**(0.5D0)*ABS(XIH) 
      Q2(24,7) = 0 
      Q2(24,8) = 0 
      Q2(24,9) = 0 
      Q2(24,10) = 0 
      Q2(24,11) = 0 
      Q2(24,12) = 0 
      Q2(24,13) = 0 
      Q2(24,14) = 0 
      Q2(24,15) = -2*15**(0.5D0)*XIH 
      Q2(24,16) = 0 
      Q2(24,17) = 0 
      Q2(24,18) = 0 
      Q2(24,19) = 0 
      Q2(24,20) = 0 
      Q2(24,21) = 0 
      Q2(24,22) = 2*5**(0.5D0)*ABS(MUH) 
      Q2(24,23) = -2*15**(0.5D0)*MUH 
      Q2(24,24) = VT + 10*ABS(MUH) + 10*ABS(XIH
     >   ) 
      Q2(24,25) = 0 
      Q2(24,26) = 0 
      Q2(24,27) = 2*15**(0.5D0)*ETAH 
      Q2(25,1) = 0 
      Q2(25,2) = 0 
      Q2(25,3) = 0 
      Q2(25,4) = 0 
      Q2(25,5) = 0 
      Q2(25,6) = 0 
      Q2(25,7) = -2*5**(0.5D0)*ABS(XIH) 
      Q2(25,8) = 0 
      Q2(25,9) = 0 
      Q2(25,10) = 0 
      Q2(25,11) = 0 
      Q2(25,12) = 0 
      Q2(25,13) = 0 
      Q2(25,14) = 0 
      Q2(25,15) = 0 
      Q2(25,16) = 2*15**(0.5D0)*XIH 
      Q2(25,17) = 0 
      Q2(25,18) = 0 
      Q2(25,19) = -2*5**(0.5D0)*ABS(ETAH) 
      Q2(25,20) = 0 
      Q2(25,21) = 0 
      Q2(25,22) = 2*15**(0.5D0)*ETAH 
      Q2(25,23) = 0 
      Q2(25,24) = 0 
      Q2(25,25) = - VT - 10*ABS(ETAH) - 2*ABS(M
     >   UH) - 10*ABS(XIH) 
      Q2(25,26) = 0 
      Q2(25,27) = -2*5**(0.5D0)*ABS(MUH) 
      Q2(26,1) = 0 
      Q2(26,2) = 0 
      Q2(26,3) = 0 
      Q2(26,4) = 0 
      Q2(26,5) = 0 
      Q2(26,6) = 0 
      Q2(26,7) = 0 
      Q2(26,8) = 2*5**(0.5D0)*ABS(XIH) 
      Q2(26,9) = 0 
      Q2(26,10) = 0 
      Q2(26,11) = 0 
      Q2(26,12) = 0 
      Q2(26,13) = 0 
      Q2(26,14) = 0 
      Q2(26,15) = 0 
      Q2(26,16) = 0 
      Q2(26,17) = -2*15**(0.5D0)*XIH 
      Q2(26,18) = 0 
      Q2(26,19) = 0 
      Q2(26,20) = 2*5**(0.5D0)*ABS(ETAH) 
      Q2(26,21) = 0 
      Q2(26,22) = 0 
      Q2(26,23) = -2*15**(0.5D0)*ETAH 
      Q2(26,24) = 0 
      Q2(26,25) = 0 
      Q2(26,26) = VT + 10*ABS(ETAH) + 10*ABS(XI
     >   H) 
      Q2(26,27) = 2*15**(0.5D0)*MUH 
      Q2(27,1) = 0 
      Q2(27,2) = 0 
      Q2(27,3) = 0 
      Q2(27,4) = 0 
      Q2(27,5) = 0 
      Q2(27,6) = 0 
      Q2(27,7) = 0 
      Q2(27,8) = 0 
      Q2(27,9) = -2*5**(0.5D0)*ABS(XIH) 
      Q2(27,10) = 0 
      Q2(27,11) = 0 
      Q2(27,12) = 0 
      Q2(27,13) = 0 
      Q2(27,14) = 0 
      Q2(27,15) = 0 
      Q2(27,16) = 0 
      Q2(27,17) = 0 
      Q2(27,18) = 2*15**(0.5D0)*XIH 
      Q2(27,19) = 0 
      Q2(27,20) = 0 
      Q2(27,21) = -2*5**(0.5D0)*ABS(ETAH) 
      Q2(27,22) = 0 
      Q2(27,23) = 0 
      Q2(27,24) = 2*15**(0.5D0)*ETAH 
      Q2(27,25) = -2*5**(0.5D0)*ABS(MUH) 
      Q2(27,26) = 2*15**(0.5D0)*MUH 
      Q2(27,27) = - VT - 10*ABS(ETAH) - 10*ABS(
     >   MUH) - 10*ABS(XIH) 
      Q2(1,28) = -Q(01)*VOL(I2,J2,J,I,K) 
      Q2(2,28) = -Q(02)*VOL(I2,J2,J,I,K) 
      Q2(3,28) = Q(03)*VOL(I2,J2,J,I,K) 
      Q2(4,28) = -Q(04)*VOL(I2,J2,J,I,K) 
      Q2(5,28) = -Q(05)*VOL(I2,J2,J,I,K) 
      Q2(6,28) = Q(06)*VOL(I2,J2,J,I,K) 
      Q2(7,28) = Q(07)*VOL(I2,J2,J,I,K) 
      Q2(8,28) = Q(08)*VOL(I2,J2,J,I,K) 
      Q2(9,28) = -Q(09)*VOL(I2,J2,J,I,K) 
      Q2(10,28) = Q(10)*VOL(I2,J2,J,I,K) 
      Q2(11,28) = Q(11)*VOL(I2,J2,J,I,K) 
      Q2(12,28) = -Q(12)*VOL(I2,J2,J,I,K) 
      Q2(13,28) = Q(13)*VOL(I2,J2,J,I,K) 
      Q2(14,28) = Q(14)*VOL(I2,J2,J,I,K) 
      Q2(15,28) = -Q(15)*VOL(I2,J2,J,I,K) 
      Q2(16,28) = -Q(16)*VOL(I2,J2,J,I,K) 
      Q2(17,28) = -Q(17)*VOL(I2,J2,J,I,K) 
      Q2(18,28) = Q(18)*VOL(I2,J2,J,I,K) 
      Q2(19,28) = Q(19)*VOL(I2,J2,J,I,K) 
      Q2(20,28) = Q(20)*VOL(I2,J2,J,I,K) 
      Q2(21,28) = -Q(21)*VOL(I2,J2,J,I,K) 
      Q2(22,28) = -Q(22)*VOL(I2,J2,J,I,K) 
      Q2(23,28) = -Q(23)*VOL(I2,J2,J,I,K) 
      Q2(24,28) = Q(24)*VOL(I2,J2,J,I,K) 
      Q2(25,28) = -Q(25)*VOL(I2,J2,J,I,K) 
      Q2(26,28) = Q(26)*VOL(I2,J2,J,I,K) 
      Q2(27,28) = -Q(27)*VOL(I2,J2,J,I,K) 
      Q2(1,28) = Q2(1,28) - 2*XNK(1,1)*ABS(XI
     >   H) - 2*ABS(ETAH)*XNJ(1,1,K) - 2*ABS(MUH)*
     >   XNI(1,1,J2,K) 
      Q2(2,28) = Q2(2,28) - 2*XNK(2,1)*ABS(XI
     >   H) - 2*ABS(ETAH)*XNJ(2,1,K) 
      Q2(3,28) = Q2(3,28) + 2*XNK(3,1)*ABS(XIH)
     >    + 2*ABS(ETAH)*XNJ(3,1,K) + 2*5**(0.5D0)*
     >   ABS(MUH)*XNI(1,1,J2,K) 
      Q2(4,28) = Q2(4,28) - 2*XNK(1,2)*ABS(XI
     >   H) - 2*ABS(MUH)*XNI(2,1,J2,K) 
      Q2(5,28) = Q2(5,28) -2*XNK(2,2)*ABS(XIH)
      Q2(6,28) = Q2(6,28) + 2*XNK(3,2)*ABS(XIH)
     >    + 2*5**(0.5D0)*ABS(MUH)*XNI(2,1,J2,K) 
      Q2(7,28) = Q2(7,28) + 2*XNK(1,3)*ABS(XIH)
     >    + 2*ABS(MUH)*XNI(3,1,J2,K) + 2*5**(0.5D0
     >   )*ABS(ETAH)*XNJ(1,1,K) 
      Q2(8,28) = Q2(8,28) + 2*XNK(2,3)*ABS(XIH)
     >    + 2*5**(0.5D0)*ABS(ETAH)*XNJ(2,1,K) 
      Q2(9,28) = Q2(9,28) - 2*XNK(3,3)*ABS(XI
     >   H) - 2*5**(0.5D0)*ABS(ETAH)*XNJ(3,1,K) - 
     >   2*5**(0.5D0)*ABS(MUH)*XNI(3,1,J2,K) 
      Q2(10,28) = Q2(10,28) + 2*ABS(ETAH)*XNJ(1
     >   ,2,K) + 2*ABS(MUH)*XNI(1,2,J2,K) 
      Q2(11,28) = Q2(11,28) + 2*ABS(ETAH)*XNJ(2
     >   ,2,K) 
      Q2(12,28) = Q2(12,28) - 2*ABS(ETAH)*XNJ
     >   (3,2,K) - 2*5**(0.5D0)*ABS(MUH)*XNI(1,2,J
     >   2,K) 
      Q2(13,28) = Q2(13,28) + 2*ABS(MUH)*XNI(2,
     >   2,J2,K) 
      Q2(14,28) = Q2(14,28) + 0 
      Q2(15,28) = Q2(15,28) - 2*5**(0.5D0)*ABS
     >   (MUH)*XNI(2,2,J2,K) 
      Q2(16,28) = Q2(16,28) -  2*ABS(MUH)*XNI(
     >   3,2,J2,K) - 2*5**(0.5D0)*ABS(ETAH)*XNJ(1,
     >   2,K) 
      Q2(17,28) = Q2(17,28) - 2*5**(0.5D0)*ABS
     >   (ETAH)*XNJ(2,2,K) 
      Q2(18,28) = Q2(18,28) + 2*5**(0.5D0)*ABS(
     >   ETAH)*XNJ(3,2,K) + 2*5**(0.5D0)*ABS(MUH)*
     >   XNI(3,2,J2,K) 
      Q2(19,28) = Q2(19,28) + 2*ABS(ETAH)*XNJ(1
     >   ,3,K) + 2*ABS(MUH)*XNI(1,3,J2,K) + 2*5**(
     >   0.5D0)*XNK(1,1)*ABS(XIH)
      Q2(20,28) = Q2(20,28) + 2*ABS(ETAH)*XNJ(2
     >   ,3,K) + 2*5**(0.5D0)*XNK(2,1)*ABS(XIH)
      Q2(21,28) = Q2(21,28) -  2*ABS(ETAH)*XNJ
     >   (3,3,K) - 2*5**(0.5D0)*ABS(MUH)*XNI(1,3,J
     >   2,K) - 2*5**(0.5D0)*XNK(3,1)*ABS(XIH) 
      Q2(22,28) = Q2(22,28) -  2*ABS(MUH)*XNI(
     >   2,3,J2,K) - 2*5**(0.5D0)*XNK(1,2)*ABS(XIH
     >   ) 
      Q2(23,28) = Q2(23,28) - 2*5**(0.5D0)*XNK
     >   (2,2)*ABS(XIH) 
      Q2(24,28) = Q2(24,28) + 2*5**(0.5D0)*ABS(
     >   MUH)*XNI(2,3,J2,K) + 2*5**(0.5D0)*XNK(3,2
     >   )*ABS(XIH) 
      Q2(25,28) = Q2(25,28) -  2*ABS(MUH)*XNI(
     >   3,3,J2,K) - 2*5**(0.5D0)*ABS(ETAH)*XNJ(1,
     >   3,K) - 2*5**(0.5D0)*XNK(1,3)*ABS(XIH) 
      Q2(26,28) = Q2(26,28) + 2*5**(0.5D0)*ABS(
     >   ETAH)*XNJ(2,3,K) + 2*5**(0.5D0)*XNK(2,3)*
     >   ABS(XIH) 
      Q2(27,28) = Q2(27,28) -  2*5**(0.5D0)*AB
     >   S(ETAH)*XNJ(3,3,K) - 2*5**(0.5D0)*ABS(MUH
     >   )*XNI(3,3,J2,K) - 2*5**(0.5D0)*XNK(3,3)*A
     >   BS(XIH) 

*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
*** ---------------------------------------------------------------- ***
      ENDIF

      CALL ALSBD(IELEM**3,1,Q2,IER,IELEM**3)
      IF(IER.NE.0) CALL XABORT('SNFDH3: SINGULAR MATRIX.')
*
************************************************************************
************************************************************************
      IF(ICELZ.LT.NCELZ)THEN
      IF(IELEM.EQ.1) THEN
      BFLUZ(J,I,I2,J2,1,1)=2.0D0*Q2(1,2)-XNK(1,1)
      XNK(1,1) = BFLUZ(J,I,I2,J2,1,1)
      ELSE IF(IELEM.EQ.2) THEN
      BFLUZ(J,I,I2,J2,1,1) = XNK(1,1) + SIGN(1.0,XIH)*2*3
     >   **(0.5D0)*Q2(5,9) 
      BFLUZ(J,I,I2,J2,2,1) = XNK(2,1) + SIGN(1.0,XIH)*2*3
     >   **(0.5D0)*Q2(6,9) 
      BFLUZ(J,I,I2,J2,1,2) = XNK(1,2) + SIGN(1.0,XIH)*2*3
     >   **(0.5D0)*Q2(7,9) 
      BFLUZ(J,I,I2,J2,2,2) = XNK(2,2) + SIGN(1.0,XIH)*2*3
     >   **(0.5D0)*Q2(8,9) 
      XNK(:,:) = BFLUZ(J,I,I2,J2,:,:)
      ELSE IF(IELEM.EQ.3) THEN
      BFLUZ(J,I,I2,J2,1,1) = 2*Q2(01,28) - XNK(1,1) + 2*5**
     >   (0.5D0)*Q2(19,28) 
      BFLUZ(J,I,I2,J2,2,1) = 2*Q2(02,28) - XNK(2,1) + 2*5**
     >   (0.5D0)*Q2(20,28) 
      BFLUZ(J,I,I2,J2,3,1) = 2*Q2(03,28) - XNK(3,1) + 2*5**
     >   (0.5D0)*Q2(21,28) 
      BFLUZ(J,I,I2,J2,1,2) = 2*Q2(04,28) - XNK(1,2) + 2*5**
     >   (0.5D0)*Q2(22,28) 
      BFLUZ(J,I,I2,J2,2,2) = 2*Q2(05,28) - XNK(2,2) + 2*5**
     >   (0.5D0)*Q2(23,28) 
      BFLUZ(J,I,I2,J2,3,2) = 2*Q2(06,28) - XNK(3,2) + 2*5**
     >   (0.5D0)*Q2(24,28) 
      BFLUZ(J,I,I2,J2,1,3) = 2*Q2(07,28) - XNK(1,3) + 2*5**
     >   (0.5D0)*Q2(25,28) 
      BFLUZ(J,I,I2,J2,2,3) = 2*Q2(08,28) - XNK(2,3) + 2*5**
     >   (0.5D0)*Q2(26,28) 
      BFLUZ(J,I,I2,J2,3,3) = 2*Q2(09,28) - XNK(3,3) + 2*5**
     >   (0.5D0)*Q2(27,28) 
      XNK(:,:) = BFLUZ(J,I,I2,J2,:,:)
      ENDIF
      ELSEIF((ICELZ.EQ.NCELZ).OR.(K2.EQ.LZ))THEN
      IF(IELEM.EQ.1) THEN
      BFLUZ(J,I,I2,J2,1,1)=2.0D0*Q2(1,2)-XNK(1,1)
      ELSE IF(IELEM.EQ.2) THEN
      BFLUZ(J,I,I2,J2,1,1) = XNK(1,1) + SIGN(1.0,XIH)*2*3
     >   **(0.5D0)*Q2(5,9) 
      BFLUZ(J,I,I2,J2,2,1) = XNK(2,1) + SIGN(1.0,XIH)*2*3
     >   **(0.5D0)*Q2(6,9) 
      BFLUZ(J,I,I2,J2,1,2) = XNK(1,2) + SIGN(1.0,XIH)*2*3
     >   **(0.5D0)*Q2(7,9) 
      BFLUZ(J,I,I2,J2,2,2) = XNK(2,2) + SIGN(1.0,XIH)*2*3
     >   **(0.5D0)*Q2(8,9) 
      ELSE IF(IELEM.EQ.3) THEN
      BFLUZ(J,I,I2,J2,1,1) = 2*Q2(01,28) - XNK(1,1) + 2*5**
     >   (0.5D0)*Q2(19,28) 
      BFLUZ(J,I,I2,J2,2,1) = 2*Q2(02,28) - XNK(2,1) + 2*5**
     >   (0.5D0)*Q2(20,28) 
      BFLUZ(J,I,I2,J2,3,1) = 2*Q2(03,28) - XNK(3,1) + 2*5**
     >   (0.5D0)*Q2(21,28) 
      BFLUZ(J,I,I2,J2,1,2) = 2*Q2(04,28) - XNK(1,2) + 2*5**
     >   (0.5D0)*Q2(22,28) 
      BFLUZ(J,I,I2,J2,2,2) = 2*Q2(05,28) - XNK(2,2) + 2*5**
     >   (0.5D0)*Q2(23,28) 
      BFLUZ(J,I,I2,J2,3,2) = 2*Q2(06,28) - XNK(3,2) + 2*5**
     >   (0.5D0)*Q2(24,28) 
      BFLUZ(J,I,I2,J2,1,3) = 2*Q2(07,28) - XNK(1,3) + 2*5**
     >   (0.5D0)*Q2(25,28) 
      BFLUZ(J,I,I2,J2,2,3) = 2*Q2(08,28) - XNK(2,3) + 2*5**
     >   (0.5D0)*Q2(26,28) 
      BFLUZ(J,I,I2,J2,3,3) = 2*Q2(09,28) - XNK(3,3) + 2*5**
     >   (0.5D0)*Q2(27,28) 
      ENDIF
      INDXY=ICEL
      TMPBZ(J,INDXY,I2,J2,:,:)=BFLUZ(J,I,I2,J2,:,:)
      ENDIF
* -------------------------------------------------------------------- *
      IF(JL.LT.ISPLH)THEN
      IF(IELEM.EQ.1) THEN
      XNJ(1,1,K)=2.0D0*Q2(1,2)-XNJ(1,1,K)
      ELSEIF(IELEM.EQ.2) THEN
      XNJ(1,1,K) = XNJ(1,1,K) + SIGN(1.0,ETAH)*2*3**
     >   (0.5D0)*Q2(3,9) 
      XNJ(2,1,K) = XNJ(2,1,K) + SIGN(1.0,ETAH)*2*3**
     >   (0.5D0)*Q2(4,9) 
      XNJ(1,2,K) = XNJ(1,2,K) + SIGN(1.0,ETAH)*2*3**
     >   (0.5D0)*Q2(7,9) 
      XNJ(2,2,K) = XNJ(2,2,K) + SIGN(1.0,ETAH)*2*3**
     >   (0.5D0)*Q2(8,9) 
      ELSEIF(IELEM.EQ.3) THEN
      XNJ(1,1,K) = 2*Q2(01,28) - XNJ(1,1,K) + 2
     >   *5**(0.5D0)*Q2(07,28) 
      XNJ(2,1,K) = 2*Q2(02,28) - XNJ(2,1,K) + 2
     >   *5**(0.5D0)*Q2(08,28) 
      XNJ(3,1,K) = 2*Q2(03,28) - XNJ(3,1,K) + 2
     >   *5**(0.5D0)*Q2(09,28) 
      XNJ(1,2,K) = 2*Q2(10,28) - XNJ(1,2,K) + 2
     >   *5**(0.5D0)*Q2(16,28) 
      XNJ(2,2,K) = 2*Q2(11,28) - XNJ(2,2,K) + 2
     >   *5**(0.5D0)*Q2(17,28) 
      XNJ(3,2,K) = 2*Q2(12,28) - XNJ(3,2,K) + 2
     >   *5**(0.5D0)*Q2(18,28) 
      XNJ(1,3,K) = 2*Q2(19,28) - XNJ(1,3,K) + 2
     >   *5**(0.5D0)*Q2(25,28) 
      XNJ(2,3,K) = 2*Q2(20,28) - XNJ(2,3,K) + 2
     >   *5**(0.5D0)*Q2(26,28) 
      XNJ(3,3,K) = 2*Q2(21,28) - XNJ(3,3,K) + 2
     >   *5**(0.5D0)*Q2(27,28) 
      ENDIF
      ELSEIF((JL.EQ.ISPLH).AND.(IHEXJ.LE.NHEX))THEN
      I3=I2
      C1=1.0D0
      IF((J.EQ.1).AND.(ILOZJ.EQ.3)) THEN
         I3=ISPLH+1 -I2
         C1=-1.0D0
      ENDIF
      IF(IELEM.EQ.1) THEN
      BFLUX(CSIDEy,CHEXy,I3,1,1,K) = 2.0D0*Q2(1,2)-XNJ(1,1,K)
      ELSEIF(IELEM.EQ.2) THEN
      BFLUX(CSIDEy,CHEXy,I3,1,1,K) = XNJ(1,1,K) + SIGN(1.0,ETAH)
     >   *2*3**(0.5D0)*Q2(3,9) 
      BFLUX(CSIDEy,CHEXy,I3,2,1,K) = (XNJ(2,1,K) + SIGN(1.0,ETAH)
     >   *2*3**(0.5D0)*Q2(4,9))*C1
      BFLUX(CSIDEy,CHEXy,I3,1,2,K) = XNJ(1,2,K) + SIGN(1.0,ETAH)
     >   *2*3**(0.5D0)*Q2(7,9)
      BFLUX(CSIDEy,CHEXy,I3,2,2,K) = (XNJ(2,2,K) + SIGN(1.0,ETAH)
     >   *2*3**(0.5D0)*Q2(8,9))*C1
      ELSEIF(IELEM.EQ.3) THEN
      BFLUX(CSIDEy,CHEXy,I3,1,1,K) = 2*Q2(01,28) - XNJ(1,1,K) + 2
     >   *5**(0.5D0)*Q2(07,28) 
      BFLUX(CSIDEy,CHEXy,I3,2,1,K) = (2*Q2(02,28) - XNJ(2,1,K) + 2
     >   *5**(0.5D0)*Q2(08,28))*C1
      BFLUX(CSIDEy,CHEXy,I3,3,1,K) = 2*Q2(03,28) - XNJ(3,1,K) + 2
     >   *5**(0.5D0)*Q2(09,28)
      BFLUX(CSIDEy,CHEXy,I3,1,2,K) = 2*Q2(10,28) - XNJ(1,2,K) + 2
     >   *5**(0.5D0)*Q2(16,28)
      BFLUX(CSIDEy,CHEXy,I3,2,2,K) = (2*Q2(11,28) - XNJ(2,2,K) + 2
     >   *5**(0.5D0)*Q2(17,28))*C1
      BFLUX(CSIDEy,CHEXy,I3,3,2,K) = 2*Q2(12,28) - XNJ(3,2,K) + 2
     >   *5**(0.5D0)*Q2(18,28)
      BFLUX(CSIDEy,CHEXy,I3,1,3,K) = 2*Q2(19,28) - XNJ(1,3,K) + 2
     >   *5**(0.5D0)*Q2(25,28)
      BFLUX(CSIDEy,CHEXy,I3,2,3,K) = (2*Q2(20,28) - XNJ(2,3,K) + 2
     >   *5**(0.5D0)*Q2(26,28))*C1
      BFLUX(CSIDEy,CHEXy,I3,3,3,K) = 2*Q2(21,28) - XNJ(3,3,K) + 2
     >   *5**(0.5D0)*Q2(27,28) 
      ENDIF
      ENDIF
* -------------------------------------------------------------------- *
      IF(IL.LT.ISPLH)THEN
      IF(IELEM.EQ.1) THEN
      XNI(1,1,J2,K) = 2.0D0*Q2(1,2)-XNI(1,1,J2,K)
      ELSEIF(IELEM.EQ.2) THEN
      XNI(1,1,J2,K) = XNI(1,1,J2,K) + SIGN(1.0,MUH)*2*3**(0.5D0
     >   )*Q2(2,9) 
      XNI(2,1,J2,K) = XNI(2,1,J2,K) + SIGN(1.0,MUH)*2*3**(0.5D0
     >   )*Q2(4,9) 
      XNI(1,2,J2,K) = XNI(1,2,J2,K) + SIGN(1.0,MUH)*2*3**(0.5D0
     >   )*Q2(6,9) 
      XNI(2,2,J2,K) = XNI(2,2,J2,K) + SIGN(1.0,MUH)*2*3**(0.5D0
     >   )*Q2(8,9) 
      ELSEIF(IELEM.EQ.3) THEN
      XNI(1,1,J2,K) = 2*Q2(01,28) - XNI(1,1,J2,K)
     >    + 2*5**(0.5D0)*Q2(03,28) 
      XNI(2,1,J2,K) = 2*Q2(04,28) - XNI(2,1,J2,K)
     >    + 2*5**(0.5D0)*Q2(06,28) 
      XNI(3,1,J2,K) = 2*Q2(07,28) - XNI(3,1,J2,K)
     >    + 2*5**(0.5D0)*Q2(09,28) 
      XNI(1,2,J2,K) = 2*Q2(10,28) - XNI(1,2,J2,K)
     >    + 2*5**(0.5D0)*Q2(12,28) 
      XNI(2,2,J2,K) = 2*Q2(13,28) - XNI(2,2,J2,K)
     >    + 2*5**(0.5D0)*Q2(15,28) 
      XNI(3,2,J2,K) = 2*Q2(16,28) - XNI(3,2,J2,K)
     >    + 2*5**(0.5D0)*Q2(18,28) 
      XNI(1,3,J2,K) = 2*Q2(19,28) - XNI(1,3,J2,K)
     >    + 2*5**(0.5D0)*Q2(21,28) 
      XNI(2,3,J2,K) = 2*Q2(22,28) - XNI(2,3,J2,K)
     >    + 2*5**(0.5D0)*Q2(24,28) 
      XNI(3,3,J2,K) = 2*Q2(25,28) - XNI(3,3,J2,K)
     >    + 2*5**(0.5D0)*Q2(27,28) 
      ENDIF
      ELSEIF((IL.EQ.ISPLH).AND.(IHEXI.LE.NHEX))THEN
      J3=J2
      C1=1.0D0
      IF((J.EQ.3).AND.(ILOZI.EQ.1)) THEN
C          J3=ISPLH+1-J2
         C1=-1.0D0
      ENDIF
      IF(IELEM.EQ.1) THEN
      BFLUX(CSIDEx,CHEXx,J3,1,1,K) = 2.0D0*Q2(1,2)-XNI(1,1,J2,K)
      ELSEIF(IELEM.EQ.2) THEN
      BFLUX(CSIDEx,CHEXx,J3,1,1,K) = XNI(1,1,J2,K) + SIGN(1.0,
     >   MUH)*2*3**(0.5D0)*Q2(2,9) 
      BFLUX(CSIDEx,CHEXx,J3,2,1,K) = (XNI(2,1,J2,K) + SIGN(1.0,
     >   MUH)*2*3**(0.5D0)*Q2(4,9))*C1
      BFLUX(CSIDEx,CHEXx,J3,1,2,K) = XNI(1,2,J2,K) + SIGN(1.0,
     >   MUH)*2*3**(0.5D0)*Q2(6,9)
      BFLUX(CSIDEx,CHEXx,J3,2,2,K) = (XNI(2,2,J2,K) + SIGN(1.0,
     >   MUH)*2*3**(0.5D0)*Q2(8,9))*C1
      ELSEIF(IELEM.EQ.3) THEN
      BFLUX(CSIDEx,CHEXx,J3,1,1,K) = 2*Q2(01,28) - XNI(1,1,J2,K)
     >    + 2*5**(0.5D0)*Q2(03,28) 
      BFLUX(CSIDEx,CHEXx,J3,2,1,K) = (2*Q2(04,28) - XNI(2,1,J2,K)
     >    + 2*5**(0.5D0)*Q2(06,28))*C1
      BFLUX(CSIDEx,CHEXx,J3,3,1,K) = 2*Q2(07,28) - XNI(3,1,J2,K)
     >    + 2*5**(0.5D0)*Q2(09,28)
      BFLUX(CSIDEx,CHEXx,J3,1,2,K) = 2*Q2(10,28) - XNI(1,2,J2,K)
     >    + 2*5**(0.5D0)*Q2(12,28)
      BFLUX(CSIDEx,CHEXx,J3,2,2,K) = (2*Q2(13,28) - XNI(2,2,J2,K)
     >    + 2*5**(0.5D0)*Q2(15,28))*C1
      BFLUX(CSIDEx,CHEXx,J3,3,2,K) = 2*Q2(16,28) - XNI(3,2,J2,K)
     >    + 2*5**(0.5D0)*Q2(18,28)
      BFLUX(CSIDEx,CHEXx,J3,1,3,K) = 2*Q2(19,28) - XNI(1,3,J2,K)
     >    + 2*5**(0.5D0)*Q2(21,28)
      BFLUX(CSIDEx,CHEXx,J3,2,3,K) = (2*Q2(22,28) - XNI(2,3,J2,K)
     >    + 2*5**(0.5D0)*Q2(24,28))*C1
      BFLUX(CSIDEx,CHEXx,J3,3,3,K) = 2*Q2(25,28) - XNI(3,3,J2,K)
     >    + 2*5**(0.5D0)*Q2(27,28) 
      ENDIF
      ENDIF
************************************************************************
************************************************************************
      DO 210 ISCT=1,NSCT
      DO 200 IEL=1,IELEM**3
      FLUX(IEL,ISCT)=FLUX(IEL,ISCT)+W(M)*Q2(IEL,IELEM**3+1)*PL(ISCT,M)
  200 CONTINUE
  210 CONTINUE
      
      ILOZ=((J-1)*ISPLH+(J2-1))*ISPLH+I2
      FLUX_G(:,:,ILOZ,I,K,IG)=FLUX_G(:,:,ILOZ,I,K,IG)+FLUX(:,:)
*
      IF((NCODE(5).NE.1).or.(NCODE(6).NE.1))THEN
      IF(K2.EQ.LZ)THEN
         DO 240 IEL=1,IELEM
         DO 230 JEL=1,IELEM
         IOF=((((((M-1)*NHEX+(I-1))*3+(J-1))*ISPLH+(J2-1))*ISPLH+
     1      (I2-1))*IELEM+(JEL-1))*IELEM+IEL
         FUNKNO(L4+IOF,IG)=REAL(XNK(IEL,JEL))
         TMP(IOF,IG)=REAL(XNK(IEL,JEL))
  230    CONTINUE
  240    CONTINUE
      ENDIF
      ENDIF
*--
  250 CONTINUE ! JL 
  260 CONTINUE ! IL

  270 CONTINUE ! lozenges

  280 CONTINUE ! zcell

  290 CONTINUE ! directions
  300 CONTINUE ! energy
*$OMP END PARALLEL DO

  304 CONTINUE ! macrocells

      ISTTZ=TASKLST(2,1,IDI,IND2)
      IENDZ=TASKLST(2,TASKSINWAVE(IDI),IDI,IND2)
      ISTTZ = (ISTTZ-1)*NCELZ +1
      IENDZ = IENDZ*NCELZ 
      IF(IENDZ.GT.LZ) IENDZ=LZ

      IF(IND.GE.7) THEN
         ISTTZ = TASKLST(2,TASKSINWAVE(IDI),IDI,IND2)
         ISTTZ = ISTTZ*NCELZ
         ISTTZ = LZ-ISTTZ+1
         IF(ISTTZ.LE.0) ISTTZ=1

         IENDZ = TASKLST(2,1,IDI,IND2)
         IENDZ = NMZ+1 -IENDZ
         IENDZ = IENDZ*NCELZ 
         IF(IENDZ.GT.LZ) IENDZ=LZ
      ENDIF

      ILENZ = (IENDZ-ISTTZ)+1

      NHEXTOPASS = TASKSINWAVE(IDI)
      DO ICEL=1,TASKSINWAVE(IDI)
         I=TASKLST(1,ICEL,IDI,IND2)
         BFLUZ(:,I,:,:,:,:)=TMPBZ(:,ICEL,:,:,:,:)
      ENDDO

  305 CONTINUE ! wavefront number

      ! This is necessary so that all octants use the boundary 
      ! fluxes of the previous iteration. Otherwise, some processes
      ! might use those of the previous iteration while others
      ! might use those of the current iteration. This prevents
      ! proper convergence. This is especially an issue when the
      ! number of processes is between 1 < n < (max. effective).
      IF((NCODE(5).NE.1).or.(NCODE(6).NE.1))THEN
      DO 778 IG=1,NGEFF
        IF(.NOT.INCONV(IG)) GO TO 778
        FUNKNO(L4+1:L4+L6,IG)=REAL(TMP2(1:L6,IG))
  778 CONTINUE
      ENDIF

  310 CONTINUE ! dodecant
************************************************************************
************************************************************************
      DO 777 IG=1,NGEFF
        IF(.NOT.INCONV(IG)) GO TO 777
        FUNKNO(:L4,IG)=
     1  RESHAPE(REAL(FLUX_G(:IELEM**3,:NSCT,:3*ISPLH**2,:NHEX,:LZ,
     2  IG)), (/L4/))
      IF((NCODE(5).NE.1).or.(NCODE(6).NE.1))THEN
         FUNKNO(L4+1:L4+L6,IG)=REAL(TMP(1:L6,IG))
      ENDIF
  777 CONTINUE
************************************************************************
************************************************************************
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(FLUX_G,FLUX,INDANG)
      DEALLOCATE(BFLUX,BFLUZ,TMPBZ)
      DEALLOCATE(TMP,TMP2)
      RETURN
*
  500 FORMAT(16H SNFDH3: thread=,I8,12H --->(group=,I4,7H angle=,I4,
     1 1H))
      END
