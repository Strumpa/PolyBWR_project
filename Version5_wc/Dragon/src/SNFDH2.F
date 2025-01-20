*DECK SNFDH2
      SUBROUTINE SNFDH2(NUN,NGEFF,IMPX,INCONV,NGIND,NHEX,NBC,IELEM,
     1 ISPLH,NMAT,NPQ,NSCT,MAT,VOL,TOTAL,QEXT,DU,DE,W,DB,DA,CONNEC,
     2 SIDE,CONFROM,TASKLST,TASKSINWAVE,CON_WHOAMI,CON_TOWHOM,FUNKNO,
     3 MN,DN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform one inner iteration for solving SN equations in 2D hexagonal
* geometry for the DISCONTINUOUS GALERKIN method. KBA algorithm and 
* angle parallelisation. VOID boundary conditions.
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
* NUN     total number of unknowns in vector FUNKNO.
* NGEFF   number of energy groups processed in parallel.
* IMPX    print flag (equal to zero for no print).
* INCONV  energy group convergence flag (set to .FALSE. if converged).
* NGIND   energy group indices assign to the NGEFF set.
* NHEX    number of hexagons in X-Y plane.
* IELEM   measure of order of the spatial approximation polynomial:
*         =1 constant (classical diamond scheme - default for HODD);
*         =2 linear (default for DG);
*         =3 parabolic;
*         =4 cubic (only for DG).
* NMAT    number of material mixtures.
* ISPLH   splitting option for hexagons.
* NPQ     number of SN directions in 6 dodecants (including zero-weight
*         directions).
* NSCT    maximum number of spherical harmonics moments of the flux.
* MAT     material mixture index in each region.
* VOL     volumes of each region.
* TOTAL   macroscopic total cross sections.
* QEXT    Legendre components of the fixed source.
* DU      first direction cosines ($\\mu$).
* DE      second direction cosines ($\\eta$).
* W       weights.
* DB      diamond-scheme parameter.
* DA      diamond-scheme parameter.
* MN      moment-to-discrete matrix.
* DN      discrete-to-moment matrix.
* SIDE    side of an hexagon.
* CONNEC  connectivity matrix for flux swapping -- which lozenges is the
*         lozenge under consideration connected to; in order to pass the
*         flux along. This is dependent on direction.
* CONFROM matrix for incoming flux -- which lozenges are feeding into
*         the lozenge under consideration. This is dependent on
*         direction.
*
*Parameters: output
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
      INTEGER NUN,NGEFF,IMPX,NGIND(NGEFF),NHEX,IELEM,ISPLH,NMAT,NPQ,
     1 NSCT,MAT(ISPLH,ISPLH,3,NHEX),NBC,CONNEC(3,(NHEX*3)*2,6),
     2 CONFROM(2,3,6),TASKLST(2,NBC,((NBC-1)*4 +1),6), 
     3 TASKSINWAVE(((NBC-1)*4 +1))
      LOGICAL INCONV(NGEFF)
      REAL SIDE,VOL(ISPLH,ISPLH,3,NHEX),TOTAL(0:NMAT,NGEFF),
     1 QEXT(NUN,NGEFF),DU(NPQ),DB(ISPLH,ISPLH,3,NHEX,NPQ),
     2 DA(ISPLH,ISPLH,3,NHEX,NPQ),DE(NPQ),W(NPQ),
     3 FUNKNO(NUN,NGEFF),MN(NPQ,NSCT),DN(NSCT,NPQ)
*----
*  LOCAL VARIABLES
*----
      INTEGER NPQD(6),IIND(6),P
      DOUBLE PRECISION :: Q(IELEM**2), Q2(IELEM**2,(IELEM**2)+1),
     1 VT, CONST0, CONST1, CONST2
      PARAMETER(IUNOUT=6,PI=3.141592654)
      INTEGER :: ILOZSWP(3,6), IFROMI, IFROMJ
      REAL :: JAC(2,2,3), MUH, ETAH, AAA, BBB, CCC, DDD, MUHTEMP,
     1 ETAHTEMP
      DOUBLE PRECISION :: THETA, XNI(IELEM,ISPLH), XNJ(IELEM)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: INDANG
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: FLUX
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: FLUX_G
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: BFLUX
      INTEGER MAXWAVECNT
      INTEGER :: CON_WHOAMI(2,(NHEX*3)*2,6), CON_TOWHOM(2,(NHEX*3)*2,6)
      INTEGER CHEXx,CSIDEx,CHEXy,CSIDEy
*----
*  DEFINITION OF CONSTANTS.
*----
      L4=3*(IELEM**2)*(ISPLH**2)*NHEX*NSCT
      CONST0=2.0D0*DSQRT(3.0D0)
      CONST1=2.0D0*DSQRT(5.0D0)
      CONST2=2.0D0*DSQRT(15.0D0)
      MAXWAVECNT = ((NBC-1)*4 +1)
      NCOL=2*NBC -1
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(INDANG(NPQ,6))
      ALLOCATE(FLUX(IELEM**2,NSCT))
      ALLOCATE(FLUX_G(IELEM**2,NSCT,3*ISPLH**2,NHEX,NGEFF))
      ALLOCATE(BFLUX(3,NCOL,ISPLH,IELEM,NPQ/6))

      FLUX_G(:IELEM**2,:NSCT,:3*ISPLH**2,:NHEX,:NGEFF)=0.0D0
*----
*  CONSTRUCT JACOBIAN MATRIX FOR EACH LOZENGE
*----
      ILOZSWP = RESHAPE((/ 3, 2, 1, 3, 1, 2, 1, 3, 2, 1, 2, 3, 2, 1,
     1   3, 2, 3, 1 /), SHAPE(ILOZSWP))
      JAC = RESHAPE((/ 1., -SQRT(3.), 1., SQRT(3.), 2., 0., 1.,
     1    SQRT(3.), 2., 0., -1., SQRT(3.) /), SHAPE(JAC))
      JAC = (SIDE/2.)*JAC

************************************************************************
************************************************************************
************************************************************************

*----
*  PARAMETER VALIDATION.
*----
      IF((IELEM.LT.1).OR.(IELEM.GT.3))
     1   CALL XABORT('SNFDH2: INVALID IELEM (DIAM) VALUE. '
     2   //'CHECK INPUT DATA FILE.')
*----
*  SET DODECANT SWAPPING ORDER.
*----
      NPQD(:6)=0
      INDANG(:NPQ,:6)=0
      IIND(:6)=0
      DO 10 M=1,NPQ
        VU=DU(M)
        VE=DE(M)
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
        IND=0
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
        IIND(IND)=IND ! assume IIND(I)=I in hexagonal geometry
        NPQD(IND)=NPQD(IND)+1
        INDANG(NPQD(IND),IND)=M
   10 CONTINUE
*----
*  MAIN LOOP OVER DODECANTS.
*----
      DO 240 JND=1,6
      IND=IIND(JND)
      IF(IND.EQ.0) GO TO 240 ! needed because of S2 LS (4 dir. for 6 sextants)
*----
*  KBA DIAGONAL LOOP
*----
      BFLUX(:3,:NCOL,:ISPLH,:IELEM,:)=0.0D0
      DO 230 IDI=1,MAXWAVECNT

* Loop over hexes in same sweep front
      DO 220 ICEL=1,TASKSINWAVE(IDI)
      I=TASKLST(1,ICEL,IDI,IND)
*----
*  MAIN SWAPPING LOOPS FOR SN FLUX CALCULATION
*----
*$OMP  PARALLEL DO
      DO 210 IG=1,NGEFF
      DO 200 IPQD=1,NPQD(IND)
#if defined(_OPENMP)
        ITID=omp_get_thread_num()
#else
        ITID=0
#endif
      IF(IMPX.GT.5) WRITE(IUNOUT,400) ITID,NGIND(IG),IPQD
      IF(.NOT.INCONV(IG)) GO TO 200
      M=INDANG(IPQD,IND)
      IF(W(M).EQ.0.0) GO TO 200

* Loop over lozenges
      DO 180 JJ=1,3
      J=ILOZSWP(JJ,IND)
*
      AAA = JAC(1,1,J)
      BBB = JAC(1,2,J)
      CCC = JAC(2,1,J)
      DDD = JAC(2,2,J)

* Identify DOMAIN-WISE indices of connected sides 
      IHEXI  = CONNEC(1,((I-1)*3*2) + ((J-1)*2) +1,IND)
      ILOZI  = CONNEC(2,((I-1)*3*2) + ((J-1)*2) +1,IND)
      IHEXJ  = CONNEC(1,((I-1)*3*2) + ((J-1)*2) +2,IND)
      ILOZJ  = CONNEC(2,((I-1)*3*2) + ((J-1)*2) +2,IND)
      IFROMI = CONFROM(1,J,IND)
      IFROMJ = CONFROM(2,J,IND)

* Identify BF-ARRAY-WISE indices of sides of current lozenge
      MHEXx  = CON_WHOAMI(1,((I-1)*3*2) + ((J-1)*2) +1,IND)
      MSIDEx = CON_WHOAMI(2,((I-1)*3*2) + ((J-1)*2) +1,IND)
      MHEXy  = CON_WHOAMI(1,((I-1)*3*2) + ((J-1)*2) +2,IND)
      MSIDEy = CON_WHOAMI(2,((I-1)*3*2) + ((J-1)*2) +2,IND)

* Identify BF-ARRAY-WISE indices of connected sides
      CHEXx  = CON_TOWHOM(1,((I-1)*3*2) + ((J-1)*2) +1,IND)
      CSIDEx = CON_TOWHOM(2,((I-1)*3*2) + ((J-1)*2) +1,IND)
      CHEXy  = CON_TOWHOM(1,((I-1)*3*2) + ((J-1)*2) +2,IND)
      CSIDEy = CON_TOWHOM(2,((I-1)*3*2) + ((J-1)*2) +2,IND)

* Loop over internal split of lozenges
      DO 170 IL=1,ISPLH
      I2=IL
      IF(CONFROM(1,J,IND).EQ.3) I2=ISPLH+1-I2
      I4=I2

* It is important that when going from lozenge 3 to 1, the Boundary 
* Surface  Fluxes (BSF) are not stored in "reverse" order; as it is 
* done when going from loz. 1 to 3 (this happens when storing the BSF 
* after flux computation below). This is because the BSF in this KBA 
* sweep are stored in a "wavefront" array, where the different x and y 
* boundaries of the lozenges are not tracked. Therefore, for lozenge 1, 
* the BSF are instead  picked up in "reverse" order (line below).
      IF(J.EQ.1) I4=ISPLH+1-I2

      DO IEL=1,IELEM
         XNJ(IEL) = BFLUX(MSIDEy,MHEXy,I4,IEL,IPQD)
         BFLUX(MSIDEy,MHEXy,I4,IEL,IPQD) = 0
      ENDDO

      DO 160 JL=1,ISPLH
      J2=JL
      IF(CONFROM(2,J,IND).EQ.4) J2=ISPLH+1-J2
      J4=J2

      IF(IL.EQ.1)THEN
         DO IEL=1,IELEM
            XNI(IEL,J2) = BFLUX(MSIDEx,MHEXx,J4,IEL,IPQD)
            BFLUX(MSIDEx,MHEXx,J4,IEL,IPQD) = 0
         ENDDO
      ENDIF

      FLUX(:IELEM**2,:NSCT)=0.0D0
*
      MUHTEMP  =  DA(I2,J2,J,I,M)
      ETAHTEMP =  DB(I2,J2,J,I,M)
      MUH = (MUHTEMP*DDD) - (ETAHTEMP*BBB)
      ETAH = (-MUHTEMP*CCC) + (ETAHTEMP*AAA)
*
      IF(MAT(I2,J2,J,I).EQ.0) GO TO 180

*     -----------------------------------------------------
      DO 110 IEL=1,IELEM**2
      Q(IEL)=0.0D0
      DO 100 P=1,NSCT
      IOF=((((((I-1)*3+(J-1))*ISPLH+(J2-1))*ISPLH+(I2-1))*NSCT+
     1 (P-1))*IELEM**2)+IEL
      Q(IEL)=Q(IEL)+QEXT(IOF,IG)*MN(M,P)
  100 CONTINUE
  110 CONTINUE
*     -----------------------------------------------------
*
      VT=VOL(I2,J2,J,I)*TOTAL(MAT(I2,J2,J,I),IG)
      Q2(:IELEM**2,:(IELEM**2)+1)=0.0D0
*
*----------
* BUILD MATRIX FOR THE RELEVANT ORDER
* Note: the matrix does seem to be skew-symmetric (or antisymmetric), 
* i.e. A^T = -A. However, from a numerical perspective, it might be 
* more efficient to assign values than to sweep over the half matrix.
*     -----------------------------------------------------
      IF(IELEM.EQ.1) THEN
      Q2(1,1) = VT + ABS(ETAH) + ABS(MUH) 
      Q2(1,2) = Q(1)*VOL(I2,J2,J,I) 

      Q2(1,2) = Q2(1,2) + ETAH*SIGN(1.0,ETAH)*XNJ(1) + 
     >   MUH*SIGN(1.0,MUH)*XNI(1,J2) 
      ELSE IF(IELEM.EQ.2) THEN
*     alpha=1 , beta=1
      Q2(1,1) = VT + ABS(ETAH) + ABS(MUH) 
      Q2(1,2) = 3**(0.5D0)*MUH 
      Q2(1,3) = 3**(0.5D0)*ETAH 
*     alpha=2 , beta=1
      Q2(2,1) = -3**(0.5D0)*MUH 
      Q2(2,2) = VT + ABS(ETAH) + 3*ABS(MUH) 
      Q2(2,4) = 3**(0.5D0)*ETAH 
*     alpha=1 , beta=2
      Q2(3,1) = -3**(0.5D0)*ETAH 
      Q2(3,3) = VT + 3*ABS(ETAH) + ABS(MUH) 
      Q2(3,4) = 3**(0.5D0)*MUH 
*     alpha=2 , beta=2
      Q2(4,2) = -3**(0.5D0)*ETAH 
      Q2(4,3) = -3**(0.5D0)*MUH 
      Q2(4,4) = VT + 3*ABS(ETAH) + 3*ABS(MUH) 
*     source
      Q2(1,5) = Q(1)*VOL(I2,J2,J,I)
      Q2(2,5) = Q(2)*VOL(I2,J2,J,I)
      Q2(3,5) = Q(3)*VOL(I2,J2,J,I)
      Q2(4,5) = Q(4)*VOL(I2,J2,J,I)
*     jump terms (dependent on direction)
      Q2(1,5) = Q2(1,5) + ETAH*SIGN(1.0,ETAH)*XNJ(1) + 
     >   MUH*SIGN(1.0,MUH)*XNI(1,J2) 
      Q2(2,5) = Q2(2,5) + ETAH*SIGN(1.0,ETAH)*XNJ(2) - 3**(0.
     >   5D0)*MUH*XNI(1,J2) 
      Q2(3,5) = Q2(3,5) + MUH*SIGN(1.0,MUH)*XNI(2,J2) - 3**(
     >   0.5D0)*ETAH*XNJ(1) 
      Q2(4,5) = Q2(4,5) - 3**(0.5D0)*ETAH*
     >   XNJ(2) - 3**(0.5D0)*MUH*XNI(2,J2) 
*** ---------------------------------------------------------------- ***
      ELSE IF(IELEM.EQ.3) THEN
      Q2(1,1) = VT + ABS(ETAH) + ABS(MUH) 
      Q2(1,2) = 3**(0.5D0)*MUH 
      Q2(1,3) = 5**(0.5D0)*ABS(MUH) 
      Q2(1,4) = 3**(0.5D0)*ETAH 
      Q2(1,7) = 5**(0.5D0)*ABS(ETAH) 
      Q2(2,1) = -3**(0.5D0)*MUH 
      Q2(2,2) = VT + ABS(ETAH) + 3*ABS(MUH) 
      Q2(2,3) = 15**(0.5D0)*MUH 
      Q2(2,5) = 3**(0.5D0)*ETAH 
      Q2(2,8) = 5**(0.5D0)*ABS(ETAH) 
      Q2(3,1) = 5**(0.5D0)*ABS(MUH) 
      Q2(3,2) = -15**(0.5D0)*MUH 
      Q2(3,3) = VT + ABS(ETAH) + 5*ABS(MUH) 
      Q2(3,6) = 3**(0.5D0)*ETAH 
      Q2(3,9) = 5**(0.5D0)*ABS(ETAH) 
      Q2(4,1) = -3**(0.5D0)*ETAH 
      Q2(4,4) = VT + 3*ABS(ETAH) + ABS(MUH) 
      Q2(4,5) = 3**(0.5D0)*MUH 
      Q2(4,6) = 5**(0.5D0)*ABS(MUH) 
      Q2(4,7) = 15**(0.5D0)*ETAH 
      Q2(5,2) = -3**(0.5D0)*ETAH 
      Q2(5,4) = -3**(0.5D0)*MUH 
      Q2(5,5) = VT + 3*ABS(ETAH) + 3*ABS(MUH) 
      Q2(5,6) = 15**(0.5D0)*MUH 
      Q2(5,8) = 15**(0.5D0)*ETAH 
      Q2(6,3) = -3**(0.5D0)*ETAH 
      Q2(6,4) = 5**(0.5D0)*ABS(MUH) 
      Q2(6,5) = -3**(0.5D0)*5**(0.5D0)*MUH 
      Q2(6,6) = VT + 3*ABS(ETAH) + 5*ABS(MUH) 
      Q2(6,9) = 2*15**(0.5D0)*ETAH - 3**(0.5D0)*5**(0.5D0)*
     >   ETAH 
      Q2(7,1) = 5**(0.5D0)*ABS(ETAH) 
      Q2(7,4) = -15**(0.5D0)*ETAH 
      Q2(7,7) = VT + 5*ABS(ETAH) + ABS(MUH) 
      Q2(7,8) = 3**(0.5D0)*MUH 
      Q2(7,9) = 5**(0.5D0)*ABS(MUH) 
      Q2(8,2) = 5**(0.5D0)*ABS(ETAH) 
      Q2(8,5) = -3**(0.5D0)*5**(0.5D0)*ETAH 
      Q2(8,7) = -3**(0.5D0)*MUH 
      Q2(8,8) = VT + 5*ABS(ETAH) + 3*ABS(MUH) 
      Q2(8,9) = 2*15**(0.5D0)*MUH - 3**(0.5D0)*5**(0.5D0)*
     >   MUH 
      Q2(9,3) = 5**(0.5D0)*ABS(ETAH) 
      Q2(9,6) = -15**(0.5D0)*ETAH 
      Q2(9,7) = 5**(0.5D0)*ABS(MUH) 
      Q2(9,8) = -15**(0.5D0)*MUH 
      Q2(9,9) = VT + 5*ABS(ETAH) + 5*ABS(MUH) 
      Q2(1,10) = Q(1)*VOL(I2,J2,J,I) 
      Q2(2,10) = Q(2)*VOL(I2,J2,J,I) 
      Q2(3,10) = Q(3)*VOL(I2,J2,J,I) 
      Q2(4,10) = Q(4)*VOL(I2,J2,J,I) 
      Q2(5,10) = Q(5)*VOL(I2,J2,J,I) 
      Q2(6,10) = Q(6)*VOL(I2,J2,J,I) 
      Q2(7,10) = Q(7)*VOL(I2,J2,J,I) 
      Q2(8,10) = Q(8)*VOL(I2,J2,J,I) 
      Q2(9,10) = Q(9)*VOL(I2,J2,J,I) 
      Q2(1,10) = Q2(1,10) + ETAH*XNJ(1)*SIGN(1.0,ETAH) + 
     >   MUH*XNI(1,J2)*SIGN(1.0,MUH) 
      Q2(2,10) = Q2(2,10) + ETAH*XNJ(2)*SIGN(1.0,ETAH) - 3**(
     >   0.5D0)*MUH*XNI(1,J2) 
      Q2(3,10) = Q2(3,10) + ETAH*XNJ(3)*SIGN(1.0,ETAH) + 5**(
     >   0.5D0)*MUH*XNI(1,J2)*SIGN(1.0,MUH) 
      Q2(4,10) = Q2(4,10) + MUH*XNI(2,J2)*SIGN(1.0,MUH) - 3*
     >   *(0.5D0)*ETAH*XNJ(1) 
      Q2(5,10) = Q2(5,10) - 3**(0.5D0)*ETAH
     >   *XNJ(2) - 3**(0.5D0)*MUH*XNI(2,J2) 
      Q2(6,10) = Q2(6,10) + 5**(0.5D0)*MUH
     >   *XNI(2,J2)*SIGN(1.0,MUH) - 3**(0.5D0)*ETAH*XNJ(3) 
      Q2(7,10) = Q2(7,10) + MUH*XNI(3,J2)*SIGN(1.0,MUH) + 5*
     >   *(0.5D0)*ETAH*XNJ(1)*SIGN(1.0,ETAH) 
      Q2(8,10) = Q2(8,10) + 5**(0.5D0)*ETAH*SIGN(1.0,ETAH)*
     >   XNJ(2) - 3**(0.5D0)*MUH*XNI(3,J2) 
      Q2(9,10) = Q2(9,10) + 5**(0.5D0)*ETAH*SIGN(1.0,ETAH)*
     >   XNJ(3) + 5**(0.5D0)*MUH*XNI(3,J2)*SIGN(1.0,MUH) 
*** ---------------------------------------------------------------- ***
      ELSEIF(IELEM.EQ.4) THEN
      Q2(1,1) = VT + ABS(ETAH) + ABS(MUH) 
      Q2(1,2) = 3**(0.5D0)*MUH 
      Q2(1,3) = 5**(0.5D0)*ABS(MUH) 
      Q2(1,4) = 7**(0.5D0)*MUH 
      Q2(1,5) = 3**(0.5D0)*ETAH 
      Q2(1,9) = 5**(0.5D0)*ABS(ETAH) 
      Q2(1,13) = 7**(0.5D0)*ETAH 

      Q2(2,1) = -3**(0.5D0)*MUH 
      Q2(2,2) = VT + ABS(ETAH) + 3*ABS(MUH) 
      Q2(2,3) = 15**(0.5D0)*MUH 
      Q2(2,4) = 21**(0.5D0)*ABS(MUH) 
      Q2(2,6) = 3**(0.5D0)*ETAH 
      Q2(2,10) = 5**(0.5D0)*ABS(ETAH) 
      Q2(2,14) = 7**(0.5D0)*ETAH 

      Q2(3,1) = 5**(0.5D0)*ABS(MUH) 
      Q2(3,2) = -15**(0.5D0)*MUH 
      Q2(3,3) = VT + ABS(ETAH) + 5*ABS(MUH) 
      Q2(3,4) = 35**(0.5D0)*MUH 
      Q2(3,7) = 3**(0.5D0)*ETAH 
      Q2(3,11) = 5**(0.5D0)*ABS(ETAH) 
      Q2(3,15) = 7**(0.5D0)*ETAH 

      Q2(4,1) = -7**(0.5D0)*MUH 
      Q2(4,2) = 21**(0.5D0)*ABS(MUH) 
      Q2(4,3) = -35**(0.5D0)*MUH 
      Q2(4,4) = VT + ABS(ETAH) + 7*ABS(MUH) 
      Q2(4,8) = 3**(0.5D0)*ETAH 
      Q2(4,12) = 5**(0.5D0)*ABS(ETAH) 
      Q2(4,16) = 7**(0.5D0)*ETAH 

      Q2(5,1) = -3**(0.5D0)*ETAH 
      Q2(5,5) = VT + 3*ABS(ETAH) + ABS(MUH) 
      Q2(5,6) = 3**(0.5D0)*MUH 
      Q2(5,7) = 5**(0.5D0)*ABS(MUH) 
      Q2(5,8) = 7**(0.5D0)*MUH 
      Q2(5,9) = 15**(0.5D0)*ETAH 
      Q2(5,13) = 21**(0.5D0)*ABS(ETAH) 

      Q2(6,2) = -3**(0.5D0)*ETAH 
      Q2(6,5) = -3**(0.5D0)*MUH 
      Q2(6,6) = VT + 3*ABS(ETAH) + 3*ABS(MUH) 
      Q2(6,7) = 15**(0.5D0)*MUH 
      Q2(6,8) = 21**(0.5D0)*ABS(MUH) 
      Q2(6,10) = 15**(0.5D0)*ETAH 
      Q2(6,14) = 21**(0.5D0)*ABS(ETAH) 

      Q2(7,3) = -3**(0.5D0)*ETAH 
      Q2(7,5) = 5**(0.5D0)*ABS(MUH) 
      Q2(7,6) = -15**(0.5D0)*MUH 
      Q2(7,7) = VT + 3*ABS(ETAH) + 5*ABS(MUH) 
      Q2(7,8) = 35**(0.5D0)*MUH 
      Q2(7,11) = 2*15**(0.5D0)*ETAH - 3**(0.5D0)*5**(0.5D0)*
     >   ETAH 
      Q2(7,15) = 3**(0.5D0)*7**(0.5D0)*ABS(ETAH) 

      Q2(8,4) = -3**(0.5D0)*ETAH 
      Q2(8,5) = -7**(0.5D0)*MUH 
      Q2(8,6) = 21**(0.5D0)*ABS(MUH) 
      Q2(8,7) = -35**(0.5D0)*MUH 
      Q2(8,8) = VT + 3*ABS(ETAH) + 7*ABS(MUH) 
      Q2(8,12) = 2*15**(0.5D0)*ETAH - 3**(0.5D0)*5**(0.5D0)*
     >   ETAH 
      Q2(8,16) = 3**(0.5D0)*7**(0.5D0)*ABS(ETAH) 

      Q2(9,1) = 5**(0.5D0)*ABS(ETAH) 
      Q2(9,5) = -15**(0.5D0)*ETAH 
      Q2(9,9) = VT + 5*ABS(ETAH) + ABS(MUH) 
      Q2(9,10) = 3**(0.5D0)*MUH 
      Q2(9,11) = 5**(0.5D0)*ABS(MUH) 
      Q2(9,12) = 7**(0.5D0)*MUH 
      Q2(9,13) = 35**(0.5D0)*ETAH 

      Q2(10,2) = 5**(0.5D0)*ABS(ETAH) 
      Q2(10,6) = -15**(0.5D0)*ETAH 
      Q2(10,9) = -3**(0.5D0)*MUH 
      Q2(10,10) = VT + 5*ABS(ETAH) + 3*ABS(MUH) 
      Q2(10,11) = 2*15**(0.5D0)*MUH - 3**(0.5D0)*5**(0.5D0)*
     >   MUH 
      Q2(10,12) = 3**(0.5D0)*7**(0.5D0)*ABS(MUH) 
      Q2(10,14) = 35**(0.5D0)*ETAH 

      Q2(11,3) = 5**(0.5D0)*ABS(ETAH) 
      Q2(11,7) = -15**(0.5D0)*ETAH 
      Q2(11,9) = 5**(0.5D0)*ABS(MUH) 
      Q2(11,10) = -15**(0.5D0)*MUH 
      Q2(11,11) = VT + 5*ABS(ETAH) + 5*ABS(MUH) 
      Q2(11,12) = 35**(0.5D0)*MUH 
      Q2(11,15) = 35**(0.5D0)*ETAH 

      Q2(12,4) = 5**(0.5D0)*ABS(ETAH) 
      Q2(12,8) = -3**(0.5D0)*5**(0.5D0)*ETAH 
      Q2(12,9) = -7**(0.5D0)*MUH 
      Q2(12,10) = 3**(0.5D0)*7**(0.5D0)*ABS(MUH) 
      Q2(12,11) = -5**(0.5D0)*7**(0.5D0)*MUH 
      Q2(12,12) = VT + 5*ABS(ETAH) + 7*ABS(MUH) 
      Q2(12,16) = 2*35**(0.5D0)*ETAH - 5**(0.5D0)*7**(0.5D0)*
     >   ETAH 

      Q2(13,1) = -7**(0.5D0)*ETAH 
      Q2(13,5) = 21**(0.5D0)*ABS(ETAH) 
      Q2(13,9) = -35**(0.5D0)*ETAH 
      Q2(13,13) = VT + 7*ABS(ETAH) + ABS(MUH) 
      Q2(13,14) = 3**(0.5D0)*MUH 
      Q2(13,15) = 5**(0.5D0)*ABS(MUH) 
      Q2(13,16) = 7**(0.5D0)*MUH 

      Q2(14,2) = -7**(0.5D0)*ETAH 
      Q2(14,6) = 21**(0.5D0)*ABS(ETAH) 
      Q2(14,10) = -35**(0.5D0)*ETAH 
      Q2(14,13) = -3**(0.5D0)*MUH 
      Q2(14,14) = VT + 7*ABS(ETAH) + 3*ABS(MUH) 
      Q2(14,15) = 2*15**(0.5D0)*MUH - 3**(0.5D0)*5**(0.5D0)*
     >   MUH 
      Q2(14,16) = 3**(0.5D0)*7**(0.5D0)*ABS(MUH) 

      Q2(15,3) = -7**(0.5D0)*ETAH 
      Q2(15,7) = 3**(0.5D0)*7**(0.5D0)*ABS(ETAH) 
      Q2(15,11) = -5**(0.5D0)*7**(0.5D0)*ETAH 
      Q2(15,13) = 5**(0.5D0)*ABS(MUH) 
      Q2(15,14) = -3**(0.5D0)*5**(0.5D0)*MUH 
      Q2(15,15) = VT + 7*ABS(ETAH) + 5*ABS(MUH) 
      Q2(15,16) = 2*35**(0.5D0)*MUH - 5**(0.5D0)*7**(0.5D0)*
     >   MUH 

      Q2(16,4) = -7**(0.5D0)*ETAH 
      Q2(16,8) = 21**(0.5D0)*ABS(ETAH) 
      Q2(16,12) = -35**(0.5D0)*ETAH 
      Q2(16,13) = -7**(0.5D0)*MUH 
      Q2(16,14) = 21**(0.5D0)*ABS(MUH) 
      Q2(16,15) = -35**(0.5D0)*MUH 
      Q2(16,16) = VT + 7*ABS(ETAH) + 7*ABS(MUH) 

      Q2(1,17)  =  Q(1)*VOL(I2,J2,J,I) 
      Q2(2,17)  =  Q(2)*VOL(I2,J2,J,I) 
      Q2(3,17)  =  Q(3)*VOL(I2,J2,J,I) 
      Q2(4,17)  =  Q(4)*VOL(I2,J2,J,I) 
      Q2(5,17)  =  Q(5)*VOL(I2,J2,J,I) 
      Q2(6,17)  =  Q(6)*VOL(I2,J2,J,I) 
      Q2(7,17)  =  Q(7)*VOL(I2,J2,J,I) 
      Q2(8,17)  =  Q(8)*VOL(I2,J2,J,I) 
      Q2(9,17)  =  Q(9)*VOL(I2,J2,J,I) 
      Q2(10,17) = Q(10)*VOL(I2,J2,J,I) 
      Q2(11,17) = Q(11)*VOL(I2,J2,J,I) 
      Q2(12,17) = Q(12)*VOL(I2,J2,J,I) 
      Q2(13,17) = Q(13)*VOL(I2,J2,J,I) 
      Q2(14,17) = Q(14)*VOL(I2,J2,J,I) 
      Q2(15,17) = Q(15)*VOL(I2,J2,J,I) 
      Q2(16,17) = Q(16)*VOL(I2,J2,J,I) 

      Q2(1,17) = Q2(1,17) + ETAH*XNJ(1)*SIGN(1.0,ETAH) + 
     >   MUH*XNI(1,J2)*SIGN(1.0,MUH)  
      Q2(2,17) = Q2(2,17) + ETAH*XNJ(2)*SIGN(1.0,ETAH) - 
     >   3**(0.5D0)*MUH*XNI(1,J2) 
      Q2(3,17) = Q2(3,17) + ETAH*XNJ(3)*SIGN(1.0,ETAH) + 
     >   5**(0.5D0)*MUH*XNI(1,J2)*SIGN(1.0,MUH)  
      Q2(4,17) = Q2(4,17) + ETAH*XNJ(4)*SIGN(1.0,ETAH) - 
     >   7**(0.5D0)*MUH*XNI(1,J2) 
      Q2(5,17) = Q2(5,17) + MUH*XNI(2,J2)*SIGN(1.0,MUH)  - 
     >   3**(0.5D0)*ETAH*XNJ(1) 
      Q2(6,17) = Q2(6,17) - 3**(0.5D0)*ETAH*XNJ(2) - 3**(0.5D0)*
     >   MUH*XNI(2,J2) 
      Q2(7,17) = Q2(7,17) + 5**(0.5D0)*MUH*XNI(2,J2)*
     >   SIGN(1.0,MUH) - 3**(0.5D0)*ETAH*XNJ(3) 
      Q2(8,17) = Q2(8,17) - 3**(0.5D0)*ETAH*XNJ(4) - 7**(0.5D0)*
     >   MUH*XNI(2,J2) 
      Q2(9,17) = Q2(9,17) + MUH*XNI(3,J2)*SIGN(1.0,MUH)  + 
     >   5**(0.5D0)*ETAH*XNJ(1)*SIGN(1.0,ETAH) 
      Q2(10,17) = Q2(10,17) + 5**(0.5D0)*ETAH*XNJ(2)*SIGN(1.0,ETAH)
     >    - 3**(0.5D0)*MUH*XNI(3,J2) 
      Q2(11,17) = Q2(11,17) + 5**(0.5D0)*ETAH*XNJ(3)*SIGN(1.0,ETAH)
     >    + 5**(0.5D0)*MUH*XNI(3,J2)*SIGN(1.0,MUH)  
      Q2(12,17) = Q2(12,17) + 5**(0.5D0)*ETAH*XNJ(4)*SIGN(1.0,ETAH)
     >    - 7**(0.5D0)*MUH*XNI(3,J2) 
      Q2(13,17) = Q2(13,17) + MUH*XNI(4,J2)*SIGN(1.0,MUH) - 
     >   7**(0.5D0)*ETAH*XNJ(1) 
      Q2(14,17) = Q2(14,17) - 7**(0.5D0)*ETAH*XNJ(2) - 3**(0.5D0
     >   )*MUH*XNI(4,J2) 
      Q2(15,17) = Q2(15,17) + 5**(0.5D0)*MUH*XNI(4,J2)*
     >   SIGN(1.0,MUH) - 7**(0.5D0)*ETAH*XNJ(3) 
      Q2(16,17) = Q2(16,17) - 7**(0.5D0)*ETAH*XNJ(4) - 7**(0.5D0
     >   )*MUH*XNI(4,J2) 
      ENDIF
*
      CALL ALSBD(IELEM**2,1,Q2,IER,IELEM**2)
      IF(IER.NE.0) CALL XABORT('SNFDH2: SINGULAR MATRIX.')
*
************************************************************************
************************************************************************
      IF(JL.LT.ISPLH)THEN
      IF(IELEM.EQ.1) THEN
         XNJ(1)  =Q2(1,2)
      ELSEIF(IELEM.EQ.2) THEN
         XNJ(1)  =Q2(1,5) + SIGN(1.0,ETAH)*Q2(3,5)*(3**(0.5D0))
         XNJ(2)  =Q2(2,5) + SIGN(1.0,ETAH)*Q2(4,5)*(3**(0.5D0))
      ELSEIF(IELEM.EQ.3) THEN
         XNJ(1)  =Q2(1,10)+SIGN(1.0,ETAH)*(3**(0.5D0))*Q2(4,10)+
     >      (5**(0.5D0))*Q2(7,10)
         XNJ(2)  =Q2(2,10)+SIGN(1.0,ETAH)*(3**(0.5D0))*Q2(5,10)+
     >      (5**(0.5D0))*Q2(8,10)
         XNJ(3)  =Q2(3,10)+SIGN(1.0,ETAH)*(3**(0.5D0))*Q2(6,10)+
     >      (5**(0.5D0))*Q2(9,10)
      ELSEIF(IELEM.EQ.4)THEN
         XNJ(1)=Q2(1,17)+SIGN(1.0,ETAH)*(3**(0.5D0))*Q2(5,17)+
     >    (5**(0.5D0))*Q2(9,17)+SIGN(1.0,ETAH)*(7**(0.5D0))*Q2(13,17)
         XNJ(2)=Q2(2,17)+SIGN(1.0,ETAH)*(3**(0.5D0))*Q2(6,17)+
     >    (5**(0.5D0))*Q2(10,17)+SIGN(1.0,ETAH)*(7**(0.5D0))*Q2(14,17)
         XNJ(3)=Q2(3,17)+SIGN(1.0,ETAH)*(3**(0.5D0))*Q2(7,17)+
     >    (5**(0.5D0))*Q2(11,17)+SIGN(1.0,ETAH)*(7**(0.5D0))*Q2(15,17)
         XNJ(4)=Q2(4,17)+SIGN(1.0,ETAH)*(3**(0.5D0))*Q2(8,17)+
     >    (5**(0.5D0))*Q2(12,17)+SIGN(1.0,ETAH)*(7**(0.5D0))*Q2(16,17)
      ENDIF
      ELSEIF((JL.EQ.ISPLH).AND.(IHEXJ.LE.NHEX))THEN
      I3=I2
      C1=1.0D0
      IF((J.EQ.1).AND.(ILOZJ.EQ.3)) THEN
         I3=ISPLH+1-I2
         C1=-1.0D0
      ENDIF
      IF((IELEM.EQ.1).AND.(CHEXy.NE.0)) THEN
      BFLUX(CSIDEy,CHEXy,I3,1,IPQD) =Q2(1,2)
         WRITE(*,*) 'SNFDH2: bflux 1', IHEXJ, NHEX!,BFLUX!(CSIDEy,CHEXy,I3,1,IPQD)
      ELSEIF(IELEM.EQ.2) THEN
      BFLUX(CSIDEy,CHEXy,I3,1,IPQD) =Q2(1,5) + 
     >      SIGN(1.0,ETAH)*Q2(3,5)*(3**(0.5D0))
      BFLUX(CSIDEy,CHEXy,I3,2,IPQD) =(Q2(2,5) + 
     >      SIGN(1.0,ETAH)*Q2(4,5)*(3**(0.5D0)))*C1
      ELSEIF(IELEM.EQ.3) THEN
      BFLUX(CSIDEy,CHEXy,I3,1,IPQD) =Q2(1,10)+
     >      SIGN(1.0,ETAH)*(3**(0.5D0))*Q2(4,10)+(5**(0.5D0))*Q2(7,10)
      BFLUX(CSIDEy,CHEXy,I3,2,IPQD) =(Q2(2,10)+
     >      SIGN(1.0,ETAH)*(3**(0.5D0))*Q2(5,10)+(5**(0.5D0))*Q2(8,10)
     >      )*C1
      BFLUX(CSIDEy,CHEXy,I3,3,IPQD) =Q2(3,10)+
     >      SIGN(1.0,ETAH)*(3**(0.5D0))*Q2(6,10)+(5**(0.5D0))*Q2(9,10)
      ELSEIF(IELEM.EQ.4)THEN
      BFLUX(CSIDEy,CHEXy,I3,1,IPQD) =Q2(1,17)+
     >      SIGN(1.0,ETAH)*(3**(0.5D0))*Q2(5,17)+(5**(0.5D0))*Q2(9,17)+
     >      SIGN(1.0,ETAH)*(7**(0.5D0))*Q2(13,17)
      BFLUX(CSIDEy,CHEXy,I3,2,IPQD) =(Q2(2,17)+
     >      SIGN(1.0,ETAH)*(3**(0.5D0))*Q2(6,17)+(5**(0.5D0))*Q2(10,17)
     >      +SIGN(1.0,ETAH)*(7**(0.5D0))*Q2(14,17))*C1
      BFLUX(CSIDEy,CHEXy,I3,3,IPQD) =Q2(3,17)+
     >      SIGN(1.0,ETAH)*(3**(0.5D0))*Q2(7,17)+(5**(0.5D0))*Q2(11,17)
     >      +SIGN(1.0,ETAH)*(7**(0.5D0))*Q2(15,17)
      BFLUX(CSIDEy,CHEXy,I3,4,IPQD) =(Q2(4,17)+
     >      SIGN(1.0,ETAH)*(3**(0.5D0))*Q2(8,17)+(5**(0.5D0))*Q2(12,17)
     >      +SIGN(1.0,ETAH)*(7**(0.5D0))*Q2(16,17))*C1
      ENDIF
      ENDIF
* -------------------------------------------------------------------- *
      IF(IL.LT.ISPLH)THEN
      IF(IELEM.EQ.1) THEN
         XNI(1,J2)=Q2(1,2)
      ELSEIF(IELEM.EQ.2) THEN
         XNI(1,J2)=Q2(1,5) + SIGN(1.0,MUH)*Q2(2,5)*(3**(0.5D0))
         XNI(2,J2)=Q2(3,5) + SIGN(1.0,MUH)*Q2(4,5)*(3**(0.5D0))
      ELSEIF(IELEM.EQ.3) THEN
         XNI(1,J2)=Q2(1,10)+SIGN(1.0,MUH)*(3**(0.5D0))*Q2(2,10)+
     >      (5**(0.5D0))*Q2(3,10)
         XNI(2,J2)=Q2(4,10)+SIGN(1.0,MUH)*(3**(0.5D0))*Q2(5,10)+
     >      (5**(0.5D0))*Q2(6,10)
         XNI(3,J2)=Q2(7,10)+SIGN(1.0,MUH)*(3**(0.5D0))*Q2(8,10)+
     >      (5**(0.5D0))*Q2(9,10)
      ELSEIF(IELEM.EQ.4)THEN
         XNI(1,J2)=Q2(1,17)+SIGN(1.0,MUH)*(3**(0.5D0))*Q2(2,17)+
     >    (5**(0.5D0))*Q2(3,17)+SIGN(1.0,MUH)*(7**(0.5D0))*Q2(4,17)
         XNI(2,J2)=Q2(5,17)+SIGN(1.0,MUH)*(3**(0.5D0))*Q2(6,17)+
     >    (5**(0.5D0))*Q2(7,17)+SIGN(1.0,MUH)*(7**(0.5D0))*Q2(8,17)
         XNI(3,J2)=Q2(9,17)+SIGN(1.0,MUH)*(3**(0.5D0))*Q2(10,17)+
     >    (5**(0.5D0))*Q2(11,17)+SIGN(1.0,MUH)*(7**(0.5D0))*Q2(12,17)
         XNI(4,J2)=Q2(13,17)+SIGN(1.0,MUH)*(3**(0.5D0))*Q2(14,17)+
     >    (5**(0.5D0))*Q2(15,17)+SIGN(1.0,MUH)*(7**(0.5D0))*Q2(16,17)
      ENDIF
      ELSEIF((IL.EQ.ISPLH).AND.(IHEXI.LE.NHEX))THEN
      J3=J2
      C1=1.0D0
      IF((J.EQ.3).AND.(ILOZI.EQ.1)) THEN
C          J3=ISPLH+1-J2
         C1=-1.0D0
      ENDIF
      IF((IELEM.EQ.1).AND.(CHEXx.NE.0)) THEN
      BFLUX(CSIDEx,CHEXx,J3,1,IPQD) = Q2(1,2)
         WRITE(*,*) 'SNFDH2: bflux 2', IHEXI, NHEX!,BFLUX!(CSIDEx,CHEXx,J3,1,IPQD)
      ELSEIF(IELEM.EQ.2) THEN
      BFLUX(CSIDEx,CHEXx,J3,1,IPQD)=Q2(1,5) + 
     >      SIGN(1.0,MUH)*Q2(2,5)*(3**(0.5D0))
      BFLUX(CSIDEx,CHEXx,J3,2,IPQD)=(Q2(3,5) + 
     >      SIGN(1.0,MUH)*Q2(4,5)*(3**(0.5D0)))*C1
      ELSEIF(IELEM.EQ.3) THEN
      BFLUX(CSIDEx,CHEXx,J3,1,IPQD)=Q2(1,10)+
     >      SIGN(1.0,MUH)*(3**(0.5D0))*Q2(2,10)+(5**(0.5D0))*Q2(3,10)
      BFLUX(CSIDEx,CHEXx,J3,2,IPQD)=(Q2(4,10)+
     >      SIGN(1.0,MUH)*(3**(0.5D0))*Q2(5,10)+(5**(0.5D0))*Q2(6,10)
     >      )*C1
      BFLUX(CSIDEx,CHEXx,J3,3,IPQD)=Q2(7,10)+
     >      SIGN(1.0,MUH)*(3**(0.5D0))*Q2(8,10)+(5**(0.5D0))*Q2(9,10)
      ELSEIF(IELEM.EQ.4)THEN
      BFLUX(CSIDEx,CHEXx,J3,1,IPQD)=Q2(01,17)+
     >      SIGN(1.0,MUH)*(3**(0.5D0))*Q2(02,17)+(5**(0.5D0))*Q2(03,17)
     >      +SIGN(1.0,MUH)*(7**(0.5D0))*Q2(04,17)
      BFLUX(CSIDEx,CHEXx,J3,2,IPQD)=(Q2(05,17)+
     >      SIGN(1.0,MUH)*(3**(0.5D0))*Q2(06,17)+(5**(0.5D0))*Q2(07,17)
     >      +SIGN(1.0,MUH)*(7**(0.5D0))*Q2(08,17))*C1
      BFLUX(CSIDEx,CHEXx,J3,3,IPQD)=Q2(09,17)+
     >      SIGN(1.0,MUH)*(3**(0.5D0))*Q2(10,17)+(5**(0.5D0))*Q2(11,17)
     >      +SIGN(1.0,MUH)*(7**(0.5D0))*Q2(12,17)
      BFLUX(CSIDEx,CHEXx,J3,4,IPQD)=(Q2(13,17)+
     >      SIGN(1.0,MUH)*(3**(0.5D0))*Q2(14,17)+(5**(0.5D0))*Q2(15,17)
     >      +SIGN(1.0,MUH)*(7**(0.5D0))*Q2(16,17))*C1
      ENDIF
      ENDIF
*
      DO 150 P=1,NSCT
      DO 140 IEL=1,IELEM**2
      FLUX(IEL,P) = FLUX(IEL,P) + Q2(IEL,IELEM**2+1)*DN(P,M)
  140 CONTINUE
  150 CONTINUE

      ILOZ=((J-1)*ISPLH+(J2-1))*ISPLH+I2
      FLUX_G(:,:,ILOZ,I,IG)=FLUX_G(:,:,ILOZ,I,IG)+FLUX(:,:)
*--
  160 CONTINUE ! JL
  170 CONTINUE ! IL

  180 CONTINUE ! lozenges

  200 CONTINUE ! directions
  210 CONTINUE ! energy
*$OMP END PARALLEL DO

  220 CONTINUE ! macrocells

  230 CONTINUE ! wavefront number

  240 CONTINUE ! dodecant
************************************************************************
************************************************************************
      DO 777 IG=1,NGEFF
        IF(.NOT.INCONV(IG)) GO TO 777
        FUNKNO(:L4,IG)=
     1  RESHAPE(REAL(FLUX_G(:IELEM**2,:NSCT,:3*ISPLH**2,:NHEX,IG)),
     2  (/ L4 /) )
  777 CONTINUE
************************************************************************
************************************************************************
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(FLUX_G,FLUX,BFLUX,INDANG)
      RETURN
  400 FORMAT(16H SNFDH2: thread=,I8,12H --->(group=,I4,7H angle=,I4,
     1 1H))
      END
