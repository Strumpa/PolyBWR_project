*DECK SNFTH2
      SUBROUTINE SNFTH2(NUN,NGEFF,IMPX,INCONV,NGIND,NHEX,IELEM,ISPLH,
     1 NMAT,NPQ,NSCT,MAT,VOL,TOTAL,QEXT,DU,DE,W,DB,DA,MN,DN,SIDE,
     2 LOZSWP,COORDMAP,FUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform one inner iteration for solving SN equations in 2D hexagonal
* geometry for the DISCONTINUOUS GALERKIN method. Energy-angle
* multithreading. VOID boundary conditions.
*
*Copyright:
* Copyright (C) 2019 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. A. Calloo
*
*Parameters: input
* NUN      total number of unknowns in vector FUNKNO.
* NGEFF    number of energy groups processed in parallel.
* IMPX     print flag (equal to zero for no print).
* INCONV   energy group convergence flag (set to .FALSE. if converged).
* NGIND    energy group indices assign to the NGEFF set.
* NHEX     number of hexagons in X-Y plane.
* IELEM    measure of order of the spatial approximation polynomial:
*          =1 constant (classical diamond scheme - default for HODD);
*          =2 linear (default for DG);
*          =3 parabolic;
*          =4 cubic (only for DG).
* NMAT     number of material mixtures.
* ISPLH    splitting option for hexagons.
* NPQ      number of SN directions in 6 dodecants (including zero-weight
*          directions).
* NSCT     maximum number of spherical harmonics moments of the flux.
* MAT      material mixture index in each region.
* VOL      volumes of each region.
* TOTAL    macroscopic total cross sections.
* QEXT     Legendre components of the fixed source.
* DU       first direction cosines ($\\mu$).
* DE       second direction cosines ($\\eta$).
* W        weights.
* DB       diamond-scheme parameter.
* DA       diamond-scheme parameter.
* MN       moment-to-discrete matrix.
* DN       discrete-to-moment matrix.
* SIDE     side of an hexagon.
* LOZSWP   lozenge sweep order depending on direction.
* COORDMAP coordinate map - mapping the hexagons from the indices 
*          within the DRAGON geometry to a Cartesian axial coordinate
*          array (see redblobgames.com website).
*
*Parameters: output
* FUNKNO   Legendre components of the flux and boundary fluxes.
*
*Comments:
* 1. The direction of the axes I, J and D for the surface boundary 
*    fluxes are shown in the diagram below. This means that 
*    i) lozenge A has I- and D-boundaries (instead of I and J)
*    i) lozenge B has I- and J-boundaries
*    i) lozenge C has D- and J-boundaries (instead of I and J)
*
*                                  ^
*                         j-axis   |
*                                  |          ^
*                              _________     /    d-axis
*                             /       / \   /
*                            /   B   /   \
*                           /       /     \
*                          (-------(   A   )
*                           \       \     /
*                            \  C    \   / 
*                             \_______\_/   \
*                                            \   i-axis
*                                             ^
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
     > NSCT,MAT(ISPLH,ISPLH,3,NHEX),LOZSWP(3,6),COORDMAP(3,NHEX)
      LOGICAL INCONV(NGEFF)
      REAL SIDE,VOL(ISPLH,ISPLH,3,NHEX),TOTAL(0:NMAT,NGEFF),
     > QEXT(NUN,NGEFF),DU(NPQ),DB(ISPLH,ISPLH,3,NHEX,NPQ),
     > DA(ISPLH,ISPLH,3,NHEX,NPQ),MN(NPQ,NSCT),DN(NSCT,NPQ),DE(NPQ),
     > W(NPQ),FUNKNO(NUN,NGEFF)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(IUNOUT=6,PI=3.141592654)
      INTEGER :: NPQD(6),IIND(6),P,DCOORD
      DOUBLE PRECISION :: Q(IELEM**2), Q2(IELEM**2,(IELEM**2)+1),
     > VT, THETA, XNI(IELEM), XNJ(IELEM)
      REAL :: JAC(2,2,3), MUH, ETAH, AAA, BBB, CCC, DDD, MUHTEMP,
     > ETAHTEMP
      LOGICAL :: LHEX(NHEX)
*----
*  LOCAL ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: INDANG
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: TMPMAT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: FLUX
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: FLUX_G
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: TMPXNI,
     > TMPXNJ, TMPXND
*----
*  MAP MATERIAL VALUES TO CARTESIAN AXIAL COORDINATE MAP
*----
      NRINGS=INT((SQRT(  REAL((4*NHEX-1)/3)  )+1.)/2.)
      NCOL=2*NRINGS -1
      ALLOCATE(TMPMAT(ISPLH,ISPLH,3,NCOL,NCOL))
      TMPMAT(:,:,:,:,:) = -1
      DO IHEX_DOM=1,NHEX
         TMPMAT(:,:,:,COORDMAP(1,IHEX_DOM),COORDMAP(2,IHEX_DOM)) = 
     >         MAT(:,:,:,IHEX_DOM)
      ENDDO
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(INDANG(NPQ,6))
      ALLOCATE(FLUX(IELEM**2,NSCT,3*ISPLH**2,NHEX))
      ALLOCATE(FLUX_G(IELEM**2,NSCT,3*ISPLH**2,NHEX,NGEFF))
      ALLOCATE(TMPXNI(IELEM,ISPLH,NCOL))
      ALLOCATE(TMPXNJ(IELEM,ISPLH,NCOL))
      ALLOCATE(TMPXND(IELEM,ISPLH,NCOL))
*----
*  CONSTRUCT JACOBIAN MATRIX FOR EACH LOZENGE
*----
      JAC = RESHAPE((/ 1., -SQRT(3.), 1., SQRT(3.), 2., 0., 1.,
     >    SQRT(3.), 2., 0., -1., SQRT(3.) /), SHAPE(JAC))
      JAC = (SIDE/2.)*JAC
*----
*  DEFINITION OF CONSTANTS
*----
      L4=3*(IELEM**2)*(ISPLH**2)*NHEX*NSCT
*----
*  PARAMETER VALIDATION
*----
      IF((IELEM.LT.1).OR.(IELEM.GT.4))
     >   CALL XABORT('SNFTH2: INVALID IELEM (DIAM) VALUE. '
     >   //'CHECK INPUT DATA FILE.')

      FLUX_G(:IELEM**2,:NSCT,:3*ISPLH**2,:NHEX,:NGEFF)=0.0D0
*----
*  SET DODECANT SWAPPING ORDER
*----
      NPQD(:6)=0
      INDANG(:NPQ,:6)=0
      IIND(:6)=0
      DO M=1,NPQ
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
        ! Assume IIND(I)=I in hexagonal geometry
        IIND(IND)=IND
        NPQD(IND)=NPQD(IND)+1
        INDANG(NPQD(IND),IND)=M
      ENDDO
*----
*  MAIN LOOP OVER DODECANTS
*----
      DO JND=1,6
      IND=IIND(JND)
      ! Needed because of S2 LS (4 dir. for 6 sextants)
      IF(IND.EQ.0) CYCLE
*----
*  MAIN SWAPPING LOOPS FOR SN FLUX CALCULATION
*  LOOP OVER ENERGY AND ANGLES
*----
*$OMP  PARALLEL DO
*$OMP1 PRIVATE(ITID,I,J,LHEX,IHEX,IIHEX,DCOORD,IOF)
*$OMP1 PRIVATE(ILOZLOOP,ILOZ,AAA,BBB,CCC,DDD,IL,I2,JL,J2,FLUX)
*$OMP1 PRIVATE(XNI,XNJ,MUHTEMP,MUH,ETAHTEMP,ETAH,IEL,Q,P)
*$OMP1 PRIVATE(VT,Q2,IER,C1,I3,I_FETCH,M,III,JJJ,IIM,JIM)
*$OMP1 PRIVATE(TMPXNI,TMPXNJ,TMPXND)
*$OMP1 SHARED(IND)
*$OMP1 REDUCTION(+:FLUX_G) COLLAPSE(2)
      DO IG=1,NGEFF
      DO IPQD=1,NPQD(IND)

#if defined(_OPENMP)
        ITID=omp_get_thread_num()
#else
        ITID=0
#endif
      IF(IMPX.GT.5) WRITE(IUNOUT,400) ITID,NGIND(IG),IPQD
      M=INDANG(IPQD,IND)
      IF(.NOT.INCONV(IG)) CYCLE
      IF(W(M).EQ.0.0) CYCLE

      FLUX(:IELEM**2,:NSCT,:3*ISPLH**2,:NHEX)=0.0D0
      TMPXNI(:IELEM,:ISPLH,:NCOL)=0.0D0
      TMPXNJ(:IELEM,:ISPLH,:NCOL)=0.0D0
      TMPXND(:IELEM,:ISPLH,:NCOL)=0.0D0
*----
*  LOOP OVER CARTESIAN MAP OF HEXAGONAL DOMAIN
*----
      DO JJJ=1,NCOL
      JIM=JJJ
      ! Account for different sweep direction depending on angle 
      IF((IND.EQ.1).OR.(IND.EQ.2).OR.(IND.EQ.3)) JIM=NCOL+1-JIM
        
      DO III=1,NCOL
      IIM=III
      ! Account for different sweep direction depending on angle 
      IF((IND.EQ.2).OR.(IND.EQ.3).OR.(IND.EQ.4)) IIM=NCOL+1-IIM

      ! For IND 3 or 6, Cartesian axial coordinate map is swept 
      ! vertically instead of horizontally. IM suffix is for 'IMmutable'
      I=IIM
      J=JIM
      IF((IND.EQ.3).OR.(IND.EQ.6))THEN
        I=JIM
        J=IIM
      ENDIF

      ! If within corners of Cartesian axial coordinate map (where
      ! there are no hexagons), skip loop
      IF(TMPMAT(1,1,1,I,J).EQ.-1) CYCLE

      ! Find DRAGON geometry hexagonal index using I and J
      LHEX=(COORDMAP(1,:).EQ.I .AND. COORDMAP(2,:).EQ.J)
      IHEX=0
      DO IIHEX=1,NHEX
        IF(LHEX(IIHEX)) THEN
          IHEX=IIHEX
          EXIT
        ENDIF
      ENDDO
      IF(IHEX.EQ.0) CALL XABORT('SNFTH2: IHEX FAILURE.')
      ! Find D coordinate
      DCOORD = ABS(COORDMAP(3,IHEX))-NRINGS
*----
*  LOOP OVER LOZENGES
*----
      DO ILOZLOOP=1,3
      ILOZ=LOZSWP(ILOZLOOP,IND)

      ! Get Jacobian elements values
      AAA = JAC(1,1,ILOZ)
      BBB = JAC(1,2,ILOZ)
      CCC = JAC(2,1,ILOZ)
      DDD = JAC(2,2,ILOZ)
*----
*  LOOP OVER SUBMESH WITHIN EACH LOZENGE
*----
      DO IL=1,ISPLH
      I2=IL
      ! Account for different sweep direction depending on angle 
      IF((ILOZ.EQ.1).OR.(ILOZ.EQ.2))THEN
        IF((IND.EQ.2).OR.(IND.EQ.3).OR.(IND.EQ.4)) I2=ISPLH+1-I2
      ELSEIF(ILOZ.EQ.3)THEN
        IF((IND.EQ.3).OR.(IND.EQ.4).OR.(IND.EQ.5)) I2=ISPLH+1-I2
      ENDIF

      DO JL=1,ISPLH
      J2=JL
      ! Account for different sweep direction depending on angle 
      IF((ILOZ.EQ.2).OR.(ILOZ.EQ.3))THEN
        IF((IND.EQ.4).OR.(IND.EQ.5).OR.(IND.EQ.6)) J2=ISPLH+1-J2
      ELSEIF(ILOZ.EQ.1)THEN
        IF((IND.EQ.3).OR.(IND.EQ.4).OR.(IND.EQ.5)) J2=ISPLH+1-J2
      ENDIF

      ! Skip loop if virtual element 
      IF(MAT(I2,J2,ILOZ,IHEX).EQ.0) CYCLE

      ! Read in XNI and XNJ depending on lozenge
      I_FETCH=0
      IF((ILOZ.EQ.1))THEN
        ! Read boundary fluxes in reverse for lozenge A since affine
        ! transformation of lozenges causes the D and I directions
        ! of lozenges C and A respectively to be reversed
        I_FETCH=ISPLH+1-I2
        XNI(:) = TMPXNI(:,J2,J)
        XNJ(:) = TMPXND(:,I_FETCH,DCOORD)
      ELSEIF((ILOZ.EQ.2))THEN
        XNI(:) = TMPXNI(:,J2,J)
        XNJ(:) = TMPXNJ(:,I2,I)
      ELSEIF((ILOZ.EQ.3))THEN
        XNI(:) = TMPXND(:,J2,DCOORD)
        XNJ(:) = TMPXNJ(:,I2,I)
      ENDIF

      ! Compute adjusted direction cosines 
      MUHTEMP  =  DA(I2,J2,ILOZ,IHEX,M)
      ETAHTEMP =  DB(I2,J2,ILOZ,IHEX,M)
      MUH = (MUHTEMP*DDD) - (ETAHTEMP*BBB)
      ETAH = (-MUHTEMP*CCC) + (ETAHTEMP*AAA)

      ! Compute source and (vol*sigma) for current element 
      DO IEL=1,IELEM**2
      Q(IEL)=0.0D0
      DO P=1,NSCT
      IOF=((((((IHEX-1)*3+(ILOZ-1))*ISPLH+(J2-1))*ISPLH+
     >  (I2-1))*NSCT+(P-1))*IELEM**2)+IEL
      Q(IEL)=Q(IEL)+QEXT(IOF,IG)*MN(M,P)
      ENDDO
      ENDDO
      VT=VOL(I2,J2,ILOZ,IHEX)*TOTAL(MAT(I2,J2,ILOZ,IHEX),IG)

      ! Initialise and compute matrix of coefficients, and source.
      ! Note: the matrix does seem to be skew-symmetric (or 
      ! antisymmetric), i.e. A^T = -A. However, from a numerical 
      ! perspective, it might be more efficient to assign values 
      ! than to sweep over the half matrix.
      Q2(:IELEM**2,:IELEM**2+1)=0.0D0
      IF(IELEM.EQ.1) THEN
      Q2(1,1) = VT + ABS(ETAH) + ABS(MUH) 
      Q2(1,2) = Q(1)*VOL(I2,J2,ILOZ,IHEX) 

      Q2(1,2) = Q2(1,2) + ETAH*SIGN(1.0,ETAH)*XNJ(1) + 
     >   MUH*SIGN(1.0,MUH)*XNI(1) 
      ELSEIF(IELEM.EQ.2) THEN
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
      Q2(1,5) = Q(1)*VOL(I2,J2,ILOZ,IHEX)
      Q2(2,5) = Q(2)*VOL(I2,J2,ILOZ,IHEX)
      Q2(3,5) = Q(3)*VOL(I2,J2,ILOZ,IHEX)
      Q2(4,5) = Q(4)*VOL(I2,J2,ILOZ,IHEX)
*     jump terms (dependent on direction)
      Q2(1,5) = Q2(1,5) + ETAH*SIGN(1.0,ETAH)*XNJ(1) + 
     >   MUH*SIGN(1.0,MUH)*XNI(1) 
      Q2(2,5) = Q2(2,5) + ETAH*SIGN(1.0,ETAH)*XNJ(2) - 3**(0.
     >   5D0)*MUH*XNI(1) 
      Q2(3,5) = Q2(3,5) + MUH*SIGN(1.0,MUH)*XNI(2) - 3**(
     >   0.5D0)*ETAH*XNJ(1) 
      Q2(4,5) = Q2(4,5) - 3**(0.5D0)*ETAH*
     >   XNJ(2) - 3**(0.5D0)*MUH*XNI(2) 
*** ---------------------------------------------------------------- ***
      ELSEIF(IELEM.EQ.3) THEN
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
      Q2(1,10) = Q(1)*VOL(I2,J2,ILOZ,IHEX) 
      Q2(2,10) = Q(2)*VOL(I2,J2,ILOZ,IHEX) 
      Q2(3,10) = Q(3)*VOL(I2,J2,ILOZ,IHEX) 
      Q2(4,10) = Q(4)*VOL(I2,J2,ILOZ,IHEX) 
      Q2(5,10) = Q(5)*VOL(I2,J2,ILOZ,IHEX) 
      Q2(6,10) = Q(6)*VOL(I2,J2,ILOZ,IHEX) 
      Q2(7,10) = Q(7)*VOL(I2,J2,ILOZ,IHEX) 
      Q2(8,10) = Q(8)*VOL(I2,J2,ILOZ,IHEX) 
      Q2(9,10) = Q(9)*VOL(I2,J2,ILOZ,IHEX) 
      Q2(1,10) = Q2(1,10) + ETAH*XNJ(1)*SIGN(1.0,ETAH) + 
     >   MUH*XNI(1)*SIGN(1.0,MUH) 
      Q2(2,10) = Q2(2,10) + ETAH*XNJ(2)*SIGN(1.0,ETAH) - 3**(
     >   0.5D0)*MUH*XNI(1) 
      Q2(3,10) = Q2(3,10) + ETAH*XNJ(3)*SIGN(1.0,ETAH) + 5**(
     >   0.5D0)*MUH*XNI(1)*SIGN(1.0,MUH) 
      Q2(4,10) = Q2(4,10) + MUH*XNI(2)*SIGN(1.0,MUH) - 3*
     >   *(0.5D0)*ETAH*XNJ(1) 
      Q2(5,10) = Q2(5,10) - 3**(0.5D0)*ETAH
     >   *XNJ(2) - 3**(0.5D0)*MUH*XNI(2) 
      Q2(6,10) = Q2(6,10) + 5**(0.5D0)*MUH
     >   *XNI(2)*SIGN(1.0,MUH) - 3**(0.5D0)*ETAH*XNJ(3) 
      Q2(7,10) = Q2(7,10) + MUH*XNI(3)*SIGN(1.0,MUH) + 5*
     >   *(0.5D0)*ETAH*XNJ(1)*SIGN(1.0,ETAH) 
      Q2(8,10) = Q2(8,10) + 5**(0.5D0)*ETAH*SIGN(1.0,ETAH)*
     >   XNJ(2) - 3**(0.5D0)*MUH*XNI(3) 
      Q2(9,10) = Q2(9,10) + 5**(0.5D0)*ETAH*SIGN(1.0,ETAH)*
     >   XNJ(3) + 5**(0.5D0)*MUH*XNI(3)*SIGN(1.0,MUH) 
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

      Q2(1,17)  =  Q(1)*VOL(I2,J2,ILOZ,IHEX) 
      Q2(2,17)  =  Q(2)*VOL(I2,J2,ILOZ,IHEX) 
      Q2(3,17)  =  Q(3)*VOL(I2,J2,ILOZ,IHEX) 
      Q2(4,17)  =  Q(4)*VOL(I2,J2,ILOZ,IHEX) 
      Q2(5,17)  =  Q(5)*VOL(I2,J2,ILOZ,IHEX) 
      Q2(6,17)  =  Q(6)*VOL(I2,J2,ILOZ,IHEX) 
      Q2(7,17)  =  Q(7)*VOL(I2,J2,ILOZ,IHEX) 
      Q2(8,17)  =  Q(8)*VOL(I2,J2,ILOZ,IHEX) 
      Q2(9,17)  =  Q(9)*VOL(I2,J2,ILOZ,IHEX) 
      Q2(10,17) = Q(10)*VOL(I2,J2,ILOZ,IHEX) 
      Q2(11,17) = Q(11)*VOL(I2,J2,ILOZ,IHEX) 
      Q2(12,17) = Q(12)*VOL(I2,J2,ILOZ,IHEX) 
      Q2(13,17) = Q(13)*VOL(I2,J2,ILOZ,IHEX) 
      Q2(14,17) = Q(14)*VOL(I2,J2,ILOZ,IHEX) 
      Q2(15,17) = Q(15)*VOL(I2,J2,ILOZ,IHEX) 
      Q2(16,17) = Q(16)*VOL(I2,J2,ILOZ,IHEX) 

      Q2(1,17) = Q2(1,17) + ETAH*XNJ(1)*SIGN(1.0,ETAH) + 
     >   MUH*XNI(1)*SIGN(1.0,MUH)  
      Q2(2,17) = Q2(2,17) + ETAH*XNJ(2)*SIGN(1.0,ETAH) - 
     >   3**(0.5D0)*MUH*XNI(1) 
      Q2(3,17) = Q2(3,17) + ETAH*XNJ(3)*SIGN(1.0,ETAH) + 
     >   5**(0.5D0)*MUH*XNI(1)*SIGN(1.0,MUH)  
      Q2(4,17) = Q2(4,17) + ETAH*XNJ(4)*SIGN(1.0,ETAH) - 
     >   7**(0.5D0)*MUH*XNI(1) 
      Q2(5,17) = Q2(5,17) + MUH*XNI(2)*SIGN(1.0,MUH)  - 
     >   3**(0.5D0)*ETAH*XNJ(1) 
      Q2(6,17) = Q2(6,17) - 3**(0.5D0)*ETAH*XNJ(2) - 3**(0.5D0)*
     >   MUH*XNI(2) 
      Q2(7,17) = Q2(7,17) + 5**(0.5D0)*MUH*XNI(2)*
     >   SIGN(1.0,MUH) - 3**(0.5D0)*ETAH*XNJ(3) 
      Q2(8,17) = Q2(8,17) - 3**(0.5D0)*ETAH*XNJ(4) - 7**(0.5D0)*
     >   MUH*XNI(2) 
      Q2(9,17) = Q2(9,17) + MUH*XNI(3)*SIGN(1.0,MUH)  + 
     >   5**(0.5D0)*ETAH*XNJ(1)*SIGN(1.0,ETAH) 
      Q2(10,17) = Q2(10,17) + 5**(0.5D0)*ETAH*XNJ(2)*SIGN(1.0,ETAH)
     >    - 3**(0.5D0)*MUH*XNI(3) 
      Q2(11,17) = Q2(11,17) + 5**(0.5D0)*ETAH*XNJ(3)*SIGN(1.0,ETAH)
     >    + 5**(0.5D0)*MUH*XNI(3)*SIGN(1.0,MUH)  
      Q2(12,17) = Q2(12,17) + 5**(0.5D0)*ETAH*XNJ(4)*SIGN(1.0,ETAH)
     >    - 7**(0.5D0)*MUH*XNI(3) 
      Q2(13,17) = Q2(13,17) + MUH*XNI(4)*SIGN(1.0,MUH) - 
     >   7**(0.5D0)*ETAH*XNJ(1) 
      Q2(14,17) = Q2(14,17) - 7**(0.5D0)*ETAH*XNJ(2) - 3**(0.5D0
     >   )*MUH*XNI(4) 
      Q2(15,17) = Q2(15,17) + 5**(0.5D0)*MUH*XNI(4)*
     >   SIGN(1.0,MUH) - 7**(0.5D0)*ETAH*XNJ(3) 
      Q2(16,17) = Q2(16,17) - 7**(0.5D0)*ETAH*XNJ(4) - 7**(0.5D0
     >   )*MUH*XNI(4) 
      ENDIF

      ! Solve system of equation given by matrix of coefficients
      ! and source
      CALL ALSBD(IELEM**2,1,Q2,IER,IELEM**2)
      IF(IER.NE.0) CALL XABORT('SNFTH2: SINGULAR MATRIX.')

      ! Assign J-boundary fluxes if lozenges B or C
      IF((ILOZ.EQ.2).OR.(ILOZ.EQ.3))THEN
        IF(IELEM.EQ.1) THEN
          TMPXNJ(1,I2,I)=Q2(1,2)
        ELSEIF(IELEM.EQ.2) THEN
          TMPXNJ(1,I2,I)=Q2(1,5) + SIGN(1.0,ETAH)*Q2(3,5)*(3**(0.5D0))
          TMPXNJ(2,I2,I)=Q2(2,5) + SIGN(1.0,ETAH)*Q2(4,5)*(3**(0.5D0))
        ELSEIF(IELEM.EQ.3) THEN
          TMPXNJ(1,I2,I)=Q2(1,10)+SIGN(1.0,ETAH)*(3**(0.5D0))*Q2(4,10)+
     >      (5**(0.5D0))*Q2(7,10)
          TMPXNJ(2,I2,I)=Q2(2,10)+SIGN(1.0,ETAH)*(3**(0.5D0))*Q2(5,10)+
     >      (5**(0.5D0))*Q2(8,10)
          TMPXNJ(3,I2,I)=Q2(3,10)+SIGN(1.0,ETAH)*(3**(0.5D0))*Q2(6,10)+
     >      (5**(0.5D0))*Q2(9,10)
        ELSEIF(IELEM.EQ.4)THEN
         TMPXNJ(1,I2,I)=Q2(1,17)+SIGN(1.0,ETAH)*(3**(0.5D0))*Q2(5,17)+
     >    (5**(0.5D0))*Q2(9,17)+SIGN(1.0,ETAH)*(7**(0.5D0))*Q2(13,17)
         TMPXNJ(2,I2,I)=Q2(2,17)+SIGN(1.0,ETAH)*(3**(0.5D0))*Q2(6,17)+
     >    (5**(0.5D0))*Q2(10,17)+SIGN(1.0,ETAH)*(7**(0.5D0))*Q2(14,17)
         TMPXNJ(3,I2,I)=Q2(3,17)+SIGN(1.0,ETAH)*(3**(0.5D0))*Q2(7,17)+
     >    (5**(0.5D0))*Q2(11,17)+SIGN(1.0,ETAH)*(7**(0.5D0))*Q2(15,17)
         TMPXNJ(4,I2,I)=Q2(4,17)+SIGN(1.0,ETAH)*(3**(0.5D0))*Q2(8,17)+
     >    (5**(0.5D0))*Q2(12,17)+SIGN(1.0,ETAH)*(7**(0.5D0))*Q2(16,17)
        ENDIF
      ENDIF
      ! Assign I-boundary fluxes if lozenges A or B
      IF((ILOZ.EQ.1).OR.(ILOZ.EQ.2))THEN
        IF(IELEM.EQ.1) THEN
          TMPXNI(1,J2,J)=Q2(1,2)
        ELSEIF(IELEM.EQ.2) THEN
          TMPXNI(1,J2,J)=Q2(1,5) + SIGN(1.0,MUH)*Q2(2,5)*(3**(0.5D0))
          TMPXNI(2,J2,J)=Q2(3,5) + SIGN(1.0,MUH)*Q2(4,5)*(3**(0.5D0))
        ELSEIF(IELEM.EQ.3) THEN
          TMPXNI(1,J2,J)=Q2(1,10)+SIGN(1.0,MUH)*(3**(0.5D0))*Q2(2,10)+
     >      (5**(0.5D0))*Q2(3,10)
          TMPXNI(2,J2,J)=Q2(4,10)+SIGN(1.0,MUH)*(3**(0.5D0))*Q2(5,10)+
     >      (5**(0.5D0))*Q2(6,10)
          TMPXNI(3,J2,J)=Q2(7,10)+SIGN(1.0,MUH)*(3**(0.5D0))*Q2(8,10)+
     >      (5**(0.5D0))*Q2(9,10)
        ELSEIF(IELEM.EQ.4)THEN
          TMPXNI(1,J2,J)=Q2(1,17)+SIGN(1.0,MUH)*(3**(0.5D0))*Q2(2,17)
     >      +(5**(0.5D0))*Q2(3,17)+SIGN(1.0,MUH)*(7**(0.5D0))*Q2(4,17)
          TMPXNI(2,J2,J)=Q2(5,17)+SIGN(1.0,MUH)*(3**(0.5D0))*Q2(6,17)
     >      +(5**(0.5D0))*Q2(7,17)+SIGN(1.0,MUH)*(7**(0.5D0))*Q2(8,17)
          TMPXNI(3,J2,J)=Q2(9,17)+SIGN(1.0,MUH)*(3**(0.5D0))*Q2(10,17)
     >      +(5**(0.5D0))*Q2(11,17)+SIGN(1.0,MUH)*(7**(0.5D0))*
     >      Q2(12,17)
          TMPXNI(4,J2,J)=Q2(13,17)+SIGN(1.0,MUH)*(3**(0.5D0))*Q2(14,17)
     >      +(5**(0.5D0))*Q2(15,17)+SIGN(1.0,MUH)*(7**(0.5D0))*
     >      Q2(16,17)
        ENDIF
      ENDIF
      ! Assign D-boundary fluxes if lozenge A using XNJ
      IF((ILOZ.EQ.1))THEN
        C1=1.0D0
        I3=I_FETCH
        IF(JL.EQ.ISPLH) C1=-1.0D0
        IF(IELEM.EQ.1) THEN
          TMPXND(1,I3,DCOORD)=Q2(1,2)
        ELSEIF(IELEM.EQ.2) THEN
          TMPXND(1,I3,DCOORD)=Q2(1,5) + SIGN(1.0,ETAH)*Q2(3,5)*
     >      (3**(0.5D0))
          TMPXND(2,I3,DCOORD)=(Q2(2,5) + SIGN(1.0,ETAH)*Q2(4,5)*
     >      (3**(0.5D0)))*C1
        ELSEIF(IELEM.EQ.3) THEN
          TMPXND(1,I3,DCOORD)=Q2(1,10)+SIGN(1.0,ETAH)*(3**(0.5D0))*
     >      Q2(4,10)+(5**(0.5D0))*Q2(7,10)
          TMPXND(2,I3,DCOORD)=(Q2(2,10)+SIGN(1.0,ETAH)*(3**(0.5D0))*
     >      Q2(5,10)+(5**(0.5D0))*Q2(8,10))*C1
          TMPXND(3,I3,DCOORD)=Q2(3,10)+SIGN(1.0,ETAH)*(3**(0.5D0))*
     >      Q2(6,10)+(5**(0.5D0))*Q2(9,10)
        ELSEIF(IELEM.EQ.4)THEN
          TMPXND(1,I3,DCOORD)=Q2(1,17)+SIGN(1.0,ETAH)*(3**(0.5D0))*
     >      Q2(5,17)+(5**(0.5D0))*Q2(9,17)+SIGN(1.0,ETAH)*(7**(0.5D0))*
     >      Q2(13,17)
          TMPXND(2,I3,DCOORD)=(Q2(2,17)+SIGN(1.0,ETAH)*(3**(0.5D0))*
     >      Q2(6,17)+(5**(0.5D0))*Q2(10,17)+SIGN(1.0,ETAH)*
     >      (7**(0.5D0))*Q2(14,17))*C1
          TMPXND(3,I3,DCOORD)=Q2(3,17)+SIGN(1.0,ETAH)*(3**(0.5D0))*
     >      Q2(7,17)+(5**(0.5D0))*Q2(11,17)+SIGN(1.0,ETAH)*
     >      (7**(0.5D0))*Q2(15,17)
          TMPXND(4,I3,DCOORD)=(Q2(4,17)+SIGN(1.0,ETAH)*(3**(0.5D0))*
     >      Q2(8,17)+(5**(0.5D0))*Q2(12,17)+SIGN(1.0,ETAH)*
     >      (7**(0.5D0))*Q2(16,17))*C1
        ENDIF
      ENDIF
      ! Assign D-boundary fluxes if lozenge C using XNI
      IF((ILOZ.EQ.3))THEN
        C1=1.0D0
        IF(IL.EQ.ISPLH) C1=-1.0D0
        IF(IELEM.EQ.1) THEN
          TMPXND(1,J2,DCOORD)=Q2(1,2)
        ELSEIF(IELEM.EQ.2) THEN
          TMPXND(1,J2,DCOORD)=Q2(1,5) + SIGN(1.0,MUH)*Q2(2,5)*
     >      (3**(0.5D0))
          TMPXND(2,J2,DCOORD)=(Q2(3,5) + SIGN(1.0,MUH)*Q2(4,5)*
     >      (3**(0.5D0)))*C1
        ELSEIF(IELEM.EQ.3) THEN
          TMPXND(1,J2,DCOORD)=Q2(1,10)+SIGN(1.0,MUH)*(3**(0.5D0))*
     >      Q2(2,10)+(5**(0.5D0))*Q2(3,10)
          TMPXND(2,J2,DCOORD)=(Q2(4,10)+SIGN(1.0,MUH)*(3**(0.5D0))*
     >      Q2(5,10)+(5**(0.5D0))*Q2(6,10))*C1
          TMPXND(3,J2,DCOORD)=Q2(7,10)+SIGN(1.0,MUH)*(3**(0.5D0))*
     >      Q2(8,10)+(5**(0.5D0))*Q2(9,10)
        ELSEIF(IELEM.EQ.4)THEN
          TMPXND(1,J2,DCOORD)=Q2(1,17)+SIGN(1.0,MUH)*(3**(0.5D0))*
     >      Q2(2,17)+(5**(0.5D0))*Q2(3,17)+SIGN(1.0,MUH)*(7**(0.5D0))*
     >      Q2(4,17)
          TMPXND(2,J2,DCOORD)=(Q2(5,17)+SIGN(1.0,MUH)*(3**(0.5D0))*
     >      Q2(6,17)+(5**(0.5D0))*Q2(7,17)+SIGN(1.0,MUH)*(7**(0.5D0))*
     >      Q2(8,17))*C1
          TMPXND(3,J2,DCOORD)=Q2(9,17)+SIGN(1.0,MUH)*(3**(0.5D0))*
     >      Q2(10,17)+(5**(0.5D0))*Q2(11,17)+SIGN(1.0,MUH)*(7**(0.5D0))*
     >      Q2(12,17)
          TMPXND(4,J2,DCOORD)=(Q2(13,17)+SIGN(1.0,MUH)*(3**(0.5D0))*
     >      Q2(14,17)+(5**(0.5D0))*Q2(15,17)+SIGN(1.0,MUH)*(7**(0.5D0))*
     >      Q2(16,17))*C1
        ENDIF
      ENDIF

      ! Compute flux
      IOF=((ILOZ-1)*ISPLH+(J2-1))*ISPLH+I2
      DO P=1,NSCT
      DO IEL=1,IELEM**2
      FLUX(IEL,P,IOF,IHEX) = FLUX(IEL,P,IOF,IHEX) +
     >   Q2(IEL,IELEM**2+1)*DN(P,M)
      ENDDO
      ENDDO

      ENDDO ! jl, isplh
      ENDDO ! il, isplh

      ENDDO ! iloz

      ENDDO ! iii, ncol
      ENDDO ! jjj, ncol

      FLUX_G(:,:,:,:,IG)=FLUX_G(:,:,:,:,IG)+FLUX(:,:,:,:)

      ENDDO ! ipqd
      ENDDO ! ig
*$OMP END PARALLEL DO

      ENDDO ! jnd

      DO IG=1,NGEFF
        IF(.NOT.INCONV(IG)) CYCLE
        FUNKNO(:L4,IG)=
     >  RESHAPE(REAL(FLUX_G(:IELEM**2,:NSCT,:3*ISPLH**2,:NHEX,IG)),
     >  (/ L4 /) )
      ENDDO
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(FLUX_G,FLUX,INDANG,TMPXNI,TMPXNJ,TMPXND,TMPMAT)
      RETURN
  400 FORMAT(16H SNFTH1: thread=,I8,12H --->(group=,I4,7H angle=,I4,1H))
      END
