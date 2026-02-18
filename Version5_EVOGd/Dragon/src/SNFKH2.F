*DECK SNFKH2
      SUBROUTINE SNFKH2(NUN,NGEFF,IMPX,INCONV,NGIND,NHEX,ISPLH,SIDE,
     1 IELEM,NM,NMX,NMY,NMAT,NPQ,NSCT,MAT,VOL,TOTAL,QEXT,LFIXUP,DU,DE,W,
     2 DB,DA,MN,DN,WX,WY,CST,ISADPT,LOZSWP,COORDMAP,FUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform one inner iteration for solving SN equations in 2D hexagonal
* geometry for the HODD/DG method. KBA-like multithreading, i.e.,
* macrocell-energy. VOID boundary conditions. Boltzmann (BTE) 
* discretization.

*
*Copyright:
* Copyright (C) 2025 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert, A. A. Calloo and C. Bienvenue
*
*Parameters: input
* NUN      total number of unknowns in vector FUNKNO.
* NGEFF    number of energy groups processed in parallel.
* IMPX     print flag (equal to zero for no print).
* INCONV   energy group convergence flag (set to .FALSE. if converged).
* NGIND    energy group indices assign to the NGEFF set.
* NHEX     number of hexagons in X-Y plane.
* ISPLH    splitting option for hexagons.
* SIDE     side of an hexagon.
* IELEM    measure of order of the spatial approximation polynomial:
*          =1 constant - default for HODD;
*          =2 linear - default for DG;
*          >3 higher orders.
* NM       number of moments in space and energy for flux components
* NMX      number of moments for X axis boundaries components
* NMY      number of moments for Y axis boundaries components
* NMAT     number of material mixtures.
* NPQ      number of SN directions in four octants (including zero-weight
*          directions).
* NSCT     maximum number of spherical harmonics moments of the flux.
* MAT      material mixture index in each region.
* VOL      volumes of each region.
* TOTAL    macroscopic total cross sections.
* QEXT     Legendre components of the fixed source.
* LFIXUP   flag to enable negative flux fixup.
* DU       first direction cosines ($\\mu$).
* DE       second direction cosines ($\\eta$).
* W        weights.
* DB       diamond-scheme parameter.
* DA       diamond-scheme parameter.
* MN       moment-to-discrete matrix.
* DN       discrete-to-moment matrix.
* WX       spatial X axis closure relation weighting factors.
* WY       spatial Y axis closure relation weighting factors.
* CST      constants for the polynomial approximations.
* ISADPT   flag to enable/disable adaptive flux calculations.
* LOZSWP   lozenge sweep order depending on direction.
* COORDMAP coordinate map - mapping the hexagons from the indices 
*          within the DRAGON geometry to a Cartesian axial coordinate
*          array (see redblobgames.com website).
*
*Parameters: input/output
* FUNKNO  Legendre components of the flux and boundary fluxes.
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
#if defined(_OPENMP)
      USE omp_lib
#endif
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NUN,NGEFF,IMPX,NGIND(NGEFF),NHEX,ISPLH,IELEM,NM,NMX,
     1 NMY,NMAT,NPQ,NSCT,MAT(ISPLH,ISPLH,3,NHEX),LOZSWP(3,6),
     2 COORDMAP(3,NHEX)
      LOGICAL INCONV(NGEFF)
      REAL SIDE,VOL(ISPLH,ISPLH,3,NHEX),TOTAL(0:NMAT,NGEFF),
     1 QEXT(NUN,NGEFF),DU(NPQ),DE(NPQ),W(NPQ),
     2 DB(ISPLH,ISPLH,3,NHEX,NPQ),DA(ISPLH,ISPLH,3,NHEX,NPQ),
     3 MN(NPQ,NSCT),DN(NSCT,NPQ),FUNKNO(NUN,NGEFF),WX(IELEM+1),
     3 WY(IELEM+1),CST(IELEM)
      LOGICAL LFIXUP,ISADPT(2)
*----
*  LOCAL VARIABLES
*----
      INTEGER :: NPQD(6),IIND(6),P,DCOORD
      REAL :: JAC(2,2,3), MUH, ETAH, AAA, BBB, CCC, DDD, MUHTEMP,
     1 ETAHTEMP, WX0(IELEM+1),WY0(IELEM+1)
      DOUBLE PRECISION :: Q(NM), Q2(NM,NM+1), V,THETA, XNI(NMX),
     1 XNJ(NMY)
      PARAMETER(IUNOUT=6,RLOG=1.0E-8,PI=3.141592654)
      LOGICAL :: LHEX(NHEX)
      LOGICAL ISFIX(2)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: INDANG
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: TMPMAT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: FLUX
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: FLUX_G
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: TMPXNI,
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
      ALLOCATE(FLUX(NM,NSCT))
      ALLOCATE(FLUX_G(NM,NSCT,3*ISPLH**2,NHEX,NGEFF))
      ALLOCATE(TMPXNI(IELEM,ISPLH,NCOL,NPQ,NGEFF))
      ALLOCATE(TMPXNJ(IELEM,ISPLH,NCOL,NPQ,NGEFF))
      ALLOCATE(TMPXND(IELEM,ISPLH,NCOL,NPQ,NGEFF))
*----
*  CONSTRUCT JACOBIAN MATRIX FOR EACH LOZENGE
*----
      JAC = RESHAPE((/ 1., -SQRT(3.), 1., SQRT(3.), 2., 0., 1.,
     >    SQRT(3.), 2., 0., -1., SQRT(3.) /), SHAPE(JAC))
      JAC = (SIDE/2.)*JAC
*----
*  LENGTH OF FUNKNO COMPONENTS (IN ORDER)
*----
      LFLX=3*NM*(ISPLH**2)*NHEX*NSCT
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

      FLUX_G(:NM,:NSCT,:3*ISPLH**2,:NHEX,:NGEFF)=0.0D0
      WX0=WX
      WY0=WY

      DO JND=1,6
      IND=IIND(JND)
      ! Needed because of S2 LS (4 dir. for 6 sextants)
      IF(IND.EQ.0) CYCLE
      TMPXNI(:IELEM,:ISPLH,:NCOL,:NPQ,:NGEFF)=0.0D0
      TMPXNJ(:IELEM,:ISPLH,:NCOL,:NPQ,:NGEFF)=0.0D0
      TMPXND(:IELEM,:ISPLH,:NCOL,:NPQ,:NGEFF)=0.0D0

*----
*  LOOP OVER WAVEFRONTS
*----
      DO IDI=1,NCOL+NCOL-1
*----
*  MAIN SWAPPING LOOPS FOR SN FLUX CALCULATION
*  LOOP OVER MACROCELLS IN WAVEFRONT AND ENERGY
*----

*$OMP  PARALLEL DO 
*$OMP+ PRIVATE(ITID,FLUX,M,IG,XNI,XNJ,Q,Q2,IOF,IER,II,JJ,IEL,I,J,P)
*$OMP+ PRIVATE(IPQD,IBM,SIGMA,V,ISFIX,IX,JX,IY,JY,AAA,BBB,CCC,DDD)
*$OMP+ PRIVATE(LHEX,IHEX,IIHEX,DCOORD,ILOZLOOP,ILOZ,IL,I2,JL,J2)
*$OMP+ PRIVATE(MUHTEMP,MUH,ETAHTEMP,ETAH,I3,I_FETCH)
*$OMP+ SHARED(TMPXNI,TMPXNJ,TMPXND,IND,IDI,FUNKNO)
*$OMP+ FIRSTPRIVATE(WX,WY,WX0,WY0)
*$OMP+ COLLAPSE(2)

      ! LOOP OVER MACROCELLS IN WAVEFRONT
      DO J_MC=MAX(1,IDI-NCOL+1),MIN(NCOL,IDI)

      ! LOOP FOR GROUPS TO EXECUTE IN PARALLEL
      DO IG=1,NGEFF

      ! LOOP OVER ALL DIRECTIONS
      DO IPQD=1,NPQD(IND)
      IF(.NOT.INCONV(IG)) CYCLE
      M=INDANG(IPQD,IND)
      IF(W(M).EQ.0.0) CYCLE

      ! GET AND PRINT THREAD NUMBER
#if defined(_OPENMP)
      ITID=omp_get_thread_num()
#else
      ITID=0
#endif
      IF(IMPX.GT.5) WRITE(IUNOUT,400) ITID,NGIND(IG),IPQD

      JIM=J_MC
      ! Account for different sweep direction depending on angle
      IF((IND.EQ.1).OR.(IND.EQ.2).OR.(IND.EQ.3)) JIM=NCOL+1-JIM

      ! Find I coordinate of macrocell based on Jth coordinate, the
      ! wavefront number and occasionally the number of rings in the
      ! domain
      IF((IND.EQ.1).OR.(IND.EQ.4)) THEN
        I_MC=IDI-J_MC+1
      ELSE
        I_MC=IDI-J_MC+(NRINGS+1-J_MC)
      ENDIF
      IIM=I_MC
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
      IF((I.GT.NCOL).OR.(I.LT.1)) CYCLE
      IF((J.GT.NCOL).OR.(J.LT.1)) CYCLE
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
      IF(IHEX.EQ.0) CALL XABORT('SNFDH2: IHEX FAILURE.')
      ! Find D coordinate
      DCOORD = ABS(COORDMAP(3,IHEX))-NRINGS

      ! LOOP OVER LOZENGES
      DO ILOZLOOP=1,3
      ILOZ=LOZSWP(ILOZLOOP,IND)

      ! Get Jacobian elements values
      AAA = JAC(1,1,ILOZ)
      BBB = JAC(1,2,ILOZ)
      CCC = JAC(2,1,ILOZ)
      DDD = JAC(2,2,ILOZ)


      ! LOOP OVER SUBMESH WITHIN EACH LOZENGE
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

      ! INITIALIZE FLUXES AND BOUNDARY FLUXES
      FLUX(:IELEM**2,:NSCT)=0.0D0

      ! READ IN XNI AND XNJ DEPENDING ON LOZENGE
      I_FETCH=0
      IF((ILOZ.EQ.1))THEN
        ! Read boundary fluxes in reverse for lozenge A since affine
        ! transformation of lozenges causes the D and I directions
        ! of lozenges C and A respectively to be reversed
        I_FETCH=ISPLH+1-I2
        XNI(:) = TMPXNI(:,J2,J,IPQD,IG)
        XNJ(:) = TMPXND(:,I_FETCH,DCOORD,IPQD,IG)
      ELSEIF((ILOZ.EQ.2))THEN
        XNI(:) = TMPXNI(:,J2,J,IPQD,IG)
        XNJ(:) = TMPXNJ(:,I2,I,IPQD,IG)
      ELSEIF((ILOZ.EQ.3))THEN
        XNI(:) = TMPXND(:,J2,DCOORD,IPQD,IG)
        XNJ(:) = TMPXNJ(:,I2,I,IPQD,IG)
      ENDIF

      ! DATA
      IBM=MAT(I2,J2,ILOZ,IHEX)
      ! Skip loop if virtual element 
      IF(IBM.EQ.0) CYCLE
      SIGMA=TOTAL(IBM,IG)
      V=VOL(I2,J2,ILOZ,IHEX)

      ! COMPUTE ADJUSTED DIRECTION COSINES
      MUHTEMP  =  DA(I2,J2,ILOZ,IHEX,M)
      ETAHTEMP =  DB(I2,J2,ILOZ,IHEX,M)
      MUH = (MUHTEMP*DDD) - (ETAHTEMP*BBB)
      ETAH = (-MUHTEMP*CCC) + (ETAHTEMP*AAA)

      ! SOURCE DENSITY TERM
      DO IEL=1,NM
      Q(IEL)=0.0D0
      DO P=1,NSCT
      IOF=((((((IHEX-1)*3+(ILOZ-1))*ISPLH+(J2-1))*ISPLH+
     1   (I2-1))*NSCT+(P-1))*NM)+IEL
      Q(IEL)=Q(IEL)+QEXT(IOF,IG)*MN(M,P)
      ENDDO
      ENDDO

      ISFIX=.FALSE.
      DO WHILE (.NOT.ALL(ISFIX)) ! LOOP FOR ADAPTIVE CALCULATION

      ! FLUX MOMENT COEFFICIENTS MATRIX
      Q2(:NM,:NM+1)=0.0D0
      
      DO IY=1,IELEM
      DO JY=1,IELEM
      DO IX=1,IELEM
      DO JX=1,IELEM
        II=IELEM*(IY-1)+IX
        JJ=IELEM*(JY-1)+JX

        ! DIAGONAL TERMS
        IF(II.EQ.JJ) THEN
          Q2(II,JJ)=SIGMA*V
     1              +CST(IX)**2*WX(JX+1)*ABS(MUH)
     2              +CST(IY)**2*WY(JY+1)*ABS(ETAH)

        ! UPPER DIAGONAL TERMS
        ELSEIF(II.LT.JJ) THEN
          ! X-SPACE TERMS
          IF(IY.EQ.JY) THEN
          IF(MOD(IX+JX,2).EQ.1) THEN
            Q2(II,JJ)=CST(IX)*CST(JX)*WX(JX+1)*MUH
          ELSE
            Q2(II,JJ)=CST(IX)*CST(JX)*WX(JX+1)*ABS(MUH)
          ENDIF
          ! Y-SPACE TERMS
          ELSEIF(IX.EQ.JX) THEN
          IF(MOD(IY+JY,2).EQ.1) THEN
            Q2(II,JJ)=CST(IY)*CST(JY)*WY(JY+1)*ETAH
          ELSE
            Q2(II,JJ)=CST(IY)*CST(JY)*WY(JY+1)*ABS(ETAH)
          ENDIF
          ENDIF

        ! UNDER DIAGONAL TERMS
        ELSE
          ! X-SPACE TERMS
          IF(IY.EQ.JY) THEN
          IF(MOD(IX+JX,2).EQ.1) THEN
            Q2(II,JJ)=CST(IX)*CST(JX)*(WX(JX+1)-2)*MUH
          ELSE
            Q2(II,JJ)=CST(IX)*CST(JX)*WX(JX+1)*ABS(MUH)
          ENDIF 
          ! Y-SPACE TERMS
          ELSEIF(IX.EQ.JX) THEN
          IF(MOD(IY+JY,2).EQ.1) THEN
            Q2(II,JJ)=CST(IY)*CST(JY)*(WY(JY+1)-2)*ETAH
          ELSE
            Q2(II,JJ)=CST(IY)*CST(JY)*WY(JY+1)*ABS(ETAH)
          ENDIF
          ENDIF

        ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      ! FLUX SOURCE VECTOR
      DO IY=1,IELEM
      DO IX=1,IELEM
        II=IELEM*(IY-1)+IX
        Q2(II,NM+1)=Q(II)*V
        ! X-SPACE TERMS
        IF(MOD(IX,2).EQ.1) THEN
          Q2(II,NM+1)=Q2(II,NM+1)+CST(IX)*(1-WX(1))
     1                *XNI(IY)*ABS(MUH)
        ELSE
          Q2(II,NM+1)=Q2(II,NM+1)-CST(IX)*(1+WX(1))
     1                *XNI(IY)*MUH
        ENDIF
        ! Y-SPACE TERMS
        IF(MOD(IY,2).EQ.1) THEN
          Q2(II,NM+1)=Q2(II,NM+1)+CST(IY)*(1-WY(1))
     1                *XNJ(IX)*ABS(ETAH)
        ELSE
          Q2(II,NM+1)=Q2(II,NM+1)-CST(IY)*(1+WY(1))
     1                *XNJ(IX)*ETAH
        ENDIF
      ENDDO
      ENDDO

      CALL ALSBD(NM,1,Q2,IER,NM)
      IF(IER.NE.0) CALL XABORT('SNFKH2: SINGULAR MATRIX.')

      ! ADAPTIVE CORRECTION OF WEIGHTING PARAMETERS
      IF(ANY(ISADPT)) THEN
        IF(ISADPT(1)) THEN
          CALL SNADPT(IELEM,NM,IELEM,Q2(1:IELEM:1,NM+1),
     1    XNI(:NMX),1.0,WX,ISFIX(1))
        ELSE
          ISFIX(1)=.TRUE.
        ENDIF
        IF(ISADPT(2)) THEN
          CALL SNADPT(IELEM,NM,IELEM,Q2(1:NM:IELEM,NM+1),
     1    XNJ(:NMY),1.0,WY,ISFIX(2))
        ELSE
          ISFIX(2)=.TRUE.
        ENDIF
      ELSE
        ISFIX=.TRUE.
      ENDIF

      END DO ! END OF ADAPTIVE LOOP

      ! CLOSURE RELATIONS
      IF(IELEM.EQ.1.AND.LFIXUP.AND.(Q2(1,2).LE.RLOG)) Q2(1,2)=0.0
      ! Read XNI/XNI into TMPXNI/J/D
      IF((ILOZ.EQ.1).OR.(ILOZ.EQ.2))THEN
        TMPXNI(:NMX,J2,J,IPQD,IG)=WX(1)*XNI(:NMX)
      ELSE
        TMPXND(:NMX,J2,DCOORD,IPQD,IG)=WX(1)*XNI(:NMX)
      ENDIF
      IF((ILOZ.EQ.2).OR.(ILOZ.EQ.3))THEN
        TMPXNJ(:NMY,I2,I,IPQD,IG)=WY(1)*XNJ(:NMY)
      ELSE
        I3=I_FETCH
        TMPXND(:NMY,I3,DCOORD,IPQD,IG)=WY(1)*XNJ(:NMY)
      ENDIF
      DO IY=1,IELEM
      DO IX=1,IELEM
        II=IELEM*(IY-1)+IX
        ! X-SPACE
        ! Assign I-boundary fluxes if lozenges A or B
        IF((ILOZ.EQ.1).OR.(ILOZ.EQ.2))THEN
        IF(MOD(IX,2).EQ.1) THEN
          TMPXNI(IY,J2,J,IPQD,IG)=TMPXNI(IY,J2,J,IPQD,IG)+CST(IX)*
     1                 WX(IX+1)*Q2(II,NM+1)
        ELSE
          TMPXNI(IY,J2,J,IPQD,IG)=TMPXNI(IY,J2,J,IPQD,IG)+CST(IX)*
     1                 WX(IX+1)*Q2(II,NM+1)*SIGN(1.0,MUH)
        ENDIF
        ENDIF
        ! Y-SPACE
        ! Assign J-boundary fluxes if lozenges B or C
        IF((ILOZ.EQ.2).OR.(ILOZ.EQ.3))THEN
        IF(MOD(IY,2).EQ.1) THEN
          TMPXNJ(IX,I2,I,IPQD,IG)=TMPXNJ(IX,I2,I,IPQD,IG)+CST(IY)*
     1                 WY(IY+1)*Q2(II,NM+1)
        ELSE
          TMPXNJ(IX,I2,I,IPQD,IG)=TMPXNJ(IX,I2,I,IPQD,IG)+CST(IY)*
     1                 WY(IY+1)*Q2(II,NM+1)*SIGN(1.0,ETAH)
        ENDIF
        ENDIF
        ! D-SPACE
        ! Assign D-boundary fluxes if lozenge A using XNJ
        IF((ILOZ.EQ.1))THEN
        I3=I_FETCH
        IF(MOD(IY,2).EQ.1) THEN
          TMPXND(IX,I3,DCOORD,IPQD,IG)=TMPXND(IX,I3,DCOORD,IPQD,IG)+
     1                 CST(IY)*WY(IY+1)*Q2(II,NM+1)
        ELSE
          TMPXND(IX,I3,DCOORD,IPQD,IG)=TMPXND(IX,I3,DCOORD,IPQD,IG)+
     1                 CST(IY)*WY(IY+1)*Q2(II,NM+1)*SIGN(1.0,ETAH)
        ENDIF
        ENDIF
        ! Assign D-boundary fluxes if lozenge C using XNI
        IF((ILOZ.EQ.3))THEN
        IF(MOD(IX,2).EQ.1) THEN
          TMPXND(IY,J2,DCOORD,IPQD,IG)=TMPXND(IY,J2,DCOORD,IPQD,IG)+
     1                 CST(IX)*WX(IX+1)*Q2(II,NM+1)
        ELSE
          TMPXND(IY,J2,DCOORD,IPQD,IG)=TMPXND(IY,J2,DCOORD,IPQD,IG)+
     1                 CST(IX)*WX(IX+1)*Q2(II,NM+1)*SIGN(1.0,MUH)
        ENDIF
        ENDIF
      ENDDO
      ENDDO
      ! FLIP GRADIENTS IF NECESSARY
      DO IY=1,IELEM
        IF((MOD(IY,2).EQ.0).AND.(ILOZ.EQ.3).AND.(IL.EQ.ISPLH))
     1    TMPXND(IY,J2,DCOORD,IPQD,IG)=TMPXND(IY,J2,DCOORD,IPQD,IG)*(-1)
      ENDDO
      I3=I_FETCH
      DO IX=1,IELEM
        IF((MOD(IX,2).EQ.0).AND.(ILOZ.EQ.1).AND.(JL.EQ.ISPLH))
     1    TMPXND(IX,I3,DCOORD,IPQD,IG)=TMPXND(IX,I3,DCOORD,IPQD,IG)*(-1)
      ENDDO
      ! LFIXUP
      IF(IELEM.EQ.1.AND.LFIXUP)THEN
        IF((ILOZ.EQ.1).OR.(ILOZ.EQ.2))THEN
          IF(TMPXNI(1,J2,J,IPQD,IG).LE.RLOG) TMPXNI(1,J2,J,IPQD,IG)=0.0
        ELSE
          IF(TMPXND(1,J2,DCOORD,IPQD,IG).LE.RLOG) 
     1      TMPXND(1,J2,DCOORD,IPQD,IG)=0.0
        ENDIF
        IF((ILOZ.EQ.2).OR.(ILOZ.EQ.3))THEN
          IF(TMPXNJ(1,I2,I,IPQD,IG).LE.RLOG) TMPXNJ(1,I2,I,IPQD,IG)=0.0
        ELSE
          I3=I_FETCH
          IF(TMPXND(1,I3,DCOORD,IPQD,IG).LE.RLOG) 
     1      TMPXND(1,I3,DCOORD,IPQD,IG)=0.0
        ENDIF
      ENDIF
      WX=WX0
      WY=WY0
      ! SAVE LEGENDRE MOMENT OF THE FLUX
      DO P=1,NSCT
      DO IEL=1,NM
      FLUX(IEL,P)=FLUX(IEL,P) + Q2(IEL,NM+1)*DN(P,M)
      ENDDO
      ENDDO

      ! SAVE FLUX INFORMATION
      IOF=((ILOZ-1)*ISPLH+(J2-1))*ISPLH+I2
      FLUX_G(:,:,IOF,IHEX,IG)=FLUX_G(:,:,IOF,IHEX,IG)+FLUX(:,:)

      ENDDO ! END OF WITHIN LOZENGE J-LOOP
      ENDDO ! END OF WITHIN LOZENGE I-LOOP

      ENDDO ! END OF LOZENGE LOOP

      ENDDO ! END OF DIRECTION LOOP
      ENDDO ! END OF ENERGY LOOP
  
      ENDDO ! END OF MACROCELL LOOP
*$OMP END PARALLEL DO
      ENDDO ! END OF WAVEFRONT LOOP

      ENDDO ! END OF OCTANT LOOP

      ! SAVE FLUX INFORMATION
      DO IG=1,NGEFF
        IF(.NOT.INCONV(IG)) CYCLE
        FUNKNO(:LFLX,IG)=
     1  RESHAPE(REAL(FLUX_G(:IELEM**2,:NSCT,:3*ISPLH**2,:NHEX,IG)),
     2   (/ LFLX /) )
      ENDDO

      ! CALL XABORT('SNFKH2: testing')

*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(FLUX_G,FLUX,INDANG,TMPXNI,TMPXNJ,TMPXND,TMPMAT)
      RETURN
  400 FORMAT(16H SNFKH2: thread=,I8,12H --->(group=,I4,7H angle=,I4,1H))
      END
