*DECK SNFDH1
      SUBROUTINE SNFDH1(NUN,NGEFF,IMPX,INCONV,NGIND,NHEX,IELEM,ISPLH,
     1 NMAT,NPQ,NSCT,MAT,VOL,TOTAL,QEXT,DU,DE,W,DB,DA,MN,DN,SIDE,
     2 LOZSWP,COORDMAP,FUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform one inner iteration for solving SN equations in 2D hexagonal
* geometry for the HODD method. KBA-like multithreading, i.e.,
* macrocell-energy. VOID boundary conditions.
*
*Copyright:
* Copyright (C) 2025 Ecole Polytechnique de Montreal
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
      DOUBLE PRECISION :: Q(IELEM**2), Q2(IELEM**2,(IELEM**2)+1), VT,
     > CONST0, CONST1, CONST2, THETA, XNI(IELEM), XNJ(IELEM), C1
      REAL :: JAC(2,2,3), MUH, ETAH, AAA, BBB, CCC, DDD, MUHTEMP,
     > ETAHTEMP
      LOGICAL :: LHEX(NHEX)
*----
*  LOCAL ALLOCATABLE ARRAYS
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
      ALLOCATE(FLUX(IELEM**2,NSCT))
      ALLOCATE(FLUX_G(IELEM**2,NSCT,3*ISPLH**2,NHEX,NGEFF))
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
*  DEFINITION OF CONSTANTS
*----
      L4=3*(IELEM**2)*(ISPLH**2)*NHEX*NSCT
      CONST0=2.0D0*DSQRT(3.0D0)
      CONST1=2.0D0*DSQRT(5.0D0)
      CONST2=2.0D0*DSQRT(15.0D0)
*----
*  PARAMETER VALIDATION
*----
      IF((IELEM.LT.1).OR.(IELEM.GT.3)) THEN
         CALL XABORT('SNFDH1: INVALID IELEM (DIAM) VALUE. '
     >   //'CHECK INPUT DATA FILE.')
      ENDIF

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
*$OMP1 PRIVATE(ITID,I,J,LHEX,IHEX,IIHEX,DCOORD,IOF)
*$OMP1 PRIVATE(ILOZLOOP,ILOZ,AAA,BBB,CCC,DDD,IL,I2,JL,J2,FLUX)
*$OMP1 PRIVATE(XNI,XNJ,MUHTEMP,MUH,ETAHTEMP,ETAH,IEL,Q,P)
*$OMP1 PRIVATE(VT,Q2,IER,C1,I3,I_FETCH,IPQD,M)
*$OMP1 SHARED(IND,TMPXNI,TMPXNJ,TMPXND,IDI)
*$OMP1 COLLAPSE(2)
      DO J_MC=MAX(1,IDI-NCOL+1),MIN(NCOL,IDI)
      DO IG=1,NGEFF
      DO IPQD=1,NPQD(IND)

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
      IF(IHEX.EQ.0) CALL XABORT('SNFDH1: IHEX FAILURE.')
      ! Find D coordinate
      DCOORD = ABS(COORDMAP(3,IHEX))-NRINGS

#if defined(_OPENMP)
        ITID=omp_get_thread_num()
#else
        ITID=0
#endif
      IF(IMPX.GT.5) WRITE(IUNOUT,400) ITID,NGIND(IG),IPQD
      M=INDANG(IPQD,IND)
      IF(.NOT.INCONV(IG)) CYCLE
      IF(W(M).EQ.0.0) CYCLE
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

      FLUX(:IELEM**2,:NSCT)=0.0D0

      ! Read in XNI and XNJ depending on lozenge
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

      ! Initialise and compute matrix of coefficients, and source
      Q2(:IELEM**2,:IELEM**2+1)=0.0D0
      IF(IELEM.EQ.1) THEN
        Q2(1,1)= 2.0D0*ABS(MUH) + 2.0D0*ABS(ETAH) + VT
        !------
        Q2(1,2)= 2.0D0*ABS(MUH)*XNI(1) + 2.0D0*ABS(ETAH)*XNJ(1)
     >    + VOL(I2,J2,J,I)*Q(1)
      ELSE IF(IELEM.EQ.2) THEN
        Q2(1,1)=VT
        Q2(2,1)=CONST0*MUH
        Q2(2,2)=-VT-6.0D0*ABS(MUH)
        Q2(3,1)=CONST0*ETAH
        Q2(3,3)=-VT-6.0D0*ABS(ETAH)
        Q2(4,2)=-CONST0*ETAH
        Q2(4,3)=-CONST0*MUH
        Q2(4,4)=VT+6.0D0*ABS(MUH)+6.0D0*ABS(ETAH)
        !------
        Q2(1,5)=VOL(I2,J2,J,I)*Q(1)
        Q2(2,5)=-VOL(I2,J2,J,I)*Q(2)+CONST0*MUH*XNI(1)
        Q2(3,5)=-VOL(I2,J2,J,I)*Q(3)+CONST0*ETAH*XNJ(1)
        Q2(4,5)=VOL(I2,J2,J,I)*Q(4)-CONST0*MUH*XNI(2)-
     >    CONST0*ETAH*XNJ(2)
      ELSE IF(IELEM.EQ.3) THEN
        Q2(1,1)=VT+2.0D0*ABS(MUH)+2.0D0*ABS(ETAH)
        Q2(2,2)=-VT-2.0D0*ABS(ETAH)
        Q2(3,1)=CONST1*ABS(MUH)
        Q2(3,2)=-CONST2*MUH
        Q2(3,3)=VT+1.0D1*ABS(MUH)+2.0D0*ABS(ETAH)
        Q2(4,4)=-VT-2.0D0*ABS(MUH)
        Q2(5,5)=VT
        Q2(6,4)=-CONST1*ABS(MUH)
        Q2(6,5)=CONST2*MUH
        Q2(6,6)=-VT-1.0D1*ABS(MUH)
        Q2(7,1)=CONST1*ABS(ETAH)
        Q2(7,4)=-CONST2*ETAH
        Q2(7,7)=VT+2.0D0*ABS(MUH)+1.0D1*ABS(ETAH)
        Q2(8,2)=-CONST1*ABS(ETAH)
        Q2(8,5)=CONST2*ETAH
        Q2(8,8)=-VT-1.0D1*ABS(ETAH)
        Q2(9,3)=CONST1*ABS(ETAH)
        Q2(9,6)=-CONST2*ETAH
        Q2(9,7)=CONST1*ABS(MUH)
        Q2(9,8)=-CONST2*MUH
        Q2(9,9)=VT+1.0D1*ABS(MUH)+1.0D1*ABS(ETAH)
        !------
        Q2(1,10)=VOL(I2,J2,J,I)*Q(1)+2.0D0*ABS(MUH)*XNI(1)+
     >    2.0D0*ABS(ETAH)*XNJ(1)
        Q2(2,10)=-VOL(I2,J2,J,I)*Q(2)-2.0D0*ABS(ETAH)*XNJ(2)
        Q2(3,10)=VOL(I2,J2,J,I)*Q(3)+CONST1*ABS(MUH)*XNI(1)+
     >    2.0D0*ABS(ETAH)*XNJ(3)
        Q2(4,10)=-VOL(I2,J2,J,I)*Q(4)-2.0D0*ABS(MUH)*XNI(2)
        Q2(5,10)=VOL(I2,J2,J,I)*Q(5)
        Q2(6,10)=-VOL(I2,J2,J,I)*Q(6)-CONST1*ABS(MUH)*XNI(2)
        Q2(7,10)=VOL(I2,J2,J,I)*Q(7)+2.0D0*ABS(MUH)*XNI(3)+
     >    CONST1*ABS(ETAH)*XNJ(1)
        Q2(8,10)=-VOL(I2,J2,J,I)*Q(8)-CONST1*ABS(ETAH)*XNJ(2)
        Q2(9,10)=VOL(I2,J2,J,I)*Q(9)+CONST1*ABS(MUH)*XNI(3)+
     >    CONST1*ABS(ETAH)*XNJ(3)
      ENDIF
      DO IEL=1,IELEM**2
      DO JEL=IEL+1,IELEM**2
      Q2(IEL,JEL)=Q2(JEL,IEL)
      ENDDO
      ENDDO

      ! Solve system of equation given by matrix of coefficients
      ! and source
      CALL ALSBD(IELEM**2,1,Q2,IER,IELEM**2)
      IF(IER.NE.0) CALL XABORT('SNFDH1: SINGULAR MATRIX.')

      ! Assign J-boundary fluxes if lozenges B or C
      IF((ILOZ.EQ.2).OR.(ILOZ.EQ.3))THEN
        IF(IELEM.EQ.1) THEN
          TMPXNJ(1,I2,I,IPQD,IG) = 2.0D0*Q2(1,2)-XNJ(1)
        ELSEIF(IELEM.EQ.2) THEN
          TMPXNJ(1,I2,I,IPQD,IG)=XNJ(1)+SIGN(1.0,ETAH)*CONST0*Q2(3,5)
          TMPXNJ(2,I2,I,IPQD,IG)=XNJ(2)+SIGN(1.0,ETAH)*CONST0*Q2(4,5)
        ELSEIF(IELEM.EQ.3) THEN
          TMPXNJ(1,I2,I,IPQD,IG)=2.0D0*Q2(1,10)+CONST1*Q2(7,10)-XNJ(1)
          TMPXNJ(2,I2,I,IPQD,IG)=2.0D0*Q2(2,10)+CONST1*Q2(8,10)-XNJ(2)
          TMPXNJ(3,I2,I,IPQD,IG)=2.0D0*Q2(3,10)+CONST1*Q2(9,10)-XNJ(3)
        ENDIF
      ENDIF
      ! Assign I-boundary fluxes if lozenges A or B
      IF((ILOZ.EQ.1).OR.(ILOZ.EQ.2))THEN
        IF(IELEM.EQ.1) THEN
          TMPXNI(1,J2,J,IPQD,IG) = 2.0D0*Q2(1,2)-XNI(1)
        ELSEIF(IELEM.EQ.2) THEN
          TMPXNI(1,J2,J,IPQD,IG)=XNI(1)+SIGN(1.0,MUH)*CONST0*Q2(2,5)
          TMPXNI(2,J2,J,IPQD,IG)=XNI(2)+SIGN(1.0,MUH)*CONST0*Q2(4,5)
        ELSEIF(IELEM.EQ.3) THEN
          TMPXNI(1,J2,J,IPQD,IG)=2.0D0*Q2(1,10)+CONST1*Q2(3,10)-XNI(1)
          TMPXNI(2,J2,J,IPQD,IG)=2.0D0*Q2(4,10)+CONST1*Q2(6,10)-XNI(2)
          TMPXNI(3,J2,J,IPQD,IG)=2.0D0*Q2(7,10)+CONST1*Q2(9,10)-XNI(3)
        ENDIF
      ENDIF
      ! Assign D-boundary fluxes if lozenge A using XNJ
      IF((ILOZ.EQ.1))THEN
        C1=1.0D0
        I3=I_FETCH
        IF(JL.EQ.ISPLH) C1=-1.0D0
        IF(IELEM.EQ.1) THEN
          TMPXND(1,I3,DCOORD,IPQD,IG) = 2.0D0*Q2(1,2)-XNJ(1)
        ELSEIF(IELEM.EQ.2) THEN
          TMPXND(1,I3,DCOORD,IPQD,IG)=XNJ(1)+
     >      SIGN(1.0,ETAH)*CONST0*Q2(3,5)
          TMPXND(2,I3,DCOORD,IPQD,IG)=(XNJ(2)+
     >      SIGN(1.0,ETAH)*CONST0*Q2(4,5))*C1
        ELSEIF(IELEM.EQ.3) THEN
          TMPXND(1,I3,DCOORD,IPQD,IG)=2.0D0*Q2(1,10)+CONST1*Q2(7,10)-
     >      XNJ(1)
          TMPXND(2,I3,DCOORD,IPQD,IG)=(2.0D0*Q2(2,10)+CONST1*Q2(8,10)-
     >      XNJ(2))*C1
          TMPXND(3,I3,DCOORD,IPQD,IG)=2.0D0*Q2(3,10)+CONST1*Q2(9,10)-
     >      XNJ(3)
        ENDIF
      ENDIF
      ! Assign D-boundary fluxes if lozenge C using XNI
      IF((ILOZ.EQ.3))THEN
        C1=1.0D0
        IF(IL.EQ.ISPLH) C1=-1.0D0
        IF(IELEM.EQ.1) THEN
          TMPXND(1,J2,DCOORD,IPQD,IG) = 2.0D0*Q2(1,2)-XNI(1)
        ELSEIF(IELEM.EQ.2) THEN
          TMPXND(1,J2,DCOORD,IPQD,IG)=XNI(1)+
     >      SIGN(1.0,MUH)*CONST0*Q2(2,5)
          TMPXND(2,J2,DCOORD,IPQD,IG)=(XNI(2)+
     >      SIGN(1.0,MUH)*CONST0*Q2(4,5))*C1
        ELSEIF(IELEM.EQ.3) THEN
          TMPXND(1,J2,DCOORD,IPQD,IG)=2.0D0*Q2(1,10)+CONST1*Q2(3,10)-
     >      XNI(1)
          TMPXND(2,J2,DCOORD,IPQD,IG)=(2.0D0*Q2(4,10)+CONST1*Q2(6,10)-
     >      XNI(2))*C1
          TMPXND(3,J2,DCOORD,IPQD,IG)=2.0D0*Q2(7,10)+CONST1*Q2(9,10)-
     >      XNI(3)
        ENDIF
      ENDIF

      ! Compute flux
      DO P=1,NSCT
      DO IEL=1,IELEM**2
      FLUX(IEL,P) = FLUX(IEL,P) + Q2(IEL,IELEM**2+1)*DN(P,M)
      ENDDO
      ENDDO

      IOF=((ILOZ-1)*ISPLH+(J2-1))*ISPLH+I2
      FLUX_G(:,:,IOF,IHEX,IG)=FLUX_G(:,:,IOF,IHEX,IG)+FLUX(:,:)

      ENDDO ! jl, isplh
      ENDDO ! il, isplh

      ENDDO ! iloz

      ENDDO ! ipqd
      ENDDO ! ig
      ENDDO ! j_mc
*$OMP END PARALLEL DO

      ENDDO ! idi

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
  400 FORMAT(16H SNFDH1: thread=,I8,12H --->(group=,I4,7H angle=,I4,1H))
      END
