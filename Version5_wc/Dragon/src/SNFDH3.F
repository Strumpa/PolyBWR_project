*DECK SNFDH3
      SUBROUTINE SNFDH3(NKBA,NUN,NGEFF,IMPX,INCONV,NGIND,NHEX,LZ,IELEM,
     1 ISPLH,SIDE,NMAT,NPQ,NSCT,MAT,VOL,NCODE,ZCODE,TOTAL,QEXT,DU,DE,
     2 DZ,W,MRMZ,DC,DB,DA,MN,DN,LOZSWP,COORDMAP,FUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform one inner iteration for solving SN equations in 3D hexagonal
* geometry for the HODD method. VOID boundary conditions on sides, and,
* albedo boundary conditions on top and bottom. KBA-like multithreading,
* i.e., macrocell-energy. 
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
* LZ       number of meshes along Z axis.
* IELEM    measure of order of the spatial approximation polynomial:
*          =1 constant (classical diamond scheme - default for HODD);
*          =2 linear (default for DG);
*          =3 parabolic;
*          =4 cubic (only for DG).
* SIDE     side of an hexagon.
* NMAT     number of material mixtures.
* ISPLH    splitting option for hexagons.
* NPQ      number of SN directions in 12 dodecants (including zero-weight
*          directions).
* NSCT     maximum number of spherical harmonics moments of the flux.
* MAT      material mixture index in each region.
* VOL      volumes of each region.
* TOTAL    macroscopic total cross sections.
* NCODE    boundary condition indices.
* ZCODE    albedos.
* QEXT     Legendre components of the fixed source.
* DU       first direction cosines ($\\mu$).
* DE       second direction cosines ($\\eta$).
* DZ       third direction cosines ($\\xi$).
* W        weights.
* MRMZ     quadrature index.
* DC       diamond-scheme parameter.
* DB       diamond-scheme parameter.
* DA       diamond-scheme parameter.
* MN       moment-to-discrete matrix.
* DN       discrete-to-moment matrix.
* LOZSWP   lozenge sweep order depending on direction.
* COORDMAP coordinate map - mapping the hexagons from the indices 
*          within the DRAGON geometry to a Cartesian axial coordinate
*          array (see redblobgames.com website).
*
*Parameters: input/output
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
      INTEGER NUN,NGEFF,IMPX,NGIND(NGEFF),NHEX,LZ,IELEM,ISPLH,NMAT,NPQ,
     > NSCT,MAT(ISPLH,ISPLH,3,NHEX,LZ),NCODE(6),MRMZ(NPQ),LOZSWP(3,6),
     > COORDMAP(3,NHEX)
      LOGICAL INCONV(NGEFF)
      REAL SIDE,VOL(ISPLH,ISPLH,3,NHEX,LZ),ZCODE(6),TOTAL(0:NMAT,NGEFF),
     > QEXT(NUN,NGEFF),DU(NPQ),DE(NPQ),DZ(NPQ),W(NPQ),
     > DC(ISPLH*ISPLH*3*NHEX,1,NPQ),DB(ISPLH*ISPLH*3*NHEX,LZ,NPQ),
     > DA(1,LZ,NPQ),MN(NPQ,NSCT),DN(NSCT,NPQ),FUNKNO(NUN,NGEFF)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(IUNOUT=6,PI=3.141592654)
      INTEGER :: NPQD(12),IIND(12),P,DCOORD
      DOUBLE PRECISION Q(IELEM**3),Q2(IELEM**3,(IELEM**3)+1),VT,CONST0,
     > CONST1,CONST2,THETA,XNI(IELEM,IELEM),XNJ(IELEM,IELEM),
     > XNK(IELEM,IELEM),C1
      REAL :: JAC(2,2,3),MUH,ETAH,XIH,AAA,BBB,CCC,DDD,MUHTEMP,ETAHTEMP
      LOGICAL :: LHEX(NHEX)
*----
*  LOCAL ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: III,JJJ,KKK
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: INDANG
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: TMPMAT
      DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:,:) :: FLUX
      DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: FLUX_G
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:) :: TMPXNI,
     > TMPXNJ, TMPXND
      TYPE MATRIX_ARRAY
        DOUBLE PRECISION, POINTER, DIMENSION(:,:,:,:,:,:) :: MATRIX
      END TYPE MATRIX_ARRAY
      TYPE(MATRIX_ARRAY), POINTER, DIMENSION(:,:) :: TMPXNK
      TYPE(MATRIX_ARRAY), POINTER :: PT
*----
*  MAP MATERIAL VALUES TO CARTESIAN AXIAL COORDINATE MAP
*----
      NRINGS=INT((SQRT(  REAL((4*NHEX-1)/3)  )+1.)/2.)
      NCOL=2*NRINGS -1
      ALLOCATE(TMPMAT(ISPLH,ISPLH,3,NCOL,NCOL,LZ))
      TMPMAT(:,:,:,:,:,:) = -1
      DO IZ=1,LZ
        DO IHEX_XY=1,NHEX
          TMPMAT(:,:,:,COORDMAP(1,IHEX_XY),COORDMAP(2,IHEX_XY),IZ) = 
     >         MAT(:,:,:,IHEX_XY,IZ)
        ENDDO
      ENDDO
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(INDANG(NPQ,12))
      ALLOCATE(FLUX(IELEM**3,NSCT))
      ALLOCATE(FLUX_G(IELEM**3,NSCT,3*ISPLH**2,NHEX,LZ,NGEFF))
      ALLOCATE(TMPXNI(IELEM,IELEM,ISPLH,NCOL,LZ,NPQ,NGEFF))
      ALLOCATE(TMPXNJ(IELEM,IELEM,ISPLH,NCOL,LZ,NPQ,NGEFF))
      ALLOCATE(TMPXND(IELEM,IELEM,ISPLH,NCOL,LZ,NPQ,NGEFF))
      ALLOCATE(TMPXNK(NPQ,NGEFF))
*----
*  CONSTRUCT JACOBIAN MATRIX FOR EACH LOZENGE
*----
      JAC = RESHAPE((/ 1., -SQRT(3.), 1., SQRT(3.), 2., 0., 1.,
     >    SQRT(3.), 2., 0., -1., SQRT(3.) /), SHAPE(JAC))
      JAC = (SIDE/2.)*JAC
*----
*  DEFINITION OF CONSTANTS
*----
      NMZ=NKBA
      NCELLZ = (LZ + NKBA - 1) / NKBA
      L4=(IELEM**3)*(ISPLH**2)*3*NHEX*LZ*NSCT
      L5=(IELEM**2)*(ISPLH**2)*3*NHEX
      CONST0=2.0D0*DSQRT(3.0D0)
      CONST1=2.0D0*DSQRT(5.0D0)
      CONST2=2.0D0*DSQRT(15.0D0)
*----
*  PARAMETER VALIDATION
*----
      IF((IELEM.LT.1).OR.(IELEM.GT.3))
     >   CALL XABORT('SNFDH3: INVALID IELEM (DIAM) VALUE. '
     >   //'CHECK INPUT DATA FILE.')

      FLUX_G(:IELEM**3,:NSCT,:3*ISPLH**2,:NHEX,:LZ,:NGEFF)=0.0D0
*----
*  SET DODECANT SWAPPING ORDER
*----
      NPQD(:12)=0
      INDANG(:NPQ,:12)=0
      IIND(:12)=0
      DO M=1,NPQ
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
        ! UNFOLD DODECANTS
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
        ! Assume IIND(I)=I in hexagonal geometry
        IIND(IND)=IND
        NPQD(IND)=NPQD(IND)+1
        INDANG(NPQD(IND),IND)=M
      ENDDO
*----
*  MAIN LOOP OVER DODECANTS
*----
      DO JND=1,12
      IND=IIND(JND)
      IND_XY=MOD(IND-1,6)+1
      ! Needed because of S2 LS (8 dir. for 12 dodecants)
      IF(IND.EQ.0) CYCLE
*----
*  PRELIMINARY LOOPS FOR SETTING BOUNDARY CONDITIONS
*----
      IF((NCODE(5).NE.1).or.(NCODE(6).NE.1))THEN
      DO IG=1,NGEFF
      DO IPQD=1,NPQD(IND)
      IF(.NOT.INCONV(IG)) CYCLE
      M=INDANG(IPQD,IND)
      VZ=DZ(M)
*
      IF(VZ.GT.0.0)THEN
        M1=MRMZ(M)
        IF(NCODE(5).NE.4)THEN
          IOF=(M-1)*(L5)
          JOF=(M1-1)*(L5)
          FUNKNO(L4+IOF+1:L4+IOF+L5,IG)=
     >      FUNKNO(L4+JOF+1:L4+JOF+L5,IG)
        ENDIF
      ELSEIF(VZ.LT.0.0)THEN
        M1=MRMZ(M)
        IF(NCODE(6).NE.4)THEN
          IOF=(M-1)*(L5)
          JOF=(M1-1)*(L5)
          FUNKNO(L4+IOF+1:L4+IOF+L5,IG)=
     >      FUNKNO(L4+JOF+1:L4+JOF+L5,IG)
        ENDIF
      ENDIF
*
      ENDDO
      ENDDO
      ENDIF

      TMPXNI(:IELEM,:IELEM,:ISPLH,:NCOL,:LZ,:NPQ,:NGEFF)=0.0D0
      TMPXNJ(:IELEM,:IELEM,:ISPLH,:NCOL,:LZ,:NPQ,:NGEFF)=0.0D0
      TMPXND(:IELEM,:IELEM,:ISPLH,:NCOL,:LZ,:NPQ,:NGEFF)=0.0D0
      DO IG=1,NGEFF
        DO IPQD=1,NPQ
          PT => TMPXNK(IPQD,IG)
          ALLOCATE(PT%MATRIX(IELEM,IELEM,ISPLH,ISPLH,3,NHEX))
          PT%MATRIX(:IELEM,:IELEM,:ISPLH,:ISPLH,:3,:NHEX)=0.0D0
        ENDDO
      ENDDO
*----
*  LOOP OVER WAVEFRONTS
*----
      DO IDI=1,NCOL+NCOL+NMZ-2
*----
*  SET SWEEP INDICES
*----
      NMAX=MIN(NCOL,IDI)*MIN(NCOL,IDI)
      ALLOCATE(III(NMAX),JJJ(NMAX),KKK(NMAX))
      NCEL=0

      J_STT=MAX(1,IDI-NCOL-NMZ+2)
      J_END=MIN(NCOL,IDI)

      DO J_MC=J_STT,J_END
      JIM=J_MC
      ! Account for different sweep direction depending on angle
      IF((IND_XY.EQ.1).OR.(IND_XY.EQ.2).OR.(IND_XY.EQ.3)) JIM=NCOL+1-JIM

      IF((IND_XY.EQ.1).OR.(IND_XY.EQ.4)) THEN
        I_STT=MAX(1,IDI-J_MC-NMZ+2)
        I_END=MIN(NCOL,IDI-J_MC+1)
      ELSE
        I_STT=MAX(1,NRINGS-J_MC+1)
        I_END=MIN(NCOL,IDI-J_MC+(NRINGS+1-J_MC))
      ENDIF

      DO I_MC=I_STT,I_END
      IIM=I_MC
      ! Account for different sweep direction depending on angle
      IF((IND_XY.EQ.2).OR.(IND_XY.EQ.3).OR.(IND_XY.EQ.4)) IIM=NCOL+1-IIM

      ! For IND_XY 3 or 6, Cartesian axial coordinate map is swept 
      ! vertically instead of horizontally. IM suffix is for 'IMmutable'
      I=IIM
      J=JIM
      IF((IND_XY.EQ.3).OR.(IND_XY.EQ.6))THEN
        I=JIM
        J=IIM
      ENDIF

      ! If within corners of Cartesian axial coordinate map (where
      ! there are no hexagons), skip loop
      IF(TMPMAT(1,1,1,I,J,1).EQ.-1) CYCLE

      ! Find I coordinate of macrocell
      IF((IND_XY.EQ.1).OR.(IND_XY.EQ.4)) THEN
        K_MC=IDI-I_MC-J_MC+2
      ELSE
        K_MC=IDI-I_MC+NRINGS-((J_MC-1)*2)
      ENDIF
      IF((K_MC.GT.NMZ)) CYCLE
      K=K_MC

      NCEL=NCEL+1
      IF(NCEL.GT.NMAX) CALL XABORT('SNFDH3: NMAX OVERFLOW.')
      III(NCEL)=I
      JJJ(NCEL)=J
      KKK(NCEL)=K
      ENDDO ! I_MC

      ENDDO ! J_MC
*
*----
*  MAIN SWAPPING LOOPS FOR SN FLUX CALCULATION
*  LOOP OVER MACROCELLS IN WAVEFRONT AND ENERGY
*----
*$OMP  PARALLEL DO
*$OMP1 PRIVATE(ITID,I,J,K_MC,LHEX,IHEX_XY,IIHEX,PT,DCOORD,IOF)
*$OMP1 PRIVATE(IMZ,K,ILOZLOOP,ILOZ,AAA,BBB,CCC,DDD,IL,I2,JL,J2,FLUX)
*$OMP1 PRIVATE(XNI,XNJ,XNK,MUHTEMP,MUH,ETAHTEMP,ETAH,XIH,IEL,Q,P)
*$OMP1 PRIVATE(VT,Q2,IER,C1,I3,I_FETCH,IPQD,M)
*$OMP1 SHARED(IND,TMPXNI,TMPXNJ,TMPXND,TMPXNK,IDI,III,JJJ,KKK)
*$OMP1 COLLAPSE(2)
      DO ICEL=1,NCEL
      DO IG=1,NGEFF
      DO IPQD=1,NPQD(IND)
      PT => TMPXNK(IPQD,IG)

#if defined(_OPENMP)
        ITID=omp_get_thread_num()
#else
        ITID=0
#endif
      I=III(ICEL)
      J=JJJ(ICEL)
      K_MC=KKK(ICEL)

      IF(IMPX.GT.5) WRITE(IUNOUT,500) ITID,NGIND(IG),IPQD
      M=INDANG(IPQD,IND)
      IF(.NOT.INCONV(IG)) CYCLE
      IF(W(M).EQ.0.0) CYCLE

      ! Find in X-Y plane DRAGON geometry hexagonal index using I and J
      LHEX=(COORDMAP(1,:).EQ.I .AND. COORDMAP(2,:).EQ.J)
      IHEX_XY=0
      DO IIHEX=1,NHEX
        IF(LHEX(IIHEX)) THEN
          IHEX_XY=IIHEX
          EXIT
        ENDIF
      ENDDO
      IF(IHEX_XY.EQ.0) CALL XABORT('SNFDH3: IHEX_XY FAILURE.')
      ! Find D coordinate
      DCOORD = ABS(COORDMAP(3,IHEX_XY))-NRINGS

      IF(IDI.EQ.1)THEN
      ! Pick up boundary elements
      IF((NCODE(5).NE.1).or.(NCODE(6).NE.1))THEN
        IOF=(M-1)*(L5) + 1
        PT%MATRIX(:IELEM,:IELEM,:ISPLH,:ISPLH,:3,:NHEX)=
     >    RESHAPE(FUNKNO(L4+IOF:L4+IOF+L5,IG),
     >    (/IELEM,IELEM,ISPLH,ISPLH,3,NHEX/))
      ENDIF
      ! Account for albedo in boundary elements
      IF(IND.LT.7) THEN
        PT%MATRIX(:IELEM,:IELEM,:ISPLH,:ISPLH,:3,:NHEX)=
     >    PT%MATRIX(:IELEM,:IELEM,:ISPLH,:ISPLH,:3,:NHEX)*ZCODE(5)
      ELSE
        PT%MATRIX(:IELEM,:IELEM,:ISPLH,:ISPLH,:3,:NHEX)=
     >    PT%MATRIX(:IELEM,:IELEM,:ISPLH,:ISPLH,:3,:NHEX)*ZCODE(6)
      ENDIF
      ENDIF
*----
*  LOOP OVER Z-AXIS PLANES IN MACROCELL
*----
      DO IMZ=1,MIN(NCELLZ,LZ-(K_MC-1)*NCELLZ)
      K=(K_MC-1)*NCELLZ+IMZ
      IF(IND.GE.7) K=LZ+1-K
*----
*  LOOP OVER LOZENGES
*----
      DO ILOZLOOP=1,3
      ILOZ=LOZSWP(ILOZLOOP,IND_XY)

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
        IF((IND_XY.EQ.2).OR.(IND_XY.EQ.3).OR.(IND_XY.EQ.4))I2=ISPLH+1-I2
      ELSEIF(ILOZ.EQ.3)THEN
        IF((IND_XY.EQ.3).OR.(IND_XY.EQ.4).OR.(IND_XY.EQ.5))I2=ISPLH+1-I2
      ENDIF

      DO JL=1,ISPLH
      J2=JL
      ! Account for different sweep direction depending on angle 
      IF((ILOZ.EQ.2).OR.(ILOZ.EQ.3))THEN
        IF((IND_XY.EQ.4).OR.(IND_XY.EQ.5).OR.(IND_XY.EQ.6))J2=ISPLH+1-J2
      ELSEIF(ILOZ.EQ.1)THEN
        IF((IND_XY.EQ.3).OR.(IND_XY.EQ.4).OR.(IND_XY.EQ.5))J2=ISPLH+1-J2
      ENDIF

      FLUX(:IELEM**3,:NSCT)=0.0D0

      ! Skip loop if virtual element 
      IF(MAT(I2,J2,ILOZ,IHEX_XY,K).EQ.0) CYCLE

      ! Read in XNI and XNJ depending on lozenge
      I_FETCH=0
      IF((ILOZ.EQ.1))THEN
        ! Read boundary fluxes in reverse for lozenge A since affine
        ! transformation of lozenges causes the D and I directions
        ! of lozenges C and A respectively to be reversed
        I_FETCH=ISPLH+1-I2
        XNI(:,:) = TMPXNI(:,:,J2,J,K,IPQD,IG)
        XNJ(:,:) = TMPXND(:,:,I_FETCH,DCOORD,K,IPQD,IG)
      ELSEIF((ILOZ.EQ.2))THEN
        XNI(:,:) = TMPXNI(:,:,J2,J,K,IPQD,IG)
        XNJ(:,:) = TMPXNJ(:,:,I2,I,K,IPQD,IG)
      ELSEIF((ILOZ.EQ.3))THEN
        XNI(:,:) = TMPXND(:,:,J2,DCOORD,K,IPQD,IG)
        XNJ(:,:) = TMPXNJ(:,:,I2,I,K,IPQD,IG)
      ENDIF
      XNK(:,:) = PT%MATRIX(:,:,I2,J2,ILOZ,IHEX_XY)

      ! Compute adjusted direction cosines 
      MUHTEMP  =  DA(1,K,M)
      ETAHTEMP =  DB(1,K,M)
      MUH = (MUHTEMP*DDD) - (ETAHTEMP*BBB)
      ETAH = (-MUHTEMP*CCC) + (ETAHTEMP*AAA)
      XIH  = DC(1,1,M)

      ! Compute source and (vol*sigma) for current element 
      DO IEL=1,IELEM**3
      Q(IEL)=0.0D0
      DO P=1,NSCT
      IOF=(((((((K-1)*NHEX+(IHEX_XY-1))*3+(ILOZ-1))*ISPLH+
     >  (J2-1))*ISPLH+(I2-1))*NSCT+(P-1))*IELEM**3)+IEL
      Q(IEL)=Q(IEL)+QEXT(IOF,IG)*MN(M,P)
      ENDDO
      ENDDO
      VT=VOL(I2,J2,ILOZ,IHEX_XY,K)*TOTAL(MAT(I2,J2,ILOZ,IHEX_XY,K),IG)

      ! Initialise and compute matrix of coefficients, and source
      Q2(:IELEM**3,:IELEM**3+1)=0.0D0
      IF(IELEM.EQ.1) THEN
      Q2(1,1)=2.0D0*ABS(MUH)+2.0D0*ABS(ETAH)+
     >        2.0D0*ABS(XIH)+VT
      Q2(1,2)=2.0D0*ABS(MUH)*XNI(1,1)+
     >        2.0D0*ABS(ETAH)*XNJ(1,1)+
     >        2.0D0*ABS(XIH)*XNK(1,1)+VOL(I2,J2,ILOZ,IHEX_XY,K)*Q(1)
      ELSE IF(IELEM.EQ.2) THEN
      Q2(1,1) = VT
      Q2(1,2) = 2.0D0*3.0D0**(0.5D0)*MUH
      Q2(1,3) = 2.0D0*3.0D0**(0.5D0)*ETAH
      Q2(1,5) = 2.0D0*3.0D0**(0.5D0)*XIH
      Q2(2,1) = 2.0D0*3.0D0**(0.5D0)*MUH
      Q2(2,2) = - VT - 6.0D0*ABS(MUH)
      Q2(2,4) = -2.0D0*3.0D0**(0.5D0)*ETAH
      Q2(2,6) = -2.0D0*3.0D0**(0.5D0)*XIH
      Q2(3,1) = 2.0D0*3.0D0**(0.5D0)*ETAH
      Q2(3,3) = - VT - 6.0D0*ABS(ETAH)
      Q2(3,4) = -2.0D0*3.0D0**(0.5D0)*MUH
      Q2(3,7) = -2.0D0*3.0D0**(0.5D0)*XIH
      Q2(4,2) = -2.0D0*3.0D0**(0.5D0)*ETAH
      Q2(4,3) = -2.0D0*3.0D0**(0.5D0)*MUH
      Q2(4,4) = VT + 6.0D0*ABS(ETAH) + 6.0D0*ABS(MUH)
      Q2(4,8) = 2.0D0*3.0D0**(0.5D0)*XIH
      Q2(5,1) = 2.0D0*3.0D0**(0.5D0)*XIH
      Q2(5,5) = - VT - 6.0D0*ABS(XIH)
      Q2(5,6) = -2.0D0*3.0D0**(0.5D0)*MUH
      Q2(5,7) = -2.0D0*3.0D0**(0.5D0)*ETAH
      Q2(6,2) = -2.0D0*3.0D0**(0.5D0)*XIH
      Q2(6,5) = -2.0D0*3.0D0**(0.5D0)*MUH
      Q2(6,6) = VT + 6.0D0*ABS(MUH) + 6.0D0*ABS(XIH)
      Q2(6,8) = 2.0D0*3.0D0**(0.5D0)*ETAH
      Q2(7,3) = -2.0D0*3.0D0**(0.5D0)*XIH
      Q2(7,5) = -2.0D0*3.0D0**(0.5D0)*ETAH
      Q2(7,7) = VT + 6.0D0*ABS(ETAH) + 6.0D0*ABS(XIH)
      Q2(7,8) = 2.0D0*3.0D0**(0.5D0)*MUH
      Q2(8,4) = 2.0D0*3.0D0**(0.5D0)*XIH
      Q2(8,6) = 2.0D0*3.0D0**(0.5D0)*ETAH
      Q2(8,7) = 2.0D0*3.0D0**(0.5D0)*MUH
      Q2(8,8) = - VT - 6.0D0*ABS(ETAH) - 6.0D0*ABS
     >   (MUH) - 6.0D0*ABS(XIH)
*---------------
      Q2(1,9) = Q(1)*VOL(I2,J2,ILOZ,IHEX_XY,K)
      Q2(2,9) = -Q(2)*VOL(I2,J2,ILOZ,IHEX_XY,K)
      Q2(3,9) = -Q(3)*VOL(I2,J2,ILOZ,IHEX_XY,K)
      Q2(4,9) = Q(4)*VOL(I2,J2,ILOZ,IHEX_XY,K)
      Q2(5,9) = -Q(5)*VOL(I2,J2,ILOZ,IHEX_XY,K)
      Q2(6,9) = Q(6)*VOL(I2,J2,ILOZ,IHEX_XY,K)
      Q2(7,9) = Q(7)*VOL(I2,J2,ILOZ,IHEX_XY,K)
      Q2(8,9) = -Q(8)*VOL(I2,J2,ILOZ,IHEX_XY,K)
      Q2(1,9) = Q2(1,9)
      Q2(2,9) = Q2(2,9) + 2.0D0*3.0D0**(0.5D0)*MUH*XNI(1,1)
      Q2(3,9) = Q2(3,9) + 2.0D0*3.0D0**(0.5D0)*ETAH*XNJ(1,1)
      Q2(4,9) = Q2(4,9) - 2.0D0*3.0D0**(0.5D0)*ETAH*XNJ(2,1) -
     >  2.0D0*3.0D0**(0.5D0)*MUH*XNI(2,1)
      Q2(5,9) = Q2(5,9) + 2.0D0*3.0D0**(0.5D0)*XNK(1,1)*XIH
      Q2(6,9) = Q2(6,9) - 2.0D0*3.0D0**(0.5D0)*MUH*XNI(1,2) - 
     >  2.0D0*3.0D0**(0.5D0)*XNK(2,1)*XIH
      Q2(7,9) = Q2(7,9) - 2.0D0*3.0D0**(0.5D0)*ETAH*XNJ(1,2) - 
     >  2.0D0*3.0D0**(0.5D0)*XNK(1,2)*XIH
      Q2(8,9) = Q2(8,9) + 2.0D0*3.0D0**(0.5D0)*ETAH*XNJ(2,2) + 
     >  2.0D0*3.0D0**(0.5D0)*MUH*XNI(2,2) + 
     >  2.0D0*3.0D0**(0.5D0)*XNK(2,2)*XIH
      ELSE IF(IELEM.EQ.3) THEN
      Q2(1,1) = - VT - 2*ABS(ETAH) - 2*ABS(MUH)
     >    - 2*ABS(XIH)
      Q2(1,3) = -2*5**(0.5D0)*ABS(MUH)
      Q2(1,7) = -2*5**(0.5D0)*ABS(ETAH)
      Q2(1,19) = -2*5**(0.5D0)*ABS(XIH)
      Q2(2,2) = - VT - 2*ABS(ETAH) - 2*ABS(XIH)
      Q2(2,3) = -2*15**(0.5D0)*MUH
      Q2(2,8) = -2*5**(0.5D0)*ABS(ETAH)
      Q2(2,20) = -2*5**(0.5D0)*ABS(XIH)
      Q2(3,1) = 2*5**(0.5D0)*ABS(MUH)
      Q2(3,2) = -2*15**(0.5D0)*MUH
      Q2(3,3) = VT + 2*ABS(ETAH) + 10*ABS(MUH)
     >   + 2*ABS(XIH)
      Q2(3,9) = 2*5**(0.5D0)*ABS(ETAH)
      Q2(3,21) = 2*5**(0.5D0)*ABS(XIH)
      Q2(4,4) = - VT - 2*ABS(MUH) - 2*ABS(XIH)
      Q2(4,6) = -2*5**(0.5D0)*ABS(MUH)
      Q2(4,7) = -2*15**(0.5D0)*ETAH
      Q2(4,22) = -2*5**(0.5D0)*ABS(XIH)
      Q2(5,5) = - VT - 2*ABS(XIH)
      Q2(5,6) = -2*15**(0.5D0)*MUH
      Q2(5,8) = -2*15**(0.5D0)*ETAH
      Q2(5,23) = -2*5**(0.5D0)*ABS(XIH)
      Q2(6,4) = 2*5**(0.5D0)*ABS(MUH)
      Q2(6,5) = -2*15**(0.5D0)*MUH
      Q2(6,6) = VT + 10*ABS(MUH) + 2*ABS(XIH)
      Q2(6,9) = 2*15**(0.5D0)*ETAH
      Q2(6,24) = 2*5**(0.5D0)*ABS(XIH)
      Q2(7,1) = 2*5**(0.5D0)*ABS(ETAH)
      Q2(7,4) = -2*15**(0.5D0)*ETAH
      Q2(7,7) = VT + 10*ABS(ETAH) + 2*ABS(MUH)
     >   + 2*ABS(XIH)
      Q2(7,9) = 2*5**(0.5D0)*ABS(MUH)
      Q2(7,25) = 2*5**(0.5D0)*ABS(XIH)
      Q2(8,2) = 2*5**(0.5D0)*ABS(ETAH)
      Q2(8,5) = -2*15**(0.5D0)*ETAH
      Q2(8,8) = VT + 10*ABS(ETAH) + 2*ABS(XIH)
      Q2(8,9) = 2*15**(0.5D0)*MUH
      Q2(8,26) = 2*5**(0.5D0)*ABS(XIH)
      Q2(9,3) = -2*5**(0.5D0)*ABS(ETAH)
      Q2(9,6) = 2*15**(0.5D0)*ETAH
      Q2(9,7) = -2*5**(0.5D0)*ABS(MUH)
      Q2(9,8) = 2*15**(0.5D0)*MUH
      Q2(9,9) = - VT - 10*ABS(ETAH) - 10*ABS(MUH) 
     >  - 2*ABS(XIH)
      Q2(9,27) = -2*5**(0.5D0)*ABS(XIH)
      Q2(10,10) = VT + 2*ABS(ETAH) + 2*ABS(MUH)
      Q2(10,12) = 2*5**(0.5D0)*ABS(MUH)
      Q2(10,16) = 2*5**(0.5D0)*ABS(ETAH)
      Q2(10,19) = 2*15**(0.5D0)*XIH
      Q2(11,11) = VT + 2*ABS(ETAH)
      Q2(11,12) = 2*15**(0.5D0)*MUH
      Q2(11,17) = 2*5**(0.5D0)*ABS(ETAH)
      Q2(11,20) = 2*15**(0.5D0)*XIH
      Q2(12,10) = -2*5**(0.5D0)*ABS(MUH)
      Q2(12,11) = 2*15**(0.5D0)*MUH
      Q2(12,12) = - VT - 2*ABS(ETAH) - 10*ABS(MUH)
      Q2(12,18) = -2*5**(0.5D0)*ABS(ETAH)
      Q2(12,21) = -2*15**(0.5D0)*XIH
      Q2(13,13) = VT + 2*ABS(MUH)
      Q2(13,15) = 2*5**(0.5D0)*ABS(MUH)
      Q2(13,16) = 2*15**(0.5D0)*ETAH
      Q2(13,22) = 2*15**(0.5D0)*XIH
      Q2(14,14) = VT
      Q2(14,15) = 2*15**(0.5D0)*MUH
      Q2(14,17) = 2*15**(0.5D0)*ETAH
      Q2(14,23) = 2*15**(0.5D0)*XIH
      Q2(15,13) = -2*5**(0.5D0)*ABS(MUH)
      Q2(15,14) = 2*15**(0.5D0)*MUH
      Q2(15,15) = - VT - 10*ABS(MUH)
      Q2(15,18) = -2*15**(0.5D0)*ETAH
      Q2(15,24) = -2*15**(0.5D0)*XIH
      Q2(16,10) = -2*5**(0.5D0)*ABS(ETAH)
      Q2(16,13) = 2*15**(0.5D0)*ETAH
      Q2(16,16) = - VT - 10*ABS(ETAH) - 2*ABS(MUH)
      Q2(16,18) = -2*5**(0.5D0)*ABS(MUH)
      Q2(16,25) = -2*15**(0.5D0)*XIH
      Q2(17,11) = -2*5**(0.5D0)*ABS(ETAH)
      Q2(17,14) = 2*15**(0.5D0)*ETAH
      Q2(17,17) = - VT - 10*ABS(ETAH)
      Q2(17,18) = -2*15**(0.5D0)*MUH
      Q2(17,26) = -2*15**(0.5D0)*XIH
      Q2(18,12) = 2*5**(0.5D0)*ABS(ETAH)
      Q2(18,15) = -2*15**(0.5D0)*ETAH
      Q2(18,16) = 2*5**(0.5D0)*ABS(MUH)
      Q2(18,17) = -2*15**(0.5D0)*MUH
      Q2(18,18) = VT + 10*ABS(ETAH) + 10*ABS(MUH)
      Q2(18,27) = 2*15**(0.5D0)*XIH
      Q2(19,1) = 2*5**(0.5D0)*ABS(XIH)
      Q2(19,10) = -2*15**(0.5D0)*XIH
      Q2(19,19) = VT + 2*ABS(ETAH) + 2*ABS(MUH)
     >  + 10*ABS(XIH)
      Q2(19,21) = 2*5**(0.5D0)*ABS(MUH)
      Q2(19,25) = 2*5**(0.5D0)*ABS(ETAH)
      Q2(20,2) = 2*5**(0.5D0)*ABS(XIH)
      Q2(20,11) = -2*15**(0.5D0)*XIH
      Q2(20,20) = VT + 2*ABS(ETAH) + 10*ABS(XIH)
      Q2(20,21) = 2*15**(0.5D0)*MUH
      Q2(20,26) = 2*5**(0.5D0)*ABS(ETAH)
      Q2(21,3) = -2*5**(0.5D0)*ABS(XIH)
      Q2(21,12) = 2*15**(0.5D0)*XIH
      Q2(21,19) = -2*5**(0.5D0)*ABS(MUH)
      Q2(21,20) = 2*15**(0.5D0)*MUH
      Q2(21,21) = - VT - 2*ABS(ETAH) - 10*ABS(MUH) - 
     >  10*ABS(XIH)
      Q2(21,27) = -2*5**(0.5D0)*ABS(ETAH)
      Q2(22,4) = -2*5**(0.5D0)*ABS(XIH)
      Q2(22,13) = 2*15**(0.5D0)*XIH
      Q2(22,22) = - VT - 2*ABS(MUH) - 10*ABS(XIH)
      Q2(22,24) = -2*5**(0.5D0)*ABS(MUH)
      Q2(22,25) = -2*15**(0.5D0)*ETAH
      Q2(23,5) = -2*5**(0.5D0)*ABS(XIH)
      Q2(23,14) = 2*15**(0.5D0)*XIH
      Q2(23,23) = - VT - 10*ABS(XIH)
      Q2(23,24) = -2*15**(0.5D0)*MUH
      Q2(23,26) = -2*15**(0.5D0)*ETAH
      Q2(24,6) = 2*5**(0.5D0)*ABS(XIH)
      Q2(24,15) = -2*15**(0.5D0)*XIH
      Q2(24,22) = 2*5**(0.5D0)*ABS(MUH)
      Q2(24,23) = -2*15**(0.5D0)*MUH
      Q2(24,24) = VT + 10*ABS(MUH) + 10*ABS(XIH)
      Q2(24,27) = 2*15**(0.5D0)*ETAH
      Q2(25,7) = -2*5**(0.5D0)*ABS(XIH)
      Q2(25,16) = 2*15**(0.5D0)*XIH
      Q2(25,19) = -2*5**(0.5D0)*ABS(ETAH)
      Q2(25,22) = 2*15**(0.5D0)*ETAH
      Q2(25,25) = - VT - 10*ABS(ETAH) - 2*ABS(MUH) - 
     >   10*ABS(XIH)
      Q2(25,27) = -2*5**(0.5D0)*ABS(MUH)
      Q2(26,8) = 2*5**(0.5D0)*ABS(XIH)
      Q2(26,17) = -2*15**(0.5D0)*XIH
      Q2(26,20) = 2*5**(0.5D0)*ABS(ETAH)
      Q2(26,23) = -2*15**(0.5D0)*ETAH
      Q2(26,26) = VT + 10*ABS(ETAH) + 10*ABS(XIH)
      Q2(26,27) = 2*15**(0.5D0)*MUH
      Q2(27,9) = -2*5**(0.5D0)*ABS(XIH)
      Q2(27,18) = 2*15**(0.5D0)*XIH
      Q2(27,21) = -2*5**(0.5D0)*ABS(ETAH)
      Q2(27,24) = 2*15**(0.5D0)*ETAH
      Q2(27,25) = -2*5**(0.5D0)*ABS(MUH)
      Q2(27,26) = 2*15**(0.5D0)*MUH
      Q2(27,27) = - VT - 10*ABS(ETAH) - 10*ABS(
     >   MUH) - 10*ABS(XIH)
      Q2(1,28) = -Q(01)*VOL(I2,J2,ILOZ,IHEX_XY,K)
      Q2(2,28) = -Q(02)*VOL(I2,J2,ILOZ,IHEX_XY,K)
      Q2(3,28) = Q(03)*VOL(I2,J2,ILOZ,IHEX_XY,K)
      Q2(4,28) = -Q(04)*VOL(I2,J2,ILOZ,IHEX_XY,K)
      Q2(5,28) = -Q(05)*VOL(I2,J2,ILOZ,IHEX_XY,K)
      Q2(6,28) = Q(06)*VOL(I2,J2,ILOZ,IHEX_XY,K)
      Q2(7,28) = Q(07)*VOL(I2,J2,ILOZ,IHEX_XY,K)
      Q2(8,28) = Q(08)*VOL(I2,J2,ILOZ,IHEX_XY,K)
      Q2(9,28) = -Q(09)*VOL(I2,J2,ILOZ,IHEX_XY,K)
      Q2(10,28) = Q(10)*VOL(I2,J2,ILOZ,IHEX_XY,K)
      Q2(11,28) = Q(11)*VOL(I2,J2,ILOZ,IHEX_XY,K)
      Q2(12,28) = -Q(12)*VOL(I2,J2,ILOZ,IHEX_XY,K)
      Q2(13,28) = Q(13)*VOL(I2,J2,ILOZ,IHEX_XY,K)
      Q2(14,28) = Q(14)*VOL(I2,J2,ILOZ,IHEX_XY,K)
      Q2(15,28) = -Q(15)*VOL(I2,J2,ILOZ,IHEX_XY,K)
      Q2(16,28) = -Q(16)*VOL(I2,J2,ILOZ,IHEX_XY,K)
      Q2(17,28) = -Q(17)*VOL(I2,J2,ILOZ,IHEX_XY,K)
      Q2(18,28) = Q(18)*VOL(I2,J2,ILOZ,IHEX_XY,K)
      Q2(19,28) = Q(19)*VOL(I2,J2,ILOZ,IHEX_XY,K)
      Q2(20,28) = Q(20)*VOL(I2,J2,ILOZ,IHEX_XY,K)
      Q2(21,28) = -Q(21)*VOL(I2,J2,ILOZ,IHEX_XY,K)
      Q2(22,28) = -Q(22)*VOL(I2,J2,ILOZ,IHEX_XY,K)
      Q2(23,28) = -Q(23)*VOL(I2,J2,ILOZ,IHEX_XY,K)
      Q2(24,28) = Q(24)*VOL(I2,J2,ILOZ,IHEX_XY,K)
      Q2(25,28) = -Q(25)*VOL(I2,J2,ILOZ,IHEX_XY,K)
      Q2(26,28) = Q(26)*VOL(I2,J2,ILOZ,IHEX_XY,K)
      Q2(27,28) = -Q(27)*VOL(I2,J2,ILOZ,IHEX_XY,K)
      Q2(1,28) = Q2(1,28) - 2*XNK(1,1)*ABS(XIH) - 
     >  2*ABS(ETAH)*XNJ(1,1) - 2*ABS(MUH)*XNI(1,1)
      Q2(2,28) = Q2(2,28) - 2*XNK(2,1)*ABS(XIH) - 
     >  2*ABS(ETAH)*XNJ(2,1)
      Q2(3,28) = Q2(3,28) + 2*XNK(3,1)*ABS(XIH)
     >  + 2*ABS(ETAH)*XNJ(3,1) + 2*5**(0.5D0)
     >  *ABS(MUH)*XNI(1,1)
      Q2(4,28) = Q2(4,28) - 2*XNK(1,2)*ABS(XIH) - 
     >  2*ABS(MUH)*XNI(2,1)
      Q2(5,28) = Q2(5,28) -2*XNK(2,2)*ABS(XIH)
      Q2(6,28) = Q2(6,28) + 2*XNK(3,2)*ABS(XIH)
     >  + 2*5**(0.5D0)*ABS(MUH)*XNI(2,1)
      Q2(7,28) = Q2(7,28) + 2*XNK(1,3)*ABS(XIH)
     >  + 2*ABS(MUH)*XNI(3,1) + 2*5**(0.5D
     >  0)*ABS(ETAH)*XNJ(1,1)
      Q2(8,28) = Q2(8,28) + 2*XNK(2,3)*ABS(XIH)
     >  + 2*5**(0.5D0)*ABS(ETAH)*XNJ(2,1)
      Q2(9,28) = Q2(9,28) - 2*XNK(3,3)*ABS(XIH) 
     >  - 2*5**(0.5D0)*ABS(ETAH)*XNJ(3,1) -
     >  2*5**(0.5D0)*ABS(MUH)*XNI(3,1)
      Q2(10,28) = Q2(10,28) + 2*ABS(ETAH)*XNJ(1,2) + 
     >  2*ABS(MUH)*XNI(1,2)
      Q2(11,28) = Q2(11,28) + 2*ABS(ETAH)*XNJ(2,2)
      Q2(12,28) = Q2(12,28) - 2*ABS(ETAH)*XNJ
     >  (3,2) - 2*5**(0.5D0)*ABS(MUH)*XNI(1,2)
      Q2(13,28) = Q2(13,28) + 2*ABS(MUH)*XNI(2,2)
      Q2(14,28) = Q2(14,28) + 0
      Q2(15,28) = Q2(15,28) -2*5**(0.5D0)*ABS
     >  (MUH)*XNI(2,2)
      Q2(16,28) = Q2(16,28) - 2*ABS(MUH)*XNI(3,2) - 
     >  2*5**(0.5D0)*ABS(ETAH)*XNJ(1,2)
      Q2(17,28) = Q2(17,28) -2*5**(0.5D0)*ABS(ETAH)*XNJ(2,2)
      Q2(18,28) = Q2(18,28) + 2*5**(0.5D0)*ABS(ETAH)*XNJ(3,2) +
     >  2*5**(0.5D0)*ABS(MUH)*XNI(3,2)
      Q2(19,28) = Q2(19,28) + 2*ABS(ETAH)*XNJ(1,3) + 
     >  2*ABS(MUH)*XNI(1,3) + 2*5**(0.5D0)*XNK(1,1)*ABS(XIH)
      Q2(20,28) = Q2(20,28) + 2*ABS(ETAH)*XNJ(2,3) + 
     >  2*5**(0.5D0)*XNK(2,1)*ABS(XIH)
      Q2(21,28) = Q2(21,28) - 2*ABS(ETAH)*XNJ(3,3) - 
     >  2*5**(0.5D0)*ABS(MUH)*XNI(1,3)
     >  - 2*5**(0.5D0)*XNK(3,1)*ABS(XIH)
      Q2(22,28) = Q2(22,28) - 2*ABS(MUH)*XNI(2,3) - 
     >  2*5**(0.5D0)*XNK(1,2)*ABS(XIH)
      Q2(23,28) = Q2(23,28) -2*5**(0.5D0)*XNK(2,2)*ABS(XIH)
      Q2(24,28) = Q2(24,28) + 2*5**(0.5D0)*ABS(MUH)*XNI(2,3) +
     >  2*5**(0.5D0)*XNK(3,2)*ABS(XIH)
      Q2(25,28) = Q2(25,28) - 2*ABS(MUH)*XNI(3,3) - 
     >  2*5**(0.5D0)*ABS(ETAH)*XNJ(1,3) - 
     >  2*5**(0.5D0)*XNK(1,3)*ABS(XIH)
      Q2(26,28) = Q2(26,28) + 2*5**(0.5D0)*ABS(ETAH)*XNJ(2,3) + 
     >  2*5**(0.5D0)*XNK(2,3)*ABS(XIH)
      Q2(27,28) = Q2(27,28) - 2*5**(0.5D0)*ABS(ETAH)*XNJ(3,3) -
     >  2*5**(0.5D0)*ABS(MUH)*XNI(3,3) - 2*5**(0.5D0)*XNK(3,3)
     >  *ABS(XIH)
      ENDIF

      ! Solve system of equation given by matrix of coefficients
      ! and source
      CALL ALSBD(IELEM**3,1,Q2,IER,IELEM**3)
      IF(IER.NE.0) CALL XABORT('SNFDH3: SINGULAR MATRIX.')

      ! Assign K-boundary fluxes
      IF(IELEM.EQ.1) THEN
      PT%MATRIX(1,1,I2,J2,ILOZ,IHEX_XY)=2.0D0*Q2(1,2)-XNK(1,1)
      ELSE IF(IELEM.EQ.2) THEN
      PT%MATRIX(1,1,I2,J2,ILOZ,IHEX_XY) = XNK(1,1) + 
     >  SIGN(1.0,XIH)*2.0D0*3.0D0**(0.5D0)*Q2(5,9)
      PT%MATRIX(2,1,I2,J2,ILOZ,IHEX_XY) = XNK(2,1) + 
     >  SIGN(1.0,XIH)*2.0D0*3.0D0**(0.5D0)*Q2(6,9)
      PT%MATRIX(1,2,I2,J2,ILOZ,IHEX_XY) = XNK(1,2) + 
     >  SIGN(1.0,XIH)*2.0D0*3.0D0**(0.5D0)*Q2(7,9)
      PT%MATRIX(2,2,I2,J2,ILOZ,IHEX_XY) = XNK(2,2) + 
     >  SIGN(1.0,XIH)*2.0D0*3.0D0**(0.5D0)*Q2(8,9)
      ELSE IF(IELEM.EQ.3) THEN
      PT%MATRIX(1,1,I2,J2,ILOZ,IHEX_XY) = 2*Q2(01,28) - 
     >  XNK(1,1) + 2*5**(0.5D0)*Q2(19,28)
      PT%MATRIX(2,1,I2,J2,ILOZ,IHEX_XY) = 2*Q2(02,28) - 
     >  XNK(2,1) + 2*5**(0.5D0)*Q2(20,28)
      PT%MATRIX(3,1,I2,J2,ILOZ,IHEX_XY) = 2*Q2(03,28) - 
     >  XNK(3,1) + 2*5**(0.5D0)*Q2(21,28)
      PT%MATRIX(1,2,I2,J2,ILOZ,IHEX_XY) = 2*Q2(04,28) - 
     >  XNK(1,2) + 2*5**(0.5D0)*Q2(22,28)
      PT%MATRIX(2,2,I2,J2,ILOZ,IHEX_XY) = 2*Q2(05,28) - 
     >  XNK(2,2) + 2*5**(0.5D0)*Q2(23,28)
      PT%MATRIX(3,2,I2,J2,ILOZ,IHEX_XY) = 2*Q2(06,28) - 
     >  XNK(3,2) + 2*5**(0.5D0)*Q2(24,28)
      PT%MATRIX(1,3,I2,J2,ILOZ,IHEX_XY) = 2*Q2(07,28) - 
     >  XNK(1,3) + 2*5**(0.5D0)*Q2(25,28)
      PT%MATRIX(2,3,I2,J2,ILOZ,IHEX_XY) = 2*Q2(08,28) - 
     >  XNK(2,3) + 2*5**(0.5D0)*Q2(26,28)
      PT%MATRIX(3,3,I2,J2,ILOZ,IHEX_XY) = 2*Q2(09,28) - 
     >  XNK(3,3) + 2*5**(0.5D0)*Q2(27,28)
      ENDIF

      ! Assign J-boundary fluxes if lozenges B or C
      IF((ILOZ.EQ.2).OR.(ILOZ.EQ.3))THEN
      IF(IELEM.EQ.1) THEN
      TMPXNJ(1,1,I2,I,K,IPQD,IG)=2.0D0*Q2(1,2)-XNJ(1,1)
      ELSEIF(IELEM.EQ.2) THEN
      TMPXNJ(1,1,I2,I,K,IPQD,IG) = XNJ(1,1) + SIGN(1.0,ETAH)*
     >  2.0D0*3.0D0**(0.5D0)*Q2(3,9)
      TMPXNJ(2,1,I2,I,K,IPQD,IG) = XNJ(2,1) + SIGN(1.0,ETAH)*
     >  2.0D0*3.0D0**(0.5D0)*Q2(4,9)
      TMPXNJ(1,2,I2,I,K,IPQD,IG) = XNJ(1,2) + SIGN(1.0,ETAH)*
     >  2.0D0*3.0D0**(0.5D0)*Q2(7,9)
      TMPXNJ(2,2,I2,I,K,IPQD,IG) = XNJ(2,2) + SIGN(1.0,ETAH)*
     >  2.0D0*3.0D0**(0.5D0)*Q2(8,9)
      ELSEIF(IELEM.EQ.3) THEN
      TMPXNJ(1,1,I2,I,K,IPQD,IG) = 2*Q2(01,28) - XNJ(1,1) +
     >  2*5**(0.5D0)*Q2(07,28)
      TMPXNJ(2,1,I2,I,K,IPQD,IG) = 2*Q2(02,28) - XNJ(2,1) +
     >  2*5**(0.5D0)*Q2(08,28)
      TMPXNJ(3,1,I2,I,K,IPQD,IG) = 2*Q2(03,28) - XNJ(3,1) +
     >  2*5**(0.5D0)*Q2(09,28)
      TMPXNJ(1,2,I2,I,K,IPQD,IG) = 2*Q2(10,28) - XNJ(1,2) +
     >  2*5**(0.5D0)*Q2(16,28)
      TMPXNJ(2,2,I2,I,K,IPQD,IG) = 2*Q2(11,28) - XNJ(2,2) +
     >  2*5**(0.5D0)*Q2(17,28)
      TMPXNJ(3,2,I2,I,K,IPQD,IG) = 2*Q2(12,28) - XNJ(3,2) +
     >  2*5**(0.5D0)*Q2(18,28)
      TMPXNJ(1,3,I2,I,K,IPQD,IG) = 2*Q2(19,28) - XNJ(1,3) +
     >  2*5**(0.5D0)*Q2(25,28)
      TMPXNJ(2,3,I2,I,K,IPQD,IG) = 2*Q2(20,28) - XNJ(2,3) +
     >  2*5**(0.5D0)*Q2(26,28)
      TMPXNJ(3,3,I2,I,K,IPQD,IG) = 2*Q2(21,28) - XNJ(3,3) +
     >  2*5**(0.5D0)*Q2(27,28)
      ENDIF
      ENDIF
      ! Assign I-boundary fluxes if lozenges A or B
      IF((ILOZ.EQ.1).OR.(ILOZ.EQ.2))THEN
      IF(IELEM.EQ.1) THEN
      TMPXNI(1,1,J2,J,K,IPQD,IG) = 2.0D0*Q2(1,2)-XNI(1,1)
      ELSEIF(IELEM.EQ.2) THEN
      TMPXNI(1,1,J2,J,K,IPQD,IG) = XNI(1,1) + SIGN(1.0,MUH)*
     >  2.0D0*3.0D0**(0.5D0)*Q2(2,9)
      TMPXNI(2,1,J2,J,K,IPQD,IG) = XNI(2,1) + SIGN(1.0,MUH)*
     >  2.0D0*3.0D0**(0.5D0)*Q2(4,9)
      TMPXNI(1,2,J2,J,K,IPQD,IG) = XNI(1,2) + SIGN(1.0,MUH)*
     >  2.0D0*3.0D0**(0.5D0)*Q2(6,9)
      TMPXNI(2,2,J2,J,K,IPQD,IG) = XNI(2,2) + SIGN(1.0,MUH)*
     >  2.0D0*3.0D0**(0.5D0)*Q2(8,9)
      ELSEIF(IELEM.EQ.3) THEN
      TMPXNI(1,1,J2,J,K,IPQD,IG) = 2*Q2(01,28) - XNI(1,1)
     >  + 2*5**(0.5D0)*Q2(03,28)
      TMPXNI(2,1,J2,J,K,IPQD,IG) = 2*Q2(04,28) - XNI(2,1)
     >  + 2*5**(0.5D0)*Q2(06,28)
      TMPXNI(3,1,J2,J,K,IPQD,IG) = 2*Q2(07,28) - XNI(3,1)
     >  + 2*5**(0.5D0)*Q2(09,28)
      TMPXNI(1,2,J2,J,K,IPQD,IG) = 2*Q2(10,28) - XNI(1,2)
     >  + 2*5**(0.5D0)*Q2(12,28)
      TMPXNI(2,2,J2,J,K,IPQD,IG) = 2*Q2(13,28) - XNI(2,2)
     >  + 2*5**(0.5D0)*Q2(15,28)
      TMPXNI(3,2,J2,J,K,IPQD,IG) = 2*Q2(16,28) - XNI(3,2)
     >  + 2*5**(0.5D0)*Q2(18,28)
      TMPXNI(1,3,J2,J,K,IPQD,IG) = 2*Q2(19,28) - XNI(1,3)
     >  + 2*5**(0.5D0)*Q2(21,28)
      TMPXNI(2,3,J2,J,K,IPQD,IG) = 2*Q2(22,28) - XNI(2,3)
     >  + 2*5**(0.5D0)*Q2(24,28)
      TMPXNI(3,3,J2,J,K,IPQD,IG) = 2*Q2(25,28) - XNI(3,3)
     >  + 2*5**(0.5D0)*Q2(27,28)
      ENDIF
      ENDIF
      ! Assign D-boundary fluxes if lozenge A using XNJ
      IF((ILOZ.EQ.1))THEN
      C1=1.0D0
      I3=I_FETCH
      IF(JL.EQ.ISPLH) C1=-1.0D0
      IF(IELEM.EQ.1) THEN
      TMPXND(1,1,I3,DCOORD,K,IPQD,IG)=2.0D0*Q2(1,2)-XNJ(1,1)
      ELSEIF(IELEM.EQ.2) THEN
      TMPXND(1,1,I3,DCOORD,K,IPQD,IG) = XNJ(1,1) + SIGN(1.0,ETAH)*
     >  2.0D0*3.0D0**(0.5D0)*Q2(3,9)
      TMPXND(2,1,I3,DCOORD,K,IPQD,IG) = (XNJ(2,1) + SIGN(1.0,ETAH)*
     >  2.0D0*3.0D0**(0.5D0)*Q2(4,9))*C1
      TMPXND(1,2,I3,DCOORD,K,IPQD,IG) = XNJ(1,2) + SIGN(1.0,ETAH)*
     >  2.0D0*3.0D0**(0.5D0)*Q2(7,9)
      TMPXND(2,2,I3,DCOORD,K,IPQD,IG) = (XNJ(2,2) + SIGN(1.0,ETAH)*
     >  2.0D0*3.0D0**(0.5D0)*Q2(8,9))*C1
      ELSEIF(IELEM.EQ.3) THEN
      TMPXND(1,1,I3,DCOORD,K,IPQD,IG) = 2*Q2(01,28) - XNJ(1,1) +
     >  2*5**(0.5D0)*Q2(07,28)
      TMPXND(2,1,I3,DCOORD,K,IPQD,IG) = (2*Q2(02,28) - XNJ(2,1) +
     >  2*5**(0.5D0)*Q2(08,28))*C1
      TMPXND(3,1,I3,DCOORD,K,IPQD,IG) = 2*Q2(03,28) - XNJ(3,1) +
     >  2*5**(0.5D0)*Q2(09,28)
      TMPXND(1,2,I3,DCOORD,K,IPQD,IG) = 2*Q2(10,28) - XNJ(1,2) +
     >  2*5**(0.5D0)*Q2(16,28)
      TMPXND(2,2,I3,DCOORD,K,IPQD,IG) = (2*Q2(11,28) - XNJ(2,2) +
     >  2*5**(0.5D0)*Q2(17,28))*C1
      TMPXND(3,2,I3,DCOORD,K,IPQD,IG) = 2*Q2(12,28) - XNJ(3,2) +
     >  2*5**(0.5D0)*Q2(18,28)
      TMPXND(1,3,I3,DCOORD,K,IPQD,IG) = 2*Q2(19,28) - XNJ(1,3) +
     >  2*5**(0.5D0)*Q2(25,28)
      TMPXND(2,3,I3,DCOORD,K,IPQD,IG) = (2*Q2(20,28) - XNJ(2,3) +
     >  2*5**(0.5D0)*Q2(26,28))*C1
      TMPXND(3,3,I3,DCOORD,K,IPQD,IG) = 2*Q2(21,28) - XNJ(3,3) +
     >  2*5**(0.5D0)*Q2(27,28)
      ENDIF
      ENDIF
      ! Assign D-boundary fluxes if lozenge C using XNI
      IF((ILOZ.EQ.3))THEN
      C1=1.0D0
      IF(IL.EQ.ISPLH) C1=-1.0D0
      IF(IELEM.EQ.1) THEN
      TMPXND(1,1,J2,DCOORD,K,IPQD,IG) = 2.0D0*Q2(1,2)-XNI(1,1)
      ELSEIF(IELEM.EQ.2) THEN
      TMPXND(1,1,J2,DCOORD,K,IPQD,IG) = XNI(1,1) + SIGN(1.0,MUH)*
     >  2.0D0*3.0D0**(0.5D0)*Q2(2,9)
      TMPXND(2,1,J2,DCOORD,K,IPQD,IG) = (XNI(2,1) + SIGN(1.0,MUH)*
     >  2.0D0*3.0D0**(0.5D0)*Q2(4,9))*C1
      TMPXND(1,2,J2,DCOORD,K,IPQD,IG) = XNI(1,2) + SIGN(1.0,MUH)*
     >  2.0D0*3.0D0**(0.5D0)*Q2(6,9)
      TMPXND(2,2,J2,DCOORD,K,IPQD,IG) = (XNI(2,2) + SIGN(1.0,MUH)*
     >  2.0D0*3.0D0**(0.5D0)*Q2(8,9))*C1
      ELSEIF(IELEM.EQ.3) THEN
      TMPXND(1,1,J2,DCOORD,K,IPQD,IG) = 2*Q2(01,28) - XNI(1,1)
     >  + 2*5**(0.5D0)*Q2(03,28)
      TMPXND(2,1,J2,DCOORD,K,IPQD,IG) = (2*Q2(04,28) - XNI(2,1)
     >  + 2*5**(0.5D0)*Q2(06,28))*C1
      TMPXND(3,1,J2,DCOORD,K,IPQD,IG) = 2*Q2(07,28) - XNI(3,1)
     >  + 2*5**(0.5D0)*Q2(09,28)
      TMPXND(1,2,J2,DCOORD,K,IPQD,IG) = 2*Q2(10,28) - XNI(1,2)
     >  + 2*5**(0.5D0)*Q2(12,28)
      TMPXND(2,2,J2,DCOORD,K,IPQD,IG) = (2*Q2(13,28) - XNI(2,2)
     >  + 2*5**(0.5D0)*Q2(15,28))*C1
      TMPXND(3,2,J2,DCOORD,K,IPQD,IG) = 2*Q2(16,28) - XNI(3,2)
     >  + 2*5**(0.5D0)*Q2(18,28)
      TMPXND(1,3,J2,DCOORD,K,IPQD,IG) = 2*Q2(19,28) - XNI(1,3)
     >  + 2*5**(0.5D0)*Q2(21,28)
      TMPXND(2,3,J2,DCOORD,K,IPQD,IG) = (2*Q2(22,28) - XNI(2,3)
     >  + 2*5**(0.5D0)*Q2(24,28))*C1
      TMPXND(3,3,J2,DCOORD,K,IPQD,IG) = 2*Q2(25,28) - XNI(3,3)
     >  + 2*5**(0.5D0)*Q2(27,28)
      ENDIF
      ENDIF

      ! Compute flux
      DO P=1,NSCT
      DO IEL=1,IELEM**3
      FLUX(IEL,P)=FLUX(IEL,P)+Q2(IEL,IELEM**3+1)*DN(P,M)
      ENDDO
      ENDDO

      IOF=((ILOZ-1)*ISPLH+(J2-1))*ISPLH+I2
      FLUX_G(:,:,IOF,IHEX_XY,K,IG)=FLUX_G(:,:,IOF,IHEX_XY,K,IG)+
     >  FLUX(:,:)

      ENDDO ! jl, isplh
      ENDDO ! il, isplh

      ENDDO ! ilozloop
      ENDDO ! lz

      ! Save K-boundary fluxes if not VOID boundary conditions
      IF(IDI.EQ.NCOL+NCOL+NMZ-2)THEN
      IF((NCODE(5).NE.1).or.(NCODE(6).NE.1))THEN
        IOF=(M-1)*(L5)
        FUNKNO(L4+IOF+1:L4+IOF+L5,IG)=
     >  RESHAPE(REAL(
     >  PT%MATRIX(:IELEM,:IELEM,:ISPLH,:ISPLH,:3,:NHEX)),(/L5/))
      ENDIF
      ENDIF

      ENDDO ! ipqd
      ENDDO ! ig
      ENDDO ! icel
*$OMP END PARALLEL DO

      DEALLOCATE(JJJ,III,KKK)
      ENDDO ! idi
      ENDDO ! jnd

      DO IG=1,NGEFF
        IF(.NOT.INCONV(IG)) CYCLE
        FUNKNO(:L4,IG)=
     >  RESHAPE(REAL(FLUX_G(:IELEM**3,:NSCT,:3*ISPLH**2,:NHEX,:LZ,IG)),
     >  (/L4/))
      ENDDO
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DO IG=1,NGEFF
        DO IPQD=1,NPQ
          DEALLOCATE(TMPXNK(IPQD,IG)%MATRIX)
        ENDDO
      ENDDO
      DEALLOCATE(TMPMAT,INDANG,FLUX,FLUX_G,TMPXNI,TMPXNJ,TMPXND,TMPXNK)
      RETURN
  500 FORMAT(16H SNFDH3: thread=,I8,12H --->(group=,I4,7H angle=,I4,1H))
      END
