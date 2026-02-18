*DECK SNFBH3
      SUBROUTINE SNFBH3(NUN,NGEFF,IMPX,INCONV,NGIND,NHEX,LZ,ISPLH,SIDE,
     1 IELEM,NM,NMX,NMY,NMZ,NMAT,NPQ,NSCT,MAT,VOL,TOTAL,NCODE,ZCODE,
     2 QEXT,LFIXUP,DU,DE,DZ,W,MRMZ,DC,DB,DA,MN,DN,WX,WY,WZ,CST,ISADPT,
     3 LOZSWP,COORDMAP,FUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform one inner iteration for solving SN equations in 3D Cartesian
* geometry for the HODD method. Energy-angle multithreading. Albedo
* boundary conditions on top/bottom, Void on sides. Boltzmann (BTE) 
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
* LZ       number of meshes along Z axis.
* IELEM    measure of order of the spatial approximation polynomial:
*          =1 constant - default for HODD;
*          =2 linear - default for DG;
*          >3 higher orders.
* NM       number of moments in space for flux components
* NMX      number of moments for X axis boundaries components
* NMY      number of moments for Y axis boundaries components
* NMZ      number of moments for Z axis boundaries components
* NMAT     number of material mixtures.
* NPQ      number of SN directions in height octants.
* NSCT     maximum number of spherical harmonics moments of the flux.
* MAT      material mixture index in each region.
* VOL      volumes of each region.
* TOTAL    macroscopic total cross sections.
* ESTOPW   stopping power.
* NCODE    boundary condition indices.
* ZCODE    albedos.
* DELTAE   energy group width in MeV.
* QEXT     Legendre components of the fixed source.
* LFIXUP   flag to enable negative flux fixup.
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
* WX       spatial X axis closure relation weighting factors.
* WY       spatial Y axis closure relation weighting factors.
* WZ       spatial Z axis closure relation weighting factors.
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
*
#if defined(_OPENMP)
      USE omp_lib
#endif
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NUN,NGEFF,IMPX,NGIND(NGEFF),NHEX,LZ,ISPLH,IELEM,NM,NMX,
     1 NMY,NMZ,NMAT,NPQ,NSCT,MAT(ISPLH,ISPLH,3,NHEX,LZ),NCODE(6),
     2 MRMZ(NPQ),LOZSWP(3,6),COORDMAP(3,NHEX)
      LOGICAL INCONV(NGEFF)
      REAL SIDE,VOL(ISPLH,ISPLH,3,NHEX,LZ),TOTAL(0:NMAT,NGEFF),ZCODE(6),
     1 QEXT(NUN,NGEFF),DU(NPQ),DE(NPQ),DZ(NPQ),W(NPQ),
     2 DC(ISPLH*ISPLH*3*NHEX,1,NPQ),DB(ISPLH*ISPLH*3*NHEX,LZ,NPQ),
     3 DA(1,LZ,NPQ),FUNKNO(NUN,NGEFF),WX(IELEM+1),WY(IELEM+1),
     4 WZ(IELEM+1),CST(IELEM),MN(NPQ,NSCT),DN(NSCT,NPQ)
      LOGICAL LFIXUP,ISADPT(3)
*----
*  LOCAL VARIABLES
*----
      INTEGER :: NPQD(12),IIND(12),P,DCOORD
      REAL :: JAC(2,2,3),MUH,ETAH,XIH,AAA,BBB,CCC,DDD,MUHTEMP,ETAHTEMP,
     1 WX0(IELEM+1),WY0(IELEM+1),WZ0(IELEM+1)
      DOUBLE PRECISION :: V,Q(NM),Q2(NM,NM+1),THETA,XNI(NMX),
     > XNJ(NMY),XNK(NMZ)
      PARAMETER(IUNOUT=6,RLOG=1.0E-8,PI=3.141592654)
      LOGICAL ISFIX(3),LHEX(NHEX)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: INDANG
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: TMPMAT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: FLUX
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: FLUX_G
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: TMPXNI,
     > TMPXNJ, TMPXND
      DOUBLE PRECISION,ALLOCATABLE, DIMENSION(:,:,:,:,:) :: TMPXNK
*----
*  MAP MATERIAL VALUES TO CARTESIAN AXIAL COORDINATE MAP
*----
      NRINGS=INT((SQRT(  REAL((4*NHEX-1)/3)  )+1.)/2.)
      NCOL=2*NRINGS -1
      ALLOCATE(TMPMAT(ISPLH,ISPLH,3,NCOL,NCOL,LZ))
      TMPMAT(:,:,:,:,:,:) = -1
      DO K=1,LZ
        DO IHEX_XY=1,NHEX
          TMPMAT(:,:,:,COORDMAP(1,IHEX_XY),COORDMAP(2,IHEX_XY),K) = 
     >         MAT(:,:,:,IHEX_XY,K)
        ENDDO
      ENDDO
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(INDANG(NPQ,12))
      ALLOCATE(FLUX(NM,NSCT,3*ISPLH**2,NHEX,LZ))
      ALLOCATE(FLUX_G(NM,NSCT,3*ISPLH**2,NHEX,LZ,NGEFF))
      ALLOCATE(TMPXNI(NMX,ISPLH,NCOL))
      ALLOCATE(TMPXNJ(NMY,ISPLH,NCOL))
      ALLOCATE(TMPXND(NMX,ISPLH,NCOL))
      ALLOCATE(TMPXNK(NMZ,ISPLH,ISPLH,3,NHEX))
*----
*  CONSTRUCT JACOBIAN MATRIX FOR EACH LOZENGE
*----
      JAC = RESHAPE((/ 1., -SQRT(3.), 1., SQRT(3.), 2., 0., 1.,
     >    SQRT(3.), 2., 0., -1., SQRT(3.) /), SHAPE(JAC))
      JAC = (SIDE/2.)*JAC
*----
*  LENGTH OF FUNKNO COMPONENTS (IN ORDER)
*----
      LFLX=3*NM*(ISPLH**2)*NHEX*LZ*NSCT
      L5=3*NMZ*(ISPLH**2)*NHEX
*----
*  SET DODECANT SWAPPING ORDER.
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
*  MAIN LOOP OVER DODECANTS.
*----

      FLUX_G(:NM,:NSCT,:3*ISPLH**2,:NHEX,:LZ,:NGEFF)=0.0D0
      WX0=WX
      WY0=WY
      WZ0=WZ

      DO JND=1,12
      IND=IIND(JND)
      IND_XY=MOD(IND-1,6)+1
      ! Needed because of S2 LS (8 dir. for 12 dodecants)
      IF(IND.EQ.0) CYCLE
*----
*  PRELIMINARY LOOPS FOR SETTING BOUNDARY CONDITIONS.
*----

      IF((NCODE(5).NE.1).or.(NCODE(6).NE.1))THEN
*$OMP  PARALLEL DO
*$OMP+ PRIVATE(M,IG,VZ,M1,IOF,JOF,IPQD)
*$OMP+ SHARED(FUNKNO) COLLAPSE(2)

      DO IG=1,NGEFF
      DO IPQD=1,NPQD(IND)
      IF(.NOT.INCONV(IG)) CYCLE
      M=INDANG(IPQD,IND)

      VZ=DZ(M)
      ! Z-BOUNDARY
      IF(VZ.GT.0.0)THEN
        M1=MRMZ(M)
        IF(NCODE(5).NE.4)THEN
          IOF=(M-1)*(L5)
          JOF=(M1-1)*(L5)
          FUNKNO(LFLX+IOF+1:LFLX+IOF+L5,IG)=
     >      FUNKNO(LFLX+JOF+1:LFLX+JOF+L5,IG)
        ENDIF
      ELSEIF(VZ.LT.0.0)THEN
        M1=MRMZ(M)
        IF(NCODE(6).NE.4)THEN
          IOF=(M-1)*(L5)
          JOF=(M1-1)*(L5)
          FUNKNO(LFLX+IOF+1:LFLX+IOF+L5,IG)=
     >      FUNKNO(LFLX+JOF+1:LFLX+JOF+L5,IG)
        ENDIF
      ENDIF
*
      ENDDO
      ENDDO

*$OMP END PARALLEL DO
      ENDIF

      ! CALL XABORT('SNFBH3: testing 1 ')
*----
*  MAIN SWAPPING LOOPS FOR SN FLUX CALCULATION
*  LOOP OVER ENERGY AND ANGLES
*----

*$OMP  PARALLEL DO 
*$OMP+ PRIVATE(ITID,FLUX,M,IG,XNI,XNJ,XNK,Q,Q2,IOF,IER,II,JJ,I,J,K,K0)
*$OMP+ PRIVATE(IPQD,IBM,SIGMA,V,ISFIX,IX,JX,IY,JY,IZ,JZ,AAA,BBB,CCC,DDD)
*$OMP+ PRIVATE(IIX,IIY,IIZ,P,IEL)
*$OMP+ PRIVATE(LHEX,IHEX_XY,IIHEX,DCOORD,ILOZLOOP,ILOZ,IL,I2,JL,J2)
*$OMP+ PRIVATE(MUHTEMP,MUH,ETAHTEMP,ETAH,XIH,I3,I_FETCH,III,JJJ,IIM,JIM)
*$OMP+ PRIVATE(TMPXNI,TMPXNJ,TMPXND,TMPXNK)
*$OMP+ FIRSTPRIVATE(WX,WY,WZ,WX0,WY0,WZ0) SHARED(FUNKNO,IND)
*$OMP+ REDUCTION(+:FLUX_G) COLLAPSE(2)      
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
      IF(IMPX.GT.5) WRITE(IUNOUT,500) ITID,NGIND(IG),IPQD

      ! INITIALIZE FLUX AND Z-BOUNDARY FLUXES
      FLUX(:NM,:NSCT,:3*ISPLH**2,:NHEX,:LZ)=0.0D0
      TMPXNK(:NMZ,:ISPLH,:ISPLH,:3,:NHEX)=0.0D0

      ! PICK UP BOUNDARY ELEMENTS
      IF((NCODE(5).NE.1).or.(NCODE(6).NE.1))THEN
        IOF=(M-1)*(L5) + 1
        TMPXNK(:NMZ,:ISPLH,:ISPLH,:3,:NHEX)=
     >    RESHAPE(FUNKNO(LFLX+IOF:LFLX+IOF+L5,IG),
     >    (/NMZ,ISPLH,ISPLH,3,NHEX/))
      ENDIF
      ! ACCOUNT FOR ALBEDO IN BOUNDARY ELEMENTS
      IF(IND.LT.7) THEN
        TMPXNK(:NMZ,:ISPLH,:ISPLH,:3,:NHEX)=
     >    TMPXNK(:NMZ,:ISPLH,:ISPLH,:3,:NHEX)*ZCODE(5)
      ELSE
        TMPXNK(:NMZ,:ISPLH,:ISPLH,:3,:NHEX)=
     >    TMPXNK(:NMZ,:ISPLH,:ISPLH,:3,:NHEX)*ZCODE(6)
      ENDIF

*----
*  LOOP OVER Z-AXIS PLANES
*----

      DO K0=1,LZ
      K=K0
      IF(IND.GE.7) K=LZ+1-K0

      ! INITIALIZE I,J,D-BOUNDARY FLUXES
      TMPXNI(:NMX,:ISPLH,:NCOL)=0.0D0
      TMPXNJ(:NMY,:ISPLH,:NCOL)=0.0D0
      TMPXND(:NMX,:ISPLH,:NCOL)=0.0D0

*----
*  LOOP OVER CARTESIAN MAP OF HEXAGONAL DOMAIN
*----

      DO JJJ=1,NCOL
      JIM=JJJ
      ! Account for different sweep direction depending on angle 
      IF((IND_XY.EQ.1).OR.(IND_XY.EQ.2).OR.(IND_XY.EQ.3)) JIM=NCOL+1-JIM
        
      DO III=1,NCOL
      IIM=III
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
      IF(TMPMAT(1,1,1,I,J,K).EQ.-1) CYCLE

      ! Find in X-Y plane DRAGON geometry hexagonal index using I and J
      LHEX=(COORDMAP(1,:).EQ.I .AND. COORDMAP(2,:).EQ.J)
      IHEX_XY=0
      DO IIHEX=1,NHEX
        IF(LHEX(IIHEX)) THEN
          IHEX_XY=IIHEX
          EXIT
        ENDIF
      ENDDO
      IF(IHEX_XY.EQ.0) CALL XABORT('SNFBH3: IHEX_XY FAILURE.')
      ! Find D coordinate
      DCOORD = ABS(COORDMAP(3,IHEX_XY))-NRINGS

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

      ! CALL XABORT('SNFBH3: testing 19 ')
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
      ! READ IN XNI AND XNJ DEPENDING ON LOZENGE
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
      XNK(:) = TMPXNK(:,I2,J2,ILOZ,IHEX_XY)

      ! DATA
      IBM=MAT(I2,J2,ILOZ,IHEX_XY,K)
      ! Skip loop if virtual element
      IF(IBM.EQ.0) CYCLE
      SIGMA=TOTAL(IBM,IG)
      V=VOL(I2,J2,ILOZ,IHEX_XY,K)

      ! COMPUTE ADJUSTED DIRECTION COSINES
      MUHTEMP  =  DA(1,K,M)
      ETAHTEMP =  DB(1,K,M)
      MUH = (MUHTEMP*DDD) - (ETAHTEMP*BBB)
      ETAH = (-MUHTEMP*CCC) + (ETAHTEMP*AAA)
      XIH  = DC(1,1,M)

      ! IF(IND.EQ.12) CALL XABORT('SNFBH3: testing 60 ')
    !   WRITE(*,*) (((((((K-1)*NHEX+(IHEX_XY-1))*3+(ILOZ-1))*ISPLH+
    !  >  (J2-1))*ISPLH+(I2-1))*NSCT+(2-1))*NM)+1
    !   WRITE(*,*) K, NHEX, IHEX_XY, ILOZ, ISPLH, J2, I2, NSCT, NM
    !   CALL XABORT('SNFBH3: testing 2')
      ! SOURCE DENSITY TERM
      DO IEL=1,NM
      Q(IEL)=0.0D0
      DO P=1,NSCT
      IOF=(((((((K-1)*NHEX+(IHEX_XY-1))*3+(ILOZ-1))*ISPLH+(J2-1))*
     >  ISPLH+(I2-1))*NSCT+(P-1))*NM)+IEL
      Q(IEL)=Q(IEL)+QEXT(IOF,IG)*MN(M,P)
      ENDDO
      ENDDO

      ! CALL XABORT('SNFBH3: testing 17 ')
      ISFIX=.FALSE.
      DO WHILE (.NOT.ALL(ISFIX)) ! LOOP FOR ADAPTIVE CALCULATION

      ! FLUX MOMENT COEFFICIENTS MATRIX
      Q2(:NM,:NM+1)=0.0D0
      DO IZ=1,IELEM
      DO JZ=1,IELEM
      DO IY=1,IELEM
      DO JY=1,IELEM
      DO IX=1,IELEM
      DO JX=1,IELEM

        II=IELEM**2*(IZ-1)+IELEM*(IY-1)+IX
        JJ=IELEM**2*(JZ-1)+IELEM*(JY-1)+JX

        ! IF(IPQD.EQ.3) CALL XABORT('SNFBH3: testing 69 ')
      ! CALL XABORT('SNFBH3: testing 17 ')
        ! DIAGONAL TERMS
        IF(II.EQ.JJ) THEN
          Q2(II,JJ)=SIGMA*V
     1              +CST(IX)**2*WX(JX+1)*ABS(MUH)
     2              +CST(IY)**2*WY(JY+1)*ABS(ETAH)
     3              +CST(IZ)**2*WZ(JZ+1)*ABS(XIH)

        ! IF(IND.EQ.12) CALL XABORT('SNFBH3: testing 70 ')
      ! CALL XABORT('SNFBH3: testing 191 ')
        ! UPPER DIAGONAL TERMS
        ELSEIF(II.LT.JJ) THEN
          ! CALL XABORT('SNFBH3: testing 1919 ')
          IF(IZ.EQ.JZ) THEN
            IF(IY.EQ.JY) THEN
              ! X-SPACE TERMS
              IF(MOD(IX+JX,2).EQ.1) THEN
                Q2(II,JJ)=CST(IX)*CST(JX)*WX(JX+1)*MUH
              ELSE
                Q2(II,JJ)=CST(IX)*CST(JX)*WX(JX+1)*ABS(MUH)
              ENDIF
            ELSEIF(IX.EQ.JX) THEN
              ! Y-SPACE TERMS
              IF(MOD(IY+JY,2).EQ.1) THEN
                Q2(II,JJ)=CST(IY)*CST(JY)*WY(JY+1)*ETAH
              ELSE
                Q2(II,JJ)=CST(IY)*CST(JY)*WY(JY+1)*ABS(ETAH)
              ENDIF
            ENDIF
          ELSEIF(IY.EQ.JY.AND.IX.EQ.JX) THEN
            ! Z-SPACE TERMS
            IF(MOD(IZ+JZ,2).EQ.1) THEN
              Q2(II,JJ)=CST(IZ)*CST(JZ)*WZ(JZ+1)*XIH
            ELSE
              Q2(II,JJ)=CST(IZ)*CST(JZ)*WZ(JZ+1)*ABS(XIH)
            ENDIF
          ENDIF
        ! IF(IND.EQ.12) CALL XABORT('SNFBH3: testing 75 ')
      ! CALL XABORT('SNFBH3: testing 19 ')

         ! UNDER DIAGONAL TERMS
        ELSE
          ! CALL XABORT('SNFBH3: testing 1920 ')
          IF(IZ.EQ.JZ) THEN
            IF(IY.EQ.JY) THEN
              ! X-SPACE TERMS
              IF(MOD(IX+JX,2).EQ.1) THEN
                Q2(II,JJ)=CST(IX)*CST(JX)*(WX(JX+1)-2)*MUH
              ELSE
                Q2(II,JJ)=CST(IX)*CST(JX)*WX(JX+1)*ABS(MUH)
              ENDIF
            ELSEIF(IX.EQ.JX) THEN
              ! Y-SPACE TERMS
              IF(MOD(IY+JY,2).EQ.1) THEN
                Q2(II,JJ)=CST(IY)*CST(JY)*(WY(JY+1)-2)*ETAH
              ELSE
                Q2(II,JJ)=CST(IY)*CST(JY)*WY(JY+1)*ABS(ETAH)
              ENDIF
            ENDIF
          ELSEIF(IY.EQ.JY.AND.IX.EQ.JX) THEN
            ! Z-SPACE TERMS
            IF(MOD(IZ+JZ,2).EQ.1) THEN
              Q2(II,JJ)=CST(IZ)*CST(JZ)*(WZ(JZ+1)-2)*XIH
            ELSE
              Q2(II,JJ)=CST(IZ)*CST(JZ)*WZ(JZ+1)*ABS(XIH)
            ENDIF
          ENDIF
        ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      ! FLUX SOURCE VECTOR
      DO IZ=1,IELEM
      DO IY=1,IELEM
      DO IX=1,IELEM
        II=IELEM**2*(IZ-1)+IELEM*(IY-1)+IX
        IIX=IELEM*(IZ-1)+IY
        IIY=IELEM*(IZ-1)+IX
        IIZ=IELEM*(IY-1)+IX
        Q2(II,NM+1)=Q(II)*V
        ! X-SPACE TERMS
        IF(MOD(IX,2).EQ.1) THEN
          Q2(II,NM+1)=Q2(II,NM+1)+CST(IX)*(1-WX(1))
     1                *XNI(IIX)*ABS(MUH)
        ELSE
          Q2(II,NM+1)=Q2(II,NM+1)-CST(IX)*(1+WX(1))
     1                *XNI(IIX)*MUH
        ENDIF
        ! Y-SPACE TERMS
        IF(MOD(IY,2).EQ.1) THEN
          Q2(II,NM+1)=Q2(II,NM+1)+CST(IY)*(1-WY(1))
     1                *XNJ(IIY)*ABS(ETAH)
        ELSE
          Q2(II,NM+1)=Q2(II,NM+1)-CST(IY)*(1+WY(1))
     1                *XNJ(IIY)*ETAH
        ENDIF
        ! Z-SPACE TERMS
        IF(MOD(IZ,2).EQ.1) THEN
          Q2(II,NM+1)=Q2(II,NM+1)+CST(IZ)*(1-WZ(1))
     1                *XNK(IIZ)*ABS(XIH)
        ELSE
          Q2(II,NM+1)=Q2(II,NM+1)-CST(IZ)*(1+WZ(1))
     1                *XNK(IIZ)*XIH
        ENDIF
      ENDDO
      ENDDO
      ENDDO

      CALL ALSBD(NM,1,Q2,IER,NM)
      IF(IER.NE.0) CALL XABORT('SNFBH3: SINGULAR MATRIX.')

      ! ADAPTIVE CORRECTION OF WEIGHTING PARAMETERS
      IF(ANY(ISADPT)) THEN
        IF(ISADPT(1)) THEN
          CALL SNADPT(IELEM,NM,IELEM**2,Q2(1:IELEM:1,NM+1),
     1    XNI(:NMX),1.0,WX,ISFIX(1))
        ELSE
          ISFIX(1)=.TRUE.
        ENDIF
        IF(ISADPT(2)) THEN
          CALL SNADPT(IELEM,NM,IELEM**2,Q2(1:IELEM**2:IELEM,NM+1),
     1    XNJ(:NMY),1.0,WY,ISFIX(2))
        ELSE
          ISFIX(2)=.TRUE.
        ENDIF
        IF(ISADPT(3)) THEN
          CALL SNADPT(IELEM,NM,IELEM**2,Q2(1:NM:IELEM**2,NM+1),
     1    XNK(:NMZ),1.0,WZ,ISFIX(3))
        ELSE
          ISFIX(3)=.TRUE.
        ENDIF
      ELSE
        ISFIX=.TRUE.
      ENDIF

      END DO ! END OF ADAPTIVE LOOP

      ! CLOSURE RELATIONS
      IF(IELEM.EQ.1.AND.LFIXUP.AND.(Q2(1,2).LE.RLOG)) Q2(1,2)=0.0
      ! Read XNI/XNI/XNK into TMPXNI/J/D/K
      IF((ILOZ.EQ.1).OR.(ILOZ.EQ.2))THEN
        TMPXNI(:NMX,J2,J)=WX(1)*XNI(:NMX)
      ELSE
        TMPXND(:NMX,J2,DCOORD)=WX(1)*XNI(:NMX)
      ENDIF
      IF((ILOZ.EQ.2).OR.(ILOZ.EQ.3))THEN
        TMPXNJ(:NMY,I2,I)=WY(1)*XNJ(:NMY)
      ELSE
        I3=I_FETCH
        TMPXND(:NMY,I3,DCOORD)=WY(1)*XNJ(:NMY)
      ENDIF
      TMPXNK(:NMZ,I2,J2,ILOZ,IHEX_XY)=WZ(1)*XNK(:NMZ)
      DO IZ=1,IELEM
      DO IY=1,IELEM
      DO IX=1,IELEM
        II=IELEM**2*(IZ-1)+IELEM*(IY-1)+IX
        IIX=IELEM*(IZ-1)+IY
        IIY=IELEM*(IZ-1)+IX
        IIZ=IELEM*(IY-1)+IX
        ! X-SPACE
        ! Assign I-boundary fluxes if lozenges A or B
        IF((ILOZ.EQ.1).OR.(ILOZ.EQ.2))THEN
        IF(MOD(IX,2).EQ.1) THEN
          TMPXNI(IIX,J2,J)=TMPXNI(IIX,J2,J)+CST(IX)*WX(IX+1)
     1                 *Q2(II,NM+1)
        ELSE
          TMPXNI(IIX,J2,J)=TMPXNI(IIX,J2,J)+CST(IX)*WX(IX+1)
     1                 *Q2(II,NM+1)*SIGN(1.0,MUH)
        ENDIF
        ENDIF
        ! Y-SPACE
        ! Assign J-boundary fluxes if lozenges B or C
        IF((ILOZ.EQ.2).OR.(ILOZ.EQ.3))THEN
        IF(MOD(IY,2).EQ.1) THEN
          TMPXNJ(IIY,I2,I)=TMPXNJ(IIY,I2,I)+CST(IY)*WY(IY+1)
     1                 *Q2(II,NM+1)
        ELSE
          TMPXNJ(IIY,I2,I)=TMPXNJ(IIY,I2,I)+CST(IY)*WY(IY+1)
     1                 *Q2(II,NM+1)*SIGN(1.0,ETAH)
        ENDIF
        ENDIF
        ! D-SPACE
        ! Assign D-boundary fluxes if lozenge A using XNJ
        IF((ILOZ.EQ.1))THEN
        I3=I_FETCH
        IF(MOD(IY,2).EQ.1) THEN
          TMPXND(IIY,I3,DCOORD)=TMPXND(IIY,I3,DCOORD)+CST(IY)*WY(IY+1)
     1                 *Q2(II,NM+1)
        ELSE
          TMPXND(IIY,I3,DCOORD)=TMPXND(IIY,I3,DCOORD)+CST(IY)*WY(IY+1)
     1                 *Q2(II,NM+1)*SIGN(1.0,ETAH)
        ENDIF
        ENDIF
        ! Assign D-boundary fluxes if lozenge C using XNI
        IF((ILOZ.EQ.3))THEN
        IF(MOD(IX,2).EQ.1) THEN
          TMPXND(IIX,J2,DCOORD)=TMPXND(IIX,J2,DCOORD)+CST(IX)*WX(IX+1)
     1                 *Q2(II,NM+1)
        ELSE
          TMPXND(IIX,J2,DCOORD)=TMPXND(IIX,J2,DCOORD)+CST(IX)*WX(IX+1)
     1                 *Q2(II,NM+1)*SIGN(1.0,MUH)
        ENDIF
        ENDIF
        ! Z-SPACE
        IF(MOD(IZ,2).EQ.1) THEN
          TMPXNK(IIZ,I2,J2,ILOZ,IHEX_XY)=TMPXNK(IIZ,I2,J2,ILOZ,IHEX_XY)
     1                 +CST(IZ)*WZ(IZ+1)*Q2(II,NM+1)
        ELSE
          TMPXNK(IIZ,I2,J2,ILOZ,IHEX_XY)=TMPXNK(IIZ,I2,J2,ILOZ,IHEX_XY)
     1                 +CST(IZ)*WZ(IZ+1)*Q2(II,NM+1)*SIGN(1.0,XIH)
        ENDIF
      ENDDO
      ENDDO
      ENDDO
      ! FLIP GRADIENTS IF NECESSARY
      DO IZ=1,IELEM
      DO IY=1,IELEM
        IIX=IELEM*(IZ-1)+IY
        IF((MOD(IY,2).EQ.0).AND.(ILOZ.EQ.3).AND.(IL.EQ.ISPLH))
     1    TMPXND(IIX,J2,DCOORD)=TMPXND(IIX,J2,DCOORD)*(-1)
      ENDDO
      ENDDO
      I3=I_FETCH
      DO IZ=1,IELEM
      DO IX=1,IELEM
        IIY=IELEM*(IZ-1)+IX
        IF((MOD(IX,2).EQ.0).AND.(ILOZ.EQ.1).AND.(JL.EQ.ISPLH))
     1    TMPXND(IIY,I3,DCOORD)=TMPXND(IIY,I3,DCOORD)*(-1)
      ENDDO
      ENDDO
      ! LFIXUP
      IF(IELEM.EQ.1.AND.LFIXUP)THEN
        IF((ILOZ.EQ.1).OR.(ILOZ.EQ.2))THEN
          IF(TMPXNI(1,J2,J).LE.RLOG) TMPXNI(1,J2,J)=0.0
        ELSE
          IF(TMPXND(1,J2,DCOORD).LE.RLOG) TMPXND(1,J2,DCOORD)=0.0
        ENDIF
        IF((ILOZ.EQ.2).OR.(ILOZ.EQ.3))THEN
          IF(TMPXNJ(1,I2,I).LE.RLOG) TMPXNJ(1,I2,I)=0.0
        ELSE
          I3=I_FETCH
          IF(TMPXND(1,I3,DCOORD).LE.RLOG) TMPXND(1,I3,DCOORD)=0.0
        ENDIF
        IF(TMPXNK(1,I2,J2,ILOZ,IHEX_XY).LE.RLOG) 
     1   TMPXNK(1,I2,J2,ILOZ,IHEX_XY)=0.0
      ENDIF
      WX=WX0
      WY=WY0
      WZ=WZ0

      ! SAVE LEGENDRE MOMENT OF THE FLUX
      IOF=((ILOZ-1)*ISPLH+(J2-1))*ISPLH+I2
      DO P=1,NSCT
      DO IEL=1,NM
      FLUX(IEL,P,IOF,IHEX_XY,K) = FLUX(IEL,P,IOF,IHEX_XY,K)
     1   +Q2(IEL,NM+1)*DN(P,M)
      ENDDO
      ENDDO

      ENDDO ! END OF WITHIN LOZENGE J-LOOP
      ENDDO ! END OF WITHIN LOZENGE I-LOOP

      ENDDO ! END OF LOZENGE LOOP

      ENDDO ! END OF I COLUMNS OF CARTESIAN MAP LOOP
      ENDDO ! END OF J COLUMNS OF CARTESIAN MAP LOOP
      ENDDO ! END OF Z-LOOP

      ! SAVE FLUX INFORMATION
      FLUX_G(:,:,:,:,:,IG)=FLUX_G(:,:,:,:,:,IG)+FLUX(:,:,:,:,:)

      ! SAVE K-BOUNDARY CONDITIONS IF NOT VOID B.C.
      IF((NCODE(5).NE.1).or.(NCODE(6).NE.1))THEN
        IOF=(M-1)*(L5)
        FUNKNO(LFLX+IOF+1:LFLX+IOF+L5,IG)=
     >  RESHAPE(REAL(TMPXNK(:NMZ,:ISPLH,:ISPLH,:3,:NHEX)),
     >  (/L5/))
      ENDIF

      ENDDO ! END OF DIRECTION LOOP
      ENDDO ! END OF ENERGY LOOP
*$OMP END PARALLEL DO
      ENDDO ! END OF OCTANT LOOP

      ! SAVE FLUX INFORMATION
      DO IG=1,NGEFF
        IF(.NOT.INCONV(IG)) CYCLE
        FUNKNO(:LFLX,IG)=
     1  RESHAPE(REAL(FLUX_G(:NM,:NSCT,:3*ISPLH**2,:NHEX,:LZ,IG)),
     2  (/LFLX/))
      ENDDO
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(FLUX_G,FLUX,INDANG,TMPXNI,TMPXNJ,TMPXND,TMPXNK,TMPMAT)
      RETURN
  500 FORMAT(16H SNFBH3: thread=,I8,12H --->(group=,I4,7H angle=,I4,1H))
      END
