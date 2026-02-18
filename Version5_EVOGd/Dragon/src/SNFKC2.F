*DECK SNFKC2
      SUBROUTINE SNFKC2(NKBA,NUN,NGEFF,IMPX,INCONV,NGIND,LX,LY,IELEM,
     1 NM,NMX,NMY,NMAT,NPQ,NSCT,MAT,VOL,TOTAL,NCODE,ZCODE,QEXT,LFIXUP,
     2 DU,DE,W,MRM,MRMY,DB,DA,MN,DN,WX,WY,CST,ISADPT,ISBS,NBS,ISBSM,BS,
     3 MAXL,FUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Perform one inner iteration for solving SN equations in 2D Cartesian
* geometry for the HODD method. KBA-like multithreading. Albedo
* boundary conditions. Boltzmann (BTE) discretization.
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
* NKBA    number of macrocells along each axis
* NUN     total number of unknowns in vector FUNKNO.
* NGEFF   number of energy groups processed in parallel.
* IMPX    print flag (equal to zero for no print).
* INCONV  energy group convergence flag (set to .FALSE. if converged).
* NGIND   energy group indices assign to the NGEFF set.
* LX      number of meshes along X axis.
* LY      number of meshes along Y axis.
* IELEM   measure of order of the spatial approximation polynomial:
*         =1 constant - default for HODD;
*         =2 linear - default for DG;
*         >3 higher orders.
* NM      number of moments in space and energy for flux components
* NMX     number of moments for X axis boundaries components
* NMY     number of moments for Y axis boundaries components
* NMAT    number of material mixtures.
* NPQ     number of SN directions in four octants (including zero-weight
*         directions).
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
* W       weights.
* MRM     quadrature index.
* MRMY    quadrature index.
* DB      diamond-scheme parameter.
* DA      diamond-scheme parameter.
* MN      moment-to-discrete matrix.
* DN      discrete-to-moment matrix.
* ISBS    flag to indicate the presence or not of boundary fixed
*         sources.
* NBS     number of boundary fixed sources.
* ISBSM   flag array to indicate the presence or not of boundary fixed
*         source in each unit surface.
* BS      boundary source array with their intensities.
* MAXL    maximum size of boundary source array.
* WX      spatial X axis closure relation weighting factors.
* WY      spatial Y axis closure relation weighting factors.
* CST     constants for the polynomial approximations.
* ISADPT  flag to enable/disable adaptive flux calculations.
*
*Parameters: input/output
* FUNKNO  Legendre components of the flux and boundary fluxes.
* FLUXC   flux at the cutoff energy.
*
*-----------------------------------------------------------------------
#if defined(_OPENMP)
      USE omp_lib
#endif
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NKBA,NUN,NGEFF,IMPX,NGIND(NGEFF),LX,LY,IELEM,NM,NMX,NMY,
     1 NMAT,NPQ,NSCT,MAT(LX,LY),NCODE(4),MRM(NPQ),MRMY(NPQ),ISBS,NBS,
     2 ISBSM(4*ISBS,NPQ*ISBS,NGEFF*ISBS),MAXL
      LOGICAL INCONV(NGEFF)
      REAL VOL(LX,LY),TOTAL(0:NMAT,NGEFF),ZCODE(4),QEXT(NUN,NGEFF),
     1 DU(NPQ),DE(NPQ),W(NPQ),DB(LX,NPQ),DA(LX,LY,NPQ),
     2 FUNKNO(NUN,NGEFF),BS(MAXL*ISBS,NBS*ISBS),WX(IELEM+1),
     3 WY(IELEM+1),CST(IELEM),MN(NPQ,NSCT),DN(NSCT,NPQ)
      LOGICAL LFIXUP,ISADPT(2)
*----
*  LOCAL VARIABLES
*----
      INTEGER NPQD(4),IIND(4),P
      REAL WX0(IELEM+1),WY0(IELEM+1)
      DOUBLE PRECISION Q(NM),Q2(NM,NM+1),V
      PARAMETER(IUNOUT=6,RLOG=1.0E-8,PI=3.141592654)
      LOGICAL ISFIX(2)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: INDANG
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: FLUX
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: FLUX_G
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: XNI,XNJ
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(INDANG(NPQ,4))
      ALLOCATE(FLUX(NM,NSCT))
      ALLOCATE(FLUX_G(NM,NSCT,LX,LY,NGEFF))
      ALLOCATE(XNI(IELEM,LY,NPQ,NGEFF),XNJ(IELEM,LX,NPQ,NGEFF))
*----
*  LENGTH OF FUNKNO COMPONENTS (IN ORDER)
*----
      LFLX=NM*LX*LY*NSCT
      LXNI=NMX*LY*NPQ
      LXNJ=NMY*LX*NPQ
*----
*  NUMBER OF MACROCELLS (MACRO*)
*  NUMBER OF LZ LAYERS IN EACH MACROCELL (NCELL*)
*----
      MACROX=NKBA
      MACROY=NKBA
      NCELLX=1+(LX-1)/MACROX
      NCELLY=1+(LY-1)/MACROY
*----
*  SET OCTANT SWAPPING ORDER.
*----
      NPQD(:4)=0
      INDANG(:NPQ,:4)=0
      IIND(:)=0
      DO M=1,NPQ
        VU=DU(M)
        VE=DE(M)
        IF(W(M).EQ.0) CYCLE
        IF((VU.GE.0.0).AND.(VE.GE.0.0)) THEN
          IND=1
          JND=4
        ELSE IF((VU.LE.0.0).AND.(VE.GE.0.0)) THEN
          IND=2
          JND=3
        ELSE IF((VU.LE.0.0).AND.(VE.LE.0.0)) THEN
          IND=3
          JND=1
        ELSE
          IND=4
          JND=2
        ENDIF
        IIND(JND)=IND
        NPQD(IND)=NPQD(IND)+1
        INDANG(NPQD(IND),IND)=M
      ENDDO
*----
*  MAIN LOOP OVER OCTANTS.
*----

      FLUX_G(:NM,:NSCT,:LX,:LY,:NGEFF)=0.0D0
      WX0=WX
      WY0=WY

      DO 240 JND=1,4
      IND=IIND(JND)
*----
*  PRELIMINARY LOOPS FOR SETTING BOUNDARY CONDITIONS.
*----

*$OMP  PARALLEL DO
*$OMP+ PRIVATE(M,IG,VU,VE,M1,IOF,JOF,IEL,I,J,IPQD)
*$OMP+ SHARED(FUNKNO) COLLAPSE(2)

      DO 70 IG=1,NGEFF
      DO 60 IPQD=1,NPQD(IND)
      IF(.NOT.INCONV(IG)) GO TO 60
      M=INDANG(IPQD,IND)
      VU=DU(M)
      VE=DE(M)
      ! X-BOUNDARY
      IF(VU.GT.0.0)THEN
         M1=MRM(M)
         IF((NCODE(1).NE.4))THEN
            DO IEL=1,NMX
               DO J=1,LY
                  IOF=((M-1)*LY+(J-1))*NMX+IEL
                  JOF=((M1-1)*LY+(J-1))*NMX+IEL
                  FUNKNO(LFLX+IOF,IG)=FUNKNO(LFLX+JOF,IG)
               ENDDO
            ENDDO
         ENDIF
      ELSEIF(VU.LT.0.0)THEN
         M1=MRM(M)
         IF((NCODE(2).NE.4))THEN
            DO IEL=1,NMX
               DO J=1,LY
                  IOF=((M-1)*LY+(J-1))*NMX+IEL
                  JOF=((M1-1)*LY+(J-1))*NMX+IEL
                  FUNKNO(LFLX+IOF,IG)=FUNKNO(LFLX+JOF,IG)
               ENDDO
            ENDDO
         ENDIF
      ENDIF
      ! Y-BOUNDARY
      IF(VE.GT.0.0)THEN
         M1=MRMY(M)
         IF((NCODE(3).NE.4))THEN
            DO IEL=1,NMY
               DO I=1,LX
                  IOF=((M-1)*LX+(I-1))*NMY+IEL
                  JOF=((M1-1)*LX+(I-1))*NMY+IEL
                  FUNKNO(LFLX+LXNI+IOF,IG)=
     >               FUNKNO(LFLX+LXNI+JOF,IG)
               ENDDO
            ENDDO
         ENDIF
      ELSEIF(VE.LT.0.0)THEN
         M1=MRMY(M)
         IF((NCODE(4).NE.4))THEN
            DO IEL=1,NMY
               DO I=1,LX
                  IOF=((M-1)*LX+(I-1))*NMY+IEL
                  JOF=((M1-1)*LX+(I-1))*NMY+IEL
                  FUNKNO(LFLX+LXNI+IOF,IG)=
     >               FUNKNO(LFLX+LXNI+JOF,IG)
               ENDDO
            ENDDO
         ENDIF
      ENDIF
   60 CONTINUE
   70 CONTINUE

*OMP END PARALLEL DO

*----
*  KBA DIAGONAL LOOP
*----
      XNI(:IELEM,:LY,:NPQ,:NGEFF)=0.0D0
      XNJ(:IELEM,:LX,:NPQ,:NGEFF)=0.0D0
      DO 230 IDI=1,MACROX+MACROY-1

*----
*  MAIN SWAPPING LOOPS FOR SN FLUX CALCULATION
*----

*$OMP  PARALLEL DO 
*$OMP+ PRIVATE(ITID,FLUX,M,IG,Q,Q2,IOF,IER,II,JJ,IEL,I,J,P)
*$OMP+ PRIVATE(IPQD,IMX,IMY,IBM,SIGMA,V,ISFIX,IX,JX,IY,JY,JCEL)
*$OMP+ FIRSTPRIVATE(WX,WY,WX0,WY0) SHARED(FUNKNO,XNI,XNJ)
*$OMP+ REDUCTION(+:FLUX_G) COLLAPSE(3)

      ! LOOP FOR GROUPS TO EXECUTE IN PARALLEL
      DO 220 IG=1,NGEFF

      ! LOOP OVER ALL DIRECTIONS
      DO 210 IPQD=1,NPQD(IND)

      ! LOOP OVER MACROCELLS IN WAVEFRONT
      DO 200 ICEL=MAX(1,IDI-MACROY+1),MIN(MACROX,IDI)

      IF(.NOT.INCONV(IG)) GO TO 200
      M=INDANG(IPQD,IND)
      IF(W(M).EQ.0.0) GO TO 200

      ! GET AND PRINT THREAD NUMBER
#if defined(_OPENMP)
      ITID=omp_get_thread_num()
#else
      ITID=0
#endif
      IF(IMPX.GT.5) WRITE(IUNOUT,400) ITID,NGIND(IG),IPQD

*----
*  LOOP OVER X- AND Y-DIRECTED AXES.
*----
      JCEL=IDI-ICEL+1

      ! X-BOUNDARIES CONDITIONS
      IF(ICEL.EQ.1) THEN
        DO IMY=1,MIN(NCELLY,LY-(JCEL-1)*NCELLY)
        J=(JCEL-1)*NCELLY+IMY
        IF((IND.EQ.3).OR.(IND.EQ.4)) J=LY+1-J
        DO IEL=1,NMX
        IOF=(M-1)*NMX*LY+(J-1)*NMX+IEL
        IF((IND.EQ.1).OR.(IND.EQ.4)) THEN
          XNI(IEL,J,IPQD,IG)=FUNKNO(LFLX+IOF,IG)*ZCODE(1)
        ELSE
          XNI(IEL,J,IPQD,IG)=FUNKNO(LFLX+IOF,IG)*ZCODE(2)
        ENDIF
        ENDDO

        ! X-BOUNDARIES FIXED SOURCES
        IF(ISBS.EQ.1) THEN
          IF((IND.EQ.2.OR.IND.EQ.3).AND.ISBSM(2,M,IG).NE.0) THEN
            XNI(1,J,IPQD,IG)=XNI(1,J,IPQD,IG)+BS(J,ISBSM(2,M,IG))
          ELSE IF((IND.EQ.1.OR.IND.EQ.4).AND.ISBSM(1,M,IG).NE.0) THEN
            XNI(1,J,IPQD,IG)=XNI(1,J,IPQD,IG)+BS(J,ISBSM(1,M,IG))
          ENDIF
        ENDIF 
        ENDDO
      ENDIF

      ! Y-BOUNDARIES CONDITIONS
      IF(JCEL.EQ.1) THEN
        DO IMX=1,MIN(NCELLX,LX-(ICEL-1)*NCELLX)
        I=(ICEL-1)*NCELLX+IMX
        IF((IND.EQ.2).OR.(IND.EQ.3)) I=LX+1-I
        DO IEL=1,NMY
        IOF=(M-1)*NMY*LX+(I-1)*NMY+IEL
        IF((IND.EQ.1).OR.(IND.EQ.2)) THEN
          XNJ(IEL,I,IPQD,IG)=FUNKNO(LFLX+LXNI+IOF,IG)*ZCODE(3)
        ELSE
          XNJ(IEL,I,IPQD,IG)=FUNKNO(LFLX+LXNI+IOF,IG)*ZCODE(4)
        ENDIF
        ENDDO

        ! Y-BOUNDARIES FIXED SOURCES
        IF(ISBS.EQ.1) THEN
          IF((IND.EQ.3.OR.IND.EQ.4).AND.ISBSM(4,M,IG).NE.0) THEN
            XNJ(1,I,IPQD,IG)=XNJ(1,I,IPQD,IG)+BS(I,ISBSM(4,M,IG))
          ELSE IF((IND.EQ.1.OR.IND.EQ.2).AND.ISBSM(3,M,IG).NE.0) THEN
            XNJ(1,I,IPQD,IG)=XNJ(1,I,IPQD,IG)+BS(I,ISBSM(3,M,IG))
          ENDIF
        ENDIF
        ENDDO
      ENDIF

      ! X-AXIS LOOP
      DO 190 IMX=1,MIN(NCELLX,LX-(ICEL-1)*NCELLX)
      I=(ICEL-1)*NCELLX+IMX
      IF((IND.EQ.2).OR.(IND.EQ.3)) I=LX+1-I

      ! Y-AXIS LOOP
      DO 180 IMY=1,MIN(NCELLY,LY-(JCEL-1)*NCELLY)
      J=(JCEL-1)*NCELLY+IMY
      IF((IND.EQ.3).OR.(IND.EQ.4)) J=LY+1-J

      ! INITIALIZE FLUXES
      FLUX(:NM,:NSCT)=0.0D0

      ! IF(ITID.EQ.1)THEN
      !    WRITE(*,*) 'I,J', JCEL, NCELLY, IMY
      !    WRITE(*,*) IND, IPQD, I, J
      !    WRITE(*,*) XNI(:,J,IPQD,IG)
      !    WRITE(*,*) XNJ(:,I,IPQD,IG)
      ! ENDIF

      ! DATA
      IBM=MAT(I,J)
      IF(IBM.EQ.0) GO TO 200
      SIGMA=TOTAL(IBM,IG)
      V=VOL(I,J)

      ! SOURCE DENSITY TERM
      DO IEL=1,NM
      Q(IEL)=0.0D0
      DO P=1,NSCT
      IOF=((J-1)*LX*NSCT+(I-1)*NSCT+(P-1))*NM+IEL
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
     1              +CST(IX)**2*WX(JX+1)*ABS(DA(I,J,M))
     2              +CST(IY)**2*WY(JY+1)*ABS(DB(I,M))

        ! UPPER DIAGONAL TERMS
        ELSEIF(II.LT.JJ) THEN
          ! X-SPACE TERMS
          IF(IY.EQ.JY) THEN
          IF(MOD(IX+JX,2).EQ.1) THEN
            Q2(II,JJ)=CST(IX)*CST(JX)*WX(JX+1)*DA(I,J,M)
          ELSE
            Q2(II,JJ)=CST(IX)*CST(JX)*WX(JX+1)*ABS(DA(I,J,M))
          ENDIF
          ! Y-SPACE TERMS
          ELSEIF(IX.EQ.JX) THEN
          IF(MOD(IY+JY,2).EQ.1) THEN
            Q2(II,JJ)=CST(IY)*CST(JY)*WY(JY+1)*DB(I,M)
          ELSE
            Q2(II,JJ)=CST(IY)*CST(JY)*WY(JY+1)*ABS(DB(I,M))
          ENDIF
          ENDIF

        ! UNDER DIAGONAL TERMS
        ELSE
          ! X-SPACE TERMS
          IF(IY.EQ.JY) THEN
          IF(MOD(IX+JX,2).EQ.1) THEN
            Q2(II,JJ)=CST(IX)*CST(JX)*(WX(JX+1)-2)*DA(I,J,M)
          ELSE
            Q2(II,JJ)=CST(IX)*CST(JX)*WX(JX+1)*ABS(DA(I,J,M))
          ENDIF 
          ! Y-SPACE TERMS
          ELSEIF(IX.EQ.JX) THEN
          IF(MOD(IY+JY,2).EQ.1) THEN
            Q2(II,JJ)=CST(IY)*CST(JY)*(WY(JY+1)-2)*DB(I,M)
          ELSE
            Q2(II,JJ)=CST(IY)*CST(JY)*WY(JY+1)*ABS(DB(I,M))
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
     1                *XNI(IY,J,IPQD,IG)*ABS(DA(I,J,M))
        ELSE
          Q2(II,NM+1)=Q2(II,NM+1)-CST(IX)*(1+WX(1))
     1                *XNI(IY,J,IPQD,IG)*DA(I,J,M)
        ENDIF
        ! Y-SPACE TERMS
        IF(MOD(IY,2).EQ.1) THEN
          Q2(II,NM+1)=Q2(II,NM+1)+CST(IY)*(1-WY(1))
     1                *XNJ(IX,I,IPQD,IG)*ABS(DB(I,M))
        ELSE
          Q2(II,NM+1)=Q2(II,NM+1)-CST(IY)*(1+WY(1))
     1                *XNJ(IX,I,IPQD,IG)*DB(I,M)
        ENDIF
      ENDDO
      ENDDO

      CALL ALSBD(NM,1,Q2,IER,NM)
      IF(IER.NE.0) CALL XABORT('SNFKC2: SINGULAR MATRIX.')

      ! IF(ITID.EQ.1)THEN
      !    WRITE(*,*) 'I,J'
      !    WRITE(*,*) Q2(:,NM+1)
      !    WRITE(*,*) ' '
      !    WRITE(*,*) ' '
      !    WRITE(*,*) ' '
      !    WRITE(*,*) ' '
      ! ENDIF

      ! ADAPTIVE CORRECTION OF WEIGHTING PARAMETERS
      IF(ANY(ISADPT)) THEN
        IF(ISADPT(1)) THEN
          CALL SNADPT(IELEM,NM,IELEM,Q2(1:IELEM:1,NM+1),
     1    XNI(:NMX,J,IPQD,IG),1.0,WX,ISFIX(1))
        ELSE
          ISFIX(1)=.TRUE.
        ENDIF
        IF(ISADPT(2)) THEN
          CALL SNADPT(IELEM,NM,IELEM,Q2(1:NM:IELEM,NM+1),
     1    XNJ(:NMY,I,IPQD,IG),1.0,WY,ISFIX(2))
        ELSE
          ISFIX(2)=.TRUE.
        ENDIF
      ELSE
        ISFIX=.TRUE.
      ENDIF

      END DO ! END OF ADAPTIVE LOOP

      ! CLOSURE RELATIONS
      IF(IELEM.EQ.1.AND.LFIXUP.AND.(Q2(1,2).LE.RLOG)) Q2(1,2)=0.0
      XNI(:NMX,J,IPQD,IG)=WX(1)*XNI(:NMX,J,IPQD,IG)
      XNJ(:NMY,I,IPQD,IG)=WY(1)*XNJ(:NMY,I,IPQD,IG)
      DO IY=1,IELEM
      DO IX=1,IELEM
        II=IELEM*(IY-1)+IX
        ! X-SPACE
        IF(MOD(IX,2).EQ.1) THEN
          XNI(IY,J,IPQD,IG)=XNI(IY,J,IPQD,IG)+CST(IX)*WX(IX+1)
     1                 *Q2(II,NM+1)
        ELSE
          XNI(IY,J,IPQD,IG)=XNI(IY,J,IPQD,IG)+CST(IX)*WX(IX+1)
     1                 *Q2(II,NM+1)*SIGN(1.0,DA(I,J,M))
        ENDIF
        ! Y-SPACE
        IF(MOD(IY,2).EQ.1) THEN
          XNJ(IX,I,IPQD,IG)=XNJ(IX,I,IPQD,IG)+CST(IY)*WY(IY+1)
     1                 *Q2(II,NM+1)
        ELSE
          XNJ(IX,I,IPQD,IG)=XNJ(IX,I,IPQD,IG)+CST(IY)*WY(IY+1)
     1                 *Q2(II,NM+1)*SIGN(1.0,DB(I,M))
        ENDIF
      ENDDO
      ENDDO
      IF(IELEM.EQ.1.AND.LFIXUP.AND.(XNI(1,J,IPQD,IG).LE.RLOG)) 
     1   XNI(1,J,IPQD,IG)=0.0
      IF(IELEM.EQ.1.AND.LFIXUP.AND.(XNJ(1,I,IPQD,IG).LE.RLOG)) 
     1   XNJ(1,I,IPQD,IG)=0.0
      WX=WX0
      WY=WY0

      ! SAVE LEGENDRE MOMENT OF THE FLUX
      DO P=1,NSCT
      DO IEL=1,NM
      FLUX(IEL,P)=FLUX(IEL,P)+Q2(IEL,NM+1)*DN(P,M)
      ENDDO
      ENDDO

      ! SAVE X-BOUNDARY CONDITIONS
      IF((ICEL-1)*NCELLX+IMX.EQ.LX) THEN
        DO IEL=1,NMX
        IOF=(M-1)*NMX*LY+(J-1)*NMX+IEL
        FUNKNO(LFLX+IOF,IG)=REAL(XNI(IEL,J,IPQD,IG))
        ENDDO
      ENDIF

      ! SAVE Y-BOUNDARY CONDITIONS
      IF((JCEL-1)*NCELLY+IMY.EQ.LY) THEN
        DO IEL=1,NMY
        IOF=(M-1)*NMY*LX+(I-1)*NMY+IEL
        FUNKNO(LFLX+LXNI+IOF,IG)=REAL(XNJ(IEL,I,IPQD,IG))
        ENDDO
      ENDIF

      ! SAVE FLUX INFORMATION
      FLUX_G(:,:,I,J,IG)=FLUX_G(:,:,I,J,IG)+FLUX(:,:)

  180 CONTINUE ! END OF Y-LOOP
  190 CONTINUE ! END OF X-LOOP

  200 CONTINUE ! END OF MACROCELL LOOP
  210 CONTINUE ! END OF DIRECTION LOOP
  220 CONTINUE ! END OF ENERGY LOOP
*$OMP END PARALLEL DO
  230 CONTINUE ! END OF WAVEFRONT LOOP
  240 CONTINUE ! END OF OCTANT LOOP

      ! SAVE FLUX INFORMATION
      DO 250 IG=1,NGEFF
        IF(.NOT.INCONV(IG)) GO TO 250
        FUNKNO(:LFLX,IG)=
     1  RESHAPE(REAL(FLUX_G(:NM,:NSCT,:LX,:LY,IG)),(/ LFLX /) )
  250 CONTINUE

      ! CALL XABORT('SNFKC2: testing')

*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(XNI,XNJ,FLUX_G,FLUX,INDANG)
      RETURN
  400 FORMAT(16H SNFKC2: thread=,I8,12H --->(group=,I4,7H angle=,I4,1H))
      END
