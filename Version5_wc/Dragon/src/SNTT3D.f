*DECK SNTT3D
      SUBROUTINE SNTT3D(IGE,IMPX,LX,LY,LZ,SIDE,IELEM,NLF,NPQ,NSCT,IQUAD,
     1 NCODE,ZCODE,MAT,XXX,YYY,ZZZ,VOL,IDL,DU,DE,DZ,W,MRMX,MRMY,MRMZ,DC,
     2 DB,DA,PL,LL4,NUN,EELEM,WX,WE,CST,IBFP,ISCHM,ESCHM,IGLK,MN,DN,IL,
     3 IM,ISCAT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Numbering corresponding to a 3-D Cartesian with discrete ordinates
* approximation of the flux.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): N. Martin and C. Bienvenue
*
*Parameters: input
* IGE     geometry type with (=0 for Cartesian, =2 for hexagonal).
* IMPX    print parameter.
* LX      number of elements along the X axis;
*         OR, number of elements in X-Y plane in hexagonal geometry,
*         including lozenges and further submeshing of lozenges
*         e.g. domain of 7 hex. by 5 levels with a submeshing of 2
*         will have LX=7*3*2*2=84 
* LY      number of elements along the Y axis.
* LZ      number of elements along the Z axis.
* SIDE    side of hexagon.
* IELEM   measure of order of the spatial approximation polynomial:
*         =1 constant - only for HODD, classical diamond scheme 
*         (default for HODD);
*         =2 linear - default for DG;
*         =3 parabolic;
*         =4 cubic - only for DG.
* NLF     SN order for the flux (even number).
* NPQ     number of SN directions in eight octants.
* NSCT    maximum number of spherical harmonics moments of the flux.
* IQUAD   type of SN quadrature (1 Level symmetric, type IQUAD;
*         4 Legendre-Chebyshev; 5 symmetric Legendre-Chebyshev;
*         6 quadruple range).
* NCODE   type of boundary condition applied on each side
*         (i=1 X-;  i=2 X+;  i=3 Y-;  i=4 Y+;  i=5 Z-;  i=6 Z+):
*         =1 VOID;   =2 REFL;  =4 TRAN.
* ZCODE   ZCODE(I) is the albedo corresponding to boundary condition
*         'VOID' on each side (ZCODE(I)=0.0 by default).
* MAT     mixture index assigned to each element.
* XXX     Cartesian coordinates along the X axis.
* YYY     Cartesian coordinates along the Y axis.
* ZZZ     Cartesian coordinates along the Z axis.
* EELEM   measure of order of the energy approximation polynomial:
*         =1 constant - default for HODD;
*         =2 linear - default for DG;
*         >3 higher orders.
* IBFP    type of energy proparation relation:
*         =0 no Fokker-Planck term;
*         =1 Galerkin type;
*         =2 heuristic Przybylski and Ligou type.
* ISCHM   method of spatial discretisation:
*         =1 High-Order Diamond Differencing (HODD) - default;
*         =2 Discontinuous Galerkin finite element method (DG);
*         =3 Adaptive weighted method (AWD).
* ESCHM   method of energy discretisation:
*         =1 High-Order Diamond Differencing (HODD) - default;
*         =2 Discontinuous Galerkin finite element method (DG);
*         =3 Adaptive weighted method (AWD).
* IGLK    angular interpolation type:
*         =0 classical SN method.
*         =1 Galerkin quadrature method (M = inv(D))
*         =2 Galerkin quadrature method (D = inv(M))
* ISCAT   maximum number of spherical harmonics moments of the flux.
*
*Parameters: output
* VOL     volume of each element.
* IDL     isotropic flux indices.
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
* PL      discrete values of the spherical harmonics corresponding
*         to the 3D SN quadrature.
* LL4     number of unknowns being solved for, over the domain. This 
*         includes the various moments of the isotropic (and if present,
*         anisotropic) flux. 
* NUN     total number of unknowns stored in the FLUX vector per group.
*         This includes LL4 (see above) as well as any surface boundary
*         fluxes, if present.
* WX      spatial closure relation weighting factors.
* WE      energy closure relation weighting factors.
* CST     constants for the polynomial approximations.
* MN      moment-to-discrete matrix.
* DN      discrete-to-moment matrix.
* IL      indexes (l) of each spherical harmonics in the
*         interpolation basis.
* IM      indexes (m) of each spherical harmonics in the
*         interpolation basis.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IMPX,LX,LY,LZ,IELEM,NLF,NPQ,NSCT,IQUAD,NCODE(6),
     1 MAT(LX,LY,LZ),IDL(LX*LY*LZ),MRMX(NPQ),MRMY(NPQ),MRMZ(NPQ),EELEM,
     2 IBFP,ISCHM,ESCHM,IL(NSCT),IM(NSCT),ISCAT,IGLK
      REAL ZCODE(6),VOL(LX,LY,LZ),XXX(LX+1),YYY(LY+1),ZZZ(LZ+1),
     1 DU(NPQ),DE(NPQ),DZ(NPQ),W(NPQ),DC(LX,LY,NPQ),DB(LX,LZ,NPQ),
     2 DA(LY,LZ,NPQ),PL(NSCT,NPQ),WX(IELEM+1),
     3 WE(EELEM+1),CST(MAX(IELEM,EELEM)),MN(NPQ,NSCT),DN(NSCT,NPQ)
*----
*  LOCAL VARIABLES
*----
      CHARACTER HSMG*131
      LOGICAL L1,L2,L3,L4,L5,L6
      PARAMETER(RLOG=1.0E-8,PI=3.141592654)
      REAL PX,PE
      DOUBLE PRECISION MND(NPQ,NPQ),NORM,IPROD
      INTEGER, ALLOCATABLE, DIMENSION(:) :: JOP
      REAL, ALLOCATABLE, DIMENSION(:) :: XX,YY,ZZ,UU,WW,TPQ,UPQ,VPQ,WPQ
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: V,V2
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: U
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: RLM
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(XX(LX),YY(LY),ZZ(LZ))
*----
*  UNFOLD HEIGHT-OCTANT QUADRATURES.
*----
      IF(MOD(NLF,2).EQ.1) CALL XABORT('SNTT3D: EVEN NLF EXPECTED.')
      IF(IQUAD.EQ.10) THEN
         NPQ0=NLF**2/4
      ELSE
         NPQ0=NLF*(NLF/2+1)/4
      ENDIF
      ALLOCATE(JOP(NLF/2),UU(NLF/2),WW(NLF/2),TPQ(NPQ0),UPQ(NPQ0),
     1 VPQ(NPQ0),WPQ(NPQ0))
      IF(IQUAD.EQ.1) THEN
!Level-symmetric quadrature of type 1
         CALL SNQU01(NLF,JOP,UU,WW,TPQ,UPQ,VPQ,WPQ)
      ELSE IF(IQUAD.EQ.2) THEN
!Level-symmetric quadrature of type 2
         CALL SNQU02(NLF,JOP,UU,WW,TPQ,UPQ,VPQ,WPQ)
      ELSE IF(IQUAD.EQ.3) THEN
!Snow Level-symmetric type quadrature 
         CALL SNQU03(NLF,JOP,UU,WW,TPQ,UPQ,VPQ,WPQ)
      ELSE IF(IQUAD.EQ.4) THEN
!Legendre-Chebyshev quadrature 
         CALL SNQU04(NLF,JOP,UU,WW,TPQ,UPQ,VPQ,WPQ)
      ELSE IF(IQUAD.EQ.5) THEN
!Symmetric Legendre-Chebyshev quadrature
         CALL SNQU05(NLF,TPQ,UPQ,VPQ,WPQ)
      ELSE IF(IQUAD.EQ.6) THEN
!Quadruple Range quadrature
         CALL SNQU06(NLF,TPQ,UPQ,VPQ,WPQ)
      ELSE IF(IQUAD.EQ.10) THEN
         CALL SNQU10(NLF,JOP,UU,WW,TPQ,UPQ,VPQ,WPQ)
      ELSE
         CALL XABORT('SNTT3D: UNKNOWN QUADRATURE TYPE.')
      ENDIF
      N=0
      IOF=0
      DO 320 I=1,NLF/2
         JOF=IOF+NLF-2*I+2
         DO 330 J=0,NLF/2-I
            IOF=IOF+1
            KOF=IOF+(NLF+2)*NLF/4
            LOF=IOF+(NLF+2)*NLF/2
            MRMX(IOF)=JOF
            MRMY(IOF)=KOF
            MRMZ(IOF)=LOF
            DU(IOF)=-UPQ(N+J+1)
            DE(IOF)=-VPQ(N+J+1)
            DZ(IOF)=-TPQ(N+J+1)
            W(IOF)=WPQ(N+J+1)
            JOF=JOF-1
 330     CONTINUE
         DO 340 J=NLF/2-I,0,-1
            IOF=IOF+1
            KOF=IOF+(NLF+2)*NLF/4
            LOF=IOF+(NLF+2)*NLF/2
            MRMX(IOF)=JOF
            MRMY(IOF)=KOF
            MRMZ(IOF)=LOF
            DU(IOF)=UPQ(N+J+1)
            DE(IOF)=-VPQ(N+J+1)
            DZ(IOF)=-TPQ(N+J+1)
            W(IOF)=WPQ(N+J+1)
            JOF=JOF-1
 340     CONTINUE
         N=N+NLF/2-I+1
 320  CONTINUE
      N=0
      DO 350 I=1,NLF/2
         JOF=IOF+NLF-2*I+2
         DO 360 J=0,NLF/2-I
            IOF=IOF+1
            KOF=IOF-(NLF+2)*NLF/4
            LOF=IOF+(NLF+2)*NLF/2
            MRMX(IOF)=JOF
            MRMY(IOF)=KOF
            MRMZ(IOF)=LOF
            DU(IOF)=-UPQ(N+J+1)
            DE(IOF)=VPQ(N+J+1)
            DZ(IOF)=-TPQ(N+J+1)
            W(IOF)=WPQ(N+J+1)
            JOF=JOF-1
 360     CONTINUE
         DO 370 J=NLF/2-I,0,-1
            IOF=IOF+1
            KOF=IOF-(NLF+2)*NLF/4
            LOF=IOF+(NLF+2)*NLF/2
            MRMX(IOF)=JOF
            MRMY(IOF)=KOF
            MRMZ(IOF)=LOF
            DU(IOF)=UPQ(N+J+1)
            DE(IOF)=VPQ(N+J+1)
            DZ(IOF)=-TPQ(N+J+1)
            W(IOF)=WPQ(N+J+1)
            JOF=JOF-1
 370     CONTINUE
         N=N+NLF/2-I+1
 350  CONTINUE
      N=0
      DO 380 I=1,NLF/2
         JOF=IOF+NLF-2*I+2
         DO 390 J=0,NLF/2-I
            IOF=IOF+1
            KOF=IOF+(NLF+2)*NLF/4
            LOF=IOF-(NLF+2)*NLF/2
            MRMX(IOF)=JOF
            MRMY(IOF)=KOF
            MRMZ(IOF)=LOF
            DU(IOF)=-UPQ(N+J+1)
            DE(IOF)=-VPQ(N+J+1)
            DZ(IOF)=TPQ(N+J+1)
            W(IOF)=WPQ(N+J+1)
            JOF=JOF-1
 390     CONTINUE
         DO 400 J=NLF/2-I,0,-1
            IOF=IOF+1
            KOF=IOF+(NLF+2)*NLF/4
            LOF=IOF-(NLF+2)*NLF/2
            MRMX(IOF)=JOF
            MRMY(IOF)=KOF
            MRMZ(IOF)=LOF
            DU(IOF)=UPQ(N+J+1)
            DE(IOF)=-VPQ(N+J+1)
            DZ(IOF)=TPQ(N+J+1)
            W(IOF)=WPQ(N+J+1)
            JOF=JOF-1
 400     CONTINUE
         N=N+NLF/2-I+1
 380  CONTINUE
      N=0
      DO 410 I=1,NLF/2
         JOF=IOF+NLF-2*I+2
         DO 420 J=0,NLF/2-I
            IOF=IOF+1
            KOF=IOF-(NLF+2)*NLF/4
            LOF=IOF-(NLF+2)*NLF/2
            MRMX(IOF)=JOF
            MRMY(IOF)=KOF
            MRMZ(IOF)=LOF
            DU(IOF)=-UPQ(N+J+1)
            DE(IOF)=VPQ(N+J+1)
            DZ(IOF)=TPQ(N+J+1)
            W(IOF)=WPQ(N+J+1)
            JOF=JOF-1
 420     CONTINUE
         DO 430 J=NLF/2-I,0,-1
            IOF=IOF+1
            KOF=IOF-(NLF+2)*NLF/4
            LOF=IOF-(NLF+2)*NLF/2
            MRMX(IOF)=JOF
            MRMY(IOF)=KOF
            MRMZ(IOF)=LOF
            DU(IOF)=UPQ(N+J+1)
            DE(IOF)=VPQ(N+J+1)
            DZ(IOF)=TPQ(N+J+1)
            W(IOF)=WPQ(N+J+1)
            JOF=JOF-1
 430     CONTINUE
         N=N+NLF/2-I+1
 410  CONTINUE
      DEALLOCATE(WPQ,VPQ,UPQ,TPQ,WW,UU,JOP)
      IF(IMPX.GE.4) THEN
         WRITE(6,'(/41H SNTT3D:HEIGHT-OCTANT ANGULAR QUADRATURES:/26X,
     1   2HMU,9X,3HETA,10X,2HXI,6X,6HWEIGHT)')
         SUM=0.0
         DO 70 N=1,NPQ
         SUM=SUM+W(N)
         WRITE(6,'(1X,4I5,1P,4E12.4)') N,MRMX(N),MRMY(N),MRMZ(N),DU(N),
     1   DE(N),DZ(N),W(N)
   70    CONTINUE
         WRITE(6,'(54X,10(1H-)/52X,1P,E12.4)') SUM
      ENDIF
*----
*  IDENTIFICATION OF THE GEOMETRY.
*----
      IF(IGE.EQ.0) THEN
* ----------
*        3D CARTESIAN
* ----------
         DO 83 N=1,NPQ
         VU=DU(N)
         VE=DE(N)
         VZ=DZ(N)
         DO 82 I=1,LX
         DO 81 J=1,LY
         DO 80 K=1,LZ
         XX(I)=XXX(I+1)-XXX(I)
         YY(J)=YYY(J+1)-YYY(J)
         ZZ(K)=ZZZ(K+1)-ZZZ(K)
         DA(J,K,N)=VU*YY(J)*ZZ(K)
         DB(I,K,N)=VE*XX(I)*ZZ(K)
         DC(I,J,N)=VZ*XX(I)*YY(J)
         VOL(I,J,K)=XX(I)*YY(J)*ZZ(K)
   80    CONTINUE
   81    CONTINUE
   82    CONTINUE
   83    CONTINUE
      ELSEIF(IGE.EQ.2) THEN
* ----------
*        3D HEXAGONAL
* ----------
         DET = SQRT(3.0)*(SIDE**2)/2.0 
         DO 93 N=1,NPQ
         VU=DU(N)
         VE=DE(N)
         VZ=DZ(N)
         DO 92 K=1,LZ
         DO 91 J=1,LY
         DO 90 I=1,LX
         ZZ(K)=ZZZ(K+1)-ZZZ(K)
         DA(J,K,N)=VU*ZZ(K)
         DB(I,K,N)=VE*ZZ(K)
         DC(I,J,N)=VZ*DET
         VOL(I,J,K)=DET*ZZ(K)
   90    CONTINUE
   91    CONTINUE
   92    CONTINUE
   93    CONTINUE
      ENDIF
*----
*  GENERATE SPHERICAL HARMONICS FOR SCATTERING SOURCE.
*----
      IOF=0
      DO 211 L=0,ISCAT-1
      DO 210 M=-L,L   
      IOF=IOF+1
      IF(IOF.GT.ISCAT) GO TO 211
      DO 200 N=1,NPQ
      PL(IOF,N)=PNSH(L,M,DU(N),DE(N),DZ(N))
  200 CONTINUE
  210 CONTINUE
  211 CONTINUE
*----
*  GENERATE MAPPING MATRIX FOR GALERKIN QUADRATURE METHOD
*----
      MN(:NPQ,:NSCT)=0.0
      DN(:NSCT,:NPQ)=0.0
      IL(:NSCT)=0
      IM(:NSCT)=0
      IF(IGLK.NE.0) THEN
         ALLOCATE(U(NPQ,NPQ),RLM(NPQ,ISCAT,2*ISCAT-1),V(NPQ),V2(NPQ))
         RLM(:NPQ,:ISCAT-1,:2*ISCAT-1)=0.0
         DO L=0,ISCAT-1
         DO M=-L,L
         DO N=1,NPQ
            RLM(N,L+1,M+L+1)=PNSH(L,M,DU(N),DE(N),DZ(N))
         ENDDO
         ENDDO
         ENDDO
         ! GRAM-SCHMIDT PROCEDURE TO FIND INDEPENDANT SET
         ! OF SPHERICAL HARMONICS WITH ANY QUADRATURE
         U(:NPQ,:NPQ)=0.0D0
         NORM=0.0D0
         DO N=1,NPQ
            NORM=NORM+RLM(N,1,1)**2
         ENDDO
         NORM=SQRT(NORM)
         DO N=1,NPQ 
            IF(IGLK.EQ.1) THEN
               MND(1,N)=W(N)*RLM(N,1,1)
            ELSEIF(IGLK.EQ.2) THEN
               MND(N,1)=(2.0*L+1.0)/(4.0*PI)*RLM(N,1,1)
            ELSE
               CALL XABORT('UNKNOWN GALERKIN QUADRATURE METHOD.')
            ENDIF
            U(N,1)=RLM(N,1,1)/NORM
         ENDDO
         IND=1
         ! ITERATE OVER THE SPHERICAL HARMONICS
         DO 212 L=1,ISCAT-1
         DO 213 M=-L,L
         V2(:NPQ)=0.0D0
         DO N=1,IND
         IPROD=0.0D0
         DO N2=1,NPQ
         IPROD=IPROD+U(N2,N)*RLM(N2,L+1,M+L+1)
         ENDDO
         DO N2=1,NPQ
         V2(N2)=V2(N2)+IPROD*U(N2,N)
         ENDDO
         ENDDO
         V(:NPQ)=0.0D0
         DO N=1,NPQ
         V(N)=RLM(N,L+1,M+L+1)-V2(N)
         ENDDO
         NORM=0.0D0
         DO N=1,NPQ
            NORM=NORM+V(N)**2
         ENDDO
         NORM=SQRT(NORM)
         ! KEEP THE SPHERICAL HARMONICS IF IT IS INDEPENDANT
         IF(NORM.GE.1.0E-5) THEN
            IND=IND+1
            DO N=1,NPQ
               U(N,IND)=V(N)/NORM
               IF(IGLK.EQ.1) THEN
                  MND(IND,N)=W(N)*RLM(N,L+1,M+L+1)
               ELSEIF(IGLK.EQ.2) THEN
                  MND(N,IND)=(2.0*L+1.0)/(4.0*PI)*RLM(N,L+1,M+L+1)
               ELSE
                  CALL XABORT('UNKNOWN GALERKIN QUADRATURE METHOD.')
               ENDIF
            ENDDO
            IL(IND)=L
            IM(IND)=M
         ENDIF
         IF(IND.EQ.NPQ) GOTO 217
  213    ENDDO
  212    ENDDO
         CALL XABORT('SNTT3D: THE'// 
     1   ' GRAM-SCHMIDTH PROCEDURE TO FIND A SUITABLE INTERPOLATION'//
     2   ' BASIS REQUIRE HIGHER LEGENDRE ORDER.')
         ! FIND INVERSE MATRIX    
  217    IF(IGLK.EQ.1) THEN    
            DN=REAL(MND)
            CALL ALINVD(NPQ,MND,NPQ,IER)
            IF(IER.NE.0) CALL XABORT('SNTT3D: SINGULAR MATRIX.')
            MN=REAL(MND)
         ELSEIF(IGLK.EQ.2) THEN
            MN=REAL(MND)
            CALL ALINVD(NPQ,MND,NPQ,IER)
            IF(IER.NE.0) CALL XABORT('SNTT3D: SINGULAR MATRIX.')
            DN=REAL(MND)
         ELSE
            CALL XABORT('UNKNOWN GALERKIN QUADRATURE METHOD.')
         ENDIF
         DEALLOCATE(U,RLM,V,V2)
      ELSE
         IND=1
         DO L=0,ISCAT-1
         DO 218 M=-L,L
         IL(IND)=L
         IM(IND)=M
         DO N=1,NPQ
         DN(IND,N)=W(N)*PNSH(L,M,DU(N),DE(N),DZ(N))
         MN(N,IND)=(2.0*L+1.0)/(4.0*PI)
     1   *PNSH(L,M,DU(N),DE(N),DZ(N))
         ENDDO
         IND=IND+1
  218    ENDDO
         ENDDO
      ENDIF
*----
* GENERATE THE WEIGHTING PARAMETERS OF THE CLOSURE RELATION.
*----
      PX=1
      PE=1
      IF(ISCHM.EQ.1.OR.ISCHM.EQ.3) THEN
        PX=1
      ELSEIF(ISCHM.EQ.2) THEN
        PX=0
      ELSE
        CALL XABORT('SNTT3D: UNKNOWN TYPE OF SPATIAL CLOSURE RELATION.')
      ENDIF
      IF(MOD(IELEM,2).EQ.1) THEN
        WX(1)=-PX
        WX(2:IELEM+1:2)=1+PX
        IF(IELEM.GE.2) WX(3:IELEM+1:2)=1-PX
      ELSE
        WX(1)=PX
        WX(2:IELEM+1:2)=1-PX
        IF(IELEM.GE.2) WX(3:IELEM+1:2)=1+PX
      ENDIF  
      IF(IBFP.NE.0) THEN
      IF(ESCHM.EQ.1.OR.ESCHM.EQ.3) THEN
        PE=1
      ELSEIF(ESCHM.EQ.2) THEN
        PE=0
      ELSE
        CALL XABORT('SNTT3D: UNKNOWN TYPE OF ENERGY CLOSURE RELATION.')
      ENDIF
      IF(MOD(EELEM,2).EQ.1) THEN
        WE(1)=-PE
        WE(2:EELEM+1:2)=1+PE
        IF(EELEM.GE.2) WE(3:EELEM+1:2)=1-PE
      ELSE
        WE(1)=PE
        WE(2:EELEM+1:2)=1-PE
        IF(EELEM.GE.2) WE(3:EELEM+1:2)=1+PE
      ENDIF
      ENDIF
      ! NORMALIZED LEGENDRE POLYNOMIAL CONSTANTS
      DO IEL=1,MAX(IELEM,EELEM)
        CST(IEL)=SQRT(2.0*IEL-1.0)
      ENDDO
*----
*  COMPUTE ISOTROPIC FLUX INDICES.
*----
      NM=IELEM**3*EELEM
      NMX=IELEM**2*EELEM
      NMY=NMX
      NMZ=NMX
      NME=IELEM**3
      LL4=LX*LY*LZ*NSCT*NM
      IF(IGE.EQ.0)THEN
        NUN=LL4+(NMX*LY*LZ+NMY*LX*LZ+NMZ*LX*LY)*NPQ
      ELSEIF(IGE.EQ.2)THEN
        IF((NCODE(5)==1).and.(NCODE(5)==NCODE(6)))THEN
           NUN=LL4
        ELSE
           NUN=LL4+(LX*LY*LZ)*NMZ*NPQ
        ENDIF
      ELSE
         CALL XABORT('SNTT3D: CHECK SPATIAL SCHEME DISCRETISATION '//
     1       'PARAMETER.')
      ENDIF
      DO I=1,LX*LY*LZ
         IDL(I)=(I-1)*NSCT*NM+1
      ENDDO
*----
*  SET BOUNDARY CONDITIONS.
*----
      DO 240 I=1,6
      IF(NCODE(I).NE.1) ZCODE(I)=1.0
      IF(NCODE(I).EQ.5) CALL XABORT('SNTT3D: SYME BC NOT ALLOWED.')
      IF(NCODE(I).EQ.7) CALL XABORT('SNTT3D: ZERO FLUX BC NOT ALLOWED.')
  240 CONTINUE
*----
*  CHECK FOR INVALID VIRTUAL ELEMENTS.
*----
      DO 292 I=2,LX-1
      DO 291 J=2,LY-1
      DO 290 K=2,LZ-1
      IF(MAT(I,J,K).EQ.0) THEN
         L1=(NCODE(1).NE.1)
         DO 251 J1=1,J-1
         DO 250 K1=1,K-1
         L1=L1.OR.(MAT(I,J1,K1).NE.0)
  250    CONTINUE
  251    CONTINUE
         L2=(NCODE(2).NE.1)
         DO 261 J1=J+1,LY
         DO 260 K1=K+1,LZ
         L2=L2.OR.(MAT(I,J1,K1).NE.0)
  260    CONTINUE
  261    CONTINUE
         L3=(NCODE(3).NE.1)
         DO 271 I1=1,I-1
         DO 270 K1=1,K-1
         L3=L3.OR.(MAT(I1,J,K1).NE.0)
  270    CONTINUE
  271    CONTINUE
         L4=(NCODE(4).NE.1)
         DO 281 I1=I+1,LX
         DO 280 K1=K+1,LZ
         L4=L4.OR.(MAT(I1,J,K1).NE.0)
  280    CONTINUE
  281    CONTINUE
         L5=(NCODE(5).NE.1)
         DO 301 I1=1,I-1
         DO 300 J1=1,J-1
         L5=L5.OR.(MAT(I1,J1,K).NE.0)
  300    CONTINUE
  301    CONTINUE
         L6=(NCODE(6).NE.0)
         DO 311 I1=I+1,LX
         DO 310 J1=I+1,LY
         L6=L6.OR.(MAT(I1,J1,K).NE.0)
  310    CONTINUE
  311    CONTINUE
         IF(L1.AND.L2.AND.L3.AND.L4.AND.L5.AND.L6) THEN
            WRITE(HSMG,'(17HSNTT3D: ELEMENT (,I3,1H,,I3,11H) CANNOT BE,
     1      9H VIRTUAL.)') I,J,K
            CALL XABORT(HSMG)
         ENDIF
      ENDIF
  290 CONTINUE
  291 CONTINUE
  292 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(ZZ,YY,XX)
      RETURN
      END
