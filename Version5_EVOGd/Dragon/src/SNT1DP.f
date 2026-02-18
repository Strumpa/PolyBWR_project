*DECK SNT1DP
      SUBROUTINE SNT1DP (IMPX,LX,IELEM,NCODE,ZCODE,XXX,NLF,NSCT,U,W,PL,
     1 VOL,IDL,LL4,NUN,LSHOOT,EELEM,WX,WE,CST,IBFP,ISCHM,ESCHM,IGLK,MN,
     2 DN,IL,IM,IQUAD)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Numbering corresponding to a 1-D slab geometry with discrete
* ordinates approximation of the flux.
*
*Copyright:
* Copyright (C) 2005 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert and C. Bienvenue
*
*Parameters: input
* IMPX    print parameter.
* LX      number of elements along the X axis.
* IELEM   measure of order of the spatial approximation polynomial:
*         =1 constant - default for HODD;
*         =2 linear - default for DG;
*         >3 higher orders.
* NCODE   type of boundary condition applied on each side
*         (i=1 X-;  i=2 X+):
*         =1 VOID;  =2 REFL;  =7 ZERO.
* ZCODE   ZCODE(I) is the albedo corresponding to boundary condition
*         'VOID' on each side (ZCODE(I)=0.0 by default).
* XXX     coordinates along the R axis.
* NLF     number of $\\mu$ levels.
* NSCT    maximum number of Legendre components in the flux.
* LSHOOT  enablig flag for the shooting method.
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
* IQUAD   quadrature type:
*         =1 Gauss-Lobatto;
*         =2 Gauss-Legendre.
*
*Parameters: output
* U       base points in $\\mu$ of the 1D quadrature.
* W       weights for the quadrature in $\\mu$.
* PL      discrete values of the Legendre polynomials corresponding
*         to the SN quadrature.
* VOL     volume of each element.
* IDL     position of integrated fluxes into unknown vector.
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
      INTEGER IMPX,LX,IELEM,NCODE(2),NLF,NSCT,IDL(LX),LL4,NUN,EELEM,
     1 IBFP,ISCHM,ESCHM,IL(NSCT),IM(NSCT),IQUAD,IGLK
      REAL ZCODE(2),XXX(LX+1),U(NLF),W(NLF),VOL(LX),WX(IELEM+1),
     1 WE(EELEM+1),CST(MAX(IELEM,EELEM)),PX,PE,MN(NLF,NSCT),
     2 DN(NSCT,NLF),PL(NSCT,NLF)
      LOGICAL LSHOOT
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION MND(NLF,NSCT)
*----
*  GENERATE QUADRATURE BASE POINTS AND CORRESPONDING WEIGHTS.
*----
      IF(MOD(NLF,2).EQ.1) CALL XABORT('SNT1DP: EVEN NLF EXPECTED.')
      IF(IQUAD.EQ.1) THEN
         CALL SNQU07(NLF,U,W)         ! GAUSS-LOBATTO
      ELSEIF(IQUAD.EQ.2) THEN 
         CALL ALGPT(NLF,-1.0,1.0,U,W) ! GAUSS-LEGENDRE
      ELSE
         CALL XABORT('SNT1DP: UNKNOWN QUADRATURE TYPE.')
      ENDIF
*----
*  PRINT COSINES AND WEIGHTS.
*----
      IF(IMPX.GT.1) THEN
         WRITE(6,60) (U(M),M=1,NLF)
   60    FORMAT(//,1X,'THE QUADRATURE COSINES FOLLOW'//
     1   (1X,5E14.6))
         WRITE(6,70) (W(M),M=1,NLF)
   70    FORMAT(//,1X,'THE CORRESPONDING QUADRATURE WEIGHTS FOLLOW'//
     1   (1X,5E14.6))
      ENDIF
*----
*  GENERATE LEGENDRE POLYNOMIALS FOR SCATTERING SOURCE.
*----
      DO 145 M=1,NLF
      PL(1,M)=1.0
      IF(NSCT.GT.1) THEN
         PL(2,M)=U(M)
         DO 140 L=2,NSCT-1
         PL(L+1,M)=((2.0*L-1.0)*U(M)*PL(L,M)-(L-1)*PL(L-1,M))/L
  140    CONTINUE
      ENDIF
  145 CONTINUE
*----
*  GENERATE MAPPING MATRIX FOR GALERKIN QUADRATURE METHOD
*----
      MN(:NLF,:NSCT)=0.0D0
      DN(:NSCT,:NLF)=0.0D0
      IL(:NSCT)=0
      IM(:NSCT)=0
      IF(IGLK.EQ.1) THEN
        DO L=0,NSCT-1
          IL(L+1)=L
          DO N=1,NLF
            DN(L+1,N)=W(N)*PL(L+1,N)
          ENDDO
        ENDDO
        MND=DN
        CALL ALINVD(NLF,MND,NLF,IER)
        IF(IER.NE.0) CALL XABORT('SNT1DP: SINGULAR MATRIX.')
        MN = REAL(MND)
      ELSEIF(IGLK.EQ.2) THEN
        DO L=0,NSCT-1
          IL(L+1)=L
          DO N=1,NLF
            MN(N,L+1)=(2.0*L+1.0)/2.0*PL(L+1,N)
          ENDDO
        ENDDO
        MND=MN
        CALL ALINVD(NLF,MND,NLF,IER)
        IF(IER.NE.0) CALL XABORT('SNT1DP: SINGULAR MATRIX.')
        DN = REAL(MND)
      ELSE
        DO L=0,NSCT-1
          IL(L+1)=L
          DO N=1,NLF
            DN(L+1,N)=W(N)*PL(L+1,N)
            MN(N,L+1)=(2.0*L+1.0)/2.0*PL(L+1,N)
          ENDDO
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
        CALL XABORT('SNT1DP: UNKNOWN TYPE OF SPATIAL CLOSURE RELATION.')
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
        CALL XABORT('SNT1DP: UNKNOWN TYPE OF ENERGY CLOSURE RELATION.')
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
*  COMPUTE VOLUMES, ISOTROPIC FLUX INDICES and UNKNOWNS. 
*----
      NM=IELEM*EELEM
      NMX=EELEM
      LL4=LX*NSCT*NM
      NUN=LL4
      DO I=1,LX
         IDL(I)=(I-1)*NSCT*NM+1
         VOL(I)=XXX(I+1)-XXX(I)
      ENDDO
      IF(.NOT.LSHOOT) NUN=NUN+NMX*NLF
*----
*  SET BOUNDARY CONDITIONS.
*----
      IF(NCODE(1).EQ.2) ZCODE(1)=1.0
      IF(NCODE(2).EQ.2) ZCODE(2)=1.0
      IF((NCODE(1).EQ.4).OR.(NCODE(2).EQ.4)) THEN
         IF((NCODE(1).NE.4).OR.(NCODE(2).NE.4)) CALL XABORT('SNT1DP: '
     1   //'INVALID TRANSLATION BOUNDARY CONDITIONS.')
         ZCODE(1)=1.0
         ZCODE(2)=1.0
      ENDIF
      IF((ZCODE(1).EQ.0.0).AND.(ZCODE(2).NE.0.0)) CALL XABORT('SNT1DP:'
     1 //' CANNOT SUPPORT LEFT VACUUM BC IF RIGHT BC IS NOT VACUUM.')
*
      RETURN
      END
