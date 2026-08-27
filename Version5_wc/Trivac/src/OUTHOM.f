*DECK OUTHOM
      SUBROUTINE OUTHOM(MAXNEL,IPGEOM,IPMTX,IMPX,NEL,IELEM,ICOL,HTRACK,
     1 MAT,NZS,IHOM)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Read an modify the merge indices.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* MAXNEL  maximum number of elements.
* IPGEOM  L_GEOM pointer to the geometry.
* IPMTX   L_MATEX pointer to matex information.
* IMPX    print parameter.
* NEL     total number of finite elements.
* IELEM   degree of the Lagrangian finite elements:
* ICOL    type of quadrature used to integrate the mass matrix
* HTRACK  type of tracking (equal to 'BIVAC' or 'TRIVAC').
* MAT     index-number of the mixture type assigned to each volume.
*
*Parameters: output
* NZS     number of merged regions.
* IHOM    merge indices.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPGEOM,IPMTX
      CHARACTER HTRACK*12
      INTEGER MAXNEL,IMPX,NEL,IELEM,ICOL,MAT(NEL),NZS,IHOM(NEL)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      INTEGER ISTATE(NSTATE),NCODE(6),ICODE(6)
      REAL ZCODE(6)
      LOGICAL ILK,LDIAG,CHEX,LFOLD1,LFOLD2,LFOLD3
      DOUBLE PRECISION DFLOTT
      CHARACTER TEXT4*4,HSMG*131
      INTEGER, DIMENSION(:), ALLOCATABLE :: MAT2,DPP,MX,XXX,YYY,ZZZ
      INTEGER, DIMENSION(:), ALLOCATABLE :: ISPLX,ISPLY,ISPLZ
      EQUIVALENCE (ITYPE,ISTATE(1)),(LR1,ISTATE(2)),(LX1,ISTATE(3)),
     1 (LY1,ISTATE(4)),(LZ1,ISTATE(5))
*----
*  DETERMINE THE MESH SPLITTING INFO FROM THE GEOMETRY.
*----
      ALLOCATE(ISPLX(MAXNEL),ISPLY(MAXNEL),ISPLZ(MAXNEL))
      ALLOCATE(XXX(MAXNEL+1),YYY(MAXNEL+1),ZZZ(MAXNEL+1))
*
      ALLOCATE(MAT2(MAXNEL))
      CALL READ3D(MAXNEL,MAXNEL,MAXNEL,MAXNEL,IPGEOM,IHEX,IR,ILK,SIDE,
     1 XXX,YYY,ZZZ,IMPX,LX,LY,LZ,MAT2,IPAS,NCODE,ICODE,ZCODE,ISPLX,
     2 ISPLY,ISPLZ,ISPLH,ISPLL)
      DEALLOCATE(MAT2)
*
      CALL LCMGET(IPGEOM,'STATE-VECTOR',ISTATE)
      LYOLD=1
      LZOLD=1
      NELOLD=0
      IF(ITYPE.EQ.2) THEN
*        1-D CARTESIAN GEOMETRY.
         LXOLD=LX1
         NELOLD=LXOLD
      ELSE IF((ITYPE.EQ.3).OR.(ITYPE.EQ.4)) THEN
*        1-D CYLINDRICAL/SPHERICAL GEOMETRY.
         LXOLD=LR1
         NELOLD=LXOLD
      ELSE IF(ITYPE.EQ.5) THEN
*        2-D CARTESIAN GEOMETRY.
         LXOLD=LX1
         LYOLD=LY1
         NELOLD=LXOLD*LYOLD
         LDIAG=.FALSE.
         DO 10 IC=1,4
         LDIAG=LDIAG.OR.(NCODE(IC).EQ.3)
   10    CONTINUE
         IF(LDIAG) NELOLD=(LXOLD+1)*LXOLD/2
      ELSE IF(ITYPE.EQ.6) THEN
*        2-D CYLINDRICAL GEOMETRY.
         LXOLD=LR1
         LZOLD=LZ1
         NELOLD=LXOLD*LZOLD
      ELSE IF(ITYPE.EQ.7) THEN
*        3-D CARTESIAN GEOMETRY.
         LXOLD=LX1
         LYOLD=LY1
         LZOLD=LZ1
         NELOLD=LXOLD*LYOLD*LZOLD
         LDIAG=.FALSE.
         DO 20 IC=1,4
         LDIAG=LDIAG.OR.(NCODE(IC).EQ.3)
   20    CONTINUE
         IF(LDIAG) NELOLD=(LXOLD+1)*LXOLD*LZOLD/2
      ELSE IF(ITYPE.EQ.8) THEN
*        2-D HEXAGONAL GEOMETRY.
         LXOLD=LX1
         NELOLD=LXOLD
      ELSE IF(ITYPE.EQ.9) THEN
*        3-D HEXAGONAL GEOMETRY.
         LXOLD=LX1
         LZOLD=LZ1
         NELOLD=LXOLD*LZOLD
      ENDIF
*----
*  READ THE MERGE INDICES.
*----
      CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
      IF(INDIC.EQ.1) THEN
        IHOM(:NELOLD)=0
        IHOM(1)=NITMA
        NZS=NITMA
        DO 30 K=2,NELOLD
        CALL REDGET(INDIC,IHOM(K),FLOTT,TEXT4,DFLOTT)
        IF(INDIC.NE.1) CALL XABORT('OUTHOM: INTEGER EXPECTED.')
        NZS=MAX(NZS,IHOM(K))
   30   CONTINUE
        IF((ITYPE.EQ.8).OR.(ITYPE.EQ.9)) CALL LCMGET(IPGEOM,'IHEX',IHEX)
      ELSE IF((INDIC.EQ.3).AND.(TEXT4.EQ.'NONE')) THEN
        NZS=NEL
        DO 40 K=1,NEL
        IHOM(K)=K
   40   CONTINUE
        GO TO 270
      ELSE IF((INDIC.EQ.3).AND.(TEXT4.EQ.'COMP')) THEN
        NZS=1
        IHOM(:NEL)=1
        GO TO 270
      ELSE IF((INDIC.EQ.3).AND.(TEXT4.EQ.'IN')) THEN
        DO 50 K=1,NELOLD
        IHOM(K)=K
   50   CONTINUE
        NZS=NELOLD
        IF((ITYPE.EQ.8).OR.(ITYPE.EQ.9)) CALL LCMGET(IPGEOM,'IHEX',IHEX)
      ELSE IF((INDIC.EQ.3).AND.(TEXT4.EQ.'MATX')) THEN
        IF(.NOT.C_ASSOCIATED(IPMTX)) THEN
          CALL XABORT('OUTHOM: MATEX LCM OBJECT IS MISSING.')
        ENDIF
        CALL LCMLEN(IPMTX,'MATEX',ILONG,ITYLCM)
        IF(ILONG.NE.NEL) THEN
          WRITE(HSMG,'(42HOUTHOM: INCONSISTENT INTG MATX OPTION (EXP,
     1    25HECTED NUMBER OF MIXTURES=,I5,23H; VALUE FOUND IN L_GEOM,
     2    8H OBJECT=,I5,2H).)') NEL,ILONG
          CALL XABORT(HSMG)
        ENDIF
        CALL LCMGET(IPMTX,'MATEX',IHOM)
        NZS=0
        DO 60 K=1,NEL
        NZS=MAX(NZS,IHOM(K))
   60   CONTINUE
        GO TO 270
      ELSE IF((INDIC.EQ.3).AND.(TEXT4.EQ.'MIX')) THEN
        CALL LCMLEN(IPGEOM,'MIX',ILONG,ITYLCM)
        IF(ILONG.NE.NELOLD) THEN
          WRITE(HSMG,'(42HOUTHOM: INCONSISTENT INTG MIX OPTION (EXPE,
     1    24HCTED NUMBER OF MIXTURES=,I5,24H; VALUE FOUND IN L_GEOM ,
     2    7HOBJECT=,I5,2H).)') NELOLD,ILONG
          CALL XABORT(HSMG)
        ENDIF
        CALL LCMGET(IPGEOM,'MIX',IHOM)
        NZS=0
        DO 70 K=1,NELOLD
        NZS=MAX(NZS,IHOM(K))
   70   CONTINUE
        IF((ITYPE.EQ.8).OR.(ITYPE.EQ.9)) CALL LCMGET(IPGEOM,'IHEX',IHEX)
      ELSE
        WRITE(HSMG,'(26HOUTHOM: INVALID KEY WORD (,A,2H).)') TEXT4
        CALL XABORT(HSMG)
      ENDIF
*----
*  UNFOLD HEXAGONAL GEOMETRY IN BIVAC AND TRIVAC CASES.
*----
      CHEX=(ITYPE.EQ.8).OR.(ITYPE.EQ.9)
      LFOLD1=CHEX.AND.(IHEX.NE.9).AND.(HTRACK.EQ.'TRIVAC')
      LFOLD2=CHEX.AND.(IHEX.NE.9).AND.(HTRACK.EQ.'BIVAC').AND.
     1       (IELEM.GT.0).AND.(ICOL.LE.3)
      LFOLD3=CHEX.AND.(IHEX.NE.9).AND.((HTRACK.EQ.'MCCG').OR.
     1       (HTRACK.EQ.'EXCELL'))
      IF(LFOLD1.OR.LFOLD2.OR.LFOLD3) THEN
         IF(NELOLD.NE.LXOLD*LZOLD) CALL XABORT('OUTHOM: HEXAGONAL SPLI'
     1   //'T ERROR.')
         ALLOCATE(DPP(MAXNEL),MX(NELOLD))
         MX(:NELOLD)=IHOM(:NELOLD)
         LXOLD=LX1
         CALL BIVALL(MAXNEL,IHEX,LXOLD,LX,DPP)
         DO 90 KZ=1,LZOLD
         DO 80 KX=1,LX
         IHOM(KX+(KZ-1)*LX)=0
         KEL=DPP(KX)+(KZ-1)*LXOLD
         IF(KEL.GT.LXOLD*LZOLD) CALL XABORT('OUTHOM: MX OVERFLOW.')
         IHOM(KX+(KZ-1)*LX)=MX(KEL)
   80    CONTINUE
   90    CONTINUE
         DEALLOCATE(MX,DPP)
         LXOLD=LX
         IHEX=9
      ENDIF
*----
*  MESH-SPLITTING FOR THE IHOM VECTOR.
*----
      IF(NZS.GT.NELOLD) CALL XABORT('OUTHOM: FAILURE 1.')
      IF(ISTATE(11).NE.0) THEN
         CALL SPLIT0(MAXNEL,ITYPE,NCODE,LXOLD,LYOLD,LZOLD,ISPLX,ISPLY,
     1   ISPLZ,0,ISPLL,NEL2,LX,LY,LZ,SIDE,XXX,YYY,ZZZ,IHOM,.FALSE.,IMPX)
      ENDIF
*----
*  FORCE DIAGONAL SYMMETRY AND UNFOLD THE IHOM VECTOR.
*----
      IF((NCODE(2).EQ.3).AND.(NCODE(3).EQ.3)) THEN
         IF(NEL.EQ.LX*LY*LZ) THEN
            K=(LX*(LX+1)/2)*LZ
            DO 232 IZ=LZ,1,-1
            IOFF=(IZ-1)*LX*LY
            DO 231 IY=LY,1,-1
            DO 220 IX=LX,IY+1,-1
            IHOM(IOFF+(IY-1)*LX+IX)=IHOM(IOFF+(IX-1)*LY+IY)
  220       CONTINUE
            DO 230 IX=IY,1,-1
            IHOM(IOFF+(IY-1)*LX+IX)=IHOM(K)
            K=K-1
  230       CONTINUE
  231       CONTINUE
  232       CONTINUE
            IF(K.NE.0) CALL XABORT('OUTHOM: FAILURE 2.')
         ENDIF
      ELSE IF((NCODE(1).EQ.3).AND.(NCODE(4).EQ.3)) THEN
         IF(NEL.EQ.LX*LY*LZ) THEN
            K=(LX*(LX+1)/2)*LZ
            DO 242 IZ=LZ,1,-1
            IOFF=(IZ-1)*LX*LY
            DO 241 IY=LY,1,-1
            DO 240 IX=LX,IY,-1
            IHOM(IOFF+(IY-1)*LX+IX)=IHOM(K)
            K=K-1
  240       CONTINUE
  241       CONTINUE
  242       CONTINUE
            DO 252 IZ=1,LZ
            IOFF=(IZ-1)*LX*LY
            DO 251 IY=1,LY
            DO 250 IX=1,IY-1
            IHOM(IOFF+(IY-1)*LX+IX)=IHOM(IOFF+(IX-1)*LY+IY)
  250       CONTINUE
  251       CONTINUE
  252       CONTINUE
            IF(K.NE.0) CALL XABORT('OUTHOM: FAILURE 3.')
         ENDIF
      ENDIF
      DEALLOCATE(ZZZ,YYY,XXX,ISPLZ,ISPLY,ISPLX)
      DO 260 K=1,NEL
      IF(MAT(K).EQ.0) IHOM(K)=0
  260 CONTINUE
  270 IF(IMPX.GT.1) THEN
        WRITE(6,'(/15H MERGING INDEX:/(1X,14I5))') (IHOM(K),K=1,NEL)
      ENDIF
      RETURN
      END
