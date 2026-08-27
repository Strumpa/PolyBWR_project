*DECK USPD3D
      SUBROUTINE USPD3D (MAXX,MAXY,MAXZ,MAXPTS,IPGEOM,IMPX,LX,LY,LZ,KSP)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover the correspondance between nodes before and after mesh
* splitting.
*
*Copyright:
* Copyright (C) 2026 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* MAXX    allocated storage for arrays of dimension LX.
* MAXY    allocated storage for arrays of dimension LY.
* MAXZ    allocated storage for arrays of dimension LZ.
* MAXPTS  allocated storage for arrays of dimension NMBLK.
* IPGEOM  L_GEOM pointer to the geometry.
* IMPX    print flag. Minimum printing if IMPX=0.
* LX      number of elements along the X-axis after mesh-splitting
*         or number of hexagons in one axial plane.
* LY      number of elements along the Y-axis.
* LZ      number of elements along the Z-axis.
*
*Parameters: output
* KSP     correspondance between nodes before and after mesh splitting.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPGEOM
      INTEGER, INTENT(IN) ::  MAXX,MAXY,MAXZ,MAXPTS,IMPX,LX,LY,LZ
      INTEGER, INTENT(OUT) :: KSP(MAXPTS)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: DPP,MX,ISPLTX,ISPLTY,ISPLTZ
      REAL, ALLOCATABLE, DIMENSION(:) :: XXX,YYY,ZZZ
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40)
      LOGICAL LL1,LL2,LCYL
      CHARACTER HSMG*131
      INTEGER ISTATE(NSTATE),NCODE(6)
      EQUIVALENCE (ITYPE,ISTATE(1)),(LR1,ISTATE(2)),(LX1,ISTATE(3)),
     1 (LY1,ISTATE(4)),(LZ1,ISTATE(5)),(NMBLK,ISTATE(6))
*
      IHEX=0
      CALL LCMGET(IPGEOM,'STATE-VECTOR',ISTATE)
      IF((ITYPE.EQ.8).OR.(ITYPE.EQ.9)) THEN
        CALL LCMLEN(IPGEOM,'IHEX',ILEN,ITYLCM)
        IF(ILEN.EQ.0) CALL LCMLIB(IPGEOM)
        IF(ILEN.EQ.0) CALL XABORT('USPD3D: MISSING IHEX RECORD.')
        CALL LCMGET(IPGEOM,'IHEX',IHEX)
      ENDIF
      IF((ISTATE(8).NE.0).OR.(ISTATE(9).NE.0).OR.(ISTATE(10).NE.0).OR.
     1 (ISTATE(13).NE.0)) CALL XABORT('USPD3D: UNABLE TO PROCESS THE G'
     2 //'EOMETRY.')
      LCYL=(ITYPE.EQ.3).OR.(ITYPE.EQ.4).OR.(ITYPE.EQ.6)
      IDIM=1
      IF((ITYPE.EQ.5).OR.(ITYPE.EQ.6).OR.(ITYPE.EQ.8)) IDIM=2
      IF((ITYPE.EQ.7).OR.(ITYPE.EQ.9)) IDIM=3
*----
*  RECOVER THE BOUNDARY CONDITIONS
*----
      CALL LCMGET(IPGEOM,'NCODE',NCODE)
      DO 10 I=1,6
      IF(NCODE(I).EQ.10) NCODE(I)=2
      IF(NCODE(I).EQ.6) NCODE(I)=1
      IF((NCODE(I).EQ.20).AND.(ITYPE.NE.5).AND.(ITYPE.NE.7)) CALL
     1 XABORT('USPD3D: CYLINDRICAL CORRECTION IS LIMITED TO CARTESIAN '
     2 //'GEOMETRIES.')
      IF((NCODE(I).GE.8).AND.(NCODE(I).NE.20)) THEN
         CALL XABORT('USPD3D: INVALID TYPE OF B.C.')
      ENDIF
   10 CONTINUE
*
      IF(ITYPE.GE.10) THEN
         CALL XABORT('USPD3D: INVALID TYPE OF GEOMETRY.')
      ELSE IF(ITYPE.GE.8) THEN
         IF((NCODE(2).NE.0).OR.(NCODE(3).NE.0).OR.(NCODE(4).NE.0))
     1   CALL XABORT('USPD3D: INVALID TYPE OF HEXAGONAL B.C.')
         IF(NCODE(1).EQ.5) THEN
            IF(IHEX.EQ.1) THEN
               IHEX=10
            ELSE IF(IHEX.EQ.2) THEN
               IHEX=11
            ELSE
               CALL XABORT('USPD3D: BOUNDARY CONDITION HBC WITH OPTION'
     1         //' SYME IS ONLY PERMITTED WITH S30 OR SA60 SYMMETRY.')
            ENDIF
         ELSE IF((NCODE(1).GT.2).AND.(NCODE(1).NE.7)) THEN
            CALL XABORT('USPD3D: BOUNDARY CONDITION HBC CAN ONLY BE US'
     1      //'ED WITH OPTIONS VOID, REFL, SYME, ALBE OR ZERO.')
         ENDIF
      ENDIF
      IF(NMBLK.GT.MAXPTS) THEN
         WRITE (HSMG,690) 'NMBLK',NMBLK,'MAXPTS',MAXPTS
         CALL XABORT(HSMG)
      ENDIF
      DO 20 I=1,NMBLK
      KSP(I)=I
   20 CONTINUE
*----
*  RECOVER THE MESH COORDINATES.
*----
      ALLOCATE(XXX(MAXX+1),YYY(MAXY+1),ZZZ(MAXZ+1))
      IF(LCYL.AND.(LR1.GT.MAXX)) THEN
         WRITE (HSMG,690) 'LXOLD',LR1,'MAXX',MAXX
         CALL XABORT(HSMG)
      ELSE IF(LX1.GT.MAXX) THEN
         WRITE (HSMG,690) 'LXOLD',LX1,'MAXX',MAXX
         CALL XABORT(HSMG)
      ENDIF
      IF(LY1.GT.MAXY) THEN
         WRITE (HSMG,690) 'LYOLD',LY1,'MAXY',MAXY
         CALL XABORT(HSMG)
      ENDIF
      IF(LZ1.GT.MAXZ) THEN
         WRITE (HSMG,690) 'LZOLD',LZ1,'MAXZ',MAXZ
         CALL XABORT(HSMG)
      ENDIF
      LL1=.FALSE.
      LL2=.FALSE.
      LYOLD=1
      YYY(1)=0.0
      YYY(2)=1.0
      LZOLD=1
      ZZZ(1)=0.0
      ZZZ(2)=1.0
      IF(ITYPE.EQ.2) THEN
*        1-D CARTESIAN GEOMETRY.
         LXOLD=LX1
         IF((NCODE(1).EQ.0).OR.(NCODE(2).EQ.0)) GO TO 610
         CALL LCMGET(IPGEOM,'MESHX',XXX)
      ELSE IF((ITYPE.EQ.3).OR.(ITYPE.EQ.4)) THEN
*        1-D CYLINDRICAL/SPHERICAL GEOMETRY.
         LXOLD=LR1
         IF(NCODE(1).NE.0) GO TO 640
         IF(NCODE(2).EQ.0) GO TO 610
         NCODE(1)=2
         CALL LCMGET(IPGEOM,'RADIUS',XXX)
      ELSE IF(ITYPE.EQ.5) THEN
*        2-D CARTESIAN GEOMETRY.
         LXOLD=LX1
         LYOLD=LY1
         I2=0
         DO 30 IC=1,4
         IF(NCODE(IC).EQ.0) GO TO 610
         IF(NCODE(IC).EQ.3) I2=I2+1
   30    CONTINUE
         IF(I2.NE.0) THEN
            IF((I2.NE.2).OR.(LXOLD.NE.LYOLD)) GO TO 630
            LL1=(NCODE(2).EQ.3).AND.(NCODE(3).EQ.3)
            LL2=(NCODE(1).EQ.3).AND.(NCODE(4).EQ.3)
            IF((.NOT.LL1).AND.(.NOT.LL2)) GO TO 620
         ENDIF
         CALL LCMGET(IPGEOM,'MESHX',XXX)
         IF(LL1.OR.LL2) THEN
            CALL LCMGET(IPGEOM,'MESHX',YYY)
         ELSE
            CALL LCMGET(IPGEOM,'MESHY',YYY)
         ENDIF
      ELSE IF(ITYPE.EQ.6) THEN
*        2-D CYLINDRICAL GEOMETRY.
         LXOLD=LR1
         LZOLD=LZ1
         IF(NCODE(1).NE.0) GO TO 650
         IF((NCODE(2).EQ.3).OR.(NCODE(3).EQ.3).OR.(NCODE(4).EQ.3))
     1   GO TO 660
         IF((NCODE(2).EQ.0).OR.(NCODE(5).EQ.0).OR.(NCODE(6).EQ.0))
     1   GO TO 610
         NCODE(1)=2
         CALL LCMGET(IPGEOM,'RADIUS',XXX)
         CALL LCMGET(IPGEOM,'MESHZ',ZZZ)
      ELSE IF(ITYPE.EQ.7) THEN
*        3-D CARTESIAN GEOMETRY.
         LXOLD=LX1
         LYOLD=LY1
         LZOLD=LZ1
         I2=0
         DO 40 IC=1,4
         IF(NCODE(IC).EQ.0) GO TO 610
         IF(NCODE(IC).EQ.3) I2=I2+1
   40    CONTINUE
         IF(I2.NE.0) THEN
            IF((I2.NE.2).OR.(LXOLD.NE.LYOLD)) GO TO 630
            LL1=(NCODE(2).EQ.3).AND.(NCODE(3).EQ.3)
            LL2=(NCODE(1).EQ.3).AND.(NCODE(4).EQ.3)
            IF((.NOT.LL1).AND.(.NOT.LL2)) GO TO 620
         ENDIF
         CALL LCMGET(IPGEOM,'MESHX',XXX)
         IF(LL1.OR.LL2) THEN
            CALL LCMGET(IPGEOM,'MESHX',YYY)
         ELSE
            CALL LCMGET(IPGEOM,'MESHY',YYY)
         ENDIF
         CALL LCMGET(IPGEOM,'MESHZ',ZZZ)
      ELSE IF(ITYPE.EQ.8) THEN
*        2-D HEXAGONAL GEOMETRY.
         LXOLD=LX1
         CALL LCMGET(IPGEOM,'SIDE',SIDE)
      ELSE IF(ITYPE.EQ.9) THEN
*        3-D HEXAGONAL GEOMETRY.
         LXOLD=LX1
         LZOLD=LZ1
         CALL LCMGET(IPGEOM,'SIDE',SIDE)
         CALL LCMGET(IPGEOM,'MESHZ',ZZZ)
      ELSE
         CALL XABORT('USPD3D: INVALID TYPE OF GEOMETRY.')
      ENDIF
*----
*  UNFOLD GEOMETRY IF HEXAGONAL IN LOZENGES
*----
      ISPLTL=0
      ISPLTH=0
      CALL LCMLEN(IPGEOM,'SPLITL',ILEN,ITYLCM)
      IF(ILEN.GT.0) CALL LCMGET(IPGEOM,'SPLITL',ISPLTL)
      CALL LCMLEN(IPGEOM,'SPLITH',ILEN,ITYLCM)
      IF(ILEN.GT.0) CALL LCMGET(IPGEOM,'SPLITH',ISPLTH)
      IF((ISPLTL.GT.0).AND.(IHEX.NE.9)) THEN
         ALLOCATE(DPP(MAXPTS),MX(LX*LZ))
         MX(:LX*LZ)=KSP(:LX*LZ)
         CALL BIVALL(MAXPTS,IHEX,LXOLD,LXOLD,DPP)
         DO 80 KZ=1,LZ
         DO 70 KX=1,LXOLD
         KEL=DPP(KX)+(KZ-1)*LXOLD
         KSP(KX+(KZ-1)*LX)=MX(KEL)
   70    CONTINUE
   80    CONTINUE
         DEALLOCATE(DPP,MX)
         IHEX=9
      ENDIF
*----
*  MESH-SPLITTING
*----
      ALLOCATE(ISPLTX(MAXX),ISPLTY(MAXY),ISPLTZ(MAXZ))
      IF(ISTATE(11).NE.0) THEN
         CALL LCMLEN(IPGEOM,'SPLITR',ILEN1,ITYLCM)
         CALL LCMLEN(IPGEOM,'SPLITX',ILEN2,ITYLCM)
         IF(LCYL.AND.(ILEN1.GT.0)) THEN
            CALL LCMGET(IPGEOM,'SPLITR',ISPLTX)
         ELSE IF(ILEN2.GT.0) THEN
            CALL LCMGET(IPGEOM,'SPLITX',ISPLTX)
         ELSE IF(ITYPE.LE.7) THEN
            ISPLTX(:LX)=1
         ENDIF
         CALL LCMLEN(IPGEOM,'SPLITY',ILEN,ITYLCM)
         IF(ILEN.GT.0) THEN
            CALL LCMGET(IPGEOM,'SPLITY',ISPLTY)
         ELSE IF(LL1.OR.LL2) THEN
            ISPLTY(:LX)=ISPLTX(:LX)
         ELSE
            ISPLTY(:LY)=1
         ENDIF
         CALL LCMLEN(IPGEOM,'SPLITZ',ILEN,ITYLCM)
         IF(ILEN.GT.0) THEN
            CALL LCMGET(IPGEOM,'SPLITZ',ISPLTZ)
         ELSE
            ISPLTZ(:LZ)=1
         ENDIF
         IF((ISPLTH.GT.0).AND.(ISPLTL.GT.0)) THEN
            CALL XABORT('USPD3D: SPLITH AND SPLITL KEYWORDS ARE EXCLUS'
     1      //'IVE.')
         ENDIF
         CALL SPLIT0(MAXPTS,ITYPE,NCODE,LXOLD,LYOLD,LZOLD,ISPLTX,ISPLTY,
     1   ISPLTZ,0,ISPLTL,NMBLK,LX,LY,LZ,SIDE,XXX,YYY,ZZZ,KSP,.FALSE.,
     2   IMPX)
         IF(NMBLK.GT.MAXPTS) THEN
            WRITE (HSMG,690) 'NMBLK',NMBLK,'MAXPTS',MAXPTS
            CALL XABORT(HSMG)
         ENDIF
      ENDIF
      DEALLOCATE(ISPLTZ,ISPLTY,ISPLTX)
      DEALLOCATE(ZZZ,YYY,XXX)
      RETURN
*
  610 CALL XABORT('USPD3D: A BOUNDARY CONDITION IS MISSING.')
  620 CALL XABORT('USPD3D: THE DIAGONAL CONDITIONS X+: DIAG Y-: DIAG A'
     1 //'ND X-: DIAG Y+: DIAG ARE THE ONLY PERMITTED.')
  630 CALL XABORT('USPD3D: LX=LY WITH A DIAGONAL SYMMETRY.')
  640 CALL XABORT('USPD3D: CYLINDRICAL GEOMETRY - ONLY THE R+: BOUNDAR'
     1 //'Y CONDITION IS REQUIRED.')
  650 CALL XABORT('USPD3D: CYLINDRICAL GEOMETRY - ONLY THE R+:, Z-: AN'
     1 //'D Z+: BOUNDARY CONDITIONS ARE REQUIRED.')
  660 CALL XABORT('USPD3D: CYLINDRICAL GEOMETRY : THE DIAG BOUNDARY CO'
     1 //'NDITION CANNOT BE USED.')
*
  690 FORMAT (29HUSPD3D: INSUFFICIENT STORAGE.,5X,A6,1H=,I7,8H ; AVAIL,
     1 13HABLE STORAGE ,A6,1H=,I7)
      END
