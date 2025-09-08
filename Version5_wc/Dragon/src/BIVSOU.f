*DECK BIVSOU
      SUBROUTINE BIVSOU(MAX1,IG,IPTRK,KPMACR,NANIS,NREG,NMAT,NUNKNO,
     > NGRP,MATCOD,VOL,FUNKNO,SUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the source for the solution of diffusion or PN equations.
* BIVAC-specific version.
*
*Copyright:
* Copyright (C) 2004 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* MAX1    first dimension of FUNKNO and SUNKNO arrays.
* IG      secondary group.
* IPTRK   pointer to the tracking LCM object.
* KPMACR  pointer to the secondary-group related macrolib information.
* NANIS   maximum cross section Legendre order.
* NREG    number of regions.
* NMAT    number of mixtures.
* NUNKNO  number of unknowns per energy group including spherical
*         harmonic terms, interface currents and fundamental
*         currents.
* NGRP    number of energy groups.
* MATCOD  mixture indices.
* VOL     volumes.
* FUNKNO  fluxes.
*
*Parameters: output
* SUNKNO  sources.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,KPMACR
      INTEGER MAX1,IG,NANIS,NREG,NMAT,NUNKNO,NGRP,MATCOD(NREG)
      REAL VOL(NREG),FUNKNO(MAX1,NGRP),SUNKNO(MAX1,NGRP)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40)
      INTEGER JPAR(NSTATE),IJ1(25),IJ2(25)
      CHARACTER CAN(0:9)*2
      LOGICAL CYLIND
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJ,NJJ,IPOS,KN,IDL
      REAL, ALLOCATABLE, DIMENSION(:) :: XSCAT,XX,DD
      REAL, ALLOCATABLE, DIMENSION(:,:) :: RR,RS
*----
*  DATA STATEMENTS
*----
      DATA CAN /'00','01','02','03','04','05','06','07','08','09'/
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IJJ(0:NMAT),NJJ(0:NMAT),IPOS(0:NMAT))
      ALLOCATE(XSCAT(0:NMAT*NGRP))
*----
*  RECOVER BIVAC SPECIFIC PARAMETERS.
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',JPAR)
      IF(JPAR(1).NE.NREG) CALL XABORT('BIVSOU: INCONSISTENT NREG.')
      IF(JPAR(2).NE.NUNKNO) CALL XABORT('BIVSOU: INCONSISTENT NUNKNO.')
      ITYPE=JPAR(6)
      IELEM=JPAR(8)
      ICOL=JPAR(9)
      L4=JPAR(11)
      LX=JPAR(12)
      NLF=JPAR(14)
      ISCAT=JPAR(16)
      CYLIND=(ITYPE.EQ.3).OR.(ITYPE.EQ.6)
      IF(ICOL.EQ.4) THEN
        CALL XABORT('BIVSOU: COLLOCATION NODAL NOT IMPLEMENTED.')
      ELSE IF((ITYPE.NE.2).AND.(ITYPE.NE.5)) THEN
        CALL XABORT('BIVSOU: CARTESIAN 1D OR 2D GEOMETRY EXPECTED.')
      ENDIF
      CALL LCMLEN(IPTRK,'KN',MAXKN,ITYLCM)
      ALLOCATE(XX(NREG),DD(NREG),KN(MAXKN),IDL(NREG))
      CALL LCMGET(IPTRK,'XX',XX)
      CALL LCMGET(IPTRK,'DD',DD)
      CALL LCMGET(IPTRK,'KN',KN)
      CALL LCMGET(IPTRK,'KEYFLX',IDL)
*----
*  RECOVER THE FINITE ELEMENT UNIT STIFFNESS MATRIX.
*----
      LL=0
      IF((NLF.GT.0).OR.(IELEM.LT.0)) THEN
         CALL LCMSIX(IPTRK,'BIVCOL',1)
         CALL LCMLEN(IPTRK,'T',LC,ITYLCM)
         ALLOCATE(RR(LC,LC),RS(LC,LC))
         CALL LCMGET(IPTRK,'R',RR)
         CALL LCMGET(IPTRK,'RS',RS)
         CALL LCMSIX(IPTRK,' ',2)
*----
*  COMPUTE VECTORS IJ1 AND IJ2
*----
         LL=LC*LC
         DO 10 I=1,LL
         IJ1(I)=1+MOD(I-1,LC)
         IJ2(I)=1+(I-IJ1(I))/LC
   10    CONTINUE
      ENDIF
*----
*  COMPUTE THE SOURCE
*----
      IF(NLF.EQ.0) THEN
*----
*  ++++ DIFFUSION THEORY ++++
*----
        CALL LCMGET(KPMACR,'NJJS00',NJJ(1))
        CALL LCMGET(KPMACR,'IJJS00',IJJ(1))
        CALL LCMGET(KPMACR,'IPOS00',IPOS(1))
        CALL LCMGET(KPMACR,'SCAT00',XSCAT(1))
        IF(IELEM.GT.0) THEN
*----
*  CARTESIAN 2D DUAL (RAVIART-THOMAS) CASE.
*----
          NUM1=0
          DO 30 IR=1,NREG
          IBM=MATCOD(IR)
          IF(IBM.LE.0) GO TO 30
          IF(VOL(IR).EQ.0.0) GO TO 26
          DO 25 I0=1,IELEM*IELEM
          IND=KN(NUM1+1)+I0-1
          JG=IJJ(IBM)
          DO 20 JND=1,NJJ(IBM)
          IF(JG.NE.IG) THEN
            SUNKNO(IND,IG)=SUNKNO(IND,IG)+FUNKNO(IND,JG)*VOL(IR)*
     >      XSCAT(IPOS(IBM)+JND-1)
          ENDIF
          JG=JG-1
   20     CONTINUE
   25     CONTINUE
   26     NUM1=NUM1+5
   30     CONTINUE
        ELSE IF(IELEM.LT.0) THEN
*----
*  CARTESIAN 2D PRIM (LAGRANGE) CASE.
*----
          NUM1=0
          DO 170 IR=1,NREG
          IBM=MATCOD(IR)
          IF(IBM.LE.0) GO TO 170
          IF(VOL(IR).EQ.0.0) GO TO 160
          DO 140 I=1,LL
          IND1=KN(NUM1+I)
          IF(IND1.EQ.0) GO TO 140
          I1=IJ1(I)
          I2=IJ2(I)
          DO 130 J=1,LL
          IND2=KN(NUM1+J)
          IF(IND2.EQ.0) GO TO 130
          IF(CYLIND) THEN
            AUXX=(RR(I1,IJ1(J))+RS(I1,IJ1(J))*XX(IR)/DD(IR))*
     >      RR(I2,IJ2(J))*VOL(IR)
          ELSE
            AUXX=RR(I1,IJ1(J))*RR(I2,IJ2(J))*VOL(IR)
          ENDIF
          JG=IJJ(IBM)
          DO 120 JND=1,NJJ(IBM)
          IF(JG.NE.IG) THEN
            SUNKNO(IND1,IG)=SUNKNO(IND1,IG)+AUXX*FUNKNO(IND2,JG)*
     >      XSCAT(IPOS(IBM)+JND-1)
          ENDIF
          JG=JG-1
  120     CONTINUE
  130     CONTINUE
  140     CONTINUE
          ! append the integrated volumic sources
          JG=IJJ(IBM)
          DO 150 JND=1,NJJ(IBM)
          SUNKNO(IDL(IR),IG)=SUNKNO(IDL(IR),IG)+FUNKNO(IDL(IR),JG)*
     >    VOL(IR)*XSCAT(IPOS(IBM)+JND-1)
          JG=JG-1
  150     CONTINUE
          !
  160     NUM1=NUM1+LL
  170     CONTINUE
        ENDIF
      ELSE
*----
*  ++++ SPN THEORY ++++
*----
        DO 330 IL=0,MIN(ABS(ISCAT)-1,NANIS)
        FACT=REAL(2*IL+1)
        CALL LCMGET(KPMACR,'NJJS'//CAN(IL),NJJ(1))
        CALL LCMGET(KPMACR,'IJJS'//CAN(IL),IJJ(1))
        CALL LCMGET(KPMACR,'IPOS'//CAN(IL),IPOS(1))
        CALL LCMGET(KPMACR,'SCAT'//CAN(IL),XSCAT(1))
        NUM1=0
        DO 320 IR=1,NREG
        IBM=MATCOD(IR)
        IF(IBM.LE.0) GO TO 320
        IF(VOL(IR).EQ.0.0) GO TO 310
        IF(MOD(IL,2).EQ.0) THEN
          DO 255 I0=1,IELEM*IELEM
          IND=(IL/2)*L4+KN(NUM1+1)+I0-1
          JG=IJJ(IBM)
          DO 250 JND=1,NJJ(IBM)
          IF(JG.NE.IG) THEN
            SUNKNO(IND,IG)=SUNKNO(IND,IG)+FACT*FUNKNO(IND,JG)*
     >      VOL(IR)*XSCAT(IPOS(IBM)+JND-1)
          ENDIF
          JG=JG-1
  250     CONTINUE
  255     CONTINUE
        ELSE
          DO 305 I0=1,IELEM
          DO 275 IC=1,2
          IIC=1+(IC-1)*IELEM
          IND1=(IL/2)*L4+ABS(KN(NUM1+1+IC))+I0-1
          S1=REAL(SIGN(1,KN(NUM1+1+IC)))
          DO 270 JC=1,2
          JJC=1+(JC-1)*IELEM
          IND2=(IL/2)*L4+ABS(KN(NUM1+1+JC))+I0-1
          IF((KN(NUM1+1+IC).NE.0).AND.(KN(NUM1+1+JC).NE.0)) THEN
            S2=REAL(SIGN(1,KN(NUM1+1+JC)))
            AUXX=S1*S2*FACT*RR(IIC,JJC)*VOL(IR)
            JG=IJJ(IBM)
            DO 260 JND=1,NJJ(IBM)
            IF(JG.NE.IG) THEN
              SUNKNO(IND1,IG)=SUNKNO(IND1,IG)-AUXX*FUNKNO(IND2,JG)*
     1        XSCAT(IPOS(IBM)+JND-1)
            ENDIF
            JG=JG-1
  260       CONTINUE
          ENDIF
  270     CONTINUE
  275     CONTINUE
          DO 300 IC=3,4
          IIC=1+(IC-3)*IELEM
          IND1=(IL/2)*L4+ABS(KN(NUM1+1+IC))+I0-1
          S1=REAL(SIGN(1,KN(NUM1+1+IC)))
          DO 290 JC=3,4
          JJC=1+(JC-3)*IELEM
          IND2=(IL/2)*L4+ABS(KN(NUM1+1+JC))+I0-1
          IF((KN(NUM1+1+IC).NE.0).AND.(KN(NUM1+1+JC).NE.0)) THEN
            S2=REAL(SIGN(1,KN(NUM1+1+JC)))
            AUXX=S1*S2*FACT*RR(IIC,JJC)*VOL(IR)
            JG=IJJ(IBM)
            DO 280 JND=1,NJJ(IBM)
            IF(JG.NE.IG) THEN
              SUNKNO(IND1,IG)=SUNKNO(IND1,IG)-AUXX*FUNKNO(IND2,JG)*
     1        XSCAT(IPOS(IBM)+JND-1)
            ENDIF
            JG=JG-1
  280       CONTINUE
          ENDIF
  290     CONTINUE
  300     CONTINUE
  305     CONTINUE
        ENDIF
  310   NUM1=NUM1+5
  320   CONTINUE
  330   CONTINUE
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      IF((NLF.GT.0).OR.(IELEM.LT.0)) DEALLOCATE(RS,RR)
      DEALLOCATE(IDL,KN,DD,XX)
      DEALLOCATE(XSCAT)
      DEALLOCATE(IPOS,NJJ,IJJ)
      RETURN
      END
