*DECK TRIVSO
      SUBROUTINE TRIVSO(MAX1,IG,IPTRK,KPMACR,NANIS,NREG,NMAT,NUNKNO,
     > NGRP,MATCOD,VOL,FUNKNO,SUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the source for the solution of diffusion or PN equations.
* TRIVAC-specific version.
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
      INTEGER JPAR(NSTATE)
      CHARACTER CAN(0:9)*2
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJ,NJJ,IPOS,KN
      REAL, ALLOCATABLE, DIMENSION(:) :: XSCAT
      REAL, ALLOCATABLE, DIMENSION(:,:) :: RR
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
*  RECOVER TRIVAC SPECIFIC PARAMETERS.
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',JPAR)
      IF(JPAR(1).NE.NREG) CALL XABORT('TRIVSO: INCONSISTENT NREG.')
      IF(JPAR(2).NE.NUNKNO) CALL XABORT('TRIVSO: INCONSISTENT NUNKNO.')
      ITYPE=JPAR(6)
      IELEM=JPAR(9)
      ICOL=JPAR(10)
      L4=JPAR(11)
      NLF=JPAR(30)
      ISCAT=JPAR(32)
      IF(ICOL.EQ.4) THEN
        CALL XABORT('TRIVSO: COLLOCATION NODAL NOT IMPLEMENTED.')
      ELSE IF((ITYPE.NE.2).AND.(ITYPE.NE.5).AND.(ITYPE.NE.7)) THEN
        CALL XABORT('TRIVSO: CARTESIAN 1D, 2D OR 3D GEOMETRY EXPECTED.')
      ELSE IF(IELEM.LT.0) THEN
        CALL XABORT('TRIVSO: RAVIART-THOMAS METHOD EXPECTED.')
      ENDIF
      CALL LCMLEN(IPTRK,'KN',MAXKN,ITYLCM)
      ALLOCATE(KN(MAXKN))
      CALL LCMGET(IPTRK,'KN',KN)
*----
*  RECOVER THE FINITE ELEMENT UNIT STIFFNESS MATRIX.
*----
      IF(NLF.GT.0) THEN
         CALL LCMSIX(IPTRK,'BIVCOL',1)
         CALL LCMLEN(IPTRK,'T',LC,ITYLCM)
         ALLOCATE(RR(LC,LC))
         CALL LCMGET(IPTRK,'R',RR)
         CALL LCMSIX(IPTRK,' ',2)
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
         IF((ITYPE.EQ.2).OR.(ITYPE.EQ.5).OR.(ITYPE.EQ.7)) THEN
*----
*  CARTESIAN DUAL (RAVIART-THOMAS) CASE.
*----
            NUM1=0
            DO 60 IR=1,NREG
            IBM=MATCOD(IR)
            IF(IBM.LE.0) GO TO 60
            IF(VOL(IR).EQ.0.0) GO TO 50
            DO 40 I0=1,IELEM**3
            IND=KN(NUM1+1)+I0-1
            JG=IJJ(IBM)
            DO 30 JND=1,NJJ(IBM)
            IF(JG.NE.IG) THEN
              SUNKNO(IND,IG)=SUNKNO(IND,IG)+FUNKNO(IND,JG)*VOL(IR)*
     >        XSCAT(IPOS(IBM)+JND-1)
            ENDIF
            JG=JG-1
   30       CONTINUE
   40       CONTINUE
   50       NUM1=NUM1+1+6*IELEM**2
   60       CONTINUE
         ELSE
            CALL XABORT('TRIVSO: DISCRETIZATION NOT IMPLEMENTED(1).')
         ENDIF
      ELSE
*----
*  ++++ SPN THEORY ++++
*----
         DO 350 IL=0,MIN(ABS(ISCAT)-1,NANIS)
         FACT=REAL(2*IL+1)
         CALL LCMGET(KPMACR,'NJJS'//CAN(IL),NJJ(1))
         CALL LCMGET(KPMACR,'IJJS'//CAN(IL),IJJ(1))
         CALL LCMGET(KPMACR,'IPOS'//CAN(IL),IPOS(1))
         CALL LCMGET(KPMACR,'SCAT'//CAN(IL),XSCAT(1))
         NUM1=0
         DO 340 IR=1,NREG
         IBM=MATCOD(IR)
         IF(IBM.LE.0) GO TO 340
         IF(MOD(IL,2).EQ.0) THEN
           DO 250 I0=1,IELEM**3
           IND=(IL/2)*L4+KN(NUM1+1)+I0-1
           JG=IJJ(IBM)
           DO 240 JND=1,NJJ(IBM)
           IF(JG.NE.IG) THEN
             SUNKNO(IND,IG)=SUNKNO(IND,IG)+FACT*FUNKNO(IND,JG)*
     >       VOL(IR)*XSCAT(IPOS(IBM)+JND-1)
           ENDIF
           JG=JG-1
  240      CONTINUE
  250      CONTINUE
         ELSE
           DO 330 I0=1,IELEM
           DO 275 IC=1,2
           IIC=1+(IC-1)*IELEM
           KN1=KN(NUM1+2+(IC-1)*IELEM**2)
           IND1=(IL/2)*L4+ABS(KN1)+I0-1
           S1=REAL(SIGN(1,KN1))
           DO 270 JC=1,2
           JJC=1+(JC-1)*IELEM
           KN2=KN(NUM1+2+(JC-1)*IELEM**2)
           IND2=(IL/2)*L4+ABS(KN2)+I0-1
           IF((KN1.NE.0).AND.(KN2.NE.0)) THEN
             S2=REAL(SIGN(1,KN2))
             AUXX=S1*S2*FACT*RR(IIC,JJC)*VOL(IR)
             JG=IJJ(IBM)
             DO 260 JND=1,NJJ(IBM)
             IF(JG.NE.IG) THEN
               SUNKNO(IND1,IG)=SUNKNO(IND1,IG)-AUXX*FUNKNO(IND2,JG)*
     1         XSCAT(IPOS(IBM)+JND-1)
             ENDIF
             JG=JG-1
  260        CONTINUE
           ENDIF
  270      CONTINUE
  275      CONTINUE
           DO 295 IC=3,4
           IIC=1+(IC-3)*IELEM
           KN1=KN(NUM1+2+(IC-1)*IELEM**2)
           IND1=(IL/2)*L4+ABS(KN1)+I0-1
           S1=REAL(SIGN(1,KN1))
           DO 290 JC=3,4
           JJC=1+(JC-3)*IELEM
           KN2=KN(NUM1+2+(JC-1)*IELEM**2)
           IND2=(IL/2)*L4+ABS(KN2)+I0-1
           IF((KN1.NE.0).AND.(KN2.NE.0)) THEN
             S2=REAL(SIGN(1,KN2))
             AUXX=S1*S2*FACT*RR(IIC,JJC)*VOL(IR)
             JG=IJJ(IBM)
             DO 280 JND=1,NJJ(IBM)
             IF(JG.NE.IG) THEN
               SUNKNO(IND1,IG)=SUNKNO(IND1,IG)-AUXX*FUNKNO(IND2,JG)*
     1         XSCAT(IPOS(IBM)+JND-1)
             ENDIF
             JG=JG-1
  280        CONTINUE
           ENDIF
  290      CONTINUE
  295      CONTINUE
           DO 320 IC=5,6
           IIC=1+(IC-5)*IELEM
           KN1=KN(NUM1+2+(IC-1)*IELEM**2)
           IND1=(IL/2)*L4+ABS(KN1)+I0-1
           S1=REAL(SIGN(1,KN1))
           DO 310 JC=5,6
           JJC=1+(JC-5)*IELEM
           KN2=KN(NUM1+2+(JC-1)*IELEM**2)
           IND2=(IL/2)*L4+ABS(KN2)+I0-1
           IF((KN1.NE.0).AND.(KN2.NE.0)) THEN
             S2=REAL(SIGN(1,KN2))
             AUXX=S1*S2*FACT*RR(IIC,JJC)*VOL(IR)
             JG=IJJ(IBM)
             DO 300 JND=1,NJJ(IBM)
             IF(JG.NE.IG) THEN
               SUNKNO(IND1,IG)=SUNKNO(IND1,IG)-AUXX*FUNKNO(IND2,JG)*
     1         XSCAT(IPOS(IBM)+JND-1)
             ENDIF
             JG=JG-1
  300        CONTINUE
           ENDIF
  310      CONTINUE
  320      CONTINUE
  330      CONTINUE
         ENDIF
         NUM1=NUM1+1+6*IELEM**2
  340    CONTINUE
  350    CONTINUE
      ENDIF
      DEALLOCATE(KN)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      IF(NLF.GT.0) DEALLOCATE(RR)
      DEALLOCATE(XSCAT)
      DEALLOCATE(IPOS,NJJ,IJJ)
      RETURN
      END
