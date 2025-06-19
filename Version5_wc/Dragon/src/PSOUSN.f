*DECK PSOUSN
      SUBROUTINE PSOUSN(NUNF,NUNS,IG,IPTRK,JPTRK,KPMACR,NANIS,NREG,NMAT,
     > NGRP1,NGRP2,MATCOD,FLUX,SOURCE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the source from companion particle for the solution of SN 
* equations.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert 
*
*Parameters: input
* NUNF    first dimension of FLUX arrays.
* NUNS    first dimension of SOURCE arrays.
* IG      secondary group.
* IPTRK   pointer to the tracking LCM object (from main particle).
* JPTRK   pointer to the tracking LCM object (from companion particle).
* KPMACR  pointer to the secondary-group related macrolib information.
* NANIS   maximum cross section Legendre order.
* NREG    number of regions.
* NMAT    number of mixtures.
* NGRP1   number of primary energy groups.
* NGRP2   number of secondary energy groups.
* MATCOD  mixture indices.
* FLUX    fluxes.
*
*Parameters: output
* SOURCE  sources.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,JPTRK,KPMACR
      INTEGER NUNF,NUNS,IG,NANIS,NREG,NMAT,NGRP1,NGRP2,MATCOD(NREG)
      REAL FLUX(NUNF,NGRP1),SOURCE(NUNS,NGRP2)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40,PI4=12.5663706144)

      INTEGER IPAR(NSTATE),JPAR(NSTATE),P,P2,IELEM,EELEM,IELEM2,EELEM2,
     1 NM,NM2,EEL,IEL,IND,IND2,EL,EL2

      CHARACTER CANIL*2
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJ,NJJ,IPOS,MAP
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: M_INDEXES
      REAL, ALLOCATABLE, DIMENSION(:) :: XSCAT

      TYPE(C_PTR) IL_PTR,IM_PTR,IL2_PTR,IM2_PTR
      INTEGER, POINTER, DIMENSION(:) :: IL,IM,IL2,IM2

*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IJJ(0:NMAT),NJJ(0:NMAT),IPOS(0:NMAT))
      ALLOCATE(XSCAT(0:NMAT*NGRP1))
*----
*  RECOVER SNT SPECIFIC PARAMETERS.
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',IPAR)
      IF(IPAR(1).NE.NREG) CALL XABORT('PSOUSN: INCONSISTENT NREG.')
      ITYPE=IPAR(6)
      NSCT=IPAR(7)
      IELEM=IPAR(8)
      ISCAT=IPAR(16)
      EELEM=IPAR(35)
      CALL LCMGPD(IPTRK,'IL',IL_PTR)
      CALL LCMGPD(IPTRK,'IM',IM_PTR)
      CALL C_F_POINTER(IL_PTR,IL,(/ NSCT /))
      CALL C_F_POINTER(IM_PTR,IM,(/ NSCT /))
      CALL LCMGET(JPTRK,'STATE-VECTOR',JPAR)
      IF(JPAR(1).NE.NREG) CALL XABORT('PSOUSN: INCONSISTENT NREG.')
      ITYPE2=IPAR(6)
      NSCT2=JPAR(7)
      IELEM2=JPAR(8)
      ISCAT2=JPAR(16)
      EELEM2=JPAR(35)
      CALL LCMGPD(JPTRK,'IL',IL2_PTR)
      CALL LCMGPD(JPTRK,'IM',IM2_PTR)
      CALL C_F_POINTER(IL2_PTR,IL2,(/ NSCT2 /))
      CALL C_F_POINTER(IM2_PTR,IM2,(/ NSCT2 /))
      IF(ITYPE.NE.ITYPE2.OR.NSCT.NE.NSCT2.OR.ISCAT.NE.ISCAT2) 
     1 CALL XABORT('PSOUSN: INCONSISTENCE OF ANGULAR DISCRETISATION'
     2 //'BETWEEN THE PARTICLE AND ITS COMPANION PARTICLE.')
*----
*  MAPPING BETWEEN SPACE-ENERGY MOMENTS
*----
      NMX=IELEM
      NMX2=IELEM2
      IF((ITYPE.EQ.5).OR.(ITYPE.EQ.6).OR.(ITYPE.EQ.8)) THEN
            NMX=IELEM**2
            NMX2=IELEM2**2
      ELSEIF((ITYPE.EQ.7).OR.(ITYPE.EQ.9)) THEN
            NMX=IELEM**3
            NMX2=IELEM2**3
      ENDIF
      ALLOCATE(M_INDEXES(MAX(NMX,NMX2),MAX(EELEM,EELEM2)))
      ALLOCATE(MAP(NMX*EELEM))
      M_INDEXES=0
      MAP=0
      DO IEL=1,NMX2
      DO EEL=1,EELEM2
            M_INDEXES(IEL,EEL)=EELEM2*(IEL-1)+EEL
      ENDDO
      ENDDO
      DO IEL=1,NMX
      DO EEL=1,EELEM
      IF(IEL.LE.NMX2.AND.EEL.LE.EELEM2) THEN
            IND=EELEM*(IEL-1)+EEL
            MAP(IND)=M_INDEXES(IEL,EEL)       
      ENDIF
      ENDDO
      ENDDO
      DEALLOCATE(M_INDEXES)

*----
*  CONSTRUCT THE SOURCE.
*----
      IJJ(0)=0
      NJJ(0)=0
      IPOS(0)=0
      XSCAT(0)=0.0
      IOF0=0
      DO 100 P=1,NSCT
      ILP = IL(P)
      IF(ILP.GT.NANIS-1) GO TO 100     
      WRITE(CANIL,'(I2.2)') ILP

      CALL LCMGET(KPMACR,'NJJS'//CANIL,NJJ(1))
      CALL LCMGET(KPMACR,'IJJS'//CANIL,IJJ(1))
      CALL LCMGET(KPMACR,'IPOS'//CANIL,IPOS(1))
      CALL LCMGET(KPMACR,'SCAT'//CANIL,XSCAT(1))
      IF((ITYPE.EQ.2).OR.(ITYPE.EQ.4)) THEN
*----
*  SLAB OR SPHERICAL 1D CASE.
*----
         NSCT=ISCAT
         NM=IELEM*EELEM
         NM2=IELEM2*EELEM2
         DO 20 IR=1,NREG
         IBM=MATCOD(IR)
         IF(IBM.LE.0) GO TO 20
         DO 15 IEL=1,NM
         IF(MAP(IEL).EQ.0) CONTINUE
         IND=(IR-1)*NSCT*NM+NM*(P-1)+IEL
         IND2=(IR-1)*NSCT*NM2+NM2*(P-1)+MAP(IEL)
         JG=IJJ(IBM)
         DO 10 JND=1,NJJ(IBM)
         SOURCE(IND,IG)=SOURCE(IND,IG)+FLUX(IND2,JG)*

     >   XSCAT(IPOS(IBM)+JND-1)
         JG=JG-1
  10     CONTINUE
  15     CONTINUE
  20     CONTINUE
      ELSE IF(ITYPE.EQ.3) THEN
*----
*  CYLINDRICAL 1D CASE.
*----
         NSCT=(ISCAT/2)*(ISCAT/2+1)+(ISCAT+1)*MOD(ISCAT,2)/2
         DO 50 IM=0,IL
         IF(MOD(IL+IM,2).EQ.1) GO TO 50
         IOF0=IOF0+1
         DO 40 IR=1,NREG
         IBM=MATCOD(IR)
         IF(IBM.LE.0) GO TO 40
         IND=(IR-1)*NSCT+IOF0
         JG=IJJ(IBM)
         DO 30 JND=1,NJJ(IBM)
         SOURCE(IND,IG)=SOURCE(IND,IG)+FACT*FLUX(IND,JG)*
     >   XSCAT(IPOS(IBM)+JND-1)
         JG=JG-1
   30    CONTINUE
   40    CONTINUE
   50    CONTINUE
      ELSE IF((ITYPE.EQ.5).OR.(ITYPE.EQ.6).OR.(ITYPE.EQ.8)) THEN
*----
*  2D CASES (CARTESIAN OR R-Z).
*----
         NSCT=ISCAT*(ISCAT+1)/2
         NM=IELEM*IELEM*EELEM
         NM2=IELEM2*IELEM2*EELEM2
         DO 70 IR=1,NREG
         IBM=MATCOD(IR)
         IF(IBM.LE.0) GO TO 70
         DO 65 IEL=1,IELEM**2
         DO 64 EEL=1,EELEM
         IF(IEL.GT.IELEM2**2.OR.EEL.GT.EELEM2) CONTINUE
         EL=EELEM*(IEL-1)+EEL
         EL2=EELEM2*(IEL-1)+EEL
         IND=(IR-1)*NSCT*NM+NM*(P-1)+EL
         IND2=(IR-1)*NSCT*NM2+NM2*(P-1)+EL2
         JG=IJJ(IBM)
         DO 60 JND=1,NJJ(IBM)
         SOURCE(IND,IG)=SOURCE(IND,IG)+FLUX(IND2,JG)*
     >   XSCAT(IPOS(IBM)+JND-1)
         JG=JG-1
   60    CONTINUE
   64    CONTINUE
   65    CONTINUE
   70    CONTINUE
   80    CONTINUE
      ELSE IF((ITYPE.EQ.7).OR.(ITYPE.EQ.9)) THEN
*----
* 3D CARTESIAN CASE
*----
         NSCT=(ISCAT)**2
         NM=IELEM*IELEM*IELEM*EELEM
         NM2=IELEM2*IELEM2*IELEM2*EELEM2
         DO 90 IR=1,NREG
         IBM=MATCOD(IR)
         IF(IBM.LE.0) GO TO 90
         DO 85 IEL=1,IELEM**3
         DO 84 EEL=1,EELEM
         IF(IEL.GT.IELEM2**3.OR.EEL.GT.EELEM2) CONTINUE
         EL=EELEM*(IEL-1)+EEL
         EL2=EELEM2*(IEL-1)+EEL
         IND=(IR-1)*NSCT*NM+NM*(P-1)+EL
         IND2=(IR-1)*NSCT*NM2+NM2*(P-1)+EL2
         JG=IJJ(IBM)
         DO 80 JND=1,NJJ(IBM)
         SOURCE(IND,IG)=SOURCE(IND,IG)+FLUX(IND2,JG)*
     >   XSCAT(IPOS(IBM)+JND-1)      
         JG=JG-1
  80     CONTINUE
  84     CONTINUE
  85     CONTINUE
  90     CONTINUE

      ELSE
         CALL XABORT('PSOUSN: TYPE OF DISCRETIZATION NOT IMPLEMENTED.')
      ENDIF
 100  CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(XSCAT,MAP)
      DEALLOCATE(IPOS,NJJ,IJJ)
      RETURN
      END
