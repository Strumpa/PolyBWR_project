*DECK SNSBFP
      SUBROUTINE SNSBFP(IG,IPTRK,KPMACR,KPSYS,NANIS,NLF,NREG,NMAT,
     1 NUNKNO,NGRP,MATCOD,FLUX,QEXT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the QEXT for the solution of SN equations with a Boltzmann-
* Fokker-Planck discretization.
*
*Copyright:
* Copyright (C) 2020 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert 
*
*Parameters: input
* IG      secondary group.
* IPTRK   pointer to the tracking LCM object.
* KPMACR  pointer to the secondary-group related macrolib information.
* KPSYS   pointer to the system matrix information.
* NANIS   maximum cross section Legendre order.
* NLF     number of Legendre components in the flux.
* NREG    number of regions.
* NMAT    number of mixtures.
* NUNKNO  number of unknowns per energy group including spherical
*         harmonic terms, interface currents, fundamental currents
*         and slowing-down angular fluxes at group boundary.
* NGRP    number of energy groups.
* MATCOD  mixture indices.
* FLUX    fluxes and slowing-down angular fluxes at group boundary.
*
*Parameters: output
* QEXT    sources and slowing-down angular fluxes at group boundary.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK,KPMACR,KPSYS
      INTEGER IG,NANIS,NLF,NREG,NMAT,NUNKNO,NGRP,MATCOD(NREG)
      REAL FLUX(NUNKNO,NGRP),QEXT(NUNKNO,NGRP)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40,PI4=12.5663706144)
      INTEGER JPAR(NSTATE),EELEM,P
      CHARACTER CAN(0:19)*2
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJ,NJJ,IPOS
      REAL, ALLOCATABLE, DIMENSION(:) :: XSCAT
      TYPE(C_PTR) IL_PTR,IM_PTR
      INTEGER, POINTER, DIMENSION(:) :: IL,IM
*----
*  DATA STATEMENTS
*----
      DATA CAN /'00','01','02','03','04','05','06','07','08','09',
     >          '10','11','12','13','14','15','16','17','18','19'/
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IJJ(0:NMAT),NJJ(0:NMAT),IPOS(0:NMAT))
      ALLOCATE(XSCAT(0:NMAT*NGRP))
*----
*  RECOVER SNT SPECIFIC PARAMETERS
*----
      MAX2=0
      CALL LCMGET(IPTRK,'STATE-VECTOR',JPAR)
      IF(JPAR(1).NE.NREG) CALL XABORT('SNSBFP: INCONSISTENT NREG.')
      IF(JPAR(2).NE.NUNKNO) CALL XABORT('SNSBFP: INCONSISTENT NUNKNO.')
      IF(JPAR(15).NE.NLF) CALL XABORT('SNSBFP: INCONSISTENT NLF.')
      ITYPE=JPAR(6)
      NSCT=JPAR(7)
      IELEM=JPAR(8)
      ISCAT=JPAR(16)
      EELEM=JPAR(35)
      CALL LCMGPD(IPTRK,'IL',IL_PTR)
      CALL LCMGPD(IPTRK,'IM',IM_PTR)
      CALL C_F_POINTER(IL_PTR,IL,(/ NSCT /))
      CALL C_F_POINTER(IM_PTR,IM,(/ NSCT /))
*----
*  CONSTRUCT THE QEXT.
*----
      IJJ(0)=0
      NJJ(0)=0
      IPOS(0)=0
      XSCAT(0)=0.0
      IOF0=0
      DO 130 P=1,NSCT
      ILP=IL(P)
      IF(ILP.GT.MIN(ISCAT-1,NANIS)) GO TO 130
      CALL LCMGET(KPMACR,'NJJS'//CAN(ILP),NJJ(1))
      CALL LCMGET(KPMACR,'IJJS'//CAN(ILP),IJJ(1))
      CALL LCMGET(KPMACR,'IPOS'//CAN(ILP),IPOS(1))
      CALL LCMGET(KPMACR,'SCAT'//CAN(ILP),XSCAT(1))
      IF((ITYPE.EQ.2).OR.(ITYPE.EQ.4)) THEN
*----
*  SLAB OR SPHERICAL 1D CASE.
*----
         NM=IELEM*EELEM
         MAX2=IELEM*NLF*NREG
         DO 20 IR=1,NREG
         IBM=MATCOD(IR)
         IF(IBM.LE.0) GO TO 20
         DO 15 IEL=1,NM
         IND=(IR-1)*NSCT*NM+(P-1)*NM+IEL
         JG=IJJ(IBM)
         DO 10 JND=1,NJJ(IBM)
         IF(JG.NE.IG) THEN
            QEXT(IND,IG)=QEXT(IND,IG)+FLUX(IND,JG)*
     1      XSCAT(IPOS(IBM)+JND-1)
         ENDIF
         JG=JG-1
   10    CONTINUE
   15    CONTINUE
   20    CONTINUE
      ELSE IF((ITYPE.EQ.5).OR.(ITYPE.EQ.6).OR.(ITYPE.EQ.8)) THEN
*----
*  2D CASES (CARTESIAN OR R-Z).
*----
        NME=IELEM*IELEM
        NM=NME*EELEM
        CALL LCMLEN(IPTRK,'DU',NPQ,ITYLCM)
        MAX2=NME*NPQ*NREG
        DO 70 IR=1,NREG
        IBM=MATCOD(IR)
        IF(IBM.LE.0) GO TO 70
        DO 65 IEL=1,NM
        IND=(IR-1)*NSCT*NM+(P-1)*NM+IEL
        JG=IJJ(IBM)
        DO 60 JND=1,NJJ(IBM)
        IF(JG.NE.IG) THEN
          QEXT(IND,IG)=QEXT(IND,IG)+FLUX(IND,JG)*
     1    XSCAT(IPOS(IBM)+JND-1)
        ENDIF
        JG=JG-1
   60   CONTINUE
   65   CONTINUE
   70   CONTINUE
      ELSE IF(ITYPE.EQ.7) THEN
*----
*  3D CASES (CARTESIAN)
*----
        NME=IELEM*IELEM*IELEM
        NM=NME*EELEM
        CALL LCMLEN(IPTRK,'DU',NPQ,ITYLCM)
        MAX2=NME*NPQ*NREG
        DO 110 IR=1,NREG
        IBM=MATCOD(IR)
        IF(IBM.LE.0) GO TO 110
        DO 100 IEL=1,NM
        IND=(IR-1)*NSCT*NM+(P-1)*NM+IEL
        JG=IJJ(IBM)
        DO 90 JND=1,NJJ(IBM)
        IF(JG.NE.IG) THEN
          QEXT(IND,IG)=QEXT(IND,IG)+FLUX(IND,JG)*
     1    XSCAT(IPOS(IBM)+JND-1)
        ENDIF
        JG=JG-1
   90   CONTINUE
  100   CONTINUE
  110   CONTINUE
      ELSE
        CALL XABORT('SNSBFP: TYPE OF DISCRETIZATION NOT IMPLEMENTED.')
      ENDIF
  130 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(XSCAT)
      DEALLOCATE(IPOS,NJJ,IJJ)
*----
*  RECOVER SLOWING-DOWN ANGULAR FLUXES
*----
      MAX1=NUNKNO-MAX2
      CALL LCMGET(KPSYS,'DRAGON-DELTE',DELTAE)
      IF(IG.EQ.1) THEN
        QEXT(MAX1+1:MAX1+MAX2,1)=0.0
      ELSE
        QEXT(MAX1+1:MAX1+MAX2,IG)=FLUX(MAX1+1:MAX1+MAX2,IG-1)*DELTAE
      ENDIF
      RETURN
      END
