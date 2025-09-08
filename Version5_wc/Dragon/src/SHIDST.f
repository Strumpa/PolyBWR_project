*DECK SHIDST
      SUBROUTINE SHIDST (IPSYS,NPSYS,IPTRK,IFTRAK,CDOOR,IMPX,NBM,NREG,
     1 NUN,NGRO,IPHASE,MAT,VOL,KEYFLX,LEAKSW,IRES,SIG0,SIG1,SIG2,TITR,
     2 FUNKNO,DILAV)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of escape probability information.
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
* IPSYS   pointer to the pij (L_PIJ signature).
* NPSYS   index array pointing to the IPSYS list component corresponding
*         to each energy group. Set to zero if a group is not to be
*         processed. Usually, NPSYS(I)=I.
* IPTRK   pointer to the tracking (L_TRACK signature).
* IFTRAK  unit number of the sequential binary tracking file.
* CDOOR   name of the geometry/solution module.
* IMPX    print flag (equal to zero for no print).
* NBM     number of mixtures.
* NREG    total number of merged blocks for which specific values
*         of the neutron flux and reactions rates are required.
* NUN     number of unknowns in the flux or source vector in one
*         energy group.
* NGRO    number of energy groups.
* IPHASE  type of flux solution (=1 use a native flux solution door;
*         =2 use collision probabilities).
* MAT     index-number of the mixture type assigned to each volume.
* VOL     volumes.
* KEYFLX  pointers of fluxes in unknown vector.
* LEAKSW  leakage flag (.TRUE. only if leakage is present on the outer
*         surface).
* IRES    resonant mixture number assigned to each mixture.
* SIG0    total macroscopic cross sections of the resonant materials
*         in each mixture.
* SIG1    total macroscopic cross sections of the light materials in
*         each mixture.
* SIG2    transport correction in each mixture.
* TITR    title.
*
*Parameters: output
* FUNKNO  information used for computing escape information for the
*         Nordheim method.
* DILAV   average dilution.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE DOORS_MOD
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPSYS,IPTRK
      CHARACTER CDOOR*12,TITR*72
      LOGICAL LEAKSW
      INTEGER NPSYS(NGRO),IFTRAK,IMPX,NBM,NREG,NUN,NGRO,IPHASE,
     1 MAT(NREG),KEYFLX(NREG),IRES(NBM)
      REAL VOL(NREG),SIG0(NBM,NGRO),SIG1(NBM,NGRO),SIG2(NBM,NGRO),
     1 FUNKNO(NUN,NGRO),DILAV(NGRO)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPSYS,KPSYS,IPMACR,IPSOU
      DOUBLE PRECISION TOT1,TOT2
      LOGICAL LNORM,LEXAC,REBFLG
      REAL, ALLOCATABLE, DIMENSION(:) :: SSIGT,SSIGW,SIGG
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SUN,FUN
      INTEGER NALBP
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(SSIGT(0:NBM),SSIGW(0:NBM),SIGG(0:NBM))
*----
*  INITIALIZATIONS.
*----
      NALBP=0
      ISTRM=1
      NANI=1
      NW=0
      IPIJK=1
      ITPIJ=1
      KNORM=1
      LNORM=.FALSE.
      IDIR=0
      LEXAC=.FALSE.
      JPSYS=LCMLID(IPSYS,'GROUP',NGRO)
*----
*  SELECT THE MACROSCOPIC CROSS SECTIONS.
*----
      SSIGT(0)=0.0
      SSIGW(0)=0.0
      DO 20 LLL=1,NGRO
      IF(NPSYS(LLL).NE.0) THEN
         DO 10 IBM=1,NBM
         SSIGT(IBM)=SIG0(IBM,LLL)+SIG1(IBM,LLL)-SIG2(IBM,LLL)
         SSIGW(IBM)=-SIG2(IBM,LLL)
   10    CONTINUE
         KPSYS=LCMDIL(JPSYS,LLL)
         CALL LCMPUT(KPSYS,'DRAGON-TXSC',NBM+1,2,SSIGT(0))
         CALL LCMPUT(KPSYS,'DRAGON-S0XSC',NBM+1,2,SSIGW(0))
      ENDIF
   20 CONTINUE
*----
*  ASSEMBLY MATRIX OR REDUCED COLLISION PROBABILITIES CALCULATION.
*----
      IF(IPHASE.EQ.1) THEN
*        USE A NATIVE DOOR.
         CALL DOORAV(CDOOR,JPSYS,NPSYS,IPTRK,IFTRAK,IMPX,NGRO,NREG,
     1   NBM,NANI,NW,MAT,VOL,KNORM,LEAKSW,TITR,NALBP,ISTRM)
      ELSE IF(IPHASE.EQ.2) THEN
*        USE A COLLISION PROBABILITY DOOR.
         CALL DOORPV(CDOOR,JPSYS,NPSYS,IPTRK,IFTRAK,IMPX,NGRO,NREG,
     1   NBM,NANI,MAT,VOL,KNORM,IPIJK,LEAKSW,ITPIJ,LNORM,TITR,NALBP)
      ENDIF
*----
*  ALLOCATE MEMORY.
*----
      ALLOCATE(SUN(NUN,NGRO),FUN(NUN,NGRO))
*----
*  SOLVE FOR THE FLUX AND SET UP VECTORS DILAV AND FUNKNO.
*----
      SUN(:NUN,:NGRO)=0.0
      DO 40 LLL=1,NGRO
      IF(NPSYS(LLL).NE.0) THEN
         SIGG(0:NBM)=0.0
         DO 30 IBM=1,NBM
         IF(IRES(IBM).GT.0) SIGG(IBM)=SIG0(IBM,LLL)
   30    CONTINUE
         CALL DOORS(CDOOR,IPTRK,NBM,0,NUN,SIGG,SUN(1,LLL))
      ENDIF
   40 CONTINUE
      CALL LCMLEN(IPSYS,'FLUX1',ILON1,ITYLCM)
      IF(ILON1.EQ.NUN*NGRO) THEN
         CALL LCMGET(IPSYS,'FLUX1',FUNKNO)
      ELSE
         FUNKNO(:NUN,:NGRO)=0.0
      ENDIF
      IPMACR=C_NULL_PTR
      IPSOU=C_NULL_PTR
      REBFLG=.FALSE.
      CALL DOORFV(CDOOR,JPSYS,NPSYS,IPTRK,IFTRAK,IMPX,NGRO,NBM,IDIR,
     1 NREG,NUN,IPHASE,LEXAC,MAT,VOL,KEYFLX,TITR,SUN,FUNKNO(1,1),IPMACR,
     2 IPSOU,REBFLG)
      CALL LCMPUT(IPSYS,'FLUX1',NUN*NGRO,2,FUNKNO)
*
      SUN(:NUN,:NGRO)=0.0
      DO 60 LLL=1,NGRO
      IF(NPSYS(LLL).NE.0) THEN
         SIGG(0:NBM)=0.0
         DO 50 IBM=1,NBM
         IF(IRES(IBM).GT.0) SIGG(IBM)=1.0
   50    CONTINUE
         CALL DOORS(CDOOR,IPTRK,NBM,0,NUN,SIGG,SUN(1,LLL))
      ENDIF
   60 CONTINUE
      CALL LCMLEN(IPSYS,'FLUX2',ILON2,ITYLCM)
      IF(ILON2.EQ.NUN*NGRO) THEN
         CALL LCMGET(IPSYS,'FLUX2',FUN)
      ELSE
         FUN(:NUN,:NGRO)=0.0
      ENDIF
      IPMACR=C_NULL_PTR
      IPSOU=C_NULL_PTR
      REBFLG=.FALSE.
      CALL DOORFV(CDOOR,JPSYS,NPSYS,IPTRK,IFTRAK,IMPX,NGRO,NBM,IDIR,
     1 NREG,NUN,IPHASE,LEXAC,MAT,VOL,KEYFLX,TITR,SUN,FUN(1,1),IPMACR,
     2 IPSOU,REBFLG)
      CALL LCMPUT(IPSYS,'FLUX2',NUN*NGRO,2,FUN)
      DO 80 LLL=1,NGRO
      IF(NPSYS(LLL).NE.0) THEN
         TOT1=0.0D0
         TOT2=0.0D0
         DO 70 I=1,NREG
         IBM=MAT(I)
         IF(IBM.EQ.0) GO TO 70
         IF(IRES(IBM).GT.0) THEN
            TOT2=TOT2+(1.0D0-FUNKNO(KEYFLX(I),LLL))*VOL(I)
            TOT1=TOT1+FUN(KEYFLX(I),LLL)*VOL(I)
         ENDIF
   70    CONTINUE
         DILAV(LLL)=REAL(TOT2/TOT1)
      ENDIF
   80 CONTINUE
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(SUN,FUN)
      DEALLOCATE(SIGG,SSIGW,SSIGT)
      RETURN
      END
