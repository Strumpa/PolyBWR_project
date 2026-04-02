*DECK TONDST
      SUBROUTINE TONDST (IPSYS,NPSYS,IPTRK,IFTRAK,CDOOR,IMPX,NBM,NBNRS,
     1 NREG,NUN,NGRO,IPHASE,MAT,VOL,KEYFLX,LEAKSW,IRES,DENM,SIGT0,SIGT2,
     2 SIGT3,TCEPS,TITR,DILAV,TK3,TK4)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of escape probability information.
*
*Copyright:
* Copyright (C) 2017 Ecole Polytechnique de Montreal
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
* IPTRK   pointer to the tracking. (L_TRACK signature).
* IFTRAK  unit number of the sequential binary tracking file.
* CDOOR   name of the geometry/solution module.
* IMPX    print flag (equal to zero for no print).
* NBM     number of mixtures.
* NBNRS   number of totaly correlated resonant regions.
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
* LEAKSW  leakage flag (=.TRUE. if leakage is present on the outer
*         surface).
* IRES    resonant mixture number assigned to each mixture.
* DENM    number density of the resonant isotope in each mixture.
* SIGT0   total macroscopic cross sections of the resonant isotope
*         in each mixture.
* SIGT2   total macroscopic cross sections of the light materials in
*         each mixture.
* SIGT3   transport correction in each mixture.
* TCEPS   relative tolerance for skipping DOORPV/AV.
* TITR    title.
*
*Parameters: output
* DILAV   average dilution.
*
*Parameters: input/output
* TK3     cpu time to compute system matrices.
* TK4     cpu time to compute fluxes.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE DOORS_MOD
      USE TONDST_CACHE_MOD
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPSYS,IPTRK
      CHARACTER CDOOR*12,TITR*72
      LOGICAL LEAKSW
      INTEGER NPSYS(NGRO),IFTRAK,IMPX,NBM,NBNRS,NREG,NUN,NGRO,IPHASE,
     1 MAT(NREG),KEYFLX(NREG),IRES(NBM)
      REAL VOL(NREG),DENM(0:NBM),SIGT0(0:NBM,NGRO),SIGT2(0:NBM,NGRO),
     1 SIGT3(0:NBM,NGRO),TCEPS,DILAV(NBNRS,NGRO),TK3,TK4
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPSYS,KPSYS,IPMACR,IPSOU
      LOGICAL LNORM,LEXAC,REBFLG
      REAL, ALLOCATABLE, DIMENSION(:) :: SSIGT,SSIGW,SSIGT_OLD
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SUN,FUN1,FUN2
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: TOT1,TOT2
      INTEGER NALBP
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NPSYS_PIJ
      REAL SMAX,SNORM
      INTEGER NSKIP
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(SSIGT(0:NBM),SSIGW(0:NBM))
*----
*  INITIALIZATIONS.
*----
      NALBP=0
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
*  PIJ/ARM INFORMATION CACHING
*----
      ALLOCATE(NPSYS_PIJ(NGRO),SSIGT_OLD(0:NBM))
*     Initialize cache if needed.
      CALL TC_INIT(NBM,NGRO)
*     Build filtered NPSYS_PIJ: skip groups with small XS change.
      NPSYS_PIJ(:NGRO)=NPSYS(:NGRO)
      NSKIP=0
*----
*  SELECT THE MACROSCOPIC CROSS SECTIONS.
*----
      SSIGT(0)=0.0
      SSIGW(0)=0.0
      DO 20 LLL=1,NGRO
      IF(NPSYS(LLL).NE.0) THEN
        DO IBM=1,NBM
          SSIGT(IBM)=SIGT0(IBM,LLL)+SIGT2(IBM,LLL)-SIGT3(IBM,LLL)
          SSIGW(IBM)=-SIGT3(IBM,LLL)
        ENDDO
        CALL LCMLEL(JPSYS,LLL,ILONG,ITYLCM)
        IF(ILONG.NE.0) THEN
          KPSYS=LCMGIL(JPSYS,LLL)
          CALL LCMGET(KPSYS,'DRAGON-TXSC',SSIGT_OLD(0))
          TC_TOTAL=TC_TOTAL+1
          IF(TC_VALID(LLL)) THEN
            SMAX=0.0
            SNORM=0.0
            DO IBM=1,NBM
              SMAX=MAX(SMAX,ABS(SSIGT(IBM)-SSIGT_OLD(IBM)))
              SNORM=MAX(SNORM,ABS(SSIGT(IBM)))
            ENDDO
            IF((SNORM.GT.0.0).AND.(SMAX/SNORM.LT.TCEPS)) THEN
              NPSYS_PIJ(LLL)=0
              TC_HITS=TC_HITS+1
              NSKIP=NSKIP+1
            ENDIF
          ENDIF
        ELSE
          KPSYS=LCMDIL(JPSYS,LLL)
        ENDIF
        IF(NPSYS_PIJ(LLL).NE.0) THEN
          CALL LCMPUT(KPSYS,'DRAGON-TXSC',NBM+1,2,SSIGT(0))
          CALL LCMPUT(KPSYS,'DRAGON-S0XSC',NBM+1,2,SSIGW(0))
        ENDIF
      ENDIF
   20 CONTINUE
      IF(IMPX.GE.5) THEN
        WRITE(6,'(33H TONDST: PIJ CACHE SKIP/TOTAL GRP,
     1  2H =,I4,1H/,I4)') NSKIP,NGRO
      ENDIF
*----
*  ASSEMBLY MATRIX OR REDUCED COLLISION PROBABILITIES CALCULATION.
*----
      CALL KDRCPU(TKA)
      ISTRM=1
      IF(IPHASE.EQ.1) THEN
*        USE A NATIVE DOOR.
         CALL DOORAV(CDOOR,JPSYS,NPSYS_PIJ,IPTRK,IFTRAK,IMPX,NGRO,
     1   NREG,NBM,NANI,NW,MAT,VOL,KNORM,LEAKSW,TITR,NALBP,ISTRM)
      ELSE IF(IPHASE.EQ.2) THEN
*        USE A COLLISION PROBABILITY DOOR.
         CALL DOORPV(CDOOR,JPSYS,NPSYS_PIJ,IPTRK,IFTRAK,IMPX,NGRO,
     1   NREG,NBM,NANI,MAT,VOL,KNORM,IPIJK,LEAKSW,ITPIJ,LNORM,TITR,
     2   NALBP)    
      ENDIF
*----
*  UPDATE CACHE FOR GROUPS THAT WERE RECOMPUTED.
*----
      DO LLL=1,NGRO
        IF(NPSYS_PIJ(LLL).NE.0) THEN
          TC_VALID(LLL)=.TRUE.
        ENDIF
      ENDDO
      DEALLOCATE(SSIGT_OLD,NPSYS_PIJ)
      CALL KDRCPU(TKB)
      TK3=TK3+(TKB-TKA)
*----
*  ALLOCATE MEMORY.
*----
      ALLOCATE(SUN(NUN,NGRO),FUN1(NUN,NGRO),FUN2(NUN,NGRO))
*----
*  SOLVE FOR THE FLUX AND SET UP VECTOR DILAV.
*----
      CALL KDRCPU(TKA)
      SUN(:NUN,:NGRO)=0.0
      DO 30 LLL=1,NGRO
      IF(NPSYS(LLL).NE.0) THEN
         CALL DOORS(CDOOR,IPTRK,NBM,0,NUN,SIGT2(0,LLL),SUN(1,LLL))
      ENDIF
   30 CONTINUE
      CALL LCMLEN(IPSYS,'FLUX1',ILON1,ITYLCM)
      IF(ILON1.EQ.NUN*NGRO) THEN
         CALL LCMGET(IPSYS,'FLUX1',FUN1)
      ELSE
         FUN1(:NUN,:NGRO)=0.0
      ENDIF
      IPMACR=C_NULL_PTR
      IPSOU=C_NULL_PTR
      REBFLG=.FALSE.
      CALL DOORFV(CDOOR,JPSYS,NPSYS,IPTRK,IFTRAK,IMPX,NGRO,NBM,IDIR,
     1 NREG,NUN,IPHASE,LEXAC,MAT,VOL,KEYFLX,TITR,SUN,FUN1,IPMACR,
     2 IPSOU,REBFLG)
      CALL LCMPUT(IPSYS,'FLUX1',NUN*NGRO,2,FUN1)
*
      SUN(:NUN,:NGRO)=0.0
      DO 40 LLL=1,NGRO
      IF(NPSYS(LLL).NE.0) THEN
         CALL DOORS(CDOOR,IPTRK,NBM,0,NUN,DENM,SUN(1,LLL))
      ENDIF
   40 CONTINUE
      CALL LCMLEN(IPSYS,'FLUX2',ILON2,ITYLCM)
      IF(ILON2.EQ.NUN*NGRO) THEN
         CALL LCMGET(IPSYS,'FLUX2',FUN2)
      ELSE
         FUN2(:NUN,:NGRO)=0.0
      ENDIF
      IPMACR=C_NULL_PTR
      REBFLG=.FALSE.
      CALL DOORFV(CDOOR,JPSYS,NPSYS,IPTRK,IFTRAK,IMPX,NGRO,NBM,IDIR,
     1 NREG,NUN,IPHASE,LEXAC,MAT,VOL,KEYFLX,TITR,SUN,FUN2,IPMACR,
     2 IPSOU,REBFLG)
      CALL LCMPUT(IPSYS,'FLUX2',NUN*NGRO,2,FUN2)
      ALLOCATE(TOT2(NBNRS),TOT1(NBNRS))
      DO 70 LLL=1,NGRO
      IF(NPSYS(LLL).NE.0) THEN
         TOT2(:)=0.0D0
         TOT1(:)=0.0D0
         DO 50 I=1,NREG
         IBM=MAT(I)
         IF(IBM.EQ.0) GO TO 50
         IRS=IRES(IBM)
         IF(IRS.GT.0) THEN
            TOT1(IRS)=TOT1(IRS)+FUN1(KEYFLX(I),LLL)*VOL(I)
            TOT2(IRS)=TOT2(IRS)+FUN2(KEYFLX(I),LLL)*VOL(I)
         ENDIF
   50    CONTINUE
         DO 60 IRS=1,NBNRS
         DILAV(IRS,LLL)=REAL(TOT1(IRS)/TOT2(IRS))
   60    CONTINUE
      ENDIF
   70 CONTINUE
      DEALLOCATE(TOT2,TOT1)
      CALL KDRCPU(TKB)
      TK4=TK4+(TKB-TKA)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(SUN,FUN2,FUN1)
      DEALLOCATE(SSIGW,SSIGT)
      RETURN
      END
