*DECK FLULPN
      SUBROUTINE FLULPN(IPMACR,NUNKNO,OPTION,TYPE,NGRP,NREG,NMAT,
     1 VOL,MATCOD,NMERG,IMERG,KEYFLX,FLUX,B2,IMPX,DIFHET,DHOM)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of heterogeneous leakage coefficients using the Todorova
* approximation.
*
*Copyright:
* Copyright (C) 2022 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPMACR  pointer to the macrolib LCM object (L_MACROLIB signature).
* NUNKNO  number of flux/current unknowns.
* OPTION  type of leakage coefficients; can be 'P0' (P-0), 'P1' (P-1),
*         'P0TR' (P-0 with transport correction).
* TYPE    type of buckling iteration.
*         Can be 'DIFF' (do a P0 calculation of DIFHET and exit);
*                 other (do another type of calculation).
* NGRP    number of groups.
* NREG    number of volumes.
* NMAT    number of material mixtures.
* VOL     volumes.
* MATCOD  mixture number of each volume.
* NMERG   number of leakage zones.
* IMERG   leakage zone index in each material mixture zone.
* KEYFLX  position of each flux in the unknown vector.
* FLUX    direct unknown vector.
* B2      buckling.
* IMPX    print flag.
*
*Parameters: output
* DIFHET  heterogeneous diffusion coefficients.
* DHOM    homogeneous diffusion coefficients.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER*4 OPTION,TYPE
      TYPE(C_PTR) IPMACR
      INTEGER NUNKNO,NGRP,NREG,NMAT,MATCOD(NREG),NMERG,IMERG(NMAT),
     1 KEYFLX(NREG),IMPX
      REAL VOL(NREG),FLUX(NUNKNO,NGRP),B2,DIFHET(NMERG,NGRP),DHOM(NGRP)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPMACR,KPMACR
      CHARACTER HSMG*131
      DOUBLE PRECISION B1GAMA,DDELN1,DDELN2,DDELD1,B2HOM,ST2,STR,GAMMA
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJ,NJJ,IPOS
      REAL, ALLOCATABLE, DIMENSION(:) :: WORK
      REAL, ALLOCATABLE, DIMENSION(:,:) :: ST,FLXIN
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: SCAT1
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: STOD
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(FLXIN(NMAT,NGRP),ST(NMAT,NGRP),SCAT1(NMAT,NGRP,NGRP))
*----
*  INITIALIZATION
*----
      IF((IMPX.GT.0).AND.(TYPE.EQ.'DIFF')) THEN
         WRITE (6,100)
      ELSE IF(IMPX.GT.0) THEN
         WRITE (6,110) OPTION,TYPE,B2
      ENDIF
      ST(:NMAT,:NGRP)=0.0
      SCAT1(:NMAT,:NGRP,:NGRP)=0.0
*----
*  RECOVER MACROSCOPIC CROSS SECTIONS
*----
      ALLOCATE(IJJ(NMAT),NJJ(NMAT),IPOS(NMAT),WORK(NMAT*NGRP))
      JPMACR=LCMGID(IPMACR,'GROUP')
      DO IGR=1,NGRP
        KPMACR=LCMGIL(JPMACR,IGR)
        CALL LCMGET(KPMACR,'NTOT0',ST(1,IGR))
        CALL LCMLEN(KPMACR,'SCAT01',ILONG,ITYLCM)
        IF((ILONG.NE.0).AND.(OPTION.NE.'P0').AND.(OPTION.NE.'B0')) THEN
          CALL LCMGET(KPMACR,'IJJS01',IJJ)
          CALL LCMGET(KPMACR,'NJJS01',NJJ)
          CALL LCMGET(KPMACR,'IPOS01',IPOS)
          CALL LCMGET(KPMACR,'SCAT01',WORK)
          DO IBM=1,NMAT
            IPO=IPOS(IBM)
            DO JGR=IJJ(IBM),IJJ(IBM)-NJJ(IBM)+1,-1
              SCAT1(IBM,IGR,JGR)=WORK(IPO) ! IGR <-- JGR
              IPO=IPO+1
            ENDDO
          ENDDO
        ENDIF
      ENDDO
      DEALLOCATE(WORK,IPOS,NJJ,IJJ)
*----
*  RECOVER INTEGRATED FLUX
*----
      FLXIN(:NMAT,:NGRP)=0.0
      DO IGR=1,NGRP
        DO IBM=1,NMAT
          DO I=1,NREG
            IND=KEYFLX(I)
            IF((MATCOD(I).EQ.IBM).AND.(IND.GT.0)) THEN
              FLXIN(IBM,IGR)=FLXIN(IBM,IGR)+FLUX(IND,IGR)*VOL(I)
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      IF((OPTION.EQ.'LKRD').OR.(OPTION.EQ.'RHS')) GO TO 10
*----
*  MAIN LOOP OVER LEAKAGE ZONES
*----
      B2HOM=DBLE(B2)
      GAMMA=1.0D0
      DO INM=1,NMERG
        IF((OPTION.EQ.'P0').OR.(OPTION.EQ.'B0')) THEN
*         P0 or B0 approximation
          DO IGR=1,NGRP
            DDELN1=0.D0
            DDELD1=0.D0
            DO IBM=1,NMAT
              IF(IMERG(IBM).EQ.INM) THEN
                DDELN1=DDELN1+ST(IBM,IGR)*FLXIN(IBM,IGR)
                DDELD1=DDELD1+FLXIN(IBM,IGR)
              ENDIF
            ENDDO
            ST2=DDELN1/DDELD1
            IF(OPTION.EQ.'B0') GAMMA=B1GAMA(2,B2HOM,ST2)
            DIFHET(INM,IGR)=REAL(1.0D0/(3.0D0*GAMMA*ST2))
          ENDDO
        ELSE IF((OPTION.EQ.'P0TR').OR.(OPTION.EQ.'B0TR').OR.
     1  (TYPE.EQ.'DIFF')) THEN
*         Outscatter approximation
          DO IGR=1,NGRP
            DDELN1=0.D0
            DDELN2=0.D0
            DDELD1=0.D0
            DO IBM=1,NMAT
              IF(IMERG(IBM).EQ.INM) THEN
                DDELN1=DDELN1+ST(IBM,IGR)*FLXIN(IBM,IGR)
                DO JGR=1,NGRP
                  DDELN2=DDELN2+SCAT1(IBM,JGR,IGR)*FLXIN(IBM,IGR)
                ENDDO
                DDELD1=DDELD1+FLXIN(IBM,IGR)
              ENDIF
            ENDDO
            ST2=DDELN1/DDELD1
            IF(OPTION.EQ.'B0TR') GAMMA=B1GAMA(2,B2HOM,ST2)
            STR=(GAMMA*DDELN1-DDELN2)/DDELD1
            DIFHET(INM,IGR)=REAL(1.0D0/(3.0D0*STR))
          ENDDO
        ELSE IF((OPTION.EQ.'P1').OR.(OPTION.EQ.'B1')) THEN
*         Inscatter approximation
          ALLOCATE(STOD(NGRP,NGRP+1))
          STOD(:NGRP,:NGRP+1)=0.0D0
          DO IGR=1,NGRP
            IF(OPTION.EQ.'B1') THEN
              DDELN1=0.D0
              DDELD1=0.D0
              DO IBM=1,NMAT
                IF(IMERG(IBM).EQ.INM) THEN
                  DDELN1=DDELN1+ST(IBM,IGR)*FLXIN(IBM,IGR)
                  DDELD1=DDELD1+FLXIN(IBM,IGR)
                ENDIF
              ENDDO
              ST2=DDELN1/DDELD1
              GAMMA=B1GAMA(2,B2HOM,ST2)
            ENDIF
            DO IBM=1,NMAT
              IF(IMERG(IBM).EQ.INM) THEN
                STOD(IGR,IGR)=STOD(IGR,IGR)+GAMMA*ST(IBM,IGR)*
     1          FLXIN(IBM,IGR)
                DO JGR=1,NGRP
                  STOD(IGR,JGR)=STOD(IGR,JGR)-SCAT1(IBM,IGR,JGR)*
     1            FLXIN(IBM,JGR)
                ENDDO
                STOD(IGR,NGRP+1)=STOD(IGR,NGRP+1)+FLXIN(IBM,IGR)/3.0D0
              ENDIF
            ENDDO
          ENDDO
          CALL ALSBD(NGRP,1,STOD,IER,NGRP)
          IF(IER.NE.0) CALL XABORT('FLULPN: SINGULAR MATRIX.')
          DO IGR=1,NGRP
            DIFHET(INM,IGR)=REAL(STOD(IGR,NGRP+1))
          ENDDO
          DEALLOCATE(STOD)
        ELSE
          WRITE(HSMG,'(15HFLULPN: OPTION ,A,23H IS INVALID WITH TODORO,
     1    17HVA APPROXIMATION.)') OPTION
          CALL XABORT(HSMG)
        ENDIF
      ENDDO
*----
*  COMPUTE THE HOMOGENEOUS LEAKAGE COEFFICIENTS
*----
   10 DO IGR=1,NGRP
        DHOM(IGR)=0.0
        FLTOT=0.0
        DO IBM=1,NMAT
          INM=IMERG(IBM)
          DHOM(IGR)=DHOM(IGR)+FLXIN(IBM,IGR)*DIFHET(INM,IGR)
          FLTOT=FLTOT+FLXIN(IBM,IGR)
        ENDDO
        DHOM(IGR)=DHOM(IGR)/FLTOT
      ENDDO
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(SCAT1,ST,FLXIN)
      RETURN
*
  100 FORMAT(/54H FLULPN: OUTSCATTER DIFFUSION COEFFICIENT CALCULATION.)
  110 FORMAT(/21H FLULPN: SOLUTION OF ,A4,21H EQUATIONS WITH TYPE ,A4,
     1 10H BUCKLING=,1P,E12.4)
      END
