*DECK PIJWIJ
      SUBROUTINE PIJWIJ(  IPTRK,   IPRT,  NSOUT,   NREG,  NBMIX,   NANI,
     >                   MATCOD, VOLUME, XSSIGT, XSSIGW, NELPIJ,  IPIJK,
     >                   LEAKSW,  N2PRO,   NSBG,  NPSYS,   NPST,  NALBP,
     >                     ALBP, MATALB, VOLSUR,  DPROB, DPROBX,    PIJ,
     >                   PROBKS )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of the scattering-reduced collision probabilities for
* EXCELL. All surfaces will disappear from the system using external
* boundary conditions.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* IPTRK   pointer to the tracking (L_TRACK signature).
* IPRT    print flag (equal to zero for no print).
* NSOUT   number of surfaces.
* NREG    total number of merged blocks for which specific values
*         of the neutron flux and reactions rates are required.
* NBMIX   number of mixtures (NBMIX=max(MATCOD(i))).
* NANI    number of Legendre orders.
* MATCOD  index number of the mixture type assigned to each volume.
* VOLUME  volumes.
* XSSIGT  total macroscopic cross sections ordered by mixture.
* XSSIGW  P0 within-group scattering macroscopic cross sections
*         ordered by mixture.
* NELPIJ  number of elements in symmetrized pij matrix.
* IPIJK   pij option (=1 pij, =4 pijk).
* LEAKSW  leakage flag (=.true. if neutron leakage through external
*         boundary is present).
* N2PRO   number of terms in collision probability matrices, including
*         surface and volume contributions.
* NSBG    number of energy groups.
* NPSYS   non-converged energy group indices.
* NPST    first dimension of matrix PROBKS.
* NALBP   number of multigroup physical albedos.
* ALBP    multigroup physical albedos.
* MATALB  global mixture/albedo identification vector.
* VOLSUR  global surface volume vector.
* DPROB   collision probabilities from EXCELP.
* DPROBX  directional collision probabilities from EXCELP.
*
*Parameters: output
* PIJ     reduced and symmetrized collision probabilities.
* PROBKS  directional collision probabilities.
*
*-----------------------------------------------------------------------
*--------+---------------- R O U T I N E S -------------+--+-----------*
*  NAME  /                  DESCRIPTION                                *
*--------+-------------------------------------------------------------*
* Boundary conditions
*   PIJABC / TO ELIMINATE SURFACES USING B.C. OF THE SYSTEM
*   PIJAAA / TO ELIMINATE SURFACES FOR PIJKS USING B.C. OF THE SYSTEM
* Various functions
*   PIJWPR / TO PRINT CP MATRICES IN SUM FORMAT
*   PIJSMD / TO EVALUATE SCATTERING-MODIFIED CP MATRIX
*   PIJCMP / COMPRESS CP MATRIX TO SYMETRIC FORMAT
*   PIJD2S / CHARGE PROBKS MATRICES IN THE DRAGON SQUARE FORMAT
*   PIJD2R / CHARGE PIJ MATRICES IN THE DRAGON SYMMETRIZED FORMAT
*   PIJKST / COMPUTE PIJK* MATRICES
*--------+-------------------------------------------------------------*
*
      USE               GANLIB
      IMPLICIT          NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      LOGICAL           LEAKSW
      TYPE(C_PTR)       IPTRK
      INTEGER           IPRT, NSOUT, NREG, NBMIX, NANI, MATCOD(NREG),
     >                  NELPIJ, IPIJK, N2PRO, NSBG, NPSYS(NSBG), NPST,
     >                  NALBP, MATALB(-NSOUT:NREG)
      REAL              VOLUME(NREG), XSSIGT(0:NBMIX,NSBG),
     >                  XSSIGW(0:NBMIX,NANI,NSBG),ALBP(NALBP,NSBG),
     >                  VOLSUR(-NSOUT:NREG,NSBG),PIJ(NELPIJ,IPIJK,NSBG),
     >                  PROBKS(NPST,NSBG)
      DOUBLE PRECISION  DPROB(N2PRO,NSBG),DPROBX(N2PRO,NSBG)
*----
*  LOCAL VARIABLES
*----
      INTEGER           IOUT, ICPALL, ICPEND, MXGAUS, NSTATE, MAXCDA
      PARAMETER       ( IOUT=6, ICPALL=4, ICPEND=3, MXGAUS=64,
     >                  NSTATE=40, MAXCDA=30 )
      CHARACTER         NAMSBR*6
      PARAMETER       ( NAMSBR='PIJWIJ')
      INTEGER           ILONG,ITYPE,NPROB,ISBG,ISTATE(NSTATE),
     >                  ICODE(MAXCDA)
      REAL              FACT,ALBEDO(MAXCDA),ALBG(MAXCDA)
      LOGICAL           LSKIP,SWNZBC,SWVOID
*
      INTEGER           MSYM,IU,IL,ISOUT,IIN,I,J,IBM,IOP,INDPIJ,IJKS,
     >                  IUN,KSPEC,LOPT,NNREG,IVV,JUN,ISA
*----
*  Variables for NXT: inline tracking
*----
      INTEGER           ILCMUP,ILCMDN
      PARAMETER        (ILCMUP=1,ILCMDN=2)
      DOUBLE PRECISION  DZERO,DONE,DTWO
      PARAMETER        (DZERO=0.0D0,DONE=1.0D0,DTWO=2.0D0)
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MATRT
      REAL, ALLOCATABLE, DIMENSION(:) :: FFACT
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SIGTAL
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: PSST,PSVT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: PCSCT
*----
*  INTRINSIC FUNCTION FOR POSITION IN CONDENSE PIJ MATRIX
*----
      INTEGER INDPOS
      INDPOS(I,J)=MAX(I,J)*(MAX(I,J)-1)/2+MIN(I,J)
*----
*  RECOVER TRACKING INFORMATION
*----
      ISTATE(:NSTATE)=0
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      KSPEC=ISTATE(10)
      CALL LCMLEN(IPTRK,'ICODE',ILONG,ITYPE)
      IF(ILONG.GT.MAXCDA) CALL XABORT('PIJWIJ: MAXCDA OVERFLOW(1).')
      ICODE(:MAXCDA)=0
      CALL LCMGET(IPTRK,'ICODE',ICODE)
      CALL LCMGET(IPTRK,'ALBEDO',ALBG)
*----
*  PREPARE FOR MULTIGROUP CALCULATION
*----
      ALLOCATE(SIGTAL(-NSOUT:NREG,NSBG))
      SWNZBC= .FALSE.
      SWVOID= .FALSE.
      DO ISBG=1,NSBG
        IF(NPSYS(ISBG).NE.0) THEN
          ALBEDO(:MAXCDA)=ALBG(:MAXCDA)
          IF(NALBP .GT. 0) THEN
            DO ISA=1,MAXCDA
              IF(ICODE(ISA).GT.0) ALBEDO(ISA)=ALBP(ICODE(ISA),ISBG)
            ENDDO
          ENDIF
          DO IUN= -NSOUT, -1
            SIGTAL(IUN,ISBG)= ALBEDO(-MATALB(IUN))
            SWNZBC= SWNZBC.OR.(SIGTAL(IUN,ISBG).NE.0.0)
          ENDDO
          IUN=0
          SIGTAL(IUN,ISBG)= 0.0
          DO IUN= 1, NREG
            SIGTAL(IUN,ISBG)= XSSIGT(MATCOD(IUN),ISBG)
            IF( SIGTAL(IUN,ISBG) .EQ. 0.0 ) SWVOID= .TRUE.
          ENDDO
        ENDIF
      ENDDO
*----
*  DOUBLE PRECISION TO REAL FOR DIRECTIONAL PIJ MATRICES
*----
      IF(IPIJK.EQ.4) THEN
         DO 2070 ISBG=1,NSBG
         IF(NPSYS(ISBG).EQ.0) GO TO 2070
         CALL PIJD2S(NREG,NSOUT,DPROBX(1,ISBG),PROBKS(1,ISBG))
 2070    CONTINUE
      ENDIF
      IF( KSPEC.EQ.0 )THEN
*----
*  ELIMINATION OF SURFACES FOR PIJ
*----
         IF( SWNZBC )THEN
            ALLOCATE(PSST(NSOUT*NSOUT),PSVT(NSOUT*NREG),MATRT(NSOUT))
            CALL LCMLEN(IPTRK,'BC-REFL+TRAN',ILONG,ITYPE)
            IF(ILONG.EQ.NSOUT) THEN
              CALL LCMGET(IPTRK,'BC-REFL+TRAN',MATRT)
            ELSE
               WRITE(IOUT,9000) NAMSBR
               DO 130 ISOUT=1,NSOUT
                 MATRT(ISOUT)=ISOUT
 130           CONTINUE
            ENDIF
            DO 2080 ISBG=1,NSBG
              IF(NPSYS(ISBG).EQ.0) GO TO 2080
              CALL PIJABC(NREG,NSOUT,NPROB,SIGTAL(-NSOUT,ISBG),MATRT,
     >                    DPROB(1,ISBG),PSST,PSVT)
*----
*  ELIMINATION OF SURFACES FOR PIJX AND CREATION OF PIJXX
*----
            IF(IPIJK.EQ.4) THEN
               CALL PIJAAA(NREG,NSOUT,SIGTAL(-NSOUT,ISBG),
     >                     DPROBX(1,ISBG),PSVT,PROBKS(1,ISBG))
               CALL PIJABC(NREG,NSOUT,NPROB,SIGTAL(-NSOUT,ISBG),MATRT,
     >                     DPROBX(1,ISBG),PSST,PSVT)
            ENDIF
 2080     CONTINUE
*
            DEALLOCATE(MATRT,PSVT,PSST)
         ENDIF
      ENDIF
*
      ALLOCATE(FFACT(NREG))
      DO 2090 ISBG=1,NSBG
      IF(NPSYS(ISBG).EQ.0) GO TO 2090
      IF( IPRT.GE.ICPEND )THEN
         LOPT= +1
         MSYM=1
         WRITE(IOUT,'(1H )')
         WRITE(IOUT,'(35H   COLLISION PROBABILITIES OUTPUT: ,
     >                35H *AFTER* ALBEDO REDUCTION          )')
         CALL PIJWPR(LOPT,NREG,NSOUT,SIGTAL(-NSOUT,ISBG),
     >               DPROB(1,ISBG),VOLSUR(1,ISBG),MSYM)
*
         IF(IPIJK.EQ.4) THEN
           WRITE(IOUT,'(35H   X-DIRECT. COLL. PROBAB. OUTPUT: ,
     >                  35H *AFTER* ALBEDO REDUCTION          )')
           CALL PIJWPR(LOPT,NREG,NSOUT,SIGTAL(-NSOUT,ISBG),
     >                 DPROBX(1,ISBG),VOLSUR(1,ISBG),MSYM)
           WRITE(IOUT,'(35H0 X-DIRECT. COLL. PROBAB." OUTPUT: ,
     >                  35H PIJX"=PIJX+PISX*(1/(1-PSS))*PSJ   )')
           MSYM=0
           CALL PIJWPR(LOPT,NREG,NSOUT,SIGTAL(-NSOUT,ISBG),
     >                 DPROBX(1,ISBG),VOLSUR(1,ISBG),MSYM)
         ENDIF
*
      ENDIF
*----
*  CHARGE PIJ MATRIX IN THE DRAGON SYMMETRIZED FORMAT
*----
      DO 160 IIN=1,NREG
         IF(SIGTAL(IIN,ISBG).EQ.0.0) THEN
            FFACT(IIN)=1.0
         ELSE
            FFACT(IIN)=1.0/SIGTAL(IIN,ISBG)
         ENDIF
  160 CONTINUE
      CALL PIJD2R(NREG,NSOUT,DPROB(1,ISBG),FFACT,.FALSE.,NELPIJ,
     >            N2PRO,PIJ(1,1,ISBG))
*----
*  CHARGE PIJX AND PIJY MATRICES IN THE DRAGON SYMMETRIZED FORMAT
*  ( PIJX=PIJY ), AND PIJZ CALCULATION ( PIJZ=3*PIJ-PIJX-PIJY )
*  AND THE SAME FOR FULL MATRICES OF PIJX", PIJY" AND PIJZ"
*----
      IF(IPIJK.EQ.4) THEN
        NNREG=NREG*NREG
        CALL PIJD2R(NREG,NSOUT,DPROBX(1,ISBG),FFACT,.TRUE.,NELPIJ,
     >              N2PRO,PIJ(1,2,ISBG))
        IVV=0
        DO 181 IUN=1,NREG
          IU=IUN
          IL=(IUN-1)*NREG+1
          DO 191 JUN=1,IUN
            IVV=IVV+1
            PROBKS(IL,ISBG)=1.5*PROBKS(IL,ISBG)*FFACT(IUN)*FFACT(JUN)
            IF(IL.NE.IU)PROBKS(IU,ISBG)=1.5*PROBKS(IU,ISBG)*
     >                      FFACT(IUN)*FFACT(JUN)
            PIJ(IVV,3,ISBG)=PIJ(IVV,2,ISBG)
            PROBKS(NNREG+IL,ISBG)=PROBKS(IL,ISBG)
            PROBKS(NNREG+IU,ISBG)=PROBKS(IU,ISBG)
            PIJ(IVV,4,ISBG)=3*PIJ(IVV,1,ISBG)-PIJ(IVV,2,ISBG)
     >                                       -PIJ(IVV,3,ISBG)
            PROBKS(2*NNREG+IL,ISBG)=3*PIJ(IVV,1,ISBG)
     >      -PROBKS(IL,ISBG)-PROBKS(NNREG+IL,ISBG)
            PROBKS(2*NNREG+IU,ISBG)=3*PIJ(IVV,1,ISBG)
     >      -PROBKS(IU,ISBG)-PROBKS(NNREG+IU,ISBG)
            IU=IUN+JUN*NREG
            IL=IL+1
  191     CONTINUE
  181   CONTINUE
*----
*  COMPUTE PIJ**(-1)*PIJK*
*----
        CALL PIJKST(IPRT,NREG,PIJ(1,1,ISBG),PROBKS(1,ISBG))
      ENDIF
 2090 CONTINUE
      DEALLOCATE(FFACT)
*
      DEALLOCATE(SIGTAL)
*----
*  CHECK IF SCATTERING REDUCTION IS REQUIRED
*----
      ALLOCATE(PCSCT(NREG,2*NREG))
      DO 3000 ISBG=1,NSBG
      IF(NPSYS(ISBG).EQ.0) GO TO 3000
      LSKIP=.TRUE.
      DO 200 IBM=1,NBMIX
        LSKIP=LSKIP.AND.(XSSIGW(IBM,1,ISBG).EQ.0.0)
  200 CONTINUE
*----
*  COMPUTE THE SCATTERING-REDUCED CP MATRICES
*----
      IOP=1
      IF(.NOT.LSKIP) THEN
        CALL PIJSMD(IPRT,NBMIX,NREG,MATCOD,VOLUME,XSSIGW(0,1,ISBG),
     >              XSSIGT(0,ISBG),LEAKSW,PIJ(1,1,ISBG),PCSCT,IOP)
        DO 220 I=1,NREG
          FACT=VOLUME(I)
          DO 210 J=1,NREG
            INDPIJ=INDPOS(I,J)
            PIJ(INDPIJ,1,ISBG)=REAL(PCSCT(I,J))*FACT
  210     CONTINUE
  220   CONTINUE
      ENDIF
*-------
      IF(IPIJK.EQ.4) THEN
        IOP=4
        IF(.NOT.LSKIP) THEN
*         P1 SCATTERING REDUCTION OF THE DIRECTIONNAL CP MATRICES.
          IF(NANI.LT.2) CALL XABORT('PIJWIJ: ANISOTROPIC SCAT MISSING.')
          DO 250 IJKS=1,3
            CALL PIJSMD(IPRT,NBMIX,NREG,MATCOD,VOLUME,XSSIGW(0,2,ISBG),
     >                  XSSIGT(0,ISBG),LEAKSW,PIJ(1,IJKS+1,ISBG),PCSCT,
     >                  IOP)
            DO 240 I=1,NREG
              FACT=VOLUME(I)
              DO 230 J=1,NREG
                INDPIJ=INDPOS(I,J)
                PIJ(INDPIJ,IJKS+1,ISBG)=REAL(PCSCT(I,J))*FACT
  230         CONTINUE
  240       CONTINUE
  250     CONTINUE
        ENDIF
      ENDIF
 3000 CONTINUE
      DEALLOCATE(PCSCT)
      RETURN
*
 9000 FORMAT(1X,A6,': *** WARNING *** '/
     >       ' REFLECTION/TRANSMISSION MATRIX MISSING'/
     >       ' USE IDENTITY REFLECTION MATRIX')
      END
