*DECK EDIRES
      SUBROUTINE EDIRES(IPEDIT,IPFLUX,IPLIB,IADJ,NL,NDEL,NBESP,NBISO,
     1 NDEPL,ISONAM,ISONRF,IPISO,MIX,TN,NED,HVECT,NOUT,HVOUT,IPRINT,
     2 NGROUP,NGCOND,NBMIX,NREGIO,NMERGE,NDFI,NDFP,ILEAKS,ILUPS,NW,
     3 MATCOD,VOLUME,KEYFLX,CURNAM,IGCOND,IMERGE,FLUXES,AFLUXE,EIGENK,
     4 EIGINF,B2,DEN,ITYPE,LSISO,EMEVF,EMEVG,DECAY,YIELD,FIPI,FIFP,
     5 PYIELD,ITRANC,LISO,NMLEAK)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of residual macroscopic cross sections.
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
* IPEDIT  pointer to the edition LCM object (L_EDIT signature).
* IPFLUX  pointer to the solution LCM object (L_FLUX signature).
* IPLIB   pointer to the reference microscopic cross section library
*         LCM object (L_LIBRARY signature).
* IADJ    type of flux weighting:
*         =0: direct flux weighting;
*         =1: direct-adjoint flux weighting.
* NL      number of Legendre orders required in the calculation
*         (NL=1 or higher).
* NDEL    number of delayed precursor groups.
* NBESP   number of energy-dependent fission spectra.
* NBISO   number of isotopes.
* NDEPL   number of depleting isotopes.
* ISONAM  local names of NBISO isotopes:
*         chars 1 to 8  is the local isotope name;
*         chars 9 to 12 is a suffix function of the mix number.
* ISONRF  library name of isotopes.
* IPISO   pointer array towards microlib isotopes.
* MIX     mixture number associated with each isotope.
* TN      absolute temperature associated with each isotope.
* NED     number of extra vector edits from MATXS.
* HVECT   MATXS names of the extra vector edits.
* NOUT    number of output cross section types (set to zero to recover
*         all cross section types).
* HVOUT   MATXS names of the output cross section types.
* IPRINT  print index.
* NGROUP  number of energy groups.
* NGCOND  number of condensed groups.
* NBMIX   number of mixtures.
* NREGIO  number of volumes.
* NMERGE  number of merged regions.
* NDFI    number of fissile isotopes.
* NDFP    number of fission products.
* ILEAKS  leakage calculation type: =0: no leakage; =1: homogeneous
*         leakage (Diffon); =2: isotropic streaming (Ecco);
*         =3: anisotropic streaming (Tibere).
* ILUPS   up-scattering removing flag (=1 to remove up-scattering from
*         output cross-sections).
* NW      type of weighting for P1 cross section info (=0: P0 ; =1: P1).
* MATCOD  mixture index per volume.
* VOLUME  volumes.
* KEYFLX  position of average fluxes.
* CURNAM  name of the LCM directory where the microscopic cross sections
*         are stored (blank name means no save).
* IGCOND  limits of condensed groups.
* IMERGE  index of merged regions.
* FLUXES  fluxes.
* AFLUXE  adjoint fluxes.
* EIGENK  effective multiplication factor.
* EIGINF  infinite multiplication factor.
* B2      bucklings.
* DEN     number density of each isotope.
* ITYPE   type of each isotope.
* LSISO   flag for isotopes saved.
* EMEVF   fission production energy.
* EMEVG   capture production energy.
* DECAY   radioactive decay constant.
* YIELD   group-ordered condensed fission product yield.
* FIPI    fissile isotope index assigned to each microlib isotope.
* FIFP    fission product index assigned to each microlib isotope.
* PYIELD  fissile isotope ordered condensed fission product yield.
* ITRANC  type of transport correction (=0: no correction).
* LISO    =.TRUE. if we want to keep all the isotopes after 
*         homogeneization.
* NMLEAK  number of leakage zones.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPEDIT,IPFLUX,IPLIB,IPISO(NBISO)
      INTEGER NL,NDEL,NBESP,NBISO,NDEPL,ISONAM(3,NBISO),ISONRF(3,NBISO),
     1 MIX(NBISO),NED,NOUT,IPRINT,NGROUP,NGCOND,NBMIX,NREGIO,NMERGE,
     2 NDFI,NDFP,ILEAKS,ILUPS,NW,MATCOD(NREGIO),KEYFLX(NREGIO),
     3 IGCOND(NGCOND),IMERGE(NREGIO),ITYPE(NBISO),LSISO(NBISO),
     4 FIPI(NBISO,NMERGE),FIFP(NBISO,NMERGE),ITRANC,NMLEAK
      REAL TN(NBISO),VOLUME(NREGIO),FLUXES(NREGIO,NGROUP,NW+1),
     1 EIGENK,EIGINF,B2(4),DEN(NBISO),EMEVF(NBISO),EMEVG(NBISO),
     2 DECAY(NBISO),YIELD(NGCOND+1,NDFP,NMERGE),PYIELD(NDFI,NDFP,NMERGE)
      CHARACTER HVECT(NED)*8,HVOUT(NOUT)*8,CURNAM*12
      LOGICAL LISO
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      TYPE(C_PTR) JPEDIT,KPEDIT,IPWORK,JPWORK,KPWORK
      CHARACTER TEXT8*8,TEXT12*12,CM*2
      LOGICAL LWD,LYIEL,LFISS
      INTEGER IPAR(NSTATE)
      INTEGER, ALLOCATABLE, DIMENSION(:) :: LSIS2,IEVOL2,ISMIX,ISTYP,
     1 ISTOD,ITYPRO,JPIFI
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IHNISO,IHNIRF
      REAL, ALLOCATABLE, DIMENSION(:) :: WDLA,SDEN,VOLISO,TNISO,WORK,
     1 WPY,DENTOT,DAWR,TNTOT,YIELD2,PYIELD2
      REAL, ALLOCATABLE, DIMENSION(:,:) :: GAS,SIGS,PNFIRA
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: WSCAT,WORK2
      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:) :: HMAKE
*----
*  SCRATCH STORAGE ALLOCATION (PART 1)
*----
      MAXH=9+NBESP+2*NDEL+NED+NL+3*NW
      ALLOCATE(LSIS2(NBISO),IEVOL2(NBISO),JPIFI(NDFI),ITYPRO(NL))
      ALLOCATE(WDLA(NDEL),WSCAT(NGCOND,NGCOND,NL),GAS(NGCOND,MAXH),
     1 WORK(NGCOND+1),WPY(NDFI),PNFIRA(0:NDEL,2),
     2 WORK2(NGCOND,NGCOND,NL),DENTOT(NMERGE),DAWR(NMERGE),
     3 TNTOT(NMERGE),YIELD2(1+NGCOND),PYIELD2(NDFI))
      ALLOCATE(HMAKE(MAXH+NL))
*----
*  RECOVER THE RADIOACTIVE DECAY CONSTANTS OF DELAYED NEUTRON GROUPS
*  FROM THE MACROLIB IF THEY EXIST.
*----
      IOF0H=8+NED+NL+3*NW
      IOF1H=8+NED+NL+3*NW+NDEL
      IF(IPRINT.GT.3) THEN
         WRITE(6,'(/36H EDIRES: COMPUTE A RESIDUAL ISOTOPE.)')
      ENDIF
      CALL LCMOP(IPWORK,'*TEMPORARY*',0,1,0)
      LWD=.FALSE.
      IF(CURNAM.EQ.' ') CALL XABORT('EDIRES: NO CURNAM DIRECTORY.')
      CALL LCMSIX(IPEDIT,CURNAM,1) ! step up CURNAM
      CALL LCMSIX(IPEDIT,'MACROLIB',1)
      CALL LCMLEN(IPEDIT,'LAMBDA-D',ILONG,ITYLCM)
      LWD=(ILONG.EQ.NDEL).AND.(NDEL.GT.0)
      IF(LWD) CALL LCMGET(IPEDIT,'LAMBDA-D',WDLA)
      CALL LCMSIX(IPEDIT,' ',2)
      CALL LCMSIX(IPEDIT,' ',2)
*
      IF(LWD) THEN
         CALL LCMSIX(IPWORK,'DEFAULT',1)
         CALL LCMSIX(IPWORK,'MACROLIB',1)
         CALL LCMPUT(IPWORK,'LAMBDA-D',NDEL,2,WDLA)
         CALL LCMSIX(IPWORK,' ',2)
         CALL LCMSIX(IPWORK,' ',2)
      ENDIF
*----
*  COMPUTE MICROSCOPIC CROSS SECTIONS OF REMAINING ISOTOPES. WE SET
*  NDFI=0 TO GET RID OF PPF YIELDS.
*----
      DO 10 ISO=1,NBISO
      LSIS2(ISO)=0
      IEVOL2(ISO)=1
      IF(LSISO(ISO).EQ.0) LSIS2(ISO)=1
   10 CONTINUE
      IPRIN2=MAX(0,IPRINT-2)
      TEXT12='DEFAULT'
      CALL EDIMIC(IPWORK,IPFLUX,IPLIB,IADJ,NL,NDEL,NBESP,NBISO,NDEPL,
     1 ISONAM,ISONRF,IPISO,MIX,TN,NED,HVECT,NOUT,HVOUT,IPRIN2,NGROUP,
     2 NGCOND,NBMIX,NREGIO,NMERGE,0,0,ILEAKS,ILUPS,NW,MATCOD,VOLUME,
     3 KEYFLX,TEXT12,IGCOND,IMERGE,FLUXES,AFLUXE,EIGENK,EIGINF,B2,DEN,
     4 ITYPE,IEVOL2,LSIS2,EMEVF,EMEVG,DECAY,YIELD,FIPI,FIFP,PYIELD,
     5 ITRANC,LISO,NMLEAK)
*
      CALL LCMSIX(IPEDIT,CURNAM,1)
      CALL LCMGET(IPEDIT,'STATE-VECTOR',IPAR)
      JJISO=IPAR(2)
      JPEDIT=LCMLID(IPEDIT,'ISOTOPESLIST',JJISO+NMERGE)
      CALL LCMSIX(IPWORK,'DEFAULT',1)
      CALL LCMGET(IPWORK,'STATE-VECTOR',IPAR)
      JJWRK=IPAR(2)
*----
*  SCRATCH STORAGE ALLOCATION (PART 2)
*----
      MAXISO=MAX(JJISO+NMERGE,JJWRK)
      ALLOCATE(IHNISO(3,MAXISO),SDEN(MAXISO),IHNIRF(3,MAXISO),
     1 ISMIX(MAXISO),ISTYP(MAXISO),ISTOD(MAXISO),VOLISO(MAXISO),
     2 TNISO(MAXISO))
*----
*  RECOVER INFORMATION FROM EDIMIC
*----
      IF(JJWRK.GT.0) THEN
         CALL LCMGET(IPWORK,'ISOTOPESUSED',IHNISO)
         CALL LCMGET(IPWORK,'ISOTOPESDENS',SDEN)
         CALL LCMGET(IPWORK,'ISOTOPESMIX',ISMIX)
         CALL LCMGET(IPWORK,'ISOTOPESTEMP',TNISO)
         JPWORK=LCMGID(IPWORK,'ISOTOPESLIST')
      ENDIF
*----
*  LOOP OVER HOMOGENEOUS MIXTURES.
*----      
      DO 240 INM=1,NMERGE
      DO 20 J=1,MAXH+NL
      HMAKE(J)=' '
   20 CONTINUE
      GAS(:NGCOND,:MAXH)=0.0
      WSCAT(:NGCOND,:NGCOND,:NL)=0.0
      PNFIRA(0:NDEL,2)=0.0
      YIELD2(:1+NGCOND)=0.0
      PYIELD2(:NDFI)=0.0
      DENTOT(INM)=0.0
      DAWR(INM)=0.0
      DECISO=0.0
      LFISS=.FALSE.
      DO 170 ISO=1,JJWRK
      IF(ISMIX(ISO).EQ.INM) THEN
         WRITE(TEXT12,'(3A4)') (IHNISO(I0,ISO),I0=1,3)
         DDEN=SDEN(ISO)
         DENTOT(INM)=DENTOT(INM)+DDEN
         KPWORK=LCMGIL(JPWORK,ISO) ! set ISO-th isotope
         CALL LCMLEN(KPWORK,'AWR',LENGTH,ITYLCM)
         IF(LENGTH.EQ.1) THEN
            CALL LCMGET(KPWORK,'AWR',FLOTT)
            DAWR(INM)=DAWR(INM)+DDEN*FLOTT
         ENDIF
         TNTOT(INM)=TNISO(ISO)
         CALL LCMLEN(KPWORK,'DECAY',LENGTH,ITYLCM)
         IF(LENGTH.EQ.1) THEN
            CALL LCMGET(KPWORK,'DECAY',FLOTT)
            DECISO=DECISO+FLOTT*DDEN
         ENDIF
         IF(NDFI.GT.0) THEN
           CALL LCMLEN(KPWORK,'YIELD',LENGTH,ITYLCM)
           IF(LENGTH.EQ.NGCOND+1) THEN
             CALL LCMGET(KPWORK,'YIELD',WORK)
             DO 30 IGR=1,NGCOND+1
             YIELD2(IGR)=YIELD2(IGR)+WORK(IGR)
   30      CONTINUE
           ENDIF
           CALL LCMLEN(KPWORK,'PYIELD',LENGTH,ITYLCM)
           IF((LENGTH.GT.0).AND.(LENGTH.EQ.NDFI)) THEN
             CALL LCMGET(KPWORK,'PIFI',JPIFI)
             CALL LCMGET(KPWORK,'PYIELD',WPY)
             DO 40 I=1,NDFI
             PYIELD2(I)=PYIELD2(I)+WPY(I)
   40        CONTINUE
           ENDIF
         ENDIF
*
*        SET ARRAY HMAKE.
         DO 45 IW=1,MIN(NW+1,10)
            WRITE(TEXT8,'(3HNWT,I1)') IW-1
            CALL LCMLEN(KPWORK,TEXT8,ILONG,ITYLCM)
            IF(ILONG.EQ.NGCOND) HMAKE(IW)=TEXT8
            WRITE(TEXT8,'(4HNWAT,I1)') IW-1
            CALL LCMLEN(KPWORK,TEXT8,ILONG,ITYLCM)
            IF(ILONG.EQ.NGCOND) HMAKE(1+NW+IW)=TEXT8
            WRITE(TEXT8,'(4HNTOT,I1)') IW-1
            CALL LCMLEN(KPWORK,TEXT8,ILONG,ITYLCM)
            IF(ILONG.EQ.NGCOND) HMAKE(2+2*NW+IW)=TEXT8
   45    CONTINUE
         IOF=3+3*NW
         DO 50 IL=0,NL-1
            IOF=IOF+1
            WRITE (CM,'(I2.2)') IL
            CALL LCMLEN(KPWORK,'SIGS'//CM,ILONG,ITYLCM)
            IF(ILONG.EQ.NGCOND) HMAKE(IOF)='SIGS'//CM
   50    CONTINUE
         IOF=IOF+1
         CALL LCMLEN(KPWORK,'NUSIGF',ILONG,ITYLCM)
         IF(ILONG.EQ.NGCOND) THEN
           LFISS=.TRUE.
           HMAKE(IOF)='NUSIGF'
         ENDIF
         DO 60 IED=1,NED
         IOF=IOF+1
         CALL LCMLEN(KPWORK,HVECT(IED),ILONG,ITYLCM)
         IF(ILONG.EQ.NGCOND) HMAKE(IOF)=HVECT(IED)
   60    CONTINUE
         CALL LCMLEN(KPWORK,'H-FACTOR',ILONG,ITYLCM)
         IF(ILONG.EQ.NGCOND) HMAKE(IOF+1)='H-FACTOR'
         CALL LCMLEN(KPWORK,'OVERV',ILONG,ITYLCM)
         IF(ILONG.EQ.NGCOND) HMAKE(IOF+2)='OVERV'
         CALL LCMLEN(KPWORK,'TRANC',ILONG,ITYLCM)
         IF(ILONG.EQ.NGCOND) HMAKE(IOF+3)='TRANC'
         CALL LCMLEN(KPWORK,'STRD',ILONG,ITYLCM)
         IF(ILONG.EQ.NGCOND) HMAKE(IOF+4)='STRD'
         IOF=IOF+4
         DO 70 IDEL=1,NDEL
         IOF=IOF+1
         WRITE(TEXT8,'(6HNUSIGF,I2.2)') IDEL
         CALL LCMLEN(KPWORK,TEXT8,ILONG,ITYLCM)
         IF(ILONG.EQ.NGCOND) HMAKE(IOF)=TEXT8
   70    CONTINUE
         IOF=IOF+1
         CALL LCMLEN(KPWORK,'CHI',ILONG,ITYLCM)
         IF(ILONG.EQ.NGCOND) HMAKE(IOF)='CHI'
         DO 80 IDEL=1,NDEL
         IOF=IOF+1
         WRITE(TEXT8,'(3HCHI,I2.2)') IDEL
         CALL LCMLEN(KPWORK,TEXT8,ILONG,ITYLCM)
         IF(ILONG.EQ.NGCOND) HMAKE(IOF)=TEXT8
   80    CONTINUE
         DO 85 ISP=1,NBESP
         IOF=IOF+1
         WRITE(TEXT8,'(5HCHI--,I2.2)') ISP
         CALL LCMLEN(KPWORK,TEXT8,ILONG,ITYLCM)
         IF(ILONG.EQ.NGCOND) HMAKE(IOF)=TEXT8
   85    CONTINUE
         IF(IOF.NE.MAXH) CALL XABORT('EDIRES: WRONG OFFSET.')
*
         PNFIRA(0:NDEL,1)=0.0
         DO 150 J=1,MAXH
         IF(HMAKE(J).NE.' ') THEN
            CALL LCMLEN(KPWORK,HMAKE(J),ILONG,ITYLCM)
            IF(ILONG.GT.0) THEN
               CALL LCMGET(KPWORK,HMAKE(J),WORK)
               IF(HMAKE(J).EQ.'NUSIGF') THEN
                  DO 90 IGR=1,NGCOND
                    DEL=WORK(IGR)*GAS(IGR,1)*MAX(DDEN,1.0E-30)
                    PNFIRA(0,1)=PNFIRA(0,1)+DEL
                    PNFIRA(0,2)=PNFIRA(0,2)+DEL
                    GAS(IGR,J)=GAS(IGR,J)+WORK(IGR)*DDEN
   90             CONTINUE
               ELSE IF(HMAKE(J)(:3).EQ.'NUS') THEN
                  IDEL=J-IOF0H
                  DO 100 IGR=1,NGCOND
                    DEL=WORK(IGR)*GAS(IGR,1)*MAX(DDEN,1.0E-30)
                    PNFIRA(IDEL,1)=PNFIRA(IDEL,1)+DEL
                    PNFIRA(IDEL,2)=PNFIRA(IDEL,2)+DEL
                    GAS(IGR,J)=GAS(IGR,J)+WORK(IGR)*DDEN
  100             CONTINUE
               ELSE IF(HMAKE(J)(:3).EQ.'NWT') THEN
                  DO 110 IGR=1,NGCOND
                    GAS(IGR,J)=WORK(IGR)
  110             CONTINUE
               ELSE IF((HMAKE(J).EQ.'CHI').OR.
     1                 (HMAKE(J)(:5).EQ.'CHI--')) THEN
                  DO 120 IGR=1,NGCOND
                    GAS(IGR,J)=GAS(IGR,J)+WORK(IGR)*PNFIRA(0,1)
  120             CONTINUE
               ELSE IF(HMAKE(J)(:3).EQ.'CHI') THEN
                  IDEL=J-IOF1H-1
                  DO 130 IGR=1,NGCOND
                    GAS(IGR,J)=GAS(IGR,J)+WORK(IGR)*PNFIRA(IDEL,1)
  130             CONTINUE
               ELSE
                  DO 140 IGR=1,NGCOND
                    GAS(IGR,J)=GAS(IGR,J)+WORK(IGR)*DDEN
  140             CONTINUE
               ENDIF
            ENDIF
         ENDIF
  150    CONTINUE
         CALL LCMLEN(KPWORK,'SCAT-SAVED',ILONG,ITYLCM)
         IF(ILONG.GT.0) THEN
            ALLOCATE(SIGS(NGCOND,NL))
            CALL XDRLGS(KPWORK,-1,IPRINT,0,NL-1,1,NGCOND,SIGS(1,1),
     1      WORK2,ITYPRO)
            DEALLOCATE(SIGS)
            DO 162 IL=1,NL
            WRITE (CM,'(I2.2)') IL-1
            IF(ITYPRO(IL).NE.0) HMAKE(MAXH+IL)='SCAT'//CM
            DO 161 JGR=1,NGCOND
            DO 160 IGR=1,NGCOND
              WSCAT(IGR,JGR,IL)=WSCAT(IGR,JGR,IL)+WORK2(IGR,JGR,IL)*DDEN
  160    CONTINUE
  161    CONTINUE
  162    CONTINUE
         ENDIF
      ENDIF
  170 CONTINUE
      IF((DENTOT(INM).GT.0.0).OR.LFISS) THEN
         JJISO=JJISO+1
         IF(JJISO.GT.MAXISO) CALL XABORT('EDIRES: MAXISO OVERFLOW(1).')
         WRITE(TEXT12,'(A8,I4.4)') '*MAC*RES',INM
         IF(IPRINT.GT.0) WRITE (6,600) TEXT12,JJISO
         KPEDIT=LCMDIL(JPEDIT,JJISO) ! set JJISO-th isotope
         CALL LCMPTC(KPEDIT,'ALIAS',12,TEXT12)
         CALL LCMPUT(KPEDIT,'AWR',1,2,DAWR(INM))
         DECISO=DECISO/DENTOT(INM)
         IF(DECISO.GT.0.0) CALL LCMPUT(KPEDIT,'DECAY',1,2,DECISO)
         IF(NDFI.GT.0) THEN
           LYIEL=.FALSE.
           DO 175 IGR=1,NGCOND+1
             LYIEL=LYIEL.OR.(YIELD2(IGR).GT.0.0)
  175      CONTINUE
           IF(LYIEL) THEN
             CALL LCMPUT(KPEDIT,'YIELD',NGCOND+1,2,YIELD2)
             CALL LCMPUT(KPEDIT,'PYIELD',NDFI,2,PYIELD2)
             CALL LCMPUT(KPEDIT,'PIFI',NDFI,1,JPIFI)
           ENDIF
         ENDIF
         IF(NOUT.GT.0) THEN
           DO J=1,MAXH+NL
             DO IOUT=1,NOUT
               IF(HMAKE(J).EQ.HVOUT(IOUT)) GO TO 180
             ENDDO
             HMAKE(J)=' '
  180        CONTINUE
           ENDDO
         ENDIF
         DO 210 J=1,MAXH
         IF(HMAKE(J).EQ.'OVERV') THEN
            DO 185 IGR=1,NGCOND
              GAS(IGR,J)=GAS(IGR,J)/DENTOT(INM)
  185       CONTINUE
         ELSE IF((HMAKE(J).EQ.'CHI').OR.(HMAKE(J)(:5).EQ.'CHI--')) THEN
            DO 190 IGR=1,NGCOND
            IF(GAS(IGR,J).NE.0.0) THEN
               GAS(IGR,J)=GAS(IGR,J)/PNFIRA(0,2)
            ENDIF
  190       CONTINUE
         ELSE IF(HMAKE(J)(:3).EQ.'CHI') THEN
            IDEL=J-IOF1H-1
            DO 200 IGR=1,NGCOND
            IF(GAS(IGR,J).NE.0.0) THEN
               GAS(IGR,J)=GAS(IGR,J)/PNFIRA(IDEL,2)
            ENDIF
  200       CONTINUE
         ENDIF
         IF((HMAKE(J).NE.' ').AND.(HMAKE(J)(:4).NE.'SIGS')) THEN
            CALL LCMPUT(KPEDIT,HMAKE(J),NGCOND,2,GAS(1,J))
         ENDIF
  210    CONTINUE
         DO 220 IL=1,NL
         ITYPRO(IL)=0
         IF(HMAKE(MAXH+IL).NE.' ') ITYPRO(IL)=1
  220    CONTINUE
         IF(ITYPRO(1).GT.0) THEN
            CALL XDRLGS(KPEDIT,1,IPRINT,0,NL-1,1,NGCOND,GAS(1,4+3*NW),
     1      WSCAT,ITYPRO)
         ENDIF
         IF(LWD) CALL LCMPUT(KPEDIT,'LAMBDA-D',NDEL,2,WDLA)
*
         IF(IPRINT.GT.3) THEN
            WRITE(6,'(/17H NUMBER DENSITY =,1P,E12.4)') 1.0
            WRITE(6,'(23H WEIGHTED ATOMIC MASS =,1P,E13.5)') DAWR(INM)
            DO 230 J=1,MAXH
            IF(HMAKE(J).NE.' ') THEN
               WRITE (6,610) HMAKE(J),(GAS(I,J),I=1,NGCOND)
            ENDIF
  230       CONTINUE
            WRITE (6,610) 'SIGA    ',(GAS(I,3+2*NW)-GAS(I,4+3*NW),
     >      I=1,NGCOND)
            WRITE (6,610) 'SIGW00  ',(WSCAT(I,I,1),I=1,NGCOND)
            IF(NL.GT.1) THEN
               WRITE (6,610) 'SIGW01  ',(WSCAT(I,I,2),I=1,NGCOND)
            ENDIF
            IF(LWD) WRITE (6,610) 'LAMBDA-D',(WDLA(I),I=1,NDEL)
         ENDIF
      ENDIF
  240 CONTINUE
      CALL LCMSIX(IPWORK,' ',2)
      CALL LCMCL(IPWORK,2)
*----
*  UPDATE RECORDS ISOTOPESUSED, ISOTOPERNAME, ISOTOPESMIX, ETC.
*----
      CALL LCMGET(IPEDIT,'STATE-VECTOR',IPAR)
      JJISO=IPAR(2)
      IF(JJISO.GT.MAXISO) CALL XABORT('EDIRES: MAXISO OVERFLOW(2).')
      IF(JJISO.GT.0) THEN
         CALL LCMGET(IPEDIT,'ISOTOPESUSED',IHNISO)
         CALL LCMGET(IPEDIT,'ISOTOPERNAME',IHNIRF)
         CALL LCMGET(IPEDIT,'ISOTOPESDENS',SDEN)
         CALL LCMGET(IPEDIT,'ISOTOPESMIX',ISMIX)
         CALL LCMGET(IPEDIT,'ISOTOPESTYPE',ISTYP)
         CALL LCMGET(IPEDIT,'ISOTOPESTODO',ISTOD)
         CALL LCMGET(IPEDIT,'ISOTOPESVOL',VOLISO)
         CALL LCMGET(IPEDIT,'ISOTOPESTEMP',TNISO)
      ENDIF
      DO 260 INM=1,NMERGE
      IF(DENTOT(INM).GT.0.0) THEN
         JJISO=JJISO+1
         IF(JJISO.GT.MAXISO) CALL XABORT('EDIRES: MAXISO OVERFLOW(3).')
         WRITE(TEXT12,'(A8,I4.4)') '*MAC*RES',INM
         READ(TEXT12,'(3A4)') (IHNISO(I0,JJISO),I0=1,3)
         WRITE(TEXT12,'(A12)') '*MAC*RES    '
         READ(TEXT12,'(3A4)') (IHNIRF(I0,JJISO),I0=1,3)
         SDEN(JJISO)=1.0
         ISMIX(JJISO)=INM
         ISTYP(JJISO)=1
         ISTOD(JJISO)=1
         DVOL=0.0
         DO 250 IREGIO=1,NREGIO
         IF(IMERGE(IREGIO).EQ.INM) DVOL=DVOL+VOLUME(IREGIO)
  250    CONTINUE
         VOLISO(JJISO)=DVOL
         TNISO(JJISO)=TNTOT(INM)
      ENDIF
  260 CONTINUE
      IPAR(2)=JJISO
      IPAR(22)=IPAR(22)+1
      CALL LCMPUT(IPEDIT,'STATE-VECTOR',NSTATE,1,IPAR)
      CALL LCMPUT(IPEDIT,'ISOTOPESUSED',3*JJISO,3,IHNISO)
      CALL LCMPUT(IPEDIT,'ISOTOPERNAME',3*JJISO,3,IHNIRF)
      CALL LCMPUT(IPEDIT,'ISOTOPESDENS',JJISO,2,SDEN)
      CALL LCMPUT(IPEDIT,'ISOTOPESMIX',JJISO,1,ISMIX)
      CALL LCMPUT(IPEDIT,'ISOTOPESTYPE',JJISO,1,ISTYP)
      CALL LCMPUT(IPEDIT,'ISOTOPESTODO',JJISO,1,ISTOD)
      CALL LCMPUT(IPEDIT,'ISOTOPESVOL',JJISO,2,VOLISO)
      CALL LCMPUT(IPEDIT,'ISOTOPESTEMP',JJISO,2,TNISO)
      CALL LCMSIX(IPEDIT,' ',2) ! step down CURNAM
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(TNISO,VOLISO,ISTOD,ISTYP,ISMIX,IHNIRF,SDEN,IHNISO)
      DEALLOCATE(HMAKE)
      DEALLOCATE(PYIELD2,YIELD2,TNTOT,DAWR,DENTOT,WORK2,PNFIRA,WPY,
     1 WORK,GAS,WSCAT,WDLA)
      DEALLOCATE(ITYPRO,JPIFI,IEVOL2,LSIS2)
      RETURN
*
  600 FORMAT (//44H CROSS SECTION OF MERGED/CONDENSED ISOTOPE ',A12,
     1 7H' (ISO=,I8,2H):)
  610 FORMAT (/11H REACTION ',A12,2H':/(1X,1P,10E12.4))
      END
