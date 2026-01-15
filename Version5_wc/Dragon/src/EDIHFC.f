*DECK EDIHFC
      SUBROUTINE EDIHFC(IPEDIT,NGROUP,NGCOND,NREGIO,NMERGE,NBISO,
     >                  MATCOD,VOLUME,ISONAM,IPISO,MIX,FLUXES,DEN,
     >                  IGCOND,IMERGE,VOLME,IPRINT,EMEVF2)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover H-factors and normalize the flux.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): G. Marleau
*
*Parameters: input
* IPEDIT  pointer to the edition LCM object.
* NGROUP  number of groups.
* NGCOND  number of condensed groups.
* NREGIO  number of regions.
* NMERGE  number of merged regions.
* NBISO   number of isotopes.
* MATCOD  material per region.
* VOLUME  volume of region.
* ISONAM  isotopes names.
* IPISO   pointer array towards microlib isotopes.
* MIX     mixture associated with isotopes.
* FLUXES  multigroup fluxes.
* DEN     isotope density.
* IGCOND  limits of condensed groups.
* IMERGE  index of merged region.
* VOLME   merged volume.
* IPRINT  print level.
*
*Parameters: output
* EMEVF2  equivalent fission production energy by isotope.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPEDIT,IPISO(NBISO)
      INTEGER     IUNOUT
      INTEGER     NGROUP,NGCOND,NREGIO,NMERGE,NBISO,MATCOD(NREGIO),
     >            ISONAM(3,NBISO),MIX(NBISO),IGCOND(NGCOND),
     >            IMERGE(NREGIO)
      REAL        VOLUME(NREGIO),FLUXES(NREGIO,NGROUP),DEN(NBISO),
     >            EMEVF2(NBISO),VOLME(NMERGE)
      INTEGER     IPRINT
      DOUBLE PRECISION TOTPOW,POWF
      REAL, ALLOCATABLE, DIMENSION(:) :: SIG,HFACT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: FLXMER,WORK
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPEDIT,KPEDIT,KPLIB
      PARAMETER  (IUNOUT=6)
      CHARACTER   HSMG*131
      LOGICAL     LH
      DOUBLE PRECISION GAR,CONV,XDRCST,Z1,Z2
*----
*  SCRATCH STORAGE ALLOCATION
*   SIG     fission/capture cross sections.
*   HFACT   H-factor in a macrogroup.
*   FLXMER  merged and condensed flux.
*   WORK    H-factors.
*----
      ALLOCATE(SIG(NGROUP),HFACT(NMERGE))
      ALLOCATE(FLXMER(NMERGE,NGCOND),WORK(NMERGE,NGCOND))
*----
*  COMPUTE H-FACTOR
*----
      CONV=1.0D6 ! convert MeV to eV
      FLXMER(:NMERGE,:NGCOND)=0.0D0
      WORK(:NMERGE,:NGCOND)=0.0D0
      LH=.FALSE.
      DO 160 ISO=1,NBISO
        KPLIB=IPISO(ISO) ! set ISO-th isotope
        IF(.NOT.C_ASSOCIATED(KPLIB)) THEN
          WRITE(HSMG,'(17HEDIHFC: ISOTOPE '',3A4,16H'' IS NOT AVAILAB,
     >    19HLE IN THE MICROLIB.)') (ISONAM(I0,ISO),I0=1,3)
          CALL XABORT(HSMG)
        ENDIF
        Z1=0.0D0
        Z2=0.0D0
        EMEVF2(ISO)=0.0
*----
*  RECOVER H-FACTOR INFORMATION IF AVAILABLE
*----
        CALL LCMLEN(KPLIB,'H-FACTOR',ILLCM,ITLCM)
        IF(ILLCM.EQ.0) GO TO 160
        LH=.TRUE.
        CALL LCMGET(KPLIB,'H-FACTOR',SIG)
        DO 90 IREG=1,NREGIO
          IMR=IMERGE(IREG)
          IF((IMR.GT.0).AND.(MATCOD(IREG).EQ.MIX(ISO))) THEN
            IGRFIN=0
            DO 80 IGC=1,NGCOND
              IGRDEB=IGRFIN+1
              IGRFIN=IGCOND(IGC)
              GAR=0.0D0
              DO 70 IGR=IGRDEB,IGRFIN
                GAR=GAR+FLUXES(IREG,IGR)*DEN(ISO)*VOLUME(IREG)*SIG(IGR)
  70          CONTINUE
              WORK(IMR,IGC)=WORK(IMR,IGC)+GAR
              Z1=Z1+GAR
  80        CONTINUE
          ENDIF
  90    CONTINUE
*----
*  COMPUTE FISSION ENERGY
*----
        CALL LCMLEN(KPLIB,'NFTOT',ILLCM,ITLCM)
        IF(ILLCM.EQ.NGROUP) THEN
          CALL LCMGET(KPLIB,'NFTOT',SIG)
          DO 120 IREG=1,NREGIO
            IMR=IMERGE(IREG)
            IF((IMR.GT.0).AND.(MATCOD(IREG).EQ.MIX(ISO))) THEN
              IGRFIN=0
              DO 110 IGC=1,NGCOND
                IGRDEB=IGRFIN+1
                IGRFIN=IGCOND(IGC)
                DO 100 IGR=IGRDEB,IGRFIN
                  Z2=Z2+FLUXES(IREG,IGR)*DEN(ISO)*VOLUME(IREG)*SIG(IGR)
 100            CONTINUE
 110          CONTINUE
            ENDIF
 120      CONTINUE
          IF(Z2.NE.0.0) EMEVF2(ISO)=REAL(Z1/Z2)
        ENDIF
 160  CONTINUE
*----
*  Normalize total power to 1 W
*  Print fission, capture and total power density
*----
      TOTPOW=0.0D0
      DO IGC=1,NGCOND
        DO IMR=1,NMERGE
          TOTPOW=TOTPOW+WORK(IMR,IGC)*XDRCST('eV','J')
        ENDDO
      ENDDO
      IF(TOTPOW.GT.0.0D0) THEN
        IF(ABS(IPRINT).GE.2) THEN
          WRITE(IUNOUT,6000)
          DO IMR=1,NMERGE
            POWF=0.0D0
            DO IGC=1,NGCOND
              POWF=POWF+WORK(IMR,IGC)
            ENDDO
            IF(VOLME(IMR).NE.0.0) THEN
              POWF=POWF/(TOTPOW*VOLME(IMR))
              WRITE(IUNOUT,6001) IMR,VOLME(IMR),POWF
            ENDIF
          ENDDO
        ENDIF
      ENDIF
*----
*  COMPUTE THE HOMOGENIZED/CONDENSED FLUX
*----
      IF(LH) THEN
        DO 190 IREG=1,NREGIO
          IMR=IMERGE(IREG)
          IF(IMR.GT.0) THEN
            IGRFIN=0
            DO 180 IGC=1,NGCOND
              IGRDEB=IGRFIN+1
              IGRFIN=IGCOND(IGC)
              GAR=0.0D0
              DO 170 IGR=IGRDEB,IGRFIN
                GAR=GAR+FLUXES(IREG,IGR)*VOLUME(IREG)
 170          CONTINUE
              FLXMER(IMR,IGC)=FLXMER(IMR,IGC)+GAR
 180        CONTINUE
          ENDIF
 190    CONTINUE
        DO 210 IGC=1,NGCOND
          DO 200 IMR=1,NMERGE
            IF(FLXMER(IMR,IGC).GT.0.0) THEN
              WORK(IMR,IGC)=WORK(IMR,IGC)/FLXMER(IMR,IGC)
            ENDIF
 200      CONTINUE
 210    CONTINUE
*----
*  SAVE ON LCM
*----
        CALL LCMSIX(IPEDIT,'MACROLIB',1)
        JPEDIT=LCMLID(IPEDIT,'GROUP',NGCOND)
        DO 230 IGC=1,NGCOND
          DO 220 IMR=1,NMERGE
            HFACT(IMR)=REAL(WORK(IMR,IGC))
 220      CONTINUE
          KPEDIT=LCMDIL(JPEDIT,IGC)
          CALL LCMPUT(KPEDIT,'H-FACTOR',NMERGE,2,HFACT)
 230    CONTINUE
        CALL LCMSIX(IPEDIT,' ',2)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(WORK,FLXMER)
      DEALLOCATE(HFACT,SIG)
      RETURN
*----
*  FORMAT
*----
 6000 FORMAT(/' EDIHFC: POWER DENSITY (W/cc) NORMALIZED TO 1 W TOTAL ',
     > 'POWER '/' REGION',6X,'VOLUME',7X,'FISSION')
 6001 FORMAT(1X,I4,1P,2E14.5)
      END
