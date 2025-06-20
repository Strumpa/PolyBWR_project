*DECK ACRLIB
      SUBROUTINE ACRLIB(MAXNIS,MAXISO,IPLIB,IPAPX,IACCS,NMIX,NGRP,IMPX,
     1 HEQUI,NCAL,ITER,MY1,MY2,MD1,MD2,TERP,NISO,LISO,HISO,CONC,ITODO,
     2 MIXC,LRES,LPURE,LTOTAL,ILUPS,B2,LFROM,VTOT,YLDS,DECAYC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Build the Microlib by scanning the NCAL elementary calculations in
* a Apex file and weighting them with TERP factors.
*
*Copyright:
* Copyright (C) 2021 Ecole Polytechnique de Montreal
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* MAXNIS  maximum value of NISO(I) in user data.
* MAXISO  maximum allocated space for output Microlib TOC information.
* IPLIB   address of the output Microlib LCM object.
* IPAPX   pointer to the Apex file.
* IACCS   =0 Microlib is created; =1 ... is updated.
* NMIX    maximum number of material mixtures in the Microlib.
* NGRP    number of energy groups.
* IMPX    print parameter (equal to zero for no print).
* HEQUI   keyword of SPH-factor set to be recovered.
* NCAL    number of elementary calculations in the Apex file.
* ITER    completion flag (=0: compute the macrolib).
* MY1     number of fissile isotopes including macroscopic sets.
* MY2     number of fission fragment.
* MD1     number of types of radioactive decay reactions.
* MD2     number of particularized isotopes including macro.
* TERP    interpolation factors.
* NISO    number of user-selected isotopes.
* LISO    type of treatment (=.true.: ALL; =.false.: ONLY).
* HISO    name of the user-selected isotopes.
* CONC    user-defined number density of the user-selected isotopes. A
*         value of -99.99 is set to indicate that the Apex file value is
*         used.
* ITODO   non-depletion mask (=1 to force a user-selected isotope to be
*         non-depleting)
* MIXC    mixture index in the Apex file corresponding to each Microlib
*         mixture. Equal to zero if a Microlib mixture is not updated.
* LRES    =.true. if the interpolation is done without updating isotopic 
*         densities
* LPURE   =.true. if the interpolation is a pure linear interpolation 
*         with TERP factors.
* LTOTAL =.true. to use the mac/TOTAL macroscopic set.
* ILUPS   up-scattering removing flag (=1 to remove up-scattering from
*         output cross-sections).
* B2      buckling
* LFROM   macroregion flag (=.true. if 'xs       n' groups are set).
* VTOT    volume of updated core.
* YLDS    fission yields.
* DECAYC  radioactive decay constants.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE hdf5_wrap
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB,IPAPX
      INTEGER MAXNIS,MAXISO,IACCS,NMIX,NGRP,IMPX,NCAL,ITER,MY1,MY2,MD1,
     1 MD2,NISO(NMIX),ITODO(NMIX,MAXNIS),MIXC(NMIX),ILUPS
      REAL TERP(NCAL,NMIX),CONC(NMIX,MAXNIS),B2
      DOUBLE PRECISION VTOT,YLDS(MY1,MY2),DECAYC(MD1,MD2)
      LOGICAL LISO(NMIX),LRES,LPURE,LTOTAL,LFROM
      CHARACTER(LEN=80) HEQUI
      CHARACTER(LEN=8) HISO(NMIX,MD2)
*----
*  LOCAL VARIABLES
*----
      INTEGER, PARAMETER::IOUT=6
      INTEGER, PARAMETER::MAXMAC=2
      INTEGER, PARAMETER::MAXREA=50
      INTEGER, PARAMETER::NSTATE=40
      INTEGER, PARAMETER::MAXFRD=4
      TYPE(C_PTR) JPLIB,KPLIB
      REAL B2APEX, FACT0, WEIGHT
      INTEGER I, I0, K, IBM, IBMOLD, ICAL, ID1, IED2, IFISS, IGR,
     & ILONG, IMAC, IOF, IPRC, IREA, IREAF, IRES, ISO, ITRANC, ITSTMP,
     & ITYLCM, IY1, IY2, JSO,  KSO, KSO1, LMY1, LSO, MAXMIX, NBISO,
     & NBISO1, NBISO2, NBISO2I, NBS1, NCALS, NED2, NL, NLAM, NBMAC,
     & NMIL, NPAR, NPRC, NREA, NSURFD, NISOF, NISOP, NISOS, NISOTS,
     & NVP, RANK, NBYTE, TYPE, ISURF, DIMSR(5)
      CHARACTER RECNAM*80,TEXT8*8, TEXT12*12,HSMG*131,HVECT2(MAXREA)*8,
     1 HRESID*8,HHAD(MAXFRD)*16
      INTEGER ISTATE(NSTATE),INAME(2),IHRES(2)
      REAL TMPDAY(3)
      LOGICAL LUSER,LSTRD
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IMIX2,ITOTM,IRESM,ISONA,
     1 ISOMI,ITOD2,ISTY1,ISTY2,IPIFI,IMICR,ITOD1,JJSO,IPYMIX,DIMS_APX
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: HUSE2,HNAM2
      REAL, ALLOCATABLE, DIMENSION(:) :: DENS2,DENS3,VOL2,VOLMI2,SPH,
     1 ENER,XVOLM,CONCE,TAUXFI,NWT0,FLUXS,DENIS,GAR1,GAR2,LAMB,BETAR,
     2 INVELS,BETARB,INVELSB
      REAL, ALLOCATABLE, DIMENSION(:,:) :: ADF,DENS1,FACT,DECAY2,
     1 CHIRS,CHIRSB
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: XS,SIGS,DENS0,FLUX,ADF2,
     1 YLDS2
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: SS2D
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: YLDSM
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: LXS,MASK,MASKL
      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:) :: HADF,HNOMIS,NOMISO,
     1 NOMMAC,HPYNAM
      CHARACTER(LEN=12), ALLOCATABLE, DIMENSION(:) :: NOMREA
*----
*  RECOVER APEX FILE CHARACTERISTICS
*----
      I=0
      CALL APXTOC(IPAPX,0,NLAM,NREA,NBISO,NBMAC,NMIL,NPAR,NVP,NISOF,
     1 NISOP,NISOS,NCALS,I,NISOTS,NSURFD,NPRC)
      IF(NGRP.NE.I) CALL XABORT('ACRLIB: INVALID VALUE OF NGRP.')
      IF(NREA.GT.MAXREA) CALL XABORT('ACRLIB: MAXREA OVERFLOW')
      IF(NBMAC.GT.MAXMAC) CALL XABORT('ACRLIB: MAXMAC OVERFLOW')
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(IMIX2(MAXISO),ITOD2(MAXISO),ISTY1(MAXISO),ISTY2(MAXISO),
     1 HUSE2(3,MAXISO),HNAM2(3,MAXISO))
      ALLOCATE(DENS2(MAXISO),DENS3(MAXISO),VOL2(MAXISO),VOLMI2(NMIX),
     1 FLUX(NMIX,NGRP,2),SPH(NGRP))
      ALLOCATE(HADF(NSURFD),ADF2(NMIX,NGRP,NSURFD))
*----
*  MICROLIB INITIALIZATION
*----
      VOLMI2(:NMIX)=0.0
      DENS2(:MAXISO)=0.0
      VOL2(:MAXISO)=0.0
      IMIX2(:MAXISO)=0
      ITOD2(:MAXISO)=0
      ISTY2(:MAXISO)=0
      IF(NSURFD.GT.0) ADF2(:NMIX,:NGRP,:NSURFD)=0.0
      IF(IACCS.EQ.0) THEN
         IF(LRES) CALL XABORT('ACRLIB: RES OPTION IS INVALID.')
         NBISO2=0
         NED2=0
         TEXT12='L_LIBRARY'
         CALL LCMPTC(IPLIB,'SIGNATURE',12,TEXT12)
      ELSE
         CALL LCMGET(IPLIB,'STATE-VECTOR',ISTATE)
         IF(ISTATE(1).NE.NMIX) CALL XABORT('ACRLIB: INVALID NUMBER OF '
     1   //'MATERIAL MIXTURES IN THE MICROLIB.')
         IF(ISTATE(3).NE.NGRP) CALL XABORT('ACRLIB: INVALID NUMBER OF '
     1   //'ENERGY GROUPS IN THE MICROLIB.')
         NBISO2=ISTATE(2)
         IF(NBISO2.GT.MAXISO) CALL XABORT('ACRLIB: MAXISO OVERFLOW(1).')
         NED2=ISTATE(13)
         IF(NED2.GT.MAXREA) CALL XABORT('ACRLIB: MAXREA OVERFLOW.')
         CALL LCMLEN(IPLIB,'MIXTURESVOL',ILONG,ITYLCM)
         IF(ILONG.GT.0) THEN
           CALL LCMGET(IPLIB,'MIXTURESVOL',VOLMI2)
         ELSE
           VOLMI2(:NMIX)=0.0
         ENDIF
         CALL LCMGET(IPLIB,'ISOTOPESUSED',HUSE2)
         CALL LCMGET(IPLIB,'ISOTOPERNAME',HNAM2)
         CALL LCMGET(IPLIB,'ISOTOPESDENS',DENS2)
         CALL LCMGET(IPLIB,'ISOTOPESVOL',VOL2)
         CALL LCMGET(IPLIB,'ISOTOPESMIX',IMIX2)
         CALL LCMGET(IPLIB,'ISOTOPESTODO',ITOD2)
         CALL LCMGET(IPLIB,'ISOTOPESTYPE',ISTY2)
         IF(NED2.GT.0) CALL LCMGTC(IPLIB,'ADDXSNAME-P0',8,NED2,HVECT2)
         IF(NSURFD.GT.0) THEN
           CALL LCMSIX(IPLIB,'MACROLIB',1)
           CALL LCMLEN(IPLIB,'ADF',ILONG,ITYLCM)
           IF(ILONG.EQ.0) THEN
             CALL LCMLIB(IPLIB)
             CALL XABORT('ACRLIB: UNABLE TO FIND DIRECTORY ADF.')
           ENDIF
           CALL LCMSIX(IPLIB,'ADF',1)
           CALL LCMGTC(IPLIB,'HADF',8,NSURFD,HADF)
           DO I=1,NSURFD
             CALL LCMGET(IPLIB,HADF(I),ADF2(1,1,I))
           ENDDO
           CALL LCMSIX(IPLIB,' ',2)
           CALL LCMSIX(IPLIB,' ',2)
         ENDIF
      ENDIF
*----
*  RECOVER INFORMATION FROM physconst GROUP.
*----
      IF(hdf5_group_exists(IPAPX,"/physconst/")) THEN
        CALL hdf5_read_data(IPAPX,"/physconst/ENRGS",ENER)
      ELSE IF(hdf5_group_exists(IPAPX,"/physco001/")) THEN
        CALL hdf5_read_data(IPAPX,"/physco001/ENRGS",ENER)
      ELSE
        CALL XABORT('ACRLIB: GROUP physconst NOT FOUND IN HDF5 FILE.')
      ENDIF
      DO IGR=1,NGRP+1
        ENER(IGR)=ENER(IGR)/1.0E-6
      ENDDO
      CALL LCMPUT(IPLIB,'ENERGY',NGRP+1,2,ENER)
      DO IGR=1,NGRP
        ENER(IGR)=LOG(ENER(IGR)/ENER(IGR+1))
      ENDDO
      CALL LCMPUT(IPLIB,'DELTAU',NGRP,2,ENER)
      DEALLOCATE(ENER)
*----
*  RECOVER INFORMATION FROM explicit GROUP.
*----
      ALLOCATE(ITOTM(NMIL),IRESM(NMIL))
      ITOTM(:)=0
      IRESM(:)=0
      IREAF=0
      IF(NREA.GT.0) THEN
        IF(hdf5_group_exists(IPAPX,"/explicit/")) THEN
          CALL hdf5_read_data(IPAPX,"/explicit/REANAME",NOMREA)
        ELSE IF(hdf5_group_exists(IPAPX,"/expli001/")) THEN
          CALL hdf5_read_data(IPAPX,"/expli001/REANAME",NOMREA)
        ELSE
          CALL XABORT('ACRLIB: GROUP explicit NOT FOUND IN HDF5 FILE.')
        ENDIF
        IF(IMPX.GT.1) THEN
          WRITE(IOUT,'(29H ACRLIB: Available reactions:/(1X,10A13))')
     1    (NOMREA(I),I=1,NREA)
        ENDIF
        DO IREA=1,NREA
          IF(NOMREA(IREA).EQ.'NUFI') THEN
            IREAF=IREA
            EXIT
          ENDIF
        ENDDO
      ENDIF
      IF(NBISO.GT.0) THEN
        IF(hdf5_group_exists(IPAPX,"/explicit/")) THEN
          CALL hdf5_read_data(IPAPX,"/explicit/ISONAME",NOMISO)
        ELSE IF(hdf5_group_exists(IPAPX,"/expli001/")) THEN
          CALL hdf5_read_data(IPAPX,"/expli001/ISONAME",NOMISO)
        ELSE
          CALL XABORT('ACRLIB: GROUP explicit NOT FOUND IN HDF5 FILE.')
        ENDIF
      ENDIF
      IF(LTOTAL.AND.(NBMAC.EQ.0)) CALL XABORT('ACRLIB: NBMAC=0.')
      IF(NBMAC.GT.0) THEN
        IF(hdf5_group_exists(IPAPX,"/explicit/")) THEN
          CALL hdf5_read_data(IPAPX,"/explicit/MACNAME",NOMMAC)
        ELSE IF(hdf5_group_exists(IPAPX,"/expli001/")) THEN
          CALL hdf5_read_data(IPAPX,"/expli001/MACNAME",NOMMAC)
        ELSE
          CALL XABORT('ACRLIB: GROUP explicit NOT FOUND IN HDF5 FILE.')
        ENDIF
        DO I=1,NBMAC
          IF(NOMMAC(I).EQ.'TOTAL') ITOTM(:)=I
          IF(NOMMAC(I).EQ.'RESIDUAL') IRESM(:)=I
        ENDDO
        NBISO1=NBISO+NBMAC
        ALLOCATE(HNOMIS(NBISO1))
        IF(NBISO.GT.0) HNOMIS(:NBISO)=NOMISO(:NBISO)
        IF(LTOTAL) THEN
          ! use the mac/TOTAL macroscopic set
          HNOMIS(:NBISO)=' '
          DO I=1,NBMAC
            IF(NOMMAC(I).EQ.'TOTAL') THEN
              HNOMIS(NBISO+I)='*MAC*RES'
            ELSE IF(NOMMAC(I).EQ.'RESIDUAL') THEN
              HNOMIS(NBISO+I)=' '
            ENDIF
          ENDDO
        ELSE
          ! use the mac/RESIDUAL macroscopic set
          DO I=1,NBMAC
            IF(NOMMAC(I).EQ.'TOTAL') THEN
              HNOMIS(NBISO+I)=' '
            ELSE IF(NOMMAC(I).EQ.'RESIDUAL') THEN
              HNOMIS(NBISO+I)='*MAC*RES'
            ENDIF
          ENDDO
        ENDIF
      ELSE
        NBISO1=NBISO
        ALLOCATE(HNOMIS(NBISO1))
        IF(NBISO.GT.0) HNOMIS(:NBISO)=NOMISO(:NBISO)
      ENDIF
*----
*  RECOVER VOLUMES.
*----
      ALLOCATE(XVOLM(NMIL))
      RECNAM='calc       1/xs/'
      DO IBMOLD=1,NMIL
        IF(LFROM) WRITE(RECNAM,'(4Hcalc,I8,3H/xs,I8,1H/)') ICAL,IBMOLD
        CALL hdf5_info(IPAPX,TRIM(RECNAM)//"MEDIA_VOLUME",RANK,
     1  TYPE,NBYTE,DIMSR)
        IF(TYPE.EQ.99) THEN
          XVOLM(IBMOLD)=1.0
          WRITE(IOUT,'(44H ACRLIB: WARNING -- Record MEDIA_VOLUME is m,
     1    42Hissing in the Apex file. Volume set to 1.0)')
        ELSE
          CALL hdf5_read_data(IPAPX,TRIM(RECNAM)//"MEDIA_VOLUME",
     1    XVOLM(IBMOLD))
        ENDIF
      ENDDO
*----
*  FIND SCATTERING ANISOTROPY.
*----
      CALL hdf5_info(IPAPX,TRIM(RECNAM)//"mac/TOTAL/DIFF",RANK,
     1 TYPE,NBYTE,DIMSR)
      IF(TYPE.EQ.99) CALL XABORT('ACRLIB: MISSING SCATTERING INFO.')
      NL=DIMSR(2)
      IF(IMPX.GT.1) THEN
        WRITE(IOUT,'(36H ACRLIB: number of Legendre orders =,I4)') NL
      ENDIF
*----
*  LOOP OVER APEX MIXTURES TO COMPUTE DENS0(NMIL,NCAL,NBISO1)
*----
      ALLOCATE(DENS0(NMIL,NCAL,NBISO1))
      DENS0(:NMIL,:NCAL,:NBISO1)=0.0
      DO 30 IBMOLD=1,NMIL
      DO ICAL=1,NCAL
        DO IBM=1,NMIX
          IF((TERP(ICAL,IBM).NE.0.0).AND.(MIXC(IBM).EQ.IBMOLD)) GO TO 10
        ENDDO
        CYCLE
   10   WRITE(RECNAM,'(4Hcalc,I8,4H/xs/)') ICAL
        IF(LFROM) WRITE(RECNAM,'(4Hcalc,I8,3H/xs,I8,1H/)') ICAL,IBMOLD
        IF(NBISO.GT.0) THEN
          CALL hdf5_read_data(IPAPX,TRIM(RECNAM)//"mic/CONC",CONCE)
          DO 20 ISO=1,NBISO
            DENS0(IBMOLD,ICAL,ISO)=CONCE(ISO)
   20     CONTINUE
        ENDIF
      ENDDO
   30 CONTINUE
      IF(NBISO.GT.0) DEALLOCATE(CONCE)
*----
*  LOOP OVER MICROLIB MIXTURES
*----
      YLDS(:MY1,:MY2)=0.0D0
      DECAYC(:MD1,:MD2)=0.0D0
      VTOT=0.0D0
      DO 40 IBM=1,NMIX
      IBMOLD=MIXC(IBM)
      IF(IBMOLD.NE.0) VTOT=VTOT+XVOLM(IBMOLD)
   40 CONTINUE
      ALLOCATE(JJSO(NBISO+NBMAC),YLDSM(MY1,MY2),ITOD1(NBISO1))
      ALLOCATE(TAUXFI(NBISO+NBMAC),NWT0(NGRP),
     1 SIGS(NGRP,NL,NBISO+NBMAC),SS2D(NGRP,NGRP,NL,NBISO+NBMAC),
     2 XS(NGRP,NREA,NBISO+NBMAC))
      ALLOCATE(LXS(NREA))
      ALLOCATE(CHIRS(NGRP,NPRC),BETAR(NPRC),INVELS(NGRP))
      CHIRS(:NGRP,:NPRC)=0.0
      BETAR(:NPRC)=0.0
      INVELS(:NGRP)=0.0
      ALLOCATE(BETARB(NPRC),INVELSB(NGRP))
      ALLOCATE(DENS1(NBISO1,NCAL),FACT(NBISO1,NCAL))
      JPLIB=LCMLID(IPLIB,'ISOTOPESLIST',(NBISO+NBMAC)*NMIX)
*
      DO 180 IBM=1,NMIX
      IBMOLD=MIXC(IBM)
      IF(IBMOLD.EQ.0) GO TO 180
      IF(NISO(IBM).GT.MAXNIS) CALL XABORT('ACRLIB: MAXNIS OVERFLOW.')
      VOLMI2(IBM)=XVOLM(IBMOLD)
      IMAC=ITOTM(IBMOLD)
      IRES=IRESM(IBMOLD)
*----
*  RECOVER ITOD1(NBISO1) INDICES.
*----
      ITOD1(:NBISO1)=0
      DO 50 ISO=1,NBISO1 ! Apex file isotope
        DO KSO=1,NISO(IBM) ! user-selected isotope
          IF(NOMISO(ISO).EQ.HISO(IBM,KSO)) THEN
            ITOD1(ISO)=ITODO(IBM,KSO)
            GO TO 50
          ENDIF
        ENDDO
   50 CONTINUE
*----
*  COMPUTE THE NUMBER DENSITIES OF EACH ELEMENTARY CALCULATION.
*----
      DENS1(:NBISO1,:NCAL)=0.0
      DENS3(:NBISO1)=0.0
      DO ICAL=1,NCAL
        WEIGHT=TERP(ICAL,IBM)
        IF(WEIGHT.EQ.0.0) CYCLE
        DO ISO=1,NBISO
          LUSER=.FALSE.
          KSO1=0
          DO KSO=1,NISO(IBM) ! user-selected isotope
            IF(NOMISO(ISO).EQ.HISO(IBM,KSO)) THEN
              KSO1=KSO
              LUSER=(CONC(IBM,KSO1).NE.-99.99)
              GO TO 60
            ENDIF
          ENDDO
   60     IF(LUSER) THEN
            DENS1(ISO,ICAL)=CONC(IBM,KSO1)
            CYCLE
          ENDIF
          IF(.NOT.LISO(IBM)) CYCLE
          DENS1(ISO,ICAL)=DENS0(IBMOLD,ICAL,ISO)
        ENDDO
        IF((NBISO.NE.0).AND.(.NOT.LTOTAL)) THEN
          DENS1(NBISO+IRES,ICAL)=1.0
        ELSE IF(IMAC.NE.0) THEN
          DENS1(NBISO+IMAC,ICAL)=1.0
        ENDIF
        DO ISO=1,NBISO1
          DENS3(ISO)=DENS3(ISO)+WEIGHT*DENS1(ISO,ICAL)
        ENDDO
      ENDDO
      FACT(:NBISO1,:NCAL)=1.0
      IF(.NOT.LPURE) THEN
        DO ICAL=1,NCAL
          IF(TERP(ICAL,IBM).EQ.0.0) CYCLE
          DO ISO=1,NBISO1
            IF(DENS3(ISO).GT.DENS1(ISO,ICAL)*1.0E-9) THEN
              FACT(ISO,ICAL)=DENS1(ISO,ICAL)/DENS3(ISO)
            ENDIF
          ENDDO
        ENDDO
      ENDIF
*----
*  INITIALIZE WORKING ARRAYS.
*----
      TAUXFI(:NBISO1)=0.0
      NWT0(:NGRP)=0.0
      SIGS(:NGRP,:NL,:NBISO1)=0.0
      SS2D(:NGRP,:NGRP,:NL,:NBISO1)=0.0
      XS(:NGRP,:NREA,:NBISO1)=0.0
      LXS(:NREA)=.FALSE.
      YLDSM(:MY1,:MY2)=0.0D0
*----
*  MAIN LOOP OVER ELEMENTARY CALCULATIONS
*----
      TEXT12='*MAC*RES'
      READ(TEXT12,'(2A4)') IHRES(1),IHRES(2)
      LSTRD=.FALSE.
      B2APEX=B2
      DO 80 ICAL=1,NCAL
      WEIGHT=TERP(ICAL,IBM)
      IF(WEIGHT.EQ.0.0) GO TO 80
*----
*  RECOVER INFORMATION FROM caldir GROUP.
*----
      WRITE(RECNAM,'(4Hcalc,I8,10H/kinetics/)') ICAL
      CALL hdf5_info(IPAPX,TRIM(RECNAM)//"LAMBDA",RANK,TYPE,NBYTE,DIMSR)
      NPRC=0
      IF(TYPE.NE.99) THEN
        NPRC=DIMSR(1)
        CALL hdf5_read_data(IPAPX,TRIM(RECNAM)//"LAMBDA",LAMB)
        CALL hdf5_read_data(IPAPX,TRIM(RECNAM)//"CHIDA",CHIRSB)
        CALL hdf5_read_data(IPAPX,TRIM(RECNAM)//"BETADA",BETARB)
        CALL hdf5_read_data(IPAPX,TRIM(RECNAM)//"INVELA",INVELSB)
      ENDIF
*----
*  SELECT APEX MIXTURE IBMOLD.
*----
      WRITE(RECNAM,'(4Hcalc,I8,4H/xs/)') ICAL
      IF(LFROM) WRITE(RECNAM,'(4Hcalc,I8,3H/xs,I8,1H/)') ICAL,IBMOLD
      IF(HEQUI.NE.' ') THEN
        CALL hdf5_read_data(IPAPX,TRIM(RECNAM)//"MEDIA_SPH/"//HEQUI,SPH)
      ELSE
        SPH(:NGRP)=1.0
      ENDIF
      CALL hdf5_read_data(IPAPX,TRIM(RECNAM)//"FLUX",FLUXS)
      DO I=1,NGRP
        FLUXS(I)=FLUXS(I)/XVOLM(IBMOLD)
        NWT0(I)=NWT0(I)+WEIGHT*FLUXS(I)/SPH(I)
      ENDDO
      IF((NBISO.NE.0).AND.(.NOT.LTOTAL)) THEN
        DO ISO=1,NBISO
          FACT0=FACT(ISO,ICAL)
          CALL ACRSX2(IPAPX,RECNAM,NREA,NGRP,NISOF,NISOP,NL,ISO,
     1    NOMREA,B2APEX,FACT0,WEIGHT,SPH,FLUXS,IREAF,LPURE,LXS,
     2    XS(1,1,ISO),SIGS(1,1,ISO),SS2D(1,1,1,ISO),TAUXFI(ISO))
        ENDDO
        IF(IRES.NE.0) THEN
          FACT0=1.0
          CALL ACRSX2(IPAPX,RECNAM,NREA,NGRP,NISOF,NISOP,NL,-2,NOMREA,
     1    B2APEX,FACT0,WEIGHT,SPH,FLUXS,IREAF,LPURE,LXS,XS(1,1,NBISO+1),
     2    SIGS(1,1,NBISO+1),SS2D(1,1,1,NBISO+1),TAUXFI(NBISO+IRES))
        ENDIF
      ELSE IF(IMAC.NE.0) THEN
        FACT0=1.0
        CALL ACRSX2(IPAPX,RECNAM,NREA,NGRP,NISOF,NISOP,NL,-1,NOMREA,
     1  B2APEX,FACT0,WEIGHT,SPH,FLUXS,IREAF,LPURE,LXS,XS(1,1,NBISO+1),
     2  SIGS(1,1,NBISO+1),SS2D(1,1,1,NBISO+1),TAUXFI(NBISO+IMAC))
      ELSE
        CALL XABORT('ACRLIB: NO MACROSCOPIC SET.')
      ENDIF
      DEALLOCATE(FLUXS)
*
      IF(NPRC.GT.0) THEN
        DO IGR=1,NGRP
          INVELS(IGR)=INVELS(IGR)+SPH(IGR)*WEIGHT*INVELSB(IGR)
          DO IPRC=1,NPRC
            CHIRS(IGR,IPRC)=CHIRS(IGR,IPRC)+WEIGHT*CHIRSB(IGR,IPRC)
          ENDDO
        ENDDO
        DO IPRC=1,NPRC
          BETAR(IPRC)=BETAR(IPRC)+WEIGHT*BETARB(IPRC)
        ENDDO
      ENDIF
*----
*  COMPUTE DEPLETION CHAIN DATA
*----
      IF(NISOF*NISOP.GT.0) THEN
        IF(hdf5_group_exists(IPAPX,"/physconst/")) THEN
          CALL hdf5_read_data(IPAPX,"/physconst/FYIELDS",YLDS2)
        ELSE IF(hdf5_group_exists(IPAPX,"/physco001/")) THEN
          CALL hdf5_read_data(IPAPX,"/physco001/FYIELDS",YLDS2)
        ELSE
          CALL XABORT('ACRLIB: GROUP physconst NOT FOUND IN HDF5 FILE.')
        ENDIF
        DO IY2=1,NISOP
          DO IY1=1,NISOF
            YLDSM(IY1,IY2)=YLDSM(IY1,IY2)+WEIGHT*YLDS2(IY1,IY2,1)
            YLDS(IY1,IY2)=YLDS(IY1,IY2)+WEIGHT*YLDS2(IY1,IY2,1)*
     >      VOLMI2(IBM)/VTOT
          ENDDO
        ENDDO
        DEALLOCATE(YLDS2)
      ENDIF
      IF((MD1*MD2.GT.0).AND.(NBISO.GT.0)) THEN
        IF(hdf5_group_exists(IPAPX,"/physconst/")) THEN
          CALL hdf5_read_data(IPAPX,"/physconst/DECAYC",DECAY2)
        ELSE IF(hdf5_group_exists(IPAPX,"/physco001/")) THEN
          CALL hdf5_read_data(IPAPX,"/physco001/DECAYC",DECAY2)
        ELSE
          CALL XABORT('ACRLIB: GROUP physconst NOT FOUND IN HDF5 FILE.')
        ENDIF
        DO ISO=1,NBISO
          DO ID1=1,NLAM
            DECAYC(ID1,ISO)=DECAYC(ID1,ISO)+WEIGHT*DECAY2(ID1,ISO)*
     >      VOLMI2(IBM)/VTOT
          ENDDO
        ENDDO
        DEALLOCATE(DECAY2)
      ENDIF
   80 CONTINUE ! end of loop over elementary calculations.
*----
*  IDENTIFY SPECIAL FLUX EDITS
*----
      DO IREA=1,NREA
        IF(NOMREA(IREA).EQ.'ABSO') THEN
          DO 90 IED2=1,NED2
          IF(HVECT2(IED2).EQ.'ABSO') GO TO 100
   90     CONTINUE
          NED2=NED2+1
          IF(NED2.GT.MAXREA) CALL XABORT('ACRLIB: MAXREA OVERFLOW(1).')
          HVECT2(NED2)='ABSO'
        ELSE IF(NOMREA(IREA).EQ.'FISS') THEN
          DO 95 IED2=1,NED2
          IF(HVECT2(IED2).EQ.'NFTOT') GO TO 100
   95     CONTINUE
          NED2=NED2+1
          IF(NED2.GT.MAXREA) CALL XABORT('ACRLIB: MAXREA OVERFLOW(2).')
          HVECT2(NED2)='NFTOT'
        ENDIF
  100   CONTINUE
      ENDDO
*----
*  SET FLAG LSTRD
*----
      LSTRD=.TRUE.
      DO IREA=1,NREA
        IF(NOMREA(IREA).EQ.'LEAK') THEN
          IF(LXS(IREA).AND.(B2APEX.NE.0.0)) LSTRD=.FALSE.
          EXIT
        ENDIF
      ENDDO
*----
*  SAVE CROSS SECTIONS IN MICROLIB FOR MIXTURE IBM
*----
      ISTY1(:NBISO1)=0
      JJSO(:NBISO1)=0
      NBISO2I=NBISO2
      IF((NBISO.NE.0).AND.(.NOT.LTOTAL)) THEN
        HRESID=' '
        DO ISO=1,NBISO
          READ(HNOMIS(ISO),'(2A4)') INAME(:2)
          CALL SCRFND(MAXISO,NBISO2I,NBISO2,INAME,IBM,HRESID,HUSE2,
     1    HNAM2,IMIX2,JJSO(ISO))
          KPLIB=LCMDIL(JPLIB,JJSO(ISO)) ! step up isot JJSO(ISO)
          CALL ACRISO(KPLIB,NREA,NGRP,NL,NPRC,NOMREA,NWT0,XS(1,1,ISO),
     1    SIGS(1,1,ISO),SS2D(1,1,1,ISO),TAUXFI(ISO),LXS,LAMB,CHIRS,
     2    BETAR,INVELS,INAME,LSTRD,LPURE,ILUPS,ITRANC,IFISS)
          IF(MY1*MY2.GT.0) CALL ACRNDF(IMPX,NBISO+NBMAC,ISO,IBM,HNOMIS,
     1    IPAPX,KPLIB,MY1,MY2,YLDSM,ISTY1(ISO))
        ENDDO
        IF(IRES.NE.0) THEN
          HRESID=NOMMAC(IRES)
          CALL SCRFND(MAXISO,NBISO2I,NBISO2,IHRES,IBM,HRESID,HUSE2,
     1    HNAM2,IMIX2,JJSO(NBISO+IRES))
          KPLIB=LCMDIL(JPLIB,JJSO(NBISO+IRES)) ! step up isot JJSO(NBISO+IRES)
          CALL ACRISO(KPLIB,NREA,NGRP,NL,NPRC,NOMREA,NWT0,
     1    XS(1,1,NBISO+1),SIGS(1,1,NBISO+1),SS2D(1,1,1,NBISO+1),
     2    TAUXFI(NBISO+IRES),LXS,LAMB,CHIRS,BETAR,INVELS,IHRES,
     3    LSTRD,LPURE,ILUPS,ITRANC,IFISS)
          IF(MY1*MY2.GT.0) CALL ACRNDF(IMPX,NBISO+NBMAC,NBISO+IRES,
     1    IBM,HNOMIS,IPAPX,KPLIB,MY1,MY2,YLDSM,ISTY1(NBISO+IRES))
        ENDIF
        DEALLOCATE(NOMMAC)
      ELSE IF(IMAC.NE.0) THEN
        HRESID=NOMMAC(IMAC)
        CALL SCRFND(MAXISO,NBISO2I,NBISO2,IHRES,IBM,HRESID,HUSE2,HNAM2,
     1  IMIX2,JJSO(NBISO+IMAC))
        KPLIB=LCMDIL(JPLIB,JJSO(NBISO+IMAC)) ! step up isot JJSO(NBISO+IMAC)
        CALL ACRISO(KPLIB,NREA,NGRP,NL,NPRC,NOMREA,NWT0,XS(1,1,NBISO+1),
     1  SIGS(1,1,NBISO+1),SS2D(1,1,1,NBISO+1),TAUXFI(NBISO+IMAC),LXS,
     2  LAMB,CHIRS,BETAR,INVELS,IHRES,LSTRD,LPURE,ILUPS,ITRANC,IFISS)
        DEALLOCATE(NOMMAC)
      ENDIF
*----
*  SET NUMBER DENSITIES AND VOLUMES IN OUTPUT MICROLIB
*----
      IF(LRES) THEN
*       -- Number densities are left unchanged except if they are
*       -- listed in HISO array.
        DO 110 KSO=1,NISO(IBM) ! user-selected isotope
          DO JSO=1,NBISO2 ! microlib isotope
            IF(IMIX2(JSO).NE.IBM) CYCLE
            WRITE(TEXT8,'(2A4)') HUSE2(1,JSO),HUSE2(2,JSO)
            IF(HISO(IBM,KSO).EQ.TEXT8) THEN
              ITOD2(JSO)=ITODO(IBM,KSO)
              IF(CONC(IBM,KSO).EQ.-99.99) THEN
*               -- Only number densities of isotopes set with "MICR" and
*               -- "*" keywords are interpolated
                DENS2(JSO)=0.0
                DO ISO=1,NBISO1 ! Apex file isotope
                  IF(JJSO(ISO).EQ.JSO) DENS2(JSO)=DENS2(JSO)+DENS3(ISO)
                ENDDO
              ELSE IF(CONC(IBM,KSO).NE.-99.99) THEN
*               -- Number densities of isotopes set with "MICR" and
*               -- fixed value are forced to this value
                DENS2(JSO)=CONC(IBM,KSO)
              ENDIF
              GO TO 110
            ENDIF
          ENDDO
          WRITE(HSMG,'(31HACRLIB: UNABLE TO FIND ISOTOPE ,A8,6H IN MI,
     1    5HXTURE,I8,1H.)') HISO(IBM,KSO),IBM
          CALL XABORT(HSMG)
  110   CONTINUE
      ELSE
*       -- Number densities are interpolated or not according to
*       -- ALL/ONLY option
        DO JSO=1,NBISO2 ! microlib isotope
          WRITE(TEXT8,'(2A4)') HUSE2(1,JSO),HUSE2(2,JSO)
          IF(IBM.EQ.IMIX2(JSO)) THEN
            DO ISO=1,NBISO1 ! Apex file isotope
              IF(HNOMIS(ISO).EQ.TEXT8) THEN
                DENS2(JSO)=0.0
                VOL2(JSO)=0.0
                CYCLE
              ENDIF
            ENDDO
          ENDIF
        ENDDO
        DO 130 ISO=1,NBISO1 ! Apex file isotope
        IF(.NOT.LISO(IBM)) THEN
*         --ONLY option
          DO KSO=1,NISO(IBM) ! user-selected isotope
            IF(HNOMIS(ISO).EQ.HISO(IBM,KSO)) GO TO 120
          ENDDO
          GO TO 130
        ENDIF
  120   JSO=JJSO(ISO)
        IF(JSO.GT.0) THEN
          ITOD2(JSO)=ITOD1(ISO)
          ISTY2(JSO)=ISTY1(ISO)
          DENS2(JSO)=DENS2(JSO)+DENS3(ISO)
          VOL2(JSO)=VOL2(JSO)+XVOLM(IBMOLD)
        ENDIF
  130   CONTINUE
      ENDIF
*----
*  SET PIFI INFORMATION
*----
      ALLOCATE(IMICR(NBISO1))
      IMICR(:NBISO1)=0
      NBS1=0
      DO 140 JSO=1,NBISO2 ! microlib isotope
      IF(IMIX2(JSO).EQ.IBM) THEN
        NBS1=NBS1+1
        IF(NBS1.GT.NBISO1) CALL XABORT('ACRLIB: NBISO1 OVERFLOW.')
        IMICR(NBS1)=JSO
      ENDIF
  140 CONTINUE
      DO 170 ISO=1,NBS1 ! Apex file isotope
      JSO=IMICR(ISO)
      KPLIB=LCMDIL(JPLIB,JSO) ! step up isot JSO
      CALL LCMLEN(KPLIB,'PYIELD',LMY1,ITYLCM)
      IF(LMY1.GT.0) THEN
        ALLOCATE(HPYNAM(LMY1),IPYMIX(LMY1),IPIFI(LMY1))
        IPIFI(:LMY1)=0
        CALL LCMGTC(KPLIB,'PYNAM',8,LMY1,HPYNAM)
        CALL LCMGET(KPLIB,'PYMIX',IPYMIX)
        DO 160 IY1=1,LMY1
        IF(HPYNAM(IY1).NE.' ') THEN
          DO 150 KSO=1,NBS1
          LSO=IMICR(KSO)
          WRITE(TEXT8,'(2A4)') HUSE2(:2,LSO)
          IF((HPYNAM(IY1).EQ.TEXT8).AND.(IPYMIX(IY1).EQ.IMIX2(LSO)))
     1    THEN
            IPIFI(IY1)=LSO
            GO TO 160
          ENDIF
  150     CONTINUE
          IF(IPIFI(IY1).EQ.0) THEN
            WRITE(HSMG,'(40HACRLIB: FAILURE TO FIND FISSILE ISOTOPE ,
     1      A12,25H AMONG MICROLIB ISOTOPES.)') HPYNAM(IY1)
            CALL XABORT(HSMG)
          ENDIF
        ENDIF
  160   CONTINUE
        CALL LCMPUT(KPLIB,'PIFI',LMY1,1,IPIFI)
        DEALLOCATE(IPIFI,IPYMIX,HPYNAM)
      ENDIF
  170 CONTINUE
      DEALLOCATE(IMICR)
  180 CONTINUE ! end of loop over microlib mixtures.
*----
*  RELEASE MEMORY
*----
      DEALLOCATE(FACT,DENS1)
      IF(NPRC.GT.0) DEALLOCATE(INVELSB,BETARB,CHIRSB,INVELS,BETAR,
     1 CHIRS,LAMB)
      DEALLOCATE(LXS,XS,SS2D,SIGS,NWT0,TAUXFI)
      DEALLOCATE(ITOD1,YLDSM)
      IF(NBISO.GT.0) DEALLOCATE(NOMISO)
      DEALLOCATE(JJSO,DENS0,XVOLM,HNOMIS,IRESM,ITOTM)
*----
*  MICROLIB FINALIZATION
*----
      IF(.NOT.LRES) THEN
        ISTATE(:NSTATE)=0
        ISTATE(1)=NMIX
        ISTATE(2)=NBISO2
        ISTATE(3)=NGRP
        ISTATE(4)=NL
        ISTATE(5)=ITRANC
        ISTATE(7)=1
        IF(ITER.EQ.3) ISTATE(12)=NMIX
        ISTATE(13)=NED2
        ISTATE(14)=NMIX
        ISTATE(18)=1
        ISTATE(19)=NPRC
        ISTATE(20)=MY1
        ISTATE(22)=MAXISO/NMIX
        IF(NSURFD.GT.0) ISTATE(24)=3 ! ADF/CPDF information
        IF(NBISO2.EQ.0) CALL XABORT('ACRLIB: NBISO2=0.')
        CALL LCMPUT(IPLIB,'STATE-VECTOR',NSTATE,1,ISTATE)
        CALL LCMPUT(IPLIB,'MIXTURESVOL',NMIX,2,VOLMI2)
        CALL LCMPUT(IPLIB,'ISOTOPESUSED',3*NBISO2,3,HUSE2)
        CALL LCMPUT(IPLIB,'ISOTOPERNAME',3*NBISO2,3,HNAM2)
        CALL LCMPUT(IPLIB,'ISOTOPESDENS',NBISO2,2,DENS2)
        CALL LCMPUT(IPLIB,'ISOTOPESMIX',NBISO2,1,IMIX2)
        CALL LCMPUT(IPLIB,'ISOTOPESVOL',NBISO2,2,VOL2)
        IF(NED2.GT.0) CALL LCMPTC(IPLIB,'ADDXSNAME-P0',8,NED2,HVECT2)
        CALL LCMPUT(IPLIB,'ISOTOPESTODO',NBISO2,1,ITOD2)
        CALL LCMPUT(IPLIB,'ISOTOPESTYPE',NBISO2,1,ISTY2)
      ELSE IF(LRES.AND.(NBISO.GT.0)) THEN
        CALL LCMPUT(IPLIB,'ISOTOPESDENS',NBISO2,2,DENS2)
        CALL LCMPUT(IPLIB,'ISOTOPESVOL',NBISO2,2,VOL2)
      ENDIF
      IF(IMPX.GT.5) CALL LCMLIB(IPLIB)
      IACCS=1
*----
*  COMPUTE THE MACROSCOPIC X-SECTIONS
*----
      IF((ITER.NE.0).AND.(ITER.NE.3)) GO TO 280
      CALL LCMGET(IPLIB,'STATE-VECTOR',ISTATE)
      MAXMIX=ISTATE(1)
      IF(MAXMIX.NE.NMIX) CALL XABORT('ACRLIB: INVALID NMIX.')
      NBISO=ISTATE(2)
      ALLOCATE(MASK(MAXMIX),MASKL(NGRP))
      ALLOCATE(ISONA(3*NBISO),ISOMI(NBISO),DENIS(NBISO))
      CALL LCMGET(IPLIB,'ISOTOPESUSED',ISONA)
      CALL LCMGET(IPLIB,'ISOTOPESMIX',ISOMI)
      CALL LCMGET(IPLIB,'ISOTOPESDENS',DENIS)
      MASK(:MAXMIX)=.TRUE.
      MASKL(:NGRP)=.TRUE.
      ITSTMP=0
      TMPDAY(1)=0.0
      TMPDAY(2)=0.0
      TMPDAY(3)=0.0
      CALL LCMLEN(IPLIB,'MACROLIB',ILONG,ITYLCM)
      IF(ILONG.NE.0) CALL LCMDEL(IPLIB,'MACROLIB')
      CALL LIBMIX(IPLIB,MAXMIX,NGRP,NBISO,ISONA,ISOMI,DENIS,MASK,MASKL,
     1 ITSTMP,TMPDAY)
      DEALLOCATE(MASKL,MASK)
      DEALLOCATE(DENIS,ISOMI,ISONA)
      IF(NSURFD.GT.0) THEN
        CALL LCMSIX(IPLIB,'MACROLIB',1)
        CALL LCMGET(IPLIB,'STATE-VECTOR',ISTATE)
        ISTATE(12)=3 ! ADF information
        CALL LCMPUT(IPLIB,'STATE-VECTOR',NSTATE,1,ISTATE)
        CALL LCMSIX(IPLIB,' ',2)
      ENDIF
*----
*  INCLUDE LEAKAGE IN THE MACROLIB (USED ONLY FOR NON-REGRESSION TESTS)
*----
      IF(B2.NE.0.0) THEN
        IF(IMPX.GT.0) WRITE(IOUT,'(/31H ACRLIB: INCLUDE LEAKAGE IN THE,
     1  14H MACROLIB (B2=,1P,E12.5,2H).)') B2
        CALL LCMSIX(IPLIB,'MACROLIB',1)
        JPLIB=LCMGID(IPLIB,'GROUP')
        ALLOCATE(GAR1(NMIX),GAR2(NMIX))
        DO 270 IGR=1,NGRP
          KPLIB=LCMGIL(JPLIB,IGR)
          CALL LCMGET(KPLIB,'NTOT0',GAR1)
          CALL LCMGET(KPLIB,'DIFF',GAR2)
          DO 260 IBM=1,NMIX
            IF(MIXC(IBM).NE.0) GAR1(IBM)=GAR1(IBM)+B2*GAR2(IBM)
  260     CONTINUE
          CALL LCMPUT(KPLIB,'NTOT0',NMIX,2,GAR1)
  270   CONTINUE
        DEALLOCATE(GAR2,GAR1)
        CALL LCMSIX(IPLIB,' ',2)
      ENDIF
*----
*  PROCESS ADF INFORMATION
*----
  280 IF(NSURFD.GT.0) THEN
        DO 285 IBM=1,NMIX ! mixtures in Macrolib
          IF(MIXC(IBM).NE.0) ADF2(IBM,:NGRP,:NSURFD)=0.0
  285   CONTINUE
        DO 300 ICAL=1,NCAL
        DO 290 IBM=1,NMIX ! mixtures in Macrolib
        IBMOLD=MIXC(IBM)
        IF(IBMOLD.EQ.0) GO TO 290
        WEIGHT=TERP(ICAL,IBM)
        IF(WEIGHT.EQ.0.0) GO TO 290
        WRITE(RECNAM,'(4Hcalc,I8,15H/miscellaneous/)') ICAL
        K=0
        CALL hdf5_info(IPAPX,TRIM(RECNAM)//"ADF",RANK,TYPE,NBYTE,DIMSR)
        IF(TYPE.NE.99) THEN
          HHAD(K+1)='ADF'
          K=K+1
        ENDIF
        CALL hdf5_info(IPAPX,TRIM(RECNAM)//"CPDF",RANK,TYPE,NBYTE,DIMSR)
        IF(TYPE.NE.99) THEN
          HHAD(K+1)='CPDF'
          K=K+1
        ENDIF
        CALL hdf5_info(IPAPX,TRIM(RECNAM)//"INTERNAL_ADF",RANK,TYPE,
     1  NBYTE,DIMSR)
        IF(TYPE.NE.99) THEN
          HHAD(K+1)='INTERNAL_ADF'
          K=K+1
        ENDIF
        CALL hdf5_info(IPAPX,TRIM(RECNAM)//"INTERNAL_CPDF",RANK,TYPE,
     1  NBYTE,DIMSR)
        IF(TYPE.NE.99) THEN
          HHAD(K+1)='INTERNAL_CPDF'
          K=K+1
        ENDIF
        IF(4*K.NE.NSURFD) CALL XABORT('ACRLIB: INVALID ADF COUNT.')
        DO I=1,K
          CALL hdf5_get_shape(IPAPX,TRIM(RECNAM)//HHAD(I),DIMS_APX)
          ISURF=DIMS_APX(1)
          DEALLOCATE(DIMS_APX)
          CALL hdf5_read_data(IPAPX,TRIM(RECNAM)//HHAD(I),ADF)
          DO I0=1,ISURF
            IF(HHAD(I).EQ.'ADF') THEN
              WRITE(TEXT8,'(3HADF,I1)') I0
            ELSE IF(HHAD(I).EQ.'CPDF') THEN
              WRITE(TEXT8,'(4HCPDF,I1)') I0
            ELSE IF(HHAD(I).EQ.'INTERNAL_ADF') THEN
              WRITE(TEXT8,'(6HIN_ADF,I1)') I0
            ELSE IF(HHAD(I).EQ.'INTERNAL_CPDF') THEN
              WRITE(TEXT8,'(7HIN_CPDF,I1)') I0
            ENDIF
            IOF=(I-1)*ISURF+I0
            HADF(IOF)=TEXT8
            DO IGR=1,NGRP
              ADF2(IBM,IGR,IOF)=ADF2(IBM,IGR,IOF)+WEIGHT*ADF(I0,IGR)
            ENDDO
          ENDDO
          DEALLOCATE(ADF)
        ENDDO
  290   CONTINUE
  300   CONTINUE
        CALL LCMSIX(IPLIB,'MACROLIB',1)
        CALL LCMSIX(IPLIB,'ADF',1)
        CALL LCMPUT(IPLIB,'NTYPE',1,1,NSURFD)
        CALL LCMPTC(IPLIB,'HADF',8,NSURFD,HADF)
        DO I=1,NSURFD
          CALL LCMPUT(IPLIB,HADF(I),NMIX*NGRP,2,ADF2(1,1,I))
        ENDDO
        CALL LCMSIX(IPLIB,' ',2)
        CALL LCMSIX(IPLIB,' ',2)
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(ADF2,HADF)
      DEALLOCATE(SPH,FLUX,VOLMI2,VOL2,DENS3,DENS2)
      DEALLOCATE(HNAM2,HUSE2,ISTY2,ISTY1,ITOD2,IMIX2)
      RETURN
      END
