*DECK LIBCTL
      SUBROUTINE LIBCTL (MAXMIX,MAXISO,IPLIB,IPLIBX,INDREC,IMAC,ISOADD,
     1 NDEPL,IMPX,NBISO,NBMIX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Catenate two microscopic cross section libraries.
*
*Copyright:
* Copyright (C) 2024 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): A. Hebert
*
*Parameters: input
* MAXMIX  maximum value of NBMIX.
* MAXISO  maximum number of isotopes permitted.
* IPLIB   pointer to the lattice microscopic cross section library
*         (L_LIBRARY signature).
* IPLIBX  pointer to a read-only microscopic cross section library
*         (L_LIBRARY signature) to catenate.
* INDREC  type of action:
*         =1 a new microlib is created; =2 the microlib is updated.
* IMAC    macrolib construction flag:
*         =0 do not compute an embedded macrolib;
*         =1 compute an embedded macrolib.
* ISOADD  flag to complete the depletion chain:
*         =0 complete; =1 do not complete.
* NDEPL   number of depleting isotopes (used by EVO:).
* IMPX    print flag.
*
*Parameters: output
* NBISO   number of isotopes present in the calculation domain.
* NBMIX   number of mixtures defined in the library.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB,IPLIBX
      INTEGER MAXMIX,MAXED,MAXISO,INDREC,IMAC,ISOADD,NDEPL,IMPX,NBISO
*----
*  LOCAL PARAMETERS
*----
      PARAMETER (IOUT=6,NSTATE=40,MAXED=50,MAXLIB=20)
      INTEGER ISTATE(NSTATE)
      REAL TMPDAY(3)
      DOUBLE PRECISION DBLLIR
      LOGICAL EMPTY,LCM
      CHARACTER HSMG*131,CARLIR*12,HSIGN*12,TEXT12*12
      TYPE(C_PTR) JPLIB,KPLIB,KPLIB2
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ISOMIX,NTFG,LSHI,NIR,
     1 IEVOL,ITYP,ILLIB,KGAS,ISOMIX2,NTFG2,LSHI2,NIR2,IEVOL2,ITYP2,
     2 ILLIB2,KGAS2,LOCUPD
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ISONAM,ISONRF,ISONAM2,
     > ISONRF2
      REAL, ALLOCATABLE, DIMENSION(:) :: DENISO,DENMIX,TMPISO,SNISO,
     > SBISO,GIR,DENISO2,TMPISO2,SNISO2,SBISO2,GIR2,DENMIX2,ENERGY
      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:,:) :: HLIB,HLIB2
      CHARACTER(LEN=12), ALLOCATABLE, DIMENSION(:) :: SHINA,SHINA2
      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:) :: HVECT,HVECT2
      CHARACTER(LEN=64), ALLOCATABLE, DIMENSION(:) :: HNAME,HNAME2
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: MASK,MASKI,MASKL,MASKI2
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: JPISO2
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(ISONAM(3,MAXISO),ISONRF(3,MAXISO),ISOMIX(MAXISO),
     > NTFG(MAXISO),LSHI(MAXISO),NIR(MAXISO),IEVOL(MAXISO),
     > ITYP(MAXISO),ILLIB(MAXISO),KGAS(MAXMIX),LOCUPD(MAXMIX))
      ALLOCATE(SHINA(MAXISO),HLIB(MAXISO,4),HVECT(MAXED),
     > HNAME(MAXLIB))
      ALLOCATE(DENISO(MAXISO),TMPISO(MAXISO),SNISO(MAXISO),
     > SBISO(MAXISO),GIR(MAXISO),DENMIX(MAXMIX))
       ALLOCATE(MASKI(MAXISO),MASK(MAXMIX))
*----
*  RECOVER FIRST MICROLIB
*----
      IF((INDREC.EQ.2).AND.(NBISO.GT.0)) THEN
*       THE LIBRARY IS UPDATED. READ OLD LIBRARY INFORMATION.
        CALL LIBINF(IPLIB,MAXISO,MAXLIB,MAXED,MAXMIX,NBISO,NGRO,NL,
     1  ITRANC,NLIB,NCOMB,NEDMAC,NBMIX,ISONAM,ISONRF,ISOMIX,DENISO,
     2  TMPISO,SHINA,SNISO,SBISO,NTFG,LSHI,GIR,NIR,MASKI,HLIB,IEVOL,
     3  ITYP,ILLIB,KGAS,DENMIX,HVECT,HNAME)
      ELSE
        NBISO=0
        NGRO=0
        NL=0
        ITRANC=-999
        NLIB=0
        NEDMAC=0
        NCOMB=0
        NEDMAC=0
        NBMIX=0
        DENMIX(:MAXMIX)=-1.0
      ENDIF
*----
*  RECOVER A READ-ONLY MICROLIB
*----
      IF(C_ASSOCIATED(IPLIBX)) THEN
        CALL LCMINF(IPLIBX,CARLIR,TEXT12,EMPTY,ILCMLN,LCM)
        CALL LCMGTC(IPLIBX,'SIGNATURE',12,HSIGN)
        IF(HSIGN.NE.'L_LIBRARY') THEN
          CALL XABORT('LIBCTL: SIGNATURE OF '//CARLIR//' IS '//HSIGN//
     1    '. L_LIBRARY EXPECTED.')
        ENDIF
        CALL LCMGET(IPLIBX,'STATE-VECTOR',ISTATE)
        NBMIX2=ISTATE(1)
        NBISO2=ISTATE(2)
        NLIB2=ISTATE(8)
        NEDMAC2=ISTATE(13)
        IF((NGRO.GT.0).AND.(ISTATE(3).NE.NGRO)) THEN
          WRITE(HSMG,'(44HLIBCTL: THE RHS MICROLIB HAS AN INCOMPATIBLE,
     1    26H NUMBER OF ENERGY GROUPS (,I5,4H VS ,I5,2H).)') ISTATE(3),
     2    NGRO
          CALL XABORT(HSMG)
        ENDIF
        NGRO=ISTATE(3)
        IF((ITRANC.GT.0).AND.(ISTATE(5).NE.ITRANC)) THEN
          WRITE(HSMG,'(44HLIBCTL: THE RHS MICROLIB HAS AN INCOMPATIBLE,
     1    23H TRANSPORT CORRECTION (,I5,4H VS ,I5,2H).)') ISTATE(5),
     2    ITRANC
          CALL XABORT(HSMG)
        ENDIF
        ITRANC=ISTATE(5)
        ALLOCATE(ISONAM2(3,NBISO2),ISONRF2(3,NBISO2),ISOMIX2(NBISO2),
     >  NTFG2(NBISO2),LSHI2(NBISO2),NIR2(NBISO2),IEVOL2(NBISO2),
     >  ITYP2(NBISO2),ILLIB2(NBISO2),KGAS2(NBMIX2))
        ALLOCATE(SHINA2(NBISO2),HLIB2(NBISO2,4),HVECT2(MAXED),
     >  HNAME2(MAXLIB))
        ALLOCATE(DENISO2(NBISO2),TMPISO2(NBISO2),SNISO2(NBISO2),
     >  SBISO2(NBISO2),GIR2(NBISO2),DENMIX2(MAXMIX),ENERGY(NGRO+1))
        ALLOCATE(MASKI2(NBISO2))
        ALLOCATE(JPISO2(NBISO2))
        CALL LCMGET(IPLIBX,'ENERGY',ENERGY)
        CALL LCMPUT(IPLIB,'ENERGY',NGRO+1,2,ENERGY)
        CALL LIBINF(IPLIBX,NBISO2,MAXLIB,MAXED,MAXMIX,NBISO2,NGRO,NL2,
     1  ITRANC,NLIB2,NCOMB2,NEDMAC2,NBMIX2,ISONAM2,ISONRF2,ISOMIX2,
     2  DENISO2,TMPISO2,SHINA2,SNISO2,SBISO2,NTFG2,LSHI2,GIR2,NIR2,
     3  MASKI2,HLIB2,IEVOL2,ITYP2,ILLIB2,KGAS2,DENMIX2,HVECT2,HNAME2)
        JPLIB=LCMLID(IPLIB,'ISOTOPESLIST',MAXISO)
        CALL LIBIPS(IPLIBX,NBISO2,JPISO2)
*----
*  ADD ENTRIES INTO HVECT2 AND HNAME
*----
        DO IED=1,NEDMAC2
          DO JED=1,NEDMAC
            IF(HVECT2(IED).EQ.HVECT(JED)) GO TO 40
          ENDDO
          NEDMAC=NEDMAC+1
          IF(NEDMAC.GT.MAXED) CALL XABORT('LIBCTL: MAXED OVERFLOW(2).')
          HVECT(NEDMAC)=HVECT2(IED)
   40     CONTINUE
        ENDDO
        DO ILIB=1,NLIB2
          DO JLIB=1,NLIB
            IF(HNAME2(ILIB).EQ.HNAME(JLIB)) GO TO 45
          ENDDO
          NLIB=NLIB+1
          IF(NLIB.GT.MAXLIB) CALL XABORT('LIBCTL: MAXLIB OVERFLOW(2).')
          HNAME(NLIB)=HNAME2(ILIB)
   45     CONTINUE
        ENDDO
*----
*  READ GROUP (descmix2)
*----
        LOCUPD(:MAXMIX)=0
        NNMIX=0
   60   CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
   70   IF(ITYPLU.NE.3) CALL XABORT('LIBCTL: CHARACTER DATA EXPECTED.')
        IF(CARLIR.EQ.';') THEN
          GO TO 80
        ELSE IF(CARLIR(1:3).EQ.'MIX') THEN
          CALL REDGET(ITYPLU,NNMIX,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU.NE. 1) CALL XABORT('LIBCTL: MIXTURE TO UPDATE MUST'
     >    //' BE GIVEN.')
          IF(NNMIX.GT.MAXMIX) THEN
            CALL XABORT('LIBCTL: MAXMIX OVERFLOW.')
          ELSE IF(NNMIX.LE.0) THEN
            CALL XABORT('LIBCTL: MIX NUMBER.LE.0.')
          ENDIF
          NBMIX=MAX(NBMIX,NNMIX)
          LOCUPD(NNMIX)=NNMIX
          CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          IF(ITYPLU.EQ.1) THEN
            IF(INTLIR.LE.0 .OR.INTLIR.GT.NBMIX2) THEN
              CALL XABORT('LIBCTL: RHS NBMIX OVERFLOW.')
            ENDIF
            LOCUPD(NNMIX)=INTLIR
            CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          ENDIF
          DENMIX(NNMIX)=DENMIX2(LOCUPD(NNMIX))
          IF(ITYPLU.EQ.2) THEN
            IF(DENMIX(NNMIX).EQ.-1.0) THEN
              CALL LIBCON(IPLIBX,NNMIX,NBISO2,ISOMIX2,DENISO2,
     >        DENMIX(NNMIX),2)
            ENDIF
            DENMIX(NNMIX)=REALIR
            CALL REDGET(ITYPLU,INTLIR,REALIR,CARLIR,DBLLIR)
          ENDIF
          GO TO 70
        ELSE IF(CARLIR(1:3).EQ.'ALL') THEN
          DO ISO=1,NBISO2
            NNMIX=ISOMIX2(ISO)
            NBMIX=MAX(NBMIX,NNMIX)
            LOCUPD(NNMIX)=NNMIX
          ENDDO
        ELSE
          WRITE(HSMG,'(16HLIBCTL: KEYWORD ,A,20H IS NOT IMPLEMENTED.)')
     >    TRIM(CARLIR)
          CALL XABORT(HSMG)
        ENDIF
        GO TO 60
*----
*  CATENATE LIBRARIES
*----
   80   NL=MAX(NL,NL2)
        MIXTURE: DO NNMIX=1,NBMIX
          IF(LOCUPD(NNMIX).EQ.0) CYCLE
          DO ISO=1,NBISO
            IF(ISOMIX(ISO).EQ.NNMIX) THEN
              WRITE(HSMG,'(15HLIBCTL: MIXTURE,I8,18H IS ALREADY DEFINE,
     >        14HD FOR ISOTOPE ,3A4,1H.)') NNMIX,(ISONAM(I,ISO),I=1,3)
              WRITE(6,'(1X,A)') TRIM(HSMG)
              CYCLE MIXTURE
            ENDIF
          ENDDO
          KGAS(NNMIX)=KGAS2(LOCUPD(NNMIX))
          OUTER: DO ISO2=1,NBISO2
            IF(ISOMIX2(ISO2).NE.LOCUPD(NNMIX)) CYCLE
            NBISO=NBISO+1
            IF(NBISO.GT.MAXISO) THEN
              WRITE(6,'(33H LIBCTL: MAXISO OVERFLOW (MAXISO=,I7,4H ISO,
     >        18HTOPES PER MIXTURE=,I7,2H).)') MAXISO,1+MAXISO/NBMIX
              CALL XABORT(HSMG)
            ENDIF
            ISOMIX(NBISO)=NNMIX
            DENISO(NBISO)=DENISO2(ISO2)
            TMPISO(NBISO)=TMPISO2(ISO2)
            NIR(NBISO)=NIR2(ISO2)
            GIR(NBISO)=GIR2(ISO2)
            IEVOL(NBISO)=IEVOL2(ISO2)
            ITYP(NBISO)=ITYP2(ISO2)
            SNISO(NBISO)=SNISO2(ISO2)
            SBISO(NBISO)=SBISO2(ISO2)
            ISONAM(:3,NBISO)=ISONAM2(:3,ISO2)
            ISONRF(:3,NBISO)=ISONRF2(:3,ISO2)
            HLIB(NBISO,:4)=HLIB2(ISO2,:4)
            ILLIB(NBISO)=ILLIB2(ISO2)
            LSHI(NBISO)=LSHI2(ISO2)
            SHINA(NBISO)=SHINA2(ISO2)
            NTFG(NBISO)=NTFG2(ISO2)
            MASKI(NBISO)=.FALSE.
            DO JSO=1,NBISO-1
              IF((ISONAM(1,JSO).EQ.ISONAM(1,NBISO)).AND.
     1           (ISONAM(2,JSO).EQ.ISONAM(2,NBISO)).AND.
     2           (ISONAM(3,JSO).EQ.ISONAM(3,NBISO))) CYCLE OUTER
            ENDDO
            MASKI(NBISO)=.TRUE.
            KPLIB2=JPISO2(ISO2) ! set ISO2-th isotope
            IF(.NOT.C_ASSOCIATED(KPLIB2)) THEN
              WRITE(HSMG,'(17HLIBCTL: ISOTOPE '',3A4,7H'' (ISO=,I8,
     1        39H) IS NOT AVAILABLE IN THE RHS MICROLIB.)')
     2        (ISONAM2(I0,ISO2),I0=1,3),ISO2
              CALL XABORT(HSMG)
            ENDIF
            JPLIB=LCMGID(IPLIB,'ISOTOPESLIST')
            KPLIB=LCMDIL(JPLIB,NBISO) ! create a new list entry
            CALL LCMEQU(KPLIB2,KPLIB) ! KPLIB2 --> KPLIB
          ENDDO OUTER
          IF(DENMIX(NNMIX).GE.0.0) THEN
            CALL LIBCON(IPLIB,NNMIX,NBISO,ISOMIX,DENISO,DENMIX(NNMIX),1)
          ENDIF
        ENDDO MIXTURE
        DEALLOCATE(JPISO2)
        DEALLOCATE(MASKI2)
        DEALLOCATE(ENERGY,DENMIX2,GIR2,SBISO2,SNISO2,TMPISO2,DENISO2)
        DEALLOCATE(HNAME2,HVECT2,HLIB2,SHINA2)
        DEALLOCATE(KGAS2,ILLIB2,ITYP2,IEVOL2,NIR2,LSHI2,NTFG2,ISOMIX2,
     1  ISONRF2,ISONAM2)
      ENDIF
*----
*  ADD THE MISSING ISOTOPES FROM THE DEPLETION CHAIN.
*----
      IF((NDEPL.NE.0).AND.(ISOADD.EQ.0)) THEN
        NBISOL=NBISO
        CALL LCMSIX(IPLIB,'DEPL-CHAIN',1)
        CALL LCMGET(IPLIB,'STATE-VECTOR',ISTATE)
        IF(ISTATE(1).NE.NDEPL) CALL XABORT('LIBCTL: INVALID NUMBER OF '
     >  //'DEPLETING ISOTOPES.')
        NFISS=ISTATE(2)
        NSUPF=ISTATE(5)
        NSUPS=ISTATE(7)
        NREAC=ISTATE(8)
        NPAR=ISTATE(9)
        CALL LIBEAD(IPLIB,MAXISO,MAXMIX,IMPX,NDEPL,NFISS,NSUPS,
     1  NREAC,NPAR,NBISO,ISONAM,ISONRF,HLIB,ILLIB,ISOMIX,TMPISO,
     2  IEVOL,ITYP,NCOMB)
        CALL LCMSIX(IPLIB,' ',2)
*
        DO ISOT=NBISOL+1,NBISO
          SNISO(ISOT)=1.0E10
          SBISO(ISOT)=1.0E10
          DENISO(ISOT)=0.0
          NTFG(ISOT)=0
          HLIB(ISOT,2)=SHINA(ISOT)(:8)
          HLIB(ISOT,3)=SHINA(ISOT)(:8)
          HLIB(ISOT,4)=SHINA(ISOT)(:8)
          LSHI(ISOT)=0
          GIR(ISOT)=1.0
          NIR(ISOT)=0
          MASKI(ISOT)=.TRUE.
        ENDDO
      ENDIF
*----
*  SET THE MIXTURE MASKS.
*----
      DO 110 I=1,NBMIX
        MASK(I)=.FALSE.
        DO 90 JJ=1,NBISO
          IF(ISOMIX(JJ).EQ.I) THEN
            MASK(I)=.TRUE.
            GO TO 100
          ENDIF
   90   CONTINUE
  100   CONTINUE
  110 CONTINUE
*----
*  SAVE THE LIBRARY SPECIFIC INFORMATION.
*----
      NBMIX=0
      DO I=1,NBISO
        NBMIX=MAX(NBMIX,ISOMIX(I))
      ENDDO
      IF(NBMIX.GT.MAXMIX) CALL XABORT('LIBCTL: MAXMIX TOO SMALL.')
      TEXT12='L_LIBRARY'
      CALL LCMPTC(IPLIB,'SIGNATURE',12,TEXT12)
      ISTATE(:NSTATE)=0
      ISTATE(1)=MAXMIX
      ISTATE(2)=NBISO
      ISTATE(3)=NGRO
      ISTATE(4)=NL
      ISTATE(5)=ITRANC
      ISTATE(8)=NLIB
      ISTATE(11)=NDEPL
      ISTATE(12)=NCOMB
      ISTATE(13)=NEDMAC
      ISTATE(14)=NBMIX
      ISTATE(17)=-1
      ISTATE(18)=IMAC
      ISTATE(21)=ISOADD
      CALL LCMPUT(IPLIB,'STATE-VECTOR',NSTATE,1,ISTATE)
      CALL LCMPUT(IPLIB,'ISOTOPESUSED',3*NBISO,3,ISONAM)
      CALL LCMPUT(IPLIB,'ISOTOPERNAME',3*NBISO,3,ISONRF)
      CALL LCMPUT(IPLIB,'ISOTOPESMIX',NBISO,1,ISOMIX)
      CALL LCMPUT(IPLIB,'ISOTOPESTODO',NBISO,1,IEVOL)
      CALL LCMPUT(IPLIB,'ISOTOPESTYPE',NBISO,1,ITYP)
      CALL LCMPTC(IPLIB,'ILIBRARYTYPE',8,NBISO,HLIB(:NBISO,1))
      CALL LCMPUT(IPLIB,'ILIBRARYINDX',NBISO,1,ILLIB)
      CALL LCMPTC(IPLIB,'ISOTOPESCOH',8,NBISO,HLIB(:NBISO,2))
      CALL LCMPTC(IPLIB,'ISOTOPESINC',8,NBISO,HLIB(:NBISO,3))
      CALL LCMPTC(IPLIB,'ISOTOPESRESK',8,NBISO,HLIB(:NBISO,4))
      CALL LCMPUT(IPLIB,'ISOTOPESNTFG',NBISO,1,NTFG)
      CALL LCMPTC(IPLIB,'ISOTOPESHIN',12,NBISO,SHINA)
      CALL LCMPUT(IPLIB,'ISOTOPESSHI',NBISO,1,LSHI)
      CALL LCMPUT(IPLIB,'ISOTOPESGIR',NBISO,2,GIR)
      CALL LCMPUT(IPLIB,'ISOTOPESNIR',NBISO,1,NIR)
      CALL LCMPUT(IPLIB,'ISOTOPESTEMP',NBISO,2,TMPISO)
      CALL LCMPUT(IPLIB,'ISOTOPESDENS',NBISO,2,DENISO)
      IF(NEDMAC.GT.0) THEN
        CALL LCMPTC(IPLIB,'ADDXSNAME-P0',8,NEDMAC,HVECT)
      ENDIF
      IF(NLIB.GT.0) THEN
        CALL LCMPTC(IPLIB,'ILIBRARYNAME',64,NLIB,HNAME)
      ENDIF
      CALL LCMPUT(IPLIB,'MIXTUREGAS',NBMIX,1,KGAS)
*----
*  CHECK FOR DUPLICATE ALIAS.
*----
      DO 130 I=1,NBISO
      IF(ISOMIX(I).EQ.0) GO TO 130
      DO 120 J=I+1,NBISO
      IF((ISOMIX(I).EQ.ISOMIX(J)).AND.(ISONRF(1,I).EQ.ISONRF(1,J))
     1   .AND.(ISONRF(2,I).EQ.ISONRF(2,J))
     2   .AND.(ISONRF(3,I).EQ.ISONRF(3,J)).AND.(LSHI(I).NE.0)) THEN
         WRITE(HSMG,200) (ISONAM(I1,I),I1=1,3),(ISONAM(I1,J),I1=1,3),
     >   (ISONRF(I1,I),I1=1,3),ISOMIX(I)
         CALL XABORT(HSMG)
      ENDIF
  120 CONTINUE
  130 CONTINUE
*----
*  PRINT TABLE OF CONTENTS.
*----
      IF(IMPX.GT.1) THEN
         IF(NLIB.GT.0) THEN
            WRITE(IOUT,150)
            DO ILIB=1,NLIB
              WRITE(IOUT,'(1X,I4,4H -- ,A)') ILIB,HNAME(ILIB)
            ENDDO
         ENDIF
         WRITE(IOUT,160)
         DO I=1,NBISO
           IF(ISOMIX(I).EQ.0) CYCLE
           IF(MASK(ISOMIX(I))) THEN
              WRITE(IOUT,170) I,(ISONAM(I0,I),I0=1,3),(ISONRF(I0,I),
     >        I0=1,3),HLIB(I,1),ILLIB(I),ISOMIX(I),DENISO(I),
     >        TMPISO(I),SNISO(I),LSHI(I),SHINA(I),NTFG(I),HLIB(I,3),
     >        HLIB(I,4),HLIB(I,2)
           ENDIF
        ENDDO
      ENDIF
*----
*  COMPUTE THE MACROSCOPIC X-SECTIONS.
*----
      IF(IMAC.EQ.1) THEN
         ALLOCATE(MASKL(NGRO))
         MASKL(:NGRO)=.TRUE.
         ITSTMP=0
         TMPDAY(1)=0.0
         TMPDAY(2)=0.0
         TMPDAY(3)=0.0
         CALL LIBMIX(IPLIB,NBMIX,NGRO,NBISO,ISONAM,ISOMIX,DENISO,MASK,
     >   MASKL,ITSTMP,TMPDAY)
         DEALLOCATE(MASKL)
      ENDIF
*----
*  PRINT STATE VECTOR
*----
      IF(IMPX.GT.0) WRITE(IOUT,180) IMPX,NLIB,NBISO,NBMIX,NDEPL,NCOMB,
     > NEDMAC,NGRO,NL,ITRANC,IMAC,ISOADD
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(MASK,MASKI)
      DEALLOCATE(DENMIX,GIR,SBISO,SNISO,TMPISO,DENISO)
      DEALLOCATE(HNAME,HVECT,HLIB,SHINA)
      DEALLOCATE(LOCUPD,KGAS,ILLIB,ITYP,IEVOL,NIR,LSHI,NTFG,ISOMIX,
     1 ISONRF,ISONAM)
      RETURN
*
  150 FORMAT(/35H AVAILABLE CROSS-SECTION LIBRARIES:)
  160 FORMAT(58X,'NUMBER'/' SPEC     LOCAL NAME    ISOTOPE       FRO',
     1 'M LIBRARY  MIX   DENSITY     TEMP(K)    SIGZERO    SELF-SHIEL',
     2 '  THERMAL CORRECTION'/' -------  ------------  ------------  ',
     3 '------------  ----  ----------  ---------  ---------  -------',
     4 '---  ------------------')
  170 FORMAT(1X,I7,2X,3A4,2X,3A4,2X,A8,I4,2X,I4,1P,E12.4,2E11.3,I4,2X,
     1 A8,I4,1X,3A8)
  180 FORMAT(/16H LIBCTL: OPTIONS/8H -------/
     1 7H IMPX  ,I8,30H   (0=NO PRINT/1=SHORT/2=MORE)/
     2 7H NLIB  ,I8,32H   (NUMBER OF SETS OF LIBRARIES)/
     3 7H NBISO ,I8,36H   (NUMBER OF ISOTOPES OR MATERIALS)/
     4 7H NBMIX ,I8,23H   (NUMBER OF MIXTURES)/
     5 7H NDEPL ,I8,33H   (NUMBER OF DEPLETING ISOTOPES)/
     6 7H NCOMB ,I8,33H   (NUMBER OF DEPLETING MIXTURES)/
     7 7H NEDMAC,I8,34H   (NUMBER OF CROSS SECTION EDITS)/
     8 7H NGRO  ,I8,28H   (NUMBER OF ENERGY GROUPS)/
     9 7H NL    ,I8,30H   (NUMBER OF LEGENDRE ORDERS)/
     1 7H ITRANC,I8,45H   (0=NO TRANSPORT CORRECTION/1=APOLLO TYPE/2,
     2 57H=RECOVER FROM LIBRARY/3=WIMS-D TYPE/4=LEAKAGE CORRECTION)/
     3 7H IMAC  ,I8,45H   (0=DO NOT/1=DO BUILD AN EMBEDDED MACROLIB)/
     4 7H ISOADD,I8,37H   (0=COMPLETE BURNUP CHAIN/1=DO NOT))
  200 FORMAT(9HLIBCTL: ',3A4,7H' AND ',3A4,24H' ARE BOTH ALIAS FOR THE,
     1 23H SAME LIBRARY ISOTOPE ',3A4,12H' IN MIXTURE,I5,1H.)
      END
