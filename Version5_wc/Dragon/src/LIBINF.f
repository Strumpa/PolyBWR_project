*DECK LIBINF
      SUBROUTINE LIBINF (IPLIB,MAXISO,MAXLIB,MAXED,MAXMIX,NBISO,NGRO,
     1 NL,ITRANC,NLIB,NCOMB,NEDMAC,NBMIX,ISONAM,ISONRF,ISOMIX,DENISO,
     2 TMPISO,SHINA,SNISO,SBISO,NTFG,LSHI,GIR,NIR,MASKI,HLIB,IEVOL,
     3 ITYP,ILLIB,KGAS,DENMIX,HVECT,HNAME)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover general information from a microlib.
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
* IPLIB   pointer to the lattice microscopic cross section library
*         (L_LIBRARY signature).
* MAXISO  maximum number of isotopes permitted.
* MAXLIB  maximum number of external cross-section libraries.
* MAXED   maximum number of extra vector edits.
* MAXMIX  maximum number of material mixtures.
*
*Parameters: output
* NBISO   number of isotopes present in the microlib.
* NGRO    number of energy groups.
* NL      anisotropy order in the microlib.
* ITRANC  type of transport correction: =0 no transport correction
*         =1 Apollo type transport correction; =2 recover from
*         library; =3 WIMS-D type; =4 leakage correction alone.
* NLIB    number of cross-section libraries.
* NCOMB   number of depleting mixtures (used by EVO:).
* NEDMAC  number of extra vector edits.
* NBMIX   number of mixtures defined in the microlib.
* ISONAM  alias name of each isotope.
* ISONRF  library name of each isotope.
* ISOMIX  mix number of each isotope.
* DENISO  density of each isotope.
* MASK    mixture masks.* TMPISO  temperature of each isotope.
* SHINA   self-shielding name of each isotope.
* SNISO   dilution cross section of each isotope. A value of 1.0E10
*         is used for infinite dilution.
* SBISO   dilution cross section of each isotope used with Livolant-
*         Jeanpierre normalization.
* NTFG    number of thermal groups where the thermal inelastic
*         correction is applied.
* LSHI    resonant region number associated with i-th isotope.
*         Infinite dilution will be assumed if LSHI(I)=0. A negative
*         value is indicating correlation of cross sections with all
*         isotopes sharing the same LSHI value.
* GIR     Goldstein-Cohen IR parameter of each isotope.
* NIR     Goldstein-Cohen IR cutoff energy index. Use IR approximation.
*         for groups with index.ge.nir; Use library value if NIR=0.
* MASKI   isotope masks.
* HLIB    isotope options.
* IEVOL   flag making an isotope non-depleting:
*         =1 to force an isotope to be non-depleting;
*         =2 to force an isotope to be depleting;
*         =3 to force an isotope to be at saturation.
* ITYP    isotopic type:
*         =1: the isotope is not fissile and not a fission product;
*         =2: the isotope is fissile; =3: is a fission product.
* ILLIB   xs library index for each isotope (.le.NLIB).
* KGAS    state of mixture (used for stopping power correction):
*         =0: solid or liquid;
*         =1: gas.
* DENMIX  mixture density (set to -1.0 to avoid using them).
* HVECT   extra vector edits names.
* HNAME   external cross-section libraries names.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPLIB
      INTEGER MAXISO,MAXLIB,MAXED,MAXMIX,NBISO,NGRO,NL,ITRANC,NLIB,
     1 NCOMB,NEDMAC,NBMIX,ISONAM(3,MAXISO),ISONRF(3,MAXISO),
     2 ISOMIX(MAXISO),NTFG(MAXISO),LSHI(MAXISO),NIR(MAXISO),
     3 IEVOL(MAXISO),ITYP(MAXISO),ILLIB(MAXISO),KGAS(MAXMIX)
      REAL DENISO(MAXISO),TMPISO(MAXISO),SNISO(MAXISO),SBISO(MAXISO),
     2 GIR(MAXISO),DENMIX(MAXMIX)
      LOGICAL MASKI(MAXISO)
      CHARACTER(LEN=12) SHINA(MAXISO)
      CHARACTER(LEN=8) HLIB(MAXISO,4),HVECT(MAXED)
      CHARACTER(LEN=64) HNAME(MAXLIB)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      TYPE(C_PTR) JPLIB
      INTEGER ISTATE(NSTATE)
*----
*  RECOVER STATE-VECTOR INFORMATION
*----
      CALL LCMGET(IPLIB,'STATE-VECTOR',ISTATE)
      NBISO=ISTATE(2)
      NGRO=ISTATE(3)
      NL=ISTATE(4)
      ITRANC=ISTATE(5)
      NLIB=ISTATE(8)
      NCOMB=ISTATE(12)
      NEDMAC=ISTATE(13)
      NBMIX=ISTATE(14)
      IF(NBISO.GT.MAXISO) CALL XABORT('LIBINF: MAXISO OVERFLOW.')
      IF(NLIB.GT.MAXLIB) CALL XABORT('LIBINF: MAXLIB OVERFLOW(1).')
      IF(NEDMAC.GT.MAXED) CALL XABORT('LIBINF: MAXED OVERFLOW(1).')
      IF(NBMIX.GT.MAXMIX) CALL XABORT('LIBINF: MAXMIX OVERFLOW.')
*----
*  RECOVER ISOTOPIC INFORMATION
*----
      CALL LCMGET(IPLIB,'ISOTOPESUSED',ISONAM)
      CALL LCMLEN(IPLIB,'ISOTOPERNAME',ILENG,ITYLCM)
      IF(ILENG.GT.0) THEN
        CALL LCMGET(IPLIB,'ISOTOPERNAME',ISONRF)
      ELSE
        CALL LCMGET(IPLIB,'ISOTOPESUSED',ISONRF)
      ENDIF
      HLIB(NBISO,:4)=' '
      ILLIB(:NBISO)=0
      CALL LCMLEN(IPLIB,'ILIBRARYTYPE',ILENG,ITYLCM)
      IF(ILENG.GT.0) THEN
        CALL LCMGTC(IPLIB,'ILIBRARYTYPE',8,NBISO,HLIB(1,1))
        CALL LCMGET(IPLIB,'ILIBRARYINDX',ILLIB)
      ENDIF
      CALL LCMLEN(IPLIB,'ISOTOPESNTFG',ILENG,ITYLCM)
      IF(ILENG.GT.0) THEN
        CALL LCMGET(IPLIB,'ISOTOPESNTFG',NTFG)
        CALL LCMGTC(IPLIB,'ISOTOPESCOH',8,NBISO,HLIB(1,2))
        CALL LCMGTC(IPLIB,'ISOTOPESINC',8,NBISO,HLIB(1,3))
      ELSE
        NTFG(:NBISO)=0
        HLIB(:NBISO,2)=' '
        HLIB(:NBISO,3)=' '
      ENDIF
      CALL LCMLEN(IPLIB,'ISOTOPESRESK',ILENG,ITYLCM)
      IF(ILENG.GT.0) THEN
        CALL LCMGTC(IPLIB,'ISOTOPESRESK',8,NBISO,HLIB(1,4))
      ELSE
        HLIB(:NBISO,4)=' '
      ENDIF
      CALL LCMLEN(IPLIB,'ISOTOPESHIN',ILENG,ITYLCM)
      IF(ILENG.GT.0) THEN
        CALL LCMGTC(IPLIB,'ISOTOPESHIN',12,NBISO,SHINA)
      ELSE
        SHINA(:NBISO)=' '
      ENDIF
      CALL LCMLEN(IPLIB,'ISOTOPESSHI',ILENG,ITYLCM)
      IF(ILENG.GT.0) THEN
        CALL LCMGET(IPLIB,'ISOTOPESSHI',LSHI)
      ELSE
        LSHI(:NBISO)=0
      ENDIF
      CALL LCMLEN(IPLIB,'ISOTOPESNIR',ILENG,ITYLCM)
      IF(ILENG.GT.0) THEN
        CALL LCMGET(IPLIB,'ISOTOPESGIR',GIR)
        CALL LCMGET(IPLIB,'ISOTOPESNIR',NIR)
      ELSE
        GIR(:NBISO)=1.0
        NIR(:NBISO)=0
      ENDIF
      CALL LCMGET(IPLIB,'ISOTOPESDENS',DENISO)
      CALL LCMGET(IPLIB,'ISOTOPESMIX',ISOMIX)
      CALL LCMLEN(IPLIB,'ISOTOPESTEMP',ILENG,ITYLCM)
      IF(ILENG.GT.0) THEN
        CALL LCMGET(IPLIB,'ISOTOPESTEMP',TMPISO)
      ELSE
        TMPISO(:NBISO)=0.0
      ENDIF
      CALL LCMLEN(IPLIB,'ISOTOPESTODO',ILENG,ITYLCM)
      IF(ILENG.GT.0) THEN
        CALL LCMGET(IPLIB,'ISOTOPESTODO',IEVOL)
      ELSE
        IEVOL(:NBISO)=0
      ENDIF
      CALL LCMLEN(IPLIB,'ISOTOPESTYPE',ILENG,ITYLCM)
      IF(ILENG.GT.0) THEN
        CALL LCMGET(IPLIB,'ISOTOPESTYPE',ITYP)
      ELSE
        ITYP(:NBISO)=1
      ENDIF
      SNISO(:NBISO)=0.0
      SBISO(:NBISO)=0.0
      JPLIB=LCMGID(IPLIB,'ISOTOPESLIST')
      DO IIISO=1,NBISO
        CALL LCMLEL(JPLIB,IIISO,ILONG,ITYLCM)
        MASKI(IIISO)=ILONG.NE.0
      ENDDO
*----
*  RECOVER MIXTURES STATES
*----
      CALL LCMLEN(IPLIB,'MIXTUREGAS',ILENG,ITYLCM)
      IF(ILENG.EQ.NBMIX) THEN
        CALL LCMGET(IPLIB,'MIXTUREGAS',KGAS)
      ELSE
        CALL XDISET(KGAS,NBMIX,0)
      ENDIF
*----
*  UNSET MIXTURES DENSITIES
*----
      DENMIX(:MAXMIX)=-1.0
*----
*  RECOVER EXTRA VECTOR EDIT NAMES
*----
      IF(NEDMAC.GT.0) THEN
        CALL LCMGTC(IPLIB,'ADDXSNAME-P0',8,NEDMAC,HVECT)
      ENDIF
*----
*  RECOVER EXTERNAL CROSS-SECTION LIBRARY NAMES
*----
      IF(NLIB.GT.0) THEN
        CALL LCMGTC(IPLIB,'ILIBRARYNAME',64,NLIB,HNAME)
      ENDIF
      RETURN
      END
