*DECK MPOCA2
      SUBROUTINE MPOCA2(IPMPO,IPEDIT,HEDIT,NREA,NISO,NADRX,NED,NPRC,
     1 ILEAK,NG,NMIL,NL,ITRANC,NALBP,IMC,NBISO,ICAL,MAXRDA,MAXIDA,
     2 FNORM,IMPX,NISOTS,NISFS,NISPS,VOLMIL,FLXMIL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover the cross sections of an elementary calculation.
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
* IPMPO   pointer to the MPO file.
* IPEDIT  pointer to the edition object (L_EDIT signature).
* HEDIT   name of output group for a (multigroup mesh, output geometry)
*         couple (generally equal to 'output_0').
* NREA    number of requested reactions.
* NISO    number of particularized isotopes.
* NADRX   total number of ADRX sets.
* NED     number of additional edition cross sections.
* NPRC    number of delayed neutron precursors.
* ILEAK   type of leakage (=0/1: off/diffusion coefficients).
* NG      number of condensed energy groups.
* NMIL    number of mixtures in the MPO file.
* NL      number of Legendre orders.
* ITRANC  type of transport correction.
* NALBP   number of physical albedos per energy group.
* IMC     type of macro-calculation (1 for diffusion or SPN;
*         2 other method).
* NBISO   number of isotopes in the condensed microlib of the edition
*         object. A given isotope may appear in many mixtures.
* ICAL    index of the current elementary calculation.
* MAXRDA  dimension of RDATAX array.
* MAXIDA  dimension of IDATAP array.
* FNORM   flux normalization factor.
* IMPX    print parameter.
*
*Parameters: output
* NISOTS  number of distinct isotopes.
* NISFS   number of particularized fissile isotopes.
* NISPS   number of particularized fission products.
* VOLMIL  mixture volumes.
* FLXMIL  averaged flux of mixtures.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE hdf5_wrap
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMPO,IPEDIT
      INTEGER NREA,NISO,NADRX,NED,NPRC,ILEAK,NG,NMIL,NL,ITRANC,NALBP,
     1 IMC,NBISO,ICAL,MAXRDA,MAXIDA,IMPX,NISOTS,NISFS,NISPS
      REAL FNORM,VOLMIL(NMIL),FLXMIL(NMIL,NG)
      CHARACTER(LEN=12) HEDIT
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NREAK=50,MAXISO=800)
      TYPE(C_PTR) JPEDIT,KPEDIT,IPTEMP,KPTEMP
      INTEGER FGYS(2),RANK,TYPE,NBYTE,DIMSR(5),ADDRZI
      CHARACTER ISOTS(MAXISO)*8,CM*2,TEXT8*8,TEXT12*12,HSMG*131,
     1 RECNAM*80
      LOGICAL EXIST,LSPH
      DOUBLE PRECISION CONV,XDRCST
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IDATAP,IFD1,IAD1,IFD2,
     1 IAD2,IJJ1,NJJ1,IPOS,IJJ2,NJJ2,MIX,ITYPE,IDATAP_MIL,VINTE1D
      INTEGER, ALLOCATABLE, DIMENSION(:) :: REACTION,ISOTOPE
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ISONAM,OUPUTID
      INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: ADRX,VINTE3D
      REAL, ALLOCATABLE, DIMENSION(:) :: RDATAX,OVERV,WORKD,WORK1,
     1 WORK2,DEN,DENISO,CONCES,DECAYC,ENERG,VREAL
      REAL, ALLOCATABLE, DIMENSION(:,:) :: DNUSIG,DCHI,DATA1,DATA2,
     1 DATA4,SPH
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: DATA3
      CHARACTER(LEN=24), ALLOCATABLE, DIMENSION(:) :: TEXT24,NOMREA,
     1 NOMISO
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: IPISO
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(ADRX(NREA+3,NISO,NADRX+NMIL),IDATAP(MAXIDA),IFD1(NG),
     1 IAD1(NG+1),IFD2(NG),IAD2(NG+1),IJJ1(NMIL),NJJ1(NMIL),
     2 IPOS(NMIL),IJJ2(NG),NJJ2(NG),ISONAM(3,NBISO),MIX(NBISO),
     3 ITYPE(NBISO),IDATAP_MIL((2*NG+1)*NISO))
      ALLOCATE(RDATAX(MAXRDA),OVERV(NG),DNUSIG(NG,NPRC+1),
     1 DCHI(NG,NPRC),WORKD(NPRC),WORK1(NG*NMIL+1),WORK2(NG),
     2 DATA1(NG,NREA),DATA2(NG,NL),DATA3(NG,NG,NL),DATA4(NG,NG),
     3 DEN(NBISO),DENISO(NISO),CONCES(NBISO),DECAYC(NBISO))
*
      CONV=1.0D6*XDRCST('eV','J')
      IF(NREA.GT.NREAK) CALL XABORT('MPOCA2: NOMREA OVERFLOW.')
*----
*  SET ENERGY MESH AND ZONE VOLUMES
*----
      CALL hdf5_read_data(IPMPO,"/energymesh/NENERGYMESH",NENERG)
      CALL hdf5_read_data(IPMPO,"/geometry/NGEOMETRY",NGEOME)
      CALL hdf5_read_data(IPMPO,"/output/OUPUTID",OUPUTID)
      READ(HEDIT,'(7X,I2)') ID
      ID_G=-1
      ID_E=-1
      DO I=1,NGEOME
        DO J=1,NENERG
          IF(OUPUTID(J,I).EQ.ID) THEN
            ID_G=I-1
            ID_E=J-1
            GO TO 10
          ENDIF
        ENDDO
      ENDDO
      CALL XABORT('MPOCA2: no ID found in /output/OUPUTID.')
   10 WRITE(RECNAM,'(23H/energymesh/energymesh_,I0,1H/)') ID_E
      CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"NG",NG2)
      CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"ENERGY",ENERG)
      IF(SIZE(ENERG,1)-1.NE.NG) CALL XABORT('MPOCA2: INVALID NG VALUE.')
      DO 20 IGR=1,NG+1
        ENERG(IGR)=ENERG(IGR)/1.0E-6
   20 CONTINUE
      WRITE(RECNAM,'(19H/geometry/geometry_,I0,1H/)') ID_G
      CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"ZONEVOLUME",VREAL)
      VOLMIL(:)=VREAL(:)
      DEALLOCATE(VREAL)
      CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"NZONE",NMIL2)
      IF(NMIL.NE.NMIL2) THEN
         WRITE(HSMG,'(42HMPOCA2: ELEMENTARY CALCULATION WITH AN INV,
     1   22HALIB NB. OF MIXTURES =,I7,3H NE,I7,1H.)') NMIL2,NMIL
         CALL XABORT(HSMG)
      ELSE IF(NG.NE.NG2) THEN
         WRITE(HSMG,'(42HMPOCA2: ELEMENTARY CALCULATION WITH AN INV,
     1   20HALIB NB. OF GROUPS =,I7,3H NE,I7,1H.)') NG2,NG
         CALL XABORT(HSMG)
      ENDIF
*----
*  CREATE DUMMY DAYASETS REACTION AND ISOTOPE
*----
      WRITE(RECNAM,'(8H/output/,A,6H/info/)') TRIM(HEDIT)
      CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"REACTION",REACTION)
      CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"ISOTOPE",ISOTOPE)
*----
*  RECOVER INFORMATION FROM THE info and contents GROUPS.
*----
      ALLOCATE(NOMREA(NREA+2),NOMISO(NISO))
      IF(NREA.GT.0) THEN
        CALL hdf5_read_data(IPMPO,"/contents/reactions/REACTIONAME",
     >  TEXT24)
        DO 30 I=1,NREA
          NOMREA(I)=TEXT24(REACTION(I)+1)
   30   continue
        DEALLOCATE(TEXT24,REACTION)
      ENDIF
      CALL hdf5_read_data(IPMPO,"/contents/isotopes/ISOTOPENAME",TEXT24)
      DO 40 I=1,NISO
        NOMISO(I)=TEXT24(ISOTOPE(I)+1)
   40 CONTINUE
      DEALLOCATE(TEXT24,ISOTOPE)
      IF(IMPX.GT.2) THEN
        WRITE(6,'(/24H MPOCA2: reaction names:)')
        DO 50 I=1,NREA
          WRITE(6,'(5X,7HNOMREA(,I3,2H)=,A)') I,TRIM(NOMREA(I))
   50   CONTINUE
        WRITE(6,'(/23H MPOCA2: isotope names:)')
        DO 60 I=1,NISO
          WRITE(6,'(5X,7HNOMISO(,I3,2H)=,A)') I,TRIM(NOMISO(I))
   60   CONTINUE
      ENDIF
*----
*  RECOVER NADRI AND IDATAP.
*  NADRI IS THE TOTAL NUMBER OF TRANSPROFILE SETS.
*----
      WRITE(RECNAM,'(8H/output/,A,6H/info/)') TRIM(HEDIT)
      NADRI=0
      CALL hdf5_info(IPMPO,TRIM(RECNAM)//"TRANSPROFILE",RANK,TYPE,NBYTE,
     1 DIMSR)
      IF(TYPE.NE.99) THEN
        NADRI=DIMSR(1)/(2*NG+1)
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"TRANSPROFILE",VINTE1D)
        IDATAP(:DIMSR(1))=VINTE1D(:DIMSR(1))
        DEALLOCATE(VINTE1D)
      ENDIF
*----
*  RECOVER INFORMATION FROM THE output_id/info GROUP.
*----
      CALL hdf5_info(IPMPO,TRIM(RECNAM)//"ADDRXS",RANK,TYPE,NBYTE,DIMSR)
      IF(TYPE.NE.99) THEN
        CALL hdf5_read_data(IPMPO,TRIM(RECNAM)//"ADDRXS",VINTE3D)
        IF(NADRX.NE.DIMSR(3)) CALL XABORT('MPOCA2: INVALID NADRX.')
        ADRX(:,:,:NADRX)=VINTE3D(:,:,:NADRX)
        DEALLOCATE(VINTE3D)
      ENDIF
*----
*  SAVE INFORMATION TO THE /output/output_id/statept_id/zone_id/yields/
*  GROUP.
*----
      WRITE(RECNAM,'(8H/output/,A,9H/statept_,I0)') TRIM(HEDIT),ICAL-1
      CALL hdf5_create_group(IPMPO,TRIM(RECNAM))
      DO 70 IMIL=1,NMIL
        WRITE(RECNAM,'(8H/output/,A,9H/statept_,I0,6H/zone_,I0,1H/)')
     >  TRIM(HEDIT),ICAL-1,IMIL-1
        NMGF=1
        CALL hdf5_create_group(IPMPO,TRIM(RECNAM))
        IF(NBISO.GT.0) THEN
          FGYS(1)=0
          FGYS(2)=1
          CALL hdf5_create_group(IPMPO,TRIM(RECNAM)//"yields")
          CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"yields/NMGF",NMGF)
          CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"yields/YIELDGROUP",
     >    FGYS)
        ENDIF
   70 CONTINUE
*----
*  FIND THE NUMBER AND NAMES OF THE ISOTOPES IN THE OUTPUT TABLES.
*----
      IF(NBISO.GT.0) THEN
         CALL LCMGET(IPEDIT,'ISOTOPESUSED',ISONAM)
         CALL LCMGET(IPEDIT,'ISOTOPESMIX',MIX)
         CALL LCMGET(IPEDIT,'ISOTOPESDENS',DEN)
         CALL LCMGET(IPEDIT,'ISOTOPESTYPE',ITYPE)
      ENDIF
      NISOTS=0
      DO 90 IBISO=1,NBISO
      IF(MIX(IBISO).EQ.0) GO TO 90
      WRITE(TEXT12,'(3A4)') (ISONAM(I0,IBISO),I0=1,3)
      DO 80 ISO=1,NISOTS
      IF(TEXT12(:8).EQ.ISOTS(ISO)) GO TO 90
   80 CONTINUE
      NISOTS=NISOTS+1
      IF(NISOTS.GT.MAXISO) CALL XABORT('MPOCA2: ISOTS OVERFLOW.')
      IF(NISOTS.GT.NBISO) CALL XABORT('MPOCA2: CONCES OVERFLOW.')
      ISOTS(NISOTS)=TEXT12(:8)
   90 CONTINUE
*----
*  RECOVER RADIOACTIVE DECAY CONSTANTS.
*----
      IF(ICAL.EQ.1) THEN
        CALL XDRSET(DECAYC,NISOTS,0.0)
        DO 120 IBISO=1,NBISO
        IF(MIX(IBISO).EQ.0) GO TO 120
        WRITE(TEXT12,'(3A4)') (ISONAM(I0,IBISO),I0=1,3)
        IISOTS=0
        DO 100 ISO=1,NISOTS
        IISOTS=ISO
        IF(TEXT12(:8).EQ.ISOTS(ISO)) GO TO 110
  100   CONTINUE
        CALL XABORT('MPOCA2: UNABLE TO FIND ISOTOPE '//TEXT12//'.')
  110   DECAYC(IISOTS)=0.0
        CALL LCMLEN(IPEDIT,'DECAY',ILONG,ITYLCM)
        IF(ILONG.EQ.1) CALL LCMGET(IPEDIT,'DECAY',DECAYC(IISOTS))
  120   CONTINUE
        DO 130 ISO=1,NISOTS
          DECAYC(ISO)=DECAYC(ISO)*1.0E-8
  130   CONTINUE
        CALL hdf5_write_data(IPMPO,"/contents/isotopes/DECAYCONST",
     1  DECAYC)
        DEALLOCATE(DECAYC)
      ENDIF
*----
*  RECOVER INVERSE OF SPH EQUIVALENCE FACTORS.
*----
      CALL LCMSIX(IPEDIT,'MACROLIB',1)
      JPEDIT=LCMGID(IPEDIT,'GROUP')
      LSPH=.FALSE.
      ALLOCATE(SPH(NMIL,NG))
      DO 160 IGR=1,NG
      KPEDIT=LCMGIL(JPEDIT,IGR)
      CALL LCMLEN(KPEDIT,'NSPH',ILONG,ITYLCM)
      IF(ILONG.GT.0) THEN
        LSPH=.TRUE.
        CALL LCMGET(KPEDIT,'NSPH',WORK1)
        DO 140 IMIL=1,NMIL
        SPH(IMIL,IGR)=1.0/WORK1(IMIL)
  140   CONTINUE
      ELSE
        DO 150 IMIL=1,NMIL
        SPH(IMIL,IGR)=1.0
  150   CONTINUE
      ENDIF
  160 CONTINUE
      CALL LCMSIX(IPEDIT,' ',2)
*----
*  CREATE A SPH-UNCORRECTED MICROLIB.
*----
      CALL LCMOP(IPTEMP,'*TEMPORARY*',0,1,0)
      ALLOCATE(IPISO(NBISO))
      CALL LCMEQU(IPEDIT,IPTEMP)
      IF(LSPH) THEN
        IF(IMC.EQ.0) CALL XABORT('MPOCA2: UNDEFINED TYPE OF SPH.')
        NW=1 ! NTOT1 cross section present
        CALL SPHCMI(IPTEMP,0,IMC,NMIL,NBISO,NG,NL,NW,NED,NPRC,NALBP,SPH)
      ENDIF
      DEALLOCATE(SPH)
*----
*  STORE INFORMATION IN THE output_id/statept_id/addons GROUP.
      WRITE(RECNAM,'(8H/output/,A,9H/statept_,I0,8H/addons/)')
     & TRIM(HEDIT),ICAL-1
      CALL hdf5_create_group(IPMPO,TRIM(RECNAM))
      CALL LCMSIX(IPTEMP,'MACROLIB',1)
      JPEDIT=LCMGID(IPTEMP,'GROUP')
      CALL LCMLEN(IPTEMP,'K-EFFECTIVE',ILONG,ITYLCM)
      IF(ILONG.EQ.1) THEN
        CALL LCMGET(IPTEMP,'K-EFFECTIVE',FLOTT)
        CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"KEFF",FLOTT)
      ENDIF
      CALL LCMLEN(IPTEMP,'K-INFINITY',ILONG,ITYLCM)
      IF(ILONG.EQ.1) THEN
        CALL LCMGET(IPTEMP,'K-INFINITY',FLOTT)
        CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"KINF",FLOTT)
      ENDIF
      CALL LCMLEN(IPTEMP,'B2  B1HOM',ILONG,ITYLCM)
      IF(ILONG.EQ.1) THEN
        CALL LCMGET(IPTEMP,'B2  B1HOM',B2)
        CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"B2",B2)
      ENDIF
      CALL LCMSIX(IPTEMP,' ',2)
*----
*  LOOP OVER MPO MIXTURES.
*----
      DO 920 IMIL=1,NMIL
      IF(NADRX+1.GT.SIZE(ADRX,3)) CALL XABORT('MPOCA2: ADRX OVERFLOW.')
      IOI=0
      IOR=0
      DO 165 IGR=1,NG
      IFD1(IGR)=NG+1
      IAD1(IGR+1)=0
  165 CONTINUE
      CALL XDRSET(DATA2,NG*NL,0.0)
      CALL XDRSET(DATA3,NG*NG*NL,0.0)
      CALL LCMSIX(IPTEMP,'MACROLIB',1)
      DO 230 IGR=1,NG
      KPEDIT=LCMGIL(JPEDIT,IGR)
*----
*  RECOVER THE NEUTRON FLUX.
*----
      CALL LCMGET(KPEDIT,'FLUX-INTG',WORK1)
      IF(FNORM.NE.1.0) THEN
        FLXMIL(IMIL,IGR)=WORK1(IMIL)*FNORM*1.0E13
      ELSE
        FLXMIL(IMIL,IGR)=WORK1(IMIL)
      ENDIF
*----
*  RECOVER DELAYED NEUTRON INFORMATION.
*----
      CALL LCMLEN(KPEDIT,'NUSIGF',ILONG,ITYLCM)
      IF((NPRC.GT.0).AND.(ILONG.NE.0)) THEN
         CALL LCMGET(KPEDIT,'NUSIGF',WORK1)
         DNUSIG(IGR,NPRC+1)=WORK1(IMIL)
         CALL LCMGET(KPEDIT,'OVERV',WORK1)
         OVERV(IGR)=WORK1(IMIL)
         DO 170 IPRC=1,NPRC
         WRITE(TEXT12,'(6HNUSIGF,I2.2)') IPRC
         CALL LCMGET(KPEDIT,TEXT12,WORK1)
         DNUSIG(IGR,IPRC)=WORK1(IMIL)
         WRITE(TEXT12,'(3HCHI,I2.2)') IPRC
         CALL LCMGET(KPEDIT,TEXT12,WORK1)
         DCHI(IGR,IPRC)=WORK1(IMIL)
  170    CONTINUE
      ELSE
         CALL XDRSET(DNUSIG,NG*(NPRC+1),0.0)
      ENDIF
*
      DO 220 IREA=1,NREA
      DATA1(IGR,IREA)=0.0
      IF(NOMREA(IREA).EQ.'Total') THEN
         CALL LCMGET(KPEDIT,'NTOT0',WORK1)
         DATA1(IGR,IREA)=WORK1(IMIL)
      ELSE IF(NOMREA(IREA).EQ.'TotalP1') THEN
         CALL LCMGET(KPEDIT,'NTOT1',WORK1)
         DATA1(IGR,IREA)=WORK1(IMIL)
      ELSE IF(NOMREA(IREA).EQ.'Absorption') THEN
         CALL LCMGET(KPEDIT,'NTOT0',WORK1)
         DATA1(IGR,IREA)=WORK1(IMIL)
         CALL LCMLEN(KPEDIT,'SIGS00',ILONG,ITYLCM)
         IF(ILONG.GT.0) THEN
            CALL LCMGET(KPEDIT,'SIGS00',WORK1)
            DATA1(IGR,IREA)=DATA1(IGR,IREA)-WORK1(IMIL)
         ENDIF
         CALL LCMLEN(KPEDIT,'Nexcess',ILONG,ITYLCM)
         IF(ILONG.GT.0) THEN
            CALL LCMGET(KPEDIT,'N2N',WORK1)
            DATA1(IGR,IREA)=DATA1(IGR,IREA)+WORK1(IMIL)
         ENDIF
         CALL LCMLEN(KPEDIT,'N3N',ILONG,ITYLCM)
         IF(ILONG.GT.0) THEN
            CALL LCMGET(KPEDIT,'N3N',WORK1)
            DATA1(IGR,IREA)=DATA1(IGR,IREA)+2.0*WORK1(IMIL)
         ENDIF
      ELSE IF(NOMREA(IREA).EQ.'Fission') THEN
         CALL LCMLEN(KPEDIT,'NFTOT',ILONG,ITYLCM)
         IF(ILONG.GT.0) THEN
            CALL LCMGET(KPEDIT,'NFTOT',WORK1)
            DATA1(IGR,IREA)=WORK1(IMIL)
         ENDIF
      ELSE IF(NOMREA(IREA).EQ.'FissionSpectrum') THEN
         CALL LCMLEN(KPEDIT,'CHI',ILONG,ITYLCM)
         IF(ILONG.GT.0) THEN
            CALL LCMGET(KPEDIT,'CHI',WORK1)
            DATA1(IGR,IREA)=WORK1(IMIL)
         ENDIF
      ELSE IF(NOMREA(IREA).EQ.'NuFission') THEN
         CALL LCMLEN(KPEDIT,'NUSIGF',ILONG,ITYLCM)
         IF(ILONG.GT.0) THEN
            CALL LCMGET(KPEDIT,'NUSIGF',WORK1)
            DATA1(IGR,IREA)=WORK1(IMIL)
         ENDIF
      ELSE IF(NOMREA(IREA).EQ.'Energy') THEN
         CALL LCMLEN(KPEDIT,'H-FACTOR',ILONG,ITYLCM)
         IF(ILONG.GT.0) THEN
            CALL LCMGET(KPEDIT,'H-FACTOR',WORK1)
            DATA1(IGR,IREA)=WORK1(IMIL)/REAL(CONV)
         ENDIF
      ELSE IF(NOMREA(IREA).EQ.'FUITES') THEN
         CALL LCMLEN(KPEDIT,'DIFF',ILONG,ITYLCM)
         IF(ILONG.GT.0) THEN
            IF(B2.EQ.0.0) B2=1.0E-10
            CALL LCMGET(KPEDIT,'DIFF',WORK1)
            DATA1(IGR,IREA)=WORK1(IMIL)*B2
         ENDIF
      ELSE IF(NOMREA(IREA).EQ.'STRD') THEN
         CALL LCMLEN(KPEDIT,'DIFF',ILONG,ITYLCM)
         IF(ILONG.GT.0) THEN
            CALL LCMGET(KPEDIT,'DIFF',WORK1)
            DATA1(IGR,IREA)=1.0/(3.0*WORK1(IMIL))
         ENDIF
      ELSE IF(NOMREA(IREA).EQ.'Diffusion') THEN
         DO 180 IL=1,NL
         WRITE (CM,'(I2.2)') IL-1
         CALL LCMGET(KPEDIT,'SIGS'//CM,WORK1)
         DATA2(IGR,IL)=WORK1(IMIL)
  180    CONTINUE
         CALL LCMLEN(KPEDIT,'N2N',ILONG,ITYLCM)
         IF(ILONG.GT.0) THEN
            CALL LCMGET(KPEDIT,'N2N',WORK1)
            DATA2(IGR,1)=DATA2(IGR,1)-WORK1(IMIL)
         ENDIF
         CALL LCMLEN(KPEDIT,'N3N',ILONG,ITYLCM)
         IF(ILONG.GT.0) THEN
            CALL LCMGET(KPEDIT,'N3N',WORK1)
            DATA2(IGR,1)=DATA2(IGR,1)-2.0*WORK1(IMIL)
         ENDIF
      ELSE IF(NOMREA(IREA).EQ.'Transport') THEN
         IF((ITRANC.EQ.1).AND.(NL.GE.2)) THEN
            CALL LCMGET(KPEDIT,'SIGS01',WORK1)
            DATA1(IGR,IREA)=WORK1(IMIL)
         ELSE IF(ITRANC.EQ.2) THEN
            CALL LCMGET(KPEDIT,'TRANC',WORK1)
            DATA1(IGR,IREA)=WORK1(IMIL)
         ENDIF
      ELSE IF(NOMREA(IREA).EQ.'Scattering') THEN
         DO 190 IL=1,NL
         WRITE (CM,'(I2.2)') IL-1
         CALL LCMLEN(KPEDIT,'IJJS'//CM,ILONG,ITYLCM)
         IF(ILONG.EQ.0) GO TO 190
         CALL LCMGET(KPEDIT,'IJJS'//CM,IJJ1)
         CALL LCMGET(KPEDIT,'NJJS'//CM,NJJ1)
         DO 185 JGR=IJJ1(IMIL)-NJJ1(IMIL)+1,IJJ1(IMIL) ! IGR <-- JGR
         IFD1(JGR)=MIN(IFD1(JGR),IGR)
         IAD1(JGR+1)=MAX(IAD1(JGR+1),IGR)
  185    CONTINUE
  190    CONTINUE
         DO 210 IL=1,NL
         WRITE (CM,'(I2.2)') IL-1
         CALL LCMGET(KPEDIT,'IJJS'//CM,IJJ1)
         CALL LCMGET(KPEDIT,'NJJS'//CM,NJJ1)
         CALL LCMGET(KPEDIT,'IPOS'//CM,IPOS)
         CALL LCMGET(KPEDIT,'SCAT'//CM,WORK1)
         IPO=IPOS(IMIL)
         J2=IJJ1(IMIL)
         J1=IJJ1(IMIL)-NJJ1(IMIL)+1
         DO 200 JGR=J2,J1,-1
         DATA3(IGR,JGR,IL)=WORK1(IPO)*REAL(2*IL-1)
         IPO=IPO+1
  200    CONTINUE
  210    CONTINUE
      ELSE
         CALL LCMLEN(KPEDIT,NOMREA(IREA)(:12),ILONG,ITYLCM)
         IF(ILONG.GT.0) THEN
            CALL LCMGET(KPEDIT,NOMREA(IREA),WORK1)
            DATA1(IGR,IREA)=WORK1(IMIL)
         ENDIF
      ENDIF
  220 CONTINUE
  230 CONTINUE
      IAD1(1)=0
      DO 235 IGR=1,NG
      IAD1(IGR+1)=IAD1(IGR)+(IAD1(IGR+1)-IFD1(IGR)+1)
  235 CONTINUE
      CALL LCMSIX(IPTEMP,' ',2)
*----
*  FIND ISOTOPE POINTERS IN INPUT MICROLIB
*----
      IF(NBISO.GT.0) THEN
        CALL LIBIPS(IPTEMP,NBISO,IPISO)
*----
*  PROCESS PARTICULARIZED ISOTOPES
*----
        DO 250 IISO=1,NISO
        DO 240 IREA=1,NREA+3
        ADRX(IREA,IISO,NADRX+1)=-1
  240   CONTINUE
  250   CONTINUE
        CALL XDRSET(CONCES,NISOTS,0.0)
        DO 540 IBISO=1,NBISO
        IF(MIX(IBISO).EQ.IMIL) THEN
          WRITE(TEXT12,'(3A4)') (ISONAM(I0,IBISO),I0=1,3)
          DO 260 ISO=1,NISO
          IISO=ISO
          IF(NOMISO(ISO).EQ.TEXT12(:8)) GO TO 270
  260     CONTINUE
          GO TO 540
  270     IF(IISO.GT.NISO-1) CALL XABORT('MPOCA2: NISO OVERFLOW.')
          KPTEMP=IPISO(IBISO) ! set IBISO-th isotope
          IF(.NOT.C_ASSOCIATED(KPTEMP)) THEN
            WRITE(HSMG,'(17HMPOCA2: ISOTOPE '',A12,7H'' (ISO=,I8,3H) I,
     1      32HS NOT AVAILABLE IN THE MICROLIB.)') TEXT12,IBISO
            CALL XABORT(HSMG)
          ENDIF
          IISOTS=0
          DO 280 ISO=1,NISOTS
          IISOTS=ISO
          IF(ISOTS(ISO).EQ.TEXT12(:8)) GO TO 290
  280     CONTINUE
          CALL XABORT('MPOCA2: UNABLE TO FIND ISOTOPE '//TEXT12//'.')
  290     CONCES(IISOTS)=DEN(IBISO)
          DENISO(IISO)=DEN(IBISO)
          DO 530 IREA=1,NREA
          CALL XDRSET(WORK2,NG,0.0)
          IF(NOMREA(IREA).EQ.'Total') THEN
             CALL LCMGET(KPTEMP,'NTOT0',WORK2)
          ELSE IF(NOMREA(IREA).EQ.'TotalP1') THEN
             CALL LCMGET(KPTEMP,'NTOT1',WORK2)
          ELSE IF(NOMREA(IREA).EQ.'Absorption') THEN
             CALL LCMGET(KPTEMP,'NTOT0',WORK2)
             CALL LCMLEN(KPTEMP,'SIGS00',ILONG,ITYLCM)
             IF(ILONG.GT.0) THEN
               CALL LCMGET(KPTEMP,'SIGS00',WORK1)
               DO 300 IGR=1,NG
               WORK2(IGR)=WORK2(IGR)-WORK1(IGR)
  300          CONTINUE
             ENDIF
             CALL LCMLEN(KPTEMP,'N2N',ILONG,ITYLCM)
             IF(ILONG.GT.0) THEN
               CALL LCMGET(KPTEMP,'N2N',WORK1)
               DO 310 IGR=1,NG
               WORK2(IGR)=WORK2(IGR)+WORK1(IGR)
  310          CONTINUE
             ENDIF
             CALL LCMLEN(KPTEMP,'N3N',ILONG,ITYLCM)
             IF(ILONG.GT.0) THEN
               CALL LCMGET(KPTEMP,'N3N',WORK1)
               DO 320 IGR=1,NG
               WORK2(IGR)=WORK2(IGR)+2.0*WORK1(IGR)
  320          CONTINUE
             ENDIF
           ELSE IF(NOMREA(IREA).EQ.'Nexcess') THEN
             CALL LCMLEN(KPTEMP,'N2N',ILONG,ITYLCM)
             IF(ILONG.GT.0) CALL LCMGET(KPTEMP,'N2N',WORK2)
             CALL LCMLEN(KPTEMP,'N3N',ILONG,ITYLCM)
             IF(ILONG.GT.0) THEN
               CALL LCMGET(KPTEMP,'N3N',WORK1)
               DO 330 IGR=1,NG
               WORK2(IGR)=WORK2(IGR)+2.0*WORK1(IGR)
  330          CONTINUE
             ENDIF
           ELSE IF(NOMREA(IREA).EQ.'Fission') THEN
             CALL LCMLEN(KPTEMP,'NFTOT',ILONG,ITYLCM)
             IF(ILONG.GT.0) CALL LCMGET(KPTEMP,'NFTOT',WORK2)
           ELSE IF(NOMREA(IREA).EQ.'FissionSpectrum') THEN
             CALL LCMLEN(KPTEMP,'CHI',ILONG,ITYLCM)
             IF(ILONG.GT.0) CALL LCMGET(KPTEMP,'CHI',WORK2)
           ELSE IF(NOMREA(IREA).EQ.'NuFission') THEN
             CALL LCMLEN(KPTEMP,'NUSIGF',ILONG,ITYLCM)
             IF(ILONG.GT.0) CALL LCMGET(KPTEMP,'NUSIGF',WORK2)
           ELSE IF(NOMREA(IREA).EQ.'Energy') THEN
             CALL LCMLEN(KPTEMP,'MEVF',ILONG,ITYLCM)
             IF(ILONG.GT.0) THEN
               CALL LCMGET(KPTEMP,'NFTOT',WORK2)
               CALL LCMGET(KPTEMP,'MEVF',FLOTT)
               DO340 IGR=1,NG
               WORK2(IGR)=WORK2(IGR)*FLOTT
  340          CONTINUE
             ENDIF
             CALL LCMLEN(KPTEMP,'MEVG',ILONG,ITYLCM)
             IF(ILONG.GT.0) THEN
               CALL LCMGET(KPTEMP,'NG',WORK1)
               CALL LCMGET(KPTEMP,'MEVG',FLOTT)
               DO 350 IGR=1,NG
               WORK2(IGR)=WORK2(IGR)+WORK1(IGR)*FLOTT
  350          CONTINUE
             ENDIF
           ELSE IF(NOMREA(IREA).EQ.'FissionEnergyFission') THEN
             CALL LCMLEN(KPTEMP,'MEVF',ILONG,ITYLCM)
             IF(ILONG.GT.0) THEN
               CALL LCMGET(KPTEMP,'NFTOT',WORK2)
               CALL LCMGET(KPTEMP,'MEVF',FLOTT)
               DO 360 IGR=1,NG
               WORK2(IGR)=WORK2(IGR)*FLOTT
  360          CONTINUE
             ENDIF
           ELSE IF(NOMREA(IREA).EQ.'CaptureEnergyCapture') THEN
             CALL LCMLEN(KPTEMP,'MEVG',ILONG,ITYLCM)
             IF(ILONG.GT.0) THEN
               CALL LCMGET(KPTEMP,'NG',WORK2)
               CALL LCMGET(KPTEMP,'MEVG',FLOTT)
               DO 370 IGR=1,NG
               WORK2(IGR)=WORK2(IGR)*FLOTT
  370          CONTINUE
             ENDIF
           ELSE IF(NOMREA(IREA).EQ.'STRD') THEN
             CALL LCMLEN(KPTEMP,'STRD',ILONG,ITYLCM)
             IF(ILONG.GT.0) CALL LCMGET(KPTEMP,'STRD',WORK2)
           ELSE IF(NOMREA(IREA).EQ.'Diffusion') THEN
             ADRX(IREA,IISO,NADRX+1)=IOR
             ADRX(NREA+1,IISO,NADRX+1)=NL
             IOR=IOR+NG*NL
             IF(IOR.GT.MAXRDA) CALL XABORT('MPOCA2: RDATAX OVERFLOW(1)')
             DO 420 IL=1,NL
             WRITE (CM,'(I2.2)') IL-1
             CALL LCMLEN(KPTEMP,'SIGS'//CM,ILONG,ITYLCM)
             IF(ILONG.GT.0) THEN
               CALL LCMGET(KPTEMP,'SIGS'//CM,WORK2)
             ELSE
               CALL XDRSET(WORK2,NG,0.0)
             ENDIF
             CALL LCMLEN(KPTEMP,'N2N',ILONG,ITYLCM)
             IF((IL.EQ.1).AND.(ILONG.GT.0)) THEN
               CALL LCMGET(KPTEMP,'N2N',WORK1)
               DO 390 IGR=1,NG
               WORK2(IGR)=WORK2(IGR)-WORK1(IGR)
  390          CONTINUE
             ENDIF
             CALL LCMLEN(KPTEMP,'N3N',ILONG,ITYLCM)
             IF((IL.EQ.1).AND.(ILONG.GT.0)) THEN
               CALL LCMGET(KPTEMP,'N3N',WORK1)
               DO 400 IGR=1,NG
               WORK2(IGR)=WORK2(IGR)-2.0*WORK1(IGR)
  400          CONTINUE
             ENDIF
             DO 410 IGR=1,NG
             RDATAX(ADRX(IREA,IISO,NADRX+1)+(IL-1)*NG+IGR-1)=WORK2(IGR)
  410        CONTINUE
  420        CONTINUE
             GO TO 530
           ELSE IF(NOMREA(IREA).EQ.'Transport') THEN
             IF((ITRANC.EQ.1).AND.(NL.GE.2)) THEN
               CALL LCMGET(KPTEMP,'SIGS01',WORK2)
             ELSE IF(ITRANC.EQ.2) THEN
               CALL LCMGET(KPTEMP,'TRANC',WORK2)
             ENDIF
           ELSE IF(NOMREA(IREA).EQ.'Scattering') THEN
             DO 430 IGR=1,NG
             IFD2(IGR)=NG+1
             IAD2(IGR+1)=0
  430        CONTINUE
             DO 450 IL=1,NL
             WRITE (CM,'(I2.2)') IL-1
             CALL LCMLEN(KPTEMP,'IJJS'//CM,ILONG,ITYLCM)
             IF(ILONG.EQ.0) GO TO 450
             CALL LCMGET(KPTEMP,'IJJS'//CM,IJJ2)
             CALL LCMGET(KPTEMP,'NJJS'//CM,NJJ2)
             DO 445 JGR=1,NG
             DO 440 IGR=IJJ2(JGR)-NJJ2(JGR)+1,IJJ2(JGR) ! JGR <-- IGR
             IFD2(IGR)=MIN(IFD2(IGR),JGR)
             IAD2(IGR+1)=MAX(IAD2(IGR+1),JGR)
  440        CONTINUE
  445        CONTINUE
  450        CONTINUE
             IAD2(1)=0
             DO 460 IGR=1,NG
             IAD2(IGR+1)=IAD2(IGR)+(IAD2(IGR+1)-IFD2(IGR)+1)
  460        CONTINUE
             ADRX(NREA+1,IISO,NADRX+1)=NL
             ADRX(NREA+2,IISO,NADRX+1)=NL
             ADRX(NREA+3,IISO,NADRX+1)=IOI
             IF(IOI+2*NG+1.GT.(2*NG+1)*NISO) THEN
               CALL XABORT('MPOCA2: IDATAP_MIL OVERFLOW(1).')
             ENDIF
             DO 470 IGR=1,NG
             IDATAP_MIL(IOI+IGR)=IFD2(IGR)-1
             IDATAP_MIL(IOI+NG+IGR)=IAD2(IGR)
  470        CONTINUE
             IDATAP_MIL(IOI+2*NG+1)=IAD2(NG+1)
             ADRX(NREA+3,IISO,NADRX+1)=IOI
             IOI=IOI+2*NG+1
*
             ADRX(IREA,IISO,NADRX+1)=IOR
             IOR=IOR+IAD2(NG+1)*NL
             IF(IOR.GT.MAXRDA) CALL XABORT('MPOCA2: RDATAX OVERFLOW(2)')
             JOFS=0
             DO 500 IL=1,NL
             IMPX=0
             CALL XDRLGS(KPTEMP,-1,IMPX,IL-1,IL-1,1,NG,WORK2,DATA4,
     1       ITYPRO)
             ZIL=REAL(2*IL-1)
             DO 490 IGR=1,NG
             DO 480 JGR=IFD2(IGR),IFD2(IGR)+(IAD2(IGR+1)-IAD2(IGR))-1 ! JGR <-- IGR
             JOFS=JOFS+1
             RDATAX(ADRX(IREA,IISO,NADRX+1)+JOFS-1)=DATA4(JGR,IGR)*ZIL
  480        CONTINUE
  490        CONTINUE
  500        CONTINUE
             GO TO 530
           ELSE
             CALL LCMLEN(KPTEMP,NOMREA(IREA),ILONG,ITYLCM)
             IF(ILONG.GT.0) CALL LCMGET(KPTEMP,NOMREA(IREA),WORK2)
           ENDIF
*
           EXIST=.FALSE.
           DO 510 IGR=1,NG
           EXIST=EXIST.OR.(WORK2(IGR).NE.0.0)
  510      CONTINUE
           IF(EXIST) THEN
             ADRX(IREA,IISO,NADRX+1)=IOR
             IOR=IOR+NG
             IF(IOR.GT.MAXRDA) CALL XABORT('MPOCA2: RDATAX OVERFLOW(3)')
             DO 520 IGR=1,NG
             RDATAX(ADRX(IREA,IISO,NADRX+1)+IGR)=WORK2(IGR)
  520        CONTINUE
           ELSE
             ADRX(IREA,IISO,NADRX+1)=-1
           ENDIF
  530      CONTINUE
        ENDIF
  540   CONTINUE
      ENDIF
*----
*  STORE MACROSCOPIC RESIDUAL (ISOTOPE NISO) CROSS SECTIONS IN RDATAX.
*----
      ADRX(NREA+1,NISO,NADRX+1)=0
      ADRX(NREA+2,NISO,NADRX+1)=0
      ADRX(NREA+3,NISO,NADRX+1)=0
      DO 680 IREA=1,NREA
      IF(NOMREA(IREA).EQ.'Diffusion') THEN
         ADRX(IREA,NISO,NADRX+1)=IOR
         ADRX(NREA+1,NISO,NADRX+1)=NL
         IOR=IOR+NG*NL
         IF(IOR.GT.MAXRDA) CALL XABORT('MPOCA2: RDATAX OVERFLOW(4)')
         JOFS=0
         DO 570 IL=1,NL
         DO 560 IGR=1,NG
         JOFS=JOFS+1
         RDATAX(ADRX(IREA,NISO,NADRX+1)+JOFS)=DATA2(IGR,IL)
  560    CONTINUE
  570    CONTINUE
      ELSE IF(NOMREA(IREA).EQ.'Scattering') THEN
         ADRX(NREA+2,NISO,NADRX+1)=NL
         ADRX(NREA+3,NISO,NADRX+1)=IOI
         IF(IOI+2*NG+1.GT.(2*NG+1)*NISO) THEN
           CALL XABORT('MPOCA2: IDATAP_MIL OVERFLOW(2).')
         ENDIF
         DO 590 IGR=1,NG
         IDATAP_MIL(IOI+IGR)=IFD1(IGR)-1
         IDATAP_MIL(IOI+NG+IGR)=IAD1(IGR)
  590    CONTINUE
         IDATAP_MIL(IOI+2*NG+1)=IAD1(NG+1)
         ADRX(NREA+3,NISO,NADRX+1)=IOI
         IOI=IOI+2*NG+1
*
         ADRX(IREA,NISO,NADRX+1)=IOR
         IOR=IOR+IAD1(NG+1)*NL
         IF(IOR.GT.MAXRDA) CALL XABORT('MPOCA2: RDATAX OVERFLOW(5)')
         JOFS=0
         DO 630 IL=1,NL
         DO 620 IGR=1,NG
         DO 610 JGR=IFD1(IGR),IFD1(IGR)+(IAD1(IGR+1)-IAD1(IGR))-1 ! JGR <-- IGR
         JOFS=JOFS+1
         RDATAX(ADRX(IREA,NISO,NADRX+1)+JOFS)=DATA3(JGR,IGR,IL)
  610    CONTINUE
  620    CONTINUE
  630    CONTINUE
      ELSE
         EXIST=.FALSE.
         DO 650 IGR=1,NG
         EXIST=EXIST.OR.(DATA1(IGR,IREA).NE.0.0)
  650    CONTINUE
         IF(EXIST) THEN
            ADRX(IREA,NISO,NADRX+1)=IOR
            IOR=IOR+NG
            IF(IOR.GT.MAXRDA) CALL XABORT('MPOCA2: RDATAX OVERFLOW(6)')
            DO 660 IGR=1,NG
            RDATAX(ADRX(IREA,NISO,NADRX+1)+IGR)=DATA1(IGR,IREA)
  660       CONTINUE
          ELSE
            ADRX(IREA,NISO,NADRX+1)=-1
         ENDIF
      ENDIF
  680 CONTINUE
*----
*  REMOVE PARTICULARIZED ISOTOPIC CONTRIBUTIONS FROM MACROS.
*  ISOTOPE NISO IS THE MACROSCOPIC RESIDUAL.
*----
      IF(NBISO.GT.0) THEN
        DO 750 IREA=1,NREA
        IMACR=ADRX(IREA,NISO,NADRX+1)
        IF(IMACR+(IAD1(NG+1)-1)*NL-1.GT.MAXRDA) THEN
          CALL XABORT('MPOCA2: RDATAX OVERFLOW(6).')
        ENDIF
        IF(IMACR.EQ.-1) GO TO 750
        IGRTOT=NG
        IF(NOMREA(IREA).EQ.'Diffusion') IGRTOT=NG*NL
        IF(NOMREA(IREA).EQ.'FissionSpectrum') GO TO 750
        DO 740 IISO=1,NISO-1
        IF(DENISO(IISO).EQ.0.0) GO TO 740
        JMACR=ADRX(IREA,IISO,NADRX+1)
        IF(JMACR.EQ.-1) GO TO 740
        IF(NOMREA(IREA).EQ.'Scattering') THEN
          IOI=ADRX(NREA+3,IISO,NADRX+1)
          DO 690 IGR=1,NG
          IFD2(IGR)=IDATAP_MIL(IOI+IGR)+1
          IAD2(IGR)=IDATAP_MIL(IOI+NG+IGR)
  690     CONTINUE
          IAD2(NG+1)=IDATAP_MIL(IOI+2*NG+1)
          JOFS=0
          DO 720 IL=1,NL
          DO 710 IGR=1,NG
          DO 700 JGR=IFD2(IGR),IFD2(IGR)+(IAD2(IGR+1)-IAD2(IGR)) ! JGR <-- IGR
          I=(IL-1)*(IAD1(NG+1)-1)+IAD1(IGR)+JGR-IFD1(IGR)
          JOFS=JOFS+1
          RDATAX(IMACR+I-1)=RDATAX(IMACR+I-1)-DENISO(IISO)*
     1                      RDATAX(JMACR+JOFS-1)
  700     CONTINUE
  710     CONTINUE
  720     CONTINUE
        ELSE
          DO 730 IGR=1,IGRTOT
          RDATAX(IMACR+IGR-1)=RDATAX(IMACR+IGR-1)-DENISO(IISO)*
     1                        RDATAX(JMACR+IGR-1)
  730     CONTINUE
        ENDIF
  740   CONTINUE
  750   CONTINUE
      ENDIF
      DENISO(NISO)=1.0
*----
*  TRY TO FIND AN EXISTING IDATAP SET. OTHERWISE, CREATE A NEW ONE.
*  STORE INFORMATION IN THE ADRX(NREA+3,IISO,NADRX+1) DATASET.
*  NADRI IS THE TOTAL NUMBER OF TRANSPROFILE SETS.
*----
      DO 780 IISO=1,NISO
        IOI=ADRX(NREA+3,IISO,NADRX+1)
        DO 770 IAD1X=0,NADRI-1
          DO 760 I=1,2*NG+1
          IF(IDATAP_MIL(IOI+I).NE.IDATAP(IAD1X*(2*NG+1)+I)) GO TO 770
  760     CONTINUE
          ADRX(NREA+3,IISO,NADRX+1)=IAD1X*(2*NG+1)
          GO TO 780
  770   CONTINUE
        IF((NADRI+1)*(2*NG+1).GT.MAXIDA) THEN
          CALL XABORT('MPOCA2: IDATAP OVERFLOW.')
        ENDIF
        DO I=1,2*NG+1
          IDATAP(NADRI*(2*NG+1)+I)=IDATAP_MIL(IOI+I)
        ENDDO
        ADRX(NREA+3,IISO,NADRX+1)=NADRI*(2*NG+1)
        NADRI=NADRI+1
  780 CONTINUE
*----
*  TRY TO FIND AN EXISTING ADRX SET. OTHERWISE, CREATE A NEW ONE.
*  STORE INFORMATION IN THE output_id/statept_id/zone_id GROUP.
*  "ADDRZI" is the index in ADDRISO[NADDRISO+1]-->ISOTOPE
*  "ADDRZX" is the index in ADDRXS[NREA+3,NISO,NADRX+1]-->CROSSEXTION
*----
      WRITE(RECNAM,'(8H/output/,A,9H/statept_,I0,6H/zone_,I0,1H/)')
     1 TRIM(HEDIT),ICAL-1,IMIL-1
      DO 810 IAD1X=1,NADRX
      DO 800 I=1,NREA+3
      DO 790 J=1,NISO
      IF(ADRX(I,J,NADRX+1).NE.ADRX(I,J,IAD1X)) GO TO 810
  790 CONTINUE
  800 CONTINUE
      CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"ADDRZX",IAD1X-1)
      GO TO 820
  810 CONTINUE
      NADRX=NADRX+1
      CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"ADDRZX",NADRX-1)
  820 ADDRZI=0
      CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"ADDRZI",ADDRZI)
*----
*  STORE FLUX, CROSS SECTIONS AND NUMBER DENSITIES.
*----
      WORK2(:)=FLXMIL(IMIL,:)
      CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"ZONEFLUX",WORK2)
      IF(IOR.GT.0) THEN
        CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"CROSSECTION",
     1  RDATAX(:IOR))
      ENDIF
      CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"CONCENTRATION",DENISO)
*----
*  STORE INFORMATION IN THE output_id/statept_id/zone_id/leakage GROUP.
*----
      IF(ILEAK.EQ.1) THEN
         WRITE(RECNAM,'(8H/output/,A,9H/statept_,I0,6H/zone_,I0,
     1   9H/leakage/)') TRIM(HEDIT),ICAL-1,IMIL-1
        DO 830 IGR=1,NG
        KPEDIT=LCMGIL(JPEDIT,IGR)
        CALL LCMLEN(KPEDIT,'DIFF',ILONG,ITYLCM)
        IF(ILONG.EQ.0) CALL XABORT('MPOCA2: MISSING DIFF INFO.')
        CALL LCMGET(KPEDIT,'DIFF',WORK1)
        WORK2(IGR)=WORK1(IMIL)
  830   CONTINUE
        CALL hdf5_create_group(IPMPO,TRIM(RECNAM))
        CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"BUCKLING",B2)
        CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"DIFFCOEF",WORK2)
        WORK2(:)=WORK2(:)*B2
        CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"DB2",WORK2)
      ENDIF
*----
*  STORE INFORMATION IN THE output_id/statept_id/zone_id/kinetics GROUP.
*----
      IF(NPRC.GT.0) THEN
        EXIST=.FALSE.
        DO 850 IPRC=1,NPRC
        DO 840 IGR=1,NG
        EXIST=EXIST.OR.(DNUSIG(IGR,IPRC).NE.0.0)
  840   CONTINUE
  850   CONTINUE
        WRITE(RECNAM,'(8H/output/,A,9H/statept_,I0,6H/zone_,I0,
     1  10H/kinetics/)') TRIM(HEDIT),ICAL-1,IMIL-1
        IF(EXIST) THEN
          CALL LCMSIX(IPTEMP,'MACROLIB',1)
          CALL LCMGET(IPTEMP,'LAMBDA-D',WORKD)
          CALL LCMSIX(IPTEMP,' ',2)
          CALL hdf5_create_group(IPMPO,TRIM(RECNAM))
          CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"LAMBDAD",WORKD)
          CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"CHID",DCHI)
          CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"INVERSESPEED",
     1    OVERV)
          TGENRS=0.0
          DENOM=0.0
          DO 860 IGR=1,NG
          TGENRS=TGENRS+OVERV(IGR)*FLXMIL(IMIL,IGR)
          DENOM=DENOM+DNUSIG(IGR,NPRC+1)*FLXMIL(IMIL,IGR)
  860     CONTINUE
          TGENRS=TGENRS/DENOM
          DO 880 IPRC=1,NPRC
          WORKD(IPRC)=0.0
          DO 870 IGR=1,NG
          WORKD(IPRC)=WORKD(IPRC)+DNUSIG(IGR,IPRC)*FLXMIL(IMIL,IGR)
  870     CONTINUE
          WORKD(IPRC)=WORKD(IPRC)/DENOM
  880     CONTINUE
          CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"BETADF",WORKD)
          CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"GENERATIONTIME",
     1    TGENRS)
        ENDIF
      ENDIF
*----
*  STORE INFORMATION IN THE output_id/statept_id/zone_id/yields GROUP.
*----
      NISFS=0
      NISPS=0
      IF(NBISO.GT.0) THEN
        DO 910 ISO=1,NISO-1
        DO 890 IBISO=1,NBISO
        WRITE(TEXT8,'(2A4)') (ISONAM(I0,IBISO),I0=1,2)
        IF(NOMISO(ISO).EQ.TEXT8) THEN
          ITY=ITYPE(IBISO)
          GO TO 900
        ENDIF
  890   CONTINUE
        GO TO 910
  900   IF(ITY.EQ.2) THEN
          NISFS=NISFS+1
        ELSE IF(ITY.EQ.3) THEN
          NISPS=NISPS+1
        ENDIF
  910   CONTINUE
        NISFS=NISFS+1 ! declare the residual as fissile
        WRITE(RECNAM,'(8H/output/,A,9H/statept_,I0,6H/zone_,I0,1H/)')
     >  TRIM(HEDIT),ICAL-1,IMIL-1
        CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"yields/NISF",NISFS)
        CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"yields/NISP",NISPS)
      ENDIF
*----
*  END OF LOOP OVER MPO MIXTURES.
*----
  920 CONTINUE
      DEALLOCATE(IPISO)
      CALL LCMCL(IPTEMP,2)
*----
*  STORE INFORMATION IN THE output_id/info GROUP.
*----
      WRITE(RECNAM,'(8H/output/,A,6H/info/)') TRIM(HEDIT)
      IF((ICAL.EQ.1).AND.(NADRI.GT.0)) THEN
        CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"TRANSPROFILE",
     1  IDATAP(:NADRI*(2*NG+1)))
      ELSE IF(ICAL.GT.1) THEN
        CALL hdf5_delete(IPMPO,TRIM(RECNAM)//"NADDRXS")
        CALL hdf5_delete(IPMPO,TRIM(RECNAM)//"ADDRXS")
      ENDIF
      CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"NADDRXS",NADRX)
      CALL hdf5_write_data(IPMPO,TRIM(RECNAM)//"ADDRXS",
     1 ADRX(:,:,:NADRX))
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(NOMISO,NOMREA)
      DEALLOCATE(CONCES,DENISO,DEN,DATA4,DATA3,DATA2,DATA1,WORK2,WORK1,
     1 WORKD,DCHI,DNUSIG,OVERV,RDATAX)
      DEALLOCATE(IDATAP_MIL,ITYPE,MIX,ISONAM,NJJ2,IJJ2,IPOS,NJJ1,IJJ1,
     1 IAD2,IFD2,IAD1,IFD1,IDATAP,ADRX)
      RETURN
      END