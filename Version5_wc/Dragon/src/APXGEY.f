*DECK APXGEY
      SUBROUTINE APXGEY(IPAPX,IPEDIT,NISO,NG,NMIL,NBISO,NDFI,NISFS,
     1 NISPS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To recover the fission yields of an elementary calculation.
*
*Copyright:
* Copyright (C) 2025 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPAPX   pointer to the Apex file.
* IPEDIT  pointer to the edition object (L_EDIT signature).
* NISO    number of particularized isotopes.
* NG      number of condensed energy groups.
* NMIL    number of mixtures in the MPO file.
* NBISO   number of isotopes in the condensed microlib of the edition
*         object. A given isotope may appear in many mixtures.
* NDFI    number of fissile isotopes producing fission products in
*         the edition object.
* NISFS   number of particularized fissile isotopes.
* NISPS   number of particularized fission products.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE hdf5_wrap
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPAPX,IPEDIT
      INTEGER NISO,NG,NMIL,NBISO,NDFI,NISFS,NISPS
*----
*  LOCAL VARIABLES
*----
      PARAMETER (MAXISO=800)
      TYPE(C_PTR) JPEDIT,KPEDIT
      CHARACTER TEXT8*8,TEXT12*12,RECNAM*80
      LOGICAL LGIMF
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: MIX,PIFI,ADRY
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: ISONAM
      REAL, ALLOCATABLE, DIMENSION(:) :: DEN,PYIELD,SIG,PFIRA
      REAL, ALLOCATABLE, DIMENSION(:,:) :: FLUXES
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: YLDS
      CHARACTER(LEN=4), ALLOCATABLE, DIMENSION(:) :: TYPISO
      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:) :: NOMISO
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) ::  IPISO
*----
*  SCRATCH STORAGE ALLOCATION
*   PFIRA   fission rate.
*   ADRY    offset in YLDS array for fissile isotopes (positive) and
*           fission products (negative).
*----
      ALLOCATE(ISONAM(3,NBISO),MIX(NBISO),PIFI(NDFI))
      ALLOCATE(YLDS(NISFS,NISPS,1),DEN(NBISO),PYIELD(NDFI),
     1 FLUXES(NMIL,NG),SIG(NG),PFIRA(NBISO),ADRY(NISO))
      ALLOCATE(IPISO(NBISO))
*----
*  RECOVER INFORMATION FROM THE /contents/isotopes GROUP.
*----
      IF(NISO.GT.0) THEN
        CALL hdf5_read_data(IPAPX,"/physconst/ISOTA",NOMISO)
        CALL hdf5_read_data(IPAPX,"/physconst/ISOTYP",TYPISO)
      ENDIF
*
      CALL LCMGET(IPEDIT,'ISOTOPESUSED',ISONAM)
      CALL LCMGET(IPEDIT,'ISOTOPESMIX',MIX)
      CALL LCMGET(IPEDIT,'ISOTOPESDENS',DEN)
      CALL LIBIPS(IPEDIT,NBISO,IPISO)
*----
*  COMPUTE ARRAY ADRY.
*----
      ISF=0
      ISP=0
      ADRY(:NISO)=0
      DO 30 ISO=1,NISO
      DO 10 IBISO=1,NBISO
      WRITE(TEXT8,'(2A4)') (ISONAM(I0,IBISO),I0=1,2)
      IF(NOMISO(ISO).EQ.TEXT8) GO TO 20
   10 CONTINUE
      GO TO 30
   20 IF(TYPISO(ISO).EQ.'FISS') THEN
         ISF=ISF+1
         ADRY(ISO)=ISF
      ELSEIF(TYPISO(ISO).EQ.'F.P.') THEN
         ISP=ISP+1
         ADRY(ISO)=-ISP
      ENDIF
   30 CONTINUE
      LGIMF=NISFS.GT.0
      IMF=0
      IF(LGIMF) IMF=ADRY(NISO)
*----
*  RECOVER THE NEUTRON FLUX.
*----
      CALL LCMSIX(IPEDIT,'MACROLIB',1)
      JPEDIT=LCMGID(IPEDIT,'GROUP')
      DO 40 IGR=1,NG
      KPEDIT=LCMGIL(JPEDIT,IGR)
      CALL LCMGET(KPEDIT,'FLUX-INTG',FLUXES(1,IGR))
   40 CONTINUE
      CALL LCMSIX(IPEDIT,' ',2)
*----
*  RECOVER THE FISSION RATES.
*----
      DO 65 IBISO=1,NBISO
      GAR=0.0
      IF(MIX(IBISO).EQ.0) GO TO 60
      KPEDIT=IPISO(IBISO)
      CALL LCMLEN(KPEDIT,'NFTOT',ILONG,ITYLCM)
      IF(ILONG.GT.0) THEN
         CALL LCMGET(KPEDIT,'NFTOT',SIG)
         DO 50 IGR=1,NG
         GAR=GAR+FLUXES(MIX(IBISO),IGR)*DEN(IBISO)*SIG(IGR)
   50    CONTINUE
      ENDIF
   60 PFIRA(IBISO)=GAR
   65 CONTINUE
*----
*  LOOP OVER MPO MIXTURES TO RECOVER THE FISSION YIELDS.
*----
      DO 140 IMIL=1,NMIL
      YLDS(:NISFS,:NISPS,1)=0.0
      DO 130 IBISO=1,NBISO
      IF(MIX(IBISO).EQ.IMIL) THEN
         WRITE(TEXT12,'(3A4)') (ISONAM(I0,IBISO),I0=1,3)
         DO 80 ISO=1,NISO
         IISO=ISO
         IF(NOMISO(ISO).EQ.TEXT12(:8)) GO TO 90
   80    CONTINUE
         GO TO 130
   90    KPEDIT=IPISO(IBISO)
*
*        RECOVER THE FISSION YIELDS.
         CALL LCMLEN(KPEDIT,'PYIELD',ILONG,ITYLCM)
         IF((ILONG.GT.0).AND.(ILONG.EQ.NDFI)) THEN
            CALL LCMGET(KPEDIT,'PIFI',PIFI)
            CALL LCMGET(KPEDIT,'PYIELD',PYIELD)
         ELSE
            GO TO 130
         ENDIF
         IFP=-ADRY(IISO)
         IF(IFP.GT.0) THEN
*           Particular fission product found.
*           If exists in medium, find position in microlib
*           and search all fissiles.
            YLDW=0.0
            DO 120 IDFI=1,NDFI
            JBISO=PIFI(IDFI)
            IF(JBISO.GT.NBISO) CALL XABORT('APXGEY: MIX OVERFLOW.')
            IF(JBISO.EQ.0) GO TO 120
            IF(MIX(JBISO).NE.IMIL) GO TO 120
            WRITE(TEXT8,'(3A4)') (ISONAM(I0,JBISO),I0=1,2)
            DO 100 JSO=1,NISO
            JISO=JSO
            IF(NOMISO(JSO).EQ.TEXT8) GO TO 110
  100       CONTINUE
*           Mother isotope is in residual macro.
            YLDW=YLDW+PFIRA(JBISO)
            IF(IMF.EQ.0) CALL XABORT('APXGEY: LGIMF IS FALSE.')
            YLDS(IMF,IFP,1)=YLDS(IMF,IFP,1)+PYIELD(IDFI)*PFIRA(JBISO)
            GO TO 120
*
*           Yield for selected isotopes.
  110       IFI=ADRY(JISO)
            IF(IFI.LE.0) CALL XABORT('APXGEY: BAD ADRY.')
            YLDS(IFI,IFP,1)=PYIELD(IDFI)
  120       CONTINUE
            IF(LGIMF) THEN
               IF(YLDW.NE.0.0) YLDS(IMF,IFP,1)=YLDS(IMF,IFP,1)/YLDW
            ENDIF
         ENDIF
      ENDIF
  130 CONTINUE
*----
*  STORE INFORMATION IN THE physconst GROUP.
*----
      IF(NMIL.EQ.1) THEN
        CALL hdf5_write_data(IPAPX,"/physconst/FYIELDS",YLDS)
      ELSE
        WRITE(RECNAM,'(18H/physconst/FYIELDS,I8)') IMIL
        CALL hdf5_write_data(IPAPX,TRIM(RECNAM),YLDS)
      ENDIF
  140 CONTINUE
      IF(NISO.GT.0) DEALLOCATE(NOMISO,TYPISO)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(IPISO)
      DEALLOCATE(ADRY)
      DEALLOCATE(PFIRA,SIG,FLUXES,PYIELD,DEN,YLDS)
      DEALLOCATE(PIFI,MIX,ISONAM)
      RETURN
      END
