*DECK APXCAL
      SUBROUTINE APXCAL(IMPX,IPAPX,IPDEPL,IPMICR,HEQUI)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Store the results of an elementary calculation in the Apex file
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
* IMPX    print parameter.
* IPAPX   pointer to the Apex file.
* IPDEPL  pointer to the burnup object (L_BURNUP signature).
* IPMICR  pointer to the microlib to include (L_LIBRARY signature).
* HEQUI   keyword of SPH-factor set in the Apex file.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE hdf5_wrap
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPAPX,IPDEPL,IPMICR,IPSPH
      INTEGER IMPX
      CHARACTER(LEN=80) HEQUI
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      INTEGER IPAR(NSTATE)
      REAL BIRRAD(2)
      INTEGER RANK,TYPE,NBYTE,DIMSR(5)
      CHARACTER RECNAM*80,RECNAM2*80,HSMG*131
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: VOLMIL,WORK1
      REAL, ALLOCATABLE, DIMENSION(:,:) :: FLXMIL,RVAL0
*
      CALL LCMLEN(IPMICR,'STATE-VECTOR',ILONG,ITYLCM)
      IF(ILONG.NE.0) THEN
         CALL LCMGET(IPMICR,'STATE-VECTOR',IPAR)
         NBISO=IPAR(2)
         NED=IPAR(13)
         NPRC=IPAR(19)
         NDFI=IPAR(20)
      ELSE
         NBISO=0
         NDFI=0
      ENDIF
      CALL LCMSIX(IPMICR,'MACROLIB',1)
      CALL LCMGET(IPMICR,'STATE-VECTOR',IPAR)
      NG=IPAR(1)
      NMIL=IPAR(2)
      NL=IPAR(3)
      IF(IPAR(4).GT.1) CALL XABORT('APXCAL: CANNOT PROCESS MULTIPLE FI'
     1 //'SSION SPECTRA.')
      NED=IPAR(5)
      ITRANC=IPAR(6)
      NPRC=IPAR(7)
      NALBP=IPAR(8)
      IDF=IPAR(12)
      CALL LCMLEN(IPMICR,'SPH',ILEN,ITYLCM)
      IF(ILEN.NE.0) THEN
         IPSPH=LCMGID(IPMICR,'SPH')
         CALL LCMGET(IPSPH,'STATE-VECTOR',IPAR)
         IMC=IPAR(6)
      ELSE
         IMC=0
      ENDIF
      CALL hdf5_info(IPAPX,"/NCALS",RANK,TYPE,NBYTE,DIMSR)
      IF(RANK.EQ.99) THEN
         NCALS=0
      ELSE
         CALL hdf5_read_data(IPAPX,"/NCALS",NCALS)
      ENDIF
      ICAL=NCALS+1
      CALL hdf5_write_data(IPAPX,"/NCALS",ICAL)
      CALL LCMSIX(IPMICR,' ',2)
      WRITE(RECNAM,'(4Hcalc,I8,1H/)') ICAL
      IF(IMPX.GT.0) WRITE(6,'(/19H APXCAL: NEW GROUP ,A)') TRIM(RECNAM)
      CALL hdf5_create_group(IPAPX,RECNAM)
      CALL hdf5_create_group(IPAPX,TRIM(RECNAM)//"miscellaneous/")
*----
*  RECOVER THE FLUX NORMALIZATION FACTOR.
*----
      IF(C_ASSOCIATED(IPDEPL)) THEN
         CALL LCMGET(IPDEPL,'BURNUP-IRRAD',BIRRAD)
         BURN=BIRRAD(1)
         CALL LCMLEN(IPDEPL,'FLUX-NORM',ILONG,ITYLCM)
         IF(ILONG.EQ.0) THEN
            WRITE(HSMG,'(40HAPXCAL: THE ''FLUX-NORM'' RECORD IS NOT SE,
     1      20HT FOR BURNUP STEP AT,E12.5,14H MW-DAY/TONNE.)') BURN
            CALL XABORT(HSMG)
         ENDIF
         CALL LCMGET(IPDEPL,'FLUX-NORM',FNORM)
         IF(IMPX.GT.0) WRITE(6,100) FNORM,BURN
      ELSE
         FNORM=1.0
         IF(IMPX.GT.0) WRITE(6,110)
      ENDIF
*----
*  RECOVER THE CROSS SECTIONS.
*----
      NISO=0
      CALL hdf5_info(IPAPX,"/explicit/ISONAME",RANK,TYPE,NBYTE,DIMSR)
      IF(RANK.NE.99) NISO=DIMSR(1)
      NMAC=0
      CALL hdf5_info(IPAPX,"/explicit/MACNAME",RANK,TYPE,NBYTE,DIMSR)
      IF(RANK.NE.99) NMAC=DIMSR(1)
      NREA=0
      CALL hdf5_info(IPAPX,"/explicit/REANAME",RANK,TYPE,NBYTE,DIMSR)
      IF(RANK.NE.99) NREA=DIMSR(1)
      ALLOCATE(VOLMIL(NMIL),FLXMIL(NMIL,NG))
      CALL APXCA2(IPAPX,IPMICR,NREA,NISO,NMAC,NED,NPRC,NG,NL,ITRANC,
     1 NALBP,IMC,NMIL,NBISO,ICAL,IMPX,FNORM,NMILNR,NISFS,NISPS,VOLMIL,
     2 FLXMIL)
*----
*  RECOVER DISCONTINUITY FACTOR INFORMATION.
*----
      IF((IDF.EQ.2).OR.(IDF.EQ.3).OR.(NALBP.GT.0)) THEN
        CALL APXIDF(IPAPX,IPMICR,NG,NMIL,ICAL,IDF,NALBP,FNORM,VOLMIL,
     1  FLXMIL)
      ENDIF
*----
*  RECOVER THE FISSION YIELDS.
*----
      IF((ICAL.EQ.1).AND.(NISFS*NISPS.GT.0)) THEN
        CALL APXGEY(IPAPX,IPMICR,NISO,NG,NMIL,NBISO,NDFI,NISFS,NISPS)
      ENDIF
*----
*  RECOVER SPH FACTOR INFORMATION.
*----
      CALL LCMSIX(IPMICR,'MACROLIB',1)
      CALL LCMLEN(IPMICR,'SPH',ILEN,ITYLCM)
      CALL LCMSIX(IPMICR,' ',2)
      IF(ILEN.NE.0) THEN
        IF(HEQUI.EQ.' ') HEQUI='default'
        ALLOCATE(WORK1(NG),RVAL0(NG,NMIL))
        CALL SAPSPH(IPMICR,NG,NMIL,1,NG,RVAL0)
        IF(NMIL.EQ.1) THEN
          WORK1(:NG)=RVAL0(:NG,1)
          WRITE(RECNAM,'(4Hcalc,I8,14H/xs/MEDIA_SPH/)') ICAL
          CALL hdf5_create_group(IPAPX,TRIM(RECNAM))
          WRITE(RECNAM2,'(A,A)') TRIM(RECNAM),TRIM(HEQUI)
          CALL hdf5_write_data(IPAPX,TRIM(RECNAM2),WORK1)
        ELSE
          DO IBM=1,NMIL
            WORK1(:NG)=RVAL0(:NG,IBM)
            WRITE(RECNAM,'(4Hcalc,I8,3H/xs,I8,11H/MEDIA_SPH/)') ICAL,IBM
            CALL hdf5_create_group(IPAPX,TRIM(RECNAM))
            WRITE(RECNAM2,'(A,A)') TRIM(RECNAM),TRIM(HEQUI)
            CALL hdf5_write_data(IPAPX,TRIM(RECNAM2),WORK1)
          ENDDO
        ENDIF
        DEALLOCATE(RVAL0,WORK1)
      ENDIF
      DEALLOCATE(FLXMIL,VOLMIL)
      RETURN
*
  100 FORMAT(45H APXCAL: NORMALIZE THE FLUX WITH THE FACTOR =,1P,E12.5,
     1 26H TAKEN FROM BURNUP STEP AT,E12.5,14H MW-DAY/TONNE.)
  110 FORMAT(36H APXCAL: THE FLUX IS NOT NORMALIZED.)
      END
