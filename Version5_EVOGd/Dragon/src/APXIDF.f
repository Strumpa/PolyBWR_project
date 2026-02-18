*DECK APXIDF
      SUBROUTINE APXIDF(IPAPX,IPEDIT,NG,NMIL,ICAL,IDF,NALBP,FNORM,
     1 VOLMIL,FLXMIL)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To store discontinuity factor and albedo information in the Apex file.
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
* NG      number of condensed energy groups.
* NMIL    number of mixtures.
* ICAL    index of the current elementary calculation.
* IDF     type of surfacic information (2/3: boundary flux/DF).
* NALBP   number of physical albedos per energy group.
* FNORM   flux normalization factor.
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
      TYPE(C_PTR) IPAPX,IPEDIT
      INTEGER NG,NMIL,ICAL,IDF,NALBP
      REAL FNORM,VOLMIL(NMIL),FLXMIL(NMIL,NG)
*----
*  LOCAL VARIABLES
*----
      CHARACTER HSMG*131,RECNAM*80,RECNAM2*80
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: SURF
      REAL, ALLOCATABLE, DIMENSION(:,:) :: VREAL,ALBP
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: DISFAC
      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:) :: HADF
*----
*  RECOVER DISCONTINUITY FACTOR INFORMATION FROM MACROLIB
*----
      CALL LCMSIX(IPEDIT,'MACROLIB',1)
      CALL LCMLEN(IPEDIT,'ADF',ILONG,ITYLCM)
      IF(ILONG.NE.0) THEN
        CALL LCMSIX(IPEDIT,'ADF',1)
        CALL LCMGET(IPEDIT,'NTYPE',NSURFD)
        NGG=0
        IF((IDF.EQ.2).OR.(IDF.EQ.3)) THEN
          NGG=NG
        ELSE
          CALL XABORT('APXIDF: INVALID ADF OPTION.')
        ENDIF
        ALLOCATE(DISFAC(NSURFD,NGG,NMIL),SURF(NMIL*NGG),HADF(NSURFD))
        CALL LCMGTC(IPEDIT,'HADF',8,NSURFD,HADF)
        DO I=1,NSURFD
          CALL LCMLEN(IPEDIT,HADF(I),ILONG,ITYLCM)
          IF(IDF.EQ.2) THEN
*           boundary flux information
            IF(ILONG.NE.NMIL*NG) THEN
              WRITE(HSMG,'(16HAPXIDF: INVALID ,A,8H LENGTH=,I5,
     1        10H EXPECTED=,I5,4H.(1))') HADF(I),ILONG,NMIL*NG
              CALL XABORT(HSMG)
            ENDIF
            CALL LCMGET(IPEDIT,HADF(I),SURF)
            DO IMIL=1,NMIL
              DO IGR=1,NG
                IF(FNORM.NE.1.0) THEN
                  DISFAC(I,IGR,IMIL)=SURF((IGR-1)*NMIL+IMIL)*
     1            FNORM*1.0E13*VOLMIL(IMIL)/FLXMIL(IMIL,IGR)
                ELSE
                  DISFAC(I,IGR,IMIL)=SURF((IGR-1)*NMIL+IMIL)*
     1            VOLMIL(IMIL)/FLXMIL(IMIL,IGR)
                ENDIF
              ENDDO
            ENDDO
          ELSE IF(IDF.EQ.3) THEN
*           discontinuity factor information
            IF(ILONG.NE.NMIL*NG) THEN
              WRITE(HSMG,'(16HAPXIDF: INVALID ,A,8H LENGTH=,I5,
     1        10H EXPECTED=,I5,4H.(2))') HADF(I),ILONG,NMIL*NG
              CALL XABORT(HSMG)
            ENDIF
            CALL LCMGET(IPEDIT,HADF(I),SURF)
            DO IMIL=1,NMIL
              DO IGR=1,NG
                IOF=(IGR-1)*NMIL+IMIL
                DISFAC(I,IGR,IMIL)=SURF(IOF)
              ENDDO
            ENDDO
          ENDIF
        ENDDO
        DEALLOCATE(HADF,SURF)
        CALL LCMSIX(IPEDIT,' ',2)
*----
*  MOVE TO THE /calc_id/miscellaneous/ GROUP.
*----
        WRITE(RECNAM,'(4Hcalc,I8,15H/miscellaneous/)') ICAL
        IF((IDF.EQ.2).OR.(IDF.EQ.3)) THEN
          IF(NMIL.EQ.1) THEN
            ALLOCATE(VREAL(NSURFD,NG))
            VREAL(:NSURFD,:NG)=DISFAC(:NSURFD,:NG,1)
            CALL hdf5_write_data(IPAPX,TRIM(RECNAM)//"ADF",VREAL)
            DEALLOCATE(VREAL)
          ELSE
            DO IMIL=1,NMIL
              WRITE(RECNAM2,'(A,3HADF,I8)') TRIM(RECNAM),IMIL
              ALLOCATE(VREAL(NSURFD,NG))
              VREAL(:NSURFD,:NG)=DISFAC(:NSURFD,:NG,IMIL)
              CALL hdf5_write_data(IPAPX,TRIM(RECNAM2),VREAL)
              DEALLOCATE(VREAL)
            ENDDO
          ENDIF
        ENDIF
        DEALLOCATE(DISFAC)
      ENDIF
*----
*  RECOVER AND SAVE ALBEDO INFORMATION
*----
      IF(NALBP.NE.0) THEN
        WRITE(RECNAM,'(4Hcalc,I8,15H/miscellaneous/)') ICAL
        CALL LCMLEN(IPEDIT,'ALBEDO',ILONG,ITYLCM)
        IF(ILONG.EQ.NALBP*NG) THEN
          ALLOCATE(ALBP(NALBP,NG))
          CALL LCMGET(IPEDIT,'ALBEDO',ALBP)
          CALL hdf5_write_data(IPAPX,TRIM(RECNAM)//"ALBEDO",ALBP)
          DEALLOCATE(ALBP)
        ELSE
          CALL XABORT('APXIDF: INCONSISTENT ALBEDO INFORMATION.')
        ENDIF
      ENDIF
      CALL LCMSIX(IPEDIT,' ',2)
      RETURN
      END
