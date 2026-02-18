MODULE DOORS_MOD
  USE GANLIB
CONTAINS
  SUBROUTINE DOORS(CDOOR,IPTRK,NMAT,NANIS,NUN,SIGG,SUNKNO,FUNKNO)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! compute the product of a cross section times a flux unknow vector.
    !
    !Copyright:
    ! Copyright (C) 2025 Ecole Polytechnique de Montreal
    ! This library is free software; you can redistribute it and/or
    ! modify it under the terms of the GNU Lesser General Public
    ! License as published by the Free Software Foundation; either
    ! version 2.1 of the License, or (at your option) any later version
    !
    !Author(s): A. Hebert
    !
    !Parameters: input
    ! CDOOR   name of the geometry/solution operator.
    ! IPTRK   pointer to the tracking (L_TRACK signature).
    ! NMAT    number of mixtures in the macrolib.
    ! NANIS   number of Legendre components in the macrolib (=0: isotropic).
    ! NUN     total number of unknowns in vectors SUNKNO and FUNKNO.
    ! SIGG    cross section.
    ! FUNKNO  optional unknown vector. If not present, a flat flux
    !         approximation is assumed.
    !
    !Parameters: output
    ! SUNKNO  source vector. Volumes are included with BIVAC and TRIVAC
    !         trackings.
    !
    !---------------------------------------------------------------------
    !
    !----
    !  SUBROUTINE ARGUMENTS
    !----
    CHARACTER(LEN=12), INTENT(IN) :: CDOOR
    TYPE(C_PTR), INTENT(IN) :: IPTRK
    INTEGER, INTENT(IN) :: NMAT,NANIS,NUN
    REAL, DIMENSION(0:NMAT,NANIS+1), INTENT(IN) :: SIGG
    REAL, DIMENSION(NUN), INTENT(IN), OPTIONAL :: FUNKNO
    REAL, DIMENSION(NUN), INTENT(INOUT) :: SUNKNO
    !----
    !  LOCAL VARIABLES
    !----
    INTEGER, PARAMETER :: NSTATE=40
    INTEGER, DIMENSION(NSTATE) :: ISTATE
    !----
    !  ALLOCATABLE ARRAYS
    !----
    INTEGER, ALLOCATABLE, DIMENSION(:) :: MATCOD
    INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: KEYFLX
    REAL, ALLOCATABLE, DIMENSION(:) :: VOL
    !----
    !  RECOVER TRACKING PARAMETERS
    !  NFUNL: number of spherical harmonics components used to expand the
    !         flux and the sources.
    !  NANIS_TRK: number of components in the angular expansion of the flux
    !----
    CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
    NREG=ISTATE(1)
    IF(ISTATE(2).GT.NUN) CALL XABORT('DOORS: WRONG NUN.')
    IF(ISTATE(4).GT.NMAT) CALL XABORT('DOORS: WRONG NMAT.')
    NDIM=0
    NLIN=1
    NFUNL=1
    NANIS_TRK=1
    IF(CDOOR.EQ.'MCCG') THEN
      NANIS_TRK=ISTATE(6)
      NDIM=ISTATE(16)
      CALL LCMGET(IPTRK,'MCCG-STATE',ISTATE)
      NFUNL=ISTATE(19)
      NLIN=ISTATE(20)
    ELSE IF(CDOOR.EQ.'SN') THEN
      NFUNL=ISTATE(7)
      NLIN=ISTATE(8)
      NDIM=ISTATE(9)
      NLIN=NLIN**NDIM
      NLIN=NLIN*ISTATE(35)
      NANIS_TRK=ISTATE(16)
    ELSE IF(CDOOR.EQ.'BIVAC') THEN
      NLIN=ABS(ISTATE(8)) ! order of finite elements
      NFUNL=MAX(1,ISTATE(14))
      NANIS_TRK=ABS(ISTATE(16))
    ELSE IF(CDOOR.EQ.'TRIVAC') THEN
      NLIN=ABS(ISTATE(9)) ! order of finite elements
      NLFUNL=MAX(1,ISTATE(30))
      NANIS_TRK=ABS(ISTATE(32))
    ENDIF
    ALLOCATE(MATCOD(NREG),VOL(NREG),KEYFLX(NREG,NLIN,NFUNL))
    KEYFLX(:NREG,:NLIN,:NFUNL)=0
    CALL LCMLEN(IPTRK,'MATCOD',ILNLCM,ITYLCM)
    IF(ILNLCM.NE.NREG) CALL XABORT('DOORS: INCOMPATIBLE NUMBER OF REGIONS.')
    CALL LCMGET(IPTRK,'MATCOD',MATCOD)
    CALL LCMGET(IPTRK,'VOLUME',VOL)
    IF((CDOOR.EQ.'MCCG').OR.(CDOOR.EQ.'SN')) THEN
      CALL LCMGET(IPTRK,'KEYFLX$ANIS',KEYFLX)
    ELSE
      CALL LCMGET(IPTRK,'KEYFLX',KEYFLX)
    ENDIF
    !----
    !  PERFORM SIGG*FUNKNO MULTIPLICATION
    !----
    IF(CDOOR.EQ.'SN') THEN
      CALL DOORS_SN(IPTRK,NANIS,NREG,NMAT,NUN,MATCOD,SIGG,SUNKNO,FUNKNO)
    ELSE IF(CDOOR.EQ.'BIVAC') THEN
      CALL DOORS_BIV(IPTRK,NANIS,NREG,NMAT,NUN,MATCOD,VOL,SIGG,SUNKNO,FUNKNO)
    ELSE IF(CDOOR.EQ.'TRIVAC') THEN
      CALL DOORS_TRI(IPTRK,NANIS,NREG,NMAT,NUN,MATCOD,VOL,SIGG,SUNKNO,FUNKNO)
    ELSE IF(PRESENT(FUNKNO)) THEN
      ! general case
      IF((NANIS.EQ.0).OR.(NFUNL.EQ.1).OR.(NANIS_TRK.EQ.1)) THEN
        ! LAB isotropy or transport correction
        DO IR=1,NREG
          IBM=MATCOD(IR)
          IF(IBM.LE.0) CYCLE
          DO IE=1,NLIN
            IND=KEYFLX(IR,IE,1)
            SUNKNO(IND)=SUNKNO(IND)+SIGG(IBM,1)*FUNKNO(IND)
          ENDDO
        ENDDO ! IR
      ELSE
        ! spherical harmonics expansion of the flux and source
        DO IR=1,NREG
          IBM=MATCOD(IR)
          IF(IBM.LE.0) CYCLE
          DO IAL=0,MIN(NFUNL-1,NANIS,NANIS_TRK-1)
            FACT=REAL(2*IAL+1)
            DO IAM=0,MIN(NFUNL-1,NANIS,NANIS_TRK-1)
              DO IE=1,NLIN
                IND=0
                IF(NDIM.EQ.3) THEN
                  IND=KEYFLX(IR,IE,1+IAL*NANIS_TRK+IAM)
                ELSE IF((NDIM.EQ.2).AND.(IAM.LE.IAL)) THEN
                  IND=KEYFLX(IR,IE,1+IAL*(IAL+1)/2+IAM)
                ELSE IF(IAM.EQ.IAL) THEN
                  IND=KEYFLX(IR,IE,1+IAL)
                ENDIF
                IF(IND.EQ.0) THEN
                  CYCLE
                ELSE IF(IND.GT.NUN) THEN
                  CALL XABORT('DOORS: NUN OVERFLOW.')
                ENDIF
                SUNKNO(IND)=SUNKNO(IND)+FACT*SIGG(IBM,IAL+1)*FUNKNO(IND)
              ENDDO ! IE
            ENDDO ! IAM
          ENDDO ! IAL
        ENDDO ! IR
      ENDIF
    ELSE
      ! general case (flat flux)
      DO IR=1,NREG
        IND=KEYFLX(IR,1,1)
        SUNKNO(IND)=SUNKNO(IND)+SIGG(MATCOD(IR),1)
      ENDDO
    ENDIF
    DEALLOCATE(KEYFLX,VOL,MATCOD)
  END SUBROUTINE DOORS
END MODULE DOORS_MOD

