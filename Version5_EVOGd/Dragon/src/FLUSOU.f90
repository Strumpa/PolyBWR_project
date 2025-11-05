SUBROUTINE FLUSOU(CDOOR,HLEAK,MAX1,IG,IPTRK,KPMACR,NMAT,NANIS,NUN,NGRP, &
  & FUNKNO,SUNKNO)
  !
  !---------------------------------------------------------------------
  !
  !Purpose:
  ! compute the out-of-group scattering source in general cases.
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
  ! HLEAK   type of model (=' ': general; ='ECCO'; ='TIBERE').
  ! MAX1    first dimension of FUNKNO and SOURCE arrays.
  ! IG      secondary group.
  ! IPTRK   pointer to the tracking (L_TRACK signature).
  ! KPMACR  pointer to the secondary-group related macrolib information.
  ! NMAT    number of mixtures in the macrolib.
  ! NANIS   number of Legendre components in the macrolib.
  ! NUN     total number of flux or source unknowns.
  ! NGRP    number of energy groups.
  ! FUNKNO  unknown vector.
  !
  !Parameters: output
  ! SUNKNO  source vector.
  !---------------------------------------------------------------------
  !
  USE GANLIB
  !----
  !  SUBROUTINE ARGUMENTS
  !----
  CHARACTER(LEN=12), INTENT(IN) :: CDOOR
  CHARACTER(LEN=6), INTENT(IN) :: HLEAK
  TYPE(C_PTR), INTENT(IN) :: IPTRK,KPMACR
  INTEGER, INTENT(IN) :: MAX1,IG,NMAT,NANIS,NUN,NGRP
  REAL, DIMENSION(MAX1,NGRP), INTENT(IN) :: FUNKNO
  REAL, DIMENSION(MAX1,NGRP), INTENT(INOUT) :: SUNKNO
  !----
  !  LOCAL VARIABLES
  !----
  INTEGER, PARAMETER :: NSTATE=40
  INTEGER, DIMENSION(NSTATE) :: ISTATE
  INTEGER, DIMENSION(3) :: INDD
  CHARACTER CAN(0:19)*2
  !----
  !  ALLOCATABLE ARRAYS
  !----
  INTEGER, ALLOCATABLE, DIMENSION(:) :: MATCOD
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) :: KEYFLX
  INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJ,NJJ,IPOS
  REAL, ALLOCATABLE, DIMENSION(:) :: XSCAT
  !----
  !  DATA STATEMENTS
  !----
  DATA CAN /'00','01','02','03','04','05','06','07','08','09', &
          & '10','11','12','13','14','15','16','17','18','19'/
  !----
  !  SCRATCH STORAGE ALLOCATION
  !----
  ALLOCATE(IJJ(0:NMAT),NJJ(0:NMAT),IPOS(0:NMAT))
  ALLOCATE(XSCAT(0:NMAT*NGRP))
  !----
  !  RECOVER TRACKING PARAMETERS
  !  NFUNL: number of spherical harmonics components used to expand the
  !     flux and the sources.
  !  NANIS_TRK: number of components in the angular expansion of the flux
  !----
  CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
  NREG=ISTATE(1)
  IF(ISTATE(2).GT.NUN) CALL XABORT('FLUSOU: WRONG NUN.')
  IF(ISTATE(4).GT.NMAT) CALL XABORT('FLUSOU: WRONG NMAT.')
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
  ENDIF
  ALLOCATE(MATCOD(NREG),KEYFLX(NREG,NLIN,NFUNL))
  KEYFLX(:NREG,:NLIN,:NFUNL)=0
  CALL LCMLEN(IPTRK,'MATCOD',ILNLCM,ITYLCM)
  IF(ILNLCM.NE.NREG) CALL XABORT('FLUSOU: INCOMPATIBLE NUMBER OF REGIONS.')
  CALL LCMGET(IPTRK,'MATCOD',MATCOD)
  IF(CDOOR.EQ.'MCCG') THEN
    CALL LCMGET(IPTRK,'KEYFLX$ANIS',KEYFLX)
  ELSE
    CALL LCMGET(IPTRK,'KEYFLX',KEYFLX)
  ENDIF
  !----
  !  COMPUTE THE SCATTERING SOURCE IN THE GENERAL CASE
  !----
  IF(HLEAK.EQ.' ') THEN
    NIAL=MIN(NFUNL-1,NANIS,NANIS_TRK-1)
    DO IAL=0,NIAL
      FACT=REAL(2*IAL+1)
      CALL LCMGET(KPMACR,'NJJS'//CAN(IAL),NJJ(1))
      CALL LCMGET(KPMACR,'IJJS'//CAN(IAL),IJJ(1))
      CALL LCMGET(KPMACR,'IPOS'//CAN(IAL),IPOS(1))
      CALL LCMGET(KPMACR,'SCAT'//CAN(IAL),XSCAT(1))
      DO IR=1,NREG
        IBM=MATCOD(IR)
        IF(IBM.LE.0) CYCLE
        DO IAM=0,NIAL
          DO IE=1,NLIN
            IND=0
            IF(NDIM.EQ.3) THEN
              IF(1+IAL*NANIS_TRK+IAM.GT.NFUNL) CALL XABORT('FLUSOU: KEYFLX OVERFLOW(1)')
              IND=KEYFLX(IR,IE,1+IAL*NANIS_TRK+IAM)
            ELSE IF((NDIM.EQ.2).AND.(IAM.LE.IAL)) THEN
              IF(1+IAL*(IAL+1)/2+IAM.GT.NFUNL) CALL XABORT('FLUSOU: KEYFLX OVERFLOW(2)')
              IND=KEYFLX(IR,IE,1+IAL*(IAL+1)/2+IAM)
            ELSE IF(IAM.EQ.IAL) THEN
              IND=KEYFLX(IR,IE,1+IAL)
            ENDIF
            IF(IND.EQ.0) THEN
              CYCLE
            ELSE IF(IND.GT.NUN) THEN
              CALL XABORT('FLUSOU: NUN OVERFLOW.')
            ENDIF
            JG=IJJ(IBM)
            DO JND=1,NJJ(IBM)
              IF(JG.NE.IG) THEN
                SUNKNO(IND,IG)=SUNKNO(IND,IG)+FACT*XSCAT(IPOS(IBM)+JND-1)* &
                & FUNKNO(IND,JG)
              ENDIF
              JG=JG-1
            ENDDO ! JND
          ENDDO ! IE
        ENDDO ! IAM
      ENDDO ! IR
    ENDDO
  !----
  !  COMPUTE THE SCATTERING SOURCE WITH ECCO MODEL
  !----
  ELSE IF(HLEAK.EQ.'ECCO') THEN
    CALL LCMGET(KPMACR,'NJJS01',NJJ(1))
    CALL LCMGET(KPMACR,'IJJS01',IJJ(1))
    CALL LCMGET(KPMACR,'IPOS01',IPOS(1))
    CALL LCMGET(KPMACR,'SCAT01',XSCAT(1))
    DO IR=1,NREG
      IBM=MATCOD(IR)
      IF(IBM.LE.0) CYCLE
      DO IE=1,NLIN
        IND=MAX1/2+KEYFLX(IR,IE,1)
        JG=IJJ(IBM)
        DO JND=1,NJJ(IBM)
          IF(JG.NE.IG) THEN
            SUNKNO(IND,IG)=SUNKNO(IND,IG)+XSCAT(IPOS(IBM)+JND-1)* &
            & FUNKNO(IND,JG)
          ENDIF
          JG=JG-1
        ENDDO ! JND
      ENDDO ! IE
    ENDDO ! IR
  !----
  !  COMPUTE THE SCATTERING SOURCE WITH TIBERE MODEL
  !----
  ELSE IF(HLEAK.EQ.'TIBERE') THEN
    CALL LCMGET(KPMACR,'NJJS01',NJJ(1))
    CALL LCMGET(KPMACR,'IJJS01',IJJ(1))
    CALL LCMGET(KPMACR,'IPOS01',IPOS(1))
    CALL LCMGET(KPMACR,'SCAT01',XSCAT(1))
    DO IR=1,NREG
      IBM=MATCOD(IR)
      IF(IBM.LE.0) CYCLE
      DO IE=1,NLIN
        INDD(1)=MAX1/4+KEYFLX(IR,IE,1)
        INDD(2)=MAX1/2+KEYFLX(IR,IE,1)
        INDD(3)=3*MAX1/4+KEYFLX(IR,IE,1)
        DO IDIR=1,3
          IND=INDD(IDIR)
          JG=IJJ(IBM)
          DO JND=1,NJJ(IBM)
            IF(JG.NE.IG) THEN
              SUNKNO(IND,IG)=SUNKNO(IND,IG)+XSCAT(IPOS(IBM)+JND-1)* &
              & FUNKNO(IND,JG)
            ENDIF
            JG=JG-1
          ENDDO ! IND
        ENDDO ! IDIR
      ENDDO ! IE
    ENDDO ! IR
  ENDIF
  !----
  !  SCRATCH STORAGE DEALLOCATION
  !----
  DEALLOCATE(KEYFLX,MATCOD)
  DEALLOCATE(XSCAT)
  DEALLOCATE(IPOS,NJJ,IJJ)
  RETURN
END SUBROUTINE FLUSOU
