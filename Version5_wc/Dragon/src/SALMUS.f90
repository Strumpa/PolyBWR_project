!---------------------------------------------------------------------
!
!Purpose:
! To analyze and track a geometry data structure using the Sanchez
! algorithm for the multicell surfacic approximation.
!
!Copyright:
! Copyright (C) 2025 Ecole Polytechnique de Montreal.
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version
!
!Author(s):
! A. Hebert
!
!Parameters: input
! FGEO    unit file number of the surfacic file in read only mode.
! IPTRK   pointer to the TRACKING data structure in creation mode.
! IFTRK   pointer to the TRACKING file in creation mode.
! RCUTOF  minimum distance between two surfacic elements.
! IPRINT  print level.
! NBSLIN  maximum number of segments in a single tracking line.
!         computed by default in SALTCG but limited to 100000
!         elements. This default value can be bypassed using
!         keyword NBSLIN.
!
!Parameters: output
! GG      geometry basic information.
!
!-----------------------------------------------------------------------
!
SUBROUTINE SALMUS(FGEO ,IPTRK, IFTRK, RCUTOF, IPRINT, NBSLIN, GG)
  USE GANLIB
  USE PRECISION_AND_KINDS, ONLY : PDB
  USE SAL_GEOMETRY_TYPES,  ONLY : T_G_BASIC,LMERGM,TYPGEO,NBFOLD,F_GEO
  USE SAL_TRACKING_TYPES,  ONLY : PRTIND,EPS1
  USE SAL_GEOMETRY_MOD,    ONLY : SAL100
  IMPLICIT     NONE
  !----
  !  Subroutine arguments
  !----
  TYPE(C_PTR) IPTRK
  INTEGER  FGEO,IFTRK,IPRINT,NBSLIN
  REAL(PDB) RCUTOF
  TYPE(T_G_BASIC) :: GG
  !----
  !  Local parameters
  !----
  INTEGER      IOUT
  PARAMETER   (IOUT=6)
  INTEGER      NSTATE
  PARAMETER   (NSTATE=40)
  !----
  !  Local variables
  !----
  TYPE(C_PTR)  JPTRK,KPTRK
  INTEGER      ISTATT(NSTATE),I,J,I1,I2,J1,J2,IMACRO,J3,NMIX,IMIX,NREG,OK, &
               & NINST,ILONG,ITYLCM
  REAL         RSTATT(NSTATE)
  LOGICAL      LGINF
  CHARACTER(LEN=131) HSMG
  !----
  !  Allocatable arrays
  !----
  INTEGER, ALLOCATABLE, DIMENSION(:) :: NSURF_MACRO,NBNODE_MACRO,ITAB, &
  & NMC_NODE,NMC_SURF,IFR,MIX,ICNT,MERGE_MACRO,SURF_MACRO,PERIM,IDMA,IMAC
  REAL, ALLOCATABLE, DIMENSION(:) :: VOLUME,ALB,DVX,GALBED
  !----
  !  Read the surfacic file and fill GG object
  !----
  PRTIND=IPRINT
  F_GEO=FGEO
  EPS1=1.E-5_PDB
  IF(RCUTOF>0._PDB) THEN
    EPS1=RCUTOF
    IF(PRTIND>0) WRITE(*,*) "SALMUS: set eps1 to ",EPS1
  ENDIF
  CALL SAL100(GG)
  !------------
  IF(IPRINT > 0) WRITE(IOUT,'(/" SALMUS: TYPGEO=",I5," NBFOLD=",I5," NMACRO=",I5)') &
  & TYPGEO,NBFOLD,GG%NB_MACRO
  !----
  ! Perform optional MERGE MIX
  !----
  IF(LMERGM) THEN
    GG%NUM_MERGE(:)=GG%MED(:)
    DO I=1,MAXVAL(GG%NUM_MERGE(:))
      1000 IF(I.GT.MAXVAL(GG%NUM_MERGE(:))) EXIT
      DO J=1,GG%NB_NODE
        IF(GG%NUM_MERGE(J).EQ.I) GO TO 2000
      ENDDO
      DO J=1,GG%NB_NODE
        IF(GG%NUM_MERGE(J).GE.I) GG%NUM_MERGE(J)=GG%NUM_MERGE(J)-1
      ENDDO
      GO TO 1000
      2000 CONTINUE
    ENDDO
  ENDIF
  !----
  ! store the STATE VECTOR for the global geometry
  !----
  CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATT)
  CALL LCMGET(IPTRK,'EXCELTRACKOP',RSTATT)
  NREG=MAXVAL(GG%NUM_MERGE)
  ISTATT(1) = NREG ! number of regions
  ISTATT(2) = NREG ! number of flux unknowns
  ISTATT(7) = 5 ! set the multicell surfacic approximation
  LGINF = .TRUE.
  DO I=1,GG%NBBCDA
    LGINF = LGINF .AND. (GG%BCDATAREAD(I)%BCDATA(6) == 1._PDB)
  ENDDO
  ISTATT(3)=1
  IF(.NOT.LGINF) ISTATT(3)=0 ! reset the leakage flag
  ISTATT(4) = MAXVAL(GG%MED(1:GG%NB_NODE)) ! maximum number of mixture
  CALL LCMPUT(IPTRK,'STATE-VECTOR',NSTATE,1,ISTATT)
  !
  ! fill-in medium number per region
  ALLOCATE(ITAB(NREG),VOLUME(NREG), STAT =OK)
  IF(OK /= 0) CALL XABORT('SALMUS: failure to allocate integer ITAB')
  ! fill in MATCOD
  DO J=1,GG%NB_NODE
    ITAB(GG%NUM_MERGE(J)) = GG%MED(J)
  ENDDO
  CALL LCMPUT(IPTRK,'MATCOD',NREG,1,ITAB(1:NREG) )
  ! fill-in KEYFLX per region
  ITAB(:NREG) = (/ (I,I=1,NREG) /)
  CALL LCMPUT(IPTRK,'MERGE',NREG,1,ITAB)
  CALL LCMPUT(IPTRK,'KEYFLX',NREG,1,ITAB)
  ! fill-in volumes per region
  VOLUME(:NREG) =0.
  DO I=1,GG%NB_NODE
    VOLUME(GG%NUM_MERGE(I)) = VOLUME(GG%NUM_MERGE(I)) + REAL(GG%VOL_NODE(I))
  ENDDO
  CALL LCMPUT(IPTRK,'VOLUME',NREG,2,VOLUME)
  IF(IPRINT .GT. 5) THEN
    I1=1
    DO I=1,(NREG-1)/8+1
      I2=I1+7
      IF(I2.GT.NREG) I2=NREG
      WRITE (IOUT,20) (J,J=I1,I2)
      DO J=1,GG%NB_NODE
        ITAB(GG%NUM_MERGE(J)) = GG%MED(J)
      ENDDO
      WRITE (IOUT,30) (ITAB(J),J=I1,I2)
      WRITE (IOUT,40) (VOLUME(J),J=I1,I2)
      I1=I1+8
    ENDDO
  ENDIF
  DEALLOCATE(VOLUME,ITAB)
  !----
  ! Extract the surfacic elements belonging to each macro geometry and
  ! perform tracking
  !----
  ALLOCATE(NSURF_MACRO(GG%NB_MACRO),NBNODE_MACRO(GG%NB_MACRO))
  ALLOCATE(NMC_NODE(GG%NB_MACRO+1),NMC_SURF(GG%NB_MACRO+1))
  NSURF_MACRO(:GG%NB_MACRO) = 0
  NBNODE_MACRO(:GG%NB_MACRO) = 0
  JPTRK=LCMLID(IPTRK,'MACRO-TRACK',GG%NB_MACRO)
  NMC_NODE(1)=0
  NMC_SURF(1)=0
  DO IMACRO=1,GG%NB_MACRO
    KPTRK=LCMDIL(JPTRK,IMACRO)
    CALL LCMPUT(KPTRK,'STATE-VECTOR',NSTATE,1,ISTATT)
    CALL LCMPUT(KPTRK,'EXCELTRACKOP',NSTATE,2,RSTATT)
    CALL MUSACG(KPTRK,IFTRK,IPRINT,IMACRO,NBSLIN,RCUTOF,GG,LGINF, &
    & NBNODE_MACRO(IMACRO),NSURF_MACRO(IMACRO))
    NMC_NODE(IMACRO+1)=NMC_NODE(IMACRO)+NBNODE_MACRO(IMACRO)
    NMC_SURF(IMACRO+1)=NMC_SURF(IMACRO)+NSURF_MACRO(IMACRO)
  ENDDO
  CALL LCMPUT(IPTRK,'NMC_NODE',GG%NB_MACRO+1,1,NMC_NODE)
  CALL LCMPUT(IPTRK,'NMC_SURF',GG%NB_MACRO+1,1,NMC_SURF)
  !----
  !  Create connectivity data
  !----
  NMIX=NMC_SURF(GG%NB_MACRO+1)
  ALLOCATE(IFR(NMIX),ALB(NMIX),MIX(NMIX),DVX(NMIX),ICNT(NMIX),IDMA(NMIX))
  ALLOCATE(IMAC(NREG))
  J1=0
  NMIX=0
  IMIX=0
  DO IMACRO=1,GG%NB_MACRO
    KPTRK=LCMDIL(JPTRK,IMACRO)
    J2=NBNODE_MACRO(IMACRO)
    J3=NSURF_MACRO(IMACRO)
    ALLOCATE(MERGE_MACRO(J2),SURF_MACRO(J3))
    CALL LCMSIX(KPTRK,'SURFACIC_TMP',1)
    CALL LCMGET(KPTRK,'MERGE_MACRO',MERGE_MACRO)
    CALL LCMGET(KPTRK,'SURF_MACRO',SURF_MACRO)
    CALL LCMSIX(KPTRK,'          ',2)
    IF(J1+J2.GT.NREG) CALL XABORT('SALMUS: NREG overflow')
    IMAC(J1+1:J1+J2)=MERGE_MACRO(:J2)
    ICNT(NMIX+1:NMIX+J3)=SURF_MACRO(:J3)
    DEALLOCATE(SURF_MACRO,MERGE_MACRO)
    CALL LCMLEN(KPTRK,'ALBEDO',ILONG,ITYLCM)
    ALLOCATE(GALBED(ILONG),PERIM(J3))
    CALL LCMGET(KPTRK,'ALBEDO',GALBED)
    CALL LCMGET(KPTRK,'PERIM_SURF',PERIM)
    DO I=1,J3
      IF(PERIM(I).GT.ILONG) THEN
        WRITE(HSMG,'(51H SALMUS: inconsistent albedo info in macro geometry,I5)') IMACRO
        CALL XABORT(HSMG)
      ENDIF
      ALB(NMIX+I)=GALBED(PERIM(I))
    ENDDO
    DEALLOCATE(PERIM,GALBED)
    OUT1: DO I=NMIX+1,NMIX+J3
      DO J=NMIX+1,I-1
        IF(ICNT(I).EQ.ICNT(J)) THEN
          MIX(I)=MIX(J)
          CYCLE OUT1
        ENDIF
      ENDDO
      IMIX=IMIX+1
      MIX(I)=IMIX
    ENDDO OUT1
    DO I=NMIX+1,NMIX+J3
      NINST=COUNT(MIX(NMIX+1:NMIX+J3) == MIX(I))
      DVX(I)=1.0/REAL(NINST)
      IDMA(I)=IMACRO
    ENDDO
    J1=J1+J2
    NMIX=NMIX+J3
  ENDDO
  CALL LCMPUT(IPTRK,'MERGE_MACRO',NREG,1,IMAC)
  DEALLOCATE(IMAC)
  OUT2: DO I=1,NMIX
    DO J=1,NMIX
      IF(IDMA(J).EQ.IDMA(I)) CYCLE
      IF(ICNT(J).EQ.ICNT(I)) THEN
        IFR(I)=MIX(J)
        CYCLE OUT2
      ENDIF
    ENDDO
    IFR(I)=MIX(I)
  ENDDO OUT2
  DEALLOCATE(IDMA)
  ISTATT(24)=GG%NB_MACRO
  ISTATT(28)=MAXVAL(IFR) ! number of current unknowns
  ISTATT(29)=NMIX ! number of perimeter elements
  CALL LCMPUT(IPTRK,'STATE-VECTOR',NSTATE,1,ISTATT)
  CALL LCMPUT(IPTRK,'IFR',NMIX,1,IFR)
  CALL LCMPUT(IPTRK,'ALB',NMIX,2,ALB)
  CALL LCMPUT(IPTRK,'MIX',NMIX,1,MIX)
  CALL LCMPUT(IPTRK,'DVX',NMIX,2,DVX)
  IF(IPRINT .GT. 1) THEN
    NMIX=0
    DO IMACRO=1,GG%NB_MACRO
      J3=NSURF_MACRO(IMACRO)
      WRITE (IOUT,50) IMACRO,(ICNT(I),I=NMIX+1,NMIX+J3)
      WRITE (IOUT,60) (MIX(I),IFR(I),I=NMIX+1,NMIX+J3)
      WRITE (IOUT,70) (ALB(I),I=NMIX+1,NMIX+J3)
      WRITE (IOUT,80) (DVX(I),I=NMIX+1,NMIX+J3)
      WRITE (IOUT,90) ('----------------',I=1,MIN(8,J3))
      NMIX=NMIX+J3
    ENDDO
  ENDIF
  DEALLOCATE(ICNT,DVX,MIX,ALB,IFR)
  DEALLOCATE(NMC_SURF,NMC_NODE)
  DEALLOCATE(NBNODE_MACRO,NSURF_MACRO)
  !----
  !  Formats
  !----
   20 FORMAT (/9H  REGION:,8(I8,6X,1HI))
   30 FORMAT (9H MIXTURE:,8(I8,6X,1HI))
   40 FORMAT (9H  VOLUME:,1P,8(E13.6,2H I))
   50 FORMAT (6H MACRO,I6.6/9H ELEMENT:,8(3H  S,I6.6,6X,1HI,:)/(9X,8(3H  S,I6.6,6X,1HI,:)))
   60 FORMAT (9H  IN/OUT:,8(I6,2H /,I5,3H  I,:)/(9X,8(I6,2H /,I5,3H  I,:)))
   70 FORMAT (9H  ALBEDO:,1P,8(E13.5,3H  I,:)/(9X,8(E13.5,3H  I,:)))
   80 FORMAT (9H     DVX:,1P,8(E13.5,3H  I,:)/(9X,8(E13.5,3H  I,:)))
   90 FORMAT (9H --------,8(A16))
END SUBROUTINE SALMUS
