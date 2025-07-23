!
!-----------------------------------------------------------------------
!
!Purpose:
! Solve N-group transport equation for fluxes using the current iteration
! method for the multicell surfacic approximation.
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
! KPSYS   pointer to the assembly matrices (L_PIJ signature). KPSYS is
!         an array of directories.
! IPTRK   pointer to the tracking (L_TRACK signature).
! IFTRAK  not used.
! IMPX    print flag (equal to zero for no print).
! NGEFF   number of energy groups processed in parallel.
! NGIND   energy group indices assign to the NGEFF set.
! IDIR    not used (=0 only for SYBIL).
! NREG    total number of regions for which specific values of the
!         neutron flux and reactions rates are required.
! NUNKNO  total number of unknowns in vectors SUNKNO and FUNKNO.
! MAT     index-number of the mixture type assigned to each volume.
! VOL     volumes.
! SUNKNO  input source vector.
! TITR    title.
!
!Parameters: input/output
! FUNKNO  unknown vector.
!
!-----------------------------------------------------------------------
!
  SUBROUTINE MUSF(KPSYS,IPTRK,IFTRAK,IMPX,NGEFF,NGIND,IDIR,NREG,NUNKNO,MAT, &
  & VOL,FUNKNO,SUNKNO,TITR)
  USE GANLIB
  !----
  !  SUBROUTINE ARGUMENTS
  !----
  TYPE(C_PTR) KPSYS(NGEFF),IPTRK
  CHARACTER TITR*72
  INTEGER NGEFF,NGIND(NGEFF),IFTRAK,IMPX,IDIR,NREG,NUNKNO,MAT(NREG)
  REAL VOL(NREG),FUNKNO(NUNKNO,NGEFF),SUNKNO(NUNKNO,NGEFF)
  !----
  !  LOCAL VARIABLES
  !----
  PARAMETER (IUNOUT=6,NSTATE=40)
  CHARACTER NAMLCM*12,NAMMY*12
  INTEGER ISTATE(NSTATE)
  REAL RSTATE(NSTATE)
  LOGICAL EMPTY,LCM
  !----
  !  ALLOCATABLE ARRAYS
  !----
  TYPE(C_PTR) PIJW_PTR,PISW_PTR,PSJW_PTR,PSSW_PTR
  INTEGER, ALLOCATABLE, DIMENSION(:) :: NMC_NODE,NMC_SURF,IFR,MIX,INUM,IGEN
  REAL, ALLOCATABLE, DIMENSION(:) :: ALB,DVX
  REAL, POINTER, DIMENSION(:) :: PIJW,PISW,PSJW,PSSW
  !
  IF(IDIR.NE.0) CALL XABORT('MUSF: EXPECTING IDIR=0')
  IF(IFTRAK.NE.0) CALL XABORT('MUSF: EXPECTING IFTRAK=0')
  IF(MAT(1).LT.0) CALL XABORT('MUSF: EXPECTING MAT(1)>=0')
  IF(VOL(1).LT.0.0) CALL XABORT('MUSF: EXPECTING VOL(1)>=0')
  CALL LCMINF(KPSYS(1),NAMLCM,NAMMY,EMPTY,ILONG,LCM)
  !----
  !  RECOVER MUST SPECIFIC PARAMETERS
  !----
  IF(IMPX.GT.2) WRITE(IUNOUT,'(//7H MUSF: ,A72)') TITR
  CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
  IF(NREG.NE.ISTATE(1)) THEN
    CALL XABORT('MUSF: STATE VECTOR HAS INVALID # OF ZONES.')
  ENDIF
  NMACRO=ISTATE(24) ! NGEN
  IF(NMACRO.EQ.0) CALL XABORT('MUSF: NO MACRO GEOMETRIES DEFINED.')
  NMCEL=NMACRO
  NMERGE=NMACRO
  NGEN=NMACRO
  ALLOCATE(IGEN(NMERGE),INUM(NMCEL))
  DO IK=1,NMERGE
    IGEN(IK)=IK
  ENDDO
  DO IK=1,NMCEL
    INUM(IK)=IK
  ENDDO
  IF(NMACRO.EQ.0) CALL XABORT('MUSP: MUST MODULE TRACKING IS MANDATORY.')
  ALLOCATE(NMC_NODE(NMACRO+1),NMC_SURF(NMACRO+1))
  CALL LCMGET(IPTRK,'NMC_NODE',NMC_NODE)
  CALL LCMGET(IPTRK,'NMC_SURF',NMC_SURF)
  CALL LCMGET(IPTRK,'EXCELTRACKOP',RSTATE)
  EPSJ=RSTATE(12)
  NMIX=NMC_SURF(NMACRO+1)
  NIFR=NMC_SURF(NMACRO+1)
  ALLOCATE(IFR(NIFR),ALB(NIFR),MIX(NMIX),DVX(NMIX))
  CALL LCMGET(IPTRK,'IFR',IFR)
  CALL LCMGET(IPTRK,'ALB',ALB)
  CALL LCMGET(IPTRK,'MIX',MIX)
  CALL LCMGET(IPTRK,'DVX',DVX)
  J1=0
  IJAT=0
  NPIJ=0
  NPIS=0
  NPSS=0
  DO IMACRO=1,NMACRO
    J2=NMC_NODE(IMACRO+1)-NMC_NODE(IMACRO)
    J3=NMC_SURF(IMACRO+1)-NMC_SURF(IMACRO)
    J1=J1+J2
    IJAT=IJAT+J3
    NPIJ=NPIJ+J2*J2
    NPIS=NPIS+J2*J3
    NPSS=NPSS+J3*J3
  ENDDO
  IF(J1.NE.NREG) CALL XABORT('MUSF: INVALID NREG.')
  IF(IJAT.NE.NMC_SURF(NMACRO+1)) CALL XABORT('MUSF: INVALID IJAT.')
  !----
  !  MAIN LOOP OVER ENERGY GROUPS.
  !----
  DO II=1,NGEFF
    IF(IMPX.GT.1) WRITE(IUNOUT,'(/23H MUSF: PROCESSING GROUP,I5, &
    & 6H WITH ,A,1H.)') NGIND(II),'MUSF'
    !
    IF(LCM) THEN
      CALL LCMGPD(KPSYS(II),'PIJW$SYBIL',PIJW_PTR)
      CALL LCMGPD(KPSYS(II),'PISW$SYBIL',PISW_PTR)
      CALL LCMGPD(KPSYS(II),'PSJW$SYBIL',PSJW_PTR)
      CALL LCMGPD(KPSYS(II),'PSSW$SYBIL',PSSW_PTR)
      !
      CALL C_F_POINTER(PIJW_PTR,PIJW,(/ NPIJ /))
      CALL C_F_POINTER(PISW_PTR,PISW,(/ NPIS /))
      CALL C_F_POINTER(PSJW_PTR,PSJW,(/ NPIS /))
      CALL C_F_POINTER(PSSW_PTR,PSSW,(/ NPSS /))
    ELSE
      ALLOCATE(PIJW(NPIJ),PISW(NPIS),PSJW(NPIS),PSSW(NPSS))
      CALL LCMGET(KPSYS(II),'PIJW$SYBIL',PIJW)
      CALL LCMGET(KPSYS(II),'PISW$SYBIL',PISW)
      CALL LCMGET(KPSYS(II),'PSJW$SYBIL',PSJW)
      CALL LCMGET(KPSYS(II),'PSSW$SYBIL',PSSW)
    ENDIF
    CALL MUSJJ2(NREG,NMCEL,NMERGE,NGEN,IJAT,NPIJ,NPIS,NPSS,EPSJ,NUNKNO, &
    & NMIX,NIFR,FUNKNO(1,II),SUNKNO(1,II),IMPX,NMC_NODE,NMC_SURF,IFR,ALB, &
    & INUM,MIX,DVX,IGEN,PIJW,PISW,PSJW,PSSW)
    IF(.NOT.LCM) DEALLOCATE(PSSW,PSJW,PISW,PIJW)
    !----
    ! END OF LOOP OVER ENERGY GROUPS
    !----
  ENDDO
  DEALLOCATE(DVX,MIX,ALB,IFR)
  DEALLOCATE(NMC_SURF,NMC_NODE,INUM,IGEN)
  RETURN
  END
