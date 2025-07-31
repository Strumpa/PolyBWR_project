!
!-----------------------------------------------------------------------
!
!Purpose:
! Compute the neutron flux and interface currents using the current
! iteration method for the multicell surfacic approximation.
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
! IPAS     total number of regions.
! NMCEL    total number of cells in the domain.
! NMERGE   total number of merged cells for which specific values
!          of the neutron flux and reactions rates are required.
!          Many cells with different position in the domain can
!          be merged before the neutron flux calculation if they
!          own the same generating cell (NMERGE.le.NMCEL).
! NGEN     total number of generating cells. A generating cell is
!          defined by its material and dimensions, irrespective of
!          its position in the domain (NGEN.le.NMERGE).
! IJAT     total number of distinct out-currents.
! NPIJ     size of cellwise scattering-reduced collision probability matrices.
! NPIS     size of cellwise scattering-reduced escape probability matrices.
! NPSS     size of cellwise scattering-reduced transmission probability matrices.
! EPSJ     stopping criterion for flux-current iterations.
! NUNKNO   total number of unknowns in vectors SUNKNO and FUNKNO.
! NMIX     nmber of out-currents (dimension of arrays MIX and DVX).
! NIFR     nmber of in-currents (dimension of arrays IFR and ALB).
! SUNKNO   input source vector.
! IMPX     print flag (equal to 0 for no print).
! NMC_NODE offset of the first volume in each generating cell.
! NMC_SURF offset of the first boundary surface in each generating cell.
! IFR      index-number of in-currents.
! ALB      transmission/albedo associated with each in-current.
! INUM     index-number of the merged cell associated to each cell.
! MIX      index-number of out-currents.
! DVX      weight associated with each out-current.
!          Note: IFR, ALB, MIX and DVX contains information to rebuild
!          the geometrical 'A' matrix.
! IGEN     index-number of the generating cell associated with each
!          merged cell.
! PIJW     cellwise scattering-reduced collision probability matrices.
! PISW     cellwise scattering-reduced escape probability matrices.
! PSJW     cellwise scattering-reduced collision probability matrices
!          for incoming neutrons.
! PSSW     cellwise scattering-reduced transmission probability matrices.
! 
!Parameters: input/output
! FUNKNO   unknown vector.
!
!-----------------------------------------------------------------------
!
SUBROUTINE MUSJJ2(IPAS,NMCEL,NMERGE,NGEN,IJAT,NPIJ,NPIS,NPSS,EPSJ,NUNKNO, &
  & NMIX,NIFR,FUNKNO,SUNKNO,IMPX,NMC_NODE,NMC_SURF,IFR,ALB,INUM,MIX,DVX, &
  & IGEN,PIJW,PISW,PSJW,PSSW)
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  !----
  !  SUBROUTINE ARGUMENTS
  !----
  INTEGER IPAS,NMCEL,NMERGE,NGEN,IJAT,NPIJ,NPIS,NUNKNO,NMIX,NIFR,IMPX, &
  & NMC_NODE(NGEN+1),NMC_SURF(NGEN+1),IFR(NIFR),INUM(NMCEL),MIX(NMIX), &
  & IGEN(NMERGE)
  REAL EPSJ,FUNKNO(NUNKNO),SUNKNO(NUNKNO),ALB(NIFR),DVX(NMIX),PIJW(NPIJ), &
  & PISW(NPIS),PSJW(NPIS),PSSW(NPSS)
  !----
  !  LOCAL VARIABLES
  !----
  REAL PIJ,PIS
  LOGICAL LOGTES
  PARAMETER (MAXIT=400,LACCFC=2,ICL1=3,ICL2=3)
  !----
  !  ALLOCATABLE ARRAYS
  !----
  INTEGER, DIMENSION(:), POINTER :: INDNMC
  DOUBLE PRECISION, DIMENSION(:), POINTER :: CIT0
  DOUBLE PRECISION, DIMENSION(:,:), POINTER :: CITR,AITR
  DOUBLE PRECISION, DIMENSION(:), POINTER :: WCURR
  !----
  !  SCRATCH STORAGE ALLOCATION
  !----
  ALLOCATE(INDNMC(NMERGE))
  ALLOCATE(CITR(3,IJAT),CIT0(IJAT),AITR(2,IJAT))
  ALLOCATE(WCURR(IJAT))
  !
  KNMC=0
  DO JKK=1,NMERGE
    JKG=IGEN(JKK)
    J2=NMC_NODE(JKG+1)-NMC_NODE(JKG)
    INDNMC(JKK)=KNMC
    KNMC=KNMC+J2
  ENDDO
  !
  DO I=1,IJAT
    WCURR(I)=1.0D0
    CIT0(I)=0.0D0
    CITR(1,I)=FUNKNO(IPAS+I)
  ENDDO
  !----
  !  COMPUTE PSJW * Q(*) CONTRIBUTION
  !----
  DO IKK=1,NMERGE
    IKG=IGEN(IKK)
    I2=NMC_NODE(IKG+1)-NMC_NODE(IKG)
    I3=NMC_SURF(IKG+1)-NMC_SURF(IKG)
    IT=0
    DO IK=1,IKK-1
      IT=IT+(NMC_SURF(IGEN(IK)+1)-NMC_SURF(IGEN(IK)))
    ENDDO
    IPSJ=0
    DO IK=1,IKG-1
      IPSJ=IPSJ+(NMC_NODE(IK+1)-NMC_NODE(IK))*(NMC_SURF(IK+1)-NMC_SURF(IK))
    ENDDO
    KNMC=INDNMC(IKK)
    DO I=1,I2
      DO IC=1,I3
        JCC=MIX(IT+IC)
        PBJ=PSJW(IPSJ+(I-1)*I3+IC)
        CIT0(JCC)=CIT0(JCC)+PBJ*DVX(IT+IC)*SUNKNO(KNMC+I)
      ENDDO
    ENDDO
  ENDDO
  !----
  !  COMPUTE NORMALIZATION VECTOR WCURR
  !----
  DO ICEL=1,NMCEL
    IKK=INUM(ICEL)
    IKG=IGEN(IKK)
    J3=NMC_SURF(IKG+1)-NMC_SURF(IKG)
    IT=0
    DO IK=1,IKK-1
      IT=IT+(NMC_SURF(IGEN(IK)+1)-NMC_SURF(IGEN(IK)))
    ENDDO
    IS=0
    DO IK=1,ICEL-1
      IS=IS+(NMC_SURF(IGEN(INUM(IK))+1)-NMC_SURF(IGEN(INUM(IK))))
    ENDDO
    IPSS=0
    DO IK=1,IKG-1
      IPSS=IPSS+(NMC_SURF(IK+1)-NMC_SURF(IK))**2
    ENDDO
    DO JC=1,J3
      J1=IFR(IS+JC)
      DO IC=1,J3
        PSS=PSSW(IPSS+(JC-1)*J3+IC)
        WCURR(J1)=WCURR(J1)-PSS*ALB(IS+JC)*DVX(IT+IC)
      ENDDO
    ENDDO
  ENDDO
  !
  ISTART=1
  TEST=0.0D0
  ITER=0
  10 ITER=ITER+1
  IF(ITER.GT.MAXIT) THEN
     WRITE(6,'(/47H MUSJJ2: *** WARNING *** MAXIMUM NUMBER OF ITER, &
     & 15HATIONS REACHED.)')
     GO TO 190
  ENDIF
  IT3=MOD(ITER,3)+1
  IT2=MOD(ITER-1,3)+1
  IT1=MOD(ITER-2,3)+1
  CITR(IT3,:IJAT)=CIT0(:IJAT)
  !----
  !  COMPUTE PSSW * J(-) CONTRIBUTION
  !----
  DO ICEL=1,NMCEL
    IKK=INUM(ICEL)
    IKG=IGEN(IKK)
    J3=NMC_SURF(IKG+1)-NMC_SURF(IKG)
    IT=0
    DO IK=1,IKK-1
      IT=IT+(NMC_SURF(IGEN(IK)+1)-NMC_SURF(IGEN(IK)))
    ENDDO
    IS=0
    DO IK=1,ICEL-1
      IS=IS+(NMC_SURF(IGEN(INUM(IK))+1)-NMC_SURF(IGEN(INUM(IK))))
    ENDDO
    IPSS=0
    DO IK=1,IKG-1
      IPSS=IPSS+(NMC_SURF(IK+1)-NMC_SURF(IK))**2
    ENDDO
    DO JC=1,J3
      J1=IFR(IS+JC)
      DO IC=1,J3
        J2=MIX(IT+IC)
        PSS=PSSW(IPSS+(JC-1)*J3+IC)
        CITR(IT3,J2)=CITR(IT3,J2)+PSS*ALB(IS+JC)*DVX(IT+IC)*CITR(IT2,J1)
      ENDDO
    ENDDO
  ENDDO
  !----
  !  NORMALIZATION
  !----
  S1=0.0D0
  S2=0.0D0
  DO I=1,IJAT
    S1=S1+WCURR(I)*CITR(IT3,I)
    S2=S2+CIT0(I)
  ENDDO
  ZNORM=S2/S1
  IF(ZNORM.LT.0.0D0) ZNORM=1.0D0
  CITR(IT3,:IJAT)=CITR(IT3,:IJAT)*ZNORM
  !----
  !  ONE/TWO PARAMETER ACCELERATION
  !----
  ALP=1.0D0
  BET=0.0D0
  LOGTES=(1+MOD(ITER-ISTART,ICL1+ICL2).GT.ICL1)
  IF(LOGTES) THEN
    AITR(1,:IJAT)=CITR(IT3,:IJAT)-CITR(IT2,:IJAT)
    AITR(2,:IJAT)=CITR(IT2,:IJAT)-CITR(IT1,:IJAT)
    DO ICEL=1,NMCEL
      IKK=INUM(ICEL)
      IKG=IGEN(IKK)
      J3=NMC_SURF(IKG+1)-NMC_SURF(IKG)
      IT=0
      DO IK=1,IKK-1
        IT=IT+(NMC_SURF(IGEN(IK)+1)-NMC_SURF(IGEN(IK)))
      ENDDO
      IS=0
      DO IK=1,ICEL-1
        IS=IS+(NMC_SURF(IGEN(INUM(IK))+1)-NMC_SURF(IGEN(INUM(IK))))
      ENDDO
      IPSS=0
      DO IK=1,IKG-1
        IPSS=IPSS+(NMC_SURF(IK+1)-NMC_SURF(IK))**2
      ENDDO
      DO JC=1,J3
        J1=IFR(IS+JC)
        DO IC=1,J3
          J2=MIX(IT+IC)
          PSS=PSSW(IPSS+(JC-1)*J3+IC)*ALB(IS+JC)*DVX(IT+IC)
          AITR(1,J2)=AITR(1,J2)-PSS*(CITR(IT3,J1)-CITR(IT2,J1))
          AITR(2,J2)=AITR(2,J2)-PSS*(CITR(IT2,J1)-CITR(IT1,J1))
        ENDDO
      ENDDO
    ENDDO
    IF((LACCFC.EQ.1).OR.(MOD(ITER-ISTART,ICL1+ICL2).EQ.ICL1)) THEN
      S1=0.0D0
      S2=0.0D0
      DO I=1,IJAT
        S1=S1+(CITR(IT3,I)-CITR(IT2,I))*AITR(1,I)
        S2=S2+AITR(1,I)*AITR(1,I)
      ENDDO
      IF(S2.EQ.0.0D0) THEN
        ISTART=ITER+1
      ELSE
        ALP=S1/S2
        IF(ALP.LE.0.0D0) THEN
          ISTART=ITER+1
          ALP=1.0D0
        ENDIF
      ENDIF
      DO I=1,IJAT
        CITR(IT3,I)=CITR(IT2,I)+ALP*(CITR(IT3,I)-CITR(IT2,I))
      ENDDO
    ELSE IF(LACCFC.EQ.2) THEN
      S1=0.0D0
      S2=0.0D0
      S3=0.0D0
      S4=0.0D0
      S5=0.0D0
      DO I=1,IJAT
        S1=S1+(CITR(IT3,I)-CITR(IT2,I))*AITR(1,I)
        S2=S2+AITR(1,I)*AITR(1,I)
        S3=S3+(CITR(IT3,I)-CITR(IT2,I))*AITR(2,I)
        S4=S4+AITR(1,I)*AITR(2,I)
        S5=S5+AITR(2,I)*AITR(2,I)
      ENDDO
      DET=S2*S5-S4*S4
      IF(DET.EQ.0.0D0) THEN
        ISTART=ITER+1
      ELSE
        ALP=(S5*S1-S4*S3)/DET
        BET=(S2*S3-S4*S1)/DET
        IF(ALP.LE.0.0D0) THEN
          ISTART=ITER+1
          ALP=1.0D0
          BET=0.0D0
        ENDIF
      ENDIF
      DO I=1,IJAT
        CITR(IT3,I)=CITR(IT2,I)+ALP*(CITR(IT3,I)-CITR(IT2,I))+ &
        & BET*(CITR(IT2,I)-CITR(IT1,I))
      ENDDO
    ENDIF
  ENDIF
  !----
  !  CHECK THE CONVERGENCE ERROR
  !----
  ERR1=0.0D0
  ERR2=0.0D0
  DO I=1,IJAT
    ERR1=MAX(ERR1,ABS(CITR(IT3,I)-CITR(IT2,I)))
    ERR2=MAX(ERR2,ABS(CITR(IT3,I)))
  ENDDO
  IF(IMPX.GT.3) WRITE(6,'(30H MUSJJ2: CURRENT ITERATION NB.,I4, &
  & 7H ERROR=,1P,E10.3,5H OVER,E10.3,15H NORMALIZATION=,E10.3, &
  & 14H ACCELERATION=,2E11.3,1H.)') ITER,ERR1,ERR2,ZNORM,ALP,BET/ALP
  IF(ITER.EQ.1) TEST=ERR1/ERR2
  IF((ITER.GT.20).AND.(ERR1/ERR2.GT.TEST)) THEN
    WRITE(6,'(/50H MUSJJ2: *** WARNING *** CONVERGENCE DIFFICULTIES.)')
    GO TO 190
  ENDIF
  IF(LOGTES.OR.(ERR1.GT.EPSJ*ERR2)) GO TO 10
  IF(IMPX.GT.2) WRITE(6,'(40H MUSJJ2: CURRENT CONVERGENCE AT ITERATIO, &
  & 5HN NB.,I4,7H ERROR=,1P,E10.3,5H OVER,E10.3,1H.)') ITER,ERR1,ERR2
  !
  190 FUNKNO(:IPAS)=0.0
  DO I=1,IJAT
    FUNKNO(IPAS+I)=REAL(CITR(IT3,I))
  ENDDO
  !----
  !  COMPUTE PISW * J(-) CONTRIBUTION
  !----
  DO ICEL=1,NMCEL
    IKK=INUM(ICEL)
    IKG=IGEN(IKK)
    I2=NMC_NODE(IKG+1)-NMC_NODE(IKG)
    I3=NMC_SURF(IKG+1)-NMC_SURF(IKG)
    IS=0
    DO IK=1,ICEL-1
      IS=IS+(NMC_SURF(IGEN(INUM(IK))+1)-NMC_SURF(IGEN(INUM(IK))))
    ENDDO
    IPIS=0
    DO IK=1,IKG-1
      IPIS=IPIS+(NMC_NODE(IK+1)-NMC_NODE(IK))*(NMC_SURF(IK+1)-NMC_SURF(IK))
    ENDDO
    KNMC=INDNMC(IKK)
    DO J=1,I2
      DO JC=1,I3
        J1=IFR(IS+JC)
        PIS=PISW(IPIS+(JC-1)*I2+J)
        FUNKNO(KNMC+J)=FUNKNO(KNMC+J)+PIS*ALB(IS+JC)*FUNKNO(IPAS+J1)
      ENDDO
    ENDDO
  ENDDO
  !----
  !  COMPUTE PIJW * Q(*) CONTRIBUTION
  !----
  DO IKK=1,NMERGE
    IKG=IGEN(IKK)
    I2=NMC_NODE(IKG+1)-NMC_NODE(IKG)
    IPIJ=0
    DO IK=1,IKG-1
      IPIJ=IPIJ+(NMC_NODE(IK+1)-NMC_NODE(IK))**2
    ENDDO
    KNMC=INDNMC(IKK)
    DO I=1,I2
      DO J=1,I2
        PIJ=PIJW(IPIJ+(I-1)*I2+J)
        FUNKNO(KNMC+J)=FUNKNO(KNMC+J)+PIJ*SUNKNO(KNMC+I)
      ENDDO
    ENDDO
  ENDDO
!----
!  SCRATCH STORAGE DEALLOCATION
!----
  DEALLOCATE(WCURR)
  DEALLOCATE(AITR,CIT0,CITR)
  DEALLOCATE(INDNMC)
  RETURN
  END
