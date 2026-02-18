!-----------------------------------------------------------------------
!
! Purpose:
!   Compute the flux solution along a direction, over the spatial domain
!   in 1D/2D/3D Cartesian geometry, within the discrete ordinates (SN)
!   framework. Boltzmann and Boltzmann Fokker-Planck solvers. Albedo 
!   boundary conditions.
!
! Copyright:
!  Copyright (C) 2025 Ecole Polytechnique de Montreal
!  This library is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License as published by the Free Software Foundation; either
!  version 2.1 of the License, or (at your option) any later version
!
! Author(s): A. Hebert, A. A. Calloo, N. Martin and C. Bienvenue
!
!-----------------------------------------------------------------------
!
MODULE SNFCD_MOD

  USE SNSWC_MOD

CONTAINS

  SUBROUTINE SNFC1D(NUN,NGEFF,INCONV,LX,ISCHM,IELEM,NM,NMX,NMAT,NPQ,NSCT,DX,MAT,VOL,TOTAL,NCODE,ZCODE,QEXT,LFIXUP,LSHOOT, &
  & MN,DN,WX,CST,ISBS,NBS,ISBSM,BS,MAXL,FUNKNO,IBFP,EELEM,NME,ESCHM,ESTOPW,DELTAE,WE,ISLG,FLUXC)
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    !   Compute the flux solution along a direction, over the spatial domain
    !   in 1D Cartesian geometry, within the discrete ordinates (SN)
    !   framework. Boltzmann and Boltzmann Fokker-Planck solvers. Albedo 
    !   boundary conditions.
    !
    ! Parameters: input
    !   NUN     total number of unknowns.
    !   NGEFF   total number of energy groups processed in parallel.
    !   INCONV  energy group convergence flag (.FALSE. if converged).
    !   LX      total number of spatial cells.
    !   ISCHM   spatial discretization scheme index.
    !   IELEM   (order+1) of the spatial approximation polynomial.
    !   NM      total number of moments of the flux in space and in energy.
    !   NMX     number of incoming/outgoing boundary flux moments in space.
    !   NMAT    total number of materials.
    !   NPQ     total number of discrete angles.
    !   NSCT    total number of scattering cross-section moments.
    !   DX      factor containing first direction cosines (mu).
    !   MAT     material index in each spatial cell.
    !   VOL     cell volume.
    !   TOTAL   macroscopic total cross section in each material.
    !   NCODE   boundary condition indices.
    !   ZCODE   albedos.
    !   QEXT    angular external source term.
    !   LFIXUP  flag to enable negative flux fixup.
    !   LSHOOT  flag to enable shooting method.
    !   MN      Moment-to-discrete matrix.
    !   DN      Discrete-to-moment matrix.
    !   WX      spatial closure relation weighting factors.
    !   CST     Legendre coefficients for the polynomial approximations.
    !   ISBS    flag for boundary sources sources.
    !   NBS     total number of boundary sources.
    !   ISBSM   flag array for boundary fixed sources in each unit surface.
    !   BS      intensities of boundary fixed sources.
    !   MAXL    maximum number of unit surface in the boundary source.
    !   IBFP    type of energy proparation relation.
    !   EELEM   (order+1) of the energy approximation polynomial.
    !   NME     number of incoming/outgoing boundary flux moments in energy.
    !   ESCHM   energy discretization scheme index.
    !   ESTOPW  stopping powers at the upper and lower group boundaries.
    !   DELTAE  energy group width.
    !   WE      energy closure relation weighting factors.
    !   ISLG    flags indicating if an energy group is the lowest group.
    !   FLUXC   flux moments at the energy cutoff.
    !
    ! Parameters: input/output
    !   FUNKNO  Legendre components of the flux and boundary fluxes.
    !
    !-----------------------------------------------------------------------

    !  SUBROUTINE ARGUMENTS
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NUN,NGEFF,LX,ISCHM,IELEM,NM,NMX,NMAT,NPQ,NSCT,MAT(LX),NCODE(6),ISBS,NBS,ISBSM(2,NPQ,NGEFF), &
    & MAXL,IBFP,EELEM,NME,ESCHM
    LOGICAL, INTENT(IN) :: INCONV(NGEFF),LFIXUP,LSHOOT
    REAL, INTENT(IN) :: VOL(LX),TOTAL(0:NMAT,NGEFF),ZCODE(6),QEXT(NUN,NGEFF),BS(MAXL,NBS),WX(IELEM+1),CST(IELEM),MN(NPQ,NSCT), &
    & DN(NSCT,NPQ),DX(NPQ)
    INTEGER, INTENT(IN), OPTIONAL :: ISLG(NGEFF)
    REAL, INTENT(IN), OPTIONAL :: ESTOPW(0:NMAT,2,NGEFF),DELTAE(NGEFF),WE(EELEM+1)
    REAL, INTENT(INOUT) :: FUNKNO(NUN,NGEFF)
    REAL, INTENT(INOUT), OPTIONAL :: FLUXC(LX)

    ! LOCAL VARIABLES
    LOGICAL :: ISSHOOT
    INTEGER :: LFLX,LXNI,LFEP,M,IG,NS,IS,M0,IEL,IOF,IX
    DOUBLE PRECISION :: XNI(NMX),XNI1(NMX),XNI2(NMX),XNIA(NMX),XNIA1(NMX),XNIA2(NMX),XNIB(NMX),XNIB1(NMX),XNIB2(NMX)
    REAL :: MNT(NSCT,NPQ)
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:) :: FEP
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: FLUX
    
    !----
    ! INDEX VARIABLES
    !----
    LFLX=NM*LX*NSCT
    IF(LSHOOT) THEN
        LXNI=0
    ELSE
        LXNI=NMX*NPQ
    ENDIF
    IF(IBFP.NE.0) THEN
        LFEP=NME*LX*NPQ
    ELSE
        LFEP=0
    ENDIF

    !----
    ! EXTRACT FLUXES
    !----
    ALLOCATE(FLUX(NM,NSCT,LX,NGEFF),FEP(NME,LX,NPQ,NGEFF))
    FLUX = 0.0
    IF(IBFP.NE.0) FEP = RESHAPE(FUNKNO(LFLX+LXNI+1:LFLX+LXNI+LFEP,:NGEFF),(/NME,LX,NPQ,NGEFF/)) 
    MNT=TRANSPOSE(MN)

    !----
    ! GROUP LOOP
    !----
    DO IG=1,NGEFF
    
        XNI = 0.0
        IF(.NOT.INCONV(IG)) CYCLE

        ! SHOOTING METHOD WHEN THERE IS A NON-VACUUM RIGHT BOUNDARY CONDITION.
        ISSHOOT=(ZCODE(2).NE.0.0).AND.LSHOOT
        IF(ISSHOOT) THEN
            NS=6
        ELSE
            NS=2
        ENDIF

        ! LOOP OVER ALL DIRECTIONS
        DO M0=1,NPQ/2

            ! LOOP FOR SHOOTING METHOD
            DO IS=1,NS

                ! CHOOSE DIRECTION
                IF(MOD(IS,2).EQ.0) THEN
                    M=NPQ-M0+1 ! FORWARD
                ELSE
                    M=M0       ! BACKWARD
                ENDIF

                ! SHOOTING METHOD BOUNDARY CONDITIONS.
                IF(ISSHOOT) THEN
                    ! 1ST BACKWARD SWEEP
                    IF(IS.EQ.1) THEN
                        XNI=0.0
                        XNI1=0.0
                        XNI2=0.0
                    ! 1ST FORWARD SWEEP
                    ELSEIF(IS.EQ.2) THEN
                        XNIA1=0.0
                        IF(NCODE(1).EQ.4) THEN
                            XNIA1=XNI
                            XNI=0.0
                        ELSE
                            XNI=ZCODE(1)*XNI
                        ENDIF
                    ! 2ND BACKWARD SWEEP
                    ELSEIF(IS.EQ.3) THEN
                        XNIA2=0.0
                        XNIA=0.0
                        IF(NCODE(1).EQ.4) THEN
                            XNIA2=XNI
                        ELSE
                            XNIA=XNI
                        ENDIF
                        XNI=1.0
                    ! 2ND FORWARD SWEEP
                    ELSEIF(IS.EQ.4) THEN
                        IF(NCODE(1).EQ.4) THEN
                            XNIB1=XNI
                            XNI1=XNIA1/(1.0+XNIA1-XNIB1)
                            XNI=1.0
                        ELSE
                            XNI=ZCODE(1)*XNI
                        ENDIF
                    ! 3RD BACKWARD SWEEP
                    ELSEIF(IS.EQ.5) THEN
                        IF(NCODE(1).EQ.4) THEN
                            XNIB2=XNI
                            XNI2=XNIA2/(1.0+XNIA2-XNIB2)
                        XNI=XNI1
                        ELSE
                            XNIB=XNI
                            XNI=ZCODE(2)*XNIA/(1.0+ZCODE(2)*(XNIA-XNIB))
                        ENDIF
                    ! 3RD FORWARD SWEEP
                    ELSEIF(IS.EQ.6) THEN
                        XNI=ZCODE(1)*XNI
                        IF(NCODE(1).EQ.4) XNI=XNI2
                    ENDIF
                    
                ! NO SHOOTING METHOD BOUNDARY CONDITIONS
                ELSE
                    IF(.NOT.LSHOOT) THEN
                        IF(DX(M).GT.0.0.AND.NCODE(1).NE.4) THEN
                            DO IEL=1,EELEM
                                IOF=(M-1)*EELEM+IEL
                                FUNKNO(LFLX+LXNI+IOF,IG)=FUNKNO(LFLX+LXNI-IOF+1,IG)
                            ENDDO
                        ELSE
                            IF(NCODE(2).NE.4) THEN
                                DO IEL=1,EELEM
                                    IOF=(M-1)*EELEM+IEL
                                    FUNKNO(LFLX+LXNI+IOF,IG)=FUNKNO(LFLX+LXNI-IOF+1,IG)
                                ENDDO
                            ENDIF
                        ENDIF
                        XNI=0.0D0
                    ELSE
                        IF(IS.EQ.1) THEN
                            XNI=0.0D0
                        ELSE
                            XNI=ZCODE(1)*XNI
                        ENDIF
                    ENDIF
                ENDIF

                ! X-BOUNDARIES CONDITIONS (NO SHOOTING)
                IF(.NOT.LSHOOT) THEN
                    DO IEL=1,EELEM
                        IOF=(M-1)*EELEM+IEL
                        IF(DX(M).GT.0.0) THEN
                            XNI(IEL)=FUNKNO(LFLX+IOF,IG)*ZCODE(1)
                        ELSE
                            XNI(IEL)=FUNKNO(LFLX+IOF,IG)*ZCODE(2)
                        ENDIF
                    ENDDO
                ENDIF

                ! SWEEPING ACROSS THE 1D GRID
                IF(IBFP.EQ.0) THEN
                    CALL SNSWC1(NMAT,LX,NM,NMX,NSCT,ISCHM,IELEM,MAT,DX(M),TOTAL(0:NMAT,IG),VOL,CST,WX,BS,QEXT(1:LFLX,IG), &
                    & FLUX(:,:,:,IG),XNI,MNT(:,M),DN(:,M),IS,ISBS,ISBSM(:,M,IG),MAXL,NBS,LFIXUP,ISSHOOT,IBFP,EELEM,NME,ESCHM)
                ELSE
                    CALL SNSWC1(NMAT,LX,NM,NMX,NSCT,ISCHM,IELEM,MAT,DX(M),TOTAL(0:NMAT,IG),VOL,CST,WX,BS,QEXT(1:LFLX,IG), &
                    & FLUX(:,:,:,IG),XNI,MNT(:,M),DN(:,M),IS,ISBS,ISBSM(:,M,IG),MAXL,NBS,LFIXUP,ISSHOOT,IBFP,EELEM,NME,ESCHM, &
                    & ESTOPW(0:NMAT,1:2,IG),DELTAE(IG),WE,FEP(:,:,M,IG),QEXT(LFLX+LXNI+1+(M-1)*LX*NME:LFLX+LXNI+M*LX*NME,IG))
                ENDIF

                ! SAVE BOUNDARIES FLUX
                IF(.NOT.LSHOOT) THEN
                    DO IEL=1,EELEM
                        IOF=(M-1)*EELEM+IEL
                        FUNKNO(LFLX+IOF,IG)=REAL(XNI(IEL))
                    ENDDO
                ENDIF

            ENDDO ! END OF SHOOTING LOOP
        ENDDO ! END OF DIRECTION LOOP
    ENDDO ! END OF ENERGY LOOP

    !----
    ! FLUX SOLUTION AT THE CUTOFF ENERGY
    !----
    IF(IBFP.NE.0) THEN
        IF(ISLG(NGEFF).EQ.1) THEN
            DO IX=1,LX
                FLUXC(IX)=0.0
                DO M=1,NPQ
                    FLUXC(IX)=FLUXC(IX)+REAL(FEP(1,IX,M,NGEFF))*DN(1,M)
                ENDDO
            ENDDO
        ENDIF
    ENDIF

    !----
    ! SAVE FLUXES
    !----
    FUNKNO(1:LFLX,:NGEFF) = RESHAPE(FLUX,(/LFLX,NGEFF/))
    IF(IBFP.NE.0) FUNKNO(LFLX+LXNI+1:LFLX+LXNI+LFEP,:NGEFF) = RESHAPE(REAL(FEP),(/LFEP,NGEFF/))
    DEALLOCATE(FLUX,FEP)

    RETURN
  END SUBROUTINE SNFC1D

  SUBROUTINE SNFC2D(NUN,NGEFF,IMPX,INCONV,NGIND,LX,LY,ISCHM,IELEM,NM,NMX,NMY,NMAT,NPQ,NSCT,MAT,VOL,TOTAL,NCODE,ZCODE,QEXT,LFIXUP, &
  & DU,DE,W,MRMX,MRMY,DB,DA,MN,DN,WX,WY,CST,ISBS,NBS,ISBSM,BS,MAXL,FUNKNO,NKBA,IBFP,EELEM,NME,ESCHM,ESTOPW,DELTAE,WE,ISLG,FLUXC)
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    !   Compute the flux solution along a direction, over the spatial domain
    !   in 2D Cartesian geometry, within the discrete ordinates (SN)
    !   framework. Boltzmann and Boltzmann Fokker-Planck solvers. Albedo 
    !   boundary conditions.
    !
    ! Parameters: input
    !   NUN     total number of unknowns.
    !   NGEFF   total number of energy groups processed in parallel.
    !   IMPX    print flag (equal to zero for no print).
    !   INCONV  energy group convergence flag (.FALSE. if converged).
    !   NGIND   energy group indices assign to the NGEFF set.
    !   LX      total number of spatial cells.
    !   LY      total number of spatial cells.
    !   ISCHM   spatial discretization scheme index.
    !   IELEM   (order+1) of the spatial approximation polynomial.
    !   NM      total number of moments of the flux in space and in energy.
    !   NMX     number of incoming/outgoing boundary flux moments in space.
    !   NMY     number of incoming/outgoing boundary flux moments in space.
    !   NMAT    total number of materials.
    !   NPQ     total number of discrete angles.
    !   NSCT    total number of scattering cross-section moments.
    !   MAT     material index in each spatial cell.
    !   VOL     cell volume.
    !   TOTAL   macroscopic total cross section in each material.
    !   NCODE   boundary condition indices.
    !   ZCODE   albedos.
    !   QEXT    angular external source term.
    !   LFIXUP  flag to enable negative flux fixup.
    !   DU      direction cosines in the x-axis for each discrete angle.
    !   DE      direction cosines in the y-axis for each discrete angle.
    !   W       weights for each discrete angle.
    !   MRMX    angular mapping for reflexion along x-axis.
    !   MRMY    angular mapping for reflexion along y-axis.
    !   DA      factor containing first direction cosines (mu).
    !   DB      factor containing second direction cosines (eta).
    !   MN      Moment-to-discrete matrix.
    !   DN      Discrete-to-moment matrix.
    !   WX      spatial closure relation weighting factors.
    !   WY      spatial closure relation weighting factors.
    !   CST     Legendre coefficients for the polynomial approximations.
    !   ISBS    flag for boundary sources sources.
    !   NBS     total number of boundary sources.
    !   ISBSM   flag array for boundary fixed sources in each unit surface.
    !   BS      intensities of boundary fixed sources.
    !   MAXL    maximum number of unit surface in the boundary source.
    !   NKBA    number of macrocells along each axis.
    !   IBFP    type of energy proparation relation.
    !   EELEM   (order+1) of the energy approximation polynomial.
    !   NME     number of incoming/outgoing boundary flux moments in energy.
    !   ESCHM   energy discretization scheme index.
    !   ESTOPW  stopping powers at the upper and lower group boundaries.
    !   DELTAE  energy group width.
    !   WE      energy closure relation weighting factors.
    !   ISLG    flags indicating if an energy group is the lowest group.
    !   FLUXC   flux moments at the energy cutoff.
    !
    ! Parameters: input/output
    !   FUNKNO  Legendre components of the flux and boundary fluxes.
    !
    !-----------------------------------------------------------------------

    !$ use OMP_LIB

    !----
    ! VARIABLE DECLARATION
    !----
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NUN,NGEFF,IMPX,NGIND(NGEFF),LX,LY,ISCHM,IELEM,NM,NMX,NMY,NMAT,NPQ,NSCT,MAT(LX,LY),NCODE(6),MRMX(NPQ), &
    & MRMY(NPQ),ISBS,NBS,ISBSM(4,NPQ,NGEFF),MAXL,NKBA,IBFP,EELEM,NME,ESCHM
    LOGICAL, INTENT(IN) :: INCONV(NGEFF),LFIXUP
    REAL, INTENT(IN) :: VOL(LX,LY),TOTAL(0:NMAT,NGEFF),ZCODE(6),QEXT(NUN,NGEFF),DU(NPQ),DE(NPQ),W(NPQ),DB(LX,NPQ),DA(LX,LY,NPQ), &
    & BS(MAXL,NBS),WX(IELEM+1),WY(IELEM+1),CST(IELEM),MN(NPQ,NSCT),DN(NSCT,NPQ)
    INTEGER, INTENT(IN), OPTIONAL :: ISLG(NGEFF)
    REAL, INTENT(IN), OPTIONAL :: ESTOPW(0:NMAT,2,NGEFF),DELTAE(NGEFF),WE(EELEM+1)
    REAL, INTENT(INOUT) :: FUNKNO(NUN,NGEFF)
    REAL, INTENT(INOUT), OPTIONAL :: FLUXC(LX,LY)

    !----
    ! LOCAL VARIABLES
    !----
    INTEGER, PARAMETER :: IUNOUT=6, JMAP(4) = (/4, 3, 1, 2/)
    INTEGER :: LFLX,LXNI,LXNJ,LFEP,ITID,M,IG,IND,JND,IIND(8),NPQD(8),INDANG(NPQ,8),MRM,IPQD,IXM,IYM,IXP,IYP,IDI,NMACRO,MACROX, &
    & MACROY,NCELLX,NCELLY,ICEL,JCEL,NCWAVEF,ICWAVEF,MACROMAX,IX,IY
    REAL :: MNT(NSCT,NPQ),VU,VE
    INTEGER, ALLOCATABLE, DIMENSION(:) :: III,JJJ
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: XNI,XNJ
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: FLUX
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: FEP

    !----
    ! CHECK IF OPTIONAL ARGUMENTS ARE PRESENT
    !----
    IF(IBFP.NE.0) THEN
        IF (.NOT. PRESENT(ESTOPW)) CALL XABORT('ERROR: ESTOPW missing')
        IF (.NOT. PRESENT(DELTAE)) CALL XABORT('ERROR: DELTAE missing')
        IF (.NOT. PRESENT(WE))     CALL XABORT('ERROR: WE missing')
    ENDIF

    !----
    ! INDEX VARIABLES
    !----
    LFLX=NM*LX*LY*NSCT
    LXNI=NMX*LY*NPQ
    LXNJ=NMY*LX*NPQ
    LFEP=0
    IF(IBFP.NE.0) LFEP=NME*LX*LY*NPQ
    IF(NKBA.EQ.0) THEN
        MACROX=1
        MACROY=1
        NMACRO=1
        NCELLX=1
        NCELLY=1
    ELSE
        MACROX=NKBA
        MACROY=NKBA
        NMACRO=MACROX+MACROY-1
        NCELLX=1+(LX-1)/MACROX
        NCELLY=1+(LY-1)/MACROY
    ENDIF

    !----
    ! EXTRACT FLUXES
    !----
    ALLOCATE(FLUX(NM,NSCT,LX,LY,NGEFF),XNI(NMX,LY,NPQ,NGEFF),XNJ(NMY,LX,NPQ,NGEFF),FEP(NME,LX,LY,NPQ,NGEFF))
    FLUX = 0.0
    XNI = RESHAPE(FUNKNO(LFLX+1:LFLX+LXNI,1:NGEFF),(/NMX,LY,NPQ,NGEFF/))
    XNJ = RESHAPE(FUNKNO(LFLX+LXNI+1:LFLX+LXNI+LXNJ,1:NGEFF),(/NMY,LX,NPQ,NGEFF/))
    IF(IBFP.NE.0) FEP = RESHAPE(FUNKNO(LFLX+LXNI+LXNJ+1:LFLX+LXNI+LXNJ+LFEP,1:NGEFF),(/NME,LX,LY,NPQ,NGEFF/))
    MNT=TRANSPOSE(MN)

    !----
    ! SET OCTANT SWAPPING ORDER
    !----
    NPQD(:4)=0
    INDANG(:NPQ,:4)=0
    IIND(:)=0
    DO M=1,NPQ
        IF(W(M).EQ.0) CYCLE
        VU=DU(M)
        VE=DE(M)
        IF((VU.GE.0.0).AND.(VE.GE.0.0)) THEN
            IND=1
            JND=4
        ELSE IF((VU.LE.0.0).AND.(VE.GE.0.0)) THEN
            IND=2
            JND=3
        ELSE IF((VU.LE.0.0).AND.(VE.LE.0.0)) THEN
            IND=3
            JND=1
        ELSE
            IND=4
            JND=2
        ENDIF
        IIND(JND)=IND
        NPQD(IND)=NPQD(IND)+1
        INDANG(NPQD(IND),IND)=M
    ENDDO

    ! LOOP OVER QUADRANTS
    DO JND=1,4
        IND=IIND(JND)
        
        !----
        ! PREMILINARY LOOPS FOR SETTING BOUNDARY CONDITIONS
        !----

        !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(M,IG,MRM,IPQD) SHARED(XNI,XNJ)
        DO IG=1,NGEFF
        DO IPQD=1,NPQD(IND)
            IF(.NOT.INCONV(IG)) CYCLE
            M=INDANG(IPQD,IND)
            ! X-BOUNDARY
            IF((DU(M).GT.0.0.AND.NCODE(1).NE.4).OR.(DU(M).LT.0.0.AND.NCODE(2).NE.4))THEN
                MRM=MRMX(M)
                XNI(:,:,M,IG)=XNI(:,:,MRM,IG)
            ENDIF
            ! Y-BOUNDARY
            IF((DE(M).GT.0.0.AND.NCODE(3).NE.4).OR.(DE(M).LT.0.0.AND.NCODE(4).NE.4))THEN
                MRM=MRMY(M)
                XNJ(:,:,M,IG)=XNJ(:,:,MRM,IG)
            ENDIF
        ENDDO
        ENDDO
        !$OMP END PARALLEL DO

        ! WAVEFRONT LOOP (KBA)
        DO IDI=1,NMACRO
        
            IF(NKBA.EQ.0) THEN
                ! NO KBA SWAPPING
                ICEL=1
                JCEL=1
                NCWAVEF=1
            ELSE
                ! KBA SWAPPING INDICES
                MACROMAX=MIN(MACROX,IDI)
                ALLOCATE(III(MACROMAX),JJJ(MACROMAX))
                NCWAVEF=0
                DO ICEL=MAX(1,IDI-MACROY+1),MIN(MACROX,IDI)
                    JCEL=IDI-ICEL+1
                    NCWAVEF=NCWAVEF+1
                    IF(NCWAVEF.GT.MACROMAX) CALL XABORT('SNFC2D: MACROMAX OVERFLOW.')
                    III(NCWAVEF)=ICEL
                    JJJ(NCWAVEF)=JCEL
                ENDDO
            ENDIF

            !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ITID,IG,IPQD,ICWAVEF,ICEL,JCEL,IXM,IYM,IXP,IYP,M)
            !$OMP DO SCHEDULE(STATIC) COLLAPSE(3) REDUCTION(+:FLUX)
            DO IG=1,NGEFF               ! LOOP FOR GROUPS
            DO ICWAVEF=1,NCWAVEF        ! LOOP OVER MACROCELLS IN WAVEFRONT
            DO IPQD=1,NPQD(IND)         ! LOOP FOR ANGLES (IN OCTANT)
                IF(.NOT.INCONV(IG)) CYCLE
                M=INDANG(IPQD,IND)
                IF(W(M).EQ.0.0) CYCLE

                ITID = 0
                !$ ITID = omp_get_thread_num()
                IF(IMPX.GT.5) WRITE(IUNOUT,500) ITID,NGIND(IG),IPQD,ICEL,JCEL

                ! MACROCELL DIMENSIONS
                IF(NKBA.EQ.0) THEN
                    ! NO KBA SWAPPING
                    IXM = 1
                    IYM = 1
                    IXP = LX
                    IYP = LY
                ELSE
                    ! KBA SWAPPING (2D)
                    ICEL=III(ICWAVEF)
                    JCEL=JJJ(ICWAVEF)
                    IXM=(ICEL-1)*NCELLX+1
                    IYM=(JCEL-1)*NCELLY+1
                    IXP=MIN(ICEL*NCELLX,LX)
                    IYP=MIN(JCEL*NCELLY,LY)
                ENDIF
                
                ! SWEEPING ACROSS THE 2D GRID
                IF(IBFP.EQ.0) THEN
                    CALL SNSWC2(NMAT,IXM,IYM,IXP,IYP,LX,LY,NM,NMX,NMY,NSCT,ISCHM,IELEM,MAT,DA(:,:,M),DB(:,M),TOTAL(0:NMAT,IG), &
                    & VOL,CST,WX,WY,BS,QEXT(1:LFLX,IG),FLUX(:,:,:,:,IG),XNI(:,:,M,IG),XNJ(:,:,M,IG),MNT(:,M),DN(:,M),ZCODE,ISBS, &
                    & ISBSM(:,M,IG),MAXL,NBS,LFIXUP,IBFP,EELEM,NME,ESCHM)
                ELSE
                    CALL SNSWC2(NMAT,IXM,IYM,IXP,IYP,LX,LY,NM,NMX,NMY,NSCT,ISCHM,IELEM,MAT,DA(:,:,M),DB(:,M),TOTAL(0:NMAT,IG), &
                    & VOL,CST,WX,WY,BS,QEXT(1:LFLX,IG),FLUX(:,:,:,:,IG),XNI(:,:,M,IG),XNJ(:,:,M,IG),MNT(:,M),DN(:,M),ZCODE,ISBS, &
                    & ISBSM(:,M,IG),MAXL,NBS,LFIXUP,IBFP,EELEM,NME,ESCHM,ESTOPW(0:NMAT,1:2,IG),DELTAE(IG),WE,FEP(:,:,:,M,IG), &
                    & QEXT(LFLX+LXNI+LXNJ+1+(M-1)*LX*LY*NME:LFLX+LXNI+LXNJ+M*LX*LY*NME,IG))
                ENDIF
            
            ENDDO ! END OF DIRECTION LOOP
            ENDDO ! END OF MACROCELL LOOP
            ENDDO ! END OF ENERGY LOOP
            !$OMP END DO
            !$OMP END PARALLEL
        IF(NKBA.GT.0) DEALLOCATE(III,JJJ)
        ENDDO ! END OF WAVEFRONT LOOP
    ENDDO ! END OF QUADRANT LOOP

    !----
    ! FLUX SOLUTION AT THE CUTOFF ENERGY
    !----
    IF(IBFP.NE.0) THEN
        IF(ISLG(NGEFF).EQ.1) THEN
            DO IX=1,LX
            DO IY=1,LY
                FLUXC(IX,IY)=0.0
                DO M=1,NPQ
                    FLUXC(IX,IY)=FLUXC(IX,IY)+REAL(FEP(1,IX,IY,M,NGEFF))*DN(1,M)
                ENDDO
            ENDDO
            ENDDO
        ENDIF
    ENDIF

    !----
    ! SAVE FLUXES
    !----
    FUNKNO(1:LFLX,1:NGEFF) = RESHAPE(FLUX,(/LFLX,NGEFF/))
    FUNKNO(LFLX+1:LFLX+LXNI,1:NGEFF) = RESHAPE(XNI,(/LXNI,NGEFF/))
    FUNKNO(LFLX+LXNI+1:LFLX+LXNI+LXNJ,1:NGEFF) = RESHAPE(XNJ,(/LXNJ,NGEFF/))
    IF(IBFP.NE.0) FUNKNO(LFLX+LXNI+LXNJ+1:LFLX+LXNI+LXNJ+LFEP,1:NGEFF) = RESHAPE(FEP,(/LFEP,NGEFF/))
    DEALLOCATE(FLUX,XNI,XNJ,FEP)

    RETURN 
    500 FORMAT(16H SNFC2D: thread=,I8,12H --->(group=,I4,7H angle=,I4,11H macrocell=,3I5,1H))
  END SUBROUTINE SNFC2D
  
  SUBROUTINE SNFC3D(NUN,NGEFF,IMPX,INCONV,NGIND,LX,LY,LZ,ISCHM,IELEM,NM,NMX,NMY,NMZ,NMAT,NPQ,NSCT,MAT,VOL,TOTAL,NCODE,ZCODE,QEXT, &
  & LFIXUP,DU,DE,DZ,W,MRMX,MRMY,MRMZ,DC,DB,DA,MN,DN,WX,WY,WZ,CST,ISBS,NBS,ISBSM,BS,MAXL,FUNKNO,NKBA,IBFP,EELEM,NME,ESCHM,ESTOPW, &
  & DELTAE,WE,ISLG,FLUXC)
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    !   Compute the flux solution along a direction, over the spatial domain
    !   in 3D Cartesian geometry, within the discrete ordinates (SN)
    !   framework. Boltzmann and Boltzmann Fokker-Planck solvers. Albedo 
    !   boundary conditions.
    !
    ! Parameters: input
    !   NUN     total number of unknowns.
    !   NGEFF   total number of energy groups processed in parallel.
    !   IMPX    print flag (equal to zero for no print).
    !   INCONV  energy group convergence flag (.FALSE. if converged).
    !   NGIND   energy group indices assign to the NGEFF set.
    !   LX      total number of spatial cells.
    !   LY      total number of spatial cells.
    !   LZ      total number of spatial cells.
    !   ISCHM   spatial discretization scheme index.
    !   IELEM   (order+1) of the spatial approximation polynomial.
    !   NM      total number of moments of the flux in space and in energy.
    !   NMX     number of incoming/outgoing boundary flux moments in space.
    !   NMY     number of incoming/outgoing boundary flux moments in space.
    !   NMZ     number of incoming/outgoing boundary flux moments in space.
    !   NMAT    total number of materials.
    !   NPQ     total number of discrete angles.
    !   NSCT    total number of scattering cross-section moments.
    !   MAT     material index in each spatial cell.
    !   VOL     cell volume.
    !   TOTAL   macroscopic total cross section in each material.
    !   NCODE   boundary condition indices.
    !   ZCODE   albedos.
    !   QEXT    angular external source term.
    !   LFIXUP  flag to enable negative flux fixup.
    !   DU      direction cosines in the x-axis for each discrete angle.
    !   DE      direction cosines in the y-axis for each discrete angle.
    !   DZ      direction cosines in the z-axis for each discrete angle.
    !   W       weights for each discrete angle.
    !   MRMX    angular mapping for reflexion along x-axis.
    !   MRMY    angular mapping for reflexion along y-axis.
    !   MRMZ    angular mapping for reflexion along z-axis.
    !   DA      factor containing first direction cosines (mu).
    !   DB      factor containing second direction cosines (eta).
    !   DC      factor containing third direction cosines (xi).
    !   MN      Moment-to-discrete matrix.
    !   DN      Discrete-to-moment matrix.
    !   WX      spatial closure relation weighting factors.
    !   WY      spatial closure relation weighting factors.
    !   WZ      spatial closure relation weighting factors.
    !   CST     Legendre coefficients for the polynomial approximations.
    !   ISBS    flag for boundary sources sources.
    !   NBS     total number of boundary sources.
    !   ISBSM   flag array for boundary fixed sources in each unit surface.
    !   BS      intensities of boundary fixed sources.
    !   MAXL    maximum number of unit surface in the boundary source.
    !   NKBA    number of macrocells along each axis.
    !   IBFP    type of energy proparation relation.
    !   EELEM   (order+1) of the energy approximation polynomial.
    !   NME     number of incoming/outgoing boundary flux moments in energy.
    !   ESCHM   energy discretization scheme index.
    !   ESTOPW  stopping powers at the upper and lower group boundaries.
    !   DELTAE  energy group width.
    !   WE      energy closure relation weighting factors.
    !   ISLG    flags indicating if an energy group is the lowest group.
    !   FLUXC   flux moments at the energy cutoff.
    !
    ! Parameters: input/output
    !   FUNKNO  Legendre components of the flux and boundary fluxes.
    !
    !-----------------------------------------------------------------------

    !$ use OMP_LIB

    !----
    ! VARIABLE DECLARATION
    !----
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NUN,NGEFF,IMPX,NGIND(NGEFF),LX,LY,LZ,ISCHM,IELEM,NM,NMX,NMY,NMZ,NMAT,NPQ,NSCT,MAT(LX,LY,LZ),NCODE(6), &
    & MRMX(NPQ),MRMY(NPQ),MRMZ(NPQ),ISBS,NBS,ISBSM(6,NPQ,NGEFF),MAXL,NKBA,IBFP,EELEM,NME,ESCHM
    LOGICAL, INTENT(IN) :: INCONV(NGEFF),LFIXUP
    REAL, INTENT(IN) :: VOL(LX,LY,LZ),TOTAL(0:NMAT,NGEFF),ZCODE(6),QEXT(NUN,NGEFF),DU(NPQ),DE(NPQ),DZ(NPQ),W(NPQ),DC(LX,LY,NPQ), &
    & DB(LX,LZ,NPQ),DA(LY,LZ,NPQ),BS(MAXL,NBS),WX(IELEM+1),WY(IELEM+1),WZ(IELEM+1),CST(IELEM),MN(NPQ,NSCT),DN(NSCT,NPQ)
    INTEGER, INTENT(IN), OPTIONAL :: ISLG(NGEFF)
    REAL, INTENT(IN), OPTIONAL :: ESTOPW(0:NMAT,2,NGEFF),DELTAE(NGEFF),WE(EELEM+1)
    REAL, INTENT(INOUT) :: FUNKNO(NUN,NGEFF)
    REAL, INTENT(INOUT), OPTIONAL :: FLUXC(LX,LY,LZ)

    !----
    ! LOCAL VARIABLES
    !----
    INTEGER, PARAMETER :: IUNOUT=6,JMAP(8) = (/8,7,5,6,4,3,1,2/)
    INTEGER :: LFLX,LXNI,LXNJ,LXNK,LFEP,ITID,M,IG,IND,JND,IIND(8),NPQD(8),INDANG(NPQ,8),IPQD,IXM,IYM,IZM,IXP,IYP,IZP,IDI,NMACRO, &
    & MACROX,MACROY,MACROZ,NCELLX,NCELLY,NCELLZ,ICEL,JCEL,KCEL,NCWAVEF,ICWAVEF,MACROMAX,MRM,IX,IY,IZ
    REAL :: MNT(NSCT,NPQ),VU,VE,VZ
    INTEGER, ALLOCATABLE, DIMENSION(:) :: III,JJJ
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: XNI,XNJ,XNK
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: FEP
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: FLUX

    !----
    ! CHECK IF OPTIONAL ARGUMENTS ARE PRESENT
    !----
    IF(IBFP.NE.0) THEN
        IF (.NOT. PRESENT(ESTOPW)) CALL XABORT('ERROR: ESTOPW missing')
        IF (.NOT. PRESENT(DELTAE)) CALL XABORT('ERROR: DELTAE missing')
        IF (.NOT. PRESENT(WE))     CALL XABORT('ERROR: WE missing')
    ENDIF

    !----
    ! INDEX VARIABLES
    !----
    LFLX=NM*LX*LY*LZ*NSCT
    LXNI=NMX*LY*LZ*NPQ
    LXNJ=NMY*LX*LZ*NPQ
    LXNK=NMZ*LX*LY*NPQ
    LFEP=0
    IF(IBFP.NE.0) LFEP=NME*LX*LY*LZ*NPQ
    IF(NKBA.EQ.0) THEN
        MACROX=1
        MACROY=1
        MACROZ=1
        NMACRO=1
        NCELLX=1
        NCELLY=1
        NCELLZ=1
    ELSE
        MACROX=NKBA
        MACROY=NKBA
        MACROZ=NKBA
        NMACRO=MACROX+MACROY+MACROZ-2
        NCELLX=1+(LX-1)/MACROX
        NCELLY=1+(LY-1)/MACROY
        NCELLZ=1+(LZ-1)/MACROZ
    ENDIF

    !----
    ! EXTRACT FLUXES
    !----
    ALLOCATE(XNI(NMX,LY,LZ,NPQ,NGEFF),XNJ(NMY,LX,LZ,NPQ,NGEFF),XNK(NMZ,LX,LY,NPQ,NGEFF),FEP(NME,LX,LY,LZ,NPQ,NGEFF))
    ALLOCATE(FLUX(NM,NSCT,LX,LY,LZ,NGEFF))
    FLUX = 0.0D0
    XNI = RESHAPE(FUNKNO(LFLX+1:LFLX+LXNI,:NGEFF),(/NMX,LY,LZ,NPQ,NGEFF/))
    XNJ = RESHAPE(FUNKNO(LFLX+LXNI+1:LFLX+LXNI+LXNJ,:NGEFF),(/NMY,LX,LZ,NPQ,NGEFF/))
    XNK = RESHAPE(FUNKNO(LFLX+LXNI+LXNJ+1:LFLX+LXNI+LXNJ+LXNK,:NGEFF),(/NMZ,LX,LY,NPQ,NGEFF/))
    IF(IBFP.NE.0) FEP = RESHAPE(FUNKNO(LFLX+LXNI+LXNJ+LXNK+1:LFLX+LXNI+LXNJ+LXNK+LFEP,:NGEFF),(/NME,LX,LY,LZ,NPQ,NGEFF/))
    MNT=TRANSPOSE(MN)

    !----
    ! SET OCTANT SWAPPING ORDER
    !----
    NPQD(:8)=0
    INDANG(:NPQ,:8)=0
    IIND(:)=0
    DO M=1,NPQ
        VU=DU(M)
        VE=DE(M)
        VZ=DZ(M)
        IF(W(M).EQ.0) CYCLE
        IF((VU.GE.0.0).AND.(VE.GE.0.0).AND.(VZ.GE.0.0)) THEN
          IND=1
          JND=8
        ELSE IF((VU.LE.0.0).AND.(VE.GE.0.0).AND.(VZ.GE.0.0)) THEN
          IND=2
          JND=7
        ELSE IF((VU.LE.0.0).AND.(VE.LE.0.0).AND.(VZ.GE.0.0)) THEN
          IND=3
          JND=5
        ELSE IF((VU.GE.0.0).AND.(VE.LE.0.0).AND.(VZ.GE.0.0)) THEN
          IND=4
          JND=6
        ELSE IF((VU.GE.0.0).AND.(VE.GE.0.0).AND.(VZ.LE.0.0)) THEN
          IND=5
          JND=4
        ELSE IF((VU.LE.0.0).AND.(VE.GE.0.0).AND.(VZ.LE.0.0)) THEN
          IND=6
          JND=3
        ELSE IF((VU.LE.0.0).AND.(VE.LE.0.0).AND.(VZ.LE.0.0)) THEN
          IND=7
          JND=1
        ELSE
          IND=8
          JND=2
        ENDIF
        IIND(JND)=IND
        NPQD(IND)=NPQD(IND)+1
        INDANG(NPQD(IND),IND)=M
    ENDDO

    ! LOOP OVER OCTANTS
    DO JND=1,8
        IND=IIND(JND)
        
        !----
        ! PREMILINARY LOOPS FOR SETTING BOUNDARY CONDITIONS
        !----
        !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(M,IG,MRM,IPQD) SHARED(XNI,XNJ,XNK)
        DO IG=1,NGEFF
        DO IPQD=1,NPQD(IND)
            IF(.NOT.INCONV(IG)) CYCLE
            M=INDANG(IPQD,IND)
            ! X-BOUNDARY
            IF((DU(M).GT.0.0.AND.NCODE(1).NE.4).OR.(DU(M).LT.0.0.AND.NCODE(2).NE.4))THEN
                MRM=MRMX(M)
                XNI(:,:,:,M,IG)=XNI(:,:,:,MRM,IG)
            ENDIF
            ! Y-BOUNDARY
            IF((DE(M).GT.0.0.AND.NCODE(3).NE.4).OR.(DE(M).LT.0.0.AND.NCODE(4).NE.4))THEN
                MRM=MRMY(M)
                XNJ(:,:,:,M,IG)=XNJ(:,:,:,MRM,IG)
            ENDIF
            ! Z-BOUNDARY
            IF((DZ(M).GT.0.0.AND.NCODE(5).NE.4).OR.(DZ(M).LT.0.0.AND.NCODE(6).NE.4))THEN
                MRM=MRMZ(M)
                XNK(:,:,:,M,IG)=XNK(:,:,:,MRM,IG)
            ENDIF
        ENDDO
        ENDDO
        !$OMP END PARALLEL DO

        ! WAVEFRONT LOOP (KBA)
        DO IDI=1,NMACRO
        
            IF(NKBA.EQ.0) THEN
                ! NO KBA SWAPPING
                NCWAVEF=1
            ELSE
                ! KBA SWAPPING INDICES
                MACROMAX=MIN(MACROX,IDI)*MIN(MACROY,IDI)
                ALLOCATE(III(MACROMAX),JJJ(MACROMAX))
                NCWAVEF=0
                DO ICEL=MAX(1,IDI-MACROY-MACROZ+2),MIN(MACROX,IDI)
                DO JCEL=MAX(1,IDI-ICEL-MACROZ+2),MIN(MACROY,IDI-ICEL+1)
                    NCWAVEF=NCWAVEF+1
                    IF(NCWAVEF.GT.MACROMAX) CALL XABORT('SNFBC3: MACROMAX OVERFLOW.')
                    III(NCWAVEF)=ICEL
                    JJJ(NCWAVEF)=JCEL
                ENDDO
                ENDDO
            ENDIF

            !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ITID,IG,IPQD,ICWAVEF,ICEL,JCEL,KCEL,IXM,IYM,IZM,IXP,IYP,IZP,M)
            !$OMP DO SCHEDULE(STATIC) COLLAPSE(3) REDUCTION(+:FLUX)
            DO IG=1,NGEFF               ! LOOP FOR GROUPS
            DO ICWAVEF=1,NCWAVEF        ! LOOP OVER MACROCELLS
            DO IPQD=1,NPQD(IND)         ! LOOP FOR ANGLES
                IF(.NOT.INCONV(IG)) CYCLE
                M=INDANG(IPQD,IND)
                IF(W(M).EQ.0.0) CYCLE

                ! MACROCELL DIMENSIONS
                IF(NKBA.EQ.0) THEN
                    ! NO KBA SWAPPING
                    ICEL=1
                    JCEL=1
                    KCEL=1
                    IXM=1
                    IYM=1
                    IZM=1
                    IXP=LX
                    IYP=LY
                    IZP=LZ
                ELSE
                    ! KBA SWAPPING
                    ICEL=III(ICWAVEF)
                    JCEL=JJJ(ICWAVEF)
                    KCEL=IDI-ICEL-JCEL+2
                    IXM=(ICEL-1)*NCELLX+1
                    IYM=(JCEL-1)*NCELLY+1
                    IZM=(KCEL-1)*NCELLZ+1
                    IXP=MIN(ICEL*NCELLX,LX)
                    IYP=MIN(JCEL*NCELLY,LY)
                    IZP=MIN(KCEL*NCELLZ,LZ)
                ENDIF

                ITID = 0
                !$ ITID = omp_get_thread_num()
                IF(IMPX.GT.5) WRITE(IUNOUT,500) ITID,NGIND(IG),IPQD,ICEL,JCEL,KCEL
                
                ! SWEEPING ACROSS THE 3D GRID
                IF(IBFP.EQ.0) THEN
                    CALL SNSWC3(NMAT,IXM,IYM,IZM,IXP,IYP,IZP,LX,LY,LZ,NM,NMX,NMY,NMZ,NSCT,ISCHM,IELEM,MAT,DA(:,:,M),DB(:,:,M), &
                    & DC(:,:,M),TOTAL(0:NMAT,IG),VOL,CST,WX,WY,WZ,BS,QEXT(1:LFLX,IG),FLUX(:,:,:,:,:,IG),XNI(:,:,:,M,IG), &
                    & XNJ(:,:,:,M,IG),XNK(:,:,:,M,IG),MNT(:,M),DN(:,M),ZCODE,ISBS,ISBSM(:,M,IG),MAXL,NBS,LFIXUP,IBFP,EELEM, &
                    & NME,ESCHM)
                ELSE
                    CALL SNSWC3(NMAT,IXM,IYM,IZM,IXP,IYP,IZP,LX,LY,LZ,NM,NMX,NMY,NMZ,NSCT,ISCHM,IELEM,MAT,DA(:,:,M),DB(:,:,M), &
                    & DC(:,:,M),TOTAL(0:NMAT,IG),VOL,CST,WX,WY,WZ,BS,QEXT(1:LFLX,IG),FLUX(:,:,:,:,:,IG),XNI(:,:,:,M,IG), &
                    & XNJ(:,:,:,M,IG),XNK(:,:,:,M,IG),MNT(:,M),DN(:,M),ZCODE,ISBS,ISBSM(:,M,IG),MAXL,NBS,LFIXUP,IBFP,EELEM,NME, &
                    & ESCHM,ESTOPW(0:NMAT,1:2,IG),DELTAE(IG),WE,FEP(:,:,:,:,M,IG), &
                    & QEXT(LFLX+LXNI+LXNJ+LXNK+1+(M-1)*LX*LY*LZ*NME:LFLX+LXNI+LXNJ+LXNK+M*LX*LY*LZ*NME,IG))
                ENDIF
            
            ENDDO ! END OF DIRECTION LOOP
            ENDDO ! END OF MACROCELL LOOP
            ENDDO ! END OF ENERGY LOOP
            !$OMP END DO
            !$OMP END PARALLEL
            IF(NKBA.GT.0) DEALLOCATE(III,JJJ)
        ENDDO ! END OF WAVEFRONT LOOP
    ENDDO ! END OF OCTANT LOOP

    !----
    ! FLUX SOLUTION AT THE CUTOFF ENERGY
    !----
    IF(IBFP.NE.0) THEN
        IF(ISLG(NGEFF).EQ.1) THEN
            DO IX=1,LX
            DO IY=1,LY
            DO IZ=1,LZ
                FLUXC(IX,IY,IZ)=0.0
                DO M=1,NPQ
                    FLUXC(IX,IY,IZ)=FLUXC(IX,IY,IZ)+REAL(FEP(1,IX,IY,IZ,M,NGEFF))*DN(1,M)
                ENDDO
            ENDDO
            ENDDO
            ENDDO
        ENDIF
    ENDIF

    !----
    ! SAVE FLUXES
    !----
    FUNKNO(1:LFLX,:NGEFF) = RESHAPE(FLUX,(/LFLX,NGEFF/))
    FUNKNO(LFLX+1:LFLX+LXNI,:NGEFF) = RESHAPE(XNI,(/LXNI,NGEFF/))
    FUNKNO(LFLX+LXNI+1:LFLX+LXNI+LXNJ,:NGEFF) = RESHAPE(XNJ,(/LXNJ,NGEFF/))
    FUNKNO(LFLX+LXNI+LXNJ+1:LFLX+LXNI+LXNJ+LXNK,:NGEFF) = RESHAPE(XNK,(/LXNK,NGEFF/))
    IF(IBFP.NE.0) FUNKNO(LFLX+LXNI+LXNJ+LXNK+1:LFLX+LXNI+LXNJ+LXNK+LFEP,:NGEFF) = RESHAPE(FEP,(/LFEP,NGEFF/))
    DEALLOCATE(FLUX,XNI,XNJ,XNK,FEP)

    RETURN 
    500 FORMAT(16H SNFC3D: thread=,I8,12H --->(group=,I4,7H angle=,I4,11H macrocell=,3I5,1H))
  END SUBROUTINE SNFC3D
END MODULE SNFCD_MOD
