!-----------------------------------------------------------------------
!
! Purpose:
!   Perform one inner iteration for solving SN equations in 2D/3D hexagonal
!   geometry, within the discrete ordinates (SN) framework. Boltzmann 
!   solvers. KBA-like multithreading. Void boundary conditions.
!
! Copyright:
!  Copyright (C) 2025 Ecole Polytechnique de Montreal
!  This library ifs free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License as published by the Free Software Foundation; either
!  version 2.1 of the License, or (at your option) any later version
!
! Author(s): A. Hebert, A. A. Calloo and C. Bienvenue
!
!-----------------------------------------------------------------------
!
MODULE SNFHD_MOD

  USE SNSWH_MOD

CONTAINS

  SUBROUTINE SNFH2D(NUN,NGEFF,IMPX,INCONV,NGIND,NHEX,ISPLH,SIDE,IELEM,NM,NMX,NMY,NMAT,NPQ,NSCT,MAT,VOL,TOTAL,QEXT,LFIXUP,DU,DE, &
  & W,DB,DA,MN,DN,WX,WY,CST,LOZSWP,H2CMAP,FUNKNO,ISCHM)
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    !   Perform one inner iteration for solving SN equations in 2D hexagonal
    !   geometry, within the discrete ordinates (SN) framework. Boltzmann 
    !   solvers. KBA-like multithreading. Void boundary conditions.
    !
    ! Parameters: input
    !   NUN     total number of unknowns.
    !   NGEFF   total number of energy groups processed in parallel.
    !   IMPX    print flag (equal to zero for no print).
    !   INCONV  energy group convergence flag (.FALSE. if converged).
    !   NGIND   energy group indices assign to the NGEFF set.
    !   NHEX    total number of hexagonal cells.
    !   ISPLH   number of subcells in each dimension within a lozenge.
    !   SIDE    length of one side of the hexagon.
    !   IELEM   (order+1) of the spatial approximation polynomial.
    !   NM      total number of moments of the flux in space.
    !   NMX     number of incoming/outgoing boundary flux moments in space.
    !   NMY     number of incoming/outgoing boundary flux moments in space.
    !   NMAT    total number of materials.
    !   NPQ     total number of discrete angles.
    !   NSCT    total number of scattering cross-section moments.
    !   MAT     material index in each spatial cell.
    !   VOL     cell volume.
    !   TOTAL   macroscopic total cross section in each material.
    !   QEXT    angular external source term.
    !   LFIXUP  flag to enable negative flux fixup.
    !   DU      direction cosines in the x-axis for each discrete angle.
    !   DE      direction cosines in the y-axis for each discrete angle.
    !   W       weights for each discrete angle.
    !   DA      factor containing first direction cosines (mu).
    !   DB      factor containing second direction cosines (eta).
    !   MN      Moment-to-discrete matrix.
    !   DN      Discrete-to-moment matrix.
    !   WX      spatial closure relation weighting factors.
    !   WY      spatial closure relation weighting factors.
    !   CST     Legendre coefficients for the polynomial approximations.
    !   LOZSWP  lozenge sweeping order for the given angle.
    !   H2CMAP  mapping from hexagonal mesh to Cartesian axial coordinates.
    !   ISCHM   spatial discretization scheme index.
    !
    ! Parameters: input/output
    !   FUNKNO  Legendre components of the flux and boundary fluxes.
    !
    ! Comments:
    !   1. The direction of the axes I, J and D for the surface boundary 
    !      fluxes are shown in the diagram below. This means that 
    !      a) lozenge A has I- and D-boundaries (instead of I and J)
    !      b) lozenge B has I- and J-boundaries
    !      c) lozenge C has D- and J-boundaries (instead of I and J)
    !
    !                                  ^
    !                         j-axis   |
    !      ^  y-axis                   |          ^
    !      |                       _________     /    d-axis
    !      |                      /       / \   /
    !      |                     /   B   /   \
    !      | - - - ->           /       /     \
    !           x-axis         (-------(   A   )
    !                           \       \     /
    !                            \  C    \   / 
    !                             \_______\_/   \
    !                                            \   i-axis
    !                                             ^
    !
    !-----------------------------------------------------------------------

    !$ use OMP_LIB

    ! SUBROUTINE ARGUMENTS
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NUN,NGEFF,IMPX,NGIND(NGEFF),NHEX,ISPLH,IELEM,NM,NMX,NMY,NMAT,NPQ,NSCT,MAT(ISPLH,ISPLH,3,NHEX), &
    & LOZSWP(3,6),H2CMAP(3,NHEX),ISCHM
    LOGICAL, INTENT(IN) :: INCONV(NGEFF),LFIXUP
    REAL, INTENT(IN) :: SIDE,VOL(ISPLH,ISPLH,3,NHEX),TOTAL(0:NMAT,NGEFF),QEXT(NUN,NGEFF),DU(NPQ),DE(NPQ),W(NPQ), &
    & DB(ISPLH,ISPLH,3,NHEX,NPQ),DA(ISPLH,ISPLH,3,NHEX,NPQ),MN(NPQ,NSCT),DN(NSCT,NPQ),WX(IELEM+1),WY(IELEM+1),CST(IELEM)
    REAL, INTENT(INOUT) :: FUNKNO(NUN,NGEFF)

    ! LOCAL VARIABLES
    INTEGER :: NPQD(6),IIND(6),DCOORD,I,I_MC,IDI,IG,IHEX,IHEX_DOM,IIM,IND,IPQD,ITID,J,J_MC,JIM,JND,LFLX,M,NCOL,NRINGS,INDANG(NPQ,6)
    REAL :: JAC(2,2,3),VE,VU,MNT(NSCT,NPQ)
    DOUBLE PRECISION :: THETA
    INTEGER, PARAMETER :: IUNOUT=6
    REAL, PARAMETER :: PI=3.141592654
    INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: TMPMAT
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IHMAP
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: FLUX
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:,:) :: TMPXNI,TMPXNJ,TMPXND

    !----
    ! MAP MATERIAL VALUES TO CARTESIAN AXIAL COORDINATE MAP
    !----
    NRINGS=INT((SQRT(REAL((4*NHEX-1)/3))+1.)/2.)
    NCOL=2*NRINGS -1
    ALLOCATE(TMPMAT(ISPLH,ISPLH,3,NCOL,NCOL))
    ALLOCATE(IHMAP(NCOL,NCOL))
    TMPMAT(:,:,:,:,:) = -1
    IHMAP(:,:) = 0
    DO IHEX_DOM=1,NHEX
        TMPMAT(:,:,:,H2CMAP(1,IHEX_DOM),H2CMAP(2,IHEX_DOM))=MAT(:,:,:,IHEX_DOM)
        IHMAP(H2CMAP(1,IHEX_DOM),H2CMAP(2,IHEX_DOM))=IHEX_DOM
    ENDDO

    !----
    ! SCRATCH STORAGE ALLOCATION
    !----
    ALLOCATE(FLUX(NM,NSCT,3*ISPLH**2,NHEX,NGEFF))
    ALLOCATE(TMPXNI(IELEM,ISPLH,NCOL,NPQ,NGEFF))
    ALLOCATE(TMPXNJ(IELEM,ISPLH,NCOL,NPQ,NGEFF))
    ALLOCATE(TMPXND(IELEM,ISPLH,NCOL,NPQ,NGEFF))
    MNT=TRANSPOSE(MN)

    !----
    ! CONSTRUCT JACOBIAN MATRIX FOR EACH LOZENGE
    !----
    JAC = RESHAPE((/ 1., -SQRT(3.), 1., SQRT(3.), 2., 0., 1.,SQRT(3.), 2., 0., -1., SQRT(3.) /), SHAPE(JAC))
    JAC = (SIDE/2.)*JAC

    !----
    ! LENGTH OF FUNKNO COMPONENTS (IN ORDER)
    !----
    LFLX=3*NM*(ISPLH**2)*NHEX*NSCT

    !----
    ! SET DODECANT SWAPPING ORDER
    !----
    NPQD(:6)=0
    INDANG(:NPQ,:6)=0
    IIND(:6)=0
    DO M=1,NPQ
        VU=DU(M)
        VE=DE(M)
        IF(W(M).EQ.0) CYCLE
        THETA=0.0D0
        IF(VE.GT.0.0)THEN
            IF(VU.EQ.0.0)THEN
                THETA = PI/2
            ELSEIF(VU.GT.0.0)THEN
                THETA = ATAN(ABS(VE/VU))
            ELSEIF(VU.LT.0.0)THEN
                THETA = PI - ATAN(ABS(VE/VU))
            ENDIF
        ELSEIF(VE.LT.0.0)THEN
            IF(VU.EQ.0.0)THEN
                THETA = 3*PI/2
            ELSEIF(VU.LT.0.0)THEN
                THETA = PI + ATAN(ABS(VE/VU))
            ELSEIF(VU.GT.0.0)THEN
                THETA = 2.*PI - ATAN(ABS(VE/VU))
            ENDIF
        ENDIF
        IND=0
        IF((THETA.GT.0.0).AND.(THETA.LT.(PI/3.)))THEN
            IND=1
        ELSEIF((THETA.GT.(PI/3.)).AND.(THETA.LT.(2.*PI/3.)))THEN
            IND=2
        ELSEIF((THETA.GT.(2.*PI/3.)).AND.(THETA.LT.(PI)))THEN
            IND=3
        ELSEIF((THETA.GT.(PI)).AND.(THETA.LT.(4.*PI/3.)))THEN
            IND=4
        ELSEIF((THETA.GT.(4.*PI/3.)).AND.(THETA.LT.(5.*PI/3.)))THEN
            IND=5
        ELSEIF((THETA.GT.(5.*PI/3.)).AND.(THETA.LT.(2.*PI)))THEN
            IND=6
        ENDIF
        IIND(IND)=IND ! Assume IIND(I)=I in hexagonal geometry
        NPQD(IND)=NPQD(IND)+1
        INDANG(NPQD(IND),IND)=M
    ENDDO

    !----
    ! LOOP OVER DODECANTS
    !----
    FLUX(:NM,:NSCT,:3*ISPLH**2,:NHEX,:NGEFF)=0.0D0
    DO JND=1,6
        IND=IIND(JND)
        IF(IND.EQ.0) CYCLE ! Needed because of S2 LS (4 dir. for 6 sextants)
        TMPXNI(:IELEM,:ISPLH,:NCOL,:NPQ,:NGEFF)=0.0D0
        TMPXNJ(:IELEM,:ISPLH,:NCOL,:NPQ,:NGEFF)=0.0D0
        TMPXND(:IELEM,:ISPLH,:NCOL,:NPQ,:NGEFF)=0.0D0

        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ITID,IG,IPQD,J_MC,JIM,I_MC,IIM,I,J,IHEX,DCOORD,M,IDI)
        DO IDI=1,NCOL+NCOL-1 ! LOOP OVER WAVEFRONTS

            !$OMP DO COLLAPSE(2) SCHEDULE(STATIC,1)
            DO J_MC=MAX(1,IDI-NCOL+1),MIN(NCOL,IDI) ! LOOP OVER MACROCELLS IN WAVEFRONT
            DO IG=1,NGEFF ! LOOP FOR GROUPS
            DO IPQD=1,NPQD(IND) ! LOOP OVER ALL DIRECTIONS

                IF(.NOT.INCONV(IG)) CYCLE
                M=INDANG(IPQD,IND)
                IF(W(M).EQ.0.0) CYCLE

                ITID = 0
                !$ ITID = omp_get_thread_num()
                IF(IMPX.GT.5) WRITE(IUNOUT,400) ITID,NGIND(IG),IPQD

                ! Account for different sweep direction depending on angle
                JIM=J_MC
                IF((IND.EQ.1).OR.(IND.EQ.2).OR.(IND.EQ.3)) JIM=NCOL+1-JIM

                ! Find I coordinate of macrocell based on Jth coordinate, the wavefront number
                ! and occasionally the number of rings in the domain
                IF((IND.EQ.1).OR.(IND.EQ.4)) THEN
                    I_MC=IDI-J_MC+1
                ELSE
                    I_MC=IDI-J_MC+(NRINGS+1-J_MC)
                ENDIF
                IIM=I_MC

                ! Account for different sweep direction depending on angle
                IF((IND.EQ.2).OR.(IND.EQ.3).OR.(IND.EQ.4)) IIM=NCOL+1-IIM

                ! For IND 3 or 6, Cartesian axial coordinate map is swept vertically instead of
                ! horizontally. IM suffix is for 'IMmutable'
                I=IIM
                J=JIM
                IF((IND.EQ.3).OR.(IND.EQ.6))THEN
                    I=JIM
                    J=IIM
                ENDIF
                IF((I.GT.NCOL).OR.(I.LT.1)) CYCLE
                IF((J.GT.NCOL).OR.(J.LT.1)) CYCLE

                ! If within corners of Cartesian axial coordinate map (where there are no hexagons), skip loop
                IF(TMPMAT(1,1,1,I,J).EQ.-1) CYCLE

                ! Find DRAGON geometry hexagonal index using precomputed inverse map
                IHEX = IHMAP(I,J)
                IF(IHEX.EQ.0) CALL XABORT('SNFDH2: IHEX FAILURE.')

                ! Find D coordinate
                DCOORD = ABS(H2CMAP(3,IHEX))-NRINGS

                CALL SNSWH2(IND,JAC,LOZSWP,ISPLH,TMPXNI(:,:,:,IPQD,IG),TMPXNJ(:,:,:,IPQD,IG),TMPXND(:,:,:,IPQD,IG), &
                & FLUX(:,:,:,IHEX,IG),MAT(:,:,:,IHEX),TOTAL(:,IG),VOL(:,:,:,IHEX),DA(:,:,:,IHEX,M),DB(:,:,:,IHEX,M),QEXT(:,IG), &
                & NUN,NM,NMX,NMY,NSCT,ISCHM,IELEM,CST,WX,WY,NCOL,NMAT,IHEX,I,J,DCOORD,LFIXUP,MNT(:,M),DN(:,M))

            ENDDO ! END OF DIRECTION LOOP
            ENDDO ! END OF ENERGY LOOP
            ENDDO ! END OF MACROCELL LOOP
            !$OMP END DO
        ENDDO ! END OF WAVEFRONT LOOP
        !$OMP END PARALLEL
    ENDDO ! END OF DODECANTS LOOP

    ! SAVE FLUX INFORMATION
    DO IG=1,NGEFF
        IF(.NOT.INCONV(IG)) CYCLE
        FUNKNO(:LFLX,IG)=RESHAPE(REAL(FLUX(:IELEM**2,:NSCT,:3*ISPLH**2,:NHEX,IG)),(/ LFLX /) )
    ENDDO

    DEALLOCATE(FLUX,TMPXNI,TMPXNJ,TMPXND,TMPMAT,IHMAP)
    RETURN
    400 FORMAT(16H SNFKH2: thread=,I8,12H --->(group=,I4,7H angle=,I4,1H))
  END SUBROUTINE SNFH2D

  SUBROUTINE SNFH3D(NKBA,NUN,NGEFF,IMPX,INCONV,NGIND,NHEX,LZ,ISPLH,SIDE,IELEM,NM,NMX,NMY,NMZ,NMAT,NPQ,NSCT,MAT,VOL,TOTAL,NCODE, &
  & ZCODE,QEXT,LFIXUP,DU,DE,DZ,W,MRMZ,DC,DB,DA,MN,DN,WX,WY,WZ,CST,LOZSWP,H2CMAP,FUNKNO,ISCHM)
    !-----------------------------------------------------------------------
    !
    ! Purpose:
    !   Perform one inner iteration for solving SN equations in 3D hexagonal
    !   geometry, within the discrete ordinates (SN) framework. Boltzmann 
    !   solvers. KBA-like multithreading. Albedo boundary conditions.
    !
    ! Parameters: input
    !   NKBA    number of macrocells along each axis.
    !   NUN     total number of unknowns.
    !   NGEFF   total number of energy groups processed in parallel.
    !   IMPX    print flag (equal to zero for no print).
    !   INCONV  energy group convergence flag (.FALSE. if converged).
    !   NGIND   energy group indices assign to the NGEFF set.
    !   NHEX    total number of hexagonal cells.
    !   LZ      total number of spatial cells in the z-axis.
    !   ISPLH   number of subcells in each dimension within a lozenge.
    !   SIDE    length of one side of the hexagon.
    !   IELEM   (order+1) of the spatial approximation polynomial.
    !   NM      total number of moments of the flux in space.
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
    !   LOZSWP  lozenge sweeping order for the given angle.
    !   H2CMAP  mapping from hexagonal mesh to Cartesian axial coordinates.
    !   ISCHM   spatial discretization scheme index.
    !
    ! Parameters: input/output
    !   FUNKNO  Legendre components of the flux and boundary fluxes.
    !
    ! Comments:
    !   1. The direction of the axes I, J and D for the surface boundary 
    !      fluxes are shown in the diagram below. This means that 
    !      a) lozenge A has I- and D-boundaries (instead of I and J)
    !      b) lozenge B has I- and J-boundaries
    !      c) lozenge C has D- and J-boundaries (instead of I and J)
    !
    !                                  ^
    !                         j-axis   |
    !      ^  y-axis                   |          ^
    !      |                       _________     /    d-axis
    !      |                      /       / \   /
    !      |                     /   B   /   \
    !      | - - - ->           /       /     \
    !           x-axis         (-------(   A   )
    !                           \       \     /
    !                            \  C    \   / 
    !                             \_______\_/   \
    !                                            \   i-axis
    !                                             ^
    !
    !-----------------------------------------------------------------------

    !$ use OMP_LIB

    ! SUBROUTINE ARGUMENTS
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: NKBA,NUN,NGEFF,IMPX,NGIND(NGEFF),NHEX,LZ,ISPLH,IELEM,NM,NMX,NMY,NMZ,NMAT,NPQ,NSCT, &
    MAT(ISPLH,ISPLH,3,NHEX,LZ),NCODE(6),MRMZ(NPQ),LOZSWP(3,6),H2CMAP(3,NHEX),ISCHM
    LOGICAL, INTENT(IN) :: INCONV(NGEFF),LFIXUP
    REAL, INTENT(IN) :: SIDE,VOL(ISPLH,ISPLH,3,NHEX,LZ),ZCODE(6),TOTAL(0:NMAT,NGEFF),QEXT(NUN,NGEFF),DU(NPQ),DE(NPQ),DZ(NPQ), &
    & W(NPQ),DC(ISPLH*ISPLH*3*NHEX,1,NPQ),DB(ISPLH*ISPLH*3*NHEX,LZ,NPQ),DA(1,LZ,NPQ),MN(NPQ,NSCT),DN(NSCT,NPQ),WX(IELEM+1), &
    & WY(IELEM+1),WZ(IELEM+1),CST(IELEM)
    REAL, INTENT(INOUT) :: FUNKNO(NUN,NGEFF)

    ! LOCAL VARIABLES
    INTEGER :: NPQD(12),IIND(12),DCOORD,I,I_MC,IDI,IG,IIM,IND,IPQD,ITID,J,J_MC,JIM,JND,LFLX,M,NCOL,NRINGS,INDANG(NPQ,12),I_END, &
    & I_STT,ICEL,IHEX_XY,IND_XY,IOF,IZ,J_END,J_STT,K,K_MC,LXNK,M1,MACROZ,NCEL,NCELLZ,NMAX,JOF
    REAL :: JAC(2,2,3),VE,VU,VZ,MNT(NSCT,NPQ)
    DOUBLE PRECISION :: THETA
    INTEGER, PARAMETER :: IUNOUT=6
    REAL, PARAMETER :: PI=3.141592654
    INTEGER, ALLOCATABLE, DIMENSION(:) :: III,JJJ,KKK
    INTEGER, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: TMPMAT
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IHMAP
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:,:,:) :: FLUX,TMPXNI,TMPXNJ,TMPXND
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:,:,:,:,:) :: TMPXNK

    !----
    ! MAP MATERIAL VALUES TO CARTESIAN AXIAL COORDINATE MAP
    !----
    NRINGS=INT((SQRT(  REAL((4*NHEX-1)/3)  )+1.)/2.)
    NCOL=2*NRINGS-1
    ALLOCATE(TMPMAT(ISPLH,ISPLH,3,NCOL,NCOL,LZ))
    ALLOCATE(IHMAP(NCOL,NCOL))
    TMPMAT(:,:,:,:,:,:) = -1
    IHMAP(:,:) = 0
    DO IZ=1,LZ
        DO IHEX_XY=1,NHEX
            TMPMAT(:,:,:,H2CMAP(1,IHEX_XY),H2CMAP(2,IHEX_XY),IZ)=MAT(:,:,:,IHEX_XY,IZ)
            IHMAP(H2CMAP(1,IHEX_XY),H2CMAP(2,IHEX_XY))=IHEX_XY
        ENDDO
    ENDDO

    !----
    ! SCRATCH STORAGE ALLOCATION
    !----
    ALLOCATE(FLUX(NM,NSCT,3*ISPLH**2,NHEX,LZ,NGEFF))
    ALLOCATE(TMPXNI(NMX,ISPLH,NCOL,LZ,NPQ,NGEFF))
    ALLOCATE(TMPXNJ(NMY,ISPLH,NCOL,LZ,NPQ,NGEFF))
    ALLOCATE(TMPXND(NMX,ISPLH,NCOL,LZ,NPQ,NGEFF))
    ALLOCATE(TMPXNK(NMZ,ISPLH,ISPLH,3,NHEX,NPQ,NGEFF))
    MNT=TRANSPOSE(MN)

    !----
    !  CONSTRUCT JACOBIAN MATRIX FOR EACH LOZENGE
    !----
    JAC = RESHAPE((/ 1., -SQRT(3.), 1., SQRT(3.), 2., 0., 1.,SQRT(3.), 2., 0., -1., SQRT(3.) /), SHAPE(JAC))
    JAC = (SIDE/2.)*JAC

    !----
    !  LENGTH OF FUNKNO COMPONENTS (IN ORDER)
    !----
    LFLX=3*NM*(ISPLH**2)*NHEX*LZ*NSCT
    LXNK=3*NMZ*(ISPLH**2)*NHEX

    !----
    !  NUMBER OF MACROCELLS (MACRO*)
    !  NUMBER OF LZ LAYERS IN EACH MACROCELL (NCELL*)
    !----
    MACROZ=1
    IF(NKBA.NE.0) MACROZ=NKBA
    NCELLZ =(LZ+MACROZ-1)/MACROZ
      
    !----
    !  SET DODECANT SWAPPING ORDER.
    !----
    NPQD(:12)=0
    INDANG(:NPQ,:12)=0
    IIND(:12)=0
    DO M=1,NPQ
        VU=DU(M)
        VE=DE(M)
        VZ=DZ(M)
        IF(W(M).EQ.0) CYCLE
        THETA=0.0D0
        IF(VE.GT.0.0)THEN
            IF(VU.EQ.0.0)THEN
                THETA = PI/2
            ELSEIF(VU.GT.0.0)THEN
                THETA = ATAN(ABS(VE/VU))
            ELSEIF(VU.LT.0.0)THEN
                THETA = PI - ATAN(ABS(VE/VU))
            ENDIF
        ELSEIF(VE.LT.0.0)THEN
            IF(VU.EQ.0.0)THEN
                THETA = 3*PI/2
            ELSEIF(VU.LT.0.0)THEN
                THETA = PI + ATAN(ABS(VE/VU))
            ELSEIF(VU.GT.0.0)THEN
                THETA = 2.*PI - ATAN(ABS(VE/VU))
            ENDIF
        ENDIF
        IND=0
        IF(VZ.GE.0.0)THEN
            IF((THETA.GT.0.0).AND.(THETA.LT.(PI/3.)))THEN
                IND=1
            ELSEIF((THETA.GT.(PI/3.)).AND.(THETA.LT.(2.*PI/3.)))THEN
                IND=2
            ELSEIF((THETA.GT.(2.*PI/3.)).AND.(THETA.LT.(PI)))THEN
                IND=3
            ELSEIF((THETA.GT.(PI)).AND.(THETA.LT.(4.*PI/3.)))THEN
                IND=4
            ELSEIF((THETA.GT.(4.*PI/3.)).AND.(THETA.LT.(5.*PI/3.)))THEN
                IND=5
            ELSEIF((THETA.GT.(5.*PI/3.)).AND.(THETA.LT.(2.*PI)))THEN
                IND=6
            ENDIF
        ELSEIF(VZ.LT.0.0)THEN
            IF((THETA.GT.0.0).AND.(THETA.LT.(PI/3.)))THEN
                IND=7
            ELSEIF((THETA.GT.(PI/3.)).AND.(THETA.LT.(2.*PI/3.)))THEN
                IND=8
            ELSEIF((THETA.GT.(2.*PI/3.)).AND.(THETA.LT.(PI)))THEN
                IND=9
            ELSEIF((THETA.GT.(PI)).AND.(THETA.LT.(4.*PI/3.)))THEN
                IND=10
            ELSEIF((THETA.GT.(4.*PI/3.)).AND.(THETA.LT.(5.*PI/3.)))THEN
                IND=11
            ELSEIF((THETA.GT.(5.*PI/3.)).AND.(THETA.LT.(2.*PI)))THEN
                IND=12
            ENDIF
        ENDIF
        IIND(IND)=IND ! Assume IIND(I)=I in hexagonal geometry
        NPQD(IND)=NPQD(IND)+1
        INDANG(NPQD(IND),IND)=M
    ENDDO

    !----
    ! LOOP OVER DODECANTS
    !----
    FLUX(:NM,:NSCT,:3*ISPLH**2,:NHEX,:LZ,:NGEFF)=0.0D0
    DO JND=1,12
        IND=IIND(JND)
        IND_XY=MOD(IND-1,6)+1
        IF(IND.EQ.0) CYCLE ! Needed because of S2 LS (8 dir. for 12 dodecants)
        TMPXNI(:NMX,:ISPLH,:NCOL,:LZ,:NPQ,:NGEFF)=0.0D0
        TMPXNJ(:NMY,:ISPLH,:NCOL,:LZ,:NPQ,:NGEFF)=0.0D0
        TMPXND(:NMX,:ISPLH,:NCOL,:LZ,:NPQ,:NGEFF)=0.0D0
        TMPXNK(:NMZ,:ISPLH,:ISPLH,:3,:NHEX,:NPQ,:NGEFF)=0.0D0

        ! PRELIMINARY LOOPS FOR SETTING BOUNDARY CONDITIONS.
        IF((NCODE(5).NE.1).or.(NCODE(6).NE.1)) THEN
            DO IG=1,NGEFF
            DO IPQD=1,NPQD(IND)
                IF(.NOT.INCONV(IG)) CYCLE
                M=INDANG(IPQD,IND)
                VZ=DZ(M)
                ! Z-BOUNDARY
                IF(VZ.GT.0.0) THEN
                    M1=MRMZ(M)
                    IF(NCODE(5).NE.4) THEN
                    IOF=(M-1)*(LXNK)
                    JOF=(M1-1)*(LXNK)
                    FUNKNO(LFLX+IOF+1:LFLX+IOF+LXNK,IG)=FUNKNO(LFLX+JOF+1:LFLX+JOF+LXNK,IG)
                    ENDIF
                ELSEIF(VZ.LT.0.0) THEN
                    M1=MRMZ(M)
                    IF(NCODE(6).NE.4) THEN
                    IOF=(M-1)*(LXNK)
                    JOF=(M1-1)*(LXNK)
                    FUNKNO(LFLX+IOF+1:LFLX+IOF+LXNK,IG)=FUNKNO(LFLX+JOF+1:LFLX+JOF+LXNK,IG)
                    ENDIF
                ENDIF
            ENDDO
            ENDDO
        ENDIF

        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ITID,IG,IPQD,ICEL,I,J,K_MC,IHEX_XY,DCOORD,IOF,M,IDI,NMAX,NCEL,J_MC,JIM,I_MC,IIM) &
        !$OMP PRIVATE(I_END,I_STT,J_STT,J_END,III,JJJ,KKK)
        ! LOOP OVER WAVEFRONTS 
        DO IDI=1,NCOL+NCOL+MACROZ-2

            ! SET SWEEP INDICES
            NMAX=MIN(NCOL,IDI)*MIN(NCOL,IDI)
            ALLOCATE(III(NMAX),JJJ(NMAX),KKK(NMAX))
            NCEL=0
            J_STT=MAX(1,IDI-NCOL-MACROZ+2)
            J_END=MIN(NCOL,IDI)
            DO J_MC=J_STT,J_END
            
                JIM=J_MC
                ! Account for different sweep direction depending on angle
                IF((IND_XY.EQ.1).OR.(IND_XY.EQ.2).OR.(IND_XY.EQ.3)) JIM=NCOL+1-JIM

                IF((IND_XY.EQ.1).OR.(IND_XY.EQ.4)) THEN
                    I_STT=MAX(1,IDI-J_MC-MACROZ+2)
                    I_END=MIN(NCOL,IDI-J_MC+1)
                ELSE
                    I_STT=MAX(1,NRINGS-J_MC+1)
                    I_END=MIN(NCOL,IDI-J_MC+(NRINGS+1-J_MC))
                ENDIF

                DO I_MC=I_STT,I_END
                    IIM=I_MC
                    ! Account for different sweep direction depending on angle
                    IF((IND_XY.EQ.2).OR.(IND_XY.EQ.3).OR.(IND_XY.EQ.4)) IIM=NCOL+1-IIM

                    ! For IND_XY 3 or 6, Cartesian axial coordinate map is swept 
                    ! vertically instead of horizontally. IM suffix is for 'IMmutable'
                    I=IIM
                    J=JIM
                    IF((IND_XY.EQ.3).OR.(IND_XY.EQ.6)) THEN
                        I=JIM
                        J=IIM
                    ENDIF

                    ! If within corners of Cartesian axial coordinate map (where
                    ! there are no hexagons), skip loop
                    IF(TMPMAT(1,1,1,I,J,1).EQ.-1) CYCLE

                    ! Find I coordinate of macrocell
                    IF((IND_XY.EQ.1).OR.(IND_XY.EQ.4)) THEN
                        K_MC=IDI-I_MC-J_MC+2
                    ELSE
                        K_MC=IDI-I_MC+NRINGS-((J_MC-1)*2)
                    ENDIF
                    IF((K_MC.GT.MACROZ)) CYCLE
                    K=K_MC

                    NCEL=NCEL+1
                    IF(NCEL.GT.NMAX) CALL XABORT('SNFDH3: NMAX OVERFLOW.')
                    III(NCEL)=I
                    JJJ(NCEL)=J
                    KKK(NCEL)=K
                ENDDO ! I_MC
            ENDDO ! J_MC

            !$OMP DO SCHEDULE(STATIC,1) COLLAPSE(2)
            DO ICEL=1,NCEL ! LOOP FOR MACROCELLS
            DO IG=1,NGEFF ! LOOP FOR GROUPS
            DO IPQD=1,NPQD(IND) ! LOOP OVER ALL DIRECTIONS

                IF(.NOT.INCONV(IG)) CYCLE
                M=INDANG(IPQD,IND)
                IF(W(M).EQ.0.0) CYCLE

                ITID = 0
                !$ ITID = omp_get_thread_num()
                IF(IMPX.GT.5) WRITE(IUNOUT,500) ITID,NGIND(IG),IPQD

                I=III(ICEL)
                J=JJJ(ICEL)
                K_MC=KKK(ICEL)

                ! Find in X-Y plane DRAGON geometry hexagonal index using I and J
                IHEX_XY = IHMAP(I,J)
                IF(IHEX_XY.EQ.0) CALL XABORT('SNFDH3: IHEX_XY FAILURE.')

                ! Find D coordinate
                DCOORD = ABS(H2CMAP(3,IHEX_XY))-NRINGS

                IF(IDI.EQ.1) THEN
                    ! PICK UP BOUNDARY ELEMENTS
                    IF((NCODE(5).NE.1).or.(NCODE(6).NE.1)) THEN
                        TMPXNK(:NMZ,:ISPLH,:ISPLH,:3,:NHEX,IPQD,IG)=RESHAPE(FUNKNO(LFLX+(M-1)*LXNK+1:LFLX+M*LXNK,IG), &
                        & (/NMZ,ISPLH,ISPLH,3,NHEX/))
                    ENDIF
                    ! ACCOUNT FOR ALBEDO IN BOUNDARY ELEMENTS
                    IF(IND.LT.7) THEN
                        TMPXNK(:NMZ,:ISPLH,:ISPLH,:3,:NHEX,IPQD,IG)=TMPXNK(:NMZ,:ISPLH,:ISPLH,:3,:NHEX,IPQD,IG)*ZCODE(5)
                    ELSE
                        TMPXNK(:NMZ,:ISPLH,:ISPLH,:3,:NHEX,IPQD,IG)=TMPXNK(:NMZ,:ISPLH,:ISPLH,:3,:NHEX,IPQD,IG)*ZCODE(6)
                    ENDIF
                ENDIF

                CALL SNSWH3(IND,IND_XY,JAC,LOZSWP,ISPLH,LZ,TMPXNI(:,:,:,:,IPQD,IG),TMPXNJ(:,:,:,:,IPQD,IG), &
                & TMPXND(:,:,:,:,IPQD,IG),TMPXNK(:,:,:,:,IHEX_XY,IPQD,IG),FLUX(:,:,:,:,:,IG),MAT(:,:,:,:,:),TOTAL(:,IG), &
                & VOL(:,:,:,:,:),DA(:,:,M),DB(:,:,M),DC(:,:,M),QEXT(:,IG),NUN,NM,NMX,NMY,NMZ,NSCT,ISCHM,IELEM,CST,WX,WY,WZ, &
                & NCOL,NHEX,NMAT,IHEX_XY,I,J,DCOORD,LFIXUP,MNT(:,M),DN(:,M),K_MC,NCELLZ)

                ! SAVE K-BOUNDARY CONDITIONS IF NOT VOID B.C. (end of sweep)
                IF(IDI.EQ.NCOL+NCOL+MACROZ-2) THEN
                    IF((NCODE(5).NE.1).or.(NCODE(6).NE.1)) THEN
                        IOF=(M-1)*(LXNK)
                        FUNKNO(LFLX+IOF+1:LFLX+IOF+LXNK,IG)=RESHAPE(REAL(TMPXNK(:NMZ,:ISPLH,:ISPLH,:3,:NHEX,IPQD,IG)),(/LXNK/))
                    ENDIF
                ENDIF

            ENDDO ! END OF DIRECTION LOOP
            ENDDO ! END OF ENERGY LOOP
            ENDDO ! END OF MACROCELL (WITHIN ONE WAVEFRONT) LOOP
            !$OMP END DO
            DEALLOCATE(JJJ,III,KKK)
        ENDDO ! END OF WAVEFRONT LOOP
        !$OMP END PARALLEL
    ENDDO ! END OF DODECANTS LOOP

    ! SAVE FLUX INFORMATION
    DO IG=1,NGEFF
        IF(.NOT.INCONV(IG)) CYCLE
        FUNKNO(:LFLX,IG)=RESHAPE(REAL(FLUX(:NM,:NSCT,:3*ISPLH**2,:NHEX,:LZ,IG)),(/LFLX/))
    ENDDO

    DEALLOCATE(FLUX,TMPXNI,TMPXNJ,TMPXND,TMPXNK,TMPMAT,IHMAP)
    RETURN
    500 FORMAT(16H SNFH3D: thread=,I8,12H --->(group=,I4,7H angle=,I4,1H))
  END SUBROUTINE SNFH3D
END MODULE SNFHD_MOD
