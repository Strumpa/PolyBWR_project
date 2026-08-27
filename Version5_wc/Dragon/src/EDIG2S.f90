!
!-----------------------------------------------------------------------
!
!Purpose:
! Perform an homogenization based on a surfacic file. Process RECT, TRIA,
! QUAD, HEXA and POLY data options
!
!Copyright:
! Copyright (C) 2017 Ecole Polytechnique de Montreal
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version
!
!Author(s): A. Hebert
!
!Parameters: input
! IPRINT  print flag.
! IPGEO2  pointer towards the macrogeometry.
! IFGEO   unit file number of the surfacic file.
! NREG    number of regions.
!
!Parameters: input
! NMERGE  number of merged indices in array IMERGE after REMIX.
! IMERGE  merged regions position.
!
!-----------------------------------------------------------------------
!
SUBROUTINE EDIG2S(IPRINT,IPGEO2,IFGEO,NREG,NMERGE,IMERGE)
    USE GANLIB
    USE PRECISION_AND_KINDS, ONLY : PDB
    USE SAL_NUMERIC_MOD, ONLY : DET_ROSETTA
    USE SALGET_FUNS_MOD
    !----
    !  Subroutine arguments
    !----
    TYPE(C_PTR) IPGEO2
    INTEGER IPRINT,IFGEO,NREG,NMERGE,IMERGE(NREG)
    !----
    !  Local variables
    !----
    PARAMETER(IFOUT0=0,NSTATE=40)
    INTEGER PREC,DATAIN(25),IPAR(5),ISTATE(NSTATE)
    REAL DATARE(25),DIR1(2),DIR2(2)
    REAL(PDB) CX,CY,DX,DY,SAA,SAB,ANGL,RPAR(5),ALIGN(3,3),DET1,DET2
    REAL(PDB) NODX1,NODX2,NODY1,NODY2
    REAL(PDB), PARAMETER :: CONV=3.141592654_PDB/180.0_PDB
    CHARACTER NAME_GEOM*12,HSMG*131
    LOGICAL LTEST
    !----
    !  Allocatable arrays
    !----
    INTEGER, DIMENSION(:), ALLOCATABLE :: NUM_MERGE,IFLUX,ITNODE,IREMIX
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: ICOUNT
    REAL(PDB), DIMENSION(:), ALLOCATABLE :: SURF,TLAMB1,TLAMB2
    REAL, DIMENSION(:,:), ALLOCATABLE :: NODE
    !----
    !  Recover macrogeometry
    !----
    IF(.NOT.C_ASSOCIATED(IPGEO2)) THEN
      CALL XABORT('EDIG2S: MACROGEOMETRY IS NOT DEFINED.')
    ENDIF
    CALL LCMGET(IPGEO2,'STATE-VECTOR',ISTATE)
    IF(ISTATE(1).NE.31) CALL XABORT('EDIG2S: INVALID TYPE OF MACROGEOMETRY.')
    CALL LCMSIX(IPGEO2,'MACRO-GEOM',1)
      CALL LCMGET(IPGEO2,'STATE-VECTOR',ISTATE)
      IDIM=ISTATE(1)
      NZS=ISTATE(2)
      NMERGE=ISTATE(4)
      ALLOCATE(NODE(IDIM,NZS),ITNODE(NZS))
      CALL LCMGET(IPGEO2,'NODE',NODE)
      CALL LCMGET(IPGEO2,'ITNODE',ITNODE)
      IF(IPRINT.GT.0) WRITE(6,200) ISTATE(:5)
      DO IM=1,NZS
        IF(ITNODE(IM).EQ.1) THEN
          IF(IDIM.LT.4) CALL XABORT('EDIG2S: IDIM OVERFLOW(1).')
        ELSE IF(ITNODE(IM).EQ.2) THEN
          CALL XABORT('EDIG2S: INVALID ITNODE VALUE.')
        ELSE IF(ITNODE(IM).GE.3) THEN
          IF(IDIM.LT.2*ITNODE(IM)) CALL XABORT('EDIG2S: IDIM OVERFLOW(2).')
        ENDIF
      ENDDO
      IF(ISTATE(5).EQ.3) CALL XABORT('EDIG2S: 3D NOT SUPPORTED.')
    CALL LCMSIX(IPGEO2,' ',2)
    !----
    !  Determine homogenization indices
    !----
    IF(IFGEO.EQ.0) CALL XABORT('EDIG2S: surfacic file not defined.')
    CALL SALGET(DATAIN,6,IFGEO,IFOUT0,'dimensions for geometry')
    NBNODE=DATAIN(3)
    NBELEM=DATAIN(4)
    NBFLUX=DATAIN(6)
    CALL SALGET(DATAIN,3,IFGEO,IFOUT0,'index kndex prec')
    INDEX=DATAIN(1)
    KNDEX=DATAIN(2)
    PREC=DATAIN(3)
    CALL SALGET(DATARE,1,IFGEO,IFOUT0,'eps')
    EPS=DATARE(1)
    ALLOCATE(NUM_MERGE(NBNODE))
    CALL SALGET(NUM_MERGE,NBNODE,IFGEO,IFOUT0,'FLUX INDEX PER NODE')
    IF(MAXVAL(NUM_MERGE).NE.NBFLUX) CALL XABORT('EDIG2S: inconsistent NBFLUX.')
    CALL SALGET(NAME_GEOM,IFGEO,IFOUT0,'NAMES OF MACROS')
    ALLOCATE(IFLUX(NBFLUX))
    CALL SALGET(IFLUX,NBFLUX,IFGEO,IFOUT0,'macro order number per flux region.')
    DEALLOCATE(IFLUX)
    ALLOCATE(ICOUNT(NBNODE,NZS))
    ICOUNT(:NBNODE,:NZS)=0
    DO IELEM=1,NBELEM
      IPAR(:)=0
      RPAR(:)=0.0
      CALL SALGET(IPAR,3,IFGEO,IFOUT0,'integer descriptors')
      ITYPE=IPAR(1)
      SELECT CASE (ITYPE)
        CASE (1)
        NBER=4
        CASE (2)
        NBER=3
        CASE (3)
        NBER=5
      CASE DEFAULT
        WRITE(6,'(1X,''==> SAL126: unknown type '',I3)') ITYPE
        CALL XABORT('EDIG2S: unknown element type.')
      END SELECT
      CALL SALGET(RPAR,NBER,IFGEO,IFOUT0,PREC,'real descriptors')
      IF(ITYPE.EQ.1) THEN
        CX=RPAR(1) ; CY=RPAR(2)
        DX=CX+RPAR(3) ; DY=CY+RPAR(4)
        LOOP1: DO IM=1,NZS
          IF(ITNODE(IM).EQ.1) THEN
            NODX1=NODE(1,IM) ; NODX2=NODE(2,IM)
            NODY1=NODE(3,IM) ; NODY2=NODE(4,IM)
            IF((MIN(CX,DX).GE.NODX1-EPS).AND.(MAX(CX,DX).LE.NODX2+EPS).AND. &
               (MIN(CY,DY).GE.NODY1-EPS).AND.(MAX(CY,DY).LE.NODY2+EPS)) THEN
              IF((ABS(CX-DX).LE.EPS).AND.(ABS(CX-NODX1).LE.EPS)) THEN ! left vertical side
                IF(DY.GT.CY) THEN
                  INODE=IPAR(2)
                ELSE
                  INODE=IPAR(3)
                ENDIF
                IF(INODE.EQ.0) THEN
                  WRITE(HSMG,'(38HEDIG2S: IPAR inconcistency for element,i10,4H(1).)') IELEM
                  CALL XABORT(HSMG)
                ENDIF
                ICOUNT(INODE,IM)=ICOUNT(INODE,IM)+1
              ELSE IF((ABS(CX-DX).LE.EPS).AND.(ABS(CX-NODX2).LE.EPS)) THEN ! right vertical side
                IF(DY.GT.CY) THEN
                  INODE=IPAR(3)
                ELSE
                  INODE=IPAR(2)
                ENDIF
                IF(INODE.EQ.0) THEN
                  WRITE(HSMG,'(38HEDIG2S: IPAR inconcistency for element,i10,4H(2).)') IELEM
                  CALL XABORT(HSMG)
                ENDIF
                ICOUNT(INODE,IM)=ICOUNT(INODE,IM)+1
              ELSE IF((ABS(CY-DY).LE.EPS).AND.(ABS(CY-NODY1).LE.EPS)) THEN ! lower horizontal side
                IF(DX.GT.CX) THEN
                  INODE=IPAR(3)
                ELSE
                  INODE=IPAR(2)
                ENDIF
                IF(INODE.EQ.0) THEN
                  WRITE(HSMG,'(38HEDIG2S: IPAR inconcistency for element,i10,4H(3).)') IELEM
                  CALL XABORT(HSMG)
                ENDIF
                ICOUNT(INODE,IM)=ICOUNT(INODE,IM)+1
              ELSE IF((ABS(CY-DY).LE.EPS).AND.(ABS(CY-NODY2).LE.EPS)) THEN ! upper horizontal side
                IF(DX.GT.CX) THEN
                  INODE=IPAR(2)
                ELSE
                  INODE=IPAR(3)
                ENDIF
                IF(INODE.EQ.0) THEN
                  WRITE(HSMG,'(38HEDIG2S: IPAR inconcistency for element,i10,4H(4).)') IELEM
                  CALL XABORT(HSMG)
                ENDIF
                ICOUNT(INODE,IM)=ICOUNT(INODE,IM)+1
              ELSE IF((ABS(CX-DX).LE.EPS).OR.(ABS(CY-DY).LE.EPS)) THEN
                IF(IPAR(2).GT.0) ICOUNT(IPAR(2),IM)=ICOUNT(IPAR(2),IM)+1
                IF(IPAR(3).GT.0) ICOUNT(IPAR(3),IM)=ICOUNT(IPAR(3),IM)+1
              ENDIF
            ENDIF
          ELSE
            ALLOCATE(TLAMB1(ITNODE(IM)),TLAMB2(ITNODE(IM)))
            TLAMB1=EDIBAR(ITNODE(IM),NODE(1,IM),CX,CY)
            TLAMB2=EDIBAR(ITNODE(IM),NODE(1,IM),DX,DY)
            LTEST=.TRUE.
            DO I=1,ITNODE(IM)
              LTEST=LTEST.AND.(TLAMB1(I).GE.-EPS).AND.(TLAMB2(I).GE.-EPS)
            ENDDO
            DEALLOCATE(TLAMB2,TLAMB1)
            IF(LTEST) THEN
              DO I=1,ITNODE(IM)
                IP1 = I+1; IF (IP1 > ITNODE(IM)) IP1 = 1 ! wrap around to first vertex
                ! Find if CXY and DXY are lying between vertices I and IP1
                ALIGN(:3,3)=1.0_PDB
                ALIGN(1,1)=NODE(2*I-1,IM) ; ALIGN(1,2)=NODE(2*I,IM)
                ALIGN(2,1)=NODE(2*IP1-1,IM) ; ALIGN(2,2)=NODE(2*IP1,IM)
                ALIGN(3,1)=CX; ALIGN(3,2)=CY;
                DET1 = DET_ROSETTA(ALIGN,3)
                ALIGN(3,1)=DX; ALIGN(3,2)=DY;
                DET2 = DET_ROSETTA(ALIGN,3)
                IF((ABS(DET1).LE.EPS).AND.(ABS(DET2).LE.EPS)) THEN
                  DIR1(:2)=EDIDIR(NODE(2*I-1,IM),NODE(2*I,IM),NODE(2*IP1-1,IM),NODE(2*IP1,IM))
                  DIR2(:2)=EDIDIR(REAL(CX),REAL(CY),REAL(DX),REAL(DY))
                  IF(DOT_PRODUCT(DIR1,DIR2).LT.0.0) THEN
                    IF(IPAR(2).EQ.0) THEN
                      WRITE(HSMG,'(38HEDIG2S: IPAR inconcistency for element,i10,4H(5).)') IELEM
                      CALL XABORT(HSMG)
                    ENDIF
                    ICOUNT(IPAR(2),IM)=ICOUNT(IPAR(2),IM)+1
                    CYCLE LOOP1
                  ELSE
                    IF(IPAR(3).EQ.0) THEN
                      WRITE(HSMG,'(38HEDIG2S: IPAR inconcistency for element,i10,4H(6).)') IELEM
                      CALL XABORT(HSMG)
                    ENDIF
                    ICOUNT(IPAR(3),IM)=ICOUNT(IPAR(3),IM)+1
                    CYCLE LOOP1
                  ENDIF                  
                ENDIF
              ENDDO ! I
              IF(IPAR(2).GT.0) ICOUNT(IPAR(2),IM)=ICOUNT(IPAR(2),IM)+1
              IF(IPAR(3).GT.0) ICOUNT(IPAR(3),IM)=ICOUNT(IPAR(3),IM)+1
            ENDIF
          ENDIF
        ENDDO LOOP1
      ELSE IF(ITYPE.EQ.2) THEN
        CX=RPAR(1) ; CY=RPAR(2)
        DO IM=1,NZS
          IF(ITNODE(IM).EQ.1) THEN
            NODX1=NODE(1,IM) ; NODX2=NODE(2,IM)
            NODY1=NODE(3,IM) ; NODY2=NODE(4,IM)
            IF((CX.GE.NODX1-EPS).AND.(CX.LE.NODX2+EPS).AND. &
               (CY.GE.NODY1-EPS).AND.(CY.LE.NODY2+EPS)) THEN
              IF(IPAR(2).GT.0) ICOUNT(IPAR(2),IM)=ICOUNT(IPAR(2),IM)+1
              IF(IPAR(3).GT.0) ICOUNT(IPAR(3),IM)=ICOUNT(IPAR(3),IM)+1
            ENDIF
          ELSE
            ALLOCATE(TLAMB1(ITNODE(IM)))
            TLAMB1=EDIBAR(ITNODE(IM),NODE(1,IM),CX,CY)
            LTEST=.TRUE.
            DO I=1,ITNODE(IM)
              LTEST=LTEST.AND.(TLAMB1(I).GE.-EPS)
            ENDDO
            DEALLOCATE(TLAMB1)
            IF(LTEST) THEN
              IF(IPAR(2).GT.0) ICOUNT(IPAR(2),IM)=ICOUNT(IPAR(2),IM)+1
              IF(IPAR(3).GT.0) ICOUNT(IPAR(3),IM)=ICOUNT(IPAR(3),IM)+1
            ENDIF
          ENDIF
        ENDDO ! IM
      ELSE IF(ITYPE.EQ.3) THEN
        SAA=RPAR(4) ; SAB=SAA+RPAR(5)
        IF(SAB>SAA) THEN
          ANGL=(SAB+SAA)*0.5
        ELSE
          ANGL=(SAB+SAA)*0.5+180.0
        ENDIF
        CX=RPAR(1)+COS(ANGL*CONV)*RPAR(3) ; CY=RPAR(2)+SIN(ANGL*CONV)*RPAR(3)
        DO IM=1,NZS
          IF(ITNODE(IM).EQ.1) THEN
            NODX1=NODE(1,IM) ; NODX2=NODE(2,IM)
            NODY1=NODE(3,IM) ; NODY2=NODE(4,IM)
            IF((CX.GE.NODX1-EPS).AND.(CX.LE.NODX2+EPS).AND. &
               (CY.GE.NODY1-EPS).AND.(CY.LE.NODY2+EPS)) THEN
              IF(IPAR(2).GT.0) ICOUNT(IPAR(2),IM)=ICOUNT(IPAR(2),IM)+1
              IF(IPAR(3).GT.0) ICOUNT(IPAR(3),IM)=ICOUNT(IPAR(3),IM)+1
            ENDIF
          ELSE
            ALLOCATE(TLAMB1(ITNODE(IM)))
            TLAMB1=EDIBAR(ITNODE(IM),NODE(1,IM),CX,CY)
            LTEST=.TRUE.
            DO I=1,ITNODE(IM)
              LTEST=LTEST.AND.(TLAMB1(I).GE.-EPS)
            ENDDO
            DEALLOCATE(TLAMB1)
            IF(LTEST) THEN
              IF(IPAR(2).GT.0) ICOUNT(IPAR(2),IM)=ICOUNT(IPAR(2),IM)+1
              IF(IPAR(3).GT.0) ICOUNT(IPAR(3),IM)=ICOUNT(IPAR(3),IM)+1
            ENDIF
          ENDIF
        ENDDO ! IM
      ENDIF
    ENDDO ! IELEM
    IMERGE(:NREG)=0
    ITEST=0
    DO IM=1,NZS
      DO INODE=1,NBNODE
        IF(ICOUNT(INODE,IM).GT.0) THEN
          IF(IMERGE(NUM_MERGE(INODE)).NE.0) THEN
            WRITE(HSMG,'(46HEDIG2S: inconsistent homogenization in mixture,I8, &
            & 14h, SURFIL node=,I8,29h belong to many output nodes.)') IM,INODE
            CALL XABORT(HSMG)
          ENDIF
          IMERGE(NUM_MERGE(INODE))=IM
          ITEST=ITEST+1
        ENDIF
      ENDDO
    ENDDO
    IF(IPRINT.GT.0) THEN
      WRITE(6,'(53H EDIG2S: NUMBER OF NODES PROCESSED BY HOMOGENIZATION=,I8/ &
      & 9X,32HNUMBER OF NODES IN THE GEOMETRY=,12X,I8/ &
      & 9X,31HNUMBER OF HOMOGENEOUS MIXTURES=,13X,I8)') ITEST,NBNODE,NZS
      ALLOCATE(SURF(NZS))
      DO IM=1,NZS
        SURF(IM)=0.0D0
        IF(ITNODE(IM).EQ.1) THEN
          SURF(IM)=(NODE(2,IM)-NODE(1,IM))*(NODE(4,IM)-NODE(3,IM))
        ELSE
          SURF(IM)=EDIAREA(ITNODE(IM),NODE(1,IM))
        ENDIF
      ENDDO
      WRITE(6,'(53H EDIG2S: ANALYTICAL OUTPUT NODE SURFACES BEFORE REMIX/(5X,1P,12E12.4))') &
      & SURF(:NZS)
      DEALLOCATE(SURF)
    ENDIF
    DEALLOCATE(NUM_MERGE,ICOUNT,ITNODE,NODE)
    !----
    !  Remix homogenized indices
    !----
    IF(NMERGE.GT.0) THEN
      ALLOCATE(IREMIX(NZS))
      CALL LCMGET(IPGEO2,'IREMIX',IREMIX)
      NMERGE=0
      DO IREG=1,NREG
        IM=IMERGE(IREG)
        IF(IM.GT.0) THEN
          IF(IM.GT.NZS) CALL XABORT('EDIG2S: IMERGE OVERFLOW')
          IMERGE(IREG)=IREMIX(IM)
          NMERGE=MAX(NMERGE,IMERGE(IREG))
        ENDIF
      ENDDO
      DEALLOCATE(IREMIX)
    ELSE
      NMERGE=NZS
    ENDIF
    RETURN
    !
    200 FORMAT(/22H MACROGEOMETRY OPTIONS/1X,21(1H-)/ &
    & 7H IDIM  ,I8,34H   (FIRST DIMENSION OF ARRAY NODE)/ &
    & 7H NZS   ,I8,23H   (NUMBER OF MIXTURES)/ &
    & 7H ISYM  ,I8,39H   (=0/1/2/3: UNDEFINED/NONE/REFL/ROTA)/ &
    & 7H NMERGE,I8,33H   (=0/NUMBER OF MERGED MIXTURES)/ &
    & 7H IDIM  ,I8,16H   (=2/3: 2D/3D))
    !
    CONTAINS
      FUNCTION EDIBAR(N, NODE, CX, CY) RESULT(TLAMB)
        ! Compute the barycentric coordinates of (CX,CY)
        INTEGER, INTENT(IN)   :: N             ! number of vertices
        REAL, INTENT(IN)  :: NODE(2*N)         ! polygon vertex coordinates
        REAL(PDB), INTENT(IN)  :: CX, CY       ! target point coordinates
        REAL(PDB) :: TLAMB(N)                  ! resulting normalized weights

        REAL(PDB), ALLOCATABLE :: R(:), ALPHA(:), W(:)
        REAL(PDB) :: DX, DY, DX_NEXT, DY_NEXT, DOT, CROSS, TOTAL_WEIGHT
        INTEGER  :: I, IM1, IP1

        ALLOCATE(R(N), ALPHA(N), W(N))
        TLAMB = 0.0_PDB
        W = 0.0_PDB

        ! step 1: precompute distances (R) from (CX,CY) to all vertices
        DO I = 1, N
            R(I) = SQRT((NODE(2*I-1) - CX)**2 + (NODE(2*I) - CY)**2)
            
            ! ROBUSTNESS CHECK: IF POINT IS EXACTLY ON CORNER I
            IF (R(I) < 1.0E-12_PDB) THEN
                TLAMB = 0.0_PDB
                TLAMB(I) = 1.0_PDB
                RETURN
            ENDIF
        ENDDO

        ! step 2: compute signed angles (ALPHA) between adjacent vectors using atan2
        DO I = 1, N
            IP1 = I + 1; IF (IP1 > N) IP1 = 1 ! wrap around to first vertex
            
            ! vectors from (CX,CY) to CORNER(I) and CORNER(I+1)
            DX = NODE(2*I-1) - CX
            DY = NODE(2*I) - CY
            DX_NEXT = NODE(2*IP1-1) - CX
            DY_NEXT = NODE(2*IP1) - CY

            ! dot product and 2d cross product determinant
            DOT = DX * DX_NEXT + DY * DY_NEXT
            CROSS = DX * DY_NEXT - DY * DX_NEXT

            ! ATAN2 yields the exact signed angle between vectors in [-pi, pi]
            ALPHA(I) = ATAN2(CROSS, DOT)

            ! robustness check: if point lies precisely on the segment between I and I+1
            IF (ABS(CROSS) < 1.0E-12_PDB .AND. DOT < 0.0_PDB) THEN
                TLAMB = 0.0_PDB
                TLAMB(I) = R(IP1) / (R(I) + R(IP1))
                TLAMB(IP1) = 1.0_PDB - TLAMB(I)
                RETURN
            ENDIF
        ENDDO

        ! step 3: compute unnormalized TLAMB using the tangent of half-angles
        DO I = 1, N
            IM1 = I - 1; IF (IM1 < 1) IM1 = N ! wrap around to last vertex
            
            W(I) = (TAN(ALPHA(IM1) / 2.0_PDB) + TAN(ALPHA(I) / 2.0_PDB)) / R(I)
        ENDDO

        ! step 4: normalize TLAMB to sum up to 1.0
        TOTAL_WEIGHT = SUM(W)
        IF (ABS(TOTAL_WEIGHT) > 1.0E-12_PDB) THEN
            TLAMB = W / TOTAL_WEIGHT
        ELSE
            TLAMB = 0.0_PDB
        ENDIF
      END FUNCTION EDIBAR
      !
      FUNCTION EDIAREA(N,PXY) RESULT(POLY_AREA)
        ! Compute the area of a polygon
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        REAL, DIMENSION(2*N), INTENT(IN) :: PXY
        REAL :: POLY_AREA
        !
        INTEGER :: I, J
        REAL :: SUM_VAL

        SUM_VAL = 0.0
        DO I = 1, N
            ! wrap around to the first vertex when at the last vertex
            J = MOD(I, N) + 1
            SUM_VAL = SUM_VAL + (PXY(2*I-1) * PXY(2*J)) - (PXY(2*J-1) * PXY(I))
        ENDDO
        POLY_AREA = 0.5 * ABS(SUM_VAL) 
      END FUNCTION EDIAREA
      !
      FUNCTION EDIDIR(CX,CY,DX,DY) RESULT(TDIR)
        ! Compute the direction of a segment
        REAL, INTENT(IN)  :: CX,CY,DX,DY  ! end points of a segment
        REAL :: TDIR(2)                   ! resulting direction vector
        REAL :: DSUM
        !
        TDIR(:2)=(/ DX-CX, DY-CY /)
        DSUM=SQRT(TDIR(1)**2+TDIR(2)**2)
        TDIR(:2)=TDIR(:2)/DSUM
      END FUNCTION EDIDIR
END SUBROUTINE EDIG2S
