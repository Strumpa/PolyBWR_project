!
!---------------------------------------------------------------------
!
!Purpose:
! To analyze a geometry made of surfacic element using the SALT
! tracking procedure.
!
!Copyright:
! Copyright (C) 2014 Ecole Polytechnique de Montreal.
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
! ITRACK  pointer to the TRACKING data structure in creation mode.
! NBSLIN  maximum number of segments in a single tracking line.
! RCUTOF  minimum distance between two surfacic elements.
! IPRINT  print level.
!
!Parameters: output
! GG      geometry basic information.
!
!---------------------------------------------------------------------
!
SUBROUTINE SALACG(FGEO ,ITRACK, NBSLIN, RCUTOF, IPRINT, GG)
  USE GANLIB
  USE PRECISION_AND_KINDS, ONLY : PDB
  USE SAL_GEOMETRY_TYPES,  ONLY : T_G_BASIC,TYPGEO,NBFOLD,NBMED,F_GEO,ISPEC,NANIS, &
                                  & LGSPEC,LMERGM
  USE SAL_TRACKING_TYPES,  ONLY : NMAX2,PRTIND,ITRAC2,RTRAC2,IPART,RPART,NIPART, &
                                  & NRPART,EPS1
  USE SAL_GEOMETRY_MOD,    ONLY : SAL100
  IMPLICIT NONE
  !----
  !  Subroutine arguments
  !----
  TYPE(C_PTR) ITRACK
  INTEGER  FGEO,NBSLIN,IPRINT
  REAL(PDB) RCUTOF
  TYPE(T_G_BASIC) :: GG
  !----
  !  Local variables
  !----
  INTEGER, PARAMETER :: NSTATE=40
  INTEGER, PARAMETER ::  NDIM=2 ! NUMBER OF DIMENSIONS
  INTEGER, PARAMETER ::  NALBG=6 ! NUMBER OF ALBEDOS
  LOGICAL LGINF
  INTEGER, DIMENSION(NSTATE) :: I_STATE,IEDIMG
  INTEGER OK,I,J,NREG,ELEM,NFREG,LEAK,NSOUT,ICODE(NALBG),INDEX,MMAX
  REAL GALBED(NALBG)
  CHARACTER(LEN=72) TEXT72
  REAL(PDB) :: DGMESHX(2),DGMESHY(2)
  INTEGER :: NMAX0 = 100000
  !----
  !  Allocatable arrays
  !----
  INTEGER, DIMENSION(:) , ALLOCATABLE :: ITAB ! LOCAL ARRAY
  REAL, DIMENSION(:), ALLOCATABLE :: VOLUME ! LOCAL VOLUME SINGLE PRECISION
  INTEGER, ALLOCATABLE, DIMENSION(:) :: MATALB,KEYMRG,IBC
  REAL(PDB), ALLOCATABLE, DIMENSION(:) :: VOLSUR
  !----
  !  Data statements
  !----
  INTEGER   SDIRE(6)
  DATA      SDIRE/ 5,6,3,4,1,2 /
  !----
  !  Recover options from state vector
  !----
  CALL LCMGET(ITRACK,'STATE-VECTOR',I_STATE) 
  NANIS=I_STATE(6)
  ISPEC=I_STATE(9)
  IF(IPRINT>0) THEN
     IF(ISPEC==0) THEN
        WRITE(6,*) 'SALACG: isotropic boundary conditions'
     ELSE IF(ISPEC==1) THEN
        WRITE(6,*) 'SALACG: specular boundary conditions'
     ENDIF
  ENDIF
  !----
  !  Read the surfacic file and fill GG object
  !----
  PRTIND=IPRINT
  F_GEO=FGEO
  EPS1=1.E-5_PDB
  IF(RCUTOF>0._PDB) THEN
    EPS1=RCUTOF
    IF(PRTIND>0) WRITE(*,*) "SALACG: set eps1 to ",EPS1
  ENDIF
  !------------
  call SAL100(GG)
  !------------
  IF(IPRINT > 0) WRITE(6,'(/" SALACG: TYPGEO=",I5," NBFOLD=",I5)') TYPGEO,NBFOLD
  !----
  ! Perform optional MERGE MIX
  !----
  IF(LMERGM) THEN
    GG%NUM_MERGE(:)=GG%MED(:)
    MMAX=MAXVAL(GG%NUM_MERGE(:))
    DO I=1,MMAX
      10 IF(I.GT.MAXVAL(GG%NUM_MERGE(:))) EXIT
      DO J=1,GG%NB_NODE
        IF(GG%NUM_MERGE(J).EQ.I) GO TO 20
      ENDDO
      DO J=1,GG%NB_NODE
        IF(GG%NUM_MERGE(J).GE.I) GG%NUM_MERGE(J)=GG%NUM_MERGE(J)-1
      ENDDO
      GO TO 10
      20 CONTINUE
    ENDDO
  ENDIF
  !----
  ! Store GG object in geometry directory on LCM
  !----
  CALL LCMSIX(ITRACK,'GEOMETRY    ',1)
  CALL LCMPUT(ITRACK,'NB_ELEM     ',1,1,GG%NB_ELEM)
  CALL LCMPUT(ITRACK,'NIPAR       ',1,1,SIZE(GG%IPAR,1))
  CALL LCMPUT(ITRACK,'IPAR        ',SIZE(GG%IPAR),1,GG%IPAR)
  CALL LCMPUT(ITRACK,'RPAR        ',SIZE(GG%RPAR),4,GG%RPAR)
  CALL LCMPUT(ITRACK,'ISURF2_ELEM ',SIZE(GG%ISURF2_ELEM),1,GG%ISURF2_ELEM)
  CALL LCMPUT(ITRACK,'NB_NODE     ',1,1,GG%NB_NODE)
  CALL LCMPUT(ITRACK,'VOL_NODE    ',GG%NB_NODE,4,GG%VOL_NODE)
  CALL LCMPUT(ITRACK,'NB_SURF2    ',1,1,GG%NB_SURF2)
  LGINF = .TRUE.
  IF(GG%NB_SURF2 > 0) THEN
     CALL LCMPUT(ITRACK,'IBC2_SURF2  ',SIZE(GG%IBC2_SURF2),1,GG%IBC2_SURF2)
     CALL LCMPUT(ITRACK,'IELEM_SURF2 ',SIZE(GG%IELEM_SURF2),1,GG%IELEM_SURF2)
     CALL LCMPUT(ITRACK,'SURF2       ',SIZE(GG%SURF2),4,GG%SURF2)
     DO I = 1,GG%NB_SURF2
       LGINF = LGINF .AND. (GG%BCDATA(6,GG%IDATA_BC2(GG%IBC2_SURF2(I))) == 1._PDB)
     ENDDO
  ENDIF
  CALL LCMPUT(ITRACK,'NPERIM_MAC2 ',1,1,GG%NPERIM_MAC2)
  CALL LCMPUT(ITRACK,'PERIM_MAC2  ',SIZE(GG%PERIM_MAC2),1,GG%PERIM_MAC2)
  CALL LCMPUT(ITRACK,'PPERIM_MAC2 ',SIZE(GG%PPERIM_MAC2),1,GG%PPERIM_MAC2)
  CALL LCMPUT(ITRACK,'PERIM_NODE  ',SIZE(GG%PERIM_NODE),1,GG%PERIM_NODE)
  CALL LCMPUT(ITRACK,'PPERIM_NODE ',SIZE(GG%PPERIM_NODE),1,GG%PPERIM_NODE)
  CALL LCMPUT(ITRACK,'BC_DATA_DIM2',1,1,SIZE(GG%BCDATA,2))
  CALL LCMPUT(ITRACK,'BC_DATA     ',SIZE(GG%BCDATA),4,GG%BCDATA)
  CALL LCMPUT(ITRACK,'NB_BC2      ',1,1,GG%NB_BC2)
  CALL LCMPUT(ITRACK,'TYPE_BC2    ',SIZE(GG%TYPE_BC2),1,GG%TYPE_BC2)
  CALL LCMPUT(ITRACK,'IDATA_BC2   ',SIZE(GG%IDATA_BC2),1,GG%IDATA_BC2)
  CALL LCMSIX(ITRACK,' ',2) ! come back to father directory
  !----
  ! Print tracking object directory
  !----
  IF(IPRINT > 1) THEN
    CALL LCMLIB(ITRACK)
    CALL LCMSIX(ITRACK,'GEOMETRY',1)
    CALL LCMLIB(ITRACK)
    CALL LCMSIX(ITRACK,' ',2)
  ENDIF
  !----
  ! store the STATE VECTOR
  !----
  NREG=MAXVAL(GG%NUM_MERGE)
  LEAK=1
  IF(.NOT.LGINF) LEAK=0 ! reset the leakage flag
  I_STATE(1) = NREG ! number of regions
  I_STATE(2) = NREG ! number of unknowns in DRAGON
  I_STATE(3) = LEAK ! 1 = absent leakage, 0 leakage
  I_STATE(4) = NBMED ! maximum number of mixture
  IF(ISPEC == 0) THEN
    I_STATE(5)=GG%NB_SURF2 ! number of outer surface
    NSOUT=GG%NB_SURF2
  ELSE IF((TYPGEO == 7).OR.(TYPGEO == 8).OR.(TYPGEO == 10).OR.(TYPGEO == 12)) THEN
    I_STATE(5)=3
    NSOUT=3
  ELSE IF(TYPGEO == 9) THEN
    I_STATE(5)=6
    NSOUT=6
  ELSE
    I_STATE(5)=4
    NSOUT=4
  ENDIF
  CALL LCMPUT(ITRACK,'STATE-VECTOR',NSTATE,1,I_STATE)   
  !
  ! fill-in medium number per region
  ALLOCATE(ITAB(NREG),VOLUME(NREG), STAT =OK)
  IF(OK /= 0) CALL XABORT('SALACG: failure to allocate integer ITAB')
  ! fill in MATCOD
  DO J=1,GG%NB_NODE
    ITAB(GG%NUM_MERGE(J)) = GG%MED(J)
  ENDDO
  CALL LCMPUT(ITRACK,'MATCOD',NREG,1,ITAB(1:NREG) ) 
  ! fill-in KEYFLX per region
  DO I=1,NREG
    ITAB(I) = I
  ENDDO
  CALL LCMPUT(ITRACK,'MERGE',NREG,1,ITAB)
  CALL LCMPUT(ITRACK,'KEYFLX',NREG,1,ITAB)
  ! fill-in volumes per region
  VOLUME(:NREG) =0.
  DO I=1,GG%NB_NODE
    VOLUME(GG%NUM_MERGE(I)) = VOLUME(GG%NUM_MERGE(I)) + REAL(GG%VOL_NODE(I))
  ENDDO
  CALL LCMPUT(ITRACK,'VOLUME',NREG,2,VOLUME)
  DEALLOCATE(VOLUME,ITAB)

  ! useful values in SAL_TRACKING_TYPES module
  NFREG=GG%NB_NODE
  CALL LCMSIX(ITRACK,'NXTRecords',1)
    DGMESHX=(/ 1.E10_PDB , -1.E10_PDB /)
    DGMESHY=(/ 1.E10_PDB , -1.E10_PDB /)
    DO ELEM=1,GG%NB_ELEM
      DGMESHX(1)=MIN(DGMESHX(1),GG%RPAR(1,ELEM))
      DGMESHX(2)=MAX(DGMESHX(2),GG%RPAR(1,ELEM))
      DGMESHY(1)=MIN(DGMESHY(1),GG%RPAR(2,ELEM))
      DGMESHY(2)=MAX(DGMESHY(2),GG%RPAR(2,ELEM))
    ENDDO
    CALL LCMPUT(ITRACK,'G00000001SMX',2,4,DGMESHX)
    CALL LCMPUT(ITRACK,'G00000001SMY',2,4,DGMESHY)
    IEDIMG(:NSTATE)=0
    IEDIMG(1)=NDIM
    IEDIMG(2)=0 ! Cartesian geometry
    IF(TYPGEO.EQ.8) IEDIMG(2)=2 ! Isocel geometry with specular reflection
    IF(TYPGEO.EQ.9) IEDIMG(2)=3 ! Hexagonal geometry with translation
    IF(TYPGEO.EQ.10) IEDIMG(2)=4 ! Isocel geometry with RA60 symmetry
    IF(TYPGEO.EQ.11) IEDIMG(2)=5 ! Lozenge geometry with R120 rotation
    IF(TYPGEO.EQ.12) IEDIMG(2)=6 ! S30 geometry with specular reflection
    IEDIMG(5)=1 ! 1 cellule
    IEDIMG(13)=1 ! 1 cellule
    IEDIMG(14)=1 ! 1 cellule
    IEDIMG(22)=NSOUT ! number of external surfaces for this geometry
    IEDIMG(23)=NFREG ! number of regions for this geometry
    IEDIMG(25)=GG%NB_NODE
    CALL LCMPUT(ITRACK,'G00000001DIM',NSTATE,1,IEDIMG)
  CALL LCMSIX(ITRACK,' ',2)  ! come back to father directory
  !----
  ! process boundary conditions
  !----
  IF(LGSPEC) THEN
    IF(ISPEC/=1) CALL XABORT('SALACG: the surfacic file can only be used with' &
    //' cyclic tracking')
  ENDIF
  IF(IPRINT>0) WRITE(6,*) 'number of regions,surfaces,nodes',NREG,NSOUT,NFREG
  ICODE(:NALBG)=0
  ICODE(1)=-1
  ICODE(2)=-2
  ICODE(3)=-3
  ICODE(4)=-4
  ALLOCATE(MATALB(NFREG+NSOUT+1),VOLSUR(NFREG+NSOUT+1),KEYMRG(NFREG+NSOUT+1))
  CALL LCMGET(ITRACK,'MATCOD',MATALB(NSOUT+2))
  ALLOCATE(VOLUME(NREG))
  CALL LCMGET(ITRACK,'VOLUME',VOLUME)
  DO I=1,NREG
    VOLSUR(NSOUT+1+I)=VOLUME(I)
  ENDDO
  DEALLOCATE(VOLUME)
  ! boundary conditions structures
  GALBED(:NALBG)=REAL(GG%ALBEDO)
  IF(ISPEC == 0) THEN
    IF(GG%NALBG > 8) CALL XABORT('SALACG: Albedo array overflow(1).')
    DO I=1,NSOUT
      VOLSUR(I)=GG%SURF2(I)
      INDEX=GG%IDATA_BC2(GG%IBC2_SURF2(I))
      IF((INDEX <= 0).OR.(INDEX > GG%NALBG)) THEN
        CALL XABORT('SALACG: Albedo array overflow(2).')
      ENDIF
      KEYMRG(NSOUT-I+1)=-I
      MATALB(NSOUT-I+1)=-SDIRE(INDEX)
      GALBED(SDIRE(INDEX))=REAL(GG%BCDATA(6,INDEX))
    ENDDO
  ELSE
    DO I=1,NSOUT
      VOLSUR(I)=0.0
      KEYMRG(NSOUT-I+1)=-I
      MATALB(NSOUT-I+1)=-I
      GALBED(I)=1.0
    ENDDO
  ENDIF
  MATALB(NSOUT+1)=0
  KEYMRG(NSOUT+1)=0
  VOLSUR(NSOUT+1)=0._PDB
  DO I=1,NREG
    KEYMRG(NSOUT+1+I)=I
  ENDDO
  !
  IF(IPRINT>1) THEN
     CALL PRINDM('VOLUME',VOLSUR(1),NREG+NSOUT+1)
     CALL PRINIM('MATALB',MATALB(1),NREG+NSOUT+1)
     CALL PRINIM('KEYMRG',KEYMRG(1),NREG+NSOUT+1)
  ENDIF
  IF(IPRINT>0) THEN
     CALL PRINIM('ICODE ',ICODE(1),NALBG)
     CALL PRINAM('GALBED',GALBED(1),NALBG)
  ENDIF
  !----
  ! fill in tracking LCM object in excelt format
  !----
  TEXT72='SAL TRACKING'
  CALL LCMPTC(ITRACK,'TITLE',72,TEXT72)
  CALL LCMPUT(ITRACK,'ICODE',NALBG,1,ICODE)
  CALL LCMSIX(ITRACK,'NXTRecords',1)
  CALL LCMPUT(ITRACK,'SAreaRvolume',(NREG+NSOUT+1),4,VOLSUR(1))
  CALL LCMPUT(ITRACK,'MATALB',NREG+NSOUT+1,1,MATALB(1))
  CALL LCMPUT(ITRACK,'KEYMRG',NREG+NSOUT+1,1,KEYMRG(1))
  CALL LCMSIX(ITRACK,' ',2)
  IF(NSOUT>0) THEN
     ALLOCATE(IBC(NSOUT))
     DO I=1,NSOUT
        IBC(I)=I
     ENDDO
     CALL LCMPUT(ITRACK,'BC-REFL+TRAN',NSOUT,1,IBC)
     DEALLOCATE(IBC)
  ENDIF
  CALL LCMPUT(ITRACK,'MATCOD',NREG,1,MATALB(NSOUT+2))
  CALL LCMPUT(ITRACK,'ALBEDO',NALBG,2,GALBED)
  DEALLOCATE(KEYMRG,VOLSUR,MATALB)
  !----
  ! allocate memory to hold tracking data
  !----
  !*      tracking data buffer:
  !       integers
  !       ITRAC2(NMAX OR 2*NMAX) = integer tracking array
  !
  !       *integer descriptors in itrac2:
  !           1 = address of last data
  !           2 = total number of sub-trajectories
  !           3 = 
  !           4 = phi for trajectory (2D)
  !
  !       reals
  !       RTRAC2(nmax or 2*nmax) = real tracking array
  !
  !       *real descriptors:
  !           1 =
  !           2 = cos phi entering basic
  !           3 = sin phi left surface
  !           4 = sin phi right surface
  !           5 = cos phi left surface
  !           6 = cos phi right surface
  !           7 = total weight (DELR*WPHI)
  !           8 = radial weight (DELR)
  !
  !       IPART(NIPART,MXELEM)   = to store integer intersection data
  !       RPART(NRPART,MXELEM)   = to store real intersection data
  !
  IF(NBSLIN <= 0)THEN
     NMAX2=NMAX0
  ELSE
     NMAX2=NBSLIN
  ENDIF
  IF(ISPEC == 1) NMAX2=NMAX2*100
  ALLOCATE(ITRAC2(2*NMAX2),IPART(NIPART,GG%NB_ELEM),RTRAC2(NMAX2),RPART(NRPART,GG%NB_ELEM),STAT=OK)
  IF(OK/=0) CALL XABORT('SALACG: not enough memory IRD')
  RETURN
END SUBROUTINE SALACG
