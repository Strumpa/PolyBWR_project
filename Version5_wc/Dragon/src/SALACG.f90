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
! RCUTOF  minimum distance between two surfacic elements.
! IPRINT  print level.
!
!Parameters: output
! GG      geometry basic information.
!
!---------------------------------------------------------------------
!
SUBROUTINE SALACG(FGEO ,ITRACK, RCUTOF, IPRINT, GG)
  USE GANLIB
  USE PRECISION_AND_KINDS, ONLY : PDB
  USE SAL_GEOMETRY_TYPES,  ONLY : T_G_BASIC,TYPGEO,NBFOLD,NBMED,F_GEO,ISPEC, &
                                  & LGSPEC,LMERGM
  USE SAL_TRACKING_TYPES,  ONLY : PRTIND,EPS1
  USE SAL_GEOMETRY_MOD,    ONLY : SAL100
  IMPLICIT NONE
  !----
  !  Subroutine arguments
  !----
  TYPE(C_PTR) ITRACK
  INTEGER  FGEO,IPRINT
  REAL(PDB) RCUTOF
  TYPE(T_G_BASIC) :: GG
  !----
  !  Local variables
  !----
  INTEGER, PARAMETER :: NSTATE=40
  INTEGER, PARAMETER :: NDIM=2 ! NUMBER OF DIMENSIONS
  LOGICAL LGINF
  INTEGER, DIMENSION(NSTATE) :: I_STATE,IEDIMG
  INTEGER NALBG,OK,I,J,NREG,ELEM,NFREG,LEAK,NSOUT,INDEX,MMAX,MAXCDA
  CHARACTER(LEN=72) TEXT72
  REAL(PDB) :: DGMESHX(2),DGMESHY(2)
  !----
  !  Allocatable arrays
  !----
  INTEGER, DIMENSION(:) , ALLOCATABLE :: ITAB,ICODE
  REAL, DIMENSION(:), ALLOCATABLE :: VOLUME,GALBED
  INTEGER, ALLOCATABLE, DIMENSION(:) :: MATALB,KEYMRG,IBC
  REAL(PDB), ALLOCATABLE, DIMENSION(:) :: VOLSUR
  !----
  !  Recover options from state vector
  !----
  CALL LCMGET(ITRACK,'STATE-VECTOR',I_STATE) 
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
  IF(GG%NBBCDA.GT.0) THEN
    LGINF = .TRUE.
    DO I=1, GG%NBBCDA
      LGINF = LGINF .AND. (GG%BCDATAREAD(I)%BCDATA(6) == 1._PDB)
    ENDDO
  ELSE
    LGINF = (GG%ALBEDO == 1._PDB)
  ENDIF
  IF(GG%NB_SURF2 > 0) THEN
     CALL LCMPUT(ITRACK,'IBC2_SURF2  ',SIZE(GG%IBC2_SURF2),1,GG%IBC2_SURF2)
     CALL LCMPUT(ITRACK,'IELEM_SURF2 ',SIZE(GG%IELEM_SURF2),1,GG%IELEM_SURF2)
     CALL LCMPUT(ITRACK,'SURF2       ',SIZE(GG%SURF2),4,GG%SURF2)
  ENDIF
  CALL LCMPUT(ITRACK,'NPERIM_MAC2 ',1,1,GG%NPERIM_MAC2)
  CALL LCMPUT(ITRACK,'PERIM_MAC2  ',SIZE(GG%PERIM_MAC2),1,GG%PERIM_MAC2)
  CALL LCMPUT(ITRACK,'PPERIM_MAC2 ',SIZE(GG%PPERIM_MAC2),1,GG%PPERIM_MAC2)
  CALL LCMPUT(ITRACK,'PERIM_NODE  ',SIZE(GG%PERIM_NODE),1,GG%PERIM_NODE)
  CALL LCMPUT(ITRACK,'PPERIM_NODE ',SIZE(GG%PPERIM_NODE),1,GG%PPERIM_NODE)
  CALL LCMPUT(ITRACK,'BC_DATA_DIM2',1,1,SIZE(GG%BCDATA,2))
  IF(SIZE(GG%BCDATA) > 0) THEN
    CALL LCMPUT(ITRACK,'BC_DATA     ',SIZE(GG%BCDATA),4,GG%BCDATA)
  ENDIF
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
  IF(IPRINT>0) WRITE(6,*) 'number of merged regions,surfaces,nodes',NREG,NSOUT,NFREG
  ALLOCATE(MATALB(-NSOUT:NFREG),VOLSUR(-NSOUT:NFREG),KEYMRG(-NSOUT:NFREG))
  CALL LCMGET(ITRACK,'MATCOD',MATALB(1))
  ALLOCATE(VOLUME(NREG))
  CALL LCMGET(ITRACK,'VOLUME',VOLUME)
  VOLSUR(1:NREG)=VOLUME(:NREG)
  DEALLOCATE(VOLUME)
  ! boundary conditions structures
  MAXCDA=MAX(6,GG%NALBG)
  ALLOCATE(ICODE(MAXCDA),GALBED(MAXCDA))
  ICODE(:MAXCDA)=(/ (-I,I=1,MAXCDA) /)
  GALBED(:MAXCDA)=REAL(GG%ALBEDO)
  IF(ISPEC == 0) THEN
    NALBG=GG%NALBG
    IF(TYPGEO.EQ.0) NALBG=6
    DO I=1,NSOUT
      KEYMRG(-I)=-I
      VOLSUR(-I)=GG%SURF2(I)
      INDEX=GG%IDATA_BC2(GG%IBC2_SURF2(I))
      IF(INDEX.EQ.0) THEN
        ! Use the default albedo
        MATALB(-I)=-1
        GALBED(1)=REAL(GG%ALBEDO)
      ELSE
        IF(INDEX.GT.MAXCDA) CALL XABORT('SALACG: INDEX overflow(1).')
        IF(INDEX.GT.GG%NALBG) CALL XABORT('SALACG: INDEX overflow(2).')
        MATALB(-I)=-INDEX
        IF(SIZE(GG%BCDATA) > 0) THEN
          GALBED(INDEX)=REAL(GG%BCDATA(6,INDEX))
        ELSE
          GALBED(INDEX)=REAL(GG%ALBEDO)
        ENDIF
      ENDIF
    ENDDO
  ELSE
    NALBG=6
    DO I=1,NSOUT
      VOLSUR(-I)=0.0
      KEYMRG(-I)=-I
      MATALB(-I)=-1
    ENDDO
    GALBED(:NALBG)=1.0
  ENDIF
  MATALB(0)=0
  KEYMRG(0)=0
  VOLSUR(0)=0._PDB
  DO I=1,NREG
    KEYMRG(I)=I
  ENDDO
  !
  IF(IPRINT>1) THEN
     CALL PRINDM('VOLUME',VOLSUR(-NSOUT),NREG+NSOUT+1)
     CALL PRINIM('MATALB',MATALB(-NSOUT),NREG+NSOUT+1)
     CALL PRINIM('KEYMRG',KEYMRG(-NSOUT),NREG+NSOUT+1)
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
  CALL LCMPUT(ITRACK,'SAreaRvolume',NREG+NSOUT+1,4,VOLSUR(-NSOUT))
  CALL LCMPUT(ITRACK,'MATALB',NREG+NSOUT+1,1,MATALB(-NSOUT))
  CALL LCMPUT(ITRACK,'KEYMRG',NREG+NSOUT+1,1,KEYMRG(-NSOUT))
  CALL LCMSIX(ITRACK,' ',2)
  IF(NSOUT>0) THEN
     ALLOCATE(IBC(NSOUT))
     DO I=1,NSOUT
        IBC(I)=I
     ENDDO
     CALL LCMPUT(ITRACK,'BC-REFL+TRAN',NSOUT,1,IBC)
     DEALLOCATE(IBC)
  ENDIF
  CALL LCMPUT(ITRACK,'MATCOD',NREG,1,MATALB(1))
  CALL LCMPUT(ITRACK,'ALBEDO',NALBG,2,GALBED)
  DEALLOCATE(GALBED,ICODE)
  DEALLOCATE(KEYMRG,VOLSUR,MATALB)
  RETURN
END SUBROUTINE SALACG
