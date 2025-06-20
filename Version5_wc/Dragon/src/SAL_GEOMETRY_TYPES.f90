!
!---------------------------------------------------------------------
!
!Purpose:
! To store common variables for the geometry analysis in SALT:
! module.
!
!Copyright:
! Copyright (C) 2001 Ecole Polytechnique de Montreal.
!
!Author(s):
! X. Warin
!
!---------------------------------------------------------------------
!
MODULE SAL_GEOMETRY_TYPES
  USE PRECISION_AND_KINDS, ONLY : PDB
  IMPLICIT NONE
  INTEGER :: ITYPE
  INTEGER, PARAMETER :: N_IN=29,N_LG=4,N_RE=6
  INTEGER :: FOUT0,F_GEO
  !
  !       types of elements
  INTEGER, PARAMETER, DIMENSION(1:4) :: G_ELE_TYPE=(/1,2,3,4/)
  ! 1 : segment element
  ! 2 : circle element
  ! 3 : arc of circle element
  !       length of element discription data
  INTEGER, PARAMETER, DIMENSION(1:4) :: G_ELE_LEN=(/5,5,6,6/)
  ! 1 : segment element
  ! 2 : circle element
  ! 3 : arc of circle element
  !       types of boundary conditions
  INTEGER, PARAMETER, DIMENSION(-1:7) :: G_BC_TYPE=(/-99,0,-1,-2,-3,-4,-5,-6,-7/)
  !-99 : internal surface: macro contact surface
  !  0 : external: vacuum or albedo (isotropic reflexion)
  ! -1 : external: specular reflexion
  !      (approx. specular reflexion => G_BC_TYPE(0))
  ! -2 : external: translation
  ! -3 : external: rotation
  ! -4 : external: axial symmetry
  ! -5 : external: central symetry
  !       boundary condition length definition
  INTEGER, PARAMETER, DIMENSION(-1:5) :: G_BC_LEN=(/1,1,1,2,5,5,2/)
  ! internal : albedo
  ! vacuum surface : albedo
  ! specular reflexion : none
  ! translation :     tx ty (t=translation vector)
  ! rotation :        cx cy cos(theta) sin(theta) theta
  !                   (c= center,theta= axis angle)
  ! axial symmetry :  cx cy cos(theta) sin(theta) theta
  !                   (c= center,theta= axis angle)
  ! central symetry : cx cy (c= center)
  INTEGER, PARAMETER :: G_BC_MAX_LEN=6
  ! max bc data length
  INTEGER, PARAMETER :: NRPAR=6, NIPAR=3
  !       TYPGEO              = type of geometry:
  !                             0 = geometry with vacuum or isotropic reflexion
  !                             1 = geometry with symmetries of two axis of angle pi/n,n>0
  !                             2 = geometry with rotation of angle 2*pi/n,n>1
  !                             3 = 1/4 assembly with symmetries on all sides
  !                             5 = rectangular geometry with translation on all sides
  !                             6 = rectangular geometry with symmetry on all sides
  !                             7 = 1/8 assembly with symmetries on all sides
  !                             8 = equilateral triangle geometry with symmetries on all sides
  !                             9 = hexagonal geometry with translations on all sides
  !                            10 = equilateral triangle geometry with RA60 rotation and translation
  !                            11 = lozenge geometry with R120 rotation and translation
  !                            12 = rectangular S30 triangle geometry with symmetries on all sides
  !       NBFOLD              = n in angle definition of rotation or symmetry geometry
  !       NBMED               = nber of media in library file
  !       ALLSUR              = (have not been programmed!)
  INTEGER            :: TYPGEO,NBFOLD,ALLSUR
  !       NANIS               = 0: isotropic scattering in the laboratory system/ >1: anisotropic
  INTEGER            :: NANIS
  !       ISPEC               = 0: isotropic boundary condition/ 1: specular boundary condition
  INTEGER            :: ISPEC
  !       LGALLS              = TRUE, if 'allsur' equal to 1
  LOGICAL            :: LGALLS
  !       EPS                 = max. distance for two points to be considered as one point
  !                             in the topological geometric change
  REAL               :: EPS
  !       ANGGEO              = angle of rotation or symmetry geometry
  !       LENGTHX             = length in x direction of a rectangular geometry
  !       LENGTHY             = length in y direction of a rectangular geometry
  REAL(PDB)          :: ANGGEO,LENGTHX,LENGTHY
  !       INDEX               = when >0, impression of geometry data
  !       KNDEX               = when >0, impression of changes of geometry data in topological test
  !       PREC                = 0: geometry data is in format 4e20.12,
  !                             1: geometry data is in format 5e12.6
  INTEGER            :: INDEX,KNDEX,PREC
  !       NBMED               = nber of physical media
  INTEGER            :: NBMED
  !       LBCDIAG             = detection of diagonal symmetry in Cartesian geometry cases
  !       LGSPEC              = detection of a surfacic file limited to cyclic tracking
  !       LMERGM              = flag to perform a merge mix on nodes
  LOGICAL            :: LBCDIAG,LGSPEC,LMERGM
  TYPE T_SALBCDATA
     INTEGER                      :: SALTYPE
  !         read type of bc and nber of elements affected
  !         - TYPE of bc = 0 ~ 5
  !               0 (vacuum + albedo), 1 (specular reflexion),
  !               2 (translation),     3 (rotation),
  !               4 (axial symmetry),  5 (central symmetry)
  !         - NBER = number of elements affected
  !
     INTEGER                      :: NBER
     INTEGER,DIMENSION(:),POINTER :: ELEMNB
  !    BCDATA(1:2) : perimeter origin
  !    BCDATA(3)   : perimeter cos(angle)
  !    BCDATA(4)   : perimeter sin(angle)
  !    BCDATA(5)   : perimeter angle (radians)
  !    BCDATA(6)   : perimeter albedo
     REAL(PDB),DIMENSION(G_BC_MAX_LEN) :: BCDATA
  END TYPE T_SALBCDATA
  !
  !       geometry basic
  TYPE T_G_BASIC
     INTEGER, DIMENSION (N_IN) :: V_IN
     LOGICAL, DIMENSION (N_LG) :: V_LG
     REAL,    DIMENSION (N_RE) :: V_RE
     !
     INTEGER :: NB_ELEM ! NUMBER OF SURFACIC ELEMENTS
     INTEGER :: NB_MACRO ! NUMBER OF MACROS
     INTEGER :: NB_FLUX ! NUMBER OF FLUX VALUES
     !       definition of elements in a macro
     INTEGER, POINTER, DIMENSION(:,:) &
          :: IPAR         ! descriptors for elements
     ! dim: IPAR(NIPAR, NB_ELEM)
     REAL(PDB), POINTER, DIMENSION(:,:) &
          :: RPAR         ! DESCRIPTORS FOR ELEMENTS
     ! dim: RPAR(NRPAR, NB_ELEM)
     TYPE(T_SALBCDATA), POINTER, DIMENSION(:) &
          :: BCDATAREAD   ! BC DATA RECOVERED FROM SALGET
     REAL(PDB)       :: ALBEDO
     INTEGER         :: DEFAUL,NBBCDA
     ! dim: BCDATA(NBBCDA)
     INTEGER, POINTER, DIMENSION(:) &
          :: IBC2_ELEM    ! RELATIVE 2D BC NBER
     ! 0: is not a bc
     ! else: order number of bc
     ! dim: IBC2_ELEM(NB_ELEM)
     INTEGER, POINTER, DIMENSION(:) &
          :: ISURF2_ELEM  ! RELATIVE 2D SURFACE NBER
     ! 0: is not a surface
     ! else: order number of surface
     ! dim: ISURF2_ELEM(NB_ELEM)
     ! def: 
     !       definition of nodes in a macro
     INTEGER &
          :: NB_NODE      ! NUMBER OF NODES
     REAL(PDB), POINTER, DIMENSION(:) &
          :: VOL_NODE     ! VOLUME OF NODES, ALLOCATED ONLY IN 2D
     ! dim: VOL_NODE(NB_NODE)
     INTEGER, POINTER, DIMENSION(:) &
          :: PPERIM_NODE  ! POINTER TO LIST OF ELEMENTS OF PERIMETER
     ! dim: pperim_node(nb_node+1)
     INTEGER, POINTER, DIMENSION(:) &
          :: PERIM_NODE   ! ARRAY OF ELEMENTS IN PERIMETER OF NODES
     ! dim: perim_node(pperim_node(nb_node+1)-1)
     !       definition of boundary conditions of a macro
     INTEGER &
          :: NB_BC2       ! NUMBER OF 2D BOUNDARY CONDITIONS
     INTEGER, POINTER, DIMENSION(:) &
          :: TYPE_BC2     ! 2D BOUNDARY CONDITION TYPE
     ! dim: TYPE_BC2(NB_BC2)
     INTEGER, POINTER, DIMENSION(:) &
          :: IDATA_BC2    ! POSITION OF 2D BC DATA PER 2D BC
     ! dim: IDATA_BC2(NB_BC2)
     !       definition of surfaces of a macro
     INTEGER &
          :: NB_SURF2     ! NUMBER OF 2D SURFACES
     ! macro contacting bc (type -1)
     ! or external vacuum bc (type 0)
     INTEGER, POINTER, DIMENSION(:) &
          :: IBC2_SURF2   ! RELATIVE 2D BC ORDER NBER
     ! dim: IBC2_SURF2(NB_SURF2)
     INTEGER, POINTER, DIMENSION(:) &
          :: IELEM_SURF2  ! RELATIVE ELEM ORDER NBER
     ! dim: IELEM_SURF2(NB_SURF2)
     REAL(PDB), POINTER, DIMENSION(:) &
          :: SURF2        ! LOCAL 2D AREAS (POINTER TO SURF2_TAB)
     ! dim: SURF2(NB_SURF2)
     !       definition of perimeter of a macro
     INTEGER &
          :: NPERIM_MAC2  ! NUMBER OF ELEMENTS IN PERIMETER OF 2D MACRO
     INTEGER, POINTER, DIMENSION(:) &
          :: PERIM_MAC2   ! ELEMENTS IN PERIMETER OF 2D MACRO
     ! in case of geometry with rotation or symmetry,order of elems is:
     ! (elems on axis 1)+(elems on axis 2)+(other elems)
     ! dim: PERIM_MAC2(NPERIM_MAC2)
     INTEGER, POINTER, DIMENSION(:) &
          :: PPERIM_MAC2  ! POINTER TO TABLE 'PERIM_MAC2'&'DIST_AXIS' (ONLY FOR TYPGEO=1,2,5)
     ! for TYPGEO=1,2:
     ! (1): first elem on axis 1 (2): first elem on axis 2
     ! (3): first of other elems (4): nperim_mac2 + 1
     ! for TYPGEO=5:
     ! (1): first elem on axis 1 (2): first elem on axis 2
     ! (3): first elem on axis 3 (4): first elem on axis 4
     ! (5): NPERIM_MAC2 + 1
     ! dim: PPERIM_MAC2(5)
     INTEGER, POINTER, DIMENSION(:) &
          :: IDATA_AXIS   ! POSITION OF AXIAL BC DATA IN 'BCDATA' (ONLY FOR TYPGEO=1,2)
     ! (1): data position for axis 1 (2): data position for axis 2
     REAL(PDB), POINTER, DIMENSION(:) &
          :: DIST_AXIS    ! DISTANCE OF POINTS ON AXIS TO THE CENTER (0,0) (ONLY FOR TYPGEO=1,2)
     ! order of elems is:
     ! (distances of points on axis 1)+(distances of points on axis 2)
     ! dim: DIST_AXIS(1:PPERIM_MAC2(3)-1)
     INTEGER, POINTER, DIMENSION(:) &
          :: MED          ! (:): MEDIUM PER LOCAL NODE
     ! dim: MED(NBNODE)
     REAL(PDB), POINTER, DIMENSION(:,:) &
          :: BCDATA       ! TABLE OF BC DATA
     ! dim: BCDATA(G_BC_MAX_LEN,NT_BC)
     CHARACTER(LEN=12), POINTER, DIMENSION(:) :: NAME_MACRO
     INTEGER, POINTER, DIMENSION(:) :: NUM_MERGE
     ! NUM_MERGE : merge index per node
     INTEGER, POINTER, DIMENSION(:) :: NUM_MACRO
     ! NUM_MACRO : macro index per flux region
     INTEGER            :: NALBG
     ! NALBG : number of boundary condition types in BCDATA
  END TYPE T_G_BASIC
  !
END MODULE SAL_GEOMETRY_TYPES
