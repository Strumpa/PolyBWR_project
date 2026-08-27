!
!---------------------------------------------------------------------
!
!Purpose:
! To store common constants relarive to accuracy or to mathematics.
!
!Copyright:
! Copyright (C) 2001 Ecole Polytechnique de Montreal.
!
!Author(s):
! X. Warin
!
!---------------------------------------------------------------------
!
MODULE PRECISION_AND_KINDS
  INTEGER, PARAMETER :: PDB = SELECTED_REAL_KIND(8)
  REAL(PDB), PARAMETER :: PI=3.14159265358979_PDB, TWOPI=2._PDB*PI, &
       HALFPI=PI/2._PDB, SMALL=1.E-20_PDB,INFINITY=1.E20_PDB
  REAL(PDB)          :: HUGE_PDB,LOG10_HUGE_PDB,TINY_PDB,LOG10_TINY_PDB
END MODULE PRECISION_AND_KINDS
