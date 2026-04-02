MODULE TONDST_CACHE_MOD
!
!-----------------------------------------------------------------------
!
!Purpose:
! Cache module for TONDST to avoid redundant DOORPV/AV calls when
! macroscopic cross sections have not changed significantly between
! calls. In TONE: self-shielding, DOORPV/AV is called once per INRS
! per iteration. Between INRS calls, only the previous INRS's
! self-shielding correction changes the total XS -- a small change
! that decreases with iteration convergence. By caching the TXSC
! values used in the last DOORPV/AV call for each group, we can skip
! recomputing Pij/ARM info for groups where the change is below a
! tolerance.
!
!Author(s): B. Cui (2026)
!
!-----------------------------------------------------------------------
!
  IMPLICIT NONE
  SAVE
  INTEGER :: TC_NBM = -1
  INTEGER :: TC_NGRO = -1
  INTEGER :: TC_HITS = 0
  INTEGER :: TC_TOTAL = 0
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: TC_VALID  ! (NGRO)
CONTAINS
  SUBROUTINE TC_INIT(NBM, NGRO)
    INTEGER, INTENT(IN) :: NBM, NGRO
    IF(TC_NBM .NE. NBM .OR. TC_NGRO .NE. NGRO) THEN
      IF(ALLOCATED(TC_VALID)) DEALLOCATE(TC_VALID)
      ALLOCATE(TC_VALID(NGRO))
      TC_NBM = NBM
      TC_NGRO = NGRO
      TC_VALID(:) = .FALSE.
      TC_HITS = 0
      TC_TOTAL = 0
    ENDIF
  END SUBROUTINE

  SUBROUTINE TC_RESET()
    TC_NBM = -1
    TC_NGRO = -1
    TC_HITS = 0
    TC_TOTAL = 0
    IF(ALLOCATED(TC_VALID)) DEALLOCATE(TC_VALID)
  END SUBROUTINE
END MODULE
