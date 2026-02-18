*DECK SAPSPH
      SUBROUTINE SAPSPH(IPMICR,NG,NMIL,ILOC,NLOC,RVALOC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To recover a set of SPH equivalence factors from a microlib and store
* them as local variables.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPMICR  pointer to the microlib (L_LIBRARY signature).
* NG      number of condensed energy groups.
* NMIL    number of mixtures in the Saphyb.
* ILOC    position of local parameter in RVALOC.
* NLOC    first dimension of matrix RVALOC.
*
*Parameters: output
* RVALOC  local variable values in mixtures located in RVALOC(ILOC,:).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMICR
      INTEGER NG,NMIL,ILOC,NLOC
      REAL RVALOC(NLOC,NMIL)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      TYPE(C_PTR) JPEDIT,KPEDIT
      INTEGER ISTATE(NSTATE)
      REAL, ALLOCATABLE, DIMENSION(:) :: WORK
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(WORK(NMIL))
*
      CALL LCMSIX(IPMICR,'MACROLIB',1)
      CALL LCMGET(IPMICR,'STATE-VECTOR',ISTATE)
      IF(ISTATE(1).NE.NG) CALL XABORT('SAPSPH: BAD VALUE OF NG.')
      IF(ISTATE(2).NE.NMIL) CALL XABORT('SAPSPH: BAD VALUE OF NMIL.')
*----
*  RECOVER SPH EQUIVALENCE FACTORS.
*----
      JPEDIT=LCMGID(IPMICR,'GROUP')
      DO 30 IGR=1,NG
      KPEDIT=LCMGIL(JPEDIT,IGR)
      CALL LCMLEN(KPEDIT,'NSPH',ILONG,ITYLCM)
      IF(ILONG.GT.0) THEN
         CALL LCMGET(KPEDIT,'NSPH',WORK)
         DO 10 IMIL=1,NMIL
         RVALOC(ILOC+IGR-1,IMIL)=WORK(IMIL)
   10    CONTINUE
      ELSE
         DO 20 IMIL=1,NMIL
         RVALOC(ILOC+IGR-1,IMIL)=1.0
   20    CONTINUE
      ENDIF
   30 CONTINUE
      CALL LCMSIX(IPMICR,' ',2)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(WORK)
      RETURN
      END
