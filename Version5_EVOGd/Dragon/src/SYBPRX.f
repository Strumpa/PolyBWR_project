*DECK SYBPRX
      SUBROUTINE SYBPRX (IND,NCOUR,IPAS,IKG,SIGT,SIGW,P,PIS,PSJ,PSS)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Print the cell-wise collision probabilities in SYBRX- modules.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IND     reduction flag (1/2: matrix reduction off/on).
* NCOUR   total number of surfaces.
* IPAS    total number of volumes.
* IKG     generating cell indices.
* SIGT    total macroscopic cross sections.
* SIGW    scattering macroscopic cross sections.
* P       reduced collision probabilities.
* PIS     volume to surface probabilities.
* PSJ     surface to volume probabilities.
* PSS     surface to surface probabilities.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IND,NCOUR,IPAS,IKG
      REAL SIGT(IPAS),SIGW(IPAS),P(IPAS,IPAS),PIS(IPAS,NCOUR),
     1 PSJ(NCOUR,IPAS),PSS(NCOUR,NCOUR)
*
      IF(IND.EQ.1) THEN
        WRITE (6,100) IKG
        WRITE (6,120) (SIGT(I),I=1,IPAS)
      ELSE
        WRITE (6,110) IKG
        WRITE (6,120) (SIGT(I),I=1,IPAS)
        WRITE (6,130)
        WRITE (6,120) (SIGW(I),I=1,IPAS)
      ENDIF
      WRITE (6,'(/16H P(I,J) MATRIX :/)')
      DO 10 I=1,IPAS
      WRITE (6,120) (P(I,J),J=1,IPAS)
10    CONTINUE
      WRITE (6,'(/16H PIS(I) MATRIX :/)')
      DO 20 I=1,IPAS
      WRITE (6,120) (PIS(I,J),J=1,NCOUR)
20    CONTINUE
      WRITE (6,'(/16H PSJ(I) MATRIX :/)')
      DO 30 I=1,IPAS
      WRITE (6,120) (PSJ(J,I),J=1,NCOUR)
30    CONTINUE
      WRITE (6,'(/13H PSS MATRIX :/)')
      DO 40 I=1,NCOUR
      WRITE (6,120) (PSS(I,J),J=1,NCOUR)
40    CONTINUE
      WRITE (6,'(//)')
      RETURN
100   FORMAT (/32H SYBPRX: NO SCATTERING REDUCTION/16H GENERATING CELL,
     1 3H NB,I4//35H TOTAL MACROSCOPIC CROSS SECTIONS :)
110   FORMAT (/29H SYBPRX: SCATTERING REDUCTION/19H GENERATING CELL NB,
     1 I4//35H TOTAL MACROSCOPIC CROSS SECTIONS :)
120   FORMAT (1X,1P,10E13.5)
130   FORMAT(/40H SCATTERING MACROSCOPIC CROSS SECTIONS :)
      END
