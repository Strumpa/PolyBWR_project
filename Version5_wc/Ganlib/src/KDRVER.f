*DECK KDRVER
      SUBROUTINE KDRVER(REV,DATE)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To extract CVS or SVN version and production date.
*
*Copyright:
* Copyright (C) 2006 Ecole Polytechnique de Montreal
*
*Author(s): A. Hebert
*
*Parameters: output
* REV     revision character identification
* DATE    revision character date 
*
*-----------------------------------------------------------------------
*
      CHARACTER REV*48,DATE*64
*
      REV='Version 5.0.12 ($Revision: 3956 $)'
      DATE='$Date: 2025-09-05 09:32:25 -0400 (Fri, 05 Sep 2025) $'
      IF(REV(22:).EQ.'ion$)') THEN
*        CVS or SVN keyword expansion not performed
         REV='Version 5.0.12'
         DATE='September 5, 2025'
      ENDIF
      RETURN
      END
