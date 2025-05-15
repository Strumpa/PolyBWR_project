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
      REV='Version 5.0.11 ($Revision: 3787 $)'
      DATE='$Date: 2025-05-10 13:12:46 -0400 (Sat, 10 May 2025) $'
      IF(REV(22:).EQ.'ion$)') THEN
*        CVS or SVN keyword expansion not performed
         REV='Version 5.0.11'
         DATE='April 13, 2025'
      ENDIF
      RETURN
      END
