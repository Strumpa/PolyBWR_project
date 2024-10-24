*DECK DRVREC
      SUBROUTINE DRVREC(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover one or many LCM objects.
*
*Copyright:
* Copyright (C) 1994 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file:
*         HENTRY(1): read-only or modification type(VECTOR).
* IENTRY  type of each LCM object or file:
*         =1 LCM memory object; =2 XSM file; =3 sequential binary file;
*         =4 sequential ascii file.
* JENTRY  access of each LCM object or file:
*         =0 the LCM object or file is created;
*         =1 the LCM object or file is open for modifications;
*         =2 the LCM object or file is open in read-only mode.
* KENTRY  LCM object address or file unit number.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
      TYPE(C_PTR) KENTRY(NENTRY)
      CHARACTER HENTRY(NENTRY)*12
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) IPLIST,JPLIST,KPLIST
      CHARACTER HMEDIA*12,TEXT12*12,TEXT4*4,NAMT*12
      DOUBLE PRECISION DFLOTT
*
      IF(NENTRY.LE.1) CALL XABORT('DRVREC: TWO PARAMETERS EXPECTED.')
      ITYPE=0
      JPLIST=C_NULL_PTR
      DO 10 I=1,NENTRY
      IF(JENTRY(I).EQ.2) THEN
         ITYPE=IENTRY(I)
         IPLIST=KENTRY(I)
         HMEDIA=HENTRY(I)
         IF((IENTRY(I).NE.1).AND.(IENTRY(I).NE.2)) CALL XABORT('DRVREC:'
     1   //' RHS LINKED LIST OR XSM FILE EXPECTED.')
         GO TO 20
      ENDIF
   10 CONTINUE
      CALL XABORT('DRVREC: UNABLE TO FIND A BACKUP MEDIA OPEN IN READ-O'
     1 //'NLY MODE.')
*
   20 IMPX=1
      IPOS=0
   30 CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
      IF(INDIC.EQ.10) GO TO 40
      IF(INDIC.NE.3) CALL XABORT('DRVREC: CHARACTER DATA EXPECTED.')
      IF(TEXT4.EQ.'EDIT') THEN
         CALL REDGET(INDIC,IMPX,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('DRVREC: INTEGER DATA EXPECTED.')
      ELSE IF(TEXT4.EQ.'ITEM') THEN
         CALL REDGET(INDIC,IPOS,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('DRVREC: INTEGER DATA EXPECTED.')
      ELSE IF(TEXT4.EQ.'STEP') THEN
*        CHANGE THE HIERARCHICAL LEVEL ON THE LCM OBJECT.
         IF(ITYPE.GT.2) CALL XABORT('DRVREC: UNABLE TO STEP INTO A SE'
     1   //'QUENTIAL FILE.')
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.3) CALL XABORT('DRVREC: CHARACTER DATA EXPECTED.')
         IF(TEXT4.EQ.'UP') THEN
            CALL REDGET(INDIC,NITMA,FLOTT,NAMT,DFLOTT)
            IF(INDIC.NE.3) CALL XABORT('DRVREC: CHARACTER DATA EXPECT'
     1      //'ED.')
            IF(IMPX.GT.0) WRITE (6,100) NAMT
            JPLIST=LCMGID(IPLIST,NAMT)
         ELSE IF(TEXT4.EQ.'AT') THEN
            CALL REDGET(INDIC,NITMA,FLOTT,NAMT,DFLOTT)
            IF(INDIC.NE.1) CALL XABORT('DRVREC: INTEGER EXPECTED.')
            IF(IMPX.GT.0) WRITE (6,110) NITMA
            JPLIST=LCMGIL(IPLIST,NITMA)
         ELSE
            CALL XABORT('DRVREC: UP OR AT EXPECTED.')
         ENDIF
         IPLIST=JPLIST
      ELSE IF(TEXT4.EQ.';') THEN
         GO TO 40
      ELSE
         CALL XABORT('DRVREC: '//TEXT4//' IS AN INVALID KEY WORD.')
      ENDIF
      GO TO 30
*
   40 CALL LCMGTC(IPLIST,'SIGNATURE',12,1,TEXT12)
      IF(TEXT12.NE.'L_ARCHIVE') THEN
         CALL XABORT('DRVREC: SIGNATURE OF '//HMEDIA//' IS '//TEXT12//
     1   '. L_ARCHIVE EXPECTED.')
      ENDIF
      DO 50 I=1,NENTRY-1
      IF((JENTRY(I).EQ.0).OR.(JENTRY(I).EQ.1)) THEN
         IF(IENTRY(I).GT.2) CALL XABORT('DRVREC: LHS LINKED LIST OR XSM'
     1   //' FILE EXPECTED.')
         IF(IMPX.GT.0) THEN
           IF(IPOS.EQ.0) THEN
             WRITE (6,'(/18H DRVREC: RECOVER '',A,8H'' FROM '',A,
     1       2H''.)') TRIM(HENTRY(I)),TRIM(HMEDIA)
           ELSE
             WRITE (6,'(/22H DRVREC: RECOVER ITEM=,I5,5H OF '',A,
     1       8H'' FROM '',A,2H''.)') IPOS,TRIM(HENTRY(I)),TRIM(HMEDIA)
           ENDIF
         ENDIF
         TEXT12=HENTRY(I)
         CALL LCMLEN(IPLIST,TEXT12,ILEN,ITYLCM)
         IF(ILEN.EQ.0) THEN
            CALL LCMLIB(IPLIST)
            CALL XABORT('DRVREC: UNABLE TO FIND '//TEXT12//' ON THE BA'
     1      //'CKUP MEDIA NAMED '//HMEDIA//'.')
         ELSE IF(ITYLCM.EQ.0) THEN
            IF(IPOS.NE.0) CALL XABORT('DRVREC: RECORD '//TEXT12//' ON '
     1      //'THE BACKUP MEDIA NAMED '//HMEDIA//' IS NOT A DIRECTORY.')
            CALL LCMSIX(IPLIST,HENTRY(I),1)
            CALL LCMEQU(IPLIST,KENTRY(I))
            CALL LCMSIX(IPLIST,' ',2)
         ELSE IF(ITYLCM.EQ.10) THEN
            IF(IPOS.EQ.0) CALL XABORT('DRVREC: RECORD '//TEXT12//' ON '
     1      //'THE BACKUP MEDIA NAMED '//HMEDIA//' IS NOT A LIST.')
            JPLIST=LCMGID(IPLIST,HENTRY(I))
            KPLIST=LCMGIL(JPLIST,IPOS)
            CALL LCMEQU(KPLIST,KENTRY(I))
         ELSE
            CALL LCMLIB(IPLIST)
            CALL XABORT('DRVREC: RECORD '//TEXT12//' ON THE BACKUP MED'
     1      //'IA NAMED '//HMEDIA//' CANNOT BE COPIED.')
         ENDIF
      ENDIF
   50 CONTINUE
      RETURN
*
  100 FORMAT (/27H DRVREC: STEP UP TO LEVEL ',A12,2H'.)
  110 FORMAT (/26H DRVREC: STEP AT COMPONENT,I6,1H.)
      END
