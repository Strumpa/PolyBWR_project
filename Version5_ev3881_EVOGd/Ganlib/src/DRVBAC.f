*DECK DRVBAC
      SUBROUTINE DRVBAC(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Backup one or many LCM objects.
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
      CHARACTER HENTRY(NENTRY)*12
      TYPE(C_PTR) KENTRY(NENTRY)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) IPLIST,JPLIST,KPLIST
      CHARACTER TEXT12*12,TEXT4*4,HMEDIA*12,NAMT*12
      DOUBLE PRECISION DFLOTT
*
      IF(NENTRY.LE.1) THEN
         CALL XABORT('DRVBAC: TWO PARAMETERS EXPECTED.')
      ELSE IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2)) THEN
         CALL XABORT('DRVBAC: LHS LINKED LIST OR XSM FILE EXPECTED.')
      ELSE IF(JENTRY(1).EQ.2) THEN
         CALL XABORT('DRVBAC: LHS PARAMETER IN CREATE OR MODIFICATION '
     1   //'MODE EXPECTED.')
      ENDIF
      ITYPE=IENTRY(1)
      IPLIST=KENTRY(1)
*
      IMPX=1
      IDIM=0
      IPOS=0
      JPLIST=C_NULL_PTR
   10 CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
      IF(INDIC.EQ.10) GO TO 30
      IF(INDIC.NE.3) CALL XABORT('DRVBAC: CHARACTER DATA EXPECTED.')
      IF(TEXT4.EQ.'EDIT') THEN
         CALL REDGET(INDIC,IMPX,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('DRVBAC: INTEGER DATA EXPECTED.')
      ELSE IF(TEXT4.EQ.'STEP') THEN
*        CHANGE THE HIERARCHICAL LEVEL ON THE LCM OBJECT.
         IF(ITYPE.GT.2) CALL XABORT('DRVBAC: UNABLE TO STEP INTO A SE'
     1   //'QUENTIAL FILE.')
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.3) CALL XABORT('DRVBAC: CHARACTER DATA EXPECTED.')
         IF(TEXT4.EQ.'UP') THEN
            CALL REDGET(INDIC,NITMA,FLOTT,NAMT,DFLOTT)
            IF(INDIC.NE.3) CALL XABORT('DRVBAC: CHARACTER DATA EXPECT'
     1      //'ED.')
            IF(IMPX.GT.0) WRITE (6,100) NAMT
            CALL LCMLEN(IPLIST,NAMT,ILONG,ITYLCM)
            IF(ILONG.GT.0) THEN
               JPLIST=LCMGID(IPLIST,NAMT)
            ELSE
               JPLIST=LCMDID(IPLIST,NAMT)
            ENDIF
         ELSE IF(TEXT4.EQ.'AT') THEN
            CALL REDGET(INDIC,NITMA,FLOTT,NAMT,DFLOTT)
            IF(INDIC.NE.1) CALL XABORT('DRVBAC: INTEGER EXPECTED.')
            IF(IMPX.GT.0) WRITE (6,110) NITMA
            CALL LCMLEL(IPLIST,NITMA,ILONG,ITYLCM)
            IF(ILONG.GT.0) THEN
               JPLIST=LCMGIL(IPLIST,NITMA)
            ELSE
               JPLIST=LCMDIL(IPLIST,NITMA)
            ENDIF
         ELSE
            CALL XABORT('DRVBAC: UP OR AT EXPECTED.')
         ENDIF
         IPLIST=JPLIST
      ELSE IF(TEXT4.EQ.'LIST') THEN
         CALL REDGET(INDIC,IDIM,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('DRVBAC: INTEGER DATA EXPECTED.')
         CALL LCMPUT(IPLIST,'LISTDIM',1,1,IDIM)
      ELSE IF(TEXT4.EQ.'ITEM') THEN
         CALL REDGET(INDIC,IPOS,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('DRVBAC: INTEGER DATA EXPECTED.')
      ELSE IF(TEXT4.EQ.';') THEN
         GO TO 30
      ELSE
         CALL XABORT('DRVBAC: '//TEXT4//' IS AN INVALID KEY WORD.')
      ENDIF
      GO TO 10
*
   30 CALL LCMLEN(IPLIST,'SIGNATURE',ILONG,ITYLCM)
      IF(ILONG.NE.0) THEN
         CALL LCMGTC(IPLIST,'SIGNATURE',12,TEXT12)
         IF(TEXT12.NE.'L_ARCHIVE') THEN
            HMEDIA=HENTRY(1)
            CALL XABORT('DRVBAC: SIGNATURE OF '//HMEDIA//' IS '//TEXT12
     1      //'. L_ARCHIVE EXPECTED.')
         ENDIF
      ELSE
         TEXT12='L_ARCHIVE'
         CALL LCMPTC(IPLIST,'SIGNATURE',12,TEXT12)
      ENDIF
      ISET=0
      DO 40 I=2,NENTRY
      IF((JENTRY(I).EQ.0).OR.(JENTRY(I).EQ.1)) THEN
         TEXT12=HENTRY(I)
         CALL XABORT('DRVBAC: ENTRY '//TEXT12//' IS NOT EXPECTED.')
      ELSE IF(IENTRY(I).GT.2) THEN
         CALL XABORT('DRVBAC: RHS LINKED LIST OR XSM FILE EXPECTED.')
      ENDIF
      IF(IDIM.EQ.0) THEN
        CALL LCMLEN(IPLIST,'LISTDIM',ILONG,ITYLCM)
        IF(ILONG.EQ.1) CALL LCMGET(IPLIST,'LISTDIM',IDIM)
      ENDIF
      IF(IDIM.EQ.0) THEN
        ! HENTRY(I) is stored as a directory
        IF(IMPX.GT.0) WRITE (6,'(/17H DRVBAC: BACKUP '',A12,7H'' INTO ,
     1  1H'',A,2H''.)') TRIM(HENTRY(I)),TRIM(HENTRY(1))
        CALL LCMSIX(IPLIST,HENTRY(I),1)
        CALL LCMEQU(KENTRY(I),IPLIST)
        CALL LCMSIX(IPLIST,' ',2)
      ELSE
        ! HENTRY(I) is stored as a list of directories
        IF(IPOS.EQ.0) CALL XABORT('DRVBAC: IPOS IS NOT DEFINED.')
        IF(IPOS.GT.IDIM) CALL XABORT('DRVBAC: LIST OVERFLOW FOR OBJECT'
     1  //' '//TRIM(HENTRY(I))//'.')
        JPLIST=LCMLID(IPLIST,HENTRY(I),IPOS)
        IF(IMPX.GT.0) WRITE (6,120) TRIM(HENTRY(I)),IPOS,TRIM(HENTRY(1))
        KPLIST=LCMDIL(JPLIST,IPOS)
        CALL LCMEQU(KENTRY(I),KPLIST)
      ENDIF
   40 CONTINUE
      RETURN
  100 FORMAT (/27H DRVBAC: STEP UP TO LEVEL ',A12,2H'.)
  110 FORMAT (/26H DRVBAC: STEP AT COMPONENT,I6,1H.)
  120 FORMAT (/16H DRVBAC: BACKUP ,A,13H INTO ELEMENT,I5,9H OF LIST ,A,
     1 1H.)
      END
