*DECK EXCELT
      SUBROUTINE EXCELT(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* EXCELL tracking operator.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version.
*
*Author(s): A. Hebert and R. Roy
*
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file:
*         HENTRY(1): creation or modification type(L_TRACK);
*         HENTRY(2): sequential binary tracking file;
*         HENTRY(3): read-only type(L_GEOM).
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
      USE        GANLIB
      IMPLICIT   NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER      NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
      TYPE(C_PTR)  KENTRY(NENTRY)
      CHARACTER    HENTRY(NENTRY)*12
*----
*  LOCAL VARIABLES
*----
      INTEGER    IOUT,NSTATE
      PARAMETER (IOUT=6,NSTATE=40)
      CHARACTER  TEXT4*4,TEXT12*12,TITLE*72,HSIGN*12,CFTRAK*12
      LOGICAL    LASS,LDRASS,LPRISM,LBIHET
      DOUBLE PRECISION DFLOTT
      INTEGER    ITITL(18),ISTATE(NSTATE),IZ
      REAL       EXTKOP(NSTATE), FLOTT, CUTOFX, DELU, FRTM
*
      TYPE(C_PTR) IPTRK, IPGEOM
      INTEGER    IFTRAK, IDISP, IMPX, MAXPTS, NANIS, NORE, LMERG, I,
     >           ISYMM, KSPEC, KTOPT, KMODL, INDIC, NITMA, LCACT,
     >           NMU, INSB , IQUA10, NBATCH, IBIHET, ILONG, ITYLCM
*----
*  PARAMETER VALIDATION
*----
      IF(NENTRY.LT.2) CALL XABORT
     >  ('EXCELT: AT LEAST TWO DATA STRUCTURES REQUIRED')
      IF(IENTRY(1).GT.2) CALL XABORT
     >  ('EXCELT: FIRST DATA STRUCTURE NOT A LCM OBJECT')
      IF(JENTRY(1).NE.0.AND.JENTRY(1).NE.1) CALL XABORT
     >  ('EXCELT: FIRST DATA STRUCTURE NOT IN CREATE OR MODIFY MODE')
      IPTRK=KENTRY(1)
*----
*  RECOVER GEOMETRY
*----
      IPGEOM=C_NULL_PTR
      DO 10 I=2,NENTRY
      IF((IENTRY(I).LE.2).AND.(JENTRY(I).EQ.2)) THEN
         CALL LCMGTC(KENTRY(I),'SIGNATURE',12,HSIGN)
         IF(HSIGN.NE.'L_GEOM') THEN
            TEXT12=HENTRY(I)
            CALL XABORT('EXCELT: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1      '. L_GEOM EXPECTED.')
         ENDIF
         IPGEOM=KENTRY(I)
         ISTATE(:NSTATE)=0
         CALL LCMGET(IPGEOM,'STATE-VECTOR',ISTATE)
         TEXT12=HENTRY(I)
         CALL LCMPTC(IPTRK,'LINK.GEOM',12,TEXT12)
         GO TO 20
      ENDIF
   10 CONTINUE
*----
*  RECOVER SEQUENTIAL BINARY TRACKING FILE CHARACTERISTICS
*----
   20 CFTRAK=' '
      IFTRAK=0
      IDISP=99
      DO 30 I=2,NENTRY
      IF(IENTRY(I).EQ.3) THEN
         CFTRAK=HENTRY(I)
         IF(JENTRY(I).EQ.0) IDISP=1
         IF(JENTRY(I).EQ.1) IDISP=-1
         IF(JENTRY(I).EQ.2) IDISP=0
         IFTRAK=FILUNIT(KENTRY(I))
         GO TO 35
      ENDIF
   30 CONTINUE
*
   35 IMPX=1
      LMERG=1
      TITLE=' '
      IF(IDISP.NE.0) LMERG=0
      IF(JENTRY(1).EQ.0) THEN
        HSIGN='L_TRACK'
        CALL LCMPTC(IPTRK,'SIGNATURE',12,HSIGN)
        HSIGN='EXCELL'
        CALL LCMPTC(IPTRK,'TRACK-TYPE',12,HSIGN)
        IF(C_ASSOCIATED(IPGEOM)) THEN
          MAXPTS=ISTATE(6)
        ELSE
          MAXPTS=0
        ENDIF
        LPRISM=.FALSE.
        DELU=1.0
        NANIS=1
        NORE=0
        IF(IDISP.NE.0) NORE=-1
        KSPEC=-1
        KTOPT=-1
        CUTOFX=0.0
        ISYMM=1
        INSB=0
        IF(IFTRAK.EQ.0) INSB=2
        LCACT=-1
        NMU=0
        IQUA10=5
        NBATCH=1
        IBIHET=2
        CALL LCMLEN(IPGEOM,'BIHET',ILONG,ITYLCM)
        LBIHET=(ILONG.NE.0)
        IF(LBIHET) IQUA10=5
      ELSE IF(JENTRY(1).EQ.1) THEN
        CALL LCMGTC(IPTRK,'SIGNATURE',12,HSIGN)
        IF(HSIGN.NE.'L_TRACK') THEN
          TEXT12=HENTRY(1)
          CALL XABORT('EXCELT: SIGNATURE OF '//TEXT12//' IS '
     >    //HSIGN//' L_TRACK EXPECTED.')
        ENDIF
        CALL LCMGTC(KENTRY(1),'TRACK-TYPE',12,HSIGN)
        IF(HSIGN.NE.'EXCELL') THEN
          TEXT12=HENTRY(1)
          CALL XABORT('EXCELT: TRACK-TYPE OF '//TEXT12//' IS '
     >    //HSIGN//'. EXCELL EXPECTED.')
        ENDIF
        ISTATE(:NSTATE)=0
        CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
        MAXPTS=ISTATE(1)
        NANIS=ISTATE(6)
        KMODL=ISTATE(7)
        NORE=ISTATE(8)
        KTOPT=ISTATE(9)
        KSPEC=ISTATE(10)
        ISYMM=ISTATE(12)
        LCACT=ISTATE(13)
        NMU=ISTATE(14)
        INSB=ISTATE(22)
        NBATCH=ISTATE(27)
        IZ=ISTATE(39)
        LPRISM=(IZ.NE.0)
        CALL LCMGET(IPTRK,'EXCELTRACKOP',EXTKOP)
        CUTOFX=EXTKOP(1)
        DELU=EXTKOP(40)
        LBIHET=(ISTATE(40).GT.0)
        IF(LBIHET) THEN
           CALL LCMSIX(IPTRK,'BIHET',1)
           CALL LCMGET(IPTRK,'PARAM',ISTATE)
           CALL LCMSIX(IPTRK,'BIHET',2)
           IBIHET=ISTATE(6)
           IQUA10=ISTATE(8)
         ELSE
           IBIHET=0
           IQUA10=0
        ENDIF
        CALL LCMLEN(IPTRK,'TITLE',ILONG,ITYLCM)
        IF(ILONG.GT.0) CALL LCMGTC(IPTRK,'TITLE',72,TITLE)
      ENDIF
      FRTM=0.05
   40 CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
   41 CONTINUE
      IF(INDIC.EQ.10) GO TO 50
      IF(INDIC.NE.3) CALL XABORT('EXCELT: CHARACTER DATA EXPECTED.')
      IF(TEXT4.EQ.'EDIT') THEN
         CALL REDGET(INDIC,IMPX,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('EXCELT: INTEGER DATA EXPECTED(1).')
      ELSE IF(TEXT4.EQ.'MAXR') THEN
         CALL REDGET(INDIC,MAXPTS,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('EXCELT: INTEGER DATA EXPECTED(2).')
      ELSE IF(TEXT4.EQ.'TITL') THEN
         CALL REDGET(INDIC,NITMA,FLOTT,TITLE,DFLOTT)
         IF(INDIC.NE.3) CALL XABORT('EXCELT: TITLE EXPECTED.')
      ELSE IF(TEXT4(1:3).EQ.'PRI') THEN
         IF(.NOT.C_ASSOCIATED(IPGEOM)) THEN
            CALL XABORT('EXCELT: NO GEOMETRY TO PROJECT.')
         ENDIF
         LPRISM=.TRUE.
         IF (TEXT4(4:4).EQ.'Z') THEN
            IZ=3
         ELSEIF (TEXT4(4:4).EQ.'Y') THEN
            IZ=2
         ELSEIF (TEXT4(4:4).EQ.'X') THEN
            IZ=1
         ELSE
            CALL XABORT('EXCELT: INVALID PROJECTION AXIS FOR 3D PRISM.')
         ENDIF
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.2) THEN
            CALL XABORT('EXCELT: REAL DATA EXPECTED')
         ELSE
            DELU=1.0/FLOTT
            IF (DELU.LT.0.0)
     >        CALL XABORT('EXCELT: DELU > 0.0 EXPECTED')
         ENDIF
      ELSE IF(TEXT4.EQ.'ANIS') THEN
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) THEN
           CALL XABORT('EXCELT: INTEGER MUST FOLLOW ANIS KEYWORD')
         ELSE
           NANIS=NITMA
           IF(NANIS.LT.1)
     >       CALL XABORT('EXCELT: NANIS GREATER THAN 1 PERMITTED ONLY')
         ENDIF
      ELSE IF(TEXT4.EQ.'RENO') THEN
         IF(IDISP.EQ.0) CALL XABORT('EXCELT: CANNOT NORMALIZE READ-ONL'
     >   //'Y BINARY TRACKING FILE')
         NORE=0
      ELSE IF(TEXT4.EQ.'REND') THEN
         IF(IDISP.EQ.0) CALL XABORT('EXCELT: CANNOT NORMALIZE READ-ONL'
     >   //'Y BINARY TRACKING FILE')
         NORE=-1
      ELSE IF(TEXT4.EQ.'RENM') THEN
         IF(IDISP.EQ.0) CALL XABORT('EXCELT: CANNOT NORMALIZE READ-ONL'
     >   //'Y BINARY TRACKING FILE')
         NORE=-2
      ELSE IF(TEXT4.EQ.'NORE') THEN
         IF(IDISP.EQ.0) CALL XABORT('EXCELT: CANNOT NORMALIZE READ-ONL'
     >   //'Y BINARY TRACKING FILE')
         NORE=1
      ELSE IF(TEXT4.EQ.'TREG') THEN
         IF(IDISP.EQ.0) CALL XABORT('EXCELT: CANNOT NORMALIZE READ-ONL'
     >   //'Y BINARY TRACKING FILE')
         LMERG=0
      ELSE IF(TEXT4.EQ.'TMER') THEN
         IF(IDISP.EQ.0) CALL XABORT('EXCELT: CANNOT NORMALIZE READ-ONL'
     >   //'Y BINARY TRACKING FILE')
         LMERG=1
      ELSE IF(TEXT4.EQ.'PISO') THEN
         KSPEC=0
      ELSE IF(TEXT4.EQ.'PSPC') THEN
         KSPEC=1
      ELSE IF(TEXT4.EQ.'QUAB') THEN
         CALL REDGET(INDIC,IQUA10,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('EXCELT: INTEGER DATA EXPECTED(3).')
      ELSE IF(TEXT4.EQ.'BATC') THEN
         CALL REDGET(INDIC,NBATCH,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) CALL XABORT('EXCELT: INTEGER DATA EXPECTED(4).')
      ELSE IF(TEXT4.EQ.'SAPO') THEN
         IBIHET=1
      ELSE IF(TEXT4.EQ.'HEBE') THEN
         IBIHET=2
      ELSE IF(TEXT4.EQ.'SLSI') THEN
         IBIHET=3
         FRTM=0.05
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
         IF (INDIC.NE.2) GOTO 41
         FRTM=FLOTT
      ELSE IF(TEXT4.EQ.'SLSS') THEN
         IBIHET=4
         FRTM=0.05
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
         IF (INDIC.NE.2) GOTO 41
         FRTM=FLOTT
      ELSE IF(TEXT4.EQ.'CUT') THEN
         CALL REDGET(INDIC,NITMA, FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.2.OR.FLOTT.LT.0.0) THEN
            CALL XABORT('EXCELT: CUTOFF MUST BE A POSITIVE REAL')
         ENDIF
         CUTOFX=FLOTT
      ELSE IF(TEXT4.EQ.'ONEG') THEN
         INSB=0
      ELSE IF(TEXT4.EQ.'ALLG') THEN
         INSB=1
      ELSE IF(TEXT4.EQ.'XCLL') THEN
         INSB=2
      ELSE IF(TEXT4.EQ.'GAUS') THEN
         LCACT=0
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) GO TO 41
         NMU=NITMA
      ELSE IF(TEXT4.EQ.'CACA') THEN
         LCACT=1
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) GO TO 41
         NMU=NITMA
      ELSE IF(TEXT4.EQ.'CACB') THEN
         LCACT=2
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) GO TO 41
         NMU=NITMA
      ELSE IF(TEXT4.EQ.'LCMD') THEN
         LCACT=3
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) GO TO 41
         NMU=NITMA
      ELSE IF(TEXT4.EQ.'OPP1') THEN
         LCACT=4
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) GO TO 41
         NMU=NITMA
      ELSE IF(TEXT4.EQ.'OGAU') THEN
         LCACT=5
         CALL REDGET(INDIC,NITMA,FLOTT,TEXT4,DFLOTT)
         IF(INDIC.NE.1) GO TO 41
         NMU=NITMA
      ELSE IF(TEXT4.EQ.'TRAK') THEN
         IF(IDISP.LE.0) CALL XABORT('EXCELT: TRAK KEYWORD NOT REQUIRED')
         GO TO 50
      ELSE IF(TEXT4.EQ.';') THEN
         IF(IDISP.GT.0) CALL XABORT('EXCELT: TRAK KEYWORD EXPECTED')
         GO TO 50
      ELSE
         CALL XABORT('EXCELT: '//TEXT4//' IS AN INVALID KEY WORD.')
      ENDIF
      GO TO 40
*----
*  CALL XELDRV TO PERFORM THE TRACKING
*----
   50 IF(C_ASSOCIATED(IPGEOM)) LASS=LDRASS(IPGEOM,IMPX)
*
      READ(TITLE,'(18A4)') (ITITL(I),I=1,18)
      CALL LCMPUT(IPTRK,'TITLE',18,3,ITITL)
      IF(IMPX.GT.1) WRITE(IOUT,'(1X,A72//)') TITLE
*
      IF(MAXPTS.EQ.0) CALL XABORT('EXCELT: MAXPTS NOT DEFINED.')
      CALL XELDRV(IPTRK ,IPGEOM,IMPX  ,MAXPTS,NANIS ,NORE  ,
     >            LMERG, KSPEC , KTOPT,TITLE ,CUTOFX,CFTRAK,
     >            IFTRAK,IDISP ,ISYMM ,LCACT ,NMU   ,INSB  ,
     >            NBATCH,LBIHET,LPRISM,IZ,DELU,FRTM )
*----
*  PROCESS DOUBLE HETEROGENEITY (BIHET) DATA (IF AVAILABLE)
*----
      IF(LBIHET) CALL XDRTBH(IPGEOM,IPTRK,IQUA10,IBIHET,IMPX,FRTM)
*
      RETURN
      END
