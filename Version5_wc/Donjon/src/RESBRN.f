*DECK RESBRN
      SUBROUTINE RESBRN(IPMAP,NCH,NB,NCOMB,NX,NY,NZ,LRSCH,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Initialize the axial shape and compute the first burnup limits per
* bundle for every channel (used for the time-average model).
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
*
*Author(s): 
*  D. Sekki, I. Trancart
*
*Parameters: input
* IPMAP   pointer to fuel-map information.
* NCH     number of reactor channels.
* NB      number of fuel bundles per channel.
* NCOMB   number of combustion zones.
* NX      number of elements along x-axis in fuel map.
* NY      number of elements along y-axis in fuel map.
* NZ      number of elements along z-axis in fuel map.
* LRSCH   flag for the refuelling scheme of channels:
*          =.true. it was read from the input file;
*          =.false. otherwise.
* IMPX    printing index (=0 for no print).
*
*----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAP
      INTEGER NCH,NB,NCOMB,NX,NY,NZ,IMPX
      LOGICAL LRSCH
*----
*  LOCAL VARIABLES
*----
      PARAMETER(IOUT=6)
      INTEGER IVECT(NCOMB,NB),NSCH(NCH),IZONE(NCH),MIX(NX*NY*NZ),
     1 NAMX(NX),NAMY(NY),RSCH(NX,NY),AGLIM,CHR(NB)
      REAL BVAL(NCOMB),DELT(NB),B0(NB),B1(NB),SHAP(NCH,NB),
     1 BURN0(NCH,NB),BURN1(NCH,NB)
      CHARACTER TEXT*12,CHANY*2,FORM1*14,FORM2*14,SHU*3
      LOGICAL LAXSH
*----
*  RECOVER INFORMATION
*----
      CALL LCMLEN(IPMAP,'REF-SCHEME',LENG1,ITYP)
      CALL LCMLEN(IPMAP,'BURN-AVG',LENG2,ITYP)
      IF((LENG1.EQ.0).OR.(LENG2.EQ.0))GOTO 100
      CALL LCMLEN(IPMAP,'AX-SHAPE',LENG3,ITYP)
      IF(LENG3.EQ.0) THEN
*       INITIAL FLAT AXIAL-SHAPE
        IF(IMPX.GT.0)WRITE(IOUT,1000)
        SHAP(:NCH,:NB)=1.0/NB
        CALL LCMPUT(IPMAP,'AX-SHAPE',NCH*NB,2,SHAP)
      ELSE
        CALL LCMGET(IPMAP,'AX-SHAPE',SHAP)
      ENDIF
      CALL LCMGET(IPMAP,'REF-VECTOR',IVECT)
      CALL LCMGET(IPMAP,'REF-SCHEME',NSCH)
      CALL LCMGET(IPMAP,'BURN-AVG',BVAL)
      CALL LCMGET(IPMAP,'B-ZONE',IZONE)
      CALL LCMGET(IPMAP,'BMIX',MIX)
      BURN0(:NCH,:NB)=0.0
      BURN1(:NCH,:NB)=0.0
      LAXSH=.FALSE.
      IF(IMPX.GT.2)WRITE(IOUT,1004)
*----
*  COMPUTE FIRST BURNUP LIMITS
*----
      ICH=0
      DO 70 IEL=1,NX*NY
      IF(MIX(IEL).EQ.0) GOTO 70
      ICH=ICH+1
      IBSH=ABS(NSCH(ICH))
      SHU=' NO'
      DO IB=1,NB
      DELT(IB)=IBSH*BVAL(IZONE(ICH))*SHAP(ICH,IB)
      B0(IB)=0.
      B1(IB)=0.
*     Axial Shuffling detection
      IF(IVECT(IZONE(ICH),IB).GT.IB)THEN
        LAXSH=.TRUE.
        SHU='YES'
      ENDIF
      ENDDO
*     Burnup attribution with axial Shuffling
      IF(LAXSH)THEN 
        AGLIM=INT(NB/IBSH)+1
        CHR(:NB)=AGLIM
*       Two loops on bundle cycles (IA) and number of bundles (IB)
        DO 45 IA=0,AGLIM-1
        DO 40 IB=1,NB
*       Index ordering
        IF (NSCH(ICH).LT.0) THEN
          KK=NB-IB+1
          KV=NB-IVECT(IZONE(ICH),IB)+1
        ELSE
          KK=IB
          KV=IVECT(IZONE(ICH),IB)
        ENDIF
*       New fuel
        IF(IVECT(IZONE(ICH),IB).EQ.0)THEN
            CHR(IB)=0
            B0(KK)=0.
            B1(KK)=DELT(KK)
        ELSE 
*         Compute new burnup if previous bundle cycle done
          IF(CHR(IVECT(IZONE(ICH),IB)).EQ.(IA-1))THEN
            CHR(IB)=IA
            B0(KK)=B1(KV)
            B1(KK)=DELT(KK)+B1(KV)
          ENDIF
        ENDIF
   40   CONTINUE
   45   CONTINUE
*     Burnup attribution without axial Shuffling 
*     One loop on number of bundles (IB)
      ELSE
*       NEGATIVE DIRECTION
        IF(NSCH(ICH).LT.0)THEN
          DO 50 IB=1,NB
          KK=NB-IB+1
          KA=NB-IVECT(IZONE(ICH),IB)+1
          IF(IVECT(IZONE(ICH),IB).LE.0)THEN
            B0(KK)=0.
          ELSE
            B0(KK)=B1(KA)
          ENDIF
          B1(KK)=B0(KK)+DELT(KK)
   50     CONTINUE
*       POSITIVE DIRECTION
        ELSE
          DO 60 IB=1,NB
          IF(IVECT(IZONE(ICH),IB).LE.0)THEN
            B0(IB)=0.
          ELSE
            B0(IB)=B1(IVECT(IZONE(ICH),IB))
          ENDIF
          B1(IB)=B0(IB)+DELT(IB)
   60     CONTINUE
        ENDIF
      ENDIF
      DO IB=1,NB
        BURN0(ICH,IB)=B0(IB)
        BURN1(ICH,IB)=B1(IB)
      ENDDO
      IF(IMPX.GE.3) THEN
*       CHECK BURNUP LIMITS
        WRITE(TEXT,'(A9,I3.3)')'CHANNEL #',ICH
        WRITE(IOUT,1001)TEXT,NSCH(ICH),IZONE(ICH),SHU
        WRITE(IOUT,1002)'B0',(B0(IB),IB=1,NB)
        WRITE(IOUT,1002)'B1',(B1(IB),IB=1,NB)
      ENDIF
*     Reset shuffling for next channel
      LAXSH=.FALSE.
   70 CONTINUE
      CALL LCMPUT(IPMAP,'BURN-BEG',NB*NCH,2,BURN0)
      CALL LCMPUT(IPMAP,'BURN-END',NB*NCH,2,BURN1)
      IF((.NOT.LRSCH).OR.(IMPX.LT.2))GOTO 100
*----
*  PRINT CHANNELS REFUELLING SCHEMES
*----
      WRITE(FORM1,'(A4,I2,A8)')'(A4,',NX,'(A3,1X))'
      WRITE(FORM2,'(A4,I2,A8)')'(A2,',NX,'(I3,1X))'
      CALL LCMGET(IPMAP,'XNAME',NAMX)
      CALL LCMGET(IPMAP,'YNAME',NAMY)
      RSCH(:NX,:NY)=0
      WRITE(IOUT,1003)
      IEL=0
      ICH=0
      DO 85 J=1,NY
      DO 80 I=1,NX
      IEL=IEL+1
      IF(MIX(IEL).EQ.0) GOTO 80
      ICH=ICH+1
      RSCH(I,J)=NSCH(ICH)
   80 CONTINUE
   85 CONTINUE
      WRITE(IOUT,FORM1)' ',(NAMX(I),I=1,NX)
      WRITE(IOUT,*)' '
      DO 90 J=1,NY
      WRITE(CHANY,'(A2)') (NAMY(J))
      IF(INDEX(CHANY,'-').EQ.1) GOTO 90
      WRITE(IOUT,FORM2)CHANY,(RSCH(I,J),I=1,NX)
   90 CONTINUE
  100 RETURN
*
 1000 FORMAT(/1X,'INITIALIZING THE FLAT AXIAL POWER-SHAPE'/
     1 1X,'COMPUTING THE FIRST BURNUP LIMITS PER EACH CHANNEL'/)
 1001 FORMAT(/10X,
     1 A12,10X,'REFUELLING SCHEME:',I3,10X,'ZONE-INDEX:',I3,10X,
     2 'SHUFFLING: ',A3)
 1002 FORMAT(A3,12(F8.1,1X))
 1003 FORMAT(//20X,'** CHANNELS REFUELLING SCHEMES **'/)
 1004 FORMAT(/20X,'** FIRST BURNUP LIMITS PER EACH CHANNEL **'/)
      END
