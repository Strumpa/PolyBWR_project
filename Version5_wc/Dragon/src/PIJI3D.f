*DECK PIJI3D
      SUBROUTINE PIJI3D(NREG,NSOUT,NSLINE,NCOR,
     >                  SWVOID,SIGTAL,WEIGHT,
     >                  SEGLEN,NRSEG,
     >                  STAYIN,GOSOUT,DPR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Integration for general 3D isotropic tracking.
*
*Copyright:
* Copyright (C) 1991 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* NREG    total number of regions.
* NSOUT   number of outer surface.
* NSLINE  number of segemnts on line.
* NCOR    maximum number of corners.
* SWVOID  flag to indicate if there are voids.
* SIGTAL  albedo-cross section vector.
* WEIGHT  line weight.
* SEGLEN  length of track.
* NRSEG   region crossed by track.
*
*Parameters: output
* DPR     CP matrix. 
*
*Parameters: scratch
* STAYIN  stay-in zone probability.
* GOSOUT  goes-out zone probability.
*
*-----------------------------------------------------------------------
*
      IMPLICIT         NONE
*----
* VARIABLES
*----
      INTEGER          NREG,NSOUT,NSLINE,NCOR,NRSEG(NSLINE)
      LOGICAL          SWVOID
      REAL             SIGTAL(-NSOUT:NREG)
      DOUBLE PRECISION WEIGHT,SEGLEN(NSLINE),STAYIN(NSLINE),
     >                 GOSOUT(NSLINE)
      DOUBLE PRECISION DPR(-NSOUT:NREG,-NSOUT:NREG)
*----
*  Local variables
*----
      INTEGER          IL,JL,NOIL
      REAL             ZERO, ONE, HALF
      DOUBLE PRECISION XSIL, PRODUC, DSCBEG, DSCEND, ZCOR, ZCOR2
      INTEGER          ICSEG,JCSEG,ISD,ISF
      PARAMETER       (ZERO=0.0E0, ONE=1.0E0, HALF=0.5E0 )
      REAL             SIXT,CUTEXP
      PARAMETER       (SIXT=HALF/3.0,CUTEXP=0.02)
      DOUBLE PRECISION EXSIL,XSIL2
*----
*  Process track required
*----
      IF( NCOR.EQ.1 )THEN
*
*1)   ONLY ONE EXTERNAL SURFACE AT END --------------------------------
      ISD=NRSEG(1)
      ISF=NRSEG(NSLINE)
      IF( SWVOID )THEN
         PRODUC= WEIGHT
*        PII CALCULATION AND ESCAPE
         DO 40 IL = 1,NSLINE-2
            ICSEG=IL+1
            NOIL  = NRSEG(ICSEG)
            XSIL  = SIGTAL(NOIL)*SEGLEN(ICSEG)
            IF( XSIL.EQ.ZERO )THEN
               GOSOUT(IL)= ONE
               STAYIN(IL)= SEGLEN(ICSEG)
               DPR(NOIL,NOIL)= DPR(NOIL,NOIL)
     >                       + HALF*WEIGHT*SEGLEN(ICSEG)*SEGLEN(ICSEG)
            ELSE IF(XSIL .LT. CUTEXP) THEN
               XSIL2=XSIL*XSIL
               EXSIL=XSIL2*(HALF-SIXT*XSIL)
               STAYIN(IL)=XSIL-EXSIL
               GOSOUT(IL)=ONE-STAYIN(IL)
               PRODUC= PRODUC * GOSOUT(IL)
               DPR(NOIL,NOIL)= DPR(NOIL,NOIL) + WEIGHT*EXSIL
            ELSE
               EXSIL=EXP( - XSIL )
               STAYIN(IL)= ONE - EXSIL
               GOSOUT(IL)= EXSIL
               PRODUC= PRODUC * GOSOUT(IL)
               DPR(NOIL,NOIL)= DPR(NOIL,NOIL) + WEIGHT*(XSIL-STAYIN(IL))
            ENDIF
   40    CONTINUE
*        PIJ CALCULATION
         DSCBEG= WEIGHT
         DO 60 IL = 1, NSLINE-2
            ICSEG=IL+1
            NOIL  = NRSEG(ICSEG)
            DSCEND= WEIGHT*STAYIN(IL)
            DO 50 JL  = IL+1, NSLINE-2
               JCSEG=JL+1
               DPR(NRSEG(JCSEG),NOIL)= 
     >         DPR(NRSEG(JCSEG),NOIL) + STAYIN(JL)*DSCEND
               DSCEND= DSCEND*GOSOUT(JL)
   50       CONTINUE
*           PIS CALCULATION
            DPR(ISD,NOIL)= DPR(ISD,NOIL)+DSCBEG*STAYIN(IL)
            DPR(ISF,NOIL)= DPR(ISF,NOIL)+DSCEND
            DSCBEG= DSCBEG * GOSOUT(IL)
   60    CONTINUE
*        PSS CALCULATION
         DPR(ISD,ISF)= DPR(ISD,ISF) + PRODUC
      ELSE
*
*1.2) NO VOID REGION
         PRODUC= WEIGHT
*        PII CALCULATION AND ESCAPE
         DO 140 IL = 1,NSLINE-2
            ICSEG=IL+1
            NOIL  = NRSEG(ICSEG)
            XSIL  = SIGTAL(NOIL)*SEGLEN(ICSEG)
            IF(XSIL .LT. CUTEXP) THEN
               XSIL2=XSIL*XSIL
               EXSIL=XSIL2*(HALF-SIXT*XSIL)
               STAYIN(IL)=XSIL-EXSIL
               GOSOUT(IL)=ONE-STAYIN(IL)
               PRODUC= PRODUC * GOSOUT(IL)
               DPR(NOIL,NOIL)= DPR(NOIL,NOIL) + WEIGHT*EXSIL
            ELSE
               EXSIL=EXP( - XSIL )
               STAYIN(IL)= ONE - EXSIL
               GOSOUT(IL)= EXSIL
               PRODUC= PRODUC * GOSOUT(IL)
               DPR(NOIL,NOIL)= DPR(NOIL,NOIL) + WEIGHT*(XSIL-STAYIN(IL))
            ENDIF
  140    CONTINUE
*        PIJ CALCULATION
         DSCBEG= WEIGHT
         DO 160 IL = 1, NSLINE-2
            ICSEG=IL+1
            NOIL  = NRSEG(ICSEG)
            DSCEND= WEIGHT*STAYIN(IL)
            DO 150 JL  = IL+1, NSLINE-2
               JCSEG=JL+1
               DPR(NRSEG(JCSEG),NOIL)= 
     >         DPR(NRSEG(JCSEG),NOIL)+ STAYIN(JL)*DSCEND
               DSCEND= DSCEND*GOSOUT(JL)
  150       CONTINUE
*           PIS CALCULATION
            DPR(ISD,NOIL)= DPR(ISD,NOIL)+DSCBEG*STAYIN(IL)
            DPR(ISF,NOIL)= DPR(ISF,NOIL)+DSCEND
            DSCBEG= DSCBEG * GOSOUT(IL)
  160    CONTINUE
*        PSS CALCULATION
         DPR(ISD,ISF)= DPR(ISD,ISF) + PRODUC
      ENDIF
      ELSE
*
*2)   MORE THAN ONE SURFACE PER LINE ----------------------------------
      ZCOR= 1./FLOAT(NCOR)
      ZCOR2= ZCOR*ZCOR
      IF( SWVOID )THEN
*
*2.1) VOIDS ARE POSSIBLE
         PRODUC= WEIGHT*ZCOR2
*        PII CALCULATION AND ESCAPE
         DO 240 IL = 1,NSLINE-2*NCOR
            ICSEG=IL+NCOR
            NOIL  = NRSEG(ICSEG)
            XSIL  = SIGTAL(NOIL)*SEGLEN(ICSEG)
            IF( XSIL.EQ.ZERO )THEN
               GOSOUT(IL)= ONE
               STAYIN(IL)= SEGLEN(ICSEG)
               DPR(NOIL,NOIL)= DPR(NOIL,NOIL)
     >                       + HALF*WEIGHT*SEGLEN(ICSEG)*SEGLEN(ICSEG)
            ELSE IF(XSIL .LT. CUTEXP) THEN
               XSIL2=XSIL*XSIL
               STAYIN(IL)=XSIL-XSIL2*(HALF-SIXT*XSIL)
               GOSOUT(IL)=ONE-STAYIN(IL)
               PRODUC= PRODUC * GOSOUT(IL)
               DPR(NOIL,NOIL)= DPR(NOIL,NOIL)
     >          + WEIGHT*XSIL2*(HALF-SIXT*XSIL)
            ELSE
               GOSOUT(IL)= EXP( - XSIL )
               STAYIN(IL)= (ONE - GOSOUT(IL))
               PRODUC= PRODUC * GOSOUT(IL)
               DPR(NOIL,NOIL)= DPR(NOIL,NOIL) + WEIGHT*(XSIL-STAYIN(IL))
            ENDIF
  240    CONTINUE
*        PIJ CALCULATION
         DSCBEG= WEIGHT*ZCOR
         DO 260 IL = 1, NSLINE-2*NCOR
            ICSEG=IL+NCOR
            NOIL  = NRSEG(ICSEG)
            DSCEND= WEIGHT*STAYIN(IL)
            DO 250 JL  = IL+1, NSLINE-2*NCOR
               JCSEG=JL+NCOR
               DPR(NRSEG(JCSEG),NOIL)= 
     >         DPR(NRSEG(JCSEG),NOIL)+ STAYIN(JL)*DSCEND
               DSCEND= DSCEND*GOSOUT(JL)
  250       CONTINUE
*           PIS CALCULATION
            DO 261 JL = 1, NCOR
               ISD=NRSEG(JL)
               ISF=NRSEG(NSLINE-NCOR+JL)
               DPR(ISD,NOIL)= DPR(ISD,NOIL)+DSCBEG*STAYIN(IL)
               DPR(ISF,NOIL)= DPR(ISF,NOIL)+DSCEND*ZCOR
  261       CONTINUE
            DSCBEG= DSCBEG*GOSOUT(IL)
  260    CONTINUE
*        PSS CALCULATION
         DO 270 IL = 1, NCOR
         ISD=NRSEG(IL)
         DO 265 JL = 1, NCOR
            ISF=NRSEG(NSLINE-NCOR+JL)
            DPR(ISD,ISF)= DPR(ISD,ISF) + PRODUC
  265    CONTINUE
  270    CONTINUE
      ELSE
*
*2.2) NO VOID REGION
         PRODUC= WEIGHT*ZCOR2
*        PII CALCULATION AND ESCAPE
         DO 340 IL = 1,NSLINE-2*NCOR
            ICSEG=IL+NCOR
            NOIL  = NRSEG(ICSEG)
            XSIL  = SIGTAL(NOIL)*SEGLEN(ICSEG)
            IF(XSIL .LT. CUTEXP) THEN
               XSIL2=XSIL*XSIL
               STAYIN(IL)=XSIL-XSIL2*(HALF-SIXT*XSIL)
               GOSOUT(IL)=ONE-STAYIN(IL)
               PRODUC= PRODUC * GOSOUT(IL)
               DPR(NOIL,NOIL)= DPR(NOIL,NOIL)
     >          + WEIGHT*XSIL2*(HALF-SIXT*XSIL)
            ELSE
               GOSOUT(IL)= EXP( - XSIL )
               STAYIN(IL)= (ONE - GOSOUT(IL))
               PRODUC= PRODUC * GOSOUT(IL)
               DPR(NOIL,NOIL)= DPR(NOIL,NOIL) + WEIGHT*(XSIL-STAYIN(IL))
            ENDIF
  340    CONTINUE
*        PIJ CALCULATION
         DSCBEG= WEIGHT*ZCOR
         DO 360 IL = 1, NSLINE-2*NCOR
            ICSEG=IL+NCOR
            NOIL  = NRSEG(ICSEG)
            DSCEND= WEIGHT*STAYIN(IL)
            DO 350 JL  = IL+1, NSLINE-2*NCOR
               JCSEG=JL+NCOR
               DPR(NRSEG(JCSEG),NOIL)=
     >         DPR(NRSEG(JCSEG),NOIL)+ STAYIN(JL)*DSCEND
               DSCEND= DSCEND*GOSOUT(JL)
  350       CONTINUE
*           PIS CALCULATION
            DO 361 JL = 1, NCOR
               ISD=NRSEG(JL)
               ISF=NRSEG(NSLINE-NCOR+JL)
               DPR(ISD,NOIL)= DPR(ISD,NOIL)+DSCBEG*STAYIN(IL)
               DPR(ISF,NOIL)= DPR(ISF,NOIL)+DSCEND*ZCOR
  361       CONTINUE
            DSCBEG= DSCBEG * GOSOUT(IL)
  360    CONTINUE
*        PSS CALCULATION
         DO 370 IL = 1, NCOR
         ISD=NRSEG(IL)
         DO 365 JL = 1, NCOR
            ISF=NRSEG(NSLINE-NCOR+JL)
            DPR(ISD,ISF)= DPR(ISD,ISF) + PRODUC
  365    CONTINUE
  370    CONTINUE
      ENDIF
      ENDIF
      RETURN
      END
