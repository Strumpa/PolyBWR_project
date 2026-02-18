*DECK TAVGLM
      SUBROUTINE TAVGLM(NB,SHIFT,BCHAN,PSI,BURN0,BURN1,IVECT,NSCH)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the burnup integration limits for a given channel.
*
*Copyright:
* Copyright (C) 2007 Ecole Polytechnique de Montreal
*
*Author(s): 
* D.Rozon, M.Beaudet, D.Sekki, I. Trancart
*
*Parameters: input
* NB     number of fuel bundles.
* SHIFT  number of bundles to refuel (bundle-shift).
* PSI    axial shape over each bundle.
* NSCH   refuelling scheme of a given channel.
* BCHAN  average exit burnup for a given channel.
* IVECT  refuelling pattern vector for a given channel.
*
*Parameters: output
* BURN0  lower burnup integration limit.
* BURN1  upper burnup integration limit.
*
*Parameters: scratch
* DELT   incremental burnup over each fuel bundle.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NB,SHIFT,NSCH,IVECT(NB),CHR(NB),AGLIM
      REAL BURN0(NB),BURN1(NB),PSI(NB),BCHAN,DELT(NB)
*----
*  LOCAL VARIABLES
*----         
      LOGICAL LAXSH
*----
*  SCRATCH STORAGE ALLOCATION
*----                
*----
*  COMPUTE BURNUP LIMITS
*----
      BURN0(:NB)=0.0
      BURN1(:NB)=0.0
      DELT(:NB)=0.0
      LAXSH=.FALSE.
      DO 10 IB=1,NB
      DELT(IB)=SHIFT*BCHAN*PSI(IB)
      IF(IVECT(IB).GT.IB)THEN
        LAXSH=.TRUE.
      ENDIF
   10 CONTINUE
*     Burnup attribution with axial Shuffling
      IF(LAXSH)THEN 
        AGLIM=INT(NB/SHIFT)+1
        CHR(:NB)=AGLIM
*       Two loops on bundle cycles (IA) and nmake tesumber of bundles (IB)
        DO 25 IA=0,AGLIM-1
        DO 20 IB=1,NB
*       Index ordering
        IF (NSCH.LT.0) THEN
          KK=NB-IB+1
          KV=NB-IVECT(IB)+1
        ELSE
          KK=IB
          KV=IVECT(IB)
        ENDIF
*       New fuel
        IF(IVECT(IB).EQ.0)THEN
          CHR(IB)=0
          BURN0(KK)=0.
          BURN1(KK)=DELT(KK)
        ELSE 
*         Compute new burnup if previous bundle cycle done
          IF(CHR(IVECT(IB)).EQ.(IA-1))THEN
            CHR(IB)=IA
            BURN0(KK)=BURN1(KV)
            BURN1(KK)=DELT(KK)+BURN1(KV)
          ENDIF
        ENDIF
   20   CONTINUE
   25   CONTINUE
*     Burnup attribution without axial Shuffling 
*     One loop on number of bundles (IB)
      ELSE 
*       NEGATIVE DIRECTION
        IF(NSCH.LT.0)THEN
          DO 40 IB=1,NB
          KK=NB-IB+1
          KA=NB-IVECT(IB)+1
          IF(IVECT(IB).LE.0)THEN
            BURN0(KK)=0.
          ELSE
            BURN0(KK)=BURN1(KA)
          ENDIF
          BURN1(KK)=BURN0(KK)+DELT(KK)
   40     CONTINUE
*       POSITIVE DIRECTION
        ELSE
          DO 50 IB=1,NB
          IF(IVECT(IB).LE.0)THEN
            BURN0(IB)=0.
          ELSE
            BURN0(IB)=BURN1(IVECT(IB))
          ENDIF
          BURN1(IB)=BURN0(IB)+DELT(IB)
   50     CONTINUE
        ENDIF
      ENDIF
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      RETURN
      END
