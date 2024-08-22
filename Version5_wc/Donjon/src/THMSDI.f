*DECK THMSDI
      FUNCTION THMSDI(T2K,T1K,FTP,IFRCDI,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the thermal conductivity integral of UOX or MOX fuel.
*
*Copyright:
* Copyright (C) 2024 Ecole Polytechnique de Montreal.
*
*Author(s): 
* C. Garrido
*
*Parameters: input
* T2K     final temperature in Kelvin.
* T1K     initial temperature in Kelvin.
* FTP     tpdata object with correlations to obtain properties of molten salt.
* IFRCDI  flag indicating if average approximation is forced during
*         fuel conductivity evaluation (0=default/1=average
*         approximation forced).
* IMPX    printing index (=0 for no print).
*
*Parameters: output
* THMSDI  thermal conductivity integral in Watt/m/K.
*
*Reference:
* A. Poncot, "Assimilation de donnees pour la dynamique du xenon dans
* les coeurs de centrale nucleaire", Ph.D Thesis, Universite de
* Toulouse, France, 2008.
*
*-----------------------------------------------------------------------
*
      USE t_saltdata
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(tpdata) FTP
      INTEGER IFRCDI,IMPX
      REAL T1K,T2K,THMSDI
*----
*  LOCAL VARIABLES
* NPAS    number of rectangles in the quadrature
* DT      rectangle width
* T2T1    temperature difference
* DTMIN   cutoff criterion for selecting the approximation
*----
      INTEGER NPAS,I
      REAL T1,T2,DT,TM,DTMIN,T2T1,TT
      REAL R1,R2,ZKONE,ZMUONE,CPONE,CINT
      DATA NPAS /10/
      DATA DTMIN /10./
*
      IF(MIN(T1K,T2K).LE.0.0) THEN
         CALL XABORT('@THMSDI: NEGATIVE TEMPERATURE.')
      ENDIF
      T1=T1K
      T2=T2K
*
      T2T1 = T2-T1
      DT = T2T1/NPAS
      TM = (T1+T2)/2.0
*     User-given conductivity, as a function of temperature
      IF((ABS(T2T1).LT.DTMIN).OR.(IFRCDI.EQ.1)) THEN
*        Use the average value approximation
         CALL THMSPT(FTP,TM,R1,R2,ZKONE,ZMUONE,CPONE,IMPX)
         THMSDI=ZKONE
      ELSE
*        Use the rectangle quadrature approximation
         TT=T1-DT*0.5
         CINT=0.
         DO I=1,NPAS
            TT=TT+DT
            CALL THMSPT(FTP,TT,R1,R2,ZKONE,ZMUONE,CPONE,IMPX)
            CINT=CINT + ZKONE
         ENDDO
         THMSDI=CINT/NPAS
      ENDIF
      RETURN
      END
