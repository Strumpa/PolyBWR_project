*DECK XDRTA2
      SUBROUTINE XDRTA2
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover the tabulated functions required by the flux solution and
* store them in common blocks.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy and A. Hebert
*
*-----------------------------------------------------------------------
*
***** OUTPUT:   THE FIVE COMMONS OF BICKLEY QUADRATIC TABLES ARE FILLED
*               AND SAVED WITH RESPECTIVE NAMES:
*                  /BICKL1/,/BICKL2/,/BICKL3/,/BICKL4/,/BICKL5/
*
*               A COMMON FOR LINEAR EXPONENTIAL TABLES IS FILLED
*               AND SAVED WITH NAME: /EXP1/
*----
*  BICKLEY FUNCTION COMMONS
*----
      DOUBLE PRECISION DX
      INTEGER MLOG(5)      
      PARAMETER (NBX=600,DX=0.02D0,MLOG=(/30,15,0,0,0/))
      REAL BIV(0:NBX,3,5),XLIMV(5),PASV(5)
      COMMON /BICKL1/BI1(0:NBX),BI11(0:NBX),BI12(0:NBX),PAS1,XLIM1,L1
      COMMON /BICKL2/BI2(0:NBX),BI21(0:NBX),BI22(0:NBX),PAS2,XLIM2,L2
      COMMON /BICKL3/BI3(0:NBX),BI31(0:NBX),BI32(0:NBX),PAS3,XLIM3,L3
      COMMON /BICKL4/BI4(0:NBX),BI41(0:NBX),BI42(0:NBX),PAS4,XLIM4,L4
      COMMON /BICKL5/BI5(0:NBX),BI51(0:NBX),BI52(0:NBX),PAS5,XLIM5,L5
      SAVE /BICKL1/,/BICKL2/,/BICKL3/,/BICKL4/,/BICKL5/
*----
*  EXPONENTIAL COMMONS
*----
      DOUBLE PRECISION DEX
      REAL PARAM(3)
      PARAMETER (NBEX=7936,DEX=1.D0/512.D0)
      COMMON /EXP1/ E10(0:NBEX),E11(0:NBEX),PASE1,DXE1,XLIME1
      COMMON /EXP0/ E00(0:NBEX),E01(0:NBEX),PASE0,DXE0,XLIME0
      SAVE   /EXP1/,/EXP0/
*----
*  CHARGE BICKLEY TABLES INTO COMMON
*----
      CALL XDRKIN(DX,NBX,MLOG,BIV,PASV,XLIMV)
      PAS1=PASV(1)
      PAS2=PASV(2)
      PAS3=PASV(3)
      PAS4=PASV(4)
      PAS5=PASV(5)
      XLIM1=XLIMV(1)
      XLIM2=XLIMV(2)
      XLIM3=XLIMV(3)
      XLIM4=XLIMV(4)
      XLIM5=XLIMV(5)
      L1=MLOG(1)
      L2=MLOG(2)
      L3=MLOG(3)
      L4=MLOG(4)
      L5=MLOG(5)
      BI1(0:NBX)=BIV(0:NBX,1,1)
      BI11(0:NBX)=BIV(0:NBX,2,1)
      BI12(0:NBX)=BIV(0:NBX,3,1)
*
      BI2(0:NBX)=BIV(0:NBX,1,2)
      BI21(0:NBX)=BIV(0:NBX,2,2)
      BI22(0:NBX)=BIV(0:NBX,3,2)
*
      BI3(0:NBX)=BIV(0:NBX,1,3)
      BI31(0:NBX)=BIV(0:NBX,2,3)
      BI32(0:NBX)=BIV(0:NBX,3,3)
*
      BI4(0:NBX)=BIV(0:NBX,1,4)
      BI41(0:NBX)=BIV(0:NBX,2,4)
      BI42(0:NBX)=BIV(0:NBX,3,4)
*
      BI5(0:NBX)=BIV(0:NBX,1,5)
      BI51(0:NBX)=BIV(0:NBX,2,5)
      BI52(0:NBX)=BIV(0:NBX,3,5)
*----
*  CHARGE EXPONENTIAL TABLES INTO COMMON
*----
      CALL XDREXP(DEX,NBEX,PARAM,E00,E01,E10,E11)
      PASE1=PARAM(1)
      DXE1=PARAM(2)
      XLIME1=PARAM(3)
      PASE0=PARAM(1)
      DXE0=PARAM(2)
      XLIME0=PARAM(3)
      RETURN
      END
