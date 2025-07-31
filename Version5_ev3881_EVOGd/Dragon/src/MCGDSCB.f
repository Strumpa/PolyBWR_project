*DECK MCGDSCB
      SUBROUTINE MCGDSCB(M,NSEG,NSUB,LPS,IS,JS,H,KANGL,NOM,NZON,TR,W,
     1                   NFI,NREG,PJJ,PSJ,IMU,NMU,NFUNL,NANGL,NPJJM,
     2                   TRHAR,LPJJAN,PJJIND,OMEGA2,PJJX,PJJY,PJJZ,
     3                   PJJXI,PJJYI,PJJZI,PSJX,PSJY,PSJZ)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of contribution in PJJ and PSJ coefficients on one track,
* as well as directional values for TIBERE.
* Step-Characteristics scheme with tabulated exponential calls.
*
*Copyright:
* Copyright (C) 2019 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): S. Musongela
*
*Parameters: input
* LPS     dimension of PSJX, PSJY and PSJZ.
* M       number of material mixtures.
* NSEG    number of elements for this track.
* NSUB    number of subtracks for this track.
* IS      arrays for surfaces neighbors.
* JS      JS(IS(ISOUT)+1:IS(ISOUT+1)) give the neighboring regions to
*         surface ISOUT.
* H       real tracking elements.
* KANGL   track direction indices.
* NOM     integer tracking elements.
* NZON    index-number of the mixture type assigned to each volume.
* TR      macroscopic total cross section.
* W       weight associated with this track.
* NFI     total number of volumes and surfaces for which specific values
*         of the neutron flux and reactions rates are required.
* NREG    number of volumes for which specific values
*         of the neutron flux and reactions rates are required.
* IMU     polar angle index.
* NMU     order of the polar quadrature set.
* NFUNL   number of moments of the flux (in 2D : NFUNL=NANI*(NANI+1)/2).
* NANGL   number of tracking angles in the plane.
* NPJJM   number of pjj modes to store for LPJJAN option.
* TRHAR   spherical harmonics components for each azimuthal angle in
*         the plane.
* LPJJAN  flag for the calculation of anisotropic moments of the pjj.
* PJJIND  index of the modes for LPJJAN option.
* OMEGA2  square x, y and z-component of the direction Omega for 2D
*         geometry.
*
*Parameters: input/output
* PJJ     collision probabilities.
* PJJX    collision probabilities for TIBERE.
* PJJY    collision probabilities for TIBERE.
* PJJZ    collision probabilities for TIBERE.
* PJJXI   collision probabilities for TIBERE.
* PJJYI   collision probabilities for TIBERE.
* PJJZI   collision probabilities for TIBERE.
* PSJ     escape probabilities.
* PSJX    escape probabilities for TIBERE.
* PSJY    escape probabilities for TIBERE.
* PSJZ    escape probabilities for TIBERE.
*
*-----------------------------------------------------------------------
*
      IMPLICIT NONE
*---
* SUBROUTINE ARGUMENTS
*---
      INTEGER M,NSEG,NSUB,NFI,NREG,LPS,IS(NFI-NREG+1),JS(LPS),NZON(NFI),
     1 KANGL(NSUB),NOM(NSEG),IMU,NMU,NFUNL,NANGL,NPJJM,PJJIND(NPJJM,2)
      REAL TR(0:M),PSJ(LPS),TRHAR(NMU,NFUNL,NANGL)
      REAL PSJX(LPS),PSJY(LPS),PSJZ(LPS)
      DOUBLE PRECISION W,H(NSUB),PJJ(NREG,NPJJM),OMEGA2(3)
      DOUBLE PRECISION PJJX(NREG,NPJJM),PJJY(NREG,NPJJM),
     1 PJJZ(NREG,NPJJM),PJJXI(NREG,NPJJM),PJJYI(NREG,NPJJM),
     2 PJJZI(NREG,NPJJM)
      LOGICAL LPJJAN
*---
* LOCAL VARIABLES
*---
      DOUBLE PRECISION TAUDMIN
      PARAMETER(TAUDMIN=2.D-2)
      INTEGER I,J,NOMI,IC,IC0,NZI,NOMJ,IMOD,INU,INUP,IANG,ISUB
      DOUBLE PRECISION TRI,TRJ,TAU,EXPT,HJD,HID,TAUD,TAUD3,TAUD4,TAUD5,
     1 EXPTD,TEMPD
      LOGICAL LNEW
*     tabulated exponential common block
      REAL             E0, E1, PAS1, DX1, XLIM1
      INTEGER          MEX1, LAU
      PARAMETER      ( MEX1=7936 )
      COMMON /EXP1/ E0(0:MEX1),E1(0:MEX1),PAS1,DX1,XLIM1
*
      ISUB=0
      LNEW=.TRUE.
      IANG=KANGL(1)
      DO I=1,NSEG
         NOMI=NOM(I)
         NZI=NZON(NOMI)
         IF(NZI.LT.0) THEN
*        Boundary Condition
            LNEW=.TRUE.
            IF(LPS.GT.0) THEN
*           SCR for a non-cyclic tracking
               IF(I.EQ.1) THEN
                  J=I+1
               ELSE !! I.EQ.NSEG
                  J=I-1
               ENDIF
               NOMJ=NOM(J)
               IC=0
               DO IC0=IS(NOMI-NREG)+1,IS(NOMI-NREG+1)
                  IC=IC0
                  IF(JS(IC0).EQ.NOMJ) GOTO 10
               ENDDO
               CALL XABORT('MCGDSCB: UNABLE TO SET IC.')
 10            HJD=H(J)
               TRJ=TR(NZON(NOMJ))
               TAU=HJD*TRJ
               IF(TAU.GE.XLIM1) THEN
                  EXPT=1.0D0/TRJ
               ELSE
                  LAU=INT(TAU*PAS1)
                  EXPT=HJD*(E0(LAU)+E1(LAU)*TAU)
               ENDIF
               PSJ(IC)=PSJ(IC)+REAL(W*EXPT)
               PSJX(IC)=PSJX(IC)+REAL(W*EXPT*3.0*OMEGA2(1))
               PSJY(IC)=PSJY(IC)+REAL(W*EXPT*3.0*OMEGA2(2))
               PSJZ(IC)=PSJZ(IC)+REAL(W*EXPT*3.0*OMEGA2(3))
            ENDIF
         ELSE
*        this cell is a volume
            IF(LNEW) THEN
               ISUB=ISUB+1
               IF(ISUB.GT.NSUB) CALL XABORT('MCGDSCB: NSUB OVERFLOW.')
               LNEW=.FALSE.
               IANG=KANGL(ISUB)
               IF(IANG.GT.NANGL) CALL XABORT('MCGDSCB: NANGL OVERFLOW.')
            ENDIF
            TRI=TR(NZI)
            HID=H(I)
            TAUD=HID*TRI
            TAU=REAL(TAUD)
            IF(TAUD.LE.TAUDMIN) THEN
*           expansion in Taylor serie in O(TAUD^3)
               TAUD3=TAUD/3.D0
               TAUD4=0.125D0*TAUD
               TAUD5=0.2D0*TAUD
               EXPTD=HID*(0.5D0-TAUD3*(0.5D0-TAUD4*(1.D0-TAUD5)))
            ELSE
               IF(TAU.GE.XLIM1) THEN
*              Out of the table range
                  EXPTD=(1.D0-1.D0/TAUD)/DBLE(TRI)
               ELSE
*              Linear interpolation in table of (1-exp(-x))/x
                  LAU=INT(TAU*PAS1)
                  EXPTD=(1.D0-DBLE(E0(LAU)+E1(LAU)*TAU))/DBLE(TRI)
               ENDIF
            ENDIF
            EXPTD=EXPTD*W*HID
            IF(LPJJAN) THEN
               DO IMOD=1,NPJJM
                  INU=PJJIND(IMOD,1)
                  INUP=PJJIND(IMOD,2)
                  TEMPD=DBLE(TRHAR(IMU,INU,IANG))*
     1                  DBLE(TRHAR(IMU,INUP,IANG))
                  PJJ(NOMI,IMOD)=PJJ(NOMI,IMOD)+EXPTD*TEMPD
               ENDDO
            ELSE
               PJJ(NOMI,1)=PJJ(NOMI,1)+EXPTD
               PJJX(NOMI,1)=PJJX(NOMI,1)+EXPTD*3.0*OMEGA2(1)
               PJJY(NOMI,1)=PJJY(NOMI,1)+EXPTD*3.0*OMEGA2(2)
               PJJZ(NOMI,1)=PJJZ(NOMI,1)+EXPTD*3.0*OMEGA2(3)
               PJJXI(NOMI,1)=PJJXI(NOMI,1)+EXPTD*9.0*
     1                       OMEGA2(1)*OMEGA2(1)                
               PJJYI(NOMI,1)=PJJYI(NOMI,1)+EXPTD*9.0*
     1                       OMEGA2(2)*OMEGA2(2)
               PJJZI(NOMI,1)=PJJZI(NOMI,1)+EXPTD*9.0*
     1                       OMEGA2(3)*OMEGA2(3)
            ENDIF
         ENDIF
      ENDDO  
*
      RETURN
      END
