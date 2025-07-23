*DECK MCGDDDF
      SUBROUTINE MCGDDDF(M,NSEG,NSUB,LPS,IS,JS,H,KANGL,NOM,NZON,TR,W,
     1                   NFI,NREG,PJJ,PSJ,IMU,NMU,NFUNL,NANGL,NPJJM,
     2                   TRHAR,LPJJAN,PJJIND)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of contribution in PJJ and PSJ coefficients on one track.
* Diamond-Differencing scheme.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Le Tellier
*
*Parameters: input
* LPS     dimension of PSJ.
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
*
*Parameters: input/output
* PJJ     collision probabilities.
* PSJ     escape probabilities.
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
      DOUBLE PRECISION W,H(NSEG),PJJ(NREG,NPJJM)
      LOGICAL LPJJAN
*---
* LOCAL VARIABLES
*---
      INTEGER I,J,NOMI,IC,IC0,NZI,NOMJ,IMOD,INU,INUP,IANG,ISUB
      DOUBLE PRECISION TRI,TRJ,TAU,EXPT,HJD,HID,TAUD,EXPTD,TEMPD
      LOGICAL LNEW
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
               CALL XABORT('MCGDDDF: UNABLE TO SET IC.')
 10            HJD=H(J)
               TRJ=TR(NZON(NOMJ))
               TAU=HJD*TRJ
               EXPT=2.0D0*HJD/(2.0D0+TAU)
               PSJ(IC)=PSJ(IC)+REAL(W*EXPT)
            ENDIF
         ELSE
*        this cell is a volume
            IF(LNEW) THEN
               ISUB=ISUB+1
               IF(ISUB.GT.NSUB) CALL XABORT('MCGDDDF: NSUB OVERFLOW.')
               LNEW=.FALSE.
               IANG=KANGL(ISUB)
               IF(IANG.GT.NANGL) CALL XABORT('MCGDDDF: NANGL OVERFLOW.')
            ENDIF
            TRI=TR(NZI)
            HID=H(I)
            TAUD=HID*TRI
            EXPTD=HID/(2.D0+TAUD)
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
            ENDIF
         ENDIF
      ENDDO  
*
      RETURN
      END
