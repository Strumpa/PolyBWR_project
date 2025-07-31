*DECK MCGDSD
      SUBROUTINE MCGDSD(NSEG,NSUB,NMU,LPS,NFUNL,NANGL,NGEFF,WEI2D,
     1                  TRHAR,H2D,ZMU,WZMU,KANGL,NOMCEL,NZON,NFI,NREG,
     2                  NDIM,M,IS,JS,PJJ,PSJ,LPJJAN,NPJJM,PJJIND,SIGAL,
     3                  MUST,PHI1,PHI2,PJJX,PJJY,PJJZ,PJJXI,PJJYI,
     4                  PJJZI,CAZ0,PSJX,PSJY,PSJZ)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of the PJJ and PSJ as well as directional values for
* TIBERE.
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
* NSEG    number of elements in the current track.
* NSUB    number of subtracks in the current track.
* NMU     order of the polar quadrature set.
* LPS     first dimension of PSJ.
* NFUNL   number of moments of the flux (in 2D : NFUNL=NANI*(NANI+1)/2).
* NANGL   number of tracking angles in the plane.
* NGEFF   number of energy groups to process.
* NFI     total number of volumes for which specific values
*         of the neutron flux and reactions rates are required.
* NREG    number of volumes for which specific values
*         of the neutron flux and reactions rates are required.
* NDIM    number of dimensions for the geometry.
* M       number of material mixtures.
* IS      arrays for surfaces neighbors.
* JS      JS(IS(ISOUT)+1:IS(ISOUT+1)) give the neighboring regions to
*         surface ISOUT.
* NZON    index-number of the mixture type assigned to each volume.
* TRHAR   spherical harmonics components for this angle in the plane.
* WEI2D   track weight.
* KANGL   track direction indices.
* NOMCEL  integer tracking elements.
* H2D     real tracking elements.
* ZMU     polar quadrature set.
* WZMU    polar quadrature set.
* LPJJAN  flag for the calculation of anisotropic moments of the pjj.
* NPJJM   number of pjj modes to store for LPJJAN option.
* PJJIND  index of the modes for LPJJAN option.
* SIGAL   albedos and total cross sections array.
* MUST    polar index in TRHAR for 3D geometry.
* CAZ0    cosines of the tracking polar angles in 3D.
* PHI1    first cosine of the tracking azimuthal angle.
* PHI2    second cosine of the tracking azimuthal angle.
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
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NGEFF,NSEG,NSUB,NMU,LPS,NFUNL,NANGL,KANGL(NSUB),
     1 NOMCEL(NSEG),NZON(NFI),NFI,M,NREG,NDIM,IS(NFI-NREG+1),JS(LPS),
     2 NPJJM,PJJIND(NPJJM,2),MUST
      DOUBLE PRECISION WEI2D,ZZZ,H2D(NSEG)
      REAL TRHAR(NMU,NFUNL,NANGL),ZMU(NMU),WZMU(NMU),PSJ(LPS,NGEFF),
     1 SIGAL(-6:M,NGEFF),PSJX(LPS,NGEFF),PSJY(LPS,NGEFF),
     2 PSJZ(LPS,NGEFF)
      DOUBLE PRECISION PJJ(NREG,NPJJM,NGEFF),PHI1,PHI2,CAZ0
      DOUBLE PRECISION PJJX(NREG,NPJJM,NGEFF),PJJY(NREG,NPJJM,NGEFF),
     > PJJZ(NREG,NPJJM,NGEFF),PJJXI(NREG,NPJJM,NGEFF),
     > PJJYI(NREG,NPJJM,NGEFF),PJJZI(NREG,NPJJM,NGEFF)
      LOGICAL LPJJAN
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION W,OMEGA2(3)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: HG
*----
*  CALCULATION OF COEFFICIENTS
*----
      IF(NDIM.EQ.3) THEN
*     3D calculation -> no loop over a polar angle   
         DO II=1,NGEFF
*           MCGDSCB: Step-Characteristics Scheme with Tabulated
*                    exponentials
            OMEGA2(3)=CAZ0*CAZ0
            ZZZ=1.0D0/SQRT(1.0D0-OMEGA2(3))
            OMEGA2(1)=(PHI1/ZZZ)**2
            OMEGA2(2)=(PHI2/ZZZ)**2
            W=WEI2D
            CALL MCGDSCB(M,NSEG,NSUB,LPS,IS,JS,H2D(1),KANGL,NOMCEL,NZON,
     1           SIGAL(0,II),W,NFI,NREG,PJJ(1,1,II),PSJ(1,II),MUST,NMU,
     2           NFUNL,NANGL,NPJJM,TRHAR,LPJJAN,PJJIND,OMEGA2,
     3           PJJX(1,1,II),PJJY(1,1,II),PJJZ(1,1,II),PJJXI(1,1,II),
     4           PJJYI(1,1,II),PJJZI(1,1,II),PSJX(1,II),PSJY(1,II),
     5           PSJZ(1,II))
         ENDDO
      ELSE
*     2D calculation -> loop over the polar angle
         ALLOCATE(HG(NSEG))
         DO IMU=1,NMU
            OMEGA2(1)=(PHI1/ZMU(IMU))**2
            OMEGA2(2)=(PHI2/ZMU(IMU))**2
            OMEGA2(3)=(1.0-1.0/ZMU(IMU)**2)
            ZMUI=ZMU(IMU)
            W=WEI2D*WZMU(IMU)
            DO I=1,NSEG
               IF(NZON(NOMCEL(I)).GE.0) THEN
                  HG(I)=H2D(I)*ZMUI
               ENDIF
            ENDDO            
            DO II=1,NGEFF
               CALL MCGDSCB(M,NSEG,NSUB,LPS,IS,JS,HG(1),KANGL,NOMCEL,
     1              NZON,SIGAL(0,II),W,NFI,NREG,PJJ(1,1,II),PSJ(1,II),
     2              IMU,NMU,NFUNL,NANGL,NPJJM,TRHAR,LPJJAN,PJJIND,
     3              OMEGA2,PJJX(1,1,II),PJJY(1,1,II),PJJZ(1,1,II),
     4              PJJXI(1,1,II),PJJYI(1,1,II),PJJZI(1,1,II),
     5              PSJX(1,II),PSJY(1,II),PSJZ(1,II))
            ENDDO
         ENDDO
         DEALLOCATE(HG)
      ENDIF
*
      RETURN
      END
