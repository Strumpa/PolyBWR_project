*DECK MCGDS4
      SUBROUTINE MCGDS4(SUBSCH,NSEG,NSUB,NMU,LPS,NFUNL,NANGL,NGEFF,
     1                  WEI2D,KANGL,TRHAR,H2D,ZMU,WZMU,NOMCEL,NZON,NFI,
     2                  NREG,NDIM,M,IS,JS,PJJ,PSJ,LPJJAN,NPJJM,PJJIND,
     3                  SIGAL,MUST)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of the PJJ and PSJ.
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
* SUBSCH  track coefficients calculation subroutine.
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
* IS      arrays for surfaces neighbors
* JS      JS(IS(ISOUT)+1:IS(ISOUT+1)) give the neighboring regions to
*         surface ISOUT.
* NZON    index-number of the mixture type assigned to each volume.
* TRHAR   spherical harmonics components for each angle in the plane.
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
*
*Parameters: input/output
* PJJ     collision probabilities.
* PSJ     escape probabilities.
*
*-----------------------------------------------------------------------
*
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NGEFF,NSEG,NSUB,NMU,LPS,NFUNL,NANGL,KANGL(NSEG),
     1 NOMCEL(NSEG),NZON(NFI),NFI,M,NREG,NDIM,IS(NFI-NREG+1),JS(LPS),
     2 NPJJM,PJJIND(NPJJM,2),MUST
      DOUBLE PRECISION WEI2D,H2D(NSEG)
      REAL TRHAR(NMU,NFUNL,NANGL),ZMU(NMU),WZMU(NMU),PSJ(LPS,NGEFF),
     1 SIGAL(-6:M,NGEFF)
      DOUBLE PRECISION PJJ(NREG,NPJJM,NGEFF)
      LOGICAL LPJJAN
      EXTERNAL SUBSCH
*----
*  LOCAL VARIABLES
*----
      DOUBLE PRECISION W
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: HG
*----
*  CALCULATION OF COEFFICIENTS
*----
      IF(NDIM.EQ.3) THEN
*     3D calculation -> no loop over a polar angle   
         DO II=1,NGEFF
*           MCGDSCA: Step-Characteristics Scheme with Tabulated Exponentials
*           MCGDSCE: Step-Characteristics Scheme with Exact Exponentials
*           MCGDDDF: Diamond-Differencing Scheme
            CALL SUBSCH(M,NSEG,NSUB,LPS,IS,JS,H2D,KANGL,NOMCEL,NZON,
     1           SIGAL(0,II),WEI2D,NFI,NREG,PJJ(1,1,II),PSJ(1,II),MUST,
     2           NMU,NFUNL,NANGL,NPJJM,TRHAR,LPJJAN,PJJIND)
         ENDDO
      ELSE
*     2D calculation -> loop over the polar angle
         ALLOCATE(HG(NSEG))
         DO IMU=1,NMU
            ZMUI=ZMU(IMU)
            W=WEI2D*WZMU(IMU)
            DO I=1,NSEG
               IF(NZON(NOMCEL(I)).GE.0) HG(I)=H2D(I)*ZMUI
            ENDDO
            DO II=1,NGEFF
               CALL SUBSCH(M,NSEG,NSUB,LPS,IS,JS,HG,KANGL,NOMCEL,NZON,
     1              SIGAL(0,II),W,NFI,NREG,PJJ(1,1,II),PSJ(1,II),IMU,
     2              NMU,NFUNL,NANGL,NPJJM,TRHAR,LPJJAN,PJJIND)
            ENDDO
         ENDDO
         DEALLOCATE(HG)
      ENDIF
*
      RETURN
      END
