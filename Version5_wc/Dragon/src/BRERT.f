*DECK BRERT
      SUBROUTINE BRERT(IPMAC1,IELEM,ICOL,NG,NL,LX1,NMIX1,IMIX,ICODE,
     1 ISPH,IDIFF,ZKEFF,B2,ENER,XXX1,VOL1,FLX1,DC1,TOT1,CHI1,SIGF1,
     2 SCAT1,JXM,JXP,FHETXM,FHETXP,ADF1,NGET,ADFREF,IMPX)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Implement the 1D DF-RT (Raviart-Thomas) reflector model.
*
*Copyright:
* Copyright (C) 2025 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* IPMAC1  nodal macrolib.
* IELEM   Raviart-Thomas polynomial order.
* ICOL    Raviart-Thomas polynomial integration type.
* NG      number of energy groups.
* NL      Legendre order of TOT1 and SCAT1 arrays (=1 for isotropic
*         scattering in LAB). (NL-1) is the SPN order (if IDIFF>1,
*         NL is an even integer).
* LX1     number of nodes in the reflector model.
* NMIX1   number of mixtures in the nodal calculation.
* IMIX    mix index of each node.
* ICODE   physical albedo index on each side of the domain.
* ISPH    SPH flag (=0: use discontinuity factors; =1: use SPH factors).
* IDIFF   PN calculation option (=0: diffusion theory; =1: SPN theory
*         with 'NTOT1'; =2: SPN theory with 1/(3*D)).
* ZKEFF   effective multiplication factor.
* B2      buckling.
* ENER    energy limits.
* XXX1    spatial mesh.
* VOL1    volumes.
* FLX1    averaged fluxes
* DC1     diffusion coefficients.
* TOT1    total cross sections.
* CHI1    fission spectra.
* SIGF1   nu*fission cross sections.
* SCAT1   scattering P0 cross sections.
* JXM     left boundary currents.
* JXP     right boundary currents.
* FHETXM  left boundary fluxes.
* FHETXP  right boundary fluxes.
* ADF1    assembly discontinuity factors from macrolib.
* NGET    type of NGET normalization if discontinuity factors
*         (=0: simple; =1: imposed ADF on fuel assembly; =2: recover
*         fuel assembly ADF from input macrolib).
* ADFREF  imposed ADF values on fuel assembly side.
* IMPX   edition flag.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPMAC1
      INTEGER IELEM,ICOL,NG,NL,LX1,NMIX1,IMIX(LX1),ICODE(2),ISPH,IDIFF,
     1 NGET,IMPX
      REAL ZKEFF,B2,ENER(NG+1),XXX1(LX1+1),VOL1(NMIX1),FLX1(NMIX1,NG),
     1 DC1(NMIX1,NG),TOT1(NMIX1,NG,NL),CHI1(NMIX1,NG),SIGF1(NMIX1,NG),
     2 SCAT1(NMIX1,NG,NG,NL),JXM(NMIX1,NG),JXP(NMIX1,NG),
     3 FHETXM(NMIX1,NG,NL),FHETXP(NMIX1,NG,NL),ADF1(NMIX1,NG),ADFREF(NG)
*----
*  LOCAL VARIABLES
*----
      PARAMETER (NSTATE=40)
      INTEGER ISTATE(NSTATE)
      CHARACTER CM*2,HADF*8,TEXT12*12
      TYPE(C_PTR) JPMAC1,KPMAC1
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJ,NJJ,IPOS
      REAL, ALLOCATABLE, DIMENSION(:) :: WORK,AFACTOR,BETA,WORK1,WORK2,
     1 VOLTOT
      REAL, ALLOCATABLE, DIMENSION(:,:) :: FDXM,FDXP,WORK3,WORK4,WORK5,
     1 WORK6,WORK7
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: FHOMM,FHOMP
      REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: RCAT
*----
*  SCRATCH STORAGE ALLOCATION
*----
      J_FUEL=0
      NLMAX=1
      IF(IDIFF.GT.0) THEN
        IF(NL.LT.2) CALL XABORT('BRERT: EVEN NL>=2 EXPECTED WITH SPN.')
        NLMAX=NL/2
      ENDIF
      ALLOCATE(FHOMM(NMIX1,NG,NLMAX),FHOMP(NMIX1,NG,NLMAX),
     1 FDXM(NMIX1,NG),FDXP(NMIX1,NG),AFACTOR(NG),BETA(NG),VOLTOT(NMIX1),
     2 WORK1(NG),WORK2(NG),WORK3(NG,NG),WORK6(NG,NLMAX),WORK7(NG,NLMAX))
*----
*  COMPUTE BOUNDARY FLUXES
*----
      FDXM(:NMIX1,:NG)=0.0
      FDXP(:NMIX1,:NG)=0.0
      FHOMM(:NMIX1,:NG,:NLMAX)=0.0
      FHOMP(:NMIX1,:NG,:NLMAX)=0.0
      VOLTOT(:NMIX1)=0.0
      DO I=1,LX1
        IBM=IMIX(I)
        IF(IBM.EQ.0) CYCLE
        WORK6(:NG,:NLMAX)=0.0
        WORK7(:NG,:NLMAX)=0.0
        DELX=XXX1(I+1)-XXX1(I)
        IF(IMPX.GT.0) WRITE(6,'(/15H BRERT: REGION=,I5)') I
        IF(IDIFF.EQ.0) THEN
*         diffusion theory
          ALLOCATE(WORK4(NG,1),WORK5(NG,1))
          WORK1(:NG)=DC1(IBM,:NG)
          WORK4(:NG,1)=JXM(IBM,:NG)
          WORK5(:NG,1)=JXP(IBM,:NG)
          DO IG=1,NG
            IF(SIGF1(IBM,IG).GT.0.0) J_FUEL=I
            WORK2(IG)=TOT1(IBM,IG,1)+B2*DC1(IBM,IG)-SCAT1(IBM,IG,IG,1)
            DO JG=1,NG
              WORK3(IG,JG)=CHI1(IBM,IG)*SIGF1(IBM,JG)/ZKEFF
              IF(JG.NE.IG) WORK3(IG,JG)=WORK3(IG,JG)+SCAT1(IBM,IG,JG,1)
            ENDDO
          ENDDO
          CALL BRERTD(IELEM,ICOL,NG,DELX,WORK1,WORK2,WORK3,WORK4,WORK5,
     1    IMPX,WORK6,WORK7)
          DEALLOCATE(WORK5,WORK4)
        ELSE
*         SPN theory
          ALLOCATE(WORK4(NG,NL/2),WORK5(NG,NL/2),RCAT(NG,NG,NL))
          DO IL=1,NL/2
            WORK4(:NG,IL)=FHETXM(IBM,:NG,2*IL)
            WORK5(:NG,IL)=FHETXP(IBM,:NG,2*IL)
          ENDDO
          RCAT(:NG,:NG,:NL)=0.0
          DO IG=1,NG
            IF(SIGF1(IBM,IG).GT.0.0) J_FUEL=I
            DO JG=1,NG
              RCAT(IG,JG,1)=-CHI1(IBM,IG)*SIGF1(IBM,JG)/ZKEFF
            ENDDO
            RCAT(IG,IG,1)=RCAT(IG,IG,1)+B2*DC1(IBM,IG)
            DO IL=1,NL,2
              RCAT(IG,IG,IL)=RCAT(IG,IG,IL)+TOT1(IBM,IG,IL)
              DO JG=1,NG
                RCAT(IG,JG,IL)=RCAT(IG,JG,IL)-SCAT1(IBM,IG,JG,IL)
              ENDDO
            ENDDO
            DO IL=2,NL,2
              IF(IDIFF.EQ.1) THEN
                DO JG=1,NG
                  RCAT(IG,JG,IL)=RCAT(IG,JG,IL)-SCAT1(IBM,IG,JG,IL)
                ENDDO
              ELSE
                TOT1(IBM,IG,IL)=1.0/(3.0*DC1(IBM,IG))
                SCAT1(IBM,IG,:NG,IL)=0.0
              ENDIF
              RCAT(IG,IG,IL)=RCAT(IG,IG,IL)+TOT1(IBM,IG,IL)
            ENDDO
          ENDDO
          DO IL=1,NL
            RCAT(:NG,:NG,IL)=RCAT(:NG,:NG,IL)*REAL(2*IL-1)
          ENDDO
          CALL BRERTS(IELEM,ICOL,NG,NL,DELX,RCAT,WORK4,WORK5,IMPX,
     1    WORK6,WORK7)
          DEALLOCATE(RCAT,WORK5,WORK4)
        ENDIF
        FHOMM(IBM,:NG,:NLMAX)=FHOMM(IBM,:NG,:NLMAX)+WORK6(:NG,:NLMAX)*
     1  DELX
        FHOMP(IBM,:NG,:NLMAX)=FHOMP(IBM,:NG,:NLMAX)+WORK7(:NG,:NLMAX)*
     1  DELX
        VOLTOT(IBM)=VOLTOT(IBM)+DELX
      ENDDO
      DEALLOCATE(WORK7,WORK6,WORK3,WORK2,WORK1)
      DO IBM=1,NMIX1
        DO IGR=1,NG
          IF(NL.LE.2) THEN
            FDXM(IBM,IGR)=VOLTOT(IBM)*FHETXM(IBM,IGR,1)/FHOMM(IBM,IGR,1)
            FDXP(IBM,IGR)=VOLTOT(IBM)*FHETXP(IBM,IGR,1)/FHOMP(IBM,IGR,1)
          ELSE
            ! Yamamoto formula
            FDXM(IBM,IGR)=VOLTOT(IBM)*(FHETXM(IBM,IGR,1)+2.0*
     1      FHETXM(IBM,IGR,2))/(FHOMM(IBM,IGR,1)+2.0*FHOMM(IBM,IGR,2))
            FDXP(IBM,IGR)=VOLTOT(IBM)*(FHETXP(IBM,IGR,1)+2.0*
     1      FHETXP(IBM,IGR,2))/(FHOMP(IBM,IGR,1)+2.0*FHOMP(IBM,IGR,2))
          ENDIF
        ENDDO
      ENDDO
      IF(IMPX.GT.0) THEN
        WRITE(6,'(/48H BRERT: DISCONTINUITY FACTORS BEFORE NORMALIZATI,
     1  2HON)')
        DO IBM=1,NMIX1
          WRITE(6,'(/9H MIXTURE=,I5)') IBM
          WRITE(6,20) 'FDXM',FDXM(IBM,:NG)
          WRITE(6,20) 'FDXP',FDXP(IBM,:NG)
        ENDDO
      ENDIF
*----
*  COMPUTE ALBEDOS
*----
      IF(ICODE(2).NE.0) THEN
        BETA(:)=0.0
        IBM=IMIX(LX1)
        DO IGR=1,NG
          IF(IBM.EQ.0) CYCLE
          AFACTOR(IGR)=FDXP(IBM,IGR)*JXP(IBM,IGR)/FHETXP(IBM,IGR,1)
          BETA(IGR)=(1.0-2.0*AFACTOR(IGR))/(1.0+2.0*AFACTOR(IGR))
        ENDDO
        IF(IMPX.GT.0) THEN
          WRITE(6,'(/15H BRERT: ALBEDOS)')
          WRITE(6,20) 'BETA',BETA(:NG)
        ENDIF
      ENDIF
*----
*     THE SPH PARAMETERS ARE NOT DEGENERATE IN NON-FUNDAMENTAL MODE
*     CONDITION. THE ONLY SOLUTION CORRESPONDS TO J_FUEL=1
*----
      IF(ISPH.EQ.1) J_FUEL=1
*----
*  NGET NORMALIZATION OF THE DISCONTINUITY FACTORS
*----
      IF(J_FUEL.GT.0) THEN
        IF(NGET.GT.0) THEN
          IBM=IMIX(J_FUEL)
          DO IGR=1,NG
            ! impose the adf on the fuel assembly side
            IF(IBM.EQ.0) CYCLE
            IF(NGET.EQ.1) THEN
              FNORM=ADFREF(IGR)/FDXP(IBM,IGR)
            ELSE
              FNORM=ADF1(IBM,IGR)/FDXP(IBM,IGR)
            ENDIF
            FDXP(IBM,IGR)=FDXP(IBM,IGR)*FNORM
            IF(J_FUEL<LX1) THEN
              IBMP=IMIX(J_FUEL+1)
              IF(IBMP.GT.0) FDXM(IBMP,IGR)=FDXM(IBMP,IGR)*FNORM
            ENDIF
          ENDDO
        ENDIF
        DO J=J_FUEL,1,-1
          IBM=IMIX(J)
          IF(IBM.EQ.0) CYCLE
          DO IGR=1,NG
            IF(J>1) THEN
              IBMM=IMIX(J-1)
              IF(IBMM.GT.0) FDXP(IBMM,IGR)=FDXP(IBMM,IGR)*FDXP(IBM,IGR)/
     1        FDXM(IBM,IGR)
            ENDIF
            FDXM(IBM,IGR)=FDXP(IBM,IGR)
          ENDDO
        ENDDO
      ENDIF
      DO J=J_FUEL+1,LX1
        IBM=IMIX(J)
        IF(IBM.EQ.0) CYCLE
        DO IGR=1,NG
          IF(J<LX1) THEN
            IBMP=IMIX(J+1)
            IF(IBMP.GT.0) FDXM(IBMP,IGR)=FDXM(IBMP,IGR)*FDXM(IBM,IGR)/
     1      FDXP(IBM,IGR)
          ENDIF
          FDXP(IBM,IGR)=FDXM(IBM,IGR)
        ENDDO
      ENDDO
      IF(IMPX.GT.0) THEN
        WRITE(6,'(/48H BRERT: DISCONTINUITY FACTORS AFTER NGET NORMALI,
     1  6HZATION)')
        DO IBM=1,NMIX1
          WRITE(6,'(/9H MIXTURE=,I5)') IBM
          WRITE(6,20) 'FDX',FDXM(IBM,:NG)
        ENDDO
      ENDIF
*----
*  APPLY SPH FACTORS
*----
      IF(ISPH.EQ.1) THEN
        DO IGR=1,NG
          DO IBM=1,NMIX1
            IF(FDXM(IBM,IGR).LE.0) CALL XABORT('BRERT: NEGATIVE SPH F'
     1      //'ACTOR.')
          ENDDO
          DO J=1,LX1
            IBM=IMIX(J)
            IF(IBM.EQ.0) CYCLE
            DC1(IBM,IGR)=DC1(IBM,IGR)/FDXM(IBM,IGR)
            SIGF1(IBM,IGR)=SIGF1(IBM,IGR)/FDXM(IBM,IGR)
            DO IL=1,NL,2
              TOT1(IBM,IGR,IL)=TOT1(IBM,IGR,IL)/FDXM(IBM,IGR)
              DO JGR=1,NG
               SCAT1(IBM,IGR,JGR,IL)=SCAT1(IBM,IGR,JGR,IL)/FDXM(IBM,JGR)
              ENDDO
            ENDDO
            DO IL=2,NL,2
              TOT1(IBM,IGR,IL)=TOT1(IBM,IGR,IL)*FDXM(IBM,IGR)
              DO JGR=1,NG
               SCAT1(IBM,IGR,JGR,IL)=SCAT1(IBM,IGR,JGR,IL)*FDXM(IBM,IGR)
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        IF(ICODE(2).NE.0) THEN
          BETA(:)=0.0
          IF(ICODE(2).NE.0) THEN
            IBM=IMIX(LX1)
            DO IGR=1,NG
              IF(IBM.EQ.0) CYCLE
              AFACTOR(IGR)=AFACTOR(IGR)/FDXM(IBM,IGR)
              BETA(IGR)=(1.0-2.0*AFACTOR(IGR))/(1.0+2.0*AFACTOR(IGR))
            ENDDO
          ENDIF
          IF(IMPX.GT.0) THEN
            WRITE(6,'(/29H BRERT: SPH CORRECTED ALBEDOS)')
            WRITE(6,20) 'BETA',BETA(:NG)
          ENDIF
        ENDIF
      ENDIF
      IF(IMPX.GT.0) THEN
        WRITE(6,'(/30H BRERT: DIFFUSION COEFFICIENTS)')
        DO IBM=1,NMIX1
          WRITE(6,'(/9H MIXTURE=,I5)') IBM
          WRITE(6,20) 'DIFF',DC1(IBM,:NG)
        ENDDO
      ENDIF
*----
*  SAVE THE OUTPUT NODAL MACROLIB
*----
      ALLOCATE(IJJ(NMIX1),NJJ(NMIX1),IPOS(NMIX1),WORK(NMIX1*NG))
      ISTATE(:)=0
      ISTATE(1)=NG
      ISTATE(2)=NMIX1
      ISTATE(3)=NL
      IF(J_FUEL.GT.0) ISTATE(4)=1
      IF(ICODE(2).NE.0) ISTATE(8)=1  ! physical albedo information
      ISTATE(9)=1  ! diffusion coefficient information
      IF(IDIFF.EQ.1) ISTATE(10)=NL-1 ! NTOT1 is present
      IF(ISPH.EQ.0) ISTATE(12)=3 ! discontinuity factor information
      IF(ISPH.EQ.1) ISTATE(14)=1 ! SPH factor information
      CALL LCMPUT(IPMAC1,'STATE-VECTOR',NSTATE,1,ISTATE)
      CALL LCMPUT(IPMAC1,'ENERGY',NG+1,2,ENER)
      CALL LCMPUT(IPMAC1,'VOLUME',NMIX1,2,VOL1)
      CALL LCMPUT(IPMAC1,'B2  B1HOM',1,2,B2)
      IF(ICODE(2).NE.0) CALL LCMPUT(IPMAC1,'ALBEDO',NG,2,BETA)
      IF(ISPH.EQ.0) THEN
        CALL LCMSIX(IPMAC1,'ADF',1)
          NTYPE=1
          HADF='FD_B'
          CALL LCMPUT(IPMAC1,'NTYPE',1,1,NTYPE)
          CALL LCMPTC(IPMAC1,'HADF',8,HADF)
          CALL LCMPUT(IPMAC1,HADF,NMIX1*NG,2,FDXM)
        CALL LCMSIX(IPMAC1,' ',2)
      ELSE IF(ISPH.EQ.1) THEN
        CALL LCMSIX(IPMAC1,'SPH',1)
          ISTATE(:)=0
          ISTATE(1)=4
          ISTATE(2)=1
          ISTATE(6)=1
          ISTATE(7)=1
          ISTATE(8)=NG
          CALL LCMPUT(IPMAC1,'STATE-VECTOR',NSTATE,1,ISTATE)
        CALL LCMSIX(IPMAC1,' ',2)
      ENDIF
      JPMAC1=LCMLID(IPMAC1,'GROUP',NG)
      DO IGR=1,NG
        KPMAC1=LCMDIL(JPMAC1,IGR)
        DO IBM=1,NMIX1
          WORK(IBM)=VOL1(IBM)*FLX1(IBM,IGR)
        ENDDO
        CALL LCMPUT(KPMAC1,'FLUX-INTG',NMIX1,2,WORK)
        DO IL=1,NL
          WRITE(TEXT12,'(4HNTOT,I1)') IL-1
          CALL LCMPUT(KPMAC1,TEXT12,NMIX1,2,TOT1(:,IGR,IL))
        ENDDO
        CALL LCMPUT(KPMAC1,'DIFF',NMIX1,2,DC1(:,IGR))
        CALL LCMPUT(KPMAC1,'CHI',NMIX1,2,CHI1(:,IGR))
        CALL LCMPUT(KPMAC1,'NUSIGF',NMIX1,2,SIGF1(:,IGR))
        IF(ISPH.EQ.1) THEN
          DO IBM=1,NMIX1
            WORK(IBM)=1.0/FDXM(IBM,IGR)
          ENDDO
          CALL LCMPUT(KPMAC1,'NSPH',NMIX1,2,WORK)
        ENDIF
        DO IL=1,NL
          WRITE(CM,'(I2.2)') IL-1
          WORK(:NMIX1)=SCAT1(:NMIX1,IGR,IGR,IL)
          CALL LCMPUT(KPMAC1,'SIGW'//CM,NMIX1,2,WORK)
          IPOSDE=0
          DO IBM=1,NMIX1
            J2=IGR
            J1=IGR
            DO JGR=1,NG
              IF(SCAT1(IBM,IGR,JGR,IL).NE.0.0) THEN
                J2=MAX(J2,JGR)
                J1=MIN(J1,JGR)
              ENDIF
            ENDDO
            NJJ(IBM)=J2-J1+1
            IJJ(IBM)=J2
            IPOS(IBM)=IPOSDE+1
            DO JGR=J2,J1,-1
             IPOSDE=IPOSDE+1
             IF(IPOSDE.GT.NG*NMIX1) CALL XABORT('BRERT: SCAT OVERFLOW.')
             WORK(IPOSDE)=SCAT1(IBM,IGR,JGR,IL)
            ENDDO
          ENDDO
          CALL LCMPUT(KPMAC1,'SCAT'//CM,IPOSDE,2,WORK)
          CALL LCMPUT(KPMAC1,'NJJS'//CM,NMIX1,1,NJJ)
          CALL LCMPUT(KPMAC1,'IJJS'//CM,NMIX1,1,IJJ)
          CALL LCMPUT(KPMAC1,'IPOS'//CM,NMIX1,1,IPOS)
        ENDDO
      ENDDO
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(WORK,IPOS,NJJ,IJJ)
      DEALLOCATE(VOLTOT,BETA,AFACTOR,FDXP,FDXM,FHOMP,FHOMM)
      RETURN
   20 FORMAT(1X,A9,1P,10E12.4,/(10X,10E12.4))
      END
