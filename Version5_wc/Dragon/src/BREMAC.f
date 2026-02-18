*DECK BREMAC
      SUBROUTINE BREMAC(NC,IPMAC2,NG,NL,LX1,NMIX1,NMIX2,IMIX,IMIX1,
     1 IGAP,ILEAKS,IDF,IPRINT,IH,ZKEFF,B2,VOL1,FLX1,DC1,TOT1,CHI1,
     2 SIGF1,SCAT1,HFACT1,JXM,JXP,FHETXM,FHETXP,ADF1)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover macroscopic cross sections.
*
*Copyright:
* Copyright (C) 2021 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert
*
*Parameters: input
* NC      number of sn macrolibs.
* IPMAC2  pointer to the sn macrolib.
* NG      number of energy groups.
* NL      Legendre order of TOT1 and SCAT1 arrays (=1 for isotropic
*         scattering in LAB).
* LX1     number of nodes in the reflector model.
* NMIX1   number of mixtures in the nodal calculation.
* NMIX2   number of mixtures in the sn calculation after edition.
* IMIX    mix index of each node in output data.
* IMIX1   mix index of each node in sn editions.
* IGAP    mix index of the right gap where the surface flux is
*         recovered.
* ILEAKS  type of leakage calculation (=0: no; =1: isotropic;
*         =2: anisotropic).
* IDF     discontinuity factor flag (=0: not used; =3: recovered).
* IPRINT  print parameter
*
*Parameters: output
* IH      H-FACTOR flag (=0: not used; =1: recovered).
* ZKEFF   effective multiplication factor.
* B2      buckling.
* VOL1    volumes.
* FLX1    averaged fluxes
* DC1     diffusion coefficients.
* TOT1    total cross sections.
* CHI1    fission spectra.
* SIGF1   nu*fission cross sections.
* SCAT1   scattering cross sections.
* HFACT1  H-FACTOR values.
* JXM     left boundary currents.
* JXP     right boundary currents.
* FHETXM  left boundary fluxes.
* FHETXP  right boundary fluxes.
* ADF1    assembly discontinuity factors.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER NC
      TYPE(C_PTR) IPMAC2(NC)
      INTEGER NG,NL,LX1,NMIX1,NMIX2,IMIX(LX1),IMIX1(LX1),IGAP(LX1),
     1 ILEAKS,IDF,IPRINT,IH
      REAL ZKEFF(NC),B2(NC),VOL1(NMIX1,NC),FLX1(NMIX1,NG,NC),
     1 DC1(NMIX1,NG,NC),TOT1(NMIX1,NG,NL,NC),CHI1(NMIX1,NG,NC),
     2 SIGF1(NMIX1,NG,NC),SCAT1(NMIX1,NG,NG,NL,NC),HFACT1(NMIX1,NG,NC),
     3 JXM(NMIX1,NG,NC),JXP(NMIX1,NG,NC),FHETXM(NMIX1,NG,NL,NC),
     4 FHETXP(NMIX1,NG,NL,NC),ADF1(NMIX1,NG,NC)
*----
*  LOCAL VARIABLES
*----
      TYPE(C_PTR) JPMAC2,KPMAC2
      DOUBLE PRECISION DSFIS
      CHARACTER CM*2,HADF*8,TEXT12*12
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IJJ,NJJ,IPOS,IMIX2
      REAL, ALLOCATABLE, DIMENSION(:) :: VOL,WORK,SFIS,SFIS1
      REAL, ALLOCATABLE, DIMENSION(:,:) :: DC,CHI,SIGF,HFACT
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: FLX,TOT
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: SCAT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: DCOU
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(VOL(NMIX2),FLX(NMIX2,NG,NL),TOT(NMIX2,NG,NL),
     1 DC(NMIX2,NG),CHI(NMIX2,NG),SIGF(NMIX2,NG),SCAT(NMIX2,NG,NG,NL),
     2 HFACT(NMIX2,NG),IMIX2(NMIX2),SFIS(NMIX2),SFIS1(NMIX1),
     3 DCOU(NMIX2+1,NG))
      ALLOCATE(IJJ(NMIX2),NJJ(NMIX2),IPOS(NMIX2),WORK(NG*NMIX2))
*----
*  SET IMIX2
*----
      IMIX2(:)=0
      DO J=1,LX1
        IBM=IMIX1(J)
        IBG=IGAP(J)
        IF(IBM.GT.NMIX2) CALL XABORT('BREMAC: NMIX2 OVERFLOW(1).')
        IF(IBG.GT.NMIX2) CALL XABORT('BREMAC: NMIX2 OVERFLOW(2).')
        IF(IBM.GT.0) IMIX2(IBM)=IMIX(J)
        IF((IBG.GT.0).AND.(J.LT.LX1)) IMIX2(IBG)=IMIX(J+1)
      ENDDO
*----
*  LOOP OVER SN MACROLIBS
*----
      DO IC=1,NC
*----
*  RECOVER FLUX, MACROSCOPIC CROSS SECTIONS AND DIFFUSION COEFFICIENTS
*----
        CALL LCMGET(IPMAC2(IC),'VOLUME',VOL)
        CALL LCMGET(IPMAC2(IC),'K-EFFECTIVE',ZKEFF(IC))
        B2(IC)=0.0
        IF(ILEAKS.GT.0) THEN
          CALL LCMLEN(IPMAC2(IC),'B2  B1HOM',ILCMLN,ITYLCM)
          IF(ILCMLN.EQ.1) CALL LCMGET(IPMAC2(IC),'B2  B1HOM',B2(IC))
        ENDIF
        JPMAC2=LCMGID(IPMAC2(IC),'GROUP')
        SCAT(:,:,:,:)=0.0
        IH=0
        DO IGR=1,NG
          KPMAC2=LCMGIL(JPMAC2,IGR)
          CALL LCMGET(KPMAC2,'FLUX-INTG',FLX(1,IGR,1))
          CALL LCMGET(KPMAC2,'NTOT0',TOT(1,IGR,1))
          IF(NL.GE.2) THEN
            CALL LCMLEN(KPMAC2,'FLUX-INTG-P1',ILCMLN,ITYLCM)
            IF(ILCMLN.GT.0) THEN
              CALL LCMGET(KPMAC2,'FLUX-INTG-P1',FLX(1,IGR,2))
            ELSE
              FLX(:NMIX2,IGR,2)=FLX(:NMIX2,IGR,1)
            ENDIF
            CALL LCMLEN(KPMAC2,'NTOT1',ILCMLN,ITYLCM)
            IF(ILCMLN.GT.0) THEN
              CALL LCMGET(KPMAC2,'NTOT1',TOT(1,IGR,2))
            ELSE
              TOT(:NMIX2,IGR,2)=TOT(:NMIX2,IGR,1)
            ENDIF
            DO IL=3,NL
              WRITE(TEXT12,'(11HFLUX-INTG-P,I1)') IL-1
              CALL LCMLEN(KPMAC2,TEXT12,ILCMLN,ITYLCM)
              IF(ILCMLN.GT.0) THEN
                CALL LCMGET(KPMAC2,TEXT12,FLX(1,IGR,IL))
              ELSE
                FLX(:NMIX2,IGR,IL)=FLX(:NMIX2,IGR,IL-2)
              ENDIF
              WRITE(TEXT12,'(4HNTOT,I1)') IL-1
              CALL LCMLEN(KPMAC2,TEXT12,ILCMLN,ITYLCM)
              IF(ILCMLN.GT.0) THEN
                CALL LCMGET(KPMAC2,TEXT12,TOT(1,IGR,IL))
              ELSE
                TOT(:NMIX2,IGR,IL)=TOT(:NMIX2,IGR,IL-2)
              ENDIF
            ENDDO
          ENDIF
          CALL LCMGET(KPMAC2,'DIFF',DC(1,IGR))
          CALL LCMGET(KPMAC2,'CHI',CHI(1,IGR))
          CALL LCMGET(KPMAC2,'NUSIGF',SIGF(1,IGR))
          DO IL=1,NL
            WRITE(CM,'(I2.2)') IL-1
            CALL LCMLEN(KPMAC2,'IJJS'//CM,ILCMLN,ITYLCM)
            IF(ILCMLN.EQ.0) CYCLE
            CALL LCMGET(KPMAC2,'IJJS'//CM,IJJ)
            CALL LCMGET(KPMAC2,'NJJS'//CM,NJJ)
            CALL LCMGET(KPMAC2,'IPOS'//CM,IPOS)
            CALL LCMGET(KPMAC2,'SCAT'//CM,WORK)
            DO IBM=1,NMIX2
              IPOSDE=IPOS(IBM)-1
              DO JGR=IJJ(IBM),IJJ(IBM)-NJJ(IBM)+1,-1
                IPOSDE=IPOSDE+1
                SCAT(IBM,IGR,JGR,IL)=WORK(IPOSDE) ! IGR <-- JGR
              ENDDO
            ENDDO
          ENDDO
          CALL LCMLEN(KPMAC2,'H-FACTOR',ILCMLN,ITYLCM)
          IF(ILCMLN.GT.0) THEN
            IH=1
            CALL LCMGET(KPMAC2,'H-FACTOR',HFACT(1,IGR))
          ENDIF
          DO IBM=1,NMIX2
            FLX(IBM,IGR,:NL)=FLX(IBM,IGR,:NL)/VOL(IBM)
          ENDDO
        ENDDO
*----
*  COMPUTE NET CURRENTS BETWEEN NODES USING A BALANCE RELATION
*----
        DCOU(:NMIX2+1,:NG)=0.0D0
        DO IBM=1,NMIX2
          DSFIS=0.0D0
          DO IGR=1,NG
            DSFIS=DSFIS+SIGF(IBM,IGR)*FLX(IBM,IGR,1)
          ENDDO
          DSFIS=DSFIS/ZKEFF(IC)
          DO IGR=1,NG
            DCOU(IBM+1,IGR)=DCOU(IBM,IGR)+VOL(IBM)*(CHI(IBM,IGR)*DSFIS-
     1      (TOT(IBM,IGR,1)+B2(IC)*DC(IBM,IGR))*FLX(IBM,IGR,1))
            DO JGR=1,NG
              DCOU(IBM+1,IGR)=DCOU(IBM+1,IGR)+VOL(IBM)*
     1        SCAT(IBM,IGR,JGR,1)*FLX(IBM,JGR,1)
            ENDDO
          ENDDO
        ENDDO
*----
*  NORMALIZE THE ODD FLUX UNKNOWNS BETWEEN NODES
*----
        IF(NL.GT.1) THEN
          DO J=1,LX1
            IBM=IMIX1(J)
            IBG=IGAP(J)
            IF((IBM.GT.0).AND.(IBG.GT.0)) THEN
              DO IGR=1,NG
                FACT=REAL(DCOU(IBM+1,IGR))/FLX(IBG,IGR,2)
                DO IL=2,NL,2
                  FLX(IBG,IGR,IL)=FLX(IBG,IGR,IL)*FACT
                ENDDO
              ENDDO
            ENDIF
          ENDDO
        ENDIF
*----
*  HOMOGENIZATION OVER THE GAPS AND NODES
*----
        VOL1(:,IC)=0.0
        SFIS1(:)=0.0
        FLX1(:,:,IC)=0.0
        TOT1(:,:,:,IC)=0.0
        DC1(:,:,IC)=0.0
        SIGF1(:,:,IC)=0.0
        CHI1(:,:,IC)=0.0
        SCAT1(:,:,:,:,IC)=0.0
        HFACT1(:,:,IC)=0.0
        DO IL=1,NL,2
          FHETXM(:,:,IL,IC)=1.0
          FHETXP(:,:,IL,IC)=1.0
        ENDDO
        DO IL=2,NL,2
          FHETXM(:,:,IL,IC)=0.0
          FHETXP(:,:,IL,IC)=0.0
        ENDDO
        JXM(:,:,IC)=0.0
        JXP(:,:,IC)=0.0
        ADF1(:,:,IC)=0.0
        DO J=1,LX1
          IBM=IMIX1(J)
          IBG=IGAP(J)
          IF(IBG.GT.0) THEN
            DO IL=1,NL
              IF(IMIX(J).GT.0) FHETXP(IMIX(J),:,IL,IC)=FLX(IBG,:,IL)
              IF(J.LT.LX1) THEN
                IF(IMIX(J+1).GT.0) THEN
                  FHETXM(IMIX(J+1),:,IL,IC)=FLX(IBG,:,IL)
                ENDIF
              ENDIF
            ENDDO
          ENDIF
          IF(IBM.GT.0) THEN
            IF(IMIX(J).GT.0) THEN
              JXM(IMIX(J),:NG,IC)=REAL(DCOU(MAX(1,IBM-1),:NG))
              JXP(IMIX(J),:NG,IC)=REAL(DCOU(IBM+1,:NG))
            ENDIF
          ENDIF
        ENDDO
        DO IBM2=1,NMIX2
          IBM=IMIX2(IBM2)
          IF(IBM.EQ.0) CYCLE
          VOL1(IBM,IC)=VOL1(IBM,IC)+VOL(IBM2)
          SFIS(IBM2)=0.0
          DO IGR=1,NG
            SFIS(IBM2)=SFIS(IBM2)+VOL(IBM2)*FLX(IBM2,IGR,1)*
     1      SIGF(IBM2,IGR)
          ENDDO
          SFIS1(IBM)=SFIS1(IBM)+SFIS(IBM2)
          DO IGR=1,NG
            FLX1(IBM,IGR,IC)=FLX1(IBM,IGR,IC)+VOL(IBM2)*FLX(IBM2,IGR,1)
            DO IL=1,NL
              TOT1(IBM,IGR,IL,IC)=TOT1(IBM,IGR,IL,IC)+VOL(IBM2)*
     1        FLX(IBM2,IGR,1)*TOT(IBM2,IGR,IL)
            ENDDO
            DC1(IBM,IGR,IC)=DC1(IBM,IGR,IC)+VOL(IBM2)*FLX(IBM2,IGR,1)*
     1      DC(IBM2,IGR)
            SIGF1(IBM,IGR,IC)=SIGF1(IBM,IGR,IC)+VOL(IBM2)*
     1      FLX(IBM2,IGR,1)*SIGF(IBM2,IGR)
            CHI1(IBM,IGR,IC)=CHI1(IBM,IGR,IC)+SFIS(IBM2)*CHI(IBM2,IGR)
            DO IL=1,NL
              DO JGR=1,NG
                SCAT1(IBM,IGR,JGR,IL,IC)=SCAT1(IBM,IGR,JGR,IL,IC)+
     1          VOL(IBM2)*FLX(IBM2,JGR,1)*SCAT(IBM2,IGR,JGR,IL)
              ENDDO
            ENDDO
            IF(IH.EQ.1) THEN
              HFACT1(IBM,IGR,IC)=HFACT1(IBM,IGR,IC)+VOL(IBM2)*
     1        FLX(IBM2,IGR,1)*HFACT(IBM2,IGR)
            ENDIF
          ENDDO
        ENDDO
        DO IBM=1,NMIX1
          DO IGR=1,NG
            DO IL=1,NL
              TOT1(IBM,IGR,IL,IC)=TOT1(IBM,IGR,IL,IC)/FLX1(IBM,IGR,IC)
            ENDDO
            DC1(IBM,IGR,IC)=DC1(IBM,IGR,IC)/FLX1(IBM,IGR,IC)
            IF(SFIS1(IBM).NE.0.0) CHI1(IBM,IGR,IC)=CHI1(IBM,IGR,IC)/
     1      SFIS1(IBM)
            SIGF1(IBM,IGR,IC)=SIGF1(IBM,IGR,IC)/FLX1(IBM,IGR,IC)
            DO JGR=1,NG
              DO IL=1,NL
                SCAT1(IBM,IGR,JGR,IL,IC)=SCAT1(IBM,IGR,JGR,IL,IC)/
     1          FLX1(IBM,JGR,IC)
              ENDDO
            ENDDO
            IF(IH.EQ.1) THEN
              HFACT1(IBM,IGR,IC)=HFACT1(IBM,IGR,IC)/FLX1(IBM,IGR,IC)
            ENDIF
          ENDDO
          DO IGR=1,NG
            FLX1(IBM,IGR,IC)=FLX1(IBM,IGR,IC)/VOL1(IBM,IC)
          ENDDO
        ENDDO
*----
*  RECOVER ADF
*----
        IF(IDF.EQ.3) THEN
          CALL LCMSIX(IPMAC2(IC),'ADF',1)
          IF(IPRINT.GT.5) CALL LCMLIB(IPMAC2(IC))
          CALL LCMLEN(IPMAC2(IC),'HADF',NTYPE,ITYLCM)
          IF(NTYPE/2.NE.1) CALL XABORT('BREMAC: NTYPE=1 EXPECTED.')
          CALL LCMGTC(IPMAC2(IC),'HADF',8,HADF)
          CALL LCMLEN(IPMAC2(IC),HADF,ILCMLN,ITYLCM)
          IF(ILCMLN.NE.NMIX1*NG) CALL XABORT('BREMAC: ADF OVERFLOW.')
          CALL LCMGET(IPMAC2(IC),HADF,ADF1(1,1,IC))
        ENDIF
*----
*  PRINT CROSS SECTIONS
*----
        IF(IPRINT.GT.1) THEN
          WRITE(6,'(/36H BREMAC: CROSS SECTION FOR MACROLIB=,I5)') IC
          WRITE(6,'(31H BREMAC: SCATTERING ANISOTROPY=,I5)') NL-1
          WRITE(6,'(/6H KEFF=,1P E13.5,4H B2=,E13.5/)') ZKEFF(IC),B2(IC)
          WRITE(6,10) 'IMIX',IMIX(:)
          WRITE(6,20) 'VOL1',VOL1(:,IC)
          WRITE(6,20) 'SFIS1',SFIS1(:)
          DO IGR=1,NG
            WRITE(6,'(/29H BREMAC: PROCESS ENERGY GROUP,I5)') IGR
            WRITE(6,20) 'FLX1',FLX1(:,IGR,IC)
            WRITE(6,20) 'TOT1-P0',TOT1(:,IGR,1,IC)
            IF(NL.GE.2) WRITE(6,20) 'TOT1-P1',TOT1(:,IGR,2,IC)
            WRITE(6,20) 'SIGR',TOT1(:,IGR,1,IC)-SCAT1(:,IGR,IGR,1,IC)
            WRITE(6,20) 'DC1',DC1(:,IGR,IC)
            WRITE(6,20) 'CHI1',CHI1(:,IGR,IC)
            WRITE(6,20) 'SIGF1',SIGF1(:,IGR,IC)
            IF(IH.EQ.1) WRITE(6,20) 'H-FACTOR',HFACT1(:,IGR,IC)
            DO JGR=1,NG
              IF(IGR.EQ.JGR) THEN
                WRITE(6,20) 'INSCAT1-P0',SCAT1(:,IGR,IGR,1,IC)
                IF(NL.EQ.2) THEN
                  WRITE(6,20) 'INSCAT1-P1',SCAT1(:,IGR,IGR,2,IC)
                ENDIF
              ELSE
                WRITE(6,20) 'OUTSCAT1-P0',SCAT1(:,JGR,IGR,1,IC)
                IF(NL.GE.2) THEN
                  WRITE(6,20) 'OUTSCAT1-P1',SCAT1(:,JGR,IGR,2,IC)
                ENDIF
              ENDIF
            ENDDO
            WRITE(6,20) 'JXM',JXM(:,IGR,IC)
            WRITE(6,20) 'JXP',JXP(:,IGR,IC)
            DO IL=1,NL
              WRITE(TEXT12,'(8HFHETXM-P,I1)') IL-1
              WRITE(6,20) TEXT12(:9),FHETXM(:,IGR,IL,IC)
              WRITE(TEXT12,'(8HFHETXP-P,I1)') IL-1
              WRITE(6,20) TEXT12(:9),FHETXP(:,IGR,IL,IC)
            ENDDO
            IF(IDF.EQ.3) WRITE(6,20) 'ADF1',ADF1(:,IGR,IC)
          ENDDO
        ENDIF
      ENDDO
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(WORK,IPOS,NJJ,IJJ,DCOU,SFIS1,SFIS,IMIX2,HFACT,SCAT,
     1 SIGF,CHI,DC,TOT,FLX,VOL)
      RETURN
*
   10 FORMAT(1X,A12,10I13/(12X,10I13))
   20 FORMAT(1X,A12,1P,10E13.5/(12X,10E13.5))
      END
