*DECK VAL
      SUBROUTINE VAL(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Interpolate the flux distribution.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Chambon
*
*Parameters: input/output
* NENTRY  number of LCM objects or files used by the operator.
* HENTRY  name of each LCM object or file:
*         HENTRY(1): create type(L_FVIEW);
*         HENTRY(2): read-only type(L_TRACK);
*         HENTRY(3): read-only type(L_FLUX).
*         HENTRY(4): read-only type(L_MACROLIB).
* IENTRY  type of each LCM object or file:
*         =1 LCM memory object; =2 XSM file; =3 sequential binary file;
*         =4 sequential ascii file.
* JENTRY  access of each LCM object or file:
*         =0 the LCM object or file is created;
*         =1 the LCM object or file is open DO modifications;
*         =2 the LCM object or file is open in read-only mode.
* KENTRY  LCM object address or file unit number.
*
*Comments:
* The VAL: calling specifications are:
* IFLU  := VAL: TRKNAM FLUNAM :: (descval) ; 
* where
*   IFLU   : name of the \dds{interpflux} data structure (L\_FVIEW} signature) 
*     where the interpolated flux distribution will be stored.
*   TRKNAM : name of the read-only \dds{tracking} data structure (L\_TRACK 
*     signature) containing the tracking. 
*   FLUNAM : name of the read-only \dds{fluxunk} data structure (L\_FLUX
*     signature) containing a transport solution.
*   descval : structure containing the input data to this module to compute 
*     interpolated flux
* 
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER      NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
      TYPE(C_PTR)  KENTRY(NENTRY)
      CHARACTER    HENTRY(NENTRY)*12
*----
*  LOCAL VARIABLES
*----
      INTEGER NSTATE
      PARAMETER (NSTATE=40)
      CHARACTER TEXT12*12,HSIGN*12,CMODUL*12
      INTEGER INDIC,NITMA
      DOUBLE PRECISION DFLOT,ZNORM,XDRCST,EVJ
      REAL FLOT,SIDE,VNORM,DX,DY,DZ,POWER,DELMXD,DELMYD
      LOGICAL L2D,L3D
      INTEGER IGP(NSTATE),IFL(NSTATE),IFV(NSTATE),IMV(NSTATE),NXD,NYD,
     1 NZD,IELEM,NUN,IMPX,DIM,NG,NLF,NXI,NYI,NZI,NREG,ICHX,IDIM,ITYPE,
     2 L4,MAXKN,MKN,LC,ITYLCM,IREG,IGMAX,NMIX,NBFIS,IBM,IFISS,LENGT,
     3 LL4F,LL4X,LL4Y,ITRIAL,ICORN,LXH,ISPLH,NRING,NBLOS
      INTEGER I,IG,J,K
      REAL E(25),SXYZ(3)
      TYPE(C_PTR) IPFVW,IPTRK,IPFLU,JPFLU,JPFVW,IPMAC,JPMAC,KPMAC
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, DIMENSION(:), ALLOCATABLE :: MAT,KFLX,KN,MATXYZ
      REAL, DIMENSION(:), ALLOCATABLE :: XX,YY,ZZ,MXD,MYD,MZD,MXI,MYI,
     1 MZI,FLXD,XXX,YYY,ZZZ,SGD,VOL
      REAL, DIMENSION(:,:), ALLOCATABLE :: FXYZ
      REAL, DIMENSION(:,:), ALLOCATABLE :: ZUFIS
*----
*  PARAMETER VALIDATION
*----
      IF((NENTRY.NE.3).AND.(NENTRY.NE.4)) THEN
        CALL XABORT('VAL: 3 OR 4 PARAMETERS EXPECTED.')
      ENDIF
      IPMAC=C_NULL_PTR
      IF((IENTRY(1).NE.1).AND.(IENTRY(1).NE.2)) CALL XABORT('FLD: LCM '
     1 //'OBJECT EXPECTED AT LHS.')
      IF(JENTRY(1).NE.0) CALL XABORT('VAL: ENTRY IN CREATE MODE '
     1 //'EXPECTED.')
      IPFVW=KENTRY(1)
      DO I=2,NENTRY
        IF(JENTRY(I).NE.2) CALL XABORT('VAL: LCM OBJECT IN READ-ONLY '
     1 //'MODE EXPECTED AT RHS.')
        CALL LCMGTC(KENTRY(I),'SIGNATURE',12,HSIGN)
        IF(HSIGN.EQ.'L_FLUX') THEN
           IPFLU=KENTRY(I)
        ELSEIF(HSIGN.EQ.'L_TRACK') THEN
           IPTRK=KENTRY(I)
           CALL LCMGTC(IPTRK,'TRACK-TYPE',12,CMODUL)
           CALL LCMPTC(KENTRY(1),'LINK.TRACK',12,HENTRY(I))
        ELSEIF(HSIGN.EQ.'L_MACROLIB') THEN
           IPMAC=KENTRY(I)
           CALL LCMPTC(KENTRY(1),'LINK.MACRO',12,HENTRY(I))
        ELSE
           TEXT12=HENTRY(I)
           CALL XABORT('VAL: SIGNATURE OF '//TEXT12//' IS '//HSIGN//
     1     '. L_FLUX, L_TRACK OR L_MACROLIB EXPECTED.')
        ENDIF
      ENDDO
      HSIGN='L_FVIEW'
      CALL LCMPTC(KENTRY(1),'SIGNATURE',12,HSIGN)
      L2D=.TRUE.
      L3D=.TRUE.
*
      CALL LCMGET(IPFLU,'STATE-VECTOR',IFL)
      NG=IFL(1)
*----
*  RECOVER GENERAL TRACKING INFORMATION
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',IGP)
      NREG=IGP(1)
      NUN=IGP(2)
      ITYPE=IGP(6)
      NLF=0
      ICHX=0
      IDIM=1
      LL4F=0
      LL4X=0
      LL4Y=0
      IGMAX=NG+1
      IF((ITYPE.EQ.5).OR.(ITYPE.EQ.6).OR.(ITYPE.EQ.8)) IDIM=2
      IF((ITYPE.EQ.7).OR.(ITYPE.EQ.9)) IDIM=3
      IF(CMODUL.EQ.'BIVAC') THEN
         L3D=.FALSE.
         IELEM=IGP(8)
         ISPLH=IGP(10)
         NLF=IGP(14)
         NXD=IGP(12)
         NYD=IGP(13)
         NZD=1
         IF(NYD.EQ.0) L2D=.FALSE.
         IF(IELEM.LT.0) ICHX=1
         IF(IELEM.GT.0) ICHX=2
         IF((IELEM.GT.0).AND.(IGP(9).EQ.4)) ICHX=3
      ELSE IF(CMODUL.EQ.'TRIVAC') THEN
         L3D=.TRUE.
         IELEM=IGP(9)
         L4=IGP(11)
         ICHX=IGP(12)
         ISPLH=IGP(13)
         NLF=IGP(30)
         NXD=IGP(14)
         NYD=IGP(15)
         NZD=IGP(16)
         LL4F=IGP(25)
         LL4X=IGP(27)
         LL4Y=IGP(28)
         IGMAX=IGP(39)
         IF(NYD.EQ.0) L2D=.FALSE.
         IF(NZD.EQ.0) L3D=.FALSE.
         NZD=MAX(1,NZD)
      ENDIF
      IF((ITYPE.GE.8).AND.(ICHX.NE.2)) THEN
        CALL XABORT('VAL: ONLY RAVIART-THOMAS-SCHNEIDER HEXAGONAL OPTI'
     1  //'ONS ARE SUPPORTED.')
      ENDIF
*----
*  CHECK FOR 'FLUX' OR 'MODE'
*----
      CALL LCMLEN(IPFLU,'FLUX',LENGT,ITYLCM)
      IF(LENGT.EQ.0) THEN
         CALL LCMLEN(IPFLU,'MODE',LENGT,ITYLCM)
         IF(LENGT.GT.0) THEN
            JPFLU=LCMGID(IPFLU,'MODE')
            IPFLU=LCMGIL(JPFLU,1)
         ELSE
            CALL LCMLIB(IPFLU)
            CALL XABORT('VAL: UNABLE TO RECOVER A DIRECT FLUX.')
         ENDIF
      ENDIF
*----
*  READ INPUTS
*----
      IMPX=0
      DX=1.
      DY=1.
      DZ=1.
      ZNORM=1.0D0
      ICORN=1
   10 CALL REDGET(INDIC,NITMA,FLOT,TEXT12,DFLOT)
      IF(INDIC.NE.3) CALL XABORT('VAL: character data expected.')
      IF(TEXT12.EQ.'EDIT') THEN
        CALL REDGET(INDIC,IMPX,FLOT,TEXT12,DFLOT)
        IF(INDIC.NE.1) CALL XABORT('VAL: integer data expected(1).')
      ELSE IF(TEXT12.EQ.'MODE') THEN
        CALL REDGET(INDIC,NITMA,FLOT,TEXT12,DFLOT)
        IF(INDIC.NE.1) CALL XABORT('VAL: integer data expected(2).')
        JPFLU=LCMGID(IPFLU,'MODE')
        IPFLU=LCMGIL(JPFLU,NITMA)
      ELSE IF(TEXT12.EQ.'DIM') THEN
        CALL REDGET(INDIC,DIM,FLOT,TEXT12,DFLOT)
        IF((DIM.LE.0).OR.(DIM.GE.4)) CALL XABORT('VAL: 1<=DIM<=3 expec'
     1   //'ted.')
        CALL REDGET(INDIC,NITMA,DX,TEXT12,DFLOT)
        IF(DIM.GE.2) CALL REDGET(INDIC,NITMA,DY,TEXT12,DFLOT)
        IF(DIM.EQ.3) CALL REDGET(INDIC,NITMA,DZ,TEXT12,DFLOT)
      ELSE IF(TEXT12.EQ.'POWR') THEN
*       NORMALIZATION TO A GIVEN FISSION POWER.
        IF(.NOT.C_ASSOCIATED(IPMAC)) CALL XABORT('VAL: MISSING RHS MAC'
     1  //'ROLIB.')
        CALL LCMGET(IPMAC,'STATE-VECTOR',IMV)
        NMIX=IMV(2)
        NBFIS=IMV(4)
        ALLOCATE(MAT(NREG),KFLX(NREG),VOL(NREG),FLXD(NUN),SGD(NMIX))
        CALL LCMGET(IPTRK,'MATCOD',MAT)
        CALL LCMGET(IPTRK,'KEYFLX',KFLX)
        CALL LCMGET(IPTRK,'VOLUME',VOL)
        CALL REDGET (INDIC,NITMA,POWER,TEXT12,DFLOT) ! power in MW
        IF(INDIC.NE.2) CALL XABORT('VAL: REAL DATA EXPECTED.')
*       NORMALIZATION FACTOR FOR THE DIRECT FLUX.
        EVJ=XDRCST('eV','J')
        ZNORM=0.0D0
        JPFLU=LCMGID(IPFLU,'FLUX')
        JPMAC=LCMGID(IPMAC,'GROUP')
        DO IG=1,NG
          CALL LCMGDL(JPFLU,IG,FLXD)
          KPMAC=LCMGIL(JPMAC,IG)
          CALL LCMLEN(KPMAC,'H-FACTOR',LENGT,ITYLCM)
          IF(LENGT.GT.0) THEN
            CALL LCMGET(KPMAC,'H-FACTOR',SGD)
            SGD(:NMIX)=SGD(:NMIX)*REAL(EVJ*1.0D-6) ! convert eV to MW
          ELSE
            WRITE(6,'(/44H VAL: *** WARNING *** NO H-FACTOR FOUND ON L,
     1      24HCM. USE NU*SIGF INSTEAD.)')
            ALLOCATE(ZUFIS(NMIX,NBFIS))
            CALL LCMGET(KPMAC,'NUSIGF',ZUFIS)
            SGD(:NMIX)=0.0
            DO IFISS=1,NBFIS
              SGD(:NMIX)=SGD(:NMIX)+ZUFIS(:NMIX,IFISS)
            ENDDO
            DEALLOCATE(ZUFIS)
          ENDIF
          DO 20 K=1,NREG
          IBM=MAT(K)
          IF((IBM.EQ.0).OR.(KFLX(K).EQ.0)) GO TO 20
          ZNORM=ZNORM+FLXD(KFLX(K))*VOL(K)*SGD(IBM)
   20     CONTINUE
        ENDDO
        ZNORM=POWER/ZNORM
        WRITE(6,300) ' DIRECT',ZNORM
        DEALLOCATE(SGD,FLXD,VOL,KFLX,MAT)
      ELSE IF(TEXT12.EQ.'NOCCOR') THEN
        ICORN=0
      ELSE IF(TEXT12.EQ.'CCOR') THEN
        ICORN=1
      ELSE IF(TEXT12.EQ.';') THEN
        GO TO 30
      ELSE
        CALL XABORT('VAL: unknownn keyword-->'//TEXT12)
      ENDIF
      GO TO 10
*----
*  Get Data in L_TRACK
*----
   30 ALLOCATE(MAT(NREG),KFLX(NREG))
      CALL LCMGET(IPTRK,'MATCOD',MAT)
      CALL LCMGET(IPTRK,'KEYFLX',KFLX)
      IF(ITYPE.LE.7) THEN
        ALLOCATE(MXD(NXD+1),MYD(NYD+1),MZD(NZD+1))
        ALLOCATE(XX(NREG),YY(NREG),ZZ(NREG))
        CALL LCMGET(IPTRK,'XX',XX)
        IF(L2D) CALL LCMGET(IPTRK,'YY',YY)
      ELSE
*       hexagonal geometry
        ALLOCATE(MZD(NZD+1))
        ALLOCATE(ZZ(NREG))
        CALL LCMGET(IPTRK,'SIDE',SIDE)
      ENDIF
      IF(L3D) CALL LCMGET(IPTRK,'ZZ',ZZ)
*----
*  Compute X and Y mesh from L_TRACK
*----
      ALLOCATE(XXX(NXD),YYY(NYD))
      XXX(:NXD)=0.0
      YYY(:NYD)=0.0
      IREG=0
      IF(L3D.AND.(ITYPE.LE.7)) THEN
*       3D Cartesian
        ALLOCATE(ZZZ(NZD))
        ZZZ(:NZD)=0.0
        DO K=1,NZD
          DO J=1,NYD
            DO I=1,NXD
              IREG=IREG+1
              IF(XX(IREG).NE.0.0) THEN
                IF(XXX(I).EQ.0.0) THEN
                  XXX(I)=XX(IREG)
                ELSE IF(ABS(XXX(I)-XX(IREG)).GT.1.0E-6) THEN
                  CALL XABORT('VAL: inconsistent tracking in X')
                ENDIF
              ENDIF
              IF(YY(IREG).NE.0.0) THEN
                IF(YYY(J).EQ.0.0) THEN
                  YYY(J)=YY(IREG)
                ELSE IF(ABS(YYY(J)-YY(IREG)).GT.1.0E-6) THEN
                  CALL XABORT('VAL: inconsistent tracking in Y')
                ENDIF
              ENDIF
              IF(ZZ(IREG).NE.0.0) THEN
                IF(ZZZ(K).EQ.0.0) THEN
                  ZZZ(K)=ZZ(IREG)
                ELSE IF(ABS(ZZZ(K)-ZZ(IREG)).GT.1.0E-6) THEN
                  CALL XABORT('VAL: inconsistent tracking in Z')
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        IF(IREG.NE.NREG) CALL XABORT('VAL: invalid tracking')
      ELSE IF(L3D) THEN
*       3D hexagonal
        ALLOCATE(ZZZ(NZD))
        ZZZ(:NZD)=0.0
        DO K=1,NZD
          DO I=1,NXD
            IREG=IREG+1
            IF(ZZ(IREG).NE.0.0) THEN
              IF(ZZZ(K).EQ.0.0) THEN
                ZZZ(K)=ZZ(IREG)
              ELSE IF(ABS(ZZZ(K)-ZZ(IREG)).GT.1.0E-6) THEN
                CALL XABORT('VAL: inconsistent tracking in Z')
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        IF(IREG.NE.NREG) CALL XABORT('VAL: invalid tracking')
      ELSE IF(L2D.AND.(ITYPE.LE.7)) THEN
*       2D Cartesian
        DO J=1,NYD
          DO I=1,NXD
            IREG=IREG+1
            IF(XX(IREG).NE.0.0) THEN
              IF(XXX(I).EQ.0.0) THEN
                XXX(I)=XX(IREG)
              ELSE IF(ABS(XXX(I)-XX(IREG)).GT.1.0E-6) THEN
                CALL XABORT('VAL: inconsistent tracking in X')
              ENDIF
            ENDIF
            IF(YY(IREG).NE.0.0) THEN
              IF(YYY(J).EQ.0.0) THEN
                YYY(J)=YY(IREG)
              ELSE IF(ABS(YYY(J)-YY(IREG)).GT.1.0E-6) THEN
                CALL XABORT('VAL: inconsistent tracking in Y')
              ENDIF
            ENDIF
          ENDDO
        ENDDO
        IF(IREG.NE.NREG) CALL XABORT('VAL: invalid tracking')
      ELSE IF(ITYPE.EQ.2) THEN
*       1D Cartesian
        DO I=1,NXD
          IREG=IREG+1
          IF(XX(IREG).NE.0.0) THEN
            IF(XXX(I).EQ.0.0) THEN
              XXX(I)=XX(IREG)
            ELSE IF(ABS(XXX(I)-XX(IREG)).GT.1.0E-6) THEN
              CALL XABORT('VAL: inconsistent tracking in X')
            ENDIF
          ENDIF
        ENDDO
        IF(IREG.NE.NREG) CALL XABORT('VAL: invalid tracking')
      ENDIF
      IF(ITYPE.LE.7) THEN
        MXD(1)=0.0
        MYD(1)=0.0
        DO I=1,NXD
          MXD(I+1)=MXD(I)+XXX(I)
        ENDDO
        IF(L2D) THEN
          MYD(1)=0.0
          DO I=1,NYD
            MYD(I+1)=MYD(I)+YYY(I)
          ENDDO
        ELSE
          MYD(2)=0.0
        ENDIF
        DEALLOCATE(YYY,XXX)
      ENDIF
      MZD(1)=0.0
      IF(L3D) THEN
        DO I=1,NZD
          MZD(I+1)=MZD(I)+ZZZ(I)
        ENDDO
        DEALLOCATE(ZZZ)
      ELSE
        MZD(2)=0.0
      ENDIF
*----
*  Perform interpolation
*----
*     Compute points to interpolate
      IF(ITYPE.LE.7) THEN
        DELMXD=MXD(NXD+1)-MXD(1)
        DELMYD=MYD(NYD+1)-MYD(1)
      ELSE
        NRING=INT((1+SQRT(1.0+(NXD/(3*ISPLH**2)-1)*4.0/3.0))/2.0)
        DELMXD=2.0*(1.0+1.5*(NRING-1))*SIDE*ISPLH
        DELMYD=(1.0+2.0*(NRING-1))*SIDE*SQRT(3.0)*ISPLH
      ENDIF
      NXI=INT(DELMXD/DX)+1
      NYI=INT(DELMYD/DY)+1
      NZI=INT((MZD(NZD+1)-MZD(1))/DZ)+1
      ALLOCATE(MXI(NXI),MYI(NYI),MZI(NZI))
      ALLOCATE(MATXYZ(NXI*NYI*NZI),FXYZ(NXI*NYI*NZI,NG))
      IF(NXI.LE.1) CALL XABORT('VAL: UNABLE TO INTERPOLATE IN 1D.')
      SXYZ(:3)=0.0
      DO I=1,NXI
        IF(ITYPE.LE.7) THEN
          MXI(I)=MXD(1)+DELMXD*REAL(I-1)/REAL(NXI-1)
        ELSE
          MXI(I)=-DELMXD/2.0+DELMXD*REAL(I-1)/REAL(NXI-1)
        ENDIF
      ENDDO
      SXYZ(1)=DELMXD/REAL(NXI-1)
      IF(L2D) THEN
        IF(NYI.LE.1) CALL XABORT('VAL: UNABLE TO INTERPOLATE IN 2D.')
        DO I=1,NYI
          IF(ITYPE.LE.7) THEN
            MYI(I)=MYD(1)+DELMYD*REAL(I-1)/REAL(NYI-1)
          ELSE
            MYI(I)=-DELMYD/2.0+DELMYD*REAL(I-1)/REAL(NYI-1)
          ENDIF
        ENDDO
        SXYZ(2)=DELMYD/REAL(NYI-1)
      ELSE
        SXYZ(2)=1.0
      ENDIF
      IF(L3D) THEN
        IF(NZI.LE.1) CALL XABORT('VAL: UNABLE TO INTERPOLATE IN 3D.')
        DO I=1,NZI
          MZI(I)=MZD(1)+(MZD(NZD+1)-MZD(1))*REAL(I-1)/REAL(NZI-1)
        ENDDO
        SXYZ(3)=(MZD(NZD+1)-MZD(1))/REAL(NZI-1)
      ELSE
        SXYZ(3)=1.0
      ENDIF
      JPFLU=LCMGID(IPFLU,'FLUX')
*     Get Data in L_FLUX
      ALLOCATE(FLXD(NUN))
      IF((ICHX.EQ.4).OR.(ICHX.EQ.5).OR.(ICHX.EQ.6)) THEN
*       recover removal xs and diffusion coefficients in JPMAC
        IF(.NOT.C_ASSOCIATED(IPMAC)) CALL XABORT('VAL: MISSING RHS MAC'
     1  //'ROLIB.')
        CALL LCMGET(IPMAC,'STATE-VECTOR',IMV)
        NMIX=IMV(2)
        JPMAC=LCMGID(IPMAC,'GROUP')
      ENDIF
      VNORM=1.0
      DO IG=1,NG
        CALL LCMGDL(JPFLU,IG,FLXD)
*       Perform normalization
        FLXD(:NUN)=FLXD(:NUN)*REAL(ZNORM)
*       Perform interpolation
        IF(L3D) THEN
          IF(ICHX.EQ.1) THEN
*           Variational collocation method
            CALL LCMLEN(IPTRK,'KN',MAXKN,ITYLCM)
            MKN=MAXKN/(NXD*NYD*NZD)
            ALLOCATE(KN(MAXKN))
            CALL LCMGET(IPTRK,'KN',KN)
            CALL LCMSIX(IPTRK,'BIVCOL',1)
            CALL LCMLEN(IPTRK,'T',LC,ITYLCM)
            CALL LCMGET(IPTRK,'E',E)
            CALL LCMSIX(IPTRK,' ',2)
            CALL VALUE2(LC,MKN,NXD,NYD,NZD,L4,MXI,MYI,MZI,MXD,MYD,MZD,
     1      FLXD,MAT,KN,NXI,NYI,NZI,E,MATXYZ,FXYZ(1,IG))
            DEALLOCATE(KN)
          ELSE IF(ICHX.EQ.2) THEN
*           Raviart-Thomas finite element method
            IF(ITYPE.LE.7) THEN
              CALL VALUE4(IELEM,NUN,NXD,NYD,NZD,MXI,MYI,MZI,MXD,MYD,MZD,
     1        FLXD,MAT,KFLX,NXI,NYI,NZI,MATXYZ,FXYZ(1,IG))
            ELSE
*             Hexagonal geometry (LXH hexagons)
              LXH=NXD/(3*ISPLH**2)
              NBLOS=LXH*NZD*ISPLH**2
              CALL VALUE6(IELEM,NUN,LXH,NZD,NBLOS,ISPLH,MXI,MYI,MZI,
     1        SIDE,MZD,FLXD,MAT,KFLX,NXI,NYI,NZI,SXYZ,VNORM,MATXYZ,
     2        FXYZ(1,IG))
            ENDIF
          ELSE IF(ICHX.EQ.3) THEN
*           Nodal collocation method (MCFD)
            CALL VALUE1(IDIM,NXD,NYD,NZD,L4,MXI,MYI,MZI,MXD,MYD,MZD,
     1      FLXD,MAT,IELEM,NXI,NYI,NZI,MATXYZ,FXYZ(1,IG))
          ELSE IF(ICHX.EQ.6) THEN
*           Analytic nodal method (ANM)
            IF(IMPX.GT.0) WRITE(6,320) ICORN
            CALL LCMLEN(IPTRK,'KN',MAXKN,ITYLCM)
            ALLOCATE(KN(MAXKN))
            CALL LCMGET(IPTRK,'KN',KN)
            KPMAC=LCMGIL(JPMAC,IG)
            CALL VALU5(KPMAC,NXD,NYD,NZD,LL4F,LL4X,LL4Y,NUN,NMIX,MXI,
     1      MYI,MZI,MXD,MYD,MZD,FLXD,MAT,KFLX,KN,NXI,NYI,NZI,ICORN,
     2      MATXYZ,FXYZ(1,IG))
            DEALLOCATE(KN)
          ELSE
            CALL XABORT('VAL: INTERPOLATION NOT IMPLEMENTED(1).')
          ENDIF
        ELSE IF(L2D) THEN
          IF(ICHX.EQ.1) THEN
*           Variational collocation method
            CALL LCMLEN(IPTRK,'KN',MAXKN,ITYLCM)
            MKN=MAXKN/(NXD*NYD)
            ALLOCATE(KN(MAXKN))
            CALL LCMGET(IPTRK,'KN',KN)
            CALL LCMSIX(IPTRK,'BIVCOL',1)
            CALL LCMLEN(IPTRK,'T',LC,ITYLCM)
            CALL LCMGET(IPTRK,'E',E)
            CALL LCMSIX(IPTRK,' ',2)
            CALL VALU2B(LC,MKN,NXD,NYD,L4,MXI,MYI,MXD,MYD,FLXD,MAT,KN,
     1      NXI,NYI,E,MATXYZ,FXYZ(1,IG))
          ELSE IF(ICHX.EQ.2) THEN
*           Raviart-Thomas finite element method
            IF(ITYPE.LE.7) THEN
              CALL VALU4B(IELEM,NUN,NXD,NYD,MXI,MYI,MXD,MYD,FLXD,MAT,
     1        KFLX,NXI,NYI,MATXYZ,FXYZ(1,IG))
            ELSE
*             Hexagonal geometry (LXH hexagons)
              LXH=NXD/(3*ISPLH**2)
              NBLOS=LXH*ISPLH**2
              CALL VALU6B(IELEM,NUN,LXH,NBLOS,ISPLH,MXI,MYI,SIDE,FLXD,
     1        MAT,KFLX,NXI,NYI,SXYZ,VNORM,MATXYZ,FXYZ(1,IG))
            ENDIF
          ELSE IF(ICHX.EQ.3) THEN
*           Nodal collocation method (MCFD)
            CALL VALU1B(IDIM,NXD,NYD,L4,MXI,MYI,MXD,MYD,FLXD,MAT,IELEM,
     1      NXI,NYI,MATXYZ,FXYZ(1,IG))
          ELSE IF(ICHX.EQ.6) THEN
*           Analytic nodal method (ANM)
            IF(IMPX.GT.0) WRITE(6,320) ICORN
            CALL LCMLEN(IPTRK,'KN',MAXKN,ITYLCM)
            ALLOCATE(KN(MAXKN))
            CALL LCMGET(IPTRK,'KN',KN)
            KPMAC=LCMGIL(JPMAC,IG)
            CALL VALU5B(KPMAC,NXD,NYD,LL4F,LL4X,NUN,NMIX,MXI,MYI,MXD,
     1      MYD,FLXD,MAT,KFLX,KN,NXI,NYI,ICORN,MATXYZ,FXYZ(1,IG))
            DEALLOCATE(KN)
          ELSE
            CALL XABORT('VAL: INTERPOLATION NOT IMPLEMENTED(2).')
          ENDIF
        ELSE
          IF(ICHX.EQ.4) THEN
*           Coarse mesh finite differences
            KPMAC=LCMGIL(JPMAC,IG)
            ITRIAL=0
            CALL VALU5C(KPMAC,NXD,L4,NMIX,MXI,MXD,FLXD,MAT,NXI,ITRIAL,
     1      MATXYZ,FXYZ(1,IG))
          ELSE IF((ICHX.EQ.5).OR.(ICHX.EQ.6)) THEN
*           Nodal expansion method (NEM) or analytic nodal method (ANM)
            KPMAC=LCMGIL(JPMAC,IG)
            ITRIAL=1
            IF((ICHX.EQ.5).AND.(IG.GE.IGMAX)) ITRIAL=2
            CALL VALU5C(KPMAC,NXD,NUN,NMIX,MXI,MXD,FLXD,MAT,NXI,ITRIAL,
     1      MATXYZ,FXYZ(1,IG))
          ELSE
            CALL XABORT('VAL: INTERPOLATION NOT IMPLEMENTED(3).')
          ENDIF
        ENDIF
      ENDDO
*----
*  Normalize voxel sides in hexagonal cases
*----
      IF(ITYPE.GE.8) THEN
        IF(IMPX.GT.0) WRITE(6,'(35H VAL: VOXEL SIDE NORMALIZATION FACT,
     1  3HOR=,1P,E12.4)') VNORM
        IF(L3D) THEN
          SXYZ(:3)=SXYZ(:3)*VNORM
        ELSE
          SXYZ(:2)=SXYZ(:2)*VNORM
        ENDIF
      ENDIF
*----
*  Save results
*----
      CALL LCMPUT(IPFVW,'MXI',NXI,2,MXI)
      IF(L2D) CALL LCMPUT(IPFVW,'MYI',NYI,2,MYI)
      IF(L3D) CALL LCMPUT(IPFVW,'MZI',NZI,2,MZI)
      CALL LCMPUT(IPFVW,'SXYZ',3,2,SXYZ)
      IFV(:NSTATE)=0
      IFV(1)=NG
      IFV(2)=NXI
      IFV(3)=NYI
      IFV(4)=NZI
      CALL LCMPUT(IPFVW,'STATE-VECTOR',NSTATE,1,IFV)
      CALL LCMPUT(IPFVW,'MATXYZ',NXI*NYI*NZI,1,MATXYZ)
      JPFVW=LCMLID(IPFVW,'FLUX',NG)
      DO IG=1,NG
        CALL LCMPDL(JPFVW,IG,NXI*NYI*NZI,2,FXYZ(1,IG))
      ENDDO
*----
*  Save results
*----
      IF(IMPX.GE.1)THEN
        WRITE(6,*) 'Mesh along X-direction'
        WRITE(6,310) (MXI(I),I=1,NXI)
        IF(L2D) THEN
          WRITE(6,*) 'Mesh along Y-direction'
          WRITE(6,310) (MYI(I),I=1,NYI)
        ENDIF
        IF(L3D) THEN
          WRITE(6,*) 'Mesh along Z-direction'
          WRITE(6,310) (MZI(I),I=1,NZI)
        ENDIF
        IF(IMPX.GE.2)THEN
          WRITE(6,*) 'Flux distribution:'
          DO IG=1,NG
            WRITE(6,*) 'Group',IG
            DO K=1,NZI
              WRITE(6,*) 'Plane',K
              DO J=1,NYI
                WRITE(6,310) (FXYZ(I+(J-1+(K-1)*NYI)*NXI,IG),I=1,NXI)
              ENDDO
            ENDDO
          ENDDO
        ENDIF
      ENDIF
*----
*  RELEASE GENERAL TRACKING INFORMATION
*----
      DEALLOCATE(FLXD)
      DEALLOCATE(FXYZ,MATXYZ)
      DEALLOCATE(MXI,MYI,MZI)
      IF(ITYPE.LE.7) THEN
        DEALLOCATE(MXD,MYD,MZD)
        DEALLOCATE(XX,YY,ZZ)
      ELSE
        DEALLOCATE(MZD)
        DEALLOCATE(ZZ)
      ENDIF
      DEALLOCATE(KFLX,MAT)
      RETURN
  300 FORMAT(/6H VAL: ,A7,28H FLUX NORMALIZATION FACTOR =,1P,E13.5)
  310 FORMAT(1X,1P,12E12.4)
  320 FORMAT(/43H VAL: CORNER FLUX CORRECTION (0/1: OFF/ON)=,I3)
      END
