*DECK FLUDBV
      SUBROUTINE FLUDBV(CDOOR,IPHASE,JPSYS,JPSTR,NPSYS,IPTRK,IFTRAK,
     1 IPRT,NREG,NUNKNO,NFUNL,NGRP,NMAT,NANIS,LEXAC,MATCOD,VOL,KEYFLX,
     2 TITLE,ILEAK,LEAKSW,XSTRC,XSDIA,B2,NMERG,IMERG,DIFHET,GAMMA,
     3 FLUOLD,SUNKNO,FUNKNO,IPMACR,IPSOU,REBFLG,FLUXC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Find a leakage parameter to match the input DB2 value and find the
* corresponding flux. Vectorial version.
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
*Parameters: input
* CDOOR   name of the geometry/solution operator.
* IPHASE  type of flux solution door (1 for asm 2 for pij).
* JPSYS   pointer to the system LCM list object.
* JPSTR   pointer to the system LCM list object containing isotropic
*         streaming information (=0 if not required).
* NPSYS   non-converged energy group indices.
* IPTRK   pointer to the tracking LCM object.
* IFTRAK  tracking file unit number.
* IPRT    print flag.
* NREG    number of regions.
* NFUNL   second dimension of matrix KEYFLX.
* NGRP    number of energy groups.
* NUNKNO  number of flux/sources unknowns per energy group.
* NMAT    number of mixtures in the internal library.
* NANIS   maximum cross section Legendre order.
* LEXAC   type of exponential function calculation (=.false. to compute
*         exponential functions using tables).
* MATCOD  mixture indices.
* VOL     volumes.
* KEYFLX  index of L-th order flux components in unknown vector.
* TITLE   title.
* ILEAK   method used to include db2 effect:
*         =1 the scattering modified cp matrix is multiplied by PNLR;
*         =2 the reduced cp matrix is multiplied by PNL;
*         =3 sigs0-db2 approximation;
*         =4 albedo approximation;
*         =5 Todorova-type isotropic streaming model;
*         =6 Ecco-type isotropic streaming model;
*         >6 Tibere type anisotropic streaming model.
* LEAKSW  leakage flag (=.true. if leakage is present on the outer
*         surface).
* XSTRC   transport-corrected macroscopic total cross sections.
* XSDIA   transport-corrected macroscopic within-group scattering cross
*         sections.
* B2      buckling.
* NMERG   number of leakage zones.
* IMERG   leakage zone index in each material mixture zone.
* DIFHET  heterogeneous leakage coefficients.
* GAMMA   gamma factors.
* IPMACR  pointer to the macrolib LCM object.
* IPSOU   pointer to the fixed source LCM object.
* REBFLG  ACA or SCR rebalancing flag.
* FLUOLD  flux of the previous outer iteration.
* SUNKNO  input sources.
*
*Parameters: input/output
* FUNKNO  neutron flux.
* SUNKNO  sources with additional db2 contributions.
* FLUXC   flux at the cutoff energy.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER CDOOR*12,TITLE*72
      LOGICAL LEXAC,REBFLG
      TYPE(C_PTR) JPSYS,JPSTR,IPTRK,IPMACR,IPSOU
      INTEGER IPHASE,NPSYS(NGRP),IFTRAK,IPRT,NREG,NUNKNO,NFUNL,
     1 NGRP,NMAT,NANIS,MATCOD(NREG),KEYFLX(NREG,NFUNL),ILEAK,NMERG,
     2 IMERG(NMAT)
      REAL VOL(NREG),XSTRC(0:NMAT,NGRP),XSDIA(0:NMAT,0:NANIS,NGRP),
     1 B2(4),DIFHET(NMERG,NGRP),GAMMA(NGRP),FLUXC(NREG)
      REAL, INTENT(IN) :: FLUOLD(NUNKNO,NGRP)
      REAL, INTENT(INOUT) :: SUNKNO(NUNKNO,NGRP),FUNKNO(NUNKNO,NGRP)
      LOGICAL LEAKSW
*----
*  LOCAL VARIABLES
*----
      PARAMETER(IOUT=6)
      TYPE(C_PTR) KPSYS,KPSTR
      CHARACTER   TEXT12*12
      INTEGER INDD(3),INDC(3),INDB(3)
*----
*  ALLOCATABLE ARRAYS
*----
      REAL, ALLOCATABLE, DIMENSION(:) :: SIGT0,SOURCE2,FUNKNO2,F1,F2,PP
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SOURCE
*
      ALLOCATE(SOURCE(NUNKNO,NGRP))
      SOURCE(:NUNKNO,:NGRP)=SUNKNO(:NUNKNO,:NGRP)
*
      IF(ILEAK.EQ.1) THEN
         IF(NMERG.GT.1) CALL XABORT('FLUDBV: NB. LEAKAGE ZONES > 1.(1)')
         IF(LEAKSW) CALL XABORT('FLUDBV: PNLR OPTION FORBIDDEN.')
         DO 30 IGR=1,NGRP
         IF((NPSYS(IGR).NE.0).AND.(B2(4).NE.0.0)) THEN
            ZNUM=0.0
            ZDEN=0.0
            DO 10 IR=1,NREG
            IBM=MATCOD(IR)
            IND=KEYFLX(IR,1)
            SSS=XSTRC(IBM,IGR)-XSDIA(IBM,0,IGR)
            ZNUM=ZNUM+SSS*FLUOLD(IND,IGR)*VOL(IR)
            ZDEN=ZDEN+FLUOLD(IND,IGR)*VOL(IR)
   10       CONTINUE
            ALP1=ZNUM/(ZNUM+DIFHET(1,IGR)*B2(4)*ZDEN)
            DO 20 IR=1,NREG
            IND=KEYFLX(IR,1)
            SOURCE(IND,IGR)=ALP1*SOURCE(IND,IGR)
   20       CONTINUE
         ENDIF
   30    CONTINUE
         IDIR=0
         CALL DOORFV(CDOOR,JPSYS,NPSYS,IPTRK,IFTRAK,IPRT,NGRP,NMAT,
     1        IDIR,NREG,NUNKNO,IPHASE,LEXAC,MATCOD,VOL,KEYFLX,TITLE,
     2        SOURCE(1,1),FUNKNO(1,1),IPMACR,IPSOU,REBFLG,FLUXC)
      ELSE IF(ILEAK.EQ.2) THEN
         IF(NMERG.GT.1) CALL XABORT('FLUDBV: NB. LEAKAGE ZONES > 1.(2)')
         IF(LEAKSW) CALL XABORT('FLUDBV: PNL OPTION FORBIDDEN.')
         DO 50 IGR=1,NGRP
         IF((NPSYS(IGR).NE.0).AND.(B2(4).NE.0.0)) THEN
            ZNUM=0.0
            ZDEN=0.0
            DO 40 IR=1,NREG
            IBM=MATCOD(IR)
            IND=KEYFLX(IR,1)
            ZNUM=ZNUM+XSTRC(IBM,IGR)*FLUOLD(IND,IGR)*VOL(IR)
            ZDEN=ZDEN+FLUOLD(IND,IGR)*VOL(IR)
   40       CONTINUE
            ALP1=ZNUM/(ZNUM+DIFHET(1,IGR)*B2(4)*ZDEN)
            DO 45 IR=1,NREG
            IND=KEYFLX(IR,1)
            SOURCE(IND,IGR)=ALP1*SOURCE(IND,IGR)-(1.0-ALP1)
     >                 *XSDIA(MATCOD(IR),0,IGR)*FLUOLD(IND,IGR)
   45       CONTINUE
         ENDIF
   50    CONTINUE
         IDIR=0
         CALL DOORFV(CDOOR,JPSYS,NPSYS,IPTRK,IFTRAK,IPRT,NGRP,NMAT,
     1        IDIR,NREG,NUNKNO,IPHASE,LEXAC,MATCOD,VOL,KEYFLX,TITLE,
     2        SOURCE(1,1),FUNKNO(1,1),IPMACR,IPSOU,REBFLG,FLUXC)
      ELSE IF((ILEAK.EQ.3).OR.(ILEAK.EQ.5)) THEN
         DO 70 IGR=1,NGRP
         IF(NPSYS(IGR).NE.0) THEN
            BB=B2(4)
            DO 60 IR=1,NREG
            IND=KEYFLX(IR,1)
            IF(IND.EQ.0) GO TO 60
            IBM=MATCOD(IR)
            IF(IBM.EQ.0) GO TO 60
            INM=IMERG(IBM)
            IF(INM.EQ.0) GO TO 60
            SOURCE(IND,IGR)=SOURCE(IND,IGR)-DIFHET(INM,IGR)*BB*
     1      FLUOLD(IND,IGR)
   60       CONTINUE
         ENDIF
   70    CONTINUE
         IDIR=0
         CALL DOORFV(CDOOR,JPSYS,NPSYS,IPTRK,IFTRAK,IPRT,NGRP,NMAT,
     1        IDIR,NREG,NUNKNO,IPHASE,LEXAC,MATCOD,VOL,KEYFLX,TITLE,
     2        SOURCE(1,1),FUNKNO(1,1),IPMACR,IPSOU,REBFLG,FLUXC)
      ELSE IF(ILEAK.EQ.4) THEN
         ALLOCATE(F1(NREG),F2(NREG))
         DO 80 IGR=1,NGRP
         IF(NPSYS(IGR).NE.0) THEN
           KPSYS=LCMGIL(JPSYS,IGR)
           CALL FLUALB(KPSYS,NREG,NUNKNO,ILCTXS,MATCOD,VOL,KEYFLX,
     >     FLUOLD(1,IGR),SOURCE(1,IGR),XSDIA(0,0,IGR),XSTRC(0,IGR),
     >     F1,F2)
*
           IF(IPRT.GT.2) THEN
             WRITE(IOUT,'(//33H N E U T R O N    S O U R C E S :)')
             WRITE(IOUT,'(1P,6(5X,E15.7))') (SOURCE(KEYFLX(I,1),IGR),
     >       I=1,NREG)
           ENDIF
           FUNKNO(:NUNKNO,IGR)=0.0
           DO 75 IR=1,NREG
           IBM=MATCOD(IR)
           IF(IBM.EQ.0) GO TO 75
           INM=IMERG(IBM)
           IF(INM.EQ.0) GO TO 75
           FUNKNO(KEYFLX(IR,1),IGR)=F1(IR)+DIFHET(INM,IGR)*B2(4)*F2(IR)
   75      CONTINUE
           IF(IPRT.GT.2) THEN
             WRITE(IOUT,'(//33H N E U T R O N    F L U X E S   :)')
             WRITE(IOUT,'(1P,6(5X,E15.7))') (FUNKNO(KEYFLX(I,1),IGR),
     >       I=1,NREG)
           ENDIF
         ENDIF
   80    CONTINUE
         DEALLOCATE(F2,F1)
      ELSE IF(ILEAK.EQ.6) THEN
*        ISOTROPIC STREAMING MODEL (ECCO).
         IF(.NOT.C_ASSOCIATED(JPSTR)) THEN
            CALL XABORT('FLUDBV: MISSING STREAMING INFO(1).')
         ELSE IF(LEAKSW) THEN
            CALL XABORT('FLUDBV: ECCO OPTION FORBIDDEN.')
         ENDIF
         DO 95 IGR=1,NGRP
         IF(NPSYS(IGR).NE.0) THEN
            BB=B2(4)
            DO 90 IR=1,NREG
            IND=KEYFLX(IR,1)
            SOURCE(IND,IGR)=SOURCE(IND,IGR)-BB*FLUOLD(NUNKNO/2+IND,IGR)
   90       CONTINUE
         ENDIF
   95    CONTINUE
         IF(IPRT.GE.3) WRITE(IOUT,'(28H FLUDBV: FUNDAMENTAL FLUXES.)')
         IDIR=0
         CALL DOORFV(CDOOR,JPSYS,NPSYS,IPTRK,IFTRAK,IPRT,NGRP,NMAT,
     1        IDIR,NREG,NUNKNO,IPHASE,LEXAC,MATCOD,VOL,KEYFLX,TITLE,
     2        SOURCE(1,1),FUNKNO(1,1),IPMACR,IPSOU,REBFLG,FLUXC)
         DO 130 IGR=1,NGRP
         IF(NPSYS(IGR).NE.0) THEN
            KPSTR=LCMGIL(JPSTR,IGR)
            CALL LCMLEN(KPSTR,'DRAGON-TXSC',ILCTXS,ITYLCM)
            ALLOCATE(SIGT0(0:ILCTXS-1))
            CALL LCMGET(KPSTR,'DRAGON-TXSC',SIGT0(0))
            ZNUM=0.0
            ZDEN=0.0
            DO 100 IR=1,NREG
               IBM=MATCOD(IR)
               IND=KEYFLX(IR,1)
               ZNUM=ZNUM+SIGT0(IBM)*FUNKNO(IND,IGR)*VOL(IR)
               ZDEN=ZDEN+FUNKNO(IND,IGR)*VOL(IR)
  100       CONTINUE
            DO 110 IR=1,NREG
            IBM=MATCOD(IR)
            IND=KEYFLX(IR,1)
            SOURCE(NUNKNO/2+IND,IGR)=SOURCE(NUNKNO/2+IND,IGR)+
     1      (1.0-GAMMA(IGR))*(ZNUM/ZDEN-SIGT0(IBM))*
     2      FUNKNO(NUNKNO/2+IND,IGR)
  110       CONTINUE
            DEALLOCATE(SIGT0)
            DO 120 IR=1,NREG
            IND=KEYFLX(IR,1)
            SOURCE(NUNKNO/2+IND,IGR)=(SOURCE(NUNKNO/2+IND,IGR)
     1      +FUNKNO(IND,IGR)/3.0)/GAMMA(IGR)
  120       CONTINUE
         ENDIF
  130    CONTINUE
         IF(IPRT.GE.3) WRITE(IOUT,'(30H FLUDBV: FUNDAMENTAL CURRENTS.)')
         ALLOCATE(SOURCE2((NUNKNO/2)*NGRP),FUNKNO2((NUNKNO/2)*NGRP))
         IOF=0
         DO 145 IGR=1,NGRP
         DO 140 IND=1,NUNKNO/2
         IOF=IOF+1
         SOURCE2(IOF)=SOURCE(NUNKNO/2+IND,IGR)
         FUNKNO2(IOF)=FUNKNO(NUNKNO/2+IND,IGR)
  140    CONTINUE
  145    CONTINUE
         IDIR=0
         CALL DOORFV(CDOOR,JPSTR,NPSYS,IPTRK,IFTRAK,IPRT,NGRP,NMAT,
     1        IDIR,NREG,NUNKNO/2,IPHASE,LEXAC,MATCOD,VOL,KEYFLX,TITLE,
     2        SOURCE2(1),FUNKNO2(1),IPMACR,IPSOU,REBFLG,FLUXC)
         IOF=0
         DO 155 IGR=1,NGRP
         DO 150 IND=1,NUNKNO/2
         IOF=IOF+1
         SOURCE(NUNKNO/2+IND,IGR)=SOURCE2(IOF)
         FUNKNO(NUNKNO/2+IND,IGR)=FUNKNO2(IOF)
  150    CONTINUE
  155    CONTINUE
         DEALLOCATE(FUNKNO2,SOURCE2)
      ELSE IF((MOD(ILEAK,10).EQ.7).AND.(IPHASE.EQ.1)) THEN
*        ----
*        TIBERE ANISOTROPIC STREAMING MODEL FOR MOC.
*        ----
         IF(.NOT.C_ASSOCIATED(JPSTR)) THEN
            CALL XABORT('FLUDBV: MISSING STREAMING INFO(2).')
         ELSE IF(LEAKSW) THEN
            CALL XABORT('FLUDBV: TIBERE OPTION FORBIDDEN.')
         ENDIF
* ADD SOURCES FOR FLUX EQUATION
         DO IGR=1,NGRP
           IF(NPSYS(IGR).NE.0) THEN 
             IF((B2(1).NE.0.0).AND.(B2(2).NE.0.0).AND.
     1       (B2(3).NE.0.0)) THEN
               S=0.0
               DO IR=1,NREG 
                 IND=KEYFLX(IR,1)
                 INDC(1)=3*NUNKNO/8+IND
                 INDC(2)=5*NUNKNO/8+IND
                 INDC(3)=7*NUNKNO/8+IND
                 SOURCE(IND,IGR)=SOURCE(IND,IGR)-(B2(1)
     1           *FLUOLD(INDC(1),IGR)+B2(2)*FLUOLD(INDC(2),IGR)+
     2           B2(3)*FLUOLD(INDC(3),IGR))
               ENDDO 
             ENDIF  
           ENDIF 
         ENDDO
         IF(IPRT.GE.3) WRITE(IOUT,'(28H FLUDBV: FUNDAMENTAL FLUXES.)')
           IDIR=0
           NUNKNO4=NUNKNO/4
           ALLOCATE(SOURCE2(NUNKNO4*NGRP),FUNKNO2(NUNKNO4*NGRP))
           IOF=0
           DO IGR=1,NGRP
             DO IND=1,NUNKNO4
               IOF=IOF+1
               SOURCE2(IOF)=SOURCE(IND,IGR)
               FUNKNO2(IOF)=FUNKNO(IND,IGR)
             ENDDO
           ENDDO
           CALL DOORFV(CDOOR,JPSYS,NPSYS,IPTRK,IFTRAK,IPRT,NGRP,NMAT,
     1        IDIR,NREG,NUNKNO4,IPHASE,LEXAC,MATCOD,VOL,KEYFLX,TITLE,
     2        SOURCE2(1),FUNKNO2(1),IPMACR,IPSOU,REBFLG,FLUXC)
           IOF=0
           DO IGR=1,NGRP
             DO IND=1,NUNKNO4
               IOF=IOF+1
               SOURCE(IND,IGR)=SOURCE2(IOF)
               FUNKNO(IND,IGR)=FUNKNO2(IOF) 
             ENDDO
           ENDDO
           DEALLOCATE(FUNKNO2,SOURCE2)   
* ADD SOURCES FOR CURRENT EQUATIONS
           DO IGR=1,NGRP
             IF(NPSYS(IGR).NE.0) THEN 
             KPSTR=LCMGIL(JPSTR,IGR)
             CALL LCMLEN(KPSTR,'DRAGON-TXSC',ILCTXS,ITYLCM)
             ALLOCATE(SIGT0(0:ILCTXS-1))
             CALL LCMGET(KPSTR,'DRAGON-TXSC',SIGT0(0))
             ZNUM=0.0
             ZDEN=0.0
             DO IR=1,NREG
               IBM=MATCOD(IR)
               IND=KEYFLX(IR,1)
               ZNUM=ZNUM+SIGT0(IBM)*FUNKNO(IND,IGR)*VOL(IR)
               ZDEN=ZDEN+FUNKNO(IND,IGR)*VOL(IR)
             ENDDO
             DO IR=1,NREG
               IBM=MATCOD(IR)
               INDD(1)=NUNKNO/4+KEYFLX(IR,1)
               INDD(2)=NUNKNO/2+KEYFLX(IR,1)
               INDD(3)=3*NUNKNO/4+KEYFLX(IR,1)
               DO IDIR=1,3
                 SOURCE(INDD(IDIR),IGR)=SOURCE(INDD(IDIR),IGR)+
     1           (1.0-GAMMA(IGR))*(ZNUM/ZDEN-SIGT0(IBM))*
     2            FUNKNO(INDD(IDIR),IGR)
               ENDDO
             ENDDO
             DO IR=1,NREG
               INDD(1)=NUNKNO/4+KEYFLX(IR,1)
               INDD(2)=NUNKNO/2+KEYFLX(IR,1)
               INDD(3)=3*NUNKNO/4+KEYFLX(IR,1)
               DO IDIR=1,3
                 SOURCE(INDD(IDIR),IGR)=(SOURCE(INDD(IDIR),IGR)
     1           +FUNKNO(KEYFLX(IR,1),IGR)/3.0)/GAMMA(IGR)
               ENDDO
             ENDDO
             DEALLOCATE(SIGT0)
           ENDIF
           ENDDO
         DO IDIR=1,3
           IF(IPRT.GE.3) 
     >     WRITE(IOUT,'(30H FLUDBV: FUNDAMENTAL CURRENTS.)') 
           IF(IDIR.EQ.1) WRITE(6,*)'FUNDAMENTAL CURRENT X '
           IF(IDIR.EQ.2) WRITE(6,*)'FUNDAMENTAL CURRENT Y '
           IF(IDIR.EQ.3) WRITE(6,*)'FUNDAMENTAL CURRENT Z '
           NUNKNO4=NUNKNO/4
           ALLOCATE(SOURCE2(NUNKNO4*NGRP),FUNKNO2(NUNKNO4*NGRP))
           IOF=0
           DO IGR=1,NGRP
             DO IND=1,NUNKNO4
               INDB(1)=NUNKNO/4+IND
               INDB(2)=NUNKNO/2+IND
               INDB(3)=3*NUNKNO/4+IND
               IOF=IOF+1
               SOURCE2(IOF)=SOURCE(INDB(IDIR),IGR)
               FUNKNO2(IOF)=FUNKNO(INDB(IDIR),IGR)
             ENDDO
           ENDDO
           CALL DOORFV(CDOOR,JPSTR,NPSYS,IPTRK,IFTRAK,IPRT,NGRP,NMAT,
     1       IDIR,NREG,NUNKNO4,IPHASE,LEXAC,MATCOD,VOL,KEYFLX,TITLE,
     2       SOURCE2(1),FUNKNO2(1),IPMACR,IPSOU,REBFLG,FLUXC)
           IOF=0
           DO IGR=1,NGRP
             DO IND=1,NUNKNO4
               INDB(1)=NUNKNO/4+IND
               INDB(2)=NUNKNO/2+IND
               INDB(3)=3*NUNKNO/4+IND
               IOF=IOF+1
               SOURCE(INDB(IDIR),IGR)=SOURCE2(IOF)
               FUNKNO(INDB(IDIR),IGR)=FUNKNO2(IOF)   
             ENDDO
           ENDDO
           DEALLOCATE(FUNKNO2,SOURCE2)  
         ENDDO
      ELSE IF((MOD(ILEAK,10).EQ.7).AND.(IPHASE.EQ.2)) THEN
*        ----
*        TIBERE ANISOTROPIC STREAMING MODEL FOR PIJ.
*        ----
         INDD(1)=NUNKNO/4
         INDD(2)=NUNKNO/2
         INDD(3)=3*NUNKNO/4
         NUN4=NUNKNO/4
         ALLOCATE(PP(NREG*NREG))
         DO 210 IGR=1,NGRP
         IF(NPSYS(IGR).NE.0) THEN
           KPSYS=LCMGIL(JPSYS,IGR)
           DO 200 IDIR=1,3
           IF(B2(IDIR).NE.0.0) THEN
             WRITE(TEXT12,'(6HDRAGON,I1,5HP*SCT)') IDIR
             CALL LCMGET(KPSYS,TEXT12,PP)
             DO 190 IR=1,NREG
             IND=KEYFLX(IR,1)
             S=0.0
             DO 180 JREG=1,NREG
             JND=KEYFLX(JREG,1)
             S=S+FLUOLD(INDD(IDIR)+JND,IGR)*PP((JREG-1)*NREG+IR)
 180         CONTINUE
             SOURCE(IND,IGR)=SOURCE(IND,IGR)-B2(IDIR)*S
 190         CONTINUE
           ENDIF
 200     CONTINUE
         ENDIF
 210     CONTINUE
         DEALLOCATE(PP)
         IDIR=0
         CALL DOORFV(CDOOR,JPSYS,NPSYS,IPTRK,IFTRAK,IPRT,NGRP,NMAT,
     1        IDIR,NREG,NUNKNO,IPHASE,LEXAC,MATCOD,VOL,KEYFLX,TITLE,
     2        SOURCE(1,1),FUNKNO(1,1),IPMACR,IPSOU,REBFLG,FLUXC)
         DO 260 IDIR=1,3
         DO 250 IGR=1,NGRP
         IF(NPSYS(IGR).NE.0) THEN
           ZNUM=0.0
           ZDEN=0.0
           DO 220 IR=1,NREG
             IBM=MATCOD(IR)
             IND=KEYFLX(IR,1)
             ZNUM=ZNUM+XSTRC(IBM,IGR)*FLUOLD(IND,IGR)*VOL(IR)
             ZDEN=ZDEN+FLUOLD(IND,IGR)*VOL(IR)
 220       CONTINUE
           DO 230 IR=1,NREG
           IBM=MATCOD(IR)
           IND=KEYFLX(IR,1)
           IND2=INDD(IDIR)+IND
           SOURCE(IND2,IGR)=SOURCE(IND2,IGR)+(1.0-GAMMA(IGR))*
     1     (ZNUM/ZDEN-XSTRC(IBM,IGR))*FLUOLD(IND2,IGR)
 230       CONTINUE
           DO 240 IND=1,NUN4
           IND2=INDD(IDIR)+IND
           SOURCE(IND2,IGR)=(SOURCE(IND2,IGR)+FUNKNO(IND,IGR)/3.0)/
     1     GAMMA(IGR)
 240       CONTINUE
         ENDIF
 250     CONTINUE
         CALL DOORFV(CDOOR,JPSYS,NPSYS,IPTRK,IFTRAK,IPRT,NGRP,
     1        NMAT,IDIR,NREG,NUNKNO,IPHASE,LEXAC,MATCOD,VOL,KEYFLX,
     2        TITLE,SOURCE(INDD(IDIR)+1,1),FUNKNO(INDD(IDIR)+1,1),
     3        IPMACR,IPSOU,REBFLG,FLUXC)
 260     CONTINUE
      ELSE
         CALL XABORT('FLUDBV: TYPE OF LEAKAGE NOT IMPLEMENTED.')
      ENDIF
*----
*  COMPUTE DB2 PARAMETER CORRESPONDING TO ACTUAL LEAKAGE
*----
      IF((IPRT.GT.10).AND.(.NOT.LEAKSW)) THEN
        NUN=NUNKNO
        IF(ILEAK.EQ.6) NUN=NUNKNO/2
        IF(ILEAK.GE.7) NUN=NUNKNO/4
        DO 280 IGR=1,NGRP
        IF(NPSYS(IGR).EQ.0) GO TO 280
        ZNUM=0.0
        ZDEN=0.0
        DO 265 IR=1,NREG
        IND=KEYFLX(IR,1)
        SSS=XSTRC(MATCOD(IR),IGR)-XSDIA(MATCOD(IR),0,IGR)
        ZNUM=ZNUM+VOL(IR)*(SUNKNO(IND,IGR)-SSS*FUNKNO(IND,IGR))
        ZDEN=ZDEN+VOL(IR)*FUNKNO(IND,IGR)
  265   CONTINUE
        DB2NEW=0.0
        IF(ZDEN.GT.0.0) DB2NEW=ZNUM/ZDEN
        DB2OLD=0.0
        VOLTOT=0.0
        DO 270 IR=1,NREG
        INM=IMERG(IR)
        IF(INM.EQ.0) GO TO 270
        DB2OLD=DB2OLD+DIFHET(INM,IGR)*B2(4)*VOL(IR)
        VOLTOT=VOLTOT+VOL(IR)
  270   CONTINUE
        DB2OLD=DB2OLD/VOLTOT
        WRITE(IOUT,'(15H FLUDBV: GROUP=,I5,24H DB2 LEAKAGE PARAMETER F,
     >  12HROM DIFFON =,1P,E13.4/26X,30HACTUAL DB2 LEAKAGE PARAMETER =,
     >  E13.4)') IGR,DB2OLD,DB2NEW
  280   CONTINUE
      ENDIF
      SUNKNO(:NUNKNO,:NGRP)=SOURCE(:NUNKNO,:NGRP)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(SOURCE)
      RETURN
      END
