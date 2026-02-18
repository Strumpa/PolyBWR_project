*DECK MCCGA
      SUBROUTINE MCCGA(IPSYS,NPSYS,IPTRK,IFTRAK,IMPX,NGRP,NBMIX,NANI,
     1 NALBP,ISTRM)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of PJJ for flux integration when isotropic scattering is
* considered and calculation of preconditioning matrices for 
* Algebraic Collapsing Acceleration or Self-Collision Probability 
* acceleration of inner iterations (vectorial version).
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): A. Hebert and R. Le Tellier
*
*Parameters: input/output
* IPSYS   pointer to the PIJ LCM object (L_PIJ signature). IPSYS is a
*         list of NGRP directories.
* NPSYS   index array pointing to the IPSYS list component corresponding
*         to each energy group. Set to zero if a group is not to be
*         processed. Usually, NPSYS(I)=I.
* IPTRK   pointer to the tracking (L_TRACK signature).
* IFTRAK  tracking file unit number.
* IMPX    print flag (equal to zero for no print).
* NGRP    number of energy groups.
* NBMIX   number of mixtures.
* NANI    number of Legendre orders.
* NALBP   number of physical albedos.
* ISTRM   type of streaming effect:
*         =1 no streaming effect;
*         =2 isotropic streaming effect;
*         =3 anisotropic streaming effect.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPSYS,IPTRK
      INTEGER IFTRAK,IMPX,NGRP,NBMIX,NANI,NALBP,ISTRM,NPSYS(NGRP)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40)
      INTEGER JPAR(NSTATE),TRTY,PACA,STIS,IGB(8)
      CHARACTER*4 TEXT4
      REAL ZREAL(4),DELU,FACSYM
      LOGICAL LEXA,LEXF,CYCLIC,LTMT,LACA,LPJJ,LPJJAN,LVOID,LPRISM,
     1 LBIHET
      TYPE(C_PTR) JPSYS
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NGIND
      REAL, ALLOCATABLE, DIMENSION(:) :: CPO
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SIGAL
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: CAZ0,CAZ1,CAZ2
      TYPE(C_PTR), ALLOCATABLE, DIMENSION(:) :: KPSYS
*----
*  GENERIC INTERFACES
*----
      INTERFACE
        SUBROUTINE SUBPJJ_TEMPLATE(M,NSEG,NSUB,LPS,IS,JS,H,KANGL,NOM,
     1                   NZON,TR,W,NFI,NREG,PJJ,PSJ,IMU,NMU,NFUNL,NANGL,
     2                   NPJJM,TRHAR,LPJJAN,PJJIND)
          INTEGER M,NSEG,NSUB,NFI,NREG,LPS,IS(NFI-NREG+1),JS(LPS),
     1    NZON(NFI),KANGL(NSUB),NOM(NSEG),IMU,NMU,NFUNL,NANGL,NPJJM,
     2    PJJIND(NPJJM,2)
          REAL TR(0:M),PSJ(LPS),TRHAR(NMU,NFUNL,NANGL)
          DOUBLE PRECISION W,H(NSUB),PJJ(NREG,NPJJM)
          LOGICAL LPJJAN
        END SUBROUTINE SUBPJJ_TEMPLATE
        !
        SUBROUTINE SUBDSP_TEMPLATE(N,NFI,NLONG,LC,NZON,NOM,KM,MCU,IM,
     1                   PREV,NEXT,H)
          INTEGER N,NFI,NLONG,LC,NZON(NFI),NOM(N),KM(NLONG),MCU(LC),
     1    IM(NLONG),PREV(N),NEXT(N)
          DOUBLE PRECISION, OPTIONAL :: H(N)
        END SUBROUTINE SUBDSP_TEMPLATE
        !
        SUBROUTINE SUBDSC_TEMPLATE(N,M,NFI,NOM,NZON,H,XST,XSW,DINV,B,A)
          INTEGER N,M,NFI,NOM(N),NZON(NFI)
          REAL XST(0:M),XSW(0:M)
          DOUBLE PRECISION H(N),DINV(N),B(N),A(N)
        END SUBROUTINE SUBDSC_TEMPLATE
        !
        SUBROUTINE SUBDS2_TEMPLATE(SUBDSC,LC,M,N,H,NOM,NZON,TR,SC,W,NFI,
     1                   DIAGF,DIAGQ,CA,CQ,PREV,NEXT,DINV2,A2,B2)
          INTEGER LC,M,N,NFI,NZON(NFI),NOM(N),PREV(N),NEXT(N)
          DOUBLE PRECISION W,H(N),CA(LC),DIAGF(NFI),DINV2(N),A2(N),B2(N)
          REAL TR(0:M),SC(0:M),DIAGQ(NFI),CQ(LC)
          EXTERNAL SUBDSC
        END SUBROUTINE SUBDS2_TEMPLATE
      END INTERFACE
      PROCEDURE(SUBPJJ_TEMPLATE), POINTER :: SUBPJJ
      PROCEDURE(SUBDSP_TEMPLATE), POINTER :: SUBDSP
      PROCEDURE(SUBDSC_TEMPLATE), POINTER :: SUBDSC
      PROCEDURE(SUBDS2_TEMPLATE), POINTER :: SUBDS2
      PROCEDURE(SUBPJJ_TEMPLATE) :: MCGDSCA,MCGDSCE,MCGDDDF
      PROCEDURE(SUBDSP_TEMPLATE) :: MOCDSP,MCGDSP
      PROCEDURE(SUBDSC_TEMPLATE) :: MCGDS2E,MCGDS2A
      PROCEDURE(SUBDS2_TEMPLATE) :: MOCDS2,MCGDS2
*----
*  RECOVER MCCG3D SPECIFIC PARAMETERS
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',JPAR)
      IF(JPAR(4).GT.NBMIX) CALL XABORT('MCCGA: INVALID NBMIX.')
      IF(IFTRAK.LE.0) CALL XABORT('MCCGA: INVALID TRACKING FILE.')
*     recover state-vector information
      LBIHET=JPAR(40).NE.0
      IF(LBIHET) THEN
         CALL LCMSIX(IPTRK,'BIHET',1)
         CALL LCMGET(IPTRK,'PARAM',IGB)
         NREG=IGB(3)
         CALL LCMSIX(IPTRK,' ',2)
      ELSE
         NREG=JPAR(1)
      ENDIF
      NFI=NREG+JPAR(5)
      IF(JPAR(6).NE.NANI) CALL XABORT('MCCGA: INVALID NANI.')
      TRTY=JPAR(9)
      IF(TRTY.EQ.1) THEN
         IF(JPAR(5).EQ.0) NFI=NREG+1
         CYCLIC=.TRUE.
         NLONG=NREG
      ELSE
         CYCLIC=.FALSE.
         NLONG=NFI
      ENDIF
      NBCDA=JPAR(30)
      NZP=JPAR(39)
      LPRISM=(NZP.NE.0)
      CALL LCMGET(IPTRK,'MCCG-STATE',JPAR)
      NMU=JPAR(2)
      NMAX=JPAR(5)
      IAAC=JPAR(7)
      STIS=JPAR(15)
      ISCR=JPAR(8)
      LC=JPAR(6)
      LPS=JPAR(9)
      PACA=JPAR(10)
      LC0=JPAR(17)
      LTMT=(JPAR(14).EQ.1)
      LEXA=(JPAR(11).EQ.1)
      LEXF=(JPAR(12).EQ.1)
      NPJJM=JPAR(16)
*     recover real parameters
      CALL LCMGET(IPTRK,'REAL-PARAM',ZREAL)
      HDD=ZREAL(2)
      DELU=ZREAL(3)
      FACSYM=ZREAL(4)
*     recover tracking file information
      REWIND IFTRAK
      READ(IFTRAK) TEXT4,NCOMNT,NBTR,IFMT
      DO ICOM=1,NCOMNT
         READ(IFTRAK)
      ENDDO
      READ(IFTRAK) NDIM,ISPEC,N2REG,N2SOU,NBALB,NCOR,NANGL,MXSUB,MXSEG
      IF(NCOR.NE.1) 
     1 CALL XABORT('MCCGA: INVALID TRACKING FILE: NCOR.NE.1')
      IF(NBCDA.NE.NBALB) CALL XABORT('MCCGA: INVALID NUMBER OF ALBEDOS'
     1 //' IN THE TRACKING FILE.')
      READ(IFTRAK)
      READ(IFTRAK)
      READ(IFTRAK)
      READ(IFTRAK)
      ALLOCATE(CAZ0(NANGL),CAZ1(NANGL),CAZ2(NANGL),CPO(NMU))
      IF(NDIM.EQ.2) THEN
         CALL LCMGET(IPTRK,'XMU$MCCG',CPO)
         READ(IFTRAK) (CAZ1(JJ),CAZ2(JJ),JJ=1,NANGL)
      ELSE ! NDIM.EQ.3
**        correction Sylvie Musongela, december 2019
         READ(IFTRAK) (CAZ1(JJ),CAZ2(JJ),CAZ0(JJ),JJ=1,NANGL)
         DO JJ=1,NANGL
            CAZ1(JJ)=CAZ1(JJ)/SQRT(1.0D0-CAZ0(JJ)*CAZ0(JJ))
            CAZ2(JJ)=CAZ2(JJ)/SQRT(1.0D0-CAZ0(JJ)*CAZ0(JJ))
         ENDDO
      ENDIF
*---
* DETERMINE THE NUMBER OF GROUPS TO BE PROCESSED
* RECOVER POINTERS TO EACH GROUP PROPERTIES
* CREATE AN INDEX FOR THE GROUPS TO BE PROCESSED
*---
      NGEFF=0
      DO IG=1,NGRP
         IOFSET=NPSYS(IG)
         IF(IOFSET.NE.0) NGEFF=NGEFF+1
      ENDDO
      ALLOCATE(NGIND(NGEFF),KPSYS(NGEFF))
      II=1
      DO IG=1,NGRP
         IOFSET=NPSYS(IG)
         IF(IOFSET.NE.0) THEN
            NGIND(II)=IG
            IF(LBIHET) THEN
               JPSYS=LCMGIL(IPSYS,IOFSET)
               KPSYS(II)=LCMGID(JPSYS,'BIHET')
            ELSE
               KPSYS(II)=LCMGIL(IPSYS,IOFSET)
            ENDIF
            II=II+1
         ENDIF
      ENDDO
*----
*  CONSTRUCT TOTAL CROSS SECTIONS ARRAY AND CHECK FOR ZERO CROSS SECTION
*----
      ALLOCATE(SIGAL(-NBCDA:NBMIX,NGEFF))
      CALL MCGSIG(IPTRK,NBMIX,NGEFF,NBCDA,NALBP,KPSYS,SIGAL,LVOID)
      IF((LVOID).AND.(STIS.EQ.-1)) THEN
         IF(IMPX.GT.0) 
     1       WRITE(6,*) 'VOID EXISTS -> STIS SET TO 1 INSTEAD OF -1'
         STIS=1
      ENDIF
*---
* IS THERE SOMETHING TO DO ?
*---
      LACA=(IAAC.GT.0)
      LPJJ=((STIS.EQ.1).OR.(ISCR.GT.0))
      IF(.NOT.(LACA.OR.LPJJ)) GOTO 10
      LPJJAN=(LPJJ.AND.(NANI.GT.1))
      IF(HDD.GT.0.0) THEN
         ISCH=0
      ELSEIF(LEXF) THEN
         ISCH=-1
      ELSE
         ISCH=1
      ENDIF
*----
*  PRECONDITIONING MATRICES CALCULATION
*----
      IF(ISCH.EQ.1) THEN
*     PJJ/SCR: Step-Characteristics Scheme with Tabulated Exponentials
         IF(CYCLIC) THEN
*        ACA: cyclic tracking
            SUBPJJ => MCGDSCA
            SUBDS2 => MOCDS2
            SUBDSP => MOCDSP
            IF(LEXA) THEN
*           ACA: Exact Exponentials
              SUBDSC => MCGDS2E
            ELSE
*           ACA: Tabulated Exponentials
*           ACA: Exact Exponentials
              SUBDSC => MCGDS2A
            ENDIF
         ELSE
*        ACA: non-cyclic tracking
            SUBPJJ => MCGDSCA
            SUBDS2 => MCGDS2
            SUBDSP => MCGDSP
            IF(LEXA) THEN
*           ACA: Exact Exponentials
              SUBDSC => MCGDS2E
            ELSE
*           ACA: Tabulated Exponentials
              SUBDSC => MCGDS2A
            ENDIF
         ENDIF        
      ELSEIF(ISCH.EQ.0) THEN
*     PJJ/SCR: Diamond-Differencing Scheme
         IF(CYCLIC) THEN
*        ACA: cyclic tracking
            SUBPJJ => MCGDDDF
            SUBDS2 => MOCDS2
            SUBDSP => MOCDSP
            IF(LEXA) THEN
*           ACA: Exact Exponentials
              SUBDSC => MCGDS2E
            ELSE
*           ACA: Tabulated Exponentials
              SUBDSC => MCGDS2A
            ENDIF
         ELSE
*        ACA: non-cyclic tracking
            SUBPJJ => MCGDDDF
            SUBDS2 => MCGDS2
            SUBDSP => MCGDSP
            IF(LEXA) THEN
*           ACA: Exact Exponentials
              SUBDSC => MCGDS2E
            ELSE
*           ACA: Tabulated Exponentials
              SUBDSC => MCGDS2A
            ENDIF
         ENDIF
      ELSEIF(ISCH.EQ.-1) THEN
*     PJJ/SCR: Step-Characteristics Scheme with Exact Exponentials
         IF(CYCLIC) THEN
*        ACA: cyclic tracking
            SUBPJJ => MCGDSCE
            SUBDS2 => MOCDS2
            SUBDSP => MOCDSP
            IF(LEXA) THEN
*           ACA: Exact Exponentials
              SUBDSC => MCGDS2E
            ELSE
*           ACA: Tabulated Exponentials
              SUBDSC => MCGDS2A
            ENDIF
         ELSE
*        ACA: non-cyclic tracking
            SUBPJJ => MCGDSCE
            SUBDS2 => MCGDS2
            SUBDSP => MCGDSP
            IF(LEXA) THEN
*           ACA: Exact Exponentials
              SUBDSC => MCGDS2E
            ELSE
*           ACA: Tabulated Exponentials
              SUBDSC => MCGDS2A
            ENDIF
         ENDIF
      ENDIF
      CALL MCGASM(SUBPJJ,SUBDS2,SUBDSP,SUBDSC,IPTRK,KPSYS,IMPX,IFTRAK,
     1 NANI,NGEFF,NBCDA,NFI,NREG,NLONG,NBMIX,NMU,NANGL,NMAX,LC,NDIM,
     2 NGIND,CYCLIC,ISCR,CAZ0,CAZ1,CAZ2,CPO,LC0,PACA,LPS,LTMT,NPJJM,
     3 LACA,LPJJ,LPJJAN,SIGAL,LPRISM,N2REG,N2SOU,NZP,DELU,FACSYM,ISTRM)
*
   10 DEALLOCATE(SIGAL,KPSYS,NGIND,CPO,CAZ2,CAZ1,CAZ0)
      RETURN
      END
