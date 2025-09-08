*DECK BIVFIS
      SUBROUTINE BIVFIS(IPTRK,NREG,NMAT,NIFIS,NUNKNO,NGRP,MATCOD,VOL,
     > XSCHI,XSNUF,FUNKNO,SUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the fission source for a BIVAC tracking.
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
* IPTRK   pointer to the tracking LCM object.
* NANIS   maximum cross section Legendre order.
* NREG    number of regions.
* NMAT    number of mixtures.
* NIFIS   number of fissile isotopes.
* NUNKNO  number of unknowns per energy group.
* NGRP    number of energy groups.
* MATCOD  mixture indices.
* VOL     volumes.
* XSCHI   fission spectra.
* XSNUP   nu times the fission cross sections.
* FUNKNO  fluxes.
*
*Parameters: output
* SUNKNO  sources.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE DOORS_MOD
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPTRK
      INTEGER NREG,NMAT,NIFIS,NUNKNO,NGRP,MATCOD(NREG)
      REAL VOL(NREG),XSCHI(0:NMAT,NIFIS,NGRP),XSNUF(0:NMAT,NIFIS,NGRP),
     1 FUNKNO(NUNKNO,NGRP),SUNKNO(NUNKNO,NGRP)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(NSTATE=40)
      INTEGER JPAR(NSTATE),IJ1(25),IJ2(25)
      LOGICAL CYLIND
      CHARACTER*12 CXDOOR
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KN,IDL
      REAL, ALLOCATABLE, DIMENSION(:) :: XX,DD,FXSOR
      REAL, ALLOCATABLE, DIMENSION(:,:) :: RR,RS
*----
*  RECOVER BIVAC SPECIFIC PARAMETERS.
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',JPAR)
      IF(JPAR(1).NE.NREG) CALL XABORT('BIVFIS: INCONSISTENT NREG.')
      IF(JPAR(2).NE.NUNKNO) CALL XABORT('BIVFIS: INCONSISTENT NUNKNO.')
      ITYPE=JPAR(6)
      IELEM=JPAR(8)
      ICOL=JPAR(9)
      ISPLH=JPAR(10)
      L4=JPAR(11)
      LX=JPAR(12)
      NLF=JPAR(14)
      ISCAT=JPAR(16)
      CYLIND=(ITYPE.EQ.3).OR.(ITYPE.EQ.6)
      IF(ICOL.EQ.4) THEN
        CALL XABORT('BIVFIS: COLLOCATION NODAL NOT IMPLEMENTED.')
      ELSE IF((ITYPE.NE.2).AND.(ITYPE.NE.5)) THEN
        CALL XABORT('BIVFIS: CARTESIAN 1D OR 2D GEOMETRY EXPECTED.')
      ENDIF
      CALL LCMLEN(IPTRK,'KN',MAXKN,ITYLCM)
      ALLOCATE(XX(NREG),DD(NREG),KN(MAXKN),IDL(NREG))
      CALL LCMGET(IPTRK,'XX',XX)
      CALL LCMGET(IPTRK,'DD',DD)
      CALL LCMGET(IPTRK,'KN',KN)
      CALL LCMGET(IPTRK,'KEYFLX',IDL)
      IF(IELEM.LT.0) THEN
        ! Lagrangian finite element method
*----
*  RECOVER THE FINITE ELEMENT UNIT STIFFNESS MATRIX.
*----
        CALL LCMSIX(IPTRK,'BIVCOL',1)
        CALL LCMLEN(IPTRK,'T',LC,ITYLCM)
        ALLOCATE(RR(LC,LC),RS(LC,LC))
        CALL LCMGET(IPTRK,'R',RR)
        CALL LCMGET(IPTRK,'RS',RS)
        CALL LCMSIX(IPTRK,' ',2)
*----
*  COMPUTE VECTORS IJ1 AND IJ2
*----
        LL=LC*LC
        DO I=1,LL
          IJ1(I)=1+MOD(I-1,LC)
          IJ2(I)=1+(I-IJ1(I))/LC
        ENDDO
*----
*  COMPUTE THE SOURCE
*----
        NUM1=0
        DO IR=1,NREG
          IBM=MATCOD(IR)
          IF(IBM.LE.0) CYCLE
          IF(VOL(IR).EQ.0.0) GO TO 10
          DO I=1,LL
            IND1=KN(NUM1+I)
            IF(IND1.EQ.0) CYCLE
            I1=IJ1(I)
            I2=IJ2(I)
            DO J=1,LL
              IND2=KN(NUM1+J)
              IF(IND2.EQ.0) CYCLE
              IF(CYLIND) THEN
                AUXX=(RR(I1,IJ1(J))+RS(I1,IJ1(J))*XX(IR)/DD(IR))*
     >          RR(I2,IJ2(J))*VOL(IR)
              ELSE
                AUXX=RR(I1,IJ1(J))*RR(I2,IJ2(J))*VOL(IR)
              ENDIF
              DO IG=1,NGRP
                DO JG=1,NGRP
                  DO IS=1,NIFIS
                    SIGG=XSCHI(IBM,IS,IG)*XSNUF(IBM,IS,JG)
                    SUNKNO(IND1,IG)=SUNKNO(IND1,IG)+AUXX*SIGG*
     >              FUNKNO(IND2,JG)
                  ENDDO ! IS
                ENDDO ! JG
              ENDDO ! IG
            ENDDO ! J
          ENDDO ! I
   10     NUM1=NUM1+LL
        ENDDO ! IR
        DEALLOCATE(RS,RR)
        ! append the integrated volumic sources
        ALLOCATE(FXSOR(NUNKNO))
        DO IS=1,NIFIS
          FXSOR(:NUNKNO)=0.0
          DO IR=1,NREG
            IBM=MATCOD(IR)
            IF(IBM.LE.0) CYCLE
            DO IG=1,NGRP
              FXSOR(IDL(IR))=FXSOR(IDL(IR))+VOL(IR)*XSNUF(IBM,IS,IG)*
     >        FUNKNO(IDL(IR),IG)
            ENDDO ! IG
            DO IG=1,NGRP
              SUNKNO(IDL(IR),IG)=SUNKNO(IDL(IR),IG)+XSCHI(IBM,IS,IG)*
     >        FXSOR(IDL(IR))
            ENDDO ! IG
          ENDDO ! IR
        ENDDO ! IS
        DEALLOCATE(FXSOR)
      ELSE
        ! Raviart-Thomas finite element method
        CXDOOR='BIVAC'
        ALLOCATE(FXSOR(NUNKNO))
        DO IS=1,NIFIS
          FXSOR(:NUNKNO)=0.0
          DO IG=1,NGRP
            CALL DOORS(CXDOOR,IPTRK,NMAT,0,NUNKNO,XSNUF(0,IS,IG),
     >      FXSOR,FUNKNO(1,IG))
          ENDDO
          DO IR=1,NREG
            IBM=MATCOD(IR)
            IF(IBM.EQ.0) CYCLE
            DO IE=1,IELEM**2
              IND=IDL(IR)+IE-1
              IF(IND.EQ.0) CYCLE
              DO IG=1,NGRP
                SUNKNO(IND,IG)=SUNKNO(IND,IG)+XSCHI(IBM,IS,IG)*
     >          FXSOR(IND)
              ENDDO
            ENDDO
          ENDDO
        ENDDO ! IS
        DEALLOCATE(FXSOR)
      ENDIF
      DEALLOCATE(IDL,KN,DD,XX)
      RETURN
      END
