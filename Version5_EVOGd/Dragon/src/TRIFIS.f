*DECK TRIFIS
      SUBROUTINE TRIFIS(IPTRK,NREG,NMAT,NIFIS,NUNKNO,NGRP,MATCOD,VOL,
     > XSCHI,XSNUF,FUNKNO,SUNKNO)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Compute the fission source for a TRIVAC tracking.
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
      INTEGER JPAR(NSTATE),IJ1(125),IJ2(125),IJ3(125)
      CHARACTER*12 CXDOOR
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KN,IDL
      REAL, ALLOCATABLE, DIMENSION(:) :: FXSOR
      REAL, ALLOCATABLE, DIMENSION(:,:) :: RR,RS
*----
*  RECOVER TRIVAC SPECIFIC PARAMETERS.
*----
      CALL LCMGET(IPTRK,'STATE-VECTOR',JPAR)
      IF(JPAR(1).NE.NREG) CALL XABORT('TRIFIS: INCONSISTENT NREG.')
      IF(JPAR(2).NE.NUNKNO) CALL XABORT('TRIFIS: INCONSISTENT NUNKNO.')
      ITYPE=JPAR(6)
      IELEM=JPAR(9)
      ICOL=JPAR(10)
      ISCAT=JPAR(32)
      IF(ICOL.EQ.4) THEN
        CALL XABORT('TRIFIS: COLLOCATION NODAL NOT IMPLEMENTED.')
      ELSE IF((ITYPE.NE.2).AND.(ITYPE.NE.5).AND.(ITYPE.NE.7)) THEN
        CALL XABORT('TRIFIS: CARTESIAN 1D, 2D OR 3D GEOMETRY EXPECTED.')
      ENDIF
      CALL LCMLEN(IPTRK,'KN',MAXKN,ITYLCM)
      ALLOCATE(KN(MAXKN),IDL(NREG))
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
        LL=LC*LC*LC
        DO L=1,LL
          L1=1+MOD(L-1,LC)
          L2=1+(L-L1)/LC
          L3=1+MOD(L2-1,LC)
          IJ1(L)=L1
          IJ2(L)=L3
          IJ3(L)=1+(L2-L3)/LC
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
            I3=IJ3(I)
            DO J=1,LL
              IND2=KN(NUM1+J)
              IF(IND2.EQ.0) CYCLE
              AUXX=RR(I1,IJ1(J))*RR(I2,IJ2(J))*RR(I3,IJ3(J))*VOL(IR)
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
        CXDOOR='TRIVAC'
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
            DO IE=1,IELEM**3
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
      DEALLOCATE(IDL,KN)
      RETURN
      END
