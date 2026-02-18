*DECK LIBECC
      SUBROUTINE LIBECC(IPDRL,NGRO,IL,AWR,ENER,DELTA,DELECC,IGECCO,
     > SCAT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Construct the scattering matrix using analytical scattering kernels.
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
* IPDRL   pointer to the draglib (L_DRAGLIB signature).
* NGRO    number of energy groups.
* IL      Legendre order (=0: isotropic kernel in LAB).
* AWR     mass ratio for current isotope.
* ENER    energy limits of the coarse groups.
* DELTA   lethargy widths of the coarse groups.
* DELECC  lethargy widths of eccolib libraries.
* IGECCO  number of equal-width lethargy groups with eccolib libraries.
* IMPX    print flag.
*
*Parameters: output
* SCAT    scattering matrix.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPDRL
      INTEGER NGRO,IL,IGECCO
      REAL AWR,ENER(NGRO+1),DELTA(NGRO),DELECC,SCAT(NGRO,NGRO)
*----
*  LOCAL VARIABLES
*----
      PARAMETER(MAXEDI=47,MAXTRA=10000)
      CHARACTER CM*2
      CHARACTER(LEN=8), SAVE, DIMENSION(MAXEDI) :: NAMEDI=
     >  (/ 'NELAS   ','N2N     ','N3N     ','NNP     ','N4N     ',
     >     'NINEL001','NINEL002','NINEL003','NINEL004','NINEL005',
     >     'NINEL006','NINEL007','NINEL008','NINEL009','NINEL010',
     >     'NINEL011','NINEL012','NINEL013','NINEL014','NINEL015',
     >     'NINEL016','NINEL017','NINEL018','NINEL019','NINEL020',
     >     'NINEL021','NINEL022','NINEL023','NINEL024','NINEL025',
     >     'NINEL026','NINEL027','NINEL028','NINEL029','NINEL030',
     >     'NINEL031','NINEL032','NINEL033','NINEL034','NINEL035',
     >     'NINEL036','NINEL037','NINEL038','NINEL039','NINEL040',
     >     'NINEL041','NINEL   '/)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NJJ,IJJ
      REAL, ALLOCATABLE, DIMENSION(:) :: GAR,PRI,STIS,UUU,QQ
      REAL, ALLOCATABLE, DIMENSION(:,:) :: SSS2
      LOGICAL, ALLOCATABLE, DIMENSION(:) :: LPAR
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(NJJ(NGRO),IJJ(NGRO),GAR(NGRO*NGRO))
      ALLOCATE(PRI(MAXTRA),STIS(NGRO),UUU(NGRO))
      ALLOCATE(LPAR(MAXEDI),SSS2(NGRO,MAXEDI),QQ(MAXEDI))
*----
*  RECOVER CROSS SECTIONS CONTRIBUTING TO THE SCATTERING MATRIX
*----
      SSS2(:NGRO,:MAXEDI)=0.0
      LPAR(:MAXEDI)=.FALSE.
      QQ(:MAXEDI)=0.0
      DO I=1,MAXEDI
        CALL LCMLEN(IPDRL,NAMEDI(I),LENGT,ITYLCM)
        IF(LENGT.GT.0) THEN
          LPAR(I)=.TRUE.
          CALL LCMGET(IPDRL,NAMEDI(I),SSS2(1,I))
          DO IG1=1,NGRO
            IF(NAMEDI(I).EQ.'N2N') THEN
              SSS2(IG1,I)=2.0*SSS2(IG1,I)
            ELSE IF(NAMEDI(I).EQ.'N3N') THEN
              SSS2(IG1,I)=3.0*SSS2(IG1,I)
            ELSE IF(NAMEDI(I).EQ.'N4N') THEN
              SSS2(IG1,I)=4.0*SSS2(IG1,I)
            ENDIF
          ENDDO
          DO IG1=NGRO,1,-1
            IF(SSS2(IG1,I).NE.0.0) EXIT
            QQ(I)=-ENER(IG1)
          ENDDO
        ENDIF
      ENDDO
*----
*  CONSTRUCT THE SCATTERING MATRIX
*----
      WRITE (CM,'(I2.2)') IL
      CALL LCMGET(IPDRL,'NJJS'//CM,NJJ)
      CALL LCMGET(IPDRL,'IJJS'//CM,IJJ)
      LENGT=0
      DO IG1=1,NGRO
        LENGT=LENGT+NJJ(IG1)
      ENDDO
      GAR(:LENGT)=0.0
      CALL LCMGET(IPDRL,'SCAT'//CM,GAR)
      UUU(1)=DELTA(1)
      DO IG1=2,NGRO
        UUU(IG1)=UUU(IG1-1)+DELTA(IG1)
      ENDDO
      IGAR=0
      SCAT(:NGRO,:NGRO)=0.0
      DO IG1=1,IGECCO
        DO I=1,MAXEDI
          IF(LPAR(I)) THEN
            IF(NAMEDI(I).EQ.'NELAS') THEN
              CALL LIBPRI(MAXTRA,DELECC,AWR,0,IL,NPRI,PRI)
            ELSE
              ! treshold reaction
              IF(ENER(IG1).LE.-QQ(I)*(AWR+1.0)/AWR) CYCLE
              CALL LIBPRQ(MAXTRA,DELECC,AWR,ENER(IG1),QQ(I),0,IL,
     >        NPRI,PRI)
            ENDIF
            DO IPRI=1,NPRI
              IG2=IG1+IPRI-1 ! IG2 <-- IG1
              IF(IG2.GT.IGECCO) EXIT
              SCAT(IG2,IG1)=SCAT(IG2,IG1)+PRI(IPRI)*SSS2(IG1,I)
            ENDDO
          ENDIF
        ENDDO
        IGAR=IGAR+NJJ(IG1)
      ENDDO ! IG1
      DO IG2=IGECCO+1,NGRO
        DO IG1=IJJ(IG2),IJJ(IG2)-NJJ(IG2)+1,-1
          IGAR=IGAR+1
          SCAT(IG2,IG1)=GAR(IGAR)
        ENDDO
      ENDDO ! IG2
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(QQ,SSS2,LPAR)
      DEALLOCATE(UUU,STIS,PRI)
      DEALLOCATE(GAR,IJJ,NJJ)
      RETURN
      END
