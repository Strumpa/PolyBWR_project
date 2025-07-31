*DECK APXCAT
      SUBROUTINE APXCAT(IPAPX,IPRHS,NORIG,NPAR,NCAL,MUPCPO,LGNCPO,LWARN)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To catenate a RHS Apex file into the output Apex file.
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
* IPAPX   pointer to the output Apex file.
* IPRHS   pointer to the rhs Apex file (contains the new calculations).
* NORIG   index of the elementary calculation associated to the
*         father node in the parameter tree.
* NPAR    number of global parameters in the output Apex file.
* NCAL    initial number of calculations in LHS Apex file.
* MUPCPO  tuple of the new global parameters in the output Apex file.
* LGNCPO  LGNEW value of the new global parameters in the output
*         Apex file.
* LWARN   logical used in case if an elementary calculation in the RHS
*         is already present in Apex file. If LWARN=.true. a warning is
*         send and the Apex file values are kept otherwise XABORT is
*         called (default).
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE hdf5_wrap
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPAPX,IPRHS
      INTEGER NORIG,NPAR,NCAL,MUPCPO(NPAR)
      LOGICAL LGNCPO(NPAR),LWARN
*----
*  LOCAL VARIABLES
*----
      PARAMETER (MAXPAR=50)
      INTEGER RANK,TYPE,NBYTE,DIMSR(5)
      INTEGER MUPLET(2*MAXPAR),MUPRHS(2*MAXPAR)
      CHARACTER HSMG*131,RECNAM*80,RECNA2*80,TEXT4*4,TEXT12*12
      LOGICAL COMTRE,LGERR,LGNEW(MAXPAR)
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IORRHS,JDEBAR,JARBVA,VINTE,
     1 IDEBAR,IARBVA,IORIGI
      REAL, ALLOCATABLE, DIMENSION(:) :: VREAL
      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:) :: PARFMT_RHS
      CHARACTER(LEN=12), ALLOCATABLE, DIMENSION(:) :: VCHAR
      CHARACTER(LEN=80), ALLOCATABLE, DIMENSION(:) :: PARNAM_RHS,
     1 PARNAM_LHS
*
      IF(NPAR.GT.MAXPAR) CALL XABORT('APXCAT: MAXPAR OVERFLOW.')
      NGR=0
      CALL APXTOC(IPRHS,IMPX,NLAM,NREA,NBISO,NBMAC,NMILR,NPARR,NVPR,
     1 NISOF,NISOP,NISOS,NCALR,NGR,NISOTS,NSURFD,NPRC)
      IF(NCALR.EQ.0) THEN
        CALL XABORT('APXCAT: NO CALCULATION IN RHS APEX FILE.')
      ELSE IF(NPARR.GT.NPAR) THEN
        WRITE(HSMG,'(42HAPXCAT: ELEMENTARY CALCULATION WITH AN INV,
     1  31HALIB NB. OF GLOBAL PARAMETERS =,I7,3H GT,I7,1H.)') NPARR,
     2  NPAR
        CALL XABORT(HSMG)
      ENDIF
      NVPO=0 ! initial number of nodes in LHS Apex file
      CALL hdf5_read_data(IPAPX,"NCALS",NCAL)
      IF(NCAL.GT.0) THEN
        NG=0
        CALL APXTOC(IPAPX,0,NLAM,NREA,NBISO,NBMAC,NMIL,NPAR1,NVPO,
     1  NISOF,NISOP,NISOS,NCAL,NG,NISOTS,NSURFD,NPRC)
        IF(NGR.NE.NG) THEN
          WRITE(HSMG,'(42HAPXCAT: ELEMENTARY CALCULATION WITH AN INV,
     1    20HALIB NB. OF GROUPS =,I7,3H NE,I7,1H.)') NGR,NG
          CALL XABORT(HSMG)
        ELSE IF(NMILR.NE.NMIL) THEN
          WRITE(HSMG,'(42HAPXCAT: ELEMENTARY CALCULATION WITH AN INV,
     1    22HALIB NB. OF MIXTURES =,I7,3H NE,I7,1H.)') NMILR,NMIL
          CALL XABORT(HSMG)
        ELSE IF(NPAR1.NE.NPAR) THEN
          WRITE(HSMG,'(42HAPXCAT: ELEMENTARY CALCULATION WITH AN INV,
     1    31HALIB NB. OF GLOBAL PARAMETERS =,I7,3H NE,I7,1H.)') NPAR1,
     2    NPAR
          CALL XABORT(HSMG)
        ENDIF
      ENDIF
*----
*  MAIN LOOP OVER THE NCALR ELEMENTARY CALCULATIONS OF THE RHS APEX FILE
*----
      IDEM=0
      NCALS=NCAL
      DO 170 ICAL=1,NCALR
*----
*  COMPUTE THE MUPLET VECTOR FROM THE RHS APEX FILE
*----
      CALL hdf5_read_data(IPRHS,"/paramtree/DEBTREE",JDEBAR)
      CALL hdf5_read_data(IPRHS,"/paramtree/TREEVAL",JARBVA)
      CALL hdf5_read_data(IPRHS,"/paramtree/ORIGIN",IORRHS)
      DO 30 I=NVPR-NCALR+1,NVPR
      IF(JDEBAR(I+1).EQ.ICAL) THEN
         I0=I
         GO TO 40
      ENDIF
   30 CONTINUE
      CALL XABORT('APXCAT: MUPLET ALGORITHM FAILURE 1.')
   40 MUPRHS(NPAR)=JARBVA(I0)
      DO 65 IPAR=NPAR-1,1,-1
      DO 50 I=1,NVPR-NCALR
      IF(JDEBAR(I+1).GT.I0) THEN
         I0=I
         GO TO 60
      ENDIF
   50 ENDDO
      CALL XABORT('APXCAT: MUPLET ALGORITHM FAILURE 2.')
   60 MUPRHS(IPAR)=JARBVA(I0)
   65 CONTINUE
      DEALLOCATE(JARBVA,JDEBAR)
*----
*  RECOVER THE GLOBAL PARAMETERS
*----
      DO 70 I=1,NPAR
      MUPLET(I)=MUPCPO(I)
      LGNEW(I)=LGNCPO(I)
   70 CONTINUE
      CALL hdf5_read_data(IPAPX,"/paramdescrip/PARNAM",PARNAM_LHS)
      CALL hdf5_read_data(IPRHS,"/paramdescrip/PARFMT",PARFMT_RHS)
      CALL hdf5_read_data(IPRHS,"/paramdescrip/PARNAM",PARNAM_RHS)
      DO 100 IPAR=1,NPARR
         DO 80 I0=1,NPAR
         IF(PARNAM_RHS(IPAR).EQ.PARNAM_LHS(I0)) THEN
            IPARN=I0
            GO TO 90
         ENDIF
   80    CONTINUE
         CALL XABORT('APXCAT: UNABLE TO FIND '//PARNAM_RHS(IPAR)//'.')
   90    WRITE(RECNAM,'(17H/paramvalues/PVAL,I8)') IPAR
         IVAL=MUPRHS(IPAR)
         IF(PARFMT_RHS(IPAR).EQ.'FLOTTANT') THEN
            CALL hdf5_read_data(IPRHS,TRIM(RECNAM),VREAL)
            FLOTT=VREAL(IVAL)
            DEALLOCATE(VREAL)
         ELSE IF(PARFMT_RHS(IPAR).EQ.'ENTIER') THEN
            CALL hdf5_read_data(IPRHS,TRIM(RECNAM),VINTE)
            NITMA=VINTE(IVAL)
            DEALLOCATE(VINTE)
         ELSE IF(PARFMT_RHS(IPAR).EQ.'CHAINE') THEN
            CALL hdf5_read_data(IPRHS,TRIM(RECNAM),VCHAR)
            TEXT12=VCHAR(IVAL)
            DEALLOCATE(VCHAR)
         ENDIF
         CALL APXPAV(IPAPX,IPARN,NPAR,PARFMT_RHS(IPAR),FLOTT,NITMA,
     1   TEXT12,MUPLET(IPARN),LGNEW(IPARN))
  100 CONTINUE
      DEALLOCATE(PARNAM_RHS,PARFMT_RHS,PARNAM_LHS)
*----
*  UPDATE THE PARAMETER TREE IN THE OUTPUT APEX FILE
*----
      IF(NVPO.EQ.0) THEN
         MAXNVP=20*(NPAR+1)
         ALLOCATE(IDEBAR(MAXNVP+1),IARBVA(MAXNVP))
         IDEBAR(:MAXNVP+1)=0
         IARBVA(:MAXNVP)=0
         IARBVA=0
         DO 140 I=1,NPAR
         IDEBAR(I)=I+1
         IARBVA(I+1)=1
  140    CONTINUE
         IDEBAR(NPAR+1)=NPAR+2
         IDEBAR(NPAR+2)=1
         NCALS=1
         NVPNEW=NPAR+1
      ELSE
         CALL hdf5_read_data(IPAPX,"/paramtree/DEBTREE",JDEBAR)
         CALL hdf5_read_data(IPAPX,"/paramtree/TREEVAL",JARBVA)
         DO 150 IPAR=1,NPAR
         IF(LGNEW(IPAR)) THEN
            II=IPAR
            GO TO 160
         ENDIF
  150    CONTINUE
         II=NPAR+1
  160    LGERR=COMTRE(NPAR,NVPO,JARBVA,JDEBAR,MUPLET,KK,I0,IORD,JJ,
     1   LAST)
         IF((II.GT.NPAR).AND.LGERR) THEN
            WRITE(TEXT4,'(I4)') IORD
            IF(LWARN) THEN
               WRITE(6,*)'APXCAT: ELEMENTARY CALCULATION HAS THE ',
     1         'SAME PARAMETERS AS ELEMENTARY CALCULATION NB ',TEXT4
               DEALLOCATE(JARBVA,JDEBAR,IORRHS)
               IDEM=IDEM+1
               GOTO 170
            ELSE
               CALL XABORT('APXCAT: ELEMENTARY CALCULATION HAS THE '//
     1         'SAME PARAMETERS AS ELEMENTARY CALCULATION NB '//TEXT4)
            ENDIF
         ENDIF
*
*        Size of the new tree.
*
         NVPNEW=NVPO+NPAR+1-MIN(II,KK)
         MAXNVP=NVPR
         IF(NVPNEW.GT.MAXNVP) MAXNVP=NVPNEW+MAXNVP
         ALLOCATE(IDEBAR(MAXNVP+1),IARBVA(MAXNVP))
         IDEBAR(NVPNEW+2:MAXNVP+1)=0
         IARBVA(NVPNEW+1:MAXNVP)=0
*
*        Update values and suppress old PARBRE.
*
         CALL COMARB(NPAR,NVPO,NVPNEW,JDEBAR,JARBVA,LGNEW,MUPLET,NCALS,
     1   IDEBAR,IARBVA)
         DEALLOCATE(JARBVA,JDEBAR)
      ENDIF
      IF(NCALS.NE.NCAL+ICAL-IDEM) CALL XABORT('APXCAT: INVALID NCALS.')
      NVPO=NVPNEW
      CALL hdf5_write_data(IPAPX,"/NCALS",NCALS)
      CALL hdf5_write_data(IPAPX,"/paramtree/DEBTREE",IDEBAR(:NVPNEW+1))
      CALL hdf5_write_data(IPAPX,"/paramtree/TREEVAL",IARBVA(:NVPNEW))
      DEALLOCATE(IARBVA,IDEBAR)
      IF(NCALS.EQ.1) THEN
         MAXNCA=1000
         ALLOCATE(IORIGI(MAXNCA))
         IORIGI(:MAXNCA)=0
      ELSE
         CALL hdf5_info(IPAPX,"/paramtree/ORIGIN",RANK,TYPE,NBYTE,DIMSR)
         MAXNCA=DIMSR(1)
         IF(NCALS.GT.MAXNCA) MAXNCA=NCALS+MAXNCA
         ALLOCATE(IORIGI(MAXNCA))
         IORIGI(:MAXNCA)=0
         CALL hdf5_read_data(IPAPX,"/paramtree/ORIGIN",VINTE)
         IORIGI(:DIMSR(1))=VINTE(:DIMSR(1))
         DEALLOCATE(VINTE)
      ENDIF
      IF(IORRHS(ICAL).EQ.0) THEN
         IORIGI(NCALS)=NORIG
      ELSE
         IORIGI(NCALS)=NCAL+IORRHS(ICAL)
      ENDIF
      CALL hdf5_write_data(IPAPX,"/paramtree/ORIGIN",IORIGI(:NCALS))
      DEALLOCATE(IORIGI,IORRHS)
      IF(NCALS.NE.NCAL+ICAL-IDEM) CALL XABORT('APXCAT: INVALID NCALS.')
*----
*  RECOVER THE ELEMENTARY CALCULATION
*----
      WRITE(RECNAM,'(4Hcalc,I8)') NCALS
      WRITE(RECNA2,'(4Hcalc,I8)') ICAL
      call hdf5_copy(IPRHS,RECNA2,IPAPX,RECNAM) ! IPRHS -> IPAPX
  170 CONTINUE
* END OF LOOP ON ELEMENTARY CALCULATIONS. ********************
      RETURN
      END
