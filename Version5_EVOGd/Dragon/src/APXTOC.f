*DECK APXTOC
      SUBROUTINE APXTOC(IPAPX,IMPX,NLAM,NREA,NBISO,NBMAC,NMIL,NPAR,
     1 NVP,NISOF,NISOP,NISOS,NCAL,NGRP,NISOTS,NSURFD,NPRC)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Recover the table of content of an Apex file.
*
*Copyright:
* Copyright (C) 2021 Ecole Polytechnique de Montreal
*
*Author(s): 
* A. Hebert
*
*Parameters: input
* IPAPX   address of the Apex file.
* IMPX    print parameter (equal to zero for no print).
*
*Parameters: output
* NLAM    number of types of radioactive decay reactions
* NREA    number of neutron-induced reaction
* NBISO   number of particularized isotopes
* NBMAC   number of macroscopic sets
* NMIL    number of mixtures in the APEX
* NPAR    number of parameters
* NVP     number of nodes in the global parameter tree
* NISOF   number of particularized fissile isotopes
* NISOP   number of particularized fission products
* NISOS   number of particularized stable isotopes
* NCAL    number of elementary calculations
* NGRP    number of energy groups
* NISOTS  maximum number of isotopes in output tables
* NSURFD  number of discontinuity factors values in the Apex file
* NPRC    number of precursors
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE hdf5_wrap
      IMPLICIT NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPAPX
      INTEGER IMPX,NLAM,NREA,NBISO,NBMAC,NMIL,NPAR,NVP,NISOF,NISOP,
     1 NISOS,NCAL,NGRP,NISOTS,NSURFD,NPRC
*----
*  LOCAL VARIABLES
*----
      INTEGER, PARAMETER::IOUT=6
      INTEGER I,II,RANK,TYPE,NBYTE,DIMSR(5)
      CHARACTER HSMG*131,RECNAM*80
      CHARACTER(LEN=100), ALLOCATABLE, DIMENSION(:) :: LIST
      INTEGER, ALLOCATABLE, DIMENSION(:) :: DIMS_APX
      CHARACTER(LEN=4), ALLOCATABLE, DIMENSION(:) :: TYPISO
*----
*  LIST GROUPS AND DATASETS ON THE ROOT FILE
*----
      IF(IMPX.GT.0) THEN
        CALL hdf5_list_groups(IPAPX, '/', LIST)
        WRITE(*,*)
        WRITE(*,*) 'APXTOC: GROUP TABLE OF CONTENTS'
        DO I=1,SIZE(LIST)
          WRITE(*,*) TRIM(LIST(I))
        ENDDO
        DEALLOCATE(LIST)
        CALL hdf5_list_datasets(IPAPX, '/', LIST)
        WRITE(*,*)
        WRITE(*,*) 'APXTOC: DATASET TABLE OF CONTENTS'
        DO I=1,SIZE(LIST)
          WRITE(*,*) TRIM(LIST(I))
        ENDDO
        DEALLOCATE(LIST)
      ENDIF
*----
*  RECOVER APEX PARAMETERS
*----
      NMIL=1
      NGRP=0
      CALL hdf5_read_data(IPAPX,"NCALS",NCAL)
      IF(hdf5_group_exists(IPAPX,"/physconst/")) THEN
        CALL hdf5_info(IPAPX,"/physconst/ENRGS",RANK,TYPE,NBYTE,DIMSR)
        IF(TYPE.NE.99) THEN
          CALL hdf5_get_shape(IPAPX,"/physconst/ENRGS",DIMS_APX)
        ELSE
          GO TO 10
        ENDIF
      ELSE IF(hdf5_group_exists(IPAPX,"/physco001/")) THEN
        CALL hdf5_info(IPAPX,"/physc001/ENRGS",RANK,TYPE,NBYTE,DIMSR)
        IF(TYPE.NE.99) THEN
          CALL hdf5_get_shape(IPAPX,"/physc001/ENRGS",DIMS_APX)
        ELSE
          GO TO 10
        ENDIF
      ELSE
        CALL XABORT('APXTOC: GROUP physconst NOT FOUND IN HDF5 FILE.')
      ENDIF
      IF(NGRP.EQ.0) THEN
        NGRP=DIMS_APX(1)-1
      ELSE IF(NGRP.NE.DIMS_APX(1)-1) THEN
        WRITE(HSMG,'(46H APXTOC: THE APEX FILE HAS AN INVALID NUMBER O,
     1  17HF ENERGY GROUPS (,I4,3H VS,I5,2H).)') NGRP,DIMS_APX(1)-1
        CALL XABORT(HSMG)
      ENDIF
      DEALLOCATE(DIMS_APX)
  10  NBMAC=0
      NREA=0
      IF(hdf5_group_exists(IPAPX,"/explicit/")) THEN
        NBISO=0
        CALL hdf5_info(IPAPX,"/explicit/ISONAME",RANK,TYPE,NBYTE,DIMSR)
        IF(RANK.NE.99) NBISO=DIMSR(1)
        NBMAC=0
        CALL hdf5_info(IPAPX,"/explicit/MACNAME",RANK,TYPE,NBYTE,DIMSR)
        IF(RANK.NE.99) NBMAC=DIMSR(1)
        NREA=0
        CALL hdf5_info(IPAPX,"/explicit/REANAME",RANK,TYPE,NBYTE,DIMSR)
        IF(RANK.NE.99) NREA=DIMSR(1)
      ELSE IF(hdf5_group_exists(IPAPX,"/expli001/")) THEN
        NBISO=0
        CALL hdf5_info(IPAPX,"/expli001/ISONAME",RANK,TYPE,NBYTE,DIMSR)
        IF(RANK.NE.99) NBISO=DIMSR(1)
        NBMAC=0
        CALL hdf5_info(IPAPX,"/expli001/MACNAME",RANK,TYPE,NBYTE,DIMSR)
        IF(RANK.NE.99) NBMAC=DIMSR(1)
        NREA=0
        CALL hdf5_info(IPAPX,"/expli001/REANAME",RANK,TYPE,NBYTE,DIMSR)
        IF(RANK.NE.99) NREA=DIMSR(1)
      ELSE
        CALL XABORT('APXTOC: GROUP explicit NOT FOUND IN APEX FILE.')
      ENDIF
*----
*  SET NISOF AND NISOP
*----
      NISOF=0
      NISOP=0
      NISOS=0
      NSURFD=0
      IF(NBISO.GT.0) THEN
        IF(hdf5_group_exists(IPAPX,"/physconst/")) THEN
          CALL hdf5_get_shape(IPAPX,"/physconst/ISOTA",DIMS_APX)
        ELSE IF(hdf5_group_exists(IPAPX,"/physco001/")) THEN
          CALL hdf5_get_shape(IPAPX,"/physc001/ISOTA",DIMS_APX)
        ENDIF
        IF(DIMS_APX(1).NE.NBISO) THEN
          WRITE(HSMG,'(44H APXTOC: INCONSISTENT number of ISOTOPES IN ,
     1    31Hexplicit AND physconst GROUPS (,I4,3H VS,I5,2H).)') NBISO,
     2    DIMS_APX(1)
          CALL XABORT(HSMG)
        ENDIF
        DEALLOCATE(DIMS_APX)
        IF(hdf5_group_exists(IPAPX,"/physconst/")) THEN
          CALL hdf5_read_data(IPAPX,"/physconst/ISOTYP",TYPISO)
        ELSE IF(hdf5_group_exists(IPAPX,"/physco001/")) THEN
          CALL hdf5_read_data(IPAPX,"/physc001/ISOTYP",TYPISO)
        ENDIF
        DO I=1,NBISO
          IF(TYPISO(I).EQ.'FISS') NISOF=NISOF+1
          IF(TYPISO(I).EQ.'F.P.') NISOP=NISOP+1
          IF(TYPISO(I).EQ.'OTHE') NISOS=NISOS+1
        ENDDO
        DEALLOCATE(TYPISO)
      ENDIF
      IF(NCAL.EQ.0) GO TO 20
*----
*  SET DECAYC, NVALUE AND TREEVAL
*----
      NLAM=0
      NISOTS=0
      IF(hdf5_group_exists(IPAPX,"/physconst/")) THEN
        CALL hdf5_info(IPAPX,"/physconst/DECAYC",RANK,TYPE,NBYTE,DIMSR)
      ELSE IF(hdf5_group_exists(IPAPX,"/physco001/")) THEN
        CALL hdf5_info(IPAPX,"/physc001/DECAYC",RANK,TYPE,NBYTE,DIMSR)
      ENDIF
      IF(RANK.NE.99) THEN
        NLAM=DIMSR(1)
        NISOTS=DIMSR(2)
      ENDIF
      NPAR=0
      CALL hdf5_info(IPAPX,"/paramdescrip/NVALUE",RANK,TYPE,NBYTE,DIMSR)
      IF(RANK.NE.99) NPAR=DIMSR(1)
      NVP=0
      IF(hdf5_group_exists(IPAPX,"/paramtree")) THEN
        CALL hdf5_info(IPAPX,"/paramtree/TREEVAL",RANK,TYPE,NBYTE,DIMSR)
        IF(RANK.NE.99) NVP=DIMSR(1)
      ENDIF
*----
*  SET NSURFD
*----
      RECNAM='calc       1/miscellaneous/'
      CALL hdf5_info(IPAPX,TRIM(RECNAM)//"ADF",RANK,TYPE,NBYTE,DIMSR)
      IF(TYPE.NE.99) NSURFD=NSURFD+DIMSR(1)
      CALL hdf5_info(IPAPX,TRIM(RECNAM)//"CPDF",RANK,TYPE,NBYTE,DIMSR)
      IF(TYPE.NE.99) NSURFD=NSURFD+DIMSR(1)
      CALL hdf5_info(IPAPX,TRIM(RECNAM)//"INTERNAL_ADF",RANK,TYPE,NBYTE,
     1 DIMSR)
      IF(TYPE.NE.99) NSURFD=NSURFD+DIMSR(1)
      CALL hdf5_info(IPAPX,TRIM(RECNAM)//"INTERNAL_CPDF",RANK,TYPE,
     1 NBYTE,DIMSR)
      IF(TYPE.NE.99) NSURFD=NSURFD+DIMSR(1)
*----
*  SET NPRC
*----
      NPRC=0
      CALL hdf5_list_groups(IPAPX, "calc       1", LIST)
      DO I=1,SIZE(LIST)
        IF(TRIM(LIST(I)).EQ.'kinetics') THEN
          RECNAM='calc       1/kinetics/LAMBDA'
          CALL hdf5_info(IPAPX,TRIM(RECNAM)//"INTERNAL_ADF",RANK,TYPE,
     1    NBYTE,DIMSR)
          IF(TYPE.NE.99) NPRC=DIMSR(1)
          EXIT
        ENDIF
*----
*  SET NMIL
*---- 
        IF(LIST(I)(:2).EQ.'xs') THEN
          IF(LEN(LIST(I)).EQ.2) CYCLE
          READ(LIST(I),'(2X,I8)') II
          NMIL=MAX(II,NMIL)
        ENDIF
      ENDDO
      DEALLOCATE(LIST)
*----
*  PRINT APEX PARAMETERS
*----
   20 IF(IMPX.GT.0) THEN
        WRITE(IOUT,'(/38H APXTOC: table of content information:)')
        WRITE(IOUT,'(32H   nb of radioactive reactions =,I3)') NLAM
        WRITE(IOUT,'(36H   nb of neutron-induced reactions =,I3)') NREA
        WRITE(IOUT,'(34H   nb of particularized isotopes =,I4)') NBISO
        WRITE(IOUT,'(27H   nb of macroscopic sets =,I2)') NBMAC
        WRITE(IOUT,'(19H   nb of mixtures =,I5)') NMIL
        WRITE(IOUT,'(28H   nb of global parameters =,I4)') NPAR
        WRITE(IOUT,'(38H   nb of nodes in the parameter tree =,I4)') NVP
        WRITE(IOUT,'(42H   nb of particularized fissile isotopes =,I4)')
     1  NISOF
        WRITE(IOUT,'(42H   nb of particularized fission products =,I4)')
     1  NISOP
        WRITE(IOUT,'(41H   nb of particularized stable isotopes =,I4)')
     1  NISOS
        WRITE(IOUT,'(23H   nb of calculations =,I9)') NCAL
        WRITE(IOUT,'(24H   nb of energy groups =,I4)') NGRP
        WRITE(IOUT,'(44H   maximum nb of isotopes in output tables =,
     1  I4)') NISOTS
        WRITE(IOUT,'(39H   nb of discontinuity factors values =,I4)')
     1  NSURFD
        WRITE(IOUT,'(21H   nb of precursors =,I4/)') NPRC
      ENDIF
      RETURN
      END
