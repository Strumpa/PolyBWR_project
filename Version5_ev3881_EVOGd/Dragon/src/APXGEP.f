*DECK APXGEP
      SUBROUTINE APXGEP(IPAPX,IPDEPL,IMPX,ITIM,NORIG,NPAR,MUPLET,LGNEW,
     1 NVPNEW,NCALAR)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To recover remaining global parameters. Update the parameter tree
* for a new elementary calculation in the Apex file.
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
* IPAPX   pointer to the Apex file.
* IPDEPL  pointer to the burnup object.
* IMPX    print parameter.
* ITIM    index of the current burnup step.
* NORIG   index of the elementary calculation associated to the
*         father node in the parameter tree.
* NPAR    number of global parameters.
* MUPLET  tuple of indices associated to each global parameter of the
*         elementary calculation.
* LGNEW   parameter modification flag (.TRUE. only if the I-th global
*         parameter has changed in the new elementary calculation).
* NCALAR  index of the old elementary calculation.
*
*Parameters: output
* NVPNEW  number of nodes in the global parameter tree.
* NCALAR  index of the new elementary calculation.
*
*-----------------------------------------------------------------------
*
      USE GANLIB
      USE hdf5_wrap
*----
*  SUBROUTINE ARGUMENTS
*----
      TYPE(C_PTR) IPAPX,IPDEPL
      INTEGER IMPX,ITIM,NORIG,NPAR,MUPLET(NPAR),NVPNEW,NCALAR
      LOGICAL LGNEW(NPAR)
*----
*  LOCAL VARIABLES
*----
      INTEGER RANK,TYPE,NBYTE,DIMSR(5)
      PARAMETER (NSTATE=40,MAXPAR=50)
      INTEGER IDATA(NSTATE)
      CHARACTER TEXT4*4,TEXT12*12,HSMG*131
      LOGICAL LGERR,COMTRE,LAST
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: IDEBAR,IARBVA,IORIGI,IVAL0,
     1 DIMS_APX
      INTEGER, ALLOCATABLE, DIMENSION(:) :: JDEBAR,JARBVA
      CHARACTER(LEN=8), ALLOCATABLE, DIMENSION(:) :: PARFMT
      CHARACTER(LEN=80), ALLOCATABLE, DIMENSION(:) :: PARNAM
*----
*  RECOVER INFORMATION FROM THE 'DIMSAP' PARAMETER LIST.
*----
      NVPNEW=0
      NVPO=0
      IF(hdf5_group_exists(IPAPX,"/paramtree")) THEN
        CALL hdf5_info(IPAPX,"/paramtree/TREEVAL",RANK,TYPE,NBYTE,DIMSR)
        IF(RANK.NE.99) NVPO=DIMSR(1)
      ENDIF
*----
*  RECOVER INFORMATION FROM THE 'paramdescrip' DIRECTORY.
*----
      IF(NPAR.EQ.0) RETURN
      CALL hdf5_read_data(IPAPX,"/paramdescrip/PARNAM",PARNAM)
      CALL hdf5_read_data(IPAPX,"/paramdescrip/PARFMT",PARFMT)
*----
*  RECOVER REMAINING GLOBAL PARAMETERS.
*----
      DO 10 IPAR=1,NPAR
      IF((PARNAM(IPAR).EQ.'Burnup').OR.(PARNAM(IPAR).EQ.'Time').OR.
     1   (PARNAM(IPAR).EQ.'Power').OR.(PARNAM(IPAR).EQ.'Exposure').OR.
     2   (PARNAM(IPAR).EQ.'Flux').OR.(PARNAM(IPAR).EQ.'Heavy')) THEN
*
*        RECOVER GLOBAL PARAMETER VALUES FROM THE DEPLETION OBJECT.
         IF(.NOT.C_ASSOCIATED(IPDEPL)) CALL XABORT('APXGEP: NO DEPLETI'
     1   //'ON OBJECT AVAILABLE AMONG THE RHS LCM OBJECTS.')
         CALL LCMGET(IPDEPL,'STATE-VECTOR',IDATA)
         NBURN=IDATA(3)
         NBISO=IDATA(4)
         NREAC=IDATA(6)
         NVAR=IDATA(7)
         NBMIX=IDATA(8)
         CALL APXGEM(IPDEPL,ITIM,PARNAM(IPAR),0,NBURN,NBMIX,NBISO,
     1   NREAC,NVAR,VALPAR)
      ELSE
         GO TO 10
      ENDIF
      IF(IMPX.GT.0) WRITE(6,100) TRIM(PARNAM(IPAR)),VALPAR
*
      CALL APXPAV(IPAPX,IPAR,NPAR,'FLOTTANT',VALPAR,NITMA,TEXT12,
     1 MUPLET(IPAR),LGNEW(IPAR))
   10 CONTINUE
      IF(IMPX.GT.2) THEN
         WRITE(6,110) (MUPLET(I),I=1,NPAR)
         WRITE(6,'(/)')
      ENDIF
      DO 15 I=1,NPAR
      IF(MUPLET(I).EQ.0) THEN
         WRITE(HSMG,'(33HAPXGEP: UNDEFINED MUPLET ELEMENT=,I6)') I
         CALL XABORT(HSMG)
      ENDIF
   15 CONTINUE
*----
*  INTRODUCE VALUES INTO GLOBAL PARAMETER TREE.
*----
**
** Parameter tree: this tree has a number of stages equal to the
** number of parameters. For each value of the i-th parameter, we
** find the position in the tree corresponding to the value of the
** (i+1)-th parameter.
** NCALAR  Number of elementary calculations stored in the tree.
** NVP     Number of nodes in the parameter tree, including the root.
**         The value corresponding to the root is not used.
** DEBTREE - If the node does not correspond to the last parameter:
**           index in DEBTREE of the first daughter of the node.
**         - If the node correspond to the last parameter: index in
**           DEBTREE where we recover the index of an elementary
**           calculation.
** TREVAL  Index of the corresponding parameter in the 'pval'//n
**         record.
*
**     EXEMPLE:   dn = value in DEBTREE,  (m) = value in TREVAL
**
**     Root                          *(0)
**                                     !
**     Param. Nb 1                  d2(1)
**                            -------------------
**                           !                   !
**     Param. Nb 2        d3(1)                4(2)
**                       ---------           ---------
**                      !         !         !    !    !
**     Param. Nb 3   d5(1)      6(3)     d7(1) 8(2) 9(3)   d10
**
**     Calculation Nb:  4         5         1    2    3
**
**     DEBTREE:     2  3  5  7 10  4  5  1  2  3
**     TREVAL:      0  1  1  2  1  3  1  2  3
*
      IF(.NOT.hdf5_group_exists(IPAPX,"/paramtree/")) THEN
         MAXNVP=100*(NPAR+1)
         ALLOCATE(IDEBAR(MAXNVP+1),IARBVA(MAXNVP))
         IDEBAR(:MAXNVP+1)=0
         IARBVA(:MAXNVP)=0
         IARBVA=0
         DO 20 I=1,NPAR
         IDEBAR(I)=I+1
         IARBVA(I+1)=1
   20    CONTINUE
         IDEBAR(NPAR+1)=NPAR+2
         IDEBAR(NPAR+2)=1
         NCALAR=1
         NVPNEW=NPAR+1
         CALL hdf5_create_group(IPAPX,'paramtree')
      ELSE
        CALL hdf5_info(IPAPX,"/paramtree/TREEVAL",RANK,TYPE,NBYTE,DIMSR)
        MAXNVP=DIMSR(1)
*
*        Find position of the new point and create new PARBRE.
*
*        "II" is the order number of first parameter which recives a
*        "brand new" value.
*        COMTRE returns .TRUE. if the sweep throught the tree reaches
*        its bottom, otherwise it returns "KK" value: level of the
*        first new node to be introduced.
*
         CALL hdf5_read_data(IPAPX,"/paramtree/DEBTREE",JDEBAR)
         CALL hdf5_read_data(IPAPX,"/paramtree/TREEVAL",JARBVA)
         DO 30 IPAR=1,NPAR
         IF(LGNEW(IPAR)) THEN
            II=IPAR
            GO TO 40
         ENDIF
   30    CONTINUE
         II=NPAR+1
   40    LGERR=COMTRE(NPAR,NVPO,JARBVA,JDEBAR,MUPLET,KK,I0,IORD,JJ,LAST)
         IF((II.GT.NPAR).AND.LGERR) THEN
            WRITE(TEXT4,'(I4)') IORD
            CALL XABORT('APXGEP: ELEMENTARY CALCULATION HAS THE SAME'//
     1      ' GLOBAL PARAMETERS AS ELEMENTARY CALCULATION NB '//TEXT4)
         ENDIF
*
*        Size of the new tree.
*
         NVPNEW=NVPO+NPAR+1-MIN(II,KK)
         IF(NVPNEW.GT.MAXNVP) MAXNVP=NVPNEW+MAXNVP
         ALLOCATE(IDEBAR(MAXNVP+1),IARBVA(MAXNVP))
         IDEBAR(NVPNEW+2:MAXNVP+1)=0
         IARBVA(NVPNEW+1:MAXNVP)=0
*
*        Update values and suppress old PARBRE.
*
         CALL COMARB(NPAR,NVPO,NVPNEW,JDEBAR,JARBVA,LGNEW,MUPLET,NCALAR,
     1   IDEBAR,IARBVA)
         DEALLOCATE(JARBVA,JDEBAR)
      ENDIF
      CALL hdf5_write_data(IPAPX,"/paramtree/DEBTREE",IDEBAR(:NVPNEW+1))
      CALL hdf5_write_data(IPAPX,"/paramtree/TREEVAL",IARBVA(:NVPNEW))
      DEALLOCATE(IARBVA,IDEBAR)
      IF(NCALAR.EQ.1) THEN
         MAXNCA=1000
         ALLOCATE(IORIGI(MAXNCA))
         IORIGI(:MAXNCA)=0
      ELSE
         CALL hdf5_get_shape(IPAPX,"/paramtree/ORIGIN",DIMS_APX)
         MAXNCA=DIMS_APX(1)
         DEALLOCATE(DIMS_APX)
         IF(NCALAR.GT.MAXNCA) MAXNCA=NCALAR+MAXNCA
         ALLOCATE(IORIGI(MAXNCA))
         IORIGI(:MAXNCA)=0
         CALL hdf5_read_data(IPAPX,"/paramtree/ORIGIN",IVAL0)
         IORIGI(:MAXNCA)=IVAL0(:MAXNCA)
         DEALLOCATE(IVAL0)
      ENDIF
      IORIGI(NCALAR)=NORIG
      CALL hdf5_write_data(IPAPX,"/paramtree/ORIGIN",IORIGI(:NCALAR))
      DEALLOCATE(IORIGI,PARFMT,PARNAM)
      RETURN
*
  100 FORMAT(31H APXGEP: SET GLOBAL PARAMETER ',A,3H' =,1P,E12.4)
  110 FORMAT(/16H APXGEP: MUPLET=,10I6:/(16X,10I6))
      END
