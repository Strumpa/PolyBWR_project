SUBROUTINE DOORS_SN1(IPTRK,NANIS,NREG,NMAT,NUN,MATCOD,SIGG,SUNKNO)
  !
  !-----------------------------------------------------------------------
  !
  !Purpose:
  ! Compute the source for the solution of SN equations with a flat flux
  ! approximation
  !
  !Copyright:
  ! Copyright (C) 2025 Ecole Polytechnique de Montreal
  ! This library is free software; you can redistribute it and/or
  ! modify it under the terms of the GNU Lesser General Public
  ! License as published by the Free Software Foundation; either
  ! version 2.1 of the License, or (at your option) any later version
  !
  !Author(s): A. Hebert and C. Bienvenue
  !
  !Parameters: input
  ! IPTRK   pointer to the tracking LCM object.
  ! NANIS   maximum cross section Legendre order (=0: isotropic).
  ! NREG    number of regions.
  ! NMAT    number of mixtures.
  ! NUN     number of unknowns per energy group including spherical
  !         harmonic terms and boundary SN fluxes.
  ! MATCOD  mixture indices.
  ! SIGG    cross section.
  !
  !Parameters: output
  ! SUNKNO  sources.
  !
  !-----------------------------------------------------------------------
  !
  USE GANLIB
  !----
  !  SUBROUTINE ARGUMENTS
  !----
  TYPE(C_PTR) IPTRK
  INTEGER NANIS,NREG,NMAT,NUN,MATCOD(NREG)
  REAL SIGG(0:NMAT,NANIS+1),SUNKNO(NUN)
  !----
  !  LOCAL VARIABLES
  !----
  PARAMETER(NSTATE=40)
  INTEGER JPAR(NSTATE),P,P2,ILP
  !----
  !  ALLOCATABLE ARRAYS
  !----
  TYPE(C_PTR) IL_PTR,IM_PTR
  INTEGER, POINTER, DIMENSION(:) :: IL,IM
  !----
  !  RECOVER SNT SPECIFIC PARAMETERS.
  !----
  CALL LCMGET(IPTRK,'STATE-VECTOR',JPAR)
  IF(JPAR(1).NE.NREG) CALL XABORT('DOORS_SN1: INCONSISTENT NREG.')
  IF(JPAR(2).NE.NUN) CALL XABORT('DOORS_SN1: INCONSISTENT NUN.')
  ITYPE=JPAR(6)
  NSCT=JPAR(7)
  IELEM=JPAR(8)
  ISCAT=JPAR(16)
  CALL LCMGPD(IPTRK,'IL',IL_PTR)
  CALL LCMGPD(IPTRK,'IM',IM_PTR)
  CALL C_F_POINTER(IL_PTR,IL,(/ NSCT /))
  CALL C_F_POINTER(IM_PTR,IM,(/ NSCT /))
  !----
  !  CONSTRUCT THE SOURCE. LOOP OVER LEGENDRE ORDERS.
  !----
  IOF0=0
  DO P=1,NSCT
    ILP=IL(P)
    IF(ILP.GT.MIN(ISCAT-1,NANIS)) CYCLE
    IF((ITYPE.EQ.2).OR.(ITYPE.EQ.4)) THEN
      !----
      !  SLAB OR SPHERICAL 1D CASE.
      !----
      DO IR=1,NREG
        IBM=MATCOD(IR)
        IF(IBM.LE.0) CYCLE
        IND=(IR-1)*NSCT*IELEM+IELEM*(P-1)+1
        SUNKNO(IND)=SUNKNO(IND)+SIGG(IBM,ILP+1)
      ENDDO ! IR
    ELSE IF(ITYPE.EQ.3) THEN
      !----
      !  CYLINDRICAL 1D CASE.
      !----
      DO P2=0,P-1
        IF(MOD((P-1)+P2,2).EQ.1) CYCLE
        IOF0=IOF0+1
        DO IR=1,NREG
          IBM=MATCOD(IR)
          IF(IBM.LE.0) CYCLE
          IND=(IR-1)*NSCT+IOF0
          SUNKNO(IND)=SUNKNO(IND)+SIGG(IBM,ILP+1)
        ENDDO ! IR
      ENDDO ! P2
    ELSE IF((ITYPE.EQ.5).OR.(ITYPE.EQ.6).OR.(ITYPE.EQ.8)) THEN
      !----
      !  2D CASES (CARTESIAN OR R-Z).
      !----
      NM=IELEM**2
      DO IR=1,NREG
        IBM=MATCOD(IR)
        IF(IBM.LE.0) CYCLE
        IND=(IR-1)*NSCT*NM+(P-1)*NM+1
        SUNKNO(IND)=SUNKNO(IND)+SIGG(IBM,ILP+1)
      ENDDO ! IR
    ELSE IF((ITYPE.EQ.7).OR.(ITYPE.EQ.9)) THEN
      !----
      ! 3D CARTESIAN CASE
      !----
      NM=IELEM**3
      DO IR=1,NREG
        IBM=MATCOD(IR)
        IF(IBM.LE.0) CYCLE
        IND=(IR-1)*NSCT*NM+(P-1)*NM+1
        SUNKNO(IND)=SUNKNO(IND)+SIGG(IBM,ILP+1)  
      ENDDO ! IR
    ELSE
      CALL XABORT('DOORS_SN1: TYPE OF DISCRETIZATION NOT IMPLEMENTED.')
    ENDIF
  ENDDO ! P
  RETURN
END SUBROUTINE DOORS_SN1
