SUBROUTINE SALT(NENTRY,HENTRY,IENTRY,JENTRY,KENTRY)
!
!-----------------------------------------------------------------------
!
!Purpose:
! To analyze and track a geometry data structure using the Sanchez
! algorithm for a PIJ, MOC or multicell surfacic solution of the flux.
!Copyright:
! Copyright (C) 2014 Ecole Polytechnique de Montreal
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version
!
!Author(s):
! A. Hebert
!
!Parameters: input
! NENTRY  number of data structures transfered to this module.
! HENTRY  name of the data structures.
! IENTRY  data structure type where:
!     IENTRY=1 for LCM memory object;
!     IENTRY=2 for XSM file;
!     IENTRY=3 for sequential binary file;
!     IENTRY=4 for sequential ASCII file.
! JENTRY  access permission for the data structure where:
!     =0 for a data structure in creation mode;
!     =1 for a data structure in modifications mode;
!     =2 for a data structure in read-only mode.
! KENTRY  data structure pointer.
!
!Comments:
! Instructions for the use of the SALT: module:
!   TRKFIL VOLTRK := SALT: SURFIL [ GEOMETRY ] :: (saltget) ;
!   where
!     TRKFIL   : sequential binary tracking file to be created
!     VOLTRK   : tracking data structure (signature L_TRACK)
!     SURFIL   : sequential ascii file used to store the surfacic
!    elements of the geometry.
!     GEOMETRY : optional geometry data structure used if BIHET is set
!    (signature L_GEOM)
!     (saltget): Processing options
!    (read from input using the NXTGET routine).
!
!-----------------------------------------------------------------------
!
  USE GANLIB
  USE SAL_GEOMETRY_TYPES,  ONLY : T_G_BASIC,LMERGM,IC,ISPEC,NANIS
  USE SALGET_FUNS_MOD
  IMPLICIT     NONE
  !----
  !  Subroutine arguments
  !----
  INTEGER      NENTRY,IENTRY(NENTRY),JENTRY(NENTRY)
  TYPE(C_PTR)  KENTRY(NENTRY)
  CHARACTER    HENTRY(NENTRY)*12
  !----
  !  Local parameters
  !----
  INTEGER      IOUT
  CHARACTER    NAMSBR*6
  PARAMETER   (IOUT=6,NAMSBR='SALT  ')
  INTEGER      ILCMUP,ILCMDN
  PARAMETER   (ILCMUP=1,ILCMDN=2)
  INTEGER      NSTATE,MAXENT
  PARAMETER   (NSTATE=40,MAXENT=2)
  INTEGER      N_DATAIN
  PARAMETER   (N_DATAIN=25)
  INTEGER      IUTYPE
  PARAMETER   (IUTYPE=2)
  !----
  !  Local variables
  !----
  TYPE(C_PTR)  IPGEO,IPTRK
  INTEGER      IMGEO,IFTRK,FGEO,NMACRO
  INTEGER      IGTRK
  INTEGER      IEN,ITC
  INTEGER      IQUA10,IBIHET
  CHARACTER    HSIGN*12
  INTEGER      ISTATT(NSTATE),DATAIN(N_DATAIN)
  REAL         RSTATT(NSTATE)
  CHARACTER    TITLE*72
  INTEGER      IPRINT
  INTEGER      NBSLIN
  INTEGER      ILONG,ITYLCM
  DOUBLE PRECISION RCUTOF
  INTEGER      OK
  !----
  !  Allocatable types
  !----
  TYPE(T_G_BASIC), ALLOCATABLE :: GG
  !----
  !  Validate entry parameters
  !----
  IF((NENTRY.LT.3).OR.(NENTRY.GT.4)) CALL XABORT(NAMSBR// &
  & ': Three or four data structures permitted')
  IPGEO=C_NULL_PTR
  FGEO=0
  IMGEO=0
  NBSLIN=100000
  !----
  !  Scan data structure to determine type and mode
  !----
  DO IEN=1,2
    IF(JENTRY(IEN).NE.0) CALL XABORT(NAMSBR// &
    & ': Object in creation mode expected')
    IF((IENTRY(IEN).EQ.1).OR.(IENTRY(IEN).EQ.2)) THEN
      IPTRK=KENTRY(IEN)
      HSIGN='L_TRACK     '
      CALL LCMPTC(IPTRK,'SIGNATURE',12,HSIGN)
      HSIGN='EXCELL'
      CALL LCMPTC(IPTRK,'TRACK-TYPE',12,HSIGN)
    ELSE IF(IENTRY(IEN).EQ.3) THEN
      IFTRK=FILUNIT(KENTRY(IEN))
    ENDIF
  ENDDO
  DO IEN=3,NENTRY
    IF(JENTRY(IEN).NE.2) CALL XABORT(NAMSBR// &
    & ': Object in read-only mode expected')
    IF((IENTRY(IEN).EQ.1).OR.(IENTRY(IEN).EQ.2)) THEN
      CALL LCMGTC(KENTRY(IEN),'SIGNATURE',12,HSIGN)
      IF(HSIGN.NE.'L_GEOM') CALL XABORT(NAMSBR// &
      & ': L_GEOM signature expected for '//HENTRY(IEN))
      IPGEO=KENTRY(IEN)
      IMGEO=-1
    ELSE IF(IENTRY(IEN).EQ.4) THEN
      FGEO=FILUNIT(KENTRY(IEN))
    ENDIF
  ENDDO
  IF(FGEO.EQ.0) CALL XABORT(NAMSBR// &
  &  ': The surfacic file is not defined')
  !----
  ! Initialize tracking parameters to 0
  !----
  ISTATT(:NSTATE)=0
  RSTATT(:NSTATE)=0.0
  !----
  ! Define default tracking options that are different from 0
  !----
  ISTATT(6)=1
  ISTATT(7)=4
  ISTATT(11)=1
  ISTATT(12)=-1
  ISTATT(13)=1
  ISTATT(15)=1
  ISTATT(16)=2
  ISTATT(22)=0
  ISTATT(23)=1
  IF(IMGEO .EQ. -1) THEN
    CALL LCMLEN(IPGEO,'BIHET',ILONG,ITYLCM)
    IF(ILONG.NE.0) ISTATT(40)=1
  ENDIF
  !----
  ! Recover processing method
  !----
  CALL SALGET(DATAIN,6,FGEO,IOUT,'dimensions for geometry')
  NMACRO=DATAIN(5)
  REWIND(FGEO)
  RSTATT(11)=1.0
  TITLE=' '
  IF(NMACRO.LE.1) THEN
    ISTATT(7)=4
  ELSE
    ISTATT(7)=5
  ENDIF
  CALL NXTGET(NSTATE,IPRINT,TITLE ,ISTATT,RSTATT,NBSLIN,IQUA10,IBIHET)
  IF((ISTATT(9).EQ.1).AND.(ISTATT(15).EQ.1)) THEN
    ISTATT(15)=8 ! replace EQW by EQW2
  ENDIF
  LMERGM=(ISTATT(26)==1)
  IF(IPRINT.GT.0) WRITE(IOUT,90) TITLE
  !----
  ! Save updated STATE-VECTOR, TITLE and EXCELL track options
  ! on tracking data structure
  !----
  CALL LCMPUT(IPTRK,'STATE-VECTOR',NSTATE,1,ISTATT)
  CALL LCMPUT(IPTRK,'EXCELTRACKOP',NSTATE,2,RSTATT)
  CALL LCMPTC(IPTRK,'TITLE',72,TITLE)
  !----
  ! Analyse geometry if required
  !----
  ALLOCATE(GG, STAT= OK)
  IF(OK /= 0) CALL XABORT('SALT: failure to allocate GG')
  RCUTOF=DBLE(RSTATT(3))
  !----
  !  Recover options from state vector
  !----
  NANIS=ISTATT(6)
  IC=ISTATT(7)
  ISPEC=ISTATT(9)
  IF(IPRINT>0) THEN
    IF(ISPEC==0) THEN
      WRITE(IOUT,*) 'SALT: isotropic boundary conditions'
    ELSE IF(ISPEC==1) THEN
      WRITE(IOUT,*) 'SALT: specular boundary conditions'
    ENDIF
  ENDIF
  !----
  !  Perform tracking
  !----
  IF(IC.EQ.4) THEN
    !----
    ! Track geometry for a PIJ or MOC solution
    !----
    IF(IPRINT>0) WRITE(IOUT,*) 'SALT: PIJ or MOC tracking'
    CALL SALACG(FGEO ,IPTRK, RCUTOF, IPRINT, GG)
    IF(ISTATT(9) .GE. 0 .AND. ISTATT(23) .EQ. 1) THEN
      IGTRK=1
      CALL SALTCG(IPTRK, IFTRK, IPRINT, IGTRK, NBSLIN, GG)
    ENDIF
  ELSE IF(IC.EQ.5) THEN
    !----
    ! Track geometry for a multicell surfacic solution
    !----
    IF(IPRINT>0) WRITE(IOUT,*) 'SALT: multicell surfacic tracking'
    IF(ISPEC.EQ.1) CALL XABORT('SALT: TSPC is forbidden with multicell surfacic tracking.')
    CALL SALMUS(FGEO ,IPTRK, IFTRK, RCUTOF, IPRINT, NBSLIN, GG)
  ELSE
    CALL XABORT('SALT: INVALID PROCESSING METHOD.')
  ENDIF
  !----
  ! Release allocated memory in SALT module
  !----
  CALL SALEND(GG)
  DEALLOCATE(GG, STAT= OK)
  IF(OK /= 0) CALL XABORT('SALT: failure to deallocate GG')
  !----
  !  Process double heterogeneity (BIHET) data (if available)
  !----
  IF(ISTATT(40) .NE. 0) THEN
     CALL XDRTBH(IPGEO,IPTRK,IQUA10,IBIHET,IPRINT,RSTATT(39))
  ENDIF
  !----
  ! Processing finished, return
  !----
  IF(IPRINT .GT. 1) THEN
     CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATT)
     WRITE(IOUT,100) (ISTATT(ITC),ITC=1,10)
     WRITE(IOUT,120) (ISTATT(ITC),ITC=11,22)
     WRITE(IOUT,130) ISTATT(23:24),ISTATT(26:30),ISTATT(40)
     IF(IC.EQ.5) WRITE(IOUT,140) RSTATT(12)
  ENDIF
  RETURN
  !----
  ! Formats
  !----
   90 FORMAT(/1H1,31H SSSSS    AA   LL     TTTTTTTT ,95(1H*)/ &
     & 32H SSSSSSS  AAAA  LL     TTTTTTTT ,58(1H*), &
     & 37H MULTIGROUP VERSION.  X. WARIN (2001)/ &
     & 28H SS   SS  AAAA  LL        TT/ &
     & 28H  SS     AA  AA LL        TT/ &
     & 28H    SS   AAAAAA LL        TT/ &
     & 28H SS   SS AAAAAA LL        TT/ &
     & 28H SSSSSSS AA  AA LLLLLLL   TT/ &
     & 28H  SSSSS  AA  AA LLLLLLL   TT//1X,A72/)
  100 FORMAT(/14H STATE VECTOR:/ &
     & 7H NREG  ,I9,22H   (NUMBER OF REGIONS)/ &
     & 7H KPN   ,I9,23H   (NUMBER OF UNKNOWNS)/ &
     & 7H ILK   ,I9,39H   (0=LEAKAGE PRESENT/1=LEAKAGE ABSENT)/ &
     & 7H NBMIX ,I9,36H   (MAXIMUM NUMBER OF MIXTURES USED)/ &
     & 7H NSURF ,I9,29H   (NUMBER OF OUTER SURFACES)/ &
     & 7H NANI  ,I9,48H   (1=P0 CROSS SECTIONS/2=P1 CROSS SECTIONS/...)/ &
     & 7H METHOD,I9,38H   (4=PIJ OR MOC/5=MULTICELL SURFACIC)/ &
     & 7H NORM  ,I9,48H   (NORMALIZATION OPTION 1=ABSENT/0=GLOBAL/-1=NO, &
     & 21HRMALIZATION BY ANGLE)/ &
     & 7H TRKT  ,I9,36H   (TRACKING TYPE 0=FINITE/1=CYCLIC)/ &
     & 7H BOUND ,I9,52H   (BOUNDARY CONDITIONS TYPE 0=ISOTROPIC/1=SPECULAR))
  120 FORMAT( &
     & 7H NANG  ,I9,30H   (NUMBER OF TRACKING ANGLES)/ &
     & 7H ASYM  ,I9,28H   (ANGULAR SYMMETRY FACTOR)/ &
     & 7H POLQUA,I9,32H   (POLAR ANGLE QUADRATURE TYPE)/ &
     & 7H POLOAQ,I9,33H   (POLAR ANGLE QUADRATURE ORDER)/ &
     & 7H AZMQUA,I9,47H   (AZIMUTHAL OR SOLID ANGULAR QUADRATURE TYPE)/ &
     & 7H NDIM  ,I9,25H   (NUMBER OF DIMENSIONS)/ &
     & 7H NPOINT,I9,40H   (NUMBER OF TRACKING POINTS ON A LINE)/ &
     & 7H MAXSGL,I9,30H   (MAXIMUM LENGTH OF A TRACK)/ &
     & 7H NTLINE,I9,37H   (TOTAL NUMBER OF TRACKS GENERATED)/ &
     & 7H NBTDIR,I9,47H   (TOTAL NUMBER OF TRACK DIRECTIONS PROCESSED)/ &
     & 7H NANGL ,I9,47H   (NUMBER OF TRACK DIRECTION ANGLES CONSIDERED, &
     & 20H IN THE INTEGRATION)/ &
     & 7H INSB  ,I9,25H   (VECTORIZATION OPTION))
  130 FORMAT( &
     & 7H ITRACK,I9,47H   (-1=MONTE-CARLO/0=DESACTIVATES TRACKING FILE, &
     & 39H BUILD/1=ACTIVATES TRACKING FILE BUILD)/ &
     & 7H NMACRO,I9,31H   (NUMBER OF MACRO GEOMETRIES)/ &
     & 7H MERGMX,I9,32H   (0/1=MERGMIX ACTIVATION FLAG)/ &
     & 7H NBATCH,I9,41H   (NUMBER OF TRACKS IN EACH OPENMP CORE)/ &
     & 7H IJAT  ,I9,54H   (NUMBER OF ADDITIONAL INTERFACE CURRENT COMPONENTS)/ &
     & 7H NMIX  ,I9,53H   (NUMBER OF PERIMETER ELEMENTS IN MACRO GEOMETRIES)/ &
     & 7H NBCDA ,I9,44H   (NUMBER OF PERIMETERS IN GLOBAL GEOMETRY)/ &
     & 7H IBIHET,I9,46H   (0/1=DOUBLE HETEROGENEITY IS NOT/IS ACTIVE))
  140 FORMAT(5H EPSJ,1P,E11.2,3X,32H(FLUX-CURRENT ITERATION EPSILON))
END SUBROUTINE SALT
