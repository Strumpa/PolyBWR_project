*DECK SALTCG
      SUBROUTINE SALTCG(IPTRK ,IFTRK ,IPRINT,IGTRK ,NBSLIN, GG)
*
*-----------------------------------------------------------------------
*
*Purpose:
* To track an assembly of cells containing clusters using the new SALT
* tracking procedure (based on NXTTCG.f).
*
*Copyright:
* Copyright (C) 2014 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s):
* A. Hebert and G. Marleau
*
*Parameters: input
* IPTRK   pointer to the TRACKING data structure in
*         update or creation mode.
* IFTRK   pointer to the TRACKING file in creation mode.
* IPRINT  print level.
* IGTRK   flag to generate the tracking file. In the case where
*         IGTRK=1, the tracking is performed and
*         used to evaluate the track normalisation factor and the
*         tracking file is generated. When IGTRK=0, the tracking is
*         still performed and used to evaluate the track normalisation
*         factor but the tracking file is not generated.
* NBSLIN  maximum number of segments in a single tracking line.
*         computed by default in SALTCG but limited to 100000
*         elements. This default value can be bypassed using
*         keyword NBSLIN.
* GG      geometry basic information.
*
*----------
*
      USE              GANLIB
      USE SAL_GEOMETRY_TYPES, ONLY : T_G_BASIC
      IMPLICIT         NONE
*----
*  Subroutine arguments
*----
      TYPE(C_PTR)      IPTRK
      INTEGER          IFTRK
      INTEGER          IPRINT,IGTRK
      INTEGER          NBSLIN
      TYPE(T_G_BASIC)  GG
*----
*  Local parameters
*----
      INTEGER          IOUT
      CHARACTER        NAMSBR*6
      PARAMETER       (IOUT=6,NAMSBR='SALTCG')
      INTEGER          NSTATE
      PARAMETER       (NSTATE=40)
      DOUBLE PRECISION PI,DZERO,DONE,DTWO,DSUM
      PARAMETER       (PI=3.14159265358979, DZERO=0.0D0,DONE=1.0D0,
     >                 DTWO=2.0D0)
*----
*  Functions
*----
      INTEGER          KDROPN,IFTEMP,KDRCLS,ICLS
*----
*  Local variables
*----
      INTEGER          ISTATE(NSTATE),IEDIMG(NSTATE),ICODE(6)
      REAL             RSTATT(NSTATE),ALBEDO(6)
      INTEGER          RENO,LTRK,AZMOAQ,ISYMM,POLQUA,POLOAQ,AZMQUA,
     >                 AZMNBA
      DOUBLE PRECISION DENUSR,RCUTOF,DENLIN,SPACLN,WEIGHT
      DOUBLE PRECISION RADIUS,CENTER(3)
      INTEGER          NDIM,ITYPBC,IDIRG,NBOCEL,NBUCEL,IDIAG,
     >                 ISAXIS(3),NOCELL(3),NUCELL(3),MXMSH,MAXMSH,
     >                 MAXREG,NBTCLS,MAXMSP,MAXRSP,NFSUR,NFREG,
     >                 MXGSUR,NUNK,NPLANE,NPOINT,NTLINE,NBTDIR,
     >                 MAXSUB,MAXSGL,NBDR,NUNKF,ILONG,ITYLCM
      INTEGER          IPER(3)
      INTEGER          JJ,KK,NCOR,NQUAD,NANGL,NBANGL,LINMAX
      DOUBLE PRECISION DQUAD(4),ABSC(3,2),RCIRC,SIDEH,ANGLE
      CHARACTER        CTRK*4,COMENT*80
      INTEGER          IFMT,NEREG,NESUR
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: KEYMRG,MATALB
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NBSANG
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: SURVOL,DGMESH,
     > DNSANG,DDANG,DVNOR,DSNOR
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: DDENWT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: DANGLT
*----
*  Processing starts:
*  print routine opening header if required
*  and initialize various parameters.
*----
      IF(IPRINT .GE. 1) WRITE(IOUT,6000) NAMSBR
*----
*  Open temporary tracking file if required
*----
      IF(IGTRK .EQ. 1) THEN
        IFTEMP= KDROPN('DUMMYSQ',0,2,0)
        IF(IFTEMP .LE. 0) WRITE(IOUT,9010) NAMSBR
      ENDIF
*----
*  Get state vectors
*----
      ISTATE(:NSTATE)=0
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      CALL LCMGET(IPTRK,'EXCELTRACKOP',RSTATT)
      NEREG=ISTATE(1)
      NESUR=ISTATE(5)
      NUNK=NEREG+NESUR+1
      RENO  =ISTATE(8)
      LTRK  =ISTATE(9)
      AZMOAQ=ISTATE(11)
      ISYMM =ISTATE(12)
      POLQUA=ISTATE(13)
      POLOAQ=ISTATE(14)
      AZMQUA=ISTATE(15)
      AZMNBA=ISTATE(16)
      IFMT=ISTATE(21)
      DENUSR=DBLE(RSTATT(2))
      RCUTOF=DBLE(RSTATT(3))
      DENLIN=DBLE(RSTATT(4))
      SPACLN=DBLE(RSTATT(5))
      WEIGHT=RSTATT(6)
*----
*  Get main tracking records
*----
      CALL LCMGET(IPTRK,'ICODE       ',ICODE )
      CALL LCMGET(IPTRK,'ALBEDO      ',ALBEDO)
      CALL LCMSIX(IPTRK,'NXTRecords  ',1)
*----
*  Get general dimensioning vector for geometry tracking
*----
      IEDIMG(:NSTATE)=0
      CALL LCMGET(IPTRK,'G00000001DIM',IEDIMG)
      NDIM     =IEDIMG( 1)
      ITYPBC   =IEDIMG( 2)
      IDIRG    =IEDIMG( 3)
      NBOCEL   =IEDIMG( 4)
      NBUCEL   =IEDIMG( 5)
      IDIAG    =IEDIMG( 6)
      ISAXIS(1)=IEDIMG( 7)
      ISAXIS(2)=IEDIMG( 8)
      ISAXIS(3)=IEDIMG( 9)
      NOCELL(1)=IEDIMG(10)
      NOCELL(2)=IEDIMG(11)
      NOCELL(3)=IEDIMG(12)
      NUCELL(1)=IEDIMG(13)
      NUCELL(2)=IEDIMG(14)
      NUCELL(3)=IEDIMG(15)
      MXMSH    =IEDIMG(16)
      MAXREG   =IEDIMG(17)
      NBTCLS   =IEDIMG(18)
      MAXMSP   =IEDIMG(20)
      MAXRSP   =IEDIMG(21)
      NFSUR    =IEDIMG(22)
      NFREG    =IEDIMG(23)
      MXGSUR   =IEDIMG(24)
      MAXMSH=MAX(1,MXMSH,MAXMSP,MAXREG)
      IF(IPRINT .GE. 1) THEN
        WRITE(IOUT,6011) NFREG,NEREG,NFSUR,NESUR
      ENDIF
      IF((ITYPBC .EQ. 0).OR.(ITYPBC .EQ. 2)) THEN
*----
*  Define Cell for periodicity
*  Cartesian Boundary
*----
        IPER(1)=2
        IF(ABS(ISAXIS(1)) .EQ. 3) IPER(1)=1
        IPER(2)=2
        IF(ABS(ISAXIS(2)) .EQ. 3) IPER(2)=1
        IPER(3)=2
        IF(ABS(ISAXIS(3)) .EQ. 3) IPER(3)=1
*----
*  Use intrinsic geometry symmetries
*  to simplify tracking unless
*  NOSY tracking option activated
*----
        IF(ISYMM .NE. 0) THEN
          ISYMM=0
          IF(ABS(ISAXIS(1)) .EQ. 1 .OR. ABS(ISAXIS(1)) .EQ. 2) THEN
*----
*  X SYMMETRY
*----
            ISYMM=2
          ENDIF
          IF(ABS(ISAXIS(2)) .EQ. 1 .OR. ABS(ISAXIS(2)) .EQ. 2) THEN
*----
*  Y SYMMETRY
*----
            ISYMM=4+2*ISYMM
          ENDIF
          IF(NDIM .EQ. 3) THEN
            IF(ABS(ISAXIS(3)) .EQ. 1 .OR. ABS(ISAXIS(3)) .EQ. 2) THEN
*----
*  Z SYMMETRY
*----
              ISYMM=16+ISYMM
            ENDIF
          ENDIF
          IF(ISYMM .EQ. 0) ISYMM=1
        ENDIF
      ENDIF
*----
*  Read global mesh for geometry
*----
      ALLOCATE(DGMESH((MAXMSH+2)*4))
      IF(ITYPBC.EQ.0) THEN
        CALL NXTXYZ(IPTRK ,IPRINT,NDIM  ,ITYPBC,MAXMSH,NUCELL,
     >              ABSC,DGMESH)
      ELSE IF(ITYPBC.GE.2) THEN
        ! hexagonal geometry
        CALL LCMGET(IPTRK,'G00000001SMX',DGMESH)
        SIDEH=DGMESH(2)-DGMESH(1)
        ABSC(:2,1)=SIDEH*SQRT(3.0)
        ABSC(:2,2)=ABSC(:2,1)
      ELSE
        CALL XABORT(NAMSBR//': geometry not implemented')
      ENDIF
*----
*  Verify tracking parameters and compute number of angles
*  associated with angular order and spatial quadrature parameters
*  1. Isotropic tracking
*----
      NCOR= 1
      NPLANE=1
      IF(LTRK .EQ. 0) THEN
        IF(NDIM .EQ. 3) THEN
          CALL XABORT(NAMSBR//': 3-D geometry is not allowed')
        ENDIF
        NQUAD=2
        DQUAD(1)=DONE
        DQUAD(2)=DONE
        NANGL=AZMOAQ
        NBANGL=NANGL
        IF(ISYMM .EQ. 2 .OR. ISYMM .EQ. 8) THEN
          DQUAD(1)=DONE/DTWO
          DQUAD(2)=DZERO
        ENDIF
        ALLOCATE(DANGLT(NDIM,NBANGL,NQUAD),DDENWT(NBANGL,NQUAD))
        CALL NXTQAS(IPRINT,NDIM  ,AZMQUA,NANGL ,NQUAD ,NBANGL,
     >              DQUAD ,DANGLT,DDENWT)
        LINMAX=NBUCEL*(4+MXGSUR+16)
*----
*  Select standard spatial tracking parameters
*----
        CALL NXTQSS(IPRINT,NDIM  ,ITYPBC,MAXMSH,NUCELL,DENUSR,
     >              DGMESH,NPLANE,NPOINT,DENLIN,SPACLN,WEIGHT,
     >              RADIUS,CENTER)
        CALL LCMPUT(IPTRK,'TrackingDirc',NDIM*NBANGL*NQUAD,4,DANGLT)
      ELSE
*----
*  2. Specular tracking
*----
        NPOINT=0
        NQUAD=2
        IF(NDIM .EQ. 3) CALL XABORT(NAMSBR//
     >  ': TSPC option not valid for 3-D geometries')
        IF(ITYPBC .EQ. 0) THEN
          ! Cartesian geometries
          IF(AZMOAQ .GT. 24) THEN
            WRITE(IOUT,9002) NAMSBR,AZMOAQ,24,30
            AZMOAQ=30
          ELSE IF(AZMOAQ .GT. 20) THEN
            IF(AZMOAQ .NE. 24) THEN
              WRITE(IOUT,9003) NAMSBR,AZMOAQ,20,24,24
              AZMOAQ=24
            ENDIF
          ELSE IF(AZMOAQ .GT. 18) THEN
            IF(AZMOAQ .NE. 20) THEN
              WRITE(IOUT,9003) NAMSBR,AZMOAQ,18,20,20
              AZMOAQ=20
            ENDIF
          ELSE IF(AZMOAQ .GT. 14) THEN
            IF(AZMOAQ .NE. 18) THEN
              WRITE(IOUT,9003) NAMSBR,AZMOAQ,14,18,18
              AZMOAQ=18
            ENDIF
          ELSE IF(AZMOAQ .GT. 12) THEN
            IF(AZMOAQ .NE. 14) THEN
              WRITE(IOUT,9003) NAMSBR,AZMOAQ,12,14,14
              AZMOAQ=14
            ENDIF
          ELSE IF(AZMOAQ .GT. 8) THEN
            IF(AZMOAQ .NE. 12) THEN
              WRITE(IOUT,9003) NAMSBR,AZMOAQ,8,12,12
              AZMOAQ=12
            ENDIF
          ELSE IF(AZMOAQ .GT. 6) THEN
            IF(AZMOAQ .NE. 8) THEN
              WRITE(IOUT,9003) NAMSBR,AZMOAQ,6,8,8
              AZMOAQ=8
            ENDIF
          ELSE IF(AZMOAQ .GT. 2) THEN
            IF(AZMOAQ .NE. 6) THEN
              WRITE(IOUT,9003) NAMSBR,AZMOAQ,2,6,6
              AZMOAQ=6
            ENDIF
          ELSE IF(AZMOAQ .GE. 0) THEN
            IF(AZMOAQ .NE. 2) THEN
              WRITE(IOUT,9003) NAMSBR,AZMOAQ,-1,2,2
              AZMOAQ=2
            ENDIF
          ENDIF
        ELSE IF(ITYPBC .GE.2) THEN
          ! hexagonal geometries
          IF(AZMOAQ .GT. 12) THEN
            WRITE(IOUT,9002) NAMSBR,AZMOAQ,12,18
            AZMOAQ=18
          ELSE IF(AZMOAQ .GT. 6) THEN
            IF(AZMOAQ .NE. 12) THEN
              WRITE(IOUT,9003) NAMSBR,AZMOAQ,6,12,12
              AZMOAQ=12
            ENDIF
          ELSE IF(AZMOAQ .GT. 3) THEN
            IF(AZMOAQ .NE. 6) THEN
              WRITE(IOUT,9003) NAMSBR,AZMOAQ,3,6,6
              AZMOAQ=6
            ENDIF
          ELSE IF(AZMOAQ .GT. 1) THEN
            IF(AZMOAQ .NE. 3) THEN
              WRITE(IOUT,9003) NAMSBR,AZMOAQ,1,3,3
              AZMOAQ=3
            ENDIF
          ELSE IF(AZMOAQ .GE. 0) THEN
            IF(AZMOAQ .NE. 1) THEN
              WRITE(IOUT,9003) NAMSBR,AZMOAQ,-1,1,1
              AZMOAQ=1
            ENDIF
          ENDIF
        ENDIF
        NBANGL=AZMOAQ
        NANGL =AZMOAQ
        ALLOCATE(NBSANG(5,NBANGL))
        ALLOCATE(DANGLT(NDIM,NBANGL,4),DDENWT(NBANGL,4),DNSANG(NBANGL),
     >  DDANG(NBANGL))
        LINMAX=8*NANGL*NBUCEL*(4+MXGSUR+16)
        RCIRC=SQRT(ABSC(1,1)**2+ABSC(2,1)**2)
        ABSC(1,1)= ABSC(1,1)/RCIRC
        ABSC(2,1)= ABSC(2,1)/RCIRC
        CALL NXTQAC(IPRINT,NDIM  ,NANGL ,NBANGL,ITYPBC,DENUSR,
     >              ABSC  ,RCIRC ,AZMQUA,IPER  ,DANGLT,DDENWT,
     >              DNSANG,NBSANG,DDANG)
        DEALLOCATE(DDANG)
        DO JJ=1,NBANGL
          DANGLT(1,NBANGL-JJ+1,2)=-DANGLT(1,JJ,1)
          DANGLT(2,NBANGL-JJ+1,2)=DANGLT(2,JJ,1)
          DDENWT(NBANGL-JJ+1,2)=DDENWT(JJ,1)
        ENDDO
        DO JJ=1,NBANGL
          DANGLT(1,NBANGL-JJ+1,4)=DANGLT(1,JJ,1)
          DANGLT(2,NBANGL-JJ+1,4)=-DANGLT(2,JJ,1)
          DDENWT(NBANGL-JJ+1,4)=DDENWT(JJ,1)
          DANGLT(1,NBANGL-JJ+1,3)=DANGLT(1,JJ,2)
          DANGLT(2,NBANGL-JJ+1,3)=-DANGLT(2,JJ,2)
          DDENWT(NBANGL-JJ+1,3)=DDENWT(JJ,2)
        ENDDO
        IF(IPRINT.GT.1) THEN
          WRITE(IOUT,'(/34H SALTCG: CYCLIC ANGULAR QUADRATURE/8X,
     >    5HANGLE,8X,7HCOSINES,16(1h-),3X,6HWEIGHT,11X,
     >    17HCYCLIC PARAMETERS)')
          DSUM=0.D0
          DO KK=1,4
            DO JJ=1,NANGL
              ANGLE=SIGN(ACOS(DANGLT(1,JJ,KK))/PI*180.0,DANGLT(2,JJ,KK))
              WRITE(IOUT,'(1X,I4,1P,4E13.4,5X,2I4)') (KK-1)*NANGL+JJ,
     >        ANGLE,DANGLT(:2,JJ,KK),0.5D0/DDENWT(JJ,KK),NBSANG(:2,JJ)
              DSUM=DSUM+0.5D0/DDENWT(JJ,KK)
            ENDDO
          ENDDO
          WRITE(IOUT,'(39X,5HDSUM=,1P,E13.4)') DSUM
        ENDIF
        CALL LCMPUT(IPTRK,'TrackingDirc',NDIM*NBANGL*4,4,DANGLT)
        CALL LCMPUT(IPTRK,'TrackingTrkW',NBANGL*4,4,DDENWT)
        CALL LCMPUT(IPTRK,'TrackingSpaD',NBANGL,4,DNSANG)
        CALL LCMPUT(IPTRK,'TrackingNbST',5*NBANGL,1,NBSANG)
      ENDIF
      RSTATT(4)=REAL(DENLIN)
      RSTATT(5)=REAL(SPACLN)
      RSTATT(6)=REAL(WEIGHT)
      RSTATT(7)=REAL(RADIUS)
      RSTATT(8)=REAL(CENTER(1))
      RSTATT(9)=REAL(CENTER(2))
*----
*  Track
*----
      LINMAX=MAX(LINMAX,NBSLIN)
      IF(IPRINT .GE. 10) WRITE(IOUT,6010) LINMAX
      NBDR=1
      IF(RENO .EQ. -1) THEN
        IF(LTRK .EQ. 0) THEN
          NBDR=NQUAD*NBANGL+1
        ELSE IF(LTRK .EQ. 1) THEN
          NBDR=4*NBANGL+1
        ENDIF
      ENDIF
      ALLOCATE(DVNOR(NFREG*NBDR),DSNOR(NFSUR*NQUAD*NBANGL))
      IF(LTRK .EQ. 0) THEN
*----
*  Standard (isotropic) tracking (white boundary conditions)
*----
        IF(IPRINT .GE. 1) WRITE(IOUT,6030) NBANGL*NQUAD,NPOINT
        MAXSUB=1
        CALL SALTLS(IFTEMP,IPRINT,IGTRK ,NFREG ,NBANGL, NQUAD ,
     >              RENO  ,NBDR  ,IFMT  ,DENUSR,DANGLT, DDENWT,
     >              GG    ,NBTDIR,MAXSGL,NTLINE,DVNOR )
      ELSE
*----
*  Cyclic (specular) tracking (mirror like boundary conditions)
*----
        NBTDIR=0
        CALL SALTLC(IFTEMP,IPRINT,IGTRK,NDIM,NFREG,NBANGL,RENO,NBDR,
     >              IFMT,DENUSR,DANGLT,DDENWT,NBSANG,GG,MAXSUB,MAXSGL,
     >              NTLINE,DVNOR )
      ENDIF
*----
*  Save track normalisation vector
*----
      CALL LCMPUT(IPTRK,'NumMerge    ',NFREG,1,GG%NUM_MERGE)
      CALL LCMPUT(IPTRK,'VolMerge    ',NFREG,4,GG%VOL_NODE)
      CALL LCMPUT(IPTRK,'VTNormalize ',NFREG,4,DVNOR)
      IF(NBDR .GT. 1) THEN
        CALL LCMPUT(IPTRK,'VTNormalizeD',NFREG*(NBDR-1),4,
     >              DVNOR(NFREG+1))
      ENDIF
*----
*  Get cell description of geometry
*----
      NUNKF=NFREG+NFSUR+1
      ALLOCATE(KEYMRG(NUNKF),MATALB(NUNKF),SURVOL(NUNKF))
      CALL LCMLEN(IPTRK,'KEYMRG      ',ILONG,ITYLCM)
      IF(ILONG>NUNKF) CALL XABORT('SALTCG: NUNKF OVERLOW.')
      CALL LCMGET(IPTRK,'KEYMRG      ',KEYMRG)
      CALL LCMGET(IPTRK,'MATALB      ',MATALB)
      CALL LCMGET(IPTRK,'SAreaRvolume',SURVOL)
*----
*  Build NXTRecords directory
*----
      CALL LCMSIX(IPTRK,'NXTRecords  ',2)
      ISTATE(12)=ISYMM
      ISTATE(14)=POLOAQ
      ISTATE(17)=NPOINT
      ISTATE(18)=LINMAX
      ISTATE(19)=NTLINE
      ISTATE(20)=NBTDIR
      IF(LTRK .EQ. 0) THEN
        ISTATE(21)=NQUAD*NBANGL
      ELSE IF(LTRK .EQ. 1) THEN
        ISTATE(21)=4*NBANGL
      ENDIF
      ISTATE(22)=NPLANE
      CALL LCMPUT(IPTRK,'STATE-VECTOR',NSTATE,1,ISTATE)
      CALL LCMPUT(IPTRK,'EXCELTRACKOP',NSTATE,2,RSTATT)
*----
*  Renormalize tracks if required and transfer to final tracking file
*----
      IF(IGTRK .EQ. 1) THEN
        CTRK  = '$TRK'
        WRITE(IFTRK) CTRK,5,NTLINE,IFMT
        COMENT='CREATOR     : DRAGON'
        WRITE(IFTRK) COMENT
        COMENT='MODULE      : SALTCG'
        WRITE(IFTRK) COMENT
        COMENT='TYPE        : CARTESIAN'
        WRITE(IFTRK) COMENT
        IF(RENO .EQ. -1) THEN
          COMENT='TRKNOR      : Directional '
        ELSE IF(RENO .EQ. 0) THEN
          COMENT='TRKNOR      : Global      '
        ELSE
          COMENT='TRKNOR      : Off         '
        ENDIF
        WRITE(IFTRK) COMENT
        IF(IFMT .EQ. 1) THEN
          COMENT='OPTION      : Extended    '
          WRITE(IFTRK) COMENT
        ELSE
          COMENT='OPTION      : Short       '
          WRITE(IFTRK) COMENT
        ENDIF
*----
*  Compress VOLSUR and MATALB according to KEYMRG and save on IFTRK
*----
        IF(LTRK .EQ. 0) THEN
          WRITE(IFTRK) NDIM,LTRK,NEREG,NESUR,6,NCOR,NQUAD*NBANGL,MAXSUB,
     >    MAXSGL
        ELSE IF(LTRK .EQ. 1) THEN
          WRITE(IFTRK) NDIM,LTRK,NEREG,NESUR,6,NCOR,4*NBANGL,MAXSUB,
     >    MAXSGL
        ENDIF
        KEYMRG(NFSUR+2:NUNKF)=GG%NUM_MERGE(:NFREG)
        SURVOL(NFSUR+2:NUNKF)=GG%VOL_NODE(:NFREG)
        MATALB(NFSUR+2:NUNKF)=GG%MED(:NFREG)
        CALL NXTCVM(IFTRK,IPRINT,NFREG,NFSUR,NEREG,NESUR,MATALB,SURVOL,
     >              KEYMRG)
        WRITE(IFTRK) ( ICODE(JJ),JJ=1,6)
        WRITE(IFTRK) (ALBEDO(JJ),JJ=1,6)
        IF(LTRK .EQ. 0) THEN
          CALL NXTSQD(IFTRK,IPRINT,NDIM,NQUAD,NBANGL,DANGLT,DDENWT)
        ELSE IF(LTRK .EQ. 1) THEN
          WRITE(IFTRK) ((DANGLT(1,JJ,KK),DANGLT(2,JJ,KK),JJ=1,NBANGL),
     >    KK=1,4)
          WRITE(IFTRK) ((DDENWT(JJ,KK),JJ=1,NBANGL),KK=1,4)
        ENDIF
        REWIND IFTEMP
        CALL NXTTNS(IFTRK ,IFTEMP,IPRINT,RENO  ,NFSUR ,NFREG ,NDIM  ,
     >              MAXSUB,MAXSGL,NTLINE,NBDR  ,IFMT  ,KEYMRG,DVNOR)
*----
*  Close temporary tracking file if required
*----
        ICLS=KDRCLS(IFTEMP,2)
        IF(ICLS .NE. 0) WRITE(IOUT,9011) NAMSBR
      ENDIF
*----
*  Deallocate memory
*----
      DEALLOCATE(DSNOR,DVNOR)
      IF(LTRK .EQ. 1) DEALLOCATE(NBSANG,DNSANG)
      DEALLOCATE(DDENWT,DANGLT)
      DEALLOCATE(DGMESH,SURVOL,MATALB,KEYMRG)
*----
*  Processing finished:
*  print routine closing output header if required
*  and return
*----
      IF(IPRINT .GE. 1) THEN
        WRITE(IOUT,6012)
        WRITE(IOUT,6001) NAMSBR
      ENDIF
      RETURN
*----
*  Output formats
*----
 6000 FORMAT('(* Output from --',A6,'-- follows ')
 6001 FORMAT('   Output from --',A6,'-- completed *)')
 6010 FORMAT(' Maximum length of a line =',I10)
 6011 FORMAT(' Tracking of geometry begins:'/
     >       ' Number of regions before merge =',I10/
     >       ' Number of regions after merge  =',I10/
     >       ' Number of surfaces before merge=',I10/
     >       ' Number of surfaces after merge =',I10)
 6012 FORMAT(' Tracking of geometry completed')
 6030 FORMAT(' Number of directions for tracking = ',I10/
     >       ' Number of lines per direction     = ',I10)
 9002 FORMAT(' ***** Warning in ',A6,' *****'/
     >       '       Number of specular angles requested :',I10/
     >       '       For values > ',I10,' use ',I10)
 9003 FORMAT(' ***** Warning in ',A6,' *****'/
     >       '       Number of specular angles requested :',I10/
     >       '       For values > ',I10,' and < ',I10,' use ',I10)
 9010 FORMAT(' ***** Warning in ',A6,' *****'/
     >       '       Impossible to open temporary tracking file ')
 9011 FORMAT(' ***** Warning in ',A6,' *****'/
     >       '       Impossible to close temporary tracking file ')
      END
