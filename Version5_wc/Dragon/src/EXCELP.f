*DECK EXCELP
      SUBROUTINE EXCELP(  IPTRK, IFTRAK, IPRNTP,  NSOUT,  NREG,  NBMIX,
     >                   MATCOD, NRENOR, XSSIGT,  IPIJK, N2PRO,   NSBG,
     >                    NPSYS, NBATCH, TITREC,  NALBP,  ALBP, MATALB,
     >                   VOLSUR,  DPROB, DPROBX )
*
*-----------------------------------------------------------------------
*
*Purpose:
* Calculation of the collision probabilities for EXCELL.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy
*
*Parameters: input
* IPTRK   pointer to the tracking (L_TRACK signature).
* IFTRAK  unit of the sequential binary tracking file.
* IPRNTP  print flag (equal to zero for no print).
* NSOUT   number of surfaces.
* NREG    total number of merged blocks for which specific values
*         of the neutron flux and reactions rates are required.
* NBMIX   number of mixtures (NBMIX=max(MATCOD(i))).
* MATCOD  index number of the mixture type assigned to each volume.
* NRENOR  normalization scheme for PIJ matrices.
* XSSIGT  total macroscopic cross sections ordered by mixture.
* IPIJK   pij option (=1 pij, =4 pijk).
* N2PRO   number of terms in collision probability matrices, including
*         surface and volume contributions.
* NSBG    number of energy groups.
* NPSYS   non-converged energy group indices.
* NBATCH  number of tracks processed in each OpenMP core (default: =1).
* TITREC  title.
* NALBP   number of multigroup physical albedos.
* ALBP    multigroup physical albedos.
*
*Parameters: output
* MATALB  global mixture/albedo identification vector.
* VOLSUR  global surface volume vector.
* DPROB   collision probabilities.
* DPROBX  directional collision probabilities.
*
*-----------------------------------------------------------------------
*--------+---------------- R O U T I N E S -------------+--+-----------*
*  NAME  /                  DESCRIPTION                                *
*--------+-------------------------------------------------------------*
* CP INtegration
*   PIJI2D / TO INTEGRATE CP IN 2D GEOMETRIES (ISOTROPIC B.C.)
*   PIJI3D / TO INTEGRATE CP IN 3D GEOMETRIES (ISOTROPIC B.C.)
*   PIJS2D / TO INTEGRATE CP IN 2D GEOMETRIES (SPECULAR  B.C.)
*   PIJS3D / TO INTEGRATE CP IN 3D GEOMETRIES (SPECULAR  B.C.)
* CP Normalisation
*   PIJRDG / TO RENORMALIZE CP USING DIAGONAL COEFFICIENTS
*   PIJRGL / TO RENORMALIZE CP USING GELBARD HOMOGENEOUS SCHEME
*   PIJRNL / TO RENORMALIZE CP USING NON-LINEAR FACTORS
*   PIJRHL / TO RENORMALIZE CP USING HELIOS METHOD
* Various functions
*   PIJWPR / TO PRINT CP MATRICES IN SUM FORMAT
*   PIJCMP / COMPRESS CP MATRIX TO SYMETRIC FORMAT
* Inline tracking
*   NXTTGC / TRACK CYCLIC NXT LINE IN GEOMETRY
*   NXTTGS / TRACK STANDARD NXT LINE IN GEOMETRY
*   NXTXYZ / READ GEOMETRY LIMITS
*--------+-------------------------------------------------------------*
*
      USE               GANLIB
      IMPLICIT          NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      CHARACTER         TITREC*72
      TYPE(C_PTR)       IPTRK
      INTEGER           IFTRAK,IPRNTP,NSOUT,NREG,NBMIX,MATCOD(NREG),
     >                  NRENOR,IPIJK,N2PRO,NSBG,NPSYS(NSBG),NBATCH,
     >                  NALBP,MATALB(-NSOUT:NREG)
      REAL              XSSIGT(0:NBMIX,NSBG),ALBP(NALBP,NSBG),
     >                  VOLSUR(-NSOUT:NREG,NSBG)
      LOGICAL           SWNZBC
      DOUBLE PRECISION  DPROB(N2PRO,NSBG),DPROBX(N2PRO,NSBG)
*----
*  LOCAL VARIABLES
*----
      INTEGER           IOUT, ICPALL, ICPEND, MXGAUS, NSTATE
      PARAMETER       ( IOUT=6, ICPALL=4, ICPEND=3, MXGAUS=64,
     >                  NSTATE=40)
      CHARACTER         NAMSBR*6
      PARAMETER       ( NAMSBR='EXCELP')
      INTEGER           MKI1, MKI2, MKI3, MKI4, MKI5
      PARAMETER       (MKI1=600,MKI2=600,MKI3=600,MKI4=600,MKI5=600)
      INTEGER           ISTATE(NSTATE)
      INTEGER           NPROB,ISBG,KSBG,ITYPBC,NBCDA
      REAL              EXTKOP(NSTATE),CUTOF,RCUTOF,ASCRP,YGSS,
     >                  XGSS(MXGAUS),WGSS(MXGAUS),WGSSX(MXGAUS)
      LOGICAL           SWVOID,LPIJK
      CHARACTER         CTRKT*4, COMENT*80
      DOUBLE PRECISION  DANG0,DASCRP
*
      INTEGER           JJ,MSYM,IL,NALLOC,ITRAK,IANG,IC,IPRT,ISPEC,
     >                  IUN,KSPEC,LOPT,MXSEG,NALBG,NANGL,NCOMNT,NCOR,
     >                  NCORT,NDIM,NGSS,NREG2,NSCRP,NTRK,NUNKNO,JGSS,
     >                  JUN,IFMT,MXSUB,ISA,IBATCH,IL1,III,IND,I,J,
     >                  ITYLCM
*----
*  Variables for NXT: inline tracking
*----
      INTEGER           ILCMUP,ILCMDN
      PARAMETER        (ILCMUP=1,ILCMDN=2)
      DOUBLE PRECISION  DZERO,DONE,DTWO
      PARAMETER        (DZERO=0.0D0,DONE=1.0D0,DTWO=2.0D0)
      INTEGER           IEDIMG(NSTATE),NPOINT,NBUCEL,MXMSH,MAXPIN,
     >                  MXGSUR,MXGREG,MAXMSH,NPLANE,NUCELL(3)
      CHARACTER         NAMREC*12
*----
*  Allocatable arrays
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: ICODE
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DSV
      REAL, ALLOCATABLE, DIMENSION(:) :: ALBG,ALBEDO
      REAL, ALLOCATABLE, TARGET, DIMENSION(:,:) :: SIGTAL,SIGT00
      REAL, POINTER, DIMENSION(:,:) :: SIGT
*-- NXT TRACKING
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: IUNFLD
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DVNOR,DWGTRK
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: DGMESH,DANGLT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: DORITR
*-- Temporary arrays
      REAL, ALLOCATABLE, DIMENSION(:) :: LOPATH
      REAL, ALLOCATABLE, DIMENSION(:,:,:) :: SIGANG
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: STAYIN,GOSOUT
*-- Tracking file arrays
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NCOIL1,NSUB,NBSEG
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NUMERO
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: WEIGHT
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: LENGTH
*----
*  Common blocks for Bickley functions
*----
      INTEGER          L1, L2, L3, L4, L5
      REAL             PAS1,XLIM1,PAS2,XLIM2,PAS3,XLIM3,
     >                 PAS4,XLIM4,PAS5,XLIM5,BI1,BI2,BI3,BI4,BI5
      COMMON /BICKL1/  BI1(0:MKI1,3),PAS1,XLIM1,L1
      COMMON /BICKL2/  BI2(0:MKI2,3),PAS2,XLIM2,L2
      COMMON /BICKL3/  BI3(0:MKI3,3),PAS3,XLIM3,L3
      COMMON /BICKL4/  BI4(0:MKI4,3),PAS4,XLIM4,L4
      COMMON /BICKL5/  BI5(0:MKI5,3),PAS5,XLIM5,L5
      DOUBLE PRECISION ABSC(3,2)

      III(I,J)=(J+NSOUT)*NUNKNO+I+NSOUT+1
      IND(I,J) = MAX(I+NSOUT+1,J+NSOUT+1)*(MAX(I+NSOUT+1,J+NSOUT+1)-1)/2
     1 + MIN(I+NSOUT+1,J+NSOUT+1)
*----
* RECOVER EXCELL SPECIFIC TRACKING INFORMATION.
*             ALBEDO: SURFACE ALBEDOS (REAL(NBCDA))
*             KSPEC : KIND OF PIJ INTEGRATION (0:ISOTROPE,1:SPECULAR)
*             CUTOF : MFP CUTOFF FOR SPECULAR INTEGRATION
*----
      ISTATE(:NSTATE)=0
      CALL LCMGET(IPTRK,'STATE-VECTOR',ISTATE)
      KSPEC=ISTATE(10)
      CALL LCMGET(IPTRK,'EXCELTRACKOP',EXTKOP)
      CUTOF=EXTKOP(1)
      CALL LCMLEN(IPTRK,'ALBEDO',NBCDA,ITYLCM)
      ALLOCATE(ICODE(NBCDA),ALBG(NBCDA),ALBEDO(NBCDA))
      ICODE(:NBCDA)=0
      CALL LCMGET(IPTRK,'ICODE',ICODE)
      CALL LCMGET(IPTRK,'ALBEDO',ALBG)
*
      IPRT  = IPRNTP
      IF( IPRT.GE.ICPEND ) WRITE(IOUT,'(1X,A72//)') TITREC
      NPLANE = 1
      IF(IFTRAK .NE. 0) THEN
        READ(IFTRAK) CTRKT,NCOMNT,NTRK,IFMT
        IF( CTRKT .NE.'$TRK' .OR.
     >      NCOMNT.LT.0      .OR.
     >      NTRK  .EQ.0          ) CALL XABORT(NAMSBR//
     >      ': Invalid tracking file')
        DO IC= 1,NCOMNT
           READ(IFTRAK) COMENT
        ENDDO
        READ(IFTRAK) NDIM,ISPEC,NREG2,NSOUT,NALBG,NCOR,NANGL,MXSUB,MXSEG
        IF(NREG2.NE.NREG )THEN
           CALL XABORT(NAMSBR//': TRACKING FILE HAS INVALID # OF ZONES')
        ENDIF
        NCORT=NCOR
      ELSE
        IF(ISTATE(7) .NE. 4) CALL XABORT(NAMSBR//
     >  ': Tracking file required unless NXT: tracking provided')
        NREG2=ISTATE(1)
        IF(NREG2.NE.NREG )THEN
           CALL XABORT(NAMSBR//': STATE VECTOR HAS INVALID # OF ZONES')
        ENDIF
        NSOUT=ISTATE(5)
        ISPEC=ISTATE(9)
        NPOINT=ISTATE(17)
        MXSEG=ISTATE(18)
        NANGL=ISTATE(20)
        NPLANE=ISTATE(25)
        CALL LCMSIX(IPTRK,'NXTRecords  ',ILCMUP)
        CALL LCMGET(IPTRK,'G00000001DIM',IEDIMG)
        NDIM=IEDIMG(1)
        ITYPBC=IEDIMG( 2)
        NBUCEL=IEDIMG( 5)
        NUCELL(1)=IEDIMG(13)
        NUCELL(2)=IEDIMG(14)
        NUCELL(3)=IEDIMG(15)
        MXMSH=IEDIMG(16)
        MAXPIN=IEDIMG(19)
        MXGSUR=IEDIMG(24)
        MXGREG=IEDIMG(25)
        NCOR=1
        NCORT=NCOR
        NTRK=NANGL*NPLANE*NPOINT**(NDIM-1)
        IF(MXSEG .LE. 1) THEN
          IF(ISPEC .EQ. 0) THEN
            MXSEG=NBUCEL*
     >           ((MAXPIN+1)*(2*MXGREG+2)+MXGSUR+16)
          ELSE
            MXSEG=8*NANGL*NBUCEL*
     >            ((MAXPIN+1)*(2*MXGREG+2)+MXGSUR+16)
          ENDIF
        ENDIF
        MAXMSH=MAX(MXMSH,IEDIMG(17),IEDIMG(20))
      ENDIF
      NUNKNO= NREG+NSOUT+1
      IF(IFTRAK .NE. 0) THEN
        READ(IFTRAK) (VOLSUR(JUN,1),JUN=-NSOUT,NREG)
        READ(IFTRAK) (MATALB(JUN),JUN=-NSOUT,NREG)
        READ(IFTRAK) ( NSCRP,JUN=1,NALBG)
        READ(IFTRAK) ( ASCRP,JUN=1,NALBG)
        READ(IFTRAK) DANG0,(DASCRP,IUN=2,NDIM),
     >               ((DASCRP,IUN=1,NDIM),JUN=2,NANGL)
        READ(IFTRAK) (DASCRP,JUN=1,NANGL)
      ELSE
        CALL LCMGET(IPTRK,'MATALB      ',MATALB)
        ALLOCATE(DSV(-NSOUT:NREG))
        CALL LCMGET(IPTRK,'SAreaRvolume',DSV)
        DO JJ=-NSOUT,0
          VOLSUR(JJ,1)=0.25*REAL(DSV(JJ))
        ENDDO
        DO JJ=1,NREG
          VOLSUR(JJ,1)=REAL(DSV(JJ))
        ENDDO
*----
*  Allocate memory for NXT tracking
*----
        ALLOCATE(DGMESH(-1:MAXMSH,4))
        CALL NXTXYZ(IPTRK,IPRNTP,NDIM,ITYPBC,MAXMSH,NUCELL,ABSC,DGMESH)
        ALLOCATE(IUNFLD(2,NBUCEL))
        ALLOCATE(DANGLT(NDIM,NANGL),DORITR(NDIM*(NDIM+1),NPLANE,NANGL),
     >  DWGTRK(NANGL),DVNOR(NREG))
        NAMREC='G00000001CUF'
        CALL LCMGET(IPTRK,NAMREC,IUNFLD)
        CALL LCMGET(IPTRK,'TrackingDirc',DANGLT)
        CALL LCMGET(IPTRK,'TrackingOrig',DORITR)
        CALL LCMGET(IPTRK,'TrackingWgtD',DWGTRK)
        CALL LCMGET(IPTRK,'VTNormalize ',DVNOR)
      ENDIF
      DO ISBG=2,NSBG
        DO IUN= -NSOUT, NREG
          VOLSUR(IUN,ISBG)=VOLSUR(IUN,1)
        ENDDO
      ENDDO
*----
*  PREPARE FOR MULTIGROUP CALCULATION
*----
      ALLOCATE(SIGTAL(-NSOUT:NREG,NSBG),SIGT00(-NSOUT:NREG,NSBG))
      LPIJK= IPIJK.EQ.4
      SWNZBC= .FALSE.
      SWVOID= .FALSE.
      DO ISBG=1,NSBG
        IF(NPSYS(ISBG).NE.0) THEN
          ALBEDO(:NBCDA)=ALBG(:NBCDA)
          IF(NALBP .GT. 0) THEN
            DO ISA=1,NBCDA
              IF(ICODE(ISA).GT.0) ALBEDO(ISA)=ALBP(ICODE(ISA),ISBG)
            ENDDO
          ENDIF
          DO IUN= -NSOUT, -1
            SIGT00(IUN,ISBG)= 0.0
            IF(-MATALB(IUN).GT.NBCDA) CALL XABORT('EXCELP: NBCDA OV'
     >      //'ERFLOW(2).')
            SIGTAL(IUN,ISBG)= ALBEDO(-MATALB(IUN))
            SWNZBC= SWNZBC.OR.(SIGTAL(IUN,ISBG).NE.0.0)
          ENDDO
          IUN=0
          SIGT00(IUN,ISBG)= 0.0
          SIGTAL(IUN,ISBG)= 0.0
          DO IUN= 1, NREG
            SIGT00(IUN,ISBG)= XSSIGT(MATCOD(IUN),ISBG)
            SIGTAL(IUN,ISBG)= XSSIGT(MATCOD(IUN),ISBG)
            IF( SIGTAL(IUN,ISBG) .EQ. 0.0 )THEN
              SWVOID= .TRUE.
            ELSE
              VOLSUR(IUN,ISBG)=VOLSUR(IUN,ISBG)*SIGTAL(IUN,ISBG)
            ENDIF
          ENDDO
        ENDIF
      ENDDO
*----
*  CHOOSE ISOTROPIC OR SPECULAR B.C.
*----
      IF( KSPEC.EQ.0 )THEN
         SIGT => SIGT00
      ELSE
         SIGT => SIGTAL
      ENDIF
*
      NPROB = (NUNKNO*(NUNKNO+1))/2
      N2PRO = NUNKNO*NUNKNO
      IF(IPRNTP .GT. 1) THEN
         NALLOC=(2*N2PRO*NSBG)
         IF(LPIJK) NALLOC=NALLOC+(2*N2PRO*NSBG)
         WRITE(IOUT,6000) NALLOC/256
      ENDIF
      DPROB(:N2PRO,:NSBG)=0.0D0
      IF(LPIJK) DPROBX(:N2PRO,:NSBG)=0.0D0
      IF(IPRNTP.GT.1) WRITE(IOUT,6001)
*
      IF(IPRNTP .GE. 10) WRITE(IOUT,6010) MXSEG
*----
*  BATCH TRACKING STORAGE ALLOCATION
*----
      ALLOCATE(NCOIL1(NBATCH),NSUB(NBATCH),NBSEG(NBATCH),WEIGHT(NBATCH),
     1 NUMERO(MXSEG,NBATCH),LENGTH(MXSEG,NBATCH))
      IF( ISPEC.EQ.0 )THEN
*----
*  Standard tracking
*----
         IF( NDIM.EQ.2 )THEN
           ALLOCATE(LOPATH(MXSEG))
           DO IBATCH=1,(NTRK-1)/NBATCH+1
           DO ITRAK=(IBATCH-1)*NBATCH+1,MIN(IBATCH*NBATCH,NTRK)
             IL1=ITRAK-(IBATCH-1)*NBATCH
             IF(IFTRAK .NE. 0) THEN
*----
*  Read tracks from file
*----
               READ(IFTRAK)  NSUB(IL1),NBSEG(IL1),WEIGHT(IL1),IANG,
     >         (NUMERO(IL,IL1),IL=1,NBSEG(IL1)),
     >         (LENGTH(IL,IL1),IL=1,NBSEG(IL1))
               IF(NSUB(IL1).NE.1) CALL XABORT('EXCELP: NSUB.NE.1.')
               NCOIL1(IL1)=1
             ELSE
*----
*  Generate selected track
*----
               CALL NXTTGS(IPTRK ,IPRNTP,NDIM  ,NANGL ,NPOINT,NTRK  ,
     >                     ITRAK ,MAXMSH,NSOUT ,NREG  ,NUCELL,NBUCEL,
     >                     MXGSUR,MXGREG,MAXPIN,MXSEG ,ITYPBC,IUNFLD,
     >                     MATALB,DSV   ,DGMESH,DANGLT,DVNOR ,DWGTRK,
     >                     DORITR,NBSEG(IL1) ,NCORT ,WEIGHT(IL1),
     >                     NUMERO(1,IL1),LENGTH(1,IL1))
               NCOIL1(IL1)=NCORT
             ENDIF
           ENDDO
*$OMP  PARALLEL DO
*$OMP1 PRIVATE(IL1,ITRAK,LOPATH)
          DO ISBG=1,NSBG
             IF(NPSYS(ISBG).EQ.0) CYCLE
             DO ITRAK=(IBATCH-1)*NBATCH+1,MIN(IBATCH*NBATCH,NTRK)
               IL1=ITRAK-(IBATCH-1)*NBATCH
               IF(NCOIL1(IL1).EQ.0) CYCLE
               CALL PIJI2D(NREG,NSOUT,NBSEG(IL1),NCOR,SWVOID,
     >                     SIGT(-NSOUT,ISBG),WEIGHT(IL1),LENGTH(1,IL1),
     >                     NUMERO(1,IL1),LOPATH,DPROB(1,ISBG),
     >                     MKI1,BI1,PAS1,L1,
     >                     MKI2,BI2,PAS2,XLIM2,L2,
     >                     MKI3,BI3,PAS3,XLIM3)
             ENDDO ! ITRAK
           ENDDO ! ISBG
*$OMP END PARALLEL DO
           IF(LPIJK)THEN
*$OMP  PARALLEL DO
*$OMP1 PRIVATE(IL1,ITRAK,LOPATH)
             DO ISBG=1,NSBG
               IF(NPSYS(ISBG).EQ.0) CYCLE
               DO ITRAK=(IBATCH-1)*NBATCH+1,MIN(IBATCH*NBATCH,NTRK)
                 IL1=ITRAK-(IBATCH-1)*NBATCH
                 IF(NCOIL1(IL1).EQ.0) CYCLE
                 CALL PIJI2D(NREG,NSOUT,NBSEG(IL1),NCOR,SWVOID,
     >                       SIGT(-NSOUT,ISBG),WEIGHT(IL1),
     >                       LENGTH(1,IL1),NUMERO(1,IL1),LOPATH,
     >                       DPROBX(1,ISBG),MKI3,BI3,PAS3,L3,
     >                       MKI4,BI4,PAS4,XLIM4,L4,
     >                       MKI5,BI5,PAS5,XLIM5)
               ENDDO ! ITRAK
             ENDDO ! ISBG
*$OMP END PARALLEL DO
           ENDIF
           ENDDO ! IBATCH
           DEALLOCATE(LOPATH)
         ELSE
           ALLOCATE(STAYIN(MXSEG),GOSOUT(MXSEG))
           DO IBATCH=1,(NTRK-1)/NBATCH+1
           DO ITRAK=(IBATCH-1)*NBATCH+1,MIN(IBATCH*NBATCH,NTRK)
             IL1=ITRAK-(IBATCH-1)*NBATCH
             IF(IFTRAK .NE. 0) THEN
               READ(IFTRAK)  NSUB(IL1),NBSEG(IL1),WEIGHT(IL1),IANG,
     >         (NUMERO(IL,IL1),IL=1,NBSEG(IL1)),
     >         (LENGTH(IL,IL1),IL=1,NBSEG(IL1))
               IF(NSUB(IL1).NE.1) CALL XABORT('EXCELP: NSUB.NE.1.')
               NCOIL1(IL1)=1
             ELSE
               CALL NXTTGS(IPTRK ,IPRNTP,NDIM  ,NANGL ,NPOINT,NTRK  ,
     >                     ITRAK ,MAXMSH,NSOUT ,NREG  ,NUCELL,NBUCEL,
     >                     MXGSUR,MXGREG,MAXPIN,MXSEG ,ITYPBC,IUNFLD,
     >                     MATALB,DSV   ,DGMESH,DANGLT,DVNOR ,DWGTRK,
     >                     DORITR,NBSEG(IL1) ,NCORT ,WEIGHT(IL1),
     >                     NUMERO(1,IL1),LENGTH(1,IL1))
               NCOIL1(IL1)=NCORT
             ENDIF
           ENDDO
*$OMP  PARALLEL DO
*$OMP1 PRIVATE(IL1,ITRAK,STAYIN,GOSOUT)
           DO ISBG=1,NSBG
             IF(NPSYS(ISBG).EQ.0) CYCLE
             DO ITRAK=(IBATCH-1)*NBATCH+1,MIN(IBATCH*NBATCH,NTRK)
               IL1=ITRAK-(IBATCH-1)*NBATCH
               IF(NCOIL1(IL1).EQ.0) CYCLE
               CALL PIJI3D(NREG,NSOUT,NBSEG(IL1),NCOR,SWVOID,
     >         SIGT(-NSOUT,ISBG),WEIGHT(IL1),LENGTH(1,IL1),
     >         NUMERO(1,IL1),STAYIN,GOSOUT,DPROB(1,ISBG))
             ENDDO ! ITRAK
           ENDDO ! ISBG
*$OMP END PARALLEL DO
           ENDDO ! IBATCH
           DEALLOCATE(GOSOUT,STAYIN)
           IF(LPIJK) CALL XABORT(NAMSBR//': 3D PIJK NOT SUPPORTED')
         ENDIF
      ELSEIF( ISPEC.EQ.1 )THEN
*----
*  CYCLIC TRACKING
*----
         RCUTOF= CUTOF
         IF( NDIM.EQ.2 )THEN
            IF( DANG0.EQ. 0.0D0 )THEN
               NGSS= NANGL/8
            ELSE
               NGSS= (NANGL/4+1)/2
            ENDIF
            CALL ALGPT( NGSS,0.0,1.0,XGSS,WGSS)
            ALLOCATE(SIGANG(NGSS,-NSOUT:NREG,NSBG),STAYIN(NGSS*MXSEG),
     >      GOSOUT(NGSS*MXSEG))
            DO JGSS= 1, NGSS
              YGSS= SQRT(1.0 - XGSS(JGSS)**2)
              WGSS(JGSS)= WGSS(JGSS) * YGSS
              XGSS(JGSS)= 1.0/YGSS
              WGSSX(JGSS)= WGSS(JGSS) / (XGSS(JGSS)**2)
              DO ISBG=1,NSBG
                IF(NPSYS(ISBG).NE.0) THEN
                  DO IUN= -NSOUT,NREG
                    IF( MATALB(IUN).LE.0 )THEN
                      SIGANG(JGSS,IUN,ISBG)= SIGT(IUN,ISBG)
                    ELSE
                      SIGANG(JGSS,IUN,ISBG)= SIGT(IUN,ISBG)*XGSS(JGSS)
                    ENDIF
                  ENDDO
                ENDIF
              ENDDO
            ENDDO
*----
*  Loop over tracks
*  then loop over groups
*----
            DO IBATCH=1,(NTRK-1)/NBATCH+1
            DO ITRAK=(IBATCH-1)*NBATCH+1,MIN(IBATCH*NBATCH,NTRK)
              IL1=ITRAK-(IBATCH-1)*NBATCH
              IF(IFTRAK .NE. 0) THEN
                READ(IFTRAK)  NSUB(IL1),NBSEG(IL1),WEIGHT(IL1),
     >          (IANG,IL=1,NSUB(IL1)),(NUMERO(IL,IL1),IL=1,NBSEG(IL1)),
     >          (LENGTH(IL,IL1),IL= 1,NBSEG(IL1))
                NCOIL1(IL1)=1
              ELSE
               CALL NXTTGC(IPTRK ,IPRNTP,NDIM  ,NANGL ,NPOINT,NTRK  ,
     >                     ITRAK ,MAXMSH,NSOUT ,NREG  ,NUCELL,NBUCEL,
     >                     MXGSUR,MXGREG,MAXPIN,MXSEG ,ITYPBC,IUNFLD,
     >                     MATALB,DSV   ,DGMESH,DANGLT,DVNOR ,DWGTRK,
     >                     DORITR,NBSEG(IL1) ,NCORT ,WEIGHT(IL1),
     >                     NUMERO(1,IL1),LENGTH(1,IL1))
                NCOIL1(IL1)=NCORT
              ENDIF
            ENDDO
*$OMP  PARALLEL DO
*$OMP1 PRIVATE(IL1,ITRAK,STAYIN,GOSOUT)
            DO ISBG=1,NSBG
              IF(NPSYS(ISBG).EQ.0) CYCLE
              DO ITRAK=(IBATCH-1)*NBATCH+1,MIN(IBATCH*NBATCH,NTRK)
                IL1=ITRAK-(IBATCH-1)*NBATCH
                IF(NCOIL1(IL1).EQ.0) CYCLE
                CALL PIJS2D(NREG,NSOUT,NBSEG(IL1),WEIGHT(IL1),RCUTOF,
     >          NGSS,SIGANG(1,-NSOUT,ISBG),XGSS,WGSS,LENGTH(1,IL1),
     >          NUMERO(1,IL1),STAYIN,GOSOUT,DPROB(1,ISBG))
              ENDDO ! ITRAK
            ENDDO ! ISBG
*$OMP END PARALLEL DO
            IF(LPIJK)THEN
*             X-DIRECTION  PROBABILITIES CALCULATIONS ( PX=PY )
*$OMP  PARALLEL DO
*$OMP1 PRIVATE(IL1,ITRAK,STAYIN,GOSOUT)
              DO ISBG=1,NSBG
                IF(NPSYS(ISBG).EQ.0) CYCLE
                DO ITRAK=(IBATCH-1)*NBATCH+1,MIN(IBATCH*NBATCH,NTRK)
                  IL1=ITRAK-(IBATCH-1)*NBATCH
                  IF(NCOIL1(IL1).EQ.0) CYCLE
                  CALL PIJS2D(NREG,NSOUT,NBSEG(IL1),WEIGHT(IL1),RCUTOF,
     >            NGSS,SIGANG(1,-NSOUT,ISBG),XGSS,WGSSX,LENGTH(1,IL1),
     >            NUMERO(1,IL1),STAYIN,GOSOUT,DPROBX(1,ISBG))
                ENDDO ! ITRAK
              ENDDO ! ISBG
*$OMP END PARALLEL DO
            ENDIF
            ENDDO ! IBATCH
            DEALLOCATE(GOSOUT,STAYIN,SIGANG)
         ELSE
            ALLOCATE(STAYIN(MXSEG),GOSOUT(MXSEG))
            DO IBATCH=1,(NTRK-1)/NBATCH+1
            DO ITRAK=(IBATCH-1)*NBATCH+1,MIN(IBATCH*NBATCH,NTRK)
              IL1=ITRAK-(IBATCH-1)*NBATCH
              IF(IFTRAK .NE. 0) THEN
                READ(IFTRAK)  NSUB(IL1),NBSEG(IL1),WEIGHT(IL1),
     >          (IANG,IL=1,NSUB(IL1)),(NUMERO(IL,IL1),IL=1,NBSEG(IL1)),
     >          (LENGTH(IL,IL1),IL=1,NBSEG(IL1))
                NCOIL1(IL1)=1
              ELSE
               CALL NXTTGC(IPTRK ,IPRNTP,NDIM  ,NANGL ,NPOINT,NTRK  ,
     >                     ITRAK ,MAXMSH,NSOUT ,NREG  ,NUCELL,NBUCEL,
     >                     MXGSUR,MXGREG,MAXPIN,MXSEG ,ITYPBC,IUNFLD,
     >                     MATALB,DSV   ,DGMESH,DANGLT,DVNOR ,DWGTRK,
     >                     DORITR,NBSEG(IL1) ,NCORT ,WEIGHT(IL1),
     >                     NUMERO(1,IL1),LENGTH(1,IL1))
                NCOIL1(IL1)=NCORT
              ENDIF
            ENDDO
*$OMP  PARALLEL DO
*$OMP1 PRIVATE(IL1,ITRAK,STAYIN,GOSOUT)
            DO ISBG=1,NSBG
              IF(NPSYS(ISBG).EQ.0) CYCLE
              DO ITRAK=(IBATCH-1)*NBATCH+1,MIN(IBATCH*NBATCH,NTRK)
                IL1=ITRAK-(IBATCH-1)*NBATCH
                IF(NCOIL1(IL1).EQ.0) CYCLE
                CALL PIJS3D(NREG,NSOUT,NBSEG(IL1),WEIGHT(IL1),RCUTOF,
     >          SIGT(-NSOUT,ISBG),LENGTH(1,IL1),NUMERO(1,IL1),STAYIN,
     >          GOSOUT,DPROBX(1,ISBG))
              ENDDO ! ITRAK
            ENDDO ! ISBG
*$OMP END PARALLEL DO
            ENDDO ! IBATCH
            DEALLOCATE(GOSOUT,STAYIN)
         ENDIF
      ENDIF
      IF(IFTRAK .EQ. 0) THEN
        DEALLOCATE(DVNOR,DWGTRK,DORITR,DANGLT)
        DEALLOCATE(IUNFLD)
        DEALLOCATE(DGMESH,DSV)
        CALL LCMSIX(IPTRK,'NXTRecords  ',ILCMDN)
      ENDIF
*
      DO 2050 ISBG=1,NSBG
         IF(NPSYS(ISBG).EQ.0) GO TO 2050
         KSBG=(ISBG-1)*NUNKNO
         CALL PIJCMP(NREG,NSOUT,NCOR,DPROB(1,ISBG),
     >               VOLSUR(-NSOUT,ISBG),.FALSE.,DPROB(1,ISBG))
         IF(LPIJK)THEN
           CALL PIJCMP(NREG,NSOUT,NCOR,DPROBX(1,ISBG),
     >                 VOLSUR(-NSOUT,ISBG),.TRUE.,DPROBX(1,ISBG))
         ENDIF
 2050 CONTINUE
*----
*  BATCH TRACKING STORAGE DEALLOCATION
*----
      DEALLOCATE(LENGTH,NUMERO,WEIGHT,NBSEG,NSUB,NCOIL1)
*----
*  RENORMALIZE ALL ISOTROPIC PROBS WITH VARIOUS OPTIONS
*----
      DO 2060 ISBG=1,NSBG
      IF(NPSYS(ISBG).EQ.0) GO TO 2060
      IF( KSPEC.EQ.0 )THEN
         IF( NRENOR.EQ.1 )THEN
*
*           NORMALIZATION USING GELBARD SCHEME
            CALL PIJRGL(IPRT,NREG,NSOUT,SIGTAL(-NSOUT,ISBG),
     >                  DPROB(1,ISBG))
            IF(LPIJK) CALL PIJRGL(IPRT,NREG,NSOUT,SIGTAL(-NSOUT,ISBG),
     >                            DPROBX(1,ISBG))
         ELSEIF( NRENOR.EQ.2 )THEN
*
*           NORMALIZATION WORKING ON DIAGONAL COEFFICIENTS
            CALL PIJRDG(IPRT,NREG,NSOUT,SIGTAL(-NSOUT,ISBG),
     >                  DPROB(1,ISBG))
            IF(LPIJK) CALL PIJRDG(IPRT,NREG,NSOUT,SIGTAL(-NSOUT,ISBG),
     >                            DPROBX(1,ISBG))
         ELSEIF( NRENOR.EQ.3 )THEN
*
*           NORMALIZATION WORKING ON WEIGHT FACTORS TO KEEP DIAG = 0.0
            CALL PIJRNL(IPRT,NREG,NSOUT,SIGTAL(-NSOUT,ISBG),
     >                  DPROB(1,ISBG))
            IF(LPIJK) CALL PIJRNL(IPRT,NREG,NSOUT,SIGTAL(-NSOUT,ISBG),
     >                            DPROBX(1,ISBG))
         ELSEIF( NRENOR .EQ. 4 )THEN  ! ATTENTION
*
*           NORMALIZATION WORKING ON WEIGHT FACTORS ADDITIVE (HELIOS)
            CALL PIJRHL(IPRT,NREG,NSOUT,SIGTAL(-NSOUT,ISBG),
     >                  DPROB(1,ISBG))
            IF(LPIJK) CALL PIJRHL(IPRT,NREG,NSOUT,SIGTAL(-NSOUT,ISBG),
     >                            DPROBX(1,ISBG))
         ENDIF
         IF( IPRT.GE.ICPALL )THEN
            LOPT= -1
            MSYM=1
            WRITE(IOUT,'(1H )')
            WRITE(IOUT,'(35H   COLLISION PROBABILITIES OUTPUT: ,
     >                   35H *BEFORE* ALBEDO REDUCTION          )')
            CALL PIJWPR(LOPT,NREG,NSOUT,SIGTAL(-NSOUT,ISBG),
     >                  DPROB(1,ISBG),VOLSUR(1,ISBG),MSYM)
*
         ENDIF
      ENDIF
 2060 CONTINUE
      DEALLOCATE(ALBEDO,ALBG,ICODE)
      RETURN
*
 6010 FORMAT(' Maximum length of a line =',I10)
 6000 FORMAT(' *** SPACE REQUIRED FOR CP MATRICES = ',I10,' K ***')
 6001 FORMAT(' *** CP MATRICES ALLOCATED            ',10X,'   ***')
      END
