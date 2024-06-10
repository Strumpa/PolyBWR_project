*DECK MOCFCF
      SUBROUTINE MOCFCF(SUBFFI,SUBFFA,SUBLDC,SUBSCH,IFTRAK,NBTR,MXSUB,
     1                  MXSEG,NDIM,KPN,NREG,NSOUT,NMAT,NALB,NGEFF,NPHI,
     2                  NGSS,NLF,NFUNL,NMOD,NLFX,NLIN,NFUNLX,KEYFLX,
     3                  MATALB,NCONV,SIGANG,CAZ1,CAZ2,XGSS,YGSS,WGSS,
     4                  SOUR,ISGNR,IDIR,NBATCH,PHIOUT)
*
*-----------------------------------------------------------------------
*
*Purpose:
* Flux integration upon the cyclic tracking.
*
*Copyright:
* Copyright (C) 2002 Ecole Polytechnique de Montreal
* This library is free software; you can redistribute it and/or
* modify it under the terms of the GNU Lesser General Public
* License as published by the Free Software Foundation; either
* version 2.1 of the License, or (at your option) any later version
*
*Author(s): R. Roy and R. Le Tellier
*
*Parameters: input
* SUBFFI  Isotropic flux integration subroutine.
* SUBFFA  Anisotropic flux integration subroutine.
* SUBLDC  flux integration subroutine with linear-discontinuous source.
* SUBSCH  Track coefficients calculation subroutine.
* IFTRAK  tracking file unit number.
* NGEFF   number of groups to process.
* NMAT    number of mixtures.
* NLF     number of Legendre orders for the flux.
* NREG    number of regions.
* KPN     number of unknowns per energy group including spherical
*         harmonic terms and fundamental currents.
* NGSS    number of polar angles.
* NSOUT   number of surfaces.
* MXSUB   maximun number of subtracks in a track.
* MXSEG   maximun number of segments in a track.
* NBTR    number of tracks.
* NPHI    number of angles in the plane.
* NFUNL   number of moments of the flux (in 2D: NFUNL=NLF*(NLF+1)/2).
* NMOD    first dimension of ISGNR.
* NLFX    scattering anisotropy used to compute spherical harmonics.
* NLIN    linear discontinuous flag (=1 SC/DD0; =3 LDC/DD1).
* NFUNLX  number of spherical harmonics components.
* ISGNR   array of spherical harmonics signs.
* KEYFLX  position of flux elements in PHIIN vector.
* MATALB  mixture and albedo indices.
* WGSS    polar weights.
* XGSS    polar angle cosines.
* YGSS    polar angle sines.
* NALB    number of albedos.
* SIGANG  arrays of total cross-sections and albedos.
* CAZ1    first cosines of the different tracking azimuthal angles.
* CAZ2    second cosines of the different tracking azimuthal angles.
* SOUR    total source vector components.
* NCONV   logical array of convergence status for each group (.TRUE. for
*         not converged).
* IDIR    direction of fundamental current for TIBERE with MoC 
*         (=0,1,2,3).
* NBATCH  number of tracks processed in each OpenMP core (default: =1).
*
*Parameters: output
* PHIOUT  vector containing the zonal flux moments.
*
*-----------------------------------------------------------------------
*
      IMPLICIT  NONE
*----
*  SUBROUTINE ARGUMENTS
*----
      INTEGER IFTRAK,NGEFF,NMAT,NLF,NREG,NDIM,KPN,NGSS,NSOUT,MXSUB,
     1 MXSEG,NBTR,NPHI,NFUNL,NMOD,NLFX,NLIN,NFUNLX,ISGNR(NMOD,NFUNLX),
     2 KEYFLX(NREG,NLIN,NFUNL),MATALB(-NSOUT:NREG),NALB,IDIR,NBATCH
      REAL WGSS(NGSS),XGSS(NGSS),YGSS(NGSS),
     1 SIGANG(-NALB:NMAT,NGEFF)
      DOUBLE PRECISION CAZ1(NPHI),CAZ2(NPHI),SOUR(KPN,NGEFF),
     1 PHIOUT(KPN,NGEFF)
      LOGICAL NCONV(NGEFF)
      EXTERNAL SUBFFI,SUBFFA,SUBLDC,SUBSCH
*----
*  LOCAL VARIABLES
*----
      INTEGER MXE
      PARAMETER (MXE=64)
      INTEGER ILINE,IANG,JANG,ISUB,I,IE,II,NOMI,NZI,IND,JF,INDX,
     1 INDY,IREG,I0,IL1,IBATCH,NDFUNLX
      DOUBLE PRECISION DWEIG(MXE),Q0,Q1,Q0X,Q1X,Q0Y,Q1Y
      LOGICAL LNEW
*----
*  ALLOCATABLE ARRAYS
*----
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NSUB,NSEG
      INTEGER, ALLOCATABLE, DIMENSION(:,:) :: NRSEG,KANGL
      REAL, ALLOCATABLE, DIMENSION(:,:) :: RHARM
      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: TRHAR
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: DSIG,WEIGHT,EXPT,
     1 EXP2
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: SEGLEN,COEFI,
     1 OMEGAX,OMEGAY,OMG2,FLUX,PHIU,FLM,FLP,CYM,CYP,DFLM,DFLP,CYM2,CYP2
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: PHIV,DPHIV,
     1 FLUV,DFLUV
*----
*  SCRATCH STORAGE ALLOCATION
*----
      ALLOCATE(NSUB(NBATCH),NSEG(NBATCH),WEIGHT(NBATCH),
     1 KANGL(MXSUB,NBATCH),NRSEG(MXSEG,NBATCH),SEGLEN(MXSEG,NBATCH))
      ALLOCATE(FLUX(KPN,NGEFF),PHIU(KPN,NGEFF),FLM(NGSS,MXSEG),
     1 FLP(NGSS,MXSEG),CYM(NGSS,MXSEG),CYP(NGSS,MXSEG),EXPT(NGSS*MXSEG))
      ALLOCATE(OMEGAX(NPHI,NGSS),OMEGAY(NPHI,NGSS),OMG2(NGSS,3))
*---
* Compute flux and currents for this tracking line
*---
      PHIOUT(:KPN,:NGEFF)=0.0D0
      PHIU(:KPN,:NGEFF)=0.0D0
      IF((NLF.EQ.1).AND.(NLIN.EQ.1)) THEN
*----
*  ISOTROPIC SCATTERING
*----
        ALLOCATE(EXP2(2*NGSS*MXSEG))
        DO IE=1,NGSS
          DO IANG=1,NPHI
            OMEGAX(IANG,IE)=CAZ1(IANG)/YGSS(IE)
            OMEGAY(IANG,IE)=CAZ2(IANG)/YGSS(IE)
          ENDDO   
        ENDDO
        DO IBATCH=1,(NBTR-1)/NBATCH+1
        DO ILINE=(IBATCH-1)*NBATCH+1,MIN(IBATCH*NBATCH,NBTR)
          IL1=ILINE-(IBATCH-1)*NBATCH
          READ(IFTRAK) NSUB(IL1),NSEG(IL1),WEIGHT(IL1),
     1    (KANGL(I,IL1),I=1,NSUB(IL1)),(NRSEG(I,IL1),I=1,NSEG(IL1)),
     2    (SEGLEN(I,IL1),I=1,NSEG(IL1))
          IF(NSUB(IL1).GT.MXSUB) CALL XABORT('MOCFCF: MXSUB OVERFLOW.')
        ENDDO
*$OMP  PARALLEL DO
*$OMP1 PRIVATE(IL1,IE,DWEIG,ISUB,LNEW,NOMI,NZI,IND,FLM,FLP,OMG2)
*$OMP2 PRIVATE(FLUX,I0,II,EXPT,EXP2,CYM,CYP)
*$OMP3 REDUCTION(+:PHIU)
        DO ILINE=(IBATCH-1)*NBATCH+1,MIN(IBATCH*NBATCH,NBTR)
          IL1=ILINE-(IBATCH-1)*NBATCH
          FLUX(:KPN,:NGEFF)=0.0D0
          OMG2(:NGSS,:3)=0.0D0
          DO IE=1,NGSS
            DWEIG(IE)=WEIGHT(IL1)*WGSS(IE)
          ENDDO
          DO II=1,NGEFF
            IF(NCONV(II)) THEN
              ISUB=0
              LNEW=.TRUE.
              DO I0=1,NSEG(IL1)
                NOMI=NRSEG(I0,IL1)
                NZI=MATALB(NOMI)
                IF(NZI.LE.0) THEN
                  LNEW=.TRUE.
                  DO IE=1,NGSS
                     FLM(IE,I0)=0.D0
                     FLP(IE,I0)=0.D0
                  ENDDO
                ELSE
                  IF(LNEW) THEN
                    ISUB=ISUB+1
                    LNEW=.FALSE.
                  ENDIF
                  IND=KEYFLX(NOMI,1,1)
                  IF(IDIR.EQ.0) THEN
                   DO IE=1,NGSS
                     FLM(IE,I0)=DWEIG(IE)*SOUR(IND,II)
                     FLP(IE,I0)=FLM(IE,I0)
                   ENDDO
                  ELSEIF(IDIR.EQ.1) THEN
                   DO IE=1,NGSS
                     OMG2(IE,1)=3.0D0*OMEGAX(KANGL(ISUB,IL1),IE)**2
                     FLM(IE,I0)=DWEIG(IE)*SOUR(IND,II)*OMG2(IE,1)
                     FLP(IE,I0)=FLM(IE,I0)
                   ENDDO
                  ELSEIF(IDIR.EQ.2) THEN
                   DO IE=1,NGSS
                     OMG2(IE,2)=3.0D0*OMEGAY(KANGL(ISUB,IL1),IE)**2
                     FLM(IE,I0)=DWEIG(IE)*SOUR(IND,II)*OMG2(IE,2)
                     FLP(IE,I0)=FLM(IE,I0)
                   ENDDO
                  ELSEIF(IDIR.EQ.3) THEN
                   DO IE=1,NGSS
                     OMG2(IE,3)=3.0D0*(1.0-1.0/YGSS(IE)**2)
                     FLM(IE,I0)=DWEIG(IE)*SOUR(IND,II)*OMG2(IE,3)
                     FLP(IE,I0)=FLM(IE,I0)
                   ENDDO
                  ENDIF
                ENDIF
              ENDDO
*               MOCFFIR: 'Source Term Isolation' Strategy turned on
*               MOCFFIS: 'Source Term Isolation' Strategy turned off
*               MOCFFIT: 'MOCC/MCI' Iterative Strategy
              CALL SUBFFI(SUBSCH,NREG,NSOUT,KPN,NMAT,NSEG(IL1),
     1             SEGLEN(1,IL1),NRSEG(1,IL1),NGSS,MATALB,
     2             SIGANG(-NALB,II),KEYFLX,YGSS,FLUX(1,II),EXPT,
     3             EXP2,FLM,FLP,CYM,CYP,IDIR,OMG2)
            ENDIF
          ENDDO
          PHIU(:KPN,:NGEFF)=PHIU(:KPN,:NGEFF)+FLUX(:KPN,:NGEFF)
        ENDDO ! ILINE
*$OMP END PARALLEL DO
        ENDDO ! IBATCH
        DEALLOCATE(OMG2,OMEGAY,OMEGAX)
        DEALLOCATE(EXP2)
      ELSE IF(NLIN.EQ.1) THEN
*----
*  ANISOTROPIC SCATTERING
*----
        ALLOCATE(EXP2(2*NGSS*MXSEG))
        ALLOCATE(RHARM(NGSS,NFUNL),TRHAR(NGSS,NFUNL,NPHI,2))
        DO IANG=1,NPHI
          CALL MOCCHR(2,NLF-1,NFUNL,NGSS,XGSS,CAZ1(IANG),CAZ2(IANG),
     1    RHARM)
          DO 15 JF=1,NFUNL
          DO 10 IE=1,NGSS
            TRHAR(IE,JF,IANG,1)=ISGNR(1,JF)*RHARM(IE,JF)
            TRHAR(IE,JF,IANG,2)=ISGNR(NMOD,JF)*RHARM(IE,JF)
 10       CONTINUE
 15       CONTINUE
        ENDDO
        DO IBATCH=1,(NBTR-1)/NBATCH+1
        DO ILINE=(IBATCH-1)*NBATCH+1,MIN(IBATCH*NBATCH,NBTR)
          IL1=ILINE-(IBATCH-1)*NBATCH
          READ(IFTRAK) NSUB(IL1),NSEG(IL1),WEIGHT(IL1),
     1    (KANGL(I,IL1),I=1,NSUB(IL1)),(NRSEG(I,IL1),I=1,NSEG(IL1)),
     1    (SEGLEN(I,IL1),I=1,NSEG(IL1))
          IF(NSUB(IL1).GT.MXSUB) CALL XABORT('MOCFCF: MXSUB OVERFLOW.')
        ENDDO
*$OMP  PARALLEL DO
*$OMP1 PRIVATE(IL1,IE,DWEIG,ISUB,LNEW,NOMI,NZI,Q0,Q1,FLM,FLP)
*$OMP2 PRIVATE(FLUX,I0,II,EXPT,EXP2,CYM,CYP)
*$OMP3 REDUCTION(+:PHIU)
        DO ILINE=(IBATCH-1)*NBATCH+1,MIN(IBATCH*NBATCH,NBTR)
          IL1=ILINE-(IBATCH-1)*NBATCH
          FLUX(:KPN,:NGEFF)=0.0D0
          DO IE=1,NGSS
            DWEIG(IE)=WEIGHT(IL1)*WGSS(IE)
          ENDDO
          DO II=1,NGEFF
            IF(NCONV(II)) THEN
              FLM(:NGSS,:MXSEG)=0.0D0
              FLP(:NGSS,:MXSEG)=0.0D0
              ISUB=0
              LNEW=.TRUE.
              DO I0=1,NSEG(IL1)
                NOMI=NRSEG(I0,IL1)
                NZI=MATALB(NOMI)
                IF(NZI.LE.0) THEN
                  LNEW=.TRUE.
                ELSE
                  IF(LNEW) THEN
                     ISUB=ISUB+1
                     LNEW=.FALSE.
                  ENDIF
                  DO IE=1,NGSS
                    Q0=0.D0
                    Q1=0.D0
                    DO JF=1,NFUNL
                      IND=KEYFLX(NOMI,1,JF)
                      Q0=Q0+SOUR(IND,II)*TRHAR(IE,JF,KANGL(ISUB,IL1),2)
                      Q1=Q1+SOUR(IND,II)*TRHAR(IE,JF,KANGL(ISUB,IL1),1)
                    ENDDO
                    FLM(IE,I0)=DWEIG(IE)*Q0
                    FLP(IE,I0)=DWEIG(IE)*Q1
                  ENDDO
                ENDIF
              ENDDO
              IF(ISUB.NE.NSUB(IL1)) CALL XABORT('MOCFCF: NSUB ERROR.')
*                MOCFFAR: 'Source Term Isolation' Strategy turned on
*                MOCFFAS: 'Source Term Isolation' Strategy turned off
*                MOCFFAT: 'MOCC/MCI' Iterative Strategy
              CALL SUBFFA(SUBSCH,NREG,NSOUT,KPN,NMAT,NSEG(IL1),
     1             SEGLEN(1,IL1),NRSEG(1,IL1),NGSS,NFUNL,MATALB,
     2             SIGANG(-NALB,II),KEYFLX,YGSS,FLUX(1,II),EXPT,EXP2,
     3             FLM,FLP,CYM,CYP,NPHI,NSUB(IL1),KANGL(1,IL1),TRHAR)
            ENDIF
          ENDDO
          PHIU(:KPN,:NGEFF)=PHIU(:KPN,:NGEFF)+FLUX(:KPN,:NGEFF)
        ENDDO ! ILINE
*$OMP END PARALLEL DO
        ENDDO ! IBATCH
        DEALLOCATE(EXP2)
        DEALLOCATE(TRHAR,RHARM)
      ELSE IF(NLIN.EQ.3) THEN
*----
*  LINEAR DISCONTINUOUS SOURCE APPROXIMATION
*----
        NDFUNLX=NDIM*NFUNLX
        ALLOCATE(PHIV(NFUNLX,NREG,NGEFF),DPHIV(NDFUNLX,NREG,NGEFF))
        ALLOCATE(FLUV(NFUNLX,NREG,NGEFF),DFLUV(NDFUNLX,NREG,NGEFF))
        ALLOCATE(EXP2(5*NGSS*MXSEG))
        ALLOCATE(RHARM(NGSS,NFUNLX),TRHAR(NGSS,NFUNLX,NPHI,2))
        DO IANG=1,NPHI
          CALL MOCCHR(2,NLFX-1,NFUNLX,NGSS,XGSS,CAZ1(IANG),CAZ2(IANG),
     1                RHARM)
          DO 25 JF=1,NFUNLX
          DO 20 IE=1,NGSS
            TRHAR(IE,JF,IANG,1)=ISGNR(1,JF)*RHARM(IE,JF)
            TRHAR(IE,JF,IANG,2)=ISGNR(NMOD,JF)*RHARM(IE,JF)
 20       CONTINUE
 25       CONTINUE
        ENDDO
        DO II=1,NGEFF            
          IF(NCONV(II)) THEN
            CALL XDDSET(PHIV(1,1,II),NFUNLX*NREG,0.0D0)
            CALL XDDSET(DPHIV(1,1,II),2*NFUNLX*NREG,0.0D0)
            PHIV(:NFUNLX,:NREG,II)=0.0D0
            DPHIV(:NDFUNLX,:NREG,II)=0.0D0
          ENDIF
        ENDDO
        ALLOCATE(DFLM(NGSS,MXSEG),DFLP(NGSS,MXSEG),DSIG(MXSEG),
     1           CYM2(NGSS,MXSEG),CYP2(NGSS,MXSEG))
        DO IBATCH=1,(NBTR-1)/NBATCH+1
        DO ILINE=(IBATCH-1)*NBATCH+1,MIN(IBATCH*NBATCH,NBTR)
          IL1=ILINE-(IBATCH-1)*NBATCH
          READ(IFTRAK) NSUB(IL1),NSEG(IL1),WEIGHT(IL1),
     1    (KANGL(I,IL1),I=1,NSUB(IL1)),(NRSEG(I,IL1),I=1,NSEG(IL1)),
     1    (SEGLEN(I,IL1),I=1,NSEG(IL1))
          IF(NSUB(IL1).GT.MXSUB) CALL XABORT('MOCFCF: MXSUB OVERFLOW.')
        ENDDO
*$OMP  PARALLEL DO
*$OMP1 PRIVATE(IL1,IE,DWEIG,ISUB,JANG,LNEW,NOMI,NZI,Q0,Q1,Q0X,Q1X,Q0Y)
*$OMP2 PRIVATE(Q1Y,IND,INDX,INDY,FLM,FLP,DFLM,DFLP,FLUX,FLUV,DFLUV,I0)
*$OMP3 PRIVATE(II,DSIG,EXPT,EXP2,CYM,CYP,CYM2,CYP2)
*$OMP4 REDUCTION(+:PHIU,PHIV,DPHIV)
        DO ILINE=(IBATCH-1)*NBATCH+1,MIN(IBATCH*NBATCH,NBTR)
          IL1=ILINE-(IBATCH-1)*NBATCH
          FLUX(:KPN,:NGEFF)=0.0D0
          FLUV(:NFUNLX,:NREG,:NGEFF)=0.0D0
          DFLUV(:NDFUNLX,:NREG,:NGEFF)=0.0D0
          DO IE=1,NGSS
            DWEIG(IE)=WEIGHT(IL1)*WGSS(IE)
          ENDDO
          DO II=1,NGEFF
            IF(NCONV(II)) THEN
              FLM(:NGSS,:MXSEG)=0.0D0
              FLP(:NGSS,:MXSEG)=0.0D0
              DFLM(:NGSS,:MXSEG)=0.0D0
              DFLP(:NGSS,:MXSEG)=0.0D0
              ISUB=0
              JANG=0
              LNEW=.TRUE.
              DO I0=1,NSEG(IL1)
                NOMI=NRSEG(I0,IL1)
                NZI=MATALB(NOMI)
                IF(NZI.LE.0) THEN
                  LNEW=.TRUE.
                ELSE
                  IF(LNEW) THEN
                     ISUB=ISUB+1
                     JANG=KANGL(ISUB,IL1)
                     LNEW=.FALSE.
                  ENDIF
                  DO IE=1,NGSS
                     Q0=0.D0
                     Q1=0.D0
                     Q0X=0.0D0
                     Q1X=0.0D0
                     Q0Y=0.0D0
                     Q1Y=0.0D0
                     DO JF=1,NFUNL
                        IND=KEYFLX(NOMI,1,JF)
                        INDX=KEYFLX(NOMI,2,JF)         
                        INDY=KEYFLX(NOMI,3,JF)         
                        Q0=Q0+SOUR(IND,II)*TRHAR(IE,JF,JANG,2)
                        Q1=Q1+SOUR(IND,II)*TRHAR(IE,JF,JANG,1)
                        Q0X=Q0X+SOUR(INDX,II)*TRHAR(IE,JF,JANG,2)
                        Q1X=Q1X+SOUR(INDX,II)*TRHAR(IE,JF,JANG,1)
                        Q0Y=Q0Y+SOUR(INDY,II)*TRHAR(IE,JF,JANG,2)
                        Q1Y=Q1Y+SOUR(INDY,II)*TRHAR(IE,JF,JANG,1)
                     ENDDO
                     FLM(IE,I0)=Q0
                     FLP(IE,I0)=Q1
                     DFLM(IE,I0)=-Q0X*CAZ1(JANG)-Q0Y*CAZ2(JANG)
                     DFLP(IE,I0)=Q1X*CAZ1(JANG)+Q1Y*CAZ2(JANG)
                  ENDDO
                ENDIF
              ENDDO
              IF(ISUB.NE.NSUB(IL1)) CALL XABORT('MOCFCF: NSUB ERROR.')
*                MOCFFAL: 'Source Term Isolation' Strategy turned off
              CALL SUBLDC(SUBSCH,NREG,NSOUT,NMAT,NSEG(IL1),
     1             SEGLEN(1,IL1),NRSEG(1,IL1),NGSS,NFUNLX,MATALB,
     2             DWEIG,SIGANG(-NALB,II),YGSS,FLM,FLP,DFLM,DFLP,
     3             NPHI,NSUB(IL1),KANGL(1,IL1),TRHAR,FLUV(1,1,II),
     4             DFLUV(1,1,II),DSIG,EXPT,EXP2,CYM,CYP,CYM2,CYP2)
            ENDIF
          ENDDO
          PHIU(:KPN,:NGEFF)=PHIU(:KPN,:NGEFF)+FLUX(:KPN,:NGEFF)
          PHIV(:NFUNLX,:NREG,:NGEFF)=PHIV(:NFUNLX,:NREG,:NGEFF)+
     1    FLUV(:NFUNLX,:NREG,:NGEFF)
          DPHIV(:NDFUNLX,:NREG,:NGEFF)=DPHIV(:NDFUNLX,:NREG,:NGEFF)+
     1    DFLUV(:NDFUNLX,:NREG,:NGEFF)
        ENDDO ! ILINE
*$OMP   END PARALLEL DO
        ENDDO ! IBATCH
        DEALLOCATE(CYP2,CYM2,DSIG,DFLP,DFLM)
        DEALLOCATE(EXP2)
        DEALLOCATE(TRHAR,RHARM)
        ALLOCATE(COEFI(2*NFUNLX,2*NFUNLX))
        CALL MCGCOEF(NFUNLX,NGSS,YGSS,WGSS,NPHI,CAZ1,CAZ2,COEFI)
        DO II=1,NGEFF            
        IF(NCONV(II)) THEN
          DO IREG=1,NREG
            DPHIV(:,IREG,II)=MATMUL(COEFI,DPHIV(:,IREG,II))
            DO JF=1,NFUNL
              PHIU(KEYFLX(IREG,1,JF),II)=PHIV(JF,IREG,II)
              PHIU(KEYFLX(IREG,2,JF),II)=DPHIV(JF,IREG,II)
              PHIU(KEYFLX(IREG,3,JF),II)=DPHIV(NFUNLX+JF,IREG,II)
            ENDDO
          ENDDO
        ENDIF
        ENDDO
        DEALLOCATE(COEFI,DFLUV,FLUV,DPHIV,PHIV)
      ENDIF
      PHIOUT(:KPN,:NGEFF)=PHIU(:KPN,:NGEFF)
*----
*  SCRATCH STORAGE DEALLOCATION
*----
      DEALLOCATE(EXPT,CYP,CYM,FLP,FLM,PHIU,FLUX)
      DEALLOCATE(SEGLEN,NRSEG,KANGL,WEIGHT,NSEG,NSUB)
*
      RETURN
      END
