!
!---------------------------------------------------------------------
!
!Purpose:
! Support subroutines for isotropic and specular boundary conditions cases
!
!Copyright:
! Copyright (C) 2001 Ecole Polytechnique de Montreal.
!
!Author(s):
! X. Warin
!
!---------------------------------------------------------------------
!
MODULE SAL_AUX_MOD

  USE PRECISION_AND_KINDS, ONLY : PDB,SMALL,PI,TWOPI,HALFPI,INFINITY
  USE SAL_NUMERIC_MOD,     ONLY : SAL141

CONTAINS
  SUBROUTINE SAL231(RTRACK,ITRACK,DELX,EX0,EY0,ANGLE)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! print out trajectory information
    !
    !Parameters: input
    ! RTRACK    floating point vectors to store trajectory information
    ! ITRACK    integer vectors to store trajectory information
    ! DELX      initial point of trajectory (D=0)
    ! EX0       first direction cosine
    ! EY0       second direction cosine
    ! ANGLE     track angle
    !
    !---------------------------------------------------------------------
    !
    USE PRECISION_AND_KINDS, ONLY : PDB
    USE SAL_TRACKING_TYPES,  ONLY : NNN,NMAX2
    !**
    IMPLICIT NONE
    INTEGER,     INTENT(IN), DIMENSION(:) :: ITRACK
    REAL(PDB),   INTENT(IN), DIMENSION(:) :: RTRACK
    REAL(PDB),   INTENT(IN)               :: DELX,EX0,EY0
    REAL(PDB),   INTENT(IN)               :: ANGLE
    !**
    INTEGER :: I,KM,K,II,JSURF,JPHI,JPSI,NTRACK,NBTOT,I0
    REAL(PDB) :: ANG0
    REAL(PDB), PARAMETER :: SMALLT=1.E-10
    INTEGER, PARAMETER :: FOUT =6
    !**
    ANG0=ANGLE
    NTRACK=ITRACK(1)
    NBTOT=ITRACK(2)
    WRITE(FOUT,'(//,3X,"TRAJECTORY",/,3X,"==========", &
         &      //3X,"DELX = ",1P,E12.4,6X,"EX EY = ",1P,2E12.4, &
         &      /,3X,"WITH SMALL = ",1P,E12.4,5X,"(",1P,E12.4," DEGREES )", &
         &      /,3X,"ANGLE # AND WEIGHT ",I6,1P,E12.4, &
         &      /,3X,"NBER OF SUB-TRAJ     = ",I6, &
         &      /,3X,"NBER OF ELEM IN TRAJ = ",I6,/)') &
         DELX,EX0,EY0,SMALLT,ANG0,ITRACK(7),RTRACK(7),NBTOT,NTRACK
    WRITE(FOUT,'(/,20X,"SURF",3X,"PHI",3X,"PSI",5X,"SINPHI",7X, &
         &      "COSPHI",/,20X,4("-"),2(3X,3("-")),5X,6("-"),7X,6("-"))')
    IF(ITRACK(5)/=0)THEN
       JSURF=1
       JPHI=ITRACK(5)
       JPSI=JPHI
       WRITE(FOUT,'(3X,A14,3I6,2X,1P,2E13.4)')'LEFT   SURFACE',JSURF,JPHI,JPSI, &
            RTRACK(3),RTRACK(5)
    ENDIF
    IF(ITRACK(6)/=0)THEN
       JSURF=1
       JPHI=ITRACK(6)
       JPSI=JPHI
       WRITE(FOUT,'(3X,A14,3I6,2X,1P,2E13.4)')'RIGHT  SURFACE',JSURF,JPHI,JPSI, &
            RTRACK(4),RTRACK(6)
    ENDIF
    IF(NTRACK/=0)THEN
       I0=NTRACK+NNN
       !        print sub-trajectories information
       WRITE(FOUT,'(/,"   SUB-TRAJECTORIES:")')
       WRITE(FOUT,'(/,"   NBER",2X,"NBER OF ELEM",2X," ANGLE", &
            &/,"   ----",2X,"------------",2X," -----")')
       DO I=1,NBTOT
          WRITE(FOUT,'(I6,"*",4X,I6,5X,I6)')I,ITRACK(I0+2*I-1),ITRACK(I0+2*I)
       ENDDO
       !        print trajectory
       WRITE(FOUT,'(/,"   TRAJECTORY:",/,"   NBER",3(5X,"REG",2X," LENGTH",3X,"ELEM"),/, &
       & 3X,"----",3(5X,"---",3X,"------",3X,"----"),1X,/)')
       DO I=1,NTRACK,3
          II=I+NNN
          KM=MIN(I+2,NTRACK)+NNN
          WRITE(FOUT,'(1P,I6,"*",3(I7,E10.2,I7))')I, &
               (ITRACK(K),RTRACK(K),ITRACK(K+NMAX2),K=II,KM)
       ENDDO
    ELSE
       WRITE(FOUT,'(1X,"==> Track without intersections")')
    ENDIF
    WRITE(FOUT,'(/)')
    !
  END SUBROUTINE SAL231
  !
  SUBROUTINE SAL232(ITRACK,RTRACK,FACNRM,GG,SURFN,CURRN)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! computes numerical volumes: uses local macro arrays
    !
    !Parameters: input
    ! ITRACK    integer vectors to store trajectory information
    ! RTRACK    floating point vectors to store trajectory information
    !
    !Parameters: input/output
    ! FACNRM    numerical volumes per direction
    ! SURFN     numerical areas
    ! CURRN     numerical currents
    ! GG        geometry basic information.
    !
    !---------------------------------------------------------------------
    !
    USE SAL_GEOMETRY_TYPES, ONLY : T_G_BASIC
    USE SAL_TRACKING_TYPES, ONLY : NNN
    !***
    IMPLICIT NONE
    INTEGER,   INTENT(IN),    DIMENSION(:)   :: ITRACK
    REAL(PDB), INTENT(IN),    DIMENSION(:)   :: RTRACK
    REAL(PDB), INTENT(INOUT), DIMENSION(:), OPTIONAL   :: SURFN,CURRN
    REAL(PDB), INTENT(INOUT), DIMENSION(:,:) :: FACNRM
    TYPE(T_G_BASIC) :: GG
    !     DIMENSION        ITRACK(*),RTRACK(*),FACNRM(NBREG,NPHI),
    !                      SURFN(NCURR,2),CURRN(NCURR,2)
    !***
    INTEGER  :: LASTI,II,I,ICURR,NPHI,IPHI,C,P
    REAL     :: WT,WR
    !***
    !     total weight and space weight
    NPHI=SIZE(FACNRM,2)
    WT=REAL(RTRACK(7)) ; WR=REAL(RTRACK(8))
    LASTI=ITRACK(1)
    C=0; P=LASTI+NNN; IPHI=0
    DO II=1+NNN,LASTI+NNN
       IF((II-NNN)>C) THEN
          C=C+ITRACK(P+1)
          IPHI=ITRACK(P+2)
          IF(IPHI>NPHI) IPHI=IPHI-NPHI
          P=P+2
       ENDIF
       I=ITRACK(II)
       IF(IPHI==0) CALL XABORT('SAL232: invalid IPHI')
       FACNRM(I,IPHI)=FACNRM(I,IPHI)+RTRACK(II)*WR
    ENDDO
    IF(GG%NB_SURF2/=0)THEN
       ICURR=ITRACK(5)
       IF(ICURR>0) THEN
          ! left surface: convert into 2d horizontal currents:
          SURFN(ICURR)=SURFN(ICURR)+WT/RTRACK(5)
          CURRN(ICURR)=CURRN(ICURR)+WT
       ENDIF
       ICURR=ITRACK(6)
       IF(ICURR>0) THEN
          ! right surface: convert into 2d horizontal currents:
          SURFN(ICURR)=SURFN(ICURR)+WT/RTRACK(6)
          CURRN(ICURR)=CURRN(ICURR)+WT
       ENDIF
    ENDIF
    !
  END SUBROUTINE SAL232
  !
  SUBROUTINE SAL235(NPIECE,THETA0,EX0,EY0,IPAR,RPAR,PEREXT,NPERIM)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! computes total perimeter projection on line orthogonal to
    ! trajectory and with origin the center of coordinates
    !
    !Parameters: input
    ! THETA0    THETA- and THETA+ for this trajectory
    ! EX0       first direction cosine
    ! EY0       second direction cosine
    ! IPAR      integer element data
    ! RPAR      floating point element data
    ! PEREXT    macro perimeter
    ! NPERIM    number of elements in perimeter
    !
    !Parameters: output
    ! NPIECE    number of pieces
    !
    !---------------------------------------------------------------------
    !
    USE SAL_GEOMETRY_TYPES,  ONLY : NIPAR,NRPAR,G_ELE_TYPE
    USE SAL_TRACKING_TYPES,  ONLY : DPIECE
    !***
    IMPLICIT NONE
    REAL(PDB), INTENT(IN)                 :: EX0,EY0
    REAL(PDB), INTENT(IN), DIMENSION(:)   :: THETA0
    INTEGER,   INTENT(OUT)                :: NPIECE
    INTEGER,   INTENT(IN), DIMENSION(:,:) :: IPAR
    REAL(PDB), INTENT(IN), DIMENSION(:,:) :: RPAR
    INTEGER,   INTENT(IN), DIMENSION(:)   :: PEREXT
    INTEGER,   INTENT(IN)                 :: NPERIM
    !***
    REAL(PDB) :: X,Y,RAD,DCENT,THETA1,THETA2,THETAM,THETA,DAUX,DMIN,DMAX
    INTEGER   :: L,ELEM,TYPE,IEND,ISIDE
    LOGICAL   :: LGONE
    REAL(PDB), PARAMETER, DIMENSION(2) :: SIGNV = (/-1._PDB,1._PDB/)
    INTEGER, PARAMETER :: FOUT =6
    !***
    LGONE=.TRUE.
    DMIN=0._PDB;  DMAX=0._PDB;
    DO L=1,NPERIM
       ELEM=PEREXT(L)
       ! treat element
       TYPE=IPAR(1,ELEM)
       IF(TYPE==G_ELE_TYPE(2))THEN
          ! circle:
          RAD=RPAR(3,ELEM)
          DCENT=RPAR(1,ELEM)*EY0-RPAR(2,ELEM)*EX0
          ! project tangents to circle
          IF(LGONE)THEN
             DMIN=DCENT-RAD
             DMAX=DCENT+RAD
             LGONE=.FALSE.
          ELSE
             DMIN=MIN(DMIN,DCENT-RAD)
             DMAX=MAX(DMAX,DCENT+RAD)
          ENDIF
       ELSE
          DO IEND=1,2
             CALL SAL141(TYPE,RPAR(:,ELEM),X,Y,IEND)
             ! project end of element
             DAUX=X*EY0-Y*EX0
             IF(LGONE)THEN
                DMIN=DAUX
                DMAX=DAUX
                LGONE=.FALSE.
             ELSEIF(DAUX<DMIN)THEN
                DMIN=DAUX
             ELSEIF(DAUX>DMAX)THEN
                DMAX=DAUX
             ENDIF
          ENDDO
          IF(TYPE==G_ELE_TYPE(3))THEN
             ! treat tangent to arc of circles
             RAD=RPAR(3,ELEM)
             DCENT=RPAR(1,ELEM)*EY0-RPAR(2,ELEM)*EX0
             THETA1=RPAR(4,ELEM)
             THETA2=RPAR(5,ELEM)
             THETAM=THETA2-TWOPI
             DO ISIDE=1,2
                THETA=THETA0(ISIDE)
                IF((THETA>THETA1.AND.THETA<THETA2).OR.THETA<THETAM) THEN
                   ! check projection of tangent
                   DAUX=DCENT+SIGNV(ISIDE)*RAD
                   IF(DAUX<DMIN)THEN
                      DMIN=DAUX
                   ELSEIF(DAUX>DMAX)THEN
                      DMAX=DAUX
                   ENDIF
                ENDIF
             ENDDO
          ELSE
             IF(TYPE/=G_ELE_TYPE(1)) CALL XABORT('SAL235: not implemented')
          ENDIF
       ENDIF
    ENDDO
    NPIECE=2
    DPIECE(1)=DMIN
    DPIECE(2)=DMAX
    !
  END SUBROUTINE SAL235
  !
  SUBROUTINE SAL237(EX0,EY0,MQ,NQ,PROJTAB,AXIS)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! computes geometry outline projections on the two symmetrical axis
    !
    !Parameters: input
    ! EX0,EY0   horizontal tracking angle cosines
    ! MQ,NQ     cyclic tracking: for a rectangular geometry,
    !           the length of track is SQRT((MQ*A)**2+(NQ*B)**2)
    !           where A and B are rectangular sides.
    !
    !---------------------------------------------------------------------
    !
    USE PRECISION_AND_KINDS, ONLY : PDB
    USE SAL_GEOMETRY_TYPES,  ONLY : LENGTHX,LENGTHY,TYPGEO
    USE SAL_TRACKING_TYPES,  ONLY : DELR,LENGTH_INV_CYCL
    !***
    IMPLICIT NONE
    REAL(PDB), INTENT(IN) :: EX0,EY0
    INTEGER,   INTENT(IN) :: MQ,NQ
    !***
    INTEGER   :: IAXIS,NPIECE,AXIS(2),IMQ
    REAL(PDB) :: X1,X2,DX,DR,R1,R2,PROJTAB(6),NNQ
    REAL, PARAMETER :: EPS3  = 1.0E-3
    !***
    !     minimum radial interval to contain one trajectory
    IF(TYPGEO.LE.7) THEN
      ! Cartesian geometry
      IF(NQ>0) THEN ! angle different from 0
        X2=LENGTHX/NQ
        IF(MQ>0) THEN
          X1=0.
        ELSE
          X1=LENGTHX-X2; X2=LENGTHX
        ENDIF
        R1=X1*EY0; R2=X2*EY0
        NPIECE=INT((R2-R1)/DELR)
        IF(NPIECE==0) NPIECE=1
        DR=(R2-R1)/NPIECE
        DX=DR/EY0
        IAXIS=1
      ELSE ! angle equal to 0
        X1=LENGTHY-LENGTHY/ABS(MQ)
        X2=LENGTHY
        R1=X1*EX0; R2=X2*EX0
        NPIECE=INT((R2-R1)/DELR)
        IF(NPIECE==0) NPIECE=1
        DR=(R2-R1)/NPIECE
        DX=DR/EX0
        IF(TYPGEO.EQ.7) DX=DX*SQRT(2.0)
        IAXIS=2
      ENDIF
      IF(TYPGEO.EQ.7) DR=DR/2.0
      PROJTAB(:)=(/EX0,EY0,X1,X2,DX,DR/)
      ! cyclical track length
      IF(TYPGEO.EQ.5) THEN
        ! translation
        X1=LENGTHX*ABS(MQ); X2=LENGTHY*ABS(NQ)
      ELSE
        ! specular reflexion
        X1=2.*LENGTHX*ABS(MQ); X2=2.*LENGTHY*ABS(NQ)
      ENDIF
    ELSE IF(TYPGEO.GE.8) THEN
      ! hexagonal geometry
      NNQ=(NQ-ABS(MQ))/2
      X2=3.*LENGTHX/REAL(ABS(MQ)+2*NNQ)
      IF(X2<=LENGTHX+EPS3) THEN
         IAXIS=1
         IF(MQ>0) THEN
           X1=0.
         ELSE
           X1=LENGTHX-X2; X2=LENGTHX
         ENDIF
      ELSE
        X1=0.
        IF((TYPGEO==8).OR.(TYPGEO==10).OR.(TYPGEO==12)) THEN
           ! MQ must be positive
           IAXIS=2
           X2=3.*LENGTHX/(MQ-NNQ)
        ELSE
           IF(MQ>0) THEN
              IAXIS=2
              X2=3.*LENGTHX/(2*MQ+NNQ)
           ELSE
              IAXIS=6
              X2=3.*LENGTHX/(2*ABS(MQ)+NNQ)
           ENDIF
        ENDIF
      ENDIF
      R1=X1*EY0; R2=X2*EY0
      NPIECE=INT((R2-R1)/DELR)
      IF(NPIECE==0) NPIECE=1
      DR=(R2-R1)/NPIECE
      DX=DR/EY0
      ! empirical correction of track weight in hexagonal cases (don't ask why)
      IF(TYPGEO==8) THEN
        IMQ=ABS(MQ)
        IF(((IMQ==1).AND.(NQ==15)).OR.((IMQ==7).AND.(NQ==9)).OR.((IMQ==8).AND.(NQ==6))) THEN
          DR=DR/3.0
        ELSE IF(((IMQ==1).AND.(NQ==9)).OR.((IMQ==4).AND.(NQ==6)).OR.((IMQ==5).AND.(NQ==3))) THEN
          DR=DR/3.0
        ELSE IF(((IMQ==1).AND.(NQ==7)).OR.((IMQ==3).AND.(NQ==5)).OR.((IMQ==4).AND.(NQ==2))) THEN
          DR=DR*5.0/12.0
        ELSE IF(((IMQ==1).AND.(NQ==5)).OR.((IMQ==2).AND.(NQ==4)).OR.((IMQ==3).AND.(NQ==1))) THEN
          DR=DR*4.0/9.0
        ELSE IF(((IMQ==2).AND.(NQ==8)).OR.((IMQ==3).AND.(NQ==7)).OR.((IMQ==5).AND.(NQ==1))) THEN
          DR=DR*7.0/15.0
        ELSE IF(((IMQ==4).AND.(NQ==14)).OR.((IMQ==5).AND.(NQ==13)).OR.((IMQ==9).AND.(NQ==1))) THEN
          DR=13.0*DR/27.0
        ENDIF
      ELSE IF(TYPGEO==9) THEN
        IF(ABS(MQ)/NQ > 1) DR=(0.5+1.5*ABS(MQ)/NQ)*DR
      ELSE IF(TYPGEO==10) THEN
        IMQ=ABS(MQ)
        IF(((IMQ==1).AND.(NQ==15)).OR.((IMQ==8).AND.(NQ==6))) THEN
          DR=0.25*DR
        ELSE IF(((IMQ==1).AND.(NQ==9)).OR.((IMQ==5).AND.(NQ==3))) THEN
          DR=0.25*DR
        ELSE IF(((IMQ==1).AND.(NQ==7)).OR.((IMQ==4).AND.(NQ==2))) THEN
          DR=DR*5.0/14.0
        ELSE IF(((IMQ==1).AND.(NQ==5)).OR.((IMQ==3).AND.(NQ==1))) THEN
          DR=DR*2.0/5.0
        ELSE IF(((IMQ==2).AND.(NQ==8)).OR.((IMQ==5).AND.(NQ==1))) THEN
          DR=DR*7.0/16.0
        ELSE IF(((IMQ==4).AND.(NQ==14)).OR.((IMQ==9).AND.(NQ==1))) THEN
          DR=13.0*DR/28.0
        ELSE
          DR=0.5*DR
        ENDIF
      ELSE IF(TYPGEO==11) THEN
        IMQ=ABS(MQ)
        IF(((IMQ==1).AND.(NQ==15)).OR.((IMQ==8).AND.(NQ==6))) THEN
          DR=0.5*DR
        ELSE IF(((IMQ==1).AND.(NQ==9)).OR.((IMQ==5).AND.(NQ==3))) THEN
          DR=0.5*DR
        ELSE IF(((IMQ==1).AND.(NQ==7)).OR.((IMQ==4).AND.(NQ==2))) THEN
          DR=0.7742663247*DR
        ELSE IF(((IMQ==1).AND.(NQ==5)).OR.((IMQ==3).AND.(NQ==1))) THEN
          DR=0.8257638060*DR
        ELSE IF(((IMQ==2).AND.(NQ==8)).OR.((IMQ==5).AND.(NQ==1))) THEN
          DR=0.8863607851*DR
        ELSE IF(((IMQ==4).AND.(NQ==14)).OR.((IMQ==9).AND.(NQ==1))) THEN
          DR=0.9327596923*DR
        ENDIF
      ELSE IF(TYPGEO==12) THEN
        IMQ=ABS(MQ)
        IF(((IMQ==1).AND.(NQ==15)).OR.((IMQ==7).AND.(NQ==9)).OR.((IMQ==8).AND.(NQ==6))) THEN
          DR=DR/6.0
        ELSE IF(((IMQ==1).AND.(NQ==9)).OR.((IMQ==4).AND.(NQ==6)).OR.((IMQ==5).AND.(NQ==3))) THEN
          DR=DR/6.0
        ELSE IF(((IMQ==1).AND.(NQ==7)).OR.((IMQ==3).AND.(NQ==5)).OR.((IMQ==4).AND.(NQ==2))) THEN
          DR=DR*5.0/24.0
        ELSE IF(((IMQ==1).AND.(NQ==5)).OR.((IMQ==2).AND.(NQ==4)).OR.((IMQ==3).AND.(NQ==1))) THEN
          DR=DR*2.0/9.0
        ELSE IF(((IMQ==2).AND.(NQ==8)).OR.((IMQ==3).AND.(NQ==7)).OR.((IMQ==5).AND.(NQ==1))) THEN
          DR=0.232645602188*DR
        ELSE IF(((IMQ==4).AND.(NQ==14)).OR.((IMQ==5).AND.(NQ==13)).OR.((IMQ==9).AND.(NQ==1))) THEN
          DR=13.0*DR/54.0
        ENDIF
      ENDIF
      PROJTAB(:)=(/EX0,EY0,X1,X2,DX,DR/)
      ! cyclical track length
      X1=LENGTHX*ABS(MQ)*1.5D0; X2=LENGTHX*(2*NNQ+ABS(MQ))*SQRT(3.D0)/2.D0
    ENDIF
    AXIS(:)=(/NPIECE,IAXIS/)
    LENGTH_INV_CYCL=1./SQRT(X1*X1+X2*X2)
  END SUBROUTINE SAL237
  !
  SUBROUTINE SAL220_1(ANGLE)
    !
    !---------------------------------------------------------------------
    !
    !Purpose:
    ! computes unit vectors for the two rotative axis
    !
    !Parameters: input
    ! ANGLE    angle between the two rotative axis
    !
    !---------------------------------------------------------------------
    !
    USE PRECISION_AND_KINDS, ONLY : PDB
    USE SAL_GEOMETRY_TYPES, ONLY : LENGTHX,LENGTHY,TYPGEO
    USE SAL_TRACKING_TYPES, ONLY : HX,HY,BX,BY
    !****
    IMPLICIT NONE
    REAL(PDB), INTENT(IN) :: ANGLE
    !****
    !     unit vector for the axis 1:
    HX(1)=1.; HY(1)=0.
    !     unit vector for the axis 2:
    HX(2)=COS(ANGLE); HY(2)=SIN(ANGLE)
    BX(:)=0.; BY(:)=0.
    IF((TYPGEO.EQ.5).OR.(TYPGEO.EQ.6)) THEN
      HX(3:4)=HX(1:2); HY(3:4)=HY(1:2)
      BX(4)=LENGTHX
      BY(3)=LENGTHY
    ELSE IF(TYPGEO.EQ.7) THEN
      HX(3)=0.; HY(3)=1.
      BX(3)=LENGTHX
    ELSE IF((TYPGEO.EQ.8).OR.(TYPGEO.EQ.10)) THEN
      HX(3)=COS(ANGLE*2.); HY(3)=SIN(ANGLE*2.)
      BX(3)=LENGTHX
    ELSE IF(TYPGEO.EQ.9) THEN
      HX(1)=1.; HY(1)=0.; BX(1)=-LENGTHX*0.5; BY(1)=-LENGTHY
      HX(4)=1.; HY(4)=0.; BX(4)=-LENGTHX*0.5; BY(4)=LENGTHY
      HX(2)=COS(ANGLE*2.); HY(2)=SIN(ANGLE*2.); BX(2)=BX(1); BY(2)=BY(1)
      HX(5)=HX(2); HY(5)=HY(2); BX(5)=LENGTHX; BY(5)=0.
      HX(3)=COS(ANGLE); HY(3)=SIN(ANGLE); BX(3)=-LENGTHX; BY(3)=0.
      HX(6)=HX(3); HY(6)=HY(3); BX(6)=LENGTHX*0.5; BY(6)=-LENGTHY
    ELSE IF(TYPGEO.EQ.11) THEN
      BX(3)=LENGTHX*0.5; BY(3)=LENGTHY
      HX(3)=1.; HY(3)=0.
      BX(4)=LENGTHX
      HX(4)=COS(ANGLE); HY(4)=SIN(ANGLE)
    ELSE IF(TYPGEO.EQ.12) THEN
      HX(3)=COS(PI/2.0+ANGLE); HY(3)=SIN(PI/2.0+ANGLE)
      BX(3)=LENGTHX
    ENDIF
    !
  END SUBROUTINE SAL220_1
END MODULE SAL_AUX_MOD
