SUBROUTINE VALUE6(IELEM,NUN,LXH,LZ,NBLOS,ISPLH,X,Y,Z,SIDE,ZZZ,EVECT,ISS,KFLX, &
& IXLG,IYLG,IZLG,SXYZ,VNORM,MATXYZ,VALUE)
  !
  !-----------------------------------------------------------------------
  !
  !Purpose:
  ! Interpolate the flux distribution for DUAL method in 3D hexagonal.
  !
  !Copyright:
  ! Copyright (C) 2026 Ecole Polytechnique de Montreal
  ! This library is free software; you can redistribute it and/or
  ! modify it under the terms of the GNU Lesser General Public
  ! License as published by the Free Software Foundation; either
  ! version 2.1 of the License, or (at your option) any later version
  !
  !Author(s): A. Hebert
  !
  !Parameters: input
  ! IELEM   finite element order
  !         =1 : linear Raviart-Thomas
  !         =2 : parabolic Raviart-Thomas
  !         =3 : cubic Raviart-Thomas
  !         =4 : quartic Raviart-Thomas
  ! NUN     number of unknowns.
  ! LXH     number of hexagons in the 2D plane (1, 7, 19, 37, etc.).
  ! LZ      number of elements along the Z axis.
  ! NBLOS   number of lozenges per direction in 3D with mesh-splitting.
  ! ISPLH   lozenge splitting (=1: 3 lozenges, =2: 12 lozenges, etc.).
  ! X       Cartesian coordinates along the X axis where the flux is
  !         interpolated.
  ! Y       Cartesian coordinates along the Y axis where the flux is
  !         interpolated.
  ! Z       Cartesian coordinates along the Z axis where the flux is
  !         interpolated.
  ! SIDE    side of a lozenge.
  ! ZZZ     Cartesian coordinates along the Z axis.
  ! EVECT   variational coefficients of the flux.
  ! ISS     mixture index assigned to each element.
  ! KFLX    position of fluxes in the unknown vector.
  ! IXLG    number of interpolated points according to X.
  ! IYLG    number of interpolated points according to Y.
  ! IZLG    number of interpolated points according to Z.
  ! SXYZ    voxel side dimensions.
  !
  !Parameters: output
  ! VNORM   2D voxel side normalization factor.
  ! MATXYZ  mixture index assigned to each interpolation point.
  ! VALUE   interpolated fluxes.
  !
  !-----------------------------------------------------------------------
  !
  IMPLICIT NONE
  !----
  !  SUBROUTINE ARGUMENTS
  !----
  INTEGER IELEM,NUN,LXH,LZ,NBLOS,ISPLH,IXLG,IYLG,IZLG,ISS(3,NBLOS),KFLX(3,NBLOS), &
  & MATXYZ(IXLG,IYLG,IZLG)
  REAL X(IXLG),Y(IYLG),Z(IZLG),SIDE,ZZZ(LZ+1),EVECT(NUN),SXYZ(3),VNORM,VALUE(IXLG,IYLG,IZLG)
  !----
  !  LOCAL VARIABLES
  !----
  INTEGER I,J,K,L,ISU,I1,I2,I3,IE,IC,ICMIN,ICMAX,NBC,ISPX,ISPY,IND,IUN,IEL,KS
  LOGICAL FVALUE
  DOUBLE PRECISION COEF(2,5),FLX(5),FLY(5),FLZ(5),U,V,W,AXY(2),BXY(2),CXY(2),CENTER(2), &
     & SUMMIT(2,6),NRING,ANGLE,RADIUS,SIDEH,XYC(3),XYR(2),DVOL1,DVOL2,VSIZE
  CHARACTER HSMG*131
  DOUBLE PRECISION, PARAMETER :: PI=3.141592653589793D0
  DOUBLE PRECISION, PARAMETER :: SQRT3=SQRT(3.0D0)
  DOUBLE PRECISION, PARAMETER :: EPS=1.0D-6
  !----
  !  ALLOCATABLE ARRAYS
  !----
  TYPE VECTOR_ARRAY
    DOUBLE PRECISION, DIMENSION(2) :: ORIGIN ! origin of each lozenge
    DOUBLE PRECISION, DIMENSION(2) :: AXY_S,BXY_S,CXY_S,DXY_S
    DOUBLE PRECISION :: XMIN,XMAX,YMIN,YMAX
  END TYPE VECTOR_ARRAY
  TYPE(VECTOR_ARRAY),POINTER :: IP
  TYPE(VECTOR_ARRAY), ALLOCATABLE, TARGET, DIMENSION(:,:,:,:) :: LOZ
  !----
  !  compute coefficient for Legendre polynomials
  !----
  COEF(:2,:5)=0.0
  COEF(1,1)=1.0
  COEF(1,2)=2.*3.**0.5
  DO IE=1,3
    COEF(1,IE+2)=2.0*REAL(2*IE+1)/REAL(IE+1)*(REAL(2*IE+3)/REAL(2*IE+1))**0.5
    COEF(2,IE+2)=REAL(IE)/REAL(IE+1)*(REAL(2*IE+3)/REAL(2*IE-1))**0.5
  ENDDO
  !----
  !  2D hexagonal geometry analysis
  !----
  NRING = (1.0+SQRT(REAL(1+(LXH-1)*4/3)))/2.0
  IF(INT(NRING) /= NRING) THEN
    CALL XABORT('VALUE6: INVALID NUMBER OF HEXAGONS IN THE PLANE.')
  ENDIF
  SIDEH=SIDE*ISPLH
  SUMMIT(1,1)=SIDEH/2  ; SUMMIT(2,1)=-SIDEH*SQRT3/2
  SUMMIT(1,2)=SIDEH    ; SUMMIT(2,2)=0
  SUMMIT(1,3)=SIDEH/2  ; SUMMIT(2,3)=SIDEH*SQRT3/2
  SUMMIT(1,4)=-SIDEH/2 ; SUMMIT(2,4)=SIDEH*SQRT3/2
  SUMMIT(1,5)=-SIDEH   ; SUMMIT(2,5)=0
  SUMMIT(1,6)=-SIDEH/2 ; SUMMIT(2,6)=-SIDEH*SQRT3/2
  ICMIN=0
  ALLOCATE(LOZ(LXH,3,ISPLH,ISPLH))
  DO IC=0,LXH-1
    CENTER(:)=0.0D0
    NBC = CEILING((1.0+SQRT(1.0+IC*4.0/3.0))/2.0)
    ICMAX = 3*(NBC*(NBC-1))
    IF(NBC == 2) THEN
      ANGLE = PI/6 +(IC-ICMIN)*PI/(3*(NBC-1))
      CENTER(:) = SQRT3*SIDEH* (/COS(ANGLE), SIN(ANGLE)/)
    ELSE IF(NBC > 2) THEN
      IF(IC-ICMIN > 2*(NBC-1) .AND. IC-ICMIN < 3*(NBC-1)) THEN
        IND = IC-ICMIN-2*(NBC-1)
        CENTER(:) = (/-1.5D0*(NBC-1),0.5D0*(NBC-1-2*IND)*SQRT3/)*SIDEH
      ELSE IF(IC-ICMIN > 5*(NBC-1).AND.IC-ICMIN < 6*(NBC-1)) THEN
        IND = IC-ICMIN-5*(NBC-1)
        CENTER(:) = (/1.5D0*(NBC-1),-0.5D0*(NBC-1-2*IND)*SQRT3/)*SIDEH
      ELSE
        ANGLE = PI/6+(IC-ICMIN)*PI/(3*(NBC-1))
        RADIUS = 0.0D0
        IF(ANGLE > 0 .AND. ANGLE <= PI/2) THEN
          RADIUS = 1.5*(NBC-1)*SIDEH/SIN(PI/6+ANGLE)
        ELSE IF(ANGLE > PI/2 .AND. ANGLE <= PI) THEN
          RADIUS = 1.5*(NBC-1)*SIDEH/SIN(7*PI/6-ANGLE)
        ELSE IF(ANGLE > PI .AND. ANGLE <= 3*PI/2) THEN
          RADIUS = -1.5*(NBC-1)*SIDEH/SIN(PI/6+ANGLE)
        ELSE IF(ANGLE > 3*PI/2 .AND. ANGLE <= 2*PI) THEN
          RADIUS = -1.5*(NBC-1)*SIDEH/SIN(7*PI/6-ANGLE)
        ENDIF
        CENTER(:) = RADIUS*(/ SIN(PI/2-ANGLE), COS(PI/2-ANGLE) /)
      ENDIF
    ENDIF
    IF(IC == ICMAX) ICMIN=ICMAX+1
    DO ISU=1,6,2
      AXY(:)=CENTER(:)+SUMMIT(:,ISU)
      BXY(:)=CENTER(:)+SUMMIT(:,ISU+1)
      CXY(:)=CENTER(:)+SUMMIT(:,1+MOD(ISU+1,6))
      DO ISPY=0,ISPLH-1
        DO ISPX=0,ISPLH-1
          IP => LOZ(IC+1,(ISU+1)/2,ISPX+1,ISPY+1)
          IF(ISU==1) THEN
            IP%AXY_S(:)=AXY(:)+(ISPLH-ISPX-1)*(BXY(:)-AXY(:))/ISPLH+ISPY*(CXY(:)-BXY(:))/ISPLH
            IP%BXY_S(:)=AXY(:)+(ISPLH-ISPX)*(BXY(:)-AXY(:))/ISPLH+ISPY*(CXY(:)-BXY(:))/ISPLH
            IP%CXY_S(:)=AXY(:)+(ISPLH-ISPX)*(BXY(:)-AXY(:))/ISPLH+(ISPY+1)*(CXY(:)-BXY(:))/ISPLH
            IP%DXY_S(:)=AXY(:)+(ISPLH-ISPX-1)*(BXY(:)-AXY(:))/ISPLH+(ISPY+1)*(CXY(:)-BXY(:))/ISPLH
          ELSE IF(ISU==3) THEN
            IP%AXY_S(:)=AXY(:)+ISPY*(BXY(:)-AXY(:))/ISPLH+ISPX*(CXY(:)-BXY(:))/ISPLH
            IP%BXY_S(:)=AXY(:)+(ISPY+1)*(BXY(:)-AXY(:))/ISPLH+ISPX*(CXY(:)-BXY(:))/ISPLH
            IP%CXY_S(:)=AXY(:)+(ISPY+1)*(BXY(:)-AXY(:))/ISPLH+(ISPX+1)*(CXY(:)-BXY(:))/ISPLH
            IP%DXY_S(:)=AXY(:)+ISPY*(BXY(:)-AXY(:))/ISPLH+(ISPX+1)*(CXY(:)-BXY(:))/ISPLH
          ELSE IF(ISU==5) THEN
            IP%AXY_S(:)=AXY(:)+ISPX*(BXY(:)-AXY(:))/ISPLH+(ISPLH-ISPY-1)*(CXY(:)-BXY(:))/ISPLH
            IP%BXY_S(:)=AXY(:)+(ISPX+1)*(BXY(:)-AXY(:))/ISPLH+(ISPLH-ISPY-1)*(CXY(:)-BXY(:))/ISPLH
            IP%CXY_S(:)=AXY(:)+(ISPX+1)*(BXY(:)-AXY(:))/ISPLH+(ISPLH-ISPY)*(CXY(:)-BXY(:))/ISPLH
            IP%DXY_S(:)=AXY(:)+ISPX*(BXY(:)-AXY(:))/ISPLH+(ISPLH-ISPY)*(CXY(:)-BXY(:))/ISPLH
          ENDIF
          IP%ORIGIN(:)=0.5D0*(IP%AXY_S(:)+IP%DXY_S(:))
          IP%XMIN=MIN(IP%AXY_S(1),IP%BXY_S(1),IP%CXY_S(1),IP%DXY_S(1))
          IP%XMAX=MAX(IP%AXY_S(1),IP%BXY_S(1),IP%CXY_S(1),IP%DXY_S(1))
          IP%YMIN=MIN(IP%AXY_S(2),IP%BXY_S(2),IP%CXY_S(2),IP%DXY_S(2))
          IP%YMAX=MAX(IP%AXY_S(2),IP%BXY_S(2),IP%CXY_S(2),IP%DXY_S(2))
        ENDDO ! ISPY
      ENDDO ! ISPX
    ENDDO ! ISU
  ENDDO ! IC
  !----
  !  compute the exact 3D volume
  !----
  DVOL1=0.0D0
  DO KS=1,LZ
    DO IC=1,LXH
      DO ISPY=1,ISPLH
        DO ISPX=1,ISPLH
          IEL=((KS-1)*LXH*ISPLH+((IC-1)*ISPLH+ISPY-1))*ISPLH+ISPX
          DO ISU=1,3
            IF(ISS(ISU,IEL)>0) DVOL1=DVOL1+2.0D0*SQRT(3.0D0)*(ZZZ(KS+1)-ZZZ(KS))*SIDE**2/REAL(ISPLH**2)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !----
  !  perform interpolation
  !----
  MATXYZ(:IXLG,:IYLG,:IZLG)=0
  VALUE(:IXLG,:IYLG,:IZLG)=0.0
  DVOL2=0.0D0
  DO 120 K=1,IZLG
  DO 110 J=1,IYLG
  DO 100 I=1,IXLG
  XYC(:)=(/X(I),Y(J),Z(K)/)
  !
  ! Find the finite element index containing the interpolation point
  KS=0
  DO L=1,LZ
    KS=L
    IF((Z(K).GE.ZZZ(L)).AND.(Z(K).LE.ZZZ(L+1))) GO TO 10
  ENDDO
  CALL XABORT('VALUE6: WRONG INTERPOLATION ALONG Z AXIS.')
  10 DO IC=1,LXH
    DO ISPY=1,ISPLH
      DO ISPX=1,ISPLH
        IEL=((KS-1)*LXH*ISPLH+((IC-1)*ISPLH+ISPY-1))*ISPLH+ISPX
        DO ISU=1,3
          IF(ISS(ISU,IEL) == 0) CYCLE
          IP => LOZ(IC,ISU,ISPX,ISPY)
          IF((XYC(1) < IP%XMIN).OR.(XYC(1) > IP%XMAX)) CYCLE
          IF((XYC(2) < IP%YMIN).OR.(XYC(2) > IP%YMAX)) CYCLE
          FVALUE=VAL_CONTAINS(IP%AXY_S,IP%BXY_S,IP%CXY_S,IP%DXY_S,XYC(1),XYC(2))
          IF(FVALUE) THEN
            MATXYZ(I,J,K)=ISS(1,IEL)
            VSIZE=SXYZ(1)*SXYZ(2)*SXYZ(3)
            IF((I==1).OR.(I==IXLG)) VSIZE=VSIZE/2.0D0
            IF((J==1).OR.(J==IYLG)) VSIZE=VSIZE/2.0D0
            IF((IZLG.GT.1).AND.((K==1).OR.(K==IZLG))) VSIZE=VSIZE/2.0D0
            DVOL2=DVOL2+VSIZE
            ! lozenge rotation
            ANGLE=PI/6.0D0+(ISU-1)*PI*2.0D0/3.0D0
            XYR(1)=(XYC(1)-IP%ORIGIN(1))*COS(ANGLE)+(XYC(2)-IP%ORIGIN(2))*SIN(ANGLE)
            XYR(2)=-(XYC(1)-IP%ORIGIN(1))*SIN(ANGLE)+(XYC(2)-IP%ORIGIN(2))*COS(ANGLE)
            ! piola transformation
            U=2.0D0*XYR(1)/(SQRT3*SIDE)-0.5D0
            V=XYR(1)/(SQRT3*SIDE)-XYR(2)/SIDE
            IF((U < -0.5-EPS).OR.(U > 0.5+EPS)) THEN
              WRITE(HSMG,200) 'U',U
              CALL XABORT(HSMG)
            ELSE IF((V < -0.5-EPS).OR.(V > 0.5+EPS)) THEN
              WRITE(HSMG,200) 'V',V
              CALL XABORT(HSMG)
            ENDIF
            FLX(1)=COEF(1,1)
            FLX(2)=COEF(1,2)*U
            FLY(1)=COEF(1,1)
            FLY(2)=COEF(1,2)*V
            W=(XYC(3)-0.5*(ZZZ(KS)+ZZZ(KS+1)))/(ZZZ(KS+1)-ZZZ(KS))
            FLZ(1)=COEF(1,1)
            FLZ(2)=COEF(1,2)*W
            IF(IELEM.GE.2) THEN
              DO IE=2,IELEM
                FLX(IE+1)=FLX(IE)*U*COEF(1,IE+1)-FLX(IE-1)*COEF(2,IE+1)
                FLY(IE+1)=FLY(IE)*V*COEF(1,IE+1)-FLY(IE-1)*COEF(2,IE+1)
                FLZ(IE+1)=FLZ(IE)*W*COEF(1,IE+1)-FLZ(IE-1)*COEF(2,IE+1)
              ENDDO
            ENDIF
            DO I3=1,IELEM
              DO I2=1,IELEM
                DO I1=1,IELEM
                  IF(KFLX(ISU,IEL)==0) CYCLE
                  IUN=KFLX(ISU,IEL)+((I3-1)*IELEM+I1-1)*IELEM+I2-1
                  VALUE(I,J,K)=VALUE(I,J,K)+EVECT(IUN)*REAL(FLX(I1)*FLY(I2)*FLZ(I3))
                ENDDO
              ENDDO
            ENDDO
            GO TO 100
          ENDIF
        ENDDO ! ISU
      ENDDO ! ISPX
    ENDDO ! ISPY
  ENDDO ! IC
  100 CONTINUE
  110 CONTINUE
  120 CONTINUE
  VNORM=REAL((DVOL1/DVOL2)**(1.0D0/3.0D0))
  DEALLOCATE(LOZ)
  RETURN
  200 FORMAT(8HVALUE6: ,A1,1H=,1P,E15.8,31H IS OUTSIDE (-1/2,1/2) SUPPORT.)
  !
  CONTAINS
    FUNCTION VAL_CONTAINS(AXY, BXY, CXY, DXY, X, Y) RESULT(FVALUE)
      DOUBLE PRECISION,INTENT(IN) :: AXY(2),BXY(2),CXY(2),DXY(2),X,Y
      LOGICAL FVALUE,FVALUE1,FVALUE2
      DOUBLE PRECISION DET
      !
      DET = (BXY(1) - AXY(1)) * (DXY(2) - AXY(2)) - (BXY(2) - AXY(2)) * (DXY(1) - AXY(1))
      FVALUE1 = DET * ((BXY(1) - AXY(1)) * (Y - AXY(2)) - (BXY(2) - AXY(2)) * (X - AXY(1))) >= 0.0D0 .AND. &
              & DET * ((DXY(1) - BXY(1)) * (Y - BXY(2)) - (DXY(2) - BXY(2)) * (X - BXY(1))) >= 0.0D0 .AND. &
              & DET * ((AXY(1) - DXY(1)) * (Y - DXY(2)) - (AXY(2) - DXY(2)) * (X - DXY(1))) >= 0.0D0
      DET = (CXY(1) - BXY(1)) * (DXY(2) - BXY(2)) - (CXY(2) - BXY(2)) * (DXY(1) - BXY(1))
      FVALUE2 = DET * ((CXY(1) - BXY(1)) * (Y - BXY(2)) - (CXY(2) - BXY(2)) * (X - BXY(1))) >= 0.0D0 .AND. &
              & DET * ((DXY(1) - CXY(1)) * (Y - CXY(2)) - (DXY(2) - CXY(2)) * (X - CXY(1))) >= 0.0D0 .AND. &
              & DET * ((BXY(1) - DXY(1)) * (Y - DXY(2)) - (BXY(2) - DXY(2)) * (X - DXY(1))) >= 0.0D0
      FVALUE = FVALUE1 .OR. FVALUE2
    END FUNCTION VAL_CONTAINS
END SUBROUTINE VALUE6
