SUBROUTINE FPAND(N,M,ITER,XX,GG,XXNEW)
!
!-----------------------------------------------------------------------
!
!Purpose:
! Anderson acceleration of fixed point iterations.
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
! N       number of unknowns.
! M       window size. Should be smaller or equal to N.
! ITER    iteration index. ITER=2 at the first call.
! XX      past iterate history over the window.
! GG      past gradient history over the window.
!
!Parameters: output
! XXNEW   single iterate after the window.
!
!-----------------------------------------------------------------------
  IMPLICIT NONE
  !
  !----
  !  SUBROUTINE ARGUMENTS
  !----
  INTEGER, INTENT(IN) :: N,M,ITER
  DOUBLE PRECISION, INTENT(IN) :: GG(N,M+1),XX(N,M+1)
  DOUBLE PRECISION, INTENT(OUT) :: XXNEW(N)
  !----
  !  LOCAL VARIABLES
  !----
  INTEGER :: M_K,R
  !----
  !  ALLOCATABLE ARRAYS
  !----
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: S,D,RHS
  DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: X_K,G_K,U,V
  !----
  !  SCRATCH STORAGE ALLOCATION
  !----
  M_K = MIN(ITER-1, M)
  ALLOCATE(X_K(N,M_K),G_K(N,M_K),U(N,M_K),S(M_K),V(M_K,M_K),D(M_K),RHS(M_K))
  !----
  !  SINGULAR VALUE DECOMPOSITION
  !  [U,S,V] = SVD(G_K,"ECON") ! equivalent Matlab command
  !----
  X_K = FPDIFF(N,M_K,M+1,XX) ; G_K = FPDIFF(N,M_K,M+1,GG)
  U(:N,:M_K) = G_K(:N,:M_K)
  CALL ALSVDF(U,N,M_K,N,M_K,S,V)
  !----
  !  DETERMINE THE EFFECTIVE RANK R OF G_K
  !----
  R = 0
  DO WHILE ( R < M_k .AND. S(R+1) > 1.0E-6*S(1) )
    R = R+1
  ENDDO
  !----
  !  ACCELERATION
  !----
  D = MATMUL(TRANSPOSE(U), GG(:N,ID(ITER,M+1)))
  RHS(:R) = D(:R)/S(:R); RHS(R+1:M_K) = 0.0D0
  XXNEW(:N) = XX(:N,ID(ITER,M+1)) + GG(:N,ID(ITER,M+1)) -  &
       & MATMUL((X_K(:N,:M_K) +G_K(:N,:M_K)), MATMUL(V,RHS))
  !----
  !  SCRATCH STORAGE DEALLOCATION
  !----
  DEALLOCATE(RHS,D,V,S,U,G_K,X_K)
  RETURN
CONTAINS
  !
  FUNCTION ID(I,IS) RESULT(J)
    !----
    !  Cyclic indexing
    !----
    INTEGER, INTENT(IN) :: I,IS
    INTEGER :: J
    J=1+MOD(I-1,IS)
    RETURN
  END FUNCTION ID
  !
  FUNCTION FPDIFF(N,M_K,IS,FF) RESULT(DFF)
    !----
    !  Compute the differences along the second axis
    !----
    INTEGER, INTENT(IN) :: N,M_K,IS
    DOUBLE PRECISION, INTENT(IN) :: FF(N,M_K+1)
    DOUBLE PRECISION :: DFF(N,M_K)
    INTEGER :: I
    DO I=1,M_K
      DFF(:,I)=FF(:,ID(I+1,IS))-FF(:,ID(I,IS))
    ENDDO
    RETURN
  END FUNCTION FPDIFF
END SUBROUTINE FPAND
