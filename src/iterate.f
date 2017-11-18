      SUBROUTINE ITERATE
     &           (NNW,     NGP,    NW_IPAR,
     &            CRD,     PHI,    PHN,
     &            A0,      A1,     B,
     &            ALPHA,   MAX_IT, TOL,     X)
C ---------------------------------------------------------------------------
C     Iteratively solves the following system of equations:
C
C     A0 * X (:, 1) + A1 * X (:, 2) = B
C                     Tr * X (:, 2) = D * X (:, 1)
C
C     where A0 and A1 are the Boundary Integral Equations.
C     and Tr and D are the Spline matrices.
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'chk_params.inc'
C
      INTEGER(KIND=IK) NNW
      INTEGER(KIND=IK) NGP
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)
      REAL   (KIND=RK) CRD (4, *)
      REAL   (KIND=RK) PHI (2, *)
      REAL   (KIND=RK) PHN (2, *)
      REAL   (KIND=RK) A0 (NGP, NGP)
      REAL   (KIND=RK) A1 (NGP, NGP)
      REAL   (KIND=RK) B (NGP)
      REAL   (KIND=RK) ALPHA (*)
      INTEGER(KIND=IK) MAX_IT
      REAL   (KIND=RK) TOL (*)
      REAL   (KIND=RK) X (NGP, 2)
C
      LOGICAL(KIND=LK) CONVERGED
      INTEGER(KIND=IK) COUNTER
      REAL   (KIND=RK) RES (2)
      DATA             RES /2*0.0D0/
      REAL   (KIND=RK) Y (NGP, 2)
C
      WRITE (USR_O, *) '=== ITERATE =='
C
      IF (CHK_EQNS) CALL CHECK_EQNS (NGP, A0, A1, X, B)
C
      CONVERGED = .FALSE.
      COUNTER   = 0
C     A0 <- Inverse (A0)
      CALL INVERT (NGP, A0)
C
      DO WHILE (.NOT. CONVERGED .AND. COUNTER. LT. MAX_IT)
C         Y (:, 1) <- X (:, 2)
          CALL DCOPY (NGP, X (1, 2), 1, Y, 1)
C         Y (:, 2) <- B
C         Y (:, 2) <- Y (:, 2) - A1 * Y (:, 1)
C         Y (:, 1) <- INVERSE (A0) * Y (:, 2)
          CALL SOLVE_AX_B (NGP, A0, A1, Y, B)
C         X (:, 1) <- (1 - alpha) * Y (:, 1) + alpha * X (:, 1)
          CALL RESIDUAL (NGP, 1, ALPHA, Y, X, RES)
          CALL SET_TANGENT (NNW, NW_IPAR, CRD,
     &                      PHI, PHN, X (1, 1), X (1, 2))
          CALL DERIVATIVES (1, 1, NNW, NW_IPAR, X (1,1), X (1,2))
C         CHECK CONVERGENCE
          CONVERGED = RES (1) .LT. TOL (1) .AND. RES (2) .LT. TOL (2)
          COUNTER   = COUNTER + 1
          WRITE (USR_O, 1101) RES, COUNTER
      END DO
C
      IF (.NOT. CONVERGED) THEN
         CALL ERROR ('ITERATE', 'No convergence')
      END IF
      WRITE (USR_O, 1001) RES, COUNTER
C
      RETURN
 1001 FORMAT (1X,'ITR:',2E9.1,' USING ', I4,' ITERATIONS')
 1101 FORMAT (1X,'RES:',2E9.1,' ITERATION ', I4)
      END
      SUBROUTINE SOLVE_AX_B (N, A0_INV, A1, X, B)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      INTEGER(KIND=IK) N
      REAL   (KIND=RK) A0_INV (N, N)
      REAL   (KIND=RK) A1 (N, N)
      REAL   (KIND=RK) X (N, 2)
      REAL   (KIND=RK) B (N)
C     X (:, 2) <- B
      CALL DCOPY (N, B, 1, X (1, 2), 1)
C     X (:, 2)  <- X (:, 2) - A1 * X (:, 1)
      CALL DGEMV
     &   ('N', N, N, -1.0D0, A1,     N, X (1, 1), 1, 1.0D0, X (1, 2), 1)
C     X (:, 1) <- INVERSE (A0) * X (:, 2)
      CALL DGEMV
     &   ('N', N, N,  1.0D0, A0_INV, N, X (1, 2), 1, 0.0D0, X (1, 1), 1)
C
      RETURN
      END
      SUBROUTINE INVERT (N, A)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      INTEGER(KIND=IK) N
      REAL   (KIND=RK) A (N, N)
C
      INTEGER(KIND=IK) INFO
      INTEGER(KIND=IK) IPIV (N)
      CHARACTER*48     MSG
      REAL   (KIND=RK) W (3 * N)
C
      CALL DGETRF (N, N, A, N, IPIV, INFO)
      IF (INFO .NE. 0) THEN
         IF (INFO .LT. 0) WRITE (MSG, 201) 'DGETRF', INFO
         IF (INFO .GT. 0) WRITE (MSG, 211) 'DGETRF', INFO
         CALL ERROR ('INVERT', MSG)
      END IF
C
      CALL DGETRI (N, A, N, IPIV, W, 3 * N, INFO)
      IF (INFO .NE. 0) THEN
         IF (INFO .LT. 0) WRITE (MSG, 201) 'DGETRI', INFO
         IF (INFO .GT. 0) WRITE (MSG, 211) 'DGETRI', INFO
         CALL ERROR ('INVERT', MSG)
      END IF
C
      RETURN
  201 FORMAT (1X,A,': ARGUMENT',I2,' IS INCORRECT')
  211 FORMAT (1X,A,': MATRIX BECAME SINGULAR AT ROW',I4)
      END
      SUBROUTINE RESIDUAL (N, K, ALPHA, Y, X, RES)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      INTEGER(KIND=IK) N
      INTEGER(KIND=IK) K
      REAL   (KIND=RK) ALPHA (K)
      REAL   (KIND=RK) Y (N, K)
      REAL   (KIND=RK) X (N, K)
      REAL   (KIND=RK) RES (K)
C
      INTEGER(KIND=IK) I
      INTEGER(KIND=IK) J
      REAL   (KIND=RK) DX
C
      DO J = 1, K
         RES (J) = 0.0D0
         DO I = 1, N
            DX = X (I, J) - Y (I, J)
            RES (J) = MAX (RES (J), ABS (DX))
            X (I, J) = (1.0D0 - ALPHA (J)) * Y (I, J)
     &               +          ALPHA (J)  * X (I, J)
         END DO
      END DO
C
      RETURN
      END
      SUBROUTINE CHECK_EQNS (NGP, A0, A1, X, B)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      INTEGER(KIND=IK) NGP
      REAL   (KIND=RK) A0 (NGP, *)
      REAL   (KIND=RK) A1 (NGP, *)
      REAL   (KIND=RK) X (NGP, *)
      REAL   (KIND=RK) B (*)
C
      INTEGER(KIND=IK) I, J
      REAL   (KIND=RK) S
C
      WRITE (*, *) '=== CHECK EQNS ==='
      DO I = 1, NGP
         S = 0.0D0
         DO J = 1, NGP
            S = S + A0 (I, J) * X (J, 1) + A1 (I, J) * X (J, 2)
         END DO
         WRITE (*, '(1X,I3,'':'',E10.1)') I, S - B (I)
      END DO
      WRITE (*, *) '=== CHECK EQNS ==='
C
      RETURN
      END
