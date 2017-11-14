      SUBROUTINE ITERATE_DOUBLE
     &           (NGP,
     &            A0,  A1,
     &            S0,  S1,
     &            B0,  B1,
     &            OMG, MAX_IT, TOL,
     &            X0,  X1)
C ---------------------------------------------------------------------------
C     Iteratively solves the following system of equations:
C
C     A0 * X0 + A1 * X1 = B0
C     S0 * X0 + S1 * X1 = B1
C
C     where A0 and A1 are the Boundary Integral Equations.
C     and S0 and S1 are the Spline matrices.
C
C     X0 = INVERSE (A0) * (B0 - A1 * X1)
C     X1 = INVERSE (S1) * (B1 - S0 * X0)
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'fle_params.inc'
C
      INTEGER(KIND=IK) NGP
      REAL   (KIND=RK) A0  (NGP * 2, *)
      REAL   (KIND=RK) A1  (NGP * 2, *)
      REAL   (KIND=RK) S0  (NGP * 2, *)
      REAL   (KIND=RK) S1  (NGP * 2, *)
      REAL   (KIND=RK) B0  (*)
      REAL   (KIND=RK) B1  (*)
      REAL   (KIND=RK) OMG (*)
      INTEGER(KIND=IK) MAX_IT
      REAL   (KIND=RK) TOL (*)
      REAL   (KIND=RK) X0  (*)
      REAL   (KIND=RK) X1  (*)
C
      LOGICAL(KIND=LK) CONVERGED
      INTEGER(KIND=IK) COUNTER
      INTEGER(KIND=IK) LDA
      REAL   (KIND=RK) RES (2)
      REAL   (KIND=RK) Y0 (NGP)
      REAL   (KIND=RK) Y1 (NGP)
      INTEGER(KIND=IK) I, J
C
      CALL CHECK_EQNS_DOUBLE (NGP, A0, A1, S0, S1, X0, X1, B0, B1)
C
      CONVERGED = .FALSE.
      COUNTER   = 0
      LDA       = NGP * 2
C
      CALL INVERT_DOUBLE (NGP, LDA, A0)
      CALL INVERT_DOUBLE (NGP, LDA, S1)
C
      DO WHILE (.NOT. CONVERGED .AND. COUNTER. LT. MAX_IT)
C         Y1 <- B0
          CALL DCOPY (NGP, B0, 1, Y1, 1)
C         Y1 <- Y1 - A1 * X1
          CALL DGEMV ('N', NGP,NGP, -1.0D0, A1, LDA, X1, 1, 1.0D0, Y1,1)
C         Y0 <- INVERSE (A0) * Y1
          CALL DGEMV ('N', NGP,NGP,  1.0D0, A0, LDA, Y1, 1, 0.0D0, Y0,1)
          CALL RESIDUAL_DOUBLE (NGP, OMG (1), Y0, X0, RES (1))
C         Y0 <- B1
          CALL DCOPY (NGP, B1, 1, Y0, 1)
C         Y0 <- Y0 - S0 * X0
          CALL DGEMV ('N', NGP,NGP, -1.0D0, S0, LDA, X0, 1, 1.0D0, Y0,1)
C         Y1 <- INVERSE (S1) * Y0
          CALL DGEMV ('N', NGP,NGP,  1.0D0, S1, LDA, Y0, 1, 0.0D0, Y1,1)
          CALL RESIDUAL_DOUBLE (NGP, OMG (2), Y1, X1, RES (2))
C         CHECK CONVERGENCE
          CONVERGED = RES (1) .LT. TOL (1) .AND. RES (2) .LT. TOL (2)
          COUNTER   = COUNTER + 1
          WRITE (USR_O, 1101) RES (1), RES (2), COUNTER
      END DO
C
      IF (.NOT. CONVERGED) THEN
         CALL ERROR ('ITERATE', 'No convergence')
      END IF
      WRITE (USR_O, 1001) RES, COUNTER
C
      RETURN
  201 FORMAT (1X,I3,':',4E14.6)
  101 FORMAT (1X,4E14.6)
 1001 FORMAT (1X,'ITR:',2E9.1,' USING ', I3,' ITERATIONS')
 1101 FORMAT (1X,'ITR:',2E9.1,' AT ITERATION ', I3)
      END
      SUBROUTINE SOLVE_AX_B_DOUBLE (N, LDA, AI, A1, X0, X1, B)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      INTEGER(KIND=IK) N
      INTEGER(KIND=IK) LDA
      REAL   (KIND=RK) AI (LDA, *)
      REAL   (KIND=RK) A1 (LDA, *)
      REAL   (KIND=RK) X0 (*)
      REAL   (KIND=RK) X1 (*)
      REAL   (KIND=RK) B  (*)
C     X1 <- B
      CALL DCOPY (N, B, 1, X1, 1)
C     X1  <- X1 - A1 * X0
      CALL DGEMV ('N', N, N, -1.0D0, A1, LDA, X0, 1, 1.0D0, X1, 1)
C     X0 <-  INVERSE (A0) * X1
      CALL DGEMV ('N', N, N,  1.0D0, AI, LDA, X1, 1, 0.0D0, X0, 1)
C
      RETURN
      END
      SUBROUTINE INVERT_DOUBLE (N, LDA, A)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      INTEGER(KIND=IK) N
      INTEGER(KIND=IK) LDA
      REAL   (KIND=RK) A (LDA, *)
C
      INTEGER(KIND=IK) INFO
      INTEGER(KIND=IK) IPIV (N)
      REAL   (KIND=RK) W (3 * N)
      CHARACTER*48 MSG
C
      CALL DGETRF (N, N, A, LDA, IPIV,           INFO)
      IF (INFO .NE. 0) THEN
         IF (INFO .LT. 0) WRITE (MSG, 201) 'DGETRF', INFO
         IF (INFO .GT. 0) WRITE (MSG, 211) 'DGETRF', INFO
         CALL ERROR ('INVERT_DOUBLE',MSG)
      END IF
      CALL DGETRI (N,    A, LDA, IPIV, W, 3 * N, INFO)
      IF (INFO .NE. 0) THEN
         IF (INFO .LT. 0) WRITE (MSG, 201) 'DGETRI', INFO
         IF (INFO .GT. 0) WRITE (MSG, 211) 'DGETRI', INFO
         CALL ERROR ('INVERT_DOUBLE',MSG)
      END IF
C
      RETURN
  201 FORMAT (1X,A,': ARGUMENT',I2,' IS INCORRECT')
  211 FORMAT (1X,A,': MATRIX BECAME SINGULAR AT ROW',I4)
      END
      SUBROUTINE RESIDUAL_DOUBLE (N, OMG, Y, X, RES)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      INTEGER(KIND=IK) N
      REAL   (KIND=RK) OMG
      REAL   (KIND=RK) Y (*)
      REAL   (KIND=RK) X (*)
      REAL   (KIND=RK) RES
C
      INTEGER(KIND=IK) I
      REAL   (KIND=RK) DX
C
      RES = 0.0D0
      DO I = 1, N
         DX    = X (I) - Y (I)
         RES   = MAX (RES, ABS (DX))
         X (I) = (1.0D0 - OMG) * Y (I) + OMG  * X (I)
      END DO
C
      RETURN
      END