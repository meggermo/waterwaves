      SUBROUTINE Spline (BC, N, K, LD, F, G)
C ---------------------------------------------------------------------------
C     Computes the derivative coefficients of K Hermite
C     interpolation splines.
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      INTEGER(KIND=IK) BC (2)
      INTEGER(KIND=IK) N
      INTEGER(KIND=IK) K
      INTEGER(KIND=IK) LD
      REAL   (KIND=RK) F (LD, N)
      REAL   (KIND=RK) G (LD, N)
C
      INTEGER(KIND=IK)  I, NRHS, NEQS, INFO
      REAL   (KIND=RK) W (N, K + 3)
C
      IF (BC (1) .GE. 0 .AND . BC (2) .GE. 0) THEN
        CALL SETUP_DIAGS (N, W)
        CALL RHS (BC, N, K, LD, F, W (1, 3), NRHS)
        CALL BCD (BC, N, K, LD, F, W (1, 3), W, G, I, NEQS)
        CALL DPTSV (NEQS, NRHS, W (I, 1), W (I, 2), W (I, 3), N, INFO)
        CALL BACK_SUBS (BC, N, K, LD, W (1, 3), G)
      END IF
C
      RETURN
      END