      SUBROUTINE RHS (BC, N, KK, LD, F, R, NRHS)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'spl_params.inc'
C
      INTEGER(KIND=IK) BC (2)
      INTEGER(KIND=IK) N
      INTEGER(KIND=IK) KK
      INTEGER(KIND=IK) LD
      REAL   (KIND=RK) F (LD, N)
      REAL   (KIND=RK) R (N, KK + 1)
      INTEGER(KIND=IK) NRHS
C
      INTEGER(KIND=IK)  I, J
C
      NRHS = KK
      DO J = 1, KK
         DO I = 2, N - 1
            R (I, J) = 0.75D0 * (F (J, I + 1) - F (J, I - 1))
         END DO
      END DO
      IF (BC (1) .EQ. SPL_PERIODIC) THEN
         R (1, KK + 1) = 0.25D0
         DO I = 2, N - 3
            R (I, KK + 1) = 0.0D0
         END DO
         R (N - 2, KK + 1) = 0.25D0
         NRHS = KK + 1
      END IF
C
      RETURN
      END
