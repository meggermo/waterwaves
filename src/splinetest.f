C
      PROGRAM SPLINETEST
C
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'spl_params.inc'
C
      INTEGER(KIND=IK) BC (2)
      INTEGER(KIND=IK) N
      INTEGER(KIND=IK) M
      INTEGER(KIND=IK) LD
      PARAMETER (N  = 11)
      PARAMETER (LD = 5)
      PARAMETER (M  = 1)
      REAL   (KIND=RK) F (LD, N)
      REAL   (KIND=RK) A (N)
C
      INTEGER(KIND=IK) I
      REAL   (KIND=RK) X, H

      BC (1) = SPL_NOT_A_KNOT
      BC (2) = SPL_FIRST_DERIV

      H = 1.0D0 / DBLE (N - 1)
      DO I = 1, N
        X = H * DBLE (I - 1)
        F (1, I) = X * X * X
        F (2, I) = 3.0 * H * X * X
        F (3, I) = 6.0 * H * X
        F (4, I) = F (2, I)
        F (5, I) = F (3, I)
      END DO
      CALL SPLINE (BC, N, M, LD, F (1, 1), F (2, 1))
      F (2, :) = F (2, :) - F (4, :)
      WRITE (USR_O, 101) F
C
  101 FORMAT (1X, 5E14.6)
      END
