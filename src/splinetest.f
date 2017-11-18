C
      PROGRAM SPLINETEST
C
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'spl_params.inc'
      INCLUDE 'rle_params.inc'
C
      INTEGER(KIND=IK) BC (2)
      INTEGER(KIND=IK) N
      INTEGER(KIND=IK) M
      INTEGER(KIND=IK) LD
      PARAMETER (N  = 65)
      PARAMETER (M  = 2)
      PARAMETER (LD = M * (M + 1))
C
      REAL   (KIND=RK) F (M, M + 1, N)
      INTEGER(KIND=IK) I
      REAL   (KIND=RK) X, H

C      BC (1) = SPL_NOT_A_KNOT
C      BC (2) = SPL_NOT_A_KNOT
      BC (1) = SPL_SECOND_DERIV
      BC (2) = SPL_SECOND_DERIV

      H = 1.0D0 / DBLE (N - 1)
      DO I = 1, N
        X = H * DBLE (I - 1)
        F (1, 1, I) = X * X * X
        F (2, 1, I) = SIN (TWO_PI * X)
C        SPL_FIRST_DERIV
C        F (1, 2, I) = 3.0 * H * X * X
C        F (2, 2, I) = TWO_PI * H * COS (TWO_PI * X)
C       SPL_SECOND_DERIV
        F (1, 2, I) = 6.0 * H * X
        F (2, 2, I) = TWO_PI * TWO_PI * H * SIN (TWO_PI * X)
C       Analytical 1st derivative
        F (1, 3, I) = 3.0 * H * X * X
        F (2, 3, I) = TWO_PI * H * COS (TWO_PI * X)
      END DO
      CALL SPLINE (BC, N, M, LD, F (1, 1, 1), F (1, 2, 1))
      F (1, 3, :) = F (2, 1, :) - F (1, 3, :)
      F (2, 3, :) = F (2, 2, :) - F (2, 3, :)
      WRITE (USR_O, 101) F
C
  101 FORMAT (1X, 6E14.6)
      END
