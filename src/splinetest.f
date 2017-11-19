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
      PARAMETER (N  = 17)
      PARAMETER (M  = 2)
      PARAMETER (LD = M * (M + 1))
C
      REAL   (KIND=RK) F (M, M + 1, N)
      INTEGER(KIND=IK) I
      REAL   (KIND=RK) X, H

      H = 1.0D0 / DBLE (N - 1)

      WRITE (USR_O, *) 'SPL_NOT_A_KNOT SPL_FIRST_DERIV:'
      BC (1) = SPL_NOT_A_KNOT
      BC (2) = SPL_FIRST_DERIV
      DO I = 1, N
        X = H * DBLE (I - 1)
        F (1, 1, I) = X * X * X
        F (2, 1, I) = SIN (TWO_PI * X)
        F (1, 3, I) = 3.0 * H * X * X
        F (2, 3, I) = TWO_PI * H * COS (TWO_PI * X)
      END DO
      DO I = N, N
        X = H * DBLE (I - 1)
        F (1, 2, I) = 3.0 * H * X * X
        F (2, 2, I) = TWO_PI * H * COS (TWO_PI * X)
      END DO

      CALL SPLINE (BC, N, M, LD, F (1, 1, 1), F (1, 2, 1))
      F (1, 3, :) = F (1, 2, :) - F (1, 3, :)
      F (2, 3, :) = F (2, 2, :) - F (2, 3, :)
      WRITE (USR_O, 101) F
      WRITE (USR_O, *)

      WRITE (USR_O, *) 'SPL_PERIODIC:'
      BC (1) = SPL_PERIODIC
      BC (2) = SPL_PERIODIC
      DO I = 1, N
        X = H * DBLE (I - 1)
        F (1, 1, I) = COS (TWO_PI * X - H)
        F (2, 1, I) = SIN (TWO_PI * X - H)
        F (1, 3, I) = -TWO_PI * H * F (2, 1, I)
        F (2, 3, I) =  TWO_PI * H * F (1, 1, I)
      END DO
      CALL SPLINE (BC, N, M, LD, F (1, 1, 1), F (1, 2, 1))
      F (1, 3, :) = F (1, 2, :) - F (1, 3, :)
      F (2, 3, :) = F (2, 2, :) - F (2, 3, :)
      WRITE (USR_O, 101) F
      WRITE (USR_O, *)
C
  101 FORMAT (1X, 6E14.6)
      END
