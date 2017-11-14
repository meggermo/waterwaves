      SUBROUTINE BACK_SUBS (BC, N, KK, LD, W, G)
C ---------------------------------------------------------------------------
C     Performs the backsubstitution of the computed spline coefficients
C     depending on the kind of boundary condition. When the spline is
C     periodic then we'll need to do some additional computations before
C     backsubstituting.
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'spl_params.inc'
C
      INTEGER(KIND=IK) BC (2)           ! spline boundary condition
      INTEGER(KIND=IK) N                ! number of nodes
      INTEGER(KIND=IK) KK               ! number of splines
      INTEGER(KIND=IK) LD               ! leading dimension of W
      REAL   (KIND=RK) W (N, KK + 1)    ! spline coeffs. to be copied
      REAL   (KIND=RK) G (LD, N)        ! array to copy result into
C
      IF (BC (1) .EQ. SPL_PERIODIC) THEN
         CALL BACK_SUBS_PER (N, KK, LD, W, G)
      ELSE
         CALL BACK_SUBS_REG (N, KK, LD, W, G)
      END IF
C
      RETURN
      END
      SUBROUTINE BACK_SUBS_REG (N, K, LD, W, G)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      INTEGER(KIND=IK) N
      INTEGER(KIND=IK) K
      INTEGER(KIND=IK) LD
      REAL   (KIND=RK) W (N, K)
      REAL   (KIND=RK) G (LD, N)
C
      INTEGER(KIND=IK) I, J
C
      DO J = 1, K
         DO I = 1, N
            G (J, I) = W (I, J)
         END DO
      END DO
C
      RETURN
      END
      SUBROUTINE BACK_SUBS_PER (N, K, LD, W, G)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      INTEGER(KIND=IK) N
      INTEGER(KIND=IK) K
      INTEGER(KIND=IK) LD
      REAL   (KIND=RK) G (LD, N)
      REAL   (KIND=RK) W ( N, K + 1)
C
      INTEGER(KIND=IK) I, J
      REAL   (KIND=RK) HTD (K + 1)
C
        DO I = 1, K + 1
          HTD (I) = 0.25D0 * (W (1, I) + W (N - 2, I))
        END DO
        DO J = 1, K
           G (J, N - 1) = (W (N - 1, J) - HTD (J)) / (1.0D0 - HTD (K+1))
           DO I = 1, N - 2
              G (J, I) = W (I, J) - G (J, N - 1) * W (I, K + 1)
           END DO
           G (J, N) = G (J, 1)
        END DO
C
        RETURN
        END