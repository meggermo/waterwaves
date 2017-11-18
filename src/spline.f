      SUBROUTINE Spline (BC, N, K, LD, F, G)
C ---------------------------------------------------------------------------
C     Computes the derivative coefficients of K Hermite
C     interpolation splines.
C
C     BC (2): Start and end boundary conditions
C     Available types are:
C       SPL_NOT_A_KNOT   = 0
C       SPL_FIRST_DERIV  = 1
C       SPL_SECOND_DERIV = 2
C       SPL_FIRST_D_HALF = 3
C       SPL_PERIODIC     = 4
C       SPL_NATURAL      = 5
C     If BC(*) is less than 0 spline computation is skipped.
C
C     N:      Number of nodes
C     K:      Number of splines to compute
C     LD:     Leading dimension of F and G , LD >= K
C     F (LD, N):
C     G (LD, N):
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
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
      ELSE
        WRITE (USR_O, *) 'BCs <= 0, spline computation skipped'
      END IF
C
      RETURN
      END

      SUBROUTINE SETUP_DIAGS (N, F)
C ---------------------------------------------------------------------------
C
C     Set the diagonals of the spline system
C
C     1/4 * F'   +  F'  + 1/4 * F'    =  3/4 * (F   -  F)
C            i-1     i           i+1             i+1    i-1
C
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      INTEGER(KIND=IK) N
      REAL   (KIND=RK) F (N, *)
C
      F (:, 1) = 1.00D0
      F (:, 2) = 0.25D0
C
      RETURN
      END

      SUBROUTINE RHS (BC, N, M, LD, F, R, NRHS)
C ---------------------------------------------------------------------------
C     Computes the righthand sides of the Hermite spline equations.
C     BC (2): Start and end boundary conditions
C     N:     Number of nodes
C     M:     Number of splines to compute
C     LD:    Leading dimension of F and G , LD >= K
C     F (LD, N):  Variables to compute first derivatives for
C     R (LD, KK + 1): Righthand sides
C     NRHS:   Number of righthand sides (output)
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'spl_params.inc'
C
      INTEGER(KIND=IK) BC (2)
      INTEGER(KIND=IK) N
      INTEGER(KIND=IK) M
      INTEGER(KIND=IK) LD
      REAL   (KIND=RK) F (LD, N)
      REAL   (KIND=RK) R (N, M + 1)
      INTEGER(KIND=IK) NRHS
C
      INTEGER(KIND=IK)  I, J
C
      NRHS = M
      DO J = 1, M
         DO I = 2, N - 1
            R (I, J) = 0.75D0 * (F (J, I + 1) - F (J, I - 1))
         END DO
      END DO
      IF (BC (1) .EQ. SPL_PERIODIC) THEN
         R (1, M + 1) = 0.25D0
         DO I = 2, N - 3
            R (I, M + 1) = 0.0D0
         END DO
         R (N - 2, M + 1) = 0.25D0
         NRHS = M + 1
      END IF
C
      RETURN
      END

      SUBROUTINE BCD (BC, N, KK, LD, F, RHS, D, G, FRST, NEQS)
C ---------------------------------------------------------------------------
C     Sets the (sub) diagonal coefficients according to the boundary
C     conditions (begin and end) of the spline as requested by the user.
C     The kind of bc. is given by BC:
C     if BC (1) or BC (2) is equal to ...
C
C     0: No derivative information is given by the user, therefore
C        the not-a-knot condition is used instead
C
C     1: First derivative at the edge (xi = 0 or xi = 1)
C        is given by the user
C
C     2: Second derivative at the edge (xi = 0 or xi = 1)
C        is given by the user
C
C     3: First derivative at centre (at xi = 1/2) of first or last
C        element is given by the user
C
C     4: The function to be splined is periodic and a periodicity
C        condition is enforced
C
C     These values are defined in the file spl_params.inc.
C     Any value other than the above will stop the calling program
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
      REAL   (KIND=RK) RHS (N, KK + 1)
      REAL   (KIND=RK) D (N)
      REAL   (KIND=RK) G (LD, N)
      INTEGER(KIND=IK) FRST
      INTEGER(KIND=IK) NEQS
C
      REAL   (KIND=RK) Xi
      PARAMETER       (Xi = 0.5D0)
C
      INTEGER(KIND=IK) I
      REAL   (KIND=RK) Alpha, Beta
C
      FRST = 1
      NEQS = N
C
      IF (BC (1) .EQ. SPL_NOT_A_KNOT) THEN
C        Natural F'' = 0
         D (1) = 2.0
         DO I = 1, KK
            RHS (1, I) = 3.0 * (F (I, 2) - F (I, 1))
         END DO
      ELSE IF (BC (1) .EQ. SPL_FIRST_DERIV) THEN
C        First derivative given at begin
         NEQS = NEQS - 1
         FRST = FRST + 1
         DO I = 1, KK
            RHS (1, I) = G (I, 1)
            RHS (2, I) = RHS (2, I) - 0.25D0 * RHS (1, I)
         END DO
      ELSE IF (BC (1) .EQ. SPL_SECOND_DERIV) THEN
C        Second derivative given at begin
         D (1) = 0.5D0
         DO I = 1, KK
            RHS (1, I) = 0.750D0 * (F (I, 2) - F (I, 1))
     &                 - 0.125D0 *  G (I, 1)
         END DO
      ELSE IF (BC (1) .EQ. SPL_FIRST_D_HALF) THEN
C        First derivative given at half of first element
         Alpha =  4.00D0 * Xi * (3.0D0 * Xi - 2.0D0)
         Beta  =  6.00D0 * Xi * (        Xi - 1.0D0)
         D (1) = (Xi - 1.0D0) * (3.0D0 * Xi - 1.0D0) / Alpha
         DO I = 1, KK
            RHS (1, I) =
     &             (G (I, 1) + Beta * (F (I, 2) - F (I, 1))) / Alpha
         END DO
      ELSE IF (BC (1) .EQ. SPL_PERIODIC) THEN
C        Function is periodic
         DO I = 1, KK
            RHS (1, I) = 0.75D0 * (F (I, 2) - F (I, N - 1))
         END DO
         NEQS = N - 2
      ELSE IF (BC (1) .EQ. SPL_NATURAL) THEN
C        First derivative is assumed to be 0
         NEQS = NEQS - 1
         FRST = FRST + 1
         DO I = 1, KK
            RHS (1, I) = 0.0D0
            RHS (2, I) = RHS (2, I) - 0.25D0 * RHS (1, I)
         END DO
      ELSE
         CALL ERROR ('BCD','Spline type not valid for BC (1)')
      END IF
      IF (BC (2) .EQ. SPL_NOT_A_KNOT) THEN
C        Natural F'' = 0
         D (N) = 2.0D0
         DO I = 1, KK
            RHS (N, I) = -3.0 * (F (I, N - 1) - F (I, N))
         END DO
      ELSE IF (BC (2) .EQ. SPL_FIRST_DERIV) THEN
C        First derivative is given at end
         NEQS = NEQS - 1
         DO I = 1, KK
            RHS (N,     I) = G (I, N)
            RHS (N - 1, I) = RHS (N - 1, I) - 0.25D0 * G (I, N)
         END DO
      ELSE IF (BC (2) .EQ. SPL_SECOND_DERIV) THEN
C        Second derivative is given at end
         D (N) = 0.5D0
         DO I = 1, KK
            RHS (N, I) = 0.750D0 * (F (I, N) - F (I, N - 1))
     &                 + 0.125D0 *  G (I, N)
         END DO
      ELSE IF (BC (2) .EQ. SPL_FIRST_D_HALF) THEN
C        First derivative at midpoint of last element is given
         Alpha =  4.00D0 * (Xi - 1.0D0) * (3.0D0 * Xi - 1.0D0)
         Beta  =  6.00D0 * (Xi - 1.0D0) * Xi
         D (N) =  Xi * (3.0D0 * Xi - 2.0D0) / Alpha
C        D (N) =   0.25D0
         DO I = 1, KK
            RHS (N, I) =
     &         (G (I, N) + Beta * (F (I, N) - F (I, N - 1))) / Alpha
         END DO
      ELSE IF (BC (2) .EQ. SPL_PERIODIC) THEN
C        Function is periodic
      ELSE IF (BC (2) .EQ. SPL_NATURAL) THEN
C        First derivative is assumed to be 0
         NEQS = NEQS - 1
         DO I = 1, KK
            RHS (N,     I) = 0.0D0
            RHS (N - 1, I) = RHS (N - 1, I) - 0.25D0 * RHS (N, I)
         END DO
      ELSE
         CALL ERROR ('BCD','Spline type not valid for BC (2)')
      END IF
C
      RETURN
      END

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
C      DO J = 1, K
C         DO I = 1, N
C            G (J, I) = W (I, J)
C         END DO
C      END DO
      DO J = 1, K
        CALL DCOPY (N, W (1, J), 1, G (J, 1), LD)
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
