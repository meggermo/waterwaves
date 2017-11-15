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
