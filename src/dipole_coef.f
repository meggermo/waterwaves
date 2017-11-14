      FUNCTION DIPOLE_COEF (S)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'spl_params.inc'
C
      REAL(KIND=RK) DIPOLE_COEF
      REAL(KIND=RK) S
C
      REAL(KIND=RK) W (4, 2)
      REAL(KIND=RK) R (2)
      REAL(KIND=RK) DOT_4
      EXTERNAL      DOT_4
C
      CALL WEIGHT (S, W)
      R (1) = P (1) - DOT_4 (W, EL (1, 1))
      R (2) = P (2) - DOT_4 (W, EL (1, 2))
      DIPOLE_COEF =
     &   W (K, 1) * (R (2) * DOT_4 (W (1, 2), EL (1, 1))
     &             - R (1) * DOT_4 (W (1, 2), EL (1, 2)))
     &             / (R (1) ** 2 + R (2) ** 2)
C
      RETURN
      END