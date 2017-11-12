C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      INTEGER(KIND=IK) SPL_NOT_A_KNOT
      INTEGER(KIND=IK) SPL_FIRST_DERIV
      INTEGER(KIND=IK) SPL_SECOND_DERIV
      INTEGER(KIND=IK) SPL_FIRST_D_HALF
      INTEGER(KIND=IK) SPL_PERIODIC
      INTEGER(KIND=IK) SPL_NATURAL
C
      PARAMETER       (SPL_NOT_A_KNOT   = 0)
      PARAMETER       (SPL_FIRST_DERIV  = 1)
      PARAMETER       (SPL_SECOND_DERIV = 2)
      PARAMETER       (SPL_FIRST_D_HALF = 3)
      PARAMETER       (SPL_PERIODIC     = 4)
      PARAMETER       (SPL_NATURAL      = 5)
C
      INTEGER(KIND=IK) K
      REAL   (KIND=RK) P  (2)
      REAL   (KIND=RK) EL (4, 3)
C
      COMMON /SPLINE_PARAMS/ EL, P, K
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
