      FUNCTION DOT_4 (X, Y)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      REAL(KIND=RK) DOT_4
      REAL(KIND=RK) X (*)
      REAL(KIND=RK) Y (*)
C
      DOT_4 = X (1) * Y (1)
     &      + X (2) * Y (2)
     &      + X (3) * Y (3)
     &      + X (4) * Y (4)
C
      RETURN
      END