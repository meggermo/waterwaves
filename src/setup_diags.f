      SUBROUTINE SETUP_DIAGS (N, D)
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
      REAL   (KIND=RK) D (N, 3)
C
      INTEGER(KIND=IK)  I
C
      DO I = 1, N
          D (I, 1) = 1.00D0
          D (I, 2) = 0.25D0
      END DO
C
      RETURN
      END
