      SUBROUTINE CIRCLESECT_GRID (N, PAR, X)
C ---------------------------------------------------------------------------
C
C     PAR (1, 1) - PAR (2, 1) : circle centre coordinates
C     PAR (2, 2) - PAR (2, 2) : circle radius in x- and z-direction
C     PAR (1, 3) - PAR (2, 3) : start- and end angle
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'rle_params.inc'
C
      INTEGER(KIND=IK) N
      REAL   (KIND=RK) PAR (2, *)
      REAL   (KIND=RK) X   (4, *)
C
      INTEGER(KIND=IK) I
      REAL   (KIND=RK) THETA, DTHETA
C
      WRITE (USR_O, 1001)
      WRITE (USR_O, 1101) 'centre', PAR (1, 1), PAR (2, 1)
      WRITE (USR_O, 1101) 'radii',  PAR (1, 2), PAR (2, 2)
      WRITE (USR_O, 1201) 'theta', 2.0D0 * PAR (1, 3),
     &                             2.0D0 * PAR (2, 3)
C
      DTHETA = 2.0D0 * PI * (PAR (2, 3) - PAR (1, 3)) / DBLE (N - 1)
C
      DO I = 1, N
         THETA = 2.0D0 * PI * PAR (1, 3) + DTHETA * DBLE (I - 1)
         X (1, I) = PAR (1, 1) + PAR (1, 2) * COS (THETA)
         X (2, I) = PAR (2, 1) + PAR (2, 2) * SIN (THETA)
         X (3, I) =    -DTHETA * PAR (1, 2) * SIN (THETA)
         X (4, I) =     DTHETA * PAR (2, 2) * COS (THETA)
      END DO
C
      RETURN
 1001 FORMAT (1X,'Generating circular grid:')
 1101 FORMAT (1X,A7,2E14.6)
 1201 FORMAT (1X,A7,F9.3,' * PI', F9.3, ' * PI')
      END