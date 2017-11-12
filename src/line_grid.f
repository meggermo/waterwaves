      SUBROUTINE LINE_GRID (N, PAR, X)
C ---------------------------------------------------------------------------
C     Generates the equidistant coordinates of a line defined by:
C      
C     PAR (1, 1) - PAR (2, 1) : Start coordinate
C     PAR (1, 2) - PAR (2, 2) : End coordinate
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'grd_params.inc'
      INCLUDE 'usr_params.inc'
C
      INTEGER(KIND=IK) N                ! nr. of grid points (IN)
      REAL   (KIND=RK) PAR (2, *)       ! line begin/end point see above (IN)
      REAL   (KIND=RK) X   (4, *)       ! coordinates of the line (OUT)
C
      INTEGER(KIND=IK) I
      REAL   (KIND=RK) DX (2)
C
      WRITE (USR_O, 1001) 
      WRITE (USR_O, 1101) 'from', PAR (1, 1), PAR (2, 1)
      WRITE (USR_O, 1101) 'to',   PAR (1, 2), PAR (2, 2)
C
      DX (1) = (PAR (1, 2) - PAR (1, 1)) / DBLE (N - 1)
      DX (2) = (PAR (2, 2) - PAR (2, 1)) / DBLE (N - 1)
C
      DO I = 1, N
         X (1, I) = PAR (1, 1) + DX (1) * DBLE (I - 1)
         X (2, I) = PAR (2, 1) + DX (2) * DBLE (I - 1)
         X (3, I) = DX (1)
         X (4, I) = DX (2)
      END DO
C
      RETURN
 1001 FORMAT (1X,'Generating line grid:')
 1101 FORMAT (1X,A4,2E14.6)
      END
