      SUBROUTINE CUBIC_GRID (N, PAR, X)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'rle_params.inc'
      INTEGER(KIND=IK) N
      REAL(KIND=RK) PAR (2, *)
      REAL(KIND=RK) X   (4, N)
C
      INTEGER(KIND=IK) I
      REAL(KIND=RK) XI, S, DS
      REAL(KIND=RK) A
      REAL(KIND=RK) DX (2)
C
      DX (1) = PAR (1, 2) - PAR (1, 1)
      DX (2) = PAR (2, 2) - PAR (2, 1)
C
      A = PAR (1, 3)
      IF (A .GE. 0.0D0) THEN
         A = 1.0D0 + 0.5D0 * A
      ELSE
         A = 1.0D0 + A
      END IF
C
      DO I = 1, N
         XI = DBLE (I -1) / DBLE (N - 1)
         S  = (3.0D0 - 2.0D0 * A) * XI
     &      + 6.0D0 * (A - 1.0D0) * XI ** 2
     &      - 4.0D0 * (A - 1.0D0) * XI ** 3
         DS = (3.0D0 - 2.0D0 * A)
     &      + 12.0D0 * (A - 1.0D0) * XI
     &      - 12.0D0 * (A - 1.0D0) * XI ** 2
         X (1, I) = PAR (1, 1) + DX (1) * S
         X (2, I) = PAR (2, 1) + DX (2) * S
         X (3, I) = DX (1) * DS / DBLE (N - 1)
         X (4, I) = DX (2) * DS / DBLE (N - 1)
      END DO
C
      RETURN
      END