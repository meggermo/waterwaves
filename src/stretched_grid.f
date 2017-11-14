      SUBROUTINE STRETCHED_GRID (N, PAR, X)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'rle_params.inc'
      INCLUDE 'grd_params.inc'
      INCLUDE 'usr_params.inc'
C
      INTEGER(KIND=IK) N
      REAL   (KIND=RK) PAR (2, *)
      REAL   (KIND=RK) X   (4, *)
C
      INTEGER(KIND=IK) I
      REAL   (KIND=RK) XI, S1, C1, S, DS
      REAL   (KIND=RK) K
      REAL   (KIND=RK) A
      REAL   (KIND=RK) DX (2)
C
      DX (1) = PAR (1, 2) - PAR (1, 1)
      DX (2) = PAR (2, 2) - PAR (2, 1)
C
      K = PAR (2, 3)
      A = PAR (1, 3) / K / PI
C
      WRITE (USR_O, 1001)
      WRITE (USR_O, 1101) 'from', PAR (1, 1), PAR (2, 1)
      WRITE (USR_O, 1101) 'to',   PAR (1, 2), PAR (2, 2)
      WRITE (USR_O, 1201) 'Amplitude:',       PAR (1, 3)
      WRITE (USR_O, 1201) 'Periodicity:',     PAR (2, 3)
C
      DO I = 1, N
         XI = DBLE (I -1) / DBLE (N - 1)
         S1 = A * SIN (PI * K * XI)
         C1 = A * COS (PI * K * XI)
         S  = XI - S1
         DS = (1.0D0 - PI * K * C1) / DBLE (N - 1)
         X (1, I) = PAR (1, 1) + DX (1) * S
         X (2, I) = PAR (2, 1) + DX (2) * S
         X (3, I) = DX (1) * DS
         X (4, I) = DX (2) * DS
      END DO
C
      RETURN
 1001 FORMAT (1X,'Generating stretched grid:')
 1101 FORMAT (1X,A4, 2E14.6)
 1201 FORMAT (1X,A12, F6.2)
      END