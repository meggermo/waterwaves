      SUBROUTINE WRITE_COEFS_SD
     &           (ISD, NGP, S0, S1, D0, D1)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      INTEGER(KIND=IK) ISD           ! subdomain id.  (IN)
      INTEGER(KIND=IK) NGP           ! nr. of gr. p.  (IN)
      REAL   (KIND=RK) S0 (*)        ! Source coeffs. (IN)
      REAL   (KIND=RK) S1 (*)        ! Source coeffs. (IN)
      REAL   (KIND=RK) D0 (*)        ! Dipole coeffs. (IN)
      REAL   (KIND=RK) D1 (*)        ! Dipole coeffs. (IN)
C
      CHARACTER*9 FILE_NAME
C
      WRITE (FILE_NAME, 101) 'S0', ISD
      CALL WRITE_MATRIX (FILE_NAME, NGP, NGP, S0)
C
      WRITE (FILE_NAME, 101) 'S1', ISD
      CALL WRITE_MATRIX (FILE_NAME, NGP, NGP, S1)
C
      WRITE (FILE_NAME, 101) 'D0', ISD
      CALL WRITE_MATRIX (FILE_NAME, NGP, NGP, D0)
C
      WRITE (FILE_NAME, 101) 'D1', ISD
      CALL WRITE_MATRIX (FILE_NAME, NGP, NGP, D1)
C
      RETURN
  101 FORMAT (A,'_',I1,'.out')
      END
       SUBROUTINE WRITE_MATRIX
     &           (FILE_NAME, N, M, C)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'fle_params.inc'
C
      CHARACTER*(*)    FILE_NAME     ! file name to open (IN)
      INTEGER(KIND=IK) N             ! nr. rows          (IN)
      INTEGER(KIND=IK) M             ! nr. coulumns      (IN)
      REAL   (KIND=RK) C (N, *)      ! Some matrix       (IN)
C
      INTEGER(KIND=IK) I, J
C
      OPEN (UNIT = STD_T, FILE = FILE_NAME, STATUS = 'UNKNOWN')
      DO J = 1, M
         DO I = 1, N
            WRITE (STD_T, 101) C (J, I)
         END DO
         WRITE (STD_T, *)
      END DO
      CLOSE (UNIT = STD_T)
C
      RETURN
  101 FORMAT (E16.8,$)
      END
      SUBROUTINE WRITE_VECTOR
     &           (FILE_NAME, N, V)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'fle_params.inc'
C
      CHARACTER*(*)    FILE_NAME ! file name to open (IN)
      INTEGER(KIND=IK) N         ! nr. elements      (IN)
      REAL   (KIND=RK) V (*)     ! Some vector       (IN)
C
      INTEGER(KIND=IK) I
C
      OPEN (UNIT = STD_T, FILE = FILE_NAME, STATUS = 'UNKNOWN')
      DO I = 1, N
         WRITE (STD_T, 101) V (I)
      END DO
      WRITE (STD_T, *)
      CLOSE (UNIT = STD_T)
C
      RETURN
  101 FORMAT (E16.8,$)
      END
