      SUBROUTINE GENERATE_GRID (NNW, NW_IPAR, NW_RPAR, CRD)
C ---------------------------------------------------------------------------
C     Generates the initial grids for the networks.
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'net_funcs.inc'
      INCLUDE 'grd_params.inc'
      INCLUDE 'usr_params.inc'
C
      INTEGER(KIND=IK) NNW               ! nr. of networks (IN)
      INTEGER(KIND=IK) NW_IPAR (N_IP, *) ! integer network params (IN)
      REAL   (KIND=RK) NW_RPAR (N_RP, *) ! real network params (IN)
      REAL   (KIND=RK) CRD (4, *)        ! grid point coordinates (OUT)
C
      INTEGER(KIND=IK) INW, IGP, NGP, GRT
      REAL   (KIND=RK) GRID_PAR (N_GRP)
C
      WRITE (USR_O, *) '========== GENERATING GRID ========'
C     Set the grid point pointer to the first grid point
      IGP = 1
C     Loop over all networks
      DO INW = 1, NNW
         NGP = GET_NGP (NW_IPAR (1, INW))
         GRT = GET_GRT (NW_IPAR (1, INW))
         WRITE (USR_O, 101) INW, NGP
C        Get the number of grid points and the grid type
C        Get the grid parameters
         CALL GET_GRP (NW_RPAR (1, INW), GRID_PAR)
C        Generate the grid depending on the value of GRT
         IF (GRT .EQ. GRD_LINE) THEN
            CALL LINE_GRID
     &           (NGP, GRID_PAR, CRD (1, IGP))
         ELSE IF (GRT .EQ. GRD_STRETCHED) THEN
            CALL STRETCHED_GRID
     &           (NGP, GRID_PAR, CRD (1, IGP))
         ELSE IF (GRT .EQ. GRD_CIRCLESECT) THEN
            CALL CIRCLESECT_GRID
     &           (NGP, GRID_PAR, CRD (1, IGP))
         ELSE IF (GRT .EQ. GRD_CUBIC) THEN
            CALL CUBIC_GRID
     &           (NGP, GRID_PAR, CRD (1, IGP))
         ELSE
            CALL ERROR ('GENERATE_GRID', 'Grid type not implemented')
         END IF
C        increment the grid pointer
         IGP = IGP + NGP
      END DO
      WRITE (USR_O, *) '========== GRID IS GENERATED ======'
C
      RETURN
  101 FORMAT (1X,'NETWORK',I3,', NGP =',I4)
      END

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
