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
