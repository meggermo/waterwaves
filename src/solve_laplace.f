      SUBROUTINE SOLVE_LAPLACE
     &     (BCD,     NSD,
     &      NNW_SD,  NGP_SD,
     &      NW_IPAR, NW_RPAR,
     &      CRD,     PHI,     PHN,
     &      S0,      S1,      D0,      D1)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'tme_params.inc'
      INCLUDE 'chk_params.inc'
      INCLUDE 'bct_params.inc'
C
      INTEGER(KIND=IK) BCD
      INTEGER(KIND=IK) NSD
      INTEGER(KIND=IK) NNW_SD (*)
      INTEGER(KIND=IK) NGP_SD (*)
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)
      REAL   (KIND=RK) NW_RPAR (N_RP, *)
      REAL   (KIND=RK) CRD (4, *)
      REAL   (KIND=RK) PHI (2, *)
      REAL   (KIND=RK) PHN (2, *)
      REAL   (KIND=RK) S0 (*)
      REAL   (KIND=RK) S1 (*)
      REAL   (KIND=RK) D0 (*)
      REAL   (KIND=RK) D1 (*)
C
      IF ((GRID_TIME_DEPENDENT .OR. FIRST_TIME)
     &    .AND. BCD .NE. BCD_PHI_T) THEN
C        Depending on the type of boundary condition some grids must be
C        recomputed (e.g. R&F bc's)
         CALL REGRID
     &        (NSD, NNW_SD, NGP_SD, NW_IPAR, NW_RPAR, CRD)
C        Compute the source and dipole coefficients only if the grid is
C        time-dependent or if this is the first time
         CALL COEFFICIENTS
     &        (NSD, NNW_SD, NGP_SD, NW_IPAR, CRD, S0, S1, D0, D1)
      END IF
C     Compute the boundary conditions for those networks that use an
C     analytic solution.
      CALL BOUNDARY_CONDITIONS
     &     (BCD, NSD, NNW_SD, NW_IPAR, NW_RPAR, CRD, PHI, PHN)
C     Check if the solution satisfies the BIE
      IF (CHK_COEFS_PRE) CALL CHECK_COEFS
     &     (NSD, NGP_SD, PHI, PHN, S0, S1, D0, D1)
C     Solve the Boundary Integral equations
      CALL SOLVE_BIE
     &     (NSD,
     &      NNW_SD, NGP_SD,  NW_IPAR,
     &      S0,     S1,
     &      D0,     D1,
     &      CRD,    PHI,     PHN)
C     Check if the solution satisfies the BIE
      IF (CHK_COEFS_POST) CALL CHECK_COEFS
     &     (NSD, NGP_SD, PHI, PHN, S0, S1, D0, D1)
C     The first boundary value problem has been solved
      FIRST_TIME = .FALSE.
C
      RETURN
      END
