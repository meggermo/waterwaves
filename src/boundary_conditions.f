      SUBROUTINE BOUNDARY_CONDITIONS
     &           (BCD, NSD, NNW_SD, NW_IPAR, NW_RPAR, CRD, PHI, PHN)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'anl_params.inc'
      INCLUDE 'bct_params.inc'
      INCLUDE 'tme_params.inc'
      INCLUDE 'net_funcs.inc'
      INCLUDE 'tme_funcs.inc'
C
      INTEGER(KIND=IK) BCD
      INTEGER(KIND=IK) NSD               ! nr. of subdomains (IN)
      INTEGER(KIND=IK) NNW_SD (*)        ! nr. of netw. per subd. (IN)
      INTEGER(KIND=IK) NW_IPAR (N_IP, *) ! network integer params. (IN)
      REAL   (KIND=RK) NW_RPAR (N_RP, *) ! network real    params. (IN)
      REAL   (KIND=RK) CRD (4, *)        ! grid point coords. (IN)
      REAL   (KIND=RK) PHI (2, *)        ! PHI   (OUT)
      REAL   (KIND=RK) PHN (2, *)        ! PHI_N (OUT)
C
      INTEGER(KIND=IK)  ISD, INW, JNW, IGP, NGP, BCT, ITF, GMT
      INTEGER(KIND=IK)  ANL_SAVE
      INTEGER(KIND=IK)  BC (2)
      DATA              BC /0, 0/
      REAL   (KIND=RK)  T
C
      T = GET_TIME ()
C
      IGP = 1
      INW = 1
      DO ISD = 1, NSD
         DO JNW = INW, INW + NNW_SD (ISD) - 1
            NGP = GET_NGP (NW_IPAR (1, JNW))
            BCT = GET_BCT (NW_IPAR (1, JNW))
            ITF = GET_ITF (NW_IPAR (1, JNW))
            GMT = GET_GMT (NW_IPAR (1, JNW))
            IF (BCT .EQ. BCT_WAVEMAKER) THEN
               ANL_SAVE = ANL_TYPE
               ANL_TYPE = ANL_WAVEMAKER
               CALL BOUND_COND
     &              (BCD, NGP, BCT, T, NW_RPAR (1, JNW),
     &               CRD (1, IGP), PHI (1, IGP), PHN (1, IGP))
               CALL Spline (BC, NGP, 1, 2, PHN (1, IGP), PHN (2, IGP))
               ANL_TYPE = ANL_SAVE
            ELSE IF  (ABS (BCT) .EQ. BCT_ANALYTIC .OR. FIRST_TIME)  THEN
               CALL BOUND_COND
     &              (BCD, NGP, BCT, T, ANL_PAR,
     &               CRD (1, IGP), PHI (1, IGP), PHN (1, IGP))
            END IF
            IGP = IGP + NGP
         END DO
         INW = INW + NNW_SD (ISD)
      END DO
C
      RETURN
      END

      SUBROUTINE BOUND_COND (BCD, NGP, BCT, T, PAR, CRD, PHI, PHN)
C ---------------------------------------------------------------------------
C     Computes the analytical solution for PHI if BCT < 0 and for PHI_N
C     if BCT > 0
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'bct_params.inc'
      INCLUDE 'tme_params.inc'
      INCLUDE 'rle_params.inc'
      INCLUDE 'chk_params.inc'
C
      INTEGER(KIND=IK) BCD         ! bc's for Laplace (PHI) or for PHI_T
      INTEGER(KIND=IK) NGP         ! nr. of grid points of the netw. (IN)
      INTEGER(KIND=IK) BCT         ! type of boundary condition (IN)
      REAL   (KIND=RK) T           ! time (IN)
      REAL   (KIND=RK) PAR (*)     ! the params. of anal. solution
      REAL   (KIND=RK) CRD (4, *)  ! grid points. of netw. (IN)
      REAL   (KIND=RK) PHI (2, *)  ! PHI   (OUT)
      REAL   (KIND=RK) PHN (2, *)  ! PHI_N (OUT)
C
      INTEGER(KIND=IK) I
      REAL   (KIND=RK) X_XI, Z_XI, J_INV, PROD, DIFF
      REAL   (KIND=RK) W (NGP, 7)
C     Get the analytic solution at the grid points at time = T and store
C     them in W
      IF (BCD .EQ. BCD_PHI) THEN
C        Get bc's for problem of type Laplace (Phi) = 0
         CALL ANALYTIC (NGP, T, PAR, CRD, W)
      ELSE IF (BCD .EQ. BCD_PHI_T) THEN
C        Get bc's for problem of type Laplace (Phi_t) = 0
         CALL ANALYTIC_T (NGP, T, PAR, CRD, W)
      END IF
C
      DO I = 1, NGP
         X_XI  = CRD (3, I)
         Z_XI  = CRD (4, I)
         J_INV = 1.0D0 / SQRT (X_XI ** 2 + Z_XI ** 2)
         PROD  = 2.0D0 * X_XI * Z_XI
         DIFF  = (X_XI - Z_XI) * (X_XI + Z_XI)
         IF (BCT .LT. 0 .OR. CHK_COEFS_PRE) THEN
C           A Dirichlet type of boundary, so compute phi and phi_xi
            PHI (1, I) =  W (I, 1)
            PHI (2, I) =  W (I, 3) * X_XI + W (I, 4) * Z_XI
            IF (BCT. EQ. -BCT_LINBERNOU) THEN
C              eta = -phi_t / g
               CRD (2, I) = -W (I, 2) / GRAV
            END IF
         END IF
         IF (BCT .GT. 0 .OR. CHK_COEFS_PRE) THEN
C           A Neumann type of boundary, so compute phi_n and (phi_n)_xi
            PHN (1, I) = (W (I, 4) * X_XI - W (I, 3) * Z_XI) * J_INV
            PHN (2, I) = (W (I, 6) * DIFF - W (I, 5) * PROD) * J_INV
         END IF
      END DO
C
      RETURN
      END
