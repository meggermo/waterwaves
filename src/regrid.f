      SUBROUTINE REGRID
     &           (NSD, NNW_SD, NGP_SD, NW_IPAR, NW_RPAR, CRD)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      INTEGER(KIND=IK) NSD
      INTEGER(KIND=IK) NNW_SD (*)
      INTEGER(KIND=IK) NGP_SD (*)
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)
      REAL   (KIND=RK) NW_RPAR (N_RP, *)
      REAL   (KIND=RK) CRD (4, *)
C
      INTEGER(KIND=IK) ISD, INW, IGP, NNW, NGP
C
      INW = 1
      IGP = 1
C
      DO ISD = 1, NSD
         NNW = NNW_SD (ISD)
         NGP = NGP_SD (ISD)
         CALL REGRID_SD
     &        (NNW, NW_IPAR (1, INW), NW_RPAR (1, INW), CRD (1, IGP))
         INW = INW + NNW
         IGP = IGP + NGP
      END DO
C
      RETURN
      END
      SUBROUTINE REGRID_SD
     &           (NNW, NW_IPAR, NW_RPAR, CRD)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'net_funcs.inc'
C
      INTEGER(KIND=IK) NNW
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)
      REAL   (KIND=RK) NW_RPAR (N_RP, *)
      REAL   (KIND=RK) CRD (4, *)
C
      INTEGER(KIND=IK) INW, IGP, NGP
C
      IGP = 1
C
      DO INW = 1, NNW
         NGP = GET_NGP (NW_IPAR (1, INW))
         CALL REGRID_NW
     &        (NGP, NW_IPAR (1, INW), NW_RPAR (1, INW), CRD (1, IGP))
         IGP = IGP + NGP
      END DO
C
      RETURN
      END
      SUBROUTINE REGRID_NW
     &           (NGP, NW_IPAR, NW_RPAR, CRD)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'anl_params.inc'
      INCLUDE 'tme_params.inc'
      INCLUDE 'net_funcs.inc'
C
      INTEGER(KIND=IK) NGP
      INTEGER(KIND=IK) NW_IPAR (N_IP)
      REAL   (KIND=RK) NW_RPAR (N_RP)
      REAL   (KIND=RK) CRD (4, *)
C
      INTEGER(KIND=IK) BCT
C
      IF (FIRST_TIME .AND. ANL_TYPE .EQ. ANL_RFWAVE) THEN
         BCT = GET_BCT (NW_IPAR)
         CALL REGRID_RF (BCT, NGP, ANL_PAR, ANL_PAR (67), CRD)
      END IF
C
      RETURN
      END
      SUBROUTINE REGRID_RF (BCT, N, A, B, CRD)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'bct_params.inc'
      INCLUDE 'tme_params.inc'
      INCLUDE 'tme_funcs.inc'
C
      INTEGER(KIND=IK) BCT
      INTEGER(KIND=IK) N
      REAL   (KIND=RK) A (*)
      REAL   (KIND=RK) B (*)
      REAL   (KIND=RK) CRD (4, N)
C
      INTEGER(KIND=IK) I, J, NRF
      REAL   (KIND=RK) L, H, X0, W, K
      REAL   (KIND=RK) X, Z, ETA, ETA_X, T, C
      REAL   (KIND=RK) Z_XI (N)
      INTEGER(KIND=IK) BC (2)
      DATA             BC /0, 0/
C
      T   = GET_TIME ()
      L   = A (2)
      H   = A (3)
      X0  = A (4)
      W   = A (5)
      K   = A (6)
      C   = W / K
      NRF = INT (A (8))
C
      DO I = 1, N
         X     = CRD (1, I) - X0 - C * T
         Z     = CRD (2, I)
         ETA   = 0.5D0 * B (1)
         ETA_X = 0.0D0
         DO J = 1, NRF
           ETA   = ETA    +         B (J + 1) * COS (J * K * X)
           ETA_X = ETA_X  - J * K * B (J + 1) * SIN (J * K * X)
         END DO
         CRD (2, I) = -H + (H + ETA)  * (Z + H) / H
         Z_XI   (I) = CRD (4, I)
      END DO
C
      IF (BCT .EQ. BCT_BERNOULLI) THEN
         DO I = 1, N
            CRD (4, I) = Z_XI (I)
         END DO
      ELSE
         CALL SPLINE (BC, N, 1, 4, CRD (2, 1), CRD (4, 1))
C        WRITE (*, '(1X,4E14.6)') CRD
      END IF
C
      RETURN
      END
