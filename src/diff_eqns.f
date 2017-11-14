      SUBROUTINE DIFF_EQNS
     &           (STAGE,   NSD,
     &            NNW_SD,  NGP_SD,
     &            NW_IPAR, NW_RPAR,
     &            CRD,     PHI,     PHN,
     &            DCRD_DT, DPHI_DT, DPHN_DT,
     &            S0,      S1,      D0,      D1)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'bct_params.inc'
      INCLUDE 'chk_params.inc'
      INCLUDE 'tme_params.inc'
      INCLUDE 'tme_funcs.inc'
C
      INTEGER(KIND=IK) STAGE
      INTEGER(KIND=IK) NSD
      INTEGER(KIND=IK) NNW_SD (*)
      INTEGER(KIND=IK) NGP_SD (*)
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)
      REAL   (KIND=RK) NW_RPAR (N_RP, *)
      REAL   (KIND=RK) CRD (4, *)
      REAL   (KIND=RK) PHI (2, *)
      REAL   (KIND=RK) PHN (2, *)
      REAL   (KIND=RK) DCRD_DT (4, *)
      REAL   (KIND=RK) DPHI_DT (2, *)
      REAL   (KIND=RK) DPHN_DT (2, *)
      REAL   (KIND=RK) S0 (*)
      REAL   (KIND=RK) S1 (*)
      REAL   (KIND=RK) D0 (*)
      REAL   (KIND=RK) D1 (*)
C     Solve Delta(Phi) = 0 (by means of Boundary Integral Equations)
      CALL SOLVE_LAPLACE
     &     (BCD_PHI, NSD,
     &      NNW_SD,  NGP_SD,
     &      NW_IPAR, NW_RPAR,
     &      CRD,     PHI,     PHN,
     &      S0,      S1,      D0,      D1)
C
      IF (PROB_TIME_DEPENDENT .AND. GET_TIME () .LE. GET_TEND ()) THEN
C        Compute the righthand sides the the time-dependent DE's
         CALL RIGHTHAND_SIDES
     &        (STAGE,   NSD,
     &         NNW_SD,  NGP_SD,
     &         NW_IPAR, NW_RPAR,
     &         CRD,     PHI,     PHN,
     &         DCRD_DT, DPHI_DT, DPHN_DT)
C        Solve Laplace equation for Phi_t
C        CALL SOLVE_LAPLACE
C    &     (BCD_PHI_T, NSD,
C    &      NNW_SD,    NGP_SD,
C    &      NW_IPAR,   NW_RPAR,
C    &      CRD,       DPHI_DT,   DPHN_DT,
C    &      S0,        S1,        D0,        D1)
      END IF
C
      IF (STAGE .EQ. 1 .AND. I_PLOT .EQ. 0) THEN
         CALL PLTOUT
     &        (NSD,     NNW_SD,  NW_IPAR,
     &         CRD,      PHI,     PHN,
     &         DCRD_DT, DPHI_DT, DPHN_DT)
         IF (CHK_CONTOUR) THEN
            CALL CONTOUR_INTEGRAL
     &           (NSD, NNW_SD, NGP_SD, NW_IPAR, CRD, PHN)
         END IF
      END IF
C
      RETURN
      END
      SUBROUTINE RIGHTHAND_SIDES
     &           (STAGE,   NSD,
     &            NNW_SD,  NGP_SD,
     &            NW_IPAR, NW_RPAR,
     &            CRD,     PHI,     PHN,
     &            DCRD_DT, DPHI_DT, DPHN_DT)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      INTEGER(KIND=IK) STAGE
      INTEGER(KIND=IK) NSD
      INTEGER(KIND=IK) NNW_SD (*)
      INTEGER(KIND=IK) NGP_SD (*)
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)
      REAL   (KIND=RK) NW_RPAR (N_RP, *)
      REAL   (KIND=RK) CRD (4, *)
      REAL   (KIND=RK) PHI (2, *)
      REAL   (KIND=RK) PHN (2, *)
      REAL   (KIND=RK) DCRD_DT (4, *)
      REAL   (KIND=RK) DPHI_DT (2, *)
      REAL   (KIND=RK) DPHN_DT (2, *)
C
      INTEGER(KIND=IK) ISD, INW, IGP, NNW, NGP
C
      INW = 1
      IGP = 1
C
      DO ISD = 1, NSD
         NNW = NNW_SD (ISD)
         NGP = NGP_SD (ISD)
         CALL RHS_SD
     &        (STAGE,             NNW,              NGP,
     &         NW_IPAR (1, INW),  NW_RPAR (1, INW),
     &         CRD     (1, IGP),  PHI     (1, IGP), PHN     (1, IGP),
     &         DCRD_DT (1, IGP),  DPHI_DT (1, IGP), DPHN_DT (1, IGP))
         INW = INW + NNW
         IGP = IGP + NGP
      END DO
C
      RETURN
      END
      SUBROUTINE RHS_SD
     &           (STAGE,   NNW,     NGP,
     &            NW_IPAR, NW_RPAR,
     &            CRD,     PHI,     PHN,
     &            DCRD_DT, DPHI_DT, DPHN_DT)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'tme_params.inc'
C
      INTEGER(KIND=IK) STAGE
      INTEGER(KIND=IK) NNW
      INTEGER(KIND=IK) NGP
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)
      REAL   (KIND=RK) NW_RPAR (N_RP, *)
      REAL   (KIND=RK) CRD (4, *)
      REAL   (KIND=RK) PHI (2, *)
      REAL   (KIND=RK) PHN (2, *)
      REAL   (KIND=RK) DCRD_DT (4, *)
      REAL   (KIND=RK) DPHI_DT (2, *)
      REAL   (KIND=RK) DPHN_DT (2, *)
C
      CALL PARTIAL_DERIVS
     &     (STAGE,   NNW,
     &      NW_IPAR, NW_RPAR,
     &      CRD,     PHI,     PHN,
     &      DCRD_DT, DPHI_DT, DPHN_DT)
C
      IF (GRID_TIME_DEPENDENT) THEN
         CALL MATERIAL_DERIVS
     &        (STAGE,   NNW,     NGP,
     &         NW_IPAR, NW_RPAR,
     &         CRD,     PHI,     PHN,
     &         DCRD_DT, DPHI_DT, DPHN_DT)
      END IF
C
      RETURN
      END
      SUBROUTINE PARTIAL_DERIVS
     &           (STAGE,   NNW,
     &            NW_IPAR, NW_RPAR,
     &            CRD,     PHI,     PHN,
     &            CRD_T,   PHI_T,   PHN_T)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'bct_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'net_funcs.inc'
C
      INTEGER(KIND=IK) STAGE             ! Current RK stage (IN)
      INTEGER(KIND=IK) NNW               ! nr. of networks (IN)
      INTEGER(KIND=IK) NW_IPAR (N_IP, *) ! integer netw. params. (IN)
      REAL   (KIND=RK) NW_RPAR (N_RP, *) ! real netw. params. (IN)
      REAL   (KIND=RK) CRD (4, *)        ! grid point coords. (IN)
      REAL   (KIND=RK) PHI (2, *)        ! potential (IN)
      REAL   (KIND=RK) PHN (2, *)        ! normal deriv of Phi (IN)
      REAL   (KIND=RK) CRD_T (4, *)      ! partial d/dT of CRD (OUT)
      REAL   (KIND=RK) PHI_T (2, *)      ! partial d/dT of PHI (OUT)
      REAL   (KIND=RK) PHN_T (2, *)      ! partial d/dT of PHN (OUT)
C
      INTEGER(KIND=IK) INW, IGP, JGP, NGP, BCT
      REAL   (KIND=RK) X_XI, Z_XI, J_I, U_L (2), U_G (2)
C
      IGP = 1
C
      DO INW = 1, NNW
         NGP = GET_NGP (NW_IPAR (1, INW))
         BCT = GET_BCT (NW_IPAR (1, INW))
         DO JGP = IGP, IGP + NGP - 1
            CRD_T (1, JGP) = 0.0D0
            CRD_T (2, JGP) = 0.0D0
            PHI_T (1, JGP) = 0.0D0
            PHN_T (1, JGP) = 0.0D0
         END DO
         IF (BCT .EQ. -BCT_BERNOULLI) THEN
            CALL BERNOULLI
     &           (NGP,
     &            NW_IPAR (1, INW), NW_RPAR (1, INW),
     &            CRD     (1, IGP), PHI     (1, IGP), PHN (1, IGP),
     &            CRD_T   (1, IGP), PHI_T   (1, IGP))
         ELSE IF (BCT .EQ. -BCT_LINBERNOU) THEN
            CALL LINBERNOU
     &           (NGP,
     &            NW_IPAR (1, INW), NW_RPAR (1, INW),
     &            CRD     (1, IGP), PHI     (1, IGP), PHN (1, IGP),
     &            CRD_T   (1, IGP), PHI_T   (1, IGP))
         ELSE IF (BCT .EQ. -BCT_SOMMERFELD) THEN
            CALL SOMMERFELD
     &           (NGP,
     &            NW_IPAR (1, INW), NW_RPAR (1, INW),
     &            CRD     (1, IGP), PHI     (1, IGP), PHN (1, IGP),
     &            CRD_T   (1, IGP), PHI_T   (1, IGP))
         ELSE IF (BCT .EQ. -BCT_TIME_POLY) THEN
            CALL TIME_POLY
     &           (NGP,
     &            NW_IPAR (1, INW), NW_RPAR (1, INW),
     &            CRD     (1, IGP), PHI     (1, IGP), PHN (1, IGP),
     &            CRD_T   (1, IGP), PHI_T   (1, IGP))
         END IF
         DO JGP = IGP, IGP + NGP - 1
            X_XI = CRD (3, JGP)
            Z_XI = CRD (4, JGP)
            J_I  = 1.0D0 / SQRT (X_XI ** 2 + Z_XI ** 2)
            U_L (1) = PHI (2, JGP) * J_I
            U_L (2) = PHN (1, JGP)
            U_G (1) = (X_XI * U_L (1) - Z_XI * U_L (2)) * J_I
            U_G (2) = (X_XI * U_L (2) + Z_XI * U_L (1)) * J_I
            CALL LOCAL_TO_GLOBAL (1, 2, CRD (1, JGP), U_L, U_G)
C           WRITE (*, 1001) JGP, PHI_T (1, JGP), U_G
         END DO
         IGP = IGP + NGP
      END DO
C
      RETURN
      END
      SUBROUTINE BERNOULLI
     &           (NGP,
     &            NW_IPAR, NW_RPAR,
     &            CRD,     PHI,     PHN,
     &            CRD_T,   PHI_T)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'rle_params.inc'
C
      INTEGER(KIND=IK) NGP
      INTEGER(KIND=IK) NW_IPAR (*)
      REAL   (KIND=RK) NW_RPAR (*)
      REAL   (KIND=RK) CRD (4, *)
      REAL   (KIND=RK) PHI (2, *)
      REAL   (KIND=RK) PHN (2, *)
      REAL   (KIND=RK) CRD_T (4, *)
      REAL   (KIND=RK) PHI_T (2, *)
C
      INTEGER(KIND=IK) IGP
      REAL   (KIND=RK) X_XI, Z_XI, J_I, U_L (2), Q
C
      DO IGP = 1, NGP
         X_XI = CRD (3, IGP)
         Z_XI = CRD (4, IGP)
         J_I  = 1.0D0 / SQRT (X_XI ** 2 + Z_XI ** 2)
         U_L (1) = PHI (2, IGP) * J_I
         U_L (2) = PHN (1, IGP)
         Q       = U_L (1) ** 2 + U_L (2) ** 2
C        phi_t + 1/2 grad (phi) ** 2 + g z = 0
         PHI_T (1, IGP) =  -0.5D0 * Q - GRAV * CRD (2, IGP)
      END DO
C
      RETURN
      END
      SUBROUTINE LINBERNOU
     &           (NGP,
     &            NW_IPAR, NW_RPAR,
     &            CRD,     PHI,     PHN,
     &            CRD_T,   PHI_T)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'rle_params.inc'
      INCLUDE 'tme_funcs.inc'
C
      INTEGER(KIND=IK) NGP
      INTEGER(KIND=IK) NW_IPAR (*)
      REAL   (KIND=RK) NW_RPAR (*)
      REAL   (KIND=RK) CRD (4, *)
      REAL   (KIND=RK) PHI (2, *)
      REAL   (KIND=RK) PHN (2, *)
      REAL   (KIND=RK) CRD_T (4, *)
      REAL   (KIND=RK) PHI_T (2, *)
C
      INTEGER(KIND=IK) IGP
      REAL   (KIND=RK) X_XI, Z_XI, J_I, U_L (2), U_G (2)
      REAL   (KIND=RK) ETA (2, NGP)
      INTEGER(KIND=IK) BC (2)
      DATA             BC /0, 0/
C
      CALL DCOPY (NGP, CRD (2, 1), 4, ETA (1, 1), 2)
      CALL Spline (BC, NGP, 1, 2, ETA (1, 1), ETA (2, 1))
      DO IGP = 1, NGP
         X_XI = CRD (3, IGP)
         Z_XI = CRD (4, IGP)
         J_I  = 1.0D0 / SQRT (X_XI ** 2 + Z_XI ** 2)
         U_L (1) = PHI (2, IGP) * J_I
         U_L (2) = PHN (1, IGP)
         U_G (1) = (X_XI * U_L (1) - Z_XI * U_L (2)) * J_I
         U_G (2) = (X_XI * U_L (2) + Z_XI * U_L (1)) * J_I
C        phi_t = -g * eta
         PHI_T (1, IGP) =  -GRAV * CRD (2, IGP)
C        phi_xi_t = -g * eta
         PHI_T (2, IGP) =  -GRAV * ETA (2, IGP)
C        eta_t =  Phi_z
         CRD_T (2, IGP) = U_G (2)
C        WRITE (*, *) IGP, ':', PHI_T (1, IGP), PHI_T (2, IGP)
      END DO
C
      RETURN
      END
      SUBROUTINE SOMMERFELD
     &           (NGP,
     &            NW_IPAR, NW_RPAR,
     &            CRD,     PHI,     PHN,
     &            CRD_T,   PHI_T)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      INTEGER(KIND=IK) NGP
      INTEGER(KIND=IK) NW_IPAR (*)
      REAL   (KIND=RK) NW_RPAR (*)
      REAL   (KIND=RK) CRD (4, *)
      REAL   (KIND=RK) PHI (2, *)
      REAL   (KIND=RK) PHN (2, *)
      REAL   (KIND=RK) CRD_T (4, *)
      REAL   (KIND=RK) PHI_T (2, *)
C
      INTEGER(KIND=IK) IGP
      REAL   (KIND=RK) PAR (N_BCP), C, B
C
      CALL GET_BCP (NW_RPAR, PAR)
C
      C = PAR (1)
      B = PAR (2)
C     B = g * (R- H) - 1/2 * C ** 2 (see page 48 of PdH's thesis)
      DO IGP = 1, NGP
         PHI_T (1, IGP) =  C * PHN (1, IGP) - B
      END DO
C
      RETURN
      END
      SUBROUTINE TIME_POLY
     &           (NGP,
     &            NW_IPAR, NW_RPAR,
     &            CRD,     PHI,     PHN,
     &            CRD_T,   PHI_T)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'anl_params.inc'
      INCLUDE 'tme_funcs.inc'
C
      INTEGER(KIND=IK) NGP
      INTEGER(KIND=IK) NW_IPAR (*)
      REAL   (KIND=RK) NW_RPAR (*)
      REAL   (KIND=RK) CRD (4, *)
      REAL   (KIND=RK) PHI (2, *)
      REAL   (KIND=RK) PHN (2, *)
      REAL   (KIND=RK) CRD_T (4, *)
      REAL   (KIND=RK) PHI_T (2, *)
C
      INTEGER(KIND=IK) IGP
      REAL   (KIND=RK) C, T
C
      T = GET_TIME ()
      C = ANL_PAR (5)
      DO IGP = 1, NGP
         PHI_T (1, IGP) = 3.0D0 * C * T ** 2
      END DO
C
      RETURN
      END
      SUBROUTINE MATERIAL_DERIVS
     &           (STAGE,   NNW,     NGP,
     &            NW_IPAR, NW_RPAR,
     &            CRD,     PHI,     PHN,
     &            DCRD_DT, DPHI_DT, DPHN_DT)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      INTEGER(KIND=IK) STAGE
      INTEGER(KIND=IK) NNW
      INTEGER(KIND=IK) NGP
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)
      REAL   (KIND=RK) NW_RPAR (N_RP, *)
      REAL   (KIND=RK) CRD (4, *)
      REAL   (KIND=RK) PHI (2, *)
      REAL   (KIND=RK) PHN (2, *)
      REAL   (KIND=RK) DCRD_DT (4, *)
      REAL   (KIND=RK) DPHI_DT (2, *)
      REAL   (KIND=RK) DPHN_DT (2, *)
C
      CALL GRID_VELOCITY
     &     (NNW, NGP, NW_IPAR, NW_RPAR, CRD, PHI, PHN, DCRD_DT)
C
      CALL TO_MATERIAL
     &     (NGP, CRD, PHI, PHN, DCRD_DT, DPHI_DT)
C
      RETURN
      END
      SUBROUTINE GRID_VELOCITY
     &           (NNW, NGP, NW_IPAR, NW_RPAR, CRD, PHI, PHN, V_GRID)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      INTEGER(KIND=IK) NNW
      INTEGER(KIND=IK) NGP
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)
      REAL   (KIND=RK) NW_RPAR (N_RP, *)
      REAL   (KIND=RK) CRD (4, *)
      REAL   (KIND=RK) PHI (2, *)
      REAL   (KIND=RK) PHN (2, *)
      REAL   (KIND=RK) V_GRID (4, *)
C
      CALL BASIC_VELOCITY
     &     (NGP, CRD, PHI, PHN, V_GRID)
C
      CALL CORRECT_TANGENTIAL
     &     (NNW, NW_IPAR, NW_RPAR, CRD, PHI, PHN, V_GRID)
C
      CALL ALIGN_NORMAL
     &     (NNW, NW_IPAR, NW_RPAR, CRD, PHI, PHN, V_GRID)
C
      CALL ALIGN_TANGENTIAL
     &     (NNW, NW_IPAR, NW_RPAR, CRD, PHI, PHN, V_GRID)
C
      RETURN
      END
      SUBROUTINE BASIC_VELOCITY
     &           (NGP, CRD, PHI, PHN, V_GRID)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      INTEGER(KIND=IK) NGP
      REAL   (KIND=RK) CRD (4, *)
      REAL   (KIND=RK) PHI (2, *)
      REAL   (KIND=RK) PHN (2, *)
      REAL   (KIND=RK) V_GRID (4, *)
C
      INTEGER(KIND=IK) JGP
      REAL   (KIND=RK) X_XI, Z_XI, J_INV
C
      DO JGP = 1, NGP
         X_XI  = CRD (3, JGP)
         Z_XI  = CRD (4, JGP)
         J_INV = 1.0D0 / SQRT (X_XI ** 2 + Z_XI ** 2)
C        d/dT (X) = Grad (Phi) = (phi_s, phi_n)
         V_GRID (1, JGP) = PHI (2, JGP) * J_INV
         V_GRID (2, JGP) = PHN (1, JGP)
      END DO
C
      RETURN
      END
      SUBROUTINE CORRECT_TANGENTIAL
     &           (NNW, NW_IPAR, NW_RPAR, CRD, PHI, PHN, V_GRID)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'bct_params.inc'
      INCLUDE 'net_funcs.inc'
      INCLUDE 'tme_funcs.inc'
C
      INTEGER(KIND=IK) NNW
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)
      REAL   (KIND=RK) NW_RPAR (N_RP, *)
      REAL   (KIND=RK) CRD (4, *)
      REAL   (KIND=RK) PHI (2, *)
      REAL   (KIND=RK) PHN (2, *)
      REAL   (KIND=RK) V_GRID (4, *)
C
      INTEGER(KIND=IK) INW, IGP, JGP, NGP, BCT
      REAL   (KIND=RK) PAR (N_GMP)
      REAL   (KIND=RK) ALPHA, DT
      REAL   (KIND=RK) S_XI (100), S_XIXI (100)
      INTEGER(KIND=IK) BC (2)
      DATA             BC /0, 0/
C
      IGP = 1
C
      DO INW = 1, NNW
         NGP = GET_NGP (NW_IPAR (1, INW))
         BCT = GET_BCT (NW_IPAR (1, INW))
         IF (BCT .EQ. -BCT_BERNOULLI) THEN
            CALL GET_GMP (NW_RPAR (1, INW), PAR)
            DT    = GET_DELTA_T ()
            ALPHA = PAR (1) / DT
C           only correct the tangential velocities of networks that
C           are  bernoulli kind of networks
            IF (ALPHA .NE. 0.0D0) THEN
               DO JGP = IGP, IGP + NGP - 1
                  S_XI(JGP) = SQRT (CRD (3,JGP) ** 2 + CRD (4,JGP) ** 2)
               END DO
               CALL SPLINE (BC, NGP, 1, 1, S_XI, S_XIXI)
               DO JGP = IGP + 1, IGP + NGP - 2
                  V_GRID (1,JGP) = V_GRID (1,JGP) + ALPHA * S_XIXI (JGP)
               END DO
            END IF
         END IF
         IGP = IGP + NGP
      END DO
C
      RETURN
      END
      SUBROUTINE ALIGN_NORMAL
     &           (NNW, NW_IPAR, NW_RPAR, CRD, PHI, PHN, V_GRID)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'bct_params.inc'
      INCLUDE 'net_funcs.inc'
C
      INTEGER(KIND=IK) NNW
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)
      REAL   (KIND=RK) NW_RPAR (N_RP, *)
      REAL   (KIND=RK) CRD (4, *)
      REAL   (KIND=RK) PHI (2, *)
      REAL   (KIND=RK) PHN (2, *)
      REAL   (KIND=RK) V_GRID (4, *)
C
      INTEGER(KIND=IK) INW, IGP, JGP, NGP, IED, ANW, AGP, ACT, BCT
      REAL   (KIND=RK) U_A (2)
C
      IGP = 1
C
      DO INW = 1, NNW
         NGP = GET_NGP (NW_IPAR (1, INW))
         BCT = GET_BCT (NW_IPAR (1, INW))
         DO IED = 1, 2
            ANW = GET_ANW (IED, NW_IPAR (1, INW))
            ACT = GET_BCT      (NW_IPAR (1, ANW))
            IF  (ACT .EQ. -BCT_BERNOULLI
     &     .AND. BCT .NE.  BCT_WAVEMAKER) THEN
C              only align the normal velocities of lateral networks that
C              are coupled to a network with a bernoulli kind of bc
               JGP = IGP + (IED - 1) * (NGP - 1)
               ANW = GET_ANW (IED,      NW_IPAR (1, INW))
               AGP = GET_AGP (IED, INW, NW_IPAR)
               CALL EDGE_VELOCITY
     &              (JGP, AGP, CRD, V_GRID, U_A)
C              Make the lateral network have a uniform normal velocity
               DO JGP = IGP, IGP + NGP - 1
                  V_GRID (2, JGP) = U_A (2)
C                 V_GRID (2, JGP) = PHN (1, IGP + (NGP - 1) * (IED - 1))
               END DO
            END IF
         END DO
         IGP = IGP + NGP
      END DO
C
      RETURN
      END
      SUBROUTINE ALIGN_TANGENTIAL
     &           (NNW, NW_IPAR, NW_RPAR, CRD, PHI, PHN, V_GRID)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'bct_params.inc'
      INCLUDE 'net_funcs.inc'
C
      INTEGER(KIND=IK) NNW
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)
      REAL   (KIND=RK) NW_RPAR (N_RP, *)
      REAL   (KIND=RK) CRD (4, *)
      REAL   (KIND=RK) PHI (2, *)
      REAL   (KIND=RK) PHN (2, *)
      REAL   (KIND=RK) V_GRID (4, *)
C
      INTEGER(KIND=IK) INW, IGP, JGP, NGP, IED, AGP, BCT
      REAL   (KIND=RK) U_A (2, 2), DVS
C
      IGP = 1
C
      DO INW = 1, NNW
         NGP = GET_NGP (NW_IPAR (1, INW))
         BCT = GET_BCT (NW_IPAR (1, INW))
         IF (BCT .NE. -BCT_BERNOULLI) THEN
C           only align the tangential velocities of networks that
C           are NOT bernoulli kind of networks
            DO IED = 1, 2
               JGP = IGP + (IED - 1) * (NGP - 1)
               AGP = GET_AGP (IED, INW, NW_IPAR)
               CALL EDGE_VELOCITY
     &              (JGP, AGP, CRD, V_GRID, U_A (1, IED))
            END DO
C           make the tangential velocity be linear, starting at
C           U_A (1, 1) and ending at U_A (1, 2)
            DVS = (U_A (1, 2) - U_A (1, 1)) / DBLE (NGP - 1)
            DO JGP = IGP, IGP + NGP - 1
               V_GRID (1, JGP) = U_A (1, 1) + DVS * DBLE (JGP - IGP)
            END DO
         END IF
         IGP = IGP + NGP
      END DO
C
      RETURN
      END
      SUBROUTINE EDGE_VELOCITY
     &           (IGP, AGP, CRD, VEL, U_A)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'net_funcs.inc'
C
      INTEGER(KIND=IK) IGP
      INTEGER(KIND=IK) AGP
      REAL   (KIND=RK) CRD (4, *)
      REAL   (KIND=RK) VEL (4, *)
      REAL   (KIND=RK) U_A (2)
C
      REAL   (KIND=RK) X_XI, Z_XI, J_INV, U_G (2)
C     GET Velocities of adjacent network in global coordinates
      X_XI    = CRD (3, AGP)
      Z_XI    = CRD (4, AGP)
      J_INV   = 1.0D0 / SQRT (X_XI ** 2 + Z_XI ** 2)
      U_G (1) = (X_XI * VEL (1, AGP) - Z_XI * VEL (2, AGP)) * J_INV
      U_G (2) = (X_XI * VEL (2, AGP) + Z_XI * VEL (1, AGP)) * J_INV
C     And transform them to the local coordinates
      X_XI    = CRD (3, IGP)
      Z_XI    = CRD (4, IGP)
      J_INV   = 1.0D0 / SQRT (X_XI ** 2 + Z_XI ** 2)
      U_A (1) = (X_XI * U_G (1) + Z_XI * U_G (2)) * J_INV
      U_A (2) = (X_XI * U_G (2) - Z_XI * U_G (1)) * J_INV
C
      RETURN
      END
      SUBROUTINE TO_MATERIAL
     &           (NGP, CRD, PHI, PHN, V_GRID, DPHI_DT)
C ---------------------------------------------------------------------------
C     Converts the partial derivatives to the material derivatives and
C     converts the grid velocity from local to global coordinates.
C
C     D/Dt = d/dt + dx/dt * Grad
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'tme_funcs.inc'
C
      INTEGER(KIND=IK) NGP
      REAL   (KIND=RK) CRD (4, *)
      REAL   (KIND=RK) PHI (2, *)
      REAL   (KIND=RK) PHN (2, *)
      REAL   (KIND=RK) V_GRID  (4, *)
      REAL   (KIND=RK) DPHI_DT (2, *)
C
      INTEGER(KIND=IK) IGP
      REAL   (KIND=RK) X_XI, Z_XI, J_INV, PHI_S, PHI_N, V_S, V_N
C
      DO IGP = 1, NGP
         X_XI  = CRD (3, IGP)
         Z_XI  = CRD (4, IGP)
         J_INV = 1.0D0 / SQRT (X_XI ** 2 + Z_XI ** 2)
C        Grad Phi in local coordinates
         PHI_S = PHI (2, IGP) * J_INV
         PHI_N = PHN (1, IGP)
C        D/DT (Phi)= d/dT (Phi) + V_grid * Grad (Phi)
         DPHI_DT (1, IGP) = DPHI_DT (1, IGP)
     &                    + V_GRID  (1, IGP) * PHI_S
     &                    + V_GRID  (2, IGP) * PHI_N
C        local grid velocity
         V_S = V_GRID (1, IGP)
         V_N = V_GRID (2, IGP)
C        Convert to global velocity
         V_GRID (1, IGP)  = (X_XI * V_S - Z_XI * V_N) * J_INV
         V_GRID (2, IGP)  = (X_XI * V_N + Z_XI * V_S) * J_INV
C        WRITE (*, 1001) IGP, V_GRID (1, IGP), V_GRID (2, IGP)
      END DO
C
C1001 FORMAT (1X, I3, ':', 10E14.6)
      RETURN
      END
