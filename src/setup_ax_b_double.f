      SUBROUTINE SETUP_AX_B_DOUBLE
     &   (NNW,      NGP,
     &    NW_IPAR,
     &    CRD,      PHI, PHN,
     &    S0,       S1,
     &    D0,       D1,
     &    A,        X,   B)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'net_funcs.inc'
C
      INTEGER(KIND=IK) NNW
      INTEGER(KIND=IK) NGP
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)
      REAL   (KIND=RK) CRD (4, *)
      REAL   (KIND=RK) PHI (2, *)
      REAL   (KIND=RK) PHN (2, *)
      REAL   (KIND=RK) S0  (NGP, *)
      REAL   (KIND=RK) S1  (NGP, *)
      REAL   (KIND=RK) D0  (NGP, *)
      REAL   (KIND=RK) D1  (NGP, *)
      REAL   (KIND=RK) A   (NGP * 2, *)
      REAL   (KIND=RK) X   (*)
      REAL   (KIND=RK) B   (*)
C
      INTEGER(KIND=IK) INW, IED, IGP, BCT
      INTEGER(KIND=IK) ANW, AED, AGP, ACT
      INTEGER(KIND=IK) JGP, KGP, BGP, NGP_NW
      REAL   (KIND=RK) J_I, J_I2, J_A, J_A2
      REAL   (KIND=RK) X_XIXI, Z_XIXI, P_XIXI, PHI_SS, DOT
      REAL   (KIND=RK) R1 (2), R2 (2)
      REAL   (KIND=RK) EVAL_2
      EXTERNAL         EVAL_2
C
      DO IGP = 1, NGP * 2
         B (IGP) = 0.0D0
         DO JGP = 1, NGP * 2
            A (JGP, IGP) = 0.0D0
         END DO
      END DO
C
      IGP = 1
      DO INW = 1, NNW
         JGP    = IGP + NGP
         BCT    = GET_BCT (NW_IPAR (1, INW))
         NGP_NW = GET_NGP (NW_IPAR (1, INW))
         IF (BCT .LT. 0) THEN
            CALL COPY_COLOUMNS_DOUBLE
     &           (NGP, NGP_NW,
     &            PHI (1, IGP), PHN (1, IGP),
     &            S0  (1, IGP), S1  (1, IGP),
     &            D0  (1, IGP), D1  (1, IGP),
     &            A   (1, IGP), A   (1, JGP),  B,  X (IGP), X (JGP))
         ELSE
            CALL COPY_COLOUMNS_DOUBLE
     &           (NGP, NGP_NW,
     &            PHN (1, IGP), PHI (1, IGP),
     &            D0  (1, IGP),  D1 (1, IGP),
     &            S0  (1, IGP),  S1 (1, IGP),
     &            A   (1, IGP),  A  (1, JGP), B,  X (IGP), X (JGP))
         END IF
C        SETUP THE DIAGONALS FOR THE SPLINES
         DO JGP = IGP + 1, IGP + NGP_NW - 2
            KGP = JGP + NGP
            A (KGP, JGP - 1) =  0.75D0
            A (KGP, JGP + 1) = -0.75D0
            A (KGP, KGP - 1) =  0.25D0
            A (KGP, KGP)     =  1.00D0
            A (KGP, KGP + 1) =  0.25D0
         END DO
         IGP = IGP + NGP_NW
      END DO
C
      IGP = 1
      DO INW = 1, NNW
         NGP_NW = GET_NGP (NW_IPAR (1, INW))
         BCT    = GET_BCT (NW_IPAR (1, INW))
         DO IED = 1, 2
            JGP = IGP + (IED - 1) * (NGP_NW - 1)
            ANW = GET_ANW (IED,      NW_IPAR (1, INW))
            AED = GET_AED (IED,      NW_IPAR (1, INW))
            AGP = GET_AGP (IED, INW, NW_IPAR)
            ACT = GET_BCT (NW_IPAR (1, ANW))
C
            J_I2 = 1.0D0 / (CRD (3, JGP) ** 2 + CRD (4, JGP) ** 2)
            J_A2 = 1.0D0 / (CRD (3, AGP) ** 2 + CRD (4, AGP) ** 2)
            J_I  = SQRT (J_I2)
            J_A  = SQRT (J_A2)
            R1 (1) = (CRD (3, JGP) * CRD (3, AGP)
     &             +  CRD (4, JGP) * CRD (4, AGP)) * J_I * J_A
            R1 (2) = (CRD (3, JGP) * CRD (4, AGP)
     &             -  CRD (4, JGP) * CRD (3, AGP)) * J_I * J_A
            R2 (1) = (R1 (1) + R1 (2)) * (R1 (1) - R1 (2))
            R2 (2) = 2.0D0 * R1 (1) * R1 (2)
            X_XIXI = EVAL_2 (DBLE (AED-1), 2, CRD (1, AGP + 1 - AED))
            Z_XIXI = EVAL_2 (DBLE (AED-1), 2, CRD (2, AGP + 1 - AED))
            P_XIXI = EVAL_2 (DBLE (AED-1), 1, PHI (1, AGP + 1 - AED))
            DOT    = CRD (3, AGP) * X_XIXI + CRD (4, AGP) * Z_XIXI
            PHI_SS = (P_XIXI - DOT * J_A2 * PHI (2, AGP)) * J_A2
C
            KGP = JGP + NGP
            BGP = AGP + NGP
C
            IF (BCT .GT. 0) THEN
               IF (ACT .GT. 0) THEN
C                 NN -> PHI_S SHOULD BE CONTINUOUS
                  A (KGP, KGP) =  J_I
                  A (KGP, BGP) = -R1 (1) * J_A
                  B (KGP)      = -R1 (2) * PHN (1, AGP)
               ELSE
C                 ND -> PHI SHOULD BE CONTINUOUS
C                 A (KGP, JGP) =  1.0D0
C                 A (KGP, AGP) = -1.0D0
C                 B (KGP)      =  0.0D0
C                 ND -> PHI_S SHOULD BE CONTINUOUS
                  A (KGP, KGP) =  J_I
                  A (KGP, AGP) =  R1 (2)
                  B (KGP)      =  R1 (1) * J_A * PHI (2, AGP)
               END IF
            ELSE
               IF (ACT .LT. 0) THEN
C                 DD -> PHI_N SHOULD BE CONTINUOUS
C                 A (KGP, JGP) =  1.0D0
C                 A (KGP, AGP) = -R1 (1)
C                 B (KGP)      =  R1 (2) * J_A * PHI (2, AGP)
C                 DD -> PHI_NS SHOULD BE CONTINUOUS
                  A (KGP, KGP) =  J_I
                  A (KGP, BGP) = -R2 (1) * J_A
                  B (KGP)      =  R2 (2) * PHI_SS
                  WRITE (*, *) PHN (2, JGP) * J_I
     &              - R2 (2) * PHI_SS
     &              - R2 (1) * PHN (2, AGP) * J_A,
     &                         PHI (2, JGP) * J_I
     &              - R1 (1) * PHI (2, AGP) * J_A
     &              + R1 (2) * PHN (1, AGP)
               ELSE
C                 DN -> PHI_N SHOULD BE CONTINUOUS
                  A (KGP, JGP) =  1.0D0
                  A (KGP, BGP) = -R1 (2) * J_A
                  B (KGP)      =  R1 (1) * PHN (1, AGP)
C                 DN -> PHI_NS SHOULD BE CONTINUOUS
C                 A (KGP, KGP) =  J_I
C                 B (KGP)      =  R2 (1) * PHN (2, AGP) * J_A
C                 IF (AED .EQ. 1) THEN
C                    A (KGP, AGP)     = -R2 (2) * (-6.0D0) * J_A2
C                    A (KGP, AGP + 1) = -R2 (2) * ( 6.0D0) * J_A2
C    &                                +  R2 (2) *   DOT    * J_A2 ** 2
C                    A (KGP, BGP)     = -R2 (2) * (-4.0D0) * J_A2
C                    A (KGP, BGP + 1) = -R2 (2) * (-2.0D0) * J_A2
C                    PHI_SS = A (KGP, AGP)     * PHI (1, AGP)
C    &                      + A (KGP, AGP + 1) * PHI (1, AGP + 1)
C    &                      + A (KGP, BGP)     * PHI (2, AGP)
C    &                      + A (KGP, BGP + 1) * PHI (2, AGP + 1)
C                 ELSE
C                    A (KGP, AGP - 1) = -R2 (2) * ( 6.0D0) * J_A2
C                    A (KGP, AGP)     = -R2 (2) * (-6.0D0) * J_A2
C    &                                +  R2 (2) *   DOT    * J_A2 ** 2
C                    A (KGP, BGP - 1) = -R2 (2) * ( 2.0D0) * J_A2
C                    A (KGP, BGP)     = -R2 (2) * ( 4.0D0) * J_A2
C                    PHI_SS = A (KGP, AGP)     * PHI (1, AGP)
C    &                      + A (KGP, AGP - 1) * PHI (1, AGP - 1)
C    &                      + A (KGP, BGP)     * PHI (2, AGP)
C    &                      + A (KGP, BGP - 1) * PHI (2, AGP - 1)
C                 END IF
               END IF
            END IF
         END DO
         IGP = IGP + NGP_NW
      END DO
C
      RETURN
      END
      SUBROUTINE CHECK_EQNS_DOUBLE (NGP, A0, A1, S0, S1, X0, X1, B0, B1)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      INTEGER(KIND=IK) NGP
      REAL   (KIND=RK) A0 (NGP * 2, *)
      REAL   (KIND=RK) A1 (NGP * 2, *)
      REAL   (KIND=RK) S0 (NGP * 2, *)
      REAL   (KIND=RK) S1 (NGP * 2, *)
      REAL   (KIND=RK) X0 (*)
      REAL   (KIND=RK) X1 (*)
      REAL   (KIND=RK) B0 (*)
      REAL   (KIND=RK) B1 (*)
C
      INTEGER(KIND=IK) I, J
      REAL   (KIND=RK) S
C
      WRITE (*, *) '=== CHECK EQNS ==='
      WRITE (*, *) ' BIE-EQNS'
      DO I = 1, NGP
         S = 0.0D0
         DO J = 1, NGP
            S = S + A0 (I, J) * X0 (J) + A1 (I, J) * X1 (J)
         END DO
         WRITE (*, '(1X,I3,'':'',E10.1)') I, S - B0 (I)
      END DO
      WRITE (*, *) ' SPL-EQNS'
      DO I = 1, NGP
         S = 0.0D0
         DO J = 1, NGP
            S = S + S0 (I, J) * X0 (J) + S1 (I, J) * X1 (J)
         END DO
         WRITE (*, '(1X,I3,'':'',E10.1)') I, S - B1 (I)
      END DO
      WRITE (*, *) '=== CHECK EQNS ==='
C
      RETURN
      END
      SUBROUTINE COPY_COLOUMNS_DOUBLE
     &           (NGP, N_GP, V, W, A_0, A_1, B_0, B_1, A0, A1, B, X, DX)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      INTEGER(KIND=IK) NGP
      INTEGER(KIND=IK) N_GP
      REAL   (KIND=RK) V (2, *)
      REAL   (KIND=RK) W (2, *)
      REAL   (KIND=RK) A_0 (NGP, *)
      REAL   (KIND=RK) A_1 (NGP, *)
      REAL   (KIND=RK) B_0 (NGP, *)
      REAL   (KIND=RK) B_1 (NGP, *)
      REAL   (KIND=RK) A0  (NGP * 2, *)
      REAL   (KIND=RK) A1  (NGP * 2, *)
      REAL   (KIND=RK) B   (*)
      REAL   (KIND=RK) X   (*)
      REAL   (KIND=RK) DX  (*)
C
      INTEGER(KIND=IK) I, J
C
      DO J = 1, N_GP
         X  (J) = W (1, J)
         DX (J) = W (2, J)
         DO I = 1, NGP
            A0 (I, J) = A_0 (I, J)
            A1 (I, J) = A_1 (I, J)
            B (I) = B (I)
     &            - B_0 (I, J) * V (1, J)
     &            - B_1 (I, J) * V (2, J)
         END DO
      END DO
C
      RETURN
      END

      FUNCTION EVAL_2 (T, LDF, F)
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      REAL   (KIND=RK) EVAL_2
      REAL   (KIND=RK) T
      INTEGER(KIND=IK) LDF
      REAL   (KIND=RK) F (LDF, 2, *)
C
      EVAL_2 =  6.0D0 * (2.0D0 * T - 1.0D0) * F (1, 1, 1)
     &       -  6.0D0 * (2.0D0 * T - 1.0D0) * F (1, 1, 2)
     &       + (6.0D0 *          T - 4.0D0) * F (1, 2, 1)
     &       + (6.0D0 *          T - 2.0D0) * F (1, 2, 2)
      RETURN
      END
