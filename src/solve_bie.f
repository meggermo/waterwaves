      SUBROUTINE SOLVE_BIE
     &           (NSD,    NNW_SD, NGP_SD,  NW_IPAR,
     &            SRC_0,  SRC_1,  DIP_0,   DIP_1,
     &            CRD,    PHI,    PHN)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'slv_funcs.inc'
C
      INTEGER(KIND=IK) NSD
      INTEGER(KIND=IK) NNW_SD (*)
      INTEGER(KIND=IK) NGP_SD (*)
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)
      REAL   (KIND=RK) SRC_0 (*)
      REAL   (KIND=RK) SRC_1 (*)
      REAL   (KIND=RK) DIP_0 (*)
      REAL   (KIND=RK) DIP_1 (*)
      REAL   (KIND=RK) CRD (4, *)
      REAL   (KIND=RK) PHI (2, *)
      REAL   (KIND=RK) PHN (2, *)
C
      LOGICAL(KIND=LK) CONVERGED
      INTEGER(KIND=IK) COUNTER
      INTEGER(KIND=IK) MAX_IT
      REAL   (KIND=RK) RES (2)
      REAL   (KIND=RK) TOL (2)
      REAL   (KIND=RK) OMG (2)
C
      MAX_IT  = GET_DOD_MAXIT ()
      TOL (1) = GET_DOD_TOL (1)
      TOL (2) = GET_DOD_TOL (2)
      OMG (1) = GET_DOD_OMG (1)
      OMG (2) = GET_DOD_OMG (2)
C
      CONVERGED = .FALSE.
      COUNTER   = 0
C
      DO WHILE (.NOT. CONVERGED .AND. COUNTER .LT. MAX_IT)
         COUNTER = COUNTER + 1
C        Solve the BIE's for the individual subdomains
         CALL SOLVE_SUBDOMAINS
     &        (NSD,     NNW_SD,  NGP_SD,   NW_IPAR,
     &         SRC_0,   SRC_1,   DIP_0,    DIP_1,
     &         CRD,     PHI,     PHN)
C        Exchange the interface data
         CALL RELAX
     &        (NSD, NNW_SD, NW_IPAR, OMG, RES, PHI, PHN)
C        Check for convergence
         CONVERGED =
     &   RES (1) .LE. TOL (1) .AND. RES (2) .LE. TOL (2)
         IF (NSD .GT. 1) WRITE (6, 1001) RES, COUNTER
      END DO
C
      IF (.NOT. CONVERGED) THEN
         CALL ERROR ('SOLVE_BIE', 'No convergence')
      END IF
C
      IF (NSD .GT. 1) WRITE (6, 1101) RES, COUNTER
C
      RETURN
 1001 FORMAT (1X,'DDC:',2E9.1,' AT ITERATION ', I3)
 1101 FORMAT (1X,'DDC:',2E9.1,' AFTER ', I3, ' ITERATIONS')
      END
      SUBROUTINE SOLVE_SUBDOMAINS
     &   (NSD,   NNW_SD, NGP_SD, NW_IPAR,
     &    SRC_0, SRC_1,  DIP_0,  DIP_1,
     &    CRD,   PHI,    PHN)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'slv_funcs.inc'
C
      INTEGER(KIND=IK) NSD
      INTEGER(KIND=IK) NNW_SD (*)
      INTEGER(KIND=IK) NGP_SD (*)
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)
      REAL   (KIND=RK) SRC_0 (*)
      REAL   (KIND=RK) SRC_1 (*)
      REAL   (KIND=RK) DIP_0 (*)
      REAL   (KIND=RK) DIP_1 (*)
      REAL   (KIND=RK) CRD (4, *)
      REAL   (KIND=RK) PHI (2, *)
      REAL   (KIND=RK) PHN (2, *)
C
      INTEGER(KIND=IK) ISD, INW, IGP, IAE
      INTEGER(KIND=IK) NNW, NGP
      INTEGER(KIND=IK) MAX_IT
      REAL   (KIND=RK) TOL (2)
      REAL   (KIND=RK) OMG (2)
C
      MAX_IT  = GET_ITR_MAXIT ()
      TOL (1) = GET_ITR_TOL (1)
      TOL (2) = GET_ITR_TOL (2)
      OMG (1) = GET_ITR_OMG (1)
      OMG (2) = GET_ITR_OMG (2)
C
      INW = 1
      IGP = 1
      IAE = 1
C
      DO ISD = 1, NSD
         NNW = NNW_SD (ISD)
         NGP = NGP_SD (ISD)
C        Replace a BIE by a coupling condition for each network
         CALL COUPLE_EDGES
     &        (NNW,          NGP,           NW_IPAR (1, INW),
     &         CRD (1, IGP),
     &         SRC_0  (IAE), SRC_1   (IAE),
     &         DIP_0  (IAE), DIP_1   (IAE))
C        Solve the BIE's with the couping conditions
         CALL SOLVE_SD
     &        (NNW, NGP,
     &         NW_IPAR (1, INW),
     &         SRC_0  (IAE), SRC_1  (IAE),
     &         DIP_0  (IAE), DIP_1  (IAE),
     &         OMG,          MAX_IT,       TOL,
     &         CRD (1, IGP), PHI (1, IGP), PHN  (1, IGP))
C        Increment the index pointers to the next subdomain
         INW = INW + NNW
         IGP = IGP + NGP
         IAE = IAE + NGP ** 2
      END DO
C
      RETURN
      END
      SUBROUTINE SOLVE_SD
     &   (NNW,   NGP,
     &    NW_IPAR,
     &    SRC_0, SRC_1,  DIP_0,  DIP_1,
     &    OMG,   MAX_IT, TOL,
     &    CRD,   PHI ,   PHN)
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
      REAL   (KIND=RK) SRC_0 (NGP, NGP)
      REAL   (KIND=RK) SRC_1 (NGP, NGP)
      REAL   (KIND=RK) DIP_0 (NGP, NGP)
      REAL   (KIND=RK) DIP_1 (NGP, NGP)
      REAL   (KIND=RK) OMG (*)
      INTEGER(KIND=IK) MAX_IT
      REAL   (KIND=RK) TOL (*)
      REAL   (KIND=RK) CRD (4, *)
      REAL   (KIND=RK) PHI (2, *)
      REAL   (KIND=RK) PHN (2, *)
C
      REAL   (KIND=RK) A0 (NGP, NGP)
      REAL   (KIND=RK) A1 (NGP, NGP)
      REAL   (KIND=RK) B  (NGP)
      REAL   (KIND=RK) X  (NGP, 2)
C
      CALL SETUP_AX_B
     &     (NNW,    NGP,    NW_IPAR,
     &      PHI,    PHN,
     &      SRC_0,  SRC_1,  DIP_0, DIP_1,
     &      A0,     A1,     B,     X)
C
      CALL ITERATE
     &     (NNW,    NGP,    NW_IPAR,
     &      CRD,    PHI,    PHN,
     &      A0,     A1,     B,
     &      OMG,    MAX_IT, TOL, X)
C
      CALL COPY_BACK
     &     (NNW, NGP, NW_IPAR, X, PHI, PHN)
C
      RETURN
      END
      SUBROUTINE SETUP_AX_B
     &   (NNW,  NGP,
     &    NW_IPAR,
     &    PHI,  PHN,
     &    S0,   S1,   D0, D1,
     &    A0,   A1,   B,  X)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'chk_params.inc'
      INCLUDE 'net_funcs.inc'
C
      INTEGER(KIND=IK) NNW
      INTEGER(KIND=IK) NGP
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)
      REAL   (KIND=RK) PHI (2, NGP)
      REAL   (KIND=RK) PHN (2, NGP)
      REAL   (KIND=RK) S0 (NGP, NGP)
      REAL   (KIND=RK) S1 (NGP, NGP)
      REAL   (KIND=RK) D0 (NGP, NGP)
      REAL   (KIND=RK) D1 (NGP, NGP)
      REAL   (KIND=RK) A0 (NGP, NGP)
      REAL   (KIND=RK) A1 (NGP, NGP)
      REAL   (KIND=RK) B  (NGP)
      REAL   (KIND=RK) X  (NGP, 2)
C
      INTEGER(KIND=IK) INW, IGP, BCT, NGP_NW
C
      DO IGP = 1, NGP
         B (IGP) = 0.0D0
      END DO
C
      IGP = 1
      DO INW = 1, NNW
         BCT    = GET_BCT (NW_IPAR (1, INW))
         NGP_NW = GET_NGP (NW_IPAR (1, INW))
         IF (BCT .LT. 0) THEN
            CALL COPY_COLOUMNS (NGP, NGP_NW,
     &          PHI (1, IGP), PHN (1, IGP),
     &          S0  (1, IGP), S1  (1, IGP),  D0  (1, IGP), D1  (1, IGP),
     &          A0  (1, IGP), A1  (1, IGP),  B,  X (IGP, 1), X (IGP, 2))
         ELSE
            CALL COPY_COLOUMNS
     &         (NGP, NGP_NW,
     &          PHN (1, IGP), PHI (1, IGP),
     &          D0  (1, IGP), D1  (1, IGP), S0  (1, IGP),  S1 (1, IGP),
     &          A0  (1, IGP), A1  (1, IGP), B,  X (IGP, 1), X (IGP, 2))
         END IF
         IGP = IGP + NGP_NW
      END DO
C     Write Matrices and Vectors to file, so that we can
C     evaluate the following equation: A0 * X0 + A1 * X1 = B
      IF (CHK_EQNS) THEN
         CALL WRITE_MATRIX
     &        ('A0.out', NGP, NGP, A0)
         CALL WRITE_MATRIX
     &        ('A1.out', NGP, NGP, A1)
         CALL WRITE_VECTOR
     &        ('B.out', NGP, B)
         CALL WRITE_VECTOR
     &        ('X0.out', NGP, X (1, 1))
         CALL WRITE_VECTOR
     &        ('X1.out', NGP, X (1, 2))
      END IF
C
      RETURN
      END
      SUBROUTINE COPY_BACK (NNW, NGP, NW_IPAR, X, PHI, PHN)
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
      REAL   (KIND=RK) X   (NGP, 2)
      REAL   (KIND=RK) PHI (2, NGP)
      REAL   (KIND=RK) PHN (2, NGP)
C
      INTEGER(KIND=IK)  INW, IGP, J, BCT, NGP_NW
C
      IGP = 1
      DO INW = 1, NNW
         BCT    = GET_BCT (NW_IPAR (1, INW))
         NGP_NW = GET_NGP (NW_IPAR (1, INW))
         IF (BCT .LT. 0) THEN
            DO J = IGP, IGP + NGP_NW - 1
               PHN (1, J) = X (J, 1)
               PHN (2, J) = X (J, 2)
            END DO
         ELSE
            DO J = IGP, IGP + NGP_NW - 1
               PHI (1, J) = X (J, 1)
               PHI (2, J) = X (J, 2)
            END DO
         END IF
         IGP = IGP + NGP_NW
      END DO
C
      RETURN
      END
      SUBROUTINE COPY_COLOUMNS (NGP, N_GP, V, W, A_0, A_1, B_0, B_1, A0,
     &           A1, B, X, DX)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      INTEGER(KIND=IK) NGP
      INTEGER(KIND=IK) N_GP
      REAL   (KIND=RK) V (2, N_GP)
      REAL   (KIND=RK) W (2, N_GP)
      REAL   (KIND=RK) A_0 (NGP, N_GP)
      REAL   (KIND=RK) A_1 (NGP, N_GP)
      REAL   (KIND=RK) B_0 (NGP, N_GP)
      REAL   (KIND=RK) B_1 (NGP, N_GP)
      REAL   (KIND=RK) A0  (NGP, N_GP)
      REAL   (KIND=RK) A1  (NGP, N_GP)
      REAL   (KIND=RK) B   (NGP)
      REAL   (KIND=RK) X   (N_GP)
      REAL   (KIND=RK) DX  (N_GP)
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
C-OLD -------------------------------------------------------------------------
C     SUBROUTINE COUPLE_EDGES
C    &   (NNW,   NGP,   NW_IPAR,
C    &    CRD,
C    &    SRC_0, SRC_1,
C    &    DIP_0, DIP_1,
C    &    NEU)
C ---------------------------------------------------------------------------
C       Transformation rule for the velocities from local (S,N) to the
C       global coordinates (X,Z) are:
C
C      1  | X_XI -Z_XI |   | PHI_S |   | PHI_X |
C      -  |            | * |       | = |       |
C      J  | Z_XI  X_XI |   | PHI_N |   | PHI_Z |
C
C      where J = SQRT (X_XI ** 2 + Z_XI ** 2) is the Jacobian
C
C      This can be used to demand equality of velocities over networks
C      as follows:
C                          I                A           A           I
C      1     | X_XI  Z_XI |   | X_XI -Z_XI |   | PHI_S |   | PHI_S |
C  --------- |            | * |            | * |       | = |       |
C  J_I * J_A |-Z_XI  X_XI |   | Z_XI  X_XI |   | PHI_N |   | PHI_N |
C
C      Where we have used the fact that the transformation matrix is
C      orthonormal.
C      Note that PHI_S = 1 / J * PHI_XI, so the two equations for
C      continuity of velocity are given by
C
C     -1                        -1
C    J * (PHI_XI)  - (R_11 * J * (PHI_XI) + R_12 * PHI_N) = 0
C     I          I              A          A              A
C
C                               -1
C        (PHI_N )  + (R_21 * J * (PHI_XI) - R_11 * PHI_N) = 0
C                I              A          A              A
C     where
C
C      R_11 = (X_XI * X_XI + Z_XI * Z_XI) / J / J
C                  I      A      I      A    I   A
C
C      R_12 = (Z_XI * X_XI - X_XI * Z_XI) / J / J
C                  I      A      I      A    I   A
C
C ---------------------------------------------------------------------------
C     IMPLICIT NONE
C     INCLUDE 'knd_params.inc'
C     INCLUDE 'net_params.inc'
C     INCLUDE 'net_funcs.inc'
C
C     INTEGER(KIND=IK) NNW
C     INTEGER(KIND=IK) NGP
C     INTEGER(KIND=IK) NW_IPAR (N_IP, *)
C     REAL   (KIND=RK) CRD (4, *)
C     REAL   (KIND=RK) SRC_0 (NGP, *)
C     REAL   (KIND=RK) SRC_1 (NGP, *)
C     REAL   (KIND=RK) DIP_0 (NGP, *)
C     REAL   (KIND=RK) DIP_1 (NGP, *)
C     INTEGER(KIND=IK) NEU
C
C     INTEGER(KIND=IK) INW, IGP, IED, TBC, JGP
C     INTEGER(KIND=IK) ANW, AGP, ABC
C     REAL   (KIND=RK) J_I, J_A
C     REAL   (KIND=RK) X_XII, Z_XII, X_XIA, Z_XIA
C     REAL   (KIND=RK) R_11, R_12
C
C     IED = 1
C     IGP = 1
C     NEU = 0
C     DO INW = 1, NNW
C        DO JGP = 1, NGP
C           SRC_0 (IGP, JGP) = 0.0D0
C           SRC_1 (IGP, JGP) = 0.0D0
C           DIP_0 (IGP, JGP) = 0.0D0
C           DIP_1 (IGP, JGP) = 0.0D0
C        END DO
C        ANW = GET_ANW (IED, NW_IPAR (1, INW))
C        AGP = GET_AGP (IED, INW, NW_IPAR)
C        TBC = GET_BCT (NW_IPAR (1, INW))
C        ABC = GET_BCT (NW_IPAR (1, ANW))
C        X_XII = CRD (3, IGP)
C        Z_XII = CRD (4, IGP)
C        J_I   = SQRT (X_XII ** 2 + Z_XII ** 2)
C        X_XIA = CRD (3, AGP)
C        Z_XIA = CRD (4, AGP)
C        J_A   = SQRT (X_XIA ** 2 + Z_XIA ** 2)
C        R_11  = (X_XII * X_XIA + Z_XII * Z_XIA) / J_I / J_A
C        R_12  = (Z_XII * X_XIA - X_XII * Z_XIA) / J_I / J_A
C        IF (ABC .GT. 0) THEN
C           Continuity of normal velocity for Dirichlet BC's
C           SRC_0 (IGP, IGP) =  1.0D0
C           DIP_1 (IGP, AGP) =  R_12 / J_A
C           SRC_0 (IGP, AGP) = -R_11
C        ELSE
C           NEU = NEU + 1
C           Continuity of potential for Neumann BC's
C           DIP_0 (IGP, IGP) =  1.0D0
C           DIP_0 (IGP, AGP) = -1.0D0
C        END IF
C        IGP = IGP + GET_NGP (NW_IPAR (1, INW))
C     END DO
C
C     IF (NEU .EQ. NNW) THEN
C        If all boundary conditions are of Neumann type then
C        we'll set PHI (NGP) to some constant
C        WRITE (*, *) 'ALL NEUMANN'
C        DO IGP = 1, NGP
C           SRC_0 (NGP, IGP) = 0.0D0
C           SRC_1 (NGP, IGP) = 0.0D0
C           DIP_0 (NGP, IGP) = 0.0D0
C           DIP_1 (NGP, IGP) = 0.0D0
C        END DO
C        DIP_0 (NGP, NGP) = 1.0D0
C     END IF
C
C     RETURN
C     END
C-OLD -------------------------------------------------------------------------
