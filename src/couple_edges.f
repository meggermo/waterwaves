      SUBROUTINE COUPLE_EDGES
     &   (NNW,   NGP,   NW_IPAR,
     &    CRD,
     &    SRC_0, SRC_1,
     &    DIP_0, DIP_1)
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
C      1.0    | X_XI  Z_XI |   | X_XI -Z_XI |   | PHI_S |   | PHI_S |
C   --------- |            | * |            | * |       | = |       |
C   J_I * J_A |-Z_XI  X_XI |   | Z_XI  X_XI |   | PHI_N |   | PHI_N |
C
C      Where we have used the fact that the transformation matrix is
C      orthonormal.
C      Note that PHI_S = 1 / J * PHI_XI, so the two equations for
C      continuity of velocity are given by
C
C     -1                      -1
C    J * (PHI_XI)  - (R_11 * J * (PHI_XI) + R_12 * PHI_N) = 0
C     I          I            A          A              A
C
C                             -1
C        (PHI_N )  + (R_21 * J * (PHI_XI) - R_11 * PHI_N) = 0
C                I            A          A              A
C     where
C
C      R_11 = (X_XI * X_XI + Z_XI * Z_XI) / J / J
C                  I      A      I      A    I   A
C
C      R_12 = (Z_XI * X_XI - X_XI * Z_XI) / J / J
C                  I      A      I      A    I   A
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'net_funcs.inc'
C
      INTEGER(KIND=IK) NNW
      INTEGER(KIND=IK) NGP
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)
      REAL   (KIND=RK) CRD (4, *)
      REAL   (KIND=RK) SRC_0 (NGP, *)
      REAL   (KIND=RK) SRC_1 (NGP, *)
      REAL   (KIND=RK) DIP_0 (NGP, *)
      REAL   (KIND=RK) DIP_1 (NGP, *)
C
      LOGICAL(KIND=LK) ALL_NEUMANN
      INTEGER(KIND=IK) INW, IGP, LGP, IED, TBC, JGP, KGP, SPT
      INTEGER(KIND=IK) ANW, AGP
      REAL   (KIND=RK) J_I, J_A
      REAL   (KIND=RK) X_XII, Z_XII, X_XIA, Z_XIA
      REAL   (KIND=RK) R_11, R_12
C
      IGP = 1
      IED = 1
      ALL_NEUMANN = .TRUE.
      DO INW = 1, NNW
         LGP = GET_NGP (NW_IPAR (1, INW))
         TBC = GET_BCT (NW_IPAR (1, INW))
         IF (TBC .LT. 0) THEN
            ALL_NEUMANN = .FALSE.
            DO IED = 1, 2
               SPT = GET_SPT (2 + IED, NW_IPAR (1, INW))
C TODO: WHY?   SPT = 0
               ANW = GET_ANW (IED, NW_IPAR (1, INW))
               AGP = GET_AGP (IED, INW, NW_IPAR)
               IF (SPT .GT. 0) THEN
                  KGP = IGP + (IED - 1) * (LGP - 1)
                  DO JGP = 1, NGP
                     SRC_0 (KGP, JGP) = 0.0D0
                     SRC_1 (KGP, JGP) = 0.0D0
                     DIP_0 (KGP, JGP) = 0.0D0
                     DIP_1 (KGP, JGP) = 0.0D0
                  END DO
                  X_XII = CRD (3, KGP)
                  Z_XII = CRD (4, KGP)
                  J_I   = 1.0D0 / SQRT (X_XII ** 2 + Z_XII ** 2)
                  X_XIA = CRD (3, AGP)
                  Z_XIA = CRD (4, AGP)
                  J_A   = 1.0D0 / SQRT (X_XIA ** 2 + Z_XIA ** 2)
                  R_11  = (X_XII * X_XIA + Z_XII * Z_XIA) * J_I * J_A
                  R_12  = (Z_XII * X_XIA - X_XII * Z_XIA) * J_I * J_A
C                 Continuity of normal velocity for Dirichlet BC's
C                 PHN_K - R_11 * PHN_A + R_12 * PHS_A = 0
                  SRC_0 (KGP, KGP) =  1.0D0
                  SRC_0 (KGP, AGP) = -R_11
                  DIP_1 (KGP, AGP) =  R_12 * J_A
               END IF
            END DO
         END IF
         IGP = IGP + LGP
      END DO
C
      IF (ALL_NEUMANN) THEN
C        If all boundary conditions are of Neumann type then
C        we'll set PHI (NGP) to some constant
         WRITE (*, *) 'ALL NEUMANN'
         DO IGP = 1, NGP
            SRC_0 (NGP, IGP) = 0.0D0
            SRC_1 (NGP, IGP) = 0.0D0
            DIP_0 (NGP, IGP) = 0.0D0
            DIP_1 (NGP, IGP) = 0.0D0
         END DO
         DIP_0 (NGP, NGP) = 1.0D0
      END IF
C
      RETURN
      END
