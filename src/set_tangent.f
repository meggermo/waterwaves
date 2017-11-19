      SUBROUTINE SET_TANGENT (NNW, NW_IPAR, CRD, PHI, PHN, F, G)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'net_funcs.inc'
      INCLUDE 'usr_params.inc'
C
      INTEGER(KIND=IK) NNW
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)
      REAL   (KIND=RK) CRD (2, 2, *)
      REAL   (KIND=RK) PHI (2, *)
      REAL   (KIND=RK) PHN (2, *)
      REAL   (KIND=RK) F   (*)
      REAL   (KIND=RK) G   (*)
C
      INTEGER(KIND=IK) INW, IED
      INTEGER(KIND=IK) IGP, JGP, AGP, NGP, BCT, ANW, ACT
      INTEGER(KIND=IK) SPT (2)
      DATA             SPT /1, 1/
C
      REAL   (KIND=RK) U_G (2)
      REAL   (KIND=RK) U_L (2)
      REAL   (KIND=RK) X_XI, Z_XI, J_INV
C
      IGP = 1
      DO INW = 1, NNW
         NGP = GET_NGP (NW_IPAR (1, INW))
         BCT = GET_BCT (NW_IPAR (1, INW))
C         WRITE (USR_O, *) INW, 'G BEFORE SET_TANGENT'
C         WRITE (USR_O, 111) CRD(:,:, IGP : IGP + NGP - 1)
C         WRITE (USR_O, *)
C         WRITE (USR_O, 101) G (IGP : IGP + NGP - 1)
         DO IED = 1, 2
            JGP = IGP + (IED - 1) * (NGP - 1)
            ANW = GET_ANW (IED, NW_IPAR (1, INW))
            AGP = GET_AGP (IED, INW, NW_IPAR)
            ACT = GET_BCT (NW_IPAR (1, ANW))
C           GET Velocities of adjacent network in global coordinates
            X_XI  = CRD (1, 2, AGP)
            Z_XI  = CRD (2, 2, AGP)
            J_INV = 1.0D0 / SQRT (X_XI ** 2 + Z_XI ** 2)
            IF (ACT .LT. 0) THEN
C              Adjacent network is Dirichlet, so PHN is to be solved
               U_L (1)  = PHI (2, AGP) * J_INV
               U_L (2)  = F (AGP)
            ELSE
C              Adjacent network is Neumann, so PHI is to be solved
               U_L (1) = G (AGP) * J_INV
               U_L (2) = PHN (1, AGP)
            END IF
            CALL L2G (CRD (1, 1, AGP), U_L, U_G)
            CALL G2L (CRD (1, 1, JGP), U_G, U_L)
C           Set the value of phi_xi to that of adjacent network
            X_XI  = CRD (1, 2, JGP)
            Z_XI  = CRD (2, 2, JGP)
            J_INV = 1.0D0 / SQRT (X_XI ** 2 + Z_XI ** 2)
            IF (BCT .LT. 0) THEN
               PHI (2, JGP) = U_L (1) / J_INV
            ELSE IF (BCT .GT. 0) THEN
               G (JGP) = U_L (1) / J_INV
            END IF
         END DO
C         WRITE (USR_O, *) 'G AFTER SET_TANGENT'
C         WRITE (USR_O, 101) G (IGP : IGP + NGP - 1)
         IGP = IGP + NGP
      END DO
C
      RETURN
C  101 FORMAT (1X, 1E14.6)
C  111 FORMAT (1X, 4E14.6)
      END

      SUBROUTINE SET_TANGENT_2 (NNW, NW_IPAR, CRD, PHI, PHN, F, G)
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
      REAL   (KIND=RK) CRD (2, 2, *)
      REAL   (KIND=RK) PHI (2, *)
      REAL   (KIND=RK) PHN (2, *)
      REAL   (KIND=RK) F   (*)
      REAL   (KIND=RK) G   (*)
C
      INTEGER(KIND=IK) INW, IED, BCT, IGP, NGP
      INTEGER(KIND=IK) ANW, AED, ACT, AGP, AGP1, AGP2
      INTEGER(KIND=IK) JGP
      INTEGER(KIND=IK) SPT (2)
      DATA             SPT /1, 1/
C
      REAL   (KIND=RK) J_I,   J_I2, J_XI
      REAL   (KIND=RK) X_XI,  X_XIXI
      REAL   (KIND=RK) Z_XI,  Z_XIXI
      REAL   (KIND=RK) PH_XI, PH_XIXI
      REAL   (KIND=RK) PH_N,  PH_N_XI
      REAL   (KIND=RK) PH_SS, PH_NS, PH_XX, PH_XZ
      REAL   (KIND=RK) R_11, R_12
      REAL   (KIND=RK) W (4)
C
      IGP = 1
      DO INW = 1, NNW
         NGP = GET_NGP (NW_IPAR (1, INW))
         BCT = GET_BCT (NW_IPAR (1, INW))
         DO IED = 1, 2
            JGP = IGP + (NGP - 1) * (IED - 1)
            ANW = GET_ANW (IED, NW_IPAR (1, INW))
            ACT = GET_BCT (NW_IPAR (1, ANW))
            AED = GET_AED (IED, NW_IPAR (1, INW))
            AGP = GET_AGP (IED, INW, NW_IPAR)
            X_XI = CRD (1, 2, AGP)
            Z_XI = CRD (2, 2, AGP)
            J_I2 = 1.0D0 / (X_XI ** 2 + Z_XI ** 2)
            J_I  = SQRT (J_I2)
            W (1) =  12.0D0 * DBLE (AED - 1) - 6.0D0
            W (2) = -12.0D0 * DBLE (AED - 1) + 6.0D0
            W (3) =   6.0D0 * DBLE (AED - 1) - 4.0D0
            W (4) =   6.0D0 * DBLE (AED - 1) - 2.0D0
            AGP1  = AGP + 1 - AED
            AGP2  = AGP + 2 - AED
            X_XIXI = CRD (1, 1, AGP1) * W (1) + CRD (1, 1, AGP2) * W (2)
     &             + CRD (1, 2, AGP1) * W (3) + CRD (1, 2, AGP2) * W (4)
            Z_XIXI = CRD (2, 1, AGP1) * W (1) + CRD (2, 1, AGP2) * W (2)
     &             + CRD (2, 2, AGP1) * W (3) + CRD (2, 2, AGP2) * W (4)
            J_XI   = X_XI * X_XIXI + Z_XI * Z_XIXI
            IF (ACT .LT. 0) THEN
               PH_XI   = PHI (2, AGP)
               PH_XIXI = PHI (1, AGP1) * W (1) + PHI (1, AGP2) * W (2)
     &                 + PHI (2, AGP1) * W (3) + PHI (2, AGP2) * W (4)
               PH_N    = F (AGP)
               PH_N_XI = G (AGP)
            ELSE
               PH_XI   = G (AGP)
               PH_XIXI = F (AGP1) * W (1) + F (AGP2) * W (2)
     &                 + G (AGP1) * W (3) + G (AGP2) * W (4)
               PH_N    = PHN (1, AGP)
               PH_N_XI = PHN (2, AGP)
            END IF
            PH_NS =  PH_N_XI * J_I
            PH_SS = (PH_XIXI - J_XI * J_I * PH_XI) * J_I2
            R_11  = (X_XI - Z_XI) * (X_XI + Z_XI)
            R_12  =  X_XI * Z_XI  * 2.0D0
            PH_XX = (R_11 * PH_SS - R_12 * PH_NS) * J_I2
            PH_XZ = (R_11 * PH_NS + R_12 * PH_SS) * J_I2
C
            X_XI  = CRD (1, 2, JGP)
            Z_XI  = CRD (2, 2, JGP)
            J_I2  = 1.0D0 / (X_XI ** 2 + Z_XI ** 2)
            J_I   = SQRT (J_I2)
            R_11  = (X_XI - Z_XI) * (X_XI + Z_XI)
            R_12  =  X_XI * Z_XI  * 2.0D0
            PH_SS = (R_11 * PH_XX + R_12 * PH_XZ) * J_I2
            PH_NS = (R_11 * PH_XZ - R_12 * PH_XX) * J_I2
            PH_N_XI = PH_NS / J_I
C           WRITE (*, *) JGP, ':', PH_S, PH_N, PH_SS, PH_NS
C           WRITE (*, *) JGP, ':', PH_XX, PH_XZ
C           WRITE (*, *) JGP, ':', AGP1, AGP2, '(',AGP,')'
            IF (BCT .LT. 0) THEN
               G (JGP) = PH_N_XI
            END IF
         END DO
         IGP = IGP + NGP
      END DO
C
      RETURN
      END
