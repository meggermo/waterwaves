      SUBROUTINE GLOBAL_TO_LOCAL (N, LDU, CRD, U_G, U_L)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      INTEGER(KIND=IK) N
      INTEGER(KIND=IK) LDU
      REAL   (KIND=RK) CRD (4, *)
      REAL   (KIND=RK) U_G (LDU, *)
      REAL   (KIND=RK) U_L (LDU, *)
C
      INTEGER(KIND=IK) I
      REAL   (KIND=RK) X_XI, Z_XI, J_INV, U_X, U_Z
C
      DO I = 1, N
         X_XI  = CRD (3, I)
         Z_XI  = CRD (4, I)
         J_INV = 1.0D0 / SQRT (X_XI ** 2 + Z_XI ** 2)
         U_X   = U_G (1, I)
         U_Z   = U_G (2, I)
         U_L (1, I) = (X_XI * U_X + Z_XI * U_Z) * J_INV
         U_L (2, I) = (X_XI * U_Z - Z_XI * U_X) * J_INV
      END DO
C
      RETURN
      END