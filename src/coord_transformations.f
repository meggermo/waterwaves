      SUBROUTINE GLOBAL_TO_LOCAL (N, LDU, CRD, U_G, U_L)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      INTEGER(KIND=IK) N
      INTEGER(KIND=IK) LDU
      REAL   (KIND=RK) CRD (2, 2, *)
      REAL   (KIND=RK) U_G (LDU, *)
      REAL   (KIND=RK) U_L (LDU, *)
C
      INTEGER(KIND=IK) I
C
      DO I = 1, N
        CALL G2L (CRD (1, 1, I), U_G (1, I), U_L (1, I))
      END DO
C
      END

      SUBROUTINE G2L (CRD, U_G, U_L)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "knd_params.inc"
C
      REAL(KIND=RK) CRD (2, *)
      REAL(KIND=RK) U_G (*)
      REAL(KIND=RK) U_L (*)
C
      REAL(KIND=RK) J_I
C
      J_I = 1.0D0 / SQRT (CRD (1, 2) ** 2 + CRD (2, 2) ** 2)
      U_L (1) = J_I * (CRD (1, 2) * U_G (1) + CRD (2, 2) * U_G (2))
      U_L (2) = J_I * (CRD (1, 2) * U_G (2) - CRD (2, 2) * U_G (1))
C
      END SUBROUTINE

      SUBROUTINE LOCAL_TO_GLOBAL (N, LDU, CRD, U_L, U_G)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      INTEGER(KIND=IK) N
      INTEGER(KIND=IK) LDU
      REAL   (KIND=RK) CRD (2, 2, *)
      REAL   (KIND=RK) U_L (LDU, *)
      REAL   (KIND=RK) U_G (LDU, *)
C
      INTEGER(KIND=IK) I
C
      DO I = 1, N
        CALL L2G (CRD (1, 1, I), U_L (1, I), U_G (1, I))
      END DO
C
      END

      SUBROUTINE L2G (CRD, U_L, U_G)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE "knd_params.inc"
C
      REAL(KIND=RK) CRD (2, *)
      REAL(KIND=RK) U_L (*)
      REAL(KIND=RK) U_G (*)
C
      REAL(KIND=RK) J_I
C
      J_I = 1.0D0 / SQRT (CRD (1, 2) ** 2 + CRD (2, 2) ** 2)
      U_G (1) = J_I * (CRD (1, 2) * U_L (1) - CRD (2, 2) * U_L (2))
      U_G (2) = J_I * (CRD (1, 2) * U_L (2) + CRD (2, 2) * U_L (1))
C
      END SUBROUTINE
