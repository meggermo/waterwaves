      FUNCTION EVAL_0 (T, LDF, F)
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      REAL   (KIND=RK) EVAL_0
      REAL   (KIND=RK) T
      INTEGER(KIND=IK) LDF
      REAL   (KIND=RK) F (LDF, 2, *)
C
      EVAL_0 = (T - 1.0D0) ** 2 * (1.0D0 + 2.0D0 * T) * F (1, 1, 1)
     &       +  T ** 2          * (3.0D0 - 2.0D0 * T) * F (1, 1, 2)
     &       + (T - 1.0D0) ** 2 *  T                  * F (1, 2, 1)
     &       +  T ** 2          * (T - 1.0D0)         * F (1, 2, 2)
      RETURN
      END
      FUNCTION EVAL_1 (T, LDF, F)
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      REAL   (KIND=RK) EVAL_1
      REAL   (KIND=RK) T
      INTEGER(KIND=IK) LDF
      REAL   (KIND=RK) F (LDF, 2, *)
C
      EVAL_1 =  6.0D0 * T          * (T - 1.0D0) * F (1, 1, 1)
     &       -  6.0D0 * T          * (T - 1.0D0) * F (1, 1, 2)
     &       + (3.0D0 * T - 1.0D0) * (T - 1.0D0) * F (1, 2, 1)
     &       + (3.0D0 * T - 2.0D0) *  T          * F (1, 2, 2)
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