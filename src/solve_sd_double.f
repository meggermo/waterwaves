      SUBROUTINE SOLVE_SD_DOUBLE
     &   (NNW,     NGP,    NEU,
     &    NW_IPAR,
     &    SRC_0,   SRC_1,  DIP_0,  DIP_1,
     &    OMG,     MAX_IT, TOL,
     &    CRD,     PHI ,   PHN)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'fle_params.inc'
C
      INTEGER(KIND=IK) NNW
      INTEGER(KIND=IK) NGP
      INTEGER(KIND=IK) NEU
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
      INTEGER(KIND=IK) INFO
      INTEGER(KIND=IK) IPIV (NGP * 2)
      INTEGER(KIND=IK) I, J
      REAL   (KIND=RK) A    (NGP * 2, NGP * 2)
      REAL   (KIND=RK) X    (NGP * 2)
      REAL   (KIND=RK) B    (NGP * 2)
      CHARACTER*48     MSG
C      
      CALL SETUP_AX_B_DOUBLE 
     &     (NNW,   NGP,   NW_IPAR,
     &      CRD,   PHI,   PHN,
     &      SRC_0, SRC_1, 
     &      DIP_0, DIP_1, 
     &      A,     X,     B)
C      
      CALL CHECK_EQNS_DOUBLE 
     &     (NGP, A (1,       1), A (1,       NGP + 1), 
     &           A (NGP + 1, 1), A (NGP + 1, NGP + 1), 
     &           X (1),          X (NGP + 1),
     &           B (1),          B (NGP + 1))
C      
C     CALL ITERATE_DOUBLE
C    &     (NGP,  
C    &      A (1,       1), A (1,       NGP + 1),
C    &      A (NGP + 1, 1), A (NGP + 1, NGP + 1),
C    &      B (1),          B (NGP + 1),
C    &      OMG,            MAX_IT,      TOL,
C    &      X (1),          X (NGP +1))

      OPEN (STD_T, FILE = 'A0.out', STATUS = 'UNKNOWN')
      DO I = 1, NGP
         WRITE (STD_T, 101) (A (I, J), J = 1, NGP)
      END DO
      CLOSE (STD_T)
      
      OPEN (STD_T, FILE = 'A1.out', STATUS = 'UNKNOWN')
      DO I = 1, NGP
         WRITE (STD_T, 101) (A (I, J + NGP), J = 1, NGP)
      END DO
      CLOSE (STD_T)
      
      CALL DGESV (NGP * 2, 1, A, NGP * 2, IPIV, B, NGP * 2, INFO)
C     
      IF (INFO .NE. 0) THEN
         IF (INFO .LT. 0) WRITE (MSG, 201) INFO
         IF (INFO .GT. 0) WRITE (MSG, 211) INFO
         CALL ERROR ('SOLVE_SD_DOUBLE',MSG)
      END IF
      
      CALL COPY_BACK
     &     (NNW, NGP, NW_IPAR, B, PHI, PHN)
C
      RETURN
  101 FORMAT (1X,40E10.2)
  201 FORMAT (1X,'DGESV: ARGUMENT',I2,' IS INCORRECT')
  211 FORMAT (1X,'DGESV: MATRIX BECAME SINGULAR AT ROW',I4)
      END
