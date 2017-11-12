      SUBROUTINE PLTOUT 
     &           (NSD,   NNW_SD, NW_IPAR, 
     &            CRD,   PHI,    PHN,
     &            CRD_T, PHI_T,  PHN_T)
C ---------------------------------------------------------------------------
C    
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'fle_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'net_funcs.inc'
      INCLUDE 'tme_funcs.inc'
C      
      INTEGER(KIND=IK) NSD
      INTEGER(KIND=IK) NNW_SD (*)
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)
      REAL   (KIND=RK) CRD (4, *)
      REAL   (KIND=RK) PHI (2, *)
      REAL   (KIND=RK) PHN (2, *)
      REAL   (KIND=RK) CRD_T (4, *)
      REAL   (KIND=RK) PHI_T (2, *)
      REAL   (KIND=RK) PHN_T (2, *)
C      
      INTEGER(KIND=IK) ISD, INW, JNW, NNW, IGP, NGP, PLT, BCT, CDF
      REAL   (KIND=RK) T, PERIOD
      REAL   (KIND=RK) NRM (6)
C
      PERIOD = GET_PERIOD ()
      T      = GET_TIME   ()
      WRITE (PLT_CRD, '(''# T = '',F8.4, '' * PERIOD'')') T / PERIOD
      WRITE (PLT_PHI, '(''# T = '',F8.4, '' * PERIOD'')') T / PERIOD
      WRITE (PLT_PHN, '(''# T = '',F8.4, '' * PERIOD'')') T / PERIOD
      WRITE (PLT_NRM, 101) T
C      
      INW = 1
      IGP = 1
C     
      DO ISD = 1, NSD
         NNW = NNW_SD (ISD)
         DO JNW = INW, INW + NNW - 1
            NGP = GET_NGP (NW_IPAR (1, JNW))
            BCT = GET_BCT (NW_IPAR (1, JNW))
            PLT = GET_PLT (NW_IPAR (1, JNW))
            IF (PLT .NE. 0) THEN
               CALL PLT_NW 
     &              (NGP, ABS (BCT), PLT,
     &               T, 
     &               CRD   (1, IGP), PHI   (1, IGP), PHN   (1, IGP),
     &               CRD_T (1, IGP), PHI_T (1, IGP), PHN_T (1, IGP),
     &               NRM)
               WRITE (PLT_NRM, 121) T / PERIOD, NRM
            END IF
            IGP = IGP + NGP
         END DO
         INW = INW + NNW
      END DO
      WRITE (PLT_NRM, *)
C
      RETURN
  101 FORMAT (1X, '#TIME =', E14.6)
  111 FORMAT (1X, 'INW ', 6A10)
  121 FORMAT (1X, E14.6,10E12.4)
      END

      SUBROUTINE PLT_NW 
     &           (NGP,   BCT,   PLT, 
     &            T, 
     &            CRD,   PHI,   PHN, 
     &            CRD_T, PHI_T, PHN_T, NRM)
C ---------------------------------------------------------------------------
C    
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'anl_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'bct_params.inc'
      INCLUDE 'fle_params.inc'
C      
      INTEGER(KIND=IK) NGP
      INTEGER(KIND=IK) BCT
      INTEGER(KIND=IK) PLT
      REAL   (KIND=RK) T
      REAL   (KIND=RK) CRD (4, NGP)
      REAL   (KIND=RK) PHI (2, *)
      REAL   (KIND=RK) PHN (2, *)
      REAL   (KIND=RK) CRD_T (4, *)
      REAL   (KIND=RK) PHI_T (2, *)
      REAL   (KIND=RK) PHN_T (2, *)
      REAL   (KIND=RK) NRM   (3, *)
C
      INTEGER(KIND=IK) I, J
      REAL   (KIND=RK) X_XI, Z_XI, J_I, PROD, DIFF
      REAL   (KIND=RK) W (NGP, 8)
C
      INTEGER(KIND=IK) IDAMAX
      EXTERNAL         IDAMAX
      REAL   (KIND=RK) DNRM2, DASUM
      EXTERNAL         DNRM2, DASUM
C      
      IF (BCT .EQ. BCT_LINBERNOU) THEN
         DO I = 1, NGP
            W   (I, 8) = CRD (2, I)
            CRD (2, I) = 0.0D0
         END DO
      END IF
      CALL ANALYTIC (NGP, T, ANL_PAR, CRD, W)
      IF (BCT .EQ. BCT_LINBERNOU) THEN
         DO I = 1, NGP
            CRD (2, I) = W (I, 8)
         END DO
      END IF
C      
       DO I = 1, NGP
         X_XI  = CRD (3, I)
         Z_XI  = CRD (4, I)
         J_I   = 1.0D0 / SQRT (X_XI ** 2 + Z_XI ** 2)
         PROD  = 2.0D0 * X_XI * Z_XI
         DIFF  = (X_XI + Z_XI) * (X_XI - Z_XI)
         W (I, 1) = PHI (1, I) -  W (I,1)
         W (I, 2) = PHI (2, I) - (W (I,3) * X_XI + W (I,4) * Z_XI)
         W (I, 3) = PHN (1, I) - (W (I,4) * X_XI - W (I,3) * Z_XI) * J_I
         W (I, 4) = PHN (2, I) - (W (I,6) * DIFF - W (I,5) * PROD) * J_I
         W (I, 5) = CRD (2, I) -  W (I,7)
         W (I, 6) = ABS (W (I, 1))
         W (I, 7) = ABS (W (I, 3))
         WRITE (PLT_CRD, 1001) (CRD (J,I), J = 1, 4), W (I,5)
         WRITE (PLT_PHI, 1001)  PHI (1,I), PHI (2,I), W (I,1), W (I,2)
         WRITE (PLT_PHN, 1001)  PHN (1,I), PHN (2,I), W (I,3), W (I,4)
      END DO
      WRITE (PLT_CRD, *)
      WRITE (PLT_PHI, *)
      WRITE (PLT_PHN, *)
C     
      NRM (1, 1) = W (IDAMAX (NGP, W (1, 6), 1), 6)
      NRM (2, 1) =     DASUM (NGP, W (1, 6), 1)
      NRM (3, 1) =     DNRM2 (NGP, W (1, 6), 1)
      NRM (1, 2) = W (IDAMAX (NGP, W (1, 7), 1), 7)
      NRM (2, 2) =     DASUM (NGP, W (1, 7), 1)
      NRM (3, 2) =     DNRM2 (NGP, W (1, 7), 1)
C
      RETURN
 1001 FORMAT (1X,  10E16.8)
      END
