      SUBROUTINE CHECK_INTERFACES
     &           (NSD, NNW_SD, NW_IPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'net_funcs.inc'
C
      INTEGER(KIND=IK) NSD
      INTEGER(KIND=IK) NNW_SD (*)
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)
C
      INTEGER(KIND=IK) ISD, INW, NNW, NTF
C
      NNW = 0
      DO ISD = 1, NSD
         NNW = NNW + NNW_SD (ISD)
      END DO
      NTF = 0
      DO INW = 1, NNW
         IF (GET_ITF (NW_IPAR (1, INW)) .NE. 0) THEN
            NTF = NTF + 1
         END IF
      END DO
      IF (MOD (NTF, 2) .NE. 0) THEN
         CALL ERROR ('CHECK_INTERFACES', 'Number of ITF''s is odd')
      END IF
      CALL CHECK_INTERF (NNW, NTF, NW_IPAR)
      RETURN
      END
      SUBROUTINE CHECK_INTERF
     &           (NNW, NTF, NW_IPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'net_funcs.inc'
C
      INTEGER(KIND=IK) NNW
      INTEGER(KIND=IK) NTF
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)
C
      INTEGER(KIND=IK) INW, ITF, JTF, N_ITFS, TMP
      INTEGER(KIND=IK) ITF_ID (NTF)
      INTEGER(KIND=IK) INW_ID (NTF)
C
      N_ITFS = 0
      DO INW = 1, NNW
         ITF = GET_ITF (NW_IPAR (1, INW))
         IF (ITF .NE. 0) THEN
            N_ITFS = N_ITFS + 1
            ITF_ID (N_ITFS) = ITF
            INW_ID (N_ITFS) = INW
         END IF
      END DO
      DO ITF = 1, N_ITFS
         DO JTF = ITF + 1,  N_ITFS
            IF (ITF_ID (JTF) .LT. ITF_ID (ITF)) THEN
               TMP = ITF_ID (JTF)
               ITF_ID (JTF) = ITF_ID (ITF)
               ITF_ID (ITF) = TMP
               TMP = INW_ID (JTF)
               INW_ID (JTF) = INW_ID (ITF)
               INW_ID (ITF) = TMP
            END IF
         END DO
      END DO
      DO ITF = 1, N_ITFS, 2
C        IF (ITF_ID (ITF) .NE. ITF_ID (ITF + 1)) THEN
            WRITE (USR_O, 101) ITF_ID (ITF),     INW_ID (ITF)
            WRITE (USR_O, 101) ITF_ID (ITF + 1), INW_ID (ITF + 1)
C           CALL ERROR ('CHECK_INTERF', 'WRONG INTERFACE ID''s')
C        END IF
      END DO
C
      RETURN
  101 FORMAT (1X,'AT INTERFACE',I2,' OF NETWORK',I4)
      END