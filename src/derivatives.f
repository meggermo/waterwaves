      SUBROUTINE DERIVATIVES (PH_ONE, PN_ONE, NNW, NW_IPAR, F, G)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'net_funcs.inc'
      INCLUDE 'usr_params.inc'
C
      INTEGER(KIND=IK) PH_ONE
      INTEGER(KIND=IK) PN_ONE
      INTEGER(KIND=IK) NNW
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)
      REAL   (KIND=RK) F   (*)
      REAL   (KIND=RK) G   (*)
C
      INTEGER(KIND=IK) INW, IGP, SPT (2), NGP, BCT
C
      IGP = 1
      DO INW = 1, NNW
         NGP = GET_NGP (NW_IPAR (1, INW))
         BCT = GET_BCT (NW_IPAR (1, INW))
         IF (BCT .GT. 0) THEN
            SPT (1) = PN_ONE * GET_SPT (1, NW_IPAR (1, INW))
            SPT (2) = PN_ONE * GET_SPT (2, NW_IPAR (1, INW))
         ELSE
            SPT (1) = PH_ONE * GET_SPT (1, NW_IPAR (1, INW))
            SPT (2) = PH_ONE * GET_SPT (2, NW_IPAR (1, INW))
         END IF
C         WRITE (USR_O, *) 'G BEFORE SPLINE'
C         WRITE (USR_O, 101) G (IGP : IGP + NGP)
         CALL Spline (SPT, NGP, 1, 1, F (IGP), G (IGP))
C         WRITE (USR_O, *) ' G AFTER SPLINE'
C         WRITE (USR_O, 101) G (IGP : IGP + NGP)
         IGP = IGP + NGP
      END DO
C
C  101 FORMAT (1X, 1E14.6)
      RETURN
      END
