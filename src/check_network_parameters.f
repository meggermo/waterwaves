      SUBROUTINE CHECK_NETW_PAR
     &           (ISD, NNW, NW_IPAR, NW_RPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      INTEGER(KIND=IK) ISD
      INTEGER(KIND=IK) NNW
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)
      REAL   (KIND=RK) NW_RPAR (N_RP, *)
C
      INTEGER(KIND=IK) INW, IED
C
      DO INW = 1, NNW
         DO IED = 1, 2
            CALL CHECK_CONNECTIONS
     &           (ISD, INW, IED, NNW, NW_IPAR, NW_RPAR)
         END DO
      END DO
C
      RETURN
      END
      SUBROUTINE CHECK_CONNECTIONS
     &           (ISD, INW, IED, NNW, NW_IPAR, NW_RPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'net_funcs.inc'
C
      INTEGER(KIND=IK) ISD
      INTEGER(KIND=IK) INW
      INTEGER(KIND=IK) IED
      INTEGER(KIND=IK) NNW
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)
      REAL   (KIND=RK) NW_RPAR (N_RP, *)
C
      INTEGER(KIND=IK) ANW, AED, BNW, BED
      REAL   (KIND=RK) PAR_I (N_GRP)
      REAL   (KIND=RK) PAR_A (N_GRP)
C
      ANW = GET_ANW (IED, NW_IPAR (1, INW))
      IF (ANW .GT. NNW) THEN
         WRITE (USR_O, 101) ISD
         WRITE (USR_O, 111) IED, INW, ANW
         CALL ERROR ('CHECK_NETW_IPAR', 'Network ID out of range')
      END IF
      AED = GET_AED (IED, NW_IPAR (1, INW))
      BNW = GET_ANW (AED, NW_IPAR (1, ANW))
      BED = GET_AED (AED, NW_IPAR (1, ANW))
      IF (INW .NE. BNW) THEN
         WRITE (USR_O, 101) ISD
         WRITE (USR_O, 121) IED, INW, ANW
         WRITE (USR_O, 121) AED, ANW, BNW
         CALL ERROR ('CHECK_NETW_IPAR', 'Network ID mismatch')
      END IF
      IF (IED .NE. BED) THEN
         WRITE (USR_O, 101) ISD
         WRITE (USR_O, 121) IED, INW, ANW
         WRITE (USR_O, 121) AED, ANW, BNW
         CALL ERROR ('CHECK_NETW_IPAR', 'Edge ID mismatch')
      END IF
      CALL GET_GRP (NW_RPAR (1, INW), PAR_I)
      CALL GET_GRP (NW_RPAR (1, ANW), PAR_A)
      IF (NNW .NE. 1) THEN
      IF (PAR_I (1 + 2 * (IED - 1)) .NE. PAR_A (1 + 2 * (AED -1)) .OR.
     &    PAR_I (2 + 2 * (IED - 1)) .NE. PAR_A (2 + 2 * (AED -1))) THEN
         WRITE (USR_O, 131) IED, INW, PAR_I (1 + 2 * (IED -1)),
     &                                PAR_I (2 + 2 * (IED -1))
         WRITE (USR_O, *) 'BUT'
         WRITE (USR_O, 131) AED, ANW, PAR_A (1 + 2 * (AED -1)),
     &                                PAR_A (2 + 2 * (AED -1))
         CALL ERROR ('CHECK_NETW_IPAR', 'Coordinate mismatch')
      END IF
      END IF
C
      RETURN
  101 FORMAT ('IN SUBDOMAIN',I2)
  111 FORMAT ('EDGE',I2,' OF NETW',I4,' CANNOT CONNECT TO NETW',I4)
  121 FORMAT ('EDGE',I2,' OF NETW',I4,' CONNECTS TO NETW',I4)
  131 FORMAT ('EDGE',I2,' OF NETW',I4,' HAS COORDS',10E14.6)
      END
