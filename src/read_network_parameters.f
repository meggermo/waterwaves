      SUBROUTINE READ_NETWORK_PARAMETERS
     &           (NSD, NNW_SD, PERIOD, NW_IPAR, NW_RPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'fle_params.inc'
C
      INTEGER(KIND=IK) NSD               ! nr. of subd. (IN)
      INTEGER(KIND=IK) NNW_SD (*)        ! nr. of netw. per subd. (IN)
      REAL   (KIND=RK) PERIOD
      INTEGER(KIND=IK) NW_IPAR (N_IP, *) ! integer netw. par (OUT)
      REAL   (KIND=RK) NW_RPAR (N_RP, *) ! real netw. par. (OUT)
C
      INTEGER(KIND=IK) ISD, INW, NNW
C
      OPEN (UNIT= STD_T, FILE = 'dim.out', STATUS = 'UNKNOWN')
      INW = 1
      DO ISD = 1, NSD
         WRITE (USR_O, 101) ISD
         NNW = NNW_SD (ISD)
         CALL READ_NETW_PARAMETERS
     &        (NNW, PERIOD, NW_IPAR (1, INW), NW_RPAR (1, INW))
         CALL CHECK_NETW_PAR
     &        (ISD, NNW,    NW_IPAR (1, INW), NW_RPAR (1, INW))
         INW = INW + NNW
      END DO
      WRITE (STD_T, *)
      CLOSE (UNIT = STD_T)
C
      CALL CHECK_INTERFACES
     &     (NSD, NNW_SD, NW_IPAR)
C
  101 FORMAT (1X,'========== SUBDOMAIN',I3,' ==========')
      RETURN
      END
      SUBROUTINE READ_NETW_PARAMETERS
     &           (NNW, PERIOD, NW_IPAR, NW_RPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'net_params.inc'
C
      INTEGER(KIND=IK) NNW               ! nr. of networks (IN)
      REAL   (KIND=RK) PERIOD
      INTEGER(KIND=IK) NW_IPAR (N_IP, *) ! integer netw. params. (OUT)
      REAL   (KIND=RK) NW_RPAR (N_RP, *) ! real netw. params. (OUT)
C
      INTEGER(KIND=IK) INW
C
      DO INW = 1, NNW
         WRITE (USR_O, 101) INW
         CALL READ_NETW_PAR (PERIOD, NW_IPAR (1, INW), NW_RPAR (1, INW))
      END DO
C
      RETURN
  101 FORMAT (1X,'---------- NETWORK',I5,' ----------')
      END