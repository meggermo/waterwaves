      SUBROUTINE COPY_NGP_TO_NW_IPAR
     &           (NNW, NGP_NW, NW_IPAR)
C ---------------------------------------------------------------------------
C     Copies the number of grid points per network into the integer
C     parameter array NW_IPAR
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      INTEGER(KIND=IK) NNW               ! nr. of networks (IN)
      INTEGER(KIND=IK) NGP_NW (*)        ! nr. of gridp. per netw. (IN)
      INTEGER(KIND=IK) NW_IPAR (N_IP, *) ! integer netw. params. (OUT)
C
      INTEGER(KIND=IK) INW
C
      DO INW = 1, NNW
         CALL SET_NGP (NGP_NW (INW), NW_IPAR (1, INW))
      END DO
C
      RETURN
      END