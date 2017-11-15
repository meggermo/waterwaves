      SUBROUTINE COMPUTE_NNW
     &           (NSD, NNW_SD, NNW)
C ---------------------------------------------------------------------------
C     Computes the total number of networks and checks if the input is
C     within the allowed range.
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'max_params.inc'
C
      INTEGER(KIND=IK) NSD              ! nr. of subdomains (IN)
      INTEGER(KIND=IK) NNW_SD (*)       ! nr. of netw. per  subd. (IN)
      INTEGER(KIND=IK) NNW              ! total nr. of netw. (OUT)
C
      INTEGER(KIND=IK) ISD, NNS
C
      NNW = 0
      DO ISD = 1, NSD
         NNS = NNW_SD (ISD)
         CALL I_ASSERT (1, NNS, NNW_MAX, 'NNS', 'COMPUTE_NNW')
         NNW = NNW + NNW_SD (ISD)
         CALL I_ASSERT (1, NNW, NNW_MAX, 'NNW', 'COMPUTE_NNW')
      END DO
C
      RETURN
      END
