      SUBROUTINE COMPUTE_NGP_SD
     &           (NSD, NNW_SD, NGP_NW, NGP_SD)
C ---------------------------------------------------------------------------
C     Computes the number of grid points per subdomain and checks if the
C     given values are within the valid range.
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'max_params.inc'
C
      INTEGER(KIND=IK) NSD              ! nr. of subdomains (IN)
      INTEGER(KIND=IK) NNW_SD (*)       ! nr. of netw. per subd. (IN)
      INTEGER(KIND=IK) NGP_NW (*)       ! nr. of gridp. per netw. (IN)
      INTEGER(KIND=IK) NGP_SD (*)       ! nr. of gridp. per subd. (OUT)
C
      INTEGER(KIND=IK) ISD, INW, JNW, NNW, NGP, NGN
C      
      INW = 1
      DO ISD = 1, NSD
         NNW = NNW_SD (ISD)
         NGP = 0
         DO JNW = INW, INW + NNW - 1
            NGN = NGP_NW (JNW)
            CALL I_ASSERT (1, NGN, NGP_MAX, 'NGN', 'COMPUTE_NGP_SD')
            NGP = NGP + NGN
            CALL I_ASSERT (1, NGP, NGP_MAX, 'NGP', 'COMPUTE_NGP_SD')
         END DO
         NGP_SD (ISD) = NGP
         INW = INW + NNW
      END DO
C     
      RETURN
      END
