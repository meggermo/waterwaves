      SUBROUTINE COMPUTE_NGP_NAE
     &           (NSD, NGP_SD, NGP, NAE)
C ---------------------------------------------------------------------------
C     Computes the total number of grid points and the total number of
C     array elements
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      INTEGER(KIND=IK) NSD          ! nr. of subdomains (IN)
      INTEGER(KIND=IK) NGP_SD (*)   ! nr. of grid points per subd. (IN)
      INTEGER(KIND=IK) NGP          ! total nr. of grid points (OUT)
      INTEGER(KIND=IK) NAE          ! total nr. of array elems. (OUT)
C
      INTEGER(KIND=IK) ISD
C
      NGP = 0
      NAE = 0
      DO ISD = 1, NSD
         NGP = NGP + NGP_SD (ISD)
         NAE = NAE + NGP_SD (ISD) ** 2
      END DO
C
      RETURN
      END