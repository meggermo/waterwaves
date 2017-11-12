      PROGRAM BEM
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'max_params.inc'
C
      INTEGER(KIND=IK) NSD               ! nr. of subdomains
      INTEGER(KIND=IK) NNW               ! nr. of networks
      INTEGER(KIND=IK) NGP               ! nr. of grid points
      INTEGER(KIND=IK) NAE               ! nr. of array elements
      INTEGER(KIND=IK) NNW_SD (NSD_MAX)  ! nr. of networks per subdomain
      INTEGER(KIND=IK) NGP_SD (NSD_MAX)  ! nr. of grid points per subd.
      INTEGER(KIND=IK) NGP_NW (NNW_MAX)  ! nr. of grid points per netw.
C
      WRITE(*,'(A80)') 'BEM'
      CALL INIT
     &     (NSD, NNW, NGP, NAE, NNW_SD, NGP_SD, NGP_NW)
C
      CALL MAIN
     &     (NSD, NNW, NGP, NAE, NNW_SD, NGP_SD, NGP_NW)
C
      CALL DONE
     &     (NSD, NNW_SD)
C
      END
