      SUBROUTINE INIT
     &           (NSD, NNW, NGP, NAE, NNW_SD, NGP_SD, NGP_NW)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      INTEGER(KIND=IK) NSD              ! Nr. of subdomains     (OUT)
      INTEGER(KIND=IK) NNW              ! Nr. of networks       (OUT)
      INTEGER(KIND=IK) NGP              ! Nr. of grid points    (OUT)
      INTEGER(KIND=IK) NAE              ! Nr. of array elements (OUT)
      INTEGER(KIND=IK) NNW_SD (*)       ! Nr. of networks per subdomain (OUT)
      INTEGER(KIND=IK) NGP_SD (*)       ! Nr. of grid points per subd.  (OUT)
      INTEGER(KIND=IK) NGP_NW (*)       ! Nr. of grid points per netw,  (OUT)

C     Initialize the common block variables
      CALL INIT_COMMON_BLOCKS

C     Process the command line arguments
      CALL PROG_ARGUMENTS

C     Open some files
      CALL OPEN_OUTPUT_FILES

C     Read the number of subdomains
      CALL READ_NSD
     &     (NSD)

C     Read the number of networks for each subdomain
      CALL READ_NNW_SD
     &     (NSD, NNW_SD)

C     Compute the total number of networks
      CALL COMPUTE_NNW
     &     (NSD, NNW_SD, NNW)

C     Read the number of grid points for all networks
      CALL READ_NGP_NW
     &     (NNW, NGP_NW)

C     Compute the number of grid points of each subdomain
      CALL COMPUTE_NGP_SD
     &     (NSD, NNW_SD, NGP_NW, NGP_SD)

C     Compute the total number of grid points and number of array elements
      CALL COMPUTE_NGP_NAE
     &     (NSD, NGP_SD, NGP, NAE)
C
      RETURN
      END
