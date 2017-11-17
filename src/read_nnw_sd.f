      SUBROUTINE READ_NNW_SD
     &           (NSD, NNW_SD)
C ---------------------------------------------------------------------------
C     Reads the number of networks for each subdomain.
C     Validity checking is done in the subroutine COMPUTE_NNW_SD.
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
C
      INTEGER(KIND=IK) NSD
      INTEGER(KIND=IK) NNW_SD (*)
C
      INTEGER(KIND=IK) ISD, N
C     Read the number of networks for each subdomain
      READ (USR_I, *) (NNW_SD (ISD), ISD = 1, NSD)
C
      RETURN
      END
