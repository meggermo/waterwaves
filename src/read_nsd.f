      SUBROUTINE READ_NSD
     &           (NSD)
C ---------------------------------------------------------------------------
C     Reads the number of subdomains and checks if the value is within
C     the allowed range.
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'max_params.inc'
C
      INTEGER(KIND=IK) NSD
C     Read the number of subdomains
      READ (USR_I, *) NSD
C     and check if the value is withing the allowed range
      CALL I_ASSERT (1, NSD, NSD_MAX, 'NSD', 'READ_NSD')
C
      RETURN
      END
