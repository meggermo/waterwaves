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
      CHARACTER*72     LINE
C     Read the number of subdomains
      CALL GET_TOKENS (USR_I, 1, LINE) 
      READ (LINE, *) NSD
C     and check if the value is withing the allowed range
      CALL I_ASSERT (1, NSD, NSD_MAX, 'NSD', 'READ_NSD')
C      
      RETURN
      END
