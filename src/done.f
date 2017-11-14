      SUBROUTINE DONE (NSD, NNW_SD)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      INTEGER(KIND=IK) NSD                ! nr. of subd. (IN)
      INTEGER(KIND=IK) NNW_SD (*)         ! nr. of netw. per subd. (IN)
C
      CALL CLOSE_FILES (NSD, NNW_SD)
C
      END
      SUBROUTINE CLOSE_FILES (NSD, NNW_SD)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'fle_params.inc'
C
      INTEGER(KIND=IK) NSD                ! nr. of subd. (IN)
      INTEGER(KIND=IK) NNW_SD (*)         ! nr. of netw. per subd. (IN)
C
      CLOSE (PLT_CRD)
      CLOSE (PLT_PHI)
      CLOSE (PLT_PHN)
      CLOSE (PLT_NRM)
      CLOSE (PLT_CTR)
C
      IF (USR_I .EQ. PLT_I) THEN
         CLOSE (USR_I)
      END IF
C
      IF (USR_O .EQ. PLT_O) THEN
         CLOSE (USR_O)
      END IF
C
      END