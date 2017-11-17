      SUBROUTINE READ_NGP_NW
     &           (NNW, NGP_NW)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
C
      INTEGER(KIND=IK) NNW
      INTEGER(KIND=IK) NGP_NW (*)
C
      INTEGER(KIND=IK) INW
C
      READ (USR_I, *) (NGP_NW (INW), INW = 1, NNW)
C
      RETURN
      END
