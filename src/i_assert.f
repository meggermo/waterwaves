      SUBROUTINE I_ASSERT
     &   (I_LO, I_VALUE, I_HI, VAL_STR, SUB_STR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      INTEGER(KIND=IK) I_LO
      INTEGER(KIND=IK) I_VALUE
      INTEGER(KIND=IK) I_HI
      CHARACTER*(*) VAL_STR
      CHARACTER*(*) SUB_STR
C
      IF (I_LO .GT. I_VALUE .OR. I_VALUE .GT. I_HI) THEN
         WRITE (*, *) '*** ASSERTION FAILURE ***'
         WRITE (*, *) VAL_STR, ' is out of range ',
     &      I_LO, ' .. ', I_HI
         WRITE (*, *) VAL_STR, ' = ', I_VALUE
         WRITE (*, *) '*** ASSERTION FAILURE ***'
         CALL ERROR (SUB_STR, 'Value out of range')
      END IF
C
      RETURN
      END