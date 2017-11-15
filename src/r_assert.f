      SUBROUTINE R_ASSERT
     &           (D_LO, D_VALUE, D_HI, VAL_STR, SUB_STR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      REAL(KIND=RK) D_LO
      REAL(KIND=RK) D_VALUE
      REAL(KIND=RK) D_HI
      CHARACTER*(*) VAL_STR
      CHARACTER*(*) SUB_STR
C
      IF (D_LO .GT. D_VALUE .OR. D_VALUE .GT. D_HI) THEN
         WRITE (*, *) '*** ASSERTION FAILURE ***'
         WRITE (*, *) VAL_STR, ' is out of range ',
     &      D_LO, ' .. ', D_HI
         WRITE (*, *) VAL_STR, ' = ', D_VALUE
         WRITE (*, *) '*** ASSERTION FAILURE ***'
         CALL ERROR (SUB_STR, 'Value out of range')
      END IF
C
      RETURN
      END
