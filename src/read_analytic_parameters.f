      SUBROUTINE READ_ANALYTIC_PARAMETERS
     &           (PERIOD)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'anl_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'tme_params.inc'
C
      REAL(KIND=RK) PERIOD
C
      INTEGER(KIND=IK) I
      CHARACTER*40 MSG
C     Initialize analytic solution parameters to nil
      DATA ANL_PAR /N_RPAR * 0.0D0/
C
      WRITE (USR_O, *) '========== ANALYTIC PARAMETERS ===='
      READ  (USR_I, *) ANL_TYPE
      IF (ANL_TYPE .EQ. ANL_LINWAVE) THEN
C        Linear wave
         CALL READ_LINWAVE (ANL_PAR, PERIOD)
         GRID_TIME_DEPENDENT = .FALSE.
      ELSE IF (ANL_TYPE .EQ. ANL_RFWAVE) THEN
C        Rienecker & Fenton wave
         CALL READ_RFWAVE (ANL_PAR, PERIOD)
      ELSE IF (ANL_TYPE .EQ. ANL_POLYNOMIAL) THEN
C        Polynomial of the form (X + I * Z) ** K + C * T
         WRITE (USR_O, *) '    -> Polynomial'
         READ  (USR_I, *)    (ANL_PAR (I), I = 1, 5)
         WRITE (USR_O, 101) (ANL_PAR (I), I = 1, 5)
         GRID_TIME_DEPENDENT = .FALSE.
      ELSE IF (ANL_TYPE .EQ. ANL_NOTAVAIL) THEN
C        DO NOTHING, SINCE THERE'S NO KNOWN SOLUTION
      ELSE
         WRITE (MSG,'(A,I3)')
     &   'Analytic solution not implemented:', ANL_TYPE
         CALL ERROR ('ANALYTIC', MSG)
      END IF
C
      RETURN
  101 FORMAT (8X,'PARAMETERS:',/,8E14.6)
      END
