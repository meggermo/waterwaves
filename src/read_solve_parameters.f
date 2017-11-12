      SUBROUTINE READ_SOLVE_PARAMETERS (NSD)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'slv_params.inc'
      INCLUDE 'slv_funcs.inc'
C
      INTEGER(KIND=IK) NSD
C
      INTEGER(KIND=IK) MAX_IT
      REAL   (KIND=RK) EPS_ABS, EPS_REL
      REAL   (KIND=RK) TOL (2), OMG (2)
C     CHARACTER*72 LINE
C
      WRITE (USR_O, *) '========== SOLVE PARAMETERS ======='
C
C     CALL GET_TOKENS (USR_I, 2, LINE)
      READ (USR_I, *) EPS_ABS, EPS_REL
C
      CALL SET_INT_ABS (EPS_ABS)
      CALL SET_INT_REL (EPS_REL)
      WRITE (USR_O,   *) 'INTEGRATION PARAMETERS:'
      WRITE (USR_O, 101) 'INT_ABS', GET_INT_ABS ()
      WRITE (USR_O, 101) 'INT_REL', GET_INT_REL ()
C
C     CALL GET_TOKENS (USR_I, 5, LINE)
      READ (USR_I, *) TOL, OMG, MAX_IT
      CALL SET_ITR_TOL (1, TOL (1))
      CALL SET_ITR_TOL (2, TOL (2))
      CALL SET_ITR_OMG (1, OMG (1))
      CALL SET_ITR_OMG (2, OMG (2))
      CALL SET_ITR_MAXIT (MAX_IT)
      WRITE (USR_O,   *) 'ITERATION PARAMETERS:'
      WRITE (USR_O, 101) 'ITR_TOL', GET_ITR_TOL (1), GET_ITR_TOL (2)
      WRITE (USR_O, 101) 'ITR_OMG', GET_ITR_OMG (1), GET_ITR_OMG (2)
      WRITE (USR_O, 111) 'ITR_MAX', GET_ITR_MAXIT ()
C
      IF (NSD .GT. 1) THEN
C        CALL GET_TOKENS (USR_I, 5, LINE)
         READ (USR_I, *) TOL, OMG, MAX_IT
         CALL SET_DOD_TOL (1, TOL (1))
         CALL SET_DOD_TOL (2, TOL (2))
         CALL SET_DOD_OMG (1, OMG (1))
         CALL SET_DOD_OMG (2, OMG (2))
         CALL SET_DOD_MAXIT (MAX_IT)
         WRITE (USR_O,   *) 'DOMAIN DECOMPOSITION PARAMETERS:'
         WRITE (USR_O, 101) 'DOD_TOL', GET_DOD_TOL (1), GET_DOD_TOL (2)
         WRITE (USR_O, 101) 'DOD_OMG', GET_DOD_OMG (1), GET_DOD_OMG (2)
         WRITE (USR_O, 111) 'DOD_MAX', GET_DOD_MAXIT ()
      ELSE
         CALL SET_DOD_TOL (1, 0.0D0)
         CALL SET_DOD_TOL (2, 0.0D0)
         CALL SET_DOD_OMG (1, 0.0D0)
         CALL SET_DOD_OMG (2, 0.0D0)
         CALL SET_DOD_MAXIT (1)
      END IF
C
      RETURN
  101 FORMAT (1X,A8,':',2E14.4)
  111 FORMAT (1X,A8,':',I7)
      END
