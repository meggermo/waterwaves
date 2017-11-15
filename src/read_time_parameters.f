      SUBROUTINE READ_TIME_PARAMETERS (PERIOD)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'tme_params.inc'
      INCLUDE 'tme_funcs.inc'
C
      REAL   (KIND=RK) PERIOD
C
      REAL   (KIND=RK) TB, TE, DT
C
      WRITE (USR_O, *) '========== TIME PARAMETERS ========'
      READ  (USR_I, *) TB
      READ  (USR_I, *) TE
      READ  (USR_I, *) DT
      READ  (USR_I, *) I_PLOT_MOD
C
      CALL SET_TIME    (PERIOD * TB)
      CALL SET_TBEG    (PERIOD * TB)
      CALL SET_TEND    (PERIOD * TE)
      CALL SET_DELTA_T (PERIOD * DT)
      IF (DT .EQ. 0.0D0) THEN
         PROB_TIME_DEPENDENT = .FALSE.
         CALL SET_TIME (0.0D0)
         CALL SET_TBEG (0.0D0)
         CALL SET_TEND (1.0D0)
      ELSE
         WRITE (USR_O, 101) 'T-BEGIN', GET_TBEG    ()
         WRITE (USR_O, 101) 'T-END  ', GET_TEND    ()
         WRITE (USR_O, 101) 'DELTA-T', GET_DELTA_T ()
         WRITE (USR_O, 111) 'DT_PLT ', I_PLOT_MOD
      END IF
C
      RETURN
  101 FORMAT (1X, A8, ':', E10.2, ' s.')
  111 FORMAT (1X, A8, ':', I3, ' * DELTA-T')
      END
