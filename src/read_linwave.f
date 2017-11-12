      SUBROUTINE READ_LINWAVE (PAR, PERIOD)
C ---------------------------------------------------------------------------
C     PAR (1) : Wave Amplitude
C     PAR (2) : Wave Length
C     PAR (3) : Water depth
C     PAR (4) : Wave crest offset 
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'rle_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'tme_params.inc'
C      
      REAL   (KIND=RK) PAR (*)
      REAL   (KIND=RK) PERIOD
C      
      INTEGER(KIND=IK) I
      REAL   (KIND=RK) K, H
      CHARACTER*72     LINE
C      
      CALL GET_TOKENS (USR_I, 4, LINE)
      READ (LINE, *) (PAR (I), I = 1, 4)

      K  = TWO_PI / PAR (2)
      H  = PAR (3)
      PERIOD = TWO_PI / SQRT (Grav * K * TANH (K * H))
C
      WRITE (USR_O, *) '    -> Linear wave (grid time-independent)'
      WRITE (USR_O, '(1X,A,F10.4)') '    Wave length:', PAR (2)
      WRITE (USR_O, '(1X,A,F10.4)') '    Wave period:', PERIOD
      WRITE (USR_O, '(1X,A,F10.4)') '    Wave height:', 2.0D0 * PAR (1)
      WRITE (USR_O, '(1X,A,F10.4)') '    Water depth:', H
      WRITE (USR_O, '(1X,A,F10.4)') '    Crest offs.:', PAR (4)
C      
      GRID_TIME_DEPENDENT = .FALSE.
C      
      RETURN
      END
