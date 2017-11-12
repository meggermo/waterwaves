      SUBROUTINE READ_CHK_PARAMETERS
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'chk_params.inc'
C
C     CHARACTER*72 LINE
C
C    CALL GET_TOKENS (USR_I, 1, LINE)
      READ (USR_I, *) CHK_EQNS
C     CALL GET_TOKENS (USR_I, 1, LINE)
      READ (USR_I, *) CHK_COEFS_PRE
C     CALL GET_TOKENS (USR_I, 1, LINE)
      READ (USR_I, *) CHK_COEFS_POST
C     CALL GET_TOKENS (USR_I, 1, LINE)
      READ (USR_I, *) CHK_CONTOUR
C
      RETURN
      END
