      SUBROUTINE READ_CHK_PARAMETERS
C ---------------------------------------------------------------------------
C     
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'chk_params.inc'
C
      CHARACTER*72 LINE
C
      CALL GET_TOKENS (USR_I, 1, LINE)
      READ (LINE, *) CHK_EQNS
      CALL GET_TOKENS (USR_I, 1, LINE)
      READ (LINE, *) CHK_COEFS_PRE
      CALL GET_TOKENS (USR_I, 1, LINE)
      READ (LINE, *) CHK_COEFS_POST
      CALL GET_TOKENS (USR_I, 1, LINE)
      READ (LINE, *) CHK_CONTOUR
C     
      RETURN
      END

      
