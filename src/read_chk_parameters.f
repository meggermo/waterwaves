      SUBROUTINE READ_CHK_PARAMETERS
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'chk_params.inc'
C
      READ (USR_I, *) CHK_EQNS
      READ (USR_I, *) CHK_COEFS_PRE
      READ (USR_I, *) CHK_COEFS_POST
      READ (USR_I, *) CHK_CONTOUR
C
      RETURN
      END
