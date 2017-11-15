      SUBROUTINE ERROR (SUB_STR, MSG_STR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
C
      CHARACTER*(*) SUB_STR
      CHARACTER*(*) MSG_STR
C
      WRITE (*, *) 'ERROR IN ROUTINE ', SUB_STR, ': ', MSG_STR
      WRITE (*, *) '*** PROGRAM EXITED ***'
      CALL EXIT (2)
C
      END
