      SUBROUTINE GET_TOKENS (U, N, TOKEN)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      INTEGER(KIND=IK) U
      INTEGER(KIND=IK) N
      CHARACTER*(*)    TOKEN
C
      INTEGER(KIND=IK) I, J, M
      INTEGER(KIND=IK) GET_TOKEN
      EXTERNAL         GET_TOKEN
      J = 1
      DO I = 1, N
         M = GET_TOKEN (U, TOKEN (J:J))
         J = J + M
         TOKEN (J:J) = ' '
         J = J + 1
      END DO
      END

      FUNCTION GET_TOKEN (U, TOKEN)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      INTEGER(KIND=IK) GET_TOKEN
      INTEGER(KIND=IK) U
      CHARACTER*(*)    TOKEN
C
      COMMON /SKIP/    C
      CHARACTER        C
      INTEGER(KIND=IK) I
      INTEGER(KIND=IK) INFO
C
      CALL SKIP_COMMENT (U, C)
C
      I = 0
      DO WHILE (C .NE.  ' ' .AND. C .NE. '\t'
     &    .AND. C .NE. '\n' .AND. C .NE. '#')
         I = I + 1
         TOKEN (I:I) = C
         CALL FGETC (U, C, INFO)
      END DO
C
      GET_TOKEN = I - 0
C
      END

      SUBROUTINE SKIP_LINE (U, S)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      INTEGER(KIND=IK) U
      CHARACTER        S
C
      INTEGER(KIND=IK) INFO
C
      DO WHILE (S .NE. '\n')
         CALL FGETC (U, S, INFO)
      END DO
C
      CALL FGETC (U, S, INFO)
C
      RETURN
      END

      SUBROUTINE SKIP_WHITE (U, S)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      INTEGER(KIND=IK) U
      CHARACTER        S
C
      INTEGER(KIND=IK) INFO
C
      DO WHILE (S .EQ. ' ' .OR. S .EQ. '\t' .OR. S .EQ. '\n')
         CALL FGETC (U, S, INFO)
      END DO
C
      RETURN
      END

      SUBROUTINE SKIP_COMMENT (U, S)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      INTEGER(KIND=IK) U
      CHARACTER        S
C
      CALL SKIP_WHITE (U, S)
C
      DO WHILE (S .EQ. '#')
         CALL SKIP_LINE  (U, S)
         CALL SKIP_WHITE (U, S)
      END DO
C
      RETURN
      END
