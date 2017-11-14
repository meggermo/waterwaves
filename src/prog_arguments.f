      SUBROUTINE PROG_ARGUMENTS
C ---------------------------------------------------------------------------
C     Processed the program arguments:
C
C     ARG 1: The name of the user input file, if given, then USR_I is
C     redirected to this file instead of the standard input stream.
C
C     ARG 2: The name of the user output file, if given, then USR_O is
C     redirected to this file instead of the standard output stream.
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'fle_params.inc'
      INCLUDE 'chk_params.inc'
C
      CHARACTER*128 FILE_NAME
      CHARACTER*8   CHK
C
      INTEGER(KIND=IK) COMMAND_ARGUMENT_COUNT
      INTRINSIC        COMMAND_ARGUMENT_COUNT
C
      IF (COMMAND_ARGUMENT_COUNT () .GT. 0) THEN
         USR_I = PLT_I
         CALL GETARG (1, FILE_NAME)
         OPEN (USR_I, FILE = FILE_NAME, STATUS = 'OLD')
      END IF
C
      IF (COMMAND_ARGUMENT_COUNT () .GT. 1) THEN
         CALL GETARG (2, FILE_NAME)
         IF (FILE_NAME (1:3) .NE. 'CHK') THEN
            USR_O = PLT_O
            OPEN (USR_O, FILE = FILE_NAME, STATUS = 'UNKNOWN')
         ELSE
            CHK_EQNS       = .TRUE.
            CHK_COEFS_PRE  = .TRUE.
            CHK_COEFS_POST = .TRUE.
         END IF
      END IF
C
      IF (COMMAND_ARGUMENT_COUNT () .GT. 2) THEN
         CALL GETARG (3, CHK)
         IF (CHK(1:3) .EQ. 'CHK') THEN
            CHK_EQNS       = .TRUE.
            CHK_COEFS_PRE  = .TRUE.
            CHK_COEFS_POST = .TRUE.
         END IF
      END IF
C
      END