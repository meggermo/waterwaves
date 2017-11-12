      SUBROUTINE READ_RFWAVE (PAR, PERIOD)
C-----------------------------------------------------------------------------
C     Reads a the parameters of a Rienecker & Fenton typ of wave.
C-----------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'fle_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'rle_params.inc'
C      
      REAL(KIND=RK) PAR (33, 3)
      REAL(KIND=RK) PERIOD
C
      LOGICAL(KIND=IK) ldummy
      INTEGER(KIND=IK) i, n
      CHARACTER*32     File_Name
      REAL   (KIND=RK) waterdepth, waveheigth, rho, gravity,
     &                 phasevelocity, Q, R, rdummy, omega, wavelength
      CHARACTER*72 LINE

      WRITE (USR_O, *)              '    -> R & F wave'
      CALL GET_TOKENS (USR_I, 1, LINE)
      READ  (LINE,  *)     File_Name
      OPEN  (STD_T, FILE = File_Name, STATUS= 'OLD')
      WRITE (USR_O, '(1X,A,A    )') '    Data file:  ', File_Name
      READ  (STD_T, *) waterdepth
      WRITE (USR_O, '(1X,A,F10.4)') '    Water depth:', waterdepth
      READ  (STD_T, *) PERIOD
      WRITE (USR_O, '(1X,A,F10.4)') '    Wave period:', period
      READ  (STD_T, *) waveheigth
      WRITE (USR_O, '(1X,A,F10.4)') '    Wave height:', waveheigth
      READ  (STD_T, *) ldummy
      READ  (STD_T, *) rdummy
      READ  (STD_T, *) n
      WRITE (USR_O, '(1X,A,  I10)') '    Wave  modes:', n
      READ  (STD_T, *) rho
      READ  (STD_T, *) gravity
      READ  (STD_T, *) phasevelocity
      READ  (STD_T, *) Q
      READ  (STD_T, *) R
      omega     = 2.0D0 * PI / PERIOD
      wavelength= period * phasevelocity
      par (1, 1) = waveheigth
      par (2, 1) = wavelength
      par (3, 1) = waterdepth
      CALL GET_TOKENS (USR_I, 1, LINE)
      READ (LINE, *)
     +par (4, 1)
      WRITE (USR_O, '(1X,A,F10.4)') '    Crest offs.:', par (4,1)
      par (5, 1) = omega
      par (6, 1) = TWO_PI / period / phasevelocity
      par (7, 1) = R
      par (8, 1) = n
      DO i= 1, n + 1
        READ (STD_T, *) rdummy, PAR (i, 3), PAR (i, 2)
      END DO
      CLOSE (STD_T)
C
      RETURN
      END
