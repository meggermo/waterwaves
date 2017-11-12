      SUBROUTINE INIT_COMMON_BLOCKS
C ---------------------------------------------------------------------------
C     Initializes common block variables:
C
C     USR_I: the user input stream is set to the standard input stream
C     USR_O: the user output stream is set to the standard output stream
C     
C     FIRST_TIME: set to TRUE, this means that we assume that we'll
C     start computations without any knowledge of previous solutions (so
C     we're not resuming previous computations)
C      
C     PROB_TIME_DEPENDENT: set to TRUE, this means that we'll assume
C     that (some) boundary conditions are time-dependent
C      
C     GRID_TIME_DEPENDENT: We assume that the problem involves a
C     time-dependent grid. Depending on the user input it can be changed
C     to time independent e.g. when we're solving linearized water
C     wave-problems.
C
C     CHK_EQNS: set to FALSE. When it is set to true, then the equations
C     solved in the routine iterate are checked.
C     
C     CHK_COEFS_PRE: set to FALSE. When set to TRUE the integral
C     equations are checked to see if the coefficients are correctly
C     computed.
C     
C     CHK_COEFS_POST: set to FALSE. When set to TRUE the integral
C     equations are checked to see if the coefficients are correctly
C     computed.
C     
C     CHK_PARAMS: set to TRUE. Not used yet
C     
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C     usr_params.inc contains the integers USR_I and USR_O
      INCLUDE 'usr_params.inc'          
C     tme_params.inc contains FIRST_TIME and GRID_INDEPENDENT
      INCLUDE 'tme_params.inc'          
C     fle_params.inc contains the integer parameters STD_I, STD_O
      INCLUDE 'fle_params.inc'          
C     chk_params contains CHECK parameters to control program checks      
      INCLUDE 'chk_params.inc'
C     cdf_params contains netCDF parameters for writing netCDF files
      INCLUDE 'cdf_params.inc'
C      
      CHARACTER     C
      COMMON /SKIP/ C
C      
      INTEGER I
C      
      USR_I = STD_I
      USR_O = STD_O
C
      I_PLOT     = 0
      I_PLOT_MOD = 1
      PROB_TIME_DEPENDENT = .TRUE.
      GRID_TIME_DEPENDENT = .TRUE.
      FIRST_TIME          = .TRUE.
C     
      CHK_EQNS        = .FALSE.
      CHK_COEFS_PRE   = .FALSE.
      CHK_COEFS_POST  = .FALSE.
      CHK_COEFS_WRITE = .FALSE.
      CHK_CONTOUR     = .TRUE.
      CHK_PARAMS      = .TRUE.
C      
      C = ' '
C
      DO I = 1, CDF_MAX
         CDF_ID (I) = -1
      END DO
C      
      RETURN
      END
