      SUBROUTINE READ_PARAMETERS
     &           (NSD, NNW_SD, NW_IPAR, NW_RPAR)
C ---------------------------------------------------------------------------
C     
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      INTEGER(KIND=IK) NSD                ! nr. of subd. (IN)
      INTEGER(KIND=IK) NNW_SD (*)         ! nr. of netw. per subd. (IN)
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)  ! integer netw. par. (OUT)
      REAL(KIND=RK)    NW_RPAR (N_RP, *)  ! real netw. par. (OUT)
C
      REAL(KIND=RK) PERIOD
      DATA          PERIOD /1.0D0/
C      
      CALL READ_ANALYTIC_PARAMETERS
     &     (PERIOD)
C      
      CALL READ_TIME_PARAMETERS 
     &     (PERIOD)
C
      CALL READ_SOLVE_PARAMETERS
     &     (NSD)
C      
      CALL READ_CHK_PARAMETERS
C      
      CALL READ_NETWORK_PARAMETERS
     &     (NSD, NNW_SD, PERIOD, NW_IPAR, NW_RPAR)
C     
C     CALL OPEN_CDF_FILES (NSD, NNW_SD, NW_IPAR)
C      
      RETURN
      END
