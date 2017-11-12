C ---------------------------------------------------------------------------
C     
C --------------------------------------------------------------------------- 
      INTEGER(KIND=IK) ANL_NOTAVAIL   ! -> 0
      INTEGER(KIND=IK) ANL_LINWAVE    ! -> 1
      INTEGER(KIND=IK) ANL_RFWAVE     ! -> 2
      INTEGER(KIND=IK) ANL_POLYNOMIAL ! -> 3
      INTEGER(KIND=IK) ANL_WAVEMAKER  ! -> 4
C      
      PARAMETER (ANL_NOTAVAIL   = 0)               
      PARAMETER (ANL_LINWAVE    = ANL_NOTAVAIL   + 1) 
      PARAMETER (ANL_RFWAVE     = ANL_LINWAVE    + 1)
      PARAMETER (ANL_POLYNOMIAL = ANL_RFWAVE     + 1)
      PARAMETER (ANL_WAVEMAKER  = ANL_POLYNOMIAL + 1)
C
      INTEGER(KIND=IK) N_RPAR
      PARAMETER       (N_RPAR = 100)
C
      INTEGER(KIND=IK) ANL_TYPE
      REAL   (KIND=RK) ANL_PAR (N_RPAR)
C      
      COMMON /ANL_PARAMS/ ANL_PAR, ANL_TYPE
C ---------------------------------------------------------------------------
C     
C ---------------------------------------------------------------------------
