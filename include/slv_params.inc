C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      INTEGER(KIND=IK) IS_ITR_MAXIT
      INTEGER(KIND=IK) IS_DOD_MAXIT
      PARAMETER       (IS_ITR_MAXIT = 1)
      PARAMETER       (IS_DOD_MAXIT = IS_ITR_MAXIT + 1)
C
      INTEGER(KIND=IK) IS_INT_ABS
      INTEGER(KIND=IK) IS_INT_REL
      INTEGER(KIND=IK) IS_ITR_TOL
      INTEGER(KIND=IK) IS_ITR_OMG
      INTEGER(KIND=IK) IS_DOD_TOL
      INTEGER(KIND=IK) IS_DOD_OMG
      PARAMETER       (IS_INT_ABS = 1)
      PARAMETER       (IS_INT_REL = IS_INT_ABS + 1)
      PARAMETER       (IS_ITR_TOL = IS_INT_REL + 1)
      PARAMETER       (IS_ITR_OMG = IS_ITR_TOL + 2)
      PARAMETER       (IS_DOD_TOL = IS_ITR_OMG + 2)
      PARAMETER       (IS_DOD_OMG = IS_DOD_TOL + 2)
C
      INTEGER(KIND=IK) NS_IPAR
      INTEGER(KIND=IK) NS_RPAR
      PARAMETER       (NS_IPAR = IS_DOD_MAXIT)
      PARAMETER       (NS_RPAR = IS_DOD_OMG + 2)
C
      INTEGER(KIND=IK)    SO_IPAR (NS_IPAR)
      REAL   (KIND=RK)    SO_RPAR (NS_RPAR)
      COMMON /SLV_PARAMS/ SO_RPAR, SO_IPAR
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------