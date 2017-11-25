C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      INTEGER(KIND=IK) RK45_FEHLBERG     ! -> 1
      INTEGER(KIND=IK) RK45_CASHKARP     ! -> 2
      INTEGER(KIND=IK) RK45_DORMANDPRICE ! -> 3
C
      PARAMETER   (RK45_FEHLBERG     = 1)
      PARAMETER   (RK45_CASHKARP     = RK45_FEHLBERG + 1)
      PARAMETER   (RK45_DORMANDPRICE = RK45_CASHKARP + 1)
C ---------------------------------------------------------------------------
