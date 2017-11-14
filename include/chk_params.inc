C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      LOGICAL(KIND=LK)    CHK_EQNS
      LOGICAL(KIND=LK)    CHK_COEFS_PRE
      LOGICAL(KIND=LK)    CHK_COEFS_POST
      LOGICAL(KIND=LK)    CHK_COEFS_WRITE
      LOGICAL(KIND=LK)    CHK_CONTOUR
      LOGICAL(KIND=LK)    CHK_PARAMS
C
      COMMON /CHK_PARAMS/ CHK_EQNS,
     &                    CHK_COEFS_PRE,
     &                    CHK_COEFS_POST,
     &                    CHK_COEFS_WRITE,
     &                    CHK_CONTOUR,
     &                    CHK_PARAMS
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------