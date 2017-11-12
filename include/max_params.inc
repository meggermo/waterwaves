C ---------------------------------------------------------------------------
C     
C ---------------------------------------------------------------------------
      INTEGER(KIND=IK) NSD_MAX      ! max. nr. of subdomains
      INTEGER(KIND=IK) NNW_MAX      ! max. nr. of networks
      INTEGER(KIND=IK) NGP_MIN      ! min. nr. of grid points
      INTEGER(KIND=IK) NGP_MAX      ! max. nr. of grid points
      INTEGER(KIND=IK) ITR_MAX      ! max. nr. of iterations
      INTEGER(KIND=IK) CDF_MAX      ! max. nr. of cdf data files
C     
      PARAMETER       (NSD_MAX = 1 000)
      PARAMETER       (NNW_MAX = 4 000)
      PARAMETER       (NGP_MIN = 3)
      PARAMETER       (NGP_MAX = 1 000 000)
      PARAMETER       (ITR_MAX = 2 000)
      PARAMETER       (CDF_MAX = NNW_MAX)
C ---------------------------------------------------------------------------
C     
C ---------------------------------------------------------------------------
