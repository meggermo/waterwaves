C ---------------------------------------------------------------------------
C     Kind of boundary conditions
C     ---------------------------
C     Dirichlet : PHI is given PHN is unknown
C     Mixed     : PHI and PHN are unknown
C     Neumann   : PHN is given PHI is unknown
C
C     Boundary condition types:
C     -------------------------
C     BCT_UNKNOWN   :
C     BCT_ANALYTIC  :     BC. is given analytically
C     BCT_BERNOULLI :     phi_t + 1/2 grad (phi) ** 2 + g z = 0
C     BCT_SOMMERFELD:     phi_t + c phi_x = 0
C     BCT_LINBERNOU :     phi_t + g z = 0, eta_t - phi_z = 0
C
C ---------------------------------------------------------------------------
      INTEGER(KIND=IK) BCK_DIRICHLET   ! -> -1
      INTEGER(KIND=IK) BCK_MIXED       ! ->  0
      INTEGER(KIND=IK) BCK_NEUMANN     ! ->  1
C
      PARAMETER       (BCK_DIRICHLET = -1)
      PARAMETER       (BCK_MIXED     =  0)
      PARAMETER       (BCK_NEUMANN   =  1)
C
      INTEGER(KIND=IK) BCT_UNKNOWN      ! -> 0
      INTEGER(KIND=IK) BCT_ANALYTIC     ! -> 1
      INTEGER(KIND=IK) BCT_BERNOULLI    ! -> 2
      INTEGER(KIND=IK) BCT_SOMMERFELD   ! -> 3
      INTEGER(KIND=IK) BCT_WAVEMAKER    ! -> 4
      INTEGER(KIND=IK) BCT_LINBERNOU    ! -> 5
      INTEGER(KIND=IK) BCT_TIME_POLY    ! -> 6
      INTEGER(KIND=IK) BCT_FIRST
      INTEGER(KIND=IK) BCT_LAST
C
      PARAMETER       (BCT_UNKNOWN    = 0)
      PARAMETER       (BCT_FIRST      = BCT_UNKNOWN)
      PARAMETER       (BCT_ANALYTIC   = BCT_UNKNOWN    + 1)
      PARAMETER       (BCT_BERNOULLI  = BCT_ANALYTIC   + 1)
      PARAMETER       (BCT_SOMMERFELD = BCT_BERNOULLI  + 1)
      PARAMETER       (BCT_WAVEMAKER  = BCT_SOMMERFELD + 1)
      PARAMETER       (BCT_LINBERNOU  = BCT_WAVEMAKER  + 1)
      PARAMETER       (BCT_TIME_POLY  = BCT_LINBERNOU  + 1)
      PARAMETER       (BCT_LAST       = BCT_TIME_POLY)
C
      INTEGER(KIND=IK) BCD_PHI
      INTEGER(KIND=IK) BCD_PHI_T
C
      PARAMETER       (BCD_PHI   = 1)
      PARAMETER       (BCD_PHI_T = 2)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
