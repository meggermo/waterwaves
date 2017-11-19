C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      INTEGER(KIND=IK) GRD_LINE         ! -> 1
      INTEGER(KIND=IK) GRD_STRETCHED    ! -> 2
      INTEGER(KIND=IK) GRD_CIRCLESECT   ! -> 3
      INTEGER(KIND=IK) GRD_CUBIC        ! -> 4
      INTEGER(KIND=IK) GRD_FROMFILE     ! -> 5
      INTEGER(KIND=IK) GRD_FIRST
      INTEGER(KIND=IK) GRD_LAST
C
      PARAMETER   (GRD_LINE       = 1)
      PARAMETER   (GRD_FIRST      = GRD_LINE)
      PARAMETER   (GRD_STRETCHED  = 2)
      PARAMETER   (GRD_CIRCLESECT = 3)
      PARAMETER   (GRD_CUBIC      = 4)
      PARAMETER   (GRD_FROMFILE   = 5)
      PARAMETER   (GRD_LAST       = GRD_FROMFILE)
C ---------------------------------------------------------------------------
C     KIND OF GRID:  NR. OF PARAMS.  MEANING OF PARAMS.
C     -------------  --------------  ------------------
C     GRD_LINE        4*DBLE         (x,z) start point
C                                    (x,z) end point
C
C     GRD_STRETCHED   6*DBLE         (x,z) start point
C                                    (x,z) end point
C                                    K : period
C                                    K = 1, one sided stretching/compression
C                                    K = 2, double sided str/compr
C                                    A : amplitude
C                                    A > 0 compression at begin of grid
C                                    A < 0 compression at end of grid
C
C     GRD_CIRCLESECT  6*DBLE         (x,z) centre of circle
C                                    (rx,rz) x- and z-radius
C                                    (theta_begin,theta_end)
C
C     GRD_FROMFILE    1*STRING       Filename
C
C ---------------------------------------------------------------------------