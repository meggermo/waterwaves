C ---------------------------------------------------------------------------
C     Contains the kind parameters for the base types.
C     On GNU fortran compilers the meaning is given as follows:
C
C     KIND = 1  : the default kind on the system
C     KIND = 2  : double the size of the default kind
C
C ---------------------------------------------------------------------------
      INTEGER    CK       ! -> the character kind
      INTEGER    LK       ! -> the logical kind
      INTEGER    IK       ! -> the integer kind
      INTEGER    RK       ! -> the real kind
C
      PARAMETER (CK = 1)
      PARAMETER (LK = 1)
      PARAMETER (IK = 4)
      PARAMETER (RK = 8)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
