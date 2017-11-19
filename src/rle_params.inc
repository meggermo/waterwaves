C-----------------------------------------------------------------------------
C     REAL PARAMETER DEFINITIONS
C-----------------------------------------------------------------------------
C     PHYSICAL PARAMETERS
      REAL(KIND=RK) Grav
      REAL(KIND=RK) WDens
      REAL(KIND=RK) ADens
C     gravitational accelleration in [m/s^2]
      PARAMETER (Grav  =  9.81D+0)
C     water density in [kg/m^3]
      PARAMETER (WDens =  1.00D+3)
C     air density in [kg/m^3]
      PARAMETER (ADens =  1.00D-1)
C     MATHEMATICAL PARAMETERS
      REAL(KIND=RK) PI
      REAL(KIND=RK) TWO_PI
      PARAMETER
     &(PI    = 3.14159265358979323846264338327950288419716939937510D+00)
      PARAMETER
     &(TWO_PI= 6.28318530717958647692528676655900576839433879875020D+00)
C
      REAL(KIND=RK) ONE_THIRD
      REAL(KIND=RK) ONE_SIXTH
C
      PARAMETER    (ONE_THIRD = 1.0D0 / 3.0D0)
      PARAMETER    (ONE_SIXTH = 1.0D0 / 6.0D0)
C-----------------------------------------------------------------------------
C     END OF REAL PARAMETER DEFINITIONS
C-----------------------------------------------------------------------------