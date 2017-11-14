      SUBROUTINE ANALYTIC_T
     &           (N, T, PAR, X, W)
C ---------------------------------------------------------------------------
C     Computes an analytic solution to Laplace's equation in the grid
C     points given by X and stores it in the workspace W.
C     The kind of analytic solution is determined by the parameter
C     ANL_TYPE, which is a common block variable in the file anl_params.inc
C     The analytic solution is stored in the workspace W as follows:
C
C     W (1:N, 1) : PHI_T
C     W (1:N, 3) : PHI_XT
C     W (1:N, 4) : PHI_ZT
C     W (1:N, 5) : PHI_XXT
C     W (1:N, 6) : PHI_XZT
C
C     Due to Laplace's equation we have PHI_XX + PHI_ZZ = 0, so PHI_ZZ
C     is given implicitly as -PHI_XX
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'anl_params.inc'
C
      INTEGER(KIND=IK) N         ! Nr. of grid points  (IN)
      REAL   (KIND=RK) T         ! Time (IN)
      REAL   (KIND=RK) PAR (*)
      REAL   (KIND=RK) X (4, *)  ! coordinates of grid points (IN)
      REAL   (KIND=RK) W (*)     ! workspace, which will contain
                                 ! the analytic solution (OUT)
C
      CHARACTER*35 MSG
C
      IF (ANL_TYPE .EQ. ANL_LINWAVE) THEN
         CALL LINWAVE_T
     &        (N, T, PAR, X, W)
      ELSE IF (ANL_TYPE .EQ. ANL_RFWAVE) THEN
         CALL RFWAVE_T
     &        (N, T, PAR, X, W)
      ELSE IF (ANL_TYPE .EQ. ANL_NOTAVAIL) THEN
C        THERE'S NO KNOWN SOLUTION, SO W <- 0
         CALL NOTAVAIL
     &        (N, W)
      ELSE
         WRITE (MSG,'(A,I3)')
     &   'Analytic solution not implemented:', ANL_TYPE
         CALL ERROR ('ANALYTIC', MSG)
      END IF
C
      RETURN
      END
      SUBROUTINE LINWAVE_T
     &           (N, T, ANL_PAR, CRD, ANL)
C ---------------------------------------------------------------------------
C     Analytic solution: Linear wave
C
C              A  W  COSH (K (Z + H))
C        PHI = ----  ---------------- * SIN (K (X - X0) - W T)
C                K     SINH (K H)
C
C     A : Wave amplitude    [m.]
C     K : wave number       [1/m.]
C     W : Wave frequency    [rad/s.]
C     H : Water depth       [m.]
C     X0: Wave crest offset [m.]
C
C     ANL (1:N, 2) : PHI_T
C     ANL (1:N, 3) : PHI_X
C     ANL (1:N, 4) : PHI_Z
C     ANL (1:N, 5) : PHI_XX
C     ANL (1:N, 6) : PHI_XZ
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'rle_params.inc'
C
      INTEGER(KIND=IK) N           ! Nr. of grid points (IN)
      REAL   (KIND=RK) T           ! Time (IN)
      REAL   (KIND=RK) ANL_PAR (*)
      REAL   (KIND=RK) CRD (4, *)  ! Grid point coordinates (IN)
      REAL   (KIND=RK) ANL (N, *)  ! Analytic solution (OUT)
C
      INTEGER(KIND=IK) I
      REAL   (KIND=RK) A, L, H, X0, K, W
      REAL   (KIND=RK) X, Z, SX, CX, SHZ, CHZ, SHKH, KH, WT, AW2
C
      A  = ANL_PAR (1)
      L  = ANL_PAR (2)
      H  = ANL_PAR (3)
      X0 = ANL_PAR (4)
C
      K    = TWO_PI / L
      KH   = K * H
      W    = SQRT (Grav * K * TANH (KH))
      WT   = W * T
      AW2  = A * W ** 2
      SHKH = SINH (KH)
C
      DO I = 1, N
         X   = CRD (1, I) - X0
         Z   = CRD (2, I)
         SX  = SIN  (K * X - WT)
         CX  = COS  (K * X - WT)
         SHZ = SINH (K * Z + KH) / SHKH
         CHZ = COSH (K * Z + KH) / SHKH
         ANL (I, 1) =  AW2 / K * CHZ * SX * W
         ANL (I, 3) = -AW2     * CHZ * SX
         ANL (I, 4) =  AW2     * SHZ * CX
         ANL (I, 5) = -AW2 * K * CHZ * CX
         ANL (I, 6) = -AW2 * K * SHZ * SX
      END DO
C
      RETURN
      END
      SUBROUTINE RFWAVE_T
     &           (N, T, ANL_PAR, CRD, ANL)
C ---------------------------------------------------------------------------
C     Analytic solution: Rienecker & Fenton type of wave
C
C            __ N
C            \      COSH (J K (Z + H))
C    PHI =   /    B ---------------- * SIN (J K (X - X0 - C T))
C            --    J  COSH (J K H)
C            J = 1                                         2
C                   + (C + B) (X - X0) - (G (R - H) - 1/2 C ) T
C                           0
C
C     K : wave number            [1/m]
C     W : Wave frequency         [rad/s]
C     C : Wave phace velocity    [m/s]
C     H : Water depth            [m]
C     X0: Wave crest offset      [m]
C     B : J-th Fourier component [-]
C      J
C     G : Gravitational accel.   [m/s^2]
C     R : R&F constant           [m]
C
C     ANL (1:N, 2) : PHI_T
C     ANL (1:N, 3) : PHI_XT
C     ANL (1:N, 4) : PHI_ZT
C     ANL (1:N, 5) : PHI_XXT
C     ANL (1:N, 6) : PHI_XZT
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'rle_params.inc'
C
      INTEGER(KIND=IK) N           ! Nr. of grid points (IN)
      REAL   (KIND=RK) T           ! Time (IN)
      REAL   (KIND=RK) ANL_PAR (*)
      REAL   (KIND=RK) CRD (4, *)  ! Grid point coordinates (IN)
      REAL   (KIND=RK) ANL (N, *)  ! Analytic solution (OUT)
C
      INTEGER(KIND=IK) I, J, NRF
      REAL   (KIND=RK) L, H, X0, W, K, K2, R, B, B_0, C
      REAL   (KIND=RK) X, Z, AA, SX, CX, SHZ, CHZ, KXT, KZH
C
      L   = ANL_PAR (2)
      H   = ANL_PAR (3)
      X0  = ANL_PAR (4)
      W   = ANL_PAR (5)
      K   = ANL_PAR (6)
      K2  = K ** 2
      C   = W / K
      R   = ANL_PAR (7)
      NRF = INT (ANL_PAR (8))
      B   = (GRAV * (R - H) - 0.5D0 * C ** 2)
      B_0 = ANL_PAR (34)
C
      DO I = 1, N
         X    = CRD (1, I) - X0
         Z    = CRD (2, I)
         KXT  = K * X - W * T
         KZH  = K * (Z + H)
         ANL (I, 1) = 0.0D0
         ANL (I, 3) = 0.0D0
         ANL (I, 4) = 0.0D0
         ANL (I, 5) = 0.0D0
         ANL (I, 6) = 0.0D0
         DO J = 1, NRF
           AA  = ANL_PAR (J + 34) / COSH (J * K * H)
           SX  = SIN  (J * KXT)
           CX  = COS  (J * KXT)
           SHZ = SINH (J * KZH) * AA
           CHZ = COSH (J * KZH) * AA
           ANL (I, 3) = -ANL (I, 3) + CHZ * SX * J
           ANL (I, 4) =  ANL (I, 4) + SHZ * CX * J
           ANL (I, 5) =  ANL (I, 5) - CHZ * CX * J ** 2
           ANL (I, 6) = -ANL (I, 6) + SHZ * SX * J ** 2
         END DO
C        NOTE: for all cases we'll be calculating, the constant
C              C + B_0 will be exactly 0
         ANL (I, 1) = -W *      ANL (I, 3) + B
         ANL (I, 3) =  W * K  * ANL (I, 3)
         ANL (I, 4) =  W * K  * ANL (I, 4)
         ANL (I, 5) =  W * K2 * ANL (I, 5)
         ANL (I, 6) =  W * K2 * ANL (I, 6)
      END DO
C
      RETURN
      END