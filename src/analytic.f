      SUBROUTINE ANALYTIC 
     &           (N, T, PAR, X, W)
C ---------------------------------------------------------------------------
C     Computes an analytic solution to Laplace's equation in the grid 
C     points given by X and stores it in the workspace W. 
C     The kind of analytic solution is determined by the parameter 
C     ANL_TYPE, which is a common block variable in the file anl_params.inc
C     The analytic solution is stored in the workspace W as follows:
C      
C     W (1:N, 1) : PHI
C     W (1:N, 2) : PHI_T
C     W (1:N, 3) : PHI_X
C     W (1:N, 4) : PHI_Z
C     W (1:N, 5) : PHI_XX 
C     W (1:N, 6) : PHI_XZ
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
         CALL LINWAVE
     &        (N, T, PAR, X, W)
      ELSE IF (ANL_TYPE .EQ. ANL_RFWAVE) THEN
         CALL RFWAVE
     &        (N, T, PAR, X, W)
      ELSE IF (ANL_TYPE .EQ. ANL_POLYNOMIAL) THEN
         CALL POLYNOMIAL
     &        (N, T, PAR, X, W)
      ELSE IF (ANL_TYPE .EQ. ANL_WAVEMAKER) THEN
         CALL WAVEMAKER
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

      SUBROUTINE NOTAVAIL 
     &           (N, ANL)
C ---------------------------------------------------------------------------
C     analytic solution : PHI = 0
C      
C     ANL (1:N, 1) : PHI
C     ANL (1:N, 2) : PHI_T
C     ANL (1:N, 3) : PHI_X
C     ANL (1:N, 4) : PHI_Z
C     ANL (1:N, 5) : PHI_XX (PHI_ZZ == -PHI_XX, due to Laplace's eqn)
C     ANL (1:N, 6) : PHI_XZ
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      INTEGER(KIND=IK) N           ! Nr. of grid points (IN)
      REAL   (KIND=RK) ANL (N, *)  ! Analytic solution (OUT)
C
      INTEGER(KIND=IK) I, J
C
      DO J = 1, 6
         DO I = 1, N
            ANL (I, J) = 0.0D0
         END DO
      END DO
C
      RETURN
      END
      
      SUBROUTINE LINWAVE 
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
C     ANL (1:N, 1) : PHI
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
      REAL   (KIND=RK) X, Z, SX, CX, SHZ, CHZ, SHKH, KH, WT, AW
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
      AW   = A * W
      SHKH = SINH (KH)
C
      DO I = 1, N
         X   = CRD (1, I) - X0
         Z   = CRD (2, I)
         SX  = SIN  (K * X - WT)
         CX  = COS  (K * X - WT)
         SHZ = SINH (K * Z + KH) / SHKH
         CHZ = COSH (K * Z + KH) / SHKH
         ANL (I, 1)  = AW / K * CHZ * SX
         ANL (I, 2) = -AW / K * CHZ * CX * W
         ANL (I, 3) =  AW     * CHZ * CX
         ANL (I, 4) =  AW     * SHZ * SX
         ANL (I, 5) = -AW * K * CHZ * SX
         ANL (I, 6) =  AW * K * SHZ * CX
         ANL (I, 7) =  A  * CX
      END DO
C
      RETURN
      END

      SUBROUTINE RFWAVE 
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
C     ANL (1:N, 1) : PHI
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
      B   = -(GRAV * (R - H) - 0.5D0 * C ** 2)
      B_0 = ANL_PAR (34)
C
      DO I = 1, N
         X    = CRD (1, I) - X0
         Z    = CRD (2, I)
         KXT  = K * X - W * T
         KZH  = K * (Z + H)
         ANL (I, 1) = 0.0D0
         ANL (I, 2) = 0.0D0
         ANL (I, 3) = 0.0D0
         ANL (I, 4) = 0.0D0
         ANL (I, 5) = 0.0D0
         ANL (I, 6) = 0.0D0
         ANL (I, 7) = 0.0D0
         DO J = 1, NRF
           AA  = ANL_PAR (J + 34) / COSH (J * K * H)
           SX  = SIN  (J * KXT)
           CX  = COS  (J * KXT)
           SHZ = SINH (J * KZH) * AA
           CHZ = COSH (J * KZH) * AA
           ANL (I, 1) = ANL (I, 1) + CHZ * SX
           ANL (I, 3) = ANL (I, 3) + CHZ * CX * J
           ANL (I, 4) = ANL (I, 4) + SHZ * SX * J
           ANL (I, 5) = ANL (I, 5) - CHZ * SX * J ** 2
           ANL (I, 6) = ANL (I, 6) + SHZ * CX * J ** 2
           ANL (I, 7) = ANL (I, 7) + ANL_PAR (J + 67) * CX
         END DO
C        NOTE: for all cases we'll be calculating, the constant
C              C + B_0 will be exactly 0
         ANL (I, 1) =       ANL (I, 1) + (C + B_0) * X + B * T
         ANL (I, 2) = -W  * ANL (I, 3)                 + B
         ANL (I, 3) =  K  * ANL (I, 3) + (C + B_0)
         ANL (I, 4) =  K  * ANL (I, 4)
         ANL (I, 5) =  K2 * ANL (I, 5)
         ANL (I, 6) =  K2 * ANL (I, 6)
      END DO
C
      RETURN
      END

      SUBROUTINE POLYNOMIAL 
     &           (N, T, ANL_PAR, CRD, ANL)
C ---------------------------------------------------------------------------
C     Analytic solution : Polynomial
C
C              4
C             __
C             \                     K
C     PHI  =  /    Re (A (X + I * Z) )   +   A * T
C             --         K                    5
C             K = 1
C     
C     where Re means the real part 
C     
C     ANL (1:N, 1) : PHI
C     ANL (1:N, 2) : PHI_T
C     ANL (1:N, 3) : PHI_X
C     ANL (1:N, 4) : PHI_Z
C     ANL (1:N, 5) : PHI_XX
C     ANL (1:N, 6) : PHI_XZ
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      INTEGER(KIND=IK) N           ! Nr. of grid points (IN)
      REAL   (KIND=RK) T           ! Time (IN)
      REAL   (KIND=RK) ANL_PAR (*)
      REAL   (KIND=RK) CRD (4, *)  ! Grid point coordinates (IN)
      REAL   (KIND=RK) ANL (N, *)  ! Analytic solution (OUT)
C
      INTEGER(KIND=IK) I
      DOUBLE COMPLEX Z, DZ, DDZ, Z1, Z2, Z3, Z4, IZ
C
      IZ = COMPLEX (0.0D0, 1.0D0)
C
      DO I = 1, N
         Z1 = COMPLEX (CRD (1, I), CRD (2, I))
         Z2 = Z1 * Z1
         Z3 = Z1 * Z2
         Z4 = Z2 * Z2
         Z   = ANL_PAR (1) * Z1
     &       + ANL_PAR (2) * Z2
     &       + ANL_PAR (3) * Z3
     &       + ANL_PAR (4) * Z4
     &       + ANL_PAR (5) * T ** 3
         DZ  = ANL_PAR (1)
     &       + ANL_PAR (2) * Z1 * 2.0D0
     &       + ANL_PAR (3) * Z2 * 3.0D0
     &       + ANL_PAR (4) * Z3 * 4.0D0
         DDZ = ANL_PAR (2) *      2.0D0
     &       + ANL_PAR (3) * Z1 * 6.0D0
     &       + ANL_PAR (4) * Z2 * 1.2D1
         ANL (I, 1) = DBLE (Z)
         ANL (I, 2) = 3.0D0 * ANL_PAR (5) ** 2
         ANL (I, 3) = DBLE (DZ)
         ANL (I, 4) = DBLE (DZ * IZ)
         ANL (I, 5) = DBLE (DDZ)
         ANL (I, 6) = DBLE (DDZ * IZ)
      END DO
C
      RETURN
      END

      SUBROUTINE WAVEMAKER
     &           (N, T, ANL_PAR, CRD, ANL)
C ---------------------------------------------------------------------------
C     Analytic solution : WAVEMAKER
C
C     ANL (1:N, 1) : PHI
C     ANL (1:N, 2) : PHI_T
C     ANL (1:N, 3) : PHI_X
C     ANL (1:N, 4) : PHI_Z
C     ANL (1:N, 5) : PHI_XX
C     ANL (1:N, 6) : PHI_XZ
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'fle_params.inc'
C
      INTEGER(KIND=IK) N           ! Nr. of grid points (IN)
      REAL   (KIND=RK) T           ! Time (IN)
      REAL   (KIND=RK) ANL_PAR (17)
      REAL   (KIND=RK) CRD (4, *)  ! Grid point coordinates (IN)
      REAL   (KIND=RK) ANL (N, *)  ! Analytic solution (OUT)
C
      IF (ANL_PAR (17) .GT. 0.0) THEN
         CALL INTERP_WVM (N, WVM_IN, T, ANL_PAR, CRD (1, 1), CRD, ANL)
         WRITE (*, *) 'BEG->', ANL_PAR (1), T, ANL_PAR (2)
      ELSE
         CALL INTERP_WVM (N, WVM_IN, T, ANL_PAR, CRD (1, N), CRD, ANL)
         WRITE (*, *) 'END->', ANL_PAR (1), T, ANL_PAR (2)
      END IF
C
      RETURN
 101  FORMAT (1X,4F10.4)
      END
 
