      SUBROUTINE INTERP_WVM (NGP, U_NR, T, F, CRD_0, CRD, ANL)
C ---------------------------------------------------------------------------
C      
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C      
      INTEGER(KIND=IK) NGP
      INTEGER(KIND=IK) U_NR
      REAL   (KIND=RK) F (4, 4)
      REAL   (KIND=RK) CRD_0 (2)
      REAL   (KIND=RK) CRD   (4, NGP)
      REAL   (KIND=RK) ANL   (NGP, 6)
C
      INTEGER(KIND=IK) IGP
      REAL   (KIND=RK) T
      REAL   (KIND=RK) G (2, 4)
      REAL   (KIND=RK) DX, DZ, L
      REAL   (KIND=RK) U, V, OMG
C
      CALL INTERPOLATE (4, U_NR, T, F, G)
C      
      U   = G (1, 2)
      V   = G (1, 3)
      OMG = G (1, 4)
C      
      DO IGP = 1, NGP
         DX = CRD (1, IGP) - CRD_0 (1)
         DZ = CRD (2, IGP) - CRD_0 (2)
         L  = SQRT (DX ** 2 + DZ ** 2)
         IF (L .EQ. 0.0) L = 1.0
         ANL (IGP, 1) = 0.0D0             ! : PHI
         ANL (IGP, 2) = 0.0D0             ! : PHI_T
         ANL (IGP, 3) = U - OMG * DZ / L  ! : PHI_X
         ANL (IGP, 4) = V + OMG * DX / L  ! : PHI_Z
         ANL (IGP, 5) = 0.0D0             ! : PHI_XX
         ANL (IGP, 6) = 0.0D0             ! : PHI_XZ
      END DO
C      
      END
      
      SUBROUTINE INTERPOLATE (N, U_NR, TIME, F, G)
C ---------------------------------------------------------------------------
C     Interpolates F_I (T) and DF_I/DT (T)  at T = TIME, where F_I (T)
C     are third degree Hermite splines and returns them in G.
C      
C     G (1, I) =  F_I    (T = TIME)  I = 1 .. N
C     G (2, I) = DF_I/DT (T = TIME)  I = 1 .. N
C      
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C      
      INTEGER(KIND=IK) N
      INTEGER(KIND=IK) U_NR
      REAL   (KIND=RK) TIME
      REAL   (KIND=RK) F (4, N)
      REAL   (KIND=RK) G (2, N)
C     
      REAL   (KIND=RK) XI
      REAL   (KIND=RK) GET_XI
      EXTERNAL         GET_XI
C     find the interval containing TIME 
      DO 
         XI = GET_XI (TIME, F)
         IF (0.0 .LE. XI .AND. XI .LE. 1.0) EXIT
         CALL READ_NEXT_ELEMENTS (N, U_NR, F)
      END DO
C     and evaluate F at XI 
      CALL EVAL (N, XI, F, G)
C      
      END
  
      FUNCTION GET_XI (TIME, T)
C ---------------------------------------------------------------------------
C     Solves T (XI) - TIME = 0 by means of Newton iteration,
C     where T (XI) is a 3rd degree Hermite interpolation spline.
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      REAL   (KIND=RK) GET_XI
      REAL   (KIND=RK) TIME
      REAL   (KIND=RK) T (4)
C      
      REAL   (KIND=RK) XI
      REAL   (KIND=RK) TI (2)
C     initial guess for XI (is exact if T (XI) is linear)
      XI = (TIME - T (1)) / (T (2) - T (1))
      DO 
C        interpolate T (XI) and T_XI (XI)
         CALL EVAL (1, XI, T, TI)
C        check if we're close enough
         IF (ABS (TI (1) - TIME) .LT. 1e-12) EXIT
C        improve guess by using Newton
         XI = XI - (TI (1) - TIME) / TI (2)
      END DO
      GET_XI = XI 
      END
      
      SUBROUTINE INITIALIZE (N, U_NR, F)
C ---------------------------------------------------------------------------
C     
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C      
      INTEGER(KIND=IK) N
      INTEGER(KIND=IK) U_NR
      REAL   (KIND=RK) F (4, N)
C     read two rows of data to get the first spline element
      CALL READ_NEXT_ELEMENTS (N, U_NR, F)
      CALL READ_NEXT_ELEMENTS (N, U_NR, F)
C     
      END

      SUBROUTINE READ_NEXT_ELEMENTS (N, U_NR, F)
C ---------------------------------------------------------------------------
C     Reads the next spline elements from a unit and stores them is F.
C     The data in unit U_NR should look like this:
C     
C     T   T_XI   F_1   DF_1/DT   F_2   DF_2/DT   ......   F_N   DF_N/DT 
C     T   T_XI   F_1   DF_1/DT   F_2   DF_2/DT   ......   F_N   DF_N/DT 
C     T  ........
C      
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C      
      INTEGER(KIND=IK) N
      INTEGER(KIND=IK) U_NR
      REAL   (KIND=RK) F (4, N)
C      
      INTEGER(KIND=IK) I
C     Copy the old end node data to the begin node data
      DO I = 1, N
         F (1, I) =  F (2, I)
         F (3, I) =  F (4, I)
      END DO
C     and read the new end node data
      READ (U_NR, *) (F (2, I), F (4, I), I = 1, N)
C      
      END
      
      SUBROUTINE EVAL (N, XI, F, G)
C ---------------------------------------------------------------------------
C     
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      INTEGER(KIND=IK) N
      REAL   (KIND=RK) XI
      REAL   (KIND=RK) F (4, N)
      REAL   (KIND=RK) G (2, N)
C    
      INTEGER(KIND=IK) I, J
      REAL   (KIND=RK) W (4, 2)
C     precomute the weights for function and derivative evaluation
      CALL WEIGHT (XI, W)
C     compute the function and derivative evaluations
      DO J = 1, 2
         DO I = 1, N
            G (J, I) = W (1, J) * F (1, I) + W (2, J) * F (2, I)
     &               + W (3, J) * F (3, I) + W (4, J) * F (4, I)
         END DO
      END DO
      DO I = 2, N
         G (2, I) = G (2, I) / G (2, 1)
      END DO
C      
      END
