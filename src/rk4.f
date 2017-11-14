      SUBROUTINE RK4
     &           (NSD,     NGP,      THE_END,
     &            NNW_SD,  NGP_SD,
     &            NW_IPAR, NW_RPAR,
     &            CRD,     PHI,      PHN,
     &            S_0,     S_1,      D_0,      D_1)
C ---------------------------------------------------------------------------
C     Performs a RK-4 time step by integrating the ODE's
C     from T to T + DT
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'rle_params.inc'
      INCLUDE 'tme_funcs.inc'
C
      INTEGER(KIND=IK) NSD                ! nr. of subdomains (IN)
      INTEGER(KIND=IK) NGP                ! nr. of grid points (IN)
      LOGICAL(KIND=LK) THE_END            ! set to true if  T > T_END
      INTEGER(KIND=IK) NNW_SD (*)         ! nr. of netw. per subd. (IN)
      INTEGER(KIND=IK) NGP_SD (*)         ! nr. of grid. p. per subd. (IN)
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)  ! integer netw. params. (IN)
      REAL   (KIND=RK) NW_RPAR (N_RP, *)  ! real netw. params. (IN)
      REAL   (KIND=RK) CRD (4 * NGP, *)   ! grid points (IN/OUT)
      REAL   (KIND=RK) PHI (2 * NGP, *)   ! potential (IN/OUT)
      REAL   (KIND=RK) PHN (2 * NGP, *)   ! it's normal derivative (IN/OUT)
      REAL   (KIND=RK) S_0 (*)            ! source coef. matrix
      REAL   (KIND=RK) S_1 (*)            ! idem.
      REAL   (KIND=RK) D_0 (*)            ! dipole coef. matrix
      REAL   (KIND=RK) D_1 (*)            ! idem.
C
      REAL   (KIND=RK) T                  ! time at stage IST)
      REAL   (KIND=RK) DT                 ! time step
      REAL   (KIND=RK) TE                 ! time step
      INTEGER(KIND=IK) I, IX, IZ, IP, IS
C
      T  = GET_TIME    ()
      TE = GET_TEND    ()
      DT = GET_DELTA_T ()
      WRITE(USR_O,*) T, TE, DT
C
      IF (DT .NE. 0.0D0) THEN
         WRITE (USR_O, 1001) 'TIME = ', T / TE, ' * T_END'
      END IF
C     STAGE 1
      CALL DIFF_EQNS
     &     (1,          NSD,
     &      NNW_SD,     NGP_SD,     NW_IPAR,    NW_RPAR,
     &      CRD (1, 1), PHI (1, 1), PHN (1, 1),
     &      CRD (1, 2), PHI (1, 2), PHN (1, 2),
     &      S_0,        S_1,        D_0,        D_1)
C     Check if we need to do more timestepping
      IF (DT .NE. 0.0D0 .AND. T .LT. TE) THEN
C     Add the RHS's to the time dependent ODE's
      DO I = 0, NGP - 1
         IX = 1 + 4 * I
         IZ = 2 + 4 * I
         IP = 1 + 2 * I
         IS = 2 + 2 * I
         CRD (IX, 1) = CRD (IX, 1) + 0.5D0 * DT * CRD (IX, 2)
         CRD (IZ, 1) = CRD (IZ, 1) + 0.5D0 * DT * CRD (IZ, 2)
         PHI (IP, 1) = PHI (IP, 1) + 0.5D0 * DT * PHI (IP, 2)
         PHI (IS, 1) = PHI (IS, 1) + 0.5D0 * DT * PHI (IS, 2)
      END DO
C     Compute the spline derivatives for the integrated variables
      CALL DERIVS (NSD, NNW_SD, NGP_SD, NW_IPAR, CRD, PHI, PHN)
C     Set the physical time
      CALL SET_TIME (T + 0.5 * DT)
C     STAGE 2
      CALL DIFF_EQNS
     &     (2,          NSD,
     &      NNW_SD,     NGP_SD,     NW_IPAR,    NW_RPAR,
     &      CRD (1, 1), PHI (1, 1), PHN (1, 1),
     &      CRD (1, 3), PHI (1, 3), PHN (1, 3),
     &      S_0,        S_1,        D_0,        D_1)
C     Add the RHS's to the time dependent ODE's
      DO I = 0, NGP - 1
         IX = 1 + 4 * I
         IZ = 2 + 4 * I
         IP = 1 + 2 * I
         IS = 2 + 2 * I
         CRD (IX, 1)= CRD (IX,1) + 0.5D0*DT * (CRD (IX,3) - CRD (IX,2))
         CRD (IZ, 1)= CRD (IZ,1) + 0.5D0*DT * (CRD (IZ,3) - CRD (IZ,2))
         PHI (IP, 1)= PHI (IP,1) + 0.5D0*DT * (PHI (IP,3) - PHI (IP,2))
         PHI (IS, 1)= PHI (IS,1) + 0.5D0*DT * (PHI (IS,3) - PHI (IS,2))
      END DO
C     Compute the spline derivatives for the integrated variables
      CALL DERIVS (NSD, NNW_SD, NGP_SD, NW_IPAR, CRD, PHI, PHN)
C     Set the physical time
      CALL SET_TIME (T + 0.5 * DT)
C     STAGE 3
      CALL DIFF_EQNS
     &     (3,          NSD,
     &      NNW_SD,     NGP_SD,     NW_IPAR,    NW_RPAR,
     &      CRD (1, 1), PHI (1, 1), PHN (1, 1),
     &      CRD (1, 4), PHI (1, 4), PHN (1, 4),
     &      S_0,        S_1,        D_0,        D_1)
C     Add the RHS's to the time dependent ODE's
      DO I = 0, NGP - 1
         IX = 1 + 4 * I
         IZ = 2 + 4 * I
         IP = 1 + 2 * I
         IS = 2 + 2 * I
         CRD (IX, 1) = CRD (IX,1) + DT* (CRD (IX,4) - 0.5D0* CRD (IX,3))
         CRD (IZ, 1) = CRD (IZ,1) + DT* (CRD (IZ,4) - 0.5D0* CRD (IZ,3))
         PHI (IP, 1) = PHI (IP,1) + DT* (PHI (IP,4) - 0.5d0* PHI (IP,3))
         PHI (IS, 1) = PHI (IS,1) + DT* (PHI (IS,4) - 0.5d0* PHI (IS,3))
         CRD (IX, 2) = CRD (IX,2) + 2.0D0 * CRD (IX,3)
         CRD (IZ, 2) = CRD (IZ,2) + 2.0D0 * CRD (IZ,3)
         PHI (IP, 2) = PHI (IP,2) + 2.0D0 * PHI (IP,3)
         PHI (IS, 2) = PHI (IS,2) + 2.0D0 * PHI (IS,3)
      END DO
C     Compute the spline derivatives for the integrated variables
      CALL DERIVS (NSD, NNW_SD, NGP_SD, NW_IPAR, CRD, PHI, PHN)
C     Set the physical time
      CALL SET_TIME (T + DT)
C     STAGE 4
      CALL DIFF_EQNS
     &     (4,          NSD,
     &      NNW_SD,     NGP_SD,     NW_IPAR,    NW_RPAR,
     &      CRD (1, 1), PHI (1, 1), PHN (1, 1),
     &      CRD (1, 3), PHI (1, 3), PHN (1, 3),
     &      S_0,        S_1,        D_0,        D_1)
C     Add the RHS's to the time dependent ODE's
      DO I = 0, NGP - 1
         IX = 1 + 4 * I
         IZ = 2 + 4 * I
         IP = 1 + 2 * I
         IS = 2 + 2 * I
         CRD (IX, 1) = CRD (IX, 1)
     &         + DT * (CRD (IX, 2)
     &              +  CRD (IX, 3) - 4.0D0 * CRD (IX, 4)) / 6.0D0
         CRD (IZ, 1) = CRD (IZ, 1)
     &         + DT * (CRD (IZ, 2)
     &              +  CRD (IZ, 3) - 4.0D0 * CRD (IZ, 4)) / 6.0D0
         PHI (IP, 1) = PHI (IP, 1)
     &         + DT * (PHI (IP, 2)
     &              +  PHI (IP, 3) - 4.0D0 * PHI (IP, 4)) / 6.0D0
         PHI (IS, 1) = PHI (IS, 1)
     &         + DT * (PHI (IS, 2)
     &              +  PHI (IS, 3) - 4.0D0 * PHI (IS, 4)) / 6.0D0
      END DO
C     Compute the spline derivatives for the integrated variables
      CALL DERIVS (NSD, NNW_SD, NGP_SD, NW_IPAR, CRD, PHI, PHN)
      ELSE
      THE_END = .TRUE.
      END IF
C
      RETURN
 1001 FORMAT (1X,A,F7.3,A)
      END
      SUBROUTINE DERIVS
     &           (NSD, NNW_SD, NGP_SD, NW_IPAR, CRD, PHI, PHN)
C ---------------------------------------------------------------------------
C
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      INTEGER(KIND=IK) NSD                ! nr. of subdomains (IN)
      INTEGER(KIND=IK) NNW_SD (*)         ! nr. of netw. per subd. (IN)
      INTEGER(KIND=IK) NGP_SD (*)         ! nr. of grid. p. per subd. (IN)
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)  ! integer netw. params. (IN)
      REAL   (KIND=RK) CRD (4, *)         ! grid points (IN/OUT)
      REAL   (KIND=RK) PHI (2, *)         ! potential (IN/OUT)
      REAL   (KIND=RK) PHN (2, *)         ! it's normal derivative (IN/OUT)
C
      INTEGER(KIND=IK) ISD, INW, NNW, IGP, NGP
C
      INW = 1
      IGP = 1
C
      DO ISD = 1, NSD
         NNW = NNW_SD (ISD)
         NGP = NGP_SD (ISD)
         CALL DERIVS_SD
     &        (NNW,          NW_IPAR (1, INW),
     &         CRD (1, IGP), PHI (1, IGP),     PHN (1, IGP))
         INW = INW + NNW
         IGP = IGP + NGP
      END DO
C
      RETURN
      END
      SUBROUTINE DERIVS_SD
     &           (NNW, NW_IPAR, CRD, PHI, PHN)
C ---------------------------------------------------------------------------
C
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'tme_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'net_funcs.inc'
C
      INTEGER(KIND=IK) NNW                ! nr. of networks (IN)
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)  ! integer netw. params. (IN)
      REAL   (KIND=RK) CRD (4, *)         ! grid points (IN/OUT)
      REAL   (KIND=RK) PHI (2, *)         ! potential (IN/OUT)
      REAL   (KIND=RK) PHN (2, *)         ! it's normal derivative (IN/OUT)
C
      INTEGER(KIND=IK) INW, IGP, NGP, BCT, JGP
      INTEGER(KIND=IK) BC (2)
      DATA             BC /0, 0/
C
      IGP = 1
C
      DO INW = 1, NNW
         NGP = GET_NGP (NW_IPAR (1, INW))
         BCT = GET_BCT (NW_IPAR (1, INW))
         IF (GRID_TIME_DEPENDENT) THEN
            CALL Spline (BC, NGP, 2, 4, CRD (1, IGP), CRD (3, IGP))
         END IF
         IF (BCT .LT. 0) THEN
            CALL Spline (BC, NGP, 1, 2, PHI (1, IGP), PHI (2, IGP))
            DO JGP = IGP, IGP + NGP - 1, NGP - 1
               WRITE(*,'(1X,I4,2E14.6)') JGP, CRD(1,JGP), PHI(2,JGP)
            END DO
         END  IF
         IGP = IGP + NGP
      END DO
C
      RETURN
      END