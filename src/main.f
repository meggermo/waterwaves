      SUBROUTINE MAIN
     &           (NSD, NNW, NGP, NAE, NNW_SD, NGP_SD, NGP_NW)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'tme_params.inc'
      INCLUDE 'tme_funcs.inc'
C
      INTEGER(KIND=IK) NSD                      ! Nr. of subdomains      (IN)
      INTEGER(KIND=IK) NNW                      ! Nr. of networks        (IN)
      INTEGER(KIND=IK) NGP                      ! Nr. of grid points     (IN)
      INTEGER(KIND=IK) NAE                      ! Nr. of array elements  (IN)
      INTEGER(KIND=IK) NNW_SD (*)               ! Nr. of netw.  per subd (IN)
      INTEGER(KIND=IK) NGP_SD (*)               ! Nr. of gridp. per subd (IN)
      INTEGER(KIND=IK) NGP_NW (*)               ! Nr. of gridp. per netw (IN)
C
      INTEGER(KIND=IK) NW_IPAR (N_IP, NNW)      ! Integer network parameters
      REAL   (KIND=RK) NW_RPAR (N_RP, NNW)      ! Real    network parameters
C
      REAL   (KIND=RK) CRD (4 * NGP * 4)        ! Node coordinates
      REAL   (KIND=RK) PHI (2 * NGP * 4)        ! Potential
      REAL   (KIND=RK) PHN (2 * NGP * 4)        ! Normal deriv. of potential
      REAL   (KIND=RK) S_0 (NAE)                ! source coef. matrix
      REAL   (KIND=RK) S_1 (NAE)
      REAL   (KIND=RK) D_0 (NAE)                ! Dipole coef. matrix
      REAL   (KIND=RK) D_1 (NAE)
C
      LOGICAL(KIND=LK) THE_END
      DATA             THE_END /.FALSE./
C
      WRITE (USR_O, *) '========== PROBLEM DIMENSIONS ====='
      WRITE (USR_O, *) 'NR. OF SUBDOMAINS: ', NSD
      WRITE (USR_O, *) 'NR. OF NETWORKS:   ', NNW
      WRITE (USR_O, *) 'NR. OF GRID POINTS:', NGP
      WRITE (USR_O, *) 'NR. OF ARRAY ELEMS:', NAE
C     Copy the nr. of grid points into the nw_ipar array
      CALL COPY_NGP_TO_NW_IPAR
     &     (NNW, NGP_NW, NW_IPAR)
C     Read all the input parameters
      CALL READ_PARAMETERS
     &           (NSD, NNW_SD, NW_IPAR, NW_RPAR)
C     generate the initial grid from the given input
      CALL GENERATE_GRID
     &     (NNW, NW_IPAR, NW_RPAR, CRD)
C
      DO WHILE (.NOT. THE_END)
C        integrate the ODE's with Runge-Kutta-4 until we're done
         CALL RK4
     &        (NSD,     NGP,      THE_END,
     &         NNW_SD,  NGP_SD,
     &         NW_IPAR, NW_RPAR,
     &         CRD,     PHI,      PHN,
     &         S_0,     S_1,      D_0,      D_1)
C        Increment I_PLOT
         I_PLOT = MOD (I_PLOT + 1, I_PLOT_MOD)
      END DO
C
      RETURN
      END
