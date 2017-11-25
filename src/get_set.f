      FUNCTION GET_NGP
     &         (NW_IPAR)
C ---------------------------------------------------------------------------
C    Returns the number of grid points for the network
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      INTEGER(KIND=IK) GET_NGP
      INTEGER(KIND=IK) NW_IPAR (N_IP)
C
      GET_NGP = NW_IPAR (I_NGP)
C
      RETURN
      END FUNCTION GET_NGP

      SUBROUTINE SET_NGP
     &           (NGP, NW_IPAR)
C ---------------------------------------------------------------------------
C    Sets the number of grid points to NGP for the network
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'max_params.inc'
C
      INTEGER(KIND=IK) NGP
      INTEGER(KIND=IK) NW_IPAR (N_IP)
C
      CALL I_ASSERT (NGP_MIN, NGP, NGP_MAX, 'NGP', 'SET_NGP')
      NW_IPAR (I_NGP) = NGP
C
      RETURN
      END SUBROUTINE SET_NGP

      FUNCTION GET_BCT
     &         (NW_IPAR)
C ---------------------------------------------------------------------------
C     Returns the boundary condition type of the network
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      INTEGER(KIND=IK) GET_BCT
      INTEGER(KIND=IK) NW_IPAR (N_IP)
C
      GET_BCT = NW_IPAR (I_BCT)
C
      RETURN
      END FUNCTION GET_BCT

      SUBROUTINE SET_BCT
     &           (BCT, NW_IPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'bct_params.inc'
C
      INTEGER(KIND=IK) BCT
      INTEGER(KIND=IK) NW_IPAR (N_IP)
C
      CALL I_ASSERT (BCT_FIRST, ABS (BCT), BCT_LAST, 'BCT', 'SET_BCT')
      NW_IPAR (I_BCT) = BCT
C
      RETURN
      END SUBROUTINE SET_BCT

      FUNCTION GET_PLT
     &         (NW_IPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      INTEGER(KIND=IK) GET_PLT
      INTEGER(KIND=IK) NW_IPAR (N_IP)
C
      GET_PLT = NW_IPAR (I_PLT)
C
      RETURN
      END FUNCTION GET_PLT

      SUBROUTINE SET_PLT
     &           (PLT, NW_IPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      INTEGER(KIND=IK) PLT
      INTEGER(KIND=IK) NW_IPAR (N_IP)
C
      NW_IPAR (I_PLT) = PLT
C
      RETURN
      END SUBROUTINE SET_PLT

      FUNCTION GET_GRT
     &         (NW_IPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      INTEGER(KIND=IK) GET_GRT
      INTEGER(KIND=IK) NW_IPAR (N_IP)
C
      GET_GRT = NW_IPAR (I_GRT)
C
      RETURN
      END FUNCTION GET_GRT

      SUBROUTINE SET_GRT
     &           (GRT, NW_IPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'grd_params.inc'
C
      INTEGER(KIND=IK) GRT
      INTEGER(KIND=IK) NW_IPAR (N_IP)
C
      CALL I_ASSERT (GRD_FIRST, GRT, GRD_LAST, 'GRD', 'SET_GRT')
      NW_IPAR (I_GRT) = GRT
C
      RETURN
      END SUBROUTINE SET_GRT

      FUNCTION GET_GMT
     &         (NW_IPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      INTEGER(KIND=IK) GET_GMT
      INTEGER(KIND=IK) NW_IPAR (N_IP)
C
      GET_GMT = NW_IPAR (I_GMT)
C
      RETURN
      END FUNCTION GET_GMT

      SUBROUTINE SET_GMT
     &           (GMT, NW_IPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      INTEGER(KIND=IK) GMT
      INTEGER(KIND=IK) NW_IPAR (N_IP)
C
      NW_IPAR (I_GMT) = GMT
C
      RETURN
      END SUBROUTINE SET_GMT

      FUNCTION GET_GAT
     &         (NW_IPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      INTEGER(KIND=IK) GET_GAT
      INTEGER(KIND=IK) NW_IPAR (N_IP)
C
      GET_GAT = NW_IPAR (I_GAT)
C
      RETURN
      END FUNCTION GET_GAT

      SUBROUTINE SET_GAT
     &           (GAT, NW_IPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      INTEGER(KIND=IK) GAT
      INTEGER(KIND=IK) NW_IPAR (N_IP)
C
      NW_IPAR (I_GAT) = GAT
C
      RETURN
      END SUBROUTINE SET_GAT

      FUNCTION GET_SPT
     &         (ISP, NW_IPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      INTEGER(KIND=IK) GET_SPT
      INTEGER(KIND=IK) ISP
      INTEGER(KIND=IK) NW_IPAR (N_IP)
C
      CALL I_ASSERT (1, ISP, N_SPT, 'ISP', 'GET_SPT')
      GET_SPT = NW_IPAR (I_SPT + ISP - 1)
C
      RETURN
      END FUNCTION GET_SPT

      SUBROUTINE SET_SPT
     &           (ISP, SPT, NW_IPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      INTEGER(KIND=IK) ISP
      INTEGER(KIND=IK) SPT
      INTEGER(KIND=IK) NW_IPAR (N_IP)
C
      CALL I_ASSERT ( 1, ISP, N_SPT, 'ISP', 'GET_SPT')
      CALL I_ASSERT (-1, SPT, 4,     'SPT', 'SET_SPT')
      NW_IPAR (I_SPT + ISP - 1) = SPT
C
      RETURN
      END SUBROUTINE SET_SPT

      FUNCTION GET_ITF
     &         (NW_IPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      INTEGER(KIND=IK) GET_ITF
      INTEGER(KIND=IK) NW_IPAR (N_IP)
C
      GET_ITF = NW_IPAR (I_ITF)
C
      RETURN
      END FUNCTION GET_ITF

      SUBROUTINE SET_ITF
     &           (ITF, NW_IPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'max_params.inc'
C
      INTEGER(KIND=IK) ITF
      INTEGER(KIND=IK) NW_IPAR (N_IP)
C
      CALL I_ASSERT (0, ITF, NSD_MAX - 1, 'ITF', 'SET_ITF')
      NW_IPAR (I_ITF) = ITF
C
      RETURN
      END SUBROUTINE SET_ITF

      FUNCTION GET_AED
     &         (IED, NW_IPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      INTEGER(KIND=IK) GET_AED
      INTEGER(KIND=IK) IED
      INTEGER(KIND=IK) NW_IPAR (N_IP)
C
      CALL I_ASSERT (1, IED, 2, 'IED', 'GET_AED')
C
      GET_AED = NW_IPAR (I_AED + IED - 1)
C
      RETURN
      END FUNCTION GET_AED

      SUBROUTINE SET_AED
     &           (IED, AED, NW_IPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      INTEGER(KIND=IK) IED
      INTEGER(KIND=IK) AED
      INTEGER(KIND=IK) NW_IPAR (N_IP)
C
      CALL I_ASSERT (1, IED, 2, 'IED', 'SET_AED')
      CALL I_ASSERT (1, AED, 2, 'IED', 'SET_AED')
C
      NW_IPAR (I_AED + IED - 1) = AED
C
      RETURN
      END SUBROUTINE SET_AED

      FUNCTION GET_ANW
     &         (IED, NW_IPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      INTEGER(KIND=IK) GET_ANW
      INTEGER(KIND=IK) IED
      INTEGER(KIND=IK) NW_IPAR (N_IP)
C
      CALL I_ASSERT (1, IED, 2, 'IED', 'GET_AED')
C
      GET_ANW = NW_IPAR (I_ANW + IED - 1)
C
      RETURN
      END FUNCTION GET_ANW

      SUBROUTINE SET_ANW
     &           (IED, ANW, NW_IPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      INTEGER(KIND=IK) IED
      INTEGER(KIND=IK) ANW
      INTEGER(KIND=IK) NW_IPAR (N_IP)
C
      CALL I_ASSERT (1, IED, 2, 'IED', 'SET_AED')
C
      NW_IPAR (I_ANW + IED - 1) = ANW
C
      RETURN
      END SUBROUTINE SET_ANW

      FUNCTION GET_AGP
     &         (IED, INW, NW_IPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      INTEGER(KIND=IK) GET_AGP
      INTEGER(KIND=IK) IED
      INTEGER(KIND=IK) INW
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)
C
      INTEGER(KIND=IK) ANW, AED, AGP, I
      INTEGER(KIND=IK) GET_ANW, GET_AED, GET_NGP
      EXTERNAL         GET_ANW, GET_AED, GET_NGP
C
      ANW = GET_ANW (IED, NW_IPAR (1, INW))
      AED = GET_AED (IED, NW_IPAR (1, INW))
      AGP = 1
      DO I = 1, ANW - 1
         AGP = AGP + GET_NGP (NW_IPAR (1, I))
      END DO
      GET_AGP = AGP + (AED - 1) * (GET_NGP (NW_IPAR (1, ANW)) - 1)
      RETURN
      END FUNCTION GET_AGP

      SUBROUTINE GET_BCP
     &           (NW_RPAR, PAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      REAL(KIND=RK) NW_RPAR (N_RP)
      REAL(KIND=RK) PAR (*)
C
      INTEGER(KIND=IK) IP
C
      DO IP = 1, N_BCP
         PAR (IP) = NW_RPAR (I_BCP + IP - 1)
      END DO
C
      RETURN
      END SUBROUTINE GET_BCP

      SUBROUTINE SET_BCP
     &           (PAR, NW_RPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      REAL(KIND=RK) PAR (*)
      REAL(KIND=RK) NW_RPAR (N_RP)
C
      INTEGER(KIND=IK) IP
C
      DO IP = 1, N_BCP
         NW_RPAR (I_BCP + IP - 1) = PAR (IP)
      END DO
C
      RETURN
      END SUBROUTINE SET_BCP

      SUBROUTINE GET_GRP
     &           (NW_RPAR, PAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      REAL(KIND=RK) NW_RPAR (N_RP)
      REAL(KIND=RK) PAR (*)
C
      INTEGER(KIND=IK) IP
C
      DO IP = 1, N_GRP
         PAR (IP) = NW_RPAR (I_GRP + IP - 1)
      END DO
C
      RETURN
      END SUBROUTINE GET_GRP

      SUBROUTINE SET_GRP
     &           (PAR, NW_RPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      REAL(KIND=RK) PAR (*)
      REAL(KIND=RK) NW_RPAR (N_RP)
C
      INTEGER(KIND=IK) IP
C
      DO IP = 1, N_GRP
         NW_RPAR (I_GRP + IP - 1) = PAR (IP)
      END DO
C
      RETURN
      END SUBROUTINE SET_GRP

      SUBROUTINE GET_GMP
     &           (NW_RPAR, PAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      REAL(KIND=RK) NW_RPAR (N_RP)
      REAL(KIND=RK) PAR (*)
C
      INTEGER(KIND=IK) IP
C
      DO IP = 1, N_GMP
         PAR (IP) = NW_RPAR (I_GMP + IP - 1)
      END DO
C
      RETURN
      END SUBROUTINE GET_GMP

      SUBROUTINE SET_GMP
     &           (PAR, NW_RPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      REAL(KIND=RK) PAR (*)
      REAL(KIND=RK) NW_RPAR (N_RP)
C
      INTEGER(KIND=IK) IP
C
      DO IP = 1, N_GMP
         NW_RPAR (I_GMP + IP - 1) = PAR (IP)
      END DO
C
      RETURN
      END SUBROUTINE SET_GMP

      SUBROUTINE GET_GAP
     &           (NW_RPAR, PAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      REAL(KIND=RK) NW_RPAR (N_RP)
      REAL(KIND=RK) PAR (*)
C
      INTEGER(KIND=IK) IP
C
      DO IP = 1, N_GAP
         PAR (IP) = NW_RPAR (I_GAP + IP - 1)
      END DO
C
      RETURN
      END SUBROUTINE GET_GAP

      SUBROUTINE SET_GAP
     &           (PAR, NW_RPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      REAL(KIND=RK) PAR (*)
      REAL(KIND=RK) NW_RPAR (N_RP)
C
      INTEGER(KIND=IK) IP
C
      DO IP = 1, N_GAP
         NW_RPAR (I_GAP + IP - 1) = PAR (IP)
      END DO
C
      RETURN
      END SUBROUTINE SET_GAP

      FUNCTION GET_TIME ()
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'tme_params.inc'
C
      REAL(KIND=RK) GET_TIME
C
      GET_TIME = TM_RPAR (IT_CUR)
C
      RETURN
      END FUNCTION GET_TIME

      SUBROUTINE SET_TIME (T)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'tme_params.inc'
C
      REAL(KIND=RK) T
C
      TM_RPAR (IT_CUR) = T
C
      RETURN
      END SUBROUTINE SET_TIME

      FUNCTION GET_TBEG ()
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'tme_params.inc'
C
      REAL(KIND=RK) GET_TBEG
C
      GET_TBEG = TM_RPAR (IT_BEG)
C
      RETURN
      END FUNCTION GET_TBEG

      SUBROUTINE SET_TBEG (T)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'tme_params.inc'
C
      REAL(KIND=RK) T
C
      TM_RPAR (IT_BEG) = T
C
      RETURN
      END SUBROUTINE SET_TBEG

      FUNCTION GET_TEND ()
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'tme_params.inc'
C
      REAL(KIND=RK) GET_TEND
C
      GET_TEND = TM_RPAR (IT_END)
C
      RETURN
      END FUNCTION GET_TEND

      SUBROUTINE SET_TEND (T)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'tme_params.inc'
C
      REAL(KIND=RK) T
C
      TM_RPAR (IT_END) = T
C
      RETURN
      END SUBROUTINE SET_TEND

      FUNCTION GET_DELTA_T ()
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'tme_params.inc'
C
      REAL(KIND=RK) GET_DELTA_T
C
      GET_DELTA_T = TM_RPAR (IT_STP)
C
      RETURN
      END FUNCTION GET_DELTA_T

      SUBROUTINE SET_DELTA_T (DT)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'tme_params.inc'
C
      REAL(KIND=RK) DT
C
      TM_RPAR (IT_STP) = DT
C
      RETURN
      END SUBROUTINE SET_DELTA_T

      FUNCTION GET_DELTA_T_MIN ()
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'tme_params.inc'
C
      REAL(KIND=RK) GET_DELTA_T_MIN
C
      GET_DELTA_T_MIN = TM_RPAR (IT_STPMIN)
C
      RETURN
      END FUNCTION GET_DELTA_T_MIN

      SUBROUTINE SET_DELTA_T_MIN (DT)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'tme_params.inc'
C
      REAL(KIND=RK) DT
C
      TM_RPAR (IT_STPMIN) = DT
C
      RETURN
      END SUBROUTINE SET_DELTA_T_MIN

      FUNCTION GET_DELTA_T_MAX ()
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'tme_params.inc'
C
      REAL(KIND=RK) GET_DELTA_T_MAX
C
      GET_DELTA_T_MAX = TM_RPAR (IT_STPMAX)
C
      RETURN
      END FUNCTION GET_DELTA_T_MAX

      SUBROUTINE SET_DELTA_T_MAX (DT)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'tme_params.inc'
C
      REAL(KIND=RK) DT
C
      TM_RPAR (IT_STPMAX) = DT
C
      RETURN
      END SUBROUTINE SET_DELTA_T_MAX

      FUNCTION GET_SUBD_MAX ()
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'tme_params.inc'
C
      INTEGER(KIND=IK) GET_SUBD_MAX
C
      GET_SUBD_MAX = I_SUBD_MAX
C
      RETURN
      END FUNCTION GET_SUBD_MAX

      SUBROUTINE SET_SUBD_MAX (SUBD_MAX)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'tme_params.inc'
C
      INTEGER(KIND=IK) SUBD_MAX
C
      I_SUBD_MAX = SUBD_MAX
C
      RETURN
      END SUBROUTINE SET_SUBD_MAX


      FUNCTION GET_PERIOD ()
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'anl_params.inc'
      INCLUDE 'rle_params.inc'
C
      REAL(KIND=RK) GET_PERIOD
      REAL(KIND=RK) L, K, H
C
      IF (ANL_TYPE .EQ. ANL_LINWAVE) THEN
C        LINEAR MONOCHROMATIC WAVE PERIOD
         L  = ANL_PAR (2)
         H  = ANL_PAR (3)
         K  = TWO_PI / L
         GET_PERIOD = TWO_PI / SQRT (Grav * K * TANH (K * H))
      ELSE IF (ANL_TYPE .EQ. ANL_RFWAVE) THEN
C        REINECKER & FENTON PERIOD
         GET_PERIOD = TWO_PI / ANL_PAR (5)
      ELSE
C        ANY OTHER TYPE
         GET_PERIOD = 1.0D0
      END IF
C
      RETURN
      END FUNCTION GET_PERIOD

      FUNCTION GET_INT_ABS ()
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'slv_params.inc'
C
      REAL(KIND=RK) GET_INT_ABS
C
      GET_INT_ABS = SO_RPAR (IS_INT_ABS)
C
      RETURN
      END FUNCTION GET_INT_ABS

      SUBROUTINE SET_INT_ABS (INT_ABS)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'slv_params.inc'
C
      REAL(KIND=RK) INT_ABS
C
      CALL R_ASSERT (0.0D0, INT_ABS, 1.0D0, 'INT_ABS', 'SET_INT_ABS')
      SO_RPAR (IS_INT_ABS) = INT_ABS
C
      RETURN
      END SUBROUTINE SET_INT_ABS

      FUNCTION GET_INT_REL ()
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'slv_params.inc'
C
      REAL(KIND=RK) GET_INT_REL
C
      GET_INT_REL = SO_RPAR (IS_INT_REL)
C
      RETURN
      END FUNCTION GET_INT_REL

      SUBROUTINE SET_INT_REL (INT_REL)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'slv_params.inc'
C
      REAL(KIND=RK) INT_REL
C
      CALL R_ASSERT (0.0D0, INT_REL, 1.0D0, 'INT_REL', 'SET_INT_REL')
      SO_RPAR (IS_INT_REL) = INT_REL
C
      RETURN
      END SUBROUTINE SET_INT_REL

      FUNCTION GET_ITR_TOL (I)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'slv_params.inc'
C
      INTEGER(KIND=IK) I
      REAL(KIND=RK) GET_ITR_TOL
C
      CALL I_ASSERT (1, I, 2, 'I', 'GET_ITR_TOL')
      GET_ITR_TOL = SO_RPAR (IS_ITR_TOL + I - 1)
C
      RETURN
      END FUNCTION GET_ITR_TOL

      SUBROUTINE SET_ITR_TOL (I, TOL)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'slv_params.inc'
C
      INTEGER(KIND=IK) I
      REAL(KIND=RK) TOL
C
      CALL I_ASSERT (1,     I,   2,     'I',   'SET_ITR_TOL')
      CALL R_ASSERT (0.0D0, TOL, 1.0D0, 'TOL', 'SET_ITR_TOL')
C
      SO_RPAR (IS_ITR_TOL + I - 1) = TOL
C
      RETURN
      END SUBROUTINE SET_ITR_TOL

      FUNCTION GET_ITR_OMG (I)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'slv_params.inc'
C
      INTEGER(KIND=IK) I
      REAL(KIND=RK) GET_ITR_OMG
C
      CALL I_ASSERT (1, I, 2, 'I', 'GET_ITR_TOL')
      GET_ITR_OMG = SO_RPAR (IS_ITR_OMG + I - 1)
C
      RETURN
      END FUNCTION GET_ITR_OMG

      SUBROUTINE SET_ITR_OMG (I, OMG)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'slv_params.inc'
C
      INTEGER(KIND=IK) I
      REAL(KIND=RK) OMG
C
      CALL I_ASSERT (1,     I,   2,     'I',   'SET_ITR_OMG')
      CALL R_ASSERT (0.0D0, OMG, 1.0D0, 'OMG', 'SET_ITR_OMG')
C
      SO_RPAR (IS_ITR_OMG + I - 1) = OMG
C
      RETURN
      END SUBROUTINE SET_ITR_OMG

      FUNCTION GET_ITR_MAXIT ()
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'slv_params.inc'
C
      INTEGER(KIND=IK) GET_ITR_MAXIT
C
      GET_ITR_MAXIT = SO_IPAR (IS_ITR_MAXIT)
C
      RETURN
      END FUNCTION GET_ITR_MAXIT

      SUBROUTINE SET_ITR_MAXIT (MAXIT)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'slv_params.inc'
      INCLUDE 'max_params.inc'
C
      INTEGER(KIND=IK) MAXIT
C
      CALL I_ASSERT (0, MAXIT, ITR_MAX, 'MAXIT', 'SET_ITR_MAXIT')
      SO_IPAR (IS_ITR_MAXIT) = MAXIT
C
      RETURN
      END SUBROUTINE SET_ITR_MAXIT

      FUNCTION GET_DOD_TOL (I)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'slv_params.inc'
C
      INTEGER(KIND=IK) I
      REAL(KIND=RK) GET_DOD_TOL
C
      CALL I_ASSERT (1, I, 2, 'I', 'GET_DOD_TOL')
      GET_DOD_TOL = SO_RPAR (IS_DOD_TOL + I - 1)
C
      RETURN
      END FUNCTION GET_DOD_TOL

      SUBROUTINE SET_DOD_TOL (I, TOL)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'slv_params.inc'
C
      INTEGER(KIND=IK) I
      REAL(KIND=RK) TOL
C
      CALL I_ASSERT (1,     I,   2,     'I',   'SET_DOD_TOL')
      CALL R_ASSERT (0.0D0, TOL, 1.0D0, 'TOL', 'SET_DOD_TOL')
C
      SO_RPAR (IS_DOD_TOL + I - 1) = TOL
C
      RETURN
      END SUBROUTINE SET_DOD_TOL

      FUNCTION GET_DOD_OMG (I)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'slv_params.inc'
C
      INTEGER(KIND=IK) I
      REAL(KIND=RK) GET_DOD_OMG
C
      CALL I_ASSERT (1, I, 2, 'I', 'GET_DOD_TOL')
      GET_DOD_OMG = SO_RPAR (IS_DOD_OMG + I - 1)
C
      RETURN
      END FUNCTION GET_DOD_OMG

      SUBROUTINE SET_DOD_OMG (I, OMG)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'slv_params.inc'
C
      INTEGER(KIND=IK) I
      REAL(KIND=RK) OMG
C
      CALL I_ASSERT (1,     I,   2,     'I',   'SET_DOD_TOL')
      CALL R_ASSERT (0.0D0, OMG, 1.0D0, 'TOL', 'SET_DOD_TOL')
C
      SO_RPAR (IS_DOD_OMG + I - 1) = OMG
C
      RETURN
      END SUBROUTINE SET_DOD_OMG

      FUNCTION GET_DOD_MAXIT ()
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'slv_params.inc'
C
      INTEGER(KIND=IK) GET_DOD_MAXIT
C
      GET_DOD_MAXIT = SO_IPAR (IS_DOD_MAXIT)
C
      RETURN
      END FUNCTION GET_DOD_MAXIT

      SUBROUTINE SET_DOD_MAXIT (MAXIT)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'slv_params.inc'
C
      INTEGER(KIND=IK) MAXIT
C
      CALL I_ASSERT (0, MAXIT, 1000, 'MAXIT', 'SET_DOD_MAXIT')
      SO_IPAR (IS_DOD_MAXIT) = MAXIT
C
      RETURN
      END SUBROUTINE SET_DOD_MAXIT
