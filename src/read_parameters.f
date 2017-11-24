      SUBROUTINE READ_PARAMETERS
     &           (NSD, NNW_SD, NW_IPAR, NW_RPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      INTEGER(KIND=IK) NSD                ! nr. of subd. (IN)
      INTEGER(KIND=IK) NNW_SD (*)         ! nr. of netw. per subd. (IN)
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)  ! integer netw. par. (OUT)
      REAL(KIND=RK)    NW_RPAR (N_RP, *)  ! real netw. par. (OUT)
C
      REAL(KIND=RK) PERIOD
      DATA          PERIOD /1.0D0/
C
      CALL READ_ANALYTIC_PARAMETERS
     &     (PERIOD)
C
      CALL READ_TIME_PARAMETERS
     &     (PERIOD)
C
      CALL READ_SOLVE_PARAMETERS
     &     (NSD)
C
      CALL READ_CHK_PARAMETERS
C
      CALL READ_NETWORK_PARAMETERS
     &     (NSD, NNW_SD, PERIOD, NW_IPAR, NW_RPAR)
C
C     CALL OPEN_CDF_FILES (NSD, NNW_SD, NW_IPAR)
C
      RETURN
      END

      SUBROUTINE READ_ANALYTIC_PARAMETERS
     &           (PERIOD)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'anl_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'tme_params.inc'
C
      REAL(KIND=RK) PERIOD
C
      INTEGER(KIND=IK) I
      CHARACTER*40 MSG
C     Initialize analytic solution parameters to nil
      DATA ANL_PAR /N_RPAR * 0.0D0/
C
      WRITE (USR_O, *) '========== ANALYTIC PARAMETERS ===='
      READ  (USR_I, *) ANL_TYPE
      IF (ANL_TYPE .EQ. ANL_LINWAVE) THEN
C        Linear wave
         CALL READ_LINWAVE (ANL_PAR, PERIOD)
         GRID_TIME_DEPENDENT = .FALSE.
      ELSE IF (ANL_TYPE .EQ. ANL_RFWAVE) THEN
C        Rienecker & Fenton wave
         CALL READ_RFWAVE (ANL_PAR, PERIOD)
      ELSE IF (ANL_TYPE .EQ. ANL_POLYNOMIAL) THEN
C        Polynomial of the form (X + I * Z) ** K + C * T
         WRITE (USR_O, *) '    -> Polynomial'
         READ  (USR_I, *)    (ANL_PAR (I), I = 1, 5)
         WRITE (USR_O, 101) (ANL_PAR (I), I = 1, 5)
         GRID_TIME_DEPENDENT = .FALSE.
      ELSE IF (ANL_TYPE .EQ. ANL_NOTAVAIL) THEN
C        DO NOTHING, SINCE THERE'S NO KNOWN SOLUTION
      ELSE
         WRITE (MSG,'(A,I3)')
     &   'Analytic solution not implemented:', ANL_TYPE
         CALL ERROR ('ANALYTIC', MSG)
      END IF
C
      RETURN
  101 FORMAT (8X,'PARAMETERS:',/,8E14.6)
      END

      SUBROUTINE READ_TIME_PARAMETERS (PERIOD)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'tme_params.inc'
      INCLUDE 'tme_funcs.inc'
C
      REAL   (KIND=RK) PERIOD
C
      REAL   (KIND=RK) TB, TE, DT
C
      WRITE (USR_O, *) '========== TIME PARAMETERS ========'
      READ  (USR_I, *) TB
      READ  (USR_I, *) TE
      READ  (USR_I, *) DT
      READ  (USR_I, *) I_PLOT_MOD
C
      CALL SET_TIME    (PERIOD * TB)
      CALL SET_TBEG    (PERIOD * TB)
      CALL SET_TEND    (PERIOD * TE)
      CALL SET_DELTA_T (PERIOD * DT)
      IF (DT .EQ. 0.0D0) THEN
         PROB_TIME_DEPENDENT = .FALSE.
         CALL SET_TIME (0.0D0)
         CALL SET_TBEG (0.0D0)
         CALL SET_TEND (1.0D0)
         WRITE (USR_O, 101) 'DELTA-T', GET_DELTA_T ()
         WRITE (USR_O, *) 'TIME INDEPENDENT'
      ELSE
         WRITE (USR_O, 101) 'T-BEGIN', GET_TBEG    ()
         WRITE (USR_O, 101) 'T-END  ', GET_TEND    ()
         WRITE (USR_O, 101) 'DELTA-T', GET_DELTA_T ()
         WRITE (USR_O, 111) 'DT_PLT ', I_PLOT_MOD
      END IF
C
      RETURN
  101 FORMAT (1X, A8, ':', E10.2, ' s.')
  111 FORMAT (1X, A8, ':', I3, ' * DELTA-T')
      END

      SUBROUTINE READ_SOLVE_PARAMETERS (NSD)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'slv_params.inc'
      INCLUDE 'slv_funcs.inc'
C
      INTEGER(KIND=IK) NSD
C
      INTEGER(KIND=IK) MAX_IT
      REAL   (KIND=RK) EPS_ABS, EPS_REL
      REAL   (KIND=RK) TOL (2), OMG (2)
C
      WRITE (USR_O, *) '========== SOLVE PARAMETERS ======='
C
      READ (USR_I, *) EPS_ABS, EPS_REL
C
      CALL SET_INT_ABS (EPS_ABS)
      CALL SET_INT_REL (EPS_REL)
      WRITE (USR_O,   *) 'INTEGRATION PARAMETERS:'
      WRITE (USR_O, 101) 'INT_ABS', GET_INT_ABS ()
      WRITE (USR_O, 101) 'INT_REL', GET_INT_REL ()
C
      READ (USR_I, *) TOL, OMG, MAX_IT
      CALL SET_ITR_TOL (1, TOL (1))
      CALL SET_ITR_TOL (2, TOL (2))
      CALL SET_ITR_OMG (1, OMG (1))
      CALL SET_ITR_OMG (2, OMG (2))
      CALL SET_ITR_MAXIT (MAX_IT)
      WRITE (USR_O,   *) 'ITERATION PARAMETERS:'
      WRITE (USR_O, 101) 'ITR_TOL', GET_ITR_TOL (1), GET_ITR_TOL (2)
      WRITE (USR_O, 101) 'ITR_OMG', GET_ITR_OMG (1), GET_ITR_OMG (2)
      WRITE (USR_O, 111) 'ITR_MAX', GET_ITR_MAXIT ()
C
      IF (NSD .GT. 1) THEN
         READ (USR_I, *) TOL, OMG, MAX_IT
         CALL SET_DOD_TOL (1, TOL (1))
         CALL SET_DOD_TOL (2, TOL (2))
         CALL SET_DOD_OMG (1, OMG (1))
         CALL SET_DOD_OMG (2, OMG (2))
         CALL SET_DOD_MAXIT (MAX_IT)
         WRITE (USR_O,   *) 'DOMAIN DECOMPOSITION PARAMETERS:'
         WRITE (USR_O, 101) 'DOD_TOL', GET_DOD_TOL (1), GET_DOD_TOL (2)
         WRITE (USR_O, 101) 'DOD_OMG', GET_DOD_OMG (1), GET_DOD_OMG (2)
         WRITE (USR_O, 111) 'DOD_MAX', GET_DOD_MAXIT ()
      ELSE
         CALL SET_DOD_TOL (1, 0.0D0)
         CALL SET_DOD_TOL (2, 0.0D0)
         CALL SET_DOD_OMG (1, 0.0D0)
         CALL SET_DOD_OMG (2, 0.0D0)
         CALL SET_DOD_MAXIT (1)
      END IF
C
      RETURN
  101 FORMAT (1X,A8,':',2E14.4)
  111 FORMAT (1X,A8,':',I7)
      END

      SUBROUTINE READ_CHK_PARAMETERS
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'chk_params.inc'
C
      READ (USR_I, *) CHK_EQNS
      READ (USR_I, *) CHK_COEFS_PRE
      READ (USR_I, *) CHK_COEFS_POST
      READ (USR_I, *) CHK_CONTOUR
C
      RETURN
      END

      SUBROUTINE READ_NETWORK_PARAMETERS
     &           (NSD, NNW_SD, PERIOD, NW_IPAR, NW_RPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'fle_params.inc'
C
      INTEGER(KIND=IK) NSD               ! nr. of subd. (IN)
      INTEGER(KIND=IK) NNW_SD (*)        ! nr. of netw. per subd. (IN)
      REAL   (KIND=RK) PERIOD
      INTEGER(KIND=IK) NW_IPAR (N_IP, *) ! integer netw. par (OUT)
      REAL   (KIND=RK) NW_RPAR (N_RP, *) ! real netw. par. (OUT)
C
      INTEGER(KIND=IK) ISD, INW, NNW
C
      OPEN (UNIT= STD_T, FILE = 'dim.out', STATUS = 'UNKNOWN')
      INW = 1
      DO ISD = 1, NSD
         WRITE (USR_O, 101) ISD
         NNW = NNW_SD (ISD)
         CALL READ_NETW_PARAMETERS
     &        (NNW, PERIOD, NW_IPAR (1, INW), NW_RPAR (1, INW))
         CALL CHECK_NETW_PAR
     &        (ISD, NNW,    NW_IPAR (1, INW), NW_RPAR (1, INW))
         INW = INW + NNW
      END DO
      WRITE (STD_T, *)
      CLOSE (UNIT = STD_T)
C
      CALL CHECK_INTERFACES
     &     (NSD, NNW_SD, NW_IPAR)
C
  101 FORMAT (1X,'========== SUBDOMAIN',I3,' ==========')
      RETURN
      END
      SUBROUTINE READ_NETW_PARAMETERS
     &           (NNW, PERIOD, NW_IPAR, NW_RPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'net_params.inc'
C
      INTEGER(KIND=IK) NNW               ! nr. of networks (IN)
      REAL   (KIND=RK) PERIOD
      INTEGER(KIND=IK) NW_IPAR (N_IP, *) ! integer netw. params. (OUT)
      REAL   (KIND=RK) NW_RPAR (N_RP, *) ! real netw. params. (OUT)
C
      INTEGER(KIND=IK) INW
C
      DO INW = 1, NNW
         WRITE (USR_O, 101) INW
         CALL READ_NETW_PAR (PERIOD, NW_IPAR (1, INW), NW_RPAR (1, INW))
      END DO
C
      RETURN
  101 FORMAT (1X,'---------- NETWORK',I5,' ----------')
      END

      SUBROUTINE READ_NETW_PAR
     &           (PERIOD, NW_IPAR, NW_RPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'bct_params.inc'
      INCLUDE 'anl_params.inc'
      INCLUDE 'rle_params.inc'
      INCLUDE 'fle_params.inc'
      INCLUDE 'net_funcs.inc'
C
      REAL   (KIND=RK) PERIOD
      INTEGER(KIND=IK) NW_IPAR (N_IP)
      REAL   (KIND=RK) NW_RPAR (N_RP)
C
      CHARACTER        BCP
      INTEGER(KIND=IK) IED, IRP, ISP
      INTEGER(KIND=IK) AED (2), ANW (2), ITF, PLT, BCT, GRT, GMT, GAT
      INTEGER(KIND=IK) SPT (4)
      REAL   (KIND=RK) PAR (100)
      REAL   (KIND=RK) L, H, W, K, C, R
      CHARACTER*72     LINE
C     Read the topology information
      READ (USR_I, *) ANW (1), AED (1), ANW (2), AED (2)
      DO IED = 1, 2
        CALL SET_ANW (IED, ANW (IED), NW_IPAR)
        CALL SET_AED (IED, AED (IED), NW_IPAR)
        ANW (IED) = GET_ANW (IED, NW_IPAR)
        AED (IED) = GET_AED (IED, NW_IPAR)
        WRITE (USR_O, 101) IED, ANW (IED), AED (IED)
      END DO
      WRITE (USR_O, 111) 'NR. OF GRID POINTS', GET_NGP (NW_IPAR)
C     Read the interface ID, if it is 0 then it is not an interface
      READ (USR_I, *) ITF
      CALL SET_ITF (ITF, NW_IPAR)
      ITF = GET_ITF (NW_IPAR)
      WRITE (USR_O, 111) 'INTERFACE ID', ITF
C     Read the plot type, if it is 0 then nothing will be written to file
      READ (USR_I, *) PLT
      CALL SET_PLT (PLT, NW_IPAR)
      PLT = GET_PLT (NW_IPAR)
      WRITE (USR_O, 111) 'PLOT TYPE', PLT
      IF (PLT .GT. 0) THEN
         WRITE (STD_T,'(I5,$)')GET_NGP (NW_IPAR)
      END IF
C     Read the boundary condition type and kind
      READ (USR_I, *) BCT, BCP
      IF (BCP .EQ. 'D') THEN
         CALL SET_BCT (BCK_DIRICHLET * ABS (BCT), NW_IPAR)
         WRITE (USR_O, 111) 'BOUNDARY CONDITION', BCT
      ELSE IF (BCP .EQ. 'N') THEN
         CALL SET_BCT (BCK_NEUMANN   * ABS (BCT), NW_IPAR)
      ELSE IF (BCP .EQ. 'M') THEN
         CALL SET_BCT (BCK_MIXED     * ABS (BCT), NW_IPAR)
      ELSE
         CALL ERROR ('READ_NETW_PAR', 'BC. Kind D, N or M expected')
      END IF
      BCT = GET_BCT (NW_IPAR)
      WRITE (USR_O, 141) 'BOUNDARY CONDITION', BCT, BCP
      DO IRP = 1, N_BCP
         PAR (IRP) = 0.0D0
      END DO

      IF (BCT .EQ. BCT_WAVEMAKER) THEN
         IF (BCP .NE. 'N') THEN
            CALL ERROR ('READ_NETW_PAR','wavemaker MUST be Neumann')
         END IF
C        read the filename of the wavemaker data file
         LINE = ' '
         READ (USR_I,*) LINE
         OPEN (UNIT = WVM_IN,
     &         FILE = LINE(1:LEN_TRIM(LINE)), STATUS = 'OLD')
         WRITE (USR_O, *) '      wavemaker input file: ',
     &                    LINE (1:LEN_TRIM (LINE))
C        read initial data from wavemaker file
         CALL INITIALIZE (4, WVM_IN, PAR)
C        read if it attached to begin or to end of network
         LINE = ' '
         READ (USR_I,*)  LINE
C        b-> begin, e-> end
         IF (LINE (1:1) .EQ. 'b' .OR. LINE (1:1) .EQ. 'B') THEN
            WRITE (USR_O, *) '      origin at begin'
            PAR (17) = 1.0D0
         ELSE IF (LINE (1:1) .EQ. 'e' .OR. LINE (1:1) .EQ. 'E') THEN
            WRITE (USR_O, *) '      origin at end'
            PAR (17) = -1.0D0
         ELSE
            CALL ERROR ('READ_NETW_PAR','keyword begin or end expected')
         END IF
      ELSE IF (BCT .EQ. BCT_PISTON) THEN
        IF (BCP .NE. 'N') THEN
           CALL ERROR ('READ_NETW_PAR','piston MUST be Neumann')
        END IF
        READ  (USR_I, *) PAR (1:3)
        WRITE (USR_O, *) 'PISTON AMPLITUDE:', PAR (1)
        WRITE (USR_O, *) 'PISTON FREQUENCY:', PAR (2)
        WRITE (USR_O, *) 'PISTON FACTOR:', PAR (3)
      ELSE
C        read parameters for this boundary condition
         READ (USR_I, *) (PAR (IRP), IRP = 1, 10)
C        Scale C to that of the analytic solution
         IF (ABS (BCT) .EQ. BCT_SOMMERFELD) THEN
C           This parameter will only have meaning for R&F waves
C           otherwise it will be 0.0 and will have no effect
            PAR (2) = 0.0D0
C           Scale C to the analytic solution (if known)
            IF (ANL_TYPE .EQ. ANL_RFWAVE) THEN
               H   = ANL_PAR (3)
               W   = ANL_PAR (5)
               K   = ANL_PAR (6)
               C   = W / K
               R   = ANL_PAR (7)
               PAR (1) = C * PAR (1)
C              See page 48 of PdH's thesis
               PAR (2) =  (GRAV * (R - H) - 0.5D0 * C ** 2)
            ELSE IF (ANL_TYPE .EQ. ANL_LINWAVE) THEN
               L  = ANL_PAR (2)
               H  = ANL_PAR (3)
               K  = TWO_PI / L
               W  = SQRT (Grav * K * TANH (K * H))
               C  = W / K
               PAR (1) = C * PAR (1)
            END IF
         END IF
      END IF
C     Copy the processed data to the bc parameters
      CALL SET_BCP (PAR, NW_RPAR)
C     Read how the initial grid must be constructed
      READ (USR_I, *) GRT
      CALL SET_GRT (GRT, NW_IPAR)
      GRT = GET_GRT (NW_IPAR)
      WRITE (USR_O, 111) 'INITIAL GRID', GRT
C     And read the grid parameters
      READ (USR_I, *) (PAR (IRP), IRP = 1, N_GRP)
      CALL SET_GRP (PAR, NW_RPAR)
C     Read the kind of grid motion
      READ (USR_I, *) GMT
      CALL SET_GMT (GMT, NW_IPAR)
      GMT = GET_GMT (NW_IPAR)
      WRITE (USR_O, 111) 'GRID MOTION', GMT
C     And read the grid motion parameters
      READ (USR_I, *) (PAR (IRP), IRP = 1, N_GMP)
      CALL SET_GMP (PAR, NW_RPAR)
C     Read the kind of grid adjustment
      READ (USR_I, *) GAT
      CALL SET_GAT (GAT, NW_IPAR)
      GAT = GET_GAT (NW_IPAR)
      WRITE (USR_O, 111) 'GRID ADJUSTMENT', GAT
C     Read the kind of bc's for the splines
      READ (USR_I, *) SPT
      DO ISP = 1, 4
         CALL SET_SPT (ISP, SPT (ISP), NW_IPAR)
      END DO
      WRITE (USR_O, 131) 'SPLINE-PHI BC', SPT (1), SPT (2)
      WRITE (USR_O, 131) 'SPLINE-PHN BC', SPT (3), SPT (4)
C     Write the boundary condition parameters back to stdout
      CALL GET_BCP (NW_RPAR, PAR)
      WRITE (USR_O, 121)
     &   'BOUNDARY COND. PARAMETERS', (PAR (IRP), IRP = 1, N_BCP)
C     Write the initial grid parameters back to stdout
      CALL GET_GRP (NW_RPAR, PAR)
      WRITE (USR_O, 121)
     &   'INITIAL GRID PARAMETERS', (PAR (IRP), IRP = 1, N_GRP)
C     Write the gri motion parameters back to stdout
      CALL GET_GMP (NW_RPAR, PAR)
      WRITE (USR_O, 121)
     &   'GRID MOTION PARAMETERS', (PAR (IRP), IRP = 1, N_GMP)
C
      RETURN
  101 FORMAT (1X,'EDGE',I2,':    ','NW =',I3,', EDGE =',I2)
  111 FORMAT (1X,A18,':',I5)
  121 FORMAT (1X,A,  ':',/,(5E14.6))
  131 FORMAT (1X,A18,':',2I5)
  141 FORMAT (1X,A18,':',I5,' (TYPE = ',A,')')
      END

      SUBROUTINE READ_LINWAVE (PAR, PERIOD)
C ---------------------------------------------------------------------------
C     PAR (1) : Wave Amplitude
C     PAR (2) : Wave Length
C     PAR (3) : Water depth
C     PAR (4) : Wave crest offset
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'rle_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'tme_params.inc'
C
      REAL   (KIND=RK) PAR (*)
      REAL   (KIND=RK) PERIOD
C
      INTEGER(KIND=IK) I
      REAL   (KIND=RK) K, H
C
      READ (USR_I, *) (PAR (I), I = 1, 4)
      K  = TWO_PI / PAR (2)
      H  = PAR (3)
      PERIOD = TWO_PI / SQRT (Grav * K * TANH (K * H))
C
      WRITE (USR_O, *) '    -> Linear wave (grid time-independent)'
      WRITE (USR_O, '(1X,A,F10.4)') '    Wave length:', PAR (2)
      WRITE (USR_O, '(1X,A,F10.4)') '    Wave period:', PERIOD
      WRITE (USR_O, '(1X,A,F10.4)') '    Wave height:', 2.0D0 * PAR (1)
      WRITE (USR_O, '(1X,A,F10.4)') '    Water depth:', H
      WRITE (USR_O, '(1X,A,F10.4)') '    Crest offs.:', PAR (4)
C
      GRID_TIME_DEPENDENT = .FALSE.
C
      RETURN
      END

      SUBROUTINE READ_RFWAVE (PAR, PERIOD)
C-----------------------------------------------------------------------------
C     Reads a the parameters of a Rienecker & Fenton typ of wave.
C-----------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'fle_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'rle_params.inc'
C
      REAL(KIND=RK) PAR (33, 3)
      REAL(KIND=RK) PERIOD
C
      LOGICAL(KIND=LK) ldummy
      INTEGER(KIND=IK) i, n
      CHARACTER*32     File_Name
      REAL   (KIND=RK) waterdepth, waveheigth, rho, gravity,
     &                 phasevelocity, Q, R, rdummy, omega, wavelength
C
      WRITE (USR_O, *)              '    -> R & F wave'
      READ  (USR_I,  *)     File_Name
      OPEN  (STD_T, FILE = File_Name, STATUS= 'OLD')
      WRITE (USR_O, '(1X,A,A    )') '    Data file:  ', File_Name
      READ  (STD_T, *) waterdepth
      WRITE (USR_O, '(1X,A,F10.4)') '    Water depth:', waterdepth
      READ  (STD_T, *) PERIOD
      WRITE (USR_O, '(1X,A,F10.4)') '    Wave period:', period
      READ  (STD_T, *) waveheigth
      WRITE (USR_O, '(1X,A,F10.4)') '    Wave height:', waveheigth
      READ  (STD_T, *) ldummy
      READ  (STD_T, *) rdummy
      READ  (STD_T, *) n
      WRITE (USR_O, '(1X,A,  I10)') '    Wave  modes:', n
      READ  (STD_T, *) rho
      READ  (STD_T, *) gravity
      READ  (STD_T, *) phasevelocity
      READ  (STD_T, *) Q
      READ  (STD_T, *) R
      omega     = 2.0D0 * PI / PERIOD
      wavelength= period * phasevelocity
      par (1, 1) = waveheigth
      par (2, 1) = wavelength
      par (3, 1) = waterdepth
      READ (USR_I, *)
     +par (4, 1)
      WRITE (USR_O, '(1X,A,F10.4)') '    Crest offs.:', par (4,1)
      par (5, 1) = omega
      par (6, 1) = TWO_PI / period / phasevelocity
      par (7, 1) = R
      par (8, 1) = n
      DO i= 1, n + 1
        READ (STD_T, *) rdummy, PAR (i, 3), PAR (i, 2)
      END DO
      CLOSE (STD_T)
C
      RETURN
      END
