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
C     CALL GET_TOKENS (USR_I, 4, LINE)
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
C     CALL GET_TOKENS (USR_I, 1, LINE)
      READ (USR_I, *) ITF
      CALL SET_ITF (ITF, NW_IPAR)
      ITF = GET_ITF (NW_IPAR)
      WRITE (USR_O, 111) 'INTERFACE ID', ITF
C     Read the plot type, if it is 0 then nothing will be written to file
C     CALL GET_TOKENS (USR_I, 1, LINE)
      READ (USR_I, *) PLT
      CALL SET_PLT (PLT, NW_IPAR)
      PLT = GET_PLT (NW_IPAR)
      WRITE (USR_O, 111) 'PLOT TYPE', PLT
      IF (PLT .GT. 0) THEN
         WRITE (STD_T,'(I5,$)')GET_NGP (NW_IPAR)
      END IF
C     Read the boundary condition type and kind
C     CALL GET_TOKENS (USR_I, 2, LINE)
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
C        CALL GET_TOKENS (USR_I, 1, LINE)
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
      ELSE
C        read parameters for this boundary condition
C        CALL GET_TOKENS (USR_I, 10, LINE)
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
C     CALL GET_TOKENS (USR_I, 1, LINE)
      READ (USR_I, *) GRT
      CALL SET_GRT (GRT, NW_IPAR)
      GRT = GET_GRT (NW_IPAR)
      WRITE (USR_O, 111) 'INITIAL GRID', GRT
C     And read the grid parameters
C     CALL GET_TOKENS (USR_I, N_GRP, LINE)
      READ (USR_I, *) (PAR (IRP), IRP = 1, N_GRP)
      CALL SET_GRP (PAR, NW_RPAR)
C     Read the kind of grid motion
C     CALL GET_TOKENS (USR_I, 1, LINE)
      READ (USR_I, *) GMT
      CALL SET_GMT (GMT, NW_IPAR)
      GMT = GET_GMT (NW_IPAR)
      WRITE (USR_O, 111) 'GRID MOTION', GMT
C     And read the grid motion parameters
C     CALL GET_TOKENS (USR_I, N_GMP, LINE)
      READ (USR_I, *) (PAR (IRP), IRP = 1, N_GMP)
      CALL SET_GMP (PAR, NW_RPAR)
C     Read the kind of grid adjustment
C     CALL GET_TOKENS (USR_I, 1, LINE)
      READ (USR_I, *) GAT
      CALL SET_GAT (GAT, NW_IPAR)
      GAT = GET_GAT (NW_IPAR)
      WRITE (USR_O, 111) 'GRID ADJUSTMENT', GAT
C     Read the kind of bc's for the splines
C     CALL GET_TOKENS (USR_I, 4, LINE)
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
