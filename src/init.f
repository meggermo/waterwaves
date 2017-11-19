      SUBROUTINE INIT
     &           (NSD, NNW, NGP, NAE, NNW_SD, NGP_SD, NGP_NW)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      INTEGER(KIND=IK) NSD              ! Nr. of subdomains     (OUT)
      INTEGER(KIND=IK) NNW              ! Nr. of networks       (OUT)
      INTEGER(KIND=IK) NGP              ! Nr. of grid points    (OUT)
      INTEGER(KIND=IK) NAE              ! Nr. of array elements (OUT)
      INTEGER(KIND=IK) NNW_SD (*)       ! Nr. of networks per subdomain (OUT)
      INTEGER(KIND=IK) NGP_SD (*)       ! Nr. of grid points per subd.  (OUT)
      INTEGER(KIND=IK) NGP_NW (*)       ! Nr. of grid points per netw,  (OUT)
C     Initialize the common block variables
      CALL INIT_COMMON_BLOCKS
C     Process the command line arguments
      CALL PROG_ARGUMENTS
C     Open some files
      CALL OPEN_OUTPUT_FILES
C     Read the number of subdomains
      CALL READ_NSD
     &     (NSD)
C     Read the number of networks for each subdomain
      CALL READ_NNW_SD
     &     (NSD, NNW_SD)
C     Compute the total number of networks
      CALL COMPUTE_NNW
     &     (NSD, NNW_SD, NNW)
C     Read the number of grid points for all networks
      CALL READ_NGP_NW
     &     (NNW, NGP_NW)
C     Compute the number of grid points of each subdomain
      CALL COMPUTE_NGP_SD
     &     (NSD, NNW_SD, NGP_NW, NGP_SD)
C     Compute the total number of grid points and number of array elements
      CALL COMPUTE_NGP_NAE
     &     (NSD, NGP_SD, NGP, NAE)
C
      RETURN
      END

      SUBROUTINE INIT_COMMON_BLOCKS
C ---------------------------------------------------------------------------
C     Initializes common block variables:
C
C     USR_I: the user input stream is set to the standard input stream
C     USR_O: the user output stream is set to the standard output stream
C
C     FIRST_TIME: set to TRUE, this means that we assume that we'll
C     start computations without any knowledge of previous solutions (so
C     we're not resuming previous computations)
C
C     PROB_TIME_DEPENDENT: set to TRUE, this means that we'll assume
C     that (some) boundary conditions are time-dependent
C
C     GRID_TIME_DEPENDENT: We assume that the problem involves a
C     time-dependent grid. Depending on the user input it can be changed
C     to time independent e.g. when we're solving linearized water
C     wave-problems.
C
C     CHK_EQNS: set to FALSE. When it is set to true, then the equations
C     solved in the routine iterate are checked.
C
C     CHK_COEFS_PRE: set to FALSE. When set to TRUE the integral
C     equations are checked to see if the coefficients are correctly
C     computed.
C
C     CHK_COEFS_POST: set to FALSE. When set to TRUE the integral
C     equations are checked to see if the coefficients are correctly
C     computed.
C
C     CHK_PARAMS: set to TRUE. Not used yet
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C     usr_params.inc contains the integers USR_I and USR_O
      INCLUDE 'usr_params.inc'
C     tme_params.inc contains FIRST_TIME and GRID_INDEPENDENT
      INCLUDE 'tme_params.inc'
C     fle_params.inc contains the integer parameters STD_I, STD_O
      INCLUDE 'fle_params.inc'
C     chk_params contains CHECK parameters to control program checks
      INCLUDE 'chk_params.inc'
C     cdf_params contains netCDF parameters for writing netCDF files
      INCLUDE 'cdf_params.inc'
C
      CHARACTER     C
      COMMON /SKIP/ C
C
      INTEGER I
C
      USR_I = STD_I
      USR_O = STD_O
C
      I_PLOT     = 0
      I_PLOT_MOD = 1
      PROB_TIME_DEPENDENT = .TRUE.
      GRID_TIME_DEPENDENT = .TRUE.
      FIRST_TIME          = .TRUE.
C
      CHK_EQNS        = .FALSE.
      CHK_COEFS_PRE   = .FALSE.
      CHK_COEFS_POST  = .FALSE.
      CHK_COEFS_WRITE = .FALSE.
      CHK_CONTOUR     = .TRUE.
      CHK_PARAMS      = .TRUE.
C
      C = ' '
C
      DO I = 1, CDF_MAX
         CDF_ID (I) = -1
      END DO
C
      RETURN
      END

      SUBROUTINE PROG_ARGUMENTS
C ---------------------------------------------------------------------------
C     Processed the program arguments:
C
C     ARG 1: The name of the user input file, if given, then USR_I is
C     redirected to this file instead of the standard input stream.
C
C     ARG 2: The name of the user output file, if given, then USR_O is
C     redirected to this file instead of the standard output stream.
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'fle_params.inc'
      INCLUDE 'chk_params.inc'
C
      CHARACTER*128 FILE_NAME
      CHARACTER*8   CHK
C
      IF (COMMAND_ARGUMENT_COUNT () .GT. 0) THEN
         USR_I = PLT_I
         CALL GETARG (1, FILE_NAME)
         OPEN (USR_I, FILE = FILE_NAME, STATUS = 'OLD')
      END IF
C
      IF (COMMAND_ARGUMENT_COUNT () .GT. 1) THEN
         CALL GETARG (2, FILE_NAME)
         IF (FILE_NAME (1:3) .NE. 'CHK') THEN
            USR_O = PLT_O
            OPEN (USR_O, FILE = FILE_NAME, STATUS = 'UNKNOWN')
         ELSE
            CHK_EQNS       = .TRUE.
            CHK_COEFS_PRE  = .TRUE.
            CHK_COEFS_POST = .TRUE.
         END IF
      END IF
C
      IF (COMMAND_ARGUMENT_COUNT () .GT. 2) THEN
         CALL GETARG (3, CHK)
         IF (CHK(1:3) .EQ. 'CHK') THEN
            CHK_EQNS       = .TRUE.
            CHK_COEFS_PRE  = .TRUE.
            CHK_COEFS_POST = .TRUE.
         END IF
      END IF
C
      END

      SUBROUTINE OPEN_OUTPUT_FILES
C ---------------------------------------------------------------------------
C     Opens the following output files:
C
C     plt.out: the output file containing the solution
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'fle_params.inc'
C     The following three files are the main data files
      OPEN (UNIT = PLT_CRD, FILE = 'crd.out', STATUS = 'UNKNOWN')
      OPEN (UNIT = PLT_PHI, FILE = 'phi.out', STATUS = 'UNKNOWN')
      OPEN (UNIT = PLT_PHN, FILE = 'phn.out', STATUS = 'UNKNOWN')
C     the norms are written to
      OPEN (UNIT = PLT_NRM, FILE = 'nrm.out', STATUS = 'UNKNOWN')
C     the contour integrals are written to
      OPEN (UNIT = PLT_CTR, FILE = 'ctr.out', STATUS = 'UNKNOWN')
C
      END

      SUBROUTINE READ_NSD
     &           (NSD)
C ---------------------------------------------------------------------------
C     Reads the number of subdomains and checks if the value is within
C     the allowed range.
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'max_params.inc'
C
      INTEGER(KIND=IK) NSD
C     Read the number of subdomains
      READ (USR_I, *) NSD
C     and check if the value is withing the allowed range
      CALL I_ASSERT (1, NSD, NSD_MAX, 'NSD', 'READ_NSD')
C
      RETURN
      END

      SUBROUTINE READ_NNW_SD
     &           (NSD, NNW_SD)
C ---------------------------------------------------------------------------
C     Reads the number of networks for each subdomain.
C     Validity checking is done in the subroutine COMPUTE_NNW_SD.
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
C
      INTEGER(KIND=IK) NSD
      INTEGER(KIND=IK) NNW_SD (*)
C
      INTEGER(KIND=IK) ISD, N
C     Read the number of networks for each subdomain
      READ (USR_I, *) (NNW_SD (ISD), ISD = 1, NSD)
C
      RETURN
      END

      SUBROUTINE COMPUTE_NNW
     &           (NSD, NNW_SD, NNW)
C ---------------------------------------------------------------------------
C     Computes the total number of networks and checks if the input is
C     within the allowed range.
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'max_params.inc'
C
      INTEGER(KIND=IK) NSD              ! nr. of subdomains (IN)
      INTEGER(KIND=IK) NNW_SD (*)       ! nr. of netw. per  subd. (IN)
      INTEGER(KIND=IK) NNW              ! total nr. of netw. (OUT)
C
      INTEGER(KIND=IK) ISD, NNS
C
      NNW = 0
      DO ISD = 1, NSD
         NNS = NNW_SD (ISD)
         CALL I_ASSERT (1, NNS, NNW_MAX, 'NNS', 'COMPUTE_NNW')
         NNW = NNW + NNW_SD (ISD)
         CALL I_ASSERT (1, NNW, NNW_MAX, 'NNW', 'COMPUTE_NNW')
      END DO
C
      RETURN
      END

      SUBROUTINE READ_NGP_NW
     &           (NNW, NGP_NW)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
C
      INTEGER(KIND=IK) NNW
      INTEGER(KIND=IK) NGP_NW (*)
C
      INTEGER(KIND=IK) INW
C
      READ (USR_I, *) (NGP_NW (INW), INW = 1, NNW)
C
      RETURN
      END

      SUBROUTINE COMPUTE_NGP_SD
     &           (NSD, NNW_SD, NGP_NW, NGP_SD)
C ---------------------------------------------------------------------------
C     Computes the number of grid points per subdomain and checks if the
C     given values are within the valid range.
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'max_params.inc'
C
      INTEGER(KIND=IK) NSD              ! nr. of subdomains (IN)
      INTEGER(KIND=IK) NNW_SD (*)       ! nr. of netw. per subd. (IN)
      INTEGER(KIND=IK) NGP_NW (*)       ! nr. of gridp. per netw. (IN)
      INTEGER(KIND=IK) NGP_SD (*)       ! nr. of gridp. per subd. (OUT)
C
      INTEGER(KIND=IK) ISD, INW, JNW, NNW, NGP, NGN
C
      INW = 1
      DO ISD = 1, NSD
         NNW = NNW_SD (ISD)
         NGP = 0
         DO JNW = INW, INW + NNW - 1
            NGN = NGP_NW (JNW)
            CALL I_ASSERT (1, NGN, NGP_MAX, 'NGN', 'COMPUTE_NGP_SD')
            NGP = NGP + NGN
            CALL I_ASSERT (1, NGP, NGP_MAX, 'NGP', 'COMPUTE_NGP_SD')
         END DO
         NGP_SD (ISD) = NGP
         INW = INW + NNW
      END DO
C
      RETURN
      END

      SUBROUTINE COMPUTE_NGP_NAE
     &           (NSD, NGP_SD, NGP, NAE)
C ---------------------------------------------------------------------------
C     Computes the total number of grid points and the total number of
C     array elements
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      INTEGER(KIND=IK) NSD          ! nr. of subdomains (IN)
      INTEGER(KIND=IK) NGP_SD (*)   ! nr. of grid points per subd. (IN)
      INTEGER(KIND=IK) NGP          ! total nr. of grid points (OUT)
      INTEGER(KIND=IK) NAE          ! total nr. of array elems. (OUT)
C
      INTEGER(KIND=IK) ISD
C
      NGP = 0
      NAE = 0
      DO ISD = 1, NSD
         NGP = NGP + NGP_SD (ISD)
         NAE = NAE + NGP_SD (ISD) ** 2
      END DO
C
      RETURN
      END

      SUBROUTINE CHECK_INTERFACES
     &           (NSD, NNW_SD, NW_IPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'net_funcs.inc'
C
      INTEGER(KIND=IK) NSD
      INTEGER(KIND=IK) NNW_SD (*)
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)
C
      INTEGER(KIND=IK) ISD, INW, NNW, NTF
C
      NNW = 0
      DO ISD = 1, NSD
         NNW = NNW + NNW_SD (ISD)
      END DO
      NTF = 0
      DO INW = 1, NNW
         IF (GET_ITF (NW_IPAR (1, INW)) .NE. 0) THEN
            NTF = NTF + 1
         END IF
      END DO
      IF (MOD (NTF, 2) .NE. 0) THEN
         CALL ERROR ('CHECK_INTERFACES', 'Number of ITF''s is odd')
      END IF
      CALL CHECK_INTERF (NNW, NTF, NW_IPAR)
      RETURN
      END
      SUBROUTINE CHECK_INTERF
     &           (NNW, NTF, NW_IPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'net_funcs.inc'
C
      INTEGER(KIND=IK) NNW
      INTEGER(KIND=IK) NTF
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)
C
      INTEGER(KIND=IK) INW, ITF, JTF, N_ITFS, TMP
      INTEGER(KIND=IK) ITF_ID (NTF)
      INTEGER(KIND=IK) INW_ID (NTF)
C
      N_ITFS = 0
      DO INW = 1, NNW
         ITF = GET_ITF (NW_IPAR (1, INW))
         IF (ITF .NE. 0) THEN
            N_ITFS = N_ITFS + 1
            ITF_ID (N_ITFS) = ITF
            INW_ID (N_ITFS) = INW
         END IF
      END DO
      DO ITF = 1, N_ITFS
         DO JTF = ITF + 1,  N_ITFS
            IF (ITF_ID (JTF) .LT. ITF_ID (ITF)) THEN
               TMP = ITF_ID (JTF)
               ITF_ID (JTF) = ITF_ID (ITF)
               ITF_ID (ITF) = TMP
               TMP = INW_ID (JTF)
               INW_ID (JTF) = INW_ID (ITF)
               INW_ID (ITF) = TMP
            END IF
         END DO
      END DO
      DO ITF = 1, N_ITFS, 2
C        IF (ITF_ID (ITF) .NE. ITF_ID (ITF + 1)) THEN
            WRITE (USR_O, 101) ITF_ID (ITF),     INW_ID (ITF)
            WRITE (USR_O, 101) ITF_ID (ITF + 1), INW_ID (ITF + 1)
C           CALL ERROR ('CHECK_INTERF', 'WRONG INTERFACE ID''s')
C        END IF
      END DO
C
      RETURN
  101 FORMAT (1X,'AT INTERFACE',I2,' OF NETWORK',I4)
      END

      SUBROUTINE CHECK_NETW_PAR
     &           (ISD, NNW, NW_IPAR, NW_RPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
C
      INTEGER(KIND=IK) ISD
      INTEGER(KIND=IK) NNW
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)
      REAL   (KIND=RK) NW_RPAR (N_RP, *)
C
      INTEGER(KIND=IK) INW, IED
C
      DO INW = 1, NNW
         DO IED = 1, 2
            CALL CHECK_CONNECTIONS
     &           (ISD, INW, IED, NNW, NW_IPAR, NW_RPAR)
         END DO
      END DO
C
      RETURN
      END
      SUBROUTINE CHECK_CONNECTIONS
     &           (ISD, INW, IED, NNW, NW_IPAR, NW_RPAR)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'net_funcs.inc'
C
      INTEGER(KIND=IK) ISD
      INTEGER(KIND=IK) INW
      INTEGER(KIND=IK) IED
      INTEGER(KIND=IK) NNW
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)
      REAL   (KIND=RK) NW_RPAR (N_RP, *)
C
      INTEGER(KIND=IK) ANW, AED, BNW, BED
      REAL   (KIND=RK) PAR_I (N_GRP)
      REAL   (KIND=RK) PAR_A (N_GRP)
C
      ANW = GET_ANW (IED, NW_IPAR (1, INW))
      IF (ANW .GT. NNW) THEN
         WRITE (USR_O, 101) ISD
         WRITE (USR_O, 111) IED, INW, ANW
         CALL ERROR ('CHECK_NETW_IPAR', 'Network ID out of range')
      END IF
      AED = GET_AED (IED, NW_IPAR (1, INW))
      BNW = GET_ANW (AED, NW_IPAR (1, ANW))
      BED = GET_AED (AED, NW_IPAR (1, ANW))
      IF (INW .NE. BNW) THEN
         WRITE (USR_O, 101) ISD
         WRITE (USR_O, 121) IED, INW, ANW
         WRITE (USR_O, 121) AED, ANW, BNW
         CALL ERROR ('CHECK_NETW_IPAR', 'Network ID mismatch')
      END IF
      IF (IED .NE. BED) THEN
         WRITE (USR_O, 101) ISD
         WRITE (USR_O, 121) IED, INW, ANW
         WRITE (USR_O, 121) AED, ANW, BNW
         CALL ERROR ('CHECK_NETW_IPAR', 'Edge ID mismatch')
      END IF
      CALL GET_GRP (NW_RPAR (1, INW), PAR_I)
      CALL GET_GRP (NW_RPAR (1, ANW), PAR_A)
      IF (NNW .NE. 1) THEN
      IF (PAR_I (1 + 2 * (IED - 1)) .NE. PAR_A (1 + 2 * (AED -1)) .OR.
     &    PAR_I (2 + 2 * (IED - 1)) .NE. PAR_A (2 + 2 * (AED -1))) THEN
         WRITE (USR_O, 131) IED, INW, PAR_I (1 + 2 * (IED -1)),
     &                                PAR_I (2 + 2 * (IED -1))
         WRITE (USR_O, *) 'BUT'
         WRITE (USR_O, 131) AED, ANW, PAR_A (1 + 2 * (AED -1)),
     &                                PAR_A (2 + 2 * (AED -1))
         CALL ERROR ('CHECK_NETW_IPAR', 'Coordinate mismatch')
      END IF
      END IF
C
      RETURN
  101 FORMAT ('IN SUBDOMAIN',I2)
  111 FORMAT ('EDGE',I2,' OF NETW',I4,' CANNOT CONNECT TO NETW',I4)
  121 FORMAT ('EDGE',I2,' OF NETW',I4,' CONNECTS TO NETW',I4)
  131 FORMAT ('EDGE',I2,' OF NETW',I4,' HAS COORDS',10E14.6)
      END
