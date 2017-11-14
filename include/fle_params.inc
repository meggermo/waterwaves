C ---------------------------------------------------------------------------
C     fle_params.inc:
C     Contains parameters for unit numbers of in and output files
C ---------------------------------------------------------------------------
      INTEGER(KIND=IK) STD_I    ! Standard input stream  (SYSTEM DEPENDENT)
      INTEGER(KIND=IK) STD_O    ! Standard output stream (SYSTEM DEPENDENT)
      INTEGER(KIND=IK) STD_T    ! For temporary files
      INTEGER(KIND=IK) PLT_I    ! Main input file
      INTEGER(KIND=IK) PLT_O    ! Main output file
      INTEGER(KIND=IK) PLT_CRD  ! coords  output file
      INTEGER(KIND=IK) PLT_PHI  ! Phi     output file
      INTEGER(KIND=IK) PLT_PHN  ! Phn     output file
      INTEGER(KIND=IK) PLT_NRM  ! Norms   output file
      INTEGER(KIND=IK) PLT_CTR  ! Contour output file
      INTEGER(KIND=IK) WVM_IN   ! Contour output file
C
      PARAMETER (STD_I = 5)     ! NOTE: This parameter is system dependent
      PARAMETER (STD_O = 6)     ! NOTE: This parameter is system dependent
C
      PARAMETER (STD_T = STD_O + 1)
      PARAMETER (PLT_I = STD_T + 1)
      PARAMETER (PLT_O = PLT_I + 1)
      PARAMETER (PLT_CRD = PLT_O   + 10)
      PARAMETER (PLT_PHI = PLT_CRD + 10)
      PARAMETER (PLT_PHN = PLT_PHI + 10)
      PARAMETER (PLT_NRM = PLT_PHN + 10)
      PARAMETER (PLT_CTR = PLT_NRM + 10)
      PARAMETER (WVM_IN  = PLT_CTR + 10)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------