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
