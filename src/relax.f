      SUBROUTINE RELAX (NSD, NNW_SD, NW_IPAR, ALPHA, RES,
     &           PHI, PHN)
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
      REAL   (KIND=RK) ALPHA (2)
      REAL   (KIND=RK) RES (2)
      REAL   (KIND=RK) PHI (2, *)
      REAL   (KIND=RK) PHN (2, *)
C
      INTEGER(KIND=IK) ISD, NNW
      INTEGER(KIND=IK) INW_A, ITF_A, IGP_A, NGP_A, TBC_A
      INTEGER(KIND=IK) INW_B, ITF_B, IGP_B, NGP_B, TBC_B
C
      NNW = 0
      DO ISD = 1, NSD
         NNW = NNW + NNW_SD (ISD)
      END DO
C
      RES (1) = 0.0D0
      RES (2) = 0.0D0
C
      IGP_A = 1
C
      DO INW_A = 1, NNW
         NGP_A = GET_NGP (NW_IPAR (1, INW_A))
         ITF_A = GET_ITF (NW_IPAR (1, INW_A))
         IF (ITF_A .NE. 0) THEN
            IGP_B = IGP_A
            DO INW_B = INW_A + 1, NNW
               ITF_B = GET_ITF (NW_IPAR (1, INW_B))
               IGP_B = IGP_B + GET_NGP (NW_IPAR (1, INW_B))
               IF (ITF_A .EQ. ITF_B) EXIT
            END DO
            IF (ITF_A .EQ. ITF_B) THEN
               NGP_B = GET_NGP (NW_IPAR (1, INW_B))
               TBC_A = GET_BCT (NW_IPAR (1, INW_A))
               TBC_B = GET_BCT (NW_IPAR (1, INW_B))
               IF (TBC_A .GT. 0) THEN
                  CALL RELAX_PHI (NGP_A, NGP_B, TBC_A, TBC_B,
     &                 ALPHA (1), RES, PHI (1,IGP_A), PHI (1,IGP_B))
               ELSE
                  CALL RELAX_PHN (NGP_A, NGP_B, TBC_A, TBC_B,
     &                 ALPHA (2), RES, PHN (1,IGP_A), PHN (1,IGP_B))
               END IF
               CALL SET_BCT (-TBC_A, NW_IPAR (1, INW_A))
               CALL SET_BCT (-TBC_B, NW_IPAR (1, INW_B))
            END IF
         END IF
         IGP_A = IGP_A + NGP_A
      END DO
C
      RETURN
      END
      SUBROUTINE RELAX_PHI (NGP_A, NGP_B, TBC_A, TBC_B,
     &                      ALPHA, RES, PHA, PHB)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      INTEGER(KIND=IK) NGP_A, NGP_B
      INTEGER(KIND=IK) TBC_A, TBC_B
      REAL(KIND=RK) ALPHA
      REAL(KIND=RK) RES (2)
      REAL(KIND=RK) PHA (2, NGP_A)
      REAL(KIND=RK) PHB (2, NGP_B)
C
      INTEGER(KIND=IK) I, J
      REAL(KIND=RK) PHA_N (2)
      REAL(KIND=RK) PHB_N (2)
      INTEGER(KIND=IK) BC (2)
      DAtA             BC /0, 0/
C
      RES (1) = 0.0D0
      RES (2) = 0.0D0
C
      WRITE (*, *) 'RELAX - PHI'
      DO I = 1, NGP_A
         J = NGP_B + 1 - I
         RES (1) = MAX (RES (1), ABS (PHA (1, I) - PHB (1, J)))
         RES (2) = MAX (RES (2), ABS (PHA (2, I) + PHB (2, J)))
         PHA_N (1) = (1.0D0 - ALPHA) * PHA (1, I)
     &                      + ALPHA  * PHB (1, J)
         PHB_N (1) = (1.0D0 - ALPHA) * PHB (1, J)
     &                      + ALPHA  * PHA (1, I)
         PHA_N (2) = (1.0D0 - ALPHA) * PHA (2, I)
     &                      - ALPHA  * PHB (2, J)
         PHB_N (2) = (1.0D0 - ALPHA) * PHB (2, J)
     &                      - ALPHA  * PHA (2, I)
         PHA (1, I) = PHA_N (1)
         PHB (1, J) = PHB_N (1)
         PHA (2, I) = PHA_N (2)
         PHB (2, J) = PHB_N (2)
      END DO
C     CALL SPLINE (BC, NGP_A, 1, 2, PHA (1, 1), PHA (2, 1))
C     CALL SPLINE (BC, NGP_B, 1, 2, PHB (1, 1), PHB (2, 1))
C
      RETURN
      END
      SUBROUTINE RELAX_PHN (NGP_A, NGP_B, TBC_A, TBC_B,
     &                      ALPHA, RES, PNA, PNB)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      INTEGER(KIND=IK) NGP_A, NGP_B
      INTEGER(KIND=IK) TBC_A, TBC_B
      REAL(KIND=RK) ALPHA
      REAL(KIND=RK) RES (2)
      REAL(KIND=RK) PNA (2, NGP_A)
      REAL(KIND=RK) PNB (2, NGP_B)
C
      INTEGER(KIND=IK) I, J
      REAL(KIND=RK) PNA_N (2)
      REAL(KIND=RK) PNB_N (2)
      INTEGER(KIND=IK) BC (2)
      DAtA             BC /0, 0/
C
      WRITE (*, *) 'RELAX - PHN'
      RES (1) = 0.0D0
      RES (2) = 0.0D0
C
      DO I = 1, NGP_A
         J = NGP_B + 1 - I
         RES (1) = MAX (RES (1), ABS (PNA (1, I) + PNB (1, J)))
         RES (2) = MAX (RES (2), ABS (PNA (2, I) - PNB (2, J)))
         PNA_N (1) = (1.0D0 - ALPHA) * PNA (1, I)
     &                      - ALPHA  * PNB (1, J)
         PNB_N (1) = (1.0D0 - ALPHA) * PNB (1, J)
     &                      - ALPHA  * PNA (1, I)
         PNA_N (2) = (1.0D0 - ALPHA) * PNA (2, I)
     &                      + ALPHA  * PNB (2, J)
         PNB_N (2) = (1.0D0 - ALPHA) * PNB (2, J)
     &                      + ALPHA  * PNA (2, I)
         PNA (1, I) = PNA_N (1)
         PNB (1, J) = PNB_N (1)
         PNA (2, I) = PNA_N (2)
         PNB (2, J) = PNB_N (2)
      END DO
C     CALL SPLINE (BC, 1, 2, NGP_A, PNA (1, 1), PNA (2, 1))
C     CALL SPLINE (BC, 1, 2, NGP_B, PNB (1, 1), PNB (2, 1))
C
      RETURN
      END
