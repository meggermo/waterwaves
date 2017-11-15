      SUBROUTINE CONTOUR_INTEGRAL
     &           (NSD, NNW_SD, NGP_SD, NW_IPAR, CRD, PHN)
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'fle_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'cfs_params.inc'
      INCLUDE 'slv_funcs.inc'
      INCLUDE 'tme_funcs.inc'
C
      INTEGER(KIND=IK) NSD
      INTEGER(KIND=IK) NNW_SD (*)
      INTEGER(KIND=IK) NGP_SD (*)
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)
      REAL   (KIND=RK) CRD (4, *)
      REAL   (KIND=RK) PHN (2, *)
C
      INTEGER(KIND=IK) ISD, INW, IGP
      INTEGER(KIND=IK) IWRK (LIM)
      REAL   (KIND=RK) TOL  (2)
      REAL   (KIND=RK) DWRK (LENW
      REAL   (KIND=RK) TI (NSD)
C
      WRITE (USR_O, *) '=== CONTOUR_INTEGRAL ==='
      TOL (1) = GET_INT_ABS ()
      TOL (2) = GET_INT_REL ()
C
      INW = 1
      IGP = 1
      WRITE (PLT_CTR, 1101) GET_TIME ()
      DO ISD = 1, NSD
         WRITE (USR_O, 1001) ISD
         CALL CONTOUR_INTEGRAL_SD
     &        (NNW_SD (ISD),  NW_IPAR (1, INW),
     &         TOL (1),       TOL (2),
     &         CRD (1, IGP),  PHN (1, IGP),
     &         TI (ISD),
     &         IWRK,          DWRK)
C
         INW = INW + NNW_SD (ISD)
         IGP = IGP + NGP_SD (ISD)
      END DO
      WRITE (PLT_CTR, 1101) TI
      WRITE (PLT_CTR, *)
C
      RETURN
 1001 FORMAT (1X, 'CONTOUR_INT OF SUBDOMAIN', I4, ':')
 1101 FORMAT (1X, E14.6,$)
      END
      SUBROUTINE CONTOUR_INTEGRAL_SD
     &           (NNW, NW_IPAR, EABS, EREL, CRD, PHN, TI, IWRK, DWRK)
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'fle_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'spl_params.inc'
      INCLUDE 'cfs_params.inc'
      INCLUDE 'net_funcs.inc'
C
      INTEGER(KIND=IK) NNW
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)
      REAL   (KIND=RK) EABS
      REAL   (KIND=RK) EREL
      REAL   (KIND=RK) CRD (4, *)
      REAL   (KIND=RK) PHN (2, *)
      REAL   (KIND=RK) TI
      INTEGER(KIND=IK) IWRK (*)
      REAL   (KIND=RK) DWRK (*)
C
      INTEGER(KIND=IK) INW, IGP, JGP, NGP, NEVAL, IERR, LAST, J
      REAL   (KIND=RK) CI
      REAL   (KIND=RK) ANS
      REAL   (KIND=RK) ERR
      REAL   (KIND=RK) CONTOUR
      EXTERNAL         CONTOUR
C
      IGP = 1
      TI  = 0.0D0
      DO INW = 1, NNW
         NGP = GET_NGP (NW_IPAR (1, INW))
         CI  = 0.0D0
         DO JGP = IGP, IGP + NGP - 2
            EL (1, 1) = PHN (1, JGP)
            EL (2, 1) = PHN (1, JGP + 1)
            EL (3, 1) = PHN (2, JGP)
            EL (4, 1) = PHN (2, JGP + 1)
            EL (1, 2) = CRD (1, JGP)
            EL (2, 2) = CRD (1, JGP + 1)
            EL (3, 2) = CRD (3, JGP)
            EL (4, 2) = CRD (3, JGP + 1)
            EL (1, 3) = CRD (2, JGP)
            EL (2, 3) = CRD (2, JGP + 1)
            EL (3, 3) = CRD (4, JGP)
            EL (4, 3) = CRD (4, JGP + 1)
            CALL DQAG
     &           (CONTOUR, XB, XE, EABS, EREL, KEY, ANS, ERR,
     &            NEVAL, IERR, LIM, LENW, LAST, IWRK, DWRK)
            IF (IERR .NE. 0) THEN
               WRITE (*, *) ' PHI_N:', (EL (J, 1), J = 1, 4)
               WRITE (*, *) ' X    :', (EL (J, 2), J = 1, 4)
               WRITE (*, *) ' Z    :', (EL (J, 3), J = 1, 4)
            END IF
            CI = CI + ANS
         END DO
         WRITE (USR_O,   1001) INW, CI
         WRITE (PLT_CTR, 1201) CI
         TI  = TI  + CI
         IGP = IGP + NGP
      END DO
      WRITE (USR_O, 1101) TI
      RETURN
 1001 FORMAT (1X,'Flow through network',I4,'=',E12.4)
 1101 FORMAT (1X,'Contour Integral =',E11.3)
 1201 FORMAT (1X,E14.6,$)
      END
      FUNCTION CONTOUR (S)
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'spl_params.inc'
C
      REAL   (KIND=RK) CONTOUR
      REAL   (KIND=RK) S
C
      REAL   (KIND=RK) W (4, 2)
      REAL   (KIND=RK) DOT_4
      EXTERNAL         DOT_4
C
      CALL WEIGHT (S, W)
C
      CONTOUR = SQRT (DOT_4 (W (1, 2), EL (1, 2)) ** 2
     &              + DOT_4 (W (1, 2), EL (1, 3)) ** 2)
     &              * DOT_4 (W (1, 1), EL (1, 1))
C
      RETURN
      END
