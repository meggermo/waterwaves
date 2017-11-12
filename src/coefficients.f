      SUBROUTINE COEFFICIENTS
     &           (NSD, NNW_SD, NGP_SD, NW_IPAR, CRD, S0, S1, D0, D1)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'usr_params.inc'
      INCLUDE 'chk_params.inc'
C
      INTEGER(KIND=IK) NSD                ! nr. of subdomains (IN)
      INTEGER(KIND=IK) NNW_SD (*)         ! nr. of netw. per subd. (IN)
      INTEGER(KIND=IK) NGP_SD (*)         ! nr. of gridp. per subd. (IN)
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)  ! integer netw. par. (IN)
      REAL   (KIND=RK) CRD (4, *)         ! grid points (IN)
      REAL   (KIND=RK) S0 (*)             ! Source coeffs. (OUT)
      REAL   (KIND=RK) S1 (*)             ! Source coeffs. (OUT)
      REAL   (KIND=RK) D0 (*)             ! Dipole coeffs. (OUT)
      REAL   (KIND=RK) D1 (*)             ! Dipole coeffs. (OUT)
C
      INTEGER(KIND=IK) ISD
      INTEGER(KIND=IK) INW_SD             ! network       index pointer
      INTEGER(KIND=IK) IGP_SD             ! grid point    index pointer
      INTEGER(KIND=IK) IAE_SD             ! array element index pointer
      INTEGER(KIND=IK) NNW, NGP

C     Initialize the index pointers
      INW_SD = 1
      IGP_SD = 1
      IAE_SD = 1
C
      DO ISD = 1, NSD
C        Get the number of networks and grid points of the subdomain
         NNW = NNW_SD (ISD)
         NGP = NGP_SD (ISD)
C        Compute the source and dipole coefficients of this subdomain
         CALL COMP_COEF_SD
     &        (NNW, NGP,          
     &         NW_IPAR (1, INW_SD),
     &         CRD (1, IGP_SD),
     &         S0 (IAE_SD),  S1 (IAE_SD),
     &         D0 (IAE_SD),  D1 (IAE_SD))
C        Modify the diagonal elements so that the rowsum is 0
         CALL DIAG_COEF_SD 
     &        (NGP, D0 (IAE_SD))
C        Write the source and dipole coefficients to file (if needed)
         IF (CHK_COEFS_WRITE) THEN 
            CALL WRITE_COEFS_SD (ISD, NGP, S0, S1, D0, D1)
         END IF
C        Increment the index pointers
         INW_SD = INW_SD + NNW
         IGP_SD = IGP_SD + NGP
         IAE_SD = IAE_SD + NGP ** 2
      END DO
C
      RETURN
      END

      SUBROUTINE COMP_COEF_SD
     &           (NNW, NGP, NW_IPAR, CRD, SRC_0, SRC_1, DIP_0, DIP_1)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'net_params.inc'
      INCLUDE 'cfs_params.inc'
      INCLUDE 'net_funcs.inc'
      INCLUDE 'slv_funcs.inc'
C
      INTEGER(KIND=IK) NNW                ! nr. of networks (IN)
      INTEGER(KIND=IK) NGP                ! nr. of grid points (IN)
      INTEGER(KIND=IK) NW_IPAR (N_IP, *)  ! integer netw. par. (IN)
      REAL   (KIND=RK) CRD (4, *)         ! grid points (IN)
      REAL   (KIND=RK) SRC_0 (NGP, *)     ! Source coeffs
      REAL   (KIND=RK) SRC_1 (NGP, *)     ! Source coeffs
      REAL   (KIND=RK) DIP_0 (NGP, *)     ! Dipole coeffs
      REAL   (KIND=RK) DIP_1 (NGP, *)     ! Dipole coeffs
C
      INTEGER(KIND=IK) INW
      INTEGER(KIND=IK) IGP_NW
      INTEGER(KIND=IK) NGP_NW
      REAL   (KIND=RK) TOL (2)
      REAL   (KIND=RK) DWRK (LENW)
      INTEGER(KIND=IK) IWRK (LIM)
      INTEGER(KIND=IK) BCT, SPT
      INTEGER(KIND=IK) JGP, KGP, IED
      REAL   (KIND=RK) CRD_TMP (4, NGP)
      REAL   (KIND=RK) W (8)
      REAL   (KIND=RK) ALPHA, EPS, DELTA (2), JAC, XI
      PARAMETER       (ALPHA = 5.0D-1)
      PARAMETER       (EPS   = 1.0D-3)
C
      CALL DCOPY (4 * NGP, CRD, 1, CRD_TMP, 1)
C     
      IGP_NW = 1
C     
      DO INW = 1, NNW
         NGP_NW = GET_NGP (NW_IPAR (1, INW))
         BCT    = GET_BCT (NW_IPAR (1, INW))
         IF (BCT .LT. 0) THEN
            DO IED = 1, 2
               JGP = IGP_NW + (IED - 1) * (NGP_NW - 2)
               KGP = IGP_NW + (IED - 1) * (NGP_NW - 1)
               SPT = GET_SPT (2 + IED, NW_IPAR (1, INW))
               IF (SPT .NE. 0) THEN
                  JAC = SQRT (CRD (3, KGP) ** 2 + CRD (4, KGP) ** 2)
                  DELTA (1) = 0.0D0
                  DELTA (2) = -JAC * EPS
                  CALL LOCAL_TO_GLOBAL (1, 2, CRD (1,KGP), DELTA, DELTA)
                  XI = ALPHA * DBLE (3  - 2 * IED) + DBLE (IED - 1)
                  CALL WEIGHT (XI, W)
                  CRD_TMP (1, KGP) =
     &               W (1) * CRD  (1, JGP) + W (2) * CRD (1, JGP + 1)
     &             + W (3) * CRD  (3, JGP) + W (4) * CRD (3, JGP + 1)
     &             + DELTA (1)
                  CRD_TMP (2, KGP) =
     &               W (1) * CRD  (2, JGP) + W (2) * CRD (2, JGP + 1)
     &             + W (3) * CRD  (4, JGP) + W (4) * CRD (4, JGP + 1)
     &             + DELTA (2)
                  WRITE (*, *) KGP, ':', 
     &            CRD_TMP (1, KGP), CRD_TMP (2, KGP), XI
               END IF
            END DO
         END IF
         IGP_NW = IGP_NW + NGP_NW
      END DO
C     
      TOL (1) = GET_INT_ABS ()
      TOL (2) = GET_INT_REL ()
C      
      IGP_NW = 1
C
      DO INW = 1, NNW
         NGP_NW = GET_NGP (NW_IPAR (1, INW))
         DO JGP = 1, NGP
            SRC_0 (JGP, IGP_NW) = 0.0D0
            SRC_1 (JGP, IGP_NW) = 0.0D0
            DIP_0 (JGP, IGP_NW) = 0.0D0
            DIP_1 (JGP, IGP_NW) = 0.0D0
         END DO
C        WRITE (*, *) 'NETWORK ' , INW
         CALL COMP_COEF_NW
     &        (NGP,               NGP_NW, 
     &         TOL,
     &         CRD_TMP,           CRD   (1, IGP_NW),
     &         SRC_0 (1, IGP_NW), SRC_1 (1, IGP_NW),
     &         DIP_0 (1, IGP_NW), DIP_1 (1, IGP_NW), 
     &         DWRK,              IWRK)
         IGP_NW = IGP_NW + NGP_NW
      END DO
C
      RETURN
      END

      SUBROUTINE DIAG_COEF_SD (NGP, DIP)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'rle_params.inc'
C
      INTEGER(KIND=IK) NGP
      REAL   (KIND=RK) DIP (NGP, *)
C
      INTEGER(KIND=IK) I, J
      REAL   (KIND=RK) ROW_SUM (NGP)
C
      DO I = 1, NGP
         ROW_SUM (I) = 0.0D0
      END DO
C
      DO J = 1, NGP
         DO I = 1, NGP
            ROW_SUM (I) = ROW_SUM (I) + DIP (I, J)
         END DO
      END DO
C
      DO I = 1, NGP
         DIP (I, I) = DIP (I, I) - ROW_SUM (I)
      END DO
C
      RETURN
      END

      SUBROUTINE COMP_COEF_NW 
     &   (NGP,   NGP_NW, 
     &    TOL,
     &    CRD,   CRD_NW, 
     &    SRC_0, SRC_1,
     &    DIP_0, DIP_1, 
     &    DWRK,  IWRK)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'spl_params.inc'
C
      INTEGER(KIND=IK) NGP
      INTEGER(KIND=IK) NGP_NW
      REAL   (KIND=RK) TOL    (*)
      REAL   (KIND=RK) CRD    (4,   *)
      REAL   (KIND=RK) CRD_NW (4,   *)
      REAL   (KIND=RK) SRC_0  (NGP, *)
      REAL   (KIND=RK) SRC_1  (NGP, *)
      REAL   (KIND=RK) DIP_0  (NGP, *)
      REAL   (KIND=RK) DIP_1  (NGP, *)
      REAL   (KIND=RK) DWRK   (*)
      INTEGER(KIND=IK) IWRK (*)
C
      INTEGER(KIND=IK) I, J
C
      DO I = 1, NGP_NW - 1
C        WRITE(*,*) 'ELEMENT ', I
C        WRITE(*,*) '(X_B,Y_B) =  ', CRD_NW (1, I), CRD_NW (2, I)
C        WRITE(*,*) '(X_E,Y_E) =  ', CRD_NW (1, I+1), CRD_NW (2, I+1)
C        WRITE(*,*) '(DXB,DYB) =  ', CRD_NW (3, I), CRD_NW (4, I)
C        WRITE(*,*) '(DXE,DYE) =  ', CRD_NW (3, I+1), CRD_NW (4, I+1)
         DO J = 0, 1
            EL (1 + J, 1) = CRD_NW (1, I + J)
            EL (1 + J, 2) = CRD_NW (2, I + J)
            EL (3 + J, 1) = CRD_NW (3, I + J)
            EL (3 + J, 2) = CRD_NW (4, I + J)
         END DO
         CALL COMP_COEF_EL
     &        (NGP,
     &         TOL,
     &         CRD,
     &         SRC_0 (1, I), SRC_1 (1, I),
     &         DIP_0 (1, I), DIP_1 (1, I), DWRK, IWRK)
      END DO
C
      RETURN
      END

      FUNCTION IS_FLAT (EL)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      LOGICAL(KIND=LK) IS_FLAT
      REAL   (KIND=RK) EL (4, *)
C
      REAL   (KIND=RK) EPS
      PARAMETER    (EPS = 1.0D-12)
C
      REAL   (KIND=RK) DX_ELM (2), LEN_ELM
      REAL   (KIND=RK) DX_BEG (2), LEN_BEG
      REAL   (KIND=RK) DX_END (2), LEN_END
C
      DX_ELM (1) = EL (2, 1) - EL (1, 1)
      DX_ELM (2) = EL (2, 2) - EL (1, 2)
      LEN_ELM    = SQRT (DX_ELM (1) ** 2 + DX_ELM (2) ** 2)
      DX_ELM (1) = DX_ELM (1) / LEN_ELM
      DX_ELM (2) = DX_ELM (2) / LEN_ELM

      LEN_BEG    = SQRT (EL (3, 1) ** 2 + EL (3, 2) ** 2)
      DX_BEG (1) = EL (3, 1) / LEN_BEG
      DX_BEG (2) = EL (3, 2) / LEN_BEG

      LEN_END    = SQRT (EL (4, 1) ** 2 + EL (4, 2) ** 2)
      DX_END (1) = EL (4, 1) / LEN_END
      DX_END (2) = EL (4, 2) / LEN_END

C     WRITE (*, *) DX_ELM
C     WRITE (*, *) DX_BEG
C     WRITE (*, *) DX_END

      IS_FLAT = ABS (DX_ELM (1) - DX_BEG (1)) .LT. EPS
     &    .AND. ABS (DX_ELM (2) - DX_BEG (2)) .LT. EPS
     &    .AND. ABS (DX_ELM (1) - DX_END (1)) .LT. EPS
     &    .AND. ABS (DX_ELM (2) - DX_END (2)) .LT. EPS

      RETURN
      END

      FUNCTION IS_SINGULAR (P, EL)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
C
      LOGICAL(KIND=LK) IS_SINGULAR
      REAL   (KIND=RK) P  (*)
      REAL   (KIND=RK) EL (4, *)
C
      REAL   (KIND=RK) EPS
      PARAMETER       (EPS = 1.0D-12)
C
      IS_SINGULAR =
     &  SQRT ((P (1) - EL (1, 1)) ** 2
     &      + (P (2) - EL (1, 2)) ** 2) .LT. EPS .OR.
     &  SQRT ((P (1) - EL (2, 1)) ** 2
     &      + (P (2) - EL (2, 2)) ** 2) .LT. EPS
      RETURN
      END

      SUBROUTINE COMP_COEF_EL 
     &           (NGP, 
     &            TOL,   
     &            CRD, 
     &            SRC_0, SRC_1, 
     &            DIP_0, DIP_1,
     &            DWRK,  IWRK)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'spl_params.inc'
C
      INTEGER(KIND=IK) NGP
      REAL   (KIND=RK) TOL   (*)
      REAL   (KIND=RK) CRD   (4,   *)
      REAL   (KIND=RK) SRC_0 (NGP, *)
      REAL   (KIND=RK) SRC_1 (NGP, *)
      REAL   (KIND=RK) DIP_0 (NGP, *)
      REAL   (KIND=RK) DIP_1 (NGP, *)
      REAL   (KIND=RK) DWRK  (*)
      INTEGER(KIND=IK) IWRK (*)
C
      INTEGER(KIND=IK) I
      LOGICAL(KIND=LK) FLAT_PANEL
      REAL   (KIND=RK) SRC (4)
      REAL   (KIND=RK) DIP (4)
C
      LOGICAL(KIND=LK) IS_FLAT, IS_SINGULAR
      EXTERNAL         IS_FLAT, IS_SINGULAR
C
      FLAT_PANEL = IS_FLAT (EL)
C
      DO I = 1, NGP
C        Copy the field point coordinates into P
         P (1) = CRD (1, I)
         P (2) = CRD (2, I)
C        Compute the coefficients
         IF (IS_SINGULAR (P, EL)) THEN
            CALL COEF_SINGULAR
     &         (FLAT_PANEL, TOL (1), TOL (2), SRC, DIP, DWRK, IWRK)
         ELSE
            CALL COEF_REGULAR
     &         (FLAT_PANEL, TOL (1), TOL (2), SRC, DIP, DWRK, IWRK)
         END IF
C        And store them in the Source and Dipole matrices
         SRC_0 (I, 1) = SRC_0 (I, 1) + SRC (1)
         SRC_1 (I, 1) = SRC_1 (I, 1) + SRC (1 + 2)
         DIP_0 (I, 1) = DIP_0 (I, 1) + DIP (1)
         DIP_1 (I, 1) = DIP_1 (I, 1) + DIP (1 + 2)
         SRC_0 (I, 2) = SRC (2)
         SRC_1 (I, 2) = SRC (2 + 2)
         DIP_0 (I, 2) = DIP (2)
         DIP_1 (I, 2) = DIP (2 + 2)
      END DO
C
      RETURN
      END

      SUBROUTINE COEF_SINGULAR (FLAT, EABS, EREL, SRC, DIP, DWRK, IWRK)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'cfs_params.inc'
      INCLUDE 'spl_params.inc'
C
      LOGICAL(KIND=LK) FLAT
      REAL   (KIND=RK) EABS
      REAL   (KIND=RK) EREL
      REAL   (KIND=RK) SRC  (*)
      REAL   (KIND=RK) DIP  (*)
      REAL   (KIND=RK) DWRK (*)
      INTEGER(KIND=IK) IWRK (*)
C
      INTEGER(KIND=IK) NEVAL
      INTEGER(KIND=IK) IERR
      INTEGER(KIND=IK) LAST
      REAL   (KIND=RK) ERR
      REAL   (KIND=RK) SOURCE_COEF, DIPOLE_COEF
      EXTERNAL         SOURCE_COEF, DIPOLE_COEF
C
      IF (FLAT) THEN
C        A flat panel has dipole coefficients 0.0
         DO K = 1, 4
            CALL DQAGS (SOURCE_COEF, XB, XE, EABS, EREL, SRC (K),
     &                  ERR, NEVAL, IERR, LIM, LENW, LAST, IWRK, DWRK)
            IF (IERR .NE. 0) THEN
               WRITE (*, *) 'K:', K, ' EL:', EL
            END IF
            SRC (K) = 0.5D0 * SRC (K)
            DIP (K) = 0.0D0
         END DO
      ELSE
         DO K = 1, 4
            CALL DQAGS (DIPOLE_COEF, XB, XE, EABS, EREL, DIP (K),
     &                  ERR, NEVAL, IERR, LIM, LENW, LAST, IWRK, DWRK)
            IF (IERR .NE. 0) THEN
               WRITE (*, *) 'K:', K, ' EL:', EL
            END IF
            CALL DQAGS (SOURCE_COEF, XB, XE, EABS, EREL, SRC (K),
     &                  ERR, NEVAL, IERR, LIM, LENW, LAST, IWRK, DWRK)
            IF (IERR .NE. 0) THEN
               WRITE (*, *) 'K:', K, ' EL:', EL
            END IF
            SRC (K) = 0.5D0 * SRC (K)
         END DO
      END IF
C
      RETURN
      END

      SUBROUTINE COEF_REGULAR (FLAT, EABS, EREL, SRC, DIP, DWRK, IWRK)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'cfs_params.inc'
      INCLUDE 'spl_params.inc'
C
      LOGICAL(KIND=LK) FLAT
      REAL   (KIND=RK) EABS
      REAL   (KIND=RK) EREL
      REAL   (KIND=RK) SRC  (*)
      REAL   (KIND=RK) DIP  (*)
      REAL   (KIND=RK) DWRK (*)
      INTEGER(KIND=IK) IWRK (*)
C
      INTEGER(KIND=IK) NEVAL
      INTEGER(KIND=IK) IERR
      INTEGER(KIND=IK) LAST
      REAL   (KIND=RK) ERR
      REAL   (KIND=RK) SOURCE_COEF, DIPOLE_COEF
      EXTERNAL         SOURCE_COEF, DIPOLE_COEF
C
      DO K = 1, 4
         CALL DQAG (DIPOLE_COEF, XB, XE, EABS, EREL, KEY, DIP (K),
     &              ERR, NEVAL, IERR, LIM, LENW, LAST, IWRK, DWRK)
         IF (IERR .NE. 0) THEN
            WRITE (*, *) 'K:', K, ' EL:', EL, ' P:', P
         END IF
         CALL DQAG (SOURCE_COEF, XB, XE, EABS, EREL, KEY, SRC (K),
     &              ERR, NEVAL, IERR, LIM, LENW, LAST, IWRK, DWRK)
         IF (IERR .NE. 0) THEN
            WRITE (*, *) 'K:', K, ' EL:', EL, ' P:', P
         END IF
         SRC (K) = 0.5D0 * SRC (K)
      END DO
C
      RETURN
      END
