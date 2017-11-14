      SUBROUTINE CHECK_COEFS (NSD, NGP_SD, FPH, FPN, S0, S1, D0, D1)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
C
      INTEGER(KIND=IK) NSD
      INTEGER(KIND=IK) NGP_SD (*)
      REAL   (KIND=RK) FPH (2, *)
      REAL   (KIND=RK) FPN (2, *)
      REAL   (KIND=RK) S0 (*)
      REAL   (KIND=RK) S1 (*)
      REAL   (KIND=RK) D0 (*)
      REAL   (KIND=RK) D1 (*)
C
      INTEGER(KIND=IK) ISD, IGP, IAE
C
      WRITE (USR_O, *) '=== CHECK_COEFS ==='
C
      IGP = 1
      IAE = 1
      DO ISD = 1, NSD
         WRITE (USR_O, 1001) ISD
         CALL CHECK_SUBDOMAIN
     &        (NGP_SD (ISD),
     &         FPH (1, IGP), FPN (1, IGP),
     &         S0 (IAE),     S1 (IAE),
     &         D0 (IAE),     D1 (IAE))
         IGP = IGP + NGP_SD (ISD)
         IAE = IAE + NGP_SD (ISD) ** 2
      END DO
C
      RETURN
 1001 FORMAT (1X,'CHECK SUBDOMAIN',I4,':')
      END
      SUBROUTINE CHECK_SUBDOMAIN (NGP, FPH, FPN, S0, S1, D0, D1)
C ---------------------------------------------------------------------------
C     Computes the Max, L1 and L2 norm of the integral equations.
C ---------------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'knd_params.inc'
      INCLUDE 'usr_params.inc'
C
      INTEGER(KIND=IK) NGP
      REAL   (KIND=RK) FPH (2, *)
      REAL   (KIND=RK) FPN (2, *)
      REAL   (KIND=RK) S0  (NGP, *)
      REAL   (KIND=RK) S1  (NGP, *)
      REAL   (KIND=RK) D0  (NGP, *)
      REAL   (KIND=RK) D1  (NGP, *)
C
      INTEGER(KIND=IK) IGP, JGP
      REAL   (KIND=RK) S (NGP)
C
      INTEGER(KIND=IK) IDAMAX
      REAL   (KIND=RK) DNRM2, DASUM
      EXTERNAL         IDAMAX
      EXTERNAL         DNRM2, DASUM
C
      DO IGP = 1, NGP
         S (IGP) = 0.0D0
      END DO
C
      DO IGP = 1, NGP
         DO JGP = 1, NGP
            S (JGP) =  S (JGP)
     &              + D0 (JGP, IGP) * FPH (1, IGP)
     &              + D1 (JGP, IGP) * FPH (2, IGP)
     &              + S0 (JGP, IGP) * FPN (1, IGP)
     &              + S1 (JGP, IGP) * FPN (2, IGP)
         END DO
      END DO
C
      IGP = IDAMAX (NGP, S, 1)
C
      WRITE (USR_O, 101) ABS (S (IGP)), IGP
      WRITE (USR_O, 111) DASUM (NGP, S, 1)
      WRITE (USR_O, 121) DNRM2 (NGP, S, 1)
C
      RETURN
  101 FORMAT (1X,'L_MAX:',E9.1,' at ',I3)
  111 FORMAT (1X,'L_1  :',E9.1)
  121 FORMAT (1X,'L_2  :',E9.1)
      END