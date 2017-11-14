C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      INTEGER(KIND=IK) N_NGP, N_BCT, N_PLT, N_GRT, N_GMT, N_GAT, N_SPT,
     &                 N_ITF, N_AED, N_ANW
      INTEGER(KIND=IK) I_NGP, I_BCT, I_PLT, I_GRT, I_GMT, I_GAT, I_SPT,
     &                 I_ITF, I_AED, I_ANW
      INTEGER(KIND=IK) N_IP
C
      PARAMETER (N_NGP = 1,
     &           N_BCT = 1,
     &           N_PLT = 1,
     &           N_GRT = 1,
     &           N_GMT = 1,
     &           N_GAT = 1,
     &           N_SPT = 6,
     &           N_ITF = 1,
     &           N_AED = 2,
     &           N_ANW = 2)
      PARAMETER (I_NGP = 1,
     &           I_BCT = I_NGP + N_NGP,
     &           I_PLT = I_BCT + N_BCT,
     &           I_GRT = I_PLT + N_PLT,
     &           I_GMT = I_GRT + N_GRT,
     &           I_GAT = I_GMT + N_GMT,
     &           I_SPT = I_GAT + N_GAT,
     &           I_ITF = I_SPT + N_SPT,
     &           I_AED = I_ITF + N_ITF,
     &           I_ANW = I_AED + N_AED,
     &           N_IP  = I_ANW + N_ANW - 1)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------
      INTEGER(KIND=IK) N_BCP, N_GRP, N_GMP, N_GAP
      INTEGER(KIND=IK) I_BCP, I_GRP, I_GMP, I_GAP
      INTEGER(KIND=IK) N_RP
C
      PARAMETER (N_BCP = 17,
     &           N_GRP = 10,
     &           N_GMP = 10,
     &           N_GAP = 10)
      PARAMETER (I_BCP = 1,
     &           I_GRP = I_BCP + N_BCP,
     &           I_GMP = I_GRP + N_GRP,
     &           I_GAP = I_GMP + N_GMP,
     &           N_RP  = I_GAP + N_GAP - 1)
C ---------------------------------------------------------------------------
C
C ---------------------------------------------------------------------------