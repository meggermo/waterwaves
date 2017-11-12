C ---------------------------------------------------------------------------
C     
C ---------------------------------------------------------------------------
      INTEGER(KIND=IK) IT_CUR
      INTEGER(KIND=IK) IT_BEG
      INTEGER(KIND=IK) IT_END
      INTEGER(KIND=IK) IT_STP
C
      PARAMETER       (IT_CUR = 1)
      PARAMETER       (IT_BEG = IT_CUR + 1)
      PARAMETER       (IT_END = IT_BEG + 1)
      PARAMETER       (IT_STP = IT_END + 1)
C      
      INTEGER(KIND=IK) NT_RPAR
      PARAMETER       (NT_RPAR = IT_STP)
C
      INTEGER(KIND=IK) I_PLOT
      INTEGER(KIND=IK) I_PLOT_MOD
      LOGICAL(KIND=LK) PROB_TIME_DEPENDENT
      LOGICAL(KIND=LK) GRID_TIME_DEPENDENT
      LOGICAL(KIND=LK) FIRST_TIME
      REAL   (KIND=RK) TM_RPAR (NT_RPAR)
C
      COMMON /TME_PARAMS/ TM_RPAR, 
     &                    I_PLOT,
     &                    I_PLOT_MOD,
     &                    PROB_TIME_DEPENDENT, 
     &                    GRID_TIME_DEPENDENT,
     &                    FIRST_TIME
C ---------------------------------------------------------------------------
C     
C ---------------------------------------------------------------------------
