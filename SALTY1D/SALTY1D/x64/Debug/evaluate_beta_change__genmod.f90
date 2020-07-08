        !COMPILER-GENERATED INTERFACE MODULE: Wed Jul 08 09:01:17 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE EVALUATE_BETA_CHANGE__genmod
          INTERFACE 
            SUBROUTINE EVALUATE_BETA_CHANGE(EVENT_TIME,                 &
     &EVENT_TIME_PREVIOUS,EVENT_COUNTER,EVENT_OCCURING)
              USE MATERIAL_INFO_M
              REAL(KIND=8), INTENT(IN) :: EVENT_TIME
              REAL(KIND=8), INTENT(IN) :: EVENT_TIME_PREVIOUS
              INTEGER(KIND=4), INTENT(INOUT) :: EVENT_COUNTER
              LOGICAL(KIND=4), INTENT(IN) :: EVENT_OCCURING
            END SUBROUTINE EVALUATE_BETA_CHANGE
          END INTERFACE 
        END MODULE EVALUATE_BETA_CHANGE__genmod
