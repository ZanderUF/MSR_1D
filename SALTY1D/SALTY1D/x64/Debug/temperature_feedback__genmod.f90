        !COMPILER-GENERATED INTERFACE MODULE: Wed Jul 08 09:04:31 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE TEMPERATURE_FEEDBACK__genmod
          INTERFACE 
            SUBROUTINE TEMPERATURE_FEEDBACK(N,NL_ITER,CURRENT_TIME,     &
     &TEMP_REACTIVITY_FEEDBACK)
              INTEGER(KIND=4), INTENT(IN) :: N
              INTEGER(KIND=4), INTENT(IN) :: NL_ITER
              INTEGER(KIND=4), INTENT(IN) :: CURRENT_TIME
              REAL(KIND=8), INTENT(OUT) :: TEMP_REACTIVITY_FEEDBACK
            END SUBROUTINE TEMPERATURE_FEEDBACK
          END INTERFACE 
        END MODULE TEMPERATURE_FEEDBACK__genmod
