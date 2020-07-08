        !COMPILER-GENERATED INTERFACE MODULE: Wed Jul 08 09:04:31 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DENSITY_FEEDBACK__genmod
          INTERFACE 
            SUBROUTINE DENSITY_FEEDBACK(N,J,CURRENT_TIME,DEN_FEEDBACK)
              INTEGER(KIND=4), INTENT(IN) :: N
              INTEGER(KIND=4), INTENT(IN) :: J
              INTEGER(KIND=4), INTENT(IN) :: CURRENT_TIME
              REAL(KIND=8), INTENT(INOUT) :: DEN_FEEDBACK
            END SUBROUTINE DENSITY_FEEDBACK
          END INTERFACE 
        END MODULE DENSITY_FEEDBACK__genmod
