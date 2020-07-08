        !COMPILER-GENERATED INTERFACE MODULE: Wed Jul 08 09:01:14 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WRITE_OUT_SS_SOLN__genmod
          INTERFACE 
            SUBROUTINE WRITE_OUT_SS_SOLN(FILE_UNIT,RANGE_ELEM,          &
     &TRANSIENT_SAVE)
              INTEGER(KIND=4), INTENT(IN) :: FILE_UNIT
              INTEGER(KIND=4), INTENT(IN) :: RANGE_ELEM
              LOGICAL(KIND=4), INTENT(IN) :: TRANSIENT_SAVE
            END SUBROUTINE WRITE_OUT_SS_SOLN
          END INTERFACE 
        END MODULE WRITE_OUT_SS_SOLN__genmod
