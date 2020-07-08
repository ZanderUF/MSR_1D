        !COMPILER-GENERATED INTERFACE MODULE: Wed Jul 08 09:01:16 2020
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SOLVE_POWER_EULER__genmod
          INTERFACE 
            SUBROUTINE SOLVE_POWER_EULER(NL_ITER,CURRENT_TIME)
              USE MATERIAL_INFO_M
              USE MESH_INFO_M
              INTEGER(KIND=4), INTENT(IN) :: NL_ITER
              REAL(KIND=8), INTENT(IN) :: CURRENT_TIME
            END SUBROUTINE SOLVE_POWER_EULER
          END INTERFACE 
        END MODULE SOLVE_POWER_EULER__genmod
