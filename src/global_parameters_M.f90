! Variables that all files may need access too
! For instance the outfile unit for writing to the main file
!

module global_parameters_M

    implicit none

!---Define floats 
    integer, parameter :: dp=kind(0.d0) !---double precision
    integer, parameter :: sp=kind(1.0)  !--- single precision

!---Nonlinear variables
    integer :: max_iter = 1 ! max num of nonlinear iterations to do
    integer :: max_nl_iter  ! numer of nonllinear iterations to do
    real(dp) :: tolerance = 0.001 ! prescribed tolerance
!---
    integer :: NumTimeSteps

!---Math constants    
    real(dp) :: pi
    parameter (pi = 3.1415926535897932_dp)
!---THIS ASSUMES 1 CM MESH SIZING
    real(dp), dimension(3) :: elem_vol_int_fe
    data elem_vol_int_fe / 0.16666667198838597,0.66666665602322805, &
                          0.16666667198838597 /
    
    integer :: beta_counter

end module global_parameters_M 
