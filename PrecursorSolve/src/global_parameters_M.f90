! Variables that all files may need access too
! For instance the outfile unit for writing to the main file
!

module global_parameters_M

    implicit none

!---Define floats 
    integer, parameter:: dp=kind(0.d0)

!---File names
    character(60) :: file_name

!---File units    
    integer :: outfile_unit = 15 
    integer :: soln_outfile_unit = 99
    integer :: soln_last_t_unit = 66
    integer :: power_outfile_unit = 20
    integer :: power_file_unit = 12 
    integer :: beta_special_unit = 64

!---Nonlinear variables
    integer :: max_iter = 1 ! max num of nonlinear iterations to do
    integer :: max_nl_iter  ! numer of nonllinear iterations to do
    real(dp) :: residual
    real(dp) :: tolerance = 0.001 ! prescribed tolerance

!---Math constants    
    real(dp) :: pi
    parameter (pi = 3.1415926535897932_dp)
!---THIS ASSUMES 1 CM MESH SIZING
    real(dp), dimension(3) :: elem_vol_int_fe
    data elem_vol_int_fe / 0.16666667198838597,0.66666665602322805, &
                          0.16666667198838597 /
    


end module global_parameters_M 
