module element_matrices_M

    USE global_parameters_M

    implicit none

!---Elemental matrices
    real(dp) , dimension(3,3) :: elem_matrix_A_times_W
    real(dp) , dimension(3,3) :: identity_matrix
    real(dp) , dimension(3,3) :: inverse_A_matrix
    real(dp) , dimension(3,3) :: analytic_elem_matrix_P_ss
    real(dp) , dimension(3,3) :: analytic_elem_matrix_A
    real(dp) , dimension(3,3) :: elem_matrix_U
    real(dp) , dimension(3,3) :: elem_matrix_A
    real(dp) , dimension(3,3) :: elem_matrix_G
    real(dp) , dimension(3,3) :: elem_matrix_F
    real(dp) , dimension(3,3) :: elem_matrix_H
    real(dp) , dimension(3,3) :: matrix_W_left_face
    real(dp) , dimension(3,3) :: matrix_W_right_face
    real(dp) , dimension(3)   :: interp_fcn_rhs, interp_fcn_lhs
    real(dp) , dimension(3) :: vol_int

    data interp_fcn_lhs / 1,0,0/
    data interp_fcn_rhs / 0,0,1/
    data identity_matrix /1,0,0,&
                         0,1,0,&
                         0,0,1/
    data vol_int /0.166666666667, 0.6666666666667, 0.166666666667/

end module element_matrices_M
