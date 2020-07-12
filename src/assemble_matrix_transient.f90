!*****************************************************************************80
!
!! Assemble the elemental to solve for coefficients the transient 
!  Discussion:
!             As the elemental matrices are created we add them to the global 
!
!  Input: n - element number
!
!  Output:
! 

subroutine assemble_matrix_transient (isotope,delay_group,n)
   
    USE flags_M
    USE element_matrices_M
    USE global_parameters_M
    USE material_info_M
    USE mesh_info_M
    USE solution_vectors_M
    USE mesh_info_M
    USE time_info_M 
    Use parameters_fe
    
    implicit none
!---Dummy variables
    integer,intent(in) :: isotope
    integer,intent(in) :: delay_group
    integer,intent(in) :: n

!---Local variables
    integer :: i, j, length   
    real(dp), dimension(3) ::  beta_lambda_times_q_vec, &
                              W_left_times_prev_elem_soln_vec
    ! H_times_soln_vec,

    elem_matrix_H = 0.0_dp 
!---Calculate [U - lambda*A - W_r]
    H_times_soln_vec = matmul(elem_matrix_U,precursor_soln_prev(isotope,delay_group,n,:)) - &
        matmul(lamda_i_mat(isotope,delay_group)*elem_matrix_A,precursor_soln_prev(isotope,delay_group,n,:)) - &
        matmul(matrix_W_right_face,precursor_soln_prev(isotope,delay_group,n,:))

    beta_lambda_times_q_vec = 0.0_dp
!---Calculate Beta/Gen Time * {q}
    do i = 1, nodes_per_elem
        beta_lambda_times_q_vec(i) = power_amplitude_prev*&
                                    total_power_initial*spatial_power_frac_fcn(n,i)*&
                           (beta_i_mat(isotope,delay_group)/gen_time)*elem_vec_q(i) 
    end do

    W_left_times_prev_elem_soln_vec = 0.0_dp
!---Calculate W_l*c_e-l
    if ( n > 1) then
        do i = 1, nodes_per_elem
            do j = 1, nodes_per_elem
                W_left_times_prev_elem_soln_vec(i) = W_left_times_prev_elem_soln_vec(i) + &
                                                 matrix_W_left_face(i,j)*&
                             precursor_soln_new(isotope,delay_group,n-1,3)
            end do
        end do 
    else
        do i = 1, nodes_per_elem
            do j = 1, nodes_per_elem
                W_left_times_prev_elem_soln_vec(i) = W_left_times_prev_elem_soln_vec(i) + &
                                                 matrix_W_left_face(i,j)*&
                             precursor_soln_new(isotope,delay_group,num_elem,3)
            end do
        end do
    end if

!---RHS final
    do i = 1, nodes_per_elem
        RHS_transient_final_vec(i) = H_times_soln_vec(i) + &
                                     beta_lambda_times_q_vec(i) +&
                                     W_left_times_prev_elem_soln_vec(i) 
    end do

!    DEBUG = .TRUE.
    if (DEBUG .eqv. .TRUE.) then
       
       write(outfile_unit,fmt='(a)'),' '
       write(outfile_unit,fmt='(a,4I6)'),'[U] element Matrix gaussian integration &
                                          | element --> ',n
       do j=1,nodes_per_elem 
              write(outfile_unit,fmt='(12es14.3)') &
                   (elem_matrix_U(j,i),i=1,nodes_per_elem)             
       end do
       
       write(outfile_unit,fmt='(a)'), ' '
       write(outfile_unit,fmt='(a,I6)'),'[W] right Matrix | element --> ',n
       do j=1,nodes_per_elem 
             write(outfile_unit,fmt='(12es14.3)') &
                  (matrix_W_right_face(j,i),i=1,nodes_per_elem)             
       end do
       
       write(outfile_unit,fmt='(a)'), ' '
       write(outfile_unit,fmt='(a,4I6)'),'[W] left Matrix | element --> ',n
       do j=1,nodes_per_elem 
             write(outfile_unit,fmt='(12es14.3)') &
                  (matrix_W_left_face(j,i),i=1,nodes_per_elem)             
       end do
       write(outfile_unit,fmt='(a,4I6)'),'{q} vector  --> ',n
             write(outfile_unit,fmt='(12es14.3)') &
                  (elem_vec_q(i),i=1,nodes_per_elem)             
 
        write(outfile_unit,fmt='(a,12es14.3)'), '@ time = ', t0
        
        write(outfile_unit,fmt='(a)'), ' '
        write(outfile_unit,fmt='(a,1I2)'),'H Matrix | element --> ',n
        do i=1,nodes_per_elem
              write(outfile_unit,fmt='(12es14.3)') &
                   (elem_matrix_H(i,j),j=1,nodes_per_elem)
        end do

        write(outfile_unit,fmt='(a)'), ' '
        write(outfile_unit,fmt='(a,1I3)'),' [H]*{c_e} | element --> ', n
        do j=1,nodes_per_elem
              write(outfile_unit,fmt='(12es14.3)')H_times_soln_vec(j)
        end do
       
        write(outfile_unit,fmt='(a)'), ' '
        write(outfile_unit,fmt='(a,1I3)'),' [W_l]*{c_e-1} | element --> ', n
        do j=1,nodes_per_elem
              write(outfile_unit,fmt='(12es14.3)')W_left_times_prev_elem_soln_vec(j)
        end do

        write(outfile_unit,fmt='(a)'),' '
        write(outfile_unit,fmt='(a,1I2)'),'{q} element source vector | element --> ',n
        write(outfile_unit,fmt='(12es14.3)') (beta_lambda_times_q_vec(i),i=1,nodes_per_elem)

        write(outfile_unit,fmt='(a)'),' '
        write(outfile_unit,fmt='(a,1I2)'),'RHS vector | element --> ',n
        write(outfile_unit,fmt='(12es14.3)') (RHS_transient_final_vec(i),i=1,nodes_per_elem)

        write(outfile_unit,fmt='(a)'),'*****************************************'
       
    end if
!DEBUG = .FALSE.

end 

