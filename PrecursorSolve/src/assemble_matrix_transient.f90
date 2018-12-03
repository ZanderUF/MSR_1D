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

    implicit none
!---Dummy variables
    integer,intent(in) :: isotope
    integer,intent(in) :: delay_group
    integer,intent(in) :: n
!---Local variables
    integer :: i, j, ii, jj, nr,nc, ncl,length   
    double precision, dimension(3,3) :: inverse_matrix,temp_matrix
    double precision, dimension(3) :: elem_vec_q_times_beta_lambda, U_times_soln_vec, &
                          A_times_lambda_times_soln_vec, W_left_times_upwind_soln, &
                          W_right_times_soln

!---Initialize
    elem_matrix_H                   = 0.0_dp
    A_times_W_times_upwind_elem_vec = 0.0_dp
    elem_matrix_A_times_W           = 0.0_dp
    H_times_soln_vec                = 0.0_dp
    elem_vec_A_times_q              = 0.0_dp

    elem_matrix_A_times_W = matmul(inverse_A_matrix, matrix_W_left_face)

    !---Calculate H matrix, will be inverted later on
    elem_matrix_H = matmul(inverse_A_matrix,elem_matrix_U) - &
                    lamda_i_mat(isotope,delay_group)*identity_matrix - &
                    matmul(inverse_A_matrix,matrix_W_right_face)
    
    elem_vec_A_times_q = matmul(inverse_A_matrix,elem_vec_q)

    !---Multiply H matrix by previous soln vec
    if ( n > 1) then
        do i = 1, nodes_per_elem
            do j = 1, nodes_per_elem
                H_times_soln_vec(i) = H_times_soln_vec(i) + &
                        elem_matrix_H(i,j)*&
                        precursor_soln_prev(isotope,delay_group,n,j)
                
                A_times_W_times_upwind_elem_vec(i) = &
                        A_times_W_times_upwind_elem_vec(i) + &
                        elem_matrix_A_times_W(i,j)*&
                        precursor_soln_new(isotope,delay_group,n-1,3)
        
            end do
        end do
    !---Last element -- connect last to first in next loops
    else
        do i = 1, nodes_per_elem
            do j = 1, nodes_per_elem
                H_times_soln_vec(i) = H_times_soln_vec(i) + &
                        elem_matrix_H(i,j)*&
                        precursor_soln_prev(isotope,delay_group,n,j)
                    
                A_times_W_times_upwind_elem_vec(i) = &
                        A_times_W_times_upwind_elem_vec(i) + &
                        elem_matrix_A_times_W(i,j)*&
                        precursor_soln_new(isotope,delay_group, num_elem,3)
            end do
        end do
    end if

    if (DEBUG .eqv. .TRUE.) then
        
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
        write(outfile_unit,fmt='(a,1I3)'),' [A-1]*[W_l] | element --> ', n
        do i=1,nodes_per_elem
              write(outfile_unit,fmt='(12es14.3)') &
                   (elem_matrix_A_times_W(i,j),j=1,nodes_per_elem)
        end do

        write(outfile_unit,fmt='(a)'), ' '
        write(outfile_unit,fmt='(a,1I3)'),' [A^-1][W_l]*{c_e-1} | element --> ', n
        do j=1,nodes_per_elem
              write(outfile_unit,fmt='(12es14.3)')A_times_W_times_upwind_elem_vec(j)
        end do

        write(outfile_unit,fmt='(a)'), ' '
        write(outfile_unit,fmt='(a,1I2)'),'-lambda*I Matrix | element --> ',n
        do j=1,nodes_per_elem
              write(outfile_unit,fmt='(12es14.3)') &
                   (-lamda_i_mat(1,1)*identity_matrix(j,i),i=1,nodes_per_elem)
        end do

        write(outfile_unit,fmt='(a)'),' '
        write(outfile_unit,fmt='(a,1I2)'),'{q} element source vector | element --> ',n
        write(outfile_unit,fmt='(12es14.3)') ((beta_i_mat(isotope,delay_group)/gen_time)*&
        elem_vec_A_times_q(i),i=1,nodes_per_elem)
    end if

end 

