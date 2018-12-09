!*****************************************************************************80
!
!  Calculates the numerical flux matrices 
!
!  Discussion:
!             Generates the elemental matrices K_ij, M_ij, S_ij 
!             calculates coefficients in this suboutine 
!  Input: 
!       n - element number
!       nl_iter - nonlinear iteration counter
!  Output:
!       matrix_W_left
!       matrix_W_right
!       elem_vec_q
!*****************************************************************************80

subroutine numerical_flux_matrices (n, nl_iter)
 
    USE element_matrices_M
    USE flags_M
    USE global_parameters_M
    USE mesh_info_M
    USE solution_vectors_M

    implicit none

!---Dummy variables
    integer, intent(in) :: n
    integer, intent(in) :: nl_iter

!---Local variables
    integer :: i, j

!---Initialize
    matrix_W_right_face = 0.0_dp
    matrix_W_left_face  = 0.0_dp
    !elem_vec_q          = 0.0_dp

!---Create source vector 'q', and W - 
    do i = 1, nodes_per_elem
        do j = 1, nodes_per_elem
            !elem_vec_q(i) = elem_vec_q(i) + &
            !                spatial_power_fcn(n,j)*&
            !                power_amplitude_prev!*&
            !                elem_matrix_A(i,j)

                            ! try just working with the fractional power
                            ! not the 'power' in Watts itself

            !---Applies for all elements except the first one
            if(n > 1) then !--- n - element #
                !---Grab previous precursor conc. + velocity at 
                !---rhs of previous element
                 matrix_W_right_face(i,j) = velocity_soln_prev(n,i)*&
                                            interp_fcn_rhs(i)*interp_fcn_rhs(j)  
                 matrix_W_left_face(i,j)  = velocity_soln_prev(n-1,i)*&
                                            interp_fcn_lhs(i)*interp_fcn_lhs(j)
            else!---First element case, need to connect with end element 
                 matrix_W_right_face(i,j) = velocity_soln_prev(n,i)*&
                                            interp_fcn_rhs(i)*interp_fcn_rhs(j)  
                 matrix_W_left_face(i,j)  = velocity_soln_prev(num_elem,i)*&
                                            interp_fcn_lhs(i)*interp_fcn_lhs(j)
            end if!---End flux calculation
        end do
    end do

!--------------------------------------------------------------------
    if( DEBUG .eqv. .TRUE.) then
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
        
    end if

end subroutine numerical_flux_matrices
