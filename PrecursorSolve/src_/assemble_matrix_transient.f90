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
    
    USE parameters_fe  

    implicit none
!---Dummy variables
    integer,intent(in) :: isotope
    integer,intent(in) :: delay_group
    integer,intent(in) :: n
!---Local variables
    integer :: i, j, ii, jj, nr,nc, ncl,length   
    real, dimension(3) ::  elem_prev_soln, flux_rhs, flux_lhs, &
                          temp_vec, basis_at_lhs, basis_at_rhs, rhs_final_vec
    real, dimension(3,3) :: inverse_matrix,temp_matrix
   
!---Initialize
    elem_matrix_H = 0
    elem_matrix_G = 0
    rhs_final_vec = 0
    A_times_W_times_upwind_elem_vec = 0
    elem_matrix_A_times_W = 0
    H_times_soln_vec = 0
    elem_vec_A_times_q = 0

    elem_matrix_A_times_W = matmul(inverse_A_matrix, matrix_W_left_face)
    
    !---Calculate H matrix, will be inverted later on
    elem_matrix_H = matmul(inverse_A_matrix,elem_matrix_U) - &
                    lamda_i_mat(isotope,delay_group)*identity_matrix - &
                    matmul(inverse_A_matrix,matrix_W_right_face)
    
    elem_vec_A_times_q = matmul(inverse_A_matrix,elem_vec_q)
    
    !---Multiply H matrix by previous soln vec
    do i = 1, nodes_per_elem
        do j = 1, nodes_per_elem
            H_times_soln_vec(i) = H_times_soln_vec(i) + &
                    elem_matrix_H(i,j)*&
                    precursor_soln_prev(isotope,delay_group,n,j) 
            if ( n > 1) then
                A_times_W_times_upwind_elem_vec(i) = &
                    A_times_W_times_upwind_elem_vec(i) + &
                    elem_matrix_A_times_W(i,j)*&
                    precursor_soln_prev(isotope,delay_group,n-1,3)
            else
                A_times_W_times_upwind_elem_vec(i) = &
                    A_times_W_times_upwind_elem_vec(i) + &
                    elem_matrix_A_times_W(i,j)*&
                    precursor_soln_prev(isotope,delay_group, num_elem,3)
            end if
        end do 
    end do
    
!****************************************************************
!---Write out    
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
        write(outfile_unit,fmt='(12es14.3)') (elem_vec_q(i),i=1,nodes_per_elem)             
    end if
end 
