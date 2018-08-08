!*****************************************************************************80
!
!! Assemble the elemental to solve for coefficients the problem 
!  Discussion:
!             As the elemental matrices are created we add them to the global 
!
!  Input: n - element number
!
!  Output:
! 

subroutine assemble_matrix (n)
!---
    USE parameters_fe  

    implicit none
!---Dummy variables
    integer :: n
!---Local variables
    integer :: i, j, ii, jj, nr,nc, ncl,length   
    real, dimension(3) ::  elem_prev_soln, flux_rhs, flux_lhs, &
                          temp_vec, basis_at_lhs, basis_at_rhs, rhs_final_vec
                          
    real, dimension(3,3) :: inverse_matrix,temp_matrix
    real :: temp_soln
    
!---Inversion routine parameters
    integer :: lda, info, lwork
    integer, dimension(3) :: ipiv
    real, dimension(3) :: work

    length = 3
!---Initialize
    ipiv =  0
    work =  0.0
    lda =   length
    lwork = length

    elem_matrix_H = 0
    elem_matrix_G = 0
    elem_vec_w_left_face = 0
    rhs_final_vec = 0
    
    elem_vec_Pu = 0

    do i = 1, nodes_per_elem
        !---Calculate q vector
        !elem_vec_q(i) = (beta/gen_time)*elem_vec_q(i)
        do j = 1, nodes_per_elem
            if (n > 1) then
                elem_vec_w_left_face(i) = elem_vec_w_left_face(i) + &
                                          matrix_W_left_face(i,j)*precursor_soln_new(n-1,3)
            else
                elem_vec_w_left_face(i) = elem_vec_w_left_face(i) + &
                                          matrix_W_left_face(i,j)*precursor_soln_new(num_elem,3)
            end if
            
            elem_matrix_G(i,j) = -elem_matrix_U(i,j) + &
                                 (log(2.0)/lambda)*elem_matrix_A(i,j) + &
                                 matrix_W_right_face(i,j)
        end do
    end do
    
    write(outfile_unit,fmt='(a)'), ' '
    write(outfile_unit,fmt='(a,1I2)'),'G Matrix | element --> ',n
    do i=1,nodes_per_elem 
          write(outfile_unit,fmt='(12es14.3)') &
               (elem_matrix_G(i,j),j=1,nodes_per_elem)             
    end do

!****************************************************************
!---Write out    
    write(outfile_unit,fmt='(a)'), ' ' 
    write(outfile_unit,fmt='(a,1I3)'),'left face vector | element --> ', n
    do j=1,nodes_per_elem 
          write(outfile_unit,fmt='(12es14.3)') elem_vec_w_left_face(j)             
    end do
    
    write(outfile_unit,fmt='(a)'),' '
    write(outfile_unit,fmt='(a,1I2)'),'{q} element source vector | element --> ',n
    write(outfile_unit,fmt='(12es14.3)') (elem_vec_q(i),i=1,nodes_per_elem)             

end 
