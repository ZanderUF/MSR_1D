! Solves for the new solution vector
!
! 
! Input:
! 
! Output:
!
! 
subroutine solve_soln_transient(n, nl_iter )

    USE parameters_fe

    implicit none

!---Dummy
    integer,intent(in) :: nl_iter
    integer,intent(in) :: n

!---Local
    integer :: i,j
    real, dimension(3)   :: rhs_final_vec, temp_H_vec 
    real, dimension(3,3) :: inverse_matrix

!---PRECURSOR SOLVE

    temp_H_vec = matmul(elem_matrix_H,precursor_soln_prev(n,:) )

    do i = 1, nodes_per_elem
        precursor_soln_new(n,i) = precursor_soln_prev(n,i) + &
        delta_t*(temp_H_vec(i) + elem_vec_q(i) + A_times_W_times_upwind_elem_vec(i) )
    end do
!---END PRECURSOR SOLVE    

    write(outfile_unit,fmt='(a)'), ' ' 
    write(outfile_unit,fmt='(a,1I3)'),'Solution | element --> ', n
    do j=1,nodes_per_elem 
          write(outfile_unit,fmt='(a,1I2,12es14.3)'), 'Node -->', n-1+j, precursor_soln_new(n,j)          
    end do   
    write(outfile_unit,fmt='(a)'), '********************************'

end
