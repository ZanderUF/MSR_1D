!*****************************************************************************80
!
!! Assemble the global matrix for the problem 
!  Discussion:
!             As the elemental matrices are created we add them to the global 
!
!  Input: n - element number
!
!  Output:
! 

subroutine assemble_matrix (n)
!
    USE parameters_fe  

    implicit none
!---Dummy variable
    integer :: n
!---local
    integer :: i, j, ii, jj, nr,nc, ncl   
    do i = 1, nodes_per_elem
        nr = conn_matrix(n,i)
        do j = 1, nodes_per_elem
            nc = conn_matrix(n,j) 
            global_matrix_A(nr,nc) = global_matrix_A(nr,nc) + elem_matrix_A(i,j) 
        end do
    end do

!   Assemble global vector sources + B.C.
    do i = 1, nodes_per_elem
        ii = (2*n - 1) + (i - 1)
        global_vec_q(ii) = global_vec_q(ii) + elem_vec_f(i)
    end do
    
!   Assemble response matrix for partial current evaluation
    

!!   Assemble global matrices
!    do i = 1, nodes_per_elem
!        ii = (2*n - 1) + (i - 1)
!        do j = 1, nodes_per_elem
!            jj = (2*n - 1) + (j - 1) 
!            global_matrix_A(ii,jj) = global_matrix_A(ii,jj) + elem_matrix_K(i,j)
!        end do 
!
!    end do
!    ii = 0

end 
