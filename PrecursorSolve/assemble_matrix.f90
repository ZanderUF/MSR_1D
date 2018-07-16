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
    
!---Assemble into global matrix 
    do i = 1, nodes_per_elem
        nr = conn_matrix(n,i)
        do j = 1, nodes_per_elem
            nc = conn_matrix(n,j) 
            global_matrix_A(nr,nc) = global_matrix_A(nr,nc) + elem_matrix_A(i,j) 
        end do
    end do

    if(n > 1) then 
        nr = conn_matrix(n,1) - 1
        nc = conn_matrix(n,1)
        global_matrix_A(nr,nc) = global_matrix_A(nr,nc) + -1.0
    end if

!---Assemble global vector source
    if(unit_test .eqv. .FALSE.) then 
        !---Assemble global vector sources + B.C.
        do i = 1, nodes_per_elem
            ii = (2*n - 1) + (i - 1)
            global_vec_q(ii) = global_vec_q(ii) + elem_vec_f(i)
        end do
    end if 

end 
