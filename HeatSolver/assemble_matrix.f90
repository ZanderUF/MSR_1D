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
    integer :: i, j, ii, jj   


!   Set zero for all matrix entries 
    global_matrix_K(:,:) = 0.0

    do i = 1, nodes_per_elem
        ii = (2*n - 1) + (i - 1)
        do j = 1, nodes_per_elem
            jj = (2*n - 1) + (j - 1) 
            global_matrix_K(ii,jj) = global_matrix_K(ii,jj) + heat_elem_matrix_M(i,j)
        end do 

    end do

!---No need for elemental matrices after they have been placed in the global one
    

end 
