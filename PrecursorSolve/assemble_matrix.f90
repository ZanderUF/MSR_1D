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
!
    USE parameters_fe  

    implicit none
!---Dummy variables
    integer :: n
!---Local variables
    integer :: i, j, ii, jj, nr,nc, ncl,length   
    real, dimension(3) :: elem_prev_soln, flux_rhs, flux_lhs, basis_at_lhs
    data basis_at_lhs / 1, 0, 0/
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

    elem_prev_soln = 0
    flux_rhs = 0
    flux_lhs = 0
    Pu_minus_flux_vec = 0

!---Solve for elemental coeficients, no global assembly 
    do i = 1, nodes_per_elem
        do j = 1, nodes_per_elem
            ii = j + (n-1)*nodes_per_elem
            elem_vec_Pu(i) = elem_vec_Pu(i) + &
                             elem_matrix_P(i,j)*previous_elem_soln_vec(ii)
        end do
    end do

!---All elements except the first
    if(n > 1) then
        !---Grab current solution at lhs of current element 
        ii = 1 + (n-1)*nodes_per_elem
        temp_soln = previous_elem_soln_vec(ii)
        do i = 1, nodes_per_elem
            flux_rhs(i) = temp_soln*basis_at_lhs(i)
        end do
        !---Grab previous solution at rhs of previous element
        temp_soln = previous_elem_soln_vec(ii-1)
        do i = 1, nodes_per_elem
            flux_lhs(i) = temp_soln*basis_at_lhs(i)
        end do
    else!---For periodic B.C. need to connect first to the last element
        
        !---Get current solution at lhs of current element
        ii = 1 + (n-1)*nodes_per_elem
        temp_soln = previous_elem_soln_vec(ii)
        do i = 1, nodes_per_elem
            flux_rhs(i) = temp_soln*basis_at_lhs(i)
        end do
        !---Connect the last element to first
        temp_soln = previous_elem_soln_vec(matrix_length)
        do i = 1, nodes_per_elem
            flux_lhs(i) = temp_soln*basis_at_lhs(i)
        end do

    end if
   
!---Combine (Pu - f)
    Pu_minus_flux_vec = elem_vec_Pu + flux_rhs - flux_lhs
    
!---Write out {Pu - f}
    write(outfile_unit,fmt='(a)'), ' ' 
    write(outfile_unit,fmt='(a,1I3)'),'{Pu -f} vector | element --> ', n
    do j=1,nodes_per_elem 
           write(outfile_unit,fmt='(12es14.3)') Pu_minus_flux_vec(j)             
    end do

!---Factorize A matrix
    call dgetrf ( length, length, elem_matrix_A, lda, ipiv, info )

!---Write out factorized A matrix
    write(outfile_unit,fmt='(a)'), ' ' 
    write(outfile_unit,fmt='(a,1I3)'),'LU Factorized [A] matrix | element --> ',n
    do j=1,nodes_per_elem
           write(outfile_unit,fmt='(12es14.3)') &
                ( elem_matrix_A(j,i) ,i=1,nodes_per_elem )             
    end do

!---Compute the inverse matrix.
    call dgetri ( length, elem_matrix_A, lda, ipiv, work, lwork, info ) 

!---Write out factorized A matrix
    write(outfile_unit,fmt='(a)'), ' ' 
    write(outfile_unit,fmt='(a,1I3)'),'inverse [A]^-1 matrix | element --> ',n
    do j=1,nodes_per_elem
           write(outfile_unit,fmt='(12es14.3)') &
                ( elem_matrix_A(j,i) ,i=1,nodes_per_elem )             
    end do

!---Determine vector (Pu - f)*A^-1
    elem_prev_soln = matmul(Pu_minus_flux_vec,elem_matrix_A)

!---Solve for next time step solution
    do i = 1, nodes_per_elem
        ii = i + (n-1)*nodes_per_elem  
        cur_elem_soln_vec(ii) = dt*elem_prev_soln(i) + previous_elem_soln_vec(ii) 
    end do

end 
