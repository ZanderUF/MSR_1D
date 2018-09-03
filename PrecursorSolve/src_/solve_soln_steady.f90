! Solves for the new solution vector
!
! Steady state - T = [K]^-1 * F 
! 
! Input:
! 
! Output:
!
! 
subroutine solve_soln_steady(isotope, delay_group, n, nl_iter )

    USE parameters_fe

    implicit none

!---Dummy
    integer,intent(in) :: isotope
    integer,intent(in) :: delay_group 
    integer,intent(in) :: nl_iter
    integer,intent(in) :: n ! current node

!---Local
    real, dimension(3)  :: rhs_final_vec 
    real, dimension(3,3) :: inverse_matrix
    !---Inversion routine parameters
    integer :: lda, info, lwork
    integer, dimension(3) :: ipiv
    real, dimension(3) :: work
    integer :: length, i, j, f,g

    length = 3
!---Initialize
    ipiv =  0
    work =  0.0
    lda =   length
    lwork = length

!---PRECURSOR SOLVE
    rhs_final_vec = 0
    !---Setup G matrix to be inverted later on
    do i = 1, nodes_per_elem
        do j = 1, nodes_per_elem
            inverse_matrix(i,j) = elem_matrix_G(i,j)
        end do 
    end do
    
    !---Compute RHS vector
    do i = 1, nodes_per_elem
        rhs_final_vec(i) = elem_vec_q_final(isotope,delay_group,i) + elem_vec_w_left_face(i)
    end do
    
    !---Write out
    write(outfile_unit,fmt='(a)'), ' ' 
    write(outfile_unit,fmt='(a,1I3)'),'RHS final vector | element --> ', n
    do j=1,nodes_per_elem 
          write(outfile_unit,fmt='(12es14.3)') rhs_final_vec(j)             
    end do
    
    !---Factorize G matrix matrix
    call dgetrf ( length, length, inverse_matrix, lda, ipiv, info )
    !---Compute the inverse matrix.
    call dgetri ( length, inverse_matrix, lda, ipiv, work, lwork, info )      
     
    !---CALCULATE SOLUTION for a given element
    precursor_soln_new(isotope,delay_group,n,:) = matmul(inverse_matrix,rhs_final_vec)
    
    !---Write out inverse matrix
    write(outfile_unit,fmt='(a)'), ' ' 
    write(outfile_unit,fmt='(a,1I3)'),'inverse G precursor matrix matrix | element --> ',n
    do j=1,nodes_per_elem
           write(outfile_unit,fmt='(12es14.3)') &
                ( inverse_matrix(j,i) ,i=1,nodes_per_elem )             
    end do
!---END PRECURSOR SOLVE   

    write(outfile_unit,fmt='(a)'), ' '
    write(outfile_unit,fmt='(a,1I3,a,1I3,a,1I3)'),&
        'Solution | element --> ', n, ' Isotope ',isotope, ' Delayed Group ',delay_group
        do j=1,nodes_per_elem 
            write(outfile_unit,fmt='(a,1I2,12es14.3)'), 'Node -->', n-1+j,&
              precursor_soln_new(isotope,delay_group,n,j)         
     
    end do
      
    write(outfile_unit,fmt='(a)'), '********************************'
   
end 
