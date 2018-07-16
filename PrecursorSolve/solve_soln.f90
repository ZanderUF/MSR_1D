! Solves for the current solution vector
!
! Steady state - T = [K]^-1 * F 
! 
! Input:
! 
! Output:
!
! 
subroutine solve_soln(nl_iter )

    USE parameters_fe

    implicit none

!   Dummy
    integer :: nl_iter

!   Local
    real :: test_diff
    integer :: j, i, lda, info, lwork
    integer, dimension(matrix_length) :: ipiv
    real, dimension(matrix_length) :: work
    ! matrix for testing inversion routines
    real, dimension(4,4) :: test_global_matrix_A 
    data test_global_matrix_A   / 1, 2, 0, 0,& 
                                  0, 1, 0, 0,&
                                  0, 0, 1, 0,&
                                  0, 0, 0, 1/
    real, dimension(4,4) :: invert_test_matrix_P
    data invert_test_matrix_P  / 1, -2, 0, 0,&
                                 0,  1, 0, 0,&
                                 0,  0, 1, 0,&
                                 0,  0, 0, 1 /
    real, dimension(2*num_elem+1) :: temperature
    ipiv = 0.0
    work = 0.0
    lda = matrix_length
    lwork = matrix_length
    info =0.0 

    unit_test = .TRUE.

!   Unit test for matrix inversion routines
    if(unit_test .eqv. .TRUE. ) then
        num_elem=2
        matrix_length = 4
        lda = matrix_length
        lwork = matrix_length
        global_matrix_A(:,:) = 0.0
        global_matrix_A= test_global_matrix_A(:,:) 

    end if

!   Check if we are doing steady state solve
    if (steady_state_flag .eqv. .TRUE. .and. unit_test .eqv. .FALSE.) then
        ! Factorizes matrix
            call dgetrf ( matrix_length, matrix_length, global_matrix_A, lda, ipiv, info )
            
            if (DEBUG .eqv. .TRUE.) then
                  write(outfile_unit,fmt='(a)'), ' ' 
                  write(outfile_unit,fmt='(a)'),'Factorization of Matrix A-steady state after B.C. applied: '
                  do j=1,matrix_length
                         write(outfile_unit,fmt='(12es14.3)') &
                              (global_matrix_A(j,i) ,i=1,matrix_length)             
                  end do
            end if
        
        !  Compute the inverse matrix.
        call dgetri (matrix_length, global_matrix_A, lda, ipiv, work, lwork, info ) 
        
    end if

!   Solve for current solution 
    cur_elem_soln_vec = matmul(global_vec_q,global_matrix_A)

!--------------------------------------------------------------------------------------
    if (DEBUG .eqv. .TRUE.) then
        
        if(unit_test .eqv. .TRUE.) then
            do j=1,matrix_length
                do i=1, num_elem
                    test_diff = global_matrix_A(i,j) - invert_test_matrix_P(i,j) 
                end do
            end do
            if (test_diff > 0.0) then
                write(outfile_unit, fmt='(a)') 'The matrix inverison routine is failing'
            end if
            write(outfile_unit,fmt='(a)'), ' ' 
            write(outfile_unit,fmt='(a)'),'Test Inversion of Matrix A '
            do j=1,matrix_length
               write(outfile_unit,fmt='(12es14.6)') &
                    (global_matrix_A(j,i) ,i=1,matrix_length)             
            end do
            write(outfile_unit,fmt='(a)'),'Known Soln of Inversion of Matrix A '
            do j=1,matrix_length
               write(outfile_unit,fmt='(12es14.6)') &
                    (invert_test_matrix_P(j,i) ,i=1,matrix_length)             
            end do
        
        else
            write(outfile_unit,fmt='(a)'), ' ' 
            write(outfile_unit,fmt='(a)'),'Inversion of Matrix A steady state after B.C. applied: '
            do j=1,matrix_length
               write(outfile_unit,fmt='(12es14.3)') &
                    (global_matrix_A(j,i) ,i=1,matrix_length)             
        end do
        
        end if

    end if
!--------------------------------------------------------------------------------------
    
    global_matrix_A = 0.0
    global_vec_q = 0.0

end 
