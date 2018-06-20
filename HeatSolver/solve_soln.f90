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
    integer :: j, i, matrix_length, lda, info, lwork
    integer, dimension(2*num_elem) :: ipiv
    real(kind=8), dimension(2*num_elem) :: work
    ! matrix for testing inversion routines
    real(kind=8), dimension(4,4) :: test_global_matrix_K 
    data test_global_matrix_K   / 1, 2, 0, 0,& 
                                  0, 1, 0, 0,&
                                  0, 0, 1, 0,&
                                  0, 0, 0, 1/
    real(kind=8), dimension(4,4) :: invert_test_matrix_K
    data invert_test_matrix_K  / 1, -2, 0, 0,&
                                 0,  1, 0, 0,&
                                 0,  0, 1, 0,&
                                 0,  0, 0, 1 /
    work = 0.0
    matrix_length = 2*num_elem
    lda = matrix_length
    lwork = matrix_length

    allocate(inverse_matrix_K(matrix_length,matrix_length) ) 
    inverse_matrix_K = 0.0

    unit_test = .FALSE.
!   Unit test for matrix inversion routines
    if(unit_test .eqv. .TRUE. ) then
        num_elem=2
        matrix_length = 4
        lda = matrix_length
        lwork = matrix_length
        final_global_matrix_K(:,:) = 0.0
        final_global_matrix_K= test_global_matrix_K(:,:) 

    end if

!   Check if we are doing steady state solve
    if (steady_state_flag .eqv. .TRUE. ) then
        ! only need to invert K once for steady state
        if( nl_iter .eq. 1) then
            ! Factorizes matrix
            call dgetrf ( matrix_length, matrix_length, global_matrix_K, lda, ipiv, info )
              if (DEBUG .eqv. .TRUE.) then
                  write(outfile_unit,fmt='(a)'), ' ' 
                  write(outfile_unit,fmt='(a)'),'Factorization of Matrix K-steady state after B.C. applied: '
                  do j=1,2*num_elem
                         write(outfile_unit,fmt='(12es14.6)') &
                              (global_matrix_K(j,i) ,i=1,2*num_elem)             
                  end do
              end if

        !  Compute the inverse matrix.
        call dgetri (matrix_length, global_matrix_K, lda, ipiv, work, lwork, info ) 
        
        end if 
        ! otherwise we don't need the inversion
        
    end if

    if (DEBUG .eqv. .TRUE.) then
        
        if(unit_test .eqv. .TRUE.) then
            do j=1,2*num_elem
                do i=1, num_elem
                    test_diff = global_matrix_K(i,j) - invert_test_matrix_K(i,j) 
                end do
            end do
            if (test_diff > 0.0) then
                write(outfile_unit, fmt='(a)') 'The matrix inverison routine is failing'
            end if
            write(outfile_unit,fmt='(a)'), ' ' 
            write(outfile_unit,fmt='(a)'),'Test Inversion of Matrix K '
            do j=1,2*num_elem
               write(outfile_unit,fmt='(12es14.6)') &
                    (global_matrix_K(j,i) ,i=1,2*num_elem)             
            end do
            write(outfile_unit,fmt='(a)'),'Known Soln of Inversion of Matrix K '
            do j=1,2*num_elem
               write(outfile_unit,fmt='(12es14.6)') &
                    (invert_test_matrix_K(j,i) ,i=1,2*num_elem)             
            end do
        
        else
            write(outfile_unit,fmt='(a)'), ' ' 
            write(outfile_unit,fmt='(a)'),'Inversion of Matrix K- steady state after B.C. applied: '
            do j=1,2*num_elem
               write(outfile_unit,fmt='(12es14.6)') &
                    (global_matrix_K(j,i) ,i=1,2*num_elem)             
        end do
        
        end if

    end if

end 
