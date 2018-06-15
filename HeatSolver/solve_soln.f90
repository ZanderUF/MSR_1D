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
    integer :: j, i, matrix_length, lda, info, lwork
    integer, dimension(2*num_elem) :: ipiv
    real (kind=dp), allocatable :: work(:) 
    allocate(work(2*num_elem))
    work = 0.0
    matrix_length = 2*num_elem
    lda = matrix_length
    lwork = matrix_length

    allocate(inverse_matrix_K(matrix_length,matrix_length) ) 
    inverse_matrix_K = 0.0
    

!   Check if we are doing steady state solve
    if (steady_state_flag .eqv. .TRUE. ) then
        ! only need to invert K once for steady state
        if( nl_iter .eq. 1) then
        !
            call dgetrf ( matrix_length, matrix_length, final_global_matrix_K, lda, ipiv, info )
            
            !if ( info /= 0 ) then
            !    write ( *, '(a)' ) ' '
            !    write ( *, '(a,i8)' ) '  DGETRF returned INFO = ', info
            !    write ( *, '(a)' ) '  The matrix is numerically singular.'
            !    return
            !end if
        print *,'lwork',lwork 
        !  Compute the inverse matrix.
        !
        call dgetri (matrix_length, final_global_matrix_K, lda, ipiv, work, lwork, info ) 
        
            print *,'final',final_global_matrix_K
        end if 
        ! otherwise we don't need the inversion
        
    end if

    if (DEBUG .eqv. .TRUE.) then
        write(outfile_unit,fmt='(a)'), ' ' 
        write(outfile_unit,fmt='(a)'),'Inversion of Matrix K- steady state after B.C. applied: '
        do j=1,2*num_elem
               write(outfile_unit,fmt='(12es14.6)') &
                    (final_global_matrix_K(j,i) ,i=1,2*num_elem)             
        end do
    end if

end 
