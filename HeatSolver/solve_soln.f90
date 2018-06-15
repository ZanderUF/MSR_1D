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
    integer :: j, i

    allocate(inverse_matrix_K(2*num_elem,2*num_elem) ) 
    inverse_matrix_K = 0.0
    

!   Check if we are doing steady state solve
    if (steady_state_flag .eqv. .TRUE. ) then
        ! only need to invert K once for steady state
        if( nl_iter .eq. 1) then
            call inverse(2*num_elem)    
        end if 
        ! otherwise we don't need the inversion
        
    end if

    if (DEBUG .eqv. .TRUE.) then
        write(outfile_unit,fmt='(a)'), ' ' 
        write(outfile_unit,fmt='(a)'),'Inversion of Matrix K- steady state after B.C. applied: '
        do j=1,2*num_elem
               write(outfile_unit,fmt='(12es14.6)') &
                    (inverse_matrix_K(j,i) ,i=1,2*num_elem)             
        end do
    end if

end 
