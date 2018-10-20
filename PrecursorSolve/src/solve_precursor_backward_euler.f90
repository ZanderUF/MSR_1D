! Solves for the new precursor solution vector
! with backward Euler time integration 
! Input:
! 
! Output:
! 

subroutine solve_precursor_backward_euler(isotope,delay_group,n, nl_iter )

    USE parameters_fe

    implicit none

!---Dummy
    integer,intent(in) :: isotope
    integer,intent(in) :: delay_group
    integer,intent(in) :: nl_iter
    integer,intent(in) :: n

!---Local
    integer :: i,j,f,g
    real, dimension(3) :: A_inv_times_RHS
  
    A_inv_times_RHS(:) = 0.0

!---Multiply A^{-1}*f(u)
    do i = 1, nodes_per_elem
        do j =1, nodes_per_elem
            A_inv_times_RHS(i) = A_inv_times_RHS(i) + &
                                 inverse_A_matrix(i,j)*RHS_transient_final_vec(j)
        end do
    end do 
    
!---Solve for new precursor at delta t
    if(td_method_type == 0) then
        do i = 1, nodes_per_elem
        precursor_soln_new(isotope,delay_group, n,i) = &
                precursor_soln_prev(isotope, delay_group, n,i) + &
                delta_t*(A_inv_times_RHS(i) )
        !---test to make sure values are not too small
        !if(precursor_soln_new(isotope, delay_group, n, i) < 1E-16_dp) then
        !    precursor_soln_new(isotope, delay_group, n, i) = 0.0
        !end if

        end do
    end if
    
    if(td_method_type == 1) then
        do i = 1, nodes_per_elem
        precursor_soln_new(isotope,delay_group, n,i) = &
                precursor_soln_last_time(isotope, delay_group, n,i) + &
                delta_t*(A_inv_times_RHS(i) )
        !---test to make sure values are not too small
        if(precursor_soln_new(isotope, delay_group, n, i) < 1E-16_dp) then
            precursor_soln_new(isotope, delay_group, n, i) = 0.0
        end if

        end do
    end if
!---End precursor solve    

!---Null transient unit test
    if( reactivity_input == 0.0) then
        !---Make sure change in the solution is zero for null transient
        do i = 1, nodes_per_elem
            if(A_inv_times_RHS(i) > 0.0) then
                !write(outfile_unit, fmt='(a)'), 'FAILING NULL TRANSIENT TEST'
            end if
        end do
    end if

!*************************************
    if (DEBUG .eqv. .TRUE.) then
    !---Write out solution for current element 
        write(outfile_unit,fmt='(a,1I3,a,1I3)'),'Previous Solution | element --> ',&
        n, ' Group -->',delay_group
        do j=1,nodes_per_elem 
              write(outfile_unit,fmt='(a,1I2,12es14.3)'), 'Node -->', &
              conn_matrix(n,j), precursor_soln_prev(isotope, delay_group,n,j)          
        end do
        write(outfile_unit,fmt='(a)'), ' ' 
        write(outfile_unit,fmt='(a,1I3,a,1I3)'),'New Solution | element --> ', n, ' Group -->',delay_group
        do j=1,nodes_per_elem 
              write(outfile_unit,fmt='(a,1I2,12es14.3)'), 'Node -->', &
                conn_matrix(n,j), precursor_soln_new(isotope, delay_group,n,j)          
        end do   
        write(outfile_unit,fmt='(a)'), '********************************'
    end if

end subroutine solve_precursor_backward_euler 
