! Steady state solve of the equation 
!
! Input:
!
! Output:
! 
subroutine steady_state()

    USE parameters_fe

implicit none

!   Dummy

!   Local
    
    integer :: i,jj,j,n, nl_iter, max_nl_iter
    real :: ii, elem_length
    real :: center_temp_initial, center_power_initial
    integer :: dist_num
    real :: norm_sin, pi, cosine_term
    parameter (pi = 3.1415926535897932)

    !---Simplified unit test
    !---Have sinusodial initial condition
    if(unit_test .eqv. .TRUE.) then
        nl_iter = 1
        max_nl_iter = 1 ! this unit test is linear 
        do i = 1, matrix_length
            norm_sin = real(global_coord(i))/real(global_coord(matrix_length))
            previous_elem_soln_vec(i) = sin(2.0*pi*norm_sin)
            !previous_elem_soln_vec(i) = 1.0
        end do
    end if
    
!-------------------------------------------------------------------------------
!---Write out temperature solution
        write(outfile_unit,fmt='(a)'), ' '
        write(outfile_unit,fmt='(a)'), 'Initial condition '
        write(outfile_unit,fmt='(a)'), 'Position(x) Solution'
        do j = 1, matrix_length
            write(outfile_unit, fmt='(f6.3, f10.3)')  global_coord(j), previous_elem_soln_vec(j)
        end do

!---Normal calculation flow - no need if doing unit test
    if (unit_test .eqv. .FALSE.) then
        allocate( power_initial(matrix_length) )
        allocate( elem_node_lengths(matrix_length) )
        !---Initial guesses 
        center_temp_initial = 900
        center_power_initial = 100 
        !---Sets the max number of nonlinear iterations    
        max_nl_iter = 100 
        !---Apply to every node point within an element
        dist_num = ((matrix_length-3) + 1)/2 + 1
        do i = 1, matrix_length
            ii = (real(i-dist_num)/real(dist_num))
            cosine_term = cos(ii*(pi/2.0))
            previous_elem_soln_vec(i) = center_temp_initial*cosine_term
        !---Power
            power_initial(i) = (center_power_initial*cosine_term)

            if(cosine_term < 0.0) then
                cosine_term = 0.0
                previous_elem_soln_vec(i) = 0.0
                power_initial(i) = 0.0
            end if

        end do
        !---Set initial soln guess
        previous_elem_soln_vec(:) = 800
        nl_iter = 1 
        steady_state_flag = .TRUE.
        
        write(outfile_unit, fmt='(a)'), ' ' 
        write(outfile_unit, fmt='(a)'), 'Start steady state calculation'
    
    !---Nonlinear loop
        do 
            
            do n = 1, num_elem
                !---Computer K_ij F_ij
                call element_matrix(n,nl_iter) 
                !---Assemble K, F
                call assemble_matrix(n)
            end do 
    
            !---Apply boundary conditions
            call boundary_cond 
    
            !---Solve T^r = [K(T^(r-1)]^-1 F^(r-1) 
            call solve_soln(nl_iter)
            
            !---Calculate residual
    
            ! If residual < tolerance exit loop
            !if (residual < tolerance) then
            !    exit
            !end if
    
            nl_iter = nl_iter + 1
            ! make previous = current solution vector
            previous_elem_soln_vec = cur_elem_soln_vec
    
            ! If we've gone thru too many nonlinear iterations exit
            if (nl_iter > max_nl_iter) then
                exit
            end if
        
        end do ! end nonlinear iteration loop
    
    end if !---end normal calculation if 
     
!---Set steady state flag off
    steady_state_flag = .FALSE.

   


end
