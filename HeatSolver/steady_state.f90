! Steady state solve of the heat equation 
!
! Solves d^2/dx^2 (T) = P/(kappa*V) , P - power, V - volume
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
    integer :: i,ii,jj,j,n, nl_iter, max_nl_iter
    real :: elem_length
    real :: center_temp_initial, center_power_initial
    integer :: dist_num
    real :: pi
    parameter (pi = 3.1415926535897932)

    !
    allocate( power_initial(num_elem*nodes_per_elem) )
    !   Initial conditions 
    center_temp_initial = 900
    center_power_initial = 10 
    max_nl_iter = 10

!   Apply initial guess to solution vector
    dist_num = (num_elem*nodes_per_elem)/2 +1
    do i = -dist_num , dist_num 
        ii = (real(i)/real(dist_num))
        previous_elem_soln_vec(i+dist_num) = center_temp_initial*cos( ii*(-pi/2.0) )
        ! Test as we are not solving the power equation yet
        elem_length = elem_lengths(i + dist_num + 1) - elem_lengths(i + dist_num)
        ! Power/Volume
        power_initial(i+dist_num) = (center_power_initial*cos(ii*(pi/2.0) ) )/(A*elem_length)
    end do

    nl_iter = 0
    steady_state_flag = .TRUE.

!  Nonlinear loop
    do 
        
        do n = 1, num_elem
            ! Computer K_ij F_ij
            call element_matrix_heat(n,nl_iter) 
            ! Assemble K, F
            call assemble_matrix(n)
        end do 

        ! Solve T^r = [K(T^(r-1)]^-1 F^(r-1) 

        ! Calculate residual

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
    
    end do 
    
    steady_state_flag = .FALSE.

end
