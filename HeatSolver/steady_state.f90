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
    
    integer :: i,jj,j,n, nl_iter, max_nl_iter
    real :: ii, elem_length
    real :: center_temp_initial, center_power_initial
    integer :: dist_num
    real :: pi, cosine_term
    parameter (pi = 3.1415926535897932)

    !
    allocate( power_initial(matrix_length) )
    allocate( elem_node_lengths(matrix_length) )

    !   Initial conditions 
    center_temp_initial = 900
    center_power_initial = 10 
    max_nl_iter = 1 
    Area = 5.0
!   Apply initial guess to solution vector
    !do i = 1, num_elem
    !    do j =1,  nodes_per_elem
    !        elem_node_lengths( (i-1)*nodes_per_elem + j) = elem_lengths(i)
    !    end do
    !end do

    !dist_num = (matrix_length + 1)/2 + 1

!   Apply to every node point within an element
    !do i = -dist_num , dist_num 
    !do i = 1, matrix_length 
    !    ii = (real(i-dist_num)/real(dist_num))
    !    cosine_term = cos(ii*(pi/2.0))
    !    previous_elem_soln_vec(i) = center_temp_initial*cosine_term
    !    ! Power/Volume
    !    power_initial(i) = (center_power_initial*cosine_term )/&
    !    (area*(elem_lengths( 1 )))

    !    if(cosine_term < 0.0) then
    !        cosine_term = 0.0
    !        previous_elem_soln_vec(i) = 0.0
    !        power_initial(i) = 0.0
    !    end if

    !end do

!   Test problem - set initial condition
    previous_elem_soln_vec(:) = 500
    nl_iter = 1 
    steady_state_flag = .TRUE.
    
    write(outfile_unit, fmt='(a)'), ' ' 
    write(outfile_unit, fmt='(a)'), 'Start steady state calculation'

!  Nonlinear loop
    do 
        
        do n = 1, num_elem
            ! Computer K_ij F_ij
            call element_matrix_heat(n,nl_iter) 
            ! Assemble K, F
            call assemble_matrix(n)
        end do 

        ! Apply boundary conditions
        call boundary_cond 

        ! Solve T^r = [K(T^(r-1)]^-1 F^(r-1) 
         !call solve_soln(nl_iter)
        
        ! Calculate residual

        ! If residual < tolerance exit loop
        !if (residual < tolerance) then
        !    exit
        !end if

        nl_iter = nl_iter + 1
        ! make previous = current solution vector
        !previous_elem_soln_vec = cur_elem_soln_vec

        ! If we've gone thru too many nonlinear iterations exit
        if (nl_iter > max_nl_iter) then
            exit
        end if
    
    end do 
 ! Write out temperature solution
        !write(outfile_unit,fmt='(a)'), ' '
        !write(outfile_unit,fmt='(a)'), 'Temperature solution as a function of position'
        !do j = 1, matrix_length
        !    write(outfile_unit, fmt='(f6.3, f10.3)')  global_coord(j), previous_elem_soln_vec(j)
        !end do   

!---No need for elemental matrices after they have been placed in the global one
   !if (DEBUG .eqv. .TRUE. ) then
   !     write(outfile_unit,fmt='(a)'), ' ' 
   !     write(outfile_unit,fmt='(a)'),'Global Matrix - steady state: '
   !     do j=1,2*num_elem+1 
   !            write(outfile_unit,fmt='(12es14.6)') &
   !                 (global_matrix_K(j,i) ,i=1,2*num_elem+1)             
   !     end do

   !     write(outfile_unit,fmt='(a)'), ' '        
   !     write(outfile_unit,fmt='(a)'), 'Global vector source f - steady state:  '
   !     write(outfile_unit,fmt='(12es14.6)') &
   !                 (global_vec_f(i) ,i=1,2*num_elem+1)
   !end if

    steady_state_flag = .FALSE.

end
