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
    real :: ii, elem_length, density
    real :: center_temp_initial, center_power_initial,total_power_initial
    integer :: dist_num
    real :: norm_sin,norm_cos, pi, cosine_term, x_last, x_curr
    parameter (pi = 3.1415926535897932)

!---Simplified UNIT TEST 
!---Have sinusodial initial condition
    if(unit_test .eqv. .TRUE.) then
        nl_iter = 1
        max_nl_iter = 1 ! this unit test is linear 
        do i = 1, num_elem 
            do j = 1, nodes_per_elem
                norm_sin = real(global_coord(i,j))/real(global_coord(i,j))
                previous_elem_soln_vec(i,j) = sin(2.0*pi*norm_sin)
            end do
        end do
    end if

!---Normal calculation flow - no need if doing unit test
    if (unit_test .eqv. .FALSE.) then
        previous_elem_soln_vec(:,:) = 0
        cur_elem_soln_vec(:,:)= 0
        nl_iter = 1 
        steady_state_flag = .TRUE.
        nonlinear_ss_flag = .TRUE.
        !---Sets the max number of nonlinear iterations    
        max_nl_iter = 100 
         
        allocate( power_initial(num_elem,nodes_per_elem) )
        !---Initial guesses 
        center_temp_initial  = 800
        total_power_initial  = 10
        center_power_initial = 10
        print *,'global_coord las',global_coord(non_fuel_start,3)
        !---Used to apply cosine shape over active domain 
        dist_num = ((non_fuel_start) + 1 )/2
        !---Apply to every node point within an element
        do i = 1,num_elem
            do j = 1, nodes_per_elem
                !---Apply to active fuel region
                if( i <= non_fuel_start) then
                    x_curr = real(global_coord(i,j) )
                    x_last =  real(global_coord(non_fuel_start,3))
                    norm_cos = (x_curr - x_last/2)/(x_last) !+ 
                    cosine_term = cos( (pi/2)*norm_cos )
                    power_initial(i,j) = (center_power_initial*cosine_term)
                    !---Set temperature distribution
                    temperature_vec(i,j) = (center_temp_initial*cosine_term)
                    !---Get density to set the velocity
                    call density_corr(temperature_vec(i,j),density)
                    density_vec(i,j) = density
                    !---Need to get initial velocity distribution
                    velocity_vec(i,j) = mass_flow/(area*density)
                    !velocity_vec(i,j) = 100 
                !---Inactive region assumed to have zero power 
                else
                    !---Temperature in inactive region same as end of active region ==> no loss
                    temperature_vec(i,j) = temperature_vec(non_fuel_start ,3)
                    !---Get density to set the velocity
                    call density_corr(temperature_vec(i,j),density)
                    density_vec = density
                    !---Need to get initial velocity distribution
                    velocity_vec(i,j) = mass_flow/(area*density)
                    power_initial(i,j) = 0.0
                end if
           end do 
        end do
        !-------------------------------------------------------------------------------
        !---Write out initial solution
        write(outfile_unit,fmt='(a)'), ' '
        write(outfile_unit,fmt='(a)'), 'Initial Power distribution '
        write(outfile_unit,fmt='(a)'), 'Position(x) Power [W]'
        do i = 1,num_elem
            do j = 1, nodes_per_elem
                write(outfile_unit, fmt='(f6.3, 12es14.3)')  global_coord(i,j), power_initial(i,j)
            end do
        end do
        
        !-------------------------------------------------------------------------------
        !---Write out initial solution
        write(outfile_unit,fmt='(a)'), ' '
        write(outfile_unit,fmt='(a)'), 'Initial temperature distribution '
        write(outfile_unit,fmt='(a)'), 'Position(x) Temperature [K]'
        do i = 1, num_elem 
            do j = 1, nodes_per_elem
                write(outfile_unit, fmt='(f6.3, 12es14.3)')  global_coord(i,j), temperature_vec(i,j)
            end do
        end do
        !---Write out initial solution
        write(outfile_unit,fmt='(a)'), ' '
        write(outfile_unit,fmt='(a)'), 'Initial velocity distribution '
        write(outfile_unit,fmt='(a)'), 'Position(x) Velocity [m/s]'
        do i = 1, num_elem 
            do j = 1, nodes_per_elem
                write(outfile_unit, fmt='(f6.3, 12es14.3)')  global_coord(i,j), velocity_vec(i,j)
            end do 
        end do
        !---Write out initial solution
        write(outfile_unit,fmt='(a)'), ' '
        write(outfile_unit,fmt='(a)'), 'Initial density distribition '
        write(outfile_unit,fmt='(a)'), 'Position(x) Density [kg/m^3]'
        do i = 1, num_elem 
            do j = 1, nodes_per_elem
                write(outfile_unit, fmt='(f6.3, 12es14.3)')  global_coord(i,j), density_vec(i,j)
            end do
        end do
        !-------------------------------------------------------------------------------
        !---Write out initial solution
        write(outfile_unit,fmt='(a)'), ' '
        write(outfile_unit,fmt='(a)'), 'Initial condition '
        write(outfile_unit,fmt='(a)'), 'Position(x) Precursor Concentration'
        do i = 1, num_elem 
            do j = 1, nodes_per_elem
                write(outfile_unit, fmt='(f6.3, 12es14.3)')  global_coord(i,j), cur_elem_soln_vec(i,j)
            end do
        end do
        
        !---Steady state solver for nonlinear ss problems
        if(nonlinear_ss_flag .eqv. .TRUE.) then 
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
                !---Calculate residual
                ! If residual < tolerance exit loop
                !if (residual < tolerance) then
                !    exit
                !end if
                
                !-------------------------------------------------------------------------------
                !---Write out initial solution
                write(outfile_unit,fmt='(a)'), ' '
                write(outfile_unit,fmt='(a,1I4)'), 'Steady state soln at nl iter', nl_iter
                write(outfile_unit,fmt='(a)'), 'Position(x) Precursor Concentration'
                do i = 1, num_elem 
                    do j = 1, nodes_per_elem
                        write(outfile_unit, fmt='(f6.3, 12es14.3)')  global_coord(i,j), cur_elem_soln_vec(i,j)
                    end do
                end do

                nl_iter = nl_iter + 1
                ! make previous = current solution vector
                !previous_elem_soln_vec = cur_elem_soln_vec
    
                ! If we've gone thru too many nonlinear iterations exit
                if (nl_iter > max_nl_iter) then
                    exit
                end if
            
            end do !---end nonlinear iteration loop
        end if !---end nonlinear if 
    end if !---end normal calculation if 
    
    !---write out converged solution for plotting
    write(soln_outfile_unit,fmt='(a)'), 'Position(x) | Precursor Concentration'
    do i = 1, num_elem 
        do j = 1, nodes_per_elem
            write(soln_outfile_unit, fmt='(f6.3, 12es14.3)')  global_coord(i,j), cur_elem_soln_vec(i,j)
        end do
    end do
!---Set steady state flag off
    steady_state_flag = .FALSE.

end
