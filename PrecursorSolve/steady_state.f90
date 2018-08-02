! Steady state solve of the equation 
!
! Input:
!
! Output:
! 
subroutine steady_state()

    USE parameters_fe

implicit none

!---Dummy

!---Local
    integer :: i,jj,j,n, nl_iter
    real :: ii, elem_length, density, &
             center_temp_initial,  &
             norm_sin,norm_cos, pi, cosine_term, x_last, x_curr,&
             temperature
    parameter (pi = 3.1415926535897932)

!---------------------------------------------------------------
!---Simplified UNIT TEST 
    if(unit_test .eqv. .TRUE.) then
    !---Have sinusodial initial condition
        nl_iter = 1
        max_nl_iter = 1 ! this unit test is linear 
        do i = 1, nodes_per_elem 
            do j = 1, nodes_per_elem
                norm_sin = real(global_coord(i,j))/real(global_coord(i,j))
                previous_elem_soln_vec(i,j) = sin(2.0*pi*norm_sin)
            end do
        end do
    end if

    cos_tot = 0.0
!---------------------------------------------------------------

!---Normal calculation flow - no need if doing unit test
    if (unit_test .eqv. .FALSE.) then
        precursor_soln_new(:,:) = 0 
        power_soln_new(:,:) = 0 
        !---Starting off total power
        total_power_prev = total_power_initial
        !--Calculate total power making sure not to double count boundary values, which are =
        do i = 1, non_fuel_start 
            do j = 1, nodes_per_elem - 1
                x_curr = real(global_coord(i,j) )
                x_last =  real(global_coord(non_fuel_start,3))
                norm_cos = ( (x_curr) - (x_last*0.5) )/ (0.5*x_last)
                cosine_term = cos( (pi/2)*norm_cos )
                cos_tot = cos_tot + cosine_term
            end do
        end do
        center_power_initial = total_power_initial/cos_tot
        cos_tot = 0.0
        
        nl_iter = 1 
        steady_state_flag = .TRUE.
        nonlinear_ss_flag = .TRUE.
         
        !---Initial guesses 
        center_temp_initial  = 800
        !---Apply to every node point within an element
        do i = 1,num_elem
            do j = 1, nodes_per_elem
                !---Apply to active fuel region
                if( i <= non_fuel_start ) then
                    call get_norm_coord(i,j,norm_cos) 
                    cosine_term = cos( (pi/2)*norm_cos )
                    !amplitude_fcn(i,j) = cosine_term
                    !if( i == 1) then
                    !    amplitude_fcn(i,j)  = 0
                    !    power_soln_new(i,j) = 0
                    !else
                    !    amplitude_fcn(i,j)  = 1.0 
                    !    power_soln_new(i,j) = 1.0
                    !end if
                    amplitude_fcn(i,j) = cosine_term 
                    cos_tot = cos_tot + cosine_term
                    power_soln_new(i,j) = (center_power_initial*cosine_term)
                    !---Set temperature distribution
                    temperature_soln_new(i,j) = (center_temp_initial*cosine_term)
                    temperature = temperature_soln_new(i,j)
                    !---Get density to set the velocity
                    call density_corr(temperature,density)
                    density_soln_new(i,j) = density
                    !---Need to get initial velocity distribution
                    velocity_soln_new(i,j) = mass_flow/(area*density)
                    !velocity_soln_new(i,j) = 100
                !---Inactive region assumed to have zero power 
                else
                    !---Temperature in inactive region same as end of active region ==> no loss
                    temperature_soln_new(i,j) = temperature_soln_new(non_fuel_start ,3)
                    !---Get density to set the velocity
                    call density_corr(temperature_soln_new(i,j),density)
                    density_soln_new = density
                    !---Need to get initial velocity distribution
                    velocity_soln_new(i,j) = mass_flow/(area*density)
                    !velocity_soln_new(i,j) = 100
                    power_soln_new(i,j) = 0.0
                    amplitude_fcn(i,j) = 0.0
                end if
            end do 
        end do
        
        do i = 1,num_elem 
            do j = 1,nodes_per_elem
                velocity_soln_prev(i,j) = velocity_soln_new(i,j) 
            end do
        end do 
        !-------------------------------------------------------------------------------
        !---Write out initial solution
        write(outfile_unit,fmt='(a)'), ' '
        write(outfile_unit,fmt='(a)'), 'Initial Power distribution '
        write(outfile_unit,fmt='(a)'), 'Position(x) Power [W]'
        do i = 1,num_elem
            do j = 1, nodes_per_elem
                write(outfile_unit, fmt='(f6.3, 12es14.3)')  global_coord(i,j), power_soln_new(i,j)
            end do
        end do
        
        !-------------------------------------------------------------------------------
        !---Write out initial solution
        write(outfile_unit,fmt='(a)'), ' '
        write(outfile_unit,fmt='(a)'), 'Initial temperature distribution '
        write(outfile_unit,fmt='(a)'), 'Position(x) Temperature [K]'
        do i = 1, num_elem 
            do j = 1, nodes_per_elem
                write(outfile_unit, fmt='(f6.3, 12es14.3)')  global_coord(i,j), temperature_soln_new(i,j)
            end do
        end do
        !---Write out initial solution
        write(outfile_unit,fmt='(a)'), ' '
        write(outfile_unit,fmt='(a)'), 'Initial velocity distribution '
        write(outfile_unit,fmt='(a)'), 'Position(x) Velocity [m/s]'
        do i = 1, num_elem 
            do j = 1, nodes_per_elem
                write(outfile_unit, fmt='(f6.3, 12es14.3)')  global_coord(i,j), velocity_soln_new(i,j)
            end do 
        end do
        !---Write out initial solution
        write(outfile_unit,fmt='(a)'), ' '
        write(outfile_unit,fmt='(a)'), 'Initial density distribition '
        write(outfile_unit,fmt='(a)'), 'Position(x) Density [kg/m^3]'
        do i = 1, num_elem 
            do j = 1, nodes_per_elem
                write(outfile_unit, fmt='(f6.3, 12es14.3)')  global_coord(i,j), density_soln_new(i,j)
            end do
        end do
        !-------------------------------------------------------------------------------
        !---Write out initial solution
        write(outfile_unit,fmt='(a)'), ' '
        write(outfile_unit,fmt='(a)'), 'Initial condition '
        write(outfile_unit,fmt='(a)'), 'Position(x) Precursor Concentration'
        do i = 1, num_elem 
            do j = 1, nodes_per_elem
                write(outfile_unit, fmt='(f6.3, 12es14.3)')  global_coord(i,j), precursor_soln_new(i,j)
            end do
        end do
  
        !-------------------------------------------------------------------------------
        !---Steady state solver for nonlinear ss problems
        if(nonlinear_ss_flag .eqv. .TRUE.) then 
            write(outfile_unit, fmt='(a)'), ' ' 
            write(outfile_unit, fmt='(a)'), 'Start steady state calculation'
        !---Nonlinear loop
            do 
                do n = 1, num_elem
                    !---Computer spatial matrices 
                    call element_matrix(n,nl_iter) 
                    !---Assemble K, F
                    call assemble_matrix(n)
                    !---Solve steady state system 
                    call solve_soln_steady(n,nl_iter)
                end do

                !-------------------------------------------------------------------------------
                !---Write out initial solution
                write(outfile_unit,fmt='(a)'), ' '
                write(outfile_unit,fmt='(a,1I4)'), 'Steady state soln at nl iter', nl_iter
                write(outfile_unit,fmt='(a)'), 'Position(x) Precursor Concentration'
                do i = 1, num_elem 
                    do j = 1, nodes_per_elem
                        write(outfile_unit, fmt='(f6.3, 12es14.3)')  global_coord(i,j), precursor_soln_new(i,j)
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
   
    !---Make the final calculated solution the 'previous' solution
    do i = 1, num_elem
        do j = 1, nodes_per_elem
            precursor_soln_prev(i,j) = precursor_soln_new(i,j)
        end do
    end do

    !---write out converged solution for plotting
    write(soln_outfile_unit,fmt='(a)'), 'Position(x) | Precursor Concentration'
    do i = 1, num_elem 
        do j = 1, nodes_per_elem
            write(soln_outfile_unit, fmt='(f6.3, 12es14.3)')  global_coord(i,j), precursor_soln_new(i,j)
        end do
    end do
!---Set steady state flag off
    steady_state_flag = .FALSE.

end

!------------------------------------------------------------------
subroutine get_norm_coord(i,j,norm_cos)
    USE parameters_fe

    implicit none
    
!---Dummy
    integer, intent(in) :: i
    integer, intent(in) :: j
    real,    intent(out) :: norm_cos

!---Local
    real :: x_curr, x_last

!---Get current global coordinate
    x_curr =  real(global_coord(i,j) )
    !---Last global coordinate
    x_last =  real(global_coord(non_fuel_start,3))
    !---Normalize coordinates so we go from -1 to 1
    norm_cos = ( (x_curr) - (x_last*0.5) )/ (0.5*x_last)

end subroutine get_norm_coord
