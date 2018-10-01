! Initialize the system 
!
! Input:
!
! Output:
! 
subroutine initialize()

   USE parameters_fe

implicit none

!---Dummy

!---Local
    integer :: i,jj,j,n
    real :: ii, elem_length, density, center_temp_initial,  &
             norm_sin,norm_cos, pi, cosine_term, x_last, x_curr,&
             temperature
    parameter (pi = 3.1415926535897932)
    real :: constant_velocity
    logical :: constant_flag

    !---Initialize to zero 
    precursor_soln_new(:,:,:,:) = 0 
    power_soln_new(:,:) = 0 
   
    !---Power amplitude set
    
    power_amplitude_new = 1.0
    power_amplitude_prev = power_amplitude_new 
    power_amplitude_start = power_amplitude_new 
    
    steady_state_flag = .TRUE.
    nonlinear_ss_flag = .TRUE.
     
    !---Initial guesses 
    center_temp_initial  = 800
    !---Constant velocity for testing
    constant_velocity = 1.5! [cm/s]
    
    !---Flag for testing, use a flat power function or not
    constant_flag = .TRUE.
    !---Create spatial power function
    do i = 1, num_elem
        do j = 1, nodes_per_elem
            if( i <= non_fuel_start) then
                !---Flat spatial shape 
                if(constant_flag .eqv. .TRUE.) then
                    spatial_power_fcn(i,j) = 1.0
                else
                    !---Cosine spatial shape 
                    call get_norm_coord(i,j,norm_cos) 
                    cosine_term = cos( (pi/2)*norm_cos )
                    spatial_power_fcn(i,j) = cosine_term
                end if
            else
                spatial_power_fcn(i,j) = 0.0
            end if
        end do
    end do

    print *,'totalpowerinitial', total_power_initial
    
    !---Apply to every node point within an element
    do i = 1,num_elem
        do j = 1, nodes_per_elem
            !---Apply to active fuel region
            if( i <= non_fuel_start ) then
                !print *,'total_power/fuel',total_power_initial/non_fuel_start 
                power_soln_new(i,j) = spatial_power_fcn(i,j)*power_amplitude_new
                !---Set temperature distribution
                temperature_soln_new(i,j) = (center_temp_initial*spatial_power_fcn(i,j))
                temperature = temperature_soln_new(i,j)
                !---Get density to set the velocity
                call density_corr(temperature,density)
                density_soln_new(i,j) = density
                !---Need to get initial velocity distribution
                velocity_soln_new(i,j) = mass_flow/(area*density)
                !velocity_soln_new(i,j) = constant_velocity 
            !---Inactive region assumed to have zero power 
            else
                power_soln_new(i,j) = 0.0
                !---Temperature in inactive region same as end of active region ==> no loss
                temperature_soln_new(i,j) = temperature_soln_new(non_fuel_start ,3)
                !---Get density to set the velocity
                call density_corr(temperature_soln_new(i,j),density)
                density_soln_new = density
                !---Need to get initial velocity distribution
                velocity_soln_new(i,j) = mass_flow/(area*density)
                !velocity_soln_new(i,j) = 0
                !velocity_soln_new(i,j) = constant_velocity 
            end if
        end do 
    end do
    
!-------------------------------------------------------------------------------
!---Write out initial solution
    write(outfile_unit,fmt='(a)'), ' '
    write(outfile_unit,fmt='(a)'), 'Initial Spatial Power distribution '
    write(outfile_unit,fmt='(a,12es14.3)'),'Initial power amplitude', power_amplitude_new
    write(outfile_unit,fmt='(a)'), 'Position(x) Power [n/cm^3*s]'
    write(outfile_unit,fmt='(a)'), '------------------------------------'
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
    write(outfile_unit,fmt='(a)'), '------------------------------------'
    do i = 1, num_elem 
        do j = 1, nodes_per_elem
            write(outfile_unit, fmt='(f6.3, 12es14.3)')  global_coord(i,j), temperature_soln_new(i,j)
        end do
    end do
!---Write out initial solution
    write(outfile_unit,fmt='(a)'), ' '
    write(outfile_unit,fmt='(a)'), 'Initial velocity distribution '
    write(outfile_unit,fmt='(a)'), 'Position(x) Velocity [cm/s]'
    write(outfile_unit,fmt='(a)'), '------------------------------------'
    do i = 1, num_elem 
        do j = 1, nodes_per_elem
            write(outfile_unit, fmt='(f6.3, 12es14.3)')  global_coord(i,j), velocity_soln_new(i,j)
        end do 
    end do
    !---Write out initial solution
    write(outfile_unit,fmt='(a)'), ' '
    write(outfile_unit,fmt='(a)'), 'Initial density distribition '
    write(outfile_unit,fmt='(a)'), 'Position(x) Density [g/cm^3]'
    write(outfile_unit,fmt='(a)'), '------------------------------------'
    do i = 1, num_elem 
        do j = 1, nodes_per_elem
            write(outfile_unit, fmt='(f6.3, 12es14.3)')  global_coord(i,j), density_soln_new(i,j)
        end do
    end do
    
end subroutine
