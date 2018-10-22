! Initialize the power, temperature, and velocity distributions 
! based on the inputted data
! Input:
!
! Output:
!

subroutine initialize()

    USE global_parameters_M
    USE mesh_info_M
    USE material_info_M
    USE flags_M
    USE solution_vectors_M

implicit none

!---Dummy

!---Local
    integer :: i,jj,j,n
    real(dp) :: ii, elem_length, density, temperature_initial,  &
                norm_cos, cosine_term, x_last, x_curr,&
                temperature
    real(dp) :: constant_velocity
    logical :: constant_flag
    
    !---Initialize to zero 
    precursor_soln_new(:,:,:,:) = 0.0_dp 
    power_soln_new(:,:)         = 0.0_dp
   
    !---Power amplitude set
    power_amplitude_new   = 1.0_dp
    power_amplitude_prev  = power_amplitude_new 
    power_amplitude_start = power_amplitude_new 
    
    steady_state_flag = .TRUE.
     
    !---Initial guesses 
    temperature_initial  = 800.0_dp
    !---Constant velocity for testing
    constant_velocity = 10.0_dp ! [cm/s]
    
    !---Flag for testing, use a flat power function or not
    constant_flag = .TRUE.
    !---Create spatial power function
    do i = 1, num_elem
        do j = 1, nodes_per_elem
            !---FUEL region
            if( (fuel_region_start <= i) .AND.  (i <= fuel_region_start) ) then
                !---Flat spatial shape 
                if(constant_flag .eqv. .TRUE.) then
                    spatial_power_fcn(i,j) = 1.0
                else
                    !---Cosine spatial shape 
                    call get_norm_coord(i,j,norm_cos) 
                    cosine_term = cos( (pi/2)*norm_cos )
                    spatial_power_fcn(i,j) = cosine_term
                end if
            else ! NON-FUEL region
                spatial_power_fcn(i,j) = 0.0
            end if
        end do
    end do

    !---Apply to every node point within an element
    do i = 1,num_elem
        do j = 1, nodes_per_elem
            !---Apply to active fuel region
            if( (fuel_region_start <= i) .AND.  (i <= fuel_region_end) ) then
                power_soln_new(i,j) = spatial_power_fcn(i,j)*power_amplitude_new
                !---Set temperature distribution
                temperature_soln_new(i,j) = (temperature_initial*spatial_power_fcn(i,j))
                temperature = temperature_soln_new(i,j)
                !---Get density to set the velocity
                call density_corr(temperature,density)
                density_soln_new(i,j) = density
                !---Need to get initial velocity distribution
                velocity_soln_new(i,j) = mass_flow/(area_core*density)
                !velocity_soln_new(i,j) = constant_velocity 
                area_variation(i,j) = area_core
            else
            !---Inactive region assumed to have zero power 
                power_soln_new(i,j) = 0.0
                !---Temperature in inactive region same as end of active region ==> no loss
                temperature_soln_new(i,j) = temperature_soln_new(fuel_region_start, 3)
                temperature=temperature_soln_new(i,j)
                !---Get density to set the velocity
                call density_corr(temperature,density)
                density_soln_new = density
                !---Need to get initial velocity distribution
                velocity_soln_new(i,j) = mass_flow/(area_pipe*density)
                !velocity_soln_new(i,j) = constant_velocity 
                !---Area change for the piping 
                area_variation(i,j) = area_pipe
            
            end if
        end do 
    end do

!---Get the average starting temperature
    total_temperature_initial = 0.0
    do i = fuel_region_start, fuel_region_end 
        do j = 1, nodes_per_elem
           total_temperature_initial = total_temperature_initial + temperature_soln_new(i,j) 
        end do
    end do
    avg_temperature_initial = total_temperature_initial/(fuel_region_start-fuel_region_end)

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

!------------------------------------------------------------------
subroutine get_norm_coord(i,j,norm_cos)
    
    USE global_parameters_M
    USE mesh_info_M

    implicit none
    
!---Dummy
    integer, intent(in) :: i
    integer, intent(in) :: j
    real(dp),    intent(out) :: norm_cos

!---Local
    real :: x_curr, x_last

!---Get curren global coordinate
    x_curr =   global_coord(i,j) 
    !---Last global coordinate
    x_last =  global_coord(fuel_region_end,3)
    !---Normalize coordinates so we go from -1 to 1
    norm_cos = ( (x_curr) - (x_last*0.5) )/ (fuel_region_end - fuel_region_start)

end subroutine get_norm_coord
