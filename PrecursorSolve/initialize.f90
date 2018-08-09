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

    !---Initialize
    precursor_soln_new(:,:,:,:) = 0 
    power_soln_new(:,:) = 0 
    
    !---Starting off total power
    total_power_initial = total_power_initial/global_coord(num_elem,3)
    total_power_prev = total_power_initial
    
    steady_state_flag = .TRUE.
    nonlinear_ss_flag = .TRUE.
     
    !---Initial guesses 
    center_temp_initial  = 800
    !---Constant velocity for testing
    constant_velocity = 1000.0 ! [cm/s]

    !---Apply to every node point within an element
    do i = 1,num_elem
        do j = 1, nodes_per_elem
            !---Apply to active fuel region
            if( i <= non_fuel_start ) then
                call get_norm_coord(i,j,norm_cos) 
                cosine_term = cos( (pi/2)*norm_cos )
                cosine_term = cosine_term 
                !amplitude_fcn(i,j) = cosine_term
                if( i == 1) then
                    amplitude_fcn(i,j)  = 0
                    power_soln_new(i,j) = 0
                else
                    amplitude_fcn(i,j)  = 1.0 
                    power_soln_new(i,j) = 1.0
                end if
                
                !amplitude_fcn(i,j) = cosine_term 
                !amplitude_fcn(i,j) = 1.0
                !power_soln_new(i,j) = (total_power_initial*cosine_term)
                !power_soln_new(i,j) =  1.0 
                !---Set temperature distribution
                temperature_soln_new(i,j) = (center_temp_initial*cosine_term)
                temperature = temperature_soln_new(i,j)
                !---Get density to set the velocity
                call density_corr(temperature,density)
                density_soln_new(i,j) = density
                !---Need to get initial velocity distribution
                !velocity_soln_new(i,j) = mass_flow/(area*density)
                velocity_soln_new(i,j) = constant_velocity 
                !velocity_soln_new(i,j) = 0
            !---Inactive region assumed to have zero power 
            else
                !---Temperature in inactive region same as end of active region ==> no loss
                temperature_soln_new(i,j) = temperature_soln_new(non_fuel_start ,3)
                !---Get density to set the velocity
                call density_corr(temperature_soln_new(i,j),density)
                density_soln_new = density
                !---Need to get initial velocity distribution
                !velocity_soln_new(i,j) = mass_flow/(area*density)
                !velocity_soln_new(i,j) = 0
                power_soln_new(i,j) = 0.0
                amplitude_fcn(i,j) = 0.0
                velocity_soln_new(i,j) = constant_velocity 
            end if
        end do 
    end do

!-------------------------------------------------------------------------------
    !---Write out initial solution
    write(outfile_unit,fmt='(a)'), ' '
    write(outfile_unit,fmt='(a)'), 'Initial Power distribution '
    write(outfile_unit,fmt='(a)'), 'Position(x) Power [n/cm^3*s]'
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
    write(outfile_unit,fmt='(a)'), 'Position(x) Velocity [cm/s]'
    do i = 1, num_elem 
        do j = 1, nodes_per_elem
            write(outfile_unit, fmt='(f6.3, 12es14.3)')  global_coord(i,j), velocity_soln_new(i,j)
        end do 
    end do
    !---Write out initial solution
    write(outfile_unit,fmt='(a)'), ' '
    write(outfile_unit,fmt='(a)'), 'Initial density distribition '
    write(outfile_unit,fmt='(a)'), 'Position(x) Density [g/cm^3]'
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
    !do i = 1, num_elem 
    !    do j = 1, nodes_per_elem
    !        write(outfile_unit, fmt='(f6.3, 12es14.3)')  global_coord(i,j), precursor_soln_new(i,j)
    !    end do
    !end do
end subroutine
