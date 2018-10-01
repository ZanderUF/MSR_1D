! Transient solver
! Notes: Solves precursor and power equations using foward Euler
!        Explicit method 
! Input: none
!
! Output:
! 
subroutine transient_forward_euler()

   USE parameters_fe

implicit none

!---Dummy

!---Local
    integer :: f,g,n,i,j
    real    :: t1 
    integer :: power_write_unit
    character(len=24) :: time_soln_name
    character(len=10) :: time_characters
    real(kind=4) :: temp_time

!---Start time-dependent solve
    temperature_soln_prev = temperature_soln_new

    transient = .TRUE.
    if ( transient .eqv. .TRUE. ) then
        write(outfile_unit, fmt='(a)'), ' ' 
        write(outfile_unit, fmt='(a)'), 'In transient loop'
        timeloop: do 
            !---Create element matrices and assemble
            elements_loop: do n = 1 , num_elem 
                !---Generate spatial matrices
                call spatial_matrices(n,1)
                isotope_loop: do f = 1, num_isotopes
                    delay_loop: do g = 1, num_delay_group
                        !---Assemble matrices solve elemental coefficients 
                        call assemble_matrix_transient(f,g,n) 
                        !---Solve for the elemental solution
                        call solve_precursor_forward_euler(f,g,n,1)
                    enddo delay_loop 
                enddo isotope_loop 
                
                if( mass_flow > 0.0 ) then
                    !---Solve for temperature
                    call solve_temperature(n)
                    !---Solve for velocity
                    call solve_velocity(n)
                end if

            enddo elements_loop 
            
            !---Swap solutions
            precursor_soln_prev   = precursor_soln_new
            power_amplitude_prev  = power_amplitude_new
            !---Solve for total power after spatial sweep through precursors
            call solve_power_forward_euler(1,t0) 
            
            !---Write solution to a file periodically
            transient_save_flag = .TRUE.
            if( modulo(t0,save_time_interval) < delta_t) then
                call write_out_soln(12, num_elem, transient_save_flag )
                !---Write out power solution 
                power_write_unit = 17
                temp_time=t0 
                time_soln_name = 'power_soln_at_time_step_'
                write(time_characters,'(f10.2)' ) temp_time
                time_characters = adjustl(time_characters)

                open (unit=power_write_unit, file= time_soln_name//time_characters,&
                status='unknown',form='formatted',position='asis')
 
                write(power_write_unit,fmt='(a,12es14.3)'), 'Power distribution at time:',&
                      t0  
                write(power_write_unit,fmt='(a)'), 'Position(x) Power [n/cm^3*s]'
                do i = 1,num_elem
                    do j = 1, nodes_per_elem
                        write(power_write_unit, fmt='(f6.3, 12es14.3)') &
                              global_coord(i,j), power_soln_new(i,j)
                    end do
                end do

                close(power_write_unit)

            end if !---End write out solution
            
            transient_save_flag = .FALSE.
	        !---Write power amp out @ every time step
	        if(t0 == 0.0) then
	            write(power_outfile_unit, ('(a)')), &
				'Time (s) | Power Amp | Norm Power | Reactivity | Beta Correction'
	        end if
	        write(power_outfile_unit, ('(12es14.6 ,12es14.5, 12es14.5, 12es14.5, 12es14.5)')), &
	          t0,power_amplitude_new,power_amplitude_new/power_amplitude_start,reactivity, beta_correction
            
	        !---Swap solutions
            precursor_soln_prev  = precursor_soln_new 
            power_soln_prev = power_soln_new
            power_amplitude_prev = power_amplitude_new
            if( mass_flow > 0.0 ) then
                temperature_soln_prev = temperature_soln_new
                velocity_soln_prev    = velocity_soln_new
            end if
            
            !---Stop if we've exceeded TMAX.
            if ( tmax <= t0 ) then
                exit
            end if
            
            t1 = t0 + delta_t

            t0 = t1
       
       enddo timeloop
    
    end if !---End transient if

end subroutine transient_forward_euler
