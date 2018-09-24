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
    real :: save_time_interval

!---Start time-dependent solve
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
            
                !---Solve for temperature
                call solve_temperature(n)
                !---Solve for velocity
                call solve_velocity(n)

            enddo elements_loop 
            
            !---Swap solutions
            !precursor_soln_prev   = precursor_soln_new
            !power_amplitude_prev  = power_amplitude_new
            !---Solve for total power after spatial sweep through precursors
            call solve_power_forward_euler(1,t0) 
            
            save_time_interval = 1.0 
            transient_save_flag = .TRUE.
            !---Write solution to a file periodically
            if( modulo(t0,save_time_interval) < delta_t) then
                call write_out_soln(12, num_elem, transient_save_flag )
            end if
            transient_save_flag = .FALSE.
	        !---Write power amp out @ every time step
	        if(t0 == 0.0) then
	            write(power_outfile_unit, ('(a)')), &
				'Time (s) | Power Amp | Norm Power | Reactivity'
	        end if
	        write(power_outfile_unit, ('(12es14.6 ,12es14.5, 12es14.5, 12es14.5)')), &
	          t0,power_amplitude_new,power_amplitude_new/power_amplitude_start,reactivity
            
	        !---Swap solutions
            precursor_soln_prev  = precursor_soln_new 
            power_soln_prev = power_soln_new
            power_amplitude_prev = power_amplitude_new
            temperature_soln_prev = temperature_soln_new
            velocity_soln_prev    = velocity_soln_new

            !---Stop if we've exceeded TMAX.
            if ( tmax <= t0 ) then
                exit
            end if
            
            t1 = t0 + delta_t

            t0 = t1
       
       enddo timeloop
    
    end if !---End transient if

end subroutine transient_forward_euler
