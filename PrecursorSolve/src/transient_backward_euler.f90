! Transient solver
! Notes: Solves precursor and power equations using Backward Euler
!        Implicit
! Input: none
!
! Output:
! 
subroutine transient_backward_euler()

   USE parameters_fe

implicit none

!---Dummy

!---Local
    integer :: f,g,n,i,j,nl_iter, power_write_unit
    real    :: t1
    character(len=24) :: time_soln_name
    character(len=10) :: time_characters
    real(kind=4) :: temp_time
    real, dimension(num_isotopes, num_delay_group) :: L2_norm_current, L2_norm_prev
    real :: nl_iter_tolerance, difference_L2
    integer :: difference_counter, abs_max_nl_iter    

!---Set to make sure we don't iterate forever if we are not converging
    abs_max_nl_iter = 500 
    nl_iter_tolerance = 1E-12
!---Start time-dependent solve
    transient = .TRUE.
    if ( transient .eqv. .TRUE. ) then
        write(outfile_unit, fmt='(a)'), ' ' 
        write(outfile_unit, fmt='(a)'), 'In backward Euler transient loop'
        timeloop: do!---Time loop 
            nl_iter = 1 
            L2_norm_prev = 0.0
            L2_norm_current = 0.0
            difference_counter = 0
            nonlinearloop: do  
                !---Create element matrices and assemble
                elements_loop: do n = 1 , num_elem 
                    !---Generate spatial matrices
                    call spatial_matrices(n,nl_iter)
                    call numerical_flux_matrices(n,nl_iter)
                    isotope_loop: do f = 1, num_isotopes
                        delay_loop: do g = 1, num_delay_group
                            !---Assemble matrices solve elemental coefficients 
                            call assemble_matrix_transient(f,g,n) 
                            !---Solve for the elemental solution
                            call solve_precursor_backward_euler(f,g,n,nl_iter)
                        enddo delay_loop 
                    enddo isotope_loop 
                   
                    if( mass_flow > 0.0) then
                        call solve_temperature(n)
                        call solve_velocity(n)
                    end if

                enddo elements_loop 
                
                precursor_soln_prev = precursor_soln_new
                !power_amplitude_prev = power_amplitude_new
                
                !---Solve for total power after spatial sweep through precursors
                call solve_power_backward_euler(nl_iter,t0) 
                
                !---Calculate L2 norm of precursor solution
                if(nl_iter > 1) then
                    do f = 1, num_isotopes
                        do g = 1, num_delay_group
                            L2_norm_current(f,g) = sqrt( sum( precursor_soln_new(f,g,:,:)*&
                                                   precursor_soln_new(f,g,:,:) ) ) 
                        end do
                    end do
                    
                    !---Calculate the difference in the L2 norms between iterations
                    do f = 1, num_isotopes
                        do g = 1, num_delay_group
                            difference_L2 = abs( L2_norm_prev(f,g) - L2_norm_current(f,g) )

                            if( difference_L2 < nl_iter_tolerance) then
                                difference_counter = difference_counter + 1
                            end if
                        end do
                    end do
                    !---Need to make sure the L2 norm converges for all precursor groups
                    if ( difference_counter == num_delay_group) then
                        max_nl_iter = nl_iter - 1
                    end if
                    
                    !print *,'difference_L2',difference_L2
                    !---Swap for next iteration
                    L2_norm_prev = L2_norm_current

                end if

                nl_iter = nl_iter + 1 !---Nonlinear iteration counter
                
                !---Check if too many nonlinear iterations and not converging
                if ( nl_iter > max_nl_iter) then
                    if(nl_iter > abs_max_nl_iter) then
                        write(outfile_unit,fmt=('(a)')) 'Gone past max amount of &
                                     nonlinear iterations &
                                     and might have a problem'
                    else
                        if (DEBUG .eqv. .TRUE.) then 
                            write(outfile_unit,fmt=('(a,I3,a,8es14.3)')) &
                            'Took this # of iterations to converge --> ',nl_iter,&
                            '  <-- at time step -->', t0
                        end if
                    end if
                    
                    exit
                end if 
            
            !print *,'difference_counter',difference_counter   
            enddo nonlinearloop 
            transient_save_flag = .TRUE.
            !print *,'nl_iter ',nl_iter 
            !---Write solution to a file periodically
            if( modulo(t0,save_time_interval) < delta_t) then
                
                call write_out_soln(12, num_elem, transient_save_flag )
                
                transient_save_flag = .FALSE.
                !---Write out power solution 
                power_write_unit = 17
                temp_time = t0 
                time_soln_name = 'power_soln_at_time_step_'
                write(time_characters,'(f10.2)' ) temp_time
                time_characters = adjustl(time_characters)

                open (unit=power_write_unit, file= time_soln_name//time_characters,&
                status='unknown',form='formatted',position='asis')
 
                write(power_write_unit,fmt='(a,12es14.3)'), &
                'Power distribution at time:',t0
                write(power_write_unit,fmt='(a)'), 'Position(x) Power [n/cm^3*s]'
                
                do i = 1,num_elem
                    do j = 1, nodes_per_elem
                        write(power_write_unit, fmt='(f6.3, 12es14.3)') &
                              global_coord(i,j), power_soln_new(i,j)
                    end do
                end do

                close(power_write_unit)

            end if !---End write out solution
            !---Write power amp out @ every time step
            if(t0 == 0.0) then
                write(power_outfile_unit, ('(a)')),&
                      'Time (s) | Power Amp | Norm Power | Reactivity | &
                       Beta Correction | Reactivity Feedback'
            end if

            write(power_outfile_unit, ('(12es14.6 ,12es14.5, 12es14.5, 12es14.5,&
            12es14.5,12es14.5)')), &
            t0, power_amplitude_new, power_amplitude_new,&
            reactivity, beta_correction, reactivity_feedback


            !---Swap solutions
            precursor_soln_prev       = precursor_soln_new 
            power_amplitude_prev      = power_amplitude_new
            precursor_soln_last_time  = precursor_soln_new
            power_amplitude_last_time = power_amplitude_new
            
            if(mass_flow > 0.0 ) then
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

    end if!---End transient if

end subroutine transient_backward_euler
