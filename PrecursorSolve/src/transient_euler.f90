! Transient solver
! Notes: Solves precursor and power equations using Backward Euler
!        Implicit
!
! Input: none
!
! Output:
! 
subroutine transient_euler()
    
    USE flags_M
    USE global_parameters_M
    USE solution_vectors_M
    USE time_info_M
    USE mesh_info_M
    USE material_info_M
    USE element_matrices_M

    implicit none

!---Dummy

!---Local
    integer :: f,g,n,i,j,nl_iter
    real(dp), dimension(num_isotopes, num_delay_group) :: L2_norm_current, L2_norm_prev,&
                                                          l2_residual
    real(dp) :: t1
    integer  :: difference_counter
    real(dp) :: sum_residual
   
!---Start time-dependent solve
    if ( time_solve .eqv. .TRUE. ) then
        
        write(outfile_unit, fmt='(a)'), ' ' 
        write(outfile_unit, fmt='(a)'), 'In backward Euler transient loop'
        
        timeloop: do!---Time loop 
            
            nl_iter = 1 
            L2_norm_prev = 0.0
            L2_norm_current = 0.0
            difference_counter = 0
          
            !---Vary beta with flow speed
            call beta_feedback 
           
            
            !---Decide if we are doing backward or forward euler
            nl_iter_flag = .TRUE. 
            residual(:,:,:,:) = 0.0_dp
           
            nonlinearloop: do  
                if(nl_iter_flag .eqv. .TRUE.) then 
                    !---Create element matrices and assemble
                    elements_loop: do n = 1 , num_elem 
                        !---Generate spatial matrices
                        call spatial_matrices(n,nl_iter)
                        call numerical_flux_matrices(n,nl_iter)
                         if( mass_flow > 0.0) then
                                call solve_velocity(n,nl_iter)
                                call solve_temperature_euler(n,nl_iter)
                        end if            
                        isotope_loop: do f = 1, num_isotopes
                            delay_loop: do g = 1, num_delay_group
                                !---Assemble matrices solve elemental coefficients 
                                call assemble_matrix_transient(f,g,n) 
                                !---Solve for the elemental solution
                                call solve_precursor_euler(f,g,n,nl_iter)
                            enddo delay_loop 
                        enddo isotope_loop 
                    enddo elements_loop 
                !---Solve for total power after spatial sweep through precursors
                call solve_power_euler(nl_iter,t0) 
                !---If doing forward Euler - no iteration
                    if(td_method_type == 0 ) then
                        nl_iter_flag = .FALSE. 
                    else
                        difference_counter = 0             
                        !---Calculate the l2 norm
                        call l2_norm(nl_iter, difference_counter, &
                                     L2_norm_prev,L2_norm_current)
                    end if
                if ( nl_iter == 100) then
                    exit
                end if
                
                nl_iter = nl_iter + 1 
               
                !---No nonlinear iterations
                else
                    exit
                end if
               
                !---Swap for next nonlinear iteration
                precursor_soln_prev       = precursor_soln_new 
                power_amplitude_prev      = power_amplitude_new
                
                power_soln_prev           = power_soln_new
                temperature_soln_prev     = temperature_soln_new
                velocity_soln_prev        = velocity_soln_new
                density_soln_prev         = density_soln_new

            enddo nonlinearloop 
            !print *,' nl iter ', nl_iter

            !***************************************************
            !---Only automate time stepping
            if(td_method_type == 1) then
                sum_residual = 0.0
                do f = 1, num_isotopes 
                    do g = 1, num_delay_group
                        l2_residual(f,g) = sqrt( sum(residual(f,g,:,:)*&
                                                     residual(f,g,:,:)))
                        do i = 1, num_elem
                           do j = 1, nodes_per_elem
                                !sum_residual(f,g,i) = sum_residual(f,g,i) + & 
                                !               residual(f,g,i,j)*vol_int(j)
                           end do
                        end do
                    
                    end do
                end do
            end if
           
            !write(nl_outfile_unit, fmt='(es16.6 ,1I6,16es16.6,16es16.6,&
            !                                16es16.6,16es16.6,16es16.6,16es16.6)'),&
            !                                t0,nl_iter, (l2_residual(1,g), g=1,num_delay_group)
            
            transient_save_flag = .TRUE.
            
            !---Swap solutions for next time step 
            precursor_soln_prev       = precursor_soln_new 
            power_amplitude_prev      = power_amplitude_new
            precursor_soln_last_time  = precursor_soln_new
            power_amplitude_last_time = power_amplitude_prev
            power_soln_starting       = power_soln_new
            power_soln_prev           = power_soln_new
            
            !---Only want to calculate velocity if mass flow is defined
            if(mass_flow > 0.0 ) then
                temperature_soln_starting = temperature_soln_new
                temperature_soln_prev     = temperature_soln_new
                velocity_soln_prev        = velocity_soln_new
                density_soln_prev         = density_soln_new
            end if

            !---Write power, amplitude, reacitivty to file
            call write_periodic

            !---Stop if we've exceeded TMAX.
            if ( tmax <= t0 ) then
                exit
            end if
          
           !---Increment the time step
           t1 = t0 + delta_t

           t0 = t1

       enddo timeloop

    end if !---End transient if

end subroutine transient_euler
