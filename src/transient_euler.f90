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
    use Mod_SetupOutputFiles
    
    implicit none

!---Dummy

!---Local
    integer :: f,g,n,i,j,nl_iter
    real(dp), dimension(num_isotopes, num_delay_group) :: L2_norm_current, L2_norm_prev,&
                                                          l2_residual
    real(dp) :: t1
    integer  :: difference_counter
    real(dp) :: sum_residual
    real(dp) :: node_temperature, total_temperature
    integer  :: TimeIndex


!---Start time-dependent solve
    if ( time_solve .eqv. .TRUE. ) then
        
        write(outfile_unit, fmt='(a)'), ' ' 
        write(outfile_unit, fmt='(a)'), 'In transient loop'
        !---
        TimeIndex = 1

        timeloop: do!---Time loop 
           
            !---Reset convergence criteria counters
            nl_iter = 1 
            L2_norm_prev        = 0.0
            L2_norm_current     = 0.0
            difference_counter  = 0.0
          
            !---Interpolate beta as a function of flow speed
            call beta_feedback 
           
            !---Decide if we are doing backward or forward euler
            nl_iter_flag      = .TRUE. 
            residual(:,:,:,:) = 0.0

            !---Nonlinear iteration for implicit time stepping schemes 
            nonlinearloop: do  
                if(nl_iter_flag .eqv. .TRUE.) then 
                    !---Create element matrices and assemble
                    elements_loop: do n = 1 , num_elem 
                        !---Generate spatial matrices 
                        call spatial_matrices(n,nl_iter)
                        !---Generate numerical flux matrices for boundary of elements
                        call numerical_flux_matrices(n,nl_iter)
                        !---If it is a stationary fuel test case do not need temperature 
                        if( mass_flow > 0.0) then
                                if (feedback_method >=1) then
                                    call solve_velocity(n,nl_iter)
                                    call solve_temperature_euler(n,nl_iter)
                                end if
                                call solve_velocity(n,nl_iter)
                                call solve_temperature_euler(n,nl_iter)
                        end if            
                        !---Loop over all fissile isotopes
                        isotope_loop: do f = 1, num_isotopes
                            !---Loop over precursor families
                            delay_loop: do g = 1, num_delay_group
                                !---Assemble matrices to solve for elemental coefficients 
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
                    !---Calculate the L2 norms to determine convergence on the time step
                    call l2_norm(nl_iter)
                end if
                
                !---Reached max # iterations exit and move to next time step
                if ( nl_iter == max_nl_iter) then
                    exit
                end if
                
                nl_iter = nl_iter + 1 
               
                !---No nonlinear iterations for explicit cases
                else
                    exit
                end if
               
                !---Swap solutions vectors for next nonlinear iteration
                precursor_soln_prev       = precursor_soln_new 
                power_amplitude_prev      = power_amplitude_new
                power_soln_prev           = power_soln_new
                temperature_soln_prev     = temperature_soln_new
                velocity_soln_prev        = velocity_soln_new
                density_soln_prev         = density_soln_new

            enddo nonlinearloop 
 
            !---Reset convergence flag for temperature for next time step
            temperature_converged = .FALSE.

            !---Write information about nonlinear iterations and residual 
            if(DEBUG .eqv. .TRUE.) then
                write(nl_outfile_unit, fmt='(es16.6 ,1I6,16es16.6,16es16.6,&
                                             16es16.6,16es16.6,16es16.6,16es16.6)'),&
                                             t0,nl_iter, (l2_residual(1,g), g=1,num_delay_group)
            end if

            transient_save_flag = .TRUE.
            
            !---Swap solution vectors for next time step 
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
            
            total_temperature = 0.0_dp
            peak_temperature  = 0.0_dp

            !---Calculate average temperature across core
            do i = Fuel_Inlet_Start, Fuel_Outlet_End
                node_temperature = 0.0_dp 
                !---Get value over node
                do j = 1, nodes_per_elem
                    node_temperature = node_temperature + temperature_soln_prev(i,j)*vol_int(j)
                end do
                total_temperature = total_temperature + node_temperature
                !---Find peak
                if (peak_temperature < node_temperature) then
                    peak_temperature = node_temperature
                end if
            end do

            average_temperature = total_temperature/(Fuel_Outlet_End - Fuel_Inlet_Start)
           
                      
            !---Write final values to large time arrays
            TimeAll(TimeIndex)              = t0 
            PowerAllTime(TimeIndex)         = power_amplitude_new  
            RhoInsertedAll(TimeIndex)       = reactivity
            DopplerFeedbackAll(TimeIndex)   = total_temperature_feedback
            DensityFeedbackAll(TimeIndex)   = total_density_feedback 
            MassFlowAll(TimeIndex)          = mass_flow
            TempHeatExchangerAll(TimeIndex) = temperature_reduction 
            AvgTempAll(TimeIndex)           = average_temperature 
            PeakTempAll(TimeIndex)          = peak_temperature
            NonLinearIterAll(TimeIndex)     = nl_iter 
          
           ! write out solution vectors 
           if( modulo(t0,save_time_interval) < delta_t) then
                call write_out_soln(12, num_elem, transient_save_flag )
                
                transient_save_flag = .FALSE.
            end if

            !---Stop once we have reached end of desired simulation time.
            if ( TimeIndex == NumTimeSteps ) then
                exit
            end if
           !---Increment the time step
           t1 = t0 + delta_t
           !---Reset starting time point
           t0 = t1
           
           !---Index for keeping track of time steps
           TimeIndex = TimeIndex + 1
       
       enddo timeloop

       !---Write power, amplitude, reacitivty to file
       call write_periodic
            
    end if !---End transient if

!***************************************************
! Currently not functioning 4/11/2019
!---Only automate time stepping for implicit calculations
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


end subroutine transient_euler
