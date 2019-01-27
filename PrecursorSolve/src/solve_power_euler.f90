!  
! Notes: Solves for the power amplitude with backward Euler
!  
! Input: nl_iter - current nonlinear iteration
!        current_time - current time we are solving at
!
! Output:
! 

subroutine solve_power_euler(nl_iter, current_time)

    USE global_parameters_M
    USE time_info_M
    USE solution_vectors_M
    USE mesh_info_M
    USE material_info_M
    USE flags_M
    USE element_matrices_M

    implicit none

!---Dummy
    integer,  intent(in) :: nl_iter
    real(dp), intent(in) :: current_time

!---Local
    integer :: i,j,f,g
    real(dp), dimension(num_isotopes,num_delay_group) :: precursors_lambda_vec
    real(dp), dimension(num_elem) :: precursors_elem

    real(dp), dimension(num_delay_group) :: test_fuel_prec, test_total_prec
    real(dp):: power_new_total, total_precursor_ref,&
               total_precursor_ref_sum, total_spatial_fcn,&
               total_precursors_fuel, &
               rho_initial, step_time
    real(dp) :: first_zag, second_zag, third_zag, reactivity_zag
    real(dp) :: total_power
    real(dp) :: average_temperature, total_temperature
    real(dp) :: total_density,average_density

!---Initialize to zero
    precursors_lambda_vec(:,:) = 0.0_dp
    precursors_elem(:) = 0.0_dp

!---Calculate precursors in the fuel region per isotope and group
    do f = 1, num_isotopes
       do g = 1, num_delay_group
            do i = Fuel_Inlet_Start , Fuel_Outlet_End
                do j = 1, nodes_per_elem
                   !---Precursors*lambda
                   precursors_elem(i) = precursors_elem(i) + &
                        lamda_i_mat(f,g)*&
                        vol_int(j)*&
                        precursor_soln_prev(f,g,i,j)

                   precursors_lambda_vec(f,g) = precursors_lambda_vec(f,g) + &
                                      vol_int(j)*&
                                      precursor_soln_new(f,g,i,j)
                end do
            end do
       end do
    end do   

!---Get contribution of power from assumed spatial profile 
    total_spatial_fcn = 0.0
    total_power = 0.0

    do i = 1, num_elem
        do j = 1, nodes_per_elem
            total_power = total_power + &
                          total_power_initial*power_amplitude_prev*&
                          spatial_power_frac_fcn(i,j)*vol_int(j)
            
            total_spatial_fcn = total_spatial_fcn + &
                           total_power_initial*spatial_power_frac_fcn(i,j)*vol_int(j)

        end do
    end do

!---Calc total precursors fuel 
    total_precursors_fuel = 0.0_dp
    do f = 1, num_isotopes
        do g = 1, num_delay_group
            total_precursors_fuel = total_precursors_fuel + &
                    lamda_i_mat(f,g)*precursors_lambda_vec(f,g)
        end do
    end do

    !--- Need for benchmarks
    if(Read_DIF3D .eqv. .FALSE. ) then
        if(t0 == 0.0) then
            beta_correction = gen_time*total_precursors_fuel/total_power
        end if
    end if
    
!---STEP perturbation
    if(step_flag .eqv. .TRUE.) then
        if(step_start_time < t0 .and. t0 < step_end_time ) then
            reactivity = reactivity_input 
        elseif ( t0 > step_end_time) then
            reactivity = reactivity_input 
        end if
    end if 
!---End STEP perturbation

!---RAMP perturbation 
    if(ramp_flag .eqv. .TRUE.) then
        if(ramp_start_time < t0 .and. t0 < ramp_end_time) then   
            rho_initial = 0.0 
            reactivity = rho_initial + &
              ((reactivity_input - rho_initial)*&
               (t0 - t_initial))/(ramp_end_time - t_initial)
        elseif ( t0 > ramp_end_time) then
            reactivity = reactivity_input 
        end if
    end if  
!---End RAMP perturbation

!---Zig-zag perturbation
    reactivity_zag = 7.5E-3
    first_zag  = 0.5
    second_zag = 1.0
    third_zag  = 1.5
    rho_initial = 0.0
    
    if(zag_flag .eqv. .TRUE.) then
        if( t0 < first_zag ) then
            reactivity = rho_initial + &
               ((reactivity_zag )*&
               (t0-t_initial))/(1.0_dp - t_initial)
        elseif( t0 < second_zag) then
            reactivity = reactivity_zag/2.0_dp + &
               (( - reactivity_zag/2.0_dp )*&
               (t0-first_zag) )/(second_zag - first_zag)
        elseif ( t0 < third_zag ) then
            reactivity = &  
               (( reactivity_zag/2.0_dp )*&
               (t0-second_zag))/(third_zag - second_zag )
        else
            reactivity = reactivity
        end if
    end if

!---Calculate average temperature across the core    
    total_temperature = 0.0_dp
    total_density = 0.0_dp

    do i = Fuel_Inlet_Start, Fuel_Outlet_End
        do j = 1, nodes_per_elem
            total_temperature = total_temperature + &
                          vol_int(j)*temperature_soln_new(i,j)
            total_density = total_density + &
                            vol_int(j)*density_soln_prev(i,j)
        end do
    end do

    average_temperature = total_temperature/(Fuel_Outlet_End - Fuel_Inlet_Start)
    average_density     = total_density/(Fuel_Outlet_End - Fuel_Inlet_Start)

!---Calculate temperature reactivity feedback
    total_temperature_feedback = 0.0_dp
    total_density_feedback     = 0.0_dp

    if(feedback_method == 3 ) then
        !---Total_temperature feedback
        !do i = Fuel_Inlet_Start, Fuel_Outlet_End 
        !    Density_Reactivity_Feedback(i) = (spatial_expansion_fcn(i,2) / &
        !                            0.01_dp) * &
        !                           (density_soln_ss(i,2)/density_soln_prev(i,2) - 1.0_dp)
        !    do j = 1, nodes_per_elem

        !        total_temperature_feedback = total_temperature_feedback + &
        !                                     vol_int(j)*Temperature_Reactivity_Feedback(i,j)
        !    end do
        !end do
        !total_density_feedback = sum(Density_Reactivity_Feedback)
    
    !---Calculate reactivity change based on average
    total_temperature_feedback = (total_doppler_read_in/total_temperature_change)*&
                                 (average_temperature - average_temperature_ss)
    total_density_feedback = (total_expansion_read_in/(total_density_change*2.5935_dp))*&
                             (average_density_ss - average_density)                          

    end if
   

    reactivity_feedback = 0.0_dp 
    reactivity_feedback = total_temperature_feedback + total_density_feedback
    !reactivity_feedback = 0.0_dp 
 
!---Power Solve
    if(td_method_type == 0) then ! Forward Euler
         power_amplitude_new = power_amplitude_prev + &
                          delta_t*(( reactivity + reactivity_feedback  &
                          - beta_correction )/gen_time)*&
                          power_amplitude_prev &
                          + delta_t*(1.0_dp/total_spatial_fcn)*&
                          total_precursors_fuel
    end if

    if(td_method_type == 1) then ! Backward Euler
        power_amplitude_new = power_amplitude_last_time + &
                          delta_t*(( reactivity + reactivity_feedback  &
                          - beta_correction )/gen_time)*power_amplitude_prev &
                          + delta_t*(1.0_dp/total_spatial_fcn)*&
                          total_precursors_fuel
    end if
     
!---Project power onto spatial shape
    power_soln_new(:,:) = 0.0
    do i = Fuel_Inlet_Start, Fuel_Outlet_End 
        do j = 1, nodes_per_elem
            power_soln_new(i,j) = power_amplitude_new* &
                                  spatial_power_frac_fcn(i,j)      
        end do
    end do
    
!---End power solve

end subroutine solve_power_euler
