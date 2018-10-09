!  
! Notes: Solves for the power amplitude with backward Euler
!  
! Input: nl_iter - current nonlinear iteration
!        current_time - current time we are solving at
!
! Output:
! 

subroutine solve_power_backward_euler(nl_iter, current_time)

    USE parameters_fe

    implicit none

!---Dummy
    integer, intent(in) :: nl_iter
    double precision, intent(in) :: current_time

!---Local
    integer :: i,j,f,g
    real (kind=16), dimension(num_delay_group,num_elem) :: precursors_vec
    real (kind=16), dimension(num_delay_group) :: test_fuel_prec, test_total_prec

    real (kind=16), dimension(num_elem) :: temp_vec_num_elem, &
                                           precursors_lambda_vec, &
                                           power_soln_new_temp 
    real (kind=16):: power_new_total, total_precursor_ref,&
                     total_precursor_ref_sum, total_fuel_length,&
                     total_precursors_fuel, &
                     rho_initial, step_time
    
    real :: ramp_end_time, ramp_start_time, step_end_time, step_start_time
    real :: first_zag, second_zag, third_zag, reactivity_zag
    real :: temp_reactivity_feedback, total_power

!---Initialize to zero
    precursors_lambda_vec(:) = 0.0
    temp_vec_num_elem(:) = 0.0
    total_precursor_ref = 0.0

    precursors_vec = 0.0

!---Calculate total precursor concentration*lamda over system 
    do f = 1, num_isotopes 
       do g = 1, num_delay_group
            do i = 1, num_elem 
                do j = 1, nodes_per_elem
                   precursors_lambda_vec(i) = precursors_lambda_vec(i) + &
                                      lamda_i_mat(f,g)*&
                                      elem_vol_int(i,j)*&
                                      precursor_soln_prev(f,g,i,j)
                   
                   precursors_vec(g,i) = precursors_vec(g,i) + &
                                       elem_vol_int(i,j)*  &
                                       precursor_soln_prev(f,g,i,j)
    
                   total_precursor_ref = total_precursor_ref + &
                                      lamda_i_mat(f,g)*&
                                      elem_vol_int(i,j)*&
                                      precursor_soln_prev(f,g,i,j)
                end do
            end do
       end do
    end do
   
    total_power = 0.0 
    !---Get total length of the fuel element
    total_fuel_length = 0.0
    do i = 1, num_elem
        do j = 1, nodes_per_elem
            total_power = total_power + &
                        elem_vol_int(i,j)*power_amplitude_prev*&
                        spatial_power_fcn(i,j)

            total_fuel_length = total_fuel_length + &
                                spatial_power_fcn(i,j)*elem_vol_int(i,j)
        end do
    end do
    
    total_precursor_ref_sum   = sum(precursors_lambda_vec)
    total_precursors_fuel     = sum(precursors_lambda_vec(1:non_fuel_start))
    ! Calc beta correction per delay group
    
    if(t0==0.0) then
        beta_correction = (gen_time*(total_precursors_fuel))/total_power 
    end if

!---Hardcoded times to start perturbation - should read from input
    step_start_time = 0.0 
    step_end_time = 1.0 
    
    ramp_start_time = 0.0 
    ramp_end_time = 0.5

!---STEP perturbation
    if(step_flag .eqv. .TRUE.) then
        if(step_start_time < t0 .and. t0 < step_end_time ) then
            reactivity = reactivity_input 
        elseif ( t0 > step_end_time) then
            reactivity = reactivity_input 
        end if
    end if !---End STEP perturbation

!---RAMP perturbation 
    if(ramp_flag .eqv. .TRUE.) then
        if(ramp_start_time < t0 .and. t0 < ramp_end_time) then   
            rho_initial = 0.0 
            reactivity = rho_initial + &
              ((reactivity_input - rho_initial)*&
               (t0 - t_initial))/(ramp_end_time - t_initial)
        elseif ( t0 > ramp_end_time) then
            reactivity = 0.0 
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

    temp_reactivity_feedback = 0.0

!---Calculate temperature reactivity feedback
    if(feedback_method == 1 ) then
       call temperature_feedback(temp_reactivity_feedback,t0,nl_iter)
        reactivity_feedback = temp_reactivity_feedback
    end if

!---Power Solve
    if( mass_flow <= 0.0) then
        power_amplitude_new = power_amplitude_last_time + &
                          delta_t*( (reactivity - sum(beta_i_mat) &
                           )/gen_time)*power_amplitude_prev &
                          + delta_t*(1.0_dp/total_fuel_length)*&
                          total_precursors_fuel
    else 
        power_amplitude_new = power_amplitude_last_time + &
                          delta_t*( (reactivity  &
                          - beta_correction  )/gen_time)*power_amplitude_prev &
                          + delta_t*(1.0_dp/total_fuel_length)*&
                          total_precursors_fuel
    end if

!---Project power onto spatial shape
    power_soln_new(:,:) = 0.0
    do i = 1, non_fuel_start 
        do j = 1, nodes_per_elem
            power_soln_new(i,j) = power_amplitude_new*spatial_power_fcn(i,j)      
        end do
    end do

!---End power solve

end subroutine solve_power_backward_euler
