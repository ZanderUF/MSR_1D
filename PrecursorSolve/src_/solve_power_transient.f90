!  
! Notes: Solves for the power amplitude
!  
! Input: nl_iter - current nonlinear iteration
!        current_time - current time we are solving at
!
! Output:
! 

subroutine solve_power_transient(nl_iter, current_time)

    USE parameters_fe

    implicit none

!---Dummy
    integer, intent(in) :: nl_iter
    double precision, intent(in) :: current_time

!---Local
    integer :: i,j,f,g
    real (kind=16), dimension(num_elem) :: temp_vec_num_elem, &
                                           precursors_lambda_vec, &
                                           power_soln_new_temp 
    real (kind=16):: power_new_total, total_precursor_ref,&
                     total_precursor_ref_sum, total_fuel_length,&
                     total_precursors_fuel, &
                     save_time_interval, rho_initial, step_time
    integer :: power_write_unit
    character(len=24) :: time_soln_name
    character(len=10) :: time_characters
    real(kind=4) :: temp_time
    real :: ramp_time

!---Initialize to zero
    precursors_lambda_vec(:) = 0.0
    temp_vec_num_elem(:) = 0.0
    total_precursor_ref = 0.0

!---Calculate total precursor concentration*lamda over system 
    do f = 1, num_isotopes 
       do g = 1, num_delay_group
            do i = 1, num_elem 
                do j = 1, nodes_per_elem
                   precursors_lambda_vec(i) = precursors_lambda_vec(i) + &
                                      lamda_i_mat(f,g)*&
                                      elem_vol_int(i,j)*&
                                      precursor_soln_prev(f,g,i,j)
                   
                   total_precursor_ref = total_precursor_ref + &
                                      lamda_i_mat(f,g)*&
                                      elem_vol_int(i,j)*&
                                      precursor_soln_prev(f,g,i,j)
                end do
            end do
       end do
    end do
    
    !---Get total length of the fuel element
    total_fuel_length = 0.0
    do i = 1, num_elem
        do j = 1, nodes_per_elem
            total_fuel_length = total_fuel_length + &
                                spatial_power_fcn(i,j)*elem_vol_int(i,j)
        end do
    end do
   
    total_precursor_ref_sum   = sum(precursors_lambda_vec)
    total_precursors_fuel     = sum(precursors_lambda_vec(1:non_fuel_start))
    beta_correction           = gen_time*((total_precursor_ref - &
                                total_precursors_fuel) / &
                                (power_amplitude_start*total_fuel_length))

!---Hardcoded times to start perturbation - should read from input
    step_time = 0.5
    ramp_time = 0.2 

!---STEP perturbation
    if(step_flag .eqv. .TRUE.) then
        if(t0 > step_time) then
            !write(outfile_unit, ('(a,12es14.3)')), 'STEP perturbation at: ',t0
            reactivity = reactivity_input 
        end if
    end if !---End STEP perturbation

!---RAMP perturbation 
    if(ramp_flag .eqv. .TRUE.) then
        if(t0 < ramp_time) then   
            !write(outfile_unit, ('(a,12es14.3)')),&
            !     'RAMP perturbation at time: ', t0
            rho_initial = 0.0 
            reactivity = rho_initial + &
              ((reactivity_input - rho_initial)*&
               (t0-t_initial))/(ramp_time - t_initial)
        elseif ( t0 > ramp_time) then
            reactivity = reactivity_input
        end if
    end if !---End RAMP perturbation 

!---Power Solve
    power_amplitude_new = power_amplitude_prev + &
                          delta_t*( (reactivity - ( (sum(beta_i_mat) &
                          - beta_correction) ) )/gen_time)*power_amplitude_prev &
                          + delta_t*(1.0_dp/total_fuel_length)*&
                          total_precursors_fuel

!---Write power amp out @ every time step
    if(t0 == 0.0) then
        write(power_outfile_unit, ('(a)')), 'Time (s) | Power Amp | Reactivity'
    end if

    write(power_outfile_unit, ('(12es14.6,12es14.5,12es14.5)')), &
          t0,power_amplitude_new,reactivity

!---Project power onto spatial shape
    power_soln_new(:,:) = 0.0
    do i = 1, non_fuel_start 
        do j = 1, nodes_per_elem
            power_soln_new(i,j) = power_amplitude_new*spatial_power_fcn(i,j)      
        end do
    end do

!---Write out power solution 
    save_time_interval = 0.25 
    if( modulo(t0,save_time_interval) < delta_t) then
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

!---End power solve

end subroutine solve_power_transient
