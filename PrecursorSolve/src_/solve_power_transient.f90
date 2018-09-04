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
    real, intent(in) :: current_time

!---Local
    integer :: i,j,f,g
    real (kind=16), dimension(num_elem) :: temp_vec_num_elem, precursors_lambda_vec
    real (kind=16):: power_new_total, total_precursor_ref,&
            total_precursors_fuel
    real(kind=16), dimension(num_elem) :: power_soln_new_temp 
    real (kind=16):: total_fuel_length
    real(kind=16) :: total_precursor_ref_sum
    real :: step_time
    
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
            total_fuel_length = total_fuel_length + spatial_power_fcn(i,j)*elem_vol_int(i,j)
        end do
    end do
   
    total_precursor_ref_sum   = sum(precursors_lambda_vec)
    total_precursors_fuel     = sum(precursors_lambda_vec(1:non_fuel_start))
    beta_correction           = gen_time*((total_precursor_ref - &
                            total_precursors_fuel)/(power_amplitude_start*total_fuel_length))
   
    !print *,'total precursors', total_precursor_ref
    !print *,'total precursors sum routing', total_precursor_ref_sum
    !print *,'total prec fuel transient', total_precursors_fuel
    !print *,'total beta', sum(beta_i_mat)
    !print *,'beta correction', beta_correction

    step_time = 0.1
!---STEP perturbation
    if(current_time > step_time) then
    	reactivity = 1E-4 
    end if

!---Power Solve
    power_amplitude_new = power_amplitude_prev + &
                          delta_t*( (reactivity - ( (sum(beta_i_mat) &
                          - beta_correction) ))/gen_time)*power_amplitude_prev + &
                          delta_t*(1.0_dp/total_fuel_length)*total_precursors_fuel
    !print *, ' '
    !print *,'sum_beta',sum(beta_i_mat)
    !print *,'beta cor',beta_correction
    !print *,'power_amp prev',power_amplitude_prev
    !print *,' ' 
    !print *,'power amplitude new ', power_amplitude_new
    !print *,'first part', delta_t*( (reactivity - ( (sum(beta_i_mat) &
    !                - beta_correction)))/gen_time)*power_amplitude_prev
    !print *,'second part', delta_t*(1.0/total_fuel_length)*total_precursors_fuel
    !print *,' ' 

!---Project power onto spatial shape
    power_soln_new(:,:) = 0.0
    do i = 1, non_fuel_start 
        do j = 1, nodes_per_elem
            power_soln_new(i,j) = power_amplitude_new*spatial_power_fcn(i,j)      
        end do
    end do

!---Write out power solution 
    !write(outfile_unit,fmt='(a)'), ' '
    !write(outfile_unit,fmt='(a,12es14.3)'), 'Power distribution at time:', current_time  
    !write(outfile_unit,fmt='(a)'), 'Position(x) Power [n/cm^3*s]'
    !do i = 1,num_elem
    !    do j = 1, nodes_per_elem
    !        write(outfile_unit, fmt='(f6.3, 12es14.3)')  global_coord(i,j), power_soln_new(i,j)
    !    end do
    !end do

!---End power solve

end subroutine solve_power_transient
