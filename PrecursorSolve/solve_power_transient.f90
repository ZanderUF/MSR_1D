!  
! Notes: Solves for the total power in the system
!  
! Input: nl_iter - current nonlinear iteration
! 
! Output:
! 

subroutine solve_power_transient(nl_iter)

    USE parameters_fe

    implicit none

!---Dummy
    integer,intent(in) :: nl_iter

!---Local
    integer :: i,j,f,g
    real, dimension(num_elem) :: temp_vec_num_elem, temp_vec_prec
    real :: power_new_total, power_prev_total,total_precursor_ref,&
            beta_correction, total_precursors_fuel
    real, dimension(num_elem) :: power_soln_new_temp 
    
!---Initialize to zero
    temp_vec_prec(:) = 0.0
    temp_vec_num_elem(:) = 0.0

!---Calculate total precursor concentration*lamda over system 
    do f = 1, num_isotopes 
       do g = 1, num_delay_group
            do i = 1, num_elem 
                do j = 1, nodes_per_elem
                   temp_vec_prec(i) = temp_vec_prec(i) + &
                                      lamda_i_mat(f,g)*&
                                      elem_vol_int(i,j)*&
                                      precursor_soln_prev(f,g,i,j)
                end do
            end do
       end do
    end do
!---Calculate power over the active fuel region 
    do i = 1,non_fuel_start 
        do j = 1, nodes_per_elem
           temp_vec_num_elem(i) = temp_vec_num_elem(i) + &
                                  elem_vol_int(i,j)*power_soln_new(i,j) 
        end do
    end do
    
    power_prev_total = sum(temp_vec_num_elem)
    
    total_precursor_ref   = sum(temp_vec_prec)
    total_precursors_fuel = sum(temp_vec_prec(1:non_fuel_start))
    beta_correction = gen_time*((total_precursor_ref - &
                      total_precursors_fuel)/power_prev_total)
    
    !print *,'tot prec', ((sum(beta_i_mat)-beta_correction)/gen_time)*power_prev_total

!---Power Solve
    power_new_total = power_prev_total + delta_t*( (reactivity - (sum(beta_i_mat) &
                    - beta_correction))/gen_time)*power_prev_total + &
                      delta_t*total_precursors_fuel
!---Project power onto spatial shape
    power_soln_new_temp(:) = 0.0
    !print *,'power_new_total per mesh', power_new_total/(non_fuel_start)
    do i = 1, non_fuel_start 
        do j = 1, nodes_per_elem
            power_soln_new(i,j) = (power_new_total/(non_fuel_start))*elem_vol_int(i,j)      
        end do
    end do
    !print *,'power soln new',power_soln_new
    !print *,'sum soln',sum(power_soln_new)

!---End power solve

 end subroutine solve_power_transient
