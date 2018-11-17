! Calculates feedback reactivity based on the temperature difference 
! 
! reactivity = temp_coeff * (T_current - T_initial)
! 
! Notes: Temperatures should be weighted average over the fuel region
!        T_current = Sum over [Spatial_shape(x)*T(x)]
!
! reactivity_coeff = 15.3 pcm/*C

! Input:
! 
! Output:
! 

subroutine temperature_feedback(Temp_Reactivity_Feedback, Density_Reactivity_Feedback,&
                               current_time, nl_iter)

    USE global_parameters_M
    USE solution_vectors_M
    USE mesh_info_M

    implicit none

!---Dummy
    real(dp), intent(out) :: Temp_Reactivity_Feedback
    real(dp), intent(out) :: Density_Reactivity_Feedback 
    integer, intent(in) :: current_time 
    integer, intent(in) :: nl_iter

!---Local
    real(dp) :: Temperature_Coefficient_Reactivity
    real(dp) :: Density_Coefficient_Reactivity
    real(dp) :: Avg_Density_Current, avg_temperature_current
    real(dp) :: Total_Density_Current, total_temperature_current
    real(dp) :: fuel_elem_len
    integer  :: i, j

    !---Temp Coefficient
    ! -- 
    !Temperature_Coefficient_Reactivity = -2.58908E-4_dp/407.0_dp 
    Temperature_Coefficient_Reactivity = -1.15796E-3_dp/407.0_dp 
      
    !---Density Coeffecien
    Density_Coefficient_Reactivity = 2.94877E-3_dp/&
                    (0.99_dp*Avg_Density_Initial)  
    !---Initialize
    total_temperature_current  = 0.0
    total_density_current = 0.0
    !---Calculate current average temperature across fuel
    do i = Fuel_Inlet_Start, Fuel_Outlet_End 
        do j = 1, nodes_per_elem
            total_temperature_current = total_temperature_current + &
                                       elem_vol_int_fe(j)*temperature_soln_new(i,j)
            total_density_current = total_density_current + &
                                       elem_vol_int_fe(j)*density_soln_new(i,j)
            !print *,'density ', density_soln_new(i,j)
        end do
    end do
    
    !---Total Length of the fuel
    fuel_elem_len = global_coord(Fuel_Outlet_End,3) - global_coord(Fuel_Inlet_Start,1) 

    avg_temperature_current  = total_temperature_current/fuel_elem_len
    Avg_Density_Current      = total_density_current/fuel_elem_len

    Temp_Reactivity_Feedback = Temperature_Coefficient_Reactivity * &
                               (avg_temperature_current - avg_temperature_initial)    


    Density_Reactivity_Feedback = Density_Coefficient_Reactivity * &
                               (Avg_Density_Current - Avg_Density_Initial)
                            
end subroutine temperature_feedback 
