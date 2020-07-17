!*****************************************************************************80
! Calculates feedback spatially dependent reactivity coeffecient 
! for doppler and fuel expansion 
! 
!
! Input: n - node number, nl_iter - nonlinear iteration, current_time - time 
! 
! Output: temp_reactivity_feedback - spatial value of temperature feedback
!         Density_reactivity_feedback - spatial density feedback value
! 
!*****************************************************************************80

subroutine temperature_feedback(n, nl_iter, current_time, Temp_Reactivity_Feedback)

    USE global_parameters_M
    USE solution_vectors_M
    USE mesh_info_M

    implicit none

!---Dummy
    integer, intent(in) :: n
    integer, intent(in) :: nl_iter
    integer, intent(in) :: current_time 
    real(dp), intent(out) :: Temp_Reactivity_Feedback
    
!---Local
    real(dp) :: Temperature_Coefficient_Reactivity
    real(dp) :: avg_temperature_current
    real(dp) :: total_temperature_current
    real(dp) :: fuel_elem_len
    integer  :: i, j

    !---Temp Coefficient
    Temperature_Coefficient_Reactivity = -1.15796E-3_dp/407.0_dp 
      
    !---Initialize
    total_temperature_current  = 0.0
    !---Calculate current average temperature across fuel
    do i = Fuel_Inlet_Start, Fuel_Outlet_End 
        do j = 1, nodes_per_elem
            total_temperature_current = total_temperature_current + &
                                       elem_vol_int_fe(j)*temperature_soln_new(i,j)
        end do
    end do
    
    !---Total Length of the fuel
    fuel_elem_len = global_coord(Fuel_Outlet_End,3) - global_coord(Fuel_Inlet_Start,1) 

    avg_temperature_current  = total_temperature_current/fuel_elem_len

    Temp_Reactivity_Feedback = Temperature_Coefficient_Reactivity * &
                               (avg_temperature_current - avg_temperature_initial)    

end subroutine temperature_feedback

!---Calculate density feedback <==> fuel expansion
!
!
subroutine density_feedback( n, j, current_time, Den_Feedback )

    USE global_parameters_M
    USE solution_vectors_M
    USE mesh_info_M

    implicit none

!---Dummy
    integer, intent(in)   :: n
    integer, intent(in)   :: j 
    integer, intent(in)   :: current_time 
    real(dp), intent(inout) :: Den_Feedback 
    
!---Local
    real(dp) :: Temperature_Coefficient_Reactivity
    real(dp) :: Density_Coefficient_Reactivity
    integer  :: i  
    real(dp) :: Total_Density_Current
    
    !---Density Coeffecien
    Density_Coefficient_Reactivity = 2.94877E-3_dp/&
                    (0.99_dp*Avg_Density_Initial)  
      
    !---Initialize
    total_density_current = 0.0
    
    !Den_Feedback = Density_Coefficient_Reactivity * &
    !                           (Avg_Density_Current - Avg_Density_Initial)
                            
end subroutine density_feedback
