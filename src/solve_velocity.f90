!---Solve velocity field
!
!---Input: n - element number 
!
!---Output:
!
!

subroutine solve_velocity(n,nl_iter)

    USE flags_M
    USE material_info_M
    USE mesh_info_M
    USE time_info_M
    USE global_parameters_M
    USE solution_vectors_M
    USE element_matrices_M

implicit none

!---Dummy
    integer, intent(in) :: n 
    integer, intent(in) :: nl_iter

!---Local
    integer :: i, j
    real(dp) :: beta_1, beta_2, mass_flow_1, mass_flow_2, a,b,f
    real(dp) :: temperature_eval, density_eval
    real(dp) :: previous_time
    logical :: beta_change_flag
    real(dp) :: element_ss_density, element_prev_density

    temperature_eval = 0.0_dp
    !---Loop over all nodes in element
    do j = 1, nodes_per_elem 
        !---valuate density based on temperature
        temperature_eval = temperature_eval + vol_int(j)*temperature_soln_new(n,j)
    end do

    do j = 1, nodes_per_elem
        
        call density_corr(temperature_eval, density_eval)
        density_soln_new(n,j)  = density_eval
        
        !---New velocity
        velocity_soln_new(n,j) = (mass_flow / &
                                 (spatial_area_fcn(n,j)*density_eval))
       
        !Density_Reactivity_Feedback(n,j) = 0.0
     end do
  
    element_ss_density   = 0.0_dp
    element_prev_density = 0.0_dp
    
    if( t0 > 0.0) then
            do j = 1, nodes_per_elem
                element_ss_density = element_ss_density + vol_int(j)*density_soln_ss(n,j)
              element_prev_density = element_prev_density + vol_int(j)*density_soln_prev(n,j)
            end do

            !---Evaluate feedback based on density change
            if( n > Fuel_Inlet_Start .and. n < Fuel_Outlet_End ) then
                Density_Reactivity_Feedback(n) = (spatial_expansion_fcn(n,2) / &
                    (0.01_dp)) * &
                    (element_ss_density/element_prev_density - 1.0_dp)

            end if
        
        end if
end subroutine solve_velocity
