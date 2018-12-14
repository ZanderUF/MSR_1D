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
    real(dp) :: previous_time, time_constant
    logical :: beta_change_flag

    !---Loop over all nodes in element
    do j = 1, nodes_per_elem 
        temperature_eval = temperature_soln_new(n,j)
        !---Evaluate density based on temperature
        call density_corr(temperature_eval, density_eval)
   
            density_soln_new(n,j)  = density_eval
            !---New velocity
            velocity_soln_new(n,j) = (mass_flow / &
                                 (spatial_area_fcn(n,j)*density_eval))
        
            !---Evaluate feedback based on density change
            !Density_Reactivity_Feedback(n,j) = (spatial_expansion_fcn(n,j) / &
            !            (density_soln_starting(n,j)*total_density_change)) * &
            !            ( density_soln_starting(n,j) - density_soln_new(n,j) )
            
            Density_Reactivity_Feedback(n,j) = 0.0
        end do
   
end subroutine solve_velocity
