!---Solve velocity field
!
!---Input: n - element number 
!
!---Output:
!
!

subroutine solve_velocity(n)

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

!---Local
    integer :: j
    real(dp) :: temperature_eval, density_eval

!---Loop over all nodes in element
    do j = 1, nodes_per_elem
        temperature_eval       = temperature_soln_prev(n,j)
        !---Evaluate density based on temperature
        call density_corr_msfr(temperature_eval, density_eval)
        density_soln_new(n,j)  = density_eval
        
        !velocity_soln_new(n,j) = mass_flow*exp(-0.2_dp*t0) / &
        !                         (area_variation(n,j)*density_eval)
        velocity_soln_new(n,j)  = mass_flow / &
                                 (area_variation(n,j)*density_eval)
        

    !---Evaluate feedback based on density change
        Density_Reactivity_Feedback(n,j) = (spatial_expansion_fcn(n,j) / &
                        (density_soln_starting(n,j)*total_density_change)) * &
                        (density_soln_starting(n,j) - density_soln_new(n,j) )

    end do

end subroutine solve_velocity
