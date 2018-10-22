! Solves for the temperature at a time step
!
! Input: n - element number
! 
! Output: 
!

subroutine solve_temperature(n)

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
    real(dp) :: temperature_eval, heat_capacity_eval

!---Loop over all nodes in an element
    do j = 1, nodes_per_elem
        !---Get thermal heat_capacity value
        temperature_eval = temperature_soln_prev(n,j)
        call heat_capacity_corr(temperature_eval, heat_capacity_eval)
        !call cond_corr(temperature,heat_capacity)
        !---This should work for forward Euler
        temperature_soln_new(n,j) = (( (total_power_initial/sum(spatial_power_fcn) )*&   
                                     power_soln_prev(n,j) ) &
                                    /( (mass_flow*heat_capacity_eval) )) &
                                    + temperature_soln_prev(n,j) 
    end do
    
end subroutine solve_temperature
