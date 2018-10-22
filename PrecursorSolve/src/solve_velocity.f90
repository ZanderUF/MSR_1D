! Solve velocity field
!
! Input: n - element number 
!
! Output:
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
    real :: temperature, density

!---Loop over all nodes in element
    do j = 1, nodes_per_elem
        temperature = temperature_soln_prev(n,j)
        !---Evaluate density based on temperature
        call density_corr(temperature,density)
        density_soln_new(n,j) = density
        velocity_soln_new(n,j) = mass_flow/(area_variation(n,j)*density)
    end do

end subroutine solve_velocity
