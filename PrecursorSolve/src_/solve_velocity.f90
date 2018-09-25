! Solve velocity field
!
! Input: n - element number 
!
! Output:
!
!

subroutine solve_velocity(n)

    USE parameters_fe

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
        velocity_soln_new(n,j) = mass_flow/(area*density)
    end do


end subroutine solve_velocity
