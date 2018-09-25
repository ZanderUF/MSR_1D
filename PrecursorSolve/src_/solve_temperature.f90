! Solves for the temperature at a time step
!
! Input: n - element number
! 
! Output: 
!

subroutine solve_temperature(n)

    USE parameters_fe

implicit none

!---Dummy
    integer, intent(in) :: n 

!---Local
    integer :: j
    real :: temperature, heat_capacity

!---Loop over all nodes in an element
    do j = 1, nodes_per_elem
        !---Get thermal heat_capacity value
        temperature = temperature_soln_prev(n,j)
        call heat_capacity_corr(temperature, heat_capacity)
        !call cond_corr(temperature,heat_capacity)
        !---This should work for forward Euler
        temperature_soln_new(n,j) = (power_soln_prev(n,j)/(mass_flow*heat_capacity) )&
                                    + temperature_soln_prev(n,j) 
    end do
    
end subroutine solve_temperature
