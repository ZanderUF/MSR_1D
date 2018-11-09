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

subroutine temperature_feedback(temp_reactivity_feedback,current_time,nl_iter)

    USE global_parameters_M
    USE solution_vectors_M
    USE mesh_info_M

    implicit none

!---Dummy
    real, intent(out) :: temp_reactivity_feedback
    integer, intent(in) :: current_time 
    integer, intent(in) :: nl_iter

!---Local
    real :: temperature_coefficient
    real ::  avg_temperature_current
    real :: total_temperature_current
    integer :: fuel_elem_len, i, j

    !---Temp Coefficient
    temperature_coefficient = 3.22E-5 

    !---Initialize
    total_temperature_current  = 0.0

    do i = Fuel_Inlet_Start, Fuel_Outlet_End 
        do j = 1, nodes_per_elem
            total_temperature_current = total_temperature_current + temperature_soln_new(i,j)
        end do
    end do


    avg_temperature_current  = total_temperature_current/fuel_elem_len
    temp_reactivity_feedback = temperature_coefficient * &
                               (avg_temperature_current - avg_temperature_initial)    

end subroutine temperature_feedback 
