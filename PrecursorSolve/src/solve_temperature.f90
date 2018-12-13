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
    real(dp) :: temperature_last, power_current
    real(dp) :: temperature_current
!---Loop over all nodes in an element
    do j = 1, nodes_per_elem
        !---Get thermal heat_capacity value
        temperature_eval = temperature_soln_prev(n,j)
        call heat_capacity_corr(temperature_eval, heat_capacity_eval)
       
    end do

    temperature_last = 0.0_dp
    power_current    = 0.0_dp
    !---Get last temperature value
    if( n < 1) then
        do j = 1, nodes_per_elem
            temperature_last = temperature_last + &
                               temperature_soln_new(num_elem,j)*vol_int(j)
            power_current = power_current + power_soln_new(num_elem,j)*vol_int(j)
        end do
    else
        do j = 1, nodes_per_elem
            temperature_last = temperature_last + temperature_soln_new(n-1,j)*vol_int(j)
            power_current = power_current + power_soln_new(n-1,j)*vol_int(j)
        end do
    end if

    call heat_capacity_corr(temperature_last, heat_capacity_eval)

    temperature_current = temperature_last + &
                          power_current/(mass_flow*heat_capacity_eval)
    
    do j = 1, nodes_per_elem
        temperature_soln_new(n,j) = temperature_current*vol_int(j)
    end do

   !!---This should work for forward Euler
   !temperature_soln_new(n,j) = ( (power_soln_prev(n,j) - &
   !                               power_soln_starting(n,j) ) &
   !                            /( (mass_flow*heat_capacity_eval) )) &
   !                            + temperature_soln_prev(n,j) 
    
    
    do j = 1, nodes_per_elem        
        if(n < Fuel_Inlet_Start) then
            temperature_soln_new(n,j) = temperature_soln_new(num_elem,3)
        end if
        
        if( n <= Heat_Exchanger_Start  .and.  Fuel_Core_End < n) then
            temperature_soln_new(n,j) = temperature_soln_new(Fuel_Core_End,3) 
        end if

        !---Start heat exchanger 
        if ( Heat_Exchanger_Start <= n  .and.  n < Heat_Exchanger_End) then
            temperature_soln_new(n,j) = (100.0_dp)/&
            (Heat_Exchanger_End - Heat_Exchanger_Start)* &
            (global_coord(Heat_Exchanger_End-1,3)  - global_coord(n,j) ) - 100 + &  
            temperature_soln_new(Heat_Exchanger_Start-1,3)

        !---Rest of the external piping
        else if ( n >= Heat_Exchanger_End) then 
            temperature_soln_new(n,j) = temperature_soln_new(Heat_Exchanger_End-1,3) 
        end if
 
        !---Calculate doppler feedback for current element 
        Temperature_Reactivity_Feedback(n,j) = spatial_doppler_fcn(n,j)* &
                    ( temperature_soln_prev(n,j) - temperature_soln_starting(n,j) )
    end do

end subroutine solve_temperature
