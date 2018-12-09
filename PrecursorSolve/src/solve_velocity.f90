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
    integer :: i, j
    real(dp) :: beta_1, beta_2, mass_flow_1, mass_flow_2, a,b,f
    real(dp) :: temperature_eval, density_eval
    real(dp) :: previous_time, time_constant
    logical :: beta_change_flag

    beta_change_flag = .TRUE. 

    !if(beta_change_flag .eqv. .TRUE.) then
    !!---Could have some logic to not call @ every time step
   
    !else
    !    !--- instant change 
    !    time_constant = -0.2_dp
    !    !---Evaluate pump coast down 
    !    mass_flow = mass_flow_initial*exp(time_constant*t0)
    !
    !!---Get Beta for the current mass flow
    !    do i = 1, number_entries_beta
    !        if( mass_flow <=  Beta_Fcn_Flow(i,1) ) then
    !            exit
    !        end if
    !    end do
  
    !    !---Interpolate to get beta at current mass flow
    !    beta_1 = Beta_Fcn_Flow(i-1,2)
    !    beta_2 = Beta_Fcn_Flow(i,  2)
    !    mass_flow_1 = Beta_Fcn_Flow(i-1,1)
    !    mass_flow_2 = Beta_Fcn_Flow(i,  1)
    !    a = mass_flow   - mass_flow_1
    !    b = mass_flow_2 - mass_flow
    !    f = a/(a+b)
    !    
    !    beta_correction = f*beta_2 + (1.0_dp - f)*beta_1
    !end if

    !---Loop over all nodes in element
    do j = 1, nodes_per_elem 
        temperature_eval       = temperature_soln_prev(n,j)
        !---Evaluate density based on temperature
        call density_corr(temperature_eval, density_eval)
   
            density_soln_new(n,j)  = density_eval
            !---New velocity
            velocity_soln_new(n,j) = (mass_flow / &
                                 (area_variation(n,j)*density_eval))
        
            !---Evaluate feedback based on density change
            Density_Reactivity_Feedback(n,j) = (spatial_expansion_fcn(n,j) / &
                        (density_soln_starting(n,j)*total_density_change)) * &
                        ( density_soln_starting(n,j) - density_soln_new(n,j) )
        end do
    
end subroutine solve_velocity
