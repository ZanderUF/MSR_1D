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
    real(dp) :: temperature_eval, density_eval
    real(dp) :: beta_1, beta_2, mass_flow_1, mass_flow_2, a,b,f
    
    
    !---Evaluate pump coast down 
    mass_flow = mass_flow_initial*exp(-0.2_dp*t0)
    !---Get Beta for the current mass flow
    do i = 1, number_entries_beta
        if( mass_flow <=  Beta_Fcn_Flow(i,1) ) then
            !print *,'mass flow ',mass_flow
            !print *,' beta fcn', Beta_Fcn_Flow(i,1)
            exit
        end if

    end do
    
    beta_1 = Beta_Fcn_Flow(i-1,2)
    beta_2 = Beta_Fcn_Flow(i,  2)
    mass_flow_1 = Beta_Fcn_Flow(i-1,1)
    mass_flow_2 = Beta_Fcn_Flow(i,  1)
    a = mass_flow - mass_flow_1
    b = mass_flow_2 - mass_flow
    f = a/(a+b)
    !print *,' beta 1 ' ,beta_1
    !print *,'beta 2 ' ,beta_2
    !print *,' mass flow 1', mass_flow_1
    !print *,' mass flow 2', mass_flow_2
    !print *,' f' , f

    beta_correction = beta_1**f * beta_2**(1-f)
    ! Use fit line from excel
    !beta_correction = -1.04E-4_dp*log(mass_flow) + 0.0082_dp
   
    !print *,'beta corre', beta_correction
    !print *,' ' 

!---Loop over all nodes in element
    do j = 1, nodes_per_elem
        temperature_eval       = temperature_soln_prev(n,j)
        !---Evaluate density based on temperature
        call density_corr(temperature_eval, density_eval)
        density_soln_new(n,j)  = density_eval
        
        
        ! Exponentially decrease the mass flow rate
        velocity_soln_new(n,j) = mass_flow / &
                                 (area_variation(n,j)*density_eval)
        
        !velocity_soln_new(n,j)  = mass_flow / &
        !                         (area_variation(n,j)*density_eval)
        

    !---Evaluate feedback based on density change
        Density_Reactivity_Feedback(n,j) = (spatial_expansion_fcn(n,j) / &
                        (density_soln_starting(n,j)*total_density_change)) * &
                        ( density_soln_new(n,j) - density_soln_starting(n,j) )

    end do

end subroutine solve_velocity
