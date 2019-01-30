
!----Evaluate change in beta due to changes in flow speed

subroutine evaluate_beta_change(event_time, event_time_previous, &
                                event_counter, event_occuring)

    USE flags_M
    USE material_info_M
    USE mesh_info_M
    USE time_info_M
    USE global_parameters_M
    USE solution_vectors_M
    USE element_matrices_M
    USE time_info_M

    implicit none

    !---Dummy
    real(dp), intent(in)    :: event_time
    real(dp), intent(in)    :: event_time_previous
    integer,  intent(inout) :: event_counter
    logical,  intent(in)    :: event_occuring
    !---Local
    real(dp) :: beta_1, beta_2, mass_flow_1, mass_flow_2, a,b,mid
    real(dp),dimension(num_isotopes,num_delay_group) :: beta_interp_current,& 
                      beta_change_all_current 
    real(dp) :: current_time, number_half_lifes
    integer  :: i, g,f 
    real(dp), dimension(num_isotopes,num_delay_group) :: time_cutoff
    real(dp) :: Large_Lambda
    integer :: event_counter_last
     
    number_half_lifes = log(2.0_dp) 
  
    !---Get cutoff time to stop looking at differences between 
    !   beta values.  10*half_life of each group.  
    do f = 1, num_isotopes
        do g = 1, num_delay_group
            time_cutoff(f,g) = 10.0*(log(2.0_dp)/lamda_i_mat(f,g))
        end do
    end do
     
    !---Get Beta for the current mass flow
    do i = 1, number_entries_beta
        if( mass_flow >=  Beta_Fcn_Flow(i,1) ) then
            exit
        end if
    end do
        beta_counter = i
        
    !---Interpolate to get beta for each group for the mass flow
     do f = 1, num_isotopes
        do g = 1, num_delay_group
            if(i == 1) then

                mass_flow_1 = Beta_Fcn_Flow(i,1)
                mass_flow_2 = Beta_Fcn_Flow(i+1,1)
                
                beta_1 = Beta_Fcn_Flow(i,g+1)
                beta_2 = Beta_Fcn_Flow(i+1,  g+1)

            else

                mass_flow_2 = Beta_Fcn_Flow(i-1,  1)
                mass_flow_1 = Beta_Fcn_Flow(i,1)
                
                beta_2 = Beta_Fcn_Flow(i,g+1)
                beta_1 = Beta_Fcn_Flow(i-1,  g+1)

            end if
                        
            a = mass_flow   - mass_flow_1
            b = mass_flow_2 - mass_flow
            mid = a/(a+b)
            
            beta_interp_current(f,g) = mid*beta_2 + (1.0_dp - mid)*beta_1
        end do
    end do
    
    !---Do instant change in beta as function of flow speed
    if(feedback_method==2) then
        do f = 1, num_isotopes
            do g = 1, num_delay_group
                beta_correction_vec(f,g) = beta_interp_current(f,g)
            end do
        end do
    end if

   current_time  = t0

    !---If very first event, initialize based on starting beta
    if(event_counter == 1) then
        
        if(current_time ==  event_time ) then
            do f = 1, num_isotopes 
                do g = 1, num_delay_group
                    beta_j_minus_1(f,g) = beta_initial_vec(f,g)
                    beta_j(f,g)         = beta_interp_current(f,g)
                    beta_change(f,g) = 0.0_dp
                    beta_change_all_previous(f,g) = 0.0_dp
                    beta_correction_vec(f,g) = beta_initial_vec(f,g)
                end do
            end do
        end if
        
        do f = 1, num_isotopes
            do g = 1, num_delay_group
                
                beta_j(f,g)         = beta_interp_current(f,g) 
                
                Large_Lambda = lamda_i_mat(f,g)*number_half_lifes
               
                beta_change_all_current(f,g) = &
                    (1.0_dp - exp(-Large_Lambda*(current_time - event_time)))
                
                if( event_occuring .eqv. .TRUE.) then
                    
                    beta_change(f,g)   = beta_change(f,g) + &
                     ( beta_j(f,g) - beta_initial_vec(f,g) )*&
                     ( beta_change_all_current(f,g)  ) 
                     !+ &
                     !( beta_j_minus_1(f,g) - beta_initial_vec(f,g))*&
                     !beta_change_all_previous(f,g)
                
                else
                    beta_change(f,g)   = beta_change(f,g) + &
                     ( beta_j(f,g) - beta_initial_vec(f,g) )*&
                     ( beta_change_all_current(f,g) - beta_change_all_previous(f,g) )
                end if
                
                !---Store correction per delay group and isotope 
                beta_correction_vec(f,g) = beta_initial_vec(f,g) + beta_change(f,g) 
                !---Keep a running sum of the changes 
                beta_change_all_previous(f,g) = &
                            beta_change_all_current(f,g)
                
                beta_j_minus_1(f,g) = beta_initial_vec(f,g)

            end do
        end do
        !---Count number of events
        !event_counter_last = event_counter
        !event_counter = event_counter + 1
    end if

    

    if(1==2) then !---Every event after the 'first'
        do f = 1, num_isotopes
            do g = 1, num_delay_group
                
                beta_j(f,g) = beta_interp_current(f,g)

                Large_Lambda = lamda_i_mat(f,g)*number_half_lifes
               
                beta_change_all_current(f,g) = &
                    (1.0_dp - exp(-Large_Lambda*(current_time - event_time)))
                
                beta_change(f,g)  = beta_change(f,g) + &
                      (beta_j(f,g) - beta_j_minus_1(f,g))*&
                      (beta_change_all_current(f,g) - beta_change_all_previous(f,g)) 
               
                beta_correction_vec(f,g) = beta_correction_vec(f,g) + &
                                           beta_change(f,g)
                
                beta_change_all_previous(f,g) = &
                            beta_change_all_previous(f,g)
                
                !beta_change_all_current(f,g) = &
                !    ( exp(Large_Lambda*delta_t) - 1.0_dp)* &
                !    ( exp(-Large_Lambda*(current_time - event_time))  )*&
                !    ( ( beta_j(f,g) - beta_j_minus_1(f,g)  )*&
                !     exp(-Large_Lambda*(event_time_previous-event_time)) + &
                !      (beta_j_plus_1(f,g) - beta_j(f,g) ) )
                
                !---Store correction per delay group and isotope 
                !beta_correction_vec(f,g) = beta_initial_vec(f,g) + &
                !                           beta_change_all_current(f,g) + &
                !                           beta_change_all_previous(f,g)

                !---Keep a running sum of the changes 
                !beta_change_all_previous(f,g) = &
                !            beta_change_all_previous(f,g) + &
                !            beta_change_all_current(f,g)
                
                !---Shuffle betas for next change
                beta_j_minus_1(f,g) = beta_j(f,g)
                !beta_j(f,g)         = beta_j_plus_1(f,g)

               end do
        end do
    end if

    beta_correction = sum(beta_correction_vec)

    !---Score an 'event'

end subroutine evaluate_beta_change
