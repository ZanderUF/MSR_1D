

subroutine allocate_arrays ( )

    USE flags_M
    USE global_parameters_M
    USE mesh_info_M
    USE material_info_M
    USE time_info_M
    USE solution_vectors_M

implicit none

!---Dummy

!---Calculate total number oftime steps in the problem
    NumTimeSteps = tmax*(1.0_dp/delta_t) 
!---Allocate arrays that contain all time information
    allocate(TimeAll(NumTimeSteps),&
             PowerAllTime(NumTimeSteps),&
             RhoInsertedAll(NumTimeSteps),&
             DopplerFeedbackAll(NumTimeSteps),&
             DensityFeedbackAll(NumTimeSteps),&
             MassFlowAll(NumTimeSteps),&
             TempHeatExchangerAll(NumTimeSteps),&
             AvgTempAll(NumTimeSteps),& 
             PeakTempAll(NumTimeSteps),&
             NonLinearIterAll(NumTimeSteps)&
             )
         

!---Allocate solution vector and global matrices
    allocate(precursor_soln_new(num_isotopes,num_delay_group,num_elem,nodes_per_elem),  &
             power_soln_new(num_elem,nodes_per_elem), &
             temperature_soln_new( num_elem,nodes_per_elem),  &
             density_soln_new( num_elem,nodes_per_elem),  &
             velocity_soln_new( num_elem,nodes_per_elem),  &
             precursor_soln_prev(num_isotopes,num_delay_group,num_elem,nodes_per_elem), &
             power_soln_prev(num_elem,nodes_per_elem), &
             temperature_soln_prev( num_elem,nodes_per_elem), &
             density_soln_prev( num_elem,nodes_per_elem),  &
             velocity_soln_prev( num_elem,nodes_per_elem), &
             spatial_vol_fcn(num_elem,nodes_per_elem), &
             spatial_area_fcn(num_elem,nodes_per_elem),&
             spatial_power_frac_fcn(num_elem,nodes_per_elem), &
             spatial_power_fcn( num_elem, nodes_per_elem), &
             spatial_doppler_fcn(num_elem, nodes_per_elem), &
             spatial_expansion_fcn(num_elem, nodes_per_elem), &
             elem_vec_q_final(num_isotopes,num_delay_group,nodes_per_elem),& 
             elem_vol_int(num_elem,nodes_per_elem),&
             precursor_soln_last_time(num_isotopes,num_delay_group,num_elem,nodes_per_elem),&
             power_soln_last_time(num_elem,nodes_per_elem),& 
             area_variation(num_elem,nodes_per_elem), & 
             density_soln_starting(num_elem,nodes_per_elem), &
             temperature_soln_starting(num_elem,nodes_per_elem), &
             power_soln_starting(num_elem,nodes_per_elem), &
             Density_Reactivity_Feedback(num_elem), &
             Temperature_Reactivity_Feedback(num_elem,nodes_per_elem),&
             density_soln_ss(num_elem,nodes_per_elem),&
             temperature_soln_ss(num_elem,nodes_per_elem) )

allocate(beta_initial_vec(num_isotopes,num_delay_group))
allocate(beta_j(num_isotopes,num_delay_group))
allocate(beta_j_minus_1(num_isotopes,num_delay_group) )
allocate(beta_j_plus_1(num_isotopes,num_delay_group) )
allocate(beta_change_all_previous(num_isotopes,num_delay_group))
allocate(beta_change(num_isotopes,num_delay_group))
allocate(beta_correction_vec(num_isotopes,num_delay_group))
allocate(residual(num_isotopes,num_delay_group,num_elem,nodes_per_elem))

allocate(L2_prev_precursors(num_isotopes,num_delay_group))
allocate(L2_current_precursors(num_isotopes,num_delay_group))
allocate(L2_diffs_precursors(num_isotopes,num_delay_group))


end subroutine allocate_arrays 
