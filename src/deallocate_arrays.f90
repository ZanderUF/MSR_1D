

subroutine deallocate_arrays ( )

    USE flags_M
    USE global_parameters_M
    USE datainput_fe_M
    USE mesh_info_M
    USE material_info_M
    USE time_info_M
    USE solution_vectors_M

implicit none

!---Dummy

!---Local

!---Allocate solution vector and global matrices
!   deallocate(precursor_soln_new, &
!             power_soln_new      (num_elem,nodes_per_elem), &
!            temperature_soln_new( num_elem,nodes_per_elem),  &
!            density_soln_new( num_elem,nodes_per_elem),  &
!            velocity_soln_new( num_elem,nodes_per_elem),  &
!            precursor_soln_prev(num_isotopes,num_delay_group,num_elem,nodes_per_elem), &
!            power_soln_prev(num_elem,nodes_per_elem), &
!            temperature_soln_prev( num_elem,nodes_per_elem), &
!            density_soln_prev( num_elem,nodes_per_elem),  &
!            velocity_soln_prev( num_elem,nodes_per_elem), &
!            spatial_vol_fcn(num_elem,nodes_per_elem), &
!            spatial_area_fcn(num_elem,nodes_per_elem),&
!            spatial_power_frac_fcn(num_elem,nodes_per_elem), &
!            spatial_power_fcn( num_elem, nodes_per_elem), &
!            spatial_doppler_fcn(num_elem, nodes_per_elem), &
!            spatial_expansion_fcn(num_elem, nodes_per_elem), &
!            elem_vec_q_final(num_isotopes,num_delay_group,nodes_per_elem),& 
!            elem_vol_int(num_elem,nodes_per_elem),&
!            precursor_soln_last_time(num_isotopes,num_delay_group,num_elem,nodes_per_elem),&
!            power_soln_last_time(num_elem,nodes_per_elem),& 
!            area_variation(num_elem,nodes_per_elem), & 
!            density_soln_starting(num_elem,nodes_per_elem), &
!            temperature_soln_starting(num_elem,nodes_per_elem), &
!            power_soln_starting(num_elem,nodes_per_elem), &
!            Density_Reactivity_Feedback(num_elem), &
!            Temperature_Reactivity_Feedback(num_elem,nodes_per_elem),&
!            density_soln_ss(num_elem,nodes_per_elem),&
!            temperature_soln_ss(num_elem,nodes_per_elem) )
!
!llocate(beta_initial_vec(num_isotopes,num_delay_group))
!llocate(beta_j(num_isotopes,num_delay_group))
!llocate(beta_j_minus_1(num_isotopes,num_delay_group) )
!llocate(beta_j_plus_1(num_isotopes,num_delay_group) )
!llocate(beta_change_all_previous(num_isotopes,num_delay_group))
!llocate(beta_change(num_isotopes,num_delay_group))
!llocate(beta_correction_vec(num_isotopes,num_delay_group))
!llocate(residual(num_isotopes,num_delay_group,num_elem,nodes_per_elem))
!
!llocate(L2_prev_precursors(num_isotopes,num_delay_group))
!llocate(L2_current_precursors(num_isotopes,num_delay_group))
!llocate(L2_diffs_precursors(num_isotopes,num_delay_group))


end subroutine deallocate_arrays 
