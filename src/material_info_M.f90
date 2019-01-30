! Stores information about materials and point kinetics parameters
!

module material_info_M

    USE global_parameters_M

implicit none

!---Point kinetics parameters 
    real(dp)               :: Long_Decay_Constant
    real(dp), allocatable  :: lamda_i_mat(:,:)
    real(dp), allocatable  :: beta_i_mat(:,:)
    real(dp), allocatable  :: beta_initial_vec(:,:)
     
    real(dp) ::               beta_correction
    real(dp) , allocatable :: beta_correction_vec(:,:)
    real(dp) ::               gen_time = 0.0_dp
    real(dp) , allocatable :: beta_j(:,:), beta_j_plus_1(:,:), beta_j_minus_1(:,:)
    real(dp) , allocatable :: beta_change_all_previous(:,:),beta_change(:,:) 
    real(dp) :: reactivity          = 0.0_dp
    real(dp) :: reactivity_input    = 0.0_dp
    real(dp) :: reactivity_feedback = 0.0_dp
    
    real(dp) :: total_temperature_feedback = 0.0_dp
    real(dp) :: total_density_feedback     = 0.0_dp
    real(dp) :: total_doppler_read_in, total_expansion_read_in, &
                total_temperature_change, total_density_change

!---Material information
    real(dp) ::               mass_flow
    real(dp) ::               mass_flow_initial
    integer  ::               num_isotopes
    integer  ::               num_delay_group
    
    real(dp), allocatable :: Beta_Fcn_Flow(:,:)
    integer :: number_entries_beta

    real(dp) :: temperature_reduction

end module material_info_M
