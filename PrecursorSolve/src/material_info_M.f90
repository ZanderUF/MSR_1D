! Stores information about materials and point kinetics parameters
!

module material_info_M

    USE global_parameters_M

implicit none

!---Point kinetics parameters 
    real(dp), allocatable  :: lamda_i_mat(:,:)
    real(dp), allocatable  :: beta_i_mat(:,:)
    real(dp) ::               beta_correction
    real(dp) ::               gen_time = 0.0_dp
    
    real(dp) :: reactivity          = 0.0_dp
    real(dp) :: reactivity_input    = 0.0_dp
    real(dp) :: reactivity_feedback = 0.0_dp
    real(dp) :: temp_reactivity_feedback = 0.0_dp
    real(dp) :: Density_Reactivity_Feedback = 0.0_dp
!---Material information
    real(dp) ::               mass_flow
    integer  ::               num_isotopes
    integer  ::               num_delay_group

end module material_info_M
