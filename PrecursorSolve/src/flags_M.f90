module flags_M

implicit none

!---Flags
    integer :: feedback_method = 0
    real :: save_time_interval
    integer :: td_method_type = 0
    logical :: DEBUG = .FALSE.
    logical :: steady_state_flag = .TRUE.    
    logical :: unit_test = .FALSE. 
    logical :: unit_test_2 = .FALSE. 
    logical :: nonlinear_ss_flag 
    logical :: lagrange = .FALSE.
    logical :: hermite = .TRUE.
    logical :: transient
    logical :: transient_save_flag = .FALSE. 
    logical :: ramp_flag = .FALSE.
    logical :: step_flag = .FALSE.
    logical :: zag_flag = .FALSE.
    logical :: read_power_from_file = .FALSE.

end module flags_M
