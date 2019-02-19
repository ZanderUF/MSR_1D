module time_info_M 

    USE global_parameters_M

    implicit none

!---Time information 
    logical :: time_solve     ! decide if we are doing time solve or not  
    real(dp)   alpha          ! time solve 
    real(dp)   t0             ! starting time
    real(dp)   delta_t        ! time step 
    real(dp)   tmax           ! max time 
    real(dp)   t_initial      ! starting time
    real(dp)   ramp_start_time
    real(dp)   ramp_end_time
    real(dp)   step_start_time
    real(dp)  step_end_time
!---Time constant for flow slowdown
    real(dp)  time_constant
    real(dp) flow_reduction_percent

end module time_info_M
