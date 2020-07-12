module Mod_GlobalConstants
    
    use Mod_DelayPrecursorData
    use global_parameters_M, only: dp
    
    implicit none
    
    type(delayPrecData), allocatable, dimension(:) :: allPrecursorData
    
    !---Flags
    logical :: flow_perturbation = .TRUE.
    
    real(dp) :: max_delta_t
    
end module Mod_GlobalConstants