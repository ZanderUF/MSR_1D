module Mod_GlobalConstants
    
    use Mod_DelayPrecursorData
    
    implicit none
    
    type(delayPrecData), allocatable, dimension(:) :: allPrecursorData
    
end module Mod_GlobalConstants