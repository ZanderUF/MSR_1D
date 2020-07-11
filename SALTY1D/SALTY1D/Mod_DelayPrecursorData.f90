!--------------------------------------------------------------------------
!> @details Defines type/class for delayed neutron precursor data
!
!--------------------------------------------------------------------------
module Mod_DelayPrecursorData

    use global_parameters_M, only : sp, dp
    
    implicit none

    private 
    
    public :: delayPrecData
    
    !---Define the type
    type :: delayPrecData
        
        character(len=10) :: isotopeName 
        real(sp), allocatable, dimension(:) :: decayConst
        real(dp), allocatable, dimension(:) :: groupBeta
    
    contains
    !--- Define type bound subroutine
    procedure, pass(self) :: writeOutDelay

    end type delayPrecData
    
    contains
    
    !--------------------------------------------------------------------------
    !> @details Defines type/class for delayed neutron precursor data
    !> @param[in] unit_num: unit number to write data to
    !--------------------------------------------------------------------------

    subroutine writeOutDelay(self, unit_num)
    
        implicit none 
        
        !---Dummy variables
        class(delayPrecData), intent(in) :: self
        integer, intent(in) :: unit_num
        
        !---Local variables
        integer :: i, j
        
        do i = 1, size(self % decayConst)
            
            write(unit_num, '(a, i2)')  "Group ", i
            write(unit_num, '(a, f6.3)') "Decay constant [1/s]:      ", self % decayConst(i)
            write(unit_num, '(a, f7.5)') "Delayed fraction (Beta_i): ", self % groupBeta(i)
            
        end do
    
    end subroutine writeOutDelay
    
end module Mod_DelayPrecursorData