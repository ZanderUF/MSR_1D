! Evalute density as functino of temperature based on correlation
!
! Correlation is of the form density(T) = a + b*T where a, b are constants
! Currently using FliNaK correlation from Salanne (2009) -    
! (Heat-transport properties of molten fluorides: Determination from
!  first-principles)
! Valid from 750 - 1100 K
! Input:
!       T - temperature to evaluate density at 
! Output:
!       density - evaluated conductivity at T [ kg/m^3]
! 
subroutine density_corr(T, density)

implicit none

! Dummy
    real, intent(in)  :: T
    real, intent(out) :: density 
! Local
    real :: a
    real :: b  
   
    a = 2603
    b= -6.69E-1

    density = a + b*T

end
