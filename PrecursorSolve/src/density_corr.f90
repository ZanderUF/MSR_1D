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
!       density - evaluated conductivity at T [ g/cm^3]
! 
subroutine density_corr(T, density)

    USE global_parameters_M

implicit none

! Dummy
    real(dp), intent(in)  :: T
    real(dp), intent(out) :: density 
! Local
    real(dp) :: a
    real(dp) :: b  
   
    a = 2.603_dp
    b= -6.69E-4_dp

    density = a + b*T

end
