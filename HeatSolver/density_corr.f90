! Evalute density as functino of temperature based on correlation
!
! Correlation is of the form density(T) = a + b*T where a, b are constants
! Currently using FliNaK correlation from Cherenkova (2003)  
! Input:
!       T - temperature to evaluate density at 
! Output:
!       density - evaluated conductivity at T
! 
subroutine density_corr(T, density)

implicit none

! Dummy
    real(kind=8), intent(in)  :: T
    real(kind=8), intent(out) :: density 
! Local
    real(kind=8) :: a
    real(kind=8) :: b  
    
    a = 2579.3 
    b = -0.624
    density = a + b*T

end
