! Evalute kappa (thermal conductivity) based on correlation
!
! Correlation is of the form k(T) = a + b*T where a, b are constants
! Currently using FliNaK correlation from Smirnov (1987)  
! Input:
!       T - temperature to evaluate kappa at 
! Output:
!       kappa - evaluated conductivity at T
! 
subroutine kappa_corr(T, kappa)

implicit none

! Dummy
    real(kind=8), intent(in)  :: T
    real(kind=8), intent(out) :: kappa
! Local
    real(kind=8) :: a
    real(kind=8) :: b  
    
    a = 0.36 
    b = 0.00056
    kappa = a + b*T

end
