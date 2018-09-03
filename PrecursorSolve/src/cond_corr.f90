! Evalute C_p (thermal conductivity as function of temperature based on correlation
! C_p units - cal/(mol*K) , 1 cal = 4.184 J

! Correlation is of the form C_p(T) = a + b*T where a, b are constants
! Currently using FliNaK correlation from Rogers (1982)
! Valid between 748 - 863 K
! Input:
!       T - temperature to evaluate density at in Kelvin 
! Output:
!       C_p - evaluated conductivity at T
! 
subroutine cond_corr(T, C_p)

implicit none

! Dummy
    real, intent(in)  :: T
    real, intent(out) :: C_p 
! Local
    real :: a
    real :: b  
    
    a = 9.636 
    b = 10.487
    C_p = a + b*T

!   Test nonlinear correlation

end
