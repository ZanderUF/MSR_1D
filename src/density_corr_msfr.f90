! Evalute density as functino of temperature based on correlation
!
! Correlation is of the form density(T) = a + b*T where a, b are constants
! From MSFR paper 'Modelling and analysis of the MSFR transient behaviour'
! Input:
!       T - temperature to evaluate density at 
! Output:
!       density - evaluated conductivity at T [ g/cm^3]
! 
subroutine density_corr_msfr(T, density)

    USE global_parameters_M

implicit none

!---Dummy
    real(dp), intent(in)  :: T
    real(dp), intent(out) :: density 
!---Local
    real(dp) :: a
    real(dp) :: b  
   
    a = 4094_dp
    b = -8.82E-1_dp

    density = a + b*(T-1008_dp)*1E-3_dp

end subroutine density_corr_msfr
