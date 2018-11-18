! Evalute density as functino of temperature based on correlation
!
! Correlation is of the form density(T) = a + b*T where a, b are constants
!
! G. J. Janz, R. P. T. Tomkins, Physical properties data compilations relevant to energy storage. IV. Molten salts: Data on Additional Single and Multi- Component Salt Systems. NSRDS-NBS-61(Pt.IV), U.S. National Bureau of Standards, 1981
! This was done for KCl-MgCl2 salts
!
! Valid from 1017 - 1174 K
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
   
    a = 2007.0_dp
    b= -0.4571_dp

    density = (a + b*T)*0.001_dp

end
