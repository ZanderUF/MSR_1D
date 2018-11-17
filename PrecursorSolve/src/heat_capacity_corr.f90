! Heat capacity as a function of temperature
! From:
!   Williams, D. F. (2006a). “Assessment of Candidate 
!   Molten Salt Coolants for the NGNP/NHI Heat-
!   Transfer Loop,” Oak Ridge National Laboratory Report ORNL/TM-2006/69
!
! Units: [J/g-K]
!
! Input: Temperature [K]
!
! Output: heat_capacity
!

subroutine heat_capacity_corr(T, heat_capacity)

    USE global_parameters_M

    implicit none

!---Dummy
    real(dp), intent(in) :: T
    real(dp), intent(out) :: heat_capacity
!---Local
    real(dp) :: a
    real(dp) :: b


    a = 976.78_dp
    b = 1.0634_dp

    heat_capacity = (a + b*T)/1000.0_dp


end subroutine heat_capacity_corr
