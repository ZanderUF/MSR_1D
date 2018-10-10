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

    implicit none

!---Dummy
    real, intent(in) :: T
    real, intent(out) :: heat_capacity
!---Local
    real :: a
    real :: b


    a = 976.78
    b = 1.0634

    heat_capacity = (a + b*T)/1000.0


end subroutine heat_capacity_corr
