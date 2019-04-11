! Heat capacity as a function of temperature
! From:
!   Williams, D. F. (2006a). “Assessment of Candidate 
!   Molten Salt Coolants for the NGNP/NHI Heat-
!   Transfer Loop,” Oak Ridge National Laboratory Report ORNL/TM-2006/69
!   For FLiNaK salt ***   
! 
!   There are no correlations for heat capacity of Potassium Chloride- Magnesium Chloride salt. There is only one constant value of
!   1.150 kJ/kg K given by Ambrosek [2010]. for KCl-MgCl_2
!
!   Taube and Ligou assumed 0.95 kJ/kg K for UCl_3 @ 983 C   
! 
!   Modified FLiNak constants to give us what we want @ nominal temperature
!   basd on published data.  Assume temperature dependence is somewhat similar
!   between salts. 
!
! Units: [J/g-K]
!
! Input: Temperature [K]
!
! Output: heat_capacity
!

subroutine heat_capacity_corr(T, heat_capacity)

    USE global_parameters_M
    USE flags_M

    implicit none

!---Dummy
    real(dp), intent(in) :: T
    real(dp), intent(out) :: heat_capacity
!---Local
    real(dp) :: a
    real(dp) :: b


    a = 76.78_dp
    !b = 1.0634_dp

    !heat_capacity = 1.150
    !heat_capacity = 900.0
    heat_capacity = 0.90
    
    if(MSRE_problem .eqv. .TRUE.) then

        heat_capacity = 1.74

    end if

    !heat_capacity = (a + b*T)/1000.0_dp

end subroutine heat_capacity_corr
