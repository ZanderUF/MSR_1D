    module parameters 

    implicit none
   
    integer (kind = 4) :: ndg 
    real, allocatable :: beta_i(:)
    real, allocatable :: lamda_i(:)
    real ( kind = 8 ) gen_time
    real ( kind = 8 ) beta_tot
    real ( kind = 8 ) rho ! starting value
    real ( kind = 8 ) rho_initial ! initial value for ramp
    real ( kind = 8 ) rho_final   ! final value for rho ramp
    logical :: step = .FALSE.
    logical :: ramp = .FALSE.

    real ( kind = 8 ) t0 ! starting time
    real ( kind = 8 ) dt ! time step 
    real ( kind = 8 ) tmax  
    real ( kind = 8 ) t_initial
    real ( kind = 8 ) t_final 

end module parameters 
