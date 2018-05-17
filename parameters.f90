    module parameters 

    implicit none
   
    integer (kind = 4) :: ndg 
	real, allocatable :: beta_i(:)
    real, allocatable :: lamda_i(:)
    real ( kind = 8 ) gen_time
    real ( kind = 8 ) beta_tot
    real ( kind = 8 ) rho
    real ( kind = 8 ) rho_initial   
    real ( kind = 8 ) rho_final
    logical :: step = .FALSE.
    logical :: ramp = .FALSE.

   end module parameters 
