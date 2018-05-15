program main

!*****************************************************************************80
!
!! MAIN is the main program for RK4_PRB.
!
!  Discussion:
!     Solve the point kinetics equations using RK4 methods 
!     This code makes use of the RK4 library.
!
!  Licensing:
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!   Starting 5/15/2018 
!
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RK4_PRB'
  write ( *, '(a)' ) 'Solve the point kinetics equations'
  write ( *, '(a)' ) 'Using the RK4 library'

  call rk4vec_test ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RK4_PRB'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end

!*****************************************************************************80
!
!! RK4VEC_TEST tests RK4VEC for a vector ODE.
!

subroutine rk4vec_test ( )
  implicit none

  external rk4vec_test_f
 
  integer ( kind = 4) :: i ! counter 
  integer ( kind = 4 ), parameter :: ndg = 2
  real ( kind = 8 ), parameter :: dt = 0.1D+00
  real ( kind = 8 ) t0 ! starting time
  real ( kind = 8 ) t1 
  real ( kind = 8 ), parameter :: tmax = 10.0D+00 
  
  real ( kind = 8 ) u0(ndg+1) ! initial condition
  real ( kind = 8 ) u1(ndg+1) ! solution vector
  real ( kind = 8 ) beta_i(ndg)
  real ( kind = 8 ) lambda_i(ndg)
  real ( kind = 8 ) gen_time
  real ( kind = 8 ) beta_tot

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RK4VEC_TEST'
  write ( *, '(a)' ) 'RK4VEC takes a Runge Kutta step for a vector ODE.'
  write ( *, '(a)' ) ' '

! Problem parameters - eventually have read in from file
  beta_i(1) = 0.0003
  beta_i(2) = 0.001
  lamda_i(1) = 0.0127
  lambda_i(2) = 0.0317
  gen_time = 0.00002
! Calculate beta total
  do i = 1, ndg
    beta_tot = beta_tot + beta_i(i)
  end do  

! Specify initial conditions
  t0 = 0.0D+00    ! Starting time
  u0(1) = 1.0D+00 ! Value of P(t) initially
! Set precursor density initial conditions  
  do i = 2, ndg+1
      u0(i) = ( beta_i(i-1)*u0(1) ) / ( lamda_i(i-1)*gen_time )
  end do 

! Loop over time steps until we reach tmax
  do
!
!  Print (T0,U0).
!
    write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) t0, u0(1), u0(2)
!
!  Stop if we've exceeded TMAX.
!
    if ( tmax <= t0 ) then
      exit
    end if
!
!  Otherwise, advance to time T1, and have RK4 estimate 
!  the solution U1 there.
!
    t1 = t0 + dt
    call rk4vec ( t0, ndg, u0, dt, rk4vec_test_f, u1)
!
!  Shift the data to prepare for another step.
!
    t0 = t1
    u0(1:ndg) = u1(1:ndg)

  end do

  return
end

!*****************************************************************************80
!
!! RK4VEC_TEST_F evaluates the right hand side of a vector ODE.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) T, the current time.
!
!    Input, integer ( kind = 4 ) N, the dimension of the system.
!
!    Input, real ( kind = 8 ) U(N), the current solution value.
!
!    Output, real ( kind = 8 ) UPRIME(N), the value of the derivative, dU/dT.
!

subroutine rk4vec_test_f ( t, ndg, u, uprime )
  implicit none

  integer ( kind = 4 ) i 
  integer ( kind = 4 ) ndg 

  real ( kind = 8 ) t
  real ( kind = 8 ) u(ndg + 1)
  real ( kind = 8 ) uprime(ndg +1)
  real ( kind = 8 ) sum_lambda
  real ( kind = 8 ) rho
 
! sum over delayed groups, SUM [lambda*precursor] 
  do i=1, ndg
     sum_lamda = sum_lambda + lamda_i(i)*u(i+1)
  end do
  
! Amplitude
! How do we want to represent rho?
  uprime(1) = (rho - beta_tot)/gen_time + sum_lambda 

! Precursors
  uprime(2:ndg+1) = ( beta_i(1:ndg)*u(1) )/gen_time - lambda(1:ndg)*u(2:ndg+1)   
 
  return
end
