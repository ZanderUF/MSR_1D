program main

!*****************************************************************************80
!
!! MAIN is the main program for RK4_PRB.
!
!  Discussion:
!     Solve the point kinetics equations using RK4 methods 
!     This code makes use of the RK4 library.
!
  USE parameters
  
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
  write ( *, '(a)' ) 'Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end

!*****************************************************************************80
!
!! RK4VEC_TEST tests RK4VEC for a vector ODE.

subroutine rk4vec_test ( )
  USE parameters
  USE datainput_M

  implicit none

  external rk4vec_test_f
  
  integer ( kind = 4) :: i ! counter 
  real ( kind = 8 ) t1 
  real, allocatable :: u0(:) ! initial condition
  real, allocatable :: u1(:) ! solution vector
  character(50) :: file_name
  character, dimension(:), allocatable :: precursor_txt*10
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'RK4VEC_TEST'
  write ( *, '(a)' ) 'RK4VEC takes a Runge Kutta step for a vector ODE.'
  write ( *, '(a)' ) ' '

! Read in problem parameters here
  call datainput 
  
  allocate(u0(ndg+1), u1(ndg+1))
  allocate(precursor_txt(ndg))
  u0(:)=0.0
  u1(:)=0.0 
! Write out files depending on problem type
  if(ramp .eqv. .TRUE.) then
      write(file_name,'(a,f5.4,a,f5.4,a,f8.6)'),"out_ramp_t_",t_initial,"_to_",t_final,"tstep_",dt
  else if(step .eqv. .TRUE.) then
      write(file_name,'(a,f5.4,a,f5.4)'),"out_step_t_",t_final
  else 
      write(file_name, '(a)'),"output.txt"
  end if 
! Open file for writing out solution
  open (unit=99, file=file_name,status='unknown',form='formatted',position='asis')
  do i=1,ndg
     write(unit=precursor_txt(i), fmt='(A9,I1)') 'precursor', i  
  end do  
  write(99, fmt='(A,A,A,10(A10,4X))'),'     Time(s)   ','     Amplitude     ',(precursor_txt(i), i=1,ndg) 

! Calculate beta total
  do i = 1, ndg
      beta_tot = beta_tot + beta_i(i)
  end do  

! Specify initial conditions
  t0    = 0.0D+00    ! Starting time
  u0(1) = 10.0D+00 ! Value of P(t) initially
! Set precursor density initial conditions  
  do i = 2, ndg+1
      u0(i) = ( beta_i(i-1)*u0(1) ) / ( lamda_i(i-1)*gen_time )
  end do 
  print *,'step',step 
! Loop over time steps until we reach tmax
  do
     ! Step perturbation
     if ( step .eqv. .TRUE.) then
         if (t0 > t_final) then
             rho = -0.005
         end if
     end if     
     ! Ramp perturbation
     if(ramp .eqv. .TRUE.) then
         if(t0< t_final) then
             rho = rho_initial + ((rho_final - rho_initial)*(t0-t_initial))/(t_final - t_initial)  
         end if
     end if
     write (99, '(2x,g14.6,2x,g14.6,2x,g14.7)' ) t0, u0(1), u0(2)
!
!    Stop if we've exceeded TMAX.
!
     if ( tmax <= t0 ) then
       exit
     end if
!
!    Otherwise, advance to time T1, and have RK4 estimate 
!    the solution U1 there.
!
      t1 = t0 + dt
      call rk4vec ( t0, ndg+1, u0, dt, rk4vec_test_f, u1)
!    Shift the data to prepare for another step.
!
      t0 = t1
      u0(1:ndg+1) = u1(1:ndg+1)
 end do

 deallocate(lamda_i,beta_i)
 deallocate(u0,u1)
 deallocate(precursor_txt)
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

subroutine rk4vec_test_f ( t, u, uprime )
  
  USE parameters

  implicit none

  integer ( kind = 4 ) i 

  real ( kind = 8 ) t
  real :: u(ndg + 1)
  real :: uprime(ndg +1)
  real ( kind = 8 )  sum_lamda
  
  sum_lamda = 0.0
! sum over delayed groups, SUM [lamda*precursor] 
  do i=1, ndg
     sum_lamda = sum_lamda + lamda_i(i)*u(i+1)
  end do
  
! Power
  uprime(1) = ((rho - beta_tot)/gen_time)*u(1) + sum_lamda 
! Precursors
  uprime(2:ndg+1) = (beta_i(1:ndg)/gen_time)*u(1) - lamda_i(1:ndg)*u(2:ndg+1)   

return
end
