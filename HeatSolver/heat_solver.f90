program main

!*****************************************************************************80
!
!! MAIN is the main program for solving heat equation 
!
!  Discussion:
!      This code solves the time-dependent non-linear heat equation
!      The non-linearity stemming from material evaluations that depend
!      on the temperature		
!
  USE parameters_fe
  USE timestamp_M
  
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Solve the time-dependent heat equation'
  write ( *, '(a)' ) ' '

  call heat_solve ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop
end

!*****************************************************************************80
!
!! Main subroutine for the heat solve 

subroutine heat_solve ( )
  
  USE parameters_fe
  USE datainput_fe_M

  implicit none

  integer ( kind = 4) :: i   ! counter 
  real ( kind = 8 )   :: t1  ! next time step  
  real, allocatable :: u0(:) ! initial condition
  real, allocatable :: u1(:) ! solution vector
  
! Read in problem parameters here
  call datainput_fe 
  
!  allocate(u0(ndg+1), u1(ndg+1), )
!  u0(:)=0.0
!  u1(:)=0.0 

! Name the output files something useful 
  call proper_file_namer

! Open file for writing out solution
  open (unit=99, file=file_name,status='unknown',form='formatted',position='asis')

!  Create 1D mesh
!  call mesh_creation

! Specify initial conditions
  t0    = 0.0D+00    ! Starting time
  !u0(1) = 300.0D+00   ! Value initially

! Loop over time steps until we reach tmax
!  do
!
!    Stop if we've exceeded TMAX.
!
!     if ( tmax <= t0 ) then
!         exit
!     end if
!
!     t1 = t0 + dt
!
!    Shift the data to prepare for another step.
!      t0 = t1
!      u0(1:ndg+1) = u1(1:ndg+1)
! end do
! deallocate(u0,u1)

return
end

!*****************************************************************************80
!
! proper_file_namer names the output data file based on the problem type
! 
! No input parameters.  Uses file_name and ramp from parameters module
!
 
subroutine proper_file_namer()

USE parameters_fe

! Write out files depending on problem type
  if(time_solve .eqv. .TRUE.) then
      write(file_name,'(a,f5.4,a,f5.4,a,f2.10)'),"out_t_",t_initial,"_to_",tmax,"tstep_",dt
  else 
      write(file_name, '(a)'),"stationary_output.txt"
  end if

return

end
