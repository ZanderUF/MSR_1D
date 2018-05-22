!*****************************************************************************80
! Solves heat equation using FOrward Time Centerd Space approxiamtion 
!!  

subroutine temp_solve ( )
  
  USE parameters
  USE datainput_M

  implicit none

  integer ( kind = 4 ) :: i   ! counter 
  real    ( kind = 8 ) :: t1  ! next time step  
  real, allocatable    :: u0(:) ! initial condition
  real, allocatable    :: u1(:) ! solution vector
  
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Simple Temperature solve'
  write ( *, '(a)' ) ' '

! Read in problem parameters here
  call datainput 
  
  allocate(u0(ndg+1), u1(ndg+1))
  u0(:)=0.0
  u1(:)=0.0 

! Name the output files something useful 
  call proper_file_namer

! Open file for writing out solution
  open (unit=99, file=file_name,status='unknown',form='formatted',position='asis')
  
  !write(99, fmt='(A,A,A,10(A12,3X))'),'    Time(s)    ','    Amplitude    ',(precursor_txt(i), i=1,ndg) 

! Specify initial conditions
  t0    = 0.0D+00    ! Starting time
  u0(1) = 10.0D+00   ! Value of P(t) initially

! Loop over time steps until we reach tmax
  do

     write (99, '(2x,es14.6,12es14.6)' ) t0, u0(1), (u0(i),i=2,ndg+1) 
!
!    Stop if we've exceeded TMAX.
!
     if ( tmax <= t0 ) then
         exit
     end if
!
     t1 = t0 + dt

!    Shift the data to prepare for another step.
     t0 = t1
     u0(1:ndg+1) = u1(1:ndg+1)
 end do

deallocate(u0,u1)

return
end

!*****************************************************************************80
!
! proper_file_namer names the output data file based on the problem type
! 
! No input parameters.  Uses file_name and ramp from parameters module
!
 
subroutine proper_file_namer()

USE parameters

! Write out files depending on problem type
  if(ramp .eqv. .TRUE.) then
      write(file_name,'(a,f5.4,a,f5.4,a,f2.8)'),"out_ramp_t_",t_initial,"_to_",t_final,"tstep_",dt
  else if(step .eqv. .TRUE.) then
      write(file_name,'(a, f6.4, a,f8.6)'),"out_step_t_at_",t_final,"_tstp_", dt
  else 
      write(file_name, '(a)'),"output.txt"
  end if

return

end 
