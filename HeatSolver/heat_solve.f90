!*****************************************************************************80
!
!! Main subroutine for the heat solve 

subroutine heat_solve ( )
  
USE parameters_fe
USE datainput_fe_M
  
implicit none

  integer ( kind = 4) :: i   ! counter 
  real ( kind = 8 )   :: t1  ! next time step  
  
! Read in problem parameters here
  call datainput_fe 
    
! Name the output files something useful 
  call proper_file_namer

  outfile_unit = 11
! Open file for writing out debug information
  open (unit=outfile_unit, file="outfile.txt",status='unknown',form='formatted',position='asis')
! Open file for writing out solution
  open (unit=99, file=file_name,status='unknown',form='formatted',position='asis')

! Create 1D mesh
  call mesh_creation

! General elemental matrices
 

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

