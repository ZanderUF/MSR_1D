program msr1d 

!*****************************************************************************80
!
!  Discussion:
!
  
  implicit none

  !call timestamp ( )
  write ( *, '(a)' ) 'Begin Solve'

  call driver_solve ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Normal end of execution.'
!  call timestamp ( )

  stop
end 
