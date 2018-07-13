program main 

!*****************************************************************************80
!
!  Discussion:
!
  
  implicit none

  !call timestamp ( )
  write ( *, '(a)' ) 'Test heat solve '

  call driver_solve ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'Normal end of execution.'
!  call timestamp ( )

  stop
end 
