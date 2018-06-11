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
! proper_file_namer names the output data file based on the problem type
! 
! No input parameters.  Uses file_name and ramp from parameters module
!
 
subroutine proper_file_namer()

USE parameters_fe

! Write out files depending on problem type
  if(time_solve .eqv. .TRUE.) then
      write(file_name, '(a)'),'time_dep_output.txt'
      !write(file_name,'(a,f5.4,a,f5.4,a,f2.10)'),"out_t_",t_initial,"_to_",tmax,"tstep_",dt
  else 
      write(file_name, '(a)'),"stationary_output.txt"
  end if

return

end
