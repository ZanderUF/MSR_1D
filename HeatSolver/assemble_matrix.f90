!*****************************************************************************80
!
!! Assemble the global matrix for the problem 
!  Discussion:
!             As the elemental matrices are created we add them to the global 
!  Parameters:

subroutine assemble_matrix ()
!
    USE parameters_fe  

    implicit none
!---Dummy variable

!---local
    integer :: i   

!---No need for elemental matrices after they have been placed in the global one

end 
