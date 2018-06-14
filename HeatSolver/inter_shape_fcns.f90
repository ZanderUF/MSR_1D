!*****************************************************************************80
!
!! Evaluate the shape functions at gaus pts  
!  Discussion:
!             Calculates shape fcns and derivatives at the gauss pts
!  Parameters:

subroutine inter_shape_fcns (xi,elem_coord, h)
!
    USE parameters_fe  

    implicit none
!---Dummy variable
    real :: xi 
    real, dimension(3) :: elem_coord
    real :: h

!---local
    integer :: i   
    
! Calculate shape function @ gauss pt
    shape_fcn(1) = -0.5*xi*(1 - xi) 
    shape_fcn(2) = (1 + xi)*(1 - xi)
    shape_fcn(3) = 0.5*xi*(1 + xi)

    ! Derivative of shape function @ gauss pt
    der_shape_fcn(1) = xi - 0.5 
    der_shape_fcn(2) = -2.0*xi
    der_shape_fcn(3) = xi + 0.5
    g_jacobian = h*0.5 
    !! Compute jacobian
    !do i = 1, nodes_per_elem
    !    g_jacobian = g_jacobian + der_shape_fcn(i)/elem_coord(i)
    !end do

    ! Compute global derivative
    do i = 1, nodes_per_elem
        global_der_shape_fcn(i) = der_shape_fcn(i)/g_jacobian 
    end do

end  
