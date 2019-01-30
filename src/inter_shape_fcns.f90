!*****************************************************************************80
!
!! Evaluate the shape functions at gaus pts  
!  Discussion:
!             Calculates shape fcns and derivatives at the gauss pts
!  Parameters:

subroutine inter_shape_fcns (xi, h)
!
    USE global_parameters_M
    USE gauss_integration_M
    USE flags_M
    USE mesh_info_M

    implicit none
!---Dummy variable
    real(dp) :: xi 
    !real, dimension(3) :: elem_coord
    real(dp) :: h

!---local
    integer :: i   
   
    lagrange = .TRUE.
    ! Calculate shape function @ gauss pt
    if(lagrange .eqv. .TRUE.) then 
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
    end if
    hermite = .FALSE.
!---Calculate using hermite polynomials 
    !if(hermite .eqv. .TRUE.) then
    !    
    !    shape_fcn(1) = 0.25*(2.0 - 3.0*xi + xi**3) 
    !    shape_fcn(2) = -h*(1.0 - xi)*(1.0 - xi*xi)/8.0
    !    shape_fcn(3) = 0.25*(2.0 + 3.0*xi - xi**3)
    !    shape_fcn(4) =  h*(1.0 + xi)*(1.0 - xi*xi)/8.0
    !    
    !    der_shape_fcn(1)= -0.75*(1.0-xi*xi) 
    !    der_shape_fcn(2)=  h*(1.0+2.0*xi - 3.0*xi*xi)/8.0
    !    der_shape_fcn(3)=  0.75*(1.0 - xi*xi) 
    !    der_shape_fcn(4)=  h*(1.0 - 2.0*xi - 3.0*xi*xi)/8.0 
    !    g_jacobian = h*0.5
    !    
    !    ! Compute global derivative
    !	do i = 1, nodes_per_elem
    !	    global_der_shape_fcn(i) = der_shape_fcn(i)/g_jacobian 
    !	end do
    !
    !end if

end  
