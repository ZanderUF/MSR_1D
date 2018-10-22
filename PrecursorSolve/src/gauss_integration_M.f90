module gauss_integration_M 

implicit none

!---Gauss pts
    double precision , dimension(4)  :: gauspt, gauswt
    data gauspt /-0.8611363116, -0.3399810435, 0.3399810435, 0.8611363116 /  
    data gauswt / 0.347854851 ,  0.6521451548, 0.6521451548, 0.347854851 /

    double precision :: g_jacobian 
!---Gauss integration 
    integer  :: num_gaus_pts = 4

!---Shape functions, Lagrange, quadratic order
    double precision , dimension(3) :: shape_fcn
    double precision , dimension(3) :: der_shape_fcn
    double precision , dimension(3) :: global_der_shape_fcn 

end module gauss_integration_M
