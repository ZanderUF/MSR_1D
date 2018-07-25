!*****************************************************************************80
!
!! Sets up the elemental matrices for the heat equation evaluating matrix entries
!  with Gauss Quadrature  
!
!  Discussion:
!             Generates the elemental matrices K_ij, M_ij, S_ij 
!     	      calculates coefficients in this suboutine 
!  Input: 
!       n - element number
!  Output:
!         elem_matrix_A
!         elem_matrix_U
!         elem_vec_F
!         elem_vec_Q
!*****************************************************************************80

subroutine element_matrix (n, nl_iter)
!
  USE parameters_fe  

  implicit none

!---Dummy variables
    integer  :: g, j, i,ii, ni, n, nl_iter
!---Local variables 
    real , dimension(3) :: elem_coord, velocity
    real , dimension(3, 3) :: K_integral, M_integral
    data M_integral / 4 , 2 , -1 , 2 , 16 , 2 , -1 , 2 ,4 / 
    data K_integral / 7, -8,   1, -8,  16, -8,   1, -8, 7 /
    real , dimension(4)  :: gauspt, gauswt
    data gauspt /-0.8611363116, -0.3399810435, 0.3399810435, 0.8611363116 /  
    data gauswt / 0.347854851 ,  0.6521451548, 0.6521451548, 0.347854851 /
    real  :: xi, wt, cnst, h, s, s2, T, P,  kappa, density, &
             C_p, K_material, F_material, evaluated_velocity
!---Initialize
    elem_matrix_A = 0.0
    elem_matrix_U = 0.0
    
    elem_matrix_F   = 0.0
    elem_vec_f      = 0.0
    
    kappa = 0.0
    density = 0.0 
    C_p = 0.0

!---Integrate over Gauss Pts - assembling only A matrix and source vector 'f'
    do g = 1 , num_gaus_pts 
    !---Calculate shape functions at gauss pt
        xi = gauspt(g)
        wt = gauswt(g)
        h  = global_coord(n,3) - global_coord(n,1) 
        call inter_shape_fcns(xi,h)
        cnst = g_jacobian*wt
        !---Evaluate velocity at gauss pt
        velocity = 0.0
        do i = 1, nodes_per_elem
                velocity(i) = velocity(i) + shape_fcn(i)*velocity_vec(n,i)
        end do
        evaluated_velocity = sum(velocity) 
        !do i = 1 , nodes_per_elem
        !    velocity = velocity + velocity*velocity_vec(n,i)
        !end do

        !---UNIT TEST 
        if(unit_test .eqv. .TRUE.) then
            !---Unit test solver
            do i = 1, nodes_per_elem
                do j = 1, nodes_per_elem
                    !---Assemble A matrix
                    elem_matrix_A(i,j) = elem_matrix_A(i,j) + &
                                         cnst*shape_fcn(i)*shape_fcn(j) 
                    !---Assemble U matrix
                    elem_matrix_U(i,j) = elem_matrix_U(i,j) + &
                                         cnst*shape_fcn(j)*global_der_shape_fcn(i)
                end do
            end do
        end if !---end UNIT TEST if

        !---Normal calculation flow
        if(unit_test .eqv. .FALSE.) then
            
            do i=1, nodes_per_elem
                
                do j = 1, nodes_per_elem
                    !---Assemble A matrix
                    elem_matrix_A(i,j) = elem_matrix_A(i,j) + &
                                         cnst*shape_fcn(i)*shape_fcn(j)
                    !---Assemble P matrix
                    elem_matrix_U(i,j) = elem_matrix_U(i,j) + &
                                        evaluated_velocity*cnst*shape_fcn(j)*global_der_shape_fcn(i)

                    !---Transient calculation
                    if ( steady_state_flag .eqv. .FALSE.) then 
                      
                    end if !---end transient case if 
                
                end do !---End loop over j matrix entries
            end do !---End loop over i matrix entries

        end if !---end unit test if 
    
    end do !---end do over gauss pts 

!-------------------------------------------------------------------------------
    if (DEBUG .eqv. .TRUE. ) then

        write(outfile_unit,fmt='(a)'),' '
        write(outfile_unit,fmt='(a,1I2)'),  'Nonlinear iteration ', nl_iter
        write(outfile_unit,fmt='(a)'),  'Elemental matrices '
        write(outfile_unit,fmt='(a)'),'*****************************************'
        write(outfile_unit,fmt='(a,1I2)'),'[A] element Matrix gaussian integration | element --> ',n
        do j=1,nodes_per_elem 
               write(outfile_unit,fmt='(12es14.3)') &
                    (elem_matrix_A(j,i),i=1,nodes_per_elem)             
        end do
        write(outfile_unit,fmt='(a,1I2)'),'[U] element Matrix gaussian integration | element --> ',n
        do j=1,nodes_per_elem 
               write(outfile_unit,fmt='(12es14.3)') &
                    (elem_matrix_U(j,i),i=1,nodes_per_elem)             
        end do
        
        !write(outfile_unit,fmt='(a)'),' '
        !write(outfile_unit,fmt='(a,1I2)'),'Element solution vector | element: ',n 
        !do j = 1, nodes_per_elem 
        !    ii = j + (n - 1)*nodes_per_elem 
        !    write(outfile_unit,fmt='(a,1I2,a12es14.3)')  'Node #-->', ii, 'Soln:', previous_elem_soln_vec(ii) 
        !end do

            end if !---End matrix write out
!------------------------------------------------------------------------------
 

end 
