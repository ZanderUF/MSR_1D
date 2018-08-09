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
    real , dimension(3) :: elem_coord, velocity, temp_prec,shape_int 
    real  :: xi, wt, cnst, h, s, s2, T, P,  kappa, density, &
             C_p, K_material, F_material, evaluated_amplitude, &
             evaluated_velocity
    !---Inversion routine parameters
    integer :: lda, info, lwork,length
    integer, dimension(3) :: ipiv
    real, dimension(3) :: work
    real :: prec_test
    length = 3
!---Initialize inversion routine parms
    ipiv =  0
    work =  0.0
    lda =   length
    lwork = length
!---Initialize
    matrix_W_right_face = 0
    matrix_W_left_face = 0
    elem_matrix_U = 0.0
    elem_matrix_A = 0.0 
    elem_matrix_F   = 0.0
    elem_vec_f      = 0.0
    elem_vec_q   = 0.0    
!---Material properties 
    kappa = 0.0
    density = 0.0 
    C_p = 0.0

    elem_vol_int(n,:) = 0 
!---Integrate over Gauss Pts - assembling only A matrix and source vector 'f'
    do g = 1 , num_gaus_pts 
    !---Calculate shape functions at gauss pt
        xi = gauspt(g)
        wt = gauswt(g)
        h  = global_coord(n,3) - global_coord(n,1) 
        call inter_shape_fcns(xi,h)
        cnst = g_jacobian*wt
        !---Evaluate velocity at gauss pt
        evaluated_velocity = 0.0
        evaluated_amplitude = 0.0
        
        do i = 1, nodes_per_elem
            evaluated_velocity = evaluated_velocity + &
                                 shape_fcn(i)*velocity_soln_prev(n,i)
            evaluated_amplitude = evaluated_amplitude + &
                                 shape_fcn(i)*amplitude_fcn(n,i)*total_power_prev
        end do
    
        !---Normal calculation flow
        if(unit_test .eqv. .FALSE.) then
            
            do i=1, nodes_per_elem
                elem_vol_int(n,i) = elem_vol_int(n,i) + cnst*shape_fcn(i)
                !---Assemble q vector
                !elem_vec_q(i) = elem_vec_q(i) + &
                !                cnst*shape_fcn(i)*evaluated_amplitude!*total_power_prev 
                !elem_vol_int(n,i) = elem_vol_int(n,i) +elem_vec_q(i)
                do j = 1, nodes_per_elem
                    
                    !---Assemble A matrix - only needs to be done once
                    elem_matrix_A(i,j) = elem_matrix_A(i,j) + &
                                         cnst*shape_fcn(i)*shape_fcn(j)
                    !---Assemble P matrix
                    elem_matrix_U(i,j) = elem_matrix_U(i,j) + &
                                        evaluated_velocity*cnst*shape_fcn(j)*global_der_shape_fcn(i)
                    !---Transient calculation
                    if ( steady_state_flag .eqv. .FALSE.) then 
                        
                    end if!---end transient case if 
                
                end do !---End loop over j matrix entries
            end do !---End loop over i matrix entries
        end if !---end unit test if 
    
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
    end do !---end do over gauss pts 

    do i = 1, nodes_per_elem
        do j = 1, nodes_per_elem
        
            elem_vec_q(i) = elem_vec_q(i) + elem_matrix_A(i,j)*power_soln_new(n,j) !amplitude_fcn(n,j)!*total_power_prev
           
        !---Applies for all elements except the first one
            if(n > 1) then !--- n - element #
                !---Grab previous precursor conc. + velocity at 
                !---rhs of previous element
                 matrix_W_right_face(i,j) = velocity_soln_prev(n,i)*&
                                            interp_fcn_rhs(i)*interp_fcn_rhs(j)  
                 matrix_W_left_face(i,j) = velocity_soln_prev(n,i)*&
                                           interp_fcn_lhs(i)*interp_fcn_lhs(j)
            else!---First element case, need to connect with end element 
                 matrix_W_right_face(i,j) = velocity_soln_prev(num_elem,i)*&
                                            interp_fcn_rhs(i)*interp_fcn_rhs(j)  
                 matrix_W_left_face(i,j) =  velocity_soln_prev(num_elem,i)*&
                                            interp_fcn_lhs(i)*interp_fcn_lhs(j)
            end if!---End flux calculation
        end do
    end do
    
!---Invert A matrix
    do i = 1, nodes_per_elem
        do j = 1, nodes_per_elem
            inverse_A_matrix(i,j) = elem_matrix_A(i,j)
        end do
    end do
    
    !---Factorize matrix
    call dgetrf ( length, length, inverse_A_matrix, lda, ipiv, info )
    !---Compute the inverse matrix.
    call dgetri ( length, inverse_A_matrix, lda, ipiv, work, lwork, info ) 
    
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
        write(outfile_unit,fmt='(a)'),' '
        write(outfile_unit,fmt='(a,1I2)'),'[A^-1] element Matrix  | element --> ',n
        do j=1,nodes_per_elem 
               write(outfile_unit,fmt='(12es14.3)') &
                    (inverse_A_matrix(j,i),i=1,nodes_per_elem)             
        end do
        write(outfile_unit,fmt='(a)'),' '
        write(outfile_unit,fmt='(a,1I2)'),'[U] element Matrix gaussian integration | element --> ',n
        do j=1,nodes_per_elem 
               write(outfile_unit,fmt='(12es14.3)') &
                    (elem_matrix_U(j,i),i=1,nodes_per_elem)             
        end do
        write(outfile_unit,fmt='(a)'), ' '
        write(outfile_unit,fmt='(a,1I2)'),'[W] right Matrix | element --> ',n
        do j=1,nodes_per_elem 
              write(outfile_unit,fmt='(12es14.3)') &
                   (matrix_W_right_face(j,i),i=1,nodes_per_elem)             
        end do
        
        write(outfile_unit,fmt='(a)'), ' '
        write(outfile_unit,fmt='(a,1I2)'),'[W] left Matrix | element --> ',n
        do j=1,nodes_per_elem 
              write(outfile_unit,fmt='(12es14.3)') &
                   (matrix_W_left_face(j,i),i=1,nodes_per_elem)             
        end do       
   end if !---End matrix write out
!------------------------------------------------------------------------------

end 
