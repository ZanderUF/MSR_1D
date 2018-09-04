!*****************************************************************************80
!
!  Sets up the elemental matrices 
!
!  Discussion:
!             Generates the elemental matrices K_ij, M_ij, S_ij 
!             calculates coefficients in this suboutine 
!  Input: 
!       n - element number
!       nl_iter - nonlinear iteration counter
!  Output:
!         elem_matrix_A
!         elem_matrix_U
!         elem_vec_F
!         elem_vec_Q
!*****************************************************************************80

subroutine spatial_matrices (n, nl_iter)
!
  USE parameters_fe  

  implicit none

!---Dummy variables
    integer  :: g, j, i,ii, ni, n, nl_iter
!---Local variables 
    real , dimension(3) :: elem_coord, velocity, temp_prec,shape_int 
    real  :: xi, wt, cnst, h, s, s2, T, P,  kappa, density, &
             C_p, K_material, F_material, evaluated_spatial_power, &
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
    elem_matrix_F   = 0.0
    elem_vec_f      = 0.0
    elem_vec_q   = 0.0    
!---Material properties 
    kappa = 0.0
    density = 0.0 
    C_p = 0.0

    elem_vol_int(n,:) = 0 
!---Integrate over Gauss pts 
    do g = 1 , num_gaus_pts 
    !---Get gauss pt and weight
        xi = gauspt(g)
        wt = gauswt(g)
        h  = global_coord(n,3) - global_coord(n,1) 
        !---Evaluate shape functions at gauss pt
        call inter_shape_fcns(xi,h)
        cnst = g_jacobian*wt
        !---Evaluate velocity and power at gauss pt
        evaluated_velocity = 0.0
        evaluated_spatial_power = 0.0
        do i = 1, nodes_per_elem
            evaluated_velocity = evaluated_velocity + &
                                 shape_fcn(i)*velocity_soln_prev(n,i)
            evaluated_spatial_power = evaluated_spatial_power + &
                                 shape_fcn(i)*spatial_power_fcn(n,i)*power_amplitude_prev
        end do
    
        !---Normal calculation flow
        if(unit_test .eqv. .FALSE.) then
            
            do i=1, nodes_per_elem
                elem_vol_int(n,i) = elem_vol_int(n,i) + cnst*shape_fcn(i)
                
                do j = 1, nodes_per_elem
                    if(n < 2 .and. nl_iter < 2) then 
                        !---Assemble A matrix - only needs to be done once
                        elem_matrix_A(i,j) = elem_matrix_A(i,j) + &
                                         cnst*shape_fcn(i)*shape_fcn(j)
                    end if
                    !---Assemble P matrix
                    elem_matrix_U(i,j) = elem_matrix_U(i,j) + &
                                        evaluated_velocity*cnst*shape_fcn(j)*global_der_shape_fcn(i)
                    !---Transient calculation
                    if ( steady_state_flag .eqv. .FALSE.) then 
                        
                    end if!---end transient case if 
                
                end do !---End loop over j matrix entries
            end do !---End loop over i matrix entries
        end if !---end unit test if 
    
    end do !---end do over gauss pts 
!---Create source vector 'q', and W - 
    do i = 1, nodes_per_elem
        do j = 1, nodes_per_elem
            elem_vec_q(i) = elem_vec_q(i) + &
                            elem_matrix_A(i,j)*spatial_power_fcn(n,j)*power_amplitude_prev
            !---Applies for all elements except the first one
            if(n > 1) then !--- n - element #
                !---Grab previous precursor conc. + velocity at 
                !---rhs of previous element
                 matrix_W_right_face(i,j) = velocity_soln_prev(n,i)*&
                                            interp_fcn_rhs(i)*interp_fcn_rhs(j)  
                 matrix_W_left_face(i,j)  = velocity_soln_prev(n,i)*&
                                            interp_fcn_lhs(i)*interp_fcn_lhs(j)
            else!---First element case, need to connect with end element 
                 matrix_W_right_face(i,j) = velocity_soln_prev(num_elem,i)*&
                                            interp_fcn_rhs(i)*interp_fcn_rhs(j)  
                 matrix_W_left_face(i,j)  = velocity_soln_prev(num_elem,i)*&
                                            interp_fcn_lhs(i)*interp_fcn_lhs(j)
            end if!---End flux calculation
        end do
    end do
   
    !---Invert A matrix, only needs to be done once
    if ( n < 2 .and. nl_iter < 2) then
        do i = 1, nodes_per_elem
            do j = 1, nodes_per_elem
                inverse_A_matrix(i,j) = elem_matrix_A(i,j)
            end do
        end do
        
        !---Factorize matrix
        call dgetrf ( length, length, inverse_A_matrix, lda, ipiv, info )
        !---Compute the inverse matrix.
        call dgetri ( length, inverse_A_matrix, lda, ipiv, work, lwork, info ) 
    end if

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

end subroutine spatial_matrices 

