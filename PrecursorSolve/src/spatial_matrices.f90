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

    USE element_matrices_M   
    USE flags_M
    USE gauss_integration_M
    USE global_parameters_M
    USE material_info_M
    USE mesh_info_M
    USE time_info_M
    USE solution_vectors_M

    implicit none

!---Dummy variables
    integer, intent(in) :: n
    integer, intent(in) :: nl_iter

!---Local variables 
    integer :: g, j, i,ii, ni
    real(dp) :: cnst,h,xi,wt
    real(dp) :: evaluated_spatial_power, evaluated_velocity
    
    !---Inversion routine parameters
    integer :: lda, info, lwork,length
    integer , dimension(3) :: ipiv
    real(dp), dimension(3) :: work
    
    length = 3
!---Initialize inversion routine parms
    ipiv =  0.0
    work =  0.0
    lda =   length
    lwork = length

!---Initialize
    matrix_W_right_face = 0.0_dp
    matrix_W_left_face  = 0.0_dp
    elem_matrix_U       = 0.0_dp
    elem_matrix_A       = 0.0_dp
    elem_matrix_F       = 0.0_dp
    elem_vec_f          = 0.0_dp
    elem_vec_q          = 0.0_dp    
    elem_vol_int(n,:)   = 0.0_dp 

!---Integrate over Gauss pts 
    gaussintegration: do g = 1 , num_gaus_pts 
    !---Get gauss pt and weight
        xi = gauspt(g)
        wt = gauswt(g)
        h  = global_coord(n,3) - global_coord(n,1) 
        
        !---Evaluate shape functions at gauss pt
        call inter_shape_fcns(xi,h)
        cnst = g_jacobian*wt
        !---Evaluate velocity and power at gauss pt
        evaluated_velocity      = 0.0_dp
        evaluated_spatial_power = 0.0_dp
        do i = 1, nodes_per_elem
            evaluated_velocity = evaluated_velocity + &
                                 shape_fcn(i)*velocity_soln_prev(n,i)
            
            evaluated_spatial_power = evaluated_spatial_power + &
                                      shape_fcn(i)*spatial_power_fcn(n,i)*&
                                      power_amplitude_prev
        end do

        do i=1, nodes_per_elem
            
            elem_vol_int(n,i) = elem_vol_int(n,i) + cnst*shape_fcn(i)
            
            elem_vec_q(i) = elem_vec_q(i)+ &
                    evaluated_spatial_power*cnst*shape_fcn(i)  
            
            do j = 1, nodes_per_elem
                 !---Determine A matrix - only needs to be done once
                elem_matrix_A(i,j) = elem_matrix_A(i,j) + &
                                     cnst*shape_fcn(i)*shape_fcn(j)
                !---Determine U matrix
                elem_matrix_U(i,j) = elem_matrix_U(i,j) + &
                                     evaluated_velocity*cnst*shape_fcn(j)*&
                                     global_der_shape_fcn(i)
            end do !---End loop over j matrix entries
        end do !---End loop over i matrix entries
     enddo gaussintegration!---end do over gauss pts
   
    !elem_matrix_A = elem_matrix_A*2.0_dp

    !---Invert A matrix, only needs to be done once
    if ( (n < 2) .and. (nl_iter < 2) ) then
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
        write(outfile_unit,fmt='(a,12es14.3)'), '@ time = ', t0

        write(outfile_unit,fmt='(a)'),' '
        write(outfile_unit,fmt='(a,4I6)'),  'Nonlinear iteration ', nl_iter
        write(outfile_unit,fmt='(a)'),  'Elemental matrices '
        write(outfile_unit,fmt='(a)'),'*****************************************'
        write(outfile_unit,fmt='(a,4I6)'),'[A] element Matrix gaussian integration | element --> ',n
        do j=1,nodes_per_elem 
               write(outfile_unit,fmt='(12es14.3)') &
                    (elem_matrix_A(j,i),i=1,nodes_per_elem)             
        end do
        write(outfile_unit,fmt='(a)'),' '
        write(outfile_unit,fmt='(a,4I6)'),'[A^-1] element Matrix  | element --> ',n
        do j=1,nodes_per_elem 
               write(outfile_unit,fmt='(12es14.3)') &
                    (inverse_A_matrix(j,i),i=1,nodes_per_elem)             
        end do

        write(outfile_unit,fmt='(a)'),' '
        write(outfile_unit,fmt='(a,4I6)'),'[U] element Matrix gaussian integration | element --> ',n
        do j=1,nodes_per_elem 
               write(outfile_unit,fmt='(12es14.3)') &
                    (elem_matrix_U(j,i),i=1,nodes_per_elem)             
        end do
   end if !---End matrix write out
!------------------------------------------------------------------------------

end subroutine spatial_matrices 
