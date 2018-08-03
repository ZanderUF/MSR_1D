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
!         heat_elem_matrix_A
!         heat_elem_matrix_K
!         heat_elem_vec_F
!         heat_elem_vec_Q
!*****************************************************************************80

subroutine element_matrix_heat (n, nl_iter)
!
  USE parameters_fe  

  implicit none

!---Dummy variables
    integer :: n 
    integer :: nl_iter
!---Local variables 
    integer  :: g, j, i, ni
    real, dimension(3) :: elem_coord
    real,dimension(3, 3) :: K_integral, M_integral
    data M_integral / 4 , 2 , -1 , 2 , 16 , 2 , -1 , 2 ,4 / 
    data K_integral / 7, -8,   1, -8,  16, -8,   1, -8, 7 /
    !data K_integral / 37, 44, -7, -32, 64, 32, -7, 44, 37 /
    real, dimension(4)  :: gauspt, gauswt
    data gauspt /-0.8611363116, -0.3399810435, 0.3399810435, 0.8611363116 /  
    data gauswt / 0.347854851 ,  0.6521451548, 0.6521451548, 0.347854851 /
    real  :: xi, wt, cnst, h, s, s2, T, P 
    real  :: kappa, density, C_p, K_material, F_material
!---Local element coordinate

    heat_elem_matrix_A = 0.0
    heat_elem_matrix_K = 0.0
    heat_elem_matrix_F   = 0.0
    heat_elem_vec_f      = 0.0
    heat_elem_vec_q      = 0.0
    
    heat_last_elem_D_s1 = 0.0
    heat_elem1_D_s2 = 0.0
    heat_elem1_vec_M_s1 = 0.0
    heat_last_elem_vec_M_s2 = 0.0
    
    !heat_elem1_vec_f    = 0.0
    !heat_last_elem_vec_f    = 0.0

    kappa = 0.0
    density = 0.0 
    C_p = 0.0

    elem_matrix_P_minus = 0.0
    elem_matrix_P_plus = 0.0
!---Setup local coordinates
    do i = 1, nodes_per_elem
        ni = conn_matrix(n,i) 
        elem_coord(i) = global_coord(ni)       
    end do
    
    h = elem_lengths(n)
!---Compute first interface values 
    xi = -1.0
    call inter_shape_fcns(xi,elem_coord,h) 
    do i = 1 , nodes_per_elem
        T = T + shape_fcn(i)*previous_elem_soln_vec( ((2*n-2) + i) )
    end do   
    call kappa_corr(T,kappa)
    do i = 1, nodes_per_elem
        heat_elem1_vec_M_s1(i) = kappa*g_jacobian*shape_fcn(i)
        do j = 1, nodes_per_elem
            heat_last_elem_D_s1(i,j) = heat_last_elem_D_s1(i,j) +  &
                kappa*g_jacobian*shape_fcn(i)*global_der_shape_fcn(j)
        end do
    end do

!---Compute second interface values 
    xi = 1.0
    call inter_shape_fcns(xi,elem_coord,h) 
    do i = 1 , nodes_per_elem
        T = T + shape_fcn(i)*previous_elem_soln_vec( ((2*n-2) + i) )
    end do   
    call kappa_corr(T,kappa)
    do i = 1, nodes_per_elem
        heat_last_elem_vec_M_s2(i) = kappa*g_jacobian*shape_fcn(i)  
        do j = 1, nodes_per_elem
            heat_elem1_D_s2(i,j) = heat_elem1_D_s2(i,j) + &
                kappa*g_jacobian*shape_fcn(i)*global_der_shape_fcn(j)
        end do
    end do

!---Integrate over Gauss Pts - assembling only A matrix and source vector 'f'
    do g = 1 , num_gaus_pts 
        xi = gauspt(g)
        wt = gauswt(g)
        h = elem_lengths(n)
        T = 0.0
        P = 0.0
    !---Calculate shape functions at gauss pt
        call inter_shape_fcns(xi,elem_coord,h)
    !---Evaluate material properties at gauss pt, first get temp @ gauss pt
        do i = 1 , nodes_per_elem
                T = T + shape_fcn(i)*previous_elem_soln_vec( ((2*n-2) + i) )
                P = P + shape_fcn(i)*power_initial( n +i )
        end do   
        call kappa_corr(T,kappa) 
        call density_corr(T,density)
        call cond_corr(T,C_p)
        K_material = kappa/(density*C_p)
        F_material = 1.0/(density*C_p)
        
        cnst = g_jacobian*wt
        
        do i=1, nodes_per_elem
	        !heat_elem_vec_f(i) = P*kappa*cnst*shape_fcn(i)
            do j = 1, nodes_per_elem
                !---Steady state A - no density and C_p, only need on 1st iteration 
            	if ( steady_state_flag .eqv. .TRUE.  ) then

                    heat_elem_vec_f(i) = P*kappa*cnst*shape_fcn(i)*shape_fcn(j)
                    if( n .eq. 1) then 
                         heat_elem1_vec_f(i) = heat_elem_vec_f(i)
                    end if
                    if( n .eq. num_elem) then
                         heat_last_elem_vec_f(i) = heat_elem_vec_f(i)
                    end if        
                    
                    !heat_elem_matrix_A(i,j) = heat_elem_matrix_A(i,j) + &
                    !    kappa*cnst*global_der_shape_fcn(i)*global_der_shape_fcn(j)
                     
                    heat_elem_matrix_A(i,j) = heat_elem_matrix_A(i,j) + cnst*shape_fcn(i)*global_der_shape_fcn(j) 
                    ! Source 
                    ! Save first element values needed later on for partial currents
                    if( n .eq. 1) then 
                        heat_elem1_matrix_A_s2(i,j) = heat_elem_matrix_A(i,j)
                    end if 
                    
                    ! Save last element values needed later on for partial currents
                    if( n .eq. num_elem) then
                        heat_last_elem_matrix_A_s1(i,j) = heat_elem_matrix_A(i,j)
                    end if
                    
                end if ! end steady state if 
                
                ! Transient calculation
                if ( steady_state_flag .eqv. .FALSE.) then 
                    heat_elem_matrix_K(i,j) = heat_elem_matrix_K(i,j) + s2
                    ! M_ij matrix
                    s = shape_fcn(i)*shape_fcn(j)*cnst
                    heat_elem_matrix_A(i,j) = heat_elem_matrix_A(i,j) + s 
                    ! K_ij matrix 
                    s2 = K_material*( der_shape_fcn(i)*der_shape_fcn(j) )*( 1/shape_fcn(j) )
                    ! F_ij matrix    
                    heat_elem_matrix_F(i,j) = heat_elem_matrix_F(i,j) + shape_fcn(i)*F_material         
                    ! Form f vector
                    heat_elem_vec_f(i) = heat_elem_vec_f(i) + heat_elem_matrix_F(i,j)*P*shape_fcn(i)    
                    ! Form q vector - place holder until apply B.C.
                    heat_elem_vec_q(i) = 0.0

                end if 
                 
            end do ! End loop over j matrix entries

        end do ! End loop over i matrix entries

    end do ! end do over gauss pts 

!   Account for modification to first and last elements
    if( n .eq. 1) then
        !heat_elem_matrix_A = heat_elem_matrix_A - heat_elem1_D_s2
    end if 
    if( n .eq. num_elem) then
        !heat_elem_matrix_A = heat_elem_matrix_A + heat_last_elem_D_s1
    end if

    !heat_elem_matrix_K = heat_elem_matrix_K + elem_matrix_P_minus + elem_matrix_P_plus
!-------------------------------------------------------------------------------
    if (DEBUG .eqv. .TRUE. ) then
        if(steady_state_flag .eqv. .FALSE. )  then 
    !---do K matrix analytically  
            do i=1, nodes_per_elem
                !---Setup local coordinates
                ! Get length of element
                h = elem_lengths(n)
                ! Populate Kij, Mij, Fij, Qij
                do j = 1, nodes_per_elem
                    analytic_heat_elem_matrix_A(i,j) =  (h/30.0)*M_integral(i,j)               
                end do ! End loop over j matrix entries
            end do ! End loop over i matrix entries

            write(outfile_unit,fmt='(a)'),' '
            !---Write out each element matrix to the outfile 
            write(outfile_unit,fmt="(a,f8.3)"), 'At time ',t0
            write(outfile_unit,fmt="(a,1I2)"), 'Nonlinear iter: ',nl_iter
            write(outfile_unit,fmt='(a,1I2)'),'M element Matrix analytical: ',n
            do j=1,nodes_per_elem 
                   write(outfile_unit,fmt='(12es14.6)') &
                        (analytic_heat_elem_matrix_A(j,i),i=1,nodes_per_elem)             
            end do
            write(outfile_unit,fmt='(a)'),' ' 
            write(outfile_unit,fmt='(a24,1I2)'),'M element Matrix gauss: ',n
            do j=1,nodes_per_elem 
                   write(outfile_unit,fmt='(12es14.6)') &
                        (heat_elem_matrix_A(j,i),i=1,nodes_per_elem)             
            end do
        end if

        write(outfile_unit,fmt='(a)'),' '
        write(outfile_unit,fmt='(a,1I2)'),  'Nonlinear iteration ', nl_iter
        write(outfile_unit,fmt='(a)'),'*****************************************'
        write(outfile_unit,fmt='(a,1I2)'),'(A-D) element Matrix gaussian integration element--> ',n
        do j=1,nodes_per_elem 
               write(outfile_unit,fmt='(12es14.6)') &
                    (heat_elem_matrix_A(j,i),i=1,nodes_per_elem)             
        end do
        write(outfile_unit,fmt='(a,1I2)'),'D element 1 surface 2 Matrix gaussian integration element--> ',n
        do j=1,nodes_per_elem 
               write(outfile_unit,fmt='(12es14.6)') &
                    (heat_elem1_D_s2(j,i),i=1,nodes_per_elem)             
        end do
        write(outfile_unit,fmt='(a,1I2)'),'D last element surface 1 Matrix gaussian integration element--> ',n
        do j=1,nodes_per_elem 
               write(outfile_unit,fmt='(12es14.6)') &
                    (heat_last_elem_D_s1(j,i),i=1,nodes_per_elem)             
        end do
        
        write(outfile_unit,fmt='(a)'),' '
        write(outfile_unit,fmt='(a26,1I2)'),'f element vector element: ',n
        write(outfile_unit,fmt='(12es14.6)') (heat_elem_vec_f(i),i=1,nodes_per_elem)             
        write(outfile_unit,fmt='(a)'),'*****************************************'
    
!---End matrix write out
    end if
!------------------------------------------------------------------------------
 

end 
