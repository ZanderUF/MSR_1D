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
!         heat_elem_matrix_M
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
    real, dimension(4) :: elem_coord
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

    heat_elem_matrix_M = 0.0
    heat_elem_matrix_K = 0.0
    ! only need to calculate K on the first nonlinear iter since it doesn't change
    !if(steady_state_flag .eqv. .TRUE. .and. nl_iter .eq. 1) then 
    !    heat_elem_matrix_K = 0.0   
    !    analytic_heat_elem_matrix_K_ss = 0.0
    !end if
    heat_elem_matrix_Q = 0.0
    heat_elem_matrix_F = 0.0
    heat_elem_vec_f = 0.0
    heat_elem_vec_q = 0.0
   
    kappa = 0.0
    density = 0.0 
    C_p = 0.0


    elem_matrix_P_minus = 0.0
    elem_matrix_P_plus = 0.0
!   compute interface matrix
    xi = -1.0
    !call inter_shape_fcns(xi,elem_coord,h) 
    !do i = 1, nodes_per_elem
    !    do j = 1, nodes_per_elem
    !        elem_matrix_P_minus(i,j) = elem_matrix_P_minus(i,j) + der_shape_fcn(j)*shape_fcn(i)        
    !    end do
    !end do
    !xi = 1.0
    !!call inter_shape_fcns(xi,elem_coord,h) 
    !do i = 1, nodes_per_elem
    !    do j = 1, nodes_per_elem
    !        elem_matrix_P_plus(i,j) = elem_matrix_P_plus(i,j) + der_shape_fcn(j)*shape_fcn(i)        
    !    end do
    !end do

!---Setup local coordinates
    do i = 1, nodes_per_elem
        ni = conn_matrix(n,i) 
        elem_coord(i) = global_coord(ni)       
    end do
!---Integrate over Gauss Pts
    do g = 1 , num_gaus_pts 
        xi = gauspt(g)
        wt = gauswt(g)
        h = elem_lengths(n)
        T = 0.0
        P = 0.0
        ! Calculate shape functions at gauss pt
        call inter_shape_fcns(xi,elem_coord,h)
        ! Evaluate material properties at gauss pt, first get temp @ gauss pt
        
        do i = 1 , nodes_per_elem 
                T = T + shape_fcn(i)*previous_elem_soln_vec( ((2*n-2) + i) )
                ! Get power at gauss pts
                P = P + shape_fcn(i)*power_initial( n +i )
        end do   
        print *,'T',T 
        
        call kappa_corr(T,kappa) 
        call density_corr(T,density)
        call cond_corr(T,C_p)
        K_material = kappa/(density*C_p)
        F_material = 1.0/(density*C_p)
        ! constant kappa
        !kappa = 20.0
        cnst = g_jacobian*wt
        print *,'kappa ',kappa 
        do i=1, nodes_per_elem
        !---Setup local coordinates
            ! Populate Kij, Mij, Fij, Qij
	        do j = 1, nodes_per_elem
                ! Steady state K - no density and C_p, only need on 1st iteration 
                !if ( steady_state_flag .eqv. .TRUE. .and. nl_iter .eq. 1 ) then 
            	if ( steady_state_flag .eqv. .TRUE.  ) then
                !    heat_elem_matrix_Q(1,j) = heat_elem_matrix_Q(1,j) + &
		   		!                              global_der_shape_fcn(i)*shape_fcn(j)*cnst        

                    heat_elem_matrix_K(i,j) = heat_elem_matrix_K(i,j) + &
                                          kappa*cnst*global_der_shape_fcn(i)*global_der_shape_fcn(j)
                end if 
                
                ! Only update f vector
                if ( steady_state_flag .eqv. .TRUE. ) then 
                    ! F_ij matrix   
                    !heat_elem_matrix_F(i,j) = heat_elem_matrix_F(i,j) + &
		    			      !cnst*shape_fcn(i)*shape_fcn(j)*( 1.0 / kappa ) 
                    !heat_elem_vec_f(i) = heat_elem_vec_f(i) + P*shape_fcn(i)
                    heat_elem_vec_q(i) = 0.0
                end if

                ! Transient calculation
                if ( steady_state_flag .eqv. .FALSE.) then 
                    heat_elem_matrix_K(i,j) = heat_elem_matrix_K(i,j) + s2
                    ! M_ij matrix
                    s = shape_fcn(i)*shape_fcn(j)*cnst
                    heat_elem_matrix_M(i,j) = heat_elem_matrix_M(i,j) + s 
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
   
    heat_elem_matrix_K = heat_elem_matrix_K + elem_matrix_P_minus + elem_matrix_P_plus
!-------------------------------------------------------------------------------
    if (DEBUG .eqv. .TRUE. ) then
        if(steady_state_flag .eqv. .FALSE. )  then 
    !---do K matrix analtically  
            do i=1, nodes_per_elem
                !---Setup local coordinates
                ! Get length of element
                h = elem_lengths(n)
                ! Populate Kij, Mij, Fij, Qij
                do j = 1, nodes_per_elem
                    analytic_heat_elem_matrix_M(i,j) =  (h/30.0)*M_integral(i,j)               
                end do ! End loop over j matrix entries
            end do ! End loop over i matrix entries

            write(outfile_unit,fmt='(a)'),' '
            !---Write out each element matrix to the outfile 
            write(outfile_unit,fmt="(a,f8.3)"), 'At time ',t0
            write(outfile_unit,fmt="(a,1I2)"), 'Nonlinear iter: ',nl_iter
            write(outfile_unit,fmt='(a,1I2)'),'M element Matrix analytical: ',n
            do j=1,nodes_per_elem 
                   write(outfile_unit,fmt='(12es14.6)') &
                        (analytic_heat_elem_matrix_M(j,i),i=1,nodes_per_elem)             
            end do
            write(outfile_unit,fmt='(a)'),' ' 
            write(outfile_unit,fmt='(a24,1I2)'),'M element Matrix gauss: ',n
            do j=1,nodes_per_elem 
                   write(outfile_unit,fmt='(12es14.6)') &
                        (heat_elem_matrix_M(j,i),i=1,nodes_per_elem)             
            end do
        end if

        write(outfile_unit,fmt='(a)'),' '

        !---do K matrix analtically  
        !do i=1, nodes_per_elem
        !    !---Setup local coordinates
        !    ! Get length of element
        !    h = elem_lengths(n)
        !    ! Populate Kij, Mij, Fij, Qij
        !    do j = 1, nodes_per_elem
        !        analytic_heat_elem_matrix_K_ss(i,j) = (1.0/(3.0*h))*K_integral(i,j) 
        !    end do ! End loop over j matrix entries
        !end do ! End loop over i matrix entries
       
        write(outfile_unit,fmt='(a,1I2)'),  'Nonlinear iteration ', nl_iter
        write(outfile_unit,fmt='(a)'),'*****************************************'
        write(outfile_unit,fmt='(a26,1I2)'),'K element Matrix analytically determined element--> ',n
        do j=1,nodes_per_elem 
               write(outfile_unit,fmt='(12es14.6)') &
                    (analytic_heat_elem_matrix_K_ss(j,i),i=1,nodes_per_elem)             
        end do
        write(outfile_unit,fmt='(a26,1I2)'),'K element Matrix gaussian integration element--> ',n
        do j=1,nodes_per_elem 
               write(outfile_unit,fmt='(12es14.6)') &
                    (heat_elem_matrix_K(j,i),i=1,nodes_per_elem)             
        end do

        write(outfile_unit,fmt='(a)'),' '
        write(outfile_unit,fmt='(a26,1I2)'),'Q element Matrix gaussian integration element--> ',n
        do j=1,nodes_per_elem 
               write(outfile_unit,fmt='(12es14.6)') &
                    (heat_elem_matrix_Q(j,i),i=1,nodes_per_elem)             
        end do
        write(outfile_unit,fmt='(a)'),' '
	    write(outfile_unit,fmt='(a26,1I2)'),'P minus element Matrix gaussian integration element--> ',n
        do j=1,nodes_per_elem 
               write(outfile_unit,fmt='(12es14.6)') &
                    (elem_matrix_P_minus(j,i),i=1,nodes_per_elem)             
        end do
        write(outfile_unit,fmt='(a)'),' '
        write(outfile_unit,fmt='(a26,1I2)'),'P plus element Matrix gaussian integration element--> ',n
        do j=1,nodes_per_elem 
               write(outfile_unit,fmt='(12es14.6)') &
                    (elem_matrix_P_plus(j,i),i=1,nodes_per_elem)             
        end do
        write(outfile_unit,fmt='(a)'),' '

	    write(outfile_unit,fmt='(a26,1I2)'),'F element Matrix element: ',n
        do j=1,nodes_per_elem 
               write(outfile_unit,fmt='(12es14.6)') &
                    (heat_elem_matrix_F(j,i),i=1,nodes_per_elem)             
        end do       
        write(outfile_unit,fmt='(a)'),' '
        write(outfile_unit,fmt='(a26,1I2)'),'f element vector element: ',n
        write(outfile_unit,fmt='(12es14.6)') (heat_elem_vec_f(i),i=1,nodes_per_elem)             
        write(outfile_unit,fmt='(a)'),'*****************************************'
    
!---End matrix write out
    end if
!------------------------------------------------------------------------------
 

end 
