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
    integer ( kind = 4) :: g, j, i, ni, l, li
    real, dimension(3) :: elem_coord
    real,dimension(3, 3) :: K_integral, M_integral
    data M_integral / 4 , 2 , -1 , 2 , 16 , 2 , -1 , 2 ,4 / 
    data K_integral / 37, 44, -7, -32, 64, 32, -7, 44, 37 /
    real, dimension(4)  :: gauspt, gauswt
    data gauspt /0.8611363116, 0.3399810435, -0.3399810435, -0.8611363116 /  
    data gauswt / 0.347854851 ,  0.6521451548, 0.6521451548, 0.347854851 /
    real  :: xi, wt, cnst, h, s, s2 
!---Local element coordinate

    heat_elem_matrix_M = 0.0
    heat_elem_matrix_K = 0.0   
    heat_elem_vector_F = 0.0

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
               ! Calculate shape functions at gauss pt
        call inter_shape_fcns(xi,elem_coord,h)
        ! Evaluate material properties at gauss pt, first get temp @ gauss pt
        do i=1, nodes_per_elem
            T = T + shape_fcn(i)*previous_elem_soln_vec(n)
        end do          
 
        cnst = g_jacobian*wt
        do i=1, nodes_per_elem
        !---Setup local coordinates
            ! Get length of element
            ! Populate Kij, Mij, Fij, Qij
            do j = 1, nodes_per_elem
                s = shape_fcn(i)*shape_fcn(j)*cnst
                s2 = (der_shape_fcn(i)*der_shap_fcn(j))*(1/shape_fcn(j))
                heat_elem_matrix_M(i,j) = heat_elem_matrix_M(i,j) + s 
                heat_elem_matrix_K(i,j) = heat_elem_matrix_K(i,j) +  
            end do ! End loop over j matrix entries

        end do ! End loop over i matrix entries

    end do ! end do over gauss pts 

!-------------------------------------------------------------------------------
    if (DEBUG .eqv. .TRUE. ) then
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
        write(outfile_unit,fmt='(a)'),' '
        !---End matrix write out
    end if
!------------------------------------------------------------------------------
 
end 
