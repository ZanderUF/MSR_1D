
!*****************************************************************************80
!
!! Sets up the elemental matrices for the heat equation  
!
!  Discussion:
!             Generates the elemental matrices K_ij, M_ij, S_ij 
!     	      calculates coefficients in this suboutine 
!  Parameters:

subroutine element_matrix_heat ( )
!
  USE parameters_fe  

  implicit none
  
    integer ( kind = 4) :: j, n, i, ni, l, li
    real, allocatable :: elem_coord(:)
    real,dimension(3, 3) :: K_integral, M_integral
    data K_integral / 4 , 2 , -1 , 2 , 16 , 2 , -1 , 2 ,4 / 
    data M_integral / 37, 44, -7, -32, 64, 32, -7, 44, 37 /

!---Local element coordinate
    allocate(elem_coord(nodes_per_elem) )
    allocate(heat_elem_matrix_K( nodes_per_elem, nodes_per_elem),&
             heat_elem_matrix_M( nodes_per_elem, nodes_per_elem),&
             heat_elem_vector_F( nodes_per_elem) ) 
   
    ! loop over all elements in the mesh
    do n = 1 , 2 
        l=0
        ! Setup local coordinates
        do i=1, nodes_per_elem
            ni = conn_matrix(n,i) 
            elem_coord(i) = global_coord(ni) 
            ! Include initial values for TD problems 
            !li = (ni - 1)*ndf
            !! Apply initial conditions for TD problems
            !do j=1,ndf
            !    li = li+1
            !    l  = l+1
            !    elem_soln_vec
            !end do
            
            ! Get length of element
            h = elem_lengths(n)
            ! Populate Kij, Mij, Fij, Qij, using analytic evaluations of integrals
            do j = 1, nodes_per_elem
                heat_elem_matrix_K(i,j) =  (h/30.0)*K_integral(i,j)               
                heat_elem_matrix_M(i,j) =  (h/30.0)*M_integral(i,j)
            end do          
        end do
        if (DEBUG == .TRUE. ) then
            !---Write out each element matrix to the outfile 
            write(outfile_unit,fmt='(a24,1I2)'),'K element Matrix: ',n
            do j=1,nodes_per_elem 
                   write(outfile_unit,fmt='(12es14.6)') &
                        (heat_elem_matrix_K(j,i),i=1,nodes_per_elem)             
            end do
            write(outfile_unit,fmt='(a)'),' ' 
            write(outfile_unit,fmt='(a24,1I2)'),'M element Matrix: ',n
            do j=1,nodes_per_elem 
                   write(outfile_unit,fmt='(12es14.6)') &
                        (heat_elem_matrix_M(j,i),i=1,nodes_per_elem)             
            end do
            write(outfile_unit,fmt='(a)'),' '
            !---End matrix write out
        end if
        
    end do 

end 
