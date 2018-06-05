
!*****************************************************************************80
!
!! Sets up the elemental matrices  
!
!  Discussion:
!             Generates the elemental matrices K_ij, M_ij, S_ij 
!     	      calculates coefficients in this suboutine 
!  Parameters:

subroutine element_matrix ( )
!
  USE parameters_fe  

  implicit none
  
    integer ( kind = 4) :: n, i, ni, l, li
    real, allocatable :: elem_coord(:)

!---Local element coordinate
    allocate(elem_coord(node_per_elem) )

    ! loop over all elements in the mesh
    do n=1, num_elem
        l=0
        ! Setup local coordinates
        do i=1, nodes_per_elem
            ni = conn_matrix(n,i)
            elem_coord(i) = global_coord(ni) 
        ! Include some TD stuff here
            li = (ni - 1)*ndf
            do j=1,ndf
                li = li+1
                l  = l+1
                elem_soln_vec
            end do
        end do

    end do 
end 
