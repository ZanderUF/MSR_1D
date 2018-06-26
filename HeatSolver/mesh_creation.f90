!*****************************************************************************80
!
!! Sets up 1D mesh assuming quadratic functions  
!
!  Discussion:
!             Creates connectivity mesh assuming 1D and elements in a line
!     	       
!  Parameters:
!  	      conn_matrix(max # elements, nodes per element) 
! 	      global_coord(max # elements)

subroutine mesh_creation ( )
!
USE parameters_fe  

implicit none

    integer :: i,j,ii

    ! allocate arrays

    allocate( conn_matrix(num_elem,nodes_per_elem) , global_coord(max_num_nodes) )  
    
    ! setup connectivity matrix
    do i=1, nodes_per_elem
       conn_matrix(1,i) = i  
    end do
    do j=2, num_elem
        do i=1,nodes_per_elem
            conn_matrix(j,i) = conn_matrix(j-1,i) + nodes_per_elem-1
        end do 
    end do
    ! Account for periodic boundary condition, get mesh value and assign to front
    !conn_matrix(num_elem,nodes_per_elem) = 1   
    
    ! setup global coordinate array 
    global_coord(1) = 0.0 
    do i=1, num_elem
        ii=2*i
        global_coord(ii)   = global_coord(ii-1) + 0.5*elem_lengths(i)
        global_coord(ii+1) = global_coord(ii-1) + elem_lengths(i) 
    end do
    
    ! Write to outfile
    write(outfile_unit,fmt='(a19)'),'Connectivity Matrix'
    do j=1, num_elem
           write(outfile_unit,fmt='(a,1I2,a4,4I6)') 'Element:',j,' -->', &
                (conn_matrix(j,i),i=1,nodes_per_elem)             
    end do
    write(outfile_unit,fmt='(a)'),' ' 
    write(outfile_unit,fmt='(a24)'),'Global Coordinate Matrix'
    do j=1, max_num_nodes 
           write(outfile_unit,fmt='(a5,1I2,a4,f8.3)') 'Node:',j,' -->', &
                global_coord(j)             
    end do

end
