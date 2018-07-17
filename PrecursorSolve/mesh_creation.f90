!*****************************************************************************80
!
!! Sets up 1D mesh assuming quadratic functions  
!
!  Discussion:
!          Creates connectivity mesh assuming 1D discontinuous elements in a line
!     	       
!  Parameters:
!  	      conn_matrix(max # elements, nodes per element) 
! 	      global_coord(max # elements)

subroutine mesh_creation ( )
!
USE parameters_fe  

implicit none

    integer :: i,j,ii
    real :: temp

!---allocate arrays
    allocate( conn_matrix(num_elem, nodes_per_elem) , global_coord(matrix_length) )  
    
!---setup connectivity matrix
    do i=1, num_elem
        do j=1,nodes_per_elem
            conn_matrix(i,j) = j + (i-1)*nodes_per_elem
        end do 
    end do
    
!---setup global coordinate array 
    global_coord(1) = 0.0 
    do i = 2, nodes_per_elem
        global_coord(i) = global_coord(i-1) + 0.5*elem_lengths(1)
    end do
    
    do i=2, num_elem
        temp = global_coord( (i-1)* 3)
        do j = 1, nodes_per_elem
            ii = j + (i-1)*nodes_per_elem
            if ( j .eq. 1) then
                global_coord(ii) = temp 
            else
                global_coord(ii) = global_coord(ii-1) + 0.5*elem_lengths(i)
            end if
        end do
    end do

!-------------------------------------------------------------------------------
!---Write to outfile
    write(outfile_unit,fmt='(a19)'),'Connectivity Matrix'
    do j=1, num_elem
           write(outfile_unit,fmt='(a,1I2,a4,4I6)') 'Element:',j,' -->', &
                (conn_matrix(j,i),i=1,nodes_per_elem)             
    end do
    write(outfile_unit,fmt='(a)'),' ' 
    write(outfile_unit,fmt='(a24)'),'Global Coordinate Matrix'
    do j=1, matrix_length 
           write(outfile_unit,fmt='(a5,1I2,a,f8.3)') 'Node:',j,' x-coord -->', &
                global_coord(j)             
    end do

end
