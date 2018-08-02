!*****************************************************************************80
!
!! Sets up 1D mesh assuming quadratic functions  
!
!  Discussion:
!          Creates connectivity mesh assuming 1D discontinuous elements in a line
!     	       
!  Parameters:
!  	      conn_matrix(max # elements, nodes per element) 
! 	      global_coord(# elements, nodes per element)

subroutine mesh_creation ( )
!
USE parameters_fe  

implicit none

    integer :: i,j,ii
    real :: temp

!---allocate arrays
    allocate( conn_matrix(num_elem, nodes_per_elem) , &
              global_coord(num_elem,nodes_per_elem) )  

!---setup connectivity matrix
    do i=1, num_elem
        do j=1,nodes_per_elem
            conn_matrix(i,j) = j + (i-1)*nodes_per_elem
        end do 
    end do

!---Setup element length array. assume all elements have the same size
    !elem_size = 2.0
    do i = 1, num_elem
        elem_lengths(i) = elem_size 
    end do

!---Setup global coordinate array 
    global_coord(1,1) = 0
    global_coord(1,2) = global_coord(1,1) + 0.5*elem_lengths(1)
    global_coord(1,3) = global_coord(1,2) + 0.5*elem_lengths(1)

    do i = 2, num_elem
        do j = 1, nodes_per_elem
            if ( j .eq. 1) then
                global_coord(i,j) = global_coord(i-1,3) 
            else
                global_coord(i,j) = global_coord(i,j-1) + 0.5*elem_lengths(i)
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
    do j=1, num_elem 
           write(outfile_unit,fmt='(12es14.3)')  &
                (global_coord(j,i),i=1,nodes_per_elem)              
    end do

end
