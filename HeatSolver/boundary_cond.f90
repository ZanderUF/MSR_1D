! Evaluate boundary condition 
!
! Apply periodic boundary conditions to the global system
! Will reduce the size of the matrix from (2E + 1) x (2E +1) to (2E x 2E)
!
! Input:

! Output:
! 
subroutine boundary_cond( )

    USE parameters_fe

    implicit none

! Dummy

! Local
    integer :: i,j 

    allocate(final_global_matrix_K(2*num_elem,2*num_elem), & 
             final_global_vec_f(2*num_elem) )

!   Add last row to first
    do j = 1, 2*num_elem + 1
        global_matrix_K(1, j) = global_matrix_K(1,j) + global_matrix_K(2*num_elem+1,j)
    end do
    
    global_vec_f(1) = global_vec_f(1) + global_vec_f(num_elem + 1)

!   Add last column to first - sweep out last column
    do i = 1, 2*num_elem +   1   
        global_matrix_K(i, 1) = global_matrix_K(i,1) + global_matrix_K(i,2*num_elem+1)
    end do
!   Remove last column and row - apply periodic B.C.
    do i = 1, 2*num_elem
        final_global_vec_f(i)=global_vec_f(i)
        do j = 1, 2*num_elem
            final_global_matrix_K(i,j) = global_matrix_K(i,j)
        end do
    end do

!  No need for elemental matrices after they have been placed in the global one
   if (DEBUG .eqv. .TRUE. ) then
        write(outfile_unit,fmt='(a)'), ' ' 
        write(outfile_unit,fmt='(a)'),'Global Matrix - steady state after B.C. applied: '
        do j=1,2*num_elem
               write(outfile_unit,fmt='(12es14.6)') &
                    (final_global_matrix_K(j,i) ,i=1,2*num_elem)             
        end do

        write(outfile_unit,fmt='(a)'), ' '        
        write(outfile_unit,fmt='(a)'), 'Global vector source f - steady state after B.C. applied:  '
        write(outfile_unit,fmt='(12es14.6)') &
                    (final_global_vec_f(i) ,i=1,2*num_elem)
   end if

end
