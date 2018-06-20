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
    integer :: i,j, vec_length

    vec_length = 2*num_elem 
    allocate(final_global_matrix_K(vec_length, vec_length), & 
             final_global_vec_f(vec_length) )

!   Print out matrix before B.C. applied
!    if (DEBUG .eqv. .TRUE. ) then
!        write(outfile_unit,fmt='(a)'), ' ' 
!        write(outfile_unit,fmt='(a)'),'Global Matrix - steady state before B.C. applied: '
!        do j=1, vec_length
!               write(outfile_unit,fmt='(12es14.6)') &
!                    (global_matrix_K(j,i) ,i=1,vec_length)             
!        end do
!
!        write(outfile_unit,fmt='(a)'), ' '        
!        write(outfile_unit,fmt='(a)'), 'Global vector source f - steady state before B.C. applied:  '
!        write(outfile_unit,fmt='(12es14.6)') &
!                    (global_vec_f(i) ,i=1,vec_length)
!    end if

! This is now handled in assembly

!   Add last first to last
!    do j = 1, vec_length 
!        global_matrix_K(vec_length, j) = global_matrix_K(1,j) + global_matrix_K(vec_length,j)
!    end do
!    
!    global_vec_f(num_elem + 1) = global_vec_f(1) + global_vec_f(num_elem + 1)
!
!!   Add frst column to last - sweep out first column
!    do i = 1, vec_length   
!        global_matrix_K(i, vec_length) = global_matrix_K(i,1) + global_matrix_K(i,vec_length)
!    end do
!!   Remove first column and row - apply periodic B.C.
!    do i = 1, 2*num_elem
!        final_global_vec_f(i)=global_vec_f(i+1)
!        do j = 1, 2*num_elem
!            final_global_matrix_K(i,j) = global_matrix_K(i+1,j+1)
!        end do
!    end do
!
!!  No need for elemental matrices after they have been placed in the global one
!   if (DEBUG .eqv. .TRUE. ) then
!        write(outfile_unit,fmt='(a)'), ' ' 
!        write(outfile_unit,fmt='(a)'),'Global Matrix - steady state after B.C. applied: '
!        do j=1,2*num_elem
!               write(outfile_unit,fmt='(12es14.6)') &
!                    (final_global_matrix_K(j,i) ,i=1,2*num_elem)             
!        end do
!
!        write(outfile_unit,fmt='(a)'), ' '        
!        write(outfile_unit,fmt='(a)'), 'Global vector source f - steady state after B.C. applied:  '
!        write(outfile_unit,fmt='(12es14.6)') &
!                    (final_global_vec_f(i) ,i=1,2*num_elem)
!   end if

end
