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
    real :: initial_cond 
    real :: beg_temp, end_temp, heat_flux_bc
    !allocate(final_global_matrix_K(matrix_length, matrix_length), & 
    !         final_global_vec_f(matrix_length) )
  
 
    write(outfile_unit,fmt='(a)'), ' ' 
    write(outfile_unit,fmt='(a)'),'Global Assembled K matrix before B.C. imposed  '
    do j=1,matrix_length 
           write(outfile_unit,fmt='(8es14.3)') &
                (global_matrix_K(j,i) ,i=1,matrix_length)             
    end do

    !global_vec_q(1) = global_vec_q(1) + 200
    beg_temp = 500
    end_temp = 300
    heat_flux_bc = 50 
    ! Impose boundary condition of set value at the end  
    !global_vec_q(1) = beg_temp 
    global_vec_q(matrix_length) = end_temp 
    
    ! Set main diagonal
    !global_matrix_K(1,1) = 1.0
    global_matrix_K(matrix_length,matrix_length) = 1.0
   
    ! Impose B.C. at beginning, steady heat flux
    
    !! combine first colum with beg B.C. to add to RHS vector
    global_vec_q(1) = global_vec_q(1) + heat_flux_bc
    
    ! Set all other values along last colum = to zero except the main diagonal value
    do j = 1, matrix_length-1
        global_matrix_K(j, matrix_length) = 0.0 
    end do
    
    !do i = 2, matrix_length
    !    global_matrix_K(i,1) = 0.0
    !end do
    
    write(outfile_unit,fmt='(a)'), ' ' 
    write(outfile_unit,fmt='(a)'),'Global Assembled K matrix after B.C. imposed  '
    do j=1,matrix_length 
           write(outfile_unit,fmt='(8es14.3)') &
                (global_matrix_K(j,i) ,i=1,matrix_length)             
    end do

!   Print out matrix before B.C. applied
    if (DEBUG .eqv. .TRUE. ) then
!        write(outfile_unit,fmt='(a)'), ' ' 
!        write(outfile_unit,fmt='(a)'),'Global Matrix - steady state before B.C. applied: '
!        do j=1, matrix_length
!               write(outfile_unit,fmt='(12es14.6)') &
!                    (global_matrix_K(j,i) ,i=1,matrix_length)             
!        end do
!
        write(outfile_unit,fmt='(a)'), ' '        
        write(outfile_unit,fmt='(a)'), 'Global vector source Q - steady state   '
        write(outfile_unit,fmt='(12es14.6)') &
                    (global_vec_q(i) ,i=1,matrix_length)
    end if

! This is now handled in assembly


!    
!    global_vec_f(num_elem + 1) = global_vec_f(1) + global_vec_f(num_elem + 1)
!
!!   Add frst column to last - sweep out first column
!    do i = 1, matrix_length   
!        global_matrix_K(i, matrix_length) = global_matrix_K(i,1) + global_matrix_K(i,matrix_length)
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
