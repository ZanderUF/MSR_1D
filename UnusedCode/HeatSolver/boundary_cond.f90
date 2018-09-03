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
    real :: beg_temp, end_temp, heat_flux_bc,source_first,source_last
    !allocate(final_global_matrix_A(matrix_length, matrix_length), & 
    !         final_global_vec_f(matrix_length) )
    real :: length, temp_val,temp_source, L, H, R_first, R_last, &
            temp_before_T
    real, dimension(3) :: G_last,G_first,temp_vec
    ! Inversion routine parameters
    integer :: lda, info, lwork
    integer, dimension(3) :: ipiv
    real, dimension(3) :: work
    length = 3 
    R_first = 0.0
    R_last = 0.0
    temp_vec = 0.0
    temp_val = 0.0
    temp_source = 0.0
    temp_before_T = 0.0
    G_first = 0.0
    G_last = 0.0
    ipiv =  0
    work =  0.0
    lda =   length
    lwork = length

    write(outfile_unit,fmt='(a)'), ' ' 
    write(outfile_unit,fmt='(a)'),'Global Assembled K matrix before B.C. imposed  '
    do j=1,matrix_length 
           write(outfile_unit,fmt='(12es14.3)') &
                (global_matrix_A(j,i) ,i=1,matrix_length)             
    end do

!---Invert first element A matrix
    call dgetrf ( length, length, heat_elem1_matrix_A_s2 , lda, ipiv, info )
!---Compute the inverse matrix.
    call dgetri ( length, heat_elem1_matrix_A_s2, lda, ipiv, work, lwork, info ) 

!---For the FIRST element
!---Calculate (2 + M_e,1^T A^-1 M_e,1)^-1
    do i = 1, nodes_per_elem
        do j = 1, nodes_per_elem
            temp_vec(i) = temp_vec(i) + &
                heat_elem1_vec_M_s1(j)*heat_elem1_matrix_A_s2(i,j) 
        end do
    end do
    
    do i = 1, nodes_per_elem
        temp_val      = temp_val      + temp_vec(i)*heat_elem1_vec_M_s1(i) 
        temp_source   = temp_source   + temp_vec(i)*heat_elem1_vec_f(i)
        do j = 1, nodes_per_elem
            G_first(i) = G_first(i) + temp_vec(j)*heat_elem1_D_s2(i,j)
        end do
    end do
    
    H = 1.0 / ( 2.0 + temp_val)
    R_first = H*(temp_val - 2.0)
    source_first = H*temp_source 

!   G_first is the evaluated value b4 the Temperature for the partial current equation
    do i = 1, nodes_per_elem
        G_first(i) = H*G_first(i)  
    end do

!   reinitialize for next inversion
    ipiv =  0
    work =  0.0

!   invert last element A matrix
    call dgetrf ( length, length, heat_last_elem_matrix_A_s1, lda, ipiv, info )

!   Compute the inverse matrix.
    call dgetri (length, heat_last_elem_matrix_A_s1, lda, ipiv, work, lwork, info ) 

!   reset temp values
    temp_vec(:)   = 0.0
    temp_val      = 0.0
    temp_source   = 0.0
    temp_before_T = 0.0
!   For the LAST element (last_elem)
!   Calculate (2 - M_e,2^T A^-1 M_e,2)^-1
    do i = 1, nodes_per_elem
        do j = 1, nodes_per_elem
            temp_vec(i) = temp_vec(i) + &
                heat_last_elem_vec_M_s2(j)*heat_last_elem_matrix_A_s1(i,j)
        end do
    end do
    do i = 1, nodes_per_elem
        temp_val      = temp_val     + temp_vec(i)*heat_last_elem_vec_M_s2(i)
        temp_source   = temp_source  + temp_vec(i)*heat_last_elem_vec_f(i)
        do j = 1, nodes_per_elem
            G_last(i) = G_last(i) + temp_vec(j)*heat_last_elem_D_s1(i,j)
        end do
    end do

    do i = 1, nodes_per_elem
        G_last(i)= L*G_last(i)
    end do
    L = 1.0 / (2.0 - temp_val)
    R_last = L*(temp_val + 2.0)
    source_last = L*temp_source 

!   Add connection with first and last nodal temperature values with partialcurrents

    global_matrix_A(1, matrix_length-1) =  heat_elem1_vec_M_s1(1)
    global_matrix_A(1, matrix_length)   = - heat_elem1_vec_M_s1(1)
    
    global_matrix_A(matrix_length-2, matrix_length-1) =   heat_last_elem_vec_M_s2(3)
    global_matrix_A(matrix_length-2, matrix_length)   = - heat_last_elem_vec_M_s2(3)
    
!   Construct partial current equation in global matrix
    global_matrix_A(matrix_length-1,matrix_length-1) = 1.0
    global_matrix_A(matrix_length-1,matrix_length) = -R_first

    global_matrix_A(matrix_length,matrix_length-1) = R_last
    global_matrix_A(matrix_length,matrix_length)   = 1.0
    
    global_matrix_A(matrix_length-1, 1) = -G_first(3)
    global_matrix_A(matrix_length,matrix_length-2) = G_last(1)

    global_vec_q(matrix_length - 1) = source_first
    global_vec_q(matrix_length)     = source_last

!---Put fixed loss of heat at a node spot
    global_vec_q(matrix_length-4) = global_vec_q(matrix_length-4) - 5E-3

    write(outfile_unit,fmt='(a)'), ' ' 
    write(outfile_unit,fmt='(a)'),'Global Assembled A matrix after B.C. imposed  '
    do j=1,matrix_length 
           write(outfile_unit,fmt='(12es14.3)') &
                (global_matrix_A(j,i) ,i=1,matrix_length)             
    end do

!   Print out matrix before B.C. applied
    if (DEBUG .eqv. .TRUE. ) then
        write(outfile_unit,fmt='(a)'), ' '        
        write(outfile_unit,fmt='(a)'), 'Global vector source Q - steady state   '
        write(outfile_unit,fmt='(12es14.6)') &
                    (global_vec_q(i) ,i=1,matrix_length)
    end if

end
