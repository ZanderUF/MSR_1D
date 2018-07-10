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
    real :: length, temp_val,temp_source, L, H, R_first, R_last,G_first, &
            G_last, temp_before_T
    real, dimension(4) :: temp_vec
    ! Inversion routine parameters
    integer :: lda, info, lwork
    integer, dimension(4) :: ipiv
    real, dimension(4) :: work
    length = 4
    R_first = 0.0
    R_last = 0.0
    temp_vec = 0.0
    temp_val = 0.0
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
           write(outfile_unit,fmt='(8es14.3)') &
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
            temp_vec(j) = temp_vec(j) + &
                heat_elem1_vec_M_s1(i)*heat_elem1_matrix_A_s2(i,j) 
        end do
    end do
    
    do i = 1, nodes_per_elem
        temp_val      = temp_val      + temp_vec(i)*heat_elem1_vec_M_s1(i) 
        temp_source   = temp_source   + temp_vec(i)*heat_elem1_vec_f(i)
        print *,'temp_vec_',temp_vec(i)
        print *,'heat elem',heat_elem1_D_s2
        temp_before_T = temp_before_T + temp_vec(i)*heat_elem1_D_s2
        print *,'temp_before_T',temp_before_T,' ',i
    end do


    H = 1.0 / ( 2.0 + temp_val)
    R_first = H*(temp_val - 2.0)
    source_first = H*temp_source 
!   G_first is the evaluated value b4 the Temperature for the partial current equation
    G_first = H*temp_before_T  

!   reinitialize for next inversion
    ipiv =  0
    work =  0.0

!   invert last element A matrix
    call dgetrf ( length, length, heat_elem2_matrix_A_s1, lda, ipiv, info )

!   Compute the inverse matrix.
    call dgetri (length, heat_elem2_matrix_A_s1, lda, ipiv, work, lwork, info ) 

!   reset temp values
    temp_vec(:)      = 0.0
    temp_val      = 0.0
    temp_source   = 0.0
    temp_before_T = 0.0
!   For the LAST element
!   Calculate (2 - M_e,2^T A^-1 M_e,2)^-1
    do i = 1, nodes_per_elem
        do j = 1, nodes_per_elem
            temp_vec(i) = temp_vec(i) + &
                heat_elem2_vec_M_s2(i)*heat_elem2_matrix_A_s1(i,j)
        end do
    end do
    do i = 1, nodes_per_elem
        temp_val      = temp_val      + temp_vec(i)*heat_elem2_vec_M_s2(i)
        temp_source   = temp_source   + temp_vec(i)*heat_elem2_vec_f(i)
        temp_before_T = temp_before_T + temp_vec(i)*heat_elem2_D_s1
    end do
    
    L = 1.0 / (2.0 - temp_val)
    R_last = L*(temp_val + 2.0)
    source_last = L*temp_source 
    G_last = L*temp_before_T

!   Add connection with first and last nodal temperature values with partialcurrents

    global_matrix_A(1, matrix_length-1) = - heat_elem1_vec_M_s1(1)
    global_matrix_A(1, matrix_length)   =  heat_elem1_vec_M_s1(1)
    
    global_matrix_A(matrix_length-2, matrix_length-1) =   heat_elem2_vec_M_s2(3)
    global_matrix_A(matrix_length-2, matrix_length)   = - heat_elem2_vec_M_s2(3)
    
!   Construct partial current equation in global matrix
    global_matrix_A(matrix_length-1,matrix_length-1) = 1.0
    global_matrix_A(matrix_length-1,matrix_length) = -R_first

    global_matrix_A(matrix_length,matrix_length-1) = R_last
    global_matrix_A(matrix_length,matrix_length)   = 1.0

    print *,'g first',G_first
    print *,'g last ',G_last

    global_matrix_A(matrix_length-1, 1) = G_first
    global_matrix_A(matrix_length,matrix_length-2) = G_last

!   Old - B.C. for non periodic cases

    !beg_temp = 500
    !end_temp = 300
    !heat_flux_bc = 50 
    !! Impose boundary condition of set value at the end  
    !global_vec_q(matrix_length) = end_temp 
    !
    !! Set main diagonal
    !global_matrix_A(matrix_length,matrix_length) = 1.0
   
    !! Impose B.C. at beginning, steady heat flux
    !
    !!! combine first colum with beg B.C. to add to RHS vector
    !global_vec_q(1) = global_vec_q(1) + heat_flux_bc
    !
    !! Set all other values along last colum = to zero except the main diagonal value
    !do j = 1, matrix_length-1
    !    global_matrix_A(j, matrix_length) = 0.0 
    !end do
    
    !do i = 2, matrix_length
    !    global_matrix_A(i,1) = 0.0
    !end do
    
    write(outfile_unit,fmt='(a)'), ' ' 
    write(outfile_unit,fmt='(a)'),'Global Assembled A matrix after B.C. imposed  '
    do j=1,matrix_length 
           write(outfile_unit,fmt='(8es14.3)') &
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
