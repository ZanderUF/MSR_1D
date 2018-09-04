!*****************************************************************************80
!
!! Assemble the elemental to solve for coefficients the problem 
!  Discussion:
!             As the elemental matrices are created we add them to the global 
!
!  Input: 
!           isotope - material isotop identifier
!           delay_group - delayed neutron group identifier
!           n - element number
!
!  Output:
! 

subroutine assemble_matrix_ss (isotope,delay_group,n)
!---
    USE parameters_fe  

    implicit none
!---Dummy variables
    integer,intent(in) :: isotope
    integer,intent(in) :: delay_group
    integer :: n
!---Local variables
    integer :: i, j, ii, jj, nr,nc, length   
    
!---Inversion routine parameters
    integer :: lda, info, lwork
    integer, dimension(3) :: ipiv
    real, dimension(3) :: work

    length = 3
!---Initialize
    ipiv =  0
    work =  0.0
    lda =   length
    lwork = length
    
    elem_vec_q = 0.0
    elem_matrix_G = 0.0
    elem_vec_w_left_face = 0.0
    matrix_W_right_face = 0.0 
    matrix_W_left_face = 0.0

    !---Create source vector 'q', and W - 
    do i = 1, nodes_per_elem
        do j = 1, nodes_per_elem
            elem_vec_q(i) = elem_vec_q(i) + &
                            elem_matrix_A(i,j)*spatial_power_fcn(n,j)*power_amplitude_prev
            !---Applies for all elements except the first one
            if(n > 1) then !--- n - element #
                !---Grab previous precursor conc. + velocity at 
                !---rhs of previous element
                 matrix_W_right_face(i,j) = velocity_soln_prev(n,i)*&
                                            interp_fcn_rhs(i)*interp_fcn_rhs(j)  
                 matrix_W_left_face(i,j)  = velocity_soln_prev(n,i)*&
                                            interp_fcn_lhs(i)*interp_fcn_lhs(j)
            else!---First element case, need to connect with end element 
                 matrix_W_right_face(i,j) = velocity_soln_prev(num_elem,i)*&
                                            interp_fcn_rhs(i)*interp_fcn_rhs(j)  
                 matrix_W_left_face(i,j)  = velocity_soln_prev(num_elem,i)*&
                                            interp_fcn_lhs(i)*interp_fcn_lhs(j)
            end if!---End flux calculation
        end do
    end do

    do i = 1, nodes_per_elem
        !---Calculate q vector
        elem_vec_q_final(isotope,delay_group,i) = &
                  (beta_i_mat(isotope,delay_group)/gen_time)*elem_vec_q(i) 
        do j = 1, nodes_per_elem
            if (n > 1) then
                elem_vec_w_left_face(i) = elem_vec_w_left_face(i) + &
                                          matrix_W_left_face(i,j)*&
                                          precursor_soln_new(isotope,delay_group, n-1,3)
            else
                elem_vec_w_left_face(i) = elem_vec_w_left_face(i) + &
                                          matrix_W_left_face(i,j)*&
                                          precursor_soln_new(isotope,delay_group,num_elem,3)
            end if
            elem_matrix_G(i,j) = -elem_matrix_U(i,j) + &
                                 (lamda_i_mat(isotope,delay_group))*elem_matrix_A(i,j) + &
                                 matrix_W_right_face(i,j)
        end do
    end do   

!****************************************************************

if(DEBUG .eqv. .TRUE.) then
    write(outfile_unit,fmt='(a)'), ' '
    write(outfile_unit,fmt='(a,1I2)'),'G Matrix | element --> ',n
    do i=1,nodes_per_elem 
          write(outfile_unit,fmt='(12es14.3)') &
               (elem_matrix_G(i,j),j=1,nodes_per_elem)             
    end do

    
    write(outfile_unit,fmt='(a)'), ' ' 
    write(outfile_unit,fmt='(a,1I3)'),'left face vector | element --> ', n
    do j=1,nodes_per_elem 
          write(outfile_unit,fmt='(12es14.3)') elem_vec_w_left_face(j)             
    end do
    
    write(outfile_unit,fmt='(a)'),' '
    write(outfile_unit,fmt='(a,1I2)'),'beta/gen time *{q} element source vector | element --> ',n
    write(outfile_unit,fmt='(12es14.3)') (elem_vec_q_final(isotope,delay_group,i),i=1,nodes_per_elem)             
end if

end subroutine assemble_matrix_ss
