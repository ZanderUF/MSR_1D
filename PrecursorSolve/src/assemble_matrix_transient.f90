!*****************************************************************************80
!
!! Assemble the elemental to solve for coefficients the transient 
!  Discussion:
!             As the elemental matrices are created we add them to the global 
!
!  Input: n - element number
!
!  Output:
! 

subroutine assemble_matrix_transient (isotope,delay_group,n)
   
    USE flags_M
    USE element_matrices_M
    USE global_parameters_M
    USE material_info_M
    USE mesh_info_M
    USE solution_vectors_M
    USE mesh_info_M

    implicit none
!---Dummy variables
    integer,intent(in) :: isotope
    integer,intent(in) :: delay_group
    integer,intent(in) :: n
!---Local variables
    integer :: i, j, ii, jj, nr,nc, ncl,length   
    double precision, dimension(3,3) :: inverse_matrix,temp_matrix
    double precision, dimension(3) :: elem_vec_q_times_beta_lambda, U_times_soln_vec, &
                          A_times_lambda_times_soln_vec, W_left_times_upwind_soln, &
                          W_right_times_soln

!---Initialize
    elem_vec_q_times_beta_lambda(:)  = 0.0
    U_times_soln_vec (:)             = 0.0
    A_times_lambda_times_soln_vec(:) = 0.0 
    W_left_times_upwind_soln         = 0.0
    W_right_times_soln(:)            = 0.0
    RHS_transient_final_vec          = 0.0

    do i = 1, nodes_per_elem
        elem_vec_q_times_beta_lambda(i) = &
        ( beta_i_mat(isotope,delay_group)/gen_time )*&
          elem_vec_q(i) 
        
        do j = 1, nodes_per_elem
            
            U_times_soln_vec(i) = U_times_soln_vec(i) + &
                                  elem_matrix_U(i,j)*&
                                  precursor_soln_prev(isotope,delay_group,n,j)                 
            A_times_lambda_times_soln_vec(i) = A_times_lambda_times_soln_vec(i)+&
                                         elem_matrix_A(i,j)*&
                                         lamda_i_mat(isotope,delay_group)*&
                                         precursor_soln_prev(isotope,delay_group,n,j)
            
            if ( n > 1) then !---Always get the previous element
                W_right_times_soln(i) = W_right_times_soln(i) + &
                                    matrix_W_right_face(i,j)*&
                                    precursor_soln_prev(isotope,delay_group,n,j)

                W_left_times_upwind_soln(i) = &
                    W_left_times_upwind_soln(i) + &
                    matrix_W_left_face(i,j)*&
                    precursor_soln_new(isotope,delay_group,n-1,3)
            else !--Account for connection end of domain to beginning
                W_right_times_soln(i) = W_right_times_soln(i) + & 
                                    matrix_W_right_face(i,j)*&
                                    precursor_soln_prev(isotope,delay_group,n,j)
                W_left_times_upwind_soln(i) = W_left_times_upwind_soln(i) + & 
                    matrix_W_left_face(i,j)*&
                    precursor_soln_new(isotope, delay_group, num_elem, 3)
 
            end if
        end do 
    end do
    
    !---Calcualte the RHS vector ** F(u)
    do i = 1, nodes_per_elem
        !print *,' U times',U_times_soln_vec(i)
        !print *,' A lambda',A_times_lambda_times_soln_vec(i)
        !print *,' elem vec',elem_vec_q_times_beta_lambda(i)
        !print *,' W right', W_right_times_soln(i)
        !print *,' W left', W_left_times_upwind_soln(i)

        RHS_transient_final_vec(i) = U_times_soln_vec(i) - &
                           A_times_lambda_times_soln_vec(i) + &
                           elem_vec_q_times_beta_lambda(i) -  &
                           W_right_times_soln(i) + &
                           W_left_times_upwind_soln(i) 
        !print *,'RHS vec',RHS_transient_final_vec(i)
        !print *,' ' 
    end do 
!****************************************************************
!---Write out    
    if (DEBUG .eqv. .TRUE.) then
                write(outfile_unit,fmt='(a)'), ' ' 
        write(outfile_unit,fmt='(a,1I3)'),' [U]*{c_e} | element --> ', n
        do j=1,nodes_per_elem 
              write(outfile_unit,fmt='(16es16.8)')U_times_soln_vec(j)              
        end do
        write(outfile_unit,fmt='(a)'), ' ' 
        write(outfile_unit,fmt='(a,1I3)'),' (lambda)*[A]*{c_e} | element --> ', n
        do i=1,nodes_per_elem 
              write(outfile_unit,fmt='(12es16.10)')A_times_lambda_times_soln_vec(i)
        end do
        
        write(outfile_unit,fmt='(a)'), ' ' 
        write(outfile_unit,fmt='(a,1I3)'),' (beta/Cap lambda)*{q_e} |&
                                            element --> ', n
        do j=1,nodes_per_elem 
              write(outfile_unit,fmt='(12es16.10)')elem_vec_q_times_beta_lambda(j)             
        end do
        
        write(outfile_unit,fmt='(a)'), ' '
        write(outfile_unit,fmt='(a,1I2)'),' [W_R]*{c_e} | element --> ',n
        do j=1,nodes_per_elem 
              write(outfile_unit,fmt='(12es16.10)')W_right_times_soln(j) 
        end do
           
        write(outfile_unit,fmt='(a)'),' '
        write(outfile_unit,fmt='(a,1I2)'),' [W_L]*{c_(e-1)} | element --> ',n
        do j=1,nodes_per_elem
            write(outfile_unit,fmt='(12es16.10)') W_left_times_upwind_soln(j)          
        end do
        
        write(outfile_unit,fmt='(a)'),' '
        write(outfile_unit,fmt='(a,1I2)'),' RHS F(u) | element --> ',n
        do j=1,nodes_per_elem
            write(outfile_unit,fmt='(16es16.8)') RHS_transient_final_vec(j)          
        end do
        write(outfile_unit,fmt='(a)'),' '
 
end if
end 
