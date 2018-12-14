!---Solve transient temperature equation 
!
!---Input: n - element number 
!          nl_iter - nonlinear iteration counter
!---Output:
!

subroutine solve_temperature_euler(n,nl_iter)

    USE flags_M
    USE material_info_M
    USE mesh_info_M
    USE time_info_M
    USE global_parameters_M
    USE solution_vectors_M
    USE element_matrices_M

implicit none

!---Dummy
    integer, intent(in) :: n 
    integer, intent(in) :: nl_iter

!---Local
    integer  :: i, j, length
    real(dp) :: temperature_eval
    real(dp), dimension(3)   :: rhs_final_vec,elem_vec_w_left,&
                                A_inv_times_q_wl, elem_vec_q_temp,&
                                A_inv_times_q_wl_vec,A_inv_U_W_times_T_vec,&
                                U_W_times_T
    real(dp), dimension(3,3) ::  U_minus_W_right_mat 
    real(dp) :: volume, heat_capacity_eval, density_eval,& 
                salt_mass, q_prime, length_core
!---Inversion routine parameters
    integer :: lda, info, lwork
    integer,  dimension(3) :: ipiv
    real(dp), dimension(3) :: work

!---Size of matrix 3x3
    length = 3

!---Initialize inversion routine
    ipiv   = 0
    work   = 0.0
    lda    = length
    lwork  = length

!---Temperature SOLVE
    rhs_final_vec        = 0.0_dp
    elem_vec_q_temp      = 0.0_dp
    elem_vec_w_left      = 0.0_dp
    U_minus_W_right_mat  = 0.0_dp

!---[U - W_r]
    do i = 1, nodes_per_elem
        do j = 1, nodes_per_elem
            U_minus_W_right_mat(i,j) = -1.0*elem_matrix_U(i,j) + matrix_W_right_face(i,j) 
        end do
    end do

    U_W_times_T = 0.0_dp
!---[U - W_r]*T^k-1
    do i = 1, nodes_per_elem
        do j = 1, nodes_per_elem
            U_W_times_T(i) = U_W_times_T(i) + &
                             U_minus_W_right_mat(i,j)*temperature_soln_prev(n,i)
        end do
    end do

    U_W_times_T = matmul(U_minus_W_right_mat,temperature_soln_prev(n,:))

!---A^-1*[U - W_r]*T^k-1
    A_inv_U_W_times_T_vec = matmul(inverse_A_matrix,U_W_times_T)

!---Loop over all nodes in element
    do j = 1, nodes_per_elem
        temperature_eval = temperature_soln_prev(n,j)
        call density_corr(temperature_eval, density_eval) 
        call heat_capacity_corr(temperature_eval, heat_capacity_eval)
        q_prime = power_amplitude_prev*spatial_power_frac_fcn(n,j)/spatial_area_fcn(n,j)
        elem_vec_q_temp(j) = vol_int(j)*q_prime*(1.0_dp/(density_eval*heat_capacity_eval))

    end do

!----Evaluate W_l*T_e-1
    do i = 1, nodes_per_elem
        do j = 1, nodes_per_elem
            if(n > 1) then
                elem_vec_w_left(j) = elem_vec_w_left(j) + &
                                     matrix_W_left_face(i,j)*&
                                     temperature_soln_new(n-1,3)
            else
                elem_vec_w_left(j) = elem_vec_w_left(j) + &
                                     matrix_W_left_face(i,j)*&
                                     temperature_soln_new(n,3)
            end if
        end do
    end do 

!--- {q} + {w_l}
    do i = 1, nodes_per_elem
        rhs_final_vec(i) = elem_vec_q_temp(i) + elem_vec_w_left(i)
    end do

!---A^-1 *(q + w_l)
    A_inv_times_q_wl_vec = matmul(inverse_A_matrix,rhs_final_vec)

!---Forward Euler
    if(td_method_type == 0) then

    end if

!---Backward Euler
    if(td_method_type == 1) then
        do i = 1, nodes_per_elem
            
            temperature_soln_new(n,i) = temperature_soln_starting(n,i) + &
                                        delta_t*(A_inv_U_W_times_T_vec(i) - A_inv_times_q_wl_vec(i)) 

        end do
    end if

!---Inlet boundary conditions
    !if( n == Fuel_Inlet_start) then 
    !    do i = 1, nodes_per_elem
    !        temperature_soln_new(n,i) = 850.0_dp   
    !    end do 
    !end if

!---Boundary condition after heat exchanger
    if( n == Heat_Exchanger_End) then 
        do i = 1, nodes_per_elem
            temperature_soln_new(n,i) = 850.0_dp   
        end do 
    end if

!---Fix delta T across heat exchanger
    do j = 1, nodes_per_elem    
        !---Start heat exchanger
        if ( Heat_Exchanger_Start <= n  .and.  n < Heat_Exchanger_End) then
            temperature_soln_new(n,j) = (100.0_dp)/&
            (Heat_Exchanger_End - Heat_Exchanger_Start)* &
            (global_coord(Heat_Exchanger_End-1,3)  - global_coord(n,j) ) - 100.0_dp + &
            temperature_soln_new(Heat_Exchanger_Start-1,1)
        end if
    end do 

    Temperature_Reactivity_Feedback(n,:)=0.0 
    !---Calculate doppler feedback for current element
    !Temperature_Reactivity_Feedback(n,j) = spatial_doppler_fcn(n,j)* &
    !            ( temperature_soln_prev(n,j) - temperature_soln_starting(n,j) )
DEBUG = .TRUE.
if( DEBUG .eqv. .TRUE.) then
      
       write(outfile_unit,fmt='(a,4I6)'), 'nl iter ',nl_iter
       
       write(outfile_unit,fmt='(a)'), ' '
       write(outfile_unit,fmt='(a,4I6)'),' [W_L]*{T_e-1} | element --> ',n
       do j=1,nodes_per_elem 
             write(outfile_unit,fmt='(12es14.3)') &
                    elem_vec_w_left(j)
       end do

       write(outfile_unit,fmt='(a,4I6)'),' 1/Cp*rho*vol * {q} vector  --> ',n
             write(outfile_unit,fmt='(12es14.3)') &
                  (elem_vec_q_temp(i),i=1,nodes_per_elem)             
       
       write(outfile_unit,fmt='(a,4I6)'),'RHS vector  vector  --> ',n
             write(outfile_unit,fmt='(12es14.3)') &
                 (rhs_final_vec(i),i=1,nodes_per_elem)

        write(outfile_unit,fmt='(a,4I6)'),'[A^-1]vector  --> ',n
        do j = 1, nodes_per_elem
           write(outfile_unit,fmt='(12es14.3)') &
               (inverse_A_matrix(j,i),i=1,nodes_per_elem)
        end do

        write(outfile_unit,fmt='(a,4I6)'),'[W_r]vector  --> ',n
        do j = 1, nodes_per_elem
           write(outfile_unit,fmt='(12es14.3)') &
               (matrix_W_right_face(j,i),i=1,nodes_per_elem)
        end do
       
       write(outfile_unit,fmt='(a,4I6)'),'[U]vector  --> ',n
        do j = 1, nodes_per_elem
           write(outfile_unit,fmt='(12es14.3)') &
               (elem_matrix_U(j,i),i=1,nodes_per_elem)
        end do
       write(outfile_unit,fmt='(a,4I6)'),'[U-W]vector  --> ',n
        do j = 1, nodes_per_elem
           write(outfile_unit,fmt='(12es14.3)') &
               (U_minus_W_right_mat(j,i),i=1,nodes_per_elem)
        end do

       write(outfile_unit,fmt='(a,4I6)'),'[U - W]*T^k-1  --> ',n
             write(outfile_unit,fmt='(12es14.3)') &
                 (U_W_times_T(i),i=1,nodes_per_elem)
       
       write(outfile_unit,fmt='(a,4I6)'),'A_inv_U_W_times_T_vec   vector  --> ',n
            write(outfile_unit,fmt='(12es14.3)') &
                (A_inv_U_W_times_T_vec(i),i=1,nodes_per_elem)
       

       write(outfile_unit,fmt='(a,4I6)'),'A_inv_times_q_wl_vec  vector  --> ',n
            write(outfile_unit,fmt='(12es14.3)') &
                (A_inv_times_q_wl_vec(i),i=1,nodes_per_elem)
       
       if( n > 1) then
            write(outfile_unit,fmt='(a,4I6)'),' T_e-1 vector  --> ',n
            write(outfile_unit,fmt='(12es14.3)') &
                  (temperature_soln_new(n-1,i),i=1,nodes_per_elem)             
       else
            write(outfile_unit,fmt='(a,4I6)'),' T_e-1 vector  --> ',n
            write(outfile_unit,fmt='(12es14.3)') &
                  (temperature_soln_prev(n,i),i=1,nodes_per_elem)
        end if
       
       write(outfile_unit,fmt='(a)'), ' ' 
       write(outfile_unit,fmt='(a,4I6)'),' T_e vector  --> ',n
       write(outfile_unit,fmt='(12es14.3)') &
                  (temperature_soln_new(n,i),i=1,nodes_per_elem) 
       
       write(outfile_unit,fmt='(a)'), ' '

       write(outfile_unit,fmt='(a)'),'*****************************************'
    end if
DEBUG = .FALSE.
end subroutine solve_temperature_euler
