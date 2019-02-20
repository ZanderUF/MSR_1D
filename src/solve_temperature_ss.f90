!---Solve velocity field
!
!---Input: n - element number 
!
!---Output:
!

subroutine solve_temperature_ss(n,nl_iter)

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
    integer :: i, j, length
    real(dp) :: temperature_eval
    real(dp), dimension(3)   :: rhs_final_vec,elem_vec_w_left,&
                                elem_vec_q_temp
    real(dp), dimension(3,3) :: inverse_matrix
    real(dp) :: volume, heat_capacity_eval, density_eval 
!---Inversion routine parameters
    integer :: lda, info, lwork
    integer,  dimension(3) :: ipiv
    real(dp), dimension(3) :: work
    real(dp) :: salt_mass, q_prime, length_core

!---size of matrix 3x3
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

!---Setup RHS matrix to be inverted later on
    do i = 1, nodes_per_elem
        do j = 1, nodes_per_elem
            inverse_matrix(i,j) = -1.0*elem_matrix_U(i,j) + matrix_W_right_face(i,j) 
        end do
    end do

    !write(outfile_unit,fmt='(a)'), ' '
    !   write(outfile_unit,fmt='(a,I6)'),'[W_r - U] Matrix | element --> ',n
    !   do j=1,nodes_per_elem 
    !         write(outfile_unit,fmt='(12es14.3)') &
    !              ( inverse_matrix(j,i),i=1,nodes_per_elem)             
    !   end do

!---Factorize RHS matrix matrix
    call dgetrf ( length, length, inverse_matrix, lda, ipiv, info )
!---Compute the inverse matrix.
    call dgetri ( length, inverse_matrix, lda, ipiv, work, lwork, info )

    volume = 0.0_dp
!---Calculate core volume
    do i = Fuel_Inlet_start, Fuel_Outlet_End
        do j = 1, nodes_per_elem
            volume = volume + vol_int(j)*area_variation(i,j)
        end do
    end do
   
    temperature_eval = 0.0
    !---Get temperature over element 
    do j = 1, nodes_per_elem 
        !---valuate density based on temperature
        temperature_eval = temperature_eval + vol_int(j)*temperature_soln_new(n,j)
    end do

    call density_corr(temperature_eval, density_eval)

!---Loop over all nodes in element
    do j = 1, nodes_per_elem
        
       call heat_capacity_corr(temperature_eval, heat_capacity_eval)
       length_core = Fuel_Outlet_End - Fuel_Inlet_start
   
       q_prime = (total_power_initial*power_amplitude_prev*&
                  spatial_power_frac_fcn(n,j))/spatial_area_fcn(n,j)
       
       elem_vec_q_temp(j) = elem_vec_q(j)*q_prime*(1.0_dp/(density_eval*heat_capacity_eval))
    end do
!----Evaluate W*T_e-1
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

    do i = 1, nodes_per_elem
        rhs_final_vec(i) = elem_vec_q_temp(i) + elem_vec_w_left(i)
    end do
    

!---Final calculation
    temperature_soln_new(n,:) = matmul(inverse_matrix, rhs_final_vec)

    if( n == Fuel_Inlet_start) then 
        do i = 1, nodes_per_elem
            temperature_soln_new(n,i) = 850.0_dp   
        end do 
    end if


    if( n == Heat_Exchanger_End) then 
        do i = 1, nodes_per_elem
            temperature_soln_new(n,i) = 850.0_dp   
        end do 
    end if

    
    temperature_reduction = 100.0_dp
    do j = 1, nodes_per_elem    
        !---Start heat exchanger
        if ( Heat_Exchanger_Start <= n  .and.  n < Heat_Exchanger_End) then
            temperature_soln_new(n,j) = (temperature_reduction)/&
            (Heat_Exchanger_End - Heat_Exchanger_Start)* &
            (global_coord(Heat_Exchanger_End-1,3)  - global_coord(n,j) ) - temperature_reduction + &
            temperature_soln_new(Heat_Exchanger_Start-1,1)
        end if
    end do 


!---Fix delta T across heat exchanger
    !do j = 1, nodes_per_elem    
    !    !---Start heat exchanger
    !    if ( Heat_Exchanger_Start <= n  .and.  n < Heat_Exchanger_End) then
    !        temperature_soln_new(n,j) = (100.0_dp)/&
    !        (Heat_Exchanger_End - Heat_Exchanger_Start)* &
    !        (global_coord(Heat_Exchanger_End-1,3)  - global_coord(n,3) ) - 100.0_dp + &
    !        temperature_soln_new(Heat_Exchanger_Start-1,1)
    !    end if
    !end do 

if( n > Heat_Exchanger_End) then
        do j = 1, nodes_per_elem
            temperature_soln_new(n,j) = temperature_soln_new(Heat_Exchanger_End,j) 
        end do
    end if

!DEBUG = .TRUE.
if( DEBUG .eqv. .TRUE.) then
      
       write(outfile_unit,fmt='(a,4I6)'), 'nl iter ',nl_iter

       if( n > 1) then
            write(outfile_unit,fmt='(a,4I6)'),' T_e-1 vector  --> ',n
            write(outfile_unit,fmt='(12es14.3)') &
                  (temperature_soln_new(n-1,i),i=1,nodes_per_elem)             
       else
            write(outfile_unit,fmt='(a,4I6)'),' T_e-1 vector  --> ',n
            write(outfile_unit,fmt='(12es14.3)') &
                  (temperature_soln_new(n,i),i=1,nodes_per_elem)
        end if
       
       write(outfile_unit,fmt='(a)'), ' ' 
       write(outfile_unit,fmt='(a,4I6)'),' T_e vector  --> ',n
       write(outfile_unit,fmt='(12es14.3)') &
                  (temperature_soln_new(n,i),i=1,nodes_per_elem) 
       
              write(outfile_unit,fmt='(a)'), ' '
       write(outfile_unit,fmt='(a,I6)'),'[W_r - U]^-1 Matrix | element --> ',n
       do j=1,nodes_per_elem 
             write(outfile_unit,fmt='(12es14.3)') &
                  (inverse_matrix(j,i),i=1,nodes_per_elem)             
       end do
       
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

       write(outfile_unit,fmt='(a)'),'*****************************************'
    end if

DEBUG = .FALSE.

end subroutine solve_temperature_ss
