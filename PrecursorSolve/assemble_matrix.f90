!*****************************************************************************80
!
!! Assemble the elemental to solve for coefficients the problem 
!  Discussion:
!             As the elemental matrices are created we add them to the global 
!
!  Input: n - element number
!
!  Output:
! 

subroutine assemble_matrix (n)
!---
    USE parameters_fe  

    implicit none
!---Dummy variables
    integer :: n
!---Local variables
    integer :: i, j, ii, jj, nr,nc, ncl,length   
    real, dimension(3) :: elem_vec_w_left_face, elem_prev_soln, flux_rhs, flux_lhs, &
                          temp_vec, basis_at_lhs, basis_at_rhs, rhs_final_vec, & 
                          interp_fcn_rhs, interp_fcn_lhs
    real, dimension(3,3) :: inverse_matrix,temp_matrix
    real :: temp_soln
    data interp_fcn_lhs / 1,0,0/
    data interp_fcn_rhs / 0,0,1/

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

    elem_prev_soln = 0
    elem_matrix_G = 0
    elem_vec_w_left_face = 0
    matrix_W_right_face = 0
    matrix_W_left_face = 0
    elem_vec_q = 0
    rhs_final_vec = 0
    
    elem_vec_Pu = 0
    temp_vec = 0
    temp_matrix = 0

!---Steady state case
    if(steady_state_flag .eqv. .TRUE. ) then
         write(outfile_unit,fmt='(a)'), ' '
         write(outfile_unit,fmt='(a)'), 'Start steady state assembly '
         
         !---Applies for all elements except the first one
         if(n > 1) then !--- n - element #
             !---Grab previous precursor conc. + velocity at rhs of previous element
             do i = 1, nodes_per_elem
                do j = 1, nodes_per_elem
                    matrix_W_right_face(i,j) = velocity_vec(n,i)*&
                                               interp_fcn_rhs(i)*interp_fcn_rhs(j)  
                    matrix_W_left_face(i,j) = velocity_vec(n,i)*&
                                              interp_fcn_lhs(i)*interp_fcn_lhs(j)
                    elem_vec_w_left_face(i) =  elem_vec_w_left_face(i) + &
                                matrix_W_left_face(i,j)*cur_elem_soln_vec(n-1,3)
                end do
             end do
         else!---First element case, need to connect with end element 
             do i = 1, nodes_per_elem
                do j = 1, nodes_per_elem
                    matrix_W_right_face(i,j) = velocity_vec(num_elem,i)*&
                                               interp_fcn_rhs(i)*interp_fcn_rhs(j)  
                    matrix_W_left_face(i,j) =  velocity_vec(num_elem,i)*&
                                               interp_fcn_lhs(i)*interp_fcn_lhs(j)
                    elem_vec_w_left_face(i) =  elem_vec_w_left_face(i) + &
                                matrix_W_left_face(i,j)*cur_elem_soln_vec(num_elem,3)
                end do
             end do
         end if!---End flux calculation

         !---Calculate source vector
         do i = 1, nodes_per_elem
            do j = 1, nodes_per_elem
                elem_vec_q(i) = elem_vec_q(i) + elem_matrix_A(i,j)*power_initial(n,i)    
            end do 
            elem_vec_q(i) = elem_vec_q(i)*(beta/gen_time)
         end do

         !---Calculate G matrix, will be inverted later on
         do i = 1, nodes_per_elem
             do j = 1, nodes_per_elem
                 elem_matrix_G(i,j) = (-elem_matrix_U(i,j)) + &
                                       (log(2.0)/lambda)*elem_matrix_A(i,j) + &
                                       matrix_W_right_face(i,j)
             end do 
         end do 
        
         !---RHS vector
         do i = 1, nodes_per_elem
            rhs_final_vec(i) = elem_vec_q(i) + elem_vec_w_left_face(i)
         end do

!****************************************************************
         write(outfile_unit,fmt='(a)'), ' '
         write(outfile_unit,fmt='(a,1I2)'),'G Matrix | element --> ',n
         do i=1,nodes_per_elem 
               write(outfile_unit,fmt='(12es14.3)') &
                    (elem_matrix_G(i,j),j=1,nodes_per_elem)             
         end do
         
         write(outfile_unit,fmt='(a)'), ' '
         write(outfile_unit,fmt='(a,1I2)'),'W right Matrix | element --> ',n
         do j=1,nodes_per_elem 
               write(outfile_unit,fmt='(12es14.3)') &
                    (matrix_W_right_face(j,i),i=1,nodes_per_elem)             
         end do
         
         write(outfile_unit,fmt='(a)'), ' '
         write(outfile_unit,fmt='(a,1I2)'),'W left Matrix | element --> ',n
         do j=1,nodes_per_elem 
               write(outfile_unit,fmt='(12es14.3)') &
                    (matrix_W_left_face(j,i),i=1,nodes_per_elem)             
         end do
         
         write(outfile_unit,fmt='(a)'), ' ' 
         write(outfile_unit,fmt='(a,1I3)'),'left face vector | element --> ', n
         do j=1,nodes_per_elem 
               write(outfile_unit,fmt='(12es14.3)') elem_vec_w_left_face(j)             
         end do
         
         write(outfile_unit,fmt='(a)'), ' ' 
         write(outfile_unit,fmt='(a)'),' '
         write(outfile_unit,fmt='(a,1I2)'),'{q} element source vector | element --> ',n
         write(outfile_unit,fmt='(12es14.3)') (elem_vec_q(i),i=1,nodes_per_elem)             
         write(outfile_unit,fmt='(a)'), ' ' 
         write(outfile_unit,fmt='(a,1I3)'),'RHS final vector | element --> ', n
         do j=1,nodes_per_elem 
               write(outfile_unit,fmt='(12es14.3)') rhs_final_vec(j)             
         end do
    end if !---End steady state

!---Transient case
    !if(unit_test .eqv. .FALSE. .and. steady_state_flag .eqv. .FALSE.) then
    !    write(outfile_unit,fmt='(a)'), 'Start transient calc '
    !    !---Calculate upwind flux
    !    !---Grab current precursor conc. + velocity at rhs of current element
    !    ii = 1 + (n-1)*nodes_per_elem
    !    temp_soln = velocity_vec(ii+2)*previous_elem_soln_vec(ii+2)
    !    do i = 1, nodes_per_elem
    !        flux_rhs(i) = temp_soln*basis_at_rhs(i)
    !    end do
    !    !---All elements except the first
    !    if(n > 1) then
    !        temp_soln = velocity_vec(ii-1)*previous_elem_soln_vec(ii-1)
    !        !---Grab previous precursor conc. + velocity at rhs of previous element
    !        do i = 1, nodes_per_elem
    !            flux_lhs(i) = temp_soln*basis_at_rhs(i)
    !        end do
    !    else!---First element, need to connect with last
    !        temp_soln = velocity_vec(matrix_length)*previous_elem_soln_vec(matrix_length)
    !        do i = 1, nodes_per_elem
    !            flux_lhs(i) = temp_soln*basis_at_lhs(i)
    !        end do
    !    end if!---End flux calculation 

    !    !---Calculate (P + lambda_i*A)*C_h
    !    temp_matrix = elem_matrix_U + (lambda*elem_matrix_A)
    !    ii = (n-1)*nodes_per_elem
    !    do i = 1, nodes_per_elem
    !        do j = 1, nodes_per_elem
    !            temp_vec(i) = temp_matrix(i,j)*previous_elem_soln_vec(ii+i)
    !        end do
    !    end do

    !    !---Prep A matrix for inversion
    !    inverse_matrix = elem_matrix_A
    !     
    !    !---Calculate RHS vector (change in the solution)
    !    rhs_final_vec = elem_vec_q + flux_rhs - flux_lhs + temp_vec 

    !end if!---END transient case

!---UNIT TEST
   ! if(unit_test .eqv. .TRUE.) then
   ! !---Grab current solution at rhs of current element 
   !     ii = 1 + (n-1)*nodes_per_elem
   !     temp_soln = previous_elem_soln_vec(ii+2)
   !     do i = 1, nodes_per_elem
   !         flux_rhs(i) = temp_soln*basis_at_rhs(i)
   !     end do
   ! !---All elements except the first
   !     if(n > 1) then
   !         !---Grab previous solution at rhs of previous element
   !         temp_soln = previous_elem_soln_vec(ii-1)
   !         do i = 1, nodes_per_elem
   !             flux_lhs(i) = temp_soln*basis_at_lhs(i)
   !         end do
   !     else!---For periodic B.C. need to connect first to the last element
   !         temp_soln = previous_elem_soln_vec(matrix_length)
   !         do i = 1, nodes_per_elem
   !             flux_lhs(i) = temp_soln*basis_at_lhs(i)
   !         end do
   !     end if
   ! !---Solve for elemental coeficients, no global assembly 
   !     do i = 1, nodes_per_elem
   !         do j = 1, nodes_per_elem
   !             ii = j + (n-1)*nodes_per_elem
   !             elem_vec_Pu(i) = elem_vec_Pu(i) + &
   !                              elem_matrix_U(i,j)*previous_elem_soln_vec(ii)
   !         end do
   !     end do
   ! !---Combine (Pu - f)
   !     rhs_final_vec = elem_vec_Pu - (flux_rhs - flux_lhs)
   !     inverse_matrix = elem_matrix_A
   ! end if 
   ! !---End UNIT TEST if for calc of upwind flux
    
    do i = 1, nodes_per_elem
        do j = 1, nodes_per_elem
            inverse_matrix(i,j) = elem_matrix_G(i,j)  
        end do
    end do 
    !---Factorize matrix
    call dgetrf ( length, length, inverse_matrix, lda, ipiv, info )
    !---Compute the inverse matrix.
    call dgetri ( length, inverse_matrix, lda, ipiv, work, lwork, info )      
    
!---Write out inverse matrix
    write(outfile_unit,fmt='(a)'), ' ' 
    write(outfile_unit,fmt='(a,1I3)'),'inverse [G]^-1 matrix | element --> ',n
    do j=1,nodes_per_elem
           write(outfile_unit,fmt='(12es14.3)') &
                ( inverse_matrix(j,i) ,i=1,nodes_per_elem )             
    end do
    
!---CALCULATE SOLUTION for a given element
    
    do i = 1, nodes_per_elem
        do j = 1, nodes_per_elem
            elem_prev_soln(i) = elem_prev_soln(i) + inverse_matrix(i,j)*rhs_final_vec(j)
        end do
    end do

!---Evaluate time step if transient case
    if(steady_state_flag .eqv. .FALSE.) then    
        !---Solve for next time step solution
        do i = 1, nodes_per_elem
            cur_elem_soln_vec(n,i) = dt*elem_prev_soln(i) + previous_elem_soln_vec(n,i) 
        end do
    else !---Steady state
        do i = 1, nodes_per_elem
            cur_elem_soln_vec(n,i) = elem_prev_soln(i)
        end do
    end if

    write(outfile_unit,fmt='(a)'), ' ' 
    write(outfile_unit,fmt='(a,1I3)'),'Solution | element --> ', n
    do j=1,nodes_per_elem 
          ii = j + (n-1)*nodes_per_elem
          write(outfile_unit,fmt='(a,1I2,12es14.3)'), 'Node -->', ii, elem_prev_soln(j)          
    end do   
    write(outfile_unit,fmt='(a)'), '********************************'

!---UNIT TEST write info out to file
    if(unit_test .eqv. .TRUE.) then
    !---Write out {flux_rhs}
        write(outfile_unit,fmt='(a)'), ' ' 
        write(outfile_unit,fmt='(a,1I3)'),'{f^+} vector | element --> ', n
        do j=1,nodes_per_elem 
               write(outfile_unit,fmt='(12es14.3)') flux_rhs(j)             
        end do
    
    !---Write out {flux_lhs}
        write(outfile_unit,fmt='(a)'), ' ' 
        write(outfile_unit,fmt='(a,1I3)'),'{f^-} vector | element --> ', n
        do j=1,nodes_per_elem 
               write(outfile_unit,fmt='(12es14.3)') flux_lhs(j)             
        end do
    
    !---Write out {Pu}
        write(outfile_unit,fmt='(a)'), ' ' 
        write(outfile_unit,fmt='(a,1I3)'),'{Pu} vector | element --> ', n
        do j=1,nodes_per_elem 
               write(outfile_unit,fmt='(12es14.3)') elem_vec_Pu(j)             
        end do
    !---Write out {Pu - f}
        write(outfile_unit,fmt='(a)'), ' ' 
        write(outfile_unit,fmt='(a,1I3)'),'{Pu - f} vector | element --> ', n
        do j=1,nodes_per_elem 
               write(outfile_unit,fmt='(12es14.3)') rhs_final_vec(j)             
        end do   
    end if

end 
