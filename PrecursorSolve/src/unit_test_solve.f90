! Solves for the current solution vector
! For unit test problem
! 
! Input:
! 
! Output:
! 
subroutine unit_test_solve(nl_iter )

    USE parameters_fe

    implicit none

!   Dummy
    integer :: nl_iter

!   Local

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
   
   !---Solve for next time step solution
   do i = 1, nodes_per_elem
       cur_elem_soln_vec(n,i) = delta_t*elem_prev_soln(i) + previous_elem_soln_vec(n,i) 
   end do

   !---End UNIT TEST if for calc of upwind flux

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
