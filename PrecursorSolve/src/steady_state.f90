! Steady state solve of the equation 
!
! Input:
!
! Output:
! 
subroutine steady_state()

    USE flags_M
    USE global_parameters_M
    USE mesh_info_M
    USE solution_vectors_M
    USE material_info_M
    
implicit none

!---Dummy

!---Local
    
    integer  :: f,g,i,jj,j,n, nl_iter
    real(dp), dimension(num_isotopes,num_delay_group) :: L2_norm_current,L2_norm_prev 
    real(dp) :: nl_iter_tolerance, difference_L2    
    integer :: difference_counter
    integer :: abs_max_nl_iter
    character(len=20) :: ss_file_name
!---------------------------------------------------------------

!---Set some of the convergance properties
    L2_norm_prev = 0.0
    L2_norm_current = 0.0
    nl_iter_tolerance = 1E-15_dp
    abs_max_nl_iter = 1500_dp

    !---Set starting values for power, velocity, temperature 
    call initialize
            
    !---Steady state solver for nonlinear ss problems
    write(outfile_unit, fmt='(a)'), ' ' 
    write(outfile_unit, fmt='(a)'), 'Start steady state calculation'
    
    nl_iter=1 
    
    nonlinearloop: do !---Nonlinear loop
        elements_loop: do n = 1, num_elem
            !---Computer spatial matrices 
	        call spatial_matrices(n,nl_iter)
            call numerical_flux_matrices(n,nl_iter)
            !---Solve steady state system 
            isotope_loop: do f = 1, num_isotopes
                delay_loop: do g = 1, num_delay_group
                    !---Assemble K, F
                    call assemble_matrix_ss(f,g,n)
                    call solve_precursor_ss(f,g,n,nl_iter)
                enddo delay_loop!---Over delayed groups
            enddo isotope_loop !---Over isotopes
        enddo elements_loop !---Over nodes
        
        !---This counted the convergence of each delayed group
        difference_counter = 0
        
        !---Write out the solution
        if(DEBUG .eqv. .TRUE.) then
            call write_out_soln(outfile_unit,num_elem,transient_save_flag) 
        end if

        !---Swap for for next iteration
        precursor_soln_prev = precursor_soln_new
        
        !---Calculate L2 norm to decide when enough iterations are complete 
        if(nl_iter > 1) then
             do f = 1, num_isotopes
                do g = 1, num_delay_group
                    L2_norm_current(f,g) = sqrt ( sum ( precursor_soln_new(f,g,:,:)*&
                                           precursor_soln_new(f,g,:,:) ))   ! L2 norm
                end do
            end do

        !---Calculate the difference in the L2 norms between iterations 
            do f = 1, num_isotopes
                do g = 1, num_delay_group
                    difference_L2 = abs( L2_norm_prev(f,g) - L2_norm_current(f,g) )
                    if( difference_L2 < nl_iter_tolerance) then
                        difference_counter = difference_counter + 1
                    end if
                end do
            end do
            
            !---Need to make sure the L2 norm converges for all precursor groups
            if ( difference_counter == num_delay_group) then
                max_nl_iter = nl_iter - 1 
            end if

            !---Swap for next iteration
            L2_norm_prev = L2_norm_current
        end if
        
        nl_iter = nl_iter + 1
    
        !---If we've gone thru too many nonlinear iterations exit
        if (nl_iter > max_nl_iter) then
             stop 
         !---If we've gone through a max prescribed and still not converged
        elseif( nl_iter> abs_max_nl_iter) then
             write(outfile_unit,('(a)')) 'We have gone through the max &
                   iterations andstill not converged, something may be wrong'
             stop 
        end if
    
    enddo nonlinearloop !---end nonlinear iteration loop
  
    !---Write out
    write(outfile_unit, fmt=('(a)') ) ' ' 
    write(outfile_unit, fmt=('(a, 100I4)')) 'Number of steady state nonlinear iterations: '&
                                            , nl_iter
    write(outfile_unit, fmt=('(a)') ) ' ' 
    
    !---Make the final converged solution the 'previous' solution
    do f = 1, num_isotopes 
        do g = 1, num_delay_group
            do i = 1, num_elem
                do j = 1, nodes_per_elem
                    precursor_soln_prev(f,g,i,j) = precursor_soln_new(f,g,i,j)
                    precursor_soln_last_time(f,g,i,j) = precursor_soln_new(f,g,i,j)
                end do
            end do
        end do
    end do
   
    !---Set 'previous' to new
    do i = 1, num_elem
        do j = 1,nodes_per_elem
            power_soln_prev(i,j) = power_soln_new(i,j)
            temperature_soln_prev(i,j) = temperature_soln_new(i,j)
            velocity_soln_prev(i,j) = velocity_soln_new(i,j)
        end do
    end do
    !---Set 'previous' to new
    power_amplitude_prev      = power_amplitude_new
    !---For backward euler, need to keep the last value at the time step
    power_amplitude_last_time = power_amplitude_new
    
!---Write precursor solution outfile
    call write_out_soln(outfile_unit,num_elem,transient_save_flag)

!---Set steady state flag off 
    steady_state_flag = .FALSE.

!---Write to outfile     
    write(outfile_unit, '(a)'), ' '
    write(outfile_unit, '(a)'), '*******************************'
    write(outfile_unit, '(a)'), 'End of Steady state solve'
    write(outfile_unit, '(a)'), '*******************************'

end subroutine steady_state
