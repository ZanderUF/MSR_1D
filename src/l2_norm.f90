

!---Calculate the l2 norm of the precursor soln

subroutine l2_norm(nl_iter)

    USE flags_M
    USE global_parameters_M
    USE solution_vectors_M
    USE time_info_M
    USE mesh_info_M
    USE material_info_M
    
    use Mod_SetupOutputFiles
    
    implicit none
   
    !---Dummy
    integer, intent(in)    :: nl_iter
    
    !---Local
    integer  :: f,g,n,i,j    
    real(dp) :: nl_iter_tolerance, difference_L2
    integer  :: abs_max_nl_iter
    integer  :: difference_counter
   
    !---Describe absolute default tolerances
    max_nl_iter = 500 
    abs_max_nl_iter = 600
    nl_iter_tolerance = 1E-4_dp
  
    !---Reset the difference counters
    L2_diffs_precursors(:,:) = 0.0
    L2_diffs_temperature     = 0.0
    L2_current_precursors    = 0.0 
    
    !---Calculate L2 norm of precursor solution
    if(nl_iter > 1) then
        do f = 1, num_isotopes
            do g = 1, num_delay_group
                do i = 1, num_elem
                    do j = 1, nodes_per_elem
                        !---Calculate current sum of the squares 
                        L2_current_precursors(f,g) = L2_current_precursors(f,g) + &
                            precursor_soln_new(f,g,i,j)*precursor_soln_new(f,g,i,j)
                        
                        !---Calculate SUM([Current - Previous]^2)
                        L2_diffs_precursors(f,g) = L2_diffs_precursors(f,g) + &
                            (precursor_soln_new(f,g,i,j) - precursor_soln_prev(f,g,i,j) )*&
                            (precursor_soln_new(f,g,i,j) - precursor_soln_prev(f,g,i,j) )
                    end do
                end do
                                
                !---Get the current L2
                L2_current_precursors(f,g) = sqrt(L2_current_precursors(f,g))
                !---Get sqrt of the sum of the differences
                L2_diffs_precursors(f,g)   = sqrt(L2_diffs_precursors(f,g))
            end do
        end do
        
        !---Calculate L2 norm for the temperature
        do i = 1, num_elem
            do j = 1, nodes_per_elem
                !--Calculate current sum of squares
                L2_current_temperature = L2_current_temperature + &
                         temperature_soln_new(i,j)*temperature_soln_new(i,j)
                
                !---SUM of the square differences
                L2_diffs_temperature = L2_diffs_temperature + &
                        (temperature_soln_new(i,j) - temperature_soln_prev(i,j))*&
                        (temperature_soln_new(i,j) - temperature_soln_prev(i,j))   
            end do
        end do

        !---Check if temperature has converged
        if(L2_diffs_temperature < nl_iter_tolerance*L2_current_temperature) then
            temperature_converged = .TRUE.
        end if

        !---Reset the difference counter
        difference_counter = 0
        
        !---Check if the convergence criteria has been reached for all groups
        do f = 1, num_isotopes
           do g = 1, num_delay_group
                !---Evaluate ||Current - Previous|| < tolerance * ||Current||
                if(L2_diffs_precursors(f,g) < nl_iter_tolerance*&
                   L2_current_precursors(f,g)) then
                    
                    difference_counter = difference_counter + 1
                
                end if
           end do
        end do

        !---Need to make sure the L2 norm converges for all precursor groups
        if ( difference_counter == num_delay_group) then
            max_nl_iter = nl_iter - 1
        end if
        
    !---Swap for next iteration
        L2_prev_precursors  = L2_current_precursors
        L2_prev_temperature = L2_current_temperature
    end if

    !---Check if too many nonlinear iterations and not converging
    if ( nl_iter > max_nl_iter) then
        if(nl_iter > abs_max_nl_iter) then
            write(outfile_unit,fmt=('(a)')) 'Gone past max amount of &
                         nonlinear iterations &
                         and might have a problem'
        else
            if (DEBUG .eqv. .TRUE.) then 
                write(outfile_unit,fmt=('(a,I3,a,8es14.3)')) &
                'Took this # of iterations to converge --> ',nl_iter,&
                '  <-- at time step -->', t0
            end if
        end if
        
        nl_iter_flag = .FALSE.

    !---Write to outfile
        if(DEBUG .eqv. .TRUE.) then
            do f = 1, num_isotopes
                write(nl_outfile_unit,fmt=('(I6,8es14.5,8es14.5,8es14.5,8es14.5,&
                                     8es14.5,8es14.5,8es14.5)'))&
                nl_iter,L2_current_precursors(1,1),L2_diffs_temperature,(L2_diffs_precursors(f,g), g=1,num_delay_group)
            end do
        end if
    end if 

end subroutine l2_norm
