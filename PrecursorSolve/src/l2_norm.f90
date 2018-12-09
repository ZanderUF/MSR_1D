

!---Calculate the l2 norm of the precursor soln

subroutine l2_norm(nl_iter,difference_counter,L2_norm_prev,L2_norm_current)

    USE flags_M
    USE global_parameters_M
    USE solution_vectors_M
    USE time_info_M
    USE mesh_info_M
    USE material_info_M
    
    
    implicit none
   
    !---Dummy
    integer, intent(in) :: nl_iter
    integer, intent(inout) :: difference_counter
    real(dp), dimension(num_isotopes, num_delay_group), intent(inout)::&
                    L2_norm_prev,L2_norm_current 
    
    !---Local
    integer  :: f,g,n,i,j    
    real(dp) :: nl_iter_tolerance, difference_L2
    integer  :: abs_max_nl_iter
    

    max_nl_iter = 30 
    abs_max_nl_iter = 600
    nl_iter_tolerance = 1E-12_dp
    
    !---Calculate L2 norm of precursor solution
    if(nl_iter > 1) then
        do f = 1, num_isotopes
            do g = 1, num_delay_group
                L2_norm_current(f,g) = sqrt(sum( precursor_soln_new(f,g,:,:)*&
                                       precursor_soln_new(f,g,:,:) ))   
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

    end if 

end subroutine l2_norm
